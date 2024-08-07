#!/usr/bin/env python3

##
# @file postprocess.py
# @author Viktar Kireyeu
# @version 1.0
# @date    2024
# @copyright Do not modify a single character without the author's permission.
#
# @brief  Postprocessing macro for the MpdRoot "nuclei" wagon.
#
#
# New subdirectories will be created to store the produced histograms (PDF):
# - plots/efficiency/tpc
# - plots/efficiency/tof
# - plots/efficiency/pid
# - plots/efficiency/pid_dedx
# - plots/contamination/secondaries
# - plots/contamination/tof
# - plots/contamination/pid
# - plots/contamination/pid_dedx
# - plots/results/pt
# - plots/results/dndy
# - plots/results/coal
# 
# If the PTCORR file is set, then instead of the post-analysis the new
# ROOT-file with the \f$ p_{T} \f$ corrections 1d and 2d profiles will be produced.
#
# Usage: postprocess [-h] [-i INPUT] [-e EFFICIENCIES] [-s SETTINGS] [-o OUTPUT]
#                    [-d DIR] [-r REPORT] [-p PTCORRECTIONS] [-t TABLE] [--qa]
#                    [--png] [--dedx] [--skip] [--clear]
# 
# MpdNuclei wagon data processing program
# 
# Options:
# <pre>
#   -h, --help            show this help message and exit
#   -i INPUT, --input INPUT
#                         Input ROOT file
#   -e EFFICIENCIES, --eff EFFICIENCIES
#                         Efficiencies ROOT file
#   -s SETTINGS, --settings SETTINGS
#                         Settings JSON file
#   -o OUTPUT, --output OUTPUT
#                         Output spectra ROOT-file
#   -d DIR, --dir DIR     Output directory for histograms
#   -r REPORT, --report REPORT
#                         Output report file
#   -p PTCORRECTIONS, --ptcorr PTCORRECTIONS
#                         Make file for pt corrections
#   -t TABLE, --table TABLE
#                         Output file with dndy info for fits etc (JSON)
#   --qa                  Produce extended QA plots
#   --png                 Produce PNG-plots onstead of PDF
#   --dedx                Produce histograms for the case PID = dE/dx only
#   --skip                Skip already triggered functions
#   --clear               Clear the output directory before populating
# 
# If the EFFICIENCIES file is not set, then the INPUT will be used instead.
# </pre>

import sys
try:
    import os
    import json
    import math
    import argparse
    import fileinput
    from array import array
    from pathlib import Path
    import ROOT

    from ROOT import TCanvas, TFile, TLatex, TLegend, TPad, TRatioPlot, TGraph, TGraphErrors, TF1, TMath, TObject, TPaveText
    from ROOT import gROOT, gSystem, gPad

    from pylatex import Document, Section, Subsection, Math, NoEscape, Package, Command
    import pylatex.config as cf

except ModuleNotFoundError as err:
    sys.exit(err)

parser = argparse.ArgumentParser(
                    prog        = 'postprocess',
                    description = 'MpdNuclei wagon data processing program',
                    epilog      = 'If the EFFICIENCIES file is not set, then the INPUT will be used instead.')
parser.add_argument('-i', '--input',    metavar = 'INPUT',         help = 'Input ROOT file',                 default = 'taskNuclei.root')
parser.add_argument('-e', '--eff',      metavar = 'EFFICIENCIES',  help = 'Efficiencies ROOT file')
parser.add_argument('-s', '--settings', metavar = 'SETTINGS',      help = 'Settings JSON file',              default = 'NucleiAna.json')
parser.add_argument('-o', '--output',   metavar = 'OUTPUT',        help = 'Output spectra ROOT-file',        default = 'nuclei_ptspectra.root')
parser.add_argument('-d', '--dir',      metavar = 'DIR',           help = 'Output directory for histograms', default = 'pyplots')
parser.add_argument('-r', '--report',   metavar = 'REPORT',        help = 'Output report file')
parser.add_argument('-p', '--ptcorr',   metavar = 'PTCORRECTIONS', help = 'Make file for pt corrections')
parser.add_argument('-t', '--table',    metavar = 'TABLE',         help = 'Output file with dndy info for fits etc (JSON)')
parser.add_argument('--qa',             help = 'Produce extended QA plots',                        action = 'store_true')
parser.add_argument('--png',            help = 'Produce PNG-plots onstead of PDF',                 action = 'store_true')
parser.add_argument('--dedx',           help = 'Produce histograms for the case PID = dE/dx only', action = 'store_true')
parser.add_argument('--skip',           help = 'Skip already triggered functions',                 action = 'store_true')
parser.add_argument('--clear',          help = 'Clear the output directory before populating',     action = 'store_true')


rapidity_bins = [[-1.0, -0.5], [-0.5, 0.5], [0.5, 1.0]]
# ~ rapidity_bins = [[-1.0, 1.0]]
SMALL_DELTA = 0.00001
LW       =  5 # Line width
MS_MC    = 20 # Marker style: full circles
MS_UNCOR = 24 # Marker style: open circles
MS_COR   = 22 # Marker style: full triangles
LC_MC    =  1 # Line color: black
LC_UNCOR =  2 # Line color: red
LC_COR   =  3 # Line color: green

STATUS_FILE = 'postprocess.status'



# This is the main function, all other subroutines are called from here
def main():
    global EXT          # File extension for the produced plots: 'pdf', 'png', maybe something else
    global Particles    # Global list with the particles definitions from the settings file
    global Centrality   # Global list with the centralities definitions from the settings file
    global SYS          # Global variable with the colliding system title
    global DEDX         # Global switch for the case 'PID = dE/dx only': 
                        #   1 = make such plots; 0 = make plots for the combined PID case
    global CURRENT_ST   # Global list of the current functions calls counters
    global SKIP         # Global switch for skipping already triggered (during previous runs) functions
    global INP_FILE     # Input ROOT file
    global EFF_FILE     # Efficiencies ROOT file
    global OUT_FILE     # Output spectra ROOT-file
    global TBL_FILE     # Output file (JSON) to dump the BW fits and other information

    CURRENT_ST = [0, 0, 0, 0, 0, 0, 0, 0]
    SYS = 'Bi+Bi, #sqrt{s_{NN}} = 9.2 GeV'

    # Command line arguments parsing
    args = parser.parse_args()
    input_file      = args.input           # Input file
    efficiency_file = args.eff             # Efficiencies file
    settings_file   = args.settings        # Settings file
    spectra_file    = args.output          # Output file for the spectra
    outdir          = args.dir             # Output directory to store plots
    extended        = args.qa              # Switch for the additional efficiencies plots
    report_file     = args.report          # Switch for the report file generation
    TBL_FILE        = args.table           # Output file for the BW fits info
    if args.eff is None:                   # If efficiencies file is not set
        efficiency_file = input_file       # we can use the input file

    EXT   = 'png' if args.png   else 'pdf' # pdf/png switch
    DEDX  = True  if args.dedx  else False # Switch for the 'PID = dE/dx only' case
    SKIP  = True  if args.skip  else False # Switch to skip previously triggered functions

    make_ptcorr = False
    if args.ptcorr is not None:            # If PTCOR file is set
        corrections_file = args.ptcorr     # then the program will only produce
        make_ptcorr = True                 # this file without furter analysis

    print(f"Input file: {input_file}")
    check_file_exist(input_file)                  # Check if the input file exists
    print(f"Settings file: {settings_file}")
    check_file_exist(settings_file)               # Check if the settings file exists

    if(make_ptcorr):
        print("*** pT corrections mode ***")
        print(f"pT corrections file: {corrections_file}")
    else:
        print(f"Efficiencies file: {efficiency_file}")
        check_file_exist(efficiency_file)         # Check if the efficiencies file exists
        print(f"Spectra file: {spectra_file}")
        if report_file is not None:
            print(f"Report: {report_file}")
        if TBL_FILE is not None:
            print(f"JSON table: {TBL_FILE}")
            if os.path.exists(TBL_FILE):
                os.system(f'rm -f {TBL_FILE}')    # Remove the JSON 'table' file if it exists
        if args.clear:                            # The 'dry' run case -- clear everything before proceed
            os.system(f'rm -rf {outdir}')         # Clear remnants from previous runs
            os.system(f'rm -vf {STATUS_FILE}')
            os.system(f'rm -vf {spectra_file}')
            os.system(f'mkdir -pv {outdir}/results/{{pt,dndy,coal}}')
        if extended:                              # The extended output case -- create the needed sub-directories
            os.system(f'mkdir -pv {outdir}/efficiency/{{tpc,tof,pid}}')
            os.system(f'mkdir -pv {outdir}/contamination/{{secondaries,pid,tof}}')
            if DEDX:                              # Same for the 'PID = dE/dx only' case
                os.system(f'mkdir -pv {outdir}/efficiency/pid_dedx')
                os.system(f'mkdir -pv {outdir}/contamination/pid_dedx')

    ROOT.gROOT.SetBatch(True)                     # Batch mode for ROOT
    ROOT.gStyle.SetErrorX(0)                      # Do not draw error along the X-axis
    ROOT.gStyle.SetOptStat(0)                     # Do not draw statistics box on the plots

    read_status('main')                           # Load file with the function trigger counters
    if os.path.exists(settings_file):             # Check and open the settings file
        config = open(settings_file)
    else:
        print(f'File not exist: {settings_file}') # Exit if the settings file does not exist
        exit()
    data = json.load(config)                      # Load data from the settings file

    Centrality = data['Events']['Centrality']     # Centrality bins
    n_centrality_bins = len(Centrality)           # Number of centrality bins
    Particles = data['Particles']                 # Particles information (PDG, mass, cuts -- everything)
    n_particles = len(Particles)                  # Number of particles

    INP_FILE = ROOT.TFile.Open(input_file,"READ") # Open and check (if not 'zombie') the input file
    check_root_file(INP_FILE)

    if(make_ptcorr):                              # Create file for the pT and pZ corrections
        corr_file = ROOT.TFile.Open(corrections_file,"RECREATE")
        check_root_file(corr_file)
        corr_file.cd()
    else:                                         # Open and check the efficiencies file
        EFF_FILE = ROOT.TFile.Open(efficiency_file,"READ")
        check_root_file(EFF_FILE)
        OUT_FILE = ROOT.TFile.Open(spectra_file, "UPDATE") # Create and check the output file
        check_root_file(OUT_FILE)
        OUT_FILE.cd()

    c1 = TCanvas('c_landscape', '', 800, 600)     # 'Landscape' canvas
    c1.SetLeftMargin(0.14)                        # All four margins: left, right, bottom and top
    c1.SetRightMargin(0.1)
    c1.SetBottomMargin(0.14)
    c1.SetTopMargin(0.1)
    c2 = TCanvas('c_portrait', '', 600, 800)      # 'Portrait' canvas

    for particle in Particles:
        if(not make_ptcorr):
            os.system(f'mkdir -pv {outdir}/results/{{pt,dndy}}/{particle}') # Create the separate sub-directory
            os.system(f'mkdir -pv {outdir}/results/coal')                   # for the each particle in the file-system
            OUT_FILE.mkdir(f'{particle}', '', returnExistingDirectory=True) # and in the output file
            OUT_FILE.cd()
        for cbin in range(0, n_centrality_bins):
            common_suffix = f"{particle}_centrality{cbin}"
            if(make_ptcorr): # Copy pT and pZ TProfiles to the corrections file
                pz_corr_2d = INP_FILE.Get(f"pz__corr_{common_suffix}").Clone(f"p2d__pzcorr_{common_suffix}")
                pz_corr_2d.Write()
                pt_corr_2d = INP_FILE.Get(f"pt__corr_{common_suffix}").Clone(f"p2d__ptcorr_{common_suffix}")
                pt_corr_2d.Write()
            else:            # Make and analysis
                # MC spectra
                make_pt_spectra(0, 0, particle, cbin, f"h__pty_mc_{common_suffix}", "mc")
                # Uncorrected reconstructed spectra with combined PID = dE/dx and m^2
                make_pt_spectra(0, 0, particle, cbin, f"h__pty_pid_{common_suffix}", "uncorr")
                # Corrected reconstructed spectra with combined PID = dE/dx and m^2
                make_pt_spectra(1, 0, particle, cbin, f"h__pty_pid_{common_suffix}", "corr")
                # Combined efficiency: TPC + ToF + PID
                make_overall_efficiency(0, particle, cbin)
                if DEDX:
                    # Uncorrected reconstructed spectra with PID = dE/dx only
                    make_pt_spectra(0, 1, particle, cbin, f"h__pty_pid_dedx_{common_suffix}", "uncorr")
                    # Corrected reconstructed spectra with PID = dE/dx only
                    make_pt_spectra(1, 1, particle, cbin, f"h__pty_pid_dedx_{common_suffix}", "corr")                    

    INP_FILE.Close() # At this point the input file is not needed anymore
    if(make_ptcorr):
        corr_file.Close() # Same for the pT and pZ corrections case 
    else:
        EFF_FILE.Close()  # Efficiencies were already applied, efficiencies file can be closed
        for cbin in range(0,n_centrality_bins):
            cname = f"{Centrality[cbin][0]} - {Centrality[cbin][1]}%"
            for p_index, particle in enumerate(Particles):
                common_suffix = f"{particle}_centrality{cbin}"
                system = f"{SYS}, {particle}, {cname}"
                # Draw results: pT and rapidity spectra, efficiencies
                draw_pt_spectra_single(c2, 0, p_index, cbin, f"{outdir}/results/pt/{particle}/{particle}_centrality{cbin}")
                draw_overall_efficiency(c1, p_index, cbin, f"{outdir}/results/pt/{particle}/{common_suffix}_overall_efficiency")
                if particle != "ap" and particle != "t" and particle != "He4":
                    draw_dndy(0, c1, p_index, cbin, f"{outdir}/results/dndy/{particle}/{common_suffix}.{EXT}")
                if DEDX:
                    draw_pt_spectra_single(c2, 1, p_index, cbin, f"{outdir}/results/pt/{particle}/{common_suffix}_dedx")
                if(extended):
                    # TPC and primaries efficiency
                    draw_efficiency(c1, f"{particle}/h__efficiency_tpc_{common_suffix}",
                                        f"TPC efficiency, {system}",
                                        f"{outdir}/efficiency/tpc/{common_suffix}.{EXT}")
                    # Secondaries contamination
                    draw_efficiency(c1, f"{particle}/h__contamination_secondaries_{common_suffix}",
                                        f"Secondaries contamination, {system}",
                                        f"{outdir}/contamination/secondaries/{common_suffix}.{EXT}")
                    # ToF efficiency
                    draw_efficiency(c1, f"{particle}/h__efficiency_tof_{common_suffix}",
                                        f"ToF efficiency, {system}",
                                        f"{outdir}/efficiency/tof/{common_suffix}.{EXT}")
                    # ToF wrong flag (contamination)
                    draw_efficiency(c1, f"{particle}/h__contamination_tof_{common_suffix}",
                                        f"ToF wrong flag contamination, {system}",
                                        f"{outdir}/contamination/tof/{common_suffix}.{EXT}")
                    # PID combined efficiency
                    draw_efficiency(c1, f"{particle}/h__efficiency_pid_{common_suffix}",
                                        f"PID efficiency, {system}",
                                        f"{outdir}/efficiency/pid/{common_suffix}.{EXT}")
                    # PID combined contamination
                    draw_efficiency(c1, f"{particle}/h__contamination_pid_{common_suffix}",
                                        f"PID contamination, {system}",
                                        f"{outdir}/contamination/pid/{common_suffix}.{EXT}")
                    if DEDX:
                        # PID dE/dx efficiency
                        draw_efficiency(c1, f"{particle}/h__efficiency_pid_dedx_{common_suffix}",
                                            f"PID (dE/dx only) efficiency, {system}",
                                            f"{outdir}/efficiency/pid_dedx/{common_suffix}.{EXT}")
                        # PID dE/dx contamination
                        draw_efficiency(c1, f"{particle}/h__contamination_pid_dedx_{common_suffix}",
                                            f"PID (dE/dx only) contamination, {system}",
                                            f"{outdir}/contamination/pid_dedx/{common_suffix}.{EXT}")
            draw_b2(c1, cbin, f'{outdir}/results/coal/b2_centrality{cbin}')
        OUT_FILE.Close()
    config.close()
    if not make_ptcorr and report_file is not None: # Produce the report file (both tex and pdf)
        generate_report(extended, outdir, Path(report_file).stem)



##  This subroutine draws the efficiency or contamination histograms
#  \param canvas TCanvas prepared for drawing
#  \param psname The name of the selected efficiency (or contamination) histogram in the input file
#  \param title The title for the output histogram
#  \param fname The output pdf-file name
#
def draw_efficiency(canvas, psname, title, fname):
    order = 3                                                 # These three lines are common for several functions
    CURRENT_ST[order] += 1                                    # Increment the function calls counts
    function_name = sys._getframe().f_code.co_name            # Here we get the function name, the execution order and a counter
    if skip_function(function_name, order): return            # Thus we can skip this function execution if it was already
                                                              # called during the previous program run
    canvas.cd()                                               # Go to this canvas
    hPhaseSpace = OUT_FILE.Get(psname).Clone("hPhaseSpace")   # Get the needed phase space to draw
    frame = canvas.DrawFrame(-2.99, 0, 2.99, 4.99);           # Draw frame with chosen X/Y axes limits
    make_fancy_frame(frame, title, "y", "p_{T}, GeV/c");      # Set frame title, axes labels and fonts
    hPhaseSpace.Draw('same colz')                             # Draw phase-space inside the frame
    canvas.Print(fname)                                       # The canvas content is printed to the file
    write_status(f'{function_name} {CURRENT_ST[order]}')      # Here the function execution counter is updated in the status file



##  This subroutine is a thermal fit for the \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param x Transverse momentum \f$ p_{T} \f$
#  \param par The fit parameters set:
#     - par[0] dN/dy
#     - par[1] Temperature T
#     - par[2] Particle mass (must be fixed)
#
def fit_thermal(x, par):
    dndy = par[0]
    T    = par[1]
    m0   = par[2]
    mt   = math.sqrt(x[0]*x[0] + m0*m0)
    val  = (dndy/(T*(m0 + T))) * math.exp(-(mt - m0)/T)
    return val



##  This subroutine is a Blast-Wave fit for the \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param x Transverse momentum \f$ p_{T} \f$
#  \param par The fit parameters set:
#     - par[0] Blastwave normalization constant
#     - par[1] Mass of the particle (must be fixed)
#     - par[2] Blastwave temperature
#     - par[3] Blastwave flow velocity
#
def fit_blastwave_pt_invariant(x, par):
    pt    = x[0]
    C     = par[0]
    m0    = par[1]
    T     = par[2]
    beta  = par[3]

    mt    = ROOT.TMath.Sqrt(pt*pt + m0*m0)
    R_max = 20
    steps = 40
    dr    = R_max / steps

    fitval = 0
    for i in range(steps):
        r      = i * dr
        betar  = beta * ROOT.TMath.Power(r/R_max, 1.0)
        rho    = TMath.ATanH(betar)
        b1     = ROOT.TMath.BesselI0( pt * ROOT.TMath.SinH(rho) / T )
        b2     = ROOT.TMath.BesselK1( mt * ROOT.TMath.CosH(rho) / T )
        fitval += r * mt * dr * b1 * b2
    return C * fitval



##  This subroutine is a Blast-Wave fit for the \f$ d^{2}N/dp_{T}dy \f$ spectra
#  \param x Transverse momentum $p_{T}$
#  \param par The fit parameters set:
#     - par[0] Blastwave normalization constant
#     - par[1] Mass of the particle (must be fixed)
#     - par[2] Blastwave temperature
#     - par[3] Blastwave flow velocity
#
def fit_blastwave_pt(x, par):
    pt    = x[0]
    C     = par[0]
    m0    = par[1]
    T     = par[2]
    beta  = par[3]

    mt    = ROOT.TMath.Sqrt(pt*pt + m0*m0)
    R_max = 20
    steps = 40
    dr    = R_max / steps

    fitval = 0
    for i in range(steps):
        r      = i * dr
        betar  = beta * ROOT.TMath.Power(r/R_max, 1.0)
        rho    = TMath.ATanH(betar)
        b1     = ROOT.TMath.BesselI0( pt * ROOT.TMath.SinH(rho) / T )
        b2     = ROOT.TMath.BesselK1( mt * ROOT.TMath.CosH(rho) / T )
        fitval += r * mt * dr * b1 * b2
    return C * pt * fitval



##  This subroutine is a double-exponential (thermal) fit for the \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param x Transverse momentum $p_{T}$
#  \param par The fit parameters set:
#     - par[0] Constant 1
#     - par[1] Constant 2
#     - par[2] Temperature T 1
#     - par[3] Temperature T 2
#     - par[4] Particle mass (must be fixed)
#
def fit_2exp(x, par):
    pt  = x[0]   # Transverse momentum
    m0  = par[4] # Particle mass (fixed parameter)
    mtm = ROOT.TMath.Sqrt(pt * pt + m0 * m0) - m0
    x1 = (1. + par[1])/(par[2] * (m0 + par[2])) / ROOT.TMath.Exp(mtm / par[2]);
    x2 = par[1] / (par[3] * (m0 + par[3])) / ROOT.TMath.Exp(mtm / par[3]);
    val = par[0] * (x1 - x2)
    return val


##  This subroutine is a Bose-Einstein fit for the \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param x Transverse momentum \f$ p_{T} \f$
#  \param par The fit parameters set:
#     - par[0] dN/dy
#     - par[1] Temperature T
#     - par[2] Particle mass (must be fixed)
#
def fit_bose_einstein(x, par):
    pt   = x[0]
    dndy = par[0]
    T    = par[1]
    m0   = par[2]
    mt   = math.sqrt(x[0]*x[0] + m0*m0)
    val  = dndy / (math.exp(mt/T) - 1)
    return val

##  This subroutine draws the \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra for the defined rapidity bins
#  \param canvas TCanvas prepared for drawing
#  \param do_dedx The switch for the PID methods: 
#     - 1 = make plots for the dE/dx only PID 
#     - 0 = make plots for the combined PID 
#  \param p_pos The particle index
#  \param cbin  The centrality bin index
#  \param fname The output pdf-file name
#
def draw_pt_spectra_single(canvas, do_dedx, p_pos, cbin, fname):
    order = 1
    function_name = sys._getframe().f_code.co_name
    CURRENT_ST[order] += 1
    if skip_function(function_name, order): return

    pname    = list(Particles)[p_pos]                         # Particle name
    cbins    = list(Centrality)[cbin]                         # Centrality bin limits
    cname    = f'{cbins[0]} - {cbins[1]}%'                    # Centrality bin name
    pt_start = list(Particles.values())[p_pos]["pt_bins"][1]  # pT limits: low and 
    pt_end   = list(Particles.values())[p_pos]["pt_bins"][2]  # high
    p_mass   = float(list(Particles.values())[p_pos]["Mass"]) # Particle mass
    print(f' draw pT =============== {pname},  {cname},  {pt_start} < pT < {pt_end}  =============== ')
    
    if(do_dedx):
        postfix='_dedx'                                       # Histograms name postfix for the 'PID = dE/dx only'
    else:
        postfix=''
    canvas.cd()
    vMC = []                                                  # The list of Monte-Carlo histograms
    vUncorr = []                                              # The list of uncorrected histograms
    vCorr = []                                                # The list of corrected histograms
    common_suffix = f'{pname}_centrality{cbin}'               # Histograms name suffix for each particle and centrality
    for rbin in range(len(rapidity_bins)):                    # Populate lists with histograms for the eah rapidity bin
        common_postfix = f'{postfix}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'
        hMC     = OUT_FILE.Get(f"{pname}/h__pt_{common_suffix}_mc_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}").Clone("hMC")
        hUncorr = OUT_FILE.Get(f"{pname}/h__pt_{common_suffix}_uncorr{common_postfix}").Clone("hUncorr")
        hCorr   = OUT_FILE.Get(f"{pname}/h__pt_{common_suffix}_corr{common_postfix}").Clone("hCorr")
        vMC.append(hMC)
        vUncorr.append(hUncorr)
        vCorr.append(hCorr)

    # Blast-Wave (BW) fit for the pT invariant spectra
    blastwave_inv = ROOT.TF1("blastwave_inv", fit_blastwave_pt_invariant, pt_start, pt_end, 4)
    blastwave_inv.SetParNames("C", "m0", "T", "beta")
    blastwave_inv.SetParameters(100, 0.938, 0.15, 0.6)
    blastwave_inv.SetParLimits(2, 0.0, 0.9)
    blastwave_inv.SetParLimits(3, 0.0, 0.9)

    # Blast-Wave fit for the pT spectra
    blastwave = ROOT.TF1("blastwave", fit_blastwave_pt, pt_start, pt_end, 4)

    # Thermal fit for the pT invariant spectra
    thermal = ROOT.TF1("thermal", fit_thermal, pt_start, pt_end, 3)
    thermal.SetParNames("dN/dy", "T", "m0")
    thermal.SetParameters(20, 0.1, 0.938)

    # Double exponential fit for the pT invariant spectra
    fit_dexp = ROOT.TF1("fit_dexp", fit_2exp, pt_start, pt_end, 5)
    fit_dexp.SetParameters(30., 0.4, 0.2, 0.1, 0.938)


    # Bose-Einstein fit for the pT invariant spectra
    bosein = ROOT.TF1("thermal", fit_bose_einstein, pt_start, pt_end, 3)
    bosein.SetParNames("dN/dy", "T", "m0")
    bosein.SetParameters(40, 0.1, 0.938)

    plot_range = [0, 6]                                       # Common histograms plotting range
    if (p_pos == 0 or p_pos == 1):   # pi- and pi+ parameters
        plot_range = [0, 1.59]                                # Plot range for the selected particle
        fit_range       = [0.1, 2.0]                          # BW fit range for the selected particle
        fit_range_dexp  = fit_range                           # Double exponential fit range for the selected particle
        dndy_exp_limits = [0.1, 1.5]                          # Range for the dN/dy calculation using the corrected pT spectra
        blastwave_inv.SetParameters(100, p_mass, 0.2, 0.2)    # BW fit parameters for the selected particle
        blastwave_inv.FixParameter(1, p_mass)                 # Mass parameter must be fixed
        thermal.SetParameters(20, 0.1, p_mass)                # Thermal fit parameters for the selected particle
        thermal.FixParameter(2, p_mass)                       # Mass parameter must be fixed
        fit_dexp.SetParameters(50, 5, 0.3, 0.3, p_mass)       # Double exponential fit parameters for the selected particle
        fit_dexp.FixParameter(4, p_mass)                      # Mass parameter must be fixed
        bosein.FixParameter(2, p_mass)                        # Mass parameter must be fixed
    if (p_pos == 2 or p_pos == 3):   # K- and K+ -- same as above
        plot_range = [0, 1.69]                    
        fit_range       = [0.2, 1.5]
        fit_range_dexp  = fit_range
        dndy_exp_limits = [0.2, 1.5]
        blastwave_inv.SetParameters(100, p_mass, 0.2, 0.2)
        blastwave_inv.FixParameter(1, p_mass)
        thermal.FixParameter(2, p_mass)
        fit_dexp.SetParameters(20, 2, 0.3, 0.3, p_mass)
        fit_dexp.FixParameter(4, p_mass)
        bosein.FixParameter(2, p_mass)
    if (p_pos == 4):                 # p -- same as above
        plot_range = [0, 2.19]
        fit_range       = [0.5, 1.5]
        fit_range_dexp  = fit_range
        dndy_exp_limits = [0.5, 1.5]
        blastwave_inv.SetParameters(1000, p_mass, 0.1, 0.1)
        blastwave_inv.FixParameter(1, p_mass)
        if (cbin == 2):
            thermal.SetParameters(40, 0.3, p_mass)
        fit_dexp.FixParameter(4, p_mass)
        thermal.FixParameter(2, p_mass)
        bosein.FixParameter(2, p_mass)
    elif (p_pos == 6):               # d -- same as above
        plot_range = [0, 3.09]
        fit_range       = [0.9, 2.5]
        fit_range_dexp  = [1.0, 4.0]
        dndy_exp_limits = [0.8, 2.6]
        blastwave_inv.SetParameters(1000, p_mass, 0.2, 0.2)
        blastwave_inv.FixParameter(1, p_mass)
        fit_dexp.SetParameters(1.146, 24.5, 0.26, 0.22, p_mass)
        fit_dexp.FixParameter(4, p_mass)
        if (cbin == 2):
            thermal.SetParameters(80, 0.3, p_mass)
        thermal.FixParameter(2, p_mass)
        bosein.FixParameter(2, p_mass)
    elif (p_pos == 8):               # 3He -- same as above
        plot_range = [0, 3.99]
        fit_range       = [1.4, 3.5]
        fit_range_dexp  = [1.5, 4.0]
        dndy_exp_limits = [1.2, 3.9]
        blastwave_inv.FixParameter(1, p_mass)
        thermal.FixParameter(2, p_mass)
        fit_dexp.SetParameters(0.02, 50, 0.3, 0.3, p_mass)
        fit_dexp.FixParameter(4, p_mass)
        bosein.FixParameter(2, p_mass)

    canvas.SetLogy(1);                      # Set the canvas Y-axis as log

    legend = TLegend(0.7, 0.6, 0.9, 0.9)    # Create, modify and place legend
    legend.SetNColumns(1);
    legend.SetBorderSize(0);
    legend.SetFillStyle(0);
    legend.SetTextFont(42);
    legend.SetTextSize(0.04);

    for rbin in range(len(rapidity_bins)):
        outname=f'{fname}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}.{EXT}'  # Output file name for the each plot
        rapidity_postfix = f'y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'

        ratio_cor = TRatioPlot(vCorr[rbin], vMC[rbin])                               # Corrected ratio plot
        ratio_uncor = TRatioPlot(vUncorr[rbin], vMC[rbin])                           # Uncorrected ratio plot
        ratio_cor.SetSeparationMargin(0.0)                                           # The margin between upper and lower (ratio) part
        ratio_cor.SetRightMargin(0.01)                                               # The right margin
        ratio_cor.SetLeftMargin(0.12)                                                # The left margin

        make_fancy_histogram(vMC[rbin], MS_MC, LC_MC, 1, LC_MC)                      # Markers, lines and colours for histograms:
        make_fancy_histogram(vUncorr[rbin], MS_UNCOR, LC_UNCOR, 1, LC_UNCOR)         # MC, uncorrected and corrected
        make_fancy_histogram(vCorr[rbin], MS_COR, LC_COR, 1, LC_COR)

        if(rbin == 0):                                                               # The legend must be created one
            legend.AddEntry(vMC[rbin], "MC", "P")                                    # but then it can be used for each
            # ~ legend.AddEntry(vUncorr[rbin], "Uncorrected", "PL")                  # rapidity bin
            legend.AddEntry(vCorr[rbin], "Corrected", "PL")
            if((pname != "ap" and pname != "t" and pname != "He4") and postfix==''): # Here anti-protons, t and He4
                obj = ROOT.TGraph()                                                  # are not analyzed
                make_fancy_histogram(obj, 23,  1, 2, 14, 0.5, 3)
                legend.AddEntry(obj, "Blast-Wave", "PL")                             # Add BW fit to the legend
                objt = ROOT.TGraph()
                make_fancy_histogram(objt, 25,  6, 3, 6, 0.5, 3)
                legend.AddEntry(objt, "Thermal", "PL")                               # Add thermal fit to the legend
                obje = ROOT.TGraph()
                make_fancy_histogram(obje, 24,  4, 3, 4, 0.5, 3)
                legend.AddEntry(obje, "Double exp.", "PL")                           # Add double-exponential fit to the legend
                objb = ROOT.TGraph()
                make_fancy_histogram(objb, 26,  42, 3,  42, 0.5, 3)
                legend.AddEntry(objb, "Bose-Einstein", "PL")                           # Add double-exponential fit to the legend

        ratio_uncor.Draw()                                      # Draw uncorrected plot to have this object in the future
        g_uncor = ratio_uncor.GetLowerRefGraph();
        g_uncor.SetLineColor(LC_UNCOR)                          # Set the uncorrected results line colour
        ratio_cor.Draw()                                        # Draw corrected results, all futher histograms will be
        ROOT.gPad.Update()                                      # plotted on these pads (lower and upper)
        ratio_cor.GetLowerRefGraph().SetLineColor(LC_COR)       # Set the corrected results line colour
        ratio_cor.GetLowerRefGraph().SetLineWidth(2)            # and the line width
        # ~ ROOT.gPad.Update()
        y_low = ratio_cor.GetUpperPad().GetFrame().GetY1()      # The upper pad Y-axis maximum
        y_max = ratio_cor.GetUpperPad().GetFrame().GetY2()      # and minimum
        ratio_cor.GetLowerRefGraph().SetMinimum(0)              # Set the lower pad minimum to 0
        ratio_cor.GetLowerRefGraph().SetMaximum(2)              # Set the lower par maximum to 2
        ratio_cor.GetLowerRefXaxis().SetRangeUser(*plot_range)  # Set the lower (and tus upper too) pad the defined above plot range
        if(p_pos == 7):                                  # Special case for tritons
            ratio_cor.GetLowerRefGraph().SetMaximum(4)
        if(p_pos >= 6):                                  # Special case for d
            ratio_cor.GetUpperRefYaxis().SetRangeUser(math.pow(10, y_max - 3), math.pow(10, y_max+0.2)) # 3 orders of magnitude less
        if(p_pos >= 8):                                  # Special case for He3 and He4
            ratio_cor.GetUpperRefYaxis().SetRangeUser(math.pow(10, y_max - 3.5), math.pow(10, y_max+0.5)) # 3.5 orders of magnitude less
        ratio_cor.GetLowerRefYaxis().SetNdivisions(6)
        ratio_cor.GetUpperPad().cd()                            # Go to the upper pad and draw histograms and fits there
        vMC[rbin].Draw("same hist PE")                          # Monte-Carlo histogram for the current rapidity bin
        # ~ vUncorr[rbin].Draw("same hist PE")                      # Uncorrected results histogram for the current rapidity bin
        vCorr[rbin].Draw("same hist PE")                        # Corrected results histogram for the current rapidity bin
        if((pname != "ap" and pname != "t" and pname != "He4") and postfix==''):
            vCorr[rbin].Fit(thermal, "0Q", "", *fit_range)                          # Do the thermal fit
            make_fancy_histogram(thermal, 25,  6, 3, 6, 0.5, 3)                     # Make the fit result fancy
            thermal.Draw("same PL")                                                 # Draw the fit result
            save_fit(thermal, pname, f'thermal_{common_suffix}_{rapidity_postfix}') # Save fit to the output file

            vCorr[rbin].Fit("fit_dexp", "", "", *fit_range_dexp)                    # Same as for the thermal fit
            make_fancy_histogram(fit_dexp, 24,  4, 3, 4, 0.5, 3)
            fit_dexp.Draw("same PL")
            save_fit(fit_dexp, pname, f'double_exponential_{common_suffix}_{rapidity_postfix}')

            vCorr[rbin].Fit("blastwave_inv", "0Q", "", *fit_range)                  # Same as for the thermal fit
            make_fancy_histogram(blastwave_inv, 23,  1, 2, 14, 0.5, 3)
            blastwave_inv.Draw("same PL")
            save_fit(blastwave_inv, pname, f'blastwave_invariant_{common_suffix}_{rapidity_postfix}')

            vCorr[rbin].Fit(bosein, "0Q", "", *fit_range)
            make_fancy_histogram(bosein, 26,  42, 3,  42, 0.5, 3)
            bosein.Draw("same PL")
            save_fit(bosein, pname, f'bosein_{common_suffix}_{rapidity_postfix}')

            x, y = array('d'), array('d')
            n = 0
            for hbin in range(1, vMC[rbin].GetSize() - 1):      # Here the BW fit is evaluated at the same points
                if vMC[rbin].GetBinContent(hbin):               # as the Monte-Carlo histogram
                    n += 1                                      # and then the TGraph with values
                    x.append(vMC[rbin].GetBinCenter(hbin))      # y = Evaluated / Monte-Carlo is collected
                    y.append(blastwave_inv.Eval(vMC[rbin].GetBinCenter(hbin)) / vMC[rbin].GetBinContent(hbin))
            gr = ROOT.TGraph(n, x, y)
            make_fancy_histogram(gr, 23,  1, 2, 14, 0.5, 3)
            ratio_cor.GetLowerPad().cd()                        # The ratio BW / Monte-Carlo
            gr.Draw("same LP")                                  # is drawn on the lower pad
            ratio_cor.GetUpperPad().cd()

            dndy_exp     = dndy_from_pt(vCorr[rbin], *dndy_exp_limits)             # dN/dy value from corrected spectra within defined pT range
            dndy_mc_low  = dndy_from_pt(vMC[rbin],   pt_start, dndy_exp_limits[0]) # dN/dy value from MC -- low pT part
            dndy_mc_high = dndy_from_pt(vMC[rbin],   dndy_exp_limits[1], pt_end)   # dN/dy value from MC -- high pT part

            blastwave.SetParameters(blastwave_inv.GetParameters())                         # BW fit for the d^2N/dptdy spectra
            save_fit(blastwave, pname, f'blastwave_pt_{common_suffix}_{rapidity_postfix}') # Fit parameters are taken from the 
                                                                                           # BW fit for the invariant pT spectra (above)
            tbox = ROOT.TPaveText(0.15, 0.05, 0.45, 0.5, "NDC")                     # Text box with different information: fits, dN/dy etc
            bw_integral = blastwave.Integral(pt_start, pt_end)                      # Integral of the BW fit
            bw_high = blastwave.Integral(dndy_exp_limits[1], pt_end)                # Low pT part of the BW fit
            bw_low = blastwave.Integral(pt_start, dndy_exp_limits[0])               # High pT part of the BW fit
            tbox.AddText('BW fit range = {:.1f} - {:.1f} GeV/c'.format(*fit_range))
            tbox.AddText('BW Integral (0 - 6  GeV/c) = {:.3f}'.format(bw_integral))
            tbox.AddText('BW low p_{{T}}  = {:.3f} ({:.2f}%)'.format(bw_low, bw_low/bw_integral*100))
            tbox.AddText('BW high p_{{T}} = {:.3f} ({:.2f}%)'.format(bw_high, bw_high/bw_integral*100))
            tbox.AddText('dN/dy_{{exp}} [{:.2f}..{:.2f}] = {:.3e} \\pm {:.3e}'.format(*dndy_exp_limits, *dndy_exp))
            dndy_bw_low      = blastwave.Integral(pt_start, dndy_exp_limits[0])     # dN/dy value from BW fit -- low pT part
            dndy_bw_low_err  = math.fabs(dndy_bw_low - dndy_mc_low[0])              # dN/dy error from BW fit -- low pT part
            dndy_bw_high     = blastwave.Integral(dndy_exp_limits[1], pt_end)       # dN/dy value from BW fit -- high pT part
            dndy_bw_high_err = math.fabs(dndy_bw_high - dndy_mc_high[0])            # dN/dy error from BW fit -- high pT part
            tbox.AddText('dN/dy_{{BW}} [{:.2f}..{:.2f}] = {:.3e} \\pm {:.3e}'.format(0, dndy_exp_limits[0], dndy_bw_low , dndy_bw_low_err))
            tbox.AddText('dN/dy_{{BW}} [{:.2f}..{:.2f}] = {:.3e} \\pm {:.3e}'.format(dndy_exp_limits[1], pt_end, dndy_bw_high, dndy_bw_high_err))
            dndy_final     = dndy_bw_low + dndy_exp[0] + dndy_bw_high               # Final dN/dy value -- sum of previous dN/dy components
            dndy_final_err = math.sqrt(dndy_bw_low_err**2 + dndy_exp[1]**2 + dndy_bw_high_err**2) # final dN/dy error
            tbox.AddText('dN/dy_{{Final}} = {:.3e} \\pm {:.3e}'.format(dndy_final, dndy_final_err))
            # ~ tbox.Draw()
            particle_info = {}                   # A set of particle information (Python dict)
            particle_info["name"] = pname
            particle_info["centrality"] = cname
            particle_info["y"] = f'{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'
            particle_info["bw_range"] = f'{fit_range[0]} - {fit_range[1]}'
            particle_info["bw"] = f'Integral: {bw_integral}, low pT: {bw_low}, high pT: {bw_high}'
            particle_info["dndy_bw_low"] = f'{dndy_bw_low} +/- {dndy_bw_low_err}  in  [{pt_start} - {dndy_exp_limits[0]}]'
            particle_info["dndy_exp"] = f'{dndy_exp[0]} +/- {dndy_exp[1]}  in  [{dndy_exp_limits[0]} - {dndy_exp_limits[1]}]'
            particle_info["dndy_bw_high"] = f'{dndy_bw_high} +/- {dndy_bw_high_err}  in  [{dndy_exp_limits[1]} - {pt_end}]'
            particle_info["dndy"] = f'{dndy_final} +/- {dndy_final_err}'
            if TBL_FILE is not None:
                with open(TBL_FILE, 'a') as f:                   # All information above will be written in the 
                    f.write(json.dumps(particle_info, indent=4)) # output JSON file (TBL_FILE)

        legend.Draw()
        ratio_cor.GetLowerPad().cd()
        # ~ g_uncor.Draw()
        canvas.Print(outname)
        canvas.Clear()
    canvas.SetLogy(0)
    write_status(f'{function_name} {CURRENT_ST[order]}')



##  This subroutine draws the \f$ dN/dy \f$ spectra
#  \param do_dedx The switch for the PID methods: 
#     - 1 = make plots for the dE/dx only PID 
#     - 0 = make plots for the combined PID 
#  \param canvas TCanvas prepared for drawing
#  \param p_pos The particle index
#  \param cbin  The centrality bin index
#  \param fname The output pdf-file name
#
def draw_dndy(do_dedx, canvas, p_pos, cbin, fname):
    order = 2
    CURRENT_ST[order] += 1
    function_name = sys._getframe().f_code.co_name
    if skip_function(function_name, order): return

    pname    = list(Particles)[p_pos]
    cbins    = list(Centrality)[cbin]
    cname    = f'{cbins[0]} - {cbins[1]}%'
    pt_start = list(Particles.values())[0]["pt_bins"][1]
    pt_end   = list(Particles.values())[0]["pt_bins"][2]
    print(f' draw dN/dy =============== {pname},  {cname},  {pt_start} < pT < {pt_end}  =============== ')

    if(do_dedx):
        postfix='_dedx'
    else:
        postfix=''
    canvas.cd()
    vMC   = []
    vCorr = []
    vBW   = []
    common_suffix = f'{pname}_centrality{cbin}'
    for rbin in range(len(rapidity_bins)):
        rapidity_postfix = f'y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'
        common_postfix = f'{postfix}_{rapidity_postfix}'
        mc   = OUT_FILE.Get(f"{pname}/h__pt_{common_suffix}_mc_{rapidity_postfix}").Clone("mc")
        corr = OUT_FILE.Get(f"{pname}/h__pt_{common_suffix}_corr{common_postfix}").Clone("corr")
        fit  = OUT_FILE.Get(f"{pname}/blastwave_pt_{common_suffix}_{rapidity_postfix}").Clone("fit")
        vMC.append(mc)
        vCorr.append(corr)
        vBW.append(fit)

    if (p_pos == 0 or p_pos == 1):   # pi- and pi+
        dndy_exp_limits = [0.1, 1.5]
    if (p_pos == 2):                 # K-
        dndy_exp_limits = [0.2, 1.5]
    if (p_pos == 3):                 # K+
        dndy_exp_limits = [0.2, 1.5]
    if (p_pos == 4):                 # p
        dndy_exp_limits = [0.5, 1.5]
    elif (p_pos == 6):               # d
        dndy_exp_limits = [0.8, 2.6]
    elif (p_pos == 8):               # 3He
        dndy_exp_limits = [1.2, 3.9]

    x_mc, y_mc   = array('d'), array('d')
    ex_mc, ey_mc = array('d'), array('d')
    x, y   = array('d'), array('d')
    ex, ey = array('d'), array('d')
    n = 0
    for rbin in range(len(rapidity_bins)):
        dndy_exp         = dndy_from_pt(vCorr[rbin], *dndy_exp_limits)

        dndy_mc_low      = dndy_from_pt(vMC[rbin],   pt_start, dndy_exp_limits[0])
        dndy_bw_low      = vBW[rbin].Integral(pt_start, dndy_exp_limits[0])
        dndy_bw_low_err  = math.fabs(dndy_bw_low - dndy_mc_low[0])

        dndy_mc_high     = dndy_from_pt(vMC[rbin],   dndy_exp_limits[1],   pt_end)
        dndy_bw_high     = vBW[rbin].Integral(dndy_exp_limits[1],   pt_end)
        dndy_bw_high_err = math.fabs(dndy_bw_high - dndy_mc_high[0])

        dndy_final       = dndy_bw_low + dndy_exp[0] + dndy_bw_high
        dndy_final_err   = math.sqrt(dndy_bw_low_err**2 + dndy_exp[1]**2 + dndy_bw_high_err**2)

        n += 1
        x.append((rapidity_bins[rbin][0] + rapidity_bins[rbin][1]) / 2)
        y.append(dndy_final)
        ex.append(0)
        ey.append(dndy_final_err)
        
        dndy_orig      = dndy_from_pt(vMC[rbin],   pt_start, pt_end)
        x_mc.append((rapidity_bins[rbin][0] + rapidity_bins[rbin][1]) / 2)
        y_mc.append(dndy_orig[0])
        ex_mc.append(0)
        ey_mc.append(dndy_orig[1])

    gr = ROOT.TGraphErrors(n, x, y, ex, ey)
    make_fancy_histogram(gr, 20, 1, 1, 1, msize = 2)
    hmax = ROOT.TMath.MaxElement(n,gr.GetY()) * 1.5 # add 20% on the top
    hmin = ROOT.TMath.MinElement(n,gr.GetY()) * 0.5 # add 20% on the bottom
    if hmin < 0.001: hmin = 0
    frame = canvas.DrawFrame(-1.19, hmin, 1.19, hmax)
    make_fancy_frame(frame, f'{SYS}, {pname}, {cname}', "y", "dN/dy")
    gr.Draw("same PE")

    gr_mc = ROOT.TGraphErrors(n, x_mc, y_mc, ex_mc, ey_mc)
    make_fancy_histogram(gr_mc, 24, 2, 1, 2, msize = 2)
    gr_mc.Draw("same PE")


    canvas.Print(fname)
    canvas.Clear()
    canvas.SetLogy(0)
    write_status(f'{function_name} {CURRENT_ST[order]}')



##  This subroutine draws the coalescence parameter B2
#  \param canvas TCanvas prepared for drawing
#  \param cbin  The centrality bin index
#  \param fname The output pdf-file name
#
def draw_b2(canvas, cbin, fname):
    order = 4
    CURRENT_ST[order] += 1
    function_name = sys._getframe().f_code.co_name
    if skip_function(function_name, order): return

    cbins    = list(Centrality)[cbin]
    cname    = f'{cbins[0]} - {cbins[1]}%'

    canvas.cd()
    rc_p = []
    rc_d = []

    for rbin in range(len(rapidity_bins)):
        rapidity_postfix = f'y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'
        hp   = OUT_FILE.Get(f"p/h__pt_p_centrality{cbin}_corr_{rapidity_postfix}").Clone("hp");
        hd   = OUT_FILE.Get(f"d/h__pt_d_centrality{cbin}_corr_{rapidity_postfix}").Clone("hd");
        rc_p.append(hp)
        rc_d.append(hd)

    p_exp_limits = [0.5, 1.5]
    p_pt_start = list(Particles.values())[4]["pt_bins"][1]
    p_pt_end   = list(Particles.values())[4]["pt_bins"][2]
    p_pt_bins  = list(Particles.values())[4]["pt_bins"][0]
    p_binw     = (p_pt_end - p_pt_start) / p_pt_bins
    p_exp_bins = (p_exp_limits[1] - p_exp_limits[0]) / p_binw

    d_exp_limits = [0.8, 2.6]
    d_pt_start = list(Particles.values())[6]["pt_bins"][1]
    d_pt_end   = list(Particles.values())[6]["pt_bins"][2]
    d_pt_bins  = list(Particles.values())[6]["pt_bins"][0]
    d_binw     = (d_pt_end - d_pt_start) / d_pt_bins
    d_exp_bins = (d_exp_limits[1] - d_exp_limits[0]) / d_binw

    b2_limits = [0, 0]
    if p_exp_limits[0]*2. < d_exp_limits[0]:  b2_limits[0] = d_exp_limits[0]
    else: b2_limits[0] = p_exp_limits[0]*2.
    if p_exp_limits[1]*2. > d_exp_limits[1]:  b2_limits[1] = d_exp_limits[1]
    else: b2_limits[1] = p_exp_limits[1]*2.

    b2_bins = int((b2_limits[1] - b2_limits[0]) / d_binw)
    scale_f = 1E3
    for rbin in range(len(rapidity_bins)):
        x, y   = array('d'), array('d')
        ex, ey = array('d'), array('d')
        n = 0
        for b2bin in range(b2_bins):
            d_llim = b2_limits[0] + d_binw*b2bin
            d_rlim = b2_limits[0] + d_binw + d_binw*b2bin
            p_llim = b2_limits[0]/2. + p_binw*b2bin
            p_rlim = b2_limits[0]/2. + p_binw + p_binw*b2bin
            nd = get_ptbin_content(rc_d[rbin], d_llim, d_rlim)
            np = get_ptbin_content(rc_p[rbin], p_llim, p_rlim)
            x.append(d_llim + d_binw/2.)
            y.append(nd[0] / (np[0]*np[0])  * scale_f)
            ex.append(0)
            ey.append(0)
            n += 1

        gr = ROOT.TGraphErrors(n, x, y, ex, ey)
        make_fancy_histogram(gr, 20, 1, 1, 1, msize = 2)
        hmax = ROOT.TMath.MaxElement(n,gr.GetY()) * 1.2 # add 20% on the top
        hmin = ROOT.TMath.MinElement(n,gr.GetY()) * 0.8 # add 20% on the bottom
        frame = canvas.DrawFrame(b2_limits[0]*0.9, hmin, b2_limits[1]*1.1, hmax)
        power = math.floor(math.log10(scale_f))
        make_fancy_frame(frame, f'{SYS}, {rapidity_bins[rbin][0]} < y < {rapidity_bins[rbin][1]}, {cname}', 
                        'p_{T}, GeV/c', f'B2 #times 10^{{{power}}}, GeV^{{2}}/c^{{3}}')
        gr.Draw("same PE")
        outname=f'{fname}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}.{EXT}'
        canvas.Print(outname)
        canvas.Clear()
    canvas.SetLogy(0)
    write_status(f'{function_name} {CURRENT_ST[order]}')



##  This subroutine makes the \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra for the defined rapidity bins
#  \param do_corr Switch for the efficiency and contamination corrections:
#     - 1 = apply
#     - 0 = do not apply
#  \param do_dedx The switch for the PID methods: 
#     - 1 = make plots for the dE/dx only PID 
#     - 0 = make plots for the combined PID 
#  \param particle The particle name (p, d, He4 etc)
#  \param cbin The centrality bin (0, 1, 2, etc)
#  \param psname The name of the selected phase space histogram in the input file
#  \param postfix The phase space histogram postfix, can be:
#     - 'mc' - For the Monte-Carlo phase space.
#     - 'uncor' - For the uncorrected reconstructed phase space.
#     - 'corr' - For the corrected reconstructed phase space.
#
def make_pt_spectra(do_corr, do_dedx, particle, cbin, psname, postfix):
    order = 0
    CURRENT_ST[order] += 1
    function_name = sys._getframe().f_code.co_name
    if skip_function(function_name, order): return

    cbins    = list(Centrality)[cbin]
    cname    = f'{cbins[0]} - {cbins[1]}%'

    if((postfix == "uncorr" or postfix == "corr") and do_dedx):
        postfix=f'{postfix}_dedx'
    elif ((postfix == "uncorr" or postfix == "corr") and not do_dedx):
        postfix=f'{postfix}'
    hPhaseSpace = INP_FILE.Get(psname).Clone(f"phasespace_{postfix}_{particle}_centrality{cbin}")
    hEvents = INP_FILE.Get("h__events").Clone("hEvents")

    name_ptspectra = f"h__pt_{particle}_centrality{cbin}_{postfix}"

    if(do_corr):
      # --- TPC, Secondaries
        # ------- efficiency
        hEfficiency = calculate_efficiency(f"h__eff_tpc_numerator_{particle}_centrality{cbin}",
                                           f"h__eff_tpc_denominator_{particle}_centrality{cbin}",
                                           f"h__efficiency_tpc_{particle}_centrality{cbin}", particle)
        apply_efficiency(hPhaseSpace, hEfficiency)
        # ------- contamination
        hContamination = calculate_efficiency(f"h__cont_sec_numerator_{particle}_centrality{cbin}",
                                              f"h__cont_sec_denominator_{particle}_centrality{cbin}",
                                              f"h__contamination_secondaries_{particle}_centrality{cbin}", particle)
        apply_contamination(hPhaseSpace, hContamination)

        if(do_dedx):
          # --- PID TPC only (dE/dx)
            # ------- efficiency
            hEfficiency = calculate_efficiency(f"h__eff_pid_dedx_numerator_{particle}_centrality{cbin}",
                                               f"h__eff_tof_denominator_{particle}_centrality{cbin}",
                                               f"h__efficiency_pid_dedx_{particle}_centrality{cbin}", particle)
            apply_efficiency(hPhaseSpace, hEfficiency)
            # ------- contamination
            hContamination = calculate_efficiency(f"h__cont_pid_dedx_numerator_{particle}_centrality{cbin}",
                                                  f"h__cont_pid_dedx_denominator_{particle}_centrality{cbin}",
                                                  f"h__contamination_pid_dedx_{particle}_centrality{cbin}", particle)
            apply_contamination(hPhaseSpace, hContamination)
        else:
          # --- ToF
            # ------- efficiency
            hEfficiency = calculate_efficiency(f"h__eff_tof_numerator_{particle}_centrality{cbin}",
                                               f"h__eff_tof_denominator_{particle}_centrality{cbin}",
                                               f"h__efficiency_tof_{particle}_centrality{cbin}", particle)
            apply_efficiency(hPhaseSpace, hEfficiency)
            # ------- contamination
            hContamination = calculate_efficiency(f"h__cont_tof_numerator_{particle}_centrality{cbin}",
                                                  f"h__eff_pid_denominator_{particle}_centrality{cbin}",
                                                  f"h__contamination_tof_{particle}_centrality{cbin}", particle)
            apply_contamination(hPhaseSpace, hContamination)

          # --- PID combined
            # ------- efficiency
            hEfficiency = calculate_efficiency(f"h__eff_pid_numerator_{particle}_centrality{cbin}",
                                               f"h__eff_pid_denominator_{particle}_centrality{cbin}",
                                               f"h__efficiency_pid_{particle}_centrality{cbin}", particle)
            apply_efficiency(hPhaseSpace, hEfficiency)
            # ------- contamination
            hContamination = calculate_efficiency(f"h__cont_pid_numerator_{particle}_centrality{cbin}",
                                                  f"h__cont_pid_denominator_{particle}_centrality{cbin}",
                                                  f"h__contamination_pid_{particle}_centrality{cbin}", particle)
            apply_contamination(hPhaseSpace, hContamination)

    dn = hEvents.GetBinContent(cbin+1);
    hPhaseSpace.Scale(1./ dn);
    OUT_FILE.cd(f'{particle}')
    hPhaseSpace.Write()
    for rbin in range(len(rapidity_bins)):
        rb_low = rapidity_bins[rbin][0]
        rb_high = rapidity_bins[rbin][1]
        rb_width = math.fabs(rb_high - rb_low)

        for hbin in range(1,hPhaseSpace.GetNbinsX()):
            if(math.fabs(hPhaseSpace.GetXaxis().GetBinLowEdge(hbin) - rb_low) < SMALL_DELTA):
                left_edge = hbin
                break
        right_edge = int(left_edge + rb_width / hPhaseSpace.GetXaxis().GetBinWidth(left_edge) - 1)
        res_name = f"{name_ptspectra}_y{rb_low}_{rb_high}"
        hResult = hPhaseSpace.ProjectionY(res_name, left_edge, right_edge)
        hResult.Scale(1./ rb_width);
        for k in range(1, hResult.GetSize() - 1):
            content = hResult.GetBinContent(k);  # N
            if(content <= 0):
                continue
            pt_binw = hResult.GetBinWidth(k) # dpT
            error   = hResult.GetBinError(k)
            pt_mean = hResult.GetBinCenter(k) # pT
            hResult.SetBinContent(k, content / (pt_mean * pt_binw)) # N / pt * dpt * dy
            hResult.SetBinError(k, error / (pt_mean * pt_binw))     # and bin error
        hResult.GetXaxis().SetTitle('p_{T}, GeV/c')
        hResult.GetYaxis().SetTitle('d^{2}n/p_{T}dp_{T}dy')
        hResult.SetTitle(f'{SYS}, {particle}, {cname}, {rb_low} < y < {rb_high}')
        hResult.Write()
    OUT_FILE.cd()
    write_status(f'{function_name} {CURRENT_ST[order]}')



##  This subroutine makes the overall efficiency histograms (TPC + ToF + PID)
#  \param do_dedx The switch for the PID methods: 
#     - 1 = make plots for the dE/dx only PID 
#     - 0 = make plots for the combined PID 
#  \param particle The particle name (p, d, He4 etc)
#  \param cbin The centrality bin (0, 1, 2, etc)
#
def make_overall_efficiency(do_dedx, particle, cbin):
    order = 5
    CURRENT_ST[order] += 1
    function_name = sys._getframe().f_code.co_name
    if skip_function(function_name, order): return
    cbins    = list(Centrality)[cbin]
    cname    = f'{cbins[0]} - {cbins[1]}%'

    # ------- efficiency
    eff = OUT_FILE.Get(f'{particle}/h__efficiency_tpc_{particle}_centrality{cbin}').Clone('EFF')
    overall_efficiency = eff.Clone(f'h__overall_efficiency_{particle}_centrality{cbin}')
    eff = OUT_FILE.Get(f'{particle}/h__efficiency_tof_{particle}_centrality{cbin}').Clone('EFF')
    overall_efficiency.Multiply(eff)
    eff = OUT_FILE.Get(f'{particle}/h__efficiency_pid_{particle}_centrality{cbin}').Clone('EFF')
    overall_efficiency.Multiply(eff)

    # ------- contamination
    cont  = OUT_FILE.Get(f'{particle}/h__contamination_secondaries_{particle}_centrality{cbin}').Clone('CONT')
    overall_contamination = cont.Clone(f'h__overall_contamination_{particle}_centrality{cbin}')
    cont  = OUT_FILE.Get(f'{particle}/h__contamination_tof_{particle}_centrality{cbin}').Clone('CONT')
    overall_contamination.Multiply(cont)
    cont  = OUT_FILE.Get(f'{particle}/h__contamination_pid_{particle}_centrality{cbin}').Clone('CONT')
    overall_contamination.Multiply(cont)


    OUT_FILE.cd(f'{particle}')
    overall_efficiency.Write()
    overall_contamination.Write()
    for rbin in range(len(rapidity_bins)):
        rb_low = rapidity_bins[rbin][0]
        rb_high = rapidity_bins[rbin][1]
        rb_width = math.fabs(rb_high - rb_low)

        for hbin in range(1, overall_efficiency.GetNbinsX()):
            if(math.fabs(overall_efficiency.GetXaxis().GetBinLowEdge(hbin) - rb_low) < SMALL_DELTA):
                left_edge = hbin
                break
        right_edge = int(left_edge + rb_width / overall_efficiency.GetXaxis().GetBinWidth(left_edge) - 1)
        res_name = f"overall_efficiency_{particle}_centrality{cbin}_y{rb_low}_{rb_high}"
        hResult = overall_efficiency.ProjectionY(res_name, left_edge, right_edge)
        hResult.Scale(1. / (right_edge - left_edge))
        hResult.GetXaxis().SetTitle('p_{T}, GeV/c')
        hResult.GetYaxis().SetTitle('Overall efficiency')
        hResult.SetTitle(f'{SYS}, {particle}, {cname}, {rb_low} < y < {rb_high}')
        hResult.Write()

    for rbin in range(len(rapidity_bins)):
        rb_low = rapidity_bins[rbin][0]
        rb_high = rapidity_bins[rbin][1]
        rb_width = math.fabs(rb_high - rb_low)

        for hbin in range(1, overall_contamination.GetNbinsX()):
            if(math.fabs(overall_contamination.GetXaxis().GetBinLowEdge(hbin) - rb_low) < SMALL_DELTA):
                left_edge = hbin
                break
        right_edge = int(left_edge + rb_width / overall_contamination.GetXaxis().GetBinWidth(left_edge) - 1)
        res_name = f"overall_contamination_{particle}_centrality{cbin}_y{rb_low}_{rb_high}"
        hResult = overall_contamination.ProjectionY(res_name, left_edge, right_edge)
        hResult.Scale(1. / (right_edge - left_edge))
        hResult.GetXaxis().SetTitle('p_{T}, GeV/c')
        hResult.GetYaxis().SetTitle('Overall contamination')
        hResult.SetTitle(f'{SYS}, {particle}, {cname}, {rb_low} < y < {rb_high}')
        hResult.Write()
    OUT_FILE.cd()
    write_status(f'{function_name} {CURRENT_ST[order]}')



##  This subroutine draws the overall efficiency histograms
#  \param canvas TCanvas prepared for drawing
#  \param p_pos The particle index
#  \param cbin  The centrality bin index
#  \param fname The output pdf-file name
#
def draw_overall_efficiency(canvas, p_pos, cbin, fname):
    order = 6
    CURRENT_ST[order] += 1
    function_name = sys._getframe().f_code.co_name
    if skip_function(function_name, order): return

    pname    = list(Particles)[p_pos]
    cbins    = list(Centrality)[cbin]
    cname    = f'{cbins[0]} - {cbins[1]}%'

    canvas.cd()
    vEff = []
    for rbin in range(len(rapidity_bins)):
        common_postfix = f'y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}'
        hEff     = OUT_FILE.Get(f"{pname}/overall_efficiency_{pname}_centrality{cbin}_{common_postfix}").Clone("hEff")
        vEff.append(hEff)

    plot_range = [0, 6]
    if (p_pos == 0 or p_pos == 1):   # pi- and pi+
        plot_range = [0, 1.59]
    if (p_pos == 2 or p_pos == 3):   # K- and K+
        plot_range = [0, 1.69]
    if (p_pos == 4):   # p
        plot_range = [0, 2.29]
    elif (p_pos == 6): # d
        plot_range = [0, 3.09]
    elif (p_pos == 8): # 3He
        plot_range = [0, 4.29]

    for rbin in range(len(rapidity_bins)):
        outname=f'{fname}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}.{EXT}'
        make_fancy_histogram(vEff[rbin], MS_MC, LC_MC, 1, LC_MC, msize = 2)
        frame = canvas.DrawFrame(plot_range[0], 0, plot_range[1], 1.19)
        make_fancy_frame(frame, f'{SYS}, {pname}, {cname}, {rapidity_bins[rbin][0]} < y < {rapidity_bins[rbin][1]}',
                         "p_{T}, GeV/c", "Overall efficiency")
        vEff[rbin].Draw("same hist P")
        canvas.Print(outname)
        canvas.Clear()
    canvas.SetLogy(0)
    write_status(f'{function_name} {CURRENT_ST[order]}')


##  This subroutine draws the frame for the final histograms
#  \param frame The TH1F frame
#  \param title The title of the frame (histogram)
#  \param x_title The title of x-axis
#  \param y_title The title of y-axis
#
def make_fancy_frame(frame, title, x_title, y_title):
    frame.SetTitle(title)
    xax = frame.GetXaxis()
    xax.SetTitle(x_title)
    xax.SetTitleFont(42)
    xax.SetTitleSize(0.06)
    xax.SetTitleOffset(1.05)
    xax.SetLabelFont(42)
    xax.SetLabelSize(0.06)
    xax.SetLabelOffset(0.01)
    yax = frame.GetYaxis()
    yax.SetTitle(y_title)
    yax.SetTitleFont(42)
    yax.SetTitleSize(0.06)
    yax.SetTitleOffset(1.1)
    yax.SetLabelFont(42)
    yax.SetLabelSize(0.06)
    yax.SetLabelOffset(0.01)
    yax.SetNdivisions(8)



##  This makes the histograms fancy
#  \param hist The histogram to modify
#  \param mstyle The marker style
#  \param mcolor The marker colour
#  \param lstyle The line style
#  \param lcolor The line colour
#
#  Optionally one can set
#  \param msize The size of the marker
#  \param lwidth The width of the line
#
def make_fancy_histogram(hist, mstyle, mcolor, lstyle, lcolor, msize = None, lwidth = None):
    hist.SetMarkerStyle(mstyle)
    hist.SetMarkerColor(mcolor)
    if msize:
        hist.SetMarkerSize(msize)
    hist.SetLineStyle(lstyle)
    hist.SetLineColor(lcolor)
    if lwidth:
        hist.SetLineWidth(lwidth)



##  This subroutine calculates the efficiency (contamination) correction
#  \param num The name of the selected efficiency (or contamination) numerator in the input file
#  \param denom The name of the selected efficiency (or contamination) denominator in the input file
#  \param res The resulting histogram name
#  \param particle The particle name (p, d, He4 etc)
#
def calculate_efficiency(num, denom, res, particle):
    hNumerator   = EFF_FILE.Get(num).Clone("hNumerator")
    hDenominator = EFF_FILE.Get(denom).Clone("hDenominator")
    hEfficiency  = hNumerator.Clone(res)
    hEfficiency.Divide(hDenominator)
    OUT_FILE.cd(f'{particle}')
    hEfficiency.Write()
    OUT_FILE.cd()
    return hEfficiency



##  This subroutine applies the efficiency correction: phase-space / efficiency
#  \param phasespace The phase space to correct
#  \param eff The efficiency
#
def apply_efficiency(phasespace, eff):
    phasespace.Divide(eff) # the usual ROOT way



##  This subroutine applies the contamination correction as: phase-space * (1 - contamination)
#  \param phasespace The phase space to correct
#  \param cont The contamination
#
def apply_contamination(phasespace, cont):
    for i in range(1, cont.GetNbinsX()):
        for j in range(1, cont.GetNbinsY()):
            cont_content = cont.GetBinContent(i, j)
            ps_content   = phasespace.GetBinContent(i, j)
            phasespace.SetBinContent(i, j, ps_content * (1. - cont_content))



##  This subroutine saves the fit function
#  \param fit The fit function object (e.g. TF1)
#  \param pname The particle name
#  \param name The name of the fit function
#
def save_fit(fit, pname, name):
    res = fit.Clone(name)
    OUT_FILE.cd(f'{pname}')
    res.Write()
    OUT_FILE.cd()



##  This subroutine calculates the dN/dy for the selected \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param hist The \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param pt_low The lower edge of the first \f$ p_{T} \f$ bin
#  \param pt_high The upper edge of the last \f$ p_{T} \f$ bin
#
#  This routine returns the dN/dy value (dndy[0]) and its error (dndy[1])
#
def dndy_from_pt(hist, pt_low, pt_high):
    dndy   = [0, 0]
    counts = 0
    err    = 0
    for h_bin in range(1, hist.GetNbinsX()+1):
        if(math.fabs(hist.GetBinLowEdge(h_bin) - pt_low) < SMALL_DELTA):
            left_edge = h_bin
        if(math.fabs(hist.GetBinLowEdge(h_bin) + hist.GetBinWidth(h_bin)  - pt_high) < SMALL_DELTA):
            right_edge = h_bin
            break
    for h_bin in range(left_edge, right_edge+1):
        n   = hist.GetBinContent(h_bin) # N
        pt  =  hist.GetBinCenter(h_bin) # pT
        dpt =  hist.GetBinWidth(h_bin)  # dpT
        e   =  hist.GetBinError(h_bin)  # error
        if n:
            counts += n * pt * dpt
            err += math.pow(e * pt, 2)
    dndy[0] = counts
    dndy[1] = math.sqrt(err)
    return dndy 



##  This subroutine calculates the dN/dy for the selected \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param hist The \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra
#  \param pt_low The lower edge of the \f$ p_{T} \f$ bin
#  \param pt_high The upper edge of the \f$ p_{T} \f$ bin
# 
#   The returned value is a set of the bin content (res[0]) and error (res[1])
#
def get_ptbin_content(hist, pt_low, pt_high):
    res = [0, 0]
    for h_bin in range(1, hist.GetNbinsX()+1):
        if(math.fabs(hist.GetBinLowEdge(h_bin) - pt_low) < SMALL_DELTA):
            left_edge = h_bin
        if(math.fabs(hist.GetBinLowEdge(h_bin) + hist.GetBinWidth(h_bin)  - pt_high) < SMALL_DELTA):
            right_edge = h_bin
            break
    for h_bin in range(left_edge, right_edge+1):
        n   = hist.GetBinContent(h_bin) # N
        e   = hist.GetBinError(h_bin)   # error
        if n:
            res[0] = n
            res[1] = e
    return res



##  This subroutine converts the selected \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra to \f$ d^{2}N/dp_{T}dy \f$
#  \param hist The \f$ d^{2}N/p_{T}dp_{T}dy \f$ spectra to convert
#
def convert_to_pt(hist):
    for h_bin in range(1, hist.GetNbinsX()+1):
        n   = hist.GetBinContent(h_bin) # N
        pt  =  hist.GetBinCenter(h_bin) # pT
        dpt =  hist.GetBinWidth(h_bin)  # dpT
        e   =  hist.GetBinError(h_bin)  # error
        if n:
            hist.SetBinContent(h_bin, n * pt)
            hist.SetBinError(h_bin, e * pt)



##  This subroutine checks if the file exists
#  \param fname The file name to check
#
def check_file_exist(fname):
    if(not os.path.isfile(fname)):
        print(f"--- {fname} does not exist")
        exit(-1)



##  This subroutine checks if the ROOT file can be opened and is not "zombie"
#  \param file The ROOT file to check
#
def check_root_file(file):
    if (not file or file.IsZombie()):
        print(f"--- Can not open {file}")
        exit(-1)



##  This subroutine reads the saved progress of the funtion from the ASCII file (STATUS_FILE)
#  \param func The function name
#
def read_status(func):
    res = ['Empty', 0]
    if os.path.exists(STATUS_FILE):
        with open(STATUS_FILE, "r") as f:
            for line in f:
                if f'{func}' in line:
                    res = line.rstrip('\n').split(' ')
                    break
    else:
        with open(STATUS_FILE, 'w') as document: pass
    
    return(int(res[1]))



##  This subroutine checks if the current function call counter is less or equal than the saved progress (skip) or not (do not skip)
#  \param fnc_name The function name
#  \param order The function execution order
#
def skip_function(fnc_name, order):
    saved_state = read_status(fnc_name)
    global CURRENT_ST
    global SKIP
    if SKIP and CURRENT_ST[order] <= saved_state:
        return True
    return False



##  This subroutine saves the current progress to the ASCII file (STATUS_FILE)
#  \param progress String of the format "(string)function_name (integer)counter" e.g.: ```draw_pt_something 2```
#
def write_status(progress):
    func, position = progress.rstrip('\n').split(' ')
    do_update = 0
    with open(STATUS_FILE, "r+") as file:
        for line in file:
            if f'{func}' in line:
                do_update = 1
                break
        else:
            file.write(f'{func} {position}\n')
    if do_update:
        for line in fileinput.input(STATUS_FILE, inplace=True):
            print(line.replace(f'{func} {int(position) - 1}', f'{func} {position}'), end='')



##  This subroutine generates the pdf file with the report
#  \param is_extended Switch for the extended QA plots production:
#     - 1 = on
#     - 0 = off
#  \param dirname The input directory with produced plots to include into the report
#  \param fname The output report file (without the extention)
#
def generate_report(is_extended, dirname, fname):
    geometry_options = {"tmargin": "1cm", "lmargin": "1cm", "rmargin": "1cm", "bmargin": "1.5cm"}
    doc = Document(geometry_options=geometry_options)
    doc.packages.append(Package('hyperref'))
    doc.packages.append(Package('graphicx'))
    cf.active = cf.Version1(indent=False)
    doc.append(Command('tableofcontents'))
    doc.append(NoEscape(r'\clearpage'))
    for part in Particles:
        with doc.create(Section(f'{part}')):
            if is_extended:
                list_efficiencies   = ["tpc", "tof", "pid"]
                list_contaminations = ["secondaries", "pid", "tof"]
                if DEDX:
                    list_efficiencies.append("pid_dedx")
                    list_contaminations.append("pid_dedx")
                for eff in list_efficiencies:
                    subsname=eff.upper()
                    if(eff == "pid_dedx"):
                        subsname = "PID dE/dx only"
                    with doc.create(Subsection(f'{subsname} efficiency')):
                        doc.append(NoEscape(r'\resizebox{\textwidth}{!} {'))
                        for cbin in range(0, len(Centrality)):
                            doc.append(Command('includegraphics', NoEscape(f'{dirname}/efficiency/{eff}/{part}_centrality{cbin}.{EXT}')))
                        doc.append(NoEscape(r'}'))
                doc.append(NoEscape(r'\newpage')) 
                for eff in list_contaminations:
                    subsname=eff.upper()
                    if(eff == "pid_dedx"):
                        subsname = "PID dE/dx only"
                    with doc.create(Subsection(f'{subsname} contamination')):
                        doc.append(NoEscape(r'\resizebox{\textwidth}{!} {'))
                        for cbin in range(0, len(Centrality)):
                            doc.append(Command('includegraphics', NoEscape(f'{dirname}/contamination/{eff}/{part}_centrality{cbin}.{EXT}')))
                        doc.append(NoEscape(r'}'))
                doc.append(NoEscape(r'\newpage')) 
            with doc.create(Subsection('Results')):
                doc.append(NoEscape(r'\begin{center}'))
                for cbin in range(0, len(Centrality)):
                    doc.append(NoEscape(r'\resizebox{0.94\textwidth}{!} {'))
                    for rbin in range(len(rapidity_bins)):
                        spectra = f'{dirname}/results/pt/{part}/{part}_centrality{cbin}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}.{EXT}'
                        doc.append(Command('includegraphics', NoEscape(spectra)))
                    doc.append(NoEscape(r'}\\'))
                doc.append(NoEscape(r'\end{center}'))
                if part != "ap" and part != "t" and part != "He4":
                    doc.append(NoEscape(r'\begin{center}'))
                    doc.append(NoEscape(r'\resizebox{\textwidth}{!} {'))
                    for cbin in range(0, len(Centrality)):
                        spectra = f'{dirname}/results/dndy/{part}/{part}_centrality{cbin}.{EXT}'
                        doc.append(Command('includegraphics', NoEscape(spectra)))
                    doc.append(NoEscape(r'}'))
                    doc.append(NoEscape(r'\end{center}'))
                if part == "d":
                    doc.append(NoEscape(r'\begin{center}'))
                    for cbin in range(0, len(Centrality)):
                        doc.append(NoEscape(r'\resizebox{\textwidth}{!} {'))
                        for rbin in range(len(rapidity_bins)):
                            spectra = f'{dirname}/results/coal/b2_centrality{cbin}_y{rapidity_bins[rbin][0]}_{rapidity_bins[rbin][1]}.{EXT}'
                            doc.append(Command('includegraphics', NoEscape(spectra)))
                        doc.append(NoEscape(r'}\\'))
                    doc.append(NoEscape(r'\end{center}'))
                doc.append(NoEscape(r'\newpage'))    
        doc.append(NoEscape(r'\clearpage'))
    doc.generate_pdf(fname, clean_tex=False)



if __name__ == "__main__":
    main()
