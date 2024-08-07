#ifndef ENERCLASS1_H
#define ENERCLASS1_H

#include "MpdAnalysisTask.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanFilter.h"
#include "MpdMCTrack.h"
#include "MpdVertex.h"
#include "MpdPid.h"
#include "MpdPidQA.h"
#include "MpdTofMatchingData.h"
#include "MpdTofMatching.h"
#include "MpdHelix.h"
//#include "MpdEventHeader.h"

#include "TH2.h"
#include "TString.h"

#include <string>

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

class EnerClass1: public MpdAnalysisTask, public MpdPid {
  MpdPid  *fPID;
  
public:
  EnerClass1();
  EnerClass1(const char *name, const char *outputName = "taskName", const char* settings_filename = "file.json");
  
  ~EnerClass1();  
  
  void UserInit();
  void ProcessEvent(MpdAnalysisEvent &event);
  void Finish();
  
  struct particle_info{
    TString name;                ///< Particle name
    int     pdg;                 ///< Monte Carlo particle code, [see the scheme](https://pdg.lbl.gov/2022/reviews/rpp2022-rev-monte-carlo-numbering.pdf)
    float   mass;                ///< Particle mass
    int     charge;              ///< Particle charge
    short   enum_position;       ///< Particle position in the `core/mpdBase/MpdTrack.h`
    float   pt_bins[3];          ///< Number of pT bins, left limit, right limit.
    float   rapidity_bins[3];    ///< Number of rapidity bins, left limit, right limit.
  };
  
  
  private:
  bool isInitialized = false;
  std::string mParamConfig;
  const char* settings_file;
  std::map<int,int> map_tof_matching;
  
//__________________________________________________________________________________

// Global settings
  bool        s__gl_Verbose;     ///< Switch to handle an extra output
  short       s__gl_nMpdPID;     ///< Number of particles predefined for the PID in the `core/mpdBase/MpdTrack.h`
  bool        s__gl_MC;          ///< Switch for the creation of the histograms for the particles identified by the MC information
  bool        s__gl_Efficiency;  ///< Switch for the creation of the efficiency histograms
  short       s__gl_PID;         ///< Switch for PID mode: 0 = "competitive"; 1 = "all you can take"; 2 = "MpdPid"
  short       s__gl_DCA;         ///< Switch for DCA mode: 0 = "classic"; 1 = "NSigma";
  short       s__gl_TOF;         ///< Switch for ToF matching mode: 0 = "classic" -- check for the ToF flags 2 and 6; 1 = "NSigma" -- check for GetTofDphiSigma() and GetTofDzSigma(); 2 = use MpdTofMatchingData
  bool        s__gl_ptcorr;      ///< Switch for the use of the pt corrections: 1 = 1D or 2 = 2D TPofiles
  TFile*      s__gl_ptcorr_file = nullptr;  ///< File with 2D profiles for the pT corrections

//__________________________________________________________________________________
  
 // Event settings
  //float s__ev_PrimaryVertexZ;                   ///< Cut for the z-coordinate of the the primary vertex
  //std::vector<centrality_bin> s__ev_Centrality; ///< Vector to store the centrality bins
  //int n_centrality_bins;                        // The number of centrality bins


  // Track settings
  int   s__tr_NHits;              ///< Cut for the minimum number of hits
  float s__tr_NSigmaDCAx;         ///< N-Sigma cut for the DCAx
  float s__tr_NSigmaDCAy;         ///< N-Sigma cut for the DCAy
  float s__tr_NSigmaDCAz;         ///< N-Sigma cut for the DCAz
  float s__tr_LowPtCut;           ///< Cut for the minimum transverse momentum
  float s__tr_HighPtCut;          ///< Cut for the maximal transverse momentum
// pid settings
  float s__id_TPCSigma;           ///< N-Sigma cut for the PID by the TPC dE/dx infomation
  float s__id_TOFSigma;           ///< N-Sigma cut for the PID by the ToF mass-squared information
  float s__id_TOFDphiSigma;       ///< N-Sigma cut for the matching to TOF
  float s__id_TOFDzSigma;         ///< N-Sigma cut for the matching to TOF
// MpdPid settings
  float s__mid_Energy;            ///< Collision energy
  float s__mid_Coef;              ///< MpdPid dE/dx coefficient
  std::string s__mid_Generator;   ///< MpdPid generator string
  std::string s__mid_Tracking;    ///< MpdPid tracking string
  std::string s__mid_IniString;   ///< MpdPid initialization sting
  bool s__mid_dEdx;               ///< Switch for the use of the MpdPid "dE/dx only" case in addition to the combined PID
// Particle settings
  std::vector<particle_info> s__p_List;    ///< Vector to store the particles information

//__________________________________________________________________________________

  TClonesArray *mMCTracks                       = nullptr;
  TClonesArray *mKalmanTracks                   = nullptr;
  TClonesArray *mMpdGlobalTracks                = nullptr;
  TClonesArray *tofMatches                      = nullptr;

//__________________________________________________________________________________

  TH2F *h__dedx;
  TH2F *h__dedxHe3;
  TH2F *h__dedxHe4;
  TH2F *h__dedxt;
  TH2F *h__dedxd;
  TH2F *h__dedxap;
  TH2F *h__dedxp;
  TH2F *h__dedxkp;
  TH2F *h__dedxkm;
  TH2F *h__dedxpip;
  TH2F *h__dedxpim;
  TH2F *h__dedxpiparam;
  
//_____________________________________________________________________________________

  TH2F *h__m2;
  TH2F *h__m2He3;
  TH2F *h__m2He4;
  TH2F *h__m2t;
  TH2F *h__m2d;
  TH2F *h__m2ap;
  TH2F *h__m2p;
  TH2F *h__m2kp;
  TH2F *h__m2km;
  TH2F *h__m2pip;
  TH2F *h__m2pim;

  TF1 *fParDeBB;

 // int   particle_by_pdg(const int value);
  void  read_settings_json(const char* fname);
  
  ClassDef(EnerClass1,1);
};
#endif

