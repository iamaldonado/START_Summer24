#include "EnerClass1.h"
#include "TMath.h"
ClassImp(EnerClass1);

/// Default constructor
EnerClass1::EnerClass1(){
}

EnerClass1::EnerClass1(const char *name, const char *outputName, const char* settings_filename ) : MpdAnalysisTask(name, outputName){
  mParamConfig = outputName;
  settings_file = settings_filename;
}

void EnerClass1::UserInit(){
  printf("EnerClass1::Initialization\n");
  
  read_settings_json(settings_file);    // read settings file (JSON-formatted)
  
  fOutputList = new TList();            // Standard list for the output histograms (MpdAnalysisTask)
  fOutputList -> SetOwner(kTRUE);       // Some default line (MpdAnalysisTask)
  TH1::AddDirectory(kFALSE);            // Some default line (MpdAnalysisTask)
  TH2::AddDirectory(kFALSE);            // Some default line (MpdAnalysisTask)
  
  
  //Histograms of ionization energy loss distribution_reconstructed tracks - by MC recons
  
  h__dedx = new TH2F("h__dedx", "dEdx vs P for all particles;p*q, GeV/c;dE/dx, arb. units", 300, -1.5, 1.5, 600, 0.0, 60);
        fOutputList->Add(h__dedx);
        
  h__dedxHe3 = new TH2F("h_dedx_He^{3}", "dEdx vs P for the helium 3; p*q GeV/c; dE/dx arb.units",300, -1.5, 1.5, 600, 0.0, 60);
        fOutputList->Add(h__dedxHe3);
 
  h__dedxHe4 = new TH2F("h_dedx_He^{4}", "dEdx vs P for the helium 4; p*q GeV/c; dE/dx arb.units",300, -1.5, 1.5, 600, 0.0, 60);
        fOutputList->Add(h__dedxHe4);

  h__dedxt = new TH2F("h_dedx_t", "dEdx vs P for the tritium; p*q GeV/c ; dE/dx arb.units", 300, -1.5, 1.5, 600, 0.0, 60);
        fOutputList->Add(h__dedxt);

  h__dedxd = new TH2F("h_dedx_d", "dEdx vs P for the deuterium; p*q GeV/c; dE/dx arb.units",300, -1.5, 1.5, 600, 0.0, 60);
        fOutputList->Add(h__dedxd);

  h__dedxap = new TH2F("h_dedx_ap", "dEdx vs P for the antiproton; p*q GeV/c; dE/dx arb.units", 300, -1.5, 1.5, 600, 0.0, 60);
        fOutputList->Add(h__dedxap);   
        
  h__dedxp = new TH2F("h_dedx_p", "dEdx vs P for the proton; p*q GeV/c; dE/dx arb.units",300, -1.5, 1.5, 600, 0.0, 60);
        fOutputList->Add(h__dedxp); 
        
  h__dedxkp = new TH2F("h_dedx_kp", "dEdx vs P for the kaon+; p*q GeV/c; dE/dx arb.units",300, -1.5, 1.5, 600, 0.0, 60);
        fOutputList->Add(h__dedxkp);  

  h__dedxkm = new TH2F("h_dedx_km", "dEdx vs P for the kaon-; p*q GeV/c; dE/dx arb.units", 300, -1.5, 1.5, 600, 0.0, 60);
        fOutputList->Add(h__dedxkm); 
     
  h__dedxpip = new TH2F("h_dedx_pip", "dEdx vs P for the pion+; p*q GeV/c; dE/dx arb.units", 300, -1.5, 1.5, 600, 0.0, 60);
        fOutputList->Add(h__dedxpip);
        
  h__dedxpim = new TH2F("h_dedx_pim", "dEdx vs P for the pion-; p*q GeV/c; dE/dx arb.units", 300, -1.5, 1.5, 600, 0.0, 60);
        fOutputList->Add(h__dedxpim);
     
      
//______________________________________________________________________________________________________________________________
//Histograms for m2

  h__m2 = new TH2F("h__m2", "m^{2} vs P for all particles; p*q GeV/c; m^{2} GeV^{2}/c^{4}", 200, -4.0, 4.0, 200, -1.0, 2.0);
        fOutputList->Add(h__m2);

  h__m2He3 = new TH2F("h__m2He3", "m^{2} vs P for the Helium3; p*q GeV/c; m^{2} GeV^{2}/c^{4}", 200, -4.0, 4.0, 200, -1.0, 2.0);
        fOutputList->Add(h__m2He3);

  h__m2He4 = new TH2F("h__m2He4", "m^{2} vs P for the Helium4; p*q GeV/c; m^{2} GeV^{2}/c^{4}", 200, -4.0, 4.0, 200, -1.0, 2.0);
        fOutputList->Add(h__m2He4);

  h__m2t = new TH2F("h__m2t", "m^{2} vs P for the tritium; p*q GeV/c; m^{2} GeV^{2}/c^{4}", 200, -4.0, 4.0, 200, -1.0, 2.0);
        fOutputList->Add(h__m2t);

  h__m2d = new TH2F("h__m2d", "m^{2} vs P for the deuterium; p*q GeV/c; m^{2} GeV^{2}/c^{4}", 200, -4.0, 4.0, 200, -1.0, 2.0);
        fOutputList->Add(h__m2d);

  h__m2ap = new TH2F("h__m2ap", "m^{2} vs P for the antiproton; p*q GeV/c; m^{2} GeV^{2}/c^{4}", 200, -4.0, 4.0, 200, -1.0, 2.0);
        fOutputList->Add(h__m2ap);

  h__m2p = new TH2F("h__m2p", "m^{2} vs P for the proton; p*q GeV/c; m^{2} GeV^{2}/c^{4}", 200, -4.0, 4.0, 200, -1.0, 2.0);
        fOutputList->Add(h__m2p);

  h__m2kp = new TH2F("h__m2kp", "m^{2} vs P for the kaon+; p*q GeV/c; m^{2} GeV^{2}/c^{4}", 200, -4.0, 4.0, 200, -1.0, 2.0);
        fOutputList->Add(h__m2kp);

  h__m2km = new TH2F("h__m2km", "m^{2} vs P for the kaon-; p*q GeV/c; m^{2} GeV^{2}/c^{4}", 200, -4.0, 4.0, 200, -1.0, 2.0);
        fOutputList->Add(h__m2km);

  h__m2pip = new TH2F("h__m2pip", "m^{2} vs P for the pion+; p*q GeV/c; m^{2} GeV^{2}/c^{4}", 200, -4.0, 4.0, 200, -1.0, 2.0);
        fOutputList->Add(h__m2pip);
  
  h__m2pim = new TH2F("h__m2pim", "m^{2} vs P for the pion-; p*q GeV/c; m^{2} GeV^{2}/c^{4}", 200, -4.0, 4.0, 200, -1.0, 2.0);
        fOutputList->Add(h__m2pim);
  
//________________________________________________________________________________________________________________________________________________      
/*  TF1 *fParDeBB = new TF1("fParDeBB", "([0] / TMath::Power(x / TMath::Sqrt(x * x + 3.52), [3]))*(([1] - TMath::Power(x / TMath::Sqrt(x * x + 3.52), [3]))-(TMath::Log([2] + TMath::Power(1.0 / (x / 1.876), [4]))))",0,5);
  fParDeBB->SetParameters(s__mid_Coef * 3.27e-07, 3.74, -0.23, 2.32, 0.987);
  h__dedxd->Fit(fParDeBB,"R");
  h__dedxd->Draw("same");
  fParDeBB -> Draw("same");
  fOutputList->Add(fParDeBB);  

  TF1 *fParTrBB = new TF1("fParTrBB", "([0] / TMath::Power(x / TMath::Sqrt(x * x + 7.89), [3]))*(([1] - TMath::Power(x / TMath::Sqrt(x * x + 7.89), [3]))-(TMath::Log([2] + TMath::Power(1.0 / (x / 2.81), [4]))))",0,5);
  fParTrBB->SetParameters(s__mid_Coef * 2.59e-07, 5.06, 0.0001, 2.2, 1.056);
  h__dedxt->Fit(fParTrBB,"R");
  h__dedxt->Draw("same");
  fParTrBB->Draw("same");
  fOutputList->Add(fParTrBB);
 
  TF1 *fParPrBB = new TF1("fParPrBB", "([0] / TMath::Power(x / TMath::Sqrt(x * x + 0.88), [3]))*(([1] - TMath::Power(x / TMath::Sqrt(x * x + 0.88), [3]))-(TMath::Log([2] + TMath::Power(1.0 / (x / 0.9383), [4]))))",0,5);
  fParPrBB->SetParameters(s__mid_Coef * 4.40008e-07, 2.97563, -0.192657, 2.16118, 0.61451);
  h__dedxp->Fit(fParPrBB,"R");
  h__dedxp->Draw("same");
  fParPrBB -> Draw("same");
  fOutputList->Add(fParPrBB);

  TF1 *fParHe3BB = new TF1("fParHe3BB", "[0] * ((1.0 + TMath::Power(x / 1.4047, 2.0)) / (TMath::Power(x / 1.4047, [3])) * ([1] + [2] * TMath::Log(1.0 + TMath::Power(x / 1.4047, 2.0))) - 1.0)",0,5);
  fParHe3BB->SetParameters(s__mid_Coef * 2.86201e-06, 2.10168, 2.74807e-01, 1.86774);
  h__dedxHe3->Fit(fParHe3BB,"R");
  h__dedxHe3->Draw("same");
  fParHe3BB -> Draw("same");
  fOutputList->Add(fParHe3BB);

  TF1 *fParHe4BB = new TF1("fParHe4BB", "[0] * ((1.0 + TMath::Power(x / 1.863, 2.0)) / (TMath::Power(x / 1.863, [3])) *([1] + [2] * TMath::Log(1.0 + TMath::Power(x / 1.863, 2.0))) - 1.0)",0,5);
  fParHe4BB->SetParameters(s__mid_Coef * 2.96e-06, 2.085, 0.256, 1.85);
  h__dedxHe4->Fit(fParHe4BB,"R");
  h__dedxHe4->Draw("same");
  fParHe4BB -> Draw("same");
  fOutputList->Add(fParHe4BB);

  TF1 *fParPiBB = new TF1("fParPiBB", "([0] / TMath::Power(x / TMath::Sqrt(x * x + 0.01949), [3])) * (([1] - TMath::Power(x / TMath::Sqrt(x * x + 0.01949), [3])) - (TMath::Log([2] + TMath::Power(1.0 / (x / 0.1396), [4]))))",0,5);
  fParPiBB->SetParameters(s__mid_Coef * (-1.19342e-07), -7.83114, 8.17749, 1.85775, -1.80695); 
  h__dedxpip->Fit(fParPiBB,"R");
  h__dedxpip->Draw("same");  
  fParPiBB -> Draw("same");
  fOutputList->Add(fParPiBB);

  TF1 *fParKaBB = new TF1("fParKaBB", "([0] / TMath::Power(x / TMath::Sqrt(x * x + 0.2437), [3])) * (([1] - TMath::Power(x / TMath::Sqrt(x * x + 0.2437), [3])) - (TMath::Log([2] + TMath::Power(1.0 / (x / 0.4937), [4]))))",0,5);
  fParKaBB->SetParameters(s__mid_Coef * 6.50167e-07, 1.01718, -0.795357, 1.80916, 0.0707667);
  h__dedxkp->Fit(fParKaBB,"R");
  h__dedxkp->Draw("same");
  fParKaBB->Draw("same");
  fOutputList->Add(fParKaBB);
*/
  
}

void EnerClass1::ProcessEvent(MpdAnalysisEvent &event){

  mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();                   // All needed branches are accessed in this code block:
  mKalmanTracks    =              event.fTPCKalmanTrack;                   // MPD global tracks, TPC Kalman tracks,
//--------------------------------------------------------------------------------------------------------------------
  mMCTracks        =                     event.fMCTrack;                   // Monte-Carlo tracks
//--------------------------------------------------------------------------------------------------------------------  
  tofMatches       =                 event.fTOFMatching;                   // ToF matching                                                            
//--------------------------------------------------------------------------------------------------------------------


 
  int ntrmc = mMCTracks -> GetEntries();
  cout << "N of MC tracks = " << ntrmc << endl;

  if(s__gl_MC){ 
  	
    for (long int i = 0; i < ntrmc; i++) {
    
  	MpdMCTrack* mctrack = (MpdMCTrack*) mMCTracks -> At(i); 
  	if(mctrack -> GetMotherId() != -1) continue;                          
  	// If the MpdMCTrack::GetMotherId() == -1  -- the particle is primary
 	//int current_particle_mc = particle_by_pdg(mctrack -> GetPdgCode());
 	float pt_mc               = mctrack -> GetPt();
 	float pt_tot              = mctrack -> GetPt();
    }
  }
//end of track loop monte carlo

  // MpdPid *pid = new MpdPid(s__id_TOFSigma, s__id_TPCSigma, s__mid_Energy, s__mid_Coef, s__mid_Generator, s__mid_Tracking, s__mid_IniString);
  // MpdPidQA *pidQA = new MpdPidQA(s__id_TOFSigma, s__id_TPCSigma, s__mid_Energy, s__mid_Coef, s__mid_Generator, s__mid_Tracking, s__mid_IniString);
   
   const Double_t AbsEtaMax = 1.3;
  // TString outpath = "/lhep/users/alejandrosj/new_mpdroot/mpdroot-v24.06.24/physics/simplept/macros";
   Double_t Theta,AbsEta;
   Int_t nHits;
   Int_t nMaxEvents = 0;
   const Double_t ImpParMax = -1.0;
   const Double_t VzMax = 50.0;
   const Int_t RECOmID = 0;
  

  //if ( ImpParMax != -1.0 ) { if ( fmcHeader->GetB() > ImpParMax ) continue; }
  //if ( VzMax != -1.0 ) { if ( TMath::Abs( fmcHeader->GetZ() ) > VzMax ) continue; }
 
// The main reco loop
  int ntr   = mKalmanTracks -> GetEntries();
  for (long int i = 0; i < ntr; i++) {
    MpdTrack          *mpdtrack = (MpdTrack *) mMpdGlobalTracks->UncheckedAt(i);       // Global track is accessed

//________________________________________________________________________________________________________________________________________
    MpdTpcKalmanTrack *kftrack  = (MpdTpcKalmanTrack *) mKalmanTracks->UncheckedAt(i); // The corresponding TPC Kalman track is also used
    int   kfcharge = kftrack -> Charge();                                              // for the charge
    double p       =  kftrack -> Momentum3().Mag() * kfcharge;                         // full momentum and
    double dedx    =  kftrack -> GetDedx();                                            // dE/dx information
    int mcId = kftrack -> GetTrackID();                                                // MpdTpcKalmanTrack::GetTrackID() gives the ID of the corresponding Monte Carlo track
//____________________________________________________________________________________________________________________________________________
    
    MpdMCTrack* mctrack = (MpdMCTrack*) mMCTracks -> At(mcId);                         // Monte Carlo track is open for reading  
    
    int   pdg                 = mctrack -> GetPdgCode();                               // Track PDG code
    int   prodId              = mctrack -> GetMotherId();                              // Track primacy: -1 = primary, any other = secondary
   // int   current_particle_mc = particle_by_pdg(pdg);                                  // Particle position in the particles vector (read in the settings file)
    float rapidity_mc         = mctrack -> GetRapidity();                              // Particle rapidity (CAN BE WRONG!!!)
    float pt_mc               = mctrack -> GetPt();                                    // Particle transverse momentum
    float p_mc                = mctrack -> GetP();                                     // Particle full momentum
    float pz_mc               = mctrack -> GetPz();                                    // Particle momentum z-component

//_____________________________________________________________________________________________________________________________________________
//_____________________Cuts____________________________________________________		
	Theta = TMath::PiOver2() - kftrack->GetParam(3);
	AbsEta = TMath::Abs( -TMath::Log(TMath::Tan(0.5*Theta)) );
	nHits = kftrack->GetNofHits();
//_______________________________________________________________________________
	
	if(RECOmID == 1 && mctrack -> GetMotherId() != -1) continue; //primary tracks
	if(RECOmID == 2 && mctrack -> GetMotherId() == -1) continue; //secondary tracks
	
        //if ( AbsEta > AbsEtaMax ) continue;
	//if ( nHits < s__tr_NHits ) continue;
//_______________________________________________________________________________
  
        //if ( ImpParMax != -1.0 ) { if ( fmcHeader->GetB() > ImpParMax ) continue; }
        //if ( VzMax != -1.0 ) { if ( TMath::Abs( fmcHeader->GetZ() ) > VzMax ) continue; }
    
    //pidQA->FillDedxHists(p, dedx, pdg);
    
    h__dedx -> Fill(p, dedx);
    h__m2 -> Fill(p, mpdtrack -> GetTofMass2()); 
    
    if (pdg == 1000020030)
    {
    	//He3
    	h__dedxHe3->Fill(p, dedx);
    	h__m2He3 -> Fill(p, mpdtrack -> GetTofMass2());
    }else if (pdg == 1000020040)
    {
   	//He4
   	h__dedxHe4->Fill(p, dedx); 
   	h__m2He4 -> Fill(p, mpdtrack -> GetTofMass2());
    }else if (pdg == 1000010030)
    {
   	//Tritio(t)
   	h__dedxt->Fill(p, dedx); 
   	h__m2t -> Fill(p, mpdtrack -> GetTofMass2()); 
    }else if (pdg == 1000010020)
    {
   	//Deuterion (d)
   	h__dedxd->Fill(p, dedx);
   	h__m2d -> Fill(p, mpdtrack -> GetTofMass2());  
   }else if (pdg == -2212)
    {
   	//Anti-proton (ap)
   	h__dedxap->Fill(p, dedx); 
   	h__m2ap -> Fill(p, mpdtrack -> GetTofMass2());
    }else if (pdg == 2212)
    {
   	//Proton (p)
   	h__dedxp->Fill(p, dedx); 
   	h__m2p -> Fill(p, mpdtrack -> GetTofMass2());
    }else if (pdg == 321)
    {
   	//Kaonplus (kp)
   	h__dedxkp->Fill(p, dedx);
   	h__m2kp -> Fill(p, mpdtrack -> GetTofMass2()); 
    }else if (pdg == -321)
    {
    	//Kaonminus (km)
   	h__dedxkm->Fill(p, dedx); 
   	h__m2km -> Fill(p, mpdtrack -> GetTofMass2());
    }else if (pdg == 211)
    {
    	//pionplus (pip)
   	h__dedxpip->Fill(p, dedx); 
   	h__m2pip -> Fill(p, mpdtrack -> GetTofMass2());
   	
    }else if (pdg == -211)
    {
    	//pionminus (pip)
   	h__dedxpim->Fill(p, dedx); 
   	h__m2pim -> Fill(p, mpdtrack -> GetTofMass2());
    }
    
  }
 
 // pidQA->GetDedxQA(outpath);
}

/// Default destructor
EnerClass1::~EnerClass1(){}

void EnerClass1::Finish(){
  if(s__gl_ptcorr) s__gl_ptcorr_file -> Close(); // Close the settings file
}

//__________________________________________________________________________________
    MpdHelix MakeHelix(const MpdKalmanTrack *tr) {
	const Double_t F_CUR0 = 0.3 * 0.01 * 5 / 10;
    	Double_t r = tr->GetPosNew();
    	Double_t phi = tr->GetParam(0) / r;
    	Double_t x = r * TMath::Cos(phi);
    	Double_t y = r * TMath::Sin(phi);
    	Double_t dip = tr->GetParam(3);
    	Double_t cur = F_CUR0 * TMath::Abs (tr->GetParam(4));
    	TVector3 o(x, y, tr->GetParam(1));
    	Int_t h = (Int_t) TMath::Sign(1.1,tr->GetParam(4));
    	MpdHelix helix(cur, dip, tr->GetParam(2)-TMath::PiOver2()*h, o, h);
    	return helix;
    	}

//___________________________________________________________________________________________________________
/**
  This subroutine reads the parameters from the
    JSON-formatted settings file
  \param fname Configuration file
*/
void EnerClass1::read_settings_json(const char* fname){
  namespace settings = boost::property_tree;
  settings::ptree s_tree;
// Open settings file or terminate the program if it is not found
  try{
    settings::read_json(fname, s_tree);
  }
  catch(std::exception & e){
    printf("EnerCLass1::Fail:: %s\n", e.what());
    throw;
  }

// Global settings
  s__gl_Verbose    =              s_tree.get<bool>("Verbose", 1);
  s__gl_nMpdPID    = s_tree.get<short>("N_MPD_PID_Particles", 8);
  s__gl_MC         =              s_tree.get<bool>("make_MC", 1);
  s__gl_Efficiency =      s_tree.get<bool>("make_Efficiency", 1);
  s__gl_PID        =            s_tree.get<short>("PID_mode", 2);
  s__gl_DCA        =            s_tree.get<short>("DCA_mode", 0);
  s__gl_TOF        =            s_tree.get<short>("TOF_mode", 0);
  s__gl_ptcorr     =   s_tree.get<bool>("use_pt_corrections", 0);
  if(s__gl_ptcorr){
    std::string fcorrname = s_tree.get<std::string>("pt_corrections_file", "pt_corrections.root");
    s__gl_ptcorr_file = new TFile(fcorrname.c_str(), "READ");
    if(!s__gl_ptcorr_file || s__gl_ptcorr_file -> IsZombie()){
      printf("Error opening file: %s\n", fcorrname.c_str());
      exit(-1);
    }
  }

  if(!s__gl_MC && s__gl_Efficiency){
    printf("The efficiencies calculation can be enabled only if the 'make_MC' option is enabled\n");
    printf("Disabling the efficiencies calculation\n");
    s__gl_Efficiency = 0;
  }

  if((s__gl_PID < 0 || s__gl_PID > 2) && s__gl_Efficiency){
    printf("The efficiencies calculation can be enabled only if PID mode is within 0..2\n");
    printf("Disabling the efficiencies calculation\n");
    s__gl_Efficiency = 0;
  }

// MpdPid settings
  if(s__gl_PID == 2){
    if(!s_tree.get_child_optional("MpdPid")){
      printf("MpdPid settings are not defined\n");
      exit(-1);
    }
    s__mid_Energy     = s_tree.get<float>("MpdPid.Energy");
    s__mid_Coef       = s_tree.get<float>("MpdPid.Coef");
    s__mid_Generator  = s_tree.get<std::string>("MpdPid.Generator");
    s__mid_Tracking   = s_tree.get<std::string>("MpdPid.Tracking");
    s__mid_IniString  = s_tree.get<std::string>("MpdPid.IniString");
    s__mid_dEdx       = s_tree.get<bool>("MpdPid.Add_dedx_only");
  }


// Event settings
  if(!s_tree.get_child_optional("Events")){
    printf("Event cuts are not defined\n");
    exit(-1);
  }
  //s__ev_PrimaryVertexZ = s_tree.get<double>("Events.PrimaryVertexZ");
  //for(auto& centrality_array : s_tree.get_child("Events.Centrality")){
   // float low, high;
    //float* const elements[2] = {&low, &high};
    //auto element = std::begin(elements);
    //for (auto& node : centrality_array.second){
      //**element++ = node.second.get_value<double>();
      //if (element == std::end(elements)) break;
   // }
   // s__ev_Centrality.push_back({low, high});
  //}

// Track settings
  if(!s_tree.get_child_optional("Tracks")){
    printf("Track cuts are not defined\n");
    exit(-1);
  }
  s__tr_NHits          = s_tree.get<int>("Tracks.NHits");
  s__tr_NSigmaDCAx     = s_tree.get<float>("Tracks.NSigmaDCAx");
  s__tr_NSigmaDCAy     = s_tree.get<float>("Tracks.NSigmaDCAy");
  s__tr_NSigmaDCAz     = s_tree.get<float>("Tracks.NSigmaDCAz");
  s__tr_LowPtCut       = s_tree.get<float>("Tracks.LowPtCut");
  s__tr_HighPtCut      = s_tree.get<float>("Tracks.HighPtCut");

// PID settings
  if(!s_tree.get_child_optional("PID")){
    printf("PID cuts are not defined\n");
    exit(-1);
  }
  s__id_TPCSigma     = s_tree.get<float>("PID.TPCSigma");
  s__id_TOFSigma     = s_tree.get<float>("PID.TOFSigma");
  s__id_TOFDphiSigma = s_tree.get<float>("PID.TOFDphiSigma");
  s__id_TOFDzSigma   = s_tree.get<float>("PID.TOFDzSigma");


// Particle settings
  if(!s_tree.get_child_optional("Particles")){
    printf("Particles are not defined\n");
    exit(-1);
  }
  for(auto& particles_array : s_tree.get_child("Particles")){
    particle_info p;
    p.name              = particles_array.first.c_str();
    p.pdg               = particles_array.second.get_child("PDG").get_value<int>();
    p.mass              = particles_array.second.get_child("Mass").get_value<float>();
    p.charge            = particles_array.second.get_child("Charge").get_value<int>();
    p.enum_position     = particles_array.second.get_child("Enum").get_value<short>();

    auto node_pt = particles_array.second.get_child("pt_bins");
    for (auto it = node_pt.begin(); it != node_pt.end(); ++it){
      int position = std::distance(node_pt.begin(), it);
      p.pt_bins[position] = it->second.get_value<float>();
    }

    auto node_y = particles_array.second.get_child("rapidity_bins");
    for (auto it = node_y.begin(); it != node_y.end(); ++it){
      int position = std::distance(node_y.begin(), it);
      p.rapidity_bins[position] = it->second.get_value<float>();
    }

    s__p_List.push_back(p);
  }
} // EnerClass1::read_settings_json()




