#include "lowMgF.h"

ClassImp(lowMgF);

/// Default constructor
lowMgF::lowMgF(){
}

lowMgF::lowMgF(const char *name, const char *outputName, const char* settings_filename ) : MpdAnalysisTask(name, outputName){
  mParamConfig = outputName;
  settings_file = settings_filename;
}

void lowMgF::UserInit(){
  printf("lowMgF::Initialization\n");
  
  read_settings_json(settings_file);    // read settings file (JSON-formatted)
  
  if(s__gl_PID == 2){                   // Initialization of the MpdPid
    fPID = new MpdPid(s__id_TOFSigma, s__id_TPCSigma, s__mid_Energy, s__mid_Coef, s__mid_Generator.c_str(), s__mid_Tracking.c_str(), s__mid_IniString.c_str());
  }
  
  fOutputList = new TList();            // Standard list for the output histograms (MpdAnalysisTask)
  fOutputList -> SetOwner(kTRUE);       // Some default line (MpdAnalysisTask)
  TH1::AddDirectory(kFALSE);            // Some default line (MpdAnalysisTask)
  TH2::AddDirectory(kFALSE);            // Some default line (MpdAnalysisTask)
  
   // All Histograms
 
	// Primary Vertex Position 
   VertexPosition = new TH1F("VertexPosition0","Vertex Position; Z_{Vertex} (cm); Entries", 200, -200, 200);
	fOutputList -> Add(VertexPosition);

	// Primary Vertez Resolution
   DZvsZReco	  =  new TH2F("DZvsZReco", "#Delta Z vs Z_{Reco} ; Z_{Reco} ; #Delta (cm)", 500, -200, 200, 500, 0, 250);
	fOutputList -> Add(DZvsZReco);
   TPDZvsZReco	  =  new TProfile("TPDZvsZReco", "#Delta Z vs Z_{Reco} ; Z_{Reco} ; #Delta (cm)", 500, -200, 200, 0, 250);
	fOutputList -> Add(TPDZvsZReco);


   DZNtracks	  =  new TH2F("DZNtracks", "#Delta Z vs NTracks ; NTracks; #Delta Z (cm)", 500, 0, 500, 500,0,250);
	fOutputList -> Add(DZNtracks);
   TPDZNtracks	  =  new TProfile("TPDZNtracks", "#Delta Z vs NTracks ; NTracks; #Detla Z (cm)", 500, 0, 500, 0, 250);
	fOutputList -> Add(TPDZNtracks);


   DZb		  =  new TH2F("DZb", "#Delta Z vs b; b ; #Delta Z (cm)", 100, 0, 16, 100, 0, 250);
	fOutputList -> Add(DZb);
   TPDZb		  =  new TProfile("TPDZb", "#Delta Z vs b; b ; #Detla Z (cm)", 100, 0, 16, 0, 250);
	fOutputList -> Add(TPDZb);


   TPDZNtracksW	=  new TProfile2D("TPDZNtracksW", "#Delta p_{t};#eta; #Detla p_{t} (Gev/c)", 100, -4, 4, 100,0,5,0,4);
	fOutputList -> Add(TPDZNtracksW);
   TPDZbW	=  new TProfile2D("TPDZbW", "#Delta p_{t};#eta; Z", 100, -4, 4, 100,0,5,0,4);
	fOutputList -> Add(TPDZbW);


	// Monte Carlo
   PtMCvsEta	=  new TH2F("PtMCvsEta", "p_{t}^{mc} vs #eta ; #eta ; p_{t}^{mc} (GeV/c)", 500, -3, 3, 500, 0, 4);
	fOutputList -> Add(PtMCvsEta);
 
	// Reco
   PtRECOvsEta	=  new TH2F("ptRECOvsEta", "p_{t}^{reco} vs #eta; #eta (cm) ; p_{t}^{reco} (Gev/c)", 500, -3, 3, 500, 0, 4);
	fOutputList -> Add(PtRECOvsEta);

	// Without Cuts
		// Parameters
   PtNHits	=  new TProfile("DptvsNHits", " #Delta p_{T} vs NHits ; NumHits ; #Delta p_{T}", 100,0,50,0,15);
	fOutputList -> Add(PtNHits);

   PtEta	=  new TProfile("DptvsEta", "#Delta p_{T} vs #eta ; #eta ; #Delta p_{T}", 100, -3, 3, 0, 5);
	fOutputList -> Add(PtEta);

			// DCA Global
   PtDCAGlobal	=  new TProfile("DptvsDCAGlobal", " #Delta p_{T} vs DCA Global ; DCA Global ; #Delta p_{T}", 150,0,5,0,15);
	fOutputList -> Add(PtDCAGlobal);

   PtDCAGlobalP	=  new TProfile("DptvsDCAGP", " #Delta p_{T} vs DCAGlobal Primary ; DCAGlobal ; #Delta p_{T}", 150,0,5,0,15);
	fOutputList -> Add(PtDCAGlobalP);

   PtDCAGlobalS =  new TProfile("DptvsDCAGS", " #Delta p_{T} vs DCAGlobal Secondary; DCAGlobal ; #Delta p_{T}", 150,0,5,0,15);
	fOutputList -> Add(PtDCAGlobalS);

		// Resolution 
   PtEtaDPt	=  new TProfile2D("PtEtaDPt", "#Delta p_{t};#eta; #Detla p_{t} (Gev/c)", 100, -4, 4, 100,0,5,0,4);
	fOutputList -> Add(PtEtaDPt);

	// With Cuts in NHits
		// Parameters
   PtNHitsC	=  new TProfile("DptvsNHitsC", " #Delta p_{T} vs NHits with cut NHits > 27 ; NumHits ; #Delta p_{T}", 100,0,50,0,15);
	fOutputList -> Add(PtNHitsC);

   PtEtaC	=  new TProfile("DptvsEtaC", " #Delta p_{T} vs #eta with cut NHits > 27 ; #eta ; #Delta p_{T}", 100, -3, 3, 0, 5);
	fOutputList -> Add(PtEtaC);

	// With Cuts in Eta
		// Parameters
   PtNHitsCE	=  new TProfile("DptvsNHitsCE", " #Delta p_{T} vs NHits wiht cut in #eta (-1.5, 1.5) ; NumHits ; #Delta p_{T} (GeV/c) ", 100,0,50,0,15);
	fOutputList -> Add(PtNHitsCE);

   PtEtaCE	=  new TProfile("DptvsEtaCE", " #Delta p_{T} vs #eta with cut #eta (-1.5, 1.5) ; #eta ; #Delta p_{T} (GeV/c) ", 100, -3, 3, 0, 5);
	fOutputList -> Add(PtEtaCE);

		// DCA GLobal
   PtDCAGlobalC	=  new TProfile("DptvsDCAGC", " #Delta p_{T} vs DCAGlobal with cut NHits > 16  ; DCAGlobal (cm) ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);
	fOutputList -> Add(PtDCAGlobalC);

   PtDCAGlobalPC=  new TProfile("DptvsDCAGPC", " #Delta p_{T} vs DCAGlobal Primary with cut NHits > 16 ; DCAGlobal ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);
	fOutputList -> Add(PtDCAGlobalPC);

   PtDCAGlobalSC=  new TProfile("DptvsDCAGSC", " #Delta p_{T} vs DCAGlobal Secondary with cut NHits > 27 ; DCAGlobal ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);
	fOutputList -> Add(PtDCAGlobalSC);

	// Resolution with cut in Number of Hits
   PtEtaDPtC	=  new TProfile2D("PtEtaDPtC","#Delta p_{T} with cut NHits > 16 ;#eta;p_{T} (GeV/c);#Delta p_{T}",100,-4,4,100,0,5,0,4);
	fOutputList -> Add(PtEtaDPtC);
	
 	// Resolution with cut in Eta
   PtEtaDPtCE	=  new TProfile2D("PtEtaDPtCE","#Delta p_{T} with cut #eta (-1,1) ;#eta;p_{T} (GeV/c);#Delta p_{T}",100,-4,4,100,0,5,0,4);
	fOutputList -> Add(PtEtaDPtCE);
	

	// DCA Global
   DCAGC = new TH1F("DCAGlobalC","DCA Global wiht cut in Number of Hits > 27 ; DCA Global (cm) ; Entries", 200, 0, 40);
	fOutputList -> Add(DCAGC);

   DCAGPC = new TH1F("DCAGlobalPC","DCA Global Primary wiht cut in Number of Hits > 27 ; DCA Global (cm) ; Entries", 200, 0, 40);
	fOutputList -> Add(DCAGPC);

   DCAGSC = new TH1F("DCAGlobalSC","DCA Global Secondary wiht cut in Number of Hits > 27 ; DCA Global (cm) ; Entries", 200, 0, 40);
	fOutputList -> Add(DCAGSC);

	// Track Efficiency
   PtRecoPionP	= new TH1F("PtRecoPionP"," p_{T}^{RECO} #pi Primary ", 200, 0, 5);
	fOutputList -> Add(PtRecoPionP);
   PtRecoPionS	= new TH1F("PtRecoPionS"," p_{T}^{RECO} #pi Secondary ", 200, 0, 5);
	fOutputList -> Add(PtRecoPionS);
   PtMCPionP	= new TH1F("PtMCPionP"," p_{T}^{MC} #pi Primary ", 200, 0, 5);
	fOutputList -> Add(PtMCPionP);
   PtMCPionS	= new TH1F("PtMCPionS"," p_{T}^{MC} #pi Secundary ", 200, 0, 5);
	fOutputList -> Add(PtMCPionS);
   
   PtRecoProtonP= new TH1F("PtRecoProtonP"," p_{T}^{RECO} p Primary ", 200, 0, 5);
	fOutputList -> Add(PtRecoProtonP);
   PtRecoProtonS= new TH1F("PtRecoProtonS"," p_{T}^{RECO} p Secondary ", 200, 0, 5);
	fOutputList -> Add(PtRecoProtonS);
   PtMCProtonP  = new TH1F("PtMCProtonP"," p_{T}^{MC} p Primary ", 200, 0, 5);
	fOutputList -> Add(PtMCProtonP);
   PtMCProtonS  = new TH1F("PtMCProtonS"," p_{T}^{MC} p Secundary", 200, 0, 5);
	fOutputList -> Add(PtMCProtonS);

   PtRecoKaonP	= new TH1F("PtRecoKaonP"," p_{T}^{RECO} #kappa Primary ", 200, 0, 5);
	fOutputList -> Add(PtRecoKaonP);
   PtRecoKaonS	= new TH1F("PtRecoKaonS"," p_{T}^{RECO} #kappa Secondary ", 200, 0, 5);
	fOutputList -> Add(PtRecoKaonS);
   PtMCKaonP	= new TH1F("PtMCKaonP"," p_{T}^{MC} #kappa Primary ", 200, 0, 5);
	fOutputList -> Add(PtMCKaonP);
   PtMCKaonS	= new TH1F("PtMCKaonS"," p_{T}^{MC} #kappa Secundary ", 200, 0, 5);
	fOutputList -> Add(PtMCKaonS);


}

void lowMgF::ProcessEvent(MpdAnalysisEvent &event){

  mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();                   // All needed branches are accessed in this code block:
  mKalmanTracks    =              event.fTPCKalmanTrack;                   // MPD global tracks, TPC Kalman track
--------------------------------------------------------------------------------------------------------------------
  mMCTracks        =                     event.fMCTrack;                   // Monte-Carlo tracks
--------------------------------------------------------------------------------------------------------------------  
  tofMatches       =                 event.fTOFMatching;                   // ToF matching
--------------------------------------------------------------------------------------------------------------------
  vtxs		= 		event.fVertex;
--------------------------------------------------------------------------------------------------------------------
  mMCEventHeader = event.fMCEventHeader;


    //PrimaryVtx point 
   Int_t nVert = vtxs->GetEntriesFast();
   MpdVertex *vtx = (MpdVertex*) vtxs->First();

   Double_t ZReco = vtx -> GetZ();
   Int_t nTVert = vtx->GetNTracks();
   Double_t ZMC = mMCEventHeader-> GetZ();
   Double_t b = mMCEventHeader -> GetB();

   Double_t DZ	= TMath::Abs(( ZReco - ZMC ) / ( ZMC ));

   TVector3 Prim_Vtx(vtx->GetX(),vtx->GetY(),vtx->GetZ()); 
 
	// Histograms about Primary Vertex
   VertexPosition -> Fill(ZReco);
//cout << "Dz: " << DZ << endl;
// cout << "ZMC: " << ZMC << endl;    
   DZvsZReco	->	Fill(ZReco, DZ);
   DZNtracks	->	Fill(nTVert, DZ);
   DZb	->	Fill(b, DZ);

   TPDZvsZReco	->	Fill(ZReco, DZ);
   TPDZNtracks	->	Fill(nTVert, DZ);
   TPDZb	->	Fill(b, DZ);

   TPDZNtracksW ->	Fill(ZReco, nTVert, DZ);
   TPDZbW	->	Fill(ZReco, b,DZ);

//   if(TMath::Abs(posVtxZ) > 4) continue;
//-----------------------------------------------

  int ntrmc = mMCTracks -> GetEntries();
 // cout << "N of MC tracks = " << ntrmc << endl;

  if(s__gl_MC){ 
  	
    for (long int i = 0; i < ntrmc; i++) {
    
  	MpdMCTrack* mctrack = (MpdMCTrack*) mMCTracks -> At(i); 
  	if(mctrack -> GetMotherId() != -1) continue;                          
  	// If the MpdMCTrack::GetMotherId() == -1  -- the particle is primary
 	//int current_particle_mc = particle_by_pdg(mctrack -> GetPdgCode());
 	float pt_mc               = mctrack -> GetPt();
 	float pt_tot               = mctrack -> GetPt();
    }
  }
//end of track loop monte carlo

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
    
  // MpdTrack *track = (MpdTrack*) event.fMPDEvent->GetGlobalTracks()->UncheckedAt(i)
   MpdTrack *track  =  (MpdTrack*) mMpdGlobalTracks->UncheckedAt(i);
	// Variables
        int Ntrack  	=  track -> GetID();
	float pt_reco	=  track -> GetPt();    
        Double_t NHits 	=  track -> GetNofHits();
	Double_t Eta	=  track -> GetEta();
	Double_t DCAX	=  track -> GetDCAX();
	Double_t DCAY	=  track -> GetDCAY();
	Double_t DCAZ	=  track -> GetDCAZ();


//____________________________________________________________________________________________________________________________________________
   Int_t ID = track -> GetID(); 
   MpdMCTrack *mctrack = (MpdMCTrack*)mMCTracks -> UncheckedAt(ID); 
// MpdMCTrack* mctrack = (MpdMCTrack*) mMCTracks -> At(mcId);                         // Monte Carlo track is open for reading  
   
    int   pdg		= mctrack -> GetPdgCode();                               // Track PDG code
    int   prodId	= mctrack -> GetMotherId();                              // Track primacy: -1 = primary, any other = secondary
   // int   current_particle_mc = particle_by_pdg(pdg);                                  // Particle position in the particles vector (read in the settings file)
    float rapidity_mc	= mctrack -> GetRapidity();                              // Particle rapidity (CAN BE WRONG!!!)
    float p_mc		= mctrack -> GetP();                                     // Particle full momentum
    float pt_mc		= mctrack -> GetPt();                                    // Particle transverse momentum
    float pz_mc		= mctrack -> GetPz();                                    // Particle momentum z-component

//_____________________________________________________________________________________________________________________________________________

   // Delta Pt
   Double_t DPt	= TMath::Abs( pt_reco - pt_mc ) / ( pt_mc );

   // DCA Global
   Double_t DCAG = TMath::Sqrt( pow(DCAX,2) + pow(DCAY,2) + pow(DCAZ,2));

   // Fill Histograms
 
  	// Reco
   PtRECOvsEta	->	Fill(Eta, pt_reco);

	// No se todavia
   PtEtaDPt 	->	Fill(Eta, pt_reco,DPt);
   if(NHits>27)PtEtaDPtC 	->	Fill(Eta, pt_reco,DPt);
  	// MC
   PtMCvsEta	->	Fill(Eta, pt_mc);

	// Numbers of Hits
   PtNHits	-> 	Fill(NHits, DPt);
   PtDCAGlobal	-> 	Fill(DCAG, DPt);
   PtEta	-> 	Fill(Eta, DPt);

	// Numbers of Hits with cut

   DCAGC 	->	Fill(DCAG);

   // Obtain the Mother ID 
   if (!mctrack)continue; 

   if(mctrack->GetMotherId()==-1) //Primarias
   {
	PtDCAGlobalP 	->	Fill(DCAG, DPt);
   	if(NHits > 27)
        {
		PtDCAGlobalPC 	->	Fill(DCAG, DPt);
		DCAGPC 	->	Fill(DCAG);
		if(Eta > -1.5 && Eta < 1.5)
		{
		}
        }
   }
   if(mctrack->GetMotherId()!=-1) // Secondary
   {
	PtDCAGlobalS	->	Fill(DCAG, DPt);
   	if(NHits > 27)
	{
		PtDCAGlobalSC ->	Fill(DCAG, DPt);
		DCAGSC	->	Fill(DCAG);
	}
   }
   if(NHits > 27)
   {
   	PtNHitsC	->	Fill(NHits, DPt);
   	PtDCAGlobalC	-> 	Fill(DCAG, DPt);
   	PtEtaC	-> 	Fill(Eta, DPt);
	if(Eta > -1.5 && Eta < 1.5)
	{
	   	PtNHitsCE	->	Fill(NHits, DPt);
   		PtEtaCE	-> 	Fill(Eta, DPt);
	        PtEtaDPtCE 	->	Fill(Eta, pt_reco,DPt);

	}
        
   }
	// Track Efficiency
   if(mctrack->GetMotherId()==-1) //Primarias
   {
   	if(NHits > 27)
        {
		if(Eta > -1.5 && Eta < 1.5)
		{
			if(DCAG > 1)
			{
			if(pdg == 211)  PtRecoPionP	->	Fill(pt_reco);
			if(pdg == 2212) PtRecoProtonP	->	Fill(pt_reco);
			if(pdg == 321)  PtRecoKaonP	->	Fill(pt_reco);

			if(pdg == 211)  PtMCPionP	->	Fill(pt_mc);
			if(pdg == 2212) PtMCProtonP	->	Fill(pt_mc);
			if(pdg == 321)  PtMCKaonP	->	Fill(pt_mc);
			}
		}
        }
   }
   if(mctrack->GetMotherId() !=-1) // Secondary Particles
   {
   	if(NHits > 27)
        {
		if(Eta > -1.5 && Eta < 1.5)
		{
			if(DCAG > 1)
			{
			if(pdg == 211)  PtRecoPionS	->	Fill(pt_reco);
			if(pdg == 2212) PtRecoProtonS	->	Fill(pt_reco);
			if(pdg == 321)  PtRecoKaonS	->	Fill(pt_reco);

			if(pdg == 211)  PtMCPionS	->	Fill(pt_mc);
			if(pdg == 2212)  PtMCProtonS	->	Fill(pt_mc);
			if(pdg == 321)  PtMCKaonS	->	Fill(pt_mc);
			}
		}
        }
   }


    
  }
}

/// Default destructor
lowMgF::~lowMgF(){}

void lowMgF::Finish(){
  if(s__gl_ptcorr) s__gl_ptcorr_file -> Close(); // Close the settings file
}

//_____________________________________________

//int MpdNuclei::particle_by_pdg(const int value){
 // for(int particle = 0; particle < s__p_List.size(); ++particle){
 //   if (value == s__p_List[particle].pdg) return particle;
 // }
 // return -1;
//}
//_______________________________________________

//___________________________________________________________________________________________________________
/**
  This subroutine reads the parameters from the
    JSON-formatted settings file
  \param fname Configuration file
*/
void lowMgF::read_settings_json(const char* fname){
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
} // lowMgF::read_settings_json()




