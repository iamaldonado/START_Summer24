#include "lowMgF.h"

ClassImp(lowMgF);

/// Default constructor
lowMgF::lowMgF(){
}

lowMgF::lowMgF(const char *name, const char *outputName) : MpdAnalysisTask(name, outputName){
}

void lowMgF::UserInit(){
  printf("lowMgF::Initialization\n");
  
  fOutputList = new TList();            // Standard list for the output histograms (MpdAnalysisTask)
  fOutputList -> SetOwner(kTRUE);       // Some default line (MpdAnalysisTask)
  TH1::AddDirectory(kFALSE);            // Some default line (MpdAnalysisTask)
  TH2::AddDirectory(kFALSE);            // Some default line (MpdAnalysisTask)

 
//_____All Histograms________________________________________________________________

//------Primary Vertex Position----------
 
   VertexPosition = new TH1F("VertexPosition0","Primary Vertex Position; Z_{Vertex} (cm); Entries", 200, -200, 200);

   DZvsZReco	  =  new TH2F("DZvsZReco", "#Delta Z vs Z_{Reco} ; Z_{Reco} (cm) ; #Delta Z (cm)", 500, -200, 200, 500, 0, 250);
   DZNtracks	  =  new TH2F("DZNtracks", "#Delta Z vs NTracks ; Number of Tracks; #Delta Z (cm)", 500, 0, 500, 500,0,250);
   DZb		  =  new TH2F("DZb", "#Delta Z vs b; b (fm) ; #Delta Z (cm)", 100, 0, 16, 100, 0, 250);

   TPDZvsZReco	  =  new TProfile("TPDZvsZReco", " TPtofile #Delta Z vs Z_{Reco} ; Z_{Reco} (cm) ; #Delta Z (cm)", 500, -200, 200, 0, 250);
   TPDZNtracks	  =  new TProfile("TPDZNtracks", "TProfile #Delta Z vs NTracks ; NTracks; #Detla Z (cm)", 500, 0, 500, 0, 250);
   TPDZb		  =  new TProfile("TPDZb", "TProfile #Delta Z vs b; b (fm) ; #Detla Z (cm)", 100, 0, 16, 0, 250);

   TPDZNtracksW	=  new TProfile2D("TPDZNtracksW", "Weight of Resolution of Z vs Number of Tracks; Number of Tracks; Z_{Reco} (cm) ; #Detla Z (cm)", 100, -4, 4, 100,0,5,0,4);
   TPDZbW	=  new TProfile2D("TPDZbW", " Weight of Resolution of Z vs Impact parameter; b (fm); Z_{Reco} (cm) ; #Delta Z (cm)", 100, -200, 200, 100,0, 16,0,4);

   VtxMult	=  new TProfile("VtxMult", " Vertex Resolution vs Track Multiplicity ; Track Multiplicity ; #Delta Z = |Z_{reco} - Z_{MC}| (%)  ", 150,0,800,0,15);

//------List Primary Vertex Position----- 
	
	fOutputList -> Add(VertexPosition);

	fOutputList -> Add(DZvsZReco);
	fOutputList -> Add(DZNtracks);
	fOutputList -> Add(DZb);
	
	fOutputList -> Add(TPDZvsZReco);
	fOutputList -> Add(TPDZNtracks);
	fOutputList -> Add(TPDZb);

	fOutputList -> Add(TPDZNtracksW);
	fOutputList -> Add(TPDZbW);

	fOutputList -> Add(VtxMult);

//------------Tranverse Momemtum----------------
   
   PtMCvsEta	=  new TH2F("PtMCvsEta", "p_{t}^{mc} vs #eta ; #eta ; p_{t}^{mc} (GeV/c)", 500, -3, 3, 500, 0, 4);

   PtRECOvsEta	=  new TH2F("ptRECOvsEta", "p_{t}^{reco} vs #eta; #eta ; p_{t}^{reco} (Gev/c)", 500, -3, 3, 500, 0, 4);
 
   DPtPtReco	=  new TProfile("DPtPtReco", " Transverse Momentum Resolution as a function of Transverse Momentum ; p_{T}^{Reco} (GeV/c) ; #Delta p_{T} = |#frac{p_{T}^{reco} - p_{T}^{MC}}{p_{T}^{MC}}| (%) ", 150,0,5,0,15);

	// Without Cuts
   PtNHits	=  new TProfile("DptvsNHits", " #Delta p_{T} vs NHits ; NumHits ; #Delta p_{T}", 100,0,50,0,15);
   PtEta	=  new TProfile("DptvsEta", "#Delta p_{T} vs #eta ; #eta ; #Delta p_{T}", 100, -3, 3, 0, 5);
   PtDCAGlobal	=  new TProfile("DptvsDCAGlobal", " #Delta p_{T} vs DCA Global ; DCA Global ; #Delta p_{T}", 150,0,5,0,15);
   PtDCAGlobalP	=  new TProfile("DptvsDCAGP", " #Delta p_{T} vs DCAGlobal Primary ; DCAGlobal ; #Delta p_{T}", 150,0,5,0,15);
   PtDCAGlobalS =  new TProfile("DptvsDCAGS", " #Delta p_{T} vs DCAGlobal Secondary; DCAGlobal ; #Delta p_{T}", 150,0,5,0,15);


   PtEtaDPt	=  new TProfile2D("PtEtaDPt", "#Delta p_{t};#eta; #Detla p_{t} (Gev/c)", 100, -4, 4, 100,0,5,0,4);
   
	// With Cuts in Number of Hits
   PtNHitsC	=  new TProfile("DptvsNHitsC", " #Delta p_{T} vs NHits with cut NHits > 27 ; NumHits ; #Delta p_{T}", 100,0,50,0,15);
   PtEtaC	=  new TProfile("DptvsEtaC", " #Delta p_{T} vs #eta with cut NHits > 27 ; #eta ; #Delta p_{T}", 100, -3, 3, 0, 5);
   PtDCAGlobalC	=  new TProfile("DptvsDCAGC", " #Delta p_{T} vs DCAGlobal with cut NHits > 16  ; DCAGlobal (cm) ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);
   PtDCAGlobalPC=  new TProfile("DptvsDCAGPC", " #Delta p_{T} vs DCAGlobal Primary with cut NHits > 16 ; DCAGlobal ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);
   PtDCAGlobalSC=  new TProfile("DptvsDCAGSC", " #Delta p_{T} vs DCAGlobal Secondary with cut NHits > 27 ; DCAGlobal ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);

   PtEtaDPtC	=  new TProfile2D("PtEtaDPtC","#Delta p_{T} with cut NHits > 16 ;#eta;p_{T} (GeV/c);#Delta p_{T}",100,-4,4,100,0,5,0,4);

	// With Cuts in Pseudo-Rapidity
   PtNHitsCE		=  new TProfile("DptvsNHitsCE", " #Delta p_{T} vs NHits wiht cut in #eta (-1.5, 1.5) ; NumHits ; #Delta p_{T} (GeV/c) ", 100,0,50,0,15);
   PtEtaCE		=  new TProfile("DptvsEtaCE", " #Delta p_{T} vs #eta with cut #eta (-1.5, 1.5) ; #eta ; #Delta p_{T} (GeV/c) ", 100, -3, 3, 0, 5);
   PtDCAGlobalCE	=  new TProfile("DptvsDCAGC", " #Delta p_{T} vs DCAGlobal with cut NHits > 16  ; DCAGlobal (cm) ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);
   PtDCAGlobalPCE	=  new TProfile("DptvsDCAGPC", " #Delta p_{T} vs DCAGlobal Primary with cut NHits > 16 ; DCAGlobal ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);
   PtDCAGlobalSCE	=  new TProfile("DptvsDCAGSC", " #Delta p_{T} vs DCAGlobal Secondary with cut NHits > 27 ; DCAGlobal ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);

   PtEtaDPtCE	=  new TProfile2D("PtEtaDPtCE","#Delta p_{T} with cut #eta (-1,1) ;#eta;p_{T} (GeV/c);#Delta p_{T}",100,-4,4,100,0,5,0,4);

	// With Cuts in DCA Global
   PtNHitsDCA		=  new TProfile("DptvsNHitsDCA", " #Delta p_{T} vs NHits wiht cut in #eta (-1.5, 1.5) ; NumHits ; #Delta p_{T} (GeV/c) ", 100,0,50,0,15);
   PtEtaDCA		=  new TProfile("DptvsEtaDCA", " #Delta p_{T} vs #eta with cut #eta (-1.5, 1.5) ; #eta ; #Delta p_{T} (GeV/c) ", 100, -3, 3, 0, 5);
   PtDCAGlobalDCA	=  new TProfile("DptvsDCAGDCA", " #Delta p_{T} vs DCAGlobal with cut NHits > 16  ; DCAGlobal (cm) ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);
   PtDCAGlobalPDCA	=  new TProfile("DptvsDCAGPDCA", " #Delta p_{T} vs DCAGlobal Primary with cut NHits > 16 ; DCAGlobal ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);
   PtDCAGlobalSDCA	=  new TProfile("DptvsDCAGSDCA", " #Delta p_{T} vs DCAGlobal Secondary with cut NHits > 27 ; DCAGlobal ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);

   PtEtaDPtDCA		=  new TProfile2D("PtEtaDPtDCA","#Delta p_{T} with cut #eta (-1,1) ;#eta;p_{T} (GeV/c);#Delta p_{T}",100,-4,4,100,0,5,0,4);

	// With Cuts in Transverse Momentum
   PtNHitsCPT		=  new TProfile("DptvsNHitsCPT", " #Delta p_{T} vs NHits wiht cut in p_{T} > 0.15 (Gev/c) ; NumHits ; #Delta p_{T} (GeV/c) ", 100,0,50,0,15);
   PtEtaCPT		=  new TProfile("DptvsEtaCPT", " #Delta p_{T} vs #eta with cut in p_{T} > 0.15 (GeV/c) ; #eta ; #Delta p_{T} (GeV/c) ", 100, -3, 3, 0, 5);
   PtDCAGlobalCPT	=  new TProfile("DptvsDCAGCPT", " #Delta p_{T} vs DCAGlobal with cut in p_{T} > 0.15 (GeV/c)  ; DCAGlobal (cm) ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);
   PtDCAGlobalPCPT	=  new TProfile("DptvsDCAGPCPT", " #Delta p_{T} vs DCAGlobal Primary with cut p_{T} > 0.15 (Gev/c) ; DCAGlobal (cm) ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);
   PtDCAGlobalSCPT	=  new TProfile("DptvsDCAGSCPT", " #Delta p_{T} vs DCAGlobal Secondary with cut p_{T} > 0.15 (GeV/c) ; DCAGlobal (cm) ; #Delta p_{T} (GeV/c) ", 150,0,5,0,15);

   PtEtaDPtCPT		=  new TProfile2D("PtEtaDPtCPT","#Delta p_{T} with cut p_{T} >0.15 (GeV/c) ;#eta ;p_{T} (GeV/c);#Delta p_{T}",100,-4,4,100,0,5,0,4);

//------------List Tranverse Momemtum-------------
	fOutputList -> Add(PtRECOvsEta);	// Monte-Carlo vs Pseudo-Rapidity

	fOutputList -> Add(PtMCvsEta);		// Reconstruted vs Pseudo-Rapidity

	fOutputList -> Add(DPtPtReco);		// Transverse Momentum as a function of Transverse Momentum

	// Without Cuts
	fOutputList -> Add(PtNHits);		// Number of Hits
	fOutputList -> Add(PtEta);		// Pseudo-Rapidity
	fOutputList -> Add(PtDCAGlobal);	// DCA Global
	fOutputList -> Add(PtDCAGlobalP);	// DCA Global Primary
	fOutputList -> Add(PtDCAGlobalS);	// DCA Global Secondary

	fOutputList -> Add(PtEtaDPt);		// Resolution 

	// With Cuts in NHits
	fOutputList -> Add(PtNHitsC);
	fOutputList -> Add(PtEtaC);
	fOutputList -> Add(PtDCAGlobalC);
	fOutputList -> Add(PtDCAGlobalPC);
	fOutputList -> Add(PtDCAGlobalSC);

	fOutputList -> Add(PtEtaDPtC);

	// With Cuts in Eta
	fOutputList -> Add(PtNHitsCE);
	fOutputList -> Add(PtEtaCE);
	fOutputList -> Add(PtDCAGlobalCE);	// DCA Global
	fOutputList -> Add(PtDCAGlobalPCE);	// DCA Global Primary
	fOutputList -> Add(PtDCAGlobalSCE);	// DCA Global Secondary

	fOutputList -> Add(PtEtaDPtCE);	

	// With Cuts in DCA Global
	fOutputList -> Add(PtNHitsDCA);		// Number of Hits
	fOutputList -> Add(PtEtaDCA);		// Pseudo-Rapidity
	fOutputList -> Add(PtDCAGlobalDCA);	// DCA Global
	fOutputList -> Add(PtDCAGlobalPDCA);	// DCA Global Primary
	fOutputList -> Add(PtDCAGlobalSDCA);	// DCA Global Secondary

	fOutputList -> Add(PtEtaDPtDCA);	// Resoltion

	// With Cuts in Transverse Momentum
	fOutputList -> Add(PtNHitsCPT);		// Number of Hits
	fOutputList -> Add(PtEtaCPT);		// Pseudo-Rapidity
	fOutputList -> Add(PtDCAGlobalCPT);	// DCA Global
	fOutputList -> Add(PtDCAGlobalPCPT);	// DCA Global Primary
	fOutputList -> Add(PtDCAGlobalSCPT);	// DCA Global Secondary

	fOutputList -> Add(PtEtaDPtCPT);	// Resolution

//--------DCA Global----------------------
   DCAGC = new TH1F("DCAGlobalC","DCA Global wiht cut in Number of Hits > 27 ; DCA Global (cm) ; Entries", 200, 0, 40);		// DCA Global
   DCAGPC = new TH1F("DCAGlobalPC","DCA Global Primary wiht cut in Number of Hits > 27 ; DCA Global (cm) ; Entries", 200, 0, 40);	// DCA Global Primary
   DCAGSC = new TH1F("DCAGlobalSC","DCA Global Secondary wiht cut in Number of Hits > 27 ; DCA Global (cm) ; Entries", 200, 0, 40);	// DCA Global Secondary

//-------List DCA Global------------------
	fOutputList -> Add(DCAGC);
	fOutputList -> Add(DCAGPC);
	fOutputList -> Add(DCAGSC);

//-------Track Efficiency------------------

   PtRecoPionP	= new TH1F("PtRecoPionP"," p_{T}^{RECO} #pi Primary ", 200, 0, 5);
   PtRecoPionS	= new TH1F("PtRecoPionS"," p_{T}^{RECO} #pi Secondary ", 200, 0, 5);
   PtMCPionP	= new TH1F("PtMCPionP"," p_{T}^{MC} #pi Primary ", 200, 0, 5);
   PtMCPionS	= new TH1F("PtMCPionS"," p_{T}^{MC} #pi Secundary ", 200, 0, 5);

   PtRecoProtonP= new TH1F("PtRecoProtonP"," p_{T}^{RECO} p Primary ", 200, 0, 5);
   PtRecoProtonS= new TH1F("PtRecoProtonS"," p_{T}^{RECO} p Secondary ", 200, 0, 5);
   PtMCProtonP  = new TH1F("PtMCProtonP"," p_{T}^{MC} p Primary ", 200, 0, 5);
   PtMCProtonS  = new TH1F("PtMCProtonS"," p_{T}^{MC} p Secundary", 200, 0, 5);

   PtRecoKaonP	= new TH1F("PtRecoKaonP"," p_{T}^{RECO} #kappa Primary ", 200, 0, 5);
   PtRecoKaonS	= new TH1F("PtRecoKaonS"," p_{T}^{RECO} #kappa Secondary ", 200, 0, 5);
   PtMCKaonP	= new TH1F("PtMCKaonP"," p_{T}^{MC} #kappa Primary ", 200, 0, 5);
   PtMCKaonS	= new TH1F("PtMCKaonS"," p_{T}^{MC} #kappa Secundary ", 200, 0, 5);

//-------List Track Efficiency-------------

	fOutputList -> Add(PtRecoPionP);
	fOutputList -> Add(PtRecoPionS);
	fOutputList -> Add(PtMCPionP);
	fOutputList -> Add(PtMCPionS);
   
	fOutputList -> Add(PtRecoProtonP);
	fOutputList -> Add(PtRecoProtonS);
	fOutputList -> Add(PtMCProtonP);
	fOutputList -> Add(PtMCProtonS);

	fOutputList -> Add(PtRecoKaonP);
	fOutputList -> Add(PtRecoKaonS);
	fOutputList -> Add(PtMCKaonP);
	fOutputList -> Add(PtMCKaonS);

}

void lowMgF::ProcessEvent(MpdAnalysisEvent &event){

  mMpdGlobalTracks	=	event.fMPDEvent->GetGlobalTracks();	// All needed branches are accessed in this code block:
  mKalmanTracks		=	event.fTPCKalmanTrack;			// MPD global tracks, TPC Kalman track
--------------------------------------------------------------------------------------------------------------------
  mMCTracks		=	event.fMCTrack;				// Monte-Carlo tracks
--------------------------------------------------------------------------------------------------------------------  
  tofMatches		=	event.fTOFMatching;			// ToF matching
--------------------------------------------------------------------------------------------------------------------
  vtxs			=	event.fVertex;				// Vertex Position
--------------------------------------------------------------------------------------------------------------------
  mMCEventHeader	=	event.fMCEventHeader;			// MC Event Header
   
	 //Primary Vertex point 
   Int_t nVert		= vtxs->GetEntriesFast();
   MpdVertex *vtx 	= (MpdVertex*) vtxs->First();

   Double_t ZReco	= vtx -> GetZ();
   Int_t nTVert		= vtx->GetNTracks();
   Double_t ZMC		= mMCEventHeader-> GetZ();
   Double_t b		= mMCEventHeader -> GetB();
   Double_t absZ	= TMath::Abs( ZReco );

   Double_t DZ		= TMath::Abs( ZReco - ZMC );

   TVector3 Prim_Vtx(vtx->GetX(),vtx->GetY(),vtx->GetZ()); 
 
	// Histograms about Primary Vertex
   VertexPosition -> Fill(ZReco);
// cout << "Dz: " << DZ << endl;	// This lines help us to see this variable.
// cout << "ZMC: " << ZMC << endl;    
   DZvsZReco	->	Fill(ZReco, DZ);
   DZNtracks	->	Fill(nTVert, DZ);
   DZb		->	Fill(b, DZ);

   TPDZvsZReco	->	Fill(ZReco, DZ);
   TPDZNtracks	->	Fill(nTVert, DZ);
   TPDZb	->	Fill(b, DZ);

   TPDZNtracksW ->	Fill(ZReco, nTVert, DZ);
   TPDZbW	->	Fill(ZReco, b,DZ);

//   if(TMath::Abs(ZReco) > -50) continue;

	// Multiplicity
   Int_t refMult;

   refMult = 0; 

//____Cuts______________________________________________________________________

  const Double_t CutNHits = 27;		// My cut is: 27

//-----------------------------------------------

  int ntrmc = mMCTracks -> GetEntries();
 // cout << "N of MC tracks = " << ntrmc << endl;
bool make_MC=1;
  if(make_MC){ 
  	
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
   MpdMCTrack *mctrack = (MpdMCTrack*)mMCTracks -> UncheckedAt(ID); 		// Monte Carlo track is open for reading	
   
    int   pdg		= mctrack -> GetPdgCode();                               // Track PDG code
    int   prodId	= mctrack -> GetMotherId();                              // Track primacy: -1 = primary, any other = secondary
    float pt_mc		= mctrack -> GetPt();                                    // Particle transverse momentum

//_____________________________________________________________________________________________________________________________________________

 	//Multiplicity
   refMult++;
   // Delta Pt
   Double_t DPt	= TMath::Abs( pt_reco - pt_mc ) / ( pt_mc );

   // DCA Global
   Double_t DCAG = TMath::Sqrt( pow(DCAX,2) + pow(DCAY,2) + pow(DCAZ,2));

   // Absolute PDG
   Int_t abspdg = TMath::Abs( pdg );
   
   // Fill Histograms
  	// Reco
   PtRECOvsEta	->	Fill(Eta, pt_reco);

	// No se todavia
   PtEtaDPt 	->	Fill(Eta, pt_reco,DPt);
   if(NHits < CutNHits)PtEtaDPtC 	->	Fill(Eta, pt_reco,DPt);
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
   	if(NHits < CutNHits)
        {
		PtDCAGlobalPC 	->	Fill(DCAG, DPt);
		DCAGPC 	->	Fill(DCAG);
		if(TMath::Abs(Eta) > 1.5)
		{
		}
        }
   }
   if(mctrack->GetMotherId()!=-1) // Secondary
   {
	PtDCAGlobalS	->	Fill(DCAG, DPt);
   	if(NHits < CutNHits)
	{
		PtDCAGlobalSC ->	Fill(DCAG, DPt);
		DCAGSC	->	Fill(DCAG);
	}
   }

   if(NHits < CutNHits)
   {
   	PtNHitsC	->	Fill(NHits, DPt);
   	PtDCAGlobalC	-> 	Fill(DCAG, DPt);
   	PtEtaC	-> 	Fill(Eta, DPt);
	if(TMath::Abs(Eta) > 1.5)
	{
	   	PtNHitsCE	->	Fill(NHits, DPt);
   		PtEtaCE	-> 	Fill(Eta, DPt);
	        PtEtaDPtCE 	->	Fill(Eta, pt_reco,DPt);

	}
        
   }
		// Cuts
  // if(TMath::Abs(pt_reco) < 1.5) continue; 
   if(NHits < CutNHits) continue;
   if(TMath::Abs(Eta) > 1.5) continue;
   if(DCAG > 1) continue;

   DPtPtReco	->	Fill(pt_reco,DPt);

	// Track Efficiency
   if(mctrack -> GetMotherId() ==-1 ) //Primarias
   {
	if(abspdg == 211)  PtRecoPionP	->	Fill(pt_reco);
	if(abspdg == 2212) PtRecoProtonP->	Fill(pt_reco);
	if(abspdg == 321)  PtRecoKaonP	->	Fill(pt_reco);
   }
   if(mctrack -> GetMotherId() !=-1 ) // Secondary Particles
   {
	if(pdg == 211)  PtRecoPionS	->	Fill(pt_reco);
	if(pdg == 2212) PtRecoProtonS	->	Fill(pt_reco);
	if(pdg == 321)  PtRecoKaonS	->	Fill(pt_reco);
   }


    
  } // Close the first loop.
 	// Multiplicity
  VtxMult -> Fill(refMult,DZ);
 
 
  Int_t nmctracks = mMCTracks->GetEntriesFast();

  for (int i=0; i<nmctracks; i++){
  auto mctrack = (MpdMCTrack*) mMCTracks->UncheckedAt(i); 

    int   pdg		= mctrack -> GetPdgCode();                              // Track PDG code
    float pt_mc		= mctrack -> GetPt();                                   // Particle transverse momentum
    TVector3 P(mctrack->GetPx(),mctrack->GetPy(),mctrack->GetPz());		// Definition of the vector of momentum
    Double_t Eta=0.5*TMath::Log((P.Mag() + mctrack->GetPz())/(P.Mag() - mctrack->GetPz()+1.e-13)); // Calculate of Pseudo-Rapidity
  
     // Absolute PDG
   Int_t abspdg = TMath::Abs( pdg );

  // if(TMath::Abs(pt_mc) < 1.5) continue; 
   if(TMath::Abs(Eta) > 1.5) continue;		// Cut on the Pseudo-Rapidity
 
	// Track Efficiency
   if(mctrack -> GetMotherId() ==-1 ) //Primarias
   {
	if(abspdg == 211)  PtMCPionP	->	Fill(pt_mc);
	if(abspdg == 2212) PtMCProtonP	->	Fill(pt_mc);
	if(abspdg == 321)  PtMCKaonP	->	Fill(pt_mc);
   }
   if(mctrack -> GetMotherId() !=-1 ) // Secondary Particles
   {
	if(pdg == 211)  PtMCPionS	->	Fill(pt_mc);
	if(pdg == 2212)  PtMCProtonS	->	Fill(pt_mc);
	if(pdg == 321)  PtMCKaonS	->	Fill(pt_mc);
   }
 } // Close the second loop.
}

/// Default destructor
lowMgF::~lowMgF(){}

void lowMgF::Finish(){

}
