#include "EnerClass1.h"
#include "TMath.h"
ClassImp(EnerClass1);

/// Default constructor
EnerClass1::EnerClass1(){
}

EnerClass1::EnerClass1(const char *name, const char *outputName) : MpdAnalysisTask(name, outputName){

}

void EnerClass1::UserInit(){
  printf("EnerClass1::Initialization\n");
  
  
  fOutputList = new TList();            // Standard list for the output histograms (MpdAnalysisTask)
  fOutputList -> SetOwner(kTRUE);       // Some default line (MpdAnalysisTask)
  TH1::AddDirectory(kFALSE);            // Some default line (MpdAnalysisTask)
  TH2::AddDirectory(kFALSE);            // Some default line (MpdAnalysisTask)
  
  
  //Histograms of ionization energy loss distribution_reconstructed tracks 
  
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
  
//Histograms for efficiency of pt_________________________________________________________________________________________________________________  
//for loop Montecarlo_____________________________________________
  PtProMC1 = new TH1F("PtProMC1","P_{T}^{MC}",200,0,5);
        fOutputList->Add(PtProMC1);
  PtKaonMC1 = new TH1F("PtKaonMC1","P_{T}^{MC}",200,0,5);
        fOutputList->Add(PtKaonMC1);
  PtPionMC1 = new TH1F("PtPionMC1","P_{T}^{MC}",200,0,5);
        fOutputList->Add(PtPionMC1);

//for loop reconstructed tracks___________________________________
//for protons________________________________________________
  PtProreco1 = new TH1F("PtProreco1","P_{T}^{Reco}",200,0,5); //for de/dx restrictions 
        fOutputList->Add(PtProreco1);
  PtProreco2 = new TH1F("PtProreco2","P_{T}^{Reco}",200,0,5); //for de/dx + PDG restrictions
        fOutputList->Add(PtProreco2);
//for kaons__________________________________________________
  PtKaonreco1 = new TH1F("PtKaonreco1","P_{T}^{Reco}",200,0,5); //for de/dx restrictions
        fOutputList->Add(PtKaonreco1);
  PtKaonreco2 = new TH1F("PtKaonreco2","P_{T}^{Reco}",200,0,5); //for de/dx + PDG restrictions
        fOutputList->Add(PtKaonreco2);
//for pions_________________________________________________
  PtPionreco1 = new TH1F("PtPionreco1","P_{T}^{Reco}",200,0,5); //for de/dx restrictions
        fOutputList->Add(PtPionreco1);
  PtPionreco2 = new TH1F("PtPionreco2","P_{T}^{Reco}",200,0,5); //fro de/dx + PDG restriccions
        fOutputList->Add(PtPionreco2);

//Histograms of energy loss using cuts in dedx + PDG_________________________________________________________

  h__dedxpcut = new TH2F("h_dedx_pcut", "dEdx vs P for the proton; p*q GeV/c; dE/dx arb.units",300, -1.5, 1.5, 600, 0.0, 60);
        fOutputList->Add(h__dedxpcut);

  h__dedxkpcut = new TH2F("h_dedx_kpcut", "dEdx vs P for the kaon+; p*q GeV/c; dE/dx arb.units",300, -1.5, 1.5, 600, 0.0, 60);
        fOutputList->Add(h__dedxkpcut);

  h__dedxpipcut = new TH2F("h_dedx_pipcut", "dEdx vs P for the pion+; p*q GeV/c; dE/dx arb.units", 300, -1.5, 1.5, 600, 0.0, 60);
        fOutputList->Add(h__dedxpipcut);
  
}

void EnerClass1::ProcessEvent(MpdAnalysisEvent &event){

//Initialization of the branches_________________________________________________________________________________________________
  mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();                   // All needed branches are accessed in this code block:
  mKalmanTracks    =              event.fTPCKalmanTrack;                   // MPD global tracks, TPC Kalman tracks,
//--------------------------------------------------------------------------------------------------------------------
  mMCTracks        =                     event.fMCTrack;                   // Monte-Carlo tracks
//--------------------------------------------------------------------------------------------------------------------  
  tofMatches       =                 event.fTOFMatching;                   // ToF matching                                                            
//--------------------------------------------------------------------------------------------------------------------

 

 int ntrmc = mMCTracks -> GetEntries();       //monte carlo counter limiter
 // cout << "N of MC tracks = " << ntrmc << endl;
 bool make_MC=1;
//Begin of the loop monte carlo track___________________________________________________________
  if(make_MC){ 
  	
    for (long int i = 0; i < ntrmc; i++) {
    
  	MpdMCTrack* mctrack = (MpdMCTrack*) mMCTracks -> At(i); 
  	if(mctrack -> GetMotherId() != -1) continue;//for primary particles                          
      
        int   pdg                 = mctrack -> GetPdgCode(); 
	Double_t pt_mc            = mctrack -> GetPt();
 	Double_t p_tot            = mctrack -> GetP();
        Double_t pz_mc            = mctrack -> GetPz();
        Double_t Eta_mc;

//______Cuts in monte carlo track_________________________________________
        float limptreco = 0.15;  
        const Double_t EtaMax = 1.5;     

	Eta_mc = 0.5*TMath::Log((p_tot + pz_mc)/(p_tot - pz_mc + 1e-16));

        if (TMath::Abs(pt_mc) < limptreco) continue;
        if (TMath::Abs(Eta_mc) > EtaMax) continue;     

        if (pdg == 2212){
	//proton
		PtProMC1->Fill(pt_mc);	
	}else if (pdg == 321){
	//kaon
        	PtKaonMC1->Fill(pt_mc);	
	}else if(pdg == 211){
	//pion
	        PtPionMC1->Fill(pt_mc);
	}
        
     }
  }
//end of track loop monte carlo____________________________________________________________


//Declaration of the general variables (limits for the cuts) fot the kalmantrack loop______________________________________________________________
   
   const Double_t EtaMax = 1.5;
   const Double_t EtaMin = -1.5;
   const Double_t limnHits = 27;
   const Double_t limDCAG = 1;
   float limptreco = 0.15;
   const Double_t ImpParMax = -1.0;
   const Double_t VzMax = 50.0;
   const Int_t RECOmIDP = 1;
   const Int_t RECOmIDS = 2;
  
 
//______________________________________________________________________________________________________________________________________

  int ntr   = mKalmanTracks -> GetEntries();//reconstructed tracks counter limiter
 
//Begin of the loop of kalmantrack recontructed tracks (The main reco loop)______________________________________________________________
  for (long int i = 0; i < ntr; i++) {
    
//__Access to the tracks_________________________________________________________________________________________________________________
    MpdTrack          *mpdtrack = (MpdTrack *) mMpdGlobalTracks->UncheckedAt(i);       // Global track is accessed
    MpdTpcKalmanTrack *kftrack  = (MpdTpcKalmanTrack *) mKalmanTracks->UncheckedAt(i); // The corresponding TPC Kalman track is also used 

//__specific variables for the GlobalTracks_____________________________________________________________________________________________
    
    Double_t Eta        = mpdtrack-> GetEta();
    Double_t nHits      = mpdtrack-> GetNofHits();
    Double_t DCAX	= mpdtrack -> GetDCAX();
    Double_t DCAY	= mpdtrack -> GetDCAY();
    Double_t DCAZ	= mpdtrack -> GetDCAZ();
    float pt_reco	= mpdtrack -> GetPt();

//__specific variables for the kalmantrack_______________________________________________________________________________________________
   
    int     kfcharge = kftrack -> Charge();                    // for the charge
    double  p      = kftrack -> Momentum3().Mag() * kfcharge;  // full momentum and
    double  dedx   = kftrack -> GetDedx();                     // dE/dx information
    int     mcId   = kftrack -> GetTrackID();                  // MpdTpcKalmanTrack::GetTrackID() gives the ID of the corresponding Monte Carlo track
  
//__specific variables for the montecarlotrack___________________________________________________________________________________________
    
    MpdMCTrack* mctrack = (MpdMCTrack*) mMCTracks -> At(mcId);                       // Monte Carlo track is open for reading  
    int   pdg                 = mctrack -> GetPdgCode();                               // Track PDG code
    float rapidity_mc         = mctrack -> GetRapidity();                              // Particle rapidity (CAN BE WRONG!!!)
    float pt_mc               = mctrack -> GetPt();                                    // Particle transverse momentum
    float p_mc                = mctrack -> GetP();                                     // Particle full momentum
    float pz_mc               = mctrack -> GetPz();                                    // Particle momentum z-component

//______other variables for cuts_______________________________________________________________________	
        Double_t DPt = TMath::Abs( pt_reco - pt_mc ) / ( pt_mc );
        Double_t DCAG = TMath::Sqrt( pow(DCAX,2) + pow(DCAY,2) + pow(DCAZ,2));


//______Selection of the primary or secondary particles__________________________________________	
	//if(RECOmID == 1 && mctrack -> GetMotherId() != -1) continue; //primary tracks
	//if(RECOmID == 2 && mctrack -> GetMotherId() == -1) continue; //secondary tracks
	
//_____Bethe Bloch functions for each particle____________________________________________________ 

//Proton_sigma_________________________________________

  //double dedxp1=((-2.9144) / TMath::Power(p / TMath::Sqrt(p * p + 0.88), 0.82576))*((1.22346 - TMath::Power(p / TMath::Sqrt(p * p + 0.88), 0.82576))-(TMath::Log(2.0642 + TMath::Power(1.0 / (p / 0.9383), 2.20544))));//proton+2sigma
  double dedxp1=((-2.58072) / TMath::Power(p / TMath::Sqrt(p * p + 0.88), 0.782313))*((1.16569 - TMath::Power(p / TMath::Sqrt(p * p + 0.88), 0.782313))-(TMath::Log(1.97189 + TMath::Power(1.0 / (p / 0.9383), 2.20629))));//proton+1sigma

  //double dedxp2=((-3.11222) / TMath::Power(p / TMath::Sqrt(p * p + 0.88), 0.214233))*((1.7496 - TMath::Power(p / TMath::Sqrt(p * p + 0.88), 0.214233))-(TMath::Log(2.83206 + TMath::Power(1.0 / (p / 0.9383), 2.18993))));//proton-2sigma    
  double dedxp2=((-3.70863) / TMath::Power(p / TMath::Sqrt(p * p + 0.88), 0.434546))*((1.70197 - TMath::Power(p / TMath::Sqrt(p * p + 0.88), 0.434546))-(TMath::Log(2.64031 + TMath::Power(1.0 / (p / 0.9383), 2.0007))));//proton-1sigma    

//Kaon_sigma___________________________________________

  //double dedxk1=((-0.83892) / TMath::Power(p / TMath::Sqrt(p * p + 0.2437), 2.70241))*(((-0.307397) - TMath::Power(p / TMath::Sqrt(p * p + 0.2437), 2.70241))-(TMath::Log(0.514154 + TMath::Power(1.0 / (p / 0.4937), 0.27171))));//kaon+2sigma
  double dedxk1=((-0.817989) / TMath::Power(p / TMath::Sqrt(p * p + 0.2437), 2.51922))*(((-0.332111) - TMath::Power(p / TMath::Sqrt(p * p + 0.2437), 2.51922))-(TMath::Log(0.452453 + TMath::Power(1.0 / (p / 0.4937), 0.356776))));//kaon+1sigma
  
  //double dedxk2=((-0.824286) / TMath::Power(p / TMath::Sqrt(p * p + 0.2437), 0.323779))*(((-0.222323) - TMath::Power(p / TMath::Sqrt(p * p + 0.2437), 0.323779))-(TMath::Log(0.963711 + TMath::Power(1.0 / (p / 0.4937), 2.98519))));//kaon-2sigma
  double dedxk2=((-0.564779) / TMath::Power(p / TMath::Sqrt(p * p + 0.2437), 0.808589))*(((-1.11641) - TMath::Power(p / TMath::Sqrt(p * p + 0.2437), 0.808589))-(TMath::Log(0.799115 + TMath::Power(1.0 / (p / 0.4937), 3.07986))));//kaon-1sigma

//Pion_sigma__________________________________________

  //double dedxpi1=((-0.342358) / TMath::Power(p / TMath::Sqrt(p * p + 0.01949), 4.05206))*(((-2.42135) - TMath::Power(p / TMath::Sqrt(p * p + 0.01949), 4.05206))-(TMath::Log((-0.381729) + TMath::Power(1.0 / (p /0.1396), (-0.383815)))));//pion+2sigma
  double dedxpi1=((-0.327131) / TMath::Power(p / TMath::Sqrt(p * p + 0.01949), 3.94556))*(((-2.32774) - TMath::Power(p / TMath::Sqrt(p * p + 0.01949), 3.94556))-(TMath::Log((-0.340479) + TMath::Power(1.0 / (p /0.1396), (-0.371293)))));//pion+1sigma
  
  //double dedxpi2=((-0.198648) / TMath::Power(p / TMath::Sqrt(p * p + 0.01949), 3.19561))*(((-2.88049) - TMath::Power(p / TMath::Sqrt(p * p + 0.01949), 3.19561))-(TMath::Log(0.072076 + TMath::Power(1.0 / (p /0.1396), (-0.547553)))));//pion-2sigma
  double dedxpi2=((-0.182489) / TMath::Power(p / TMath::Sqrt(p * p + 0.01949), 3.10849))*(((-3.16187) - TMath::Power(p / TMath::Sqrt(p * p + 0.01949), 3.10849))-(TMath::Log(0.899613 + TMath::Power(1.0 / (p /0.1396), (-0.836771)))));//pion-1sigma


//_____Cuts in reconstructed tracks__________________________________________________________________________________________
    
       if (TMath::Abs(pt_reco) < limptreco) continue;
       if (nHits < limnHits ) continue;
       if (TMath::Abs(Eta) > EtaMax) continue;
       if (DCAG > limDCAG) continue;
//__________________________________________________________________________________________________________________________
       
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

	if(mctrack -> GetMotherId() != -1) continue;//for primary particles

//Limits in dE/dx of the proton______________________________________________________   
	if(dedx >= dedxp2 && dedx <= dedxp1 ){
          
          PtProreco1->Fill(pt_reco);	  
          if(pdg == 2212){ //we use the pdg code to clean up more signal and get the real particles
		h__dedxpcut->Fill(p,dedx);
                PtProreco2->Fill(pt_reco);         
           }
          
	}
//Limits in dE/dx of the kaon______________________________________________________
        if(dedx >= dedxk2 && dedx <= dedxk1 ){
         
          PtKaonreco1->Fill(pt_reco);
          if(pdg == 321){//we use the pdg code to clean up more signal and get the real particles
	
		h__dedxkpcut->Fill(p,dedx);
	        PtKaonreco2->Fill(pt_reco);
           }

        }
//Limits in dE/dx of the pion______________________________________________________
        if(dedx >= dedxpi2 && dedx <= dedxpi1 ){
 
          PtPionreco1->Fill(pt_reco);
          if(pdg == 211){//we use the pdg code to clean up more signal and get the real particles

		h__dedxpipcut->Fill(p,dedx); 
                PtPionreco2->Fill(pt_reco);
           }

        }
//______________________________________________________________
  }
 
 // pidQA->GetDedxQA(outpath);
}

/// Default destructor
EnerClass1::~EnerClass1(){}

void EnerClass1::Finish(){

}

