#include "Fixed_Analysis.h"

ClassImp(Fixed_Analysis);

Fixed_Analysis::Fixed_Analysis(const char *name, const char *outputName) : MpdAnalysisTask(name,outputName)
{
   fOutputList = nullptr;
}
//_____________________________________________________________________________________

void Fixed_Analysis::init_Hist1D(){

  char name1[15];

  char hist_name[300];
  char title_name[500];

  for (Int_t i = 0; i < 3; i++)
  {
    sprintf(name1,"");
    switch (i)
    {
      case 0:
        sprintf(name1,"Primary ");
        break;
      case 1:
        sprintf(name1,"Secondary ");
        break;
      case 2:
        sprintf(name1,"All ");
        break;
      default:
        sprintf(name1,"");
    }
    for (Int_t k = 0; k < NumHistPST1; k++)
    {
      sprintf(hist_name,"");
      sprintf(title_name,"");
      strcat(hist_name,name1);
      strcat(hist_name,HistPST1_titles[k]);
      strcat(title_name,hist_name);
      strcat(title_name,";");
      strcat(title_name,HistPST1_axis[k].X);
      HistPST1[i][k] = new TH1F(hist_name , title_name , HistPST1_bins[k].X , HistPST1_inter[k].Xlow , HistPST1_inter[k].Xup );
      fOutputList->Add(HistPST1[i][k]);   
    }
  }
  
  return;
}

void Fixed_Analysis::init_Hist2D(){
  char name1[15];
  char name2[15];

  char hist_name[300];
  char title_name[400];

  sprintf(hist_name,"");
  sprintf(title_name,"");
  sprintf(name1,"");
  sprintf(name2,"");

  for (Int_t j = 0; j < 3; j++)
  {
    sprintf(name2,"");
    switch (j)
    {
      case 0:
        sprintf(name2,"Primary ");
        break;
      case 1:
        sprintf(name2,"Secondary ");
        break;
      case 2:
        sprintf(name2,"All ");
        break;
      default:
        sprintf(name2,"");
    }
    for (Int_t i = 0; i < 2; i++)
    {
      sprintf(name1,"");
      if (i%2 == 0){
        sprintf(name1,"MC ");
      }else{
        sprintf(name1, "Reco ");
      }
      //Histograms 2D
      for (Int_t k = 0; k < NumHistPST2; k++)
      {
        sprintf(hist_name,"");
        sprintf(title_name,"");     

        strcat(hist_name,name1);
        strcat(hist_name,name2);
        strcat(hist_name,HistPST2_titles[k]);
        strcat(title_name,hist_name);
        strcat(title_name,";");
        strcat(title_name,HistPST2_axis[k].X);
        strcat(title_name,";");
        strcat(title_name,HistPST2_axis[k].Y);
        HistPriSecTot2[i][j][k] = new TH2F(hist_name , title_name , HistPST2_bins[k].X , HistPST2_inter[k].Xlow , HistPST2_inter[k].Xup , HistPST2_bins[k].Y , HistPST2_inter[k].Ylow , HistPST2_inter[k].Yup );
        fOutputList->Add(HistPriSecTot2[i][j][k]);
      } 
    }
      //Profiles 1D
    for (Int_t k = 0; k < NumProf1; k++)
    {
      sprintf(hist_name,"");
      sprintf(title_name,"");

      strcat(hist_name,name2);
      strcat(hist_name,Prof1_titles[k]);
      strcat(title_name,hist_name);
      strcat(title_name,";");
      strcat(title_name,Prof1_axis[k].X);
      strcat(title_name,";");
      strcat(title_name,Prof1_axis[k].Y);
      Prof1[j][k] = new TProfile(hist_name , title_name , Prof1_bins[k].X , Prof1_inter[k].Xlow , Prof1_inter[k].Xup , Prof1_inter[k].Ylow , Prof1_inter[k].Yup );
      fOutputList->Add(Prof1[j][k]);
    }
  }
  return;
}

void Fixed_Analysis::init_Profile2D(){
  char name2[15];

  char hist_name[300];
  char title_name[500];

  sprintf(hist_name,"");
  sprintf(title_name,"");
  sprintf(name2,"");
  for (Int_t i = 0; i < 3; i++)
  {
    sprintf(name2,"");
    switch (i)
    {
      case 0:
        sprintf(name2,"Primary ");
        break;
      case 1:
        sprintf(name2,"Secondary ");
        break;
      case 2:
        sprintf(name2,"All ");
        break;
      default:
        sprintf(name2,"");
    }
    for (Int_t k = 0; k < NumProf3; k++)
      {
        sprintf(hist_name,"");
        sprintf(title_name,"");     

        strcat(hist_name,name2);
        strcat(hist_name,Prof3D_titles[k]);
        strcat(title_name,hist_name);
        strcat(title_name,";");
        strcat(title_name,Prof3D_axis[k].X);
        strcat(title_name,";");
        strcat(title_name,Prof3D_axis[k].Y);
        strcat(title_name,";");
        strcat(title_name,Prof3D_axis[k].Z);
        // cout << title_name << endl;
        Prof3D[i][k] = new TProfile2D(hist_name , title_name , Prof3D_bins[k].X , Prof3D_inter[k].Xlow , Prof3D_inter[k].Xup , Prof3D_bins[k].Y , Prof3D_inter[k].Ylow , Prof3D_inter[k].Yup, Prof3D_inter[k].Zlow , Prof3D_inter[k].Zup );
        fOutputList->Add(Prof3D[i][k]);
      }
  }
  return;
}


//_____________________________________________________________________________________
void Fixed_Analysis::UserInit()
{
  
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  TH1::AddDirectory(kFALSE);

  hRefMult = new TH1F("hRefMultSTAR","hRefMultSTAR",300,0,300);
  hBvsRefMult = new TH2F("hBvsRefMult","hBvsRefMult",300,0,300,200,0.,20.);
  fOutputList->Add(hRefMult);
  fOutputList->Add(hBvsRefMult);

  //  Histograms

  Hist_MC_ZVertex = new TH1F("z-vertex_MC","MC z-Vertex; z-Vertex (cm); Entries",500,-200,200);
	Hist_R_ZVertex = new TH1F("z-vertex_REC","Reco z-Vertex; z-Vertex (cm); Entries",500,-200,200);
  fOutputList->Add(Hist_MC_ZVertex);
  fOutputList->Add(Hist_R_ZVertex);

	Hist_MC_NTrack_ZVertex = new TH2F("num_tracks_MC", " MC z-Vertex vs NTracks;z-Vertex (cm) ; NTraks",200,-200,200,500,0,300);
	Hist_R_NTrack_ZVertex = new TH2F("num_tracks_REC", " Reco z-Vertex vs NTracks;z-Vertex (cm) ; NTraks",200,-200,200,500,0,300);
  fOutputList->Add(Hist_MC_NTrack_ZVertex);
  fOutputList->Add(Hist_R_NTrack_ZVertex);

	Hist_Chi = new TH1F("chi","#frac{#chi^{2}}{NDF}; #frac{#chi^{2}}{NDF}; Entries",300,0,12);
  fOutputList->Add(Hist_Chi);

	Hist_MC_Chi_zVertex = new TH2F("chi_zVertex_MC", "  ZvertexMC vs #frac{#chi^{2}}{NDF};z-Vertex (cm) ; #frac{#chi^{2}}{NDF}",100,-100,100,100,0,12);
	Hist_R_Chi_zVertex = new TH2F("chi_zVertex_R", " ZvertexReco vs #frac{#chi^{2}}{NDF};z-Vertex (cm) ; #frac{#chi^{2}}{NDF}",500,-150,150,100,0,12);
  fOutputList->Add(Hist_MC_Chi_zVertex);
  fOutputList->Add(Hist_R_Chi_zVertex);

  Hist_Multipl_MC = new TH1F("hMultiplicity_MC","MC hMultiplicity",1000,0,5000);
  Hist_Multipl_R = new TH1F("hMultiplicity_Reco","Reco hMultiplicity",1000,0,1000);
  fOutputList->Add(Hist_Multipl_MC);
  fOutputList->Add(Hist_Multipl_R);
  
	Hist_b_Vertex = new TH2F("b_Vertex", "Reco z-Vertex vs Impact Parameter;z-Vertex (cm) ; b (fm)",500,-150,150,100,0,16);
  Hist_b = new TH1F("Impact_b","Impact Parameter; b (fm); Entries",100,0,16);
  fOutputList->Add(Hist_b_Vertex);
  fOutputList->Add(Hist_b);

  Hist_b_Multiplicity_R = new TH2F("b_hMultiplicity_Reco", "Reco hMultiplicity vs Impact Parameter;hMultiplicity ; b (fm)",1000,0,1000,100,0,16);
  fOutputList->Add(Hist_b_Multiplicity_R);

  Hist_Diff_ZVertex_b = new TH2F("b_DiffVertex", "#Delta z-Vertex vs Impact Parameter;#Delta z-Vertex (cm) ; b (fm)",500,-200,200,100,0,16);
  Hist_Diff_ZVertex_NTracks = new TH2F("Diff_Vertex_num_tracks", " #Delta z-Vertex vs NTracks;#Delta z-Vertex; NTraks",100,-10,10,500,0,300);

  fOutputList->Add(Hist_Diff_ZVertex_b);
  fOutputList->Add(Hist_Diff_ZVertex_NTracks);


  ZVertex_NTracksMC = new TH2F("Diff_Vertex_num_tracksMC", " #Delta z-Vertex vs NTracksMC;#Delta z-Vertex; NTraksMC",100,-4,4,500,200,700);
  fOutputList->Add(ZVertex_NTracksMC);


  Prof_2D_eta_pt_DPT_Reco = new TProfile2D("PtEtaDpTReco", "Reco #eta vs pT ; #eta ; pT ; #Delta p_{T}",100,-4,4,100,0,5,0,4);
  fOutputList->Add(Prof_2D_eta_pt_DPT_Reco);
  Prof_2D_eta_pt_DPT_MC = new TProfile2D("PtEtaDpTMC", "Reco #eta vs pT ; #eta ; pT ; #Delta p_{T}",100,-4,4,100,0,5,0,4);
  fOutputList->Add(Prof_2D_eta_pt_DPT_MC);

  init_Profile2D();
  init_Hist1D();
  init_Hist2D();
  

  MC_pT = new TH1F("MC_PT_cuts","MC p_T cuts; p_T (Gev/c) ; Entries",100,0,2.5);
	Reco_pT = new TH1F("Reco_PT_cuts","Reco p_T cuts; p_T (Gev/c); Entries",100,0,2.5);
  fOutputList->Add(MC_pT);
  fOutputList->Add(Reco_pT);

}
//____________________________________________________________________________________
void Fixed_Analysis::ProcessEvent(MpdAnalysisEvent &event)
{

  // if(!selectEvent(event)){
  // return;
  // }

  //  Branches
  fVertex = event.fVertex;
  fMCEventHeader = event.fMCEventHeader;

      // First Vertex
  Vertex = (MpdVertex*) fVertex->First();

  fMCTracks = event.fMCTrack;
  
  // Variables
  Int_t NumTracksMC;
  Int_t NumTracksVertex;
  Int_t NumTracksReco;
  Double_t Chi2;
  Double_t NDF;
  Double_t Z;
  Double_t Impactb;

    //Number Tracks for MC, Reco, Vertex
  NumTracksMC = fMCTracks->GetEntriesFast();
  NumTracksVertex = Vertex->GetNTracks();
  NumTracksReco = event.fMPDEvent->GetGlobalTracks()->GetEntriesFast();

  //Vertex Variables
  Chi2 = Vertex->GetChi2();
  NDF = Vertex->GetNDF();
  Z = Vertex->GetZ();

    //Event Header variables
  Impactb = fMCEventHeader->GetB();


  Double_t ZMC = fMCEventHeader->GetZ();
  Double_t ZReco = Vertex->GetZ();
  Double_t DzVertex = TMath::Abs((ZReco - ZMC)/ZMC);

  
  Int_t CounterMC = 0;

    //MonteCarlo Loop
  for (Int_t i = 0; i < NumTracksMC; i++){
    MpdMCTrack *MCTrack = (MpdMCTrack*)fMCTracks->UncheckedAt(i);

    Int_t abspdgcode;
    Double_t pTMC;
    Double_t P;
    Double_t Pz;
    Double_t Eta;

    abspdgcode = TMath::Abs(MCTrack->GetPdgCode());

    pTMC = MCTrack->GetPt(); 
    P = MCTrack->GetP();
    Pz = MCTrack->GetPz();
    Eta = 0.5*TMath::Log((P + Pz)/(P - Pz));

    if (abspdgcode == 321 || abspdgcode == 2212 || abspdgcode == 211 || abspdgcode == 11)// kaons, protons, pions, electrons
    {
      if(pTMC <= 2.5 && MCTrack->GetMotherId() == -1 && Impactb < 13 && Eta > -1 && Eta < 2){//Matching cuts of RecoTracks
        MC_pT->Fill(pTMC);
      }
    }

    if(abspdgcode == 321 || abspdgcode == 2212 || abspdgcode == 211){// kaons, protons, pions
		  CounterMC++; //Counting the number of kaones, protons and pions
		}
  }

  //Histograms with no cuts
  Hist_Diff_ZVertex_b->Fill(DzVertex,Impactb);
  Hist_Diff_ZVertex_NTracks->Fill(DzVertex,NumTracksVertex);
  ZVertex_NTracksMC->Fill(DzVertex,CounterMC);

  if(CounterMC > 308){//Cut on number of particles > Xe(124) + W(184). Only particles with more than the initial ones

    if(NumTracksVertex > 3) { //Cut on Number of Tracks Vertex

      //Histograms for Vertex Analysis
        //Reco
      Hist_R_NTrack_ZVertex->Fill(ZReco,NumTracksVertex);
      Hist_R_ZVertex->Fill(ZReco);
      Hist_R_Chi_zVertex->Fill(ZReco,(Chi2/NDF));
      Hist_Chi->Fill((Chi2/NDF));
      Hist_b->Fill(Impactb);
      Hist_b_Vertex->Fill(ZReco,Impactb);
        //MC
      Hist_MC_NTrack_ZVertex->Fill(ZMC,NumTracksVertex);
      Hist_MC_ZVertex->Fill(ZMC);
      Hist_MC_Chi_zVertex->Fill(ZMC,(Chi2/NDF)); 
      //Multiplicity - Cuts on MC and Reco Tracks
      Hist_Multipl_MC->Fill(CounterMC);
      Hist_Multipl_R->Fill(NumTracksReco);
      Hist_b_Multiplicity_R->Fill(NumTracksReco,Impactb);

      //Reconstructed Loop
      Int_t Multiplicity = 0;
      for (Int_t j = 0; j < NumTracksReco; j++){
        MpdTrack *RecoTrack = (MpdTrack*) event.fMPDEvent->GetGlobalTracks()->UncheckedAt(j);
        // ------ Montecarlo Identification ----

        Int_t ID = RecoTrack->GetID();
        MpdMCTrack *MCTrack = (MpdMCTrack*)fMCTracks->UncheckedAt(ID);

        // -------------------------------------

        Double_t Pt_MC;
        Double_t P_MC;
        Double_t Pz_MC;
        Double_t Eta_MC;
        Double_t Pt_Reco;
        Double_t Eta_Reco;
        Double_t Diff_Pt;
        Double_t DCAX;
        Double_t DCAY;
        Double_t DCAZ;
        Double_t DCA;
        Double_t NumHits;

        Pt_MC = MCTrack->GetPt();
        P_MC = MCTrack->GetP();
        Pz_MC = MCTrack->GetPz();

        Pt_Reco = RecoTrack->GetPt();

        //Resolution of pT - Definition
        Diff_Pt = TMath::Abs(Pt_Reco - Pt_MC)/Pt_MC;

        Eta_Reco = RecoTrack->GetEta();
        Eta_MC = 0.5*TMath::Log((P_MC + Pz_MC)/(P_MC - Pz_MC + 1e-16));// 1e-16 to avoid zero


        DCAX = RecoTrack->GetDCAX();
        DCAY = RecoTrack->GetDCAY();
        DCAZ = RecoTrack->GetDCAZ();

        //DCA - Definition
        DCA = TMath::Sqrt(TMath::Power(DCAX,2) + TMath::Power(DCAY,2) + TMath::Power(DCAZ,2) );

        NumHits = RecoTrack->GetNofHits();

        //No Cuts on eta, DCA, NHits
        Prof_2D_eta_pt_DPT_MC->Fill(Eta_MC,Pt_MC,Diff_Pt);
        Prof_2D_eta_pt_DPT_Reco->Fill(Eta_Reco,Pt_Reco,Diff_Pt);

        //Primary p = 0 ; Secondary p = 1 
        Int_t p = (MCTrack->GetMotherId() == -1) ? 0 : 1; 

        for (Int_t k = 0; k < 2; k++)//MC k = 0 ; Reco k = 1 
        {
          Double_t eta = (k%2 == 0) ? Eta_MC : Eta_Reco; 
          Double_t pt = (k%2 == 0) ? Pt_MC : Pt_Reco;
          //Fill MC or Reco [k]; Primary, Secondary [p] and All [2]
            //No Cuts on eta, DCA, NHits
          HistPriSecTot2[k][p][0]->Fill(eta,pt);
          HistPriSecTot2[k][2][0]->Fill(eta,pt);
        }

        //No Cuts on eta, DCA, NHits
        // Fill Eta vs D pT
        Prof1[p][0]->Fill(Eta_Reco,Diff_Pt);
        Prof1[2][0]->Fill(Eta_Reco,Diff_Pt);
        //Fill DCA vs D pT
        Prof1[p][1]->Fill(DCA,Diff_Pt);
        Prof1[2][1]->Fill(DCA,Diff_Pt);
        //Fill NHits vs D pT
        Prof1[p][2]->Fill(NumHits,Diff_Pt);
        Prof1[2][2]->Fill(NumHits,Diff_Pt);

        if (ZReco > -100 && ZReco < -70)//Cut on Primary Vertex
        {
          
          // Fill Eta vs D pT
          Prof1[p][3]->Fill(Eta_Reco,Diff_Pt);
          Prof1[2][3]->Fill(Eta_Reco,Diff_Pt);
          //Fill DCA vs D pT
          Prof1[p][4]->Fill(DCA,Diff_Pt);
          Prof1[2][4]->Fill(DCA,Diff_Pt);
          //Fill NHits vs D pT
          Prof1[p][5]->Fill(NumHits,Diff_Pt);
          Prof1[2][5]->Fill(NumHits,Diff_Pt);

          //Fill Eta vs pT weight DpT
          Prof3D[p][0]->Fill(Eta_Reco,Pt_Reco,Diff_Pt);
          Prof3D[2][0]->Fill(Eta_Reco,Pt_Reco,Diff_Pt);

          //Fill DCA
          HistPST1[p][0]->Fill(DCA);
          HistPST1[2][0]->Fill(DCA);

          //First cut NHits then Eta
          if(NumHits >= 20){ //Cut on Number of Hits
            //Fill Eta vs pT weight DpT
            Prof3D[p][1]->Fill(Eta_Reco,Pt_Reco,Diff_Pt);
            Prof3D[2][1]->Fill(Eta_Reco,Pt_Reco,Diff_Pt);
            //Fill DCA
            HistPST1[p][1]->Fill(DCA);
            HistPST1[2][1]->Fill(DCA);
            
            if (Eta_Reco > -1 && Eta_Reco < 2) // Cut on Eta
            {
              //Fill Eta vs pT weight DpT
              Prof3D[p][2]->Fill(Eta_Reco,Pt_Reco,Diff_Pt);
              Prof3D[2][2]->Fill(Eta_Reco,Pt_Reco,Diff_Pt);
              //Fill DCA
              HistPST1[p][2]->Fill(DCA);
              HistPST1[2][2]->Fill(DCA);
            }
            
          }

          //First cut Eta then NHits
          if (Eta_Reco > -1 && Eta_Reco < 2)// Cut on Eta
          {
            //Fill DCA vs D pT
            Prof1[p][6]->Fill(DCA,Diff_Pt);
            Prof1[2][6]->Fill(DCA,Diff_Pt);
            //Fill NHits vs D pT
            Prof1[p][7]->Fill(NumHits,Diff_Pt);
            Prof1[2][7]->Fill(NumHits,Diff_Pt);

            
            if(NumHits >= 20){//Cut on Number of Hits
              //Fill DCA vs D pT
              Prof1[p][8]->Fill(DCA,Diff_Pt);
              Prof1[2][8]->Fill(DCA,Diff_Pt);
            }

          }
          
        }
        

        if (Impactb < 13 && Eta_Reco > -1 && Eta_Reco < 2 && DCA <= 2 && NumHits >= 20 ){ //Impact Parameter Cut
          Prof3D[p][3]->Fill(Eta_Reco,Pt_Reco,Diff_Pt);
          Prof3D[2][3]->Fill(Eta_Reco,Pt_Reco,Diff_Pt);
          Multiplicity++;
          if (Pt_Reco <= 2.5)
          {
            Reco_pT->Fill(Pt_Reco);
          }
        }
      }
      //Fill Multiplicity cuts on b, eta, DCA, NHits
      hRefMult->Fill(Multiplicity);
      hBvsRefMult->Fill(Multiplicity,Impactb);
    }
  }

}
//___________________________________________________________________________________
void Fixed_Analysis::Finish()
{
}



