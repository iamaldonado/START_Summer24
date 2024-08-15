#include "Fixed_Track.h"

ClassImp(Fixed_Track);

Fixed_Track::Fixed_Track(const char *name, const char *outputName) : MpdAnalysisTask(name,outputName)
{
   fOutputList = nullptr;
}
//_____________________________________________________________________________________
void Fixed_Track::init_HistPart(){

  char name1[20];
  char name2[20];
  char name3[20];
  char hist_name[300];
  char title_name[500];

  for(Int_t i = 0 ; i < 2 ; i++)
  {
    sprintf(name1,"");
    Int_t u = (i%2 == 0) ? sprintf(name1,"MC_") : sprintf(name1, "Reco_");
    for (Int_t j = 0; j < 2; j++)
    {
      sprintf(name2,"");
      Int_t y = (j%2 == 0) ? sprintf(name2,"Primary_") : sprintf(name2,"Secondary_");
      for (Int_t l = 0; l < 4; l++)
      {
        sprintf(name3,"");
        switch (l)
        {
          case 0:
            sprintf(name3,"Kaons_");
            break;
          case 1:
            sprintf(name3,"Protons_");
            break;
          case 2:
            sprintf(name3,"Pions_");
            break;
          case 3:
            sprintf(name3,"All_");
            break;
          default:
            sprintf(name3,"Error");
        }
        for (Int_t k = 0; k < NumHistMCPSKpPi1; k++)
        {
          sprintf(hist_name,"");
          sprintf(title_name,"");

          strcat(hist_name,name1);
          strcat(hist_name,name2);
          strcat(hist_name,name3);

          strcat(hist_name,HistMCPSKpPi1_titles[k]);
          strcat(title_name,hist_name);
          strcat(title_name,";");
          strcat(title_name,HistMCPSKpPi1_axis[k].X);
          HistMCPSKpPi1[i][j][l][k] = new TH1F(hist_name , title_name , HistMCPSKpPi1_bins[k].X , HistMCPSKpPi1_inter[k].Xlow , HistMCPSKpPi1_inter[k].Xup );
          fOutputList->Add(HistMCPSKpPi1[i][j][l][k]);   
        }
      } 
    } 
  }
  return;
}

void Fixed_Track::init_Hist1D(){

  char name1[15];
  char name2[15];
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

void Fixed_Track::init_Hist2D(){
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

void Fixed_Track::init_Profile2D(){
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
      Prof3D[i][k] = new TProfile2D(hist_name , title_name , Prof3D_bins[k].X , Prof3D_inter[k].Xlow , Prof3D_inter[k].Xup , Prof3D_bins[k].Y , Prof3D_inter[k].Ylow , Prof3D_inter[k].Yup, Prof3D_inter[k].Zlow , Prof3D_inter[k].Zup );
      fOutputList->Add(Prof3D[i][k]);
    }
  }
  return;
}


//_____________________________________________________________________________________
void Fixed_Track::UserInit()
{
  
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  TH1::AddDirectory(kFALSE);

  init_HistPart();

  DZVertex_Multi = new TH2F("DzVertex_Multi", " z-Vertex vs Tracks Reco ; Tracks Reco; #Delta z-Vertex",200,0,200,300,0,3);
  fOutputList->Add(DZVertex_Multi);

  DZVertex_Multi_Profile = new TProfile("DzVertex_Multi", " z-Vertex vs Tracks Reco ; Tracks Reco; #Delta z-Vertex",200,0,200,0,3);
  fOutputList->Add(DZVertex_Multi_Profile);

  Dpt_pt_Multi = new TH2F("Dpt_pt_Multi", "pT vs #Delta pT / pT ;pT ; #Delta pT / pT",200,0,5,300,0,3);
  fOutputList->Add(Dpt_pt_Multi);

  Dpt_pt_Multi_Profile = new TProfile("Dpt_pt_Multi", "pT vs #Delta pT ;pT ; #Delta pT",200,0,5,0,3);
  fOutputList->Add(Dpt_pt_Multi_Profile);

  Dpt_pt_Multi_Profile_cuts = new TProfile("Dpt_pt_Multi", "pT vs #Delta pT ;pT ; #Delta pT",200,0,5,0,3);
  fOutputList->Add(Dpt_pt_Multi_Profile_cuts);

  DZVertex_zVertex = new TH2F("DzVertex_zVertex", " z-Vertex vs #Delta z-Vertex ; z-Vertex (cm); #Delta z-Vertex",500,-200,200,200,0,4);
  fOutputList->Add(DZVertex_zVertex);

  zVertex_b = new TH2F("zVertex_b", " z-Vertex vs b ; z-Vertex (cm); b (fm)",500,-200,200,100,0,16);
  fOutputList->Add(zVertex_b);
}
//____________________________________________________________________________________
void Fixed_Track::ProcessEvent(MpdAnalysisEvent &event)
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
  Double_t DzVertex;

  DzVertex = TMath::Abs((ZReco - ZMC)/ZMC);
 
  
  
  Int_t CounterMC = 0;

    //MonteCarlo Loop
  for (Int_t i = 0; i < NumTracksMC; i++){
    MpdMCTrack *MCTrack = (MpdMCTrack*)fMCTracks->UncheckedAt(i);

    Int_t abspdgcode;
    Double_t pTMC;
    Double_t P;
    Double_t Pz;
    Double_t Eta;

    
    //Primary p = 0 ; Secondary p = 1 
    Int_t p = (MCTrack->GetMotherId() == -1) ? 0 : 1;

    //Kaons, Protons, Pions
    abspdgcode = TMath::Abs(MCTrack->GetPdgCode());
    Int_t l = 3;

    switch (abspdgcode)
    {
      case 321://Kaons
        l = 0;
        break;
      case 2212://Protons
        l = 1;
        break;
      case 211://Pions
        l = 2;
        break;    
    }

    pTMC = MCTrack->GetPt(); 
    P = MCTrack->GetP();
    Pz = MCTrack->GetPz();
    Eta = 0.5*TMath::Log((P + Pz)/(P - Pz));

    if (abspdgcode == 321 || abspdgcode == 2212 || abspdgcode == 211 || abspdgcode == 11)// kaons, protons, pions, electrons
    {
      if(Impactb < 13 && Eta > -1 && Eta < 2 && ZReco > -100 && ZReco < -70)
      {//Matching cuts of RecoTracks
        if(pTMC == 0) continue;
        HistMCPSKpPi1[0][p][l][0]->Fill(pTMC);//p primary or secondary, l - kaons, protons, pions
        HistMCPSKpPi1[0][p][3][0]->Fill(pTMC);//3-all particles
      }

    }

    if(abspdgcode == 321 || abspdgcode == 2212 || abspdgcode == 211){// kaons, protons, pions
		  CounterMC++; //Counting the number of kaones, protons and pions
		}
  }

  //Histograms with no cuts
  

  if(CounterMC > 308){//Cut on number of particles > Xe(124) + W(184). Only particles with more than the initial ones

    if(NumTracksVertex > 3) { //Cut on Number of Tracks Vertex
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

        //Primary p = 0 ; Secondary p = 1 
        Int_t p = (MCTrack->GetMotherId() == -1) ? 0 : 1;

        //Kaons, Protons, Pions
        Int_t abspgdcode = TMath::Abs(MCTrack->GetPdgCode());
        Int_t l = 3;

        switch (abspgdcode)
        {
          case 321://Kaons
            l = 0;
            break;
          case 2212://Protons
            l = 1;
            break;
          case 211://Pions
            l = 2;
            break;    
        }
        
        //No Cuts on eta, DCA, NHits
        Dpt_pt_Multi->Fill(Pt_Reco,Diff_Pt);
        Dpt_pt_Multi_Profile->Fill(Pt_Reco,Diff_Pt);
        

        if (ZReco > -100 && ZReco < -70)//Cut on Primary Vertex
        {
          if (Impactb < 13 && Eta_Reco > -1 && Eta_Reco < 2 && DCA <= 2 && NumHits >= 20 )
          { //Impact Parameter Cut
            Multiplicity++;

            if(Pt_MC == 0) continue;
            HistMCPSKpPi1[1][p][l][0]->Fill(Pt_Reco);
            HistMCPSKpPi1[1][p][3][0]->Fill(Pt_Reco);

            Dpt_pt_Multi_Profile_cuts->Fill(Pt_Reco,Diff_Pt);
          }
        } 
      }
      DZVertex_zVertex->Fill(ZReco,DzVertex);
      DZVertex_Multi->Fill(Multiplicity,DzVertex);
      DZVertex_Multi_Profile->Fill(Multiplicity,DzVertex);
      zVertex_b->Fill(ZReco,Impactb);
    }
  }
}
//___________________________________________________________________________________
void Fixed_Track::Finish()
{
}



