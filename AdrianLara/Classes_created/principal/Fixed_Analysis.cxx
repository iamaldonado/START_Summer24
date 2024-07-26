#include "Fixed_Analysis.h"

ClassImp(Fixed_Analysis);

Fixed_Analysis::Fixed_Analysis(const char *name, const char *outputName) : MpdAnalysisTask(name,outputName)
{
   fOutputList = nullptr;
}
//_____________________________________________________________________________________

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
  Hist_Diff_ZVertex_NTracks = new TH2F("Diff_Vertex_num_tracks", " #Delta z-Vertex vs NTracks;#Delta z-Vertex (cm) ; NTraks",200,-200,200,500,0,300);
  fOutputList->Add(Hist_Diff_ZVertex_b);
  fOutputList->Add(Hist_Diff_ZVertex_NTracks);

  // Profile_Diff_Vertex_b = new TProfile("DiffVertex_b_Profile","Impact Parameter vs #Delta z-Vertex",100,0,16,-200,200);
  // Profile_Diff_Vertex_NTracks = new TProfile("DiffVertex_Ntracks_Profile","NTracks vs #Deltaz-Vertex",500,0,200,-200,200);
  // fOutputList->Add(Profile_Diff_Vertex_b);
  // fOutputList->Add(Profile_Diff_Vertex_NTracks);



  Hist_Pt_Eta_Reco = new TH2F("EtavsPtReco", " Reco #eta vs pT ;#eta ; pT",100,-3,3,100,0,3);
  Hist_Pt_Eta_MC = new TH2F("EtavsPtMC", "MC #eta vs pT ; #eta ; pT ; #Delta p_{T}",100,-3,3,100,0,3);
  fOutputList->Add(Hist_Pt_Eta_Reco);
  fOutputList->Add(Hist_Pt_Eta_MC);

  Prof_2D_eta_pt_DPT_Reco = new TProfile2D("PtEtaDpTReco", "Reco #eta vs pT ; #eta ; pT ; #Delta p_{T}",100,-4,4,100,0,5,0,4);
  fOutputList->Add(Prof_2D_eta_pt_DPT_Reco);
  Prof_2D_eta_pt_DPT_MC = new TProfile2D("PtEtaDpTMC", "Reco #eta vs pT ; #eta ; pT ; #Delta p_{T}",100,-4,4,100,0,5,0,4);
  fOutputList->Add(Prof_2D_eta_pt_DPT_MC);

  // char type[30];

  // for (Int_t i = 0; i < 2; i++)
  // {
  //   if (i%2 == 0){ 
  //     sprintf(type,"Primary");
  //   }else{
  //     sprintf(type,"Secondary");
  //   }
    
  //   Profile_Diff_Pt_Eta[i] = new TProfile(Form("%s DiffPt_eta_Profile",&type),Form("%s #eta vs #Delta p_{T} ;#eta ; #frac{|pT_{Reco} - pT_{MC}|}{pT_{MC}}",&type),100,-3,3,0,3);
  //   fOutputList->Add(Profile_Diff_Pt_Eta[i]);
  // }
  
  
  // Profile_Diff_Pt_DCA = new TProfile("DiffPt_DCA_Profile","DCAZ vs #Delta p_{T} ; DCAZ ; #frac{|pT_{Reco} - pT_{MC}|}{pT_{MC}}",100,0,5,0,3);
  // Profile_Diff_Pt_NumHits = new TProfile("DiffPt_NHits_Profile","NHits vs #Delta p_{T}; NHits ; #frac{|pT_{Reco} - pT_{MC}|}{pT_{MC}}",100,0,50,0,3);
  
  // fOutputList->Add(Profile_Diff_Pt_DCA);
  // fOutputList->Add(Profile_Diff_Pt_NumHits);
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
      for (Int_t k = 0; k < NumProf3; k++)
      {
        sprintf(hist_name,"");
        sprintf(title_name,"");     

        strcat(hist_name,name1);
        strcat(hist_name,name2);
        strcat(hist_name,Prof3D_titles[k]);
        strcat(title_name,hist_name);
        strcat(title_name,";");
        strcat(title_name,Prof3D_axis[k].X);
        strcat(title_name,";");
        strcat(title_name,Prof3D_axis[k].Y);
        strcat(title_name,";");
        strcat(title_name,Prof3D_axis[k].Z);
        Prof3D[i][j][k] = new TProfile2D(hist_name , title_name , Prof3D_bins[k].X , Prof3D_inter[k].Xlow , Prof3D_inter[k].Xup , Prof3D_bins[k].Y , Prof3D_inter[k].Ylow , Prof3D_inter[k].Yup, Prof3D_inter[k].Zlow , Prof3D_inter[k].Zup );
        fOutputList->Add(Prof3D[i][j][k]);
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


  // Profile_Diff_Pt_Eta_Cuts = new TProfile("DiffPt_eta_Profile_Cuts","Cut #eta vs #Delta p_{T} ;#eta ; #frac{|pT_{Reco} - pT_{MC}|}{pT_{MC}} ",100,-3,3,0,3);
  // Profile_Diff_Pt_DCA_Cuts = new TProfile("DiffPt_DCA_Profile_Cuts","Cut DCA vs #Delta p_{T} ; NHits ; #frac{|pT_{Reco} - pT_{MC}|}{pT_{MC}}",100,0,5,0,3);
  // Profile_Diff_Pt_NumHits_Cuts = new TProfile("DiffPt_NHits_Profile_Cuts","Cut NHits vs #Delta p_{T} ; NHits ; #frac{|pT_{Reco} - pT_{MC}|}{pT_{MC}}",100,0,50,0,3);
  // fOutputList->Add(Profile_Diff_Pt_Eta_Cuts);
  // fOutputList->Add(Profile_Diff_Pt_DCA_Cuts);
  // fOutputList->Add(Profile_Diff_Pt_NumHits_Cuts);

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

    //Number Tracks for MC, Reco, Vertex
  Int_t NumTracksMC = fMCTracks->GetEntriesFast();
  Int_t NumTracks_Vertex = Vertex->GetNTracks();
  Int_t NumTracks_Reco = event.fMPDEvent->GetGlobalTracks()->GetEntriesFast();

    //Vertex Variables
  Double_t Chi2 = Vertex->GetChi2();
  Double_t NDF = Vertex->GetNDF();
  Double_t Z = Vertex->GetZ();

    //Event Header variables
  Double_t Impactb = fMCEventHeader->GetB();

  
  Int_t CounterMC = 0;

    //MonteCarlo Loop
  for (Int_t i = 0; i < NumTracksMC; i++){
    MpdMCTrack *MCTrack = (MpdMCTrack*)fMCTracks->UncheckedAt(i);

    Int_t abspdgcode = TMath::Abs(MCTrack->GetPdgCode());


    if(abspdgcode == 321 || abspdgcode == 2212 || abspdgcode == 211){
		  CounterMC++;
		}
  }

  if(CounterMC > 308){

    Double_t ZMC = fMCEventHeader->GetZ();
    Double_t ZReco = Vertex->GetZ();
    Double_t Diff = ZReco - ZMC;

    Hist_Diff_ZVertex_b->Fill(Diff,Impactb);
    Hist_Diff_ZVertex_NTracks->Fill(Diff,NumTracks_Vertex);

    if(NumTracks_Vertex > 3) { //Cut on Number of Tracks

      Hist_R_NTrack_ZVertex->Fill(Z,NumTracks_Vertex);
      Hist_R_ZVertex->Fill(Z);
      Hist_R_Chi_zVertex->Fill(Z,(Chi2/NDF));
      Hist_Chi->Fill((Chi2/NDF));
      Hist_b->Fill(Impactb);
      Hist_b_Vertex->Fill(Z,Impactb);
      //MC
      Z = fMCEventHeader->GetZ();
      Hist_MC_NTrack_ZVertex->Fill(Z,NumTracks_Vertex);
      Hist_MC_ZVertex->Fill(Z);
      Hist_MC_Chi_zVertex->Fill(Z,(Chi2/NDF)); 
      //  Multiplicity
      Hist_Multipl_MC->Fill(CounterMC);
      Hist_Multipl_R->Fill(NumTracks_Reco);
      Hist_b_Multiplicity_R->Fill(NumTracks_Reco,Impactb);




      //Reconstructed Loop
      Int_t Multiplicity = 0;
      for (Int_t j = 0; j < NumTracks_Reco; j++){
        MpdTrack *RecoTrack = (MpdTrack*) event.fMPDEvent->GetGlobalTracks()->UncheckedAt(j);
        Int_t ID = RecoTrack->GetID();
        MpdMCTrack *MCTrack = (MpdMCTrack*)fMCTracks->UncheckedAt(ID);

        

        Double_t Pt_MC = MCTrack->GetPt();
        Double_t Pt_Reco = RecoTrack->GetPt();
        Double_t Diff_Pt = TMath::Abs(Pt_Reco - Pt_MC)/Pt_MC;

        Double_t P_MC = MCTrack->GetP();
        Double_t Pz_MC = MCTrack->GetPz();

        Double_t Eta_Reco = RecoTrack->GetEta();
        Double_t Eta_MC = 0.5*TMath::Log((P_MC + Pz_MC)/(P_MC - Pz_MC + 1e-16));


        Double_t DCAX = RecoTrack->GetDCAX();
        Double_t DCAY = RecoTrack->GetDCAY();
        Double_t DCAZ = RecoTrack->GetDCAZ();
        Double_t DCA = TMath::Sqrt(TMath::Power(DCAX,2) + TMath::Power(DCAY,2) + TMath::Power(DCAZ,2) );

        Double_t NumHits = RecoTrack->GetNofHits();

          // No Cuts

        Hist_Pt_Eta_MC->Fill(Eta_MC,Pt_MC,(Pt_Reco/Pt_MC));
        Prof_2D_eta_pt_DPT_MC->Fill(Eta_MC,Pt_MC,Diff_Pt);
        Prof_2D_eta_pt_DPT_Reco->Fill(Eta_Reco,Pt_Reco,Diff_Pt);

        Hist_Pt_Eta_Reco->Draw("COLZ");
        Hist_Pt_Eta_MC->Draw("COLZ");
        
        Int_t p = (MCTrack->GetMotherId() == -1) ? 0 : 1; 

        for (Int_t k = 0; k < 2; k++)
        {
          Double_t eta = (k%2 == 0) ? Eta_MC : Eta_Reco;
          Double_t pt = (k%2 == 0) ? Pt_MC : Pt_Reco;
          HistPriSecTot2[k][p][0]->Fill(eta,pt);
          HistPriSecTot2[k][2][0]->Fill(eta,pt);
        }

        // Fill Eta vs D pT
        Prof1[p][0]->Fill(Eta_Reco,Diff_Pt);
        Prof1[2][0]->Fill(Eta_Reco,Diff_Pt);
        //Fill DCA vs D pT
        Prof1[p][1]->Fill(DCA,Diff_Pt);
        Prof1[2][1]->Fill(DCA,Diff_Pt);
        //Fill NHits vs D pT
        Prof1[p][2]->Fill(NumHits,Diff_Pt);
        Prof1[2][2]->Fill(NumHits,Diff_Pt);

        
        
        if (ZReco > -100 && ZReco < -70)
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

          if (Eta_Reco > -1 && Eta_Reco < 2)
          {
            //Fill DCA vs D pT
            Prof1[p][6]->Fill(DCA,Diff_Pt);
            Prof1[2][6]->Fill(DCA,Diff_Pt);
            //Fill NHits vs D pT
            Prof1[p][7]->Fill(NumHits,Diff_Pt);
            Prof1[2][7]->Fill(NumHits,Diff_Pt);
            
            if(NumHits >= 20){
              //Fill DCA vs D pT
              Prof1[p][8]->Fill(DCA,Diff_Pt);
              Prof1[2][8]->Fill(DCA,Diff_Pt);
            }

          }
          
        }
        

        if (Impactb < 13 && Eta_Reco > -1 && Eta_Reco < 2 && DCA <= 2 && NumHits >= 20 ){ //Impact Parameter Cut
          Multiplicity++;
        }
      }

      hRefMult->Fill(Multiplicity);
      hBvsRefMult->Fill(Multiplicity,Impactb);
    }
  }

}
//___________________________________________________________________________________
void Fixed_Analysis::Finish()
{
}



