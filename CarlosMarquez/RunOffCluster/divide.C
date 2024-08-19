void divide(){
	
	TFile *fA = TFile::Open("tasklowMgF.root");

	TTree *t2 = (TTree*)fA -> Get("demo/AnaTree");
	
	  // Histrograms
	// pT vs Eta 
   TH2F *h1	 =  (TH2F*)fA->Get("PtRecoPionP");
   TH2F *h2	 =  (TH2F*)fA->Get("PtRecoProtonP");
   TH2F *h3	 =  (TH2F*)fA->Get("PtRecoKaonP");
   TH2F *h1d	 =  (TH2F*)fA->Get("PtMCPionP");
   TH2F *h2d	 =  (TH2F*)fA->Get("PtMCProtonP");
   TH2F *h3d	 =  (TH2F*)fA->Get("PtMCKaonP");
   h1 -> Sumw2();
   h2 -> Sumw2();
   h3 -> Sumw2();
   h1d -> Sumw2();
   h2d -> Sumw2();
   h3d -> Sumw2();

   TH2F *h4	 =  (TH2F*)fA->Get("PtRecoPionS");
   TH2F *h5 	 =  (TH2F*)fA->Get("PtRecoProtonS");
   TH2F *h6 	 =  (TH2F*)fA->Get("PtRecoKaonS");
   TH2F *h4d	 =  (TH2F*)fA->Get("PtMCPionS");
   TH2F *h5d	 =  (TH2F*)fA->Get("PtMCProtonS");
   TH2F *h6d	 =  (TH2F*)fA->Get("PtMCKaonS");
   h4 -> Sumw2();
   h5 -> Sumw2();
   h6 -> Sumw2();
   h4d -> Sumw2();
   h5d -> Sumw2();
   h6d -> Sumw2();

   TH2F *hA 	 =  new TH2F("", "", 100, -3, 3, 100, 0, 4);
   TH2F *hB 	 =  new TH2F("", "", 100, -3, 3, 100, 0, 4);
   TH2F *hC 	 =  new TH2F("", "", 100, -3, 3, 100, 0, 4);
   TH2F *hD 	 =  new TH2F("", "", 100, -3, 3, 100, 0, 4);
   TH2F *hE 	 =  new TH2F("", "", 100, -3, 3, 100, 0, 4);
   TH2F *hF 	 =  new TH2F("", "", 100, -3, 3, 100, 0, 4);


	// Crear canvas
   TCanvas *c1 = new TCanvas("c1","", 1000, 800);
   c1 -> Divide(2,1); // (Columnas, Filas)

   TCanvas *c2 = new TCanvas("c2","", 1000, 800);
   c2 -> Divide(2,1); // (Columnas, Filas)

   TCanvas *c3 = new TCanvas("c3","", 1000, 800);
   c3 -> Divide(2,1); // (Columnas, Filas)

   TCanvas *cd = new TCanvas("CDivide","",1000,800);
   cd -> Divide(3,1);

	// Dibujar histogramas para Protones
   hA = (TH2F*)h1->Clone();
   hA -> Sumw2();
   hA -> Divide(h1d);
   hA -> SetOption("E1");
   hA -> SetTitle(" p_{T}^{RECO} / p_{T}^{MC} #pi Primary ");
   c1 -> cd(1);
   hA -> Draw();
   
   hB = (TH2F*)h4->Clone();
   hB -> Sumw2();
   hB -> Divide(h4d);
   hB -> SetOption("E1");
   hB -> SetTitle(" p_{T}^{RECO} / p_{T}^{MC} #pi Secondary ");
   c1 -> cd(2);
   hB -> Draw();

   c1 -> Update();
   
   hC = (TH2F*)h2->Clone();
   hC -> Sumw2();
   hC -> Divide(h2d);
   hC -> SetOption("E1");
   hC -> SetTitle(" p_{T}^{RECO} / p_{T}^{MC} p Primary ");
   c2 -> cd(1);
   hC -> Draw();
   
   hD = (TH2F*)h5->Clone();
   hD -> Sumw2();
   hD -> Divide(h5d);
   hD -> SetOption("E1");
   hD -> SetTitle(" p_{T}^{RECO} / p_{T}^{MC} p Secondary ");
   c2 -> cd(2);
   hD -> Draw();
 
   c2 -> Update();
   
   hE = (TH2F*)h3->Clone();
   hE -> Sumw2();
   hE -> Divide(h3d);
   hE -> SetOption("E1");
   hE -> SetTitle(" p_{T}^{RECO} / p_{T}^{MC} #kappa Primary ");
   c3 -> cd(1);
   hE -> Draw();
   
   hF = (TH2F*)h6->Clone();
   hF -> Sumw2();
   hF -> Divide(h6d);
   hF -> SetOption("E1");
   hF -> SetTitle(" p_{T}^{RECO} / p_{T}^{MC} #kappa Secondary ");
   c3 -> cd(2);
   hF -> Draw();

   c3 -> Update();

	// Joint the primary and secondary particles in the same histogram
   cd -> cd(1);
   hA -> SetLineColor(kRed);
   hA -> SetTitle(" Comparative p_{T}^{RECO} / p_{T}^{MC}  Primary and Secondary #pi ");
   hA -> Draw();

   hB -> SetLineColor(kBlue);
   hB -> Draw("SAME");

   auto legend1 = new TLegend(0.1,0.7,0.48,0.9);
   legend1 -> AddEntry(hA," Primary #pi ");
   legend1 -> AddEntry(hB," Secondary #pi ");
   legend1 -> Draw("SAME");

   cd -> cd(2);
   hC -> SetLineColor(kRed);
   hC -> SetTitle(" Comparative p_{T}^{RECO} / p_{T}^{MC} Primary and Secondary p ");
   hC -> Draw();

   hD -> SetLineColor(kBlue);
   hD -> Draw("SAME");

   auto legend2 = new TLegend(0.1,0.7,0.48,0.9);
   legend2 -> AddEntry(hA," Primary #p ");
   legend2 -> AddEntry(hB," Secondary #p ");
   legend2 -> Draw("SAME");

   cd -> cd(3);
   hE -> SetLineColor(kRed);
   hE -> SetTitle(" Comparative p_{T}^{RECO} / p_{T}^{MC} Primary and Secondary #kappa ");
   hE -> Draw();

   hF -> SetLineColor(kBlue);
   hF -> Draw("SAME");

   auto legend3 = new TLegend(0.1,0.7,0.48,0.9);
   legend3 -> AddEntry(hA," Primary #kappa ");
   legend3 -> AddEntry(hB," Secondary #kapp ");
   legend3 -> Draw("SAME");

   cd -> Update();
}
