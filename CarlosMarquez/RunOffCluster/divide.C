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

   TH2F *h4	 =  (TH2F*)fA->Get("PtRecoPionS");
   TH2F *h5 	 =  (TH2F*)fA->Get("PtRecoProtonS");
   TH2F *h6 	 =  (TH2F*)fA->Get("PtRecoKaonS");
   TH2F *h4d	 =  (TH2F*)fA->Get("PtMCPionS");
   TH2F *h5d	 =  (TH2F*)fA->Get("PtMCProtonS");
   TH2F *h6d	 =  (TH2F*)fA->Get("PtMCKaonS");

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

	// Dibujar histogramas para Protones
   hA = (TH2F*)h1->Clone();
   hA -> Divide(h1d);
   hA -> SetTitle(" p_{T}^{RECO}/p_{T}^{MC} #pi Secondary ");
   c1 -> cd(1);
   hA -> Draw();
   
   hB = (TH2F*)h4->Clone();
   hB -> Divide(h4d);
   hB -> SetTitle(" p_{T}^{RECO}/p_{T}^{MC} #pi Secondary ");
   c1 -> cd(2);
   hB -> Draw();

   hC = (TH2F*)h2->Clone();
   hC -> Divide(h2d);
   hC -> SetTitle(" p_{T}^{RECO}/p_{T}^{MC} p Primary ");
   c2 -> cd(1);
   hC -> Draw();
   
   hD = (TH2F*)h5->Clone();
   hD -> Divide(h5d);
   hD -> SetTitle(" p_{T}^{RECO} / p_{T}^{MC} p Secondary ");
   c2 -> cd(2);
   hD -> Draw();
 
   hE = (TH2F*)h3->Clone();
   hE -> Divide(h3d);
   hE -> SetTitle(" p_{T}^{RECO} / p_{T}^{MC} #kappa Primary ");
   c3 -> cd(1);
   hE -> Draw();
   
   hF = (TH2F*)h6->Clone();
   hF -> Divide(h6d);
   //hF -> SetOption("E1");
   hF -> SetTitle(" p_{T}^{RECO} / p_{T}^{MC} #kappa Secondary ");
   c3 -> cd(2);
   hF -> Draw();

   c1 -> Update();
}
