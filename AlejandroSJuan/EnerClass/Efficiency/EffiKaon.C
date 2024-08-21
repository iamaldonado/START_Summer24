void EffiKaon(){

  gStyle->SetPalette(1);
  gStyle->SetCanvasColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);

TFile *f1 = TFile::Open("/home/alejandro/Documentos/EnerClass/simplept/macros/taskEner.root");
//razon con el corte unicamente en dedx
TH1F *hptre = (TH1F*)f1->Get("PtKaonreco1");
TH1F *hptmc = (TH1F*)f1->Get("PtKaonMC1");

hptre->Sumw2();
hptmc->Sumw2();

TH1F *hEff = (TH1F*) hptre->Clone("hEff");
hEff->Divide(hptmc);

TCanvas* c1 = new TCanvas("c1","Pt_Efficiency_Kaon",800,800);

TPad *pad = new TPad("pad1","pad1",0.01,0.01,0.99,0.99,0);

pad->SetLeftMargin(0.15);
pad->SetRightMargin(0.015);
pad->SetTopMargin(0.01);
pad->SetBottomMargin(0.12);
pad->SetFillColor(0);
pad->SetGridy();

pad->Draw();
pad->cd();
TH1F *arm = new TH1F("arm","arm",200,0,5);

        arm->GetYaxis()->SetTitleSize(0.05);
        arm->GetXaxis()->SetTitleSize(0.05);
        arm->GetXaxis()->SetTitleOffset(1.1);
        arm->GetYaxis()->SetTitleOffset(1.2);
        arm->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        arm->GetXaxis()->SetLabelFont(42);
        arm->GetXaxis()->SetLabelSize(0.04);

        arm->GetXaxis()->SetRangeUser(0.0,2.0);
        arm->GetYaxis()->SetRangeUser(0.0,3.0);
        arm->GetYaxis()->SetTitle("Efficiency");
        arm->GetYaxis()->SetLabelFont(42);
        arm->GetYaxis()->SetLabelSize(0.04);
        arm->GetYaxis()->SetTitleFont(12);
        arm->GetXaxis()->SetTitleFont(12);
        arm->Draw();

hEff->Draw("sames");
}
