//Bethe-Bloch equations for tritium 
Double_t parT(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 7.89), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 7.89), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 2.81), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}
 void FitTrinocut(){
    //__mean values of Gaussian functions and their errors 
     const int nBins=18;
    Stat_t data1[nBins]={28.6264,25.5555,23.3364,20.95,18.652,18.6571,15.172,14.1364,12.7,11.4221,10.3551,9.58458,8.85554,8.0735,7.7786,7.4989,6.65372,6.20932};

    const int nError=18;
    Stat_t data2[nError]={0.115406,0.144772,0.0641441,0.0775749,0.142569,0.717905,0.276563,0.0485215,0.0636811,0.146642,0.119405,0.398137,0.0452,0.0660899,0.408914,0.364173,0.0518328,0.0464529};

    //Histogram dE/dx
     TFile *file1 = new TFile("/home/alejandro/Documentos/EnerClass/Fit functions/taskEnerAll.root");
     TH2F *h_dedx_t = (TH2F*)file1->Get("h_dedx_t");
     TH1F *fa1 = new TH1F("dEdx_tritium", "dE/dx parameterization for the tritium; p*q GeV/c; dE/dx arb.units", 18, 0.39, 1.28);

     for(int i=0; i<nBins; i++){

     fa1->SetBinContent(i+1,data1[i]);
       for (int j=0; j<nError; j++){
         fa1->SetBinError(j+1,data2[j]);
       }
     }

     TF1 *fitparT = new TF1("fitparT",parT,0.39,1.28,5);

     h_dedx_t->Draw("same");
     fa1->Draw("same");
     fa1->Fit("fitparT");

 }
