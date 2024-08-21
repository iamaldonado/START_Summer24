//Bethe-Bloch equations for Helium 4
Double_t parHe4(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 13.898), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 13.898), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 3.728), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}

 void FitHe4nocut(){
   //__mean values of Gaussian functions and their errors 
    const int nBins=18;
    Stat_t data1[nBins]={54.3364,51.4505,47.0809,45.0653,41.4086,40.3827,37.2398,34.9885,33.9495,32.6539,31.5679,29.5501,28.2179,27.5399,26.6678,25.1155,25.9452,22.1674};

    const int nError=18;
    Stat_t data2[nError]={0.798889,0.534871,0.352881,0.307755,0.328484,0.1532,0.143935,0.170712,0.1976,0.146428,0.156418,0.169914,0.190561,0.198644,0.430992,0.242638,1.44397,0.768671};

    //Histogram dE/dx
     TFile *file1 = new TFile("/home/alejandro/Documentos/EnerClass/Fit functions/taskEnerAll.root");
     TH2F *h_dedx_He4 = (TH2F*)file1->Get("h_dedx_He^{4}");
     TH1F *fa1 = new TH1F("dEdx_Helium4", "dE/dx parameterization for the Helium4; p*q GeV/c; dE/dx arb.units", 18, 0.350, 0.900);

     for(int i=0; i<nBins; i++){

     fa1->SetBinContent(i+1,data1[i]);
       for (int j=0; j<nError; j++){
         fa1->SetBinError(j+1,data2[j]);
       }
     }

     TF1 *fitparHe4 = new TF1("fitparHe4",parHe4,0.350,0.900,5);

     h_dedx_He4->Draw("same");
     fa1->Draw("same");
     fa1->Fit("fitparHe4");

 }
