//Bethe-Bloch equations for Deuterium
Double_t parD(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 3.52), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 3.52), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 1.876), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}
 void FitDnocut(){
   //__mean values of Gaussian functions and their errors
     const int nBins=16;
    Stat_t data1[nBins]={21.0392,17.8186,15.3921,13.2319,11.4434,9.94222,8.75991,7.81609,6.95,6.28567,5.74349,5.29279,4.84719,4.48742,4.18554,3.9182};

    const int nError=16;
    Stat_t data2[nError]={0.0737788,0.140581,0.0444482,0.0244504,0.0343747,0.0422325,0.0147567,0.0337827,0.0192728,0.0150906,0.00990263,0.0102437,0.0106329,0.00857485,0.00713764,0.00620705};

    //Histogram dE/dx
     TFile *file1 = new TFile("/home/alejandro/Documentos/EnerClass/Fit functions/taskEnerAll.root");
     TH2F *h_dedx_d = (TH2F*)file1->Get("h_dedx_d");
     TH1F *fa1 = new TH1F("dEdx_deuterium", "dE/dx parameterization for the deuterium; p*q GeV/c; dE/dx arb.units", 16, 0.32, 1.16);

     for(int i=0; i<nBins; i++){

     fa1->SetBinContent(i+1,data1[i]);
       for (int j=0; j<nError; j++){
         fa1->SetBinError(j+1,data2[j]);
       }
     }

     TF1 *fitparD = new TF1("fitparD",parD,0.32,1.16,5);

     h_dedx_d->Draw("same");
     fa1->Draw("same");
     fa1->Fit("fitparD");


 }
