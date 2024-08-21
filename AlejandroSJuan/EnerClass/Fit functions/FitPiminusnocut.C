//Bethe-Bloch equations for AntiPions
Double_t parPiminus(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power((-1)*x[0] / TMath::Sqrt(x[0] * x[0] + 0.01949), p[3]);
   x2  = p[1] - TMath::Power((-1)*x[0] / TMath::Sqrt(x[0] * x[0] + 0.01949), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / ((-1)*x[0] /0.1396), p[4]));
   ans = (-1)* x1 * (x2 - x3);

   return ans;
}
 void FitPiminusnocut(){
 //__mean values of Gaussian functions and their errors
    const int nBins=16;
    Stat_t data1[nBins]={1.19731,1.19702,1.19791,1.20075,1.20612,1.21493,1.22899,1.25753,1.28768,1.33784,1.42179,1.54381,1.74105,2.10238,2.66889,4.83726};

    const int nError=16;
    Stat_t data2[nError]={4.63751E-05,4.45443E-05,4.34844E-05,4.38774E-05,4.60078E-05,5.13052E-05,6.38093E-05,0.000029952,0.00004054,8.46954E-05,3.31638E-05,0.0000483826,0.000189425,0.000327029,0.0161612,0.00430803};


    //Histogram dE/dx
     TFile *file1 = new TFile("/home/alejandro/Documentos/EnerClass/Fit functions/taskEnerAll.root");
     TH2F *h_dedx_pim = (TH2F*)file1->Get("h_dedx_pim");
     TH1F *fa1 = new TH1F("dEdx_pionminus", "dE/dx parameterization for the pionminus; p*q GeV/c; dE/dx arb.units", 16,-0.625,-0.060);

     for(int i=0; i<nBins; i++){

     fa1->SetBinContent(i+1,data1[i]);
       for (int j=0; j<nError; j++){
         fa1->SetBinError(j+1,data2[j]);
       }
     }

     TF1 *fitparPiminus = new TF1("fitparPiminus",parPiminus,-0.625,-0.060,5);

     h_dedx_pim->Draw("same");
     fa1->Draw("same");
     fa1->Fit("fitparPiminus");

 }
