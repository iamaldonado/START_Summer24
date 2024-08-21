//Bethe-Bloch equations for Helium 3
Double_t parHe3(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 7.890), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 7.890), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 2.809), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}

 void FitHe3nocut(){
   //__mean values of Gaussian functions and their errors 
    const int nBins=17;
    Stat_t data1[nBins]={52.3567,48.3904,44.4226,41.8972,37.7458,34.9401,32.8141,30.9003,28.0807,26.9465,26.6242,25.8119,23.435,23.3145,22.0309,21.344,20.3264};

    const int nError=17;
    Stat_t data2[nError]={0.745219,0.762677,0.383833,0.181101,0.170119,0.708364,0.79278,0.28176,1.25059,0.348681,0.143361,0.131917,0.394502,0.122222,0.680441,0.681122,0.226685};

    //Histogram dE/dx
     TFile *file1 = new TFile("/home/alejandro/Documentos/EnerClass/Fit functions/taskEnerAll.root");
     TH2F *h_dedx_He3 = (TH2F*)file1->Get("h_dedx_He^{3}");
     TH1F *fa1 = new TH1F("dEdx_Helium3", "dE/dx parameterization for the Helium3; p*q GeV/c; dE/dx arb.units", 17, 0.260, 0.750);

     for(int i=0; i<nBins; i++){

     fa1->SetBinContent(i+1,data1[i]);
       for (int j=0; j<nError; j++){
         fa1->SetBinError(j+1,data2[j]);
       }
     }

     TF1 *fitparHe3 = new TF1("fitparHe3",parHe3,0.260,0.750,5);

     
     h_dedx_He3->Draw("same");
     fa1->Draw("same");
     fa1->Fit("fitparHe3");

 }
