//Bethe-Bloch equation for kaons
Double_t parK(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.2437), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.2437), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 0.4937), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}
 void FitKnocut(){
 //__mean values of Gaussian functions and their errors
     const int nBins=12;
    Stat_t data1[nBins]={4.64478,3.58656,2.92935,2.49404,2.20807,1.96815,1.80645,1.68376,1.59043,1.52035,1.46073,1.41887};

    const int nError=12;
    Stat_t data2[nError]={0.00379767,0.00426258,0.00549599,0.00113166,0.0018749,0.000359293,0.000411637,0.000266968,0.000250994,0.000126489,0.000154411,0.000115276};

    //Histogram dE/dx
     TFile *file1 = new TFile("/home/alejandro/Documentos/EnerClass/Fit functions/taskEnerAll.root");
     TH2F *h_dedx_kp = (TH2F*)file1->Get("h_dedx_kp");
     TH1F *fa1 = new TH1F("dEdx_kaon", "dE/dx parameterization for the kaon; p*q GeV/c; dE/dx arb.units", 12, 0.235, 0.83);

     for(int i=0; i<nBins; i++){

     fa1->SetBinContent(i+1,data1[i]);
       for (int j=0; j<nError; j++){
         fa1->SetBinError(j+1,data2[j]);
       }
     }

     TF1 *fitparK = new TF1("fitparK",parK,0.235,0.83,5);


     h_dedx_kp->Draw("same");
     fa1->Draw("same");
     fa1->Fit("fitparK");

 }
