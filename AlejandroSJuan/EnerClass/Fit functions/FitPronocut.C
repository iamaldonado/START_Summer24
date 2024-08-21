//Bethe-Bloch equation for protons
Double_t parPr(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.88), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.88), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 0.9383), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}
 void FitPronocut(){
 
    //Histogram dE/dx
     TFile *file1 = new TFile("/home/alejandro/Documentos/EnerClass/Fit functions/taskEnerAll.root");
     TH2F *h_dedx_p = (TH2F*)file1->Get("h_dedx_p");

 //__mean values of Gaussian functions and their errors
    const int nBins=18;
    Stat_t data1[nBins]={24.9622,17.9635,13.1672,9.88914,7.81359,6.31973,5.25655,4.47178,3.88242,3.43243,3.07454,2.78969,2.56354,2.37552,2.22018,2.09064,1.98205,1.89001};

    const int nError=18;
    Stat_t data2[nError]={0.0258993,0.0148647,0.0167922,0.00605115,0.0039364,0.00293919,0.00411013,0.00169941,0.0015333,0.000743479,0.000772767,0.000642989,0.00035937,0.000349082,0.000448932,0.000292752,0.000237986,0.000234231};

     TH1F *fa1 = new TH1F("dEdx_proton", "dE/dx parameterization for the proton; p*q GeV/c; dE/dx arb.units", 18, 0.13, 1.09);

     for(int i=0; i<nBins; i++){

     fa1->SetBinContent(i+1,data1[i]);
       for (int j=0; j<nError; j++){
         fa1->SetBinError(j+1,data2[j]);
       }
     }

     TF1 *fitparPr = new TF1("fitparPr",parPr,0.13,1.09,5);


     h_dedx_p->Draw("same");
     fa1->Draw("same");
     fa1->Fit("fitparPr");


 }
