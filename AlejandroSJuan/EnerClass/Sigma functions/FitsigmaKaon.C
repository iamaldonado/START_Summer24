//Bethe-Bloch equation for kaon
Double_t parK(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.2437), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.2437), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 0.4937), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}

void FitsigmaKaon(){

    //Histogram dE/dx
     TFile *file1 = new TFile("/home/alejandro/Documentos/EnerClass/Fit functions/taskEnerAll.root");
     TH2F *h_dedx_kp = (TH2F*)file1->Get("h_dedx_kp");

//__Kaon+1sima_Bins___________________________________________________________________
    const int nBinsKplus1=12;
    Stat_t data1Kplus1[nBinsKplus1]={5.367451,4.032883,3.251013,2.745502,2.431266,2.146166,1.968215,1.831212,1.729213,1.653182,1.585386,1.542418};

//__Kaon-1sima_Bins___________________________________________________________________
    const int nBinsKminus1=12;
    Stat_t data1Kminus1[nBinsKminus1]={3.922109,3.140237,2.607687,2.242578,1.984874,1.790134,1.644685,1.536308,1.451647,1.387518,1.336074,1.295322};

//__Kaon+2sima_Bins___________________________________________________________________
    const int nBinsKplus2=12;
    Stat_t data1Kplus2[nBinsKplus2]={6.090122,4.479206,3.572676,2.996964,2.654462,2.324182,2.12998,1.978664,1.867996,1.786014,1.710042,1.665966};

//__Kaon-2sima_Bins___________________________________________________________________
    const int nBinsKminus2=12;
    Stat_t data1Kminus2[nBinsKminus2]={3.199438,2.693914,2.286024,1.991116,1.761678,1.612118,1.48292,1.388856,1.312864,1.254686,1.211418,1.171774};

//__errors__________________________________________________________
    const int nError=12;
    Stat_t data2[nError]={0.01954147,0.01083589,0.0094293,0.00256685,0.00294743,0.000784369,0.000746148,0.000509415,0.000461023,0.000290335,0.000303722,0.000254303};

//__Arrangement_for_kaon+1sigma_parameters_______________________________________________ 
     TH1F *fa1 = new TH1F("dEdx_kaon+1sigma", "dE/dx parameterization for the kaon+1sigma; p*q GeV/c; dE/dx arb.units", 12, 0.235, 0.83);

     for(int i=0; i<nBinsKplus1; i++){

     fa1->SetBinContent(i+1,data1Kplus1[i]);
       for (int j=0; j<nError; j++){
         fa1->SetBinError(j+1,data2[j]);
       }
     }

//__Arrangement_for_kaon-1sigma_parameters_______________________________________________ 
     TH1F *fa2 = new TH1F("dEdx_kaon-1sigma", "dE/dx parameterization for the kaon-1sigma; p*q GeV/c; dE/dx arb.units", 12, 0.235, 0.83);

     for(int i=0; i<nBinsKminus1; i++){

     fa2->SetBinContent(i+1,data1Kminus1[i]);
       for (int j=0; j<nError; j++){
         fa2->SetBinError(j+1,data2[j]);
       }
     }

     TF1 *fitparK = new TF1("fitparK",parK,0.235,0.83,5);
    
     h_dedx_kp->Draw("same");
    
     fa1->Draw("same");
     fa2->Draw("same");

     fitparK->SetLineColor(2);
     fa1->Fit("fitparK");
     printf("Parameters for the Bethe-Bloch function of the Kaon+1sigma\n");
     fa2->Fit("fitparK");
     printf("Parameters for the Bethe-Bloch function of the Kaon-1sigma\n");


}
