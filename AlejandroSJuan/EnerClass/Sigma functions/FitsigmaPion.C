//Bethe-Bloch equation for Pion 
Double_t parPi(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.01949), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.01949), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] /0.1396), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}

void FitsigmaPion(){

    //Histogram dE/dx
     TFile *file1 = new TFile("/home/alejandro/Documentos/EnerClass/Fit functions/taskEnerAll.root");
     TH2F *h_dedx_pip = (TH2F*)file1->Get("h_dedx_pip");

//__Pion+1sima_Bins___________________________________________________________________
    const int nBinsPiplus1=16;
    Stat_t data1Piplus1[nBinsPiplus1]={5.637703,3.377253,2.324922,1.899806,1.67996,1.545588,1.449198,1.392983,1.361642,1.3282775,1.3128691,1.3033808,1.2976517,1.2946204,1.2937497,1.294074};

//__Pion-1sigma_Bins___________________________________________________________________
    const int nBinsPiminus1=16;
    Stat_t data1Piminus1[nBinsPiminus1]={4.036817,1.960527,1.879838,1.582294,1.40766,1.297992,1.226482,1.182377,1.153418,1.1297025,1.1169909,1.1088592,1.1038483,1.1011996,1.1002903,1.100546};

//__Pion+2sigma_Bins___________________________________________________________________
    const int nBinsPiplus2=16;
    Stat_t data1Piplus2[nBinsPiplus2]={6.438146,4.085616,2.547464,2.058562,1.81611,1.669386,1.560556,1.498286,1.465754,1.427565,1.4108082,1.4006416,1.3945534,1.3913308,1.3904794,1.390838};

//__Pion-2sigma_Bins___________________________________________________________________
    const int nBinsPiminus2=16;
    Stat_t data1Piminus2[nBinsPiminus2]={3.236374,1.252164,1.657296,1.423538,1.27151,1.174194,1.115124,1.077074,1.049306,1.030415,1.0190518,1.0115984,1.0069466,1.0044892,1.0035606,1.003782};


//__errors__________________________________________________________
    const int nError=16;
    Stat_t data2[nError]={0.01843353,0.0334273,0.000674437,0.000309777,0.0001042,7.27723E-05,0.000133965,7.15084E-05,5.62129E-05,0.000103032,8.54249E-05,0.000078003,7.51629E-05,7.49232E-05,7.68819E-05,8.00301E-05};

//__Arrangement_for_pion+1sigma_parameters_______________________________________________ 
    
    TH1F *fa1 = new TH1F("dEdx_pion+1sigma", "dE/dx parameterization for the pion+1sigma; p*q GeV/c; dE/dx arb.units", 16, 0.069, 0.625);
    
    for(int i=0; i<nBinsPiplus1; i++){

    fa1->SetBinContent(i+1,data1Piplus1[i]);
      for (int j=0; j<nError; j++){
        fa1->SetBinError(j+1,data2[j]);
      }
    }

//__Arrangement_for_pion-1sigma_parameters_______________________________________________ 
    
    TH1F *fa2 = new TH1F("dEdx_pion-1sigma", "dE/dx parameterization for the pion-1sigma; p*q GeV/c; dE/dx arb.units", 16, 0.069, 0.625);

    for(int i=0; i<nBinsPiminus1; i++){

    fa2->SetBinContent(i+1,data1Piminus1[i]);
      for (int j=0; j<nError; j++){
        fa2->SetBinError(j+1,data2[j]);
      }
    }

    TF1 *fitparPi = new TF1("fitparPi",parPi,0.069,0.625,5);
    
    h_dedx_pip->Draw("same");
    
    fa1->Draw("same");
    fa2->Draw("same");
   
    fitparPi->SetLineColor(2); 
    fa1->Fit("fitparPi");
    printf("Parameters for the Bethe-Bloch function of the Pion+1sigma\n");
    fa2->Fit("fitparPi");
    printf("Parameters for the Bethe-Bloch function of the Pion-1sigma\n");

}
