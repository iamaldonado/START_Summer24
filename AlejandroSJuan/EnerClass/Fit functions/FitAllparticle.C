
//Bethe_Bloch_functions_for_all_particles
//the mass is in units of GeV/c^{2}
Double_t parPr(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.88), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.88), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 0.9383), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}

Double_t parAntiPr(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power((-1)*x[0] / TMath::Sqrt(x[0] * x[0] + 0.88), p[3]);
   x2  = p[1] - TMath::Power((-1)*x[0] / TMath::Sqrt(x[0] * x[0] + 0.88), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / ((-1)*x[0] / 0.9383), p[4]));
   ans =(-1)* x1 * (x2 - x3);

   return ans;
}

Double_t parK(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.2437), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.2437), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 0.4937), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}

Double_t parKminus(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power((-1)*x[0] / TMath::Sqrt(x[0] * x[0] + 0.2437), p[3]);
   x2  = p[1] - TMath::Power((-1)*x[0] / TMath::Sqrt(x[0] * x[0] + 0.2437), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / ((-1)*x[0] / 0.4937), p[4]));
   ans = (-1)* x1 * (x2 - x3);

   return ans;
}

Double_t parPi(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.01949), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 0.01949), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] /0.1396), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}

Double_t parPiminus(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power((-1)*x[0] / TMath::Sqrt(x[0] * x[0] + 0.01949), p[3]);
   x2  = p[1] - TMath::Power((-1)*x[0] / TMath::Sqrt(x[0] * x[0] + 0.01949), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / ((-1)*x[0] /0.1396), p[4]));
   ans = (-1)* x1 * (x2 - x3);

   return ans;
}

Double_t parD(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 3.52), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 3.52), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 1.876), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}

Double_t parT(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 7.89), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 7.89), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 2.81), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}

Double_t parHe3(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 7.890), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 7.890), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 2.809), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}

Double_t parHe4(Double_t *x, Double_t *p)
{
   Double_t x1, x2, x3,ans;

   x1  = p[0] / TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 13.898), p[3]);
   x2  = p[1] - TMath::Power(x[0] / TMath::Sqrt(x[0] * x[0] + 13.898), p[3]);
   x3  = TMath::Log(p[2] + TMath::Power(1.0 / (x[0] / 3.728), p[4]));
   ans = x1 * (x2 - x3);

   return ans;
}

void FitAllparticle(){
//__Histogram_for_all_particles___________________________________________________
     TFile *file1 = new TFile("/home/alejandro/Documentos/EnerClass/Fit functions/taskEnerAll.root");
     TH2F *h__dedx = (TH2F*)file1->Get("h__dedx");

//__Proton_Bins___________________________________________________________________
    const int nBinsPro=18;
    Stat_t data1Pro[nBinsPro]={24.9622,17.9635,13.1672,9.88914,7.81359,6.31973,5.25655,4.47178,3.88242,3.43243,3.07454,2.78969,2.56354,2.37552,2.22018,2.09064,1.98205,1.89001};

    const int nErrorPro=18;
    Stat_t data2Pro[nErrorPro]={0.0258993,0.0148647,0.0167922,0.00605115,0.0039364,0.00293919,0.00411013,0.00169941,0.0015333,0.000743479,0.000772767,0.000642989,0.00035937,0.000349082,0.000448932,0.000292752,0.000237986,0.000234231};
//__Anti_proton_Bins______________________________________________________________
    const int nBinsAP=18;
    Stat_t data1AP[nBinsAP]={1.89001,1.98205,2.09064,2.22018,2.37552,2.56354,2.78969,3.07454,3.43243,3.88242,4.47178,5.25655,6.31973,7.81359,9.88914,13.1672,17.9635,24.9622};

    const int nErrorAP=18;
    Stat_t data2AP[nErrorAP]={0.000234231,0.000237986,0.000292752,0.000448932,0.000349082,0.00035937,0.000642989,0.000772767,0.000743479,0.0015333,0.00169941,0.00411013,0.00293919,0.0039364,0.00605115,0.0167922,0.0148647,0.0258993};
//__Kaon+_Bins______________________________________________________________________
    const int nBinsK=12;
    Stat_t data1K[nBinsK]={4.64478,3.58656,2.92935,2.49404,2.20807,1.96815,1.80645,1.68376,1.59043,1.52035,1.46073,1.41887};

    const int nErrorK=12;
    Stat_t data2K[nErrorK]={0.00379767,0.00426258,0.00549599,0.00113166,0.0018749,0.000359293,0.000411637,0.000266968,0.000250994,0.000126489,0.000154411,0.000115276};

//__Kaon-_Bins______________________________________________________________________
    const int nBinsAK=12;
    Stat_t data1AK[nBinsAK]={1.41887,1.46073,1.52035,1.59043,1.68376,1.80645,1.96815,2.20807,2.49404,2.92935,3.58656,4.64478};

    const int nErrorAK=12;
    Stat_t data2AK[nErrorAK]={0.000115276,0.000154411,0.000126489,0.000250994,0.000266968,0.000411637,0.000359293,0.0018749,0.00113166,0.00549599,0.00426258,0.00379767};
    
//__Pion+_Bins______________________________________________________________________
    const int nBinsPi=16;
    Stat_t data1Pi[nBinsPi]={4.83726,2.66889,2.10238,1.74105,1.54381,1.42179,1.33784,1.28768,1.25753,1.22899,1.21493,1.20612,1.20075,1.19791,1.19702,1.19731};

    const int nErrorPi=16;
    Stat_t data2Pi[nErrorPi]={0.00430803,0.0161612,0.000327029,0.000189425,0.0000483826,3.31638E-05,8.46954E-05,0.00004054,0.000029952,6.38093E-05,5.13052E-05,4.60078E-05,4.38774E-05,4.34844E-05,4.45443E-05,4.63751E-05};

//__Pion-_Bins______________________________________________________________________
    const int nBinsAPi=16;
    Stat_t data1APi[nBinsAPi]={1.19731,1.19702,1.19791,1.20075,1.20612,1.21493,1.22899,1.25753,1.28768,1.33784,1.42179,1.54381,1.74105,2.10238,2.66889,4.83726};

    const int nErrorAPi=16;
    Stat_t data2APi[nErrorAPi]={4.63751E-05,4.45443E-05,4.34844E-05,4.38774E-05,4.60078E-05,5.13052E-05,6.38093E-05,0.000029952,0.00004054,8.46954E-05,3.31638E-05,0.0000483826,0.000189425,0.000327029,0.0161612,0.00430803};

//__Deuterium_Bins__________________________________________________________________
    const int nBinsD=16;
    Stat_t data1D[nBinsD]={21.0392,17.8186,15.3921,13.2319,11.4434,9.94222,8.75991,7.81609,6.95,6.28567,5.74349,5.29279,4.84719,4.48742,4.18554,3.9182};

    const int nErrorD=16;
    Stat_t data2D[nErrorD]={0.0737788,0.140581,0.0444482,0.0244504,0.0343747,0.0422325,0.0147567,0.0337827,0.0192728,0.0150906,0.00990263,0.0102437,0.0106329,0.00857485,0.00713764,0.00620705};

//__Tritium_Bins____________________________________________________________________
    const int nBinsT=18;
    Stat_t data1T[nBinsT]={28.6264,25.5555,23.3364,20.95,18.652,18.6571,15.172,14.1364,12.7,11.4221,10.3551,9.58458,8.85554,8.0735,7.7786,7.4989,6.65372,6.20932};

    const int nErrorT=18;
    Stat_t data2T[nErrorT]={0.115406,0.144772,0.0641441,0.0775749,0.142569,0.717905,0.276563,0.0485215,0.0636811,0.146642,0.119405,0.398137,0.0452,0.0660899,0.408914,0.364173,0.0518328,0.0464529};

//__He3_Bins________________________________________________________________________
    const int nBinsHe3=17;
    Stat_t data1He3[nBinsHe3]={52.3567,48.3904,44.4226,41.8972,37.7458,34.9401,32.8141,30.9003,28.0807,26.9465,26.6242,25.8119,23.435,23.3145,22.0309,21.344,20.3264};

    const int nErrorHe3=17;
    Stat_t data2He3[nErrorHe3]={0.745219,0.762677,0.383833,0.181101,0.170119,0.708364,0.79278,0.28176,1.25059,0.348681,0.143361,0.131917,0.394502,0.122222,0.680441,0.681122,0.226685};

//__He4_Bins________________________________________________________________________
    const int nBinsHe4=18;
    Stat_t data1He4[nBinsHe4]={54.3364,51.4505,47.0809,45.0653,41.4086,40.3827,37.2398,34.9885,33.9495,32.6539,31.5679,29.5501,28.2179,27.5399,26.6678,25.1155,25.9452,22.1674};

    const int nErrorHe4=18;
    Stat_t data2He4[nErrorHe4]={0.798889,0.534871,0.352881,0.307755,0.328484,0.1532,0.143935,0.170712,0.1976,0.146428,0.156418,0.169914,0.190561,0.198644,0.430992,0.242638,1.44397,0.768671};



//__Arrangement_for_proton_parameters_______________________________________________ 
    TH1F *fa1 = new TH1F("dEdx_proton", "dE/dx parameterization for the proton; p*q GeV/c; dE/dx arb.units", 18, 0.13, 1.09);

    for(int i=0; i<nBinsPro; i++){

    fa1->SetBinContent(i+1,data1Pro[i]);
      for (int j=0; j<nErrorPro; j++){
        fa1->SetBinError(j+1,data2Pro[j]);
      }
    }
//__Arrangement_for_Antiproton_parameters___________________________________________     
    TH1F *fa2 = new TH1F("dEdx_Antiproton", "dE/dx parameterization for the Antiproton; p*q GeV/c; dE/dx arb.units", 18,-1.09,-0.148);

    for(int i=0; i<nBinsAP; i++){

    fa2->SetBinContent(i+1,data1AP[i]);
      for (int j=0; j<nErrorAP; j++){
        fa2->SetBinError(j+1,data2AP[j]);
      }
    }

//__Arrangement_for_Kaon+_parameters________________________________________________
    TH1F *fa3 = new TH1F("dEdx_kaon", "dE/dx parameterization for the kaon; p*q GeV/c; dE/dx arb.units", 12, 0.235, 0.83);

    for(int i=0; i<nBinsK; i++){

    fa3->SetBinContent(i+1,data1K[i]);
      for (int j=0; j<nErrorK; j++){
        fa3->SetBinError(j+1,data2K[j]);
      }
    }
//__Arrangement_for_Kaon-_parameters________________________________________________
     TH1F *fa4 = new TH1F("dEdx_kaonminus", "dE/dx parameterization for the kaonminus; p*q GeV/c; dE/dx arb.units", 12,-0.83,-0.255);

     for(int i=0; i<nBinsAK; i++){

     fa4->SetBinContent(i+1,data1AK[i]);
       for (int j=0; j<nErrorAK; j++){
         fa4->SetBinError(j+1,data2AK[j]);
       }
     }

//__Arrangement_for_Pion+_parameters________________________________________________
     TH1F *fa5 = new TH1F("dEdx_pion", "dE/dx parameterization for the pion; p*q GeV/c; dE/dx arb.units", 16, 0.069, 0.625);

     for(int i=0; i<nBinsPi; i++){

     fa5->SetBinContent(i+1,data1Pi[i]);
       for (int j=0; j<nErrorPi; j++){
         fa5->SetBinError(j+1,data2Pi[j]);
       }
     }
//__Arrangement_for_Pion-_parameters________________________________________________
     TH1F *fa6 = new TH1F("dEdx_pionminus", "dE/dx parameterization for the pionminus; p*q GeV/c; dE/dx arb.units", 16,-0.625,-0.060);

     for(int i=0; i<nBinsAPi; i++){

     fa6->SetBinContent(i+1,data1APi[i]);
       for (int j=0; j<nErrorAPi; j++){
         fa6->SetBinError(j+1,data2APi[j]);
       }
     }
//__Arrangement_for_Deuterium_parameters____________________________________________
     TH1F *fa7 = new TH1F("dEdx_deuterium", "dE/dx parameterization for the deuterium; p*q GeV/c; dE/dx arb.units", 16, 0.32, 1.16);

     for(int i=0; i<nBinsD; i++){

     fa7->SetBinContent(i+1,data1D[i]);
       for (int j=0; j<nErrorD; j++){
         fa7->SetBinError(j+1,data2D[j]);
       }
     }
//__Arrangement_for_Tritium_parameters______________________________________________
     TH1F *fa8 = new TH1F("dEdx_tritium", "dE/dx parameterization for the tritium; p*q GeV/c; dE/dx arb.units", 18, 0.39, 1.28);

     for(int i=0; i<nBinsT; i++){

     fa8->SetBinContent(i+1,data1T[i]);
       for (int j=0; j<nErrorT; j++){
         fa8->SetBinError(j+1,data2T[j]);
       }
     }
//__Arrangement_for_He3_parameters__________________________________________________
     TH1F *fa9 = new TH1F("dEdx_Helium3", "dE/dx parameterization for the Helium3; p*q GeV/c; dE/dx arb.units", 17, 0.260, 0.750);

     for(int i=0; i<nBinsHe3; i++){

     fa9->SetBinContent(i+1,data1He3[i]);
       for (int j=0; j<nErrorHe3; j++){
         fa9->SetBinError(j+1,data2He3[j]);
       }
     }
//__Arrangement_for_He4_parameters__________________________________________________
     TH1F *fa10 = new TH1F("dEdx_Helium4", "dE/dx parameterization for the Helium4; p*q GeV/c; dE/dx arb.units", 18, 0.350, 0.900);

     for(int i=0; i<nBinsHe4; i++){

     fa10->SetBinContent(i+1,data1He4[i]);
       for (int j=0; j<nErrorHe4; j++){
         fa10->SetBinError(j+1,data2He4[j]);
       }
     }

//__Fit_funtions_for_all_particles__________________________________________________
    TF1 *fitparPr = new TF1("fitparPr",parPr,0.13,1.09,5);
    TF1 *fitparAntiPr = new TF1("fitparAntiPr",parAntiPr,-1.09,-0.148,5);
    TF1 *fitparK = new TF1("fitparK",parK,0.235,0.83,5);
    TF1 *fitparKminus = new TF1("fitparKminus",parKminus,-0.83,-0.255,5);
    TF1 *fitparPi = new TF1("fitparPi",parPi,0.069,0.625,5);
    TF1 *fitparPiminus = new TF1("fitparPiminus",parPiminus,-0.625,-0.060,5);
    TF1 *fitparD = new TF1("fitparD",parD,0.32,1.16,5);
    TF1 *fitparT = new TF1("fitparT",parT,0.39,1.28,5);
    TF1 *fitparHe3 = new TF1("fitparHe3",parHe3,0.260,0.750,5);
    TF1 *fitparHe4 = new TF1("fitparHe4",parHe4,0.350,0.900,5);

//__Outline of the functions________________________________________________________
     h__dedx->Draw("same");
 
     fa1->Draw("same");
     fa2->Draw("same");
     fa3->Draw("same");
     fa4->Draw("same");
     fa5->Draw("same");
     fa6->Draw("same");
     fa7->Draw("same");
     fa8->Draw("same");
     fa9->Draw("same");
     fa10->Draw("same");

     fitparPr->SetLineColor(2);
     fa1->Fit("fitparPr");
     printf("Parameters for the Bethe-Bloch function of the Proton\n");
     fitparAntiPr->SetLineColor(2);
     fa2->Fit("fitparAntiPr");
     printf("Parameters for the Bethe-Bloch function of the Anti-proton\n");
     fitparK->SetLineColor(3);
     fa3->Fit("fitparK");
     printf("Parameters for the Bethe-Bloch function of the Kaon+\n");
     fitparKminus->SetLineColor(3);
     fa4->Fit("fitparKminus");
     printf("Parameters for the Bethe-Bloch function of the Kaon-\n");
     fitparPi->SetLineColor(1);
     fa5->Fit("fitparPi");
     printf("Parameters for the Bethe-Bloch function of the Pion+\n");
     fitparPiminus->SetLineColor(1);
     fa6->Fit("fitparPiminus");
     printf("Parameters for the Bethe-Bloch function of the Pion-\n");
     fitparD->SetLineColor(5);
     fa7->Fit("fitparD");
     printf("Parameters for the Bethe-Bloch function of the Deuterium\n");
     fitparT->SetLineColor(6);
     fa8->Fit("fitparT");
     printf("Parameters for the Bethe-Bloch function of the Tritium\n");
     fitparHe3->SetLineColor(7);
     fa9->Fit("fitparHe3");
     printf("Parameters for the Bethe-Bloch function of the Helium3\n");
     fitparHe4->SetLineColor(8);
     fa10->Fit("fitparHe4");
     printf("Parameters for the Bethe-Bloch function of the Helium4\n");

}

