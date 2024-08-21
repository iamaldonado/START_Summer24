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

void FitsigmaPro(){

    //Histogram dE/dx
     TFile *file1 = new TFile("/home/alejandro/Documentos/EnerClass/Fit functions/taskEnerAll.root");
     TH2F *h_dedx_p = (TH2F*)file1->Get("h_dedx_p");

//__Proton+1sima_Bins___________________________________________________________________
    const int nBinsPplus1=18;
    Stat_t data1Pplus1[nBinsPplus1]={30.80537,21.81227,15.77434,11.66682,8.97788,7.177236,5.939882,4.971751,4.295575,3.788829,3.375116,3.051654,2.801907,2.590339,2.420626,2.27716,2.158214,2.058432};

//__Proton-1sima_Bins___________________________________________________________________
    const int nBinsPminus1=18;
    Stat_t data1Pminus1[nBinsPminus1]={19.11903,14.11473,10.56006,8.11146,6.6493,5.462224,4.573218,3.971809,3.469265,3.076031,2.773964,2.527726,2.325173,2.160701,2.019734,1.90412,1.805886,1.721588};


//__Proton+2sima_Bins___________________________________________________________________
    const int nBinsPplus2=18;
    Stat_t data1Pplus2[nBinsPplus2]={36.64854,25.66104,18.38148,13.4445,10.14217,8.034742,6.623214,5.471722,4.70873,4.145228,3.675692,3.313618,3.040274,2.805158,2.621072,2.46368,2.334378,2.226854};

//__Proton-2sima_Bins___________________________________________________________________
    const int nBinsPminus2=18;
    Stat_t data1Pminus2[nBinsPminus2]={13.27586,10.26596,7.95292,6.33378,5.48501,4.604718,3.889886,3.471838,3.05611,2.719632,2.473388,2.265762,2.086806,1.945882,1.819288,1.7176,1.629722,1.553166};


//__Proton+3sima_Bins___________________________________________________________________
    const int nBinsPplus3=18;
    Stat_t data1Pplus3[nBinsPplus3]={42.49171,29.50981,20.98862,15.22218,11.30646,8.892248,7.306546,5.971693,5.121885,4.501627,3.976268,3.575582,3.278641,3.019977,2.821518,2.6502,2.510542,2.395276};

//__Proton-3sima_Bins___________________________________________________________________
    const int nBinsPminus3=18;
    Stat_t data1Pminus3[nBinsPminus3]={7.43269,6.41719,5.34578,4.5561,4.32072,3.747212,3.206554,2.971867,2.642955,2.363233,2.172812,2.003798,1.848439,1.731063,1.618842,1.53108,1.453558,1.384744};

//__errors__________________________________________________________
    const int nError=18;
    Stat_t data2[nError]={0.1170546,0.0797691,0.0998064,0.02884405,0.0128405,0.01112165,0.02668003,0.00357914,0.00378996,0.003144379,0.002145787,0.001520423,0.000728888,0.000801909,0.000834705,0.000589733,0.000485966,0.000455834};

//__Arrangement_for_proton+1sigma_parameters_______________________________________________ 
    TH1F *fa1 = new TH1F("dEdx_proton+1sigma", "dE/dx parameterization for the proton+1sigma; p*q GeV/c; dE/dx arb.units", 18, 0.13, 1.09);

    for(int i=0; i<nBinsPplus1; i++){

    fa1->SetBinContent(i+1,data1Pplus1[i]);
      for (int j=0; j<nError; j++){
        fa1->SetBinError(j+1,data2[j]);
      }
    }

//__Arrangement_for_proton-1sigma_parameters_______________________________________________ 
    TH1F *fa2 = new TH1F("dEdx_proton-1sigma", "dE/dx parameterization for the proton-1sigma; p*q GeV/c; dE/dx arb.units", 18, 0.13, 1.09);

    for(int i=0; i<nBinsPminus1; i++){

    fa2->SetBinContent(i+1,data1Pminus1[i]);
      for (int j=0; j<nError; j++){
        fa2->SetBinError(j+1,data2[j]);
      }
    }


    TF1 *fitparPr = new TF1("fitparPr",parPr,0.13,1.09,5);

     h_dedx_p->Draw("same");

     fa1->Draw("same");
     fa2->Draw("same");

     fitparPr->SetLineColor(2);
     fa1->Fit("fitparPr");
     printf("Parameters for the Bethe-Bloch function of the Proton+1sigma\n");
     fa2->Fit("fitparPr");
     printf("Parameters for the Bethe-Bloch function of the Proton-1sigma\n");

}
