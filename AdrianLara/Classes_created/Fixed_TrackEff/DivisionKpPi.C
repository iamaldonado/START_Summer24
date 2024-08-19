struct axis1D
{
  char X[70];
};
struct interval1D
{
  Double_t Xlow;
  Double_t Xup;
};
struct bin1D
{
  Int_t X;
};


void DivisionKpPi(){
  const static int NumHistMCPSKpPi1 = 1;
  char HistMCPSKpPi1_titles[NumHistMCPSKpPi1][300]={"pT"                                             
                                                    };

  char Division_titles[NumHistMCPSKpPi1][300]={"pT_Division"                                             
                                              };

  axis1D HistMCPSKpPi1_axis[NumHistMCPSKpPi1] = {{"pT (Gev/c)"}
                                                };
  interval1D HistMCPSKpPi1_inter[NumHistMCPSKpPi1] = {{0,5}
                                            };

  bin1D HistMCPSKpPi1_bins[NumHistMCPSKpPi1] = {{200}
                                      };
  TH1F *HistMCPSKpPi1[2][2][4][NumHistMCPSKpPi1];

  TH1F *Division[1][2][4][NumHistMCPSKpPi1];

  TFile *file= TFile::Open("Salida_Track_Eff_18.root");


  char name1[20];
  char name2[20];
  char name3[20];
  char hist_name[300];
  char title_name[500];

  for(Int_t i = 0 ; i < 2 ; i++)
  {
    sprintf(name1,"");
    Int_t u = (i%2 == 0) ? sprintf(name1,"MC_") : sprintf(name1, "Reco_");
    for (Int_t j = 0; j < 2; j++)
    {
      sprintf(name2,"");
      Int_t y = (j%2 == 0) ? sprintf(name2,"Primary_") : sprintf(name2,"Secondary_");
      for (Int_t l = 0; l < 4; l++)
      {
        sprintf(name3,"");
        switch (l)
        {
          case 0:
            sprintf(name3,"Kaons_");
            break;
          case 1:
            sprintf(name3,"Protons_");
            break;
          case 2:
            sprintf(name3,"Pions_");
            break;
          case 3:
            sprintf(name3,"All_");
            break;
          default:
            sprintf(name3,"Error");
        }
        for (Int_t k = 0; k < NumHistMCPSKpPi1; k++)
        {
          sprintf(hist_name,"");
          sprintf(title_name,"");

          strcat(hist_name,name1);
          strcat(hist_name,name2);
          strcat(hist_name,name3);
          strcat(hist_name,HistMCPSKpPi1_titles[k]);


          //HistMCPSKpPi1[i][j][l][k] = new TH1F(hist_name , title_name , HistMCPSKpPi1_bins[k].X , HistMCPSKpPi1_inter[k].Xlow , HistMCPSKpPi1_inter[k].Xup );
          HistMCPSKpPi1[i][j][l][k] = (TH1F*)file->Get(hist_name);
          HistMCPSKpPi1[i][j][l][k]->Sumw2();
        }
      } 
    } 
  }
  
  
  TFile *out = TFile::Open("PTDIV.root","RECREATE");
  
  TCanvas *c[8];
  THStack *hstack[4];
  
  for(Int_t i = 0; i < 8 ; i++)
  {
    c[i] = new TCanvas();
    if(i < 4){
      c[i]->Divide(2);
    }
    
  }

  for (Int_t j = 0; j < 2; j++)
    {
      sprintf(name2,"");
      Int_t y = (j%2 == 0) ? sprintf(name2,"Primary_") : sprintf(name2,"Secondary_");
      for (Int_t l = 0; l < 4; l++)
      {
        sprintf(name3,"");
        switch (l)
        {
          case 0:
            sprintf(name3,"Kaons_");
            break;
          case 1:
            sprintf(name3,"Protons_");
            break;
          case 2:
            sprintf(name3,"Pions_");
            break;
          case 3:
            sprintf(name3,"All_");
            break;
          default:
            sprintf(name3,"Error");
        }
        for (Int_t k = 0; k < NumHistMCPSKpPi1; k++)
        {
          sprintf(hist_name,"");
          sprintf(title_name,"");

          //strcat(hist_name,name1);
          strcat(hist_name,name2);
          strcat(hist_name,name3);
          strcat(hist_name,Division_titles[k]);
          strcat(title_name,hist_name);
          strcat(title_name,";");
          strcat(title_name,HistMCPSKpPi1_axis[k].X);
          Division[0][j][l][k] = new TH1F(hist_name , title_name , HistMCPSKpPi1_bins[k].X , HistMCPSKpPi1_inter[k].Xlow , HistMCPSKpPi1_inter[k].Xup );
          
          Division[0][j][l][k]->Divide(HistMCPSKpPi1[1][j][l][k],HistMCPSKpPi1[0][j][l][k]);
          Division[0][j][l][k]->SetOption("E1");
          
          if(j%2 == 1){
            Division[0][j][l][k]->SetLineColor(kRed);
          }
          c[l]->cd(j+1);
          Division[0][j][l][k]->Draw();
          c[l+4]->cd();
          if(j%2 == 0){
            Division[0][j][l][k]->Draw();
          }else{
            Division[0][j][l][k]->Draw("same");
          }
          
          //Division[0][j][l][k]->Draw("same");

          Division[0][j][l][k]->Write();
          //HistMCPSKpPi1[i][j][l][k] = (TH1F*)file->Get(hist_name);
        }
      } 
    }
    
    for(Int_t i = 0; i < 8 ; i++)
  {
    c[i]->Update();
    c[i]->Write();
  }
     
  return;
}
