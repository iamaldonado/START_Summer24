#ifndef Fixed_Analysis_h
#define Fixed_Analysis_h

#include "MpdMCTrack.h"
#include "MpdAnalysisTask.h"
#include "MpdVertex.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"


#include <iostream>
#include <Rtypes.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TStopwatch.h>
using std::cout;
using std::endl;



class Fixed_Analysis: public MpdAnalysisTask {
  

  public:
    Fixed_Analysis() {}
    Fixed_Analysis(const char *name, const char *outputName = "taskName");
    ~Fixed_Analysis() {} // TODO: normal descructor with cleaning off histos
    void init_Profile2D();
    void init_Hist1D();
    void init_Hist2D();
    void UserInit();
    void ProcessEvent(MpdAnalysisEvent &event);
    void Finish();
    

    TClonesArray *fMCTracks = nullptr;


  protected:
    struct axis1D
    {
      char X[70];
    };

    struct axis2D
    {
      char X[70];
      char Y[70];
    };

    struct axis3D
    {
      char X[70];
      char Y[70];
      char Z[70];
    };

    struct interval1D
    {
      Double_t Xlow;
      Double_t Xup;
    };

    struct interval2D
    {
      Double_t Xlow;
      Double_t Xup;
      Double_t Ylow;
      Double_t Yup;
    };
    struct interval3D
    {
      Double_t Xlow;
      Double_t Xup;
      Double_t Ylow;
      Double_t Yup;
      Double_t Zlow;
      Double_t Zup;
    };
    struct bin2D
    {
      Int_t X;
      Int_t Y;
    };

    struct bin1D
    {
      Int_t X;
    };

// ---- Arrays of Histograms---------
//How to use it
  //Lets use a TH1F Primary, Secondary, All particles as an example
  //For each title on HistPST1_titles it will create a Histogram of Primary, Secondary and All.
  //The number of titles should be indicated on NumHistPST1.
  //Finally for each title you should put the name, the interval and the number of bins of the x axis


  //It creates a TH1F for Primary , Secondary or All particles
    const static int NumHistPST1 = 3;
    char HistPST1_titles[NumHistPST1][300]={"DCA"
                                            ,"DCA Nhits"
                                            ,"DCA #eta"                                              
                                            };

    axis1D HistPST1_axis[NumHistPST1] = {{"DCA (cm)"}
                                        ,{"DCA (cm)"}
                                        ,{"DCA (cm)"}
                                        };
    interval1D HistPST1_inter[NumHistPST1] = {{0,5} 
                                              ,{0,5}
                                              ,{0,5}
                                              };

    bin1D HistPST1_bins[NumHistPST1] = {{100}
                                        ,{100} 
                                        ,{100}
                                        };

    TH1F *HistPST1[3][NumHistPST1];

  //It creates a TProfile for Primary , Secondary or All particles
      
    const static int NumProf1 = 9;
    char Prof1_titles[NumProf1][70]={"#eta vs #Delta p_{T}"
                                    ,"DCA vs #Delta p_{T}"
                                    ,"NHits vs #Delta p_{T}"
                                    ,"#eta vs #Delta p_{T} cut z-vertex"
                                    ,"DCA vs #Delta p_{T} cut z-vertex"
                                    ,"NHits vs #Delta p_{T} cut z-vertex"
                                    ,"DCA vs #Delta p_{T} cut eta"
                                    ,"NHits vs #Delta p_{T} cut eta"
                                    ,"DCA vs #Delta p_{T} cut NHits"
                                    };

    axis2D Prof1_axis[NumProf1] = {{"#eta","#Delta p_{T}"}
                                  ,{"DCA","#Delta p_{T}"}
                                  ,{"NHits","#Delta p_{T}"}
                                  ,{"#eta","#Delta p_{T}"}
                                  ,{"DCA","#Delta p_{T}"}
                                  ,{"NHits","#Delta p_{T}"}
                                  ,{"DCA","#Delta p_{T}"}
                                  ,{"NHits","#Delta p_{T}"}
                                  ,{"DCA","#Delta p_{T}"}
                                       };
    interval2D Prof1_inter[NumProf1] = {{-3,3,0,3}
                                       ,{0,5,0,3}
                                       ,{0,50,0,3}
                                       ,{-3,3,0,3}
                                       ,{0,5,0,3}
                                       ,{0,50,0,3}
                                       ,{0,5,0,3}
                                       ,{0,50,0,3}
                                       ,{0,5,0,3}
                                            };

    bin1D Prof1_bins[NumProf1] = {{200}
                                 ,{200}
                                 ,{200}
                                 ,{200}
                                 ,{200}
                                 ,{200}
                                 ,{200}
                                 ,{200}
                                 ,{200}
                                      };

    TProfile *Prof1[3][NumProf1];

    //It creates a TH2F for each MC Primary, Secondary and All particles, Reco Primary, Secondary and All particles
    const static int NumHistPST2 = 1;
    char HistPST2_titles[NumHistPST2][300]={"#eta vs p_{T}"
                                           };

    axis2D HistPST2_axis[NumHistPST2] = {{"#eta", "p_{T} (Gev/c)"}
                                        };


    interval2D HistPST2_inter[NumHistPST2] = {{-2,2,0,3}
                                             };

    bin2D HistPST2_bins[NumHistPST2] = {{100,100}
                                       };

    TH2F *HistPriSecTot2[2][3][NumHistPST2];

    //It creates a TProfile3D for Primary , Secondary or All particles
    const static int NumProf3 = 4;
    char Prof3D_titles[NumProf3][70]={"#eta vs pT"
                                     ,"#eta vs pT cut Nhits"
                                     ,"#eta vs pT cut #eta"
                                     ,"#eta vs pT all cuts"
                                     };

    axis3D Prof3D_axis[NumProf3] = {{"#eta","pT"," #Delta pT"}
                                   ,{"#eta","pT"," #Delta pT"}
                                   ,{"#eta","pT"," #Delta pT"}
                                   ,{"#eta","pT"," #Delta pT"}
                                   };
    interval3D Prof3D_inter[NumProf3] = {{-4,4,0,5,0,4}
                                        ,{-4,4,0,5,0,4}
                                        ,{-4,4,0,5,0,4}
                                        ,{-4,4,0,5,0,4}
                                        };

    bin2D Prof3D_bins[NumProf3] = {{100,100}
                                  ,{100,100}
                                  ,{100,100}
                                  ,{100,100}
                                  };

    TProfile2D *Prof3D[3][NumProf3];

    TH1F *Hist_MC_ZVertex = nullptr;
    TH1F *Hist_R_ZVertex = nullptr;
    TH1F *Hist_Chi = nullptr;
    TH1F *Hist_Multipl_MC = nullptr;
    TH1F *Hist_Multipl_R = nullptr;
    TH1F *Hist_b = nullptr;
    TH1F *MC_pT = nullptr;
    TH1F *Reco_pT = nullptr;
    TH1F *hRefMult = nullptr;

    TH2F *hBvsRefMult = nullptr;
    TH2F *Hist_MC_NTrack_ZVertex = nullptr;
    TH2F *Hist_R_NTrack_ZVertex = nullptr;
    TH2F *Hist_MC_Chi_zVertex = nullptr;
    TH2F *Hist_R_Chi_zVertex = nullptr;
    TH2F *Hist_b_Vertex = nullptr;
    TH2F *Hist_b_Multiplicity_R = nullptr;
    TH2F *Hist_Diff_ZVertex_b = nullptr;
    TH2F *Hist_Diff_ZVertex_NTracks = nullptr;
    TH2F *ZVertex_NTracksMC = nullptr;

    TProfile *Profile_Diff_Vertex_b = nullptr;
    TProfile *Profile_Diff_Vertex_NTracks = nullptr;

    TProfile2D *Prof_2D_eta_pt_DPT_Reco = nullptr;
    TProfile2D *Prof_2D_eta_pt_DPT_MC = nullptr;


    //  Branches
    TClonesArray *fVertex = nullptr;
    FairMCEventHeader *fMCEventHeader = nullptr;

    MpdVertex *Vertex = nullptr;

    

  private:

  ClassDef(Fixed_Analysis, 1);
};
#endif

