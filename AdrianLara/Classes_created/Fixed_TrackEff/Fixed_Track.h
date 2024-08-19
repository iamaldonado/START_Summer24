#ifndef Fixed_Track_h
#define Fixed_Track_h

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



class Fixed_Track: public MpdAnalysisTask {
  

  public:
    Fixed_Track() {}
    Fixed_Track(const char *name, const char *outputName = "taskName");
    ~Fixed_Track() {} // TODO: normal descructor with cleaning off histos
    void init_Profile2D();
    void init_Hist1D();
    void init_Hist2D();
    void init_HistPart();
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


  //MC Reco and Pri, Sec Particles, Kaons, Protons, Pions
  const static int NumHistMCPSKpPi1 = 1;
  char HistMCPSKpPi1_titles[NumHistMCPSKpPi1][300]={"pT"                                             
                                                    };

  axis1D HistMCPSKpPi1_axis[NumHistMCPSKpPi1] = {{"pT (Gev/c)"}
                                                };
  interval1D HistMCPSKpPi1_inter[NumHistMCPSKpPi1] = {{0,5}
                                                     };

  bin1D HistMCPSKpPi1_bins[NumHistMCPSKpPi1] = {{200}
                                               };
  TH1F *HistMCPSKpPi1[2][2][4][NumHistMCPSKpPi1];

    TH2F *DZVertex_zVertex = nullptr;
    TH2F *DZVertex_Multi = nullptr;
    TH2F *Dpt_pt_Multi = nullptr;
    TH2F *zVertex_b = nullptr;
    
    
    TProfile *Dpt_pt_Multi_Profile = nullptr;
    TProfile *Dpt_pt_Multi_Profile_cuts = nullptr;
    
    TProfile *DZVertex_Multi_Profile = nullptr;


    //  Branches
    TClonesArray *fVertex = nullptr;
    FairMCEventHeader *fMCEventHeader = nullptr;

    MpdVertex *Vertex = nullptr;

    

  private:

  ClassDef(Fixed_Track, 1);
};
#endif

