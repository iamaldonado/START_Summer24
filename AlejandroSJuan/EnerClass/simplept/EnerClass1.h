#ifndef ENERCLASS1_H
#define ENERCLASS1_H

#include "MpdAnalysisTask.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanFilter.h"
#include "MpdMCTrack.h"
#include "MpdVertex.h"
#include "MpdPid.h"
#include "MpdPidQA.h"
#include "MpdTofMatchingData.h"
#include "MpdTofMatching.h"
#include "MpdHelix.h"
//#include "MpdEventHeader.h"

#include "TH2.h"
#include "TString.h"

#include <string>

class EnerClass1: public MpdAnalysisTask, public MpdPid {
  MpdPid  *fPID;
  
public:
  EnerClass1();
  EnerClass1(const char *name, const char *outputName = "taskName");
  
  ~EnerClass1();  
  
  void UserInit();
  void ProcessEvent(MpdAnalysisEvent &event);
  void Finish();

//__________________________________________________________________________________

  TClonesArray *mMCTracks                       = nullptr;
  TClonesArray *mKalmanTracks                   = nullptr;
  TClonesArray *mMpdGlobalTracks                = nullptr;
  TClonesArray *tofMatches                      = nullptr;

//Histograms for dE/dx of each particle__________________________________________________________________________________

  TH2F *h__dedx;
  TH2F *h__dedxHe3;
  TH2F *h__dedxHe4;
  TH2F *h__dedxt;
  TH2F *h__dedxd;
  TH2F *h__dedxap;
  TH2F *h__dedxp;
  TH2F *h__dedxkp;
  TH2F *h__dedxkm;
  TH2F *h__dedxpip;
  TH2F *h__dedxpim;
  TH2F *h__dedxpiparam;
  
//Histograms for m² of each particle_____________________________________________________________________________________

  TH2F *h__m2;
  TH2F *h__m2He3;
  TH2F *h__m2He4;
  TH2F *h__m2t;
  TH2F *h__m2d;
  TH2F *h__m2ap;
  TH2F *h__m2p;
  TH2F *h__m2kp;
  TH2F *h__m2km;
  TH2F *h__m2pip;
  TH2F *h__m2pim;

//Histograms of transverse momentum for Monte Carlo tracks____________________________

  TH1F *PtProMC1;
  TH1F *PtKaonMC1;
  TH1F *PtPionMC1;
//Histograms of transverse momentum for reconstructed tracks__________________________
  TH1F *PtProreco1;
  TH1F *PtProreco2;
  TH1F *PtKaonreco1;
  TH1F *PtKaonreco2;
  TH1F *PtPionreco1;
  TH1F *PtPionreco2; 
//Histograms of dE/dx but with limits_________________________________________________

  TH2F *h__dedxpcut;
  TH2F *h__dedxkpcut;
  TH2F *h__dedxpipcut;
 
  ClassDef(EnerClass1,1);
};
#endif

