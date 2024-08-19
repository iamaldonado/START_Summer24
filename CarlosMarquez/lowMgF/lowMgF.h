#ifndef LOWMGF_H
#define LOWMGF_H

#include "MpdAnalysisTask.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanFilter.h"
#include "MpdMCTrack.h"
#include "MpdVertex.h"
#include "MpdPid.h"
#include "MpdTofMatchingData.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TString.h"
#include "FairMCEventHeader.h"

#include <string>

class lowMgF: public MpdAnalysisTask, public MpdPid {
  MpdPid  *fPID;
  
public:
  lowMgF();
  lowMgF(const char *name, const char *outputName = "taskName");
  
  ~lowMgF();  
  
  void UserInit();
  void ProcessEvent(MpdAnalysisEvent &event);
  void Finish();
  
//__________________________________________________________________________________

  TClonesArray *mMCTracks                       = nullptr;
  TClonesArray *mKalmanTracks                   = nullptr;
  TClonesArray *mMpdGlobalTracks                = nullptr;
  TClonesArray *tofMatches                      = nullptr;
  TClonesArray *vtxs 			        = nullptr;
  FairMCEventHeader *mMCEventHeader 		= nullptr;

//__________________________________________________________________________________

 // Histograms Transverse Momemtum vs Rapidity

	// Monte Carlo
   TH2F *PtMCvsEta;

	// Reco
   TH2F *PtRECOvsEta;
   TH1F *VertexPosition;

	// Vtx Resolution
   // TH2F
   TH2F *DZvsZReco;
   TH2F *DZNtracks;
   TH2F *DZb;

   // TProfile
   TProfile *TPDZvsZReco;
   TProfile *TPDZNtracks;
   TProfile *TPDZb;

   TProfile2D *TPDZNtracksW;
   TProfile2D *TPDZbW;

	// without Cuts
   TProfile *PtNHits;		// Number of Hits
   TProfile *PtDCAGlobal;	// DCA Global
   TProfile *PtDCAGlobalP;	// DCA Global Primary
   TProfile *PtDCAGlobalS;	// DCA Global Secondary
   TProfile *PtEta;		// PseudoRapidity

   TProfile2D *PtEtaDPt;     	// Resolution

	// With Cuts in NHits
   TProfile *PtNHitsC;		// Number of Hits
   TProfile *PtDCAGlobalC;	// DCA Global
   TProfile *PtDCAGlobalPC;	// DCA Global Primary
   TProfile *PtDCAGlobalSC;	// DCA Global Secondary
   TProfile *PtEtaC;		// PseudoRapidity 

   TProfile2D *PtEtaDPtC; 	// Resolution

	// With Cuts in Eta
   TProfile *PtNHitsCE;		// Number of Hits
   TProfile *PtDCAGlobalCE;	// DCA Global
   TProfile *PtDCAGlobalPCE;	// DCA Global Primary
   TProfile *PtDCAGlobalSCE;	// DCA Global Secondary
   TProfile *PtEtaCE;		// PseudoRapidity 

   TProfile2D *PtEtaDPtCE; 	// Resolution


	// DCA Global with cuts in Number of Hits
   TH1F *DCAGPC; 		// DCA Global Primary
   TH1F *DCAGSC; 		// DCA Global Secondary
   TH1F *DCAGC;			// DCA Global

	// Track Efficiency
   TH1F *PtRecoPionP;
   TH1F *PtRecoPionS;
   TH1F *PtMCPionP;
   TH1F *PtMCPionS;

   TH1F *PtRecoProtonP;
   TH1F *PtRecoProtonS;
   TH1F *PtMCProtonP;
   TH1F *PtMCProtonS;

   TH1F *PtRecoKaonP;
   TH1F *PtRecoKaonS;
   TH1F *PtMCKaonP;
   TH1F *PtMCKaonS;

   TProfile *DPtPtReco;

	// Multiplicity
   TProfile *VtxMult;

  ClassDef(lowMgF,1);
};
#endif

