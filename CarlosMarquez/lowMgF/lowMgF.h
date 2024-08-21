#ifndef LOWMGF_H
#define LOWMGF_H

#include "MpdAnalysisTask.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanFilter.h"
#include "MpdMCTrack.h"
#include "MpdVertex.h"
#include "MpdTofMatchingData.h"

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TString.h"
#include "FairMCEventHeader.h"

#include <string>

class lowMgF: public MpdAnalysisTask {
  
public:
  lowMgF();
  lowMgF(const char *name, const char *outputName = "taskName");
  
  ~lowMgF();  
  
  void UserInit();
  void ProcessEvent(MpdAnalysisEvent &event);
  void Finish();
  
//__________________________________________________________________________________

  TClonesArray *mMCTracks                       = nullptr;	// Monte-Carlo tracks
  TClonesArray *mKalmanTracks                   = nullptr;	// TPC Kalman tracks
  TClonesArray *mMpdGlobalTracks                = nullptr;	// MPD Global tracks
  TClonesArray *tofMatches                      = nullptr;	// ToF Matching
  TClonesArray *vtxs 			        = nullptr;	// Vertex Position
  FairMCEventHeader *mMCEventHeader 		= nullptr;	// MC Event Header track

//__________________________________________________________________________________

//____All Histograms________________________________________________________________

//------Primary Vertex Position-------
   
   TH1F *VertexPosition;	// Vtx Resolution
   
   TH2F *DZvsZReco;		// Resolution vs Z Recontruction
   TH2F *DZNtracks;		// Resolution vs Number of Track
   TH2F *DZb;      		// Resolution vs Impact Parameter

   TProfile *TPDZvsZReco;	// TProfile Resolution vs Z Recontruction
   TProfile *TPDZNtracks;	// TProfile Resolution vs Number of Track
   TProfile *TPDZb;		// TProfile Resolution vs Impact Parameter

   TProfile2D *TPDZNtracksW;	// TProfile2D/Weight Resolution vs Number of Track
   TProfile2D *TPDZbW;		// TProfile2D/Weight Resolution vs Impact Parameter

   TProfile *VtxMult;		// Transvese Momentum vs Multiplicity

//-------Tranverse Momemtum----------

   TH2F *PtMCvsEta;		// Monte-Carlos vs Psuedo-Rapidity

   TH2F *PtRECOvsEta;		// Recontructed vs Pseudo-Rapidity

   TProfile *DPtPtReco;		// Transverse Momentum as a function of Transverse Momentum

	// without Cuts
   TProfile *PtNHits;		// Number of Hits
   TProfile *PtEta;		// Pseudo-Rapidity
   TProfile *PtDCAGlobal;	// DCA Global
   TProfile *PtDCAGlobalP;	// DCA Global Primary
   TProfile *PtDCAGlobalS;	// DCA Global Secondary

   TProfile2D *PtEtaDPt;     	// Resolution of the transverse momentum

	// With Cuts in NHits
   TProfile *PtNHitsC;		// Number of Hits
   TProfile *PtEtaC;		// PseudoRapidity 
   TProfile *PtDCAGlobalC;	// DCA Global
   TProfile *PtDCAGlobalPC;	// DCA Global Primary
   TProfile *PtDCAGlobalSC;	// DCA Global Secondary

   TProfile2D *PtEtaDPtC; 	// Resolution

	// With Cuts in Eta
   TProfile *PtNHitsCE;		// Number of Hits
   TProfile *PtEtaCE;		// PseudoRapidity 
   TProfile *PtDCAGlobalCE;	// DCA Global
   TProfile *PtDCAGlobalPCE;	// DCA Global Primary
   TProfile *PtDCAGlobalSCE;	// DCA Global Secondary

   TProfile2D *PtEtaDPtCE; 	// Resolution

	// With Cuts in DCA Global
   TProfile *PtNHitsDCA;	// Number of Hits
   TProfile *PtEtaDCA;		// PseudoRapidity 
   TProfile *PtDCAGlobalDCA;	// DCA Global
   TProfile *PtDCAGlobalPDCA;	// DCA Global Primary
   TProfile *PtDCAGlobalSDCA;	// DCA Global Secondary

   TProfile2D *PtEtaDPtDCA; 	// Resolution

	// With Cuts in Transverse Momentum
   TProfile *PtNHitsCPT;	// Number of Hits
   TProfile *PtEtaCPT;		// PseudoRapidity 
   TProfile *PtDCAGlobalCPT;	// DCA Global
   TProfile *PtDCAGlobalPCPT;	// DCA Global Primary
   TProfile *PtDCAGlobalSCPT;	// DCA Global Secondary

   TProfile2D *PtEtaDPtCPT; 	// Resolution

//-------DCA Global with cuts in Number of Hits
   TH1F *DCAGC;			// DCA Global
   TH1F *DCAGPC; 		// DCA Global Primary
   TH1F *DCAGSC; 		// DCA Global Secondary

//-------Track Efficiency-----------------------

   TH1F *PtRecoPionP;		// Reconstruted Primary Pions
   TH1F *PtRecoPionS;		// Reconstruted Secondary Pions
   TH1F *PtMCPionP;		// Monte-Carlo Primary Pions
   TH1F *PtMCPionS;		// Monte-Carlo Secondary Pions

   TH1F *PtRecoProtonP;		// Reconstruted Primary Protons
   TH1F *PtRecoProtonS;		// Reconstruted Secondary Protons
   TH1F *PtMCProtonP;		// Monte-Carlo Primary Protons
   TH1F *PtMCProtonS;		// Monte-Carlo Secondary Protons

   TH1F *PtRecoKaonP;		// Reconstruted Primary Kaons
   TH1F *PtRecoKaonS;		// Reconstruted Secondary Kaons
   TH1F *PtMCKaonP;		// Monte-Carlo Primary Kaons
   TH1F *PtMCKaonS;		// Monte-Carlo Secondary Kaons

  ClassDef(lowMgF,1);
};
#endif
