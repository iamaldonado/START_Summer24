#include <iostream>
#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <vector>

#include <Rtypes.h>
#include <TChain.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TRandom3.h>

//#include <MpdEvent.h>
//#include <MpdZdcDigi.h>
//#include <MpdPid.h>
//#include <MpdTrack.h>
//#include <MpdKalmanTrack.h>
//#include <MpdVertex.h>
R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "mpdloadlibs.C"

void MpdDstReader(TChain *inChain, TString outFileName)
{
  TStopwatch timer;
  timer.Start();

  if (!inChain) return;

  const Double_t cut_pt = 0.05; // default: 0.15 GeV/c
  //const Double_t cut_eta = 0.5; // default: 0.5
  const Double_t cut_eta_min = -1.0; // default: 0.5
  const Double_t cut_eta_max = 2.0; // default: 0.5
  const Int_t cut_nhits = 16;   // default: 16
  const Double_t dca_cut = 2.0; // default: 0.5 cm

  TH1F *hRefMult = new TH1F("hRefMultSTAR","hRefMultSTAR",2500,0,2500);
  TH2F *hBvsRefMult = new TH2F("hBvsRefMult","hBvsRefMult",2500,0,2500,200,0.,20.);
  TH1F *hrap = new TH1F("hrap","hrap",200,-10,10);
  TH1F *hpt = new TH1F("hpt","hpt",200,-5,5);
  TH1F *hpttpc = new TH1F("hpttpc","hpt",200,-5,5);

  FairMCEventHeader *MCHeader;
  TClonesArray *MCTracks;
  MpdEvent *MPDEvent;
  TClonesArray *MpdGlobalTracks;
  TClonesArray *TPCEvent;

  MCHeader = nullptr;
  MCTracks = nullptr;
  MPDEvent = nullptr;
  TPCEvent = nullptr;

  inChain->SetBranchAddress("MCEventHeader.", &MCHeader);
  inChain->SetBranchAddress("MCTrack", &MCTracks);
  inChain->SetBranchAddress("MPDEvent.", &MPDEvent);
  inChain->SetBranchAddress("TpcKalmanTrack", &TPCEvent);
  
  std::vector<Long64_t> vEvents;
  std::random_device rd;
  std::mt19937 g(rd());

//  Long64_t Nentries = 5e5; // or (Long64_t) inChain->GetEntries();
  Long64_t Nentries =  inChain->GetEntries();

  // Starting event loop
  Long64_t Nevents;
  Nevents = (Nentries < 5e5) ? Nentries : 5e5;
  Int_t refMult, nhits;
  Double_t pt, eta;
  Long64_t iEvent;
  TRandom3 *rnd = new TRandom3();
  Double_t prob_skip;
  for (Long64_t jentry=0; jentry<Nevents;jentry++)
  {
    iEvent = jentry;

    inChain->GetEntry(iEvent);
    if (jentry%1000 == 0) std::cout << "Event ["
      << jentry << "/" << Nevents << "]" << std::endl;
    refMult = 0;
    MpdGlobalTracks = (TClonesArray*) MPDEvent->GetGlobalTracks();
    Int_t ntracks = MpdGlobalTracks->GetEntriesFast();
    for (int iTr=0; iTr<ntracks; iTr++)
    {
      auto mpdtrack = (MpdTrack*) MpdGlobalTracks->UncheckedAt(iTr);
//#ifdef _OLD_MCSTACK_
//      auto mctrack = (FairMCTrack*) MCTracks->UncheckedAt(mpdtrack->GetID());
//#endif
//#ifdef _NEW_MCSTACK_
      auto mctrack = (MpdMCTrack*) MCTracks->UncheckedAt(mpdtrack->GetID());
//#endif
    
      pt = mpdtrack->GetPt();
      eta = mpdtrack->GetEta();
      nhits = mpdtrack->GetNofHits();
      hrap->Fill(eta);
      hpt->Fill(TMath::Abs(pt));

      if (TMath::Abs(pt) < cut_pt) continue;
      if (eta < cut_eta_min || eta > cut_eta_max) continue;
      //if (TMath::Abs(eta) > cut_eta) continue;
      if (nhits < cut_nhits) continue;

      // Primary track selection
      if (TMath::Sqrt(TMath::Power(mpdtrack->GetDCAX(),2) + TMath::Power(mpdtrack->GetDCAY(),2) + TMath::Power(mpdtrack->GetDCAZ(),2)) > dca_cut) continue;

      refMult++;
    }
    hRefMult->Fill(refMult);
    hBvsRefMult->Fill(refMult,MCHeader->GetB());

 //   MpdGlobalTracks = (TClonesArray*) MPDEvent->GetGlobalTracks();
    Int_t ntrackstpc = TPCEvent->GetEntriesFast();
    for (int iTr=0; iTr<ntrackstpc; iTr++)
    {
      auto mpdtrack = (MpdKalmanTrack*) TPCEvent->UncheckedAt(iTr);
      hpttpc->Fill(mpdtrack->Pt());
}
  }
  TFile *fo = new TFile(outFileName.Data(),"recreate");
  fo->cd();
  hRefMult->Write();
  hBvsRefMult->Write();
  hrap->Write();
  hpt->Write();
  hpttpc->Write();
  fo->Close();

  timer.Stop();
  timer.Print();
}
