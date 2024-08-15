#include "CleanClass.h"

ClassImp(CleanClass); 

CleanClass::CleanClass(const char *name, const char *outputName) : MpdAnalysisTask(name,outputName)
{
   fOutputList = nullptr;
}
void CleanClass::UserInit()
{
  //Define a List, it will have all the Histograms, TProfile, etc
  fOutputList = new TList(); //If you dont put your Histogram on this List it will not be display
  fOutputList->SetOwner(kTRUE);

  TH1::AddDirectory(kFALSE);

	Double_t bins = 200;
	Double_t xlow = 0;
	Double_t xup = 5;

	//Declaration of Histograms should be done on .h file
	Hist_Example = new TH1F("name_hist","title_hist; x-axis; y-axis",bins,xlow,xup);
  fOutputList->Add(Hist_Example);//Adding Histogram to a List
	//Display of Histograms on a file it will be automatically done by MpdAnalysisManager
}

void CleanClass::ProcessEvent(MpdAnalysisEvent &event)
{
	//It will process Event per Event

	//To obtain Branches (wanted branches define on RunAnalysis.C) event.fbranch_name
	fVertex = event.fVertex;
  fMCEventHeader = event.fMCEventHeader;

  // First Vertex
  Vertex = (MpdVertex*) fVertex->First();

  fMCTracks = event.fMCTrack;

  // if(!selectEvent(event)){
  // return;
  // }

	Int_t NumTracksMC;
  Int_t NumTracksVertex;
  Int_t NumTracksReco;

	  //Number Tracks for MC, Reco, Vertex
  NumTracksMC = fMCTracks->GetEntriesFast();//Tracks MonteCarlo (simulated)
	NumTracksReco = event.fMPDEvent->GetGlobalTracks()->GetEntriesFast();// Tracks use for Reconstruction
  NumTracksVertex = Vertex->GetNTracks(); //Tracks used to obtain Vertex

	//MonteCarlo Loop
	for(Int_t i = 0; i < NumTracksMC; i++){
		// Getting Tracks MonteCarlo
		MpdMCTrack *MCTrack = (MpdMCTrack*)fMCTracks->UncheckedAt(i);
		if(MCTrack->GetMotherId() == -1){ //Primary Particles = -1 , Secondary Particles != -1

		}else{

		}
		
	}

	//Reco Loop
	for(Int_t j = 0; j < NumTracksReco; j++){
		//Getting Reco Tracks
		MpdTrack *RecoTrack = (MpdTrack*) event.fMPDEvent->GetGlobalTracks()->UncheckedAt(j);
		// ------ Montecarlo Identification ----

		Int_t ID = RecoTrack->GetID();
		MpdMCTrack *MCTrack = (MpdMCTrack*)fMCTracks->UncheckedAt(ID);
		if(!MCTrack) continue;

		// -------------------------------------
		Double_t pT;
		
		pT = RecoTrack->GetPt();

		Hist_Example->Fill(pT);
	}
  
}
void CleanClass::Finish()
{
}