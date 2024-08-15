#ifndef CleanClass_h
#define CleanClass_h

#include "MpdMCTrack.h"
#include "MpdAnalysisTask.h"
#include "MpdVertex.h"

//Include the Histograms 
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


//To change the name of the class, replace all "CleanClass" for your own "Class_Name"
class CleanClass: public MpdAnalysisTask {
  

  public:
    CleanClass() {}
    CleanClass(const char *name, const char *outputName = "taskName");
    ~CleanClass() {} // TODO: normal descructor with cleaning off histos
		//UserInit, ProcessEvent, Finish are default functions on classes
    void UserInit();
    void ProcessEvent(MpdAnalysisEvent &event);
    void Finish();
    


  protected:
		//Branches
		TClonesArray *fVertex = nullptr;
		FairMCEventHeader *fMCEventHeader = nullptr;
		//Vertex Object
		MpdVertex *Vertex = nullptr;

		//Histograms
		TH1F *Hist_Example = nullptr;
  private:

  ClassDef(CleanClass, 1);
};
#endif

