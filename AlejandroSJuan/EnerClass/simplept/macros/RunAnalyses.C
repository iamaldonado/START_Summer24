/*! \file RunAnalyses.C
    \brief Macro to run the MpdRoot analysis "train".

    Usage:
    root -l -b -q RunAnalyses.C
*/

void RunAnalyses(int nEvents = -1){
  gSystem->Load("libZdc.so");
  gSystem->Load("libMpdPhysics.so");


  MpdAnalysisManager man("ManagerAnal", nEvents);
  man.InputFileList("listareq28_22.txt");
  //man.ReadBranches("*");
  man.ReadBranches("MCTrack,TpcKalmanTrack,ZdcDigi,Vertex,MPDEvent,TOFMatching");
  man.SetOutput("histos.root");

  MpdCentralityAll pCentr("pCentr", "pCentr");
  man.AddTask(&pCentr);

 /* MpdPid taskMpdPid("taskMpdPid","taskMpdPid");
  man.AddTask(&taskMpdPid);
*/
/*
  MpdTrackPidMaker pPID("pPID","pPID") ;
  man.AddTask(&pPID) ;
*/

/*  MpdNuclei taskNuclei("taskNuclei", "taskNuclei", "NucleiAna.json");
  man.AddTask(&taskNuclei);
*/
  EnerClass1 taskEnerClass1("taskEner","taskEner");
  //EnerClass1 taskEnerClass1("taskEner","taskEner","NucleiAna.json");
  man.AddTask(&taskEnerClass1);

  man.Process();
}

