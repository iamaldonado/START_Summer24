void RunAnalyses(int nEvents = -1){

   gSystem->Load("libZdc.so") ;
   gSystem->Load("libMpdPhysics.so") ;

   MpdAnalysisManager man("ManagerAnal", nEvents) ;
   man.InputFileList("lista.txt") ;
   man.ReadBranches("MCTrack,MCEventHeader,Vertex,MPDEvent,TpcKalmanTrack,ZdcDigi,TOFMatching") ; 
   man.SetOutput("histos.root") ;
   
//   MpdCentralityAll pCentr("pCentr","pCentr") ;
//   man.AddTask(&pCentr) ;
   
//   MpdEventPlaneAll pEP("pEP","pEP") ;
//   man.AddTask(&pEP) ;

	MpdPtMCAnalysisTask Task("Analysis","salida1") ;
   man.AddTask(&Task) ;
   man.Process() ;

}
