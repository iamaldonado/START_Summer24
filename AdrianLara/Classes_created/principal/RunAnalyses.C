void RunAnalyses(int nEvents = -1){

   gSystem->Load("libZdc.so") ;
   gSystem->Load("libMpdPhysics.so") ;

   MpdAnalysisManager man("ManagerAnal", nEvents) ;
   man.InputFileList("lista_no.txt") ;
   man.ReadBranches("MCTrack,MCEventHeader,Vertex,MPDEvent") ; 
   man.SetOutput("histos.root") ;
   
//   MpdCentralityAll pCentr("pCentr","pCentr") ;
//   man.AddTask(&pCentr) ;
   
//   MpdEventPlaneAll pEP("pEP","pEP") ;
//   man.AddTask(&pEP) ;

	Fixed_Analysis Task("Analysis","Salida_prueba") ;
   // FixedAnalysis Task("Analysis","V20");
   man.AddTask(&Task) ;
   man.Process() ;

}
