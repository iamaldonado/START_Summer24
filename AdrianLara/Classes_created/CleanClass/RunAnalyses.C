void RunAnalyses(int nEvents = -1){//nEvents = -1 will process all events

   gSystem->Load("libZdc.so");
   gSystem->Load("libMpdPhysics.so");

   MpdAnalysisManager train("ManagerAnal", nEvents); //Declare a "Train"
   train.InputFileList("list.txt");//.txt with all directories of .root files with the events on them
   train.ReadBranches("MCTrack,MCEventHeader,Vertex,MPDEvent"); //Input all the Branches that will be used, "*" for all Branches
   train.SetOutput("histos.root");
   
   // MpdCentralityAll pCentr("pCentr", "pCentr");
   // train.AddTask(&pCentr);
   
   //"Wagon" Declaration
   CleanClass wagon("Analysis","output");//"wagon_name", "output_file_name"

   train.AddTask(&wagon);//Add wagon to train
   train.Process();//start wagon

}