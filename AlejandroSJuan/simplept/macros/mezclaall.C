void mezclaall(){

TFileMerger m;

for(Int_t jk=0; jk<20; jk++){

  char NameofFile[50];
  sprintf(NameofFile,"dir%d/taskNuclei.root",jk); //cambia polana de acuerdo al nombre y directorio de tu archivo

if(!(m.AddFile(NameofFile)))continue;

m.OutputFile("taskNuclei.root"); // nombre de tu archivo
m.Merge();

exit(0);
}

}
