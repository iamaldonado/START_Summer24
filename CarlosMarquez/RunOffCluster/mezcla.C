void mezcla(){

char NameofFile[90];

TFileMerger m;

for(Int_t i=0; i<20; ++i){
for(Int_t j=0; j<2; ++j){

sprintf(NameofFile,"dir%d/xa%d/tasklowMgF.root",i,j);
if(!(m.AddFile(NameofFile)))continue;
cout<< NameofFile << endl;
}
}

m.OutputFile("tasklowMgFAll.root");
m.Merge();
exit(0);

}
