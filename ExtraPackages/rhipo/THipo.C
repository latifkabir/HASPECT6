#include "THipo.h"
#include <TSystem.h>
#include <TApplication.h>
#include <TObjString.h>
#include <TObjArray.h>


THipo::THipo(){
  fParser=new THipoBankParser("");
  
  //Decide whether input files or dirs from command line
  Int_t argc=gApplication->Argc();
  char** argv=gApplication->Argv();
  for(Int_t i=1;i<argc;i++){
    TString dirname=argv[i];
    if(dirname.Contains("--in=")){
      fUseCLDirs=kTRUE;
      break;
    }
    else
      fUseCLDirs=kFALSE;
  }
  //Look for output dir
  for(Int_t i=1;i<argc;i++){
    TString dirname=argv[i];
    if(dirname.Contains("--out=")){
      dirname.Remove(0,6);
      fOutDirName=dirname;
      cout<<"THipo will write root files to "<<fOutDirName<<endl;
      if(gSystem->MakeDirectory(fOutDirName)==-1){
	cout<<"Warning directory "<<fOutDirName<<" already exists or its parent directory does not exist"<<endl;
      }
    }
  }
  
}
THipo::THipo(TString filename){
  fParser=new THipoBankParser(filename);
}
THipo::~THipo(){
  if(fParser) delete fParser;
  for(UInt_t ib=0;ib<fBanks.size();ib++)
    delete fBanks.at(ib);
  fBanks.clear();
  
}
void THipo::InitFile(TString filename){
  fReader.open(filename.Data());
  fReader.showInfo();
  fReader.readRecordIndex();
  fNRecords = fReader.getRecordCount();
  fRecordCount=0;
}
void THipo::AddCommandLineFiles(){
  //If using command line dirs do that then return
  if(fUseCLDirs) //This is set in the constructor
    return AddCommandLineDirs();
  //look for files from command line
  cout<<"THipo::AddCommandLineFiles() "<<" using input files from commandline"<<endl;
  Int_t argc=gApplication->Argc();
  char** argv=gApplication->Argv();
  for(Int_t i=0;i<argc;i++){
    if(TString(argv[i]).Contains(".hipo")) fInFiles.push_back(argv[i]);
  } 
}
void THipo::AddCommandLineDirs(){
  cout<<"THipo::AddCommandLineDirs() "<<" using input directory from commandline"<<endl;
  //look for files from command line
  Int_t argc=gApplication->Argc();
  char** argv=gApplication->Argv();
  for(Int_t i=1;i<argc;i++){
    TString dirname=argv[i];
    if(dirname.Contains("--in=")){
      dirname.Remove(0,5);
      cout<<"THipo::AddCommandLineDirs() checking for files in "<<dirname<<endl;
      void *dir=gSystem->OpenDirectory(dirname);
      cout<<dir<<endl;
      TString fileName;
      while( (fileName=(gSystem->GetDirEntry(dir)))){
	cout<<fileName<<endl;
 	if(fileName.Contains(TString(".hipo")))
	  fInFiles.push_back(argv[i]);
 	if(fileName==TString("")) break;
     }
    }
    
  }
  
  cout<<"THipo::AddCommandLineDirs() found "<<fInFiles.size()<<" hipo files in --in directories"<<endl;
  
 
}

void THipo::ToRoot(){
  //Use file list
  for(UInt_t i=0;i<fInFiles.size();i++)
    ToRoot(fInFiles.at(i));

}
void THipo::ToRoot(TString filename){
  cout<<"THipo::ToRoot "<<filename<<endl;
  //Convert file called filename
  InitFile(filename);
  cout<<"THipo::ToRoot Number of records "<<fNRecords<<endl;
  //if perviously defined by user configuration or previous hipo file
  if(fBanks.size()==0)ConfigAllBanks();
  filename.ReplaceAll(".hipo",".root");

  //Append the output directory name
  if(fOutDirName.Sizeof()>1){
    TString newfilename=fOutDirName+"/"+gSystem->BaseName(filename);
    filename=newfilename;
  }
  InitOutput(filename);
  //loop over events
  while(NextEvent());
  CloseOutput();
  
}
Bool_t  THipo::NextFile(){
  if(fFilei>=fInFiles.size()) return kFALSE;
  cout<<" THipo::NextFile() Analyse file "<<fInFiles.at(fFilei)<<endl;
  InitFile(fInFiles.at(fFilei++));
  return kTRUE;
}
Long64_t THipo::ProcessFile(TString filename){
  InitFile(filename);
  //cout<<"Number of records "<<fNRecords<<endl;

  //loop over records
  while(NextEvent());
  if(fOutTree) CloseOutput();
  return fNEvents;
}
Bool_t THipo::NextRecord(){//Move to next record in file
  if(fRecordCount>=fNRecords) return kFALSE;
  fReader.readRecord(fRecord,fRecordCount++);
  fNRecEvents = fRecord.getEventCount(); //get events in this record
  fRecEventCount=0;//restart record event counter
  // cout<<"Number of rec events "<<fNRecEvents<<" in record "<<fRecordCount<<endl;
  return kTRUE;
}
Bool_t THipo::NextEvent(){//Read next event in file
 if(fRecEventCount>=fNRecEvents){//Look for next record
    if(!NextRecord()) return kFALSE;//No more records
  }
  fHipoEvent = fRecord.getHipoEvent(fRecEventCount++);
  //Loop over group
  for(UInt_t ig=0;ig<fBanks.size();ig++){
    fCurBank=fBanks.at(ig);
    fCurBank->ClearEvent();
    //Loop over int items
    for(UInt_t in=0;in<fCurBank->NInt();in++){
      ReadItemI(in);
    }
    //Loop over float items
    for(UInt_t fn=0;fn<fCurBank->NFloat();fn++)
      ReadItemF(fn);
    
  }//end group loop
  if(fOutTree)fOutTree->Fill();
  fNEvents++;
  return kTRUE;
}
Bool_t THipo::ConfigBank(TString bankname){
  fBanks.push_back(new THipoBank(bankname));
  Bool_t isOK= fParser->ConfigBank(fBanks.back());
  fBankMap[fBanks[fBanks.size()-1]->GetName()]=fBanks.size()-1;//map name after configure in case items removed
  return isOK;
}
void THipo::ConfigAllBanks(){
  //Look for "bank": lines in .json file and configure
  TString bankline;
  while((bankline=fParser->FindNextLineWith("\"bank\""))!=TString("")){
    bankline.Remove(0,bankline.First(":")+1);
    bankline.Remove(0,bankline.First("\"")+1);
    bankline.Remove(bankline.First("\""),bankline.Length()-bankline.First("\"")); 
    ConfigBank(bankline.Data());
  }
  
}
void THipo::InitOutput(TString filename){
  fOutFile=new TFile(filename,"recreate");
  fOutTree=new TTree("hipo2root","CLAS12 banks in ROOT");
  //Make branches from all the banks
  for(UInt_t ib=0;ib<fBanks.size();ib++)
    fBanks.at(ib)->InitTree(fOutTree);
  // fOutTree->Print("");
}
void THipo::CloseOutput(){
  fOutTree->Write();
  fOutFile->Close();
  delete fOutFile;
}
THipoBank* THipo::GetBank(TString name){
  auto it=fBankMap.find(name);
  if(it==fBankMap.end()){
    cout<<"Error THipo* GetBank( No bank comfigured with name = "<<name<<endl;
    return nullptr;
  }
  return fBanks[fBankMap[name]];
};

////////////////////////////////////////////////THipoBank
THipoBank::~THipoBank(){
   for(UInt_t i=0;i<fVecItems.size();i++) delete fVecItems.at(i);
   for(UInt_t i=0;i<fVecI.size();i++) delete fVecI.at(i);
   for(UInt_t i=0;i<fVecF.size();i++) delete fVecF.at(i);
   fVecItems.clear();
   fVecI.clear();
   fVecF.clear();
   fItemsI.clear();
   fItemsF.clear();
}
void THipoBank::AddIntItem(TString name,Int_t id){
  fItemsI.push_back(name);
  fVecI.push_back(new vector<Int_t >);
  fItemMap[name]=fIndexI.size();
  fIndexI.push_back(id);
  fFloatMap[name]=kFALSE;
}
void THipoBank::AddFloatItem(TString name,Int_t id){
  fItemsF.push_back(name);
  fVecF.push_back(new vector<Float_t >);
  fItemMap[name]=fIndexF.size();
  fIndexF.push_back(id);
  fFloatMap[name]=kTRUE;
}
void THipoBank::InitTree(TTree* tree){
  cout<<"Init tree branch for "<<fName<<" "<<NInt()<<endl;
  //Make a branch for each integer item
  for(UInt_t in=0;in<NInt();in++){
    TString branchname=fName;
    branchname.ReplaceAll("::","_");
    branchname.ReplaceAll(" ","");
    tree->Branch(branchname+"_"+fItemsI.at(in),(fVecI.at(in)));
  }
  //Make a branch for each float item
  for(UInt_t fn=0;fn<NFloat();fn++){
    TString branchname=fName;
    branchname.ReplaceAll("::","_");
    branchname.ReplaceAll(" ","");
    tree->Branch(branchname+"_"+fItemsF.at(fn),(fVecF.at(fn)));
  }
  
}
void THipoBank::Print(){
  cout<<"Printing Hipo Bank"<<endl;
  cout<<"Name = "<<fName<<endl;
  cout<<"Group = "<<fGroup<<endl;
  cout<<"Int Items :"<<endl;
  for(UInt_t i=0;i<fItemsI.size();i++)
    cout<<fItemsI[i]<<endl;
  cout<<"Float Items :"<<endl;
  for(UInt_t i=0;i<fItemsF.size();i++)
    cout<<fItemsF[i]<<endl;
}
THipoItem* THipoBank::GetItem(TString name){
  if((std::find(fItemsI.begin(),fItemsI.end(),name)==fItemsI.end())&&(std::find(fItemsF.begin(),fItemsF.end(),name)==fItemsF.end())){
    cout<<"Error in THipoItem THipoBank::GetItem item "<<name<<" is not configured  in bank "<<GetName()<<endl;
    exit(0);
  }
  Int_t index=fItemMap[name];
  if(fFloatMap[name])
    fVecItems.push_back(new THipoItem(fVecF[index],&fEntry));
  else
    fVecItems.push_back(new THipoItem(fVecI[index],&fEntry));
  return fVecItems.back();
}
////////////////////////////////////////////////THipoBankParser
THipoBankParser::THipoBankParser(TString filename){
  //if filename="" it will just load all files in banksdefs/hipo dir
  //if filename="DATA:FT:..." will load file DATA.json from banksdefs/hipo
  //if filename="file.any" will load banks define in file.any

  if(filename.Contains("."))
    fBanksFile=new TMacro(filename);  
  if(filename.Length()==0){
    gSystem->Exec("rm banks.txt");
    gSystem->Exec(Form("cat $CLARA_HOME/plugins/clas12/etc/bankdefs/hipo/*json > banks.txt"));
    fBanksFile=new TMacro("banks.txt");
  }
  else{
    TObjArray* bankfiles=filename.Tokenize(":");
    gSystem->Exec("rm banks.txt");
    gSystem->Exec("touch banks.txt");
    for(Int_t ib=0;ib<bankfiles->GetEntries();ib++)
      gSystem->Exec(Form("cat $CLARA_HOME/plugins/clas12/etc/bankdefs/hipo/%s.json >> banks.txt",bankfiles->At(ib)->GetName()));
    fBanksFile=new TMacro("banks.txt");
    delete bankfiles;  
  }
  
}
Bool_t THipoBankParser::ConfigBank(THipoBank* bank){
  //Search for bank definition in file
  //Check if variable selection in name i.e bankname>>var1:var2:var3...
  TObjArray* selectname=bank->GetName().Tokenize(">>");
  if(selectname->GetEntries()>1)bank->SetName(((TObjString*)selectname->At(0))->GetName());
  //Now split the variable names seperated by :
  TObjArray* use_var=nullptr;
  if(selectname->GetEntries()>1)
    use_var=((TObjString*)selectname->At(1))->String().Tokenize(":");

  cout<<"THipoBankParser::ConfigBank "<<bank->GetName()<<endl;
  //Find the bank
  fPlace=0;
  FindNextLineWith(bank->GetName());

  //Find the group
  TString line=FindNextLineWith("group");
  line.ReplaceAll(",","");line.ReplaceAll(" ","");
  TObjArray* group=line.Tokenize(":");
  bank->SetGroup(((TObjString*)group->At(1))->String().Atoi());
  //Find the items
  FindNextLineWith("items");fPlace++;
  TList *lines=fBanksFile->GetListOfLines();
  while(dynamic_cast<TObjString*>(lines->At(fPlace))->String().Contains("name")){
    TObjArray* item=dynamic_cast<TObjString*>(lines->At(fPlace))->String().Tokenize(",");
    TString item_name=((TObjString*)item->At(0))->String();
    item_name.ReplaceAll("\"","");
    item_name.ReplaceAll(" ","");
    item_name.ReplaceAll("\t","");
    item_name.ReplaceAll("{name:","");
    if(use_var) //There was a variable selection given with bankname
      if(!use_var->FindObject(item_name)) {fPlace++;continue;} //this variable was not selected
    
    TString item_type=((TObjString*)item->At(2))->String();
    item_type.ReplaceAll("\"","");
    item_type.ReplaceAll("type:","");

    TString item_id=((TObjString*)item->At(1))->String();
    item_id.ReplaceAll("\"","");
    item_id.ReplaceAll("id:","");
    if(item_type.Contains("float"))
      bank->AddFloatItem(item_name,item_id.Atoi());
    else if(item_type.Contains("int"))
      bank->AddIntItem(item_name,item_id.Atoi());
    else cout<<"Warning :: THipoBankParser::ConfigBank "<<bank->GetName()<< " item of invalid type : "<<item_name<<" given as "<<item_type<<endl;
    //move to next line
    fPlace++;

    delete item; item=nullptr;
    
  }
  
  //bank->Print();

  delete group;group=nullptr;
  delete use_var;use_var=nullptr;
  delete selectname;selectname=nullptr;
  
  return kTRUE;
}
TString THipoBankParser::FindNextLineWith(TString linelike){
  TList *lines=fBanksFile->GetListOfLines();
  TObjString* thisline=0;
  Int_t count=0;
  for(count=fPlace;count<lines->GetEntries();count++){
    thisline=dynamic_cast<TObjString*>(lines->At(count));
    if(thisline->String().Contains(linelike))
      break;
  }
  if(count==lines->GetEntries()) {fPlace=count; return "";} //didn't find line go to end
  else fPlace=count; //get line number
  return thisline->String();
}
