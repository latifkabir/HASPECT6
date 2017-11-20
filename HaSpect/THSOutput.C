//--Author      DI Glazier 30/06/2014
//--Rev
//--Update
//--Description
//HASPECT Event Reconstruction
//THSOuput
//Controls the output of the TSelector analysis
//Creates THSParticle branches in ouput tree
//Writes source code to output file
//Writes either one combined file, or one file for each input file 
//fStrParticles should be passed a string of particles to be
//searched for and output in the tree, e.g. "pi+:pi-:pi+:neutron"
//the particles will be ordered in fDetParticle[i] for further analysis
//ordered by type then decreasing momentum, see output message for connecting
//branch names to fDetParticle indices
//After making a skeleton selector class you can add inheritance from this class
//to obtain its functionality

#include <TObject.h>
#include <TProofOutputFile.h>
#include <TDatabasePDG.h>
#include <TSystem.h>
#include <TEntryList.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TH3.h>
#include <TKey.h>
#include <TProof.h>
#include <TError.h>
#include <TTreeIndex.h>
#include <iostream>
#include <algorithm>
#include "THSOutput.h"

using namespace std;

THSOutput::~THSOutput(){
  if(fOutTree) SafeDelete(fOutTree);
  if(fFile){
    SafeDelete(fFile);
    //SafeDelete(fProofFile);
  }
  //if(fListOfFiles) SafeDelete(fListOfFiles);
  // if(fParentFile) {fParentFile->Close(); SafeDelete(fParentFile);fParentTree=0;}
}
void THSOutput::HSBegin(TList* input,TList* output){
  fOutName = gSystem->Getenv("HSOUT"); //outfile name must be set as an enviroment variable
  if(fOutName==TString("")){fMessages.Warning("THSOutput::HSBegin","No output destination given via HSOUT enviroment variable...");}
    
  PrepareOutDir();  //Make directory for output root files
  fSelInput=input; //connect to the selector input list
  fSelOutput=output; //connect to the selector output list
  //Copy all the source code macros
  fCodeList=new TList();
  fCodeList->SetOwner();
  fCodeList->SetName("CodeList");
  ImportSysDirtoList(gSystem->pwd(),fCodeList);//the directory the user analysis code is in
  ImportSysDirtoList(gSystem->Getenv("HSANA"),fCodeList); //the directory the HasPEct core code is
  //if(gProof)fSelInput->Add(fCodeList);//add to input so can give to Slaves, if fInput doesn't exit then not running on PROOF and this fCodeList is OK
 }
void  THSOutput::HSSlaveBegin(TList* input,TList* output){
  //need to reconnect with selector in case running on proof
  fOutName = gSystem->Getenv("HSOUT"); //outfile name must be set as an enviroment variable
  if(input) {//input only exists if proof Proof loses env variables (in root 6)
    if(input->FindObject("HSOUT")) fOutName=input->FindObject("HSOUT")->GetTitle();
  }
 if(!fSelOutput){
    fSelInput=input; //connect to the selector input list
    fSelOutput=output; //connect to the selector output list
  }
  //Make list for storing names of analysed files for use in terminate
  fListOfFiles=new TList();
  fListOfFiles->SetOwner(kTRUE);
  fListOfFiles->SetName("HSInFiles");
  fSelOutput->Add(fListOfFiles);
  
  //Start the ID from the offset if
  fgID=fgIDoff;	
};
void THSOutput::HSNotify(TTree* tree){
  cout<<"THSOutput::HSNotify(TTree* tree) "<<tree->GetDirectory()->GetName()<<endl;
  //Function that looks after the output file
  //Needs to be called in TSelector::Notify()
  fCurTree=tree;
  fListOfFiles->Add(new TNamed(gSystem->BaseName(tree->GetDirectory()->GetName()),gSystem->DirName(tree->GetDirectory()->GetName())));
  //Check if gID exists in current tree, if not will start saving it now
  if(!fCurTree->GetBranch("UID")) fSaveID=kTRUE;
  if(fOutName.EndsWith(".root")){
    //Only want to do some things once when writing 1 file only
     if(fStepName.Length()==0){
       InitOutFile(fCurTree); //initialise output file
       return;
     }
    else return; //only saving one combined file
   }
 //Case  making many output files
  FinishOutput(); //close the last file
  InitOutFile(fCurTree);   //start the new file
  return;

}
void THSOutput::HSProcessStart(Long64_t entry){
  fEntry=entry;
  //if(fSaveID) fgID++;
  if(!fSaveID){
    fCurTree->GetBranch("UID")->GetEntry(fEntry); //make sure get ID branch to write to new file
    fgID=fCurTree->GetLeaf("UID")->GetValue();
  }
}
void THSOutput::HSProcessFill(){
  //if(fSaveID) fgID=fEntry;
  if(fSaveID) fgID++;
  if(fOutTree) fOutTree->Fill();
  //if(!fSaveID)fgID+=fgIDoff; //in case multiple calls to fill change ID
}
void THSOutput::HSSlaveTerminate(){
  FinishOutput();
}
void THSOutput::HSTerminate(){

  
  //TURN OFF saving of code and entrylists
  //LEave code here as might prove useful in the future
  //Write analysis code to a file.
  //if 1 root file given as output just store there
  //else if writing to directory create a new file HSCode.root
  //fListOfFiles=dynamic_cast<TList*>(fSelOutput->FindObject("HSInFiles"));
  TString CodeFileName;
  TFile *outfile=0;
  TDirectory* savedir=gDirectory;
  if(fOutName.EndsWith(".root")) {CodeFileName=fOutName;outfile=TFile::Open(fOutName,"update");}
  else {CodeFileName="/HSCode.root"; outfile=TFile::Open(fOutName+CodeFileName,"recreate");}
  //Also want input file so can copy code from there
   TString InDirName=gSystem->Getenv("HSIN");
   TString InFileName;
   TFile* infile=0;
   //Check if HSCode.root exists in input directory
   gErrorIgnoreLevel = kFatal; //temp turn off errors as it is OK if this file does not exist
   infile=TFile::Open(InDirName+"/HSCode.root");
   gErrorIgnoreLevel = 0;
   if(!infile){//No HSCode file, assume in input file, if not will create first save dir
     infile=TFile::Open(TString(fListOfFiles->At(0)->GetTitle())+"/"+fListOfFiles->At(0)->GetName()); //Note need to take 1 as add dummy entry Start for proof workers that don't get data
     if(fListOfFiles->GetEntries()>1)Warning("HSTerminate()"," Just going to find previous code from first file, there were %d",fListOfFiles->GetEntries());
   }
   //Copy code from current directory and add to previous code in input file
   CopyCode(outfile,infile);
   outfile->cd();
   WriteListtoFile(fStepDir);
   infile->Close();
   SafeDelete(infile);

   outfile->Close();
   SafeDelete(outfile);
   //make sure new tree events are back in fgID order
   // if(gProof&&(fSort==kTRUE)){
   //   TIter next(fSelOutput);
   //   TProofOutputFile* elpofile=0;
   //   TObject* outo=0;
   //   TFile* infile=0; //pointer to input file
   //   //iterate over any proof files which are in the ouput list
   //   while((outo=dynamic_cast<TObject*>(next()))){
   //     if((elpofile=dynamic_cast<TProofOutputFile*>(outo))){
   // 	 TFile* elfile = elpofile->OpenFile("UPDATE");
   // 	 //First sort tree to regain original ordering
   // 	 TIter fnext(elfile->GetListOfKeys());
   // 	 TKey* fkey=0;
   // 	 //Look for a trees in saved file
   // 	 while ((fkey = (TKey*)fnext())){
   // 	   if(TString(fkey->GetClassName())==TString("TTree")){
   // 	     SortTree(dynamic_cast<TTree*>(elfile->Get(fkey->GetName())));
   // 	   }
   // 	 }
   //     }
   //   }
   // }

   //in case of proof must add uid at end as cannot synch fgID in proof
   if(gProof&&(fSaveID==kTRUE)){
     cout<<"THSOutput::Terminate going to add UID branch, may take a while"<<endl;
     TIter next(fSelOutput);
     TProofOutputFile* elpofile=0;
     TObject* outo=0;
     TFile* infile=0; //pointer to input file
     //iterate over any proof files which are in the ouput list
     fgID=fgIDoff;
     while((outo=dynamic_cast<TObject*>(next()))){
       if((elpofile=dynamic_cast<TProofOutputFile*>(outo))){
	 TFile* elfile = elpofile->OpenFile("UPDATE");
	 //First sort tree to regain original ordering
	 TIter fnext(elfile->GetListOfKeys());
	 TKey* fkey=0;
	 //Look for a trees in saved file
	 while ((fkey = (TKey*)fnext())){
	   if(TString(fkey->GetClassName())==TString("TTree")){
	     TTree* tree=dynamic_cast<TTree*>(elfile->Get(fkey->GetName()));
	     if(tree->GetBranch("UID")) continue; //aready exists
	     TBranch* branch=tree->Branch("UID",&fgID,"UID/D");
	     if(!branch) continue;
	     Long64_t id=(Long64_t)fgID;
	     fgID=id;
	     for(Long64_t i=0;i<tree->GetEntries();i++){
	       branch->Fill();
	       id++;
	       fgID=id;   
	     }
	     tree->Write(tree->GetName(),TObject::kOverwrite);
	     elfile->Close();
	   }
	 }
       }
     }
   }

  //   OrganiseSavedHists(fSelOutput);
   
   return;
 
}
void THSOutput::FinishOutput(){
   // Write the output to the file
  if (fFile) {
    Bool_t cleanup = kFALSE;
    TDirectory *savedir = gDirectory;
    fFile->cd();
    if (!fOutTree) {//no output tree 
      fSelOutput->Add(fProofFile); //give this proof file to output for merging in case there are histograms
      fProofFile->Print();
     }
    else  if (fOutTree->GetEntries() > 0) {//if tree has entries save it
      fOutTree->Write(0, TObject::kOverwrite);  //write ouput tree
      fOutTree->Reset();  //remove saved entries
      fSelOutput->Add(fProofFile); //give this proof file to output for merging
      fProofFile->Print();
      fOutTree->SetDirectory(0);
    } 
    else {//an empty tree so don't save it
      cleanup = kTRUE;
      fOutTree->SetDirectory(0);
    }
    OrganiseSavedHists(fFile,fSelOutput);
    gDirectory = savedir;
    fFile->Close(); //Close this file
    // Cleanup, if needed
    if (cleanup) {
      TUrl uf(*(fFile->GetEndpointUrl()));
      SafeDelete(fFile);
      gSystem->Unlink(uf.GetFile());
      SafeDelete(fProofFile);
      //	SafeDelete(fOutTree);
    }
  }
}

void THSOutput::SortTree(TTree* tree){
  //reorder events in the tree so the order is preserved
  //typically PROOF will reorder events, but order preservation is
  //useful for connecting with parent trees
  //Make sure the file you are writing to is the current directory
  
  //order the events based on the global ID variable, which came from the original tree
  fMessages.Info("THSOutput::SortTree","Sorting tree You may want to disable this?");
  if(!tree) return;
  if(!tree->GetBranch("UID")) return;
  tree->BuildIndex("(Long64_t)UID");
  TTreeIndex *index = (TTreeIndex*)tree->GetTreeIndex();
  TTree* cltree=tree->CloneTree(0); //create empty tree with branch adresses set
  cltree->SetAutoSave();
  // cltree->CopyEntries(tree,);//check if this works with index without needing the rest
  for( int i =  0; i < index->GetN() ; i++ ) {
    Long64_t local = tree->LoadTree( index->GetIndex()[i] );
    tree->GetEntry(local);
    cltree->Fill(); //fill as ordered by the build index
  }
  cltree->Write("",TObject::kOverwrite); //replace current tree
  SafeDelete(index);
}

void THSOutput::PrepareOutDir(){
  //Make the outpur directory if requested
  //If it exists exit so we do not copy over any files
  if(!fOutName.EndsWith(".root")&&gSystem->MakeDirectory(fOutName)!=0){// onlymake directory first time and only if it doesn't exist
      fMessages.Error("THSOutput::PrepareOutDir()", "Directory given  is not approriate, it either exists or its parent doesn't, I would rather not overwrite it in case it contains files you need. It is advised t set fOutput=\"outfile or outdir\" in your selector constructor, this way it can be seen by all the worker servers.");
      cout<<"Directory Name = "<<fOutName<<endl;
      gSystem->Exit(1);
    }
}
void THSOutput::InitOutFile(TTree* chain){
   //check fOutName to see if making 1 or many files
  if(fFile&&fOutName.EndsWith(".root"))return;//only saving one file

  //create a new output file to save tree to
  //can be called as part of notify to switch filename to match input file
   TString ofname;
  if(fOutName.EndsWith(".root")) ofname=fOutName;
  else  {
    ofname =gSystem->BaseName((fCurTree->GetCurrentFile()->GetName()));
    ofname.Prepend(fOutName+"/");
  }
  TDirectory* savedir=gDirectory;
  if(fFile)	SafeDelete(fFile);
   //if(fProofFile)	SafeDelete(fProofFile);
  fMessages.Info("InitOutFile","Making new proof file %s",ofname.Data());
  fProofFile = 0;
  fFile = 0;
  fProofFile = new TProofOutputFile(ofname,"M");
  fFile = fProofFile->OpenFile("RECREATE");
  if (fFile && fFile->IsZombie()) SafeDelete(fFile);
  // Cannot continue
   if (!fFile) {
     Info("SlaveBegin", "could not create '%s': instance is invalid!", fProofFile->GetName());
     return;
   }
   //if it exists give the tree to the file
   
   gDirectory = savedir;

}
void THSOutput::InitOutTree(){
  //Function to set tree file and make sure fgID branch exists
  //used to be part of InitOutTree but split because needs to be called after fOutTree create
   if(fOutTree&&fFile){
     fOutTree->SetDirectory(fFile);
     fOutTree->AutoSave();
     //if proof will add UID in terminate so it is incremented properly
     if(!fOutTree->GetBranch("UID")&&!gProof)fOutTree->Branch("UID", &fgID, "UID/D");
     //if(!fSaveID)//copy existing global ID
       // fOutTree->SetBranchAddress("UID",fCurTree->GetBranch("UID")->GetAddress());
   }
}
void THSOutput::CopyCode(TDirectory* curDir,TDirectory* prevDir){

  //Function to copy all source codes to the output file (curDir)or outputlist
  //If any code exists in the input file (prevDir) this is copied as a sublist
  //in the new list.
  //The code is saved as TMacros ans are stored in a TList (easier to use than TDirectory)
  TDirectory* savedir=gDirectory;
  TDirectory* prevStep=0;
  //Check to see if code already exists in parent file
  TIter next(prevDir->GetListOfKeys());
  TKey* key=0;
  //Look for a directory in the input root file which includes HSStep
    //it will be copied to the new output root file
  while ((key = (TKey*)next())) if(TString(key->GetName()).Contains("HSStep")){
      prevStep=dynamic_cast<TDirectory*>(key->ReadObj());
      break;
    }
  if(prevStep){
    //We have the saved last step, we will save this in the current source step
    //so as to contain the full analysis chain
    //the current step will have the HSStep_number incremented by 1
    fStepName= prevStep->GetName();
    TString prevstepi=TString(fStepName(fStepName.Index("_")+1,fStepName.Length()-fStepName.Index("_"))); //get the number as a string
    Int_t prevStepN=prevstepi.Atoi(); //convert it to an int
    TString stepi;
    stepi=stepi.Itoa(prevStepN+1,10); //add 1 and convert back to string
    fStepName.ReplaceAll(prevstepi,stepi); //now have the new step name
  }
  else fStepName="HSStep_0";

  //create list of current source, prepare to add previous code form in file
  //if(fStepDir) delete fStepDir; //cleanup previous step directory
  fStepDir=(TList*)fCodeList->Clone();
  //fStepDir->SetOwner();
  fStepDir->SetName(fStepName);
  
  //If there was a previous step copy its source to the new step list
  if(prevStep) fStepDir->Add(CopyDirtoList(prevStep));
  //Write the source code to the output file curDir
  // curDir->cd();
  // WriteListtoFile(fStepDir);
  savedir->cd();

}
void THSOutput::ImportSysDirtoList(const char *dirname,TList* list) {
  //based on $ROOTSYS/tutorials/io/importCode.C       
  //  char *slash = (char*)strrchr(dirname,'/');
  TDirectory *savdir = gDirectory;
  void *dirp = gSystem->OpenDirectory(dirname); //open directory on system
  if (!dirp) return;
  char *direntry;
  Long_t id, size,flags,modtime;
  //loop on all entries of this directory
  while ((direntry=(char*)gSystem->GetDirEntry(dirp))) {
    TString afile = Form("%s/%s",dirname,direntry);
    gSystem->GetPathInfo(afile,&id,&size,&flags,&modtime);
    if (direntry[0] == '.')             continue; //forget the "." and ".." special cases
    if (!strcmp(direntry,"CVS"))        continue; //forget some special directories
    if (!strcmp(direntry,"htmldoc"))    continue;
    if (strstr(dirname,"root/include")) continue;
    if (strstr(direntry,"G__"))         continue;
    if (strstr(direntry,".root"))         continue;
    if (strstr(direntry,"~"))         continue;
    //Copy all of these types of files
    if (strstr(direntry,".c")    ||
	strstr(direntry,".h")    ||
	strstr(direntry,".hh")    ||
	strstr(direntry,".dat")  ||
	strstr(direntry,".py")   ||
	strstr(direntry,".txt")   ||
	strstr(direntry,".cxx")   ||
	strstr(direntry,".C")) {
      TMacro *m = new TMacro(afile);
      m->SetName(gSystem->BaseName(afile));
      list->Add(m);
  
    } else {
      if (flags != 3)                     continue; //must be a directory
 	//*************remove copying of lower level directories
       //we have found a valid sub-directory. Process it
      // TList* list1=new TList(); //create new sublist
      // list1->SetName(afile);
      // list1->SetOwner();
      // ImportSysDirtoList(afile,list1); //write the sub directory to the sublist
      // list->Add(list1); //add sub list
    }
  }
  gSystem->FreeDirectory(dirp);
  savdir->cd();
}

TList* THSOutput::CopyDirtoList(TDirectory *source) {
  //copy all objects and subdirs of directory source as a tlist
  //This should be easy, but for some reason is not. 
  //This code is based on tutorials/io/copyFiles.C   
  TList* list=new TList();
  list->SetName(source->GetName());
  list->SetOwner();
  TDirectory *savdir = gDirectory;
  //loop on all entries of this directory
  TKey *key=0;
  TIter nextkey(source->GetListOfKeys());
  while ((key = (TKey*)nextkey())) {
      const char *classname = key->GetClassName();
      TClass *cl = gROOT->GetClass(classname);
      if (!cl) continue;
      if (cl->InheritsFrom(TDirectory::Class())) {
       source->cd(key->GetName());
         TDirectory *subdir = gDirectory;
         list->Add(CopyDirtoList(subdir));
      } else if (cl->InheritsFrom(TTree::Class())) {
         TTree *T = (TTree*)source->Get(key->GetName());
	 TTree *newT = T->CloneTree(-1,"fast");
         list->Add(newT);
      } else {
         source->cd();
         TObject *obj = key->ReadObj();
	 list->Add(obj);
      }
  }
  savdir->cd();
  return list;
}
void THSOutput::WriteListtoFile(TList* list0){
  //Contents of list will be written in current file  
  //loop on all entries of this list
  if(!gDirectory->IsWritable()) return;
  //make directory for this list
  TDirectory *saveDir = gDirectory;
  TDirectory *curDir = gDirectory->mkdir(list0->GetName());
  curDir->cd();
  TKey *key;
  TIter nextkey(list0);
  //loop through list and  make new directories
  while ((key = (TKey*)nextkey())) {if(TString(key->GetName()).Contains("HSStep")){//key->GetClassName() causes a crash...
      TList* list1=(TList*)list0->FindObject(key->GetName());
      WriteListtoFile(list1); //write sublist to a subdirectory
    } 
    //write the objects in this list
    else key->Write(0, TObject::kOverwrite);
  }
  // list0->Write(0, TObject::kOverwrite);
  //back to original file directory
  saveDir->cd();
}

Int_t THSOutput::CountChar(TString tests,char testc){
  //based on TString function of same name which for some reason takes an int not a char
  Int_t count = 0;
  Int_t len   = tests.Length();
  const char *data  = tests.Data();
  for (Int_t ccn = 0; ccn < len; ccn++)
    if (data[ccn] == testc) count++;
  
  return count;
}
