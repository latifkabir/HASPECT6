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
//After making a skeleton selector class you can add inheritance from this class
//to obtain its functionality

#ifndef THSOutput_h
#define THSOutput_h

#include <TFile.h>
#include <TChain.h>
#include <TList.h>
#include <TString.h>
#include <TDirectory.h>
#include <TMacro.h>
#include <vector>
#include "THSHisto.h"

class TProofOutputFile;
class THSHisto;

class THSOutput :  public THSHisto{
 protected :
  TFile            *fFile=nullptr;      //Output file
  TProofOutputFile *fProofFile=nullptr; // For optimized merging of the ntuple
  TTree            *fOutTree=nullptr;  //Output tree
  TTree            *fCurTree=nullptr;  //Current input tree
  TEntryList* fEntryList=nullptr; //Entry List prooduced during selector Process
  TString fOutName;   //name of ouput file or directory
  TString fStrParticles; //string containing detected particles name
  TString fStepName; //Name of directory in fFile to write source code to
  TList* fStepDir=nullptr; //Directory containing source code
  TList* fCodeList=nullptr;  //pointer to the list of source code
  TList* fListOfFiles=nullptr; //list of analysed file names
  Long64_t fEntry=0; //current entry number in Chain
  TObject fMessages; //For access to Info,Error,...
  Bool_t fSort=kFALSE; //option to reorder events back to inital order in PROOF (which will randomise tree order)
 
  Bool_t fSaveID=kFALSE; //Ony save ID if this is the first THSOutput tree, subsequent trees will aready contain this branch
  Double_t fgID=0; //global event ID number, should be set in first instance of THSOutput and preserved through further steps. Required to synchornise differnt trees as PROOF does not preserve event ordering
  Long64_t fgIDoff=0; //starting offset for ID number
 public :
 THSOutput() {}   
  //THSOutput() : fFile(0), fProofFile(0), fOutTree(0),fCurTree(0), fEntryList(0),fStepDir(0),fSelInput(0),fSelOutput(0), fCodeList(0),fListOfFiles(0) {fSort=kFALSE;fSaveID=kFALSE;fgID=0;fgIDoff=0;}   
  virtual ~THSOutput();

  // virtual void InitOutput(); //Configure the output file and tree 
  virtual void FinishOutput();
  //virtual void NextFile(TTree* chain);  //call function when switch files in chain
  virtual void InitOutFile(TTree* chain);
  void InitOutTree();
  void PrepareOutDir(); //create output file directory
  Int_t CountChar(TString tests,char testc);
  void CopyCode(TDirectory* curDir=0,TDirectory* prevDir=0);  //Give new file,curDir to copy to, and and previous analysis steps, prevDir
  void ImportSysDirtoList(const char *dirname,TList* list);  //look through the given OS directory and save source code files

  //Utility functions for copying objects from directories to lists 
  //and writing from lists to directories
  //these were created as it proved too much effort to use Tdirectory throughout 
  TList* CopyDirtoList(TDirectory *source); //write source to list
  void WriteListtoFile(TList* list0); //write to current directory

  //Mimic Selector functions
  void HSBegin(TList* input,TList* output);
  void HSSlaveBegin(TList* input,TList* output);
  virtual void HSProcessStart(Long64_t entry);
  virtual void HSProcessFill();
  virtual void HSSlaveTerminate();
  virtual void HSNotify(TTree* tree);
  virtual void HSTerminate();

  void SortTree(TTree* tree);
  void SetSort(Bool_t sort=kTRUE){fSort=kTRUE;}
  
  TTree* GetOutTree() {return fOutTree;}

  void SetIDoff(Long64_t off) {fgIDoff=off;}
  //Setter functions may be needed if using without TSelector
  void SetParticles(TString sp){fStrParticles=sp;}
};
#endif
