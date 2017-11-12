#ifndef THIPO_h
#define THIPO_h

#include <TString.h>
#include <TMacro.h>
#include <TTree.h>
#include <TFile.h>
#include <vector>
#include <map>
#include <algorithm>

#include "reader.h"
#include "record.h"
#include "event.h"

using namespace std;

class THipoBank;

class THipoBankParser  {

  public:
  THipoBankParser(){};
  THipoBankParser(TString filename);
  virtual ~THipoBankParser(){if(fBanksFile) delete fBanksFile;};
  
  Bool_t ConfigBank(THipoBank* bank);
  TString FindNextLineWith(TString linelike);
 protected:
  TMacro *fBanksFile=nullptr;
  Int_t fPlace=0;
  
  private:

};
class THipoItem  {

  public:
 THipoItem(vector<Int_t >* vec,Int_t *entry) : fBankEntry(*entry){fItemI=vec;fIsFloat=kFALSE;};
 THipoItem(vector<Float_t >* vec,Int_t *entry) : fBankEntry(*entry) {fItemF=vec;};
  virtual ~THipoItem(){};

  Float_t Val(UInt_t entry){
    if(fIsFloat)
      return fItemF->at(entry);
    return fItemI->at(entry);
  }
  Int_t ValI(UInt_t entry){
    return fItemI->at(entry);
  }
  Float_t Val(){return Val(fBankEntry);}
  Int_t ValI(){if(!fIsFloat)return ValI(fBankEntry);return (Int_t) Val(fBankEntry);}
  
  UInt_t Size(){if(fIsFloat) return fItemF->size();return fItemI->size();}
  Bool_t FindEntry(Float_t val);
  void Reset(){fPos=0;}
 protected:
  
  private:
  vector<Int_t >* fItemI=nullptr;
  vector<Float_t >* fItemF=nullptr;
  Bool_t fIsFloat=kTRUE;
  Int_t &fBankEntry;//Points to THipoBank fEntry
  Int_t fPos=0;

};
inline Bool_t THipoItem::FindEntry(Float_t val){
  
  if(fIsFloat){
    if(fItemF->size()==0) return kFALSE;
    auto it = std::find(fItemF->begin()+fPos,fItemF->end(),val);
    fPos= it-fItemF->begin();
    if(fPos>=(Int_t)fItemF->size())fPos=-1;
    
  }
  else{
    if(fItemI->size()==0) return kFALSE;
    auto it = std::find(fItemI->begin()+fPos,fItemI->end(),(Int_t)val);
    fPos= it-fItemI->begin();
    if(fPos>=(Int_t)fItemI->size()) fPos=-1;
  }
  if(fPos>-1) {
    fBankEntry=fPos;
    fPos++;
    return kTRUE;
  }
  else{
    fPos=0;
    return kFALSE;
  }
}
class THipoBank  {

 public:
  THipoBank(){};
  THipoBank(TString name):fName(name){};
  virtual ~THipoBank();

  void SetGroup(Int_t group){fGroup=group;}
  void AddIntItem(TString name,Int_t id);
  void AddFloatItem(TString name,Int_t id);

  void SetName(TString name){fName=name;}
  TString GetName(){return fName;}
  Int_t Group(){return fGroup;}

  UInt_t NInt(){return fVecI.size();}
  UInt_t NFloat(){return fVecF.size();}

  UInt_t IdxI(UInt_t ii){return fIndexI.at(ii);}
  UInt_t IdxF(UInt_t ii){return fIndexF.at(ii);}

  //Note in readItems hipoEvent
  void ReadItemI(hipo::event *event,UInt_t ii);
  void ReadItemF(hipo::event *event,UInt_t ii);

  void ClearEvent();
  void Print();
  void InitTree(TTree* tree);
  //Iterate over entries for each event in item name
  Bool_t NextEntry();
  Int_t GetEntry(){return fEntry;}
  THipoItem* GetItem(TString name);
  
  protected:
  TString fName;
  Int_t fGroup=0;
  vector<vector<Int_t>* > fVecI;
  vector<vector<Float_t>* > fVecF;
  vector<UInt_t> fIndexI;
  vector<UInt_t> fIndexF;
  vector<TString> fItemsI;
  vector<TString> fItemsF;
  vector<THipoItem* > fVecItems;
  
  map<TString,Int_t > fItemMap; 
  map<TString,Bool_t > fFloatMap; 

    
  Int_t fEntry=-1;
  
  private:


};
inline   void THipoBank::ClearEvent(){
  for(UInt_t i=0;i<fVecI.size();i++) fVecI.at(i)->clear();
  for(UInt_t i=0;i<fVecF.size();i++) fVecF.at(i)->clear();
  for(UInt_t i=0;i<fVecItems.size();i++) fVecItems.at(i)->Reset();
  fEntry=-1;
}
inline Bool_t THipoBank::NextEntry(){
  //try to iterate using first integer vector
  if((Int_t)fVecI[0]->size()>fEntry+1){
    fEntry++;
    return kTRUE;
  }
  if(fEntry>-1) return kFALSE;
  //No integer vectors perhaps float
  if((Int_t)fVecF[0]->size()>fEntry+1){
    fEntry++;
    return kTRUE;
  }
  return kFALSE;
}

inline  void THipoBank::ReadItemI(hipo::event *event,UInt_t ii){
  vector<Int_t> vec=event->getInt(fGroup,IdxI(ii));
  for(UInt_t iv=0;iv<vec.size();iv++)
    (fVecI.at(ii))->push_back(vec.at(iv));
}
inline  void THipoBank::ReadItemF(hipo::event *event,UInt_t ii){
  vector<Float_t> vec=event->getFloat(fGroup,IdxF(ii));
  for(UInt_t iv=0;iv<vec.size();iv++)
    (fVecF.at(ii))->push_back(vec.at(iv));
}

class THipo  {

  public:
  THipo();
  THipo(TString filename);
  virtual ~THipo();
  
  protected:
  hipo::reader fReader;
  hipo::record fRecord;
  hipo::event fHipoEvent;

  vector<TString> fInFiles;

  TTree* fOutTree=nullptr;
  TFile* fOutFile=nullptr;
  
  vector<THipoBank*> fBanks;
  THipoBankParser *fParser=nullptr;
  THipoBank* fCurBank=nullptr;

  Long64_t fNRecEvents=0;
  Long64_t fRecEventCount=0;
  
  Long64_t fNRecords=0;
  Long64_t fRecordCount=0;

  Long64_t fNEvents=0;
  UInt_t fFilei=0;

  Bool_t fUseCLDirs=kFALSE;
  TString fOutDirName;
  
  map<TString,Int_t> fBankMap;
  
 private:
 public:
  Long64_t ProcessFile(TString filename);
  Bool_t  NextFile();
  void AddCommandLineFiles();
  void AddCommandLineDirs();
  void ToRoot();
  void ToRoot(TString filename);
  void InitFile(TString filename);
  Bool_t ConfigBank(TString bankname);
  void ConfigAllBanks();
  THipoBank* GetBank(TString name);
  void ReadItemI(UInt_t ii){fCurBank->ReadItemI(&fHipoEvent,ii);}
  void ReadItemF(UInt_t ii){fCurBank->ReadItemF(&fHipoEvent,ii);}

   Bool_t NextEvent();
   Bool_t NextRecord();

   void InitOutput(TString filename);
   TString GetOutDirName(){return fOutDirName;}
   void CloseOutput();
};

#endif
