
#ifndef THSFINALSTATE_h
#define THSFINALSTATE_h

#include <TString.h>
//#include <TNamed.h>
#include <TTree.h>
#include <vector>
#include <map>

#include "THSParticle.h"
#include "THSParticleIter.h"
#include "THSDataManager.h"
#include "THSKinematics.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>


//Developments :
//Add THSParticleIter class for combitorials
//To limit number of combitorials
//First loop over beam particles, assign status of other particles
//based on timing coincidence (0=reject;1=coincidence;2=random sample)
//Run check topology based on status


class THSFinalState{

 public :
  THSFinalState();
  virtual ~THSFinalState(){fDetIter.clear();};

  virtual void GetEvent(Long64_t uid){}//interface to different readers
  virtual Bool_t WorkOnEvent(){return kFALSE;};
  virtual void FinaliseEvent(){};
  virtual void InitEvent(){fGotCorrectOne=kFALSE;  fNPerm=0; };
  
  void SetDetParts(vector<THSParticle*> * dpp){frDetParts=dpp;}
  void SetGenParts(vector<THSParticle*> * dpp){frGenParts=dpp;}

  Int_t AddTopology(TString topo);
  vector<Int_t> GetTopology(Int_t i){return fTopo.at(i);};
  Int_t FindTopology(); //For current event
  Int_t FindInclusiveTopology(Int_t incType=-1E9);//Inclusive For current event
  void SetPermutate(){fTryPerm=kTRUE;} //Turn on permuations
  void SetInclusive(){fIsInclusive=kTRUE;} //Turn on permuations
  void SetGenerated(Bool_t gen=kTRUE){fIsGenerated=gen;}
  void SetUseChargePID(Bool_t gen=kTRUE){fUseChargePID=gen;}

 protected :
  Int_t fNDet=0;
  vector<vector<Int_t> >fTopo; //vector of topologies
  Int_t fNTopo=0; //Number of added topologies
  Int_t fCurrTopo=-1;
  UInt_t fNParts=0; //Number of particles detected
  UInt_t fNGen=0; //Number of generated particles
  Bool_t fGoodEvent = kTRUE;
  Long64_t fNUsedReadEvent=0; //Number of read events with at least 1 GoodEvent
  //Special behaviour for generated MC events
  Bool_t fIsGenerated=kFALSE;
  virtual void InitGenerated(){}
 //For simulated events flag with correct permutation
  Int_t fCorrect=kTRUE;
  Bool_t fGotCorrectOne=kFALSE;
  //Flag to do PID by charge not given pdg code
  Bool_t fUseChargePID=kFALSE;
  Bool_t fIsInclusive=kFALSE;
  
  //Combitorial
  Bool_t IsMissing(THSParticle* part){if(part->P4p()->Theta()==0){return kTRUE;}return kFALSE;}
  Bool_t fTryPerm=kFALSE;
  Int_t fNPerm=0;
  void InitParticles(); //Allocate particles to vector for this event
  void InitDetParts(); //Fill vector based on pdg
  void InitDetPartsCharge(); //Fill vector based on charge
  Bool_t PermutateParticles(); //Permutate particles in vectors
  
  vector<THSParticle*> fVecProtons;
  vector<THSParticle*> fVecPiPs;
  vector<THSParticle*> fVecPiMs;
  vector<THSParticle*> fVecKPs;
  vector<THSParticle*> fVecKMs;
  vector<THSParticle*> fVecEls;
  vector<THSParticle*> fVecPos;
  vector<THSParticle*> fVecPi0s;
  vector<THSParticle*> fVecGams;
  vector<THSParticle*> fVecPlus;
  vector<THSParticle*> fVecMinus;
  vector<THSParticle*> fVec0;
  vector<THSParticle*> fVecBeams;
  vector<Int_t> fDetTypes;

  vector<THSParticleIter*> fDetIter;  //Particle iterator for each topology

  Int_t fIterLevel=0;

    
  map< Int_t , vector< THSParticle* >* > fMapPDGtoParticle;
  //map< Int_t , Int_t > fParticleCounter; //count how many of each type of particle have been asked for this permutation

  Bool_t fIsPermutating0=kFALSE;
  Bool_t fIsPermutating1=kFALSE;

  TTree* fFinalTree=nullptr;
  
  //Kinematics calculator
  THSKinematics fKine;

  //Unique event ID
  Double_t fUID=0;
  
  //Detected Particles
  vector<THSParticle* >* frDetParts=nullptr;
  vector<THSParticle* >* frGenParts=nullptr;

  //Final Particles vector for saving in output tree
  vector<THSParticle*> fFinal;
 
  //Observables

  //Discriminators
  
 public :
  Bool_t IsGoodEvent(){return fGoodEvent;}
  virtual Bool_t IsPermutating(){ return fTryPerm*(fIsPermutating0+fIsPermutating1);}

  virtual void ProcessEvent();
  virtual Bool_t FSProcess();
  virtual void UserProcess();
  
  
  virtual void MatchWithGen(THSParticle* part);  
  virtual Bool_t IsCorrectTruth(THSParticle *part);
  virtual void FinalStateOutTree(TTree* tree){tree->Branch("Final",&fFinal);};
  virtual void CheckTruth();


  THSParticleIter* AddSelectToSelected(THSParticleIter* diter,Int_t ni,Int_t nsel,THSParticle* part0=nullptr,THSParticle* part1=nullptr,THSParticle* part2=nullptr,THSParticle* part3=nullptr,THSParticle* part4=nullptr,THSParticle* part5=nullptr,THSParticle* part6=nullptr,THSParticle* part7=nullptr,THSParticle* part8=nullptr,THSParticle* part9=nullptr);
  
  THSParticleIter* AddSelectToRemainder(THSParticleIter* diter,Int_t ni,Int_t nsel,THSParticle* part0=nullptr,THSParticle* part1=nullptr,THSParticle* part2=nullptr,THSParticle* part3=nullptr,THSParticle* part4=nullptr,THSParticle* part5=nullptr,THSParticle* part6=nullptr,THSParticle* part7=nullptr,THSParticle* part8=nullptr,THSParticle* part9=nullptr);

  THSParticleIter* AddSelectXofY(Int_t ni,Int_t nsel,THSParticle* part0=nullptr,THSParticle* part1=nullptr,THSParticle* part2=nullptr,THSParticle* part3=nullptr,THSParticle* part4=nullptr,THSParticle* part5=nullptr,THSParticle* part6=nullptr,THSParticle* part7=nullptr,THSParticle* part8=nullptr,THSParticle* part9=nullptr);

  THSParticleIter* AddPermutateToSelected(THSParticleIter* diter,THSParticle* part0=nullptr,THSParticle* part1=nullptr,THSParticle* part2=nullptr,THSParticle* part3=nullptr,THSParticle* part4=nullptr,THSParticle* part5=nullptr,THSParticle* part6=nullptr,THSParticle* part7=nullptr,THSParticle* part8=nullptr,THSParticle* part9=nullptr);

  THSParticleIter* AddPermutateToRemainder(THSParticleIter* diter,THSParticle* part0=nullptr,THSParticle* part1=nullptr,THSParticle* part2=nullptr,THSParticle* part3=nullptr,THSParticle* part4=nullptr,THSParticle* part5=nullptr,THSParticle* part6=nullptr,THSParticle* part7=nullptr,THSParticle* part8=nullptr,THSParticle* part9=nullptr);

  THSParticleIter* AddPermutate(THSParticle* part0=nullptr,THSParticle* part1=nullptr,THSParticle* part2=nullptr,THSParticle* part3=nullptr,THSParticle* part4=nullptr,THSParticle* part5=nullptr,THSParticle* part6=nullptr,THSParticle* part7=nullptr,THSParticle* part8=nullptr,THSParticle* part9=nullptr);
  
  THSParticleIter* CreateParticleIter(Int_t topo,vector<THSParticle*> *parts,Int_t Nsel);
  

    
};

/* inline void THSFinalState::SetMinus(THSParticle *particle,Int_t index){ */
/*   if(index==-1)index=fParticleCounter[-1E6]++; */
/*   particle->CopyParticle(fVecMinus.at(index),kFALSE); */
/*   particle->TakePDGMass(); */
/* } */
/* inline void THSFinalState::SetPlus(THSParticle *particle,Int_t index){ */
/*  if(index==-1)index=fParticleCounter[1E6]++; */
/* /\*   fParticleCounter[1E6]++; *\/ */
/* /\* cout<<"plus "<<index<<endl; *\/ */
/*   particle->CopyParticle(fVecPlus.at(index),kFALSE); */
/*   particle->TakePDGMass(); */
/* } */
/* inline void THSFinalState::Set0(THSParticle *particle,Int_t index){ */
/*  if(index==-1)index=fParticleCounter[0]++; */
/* //  fParticleCounter[0]++; */
/*   particle->CopyParticle(fVec0.at(index),kFALSE); */
/*   particle->TakePDGMass(); */
/* } */
/* inline void THSFinalState::SetFromPDG(THSParticle *particle,Int_t index){ */
/* if(!fUseChargePID){ */
/* if(index==-1)index=fParticleCounter[particle->PDG()]++; */
/* //  fParticleCounter[particle->PDG()]++; */
/* particle->CopyParticle(fMapPDGtoParticle[particle->PDG()]->at(index),kFALSE); */
/* particle->TakePDGMass(); */
/* } */
/*  else{ */
/* Int_t charge=particle->Charge(); */
/* if(charge) */
/*   if(charge>0) */
/*     SetPlus(particle,index); */
/*   else */
/*     SetMinus(particle,index); */

/*  else */
/*    Set0(particle,index); */
   
/* } */
/* } */

#endif //ifdef THSFINALSTATE
