#ifndef THSPARTICLEITER_h
#define THSPARTICLEITER_h

#include "THSCombitorial.h"
#include "THSParticle.h"

class THSCombitorial;

class THSParticleIter{
  
 public :
  THSParticleIter(){};
  virtual ~ THSParticleIter(){delete fExtraParticle;if(fInnerIter) delete fInnerIter;}//if(fSelIter) delete fSelIter; if(fRemIter) delete fRemIter;};
  
 private:
  Int_t fPDG=0; //PDG code for particles in this iterator
  THSCombitorial fCombi; //handler for the combination/permuations
  THSParticleIter *fInnerIter=nullptr;  //recursive iterator
  THSParticleIter *fSelIter=nullptr;  //selected iterator
  THSParticleIter *fRemIter=nullptr;  //remaining iterator
  vector<THSParticle*> *fAllParticles=nullptr; //vector of particles for this event
  vector<THSParticle*>  fSelected; //combitorial for selected particles
  vector<THSParticle*>  fRemainder; //combitorial for unselected particles
  
  vector<THSParticle*>  fEvParts; //particles to be used in event
  
  UInt_t fParti=0;// index of next particle
  Bool_t fUseCombi=kFALSE;
  Bool_t   fDoneSelRem=kTRUE;
  Int_t fNSel=0;
  Int_t fNIdentical=1;

  THSParticle* fExtraParticle=new THSParticle(-1);
  
 public:
  THSCombitorial* GetCombi(){return &fCombi;}
  
  void SetNSel(Int_t n){fNSel=n;}
  void SetNIdentical(Int_t n){fNIdentical=n;}
  Int_t GetNSel(){return fNSel;}
  Int_t GetNIdentical(){return fNIdentical;}
  
  vector<THSParticle*>* GetAllParticles(){return fAllParticles;}
  // vector<THSParticle*>* GetAllParticles(){return &fAllParticles;}
  void SetParticles(vector<THSParticle*> *parts){fAllParticles=parts;}
  void SetSelParticles(vector<THSParticle*> parts){fSelected=parts;}
  
  
  void SetSelIter(THSParticleIter* com){fSelIter=com;}
  void SetRemIter(THSParticleIter* com){fRemIter=com;}

  void AddEventParticle(THSParticle* part){fEvParts.push_back(part);};
  void AddEventParticles(THSParticle* part0=nullptr,THSParticle* part1=nullptr,THSParticle* part2=nullptr,THSParticle* part3=nullptr,THSParticle* part4=nullptr,THSParticle* part5=nullptr,THSParticle* part6=nullptr,THSParticle* part7=nullptr,THSParticle* part8=nullptr,THSParticle* part9=nullptr);
  void SetEventParticles(vector<THSParticle*> parts){fEvParts=parts;};
  
  void SortEvent();
  void SetCombi(THSCombitorial combi){fCombi=combi;fUseCombi=kTRUE;}
  
  
  THSParticle* NextParticle();
  Bool_t NextCombitorial();
  
  void ConfigureIters();
  void SetNextInnerIter(THSParticleIter *_iter);
  void SelectXofY(Int_t _X, Int_t _Y, THSCombitorial* sel_iter=nullptr);
  
  Bool_t Testing();

};

#endif
