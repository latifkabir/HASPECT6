
#ifndef THSEP_OMEGAP_h
#define THSEP_OMEGAP_h

#include "THSFinalState.h"
#include "TLorentzVector.h"
#include <vector>

class THSep_omegap : public THSFinalState{

 public :
  THSep_omegap();
  ~THSep_omegap(){};

  virtual Bool_t WorkOnEvent();

  //Init functions
  void Init_Iter0();
  void Init_Iter1();
  void Topo_0();
  void Topo_1();
  void Init_Generated();
  //void Init_IterX();
  //void Topo_X();

  void Kinematics();
  protected :
  //Topology flags for this project
  Int_t fTID_0=-1;
  Int_t fTID_1=-1;
  // Int_t fTID_=-1;
 
  //Initial state
  TLorentzVector fBeam=TLorentzVector(0,0,10.6,10.6);
  TLorentzVector fTarget=TLorentzVector(0,0,0,0.938272);
  //Final Particles
  THSParticle fElectron=THSParticle("e-");
  THSParticle fProton=THSParticle("proton");
  THSParticle fPip=THSParticle("pi+");
  THSParticle fPim=THSParticle("pi-");
  THSParticle fPi0=THSParticle("pi0");
  THSParticle fGamma1=THSParticle("gamma");
  THSParticle fGamma2=THSParticle("gamma");
  THSParticle fOmega=THSParticle("omega");
 
  //Observables
  Double_t f_t=0;
  Double_t f_Q2=0;
  Double_t f_W=0;
  Double_t f_Pol=0;
  Double_t fCMCosTh=0;
  Double_t fCMPhi=0;
  Double_t fHelCosTh=0;
  Double_t fHelPhi=0;

  //Discriminators
  Double_t fMissMass=0;
  Double_t fMissMassP=0;
  Double_t fOMass=0;
  Double_t fPi0Mass=0;

   public :
  virtual void FinalStateOutTree(TTree* tree);

 

};

#endif //ifdef THSep_omegap
