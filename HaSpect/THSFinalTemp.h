
#ifndef THSPROJTEMP_h
#define THSPROJTEMP_h

#include "THSFinalState.h"
#include "TLorentzVector.h"
#include <vector>

class THSFinalTemp : public THSFinalState{

 public :
  THSFinalTemp();
  ~THSFinalTemp(){};

  virtual Bool_t WorkOnEvent();

  //Init functions
  void Init_Generated();
  //void Init_IterX();
  //void Topo_X();

  void Kinematics();
  protected :
  //Topology flags for this project
  // Int_t fTID_=-1;
 
  //Initial state
  TLorentzVector fBeam=TLorentzVector(0,0,10.6,10.6);
  TLorentzVector fTarget=TLorentzVector(0,0,0,0.938272);
  //Final Particles
 
  //Observables
  Double_t f_t;
  Double_t f_Q2;
  Double_t f_W;
  Double_t f_Pol;
  Double_t fCMCosTh;
  Double_t fCMPhi;

  //Discriminators
  Double_t fMissMass=0;

   public :
  virtual void FinalStateOutTree(TTree* tree);

 

};

#endif //ifdef THSFinalTemp
