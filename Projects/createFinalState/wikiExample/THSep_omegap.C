/**
 * \class THSep_omegap
 * 
 * Template final class .
 * 
 * Users should create their own analysis specific project classes.
 * 
 */



#include "TDatabasePDG.h"
#include "THSep_omegap.h"
#include <algorithm>


THSep_omegap::THSep_omegap(){
  SetUseChargePID();//If want to only use charge for IDing
  //include topology for analysis and get index
  fTID_0=AddTopology("pi+:pi-:gamma:gamma:proton:e-");
  fTID_1=AddTopology("pi+:pi-:proton:e-");
 
  //Set final state
  //Note for CheckTruth algorithm to work all particles added 
  //to fFinal should have been in Generated 
  fFinal.push_back(&fElectron);
  fFinal.push_back(&fProton);
  fFinal.push_back(&fPip);
  fFinal.push_back(&fPim);
  fFinal.push_back(&fPi0);

  //Initialise particle iters
  fDetIter.resize(fNTopo);
  Init_Iter0();
  Init_Iter1();
}


void THSep_omegap::Init_Generated(){
 if(!THSFinalState::frGenParts) return;
  if(THSFinalState::frGenParts->size()!=5) {fGoodEvent=kFALSE;return;}
  //Fill our data member particles
  //User is responsible for indicing right
  //comes from order in generated file (e.g LUND)
  if(THSFinalState::fIsGenerated){
    fPip=*frGenParts->at(0);
    fPim=*frGenParts->at(1);
    fPi0=*frGenParts->at(2);
    fProton=*frGenParts->at(3);
    fElectron=*frGenParts->at(4);
  }
  else{//Just assign truth values
    //fElectron.SetTruth(frGenParts->at(0));
    fPip.SetTruth(frGenParts->at(0));
    fPim.SetTruth(frGenParts->at(1));
    fPi0.SetTruth(frGenParts->at(2));
    fProton.SetTruth(frGenParts->at(3));
    fElectron.SetTruth(frGenParts->at(4));
 
  } 
}

//Define topology Iterator functions
void THSep_omegap::Init_Iter0(){
  //THSParticle iterator initialisation
  ////////////-ve particles
  THSParticleIter* diter=CreateParticleIter(fTID_0,&fVecMinus,2);
  AddPermutateToSelected(diter,&fElectron,&fPim);
  ////////////+ve particles
  diter=CreateParticleIter(fTID_0,&fVecPlus,2);
  AddPermutateToSelected(diter,&fProton,&fPip);
  ///////////neutral particles
  diter=CreateParticleIter(fTID_0,&fVec0,2);
  AddSelectToSelected(diter,1,2,&fGamma1,&fGamma2);
  
  fDetIter[fTID_0]->ConfigureIters();

}
void THSep_omegap::Init_Iter1(){
  //THSParticle iterator initialisation
  ////////////-ve particles
  THSParticleIter* diter=CreateParticleIter(fTID_1,&fVecMinus,2);
  AddPermutateToSelected(diter,&fElectron,&fPim);
  ////////////+ve particles
  diter=CreateParticleIter(fTID_1,&fVecPlus,2);
  AddPermutateToSelected(diter,&fProton,&fPip);
 
  fDetIter[fTID_1]->ConfigureIters();
 }
// void THSep_omegap::Init_IterX(){
//  THSParticleIter* diter=CreateParticleIter(fTID_X,&fVecY,NPARTICLES);
//  ...
//  fDetIter[fTID_0]->ConfigureIters();
// }
//Or if iterator the same as a previous one
// void THSep_omegap::Init_IterY(){
//   fDetIter[fTID_1]=fDetIter[fTID_0];
// }
//Define topology functions
void THSep_omegap::Topo_0(){
  //Reconstruct missing or combined particles
  //I can put a cut on the electron, for example on the Detector
  if(fElectron.Detector()>0) {fGoodEvent=kFALSE;return;}
  //Reconstruct missing or combined particles
  TLorentzVector miss=fBeam+fTarget-fElectron.P4()-fProton.P4()-fPip.P4()-fPim.P4()-fGamma1.P4()-fGamma2.P4();
  
  fPi0.Add(&fGamma1,&fGamma2,fPi0.PDG());
  
  fMissMass=miss.M();
  fPi0Mass=fPi0.MassDiff();
  //I could place a cut on the 2gamma invariant mass if I wanted
  //  if(TMath::Abs(fPi0Mass)>0.02){fGoodEvent=kFALSE;return;}
  
}
void THSep_omegap::Topo_1(){
  //Reconstruct missing or combined particles
  if(fElectron.Detector()>0) {fGoodEvent=kFALSE;return;}
  //Reconstruct missing or combined particles
  TLorentzVector miss=fBeam+fTarget-fElectron.P4()-fProton.P4()-fPip.P4()-fPim.P4();
  fPi0.SetP4(miss);
  fPi0.SetMeasMass(fPi0.P4().M());
  fPi0.TakePDGMass();
  
  fMissMass=miss.M();
  fPi0Mass=fPi0.MassDiff();
  
}
// void THSep_omegap::Topo_X(){
// }

void THSep_omegap::Kinematics(){
   if(!THSFinalState::fGoodEvent) return;//don't do calculations
   //Do calculations if Good Event
   //missing mass off proton
   TLorentzVector missP=fBeam+fTarget-fElectron.P4()-fProton.P4();
   fMissMassP=missP.M();
   //invariant mass of 3pi
   fOmega.SetP4(fPip.P4()+fPim.P4()+fPi0.P4());
   fOmega.SetMeasMass(fOmega.P4p()->M());
   fOmega.TakePDGMass();
   fOMass=fOmega.MassDiff();
   
   //Set the initial and scattered electron and target 4-vectors
   fKine.SetElecsTarget(fBeam,fElectron.P4(),fTarget);
   //set the final state meson and baryon 4-vectors
   fKine.SetMesonBaryon(fOmega.P4(),fProton.P4());
   //Get electron scattering variables
   f_W=fKine.W();
   f_Q2=fKine.Q2();
   f_Pol=fKine.GammaPol();
   //Get Omega production variables
   fKine.ElectroCMDecay();
   fCMCosTh=fKine.CosTheta();
   fCMPhi=fKine.Phi();
   f_t=fKine.t(); //t for Omega meson production
   
}

void THSep_omegap::FinalStateOutTree(TTree* tree){
  THSFinalState::fFinalTree=tree;
  //tree->Branch("Final",&fFinal);
  tree->Branch("MissMass",&fMissMass,"MissMass/D");
  tree->Branch("MissMassP",&fMissMassP,"MissMassP/D");
  tree->Branch("Topo",&fCurrTopo,"Topo/I");
  tree->Branch("NPerm",&fNPerm,"NPerm/I");
  tree->Branch("NDet",&fNDet,"NDet/I");
  tree->Branch("Correct",&fCorrect,"Correct/I");
  tree->Branch("OMass",&fOMass,"OMass/D");
  tree->Branch("Pi0Mass",&fPi0Mass,"Pi0Mass/D");
  tree->Branch("t",&f_t,"t/D");
  tree->Branch("Q2",&f_Q2,"Q2/D");
  tree->Branch("W",&f_W,"W/D");
  tree->Branch("Pol",&f_Pol,"Pol/D");
  tree->Branch("CMCosTh",&fCMCosTh,"CMCosTh/D");
  tree->Branch("CMPhi",&fCMPhi,"CMPhi/D");
}
Bool_t THSep_omegap::WorkOnEvent(){
  //Should this event be saved?
  THSFinalState::fGoodEvent=kTRUE;
  THSFinalState:: fCorrect=0; //Correct permutation? used for simulation only

  //If generated MC events
  if(THSFinalState::fIsGenerated) Init_Generated();
  else{//Look for reconstructed events
    //if reconstructed Find the detected particles in this event
    //    if(FindInclusiveTopology()==-1) {fGoodEvent=kFALSE;return fIsPermutating0=kFALSE;}
    if(FindTopology()==-1) {fGoodEvent=kFALSE;return fIsPermutating0=kFALSE;}
    //Do they correspond to a defined topology?
    else if(fCurrTopo==fTID_0) Topo_0();
    else if(fCurrTopo==fTID_1) Topo_1();
    // else if(fCurrTopo==fTID_) Topo_();
    else fGoodEvent=kFALSE;
    //Get truth values
    Init_Generated();
  }
  
  //Calc kinematics
  Kinematics();
  
  //Check if assigned vectors agree with true generated
  //Simulation only
  THSFinalState::CheckTruth();
  
  //Can do some checks if event is worth writing or not
  //if()fGoodEvent=kTRUE;
  //else() fGoodEvent =kFALSE;
  if(fIsGenerated) return kTRUE; //Generated only 1 permutation
  return kTRUE;
}
