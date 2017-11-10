/**
	\class THSParticleIter
	
	Class to loop over THSParticle vectors
	Permuation, Rotation or Combination (Selection) done
	via THSCombinitorial class
	This class contains a pointer to another THSParticleIter
	to allow for recursive combinations.
	This class is used to control the permutating of particles
	in THSFinalState
*/

#include "THSParticleIter.h"
//////////////////////////////////////////////////////////////////////////
///Add function comments like this

Bool_t THSParticleIter::Testing(){
  THSParticle ppp;
  cout<<"PPP "<<ppp.Detector()<<endl;;
  THSCombitorial combii;
  cout<<combii.GetType()<<endl;
  return kTRUE;
}
//////////////////////////////////////////////////////////////////////////
///Get the next combination of the iterator
///Returns kFALSE when all possible combinations have been given
///Returns kTRUE when still more combinations to go
///The new ordering can be accessed via NextParticle()
///Or if event particles have been set they can be updated
///by a call to SortEvent();

Bool_t THSParticleIter::NextCombitorial(){
  // cout<<"THSParticleIter::NextCombitorial() "<<fInnerIter<<" "<<fSelIter<<" "<<fAllParticles->size()<<" "<<fEvParts.size()<<endl;
  //  if(!fUseCombi) return kFALSE;
  //Make sure Iterator fully configured
  if((fSelIter||fRemIter)&&(!fInnerIter))ConfigureIters();
  
  //Look for inner iterators after we have made the first
  //combitorial of this one
  if(!fDoneSelRem){ //First combination
    if(fInnerIter) //recursively go through all linked iterators
      if(fInnerIter->NextCombitorial())
	return kTRUE;
    //no further inner combinations to go
    fDoneSelRem=kTRUE;
  }
  //else get the next combitorial for this iterator
  //Check if all combinatioon already tried
  if(fCombi.IsFinished()){return kFALSE;}
  //If more to go get the next one
  //  fSelected=fCombi.Next(&fAllParticles,fNSel,fNIdentical);
  vector<THSParticle*>  vec_combi=*fAllParticles;
  if(fCombi.GetType()==2){ //for selections with more requested than in allparticles, add extra particles (with zero energy etc)
    while(vec_combi.size()<UInt_t(fNSel)){
      vec_combi.push_back(fExtraParticle);
      //  cout<<"ADD EXTRA "<<endl;
    }
    // cout<<"back "<<vec_combi.size()<<" "<<vec_combi.back()<<endl;
    fSelected=fCombi.Next(&vec_combi,fNSel,fNIdentical);
  }
  //cout<<fSelected.size()<<endl;
  else fSelected=fCombi.Next(fAllParticles,fNSel,fNIdentical);
  //update selected and remaining particles
  if(fRemIter){
    if(fCombi.GetType()==2) fRemainder=fCombi.Remainder(&vec_combi);
    else fRemainder=fCombi.Remainder(fAllParticles);
    fRemIter->SetParticles(&fRemainder);
  }
  if(fSelIter){
    fSelIter->SetParticles(&fSelected);
  }
  fDoneSelRem=kFALSE;
  //get the first combination of inner iterators
  if(fInnerIter) fInnerIter->NextCombitorial();
  
  return kTRUE;  
  
}
//////////////////////////////////////////////////////////////////////////
///Set the order that the selected and remaining iterators will be called
///First Add selected then Remaining recursively
void THSParticleIter::ConfigureIters(){
  if(fSelIter){
    SetNextInnerIter(fSelIter);
    if(fRemIter){
      SetNextInnerIter(fRemIter);
    }
  }
  else if(fRemIter){	
    SetNextInnerIter(fRemIter);
  }
  if(fInnerIter) fInnerIter->ConfigureIters();
}
//////////////////////////////////////////////////////////////////////////
///Add _iter to the end of the current iterator sequence
void THSParticleIter::SetNextInnerIter(THSParticleIter *_iter){
  if(fInnerIter) //If inner already exists move on inside it to look for next space
    fInnerIter->SetNextInnerIter(_iter);
  else
    fInnerIter=_iter;
}
////////////////////////////////////////////////////////////////////////
///In case all particle vector changes size need to redo selection map
// ResetSelection(){
//   if(fInnerIter)
//     fInnerIter->ResetSelection();
//   else
//     if(fCombi.GetType()==2) fCombi.ResetCombi();
// }

//////////////////////////////////////////////////////////////////////////
///Get the next particle in the current combination
///Returns a pointer to the next THSParticle
///If all particles done return null

THSParticle* THSParticleIter::NextParticle(){
  
  THSParticle* particle=nullptr;
    //If we have a THSCombination get the particle from there
  if(fUseCombi){
    if(fSelIter){
      //If we have a Selected THSCombination get the particle from there
      if(fParti<fSelected.size()){
	if((particle=fSelIter->NextParticle())){
	  return particle;
	}
	else fParti=fSelected.size(); //Don't look in SelIter again until gone though remainder
      }
    }
    //else //no combitorial for selected events just return selected
    if(!fSelIter&&fParti<fSelected.size())
      return fSelected[fParti++];//Just use results of this combitorial
    
    if(fRemIter){//if we want to use the remaining particles
      //Get the next particle from there
      if((particle=fRemIter->NextParticle())){
	  return particle;
      }
      //If no remaining iterator given, throw away remaining particles
    }
    fParti=0;
    return nullptr; //end of vector
    
  }
  
  else{//No combinitorial behaviour just iterate through vector
    // if(fParti<fAllParticles.size()){
    // 	return fAllParticles[fParti++];
    // }
    if(fParti<fAllParticles->size()){
      return fAllParticles->at(fParti++);
    }
    else{
      fParti=0;
      return nullptr;     
    }
  }
  return nullptr;
}

//////////////////////////////////////////////////////////////////////////
///Automates the recursive selection of Y identical particle types  X times
///e.g. select 3pi0 from 6 photons
///== select 3 from 6; select 2 from remaining 4; use remaining 2
///example usage : piterator.SelectXofY(3,2);
///An additional level of combinations sel_iter can be given and will be
///performed on all Y selections
/// i.e.can add a permutation to the X selections
///example usage : piterator.SelectXofY(3,2, new THSPermutation());

void THSParticleIter::SelectXofY(Int_t _X, Int_t _Y, THSCombitorial* sel_iter){
  cout<<"SelectXofY "<<_X<<" "<<_Y<<endl;
  //recursively create iterators
  if(fEvParts.size()!=UInt_t(_X*_Y)) cout<<"Warning ParticleIter::SelectXofY not sufficeint Event Particles Set "<<fEvParts.size()<<" when we need "<<_X*_Y<<endl;
  SetCombi(THSSelection());
  SetNSel(_Y);
  SetNIdentical(_X);
  
  //Now sort the event particles
  //Take the first _Y for those selected here
  //if(fEvParts.size()<_Y) cout<<"Warning ParticleIter::SelectXofY insufficient particles in Event Particles to cover "<<endl;
  vector<THSParticle*>  selParts=fCombi.SubSet(&fEvParts,_Y);
  vector<THSParticle*>  remParts=fCombi.SubSet(&fEvParts,fEvParts.size()-_Y,_Y);
  
  fEvParts=selParts;
  
  if(sel_iter){
    fSelIter=new THSParticleIter();
    fSelIter->SetEventParticles(selParts);
    fSelIter->SetCombi(*sel_iter);
  }
  if(_X>1){
    fRemIter=new THSParticleIter();
    fRemIter->SetEventParticles(remParts);
    fRemIter->SelectXofY(_X-1,_Y,sel_iter);
  }
}
//////////////////////////////////////////////////////////////////////////
///If event particles have been given assign copy a HSParticle
///to them in the current combination order
///The particles have to have been assigned to the correct iterator
void THSParticleIter::SortEvent(){
  // cout<<"sort "<<fAllParticles<<" "<<fUseCombi<<" "<<fAllParticles->size()<<" "<<fEvParts.size()<<" "<<fSelected.size()<<" "<<fCombi.GetType()<<endl;
  if(fUseCombi){
    // if(fSelIter){
    //   fSelIter->SortEvent();
    // }
    if(fInnerIter){
      fInnerIter->SortEvent();
    }
    //else{
      for(UInt_t isel=0;isel<fEvParts.size();isel++){
    	//if(isel>=fSelected.size())cout<<" THSParticleIter::SortEvent() Warning too many particles "<<isel<<" "<<fEvParts.size()<<endl;
    	if(isel<fSelected.size()){
	  fEvParts[isel]->CopyParticle(fSelected[isel],kFALSE);
	  fEvParts[isel]->TakePDGMass();
	}
	else{
	  fEvParts[isel]->P4p()->SetE(0);
	  fEvParts[isel]->P4p()->SetRho(0);
	  fEvParts[isel]->P4p()->SetTheta(0);
	}
      }
    // }
    // if(fRemIter){
    //   fRemIter->SortEvent();
    // }
    
    }
 
    else{ //No combitorial just take all particles as given
    // cout<<"Sorting "<<fEvParts.size()<<" "<<fAllParticles->size()<<endl;
   for(UInt_t iall=0;iall<fAllParticles->size();iall++){
     if(iall>=fAllParticles->size())cout<<" THSParticleIter::SortEvent() Warning too may particles in AllParticles "<<iall<<" "<<fEvParts.size()<<endl;
      // fEvParts[iall]->CopyParticle(fAllParticles[iall],kFALSE);
      // fEvParts[iall]->TakePDGMass();
     if(iall>=fEvParts.size())cout<<" THSParticleIter::SortEvent() Warning too may particles in EvParticles "<<iall<<" "<<fEvParts.size()<<endl;
     //  cout<<iall<<" "<<fEvParts[iall]<<" "<<fAllParticles->at(iall)<<endl;
      fEvParts[iall]->CopyParticle(fAllParticles->at(iall),kFALSE);
      fEvParts[iall]->TakePDGMass();
    }
  }
}
////////////////////////////////////////////////////////////////////
///Set up to 10 event particle with 1 function call
///If you need more just call it again
void THSParticleIter::AddEventParticles(THSParticle* part0,THSParticle* part1,THSParticle* part2,THSParticle* part3,THSParticle* part4,THSParticle* part5,THSParticle* part6,THSParticle* part7,THSParticle* part8,THSParticle* part9){
  if(part0) AddEventParticle(part0);
  if(part1) AddEventParticle(part1);
  if(part2) AddEventParticle(part2);
  if(part3) AddEventParticle(part3);
  if(part4) AddEventParticle(part4);
  if(part5) AddEventParticle(part5);
  if(part6) AddEventParticle(part6);
  if(part7) AddEventParticle(part7);
  if(part8) AddEventParticle(part8);
  if(part9) AddEventParticle(part9);

}
