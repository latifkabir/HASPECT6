/**
	\class THSFinalState
	
	Class to control particlular analysis
	Users need to create their own specific final state class
	It should operate on 3 classes of data:
	real, mc rec, mc gen
	input should be in form of events extracted
	from trees or text files,...
	
	The main task is to take detected particles and
	convert them into information needed for subsequent
	fitting of observables, e.g TLorentzVectors for AmpTools
	Double_t s for RooFit
	
	@example THSFinalTemp.C
*/



#include "THSFinalState.h"

THSFinalState::THSFinalState(){
  
  fMapPDGtoParticle[2212]=&fVecProtons;
  fMapPDGtoParticle[211]=&fVecPiPs;
  fMapPDGtoParticle[-211]=&fVecPiMs;
  fMapPDGtoParticle[111]=&fVecPi0s;
  fMapPDGtoParticle[321]=&fVecKPs;
  fMapPDGtoParticle[-321]=&fVecKMs;
  fMapPDGtoParticle[11]=&fVecEls;
  fMapPDGtoParticle[-11]=&fVecPos;
  fMapPDGtoParticle[22]=&fVecGams;
  fMapPDGtoParticle[-22]=&fVecBeams;
  fMapPDGtoParticle[1E6]=&fVecPlus;
  fMapPDGtoParticle[-1E6]=&fVecMinus;
  fMapPDGtoParticle[0]=&fVec0;


}
//////////////////////////////////////////////////////////////////////////
///topo should be a pdg particle list seperated by : e.g. "proton:pi+:pi-"
Int_t THSFinalState::AddTopology(TString topo){
  vector<Int_t> detpart;
  TObjArray* OptList = topo.Tokenize(":");
  for(Int_t i=0;i<OptList->GetEntries();i++){
    Int_t pdg=0;
    Int_t charge=0;
    if(((TObjString*)OptList->At(i))->String()==TString("Rootino+")) {pdg=1E6;charge=1;}
    else if(((TObjString*)OptList->At(i))->String()==TString("Rootino-")) {pdg=-1E6;charge=-1;}
    else if(((TObjString*)OptList->At(i))->String()==TString("Beam")) {pdg=-22;charge=-1;}
    else if(!TDatabasePDG::Instance()->GetParticle(((TObjString*)OptList->At(i))->String()))Error("THSFinalState::AddTopology","Particle not found = %s",((TObjString*)OptList->At(i))->String().Data());
    else{
      pdg=TDatabasePDG::Instance()->GetParticle(((TObjString*)OptList->At(i))->String())->PdgCode();
      charge=TDatabasePDG::Instance()->GetParticle(pdg)->Charge();
      if(charge!=0)charge=charge/TMath::Abs(charge); //just get sign not mag.
    }
    if(fUseChargePID) detpart.push_back(charge);
    else  detpart.push_back(pdg);
    //Make list of all possible detected types to optimise particle vectors
    if(std::find(fDetTypes.begin(), fDetTypes.end(), pdg) == fDetTypes.end()) fDetTypes.push_back(pdg);
  }
  //sort into order for checking
  std::sort(detpart.begin(),detpart.begin()+detpart.size());
  cout<<"AddTopology ";
  for(UInt_t i=0;i<detpart.size();i++)
    cout<<detpart[i]<<" ";
  cout<<endl;
  fTopo.push_back(detpart);
  fNTopo++;
  return fNTopo-1;//return index of this topology
}

/////////////////////////////////////////////////////////////////////
///Note particles with TruthOnly are reserved for undetected particles
///which may be kept in simulated data, but shouldn't be included
///in final state check
Int_t THSFinalState::FindTopology(){
  if(IsPermutating()) return fCurrTopo;//Topology hasn't changed  with permutation
  vector<Int_t> thisTopo;
  for(UInt_t ip=0;ip<frDetParts->size();ip++){
    //    if((frDetParts->At(ip).TruthOnly())==kFALSE)thisTopo.push_back(frDetParts->At(ip).PDG());
    if(!fUseChargePID)thisTopo.push_back(frDetParts->at(ip)->PDG());
    else{
      thisTopo.push_back(frDetParts->at(ip)->Charge());
    }
  }
  std::sort(thisTopo.begin(),thisTopo.begin()+thisTopo.size());
  // cout<<"ThisTopology ";
  // for(UInt_t i=0;i<thisTopo.size();i++)
  //   cout<<thisTopo[i]<<" ";
  // cout<<endl;

  //now check to see if this topology matches any of our predefines ones
  fCurrTopo=-1;
  for(Int_t itopo=0;itopo<fNTopo;itopo++)
    if(thisTopo==fTopo[itopo]){
      fCurrTopo=itopo;
    }
  fNParts=thisTopo.size();
  thisTopo.clear();
  //Init particle vectors, this should only be done once for each event
  //Not for each permuation
  if(fCurrTopo>-1) InitParticles();
  
  return fCurrTopo;
}

//////////////////////////////////////////////////////////////////////
///Note particles with TruthOnly are reserved for undetected particles
///which may be kept in simulated data, but shouldn't be included
///in final state check
Int_t THSFinalState::FindInclusiveTopology(Int_t incType){
  if(IsPermutating()) return fCurrTopo;//Topology hasn't changed  with permutation
  vector<Int_t> thisTopo;
  for(UInt_t ip=0;ip<frDetParts->size();ip++){
    //    thisTopo.push_back(frDetParts->at(ip)->PDG());
    if(!fUseChargePID)thisTopo.push_back(frDetParts->at(ip)->PDG());
    else{
      thisTopo.push_back(frDetParts->at(ip)->Charge());
    }
  }
  std::sort(thisTopo.begin(),thisTopo.begin()+thisTopo.size());
  // cout<<"ThisTopology ";
  // for(UInt_t i=0;i<thisTopo.size();i++)
  //   cout<<thisTopo[i]<<" ";
  // cout<<endl;

  //now check to see if this topology matches any of our predefines ones
  fCurrTopo=-1;
  for(Int_t itopo=0;itopo<fNTopo;itopo++){
    vector<Int_t>ptypes;
    Bool_t isThis=kTRUE;//assume it is this topo
    for(UInt_t ipart=0;ipart<fTopo[itopo].size();ipart++){
      //check if ptypes already includes this type fTopo[itopo][ipart]
      if(std::find(ptypes.begin(), ptypes.end(), fTopo[itopo][ipart]) == ptypes.end()){
	ptypes.push_back(fTopo[itopo][ipart]);//Not there so add it
	//Now see if sufficient number in thisTopo
	Int_t thiscount = std::count (thisTopo.begin(), thisTopo.end(), fTopo[itopo][ipart]);
	Int_t topocount = std::count (fTopo[itopo].begin(), fTopo[itopo].end(), fTopo[itopo][ipart]);
	if(thiscount<topocount)//Not sufficient of this type, not this topo
	  isThis=kFALSE;
      }
    }
    ptypes.clear();
    if(isThis){
      //contains at least the particles we need
      if(incType!=-1E9){ //if requested only inclusive in 1 particle type check that is case
       
	UInt_t thiscount = std::count (thisTopo.begin(), thisTopo.end(), incType);
	UInt_t topocount = std::count (fTopo[itopo].begin(), fTopo[itopo].end(), incType); //how many incTypes in original topo
	if(thiscount-topocount!=thisTopo.size()-fTopo[itopo].size())
	  continue; //not just incType extra particles
      }
      fCurrTopo=itopo;
      break; //don't look for any more    
    }
  }
    
  fNParts=thisTopo.size();
  thisTopo.clear();
  //Init particle vectors, this should only be done once for each event
  //Not for each permuation
  if(fCurrTopo>-1) InitParticles();
  
  return fCurrTopo;
}

///////////////////////////////////////////////////////////////////
///Read in the particles from this input event and
///Get the first combination from the particle iterator
void THSFinalState::InitParticles(){
  //This should be optimised to run 1 loop not number of species loops
  if(fIsPermutating0==kTRUE) return;
  if(!fUseChargePID)
    //Init all particle type vectors
    InitDetParts();
  
  else//Just use the charge for PID (Need to SetUseChargePID())
    InitDetPartsCharge(); 
  
  //Get first combination and sort event particles
  if(!fDetIter[fCurrTopo]) {cout<<"ERROR THSFinalState::InitParticles() No PHSParticleITer defiend for this topology "<<fCurrTopo<<endl; exit(0);}

  //Initialise iterators for this event
  //fDetIter[fCurrTopo]->NextCombitorial();
  //fDetIter[fCurrTopo]->SortEvent();
  PermutateParticles();
}

///////////////////////////////////////////////////////////////////
///Read in the particles from this input event
///Place detected particles in array associated with their pdg type
///Loop over detected particles
void THSFinalState::InitDetParts(){
  //parts->clear();
  map< Int_t , vector< THSParticle* >* >::iterator itmap;
  for (itmap = fMapPDGtoParticle.begin(); itmap != fMapPDGtoParticle.end(); itmap++)
    {
      itmap->second->clear();
    }  
  for(UInt_t ip=0;ip<fNParts;ip++){
    THSParticle* part=frDetParts->at(ip);
    fMapPDGtoParticle[part->PDG()]->push_back(part);
   }
 
}
//////////////////////////////////////////////////////////////
///Read in the particles from this input event
///assign them in vectors corrsponing to their charge
void THSFinalState::InitDetPartsCharge(){

  fMapPDGtoParticle[Int_t(1E6)]->clear();
  fMapPDGtoParticle[Int_t(-1E6)]->clear();
  fMapPDGtoParticle[Int_t(0)]->clear();

  for(UInt_t ip=0;ip<fNParts;ip++){
    THSParticle* part=frDetParts->at(ip);
    fMapPDGtoParticle[Int_t(part->Charge()*1E6)]->push_back(part);
   }
 
}

//////////////////////////////////////////////////////////////////
///Process one event from the input tree
void THSFinalState::ProcessEvent(){
 //Process one input event
  InitEvent();
  if(frDetParts)fNDet=frDetParts->size();
  // Int_t firstuse=1;

  while(FSProcess());
  
  FinaliseEvent();
}
Bool_t THSFinalState::FSProcess(){
  
  if(!WorkOnEvent()) return kFALSE;
  
  if(IsGoodEvent()){
    UserProcess(); 
  }
  
  PermutateParticles();
  if(IsPermutating())
    return kTRUE;
  else
    return kFALSE;
}
void THSFinalState::UserProcess(){
  //fNUsedReadEvent+=firstuse;
  //if(firstuse) firstuse=0;
  if(fFinalTree) fFinalTree->Fill(); //fill for first combination
  fUID++;
  
}
// void THSFinalState::ProcessEvent(){
//  //Process one input event
//   InitEvent();
//   if(frDetParts)fNDet=frDetParts->size();
//   Int_t firstuse=1;
//   do{
//     if(!WorkOnEvent()) return;
//      if(fFinalTree)
//        if(IsGoodEvent()){
// 	 fNUsedReadEvent+=firstuse;
// 	 if(firstuse) firstuse=0;
// 	 fFinalTree->Fill(); //fill for first combination
// 	 fUID++;
//        }
//      PermutateParticles();
     
//   }
  
//   while(IsPermutating());
//   FinaliseEvent();
// }

//////////////////////////////////////////////////////////////
///returns true if another valid permuation to be tried \n
///returns false when time to move on \n
///Check through particle types,
///if more than 1 of a type permutate through all combinations
Bool_t THSFinalState::PermutateParticles(){
  if(!fTryPerm&&fNPerm>0) return kFALSE;
  if(fIsGenerated) return kFALSE;
  if(fCurrTopo<0) return kFALSE;
  fNPerm++;
  if(fIsPermutating1) return kTRUE; //so user can do their own permutation
  fIsPermutating0=kTRUE; //Will be set to false when event over

  if(fDetIter[fCurrTopo]->NextCombitorial()){fDetIter[fCurrTopo]->SortEvent();return kTRUE;}
  else fDetIter[fCurrTopo]->SortEvent();
  fIsPermutating0=kFALSE;
  return kFALSE; 
}

///////////////////////////////////////////////////////
///Loop through generated and record momentum distance
void THSFinalState::MatchWithGen(THSParticle *part){
  UInt_t match=0;
  Double_t mindist=1E10;
  for(UInt_t ip=0;ip<frGenParts->size();ip++){
    THSParticle* gen=frGenParts->at(ip);
    Double_t dist=part->p3Distance(gen->P4p()->Vect());
    if(dist<mindist){
      mindist=dist;
      match=ip;
    }
  }
  part->SetTruth(frGenParts->at(match)->P4(),frGenParts->at(match)->Vertex(),frGenParts->at(match)->PDG());
}

//////////////////////////////////////////////////////////////////////
///Loop through generated and record momentum distance 
///if the nearest generated track is the same as the TRUTH return true
Bool_t THSFinalState::IsCorrectTruth(THSParticle *part){
  if(!frGenParts) return kFALSE;
  if(!frGenParts->size())return kFALSE;
  
  UInt_t match=0;
  Double_t mindist=1E10;
  for(UInt_t ip=0;ip<frGenParts->size();ip++){
    THSParticle* gen=frGenParts->at(ip);
    Double_t dist=part->p3Distance(gen->P4p()->Vect());
    if(dist<mindist){
      mindist=dist;
      match=ip;
    }
  }
  if(part->TruthP4()==frGenParts->at(match)->P4()&&part->TruthPDG()==frGenParts->at(match)->PDG())
    return kTRUE;
  return kFALSE;
}

///////////////////////////////////////////////
///Check if our final vectors match with truth
void THSFinalState::CheckTruth(){
  if(fIsGenerated) return;
  //Already got one for this event
  //If we set this correct it would be double counting
  //This can happen when using inclusive topologies
  //And rotate particles that are not in the defined topology
  fCorrect=0;
  if(fGotCorrectOne) return;
  //Check if our final vectors match with truth
  fCorrect=1;
  for(UInt_t ip=0;ip<fFinal.size();ip++)
    fCorrect*=IsCorrectTruth(fFinal[ip]);
  
  if(fCorrect) fGotCorrectOne=kTRUE;
}

//////////////////////////////////////////////////
///Interface to THSParticleIter
///Add selected iterator to previously selected particles
THSParticleIter* THSFinalState::AddSelectToSelected(THSParticleIter* diter,Int_t ni,Int_t nsel,THSParticle* part0,THSParticle* part1,THSParticle* part2,THSParticle* part3,THSParticle* part4,THSParticle* part5,THSParticle* part6,THSParticle* part7,THSParticle* part8,THSParticle* part9){
  if(nsel*ni>diter->GetNSel()) cout<<"WARNING THSFinalState::AddSelectToSelected : trying to select more particles than exist in original iterator "<<nsel*ni<<" "<< diter->GetNSel()<<endl;
  THSParticleIter *new_iter=AddSelectXofY(ni,nsel,part0,part1,part2,part3,part4,part5,part6,part7,part8,part9);
  diter->SetSelIter(new_iter);
  return new_iter;
}
//////////////////////////////////////////////////
///Interface to THSParticleIter
///Add selected iterator to remaining particles
THSParticleIter* THSFinalState::AddSelectToRemainder(THSParticleIter* diter,Int_t ni,Int_t nsel,THSParticle* part0,THSParticle* part1,THSParticle* part2,THSParticle* part3,THSParticle* part4,THSParticle* part5,THSParticle* part6,THSParticle* part7,THSParticle* part8,THSParticle* part9){
  THSParticleIter *new_iter=AddSelectXofY(ni,nsel,part0,part1,part2,part3,part4,part5,part6,part7,part8,part9);
  diter->SetRemIter(new_iter);
  return new_iter;
}
//////////////////////////////////////////////////
///Interface to THSParticleIter
///Add selected iterator
THSParticleIter* THSFinalState::AddSelectXofY(Int_t ni,Int_t nsel,THSParticle* part0,THSParticle* part1,THSParticle* part2,THSParticle* part3,THSParticle* part4,THSParticle* part5,THSParticle* part6,THSParticle* part7,THSParticle* part8,THSParticle* part9){
  
  THSParticleIter *new_iter= new THSParticleIter();
  new_iter->AddEventParticles(part0,part1,part2,part3,part4,part5,part6,part7,part8,part9);
  new_iter->SelectXofY(ni,nsel);
  return new_iter;
}
//////////////////////////////////////////////////
///Interface to THSParticleIter
///Add permutate iterator to previously selected particles
THSParticleIter* THSFinalState::AddPermutateToSelected(THSParticleIter* diter,THSParticle* part0,THSParticle* part1,THSParticle* part2,THSParticle* part3,THSParticle* part4,THSParticle* part5,THSParticle* part6,THSParticle* part7,THSParticle* part8,THSParticle* part9){ 
  THSParticleIter *new_iter=AddPermutate(part0,part1,part2,part3,part4,part5,part6,part7,part8,part9);
  diter->SetSelIter(new_iter);
  return new_iter;
}
/////////////////////////////////////////////////////////////////
///Interface to THSParticleIter
///Add permutate iterator to remaining particles
THSParticleIter* THSFinalState::AddPermutateToRemainder(THSParticleIter* diter,THSParticle* part0,THSParticle* part1,THSParticle* part2,THSParticle* part3,THSParticle* part4,THSParticle* part5,THSParticle* part6,THSParticle* part7,THSParticle* part8,THSParticle* part9){
  THSParticleIter *new_iter=AddPermutate(part0,part1,part2,part3,part4,part5,part6,part7,part8,part9);
  diter->SetRemIter(new_iter);
  return new_iter;
}
/////////////////////////////////////////////////////////////////
///Interface to THSParticleIter
///Add selected iterator 
THSParticleIter* THSFinalState::AddPermutate(THSParticle* part0,THSParticle* part1,THSParticle* part2,THSParticle* part3,THSParticle* part4,THSParticle* part5,THSParticle* part6,THSParticle* part7,THSParticle* part8,THSParticle* part9){
  
  THSParticleIter *new_iter= new THSParticleIter();
  new_iter->SetCombi(THSPermutation());
  new_iter->AddEventParticles(part0,part1,part2,part3,part4,part5,part6,part7,part8,part9);
  return new_iter;
}

////////////////////////////////////////////////////////////////
///Interface to THSParticleIter
///Create particle iterator for a particticular species
///Given by the particle vector parts (e.g. fVecMinus)
THSParticleIter* THSFinalState::CreateParticleIter(Int_t topo,vector<THSParticle*> *parts,Int_t Nsel){

  THSParticleIter *diter= new THSParticleIter();
  diter->SetCombi(THSSelection());
  diter->SetParticles(parts);
  diter->SetNSel(Nsel);
   if(!fDetIter[topo])
    fDetIter[topo]=diter;
  else
    fDetIter[topo]->SetNextInnerIter(diter);
    
  return diter;
}
