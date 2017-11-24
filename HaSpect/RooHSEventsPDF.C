#include "RooHSEventsPDF.h"
#include "Riostream.h" 
#include "RooRealVar.h" 
#include "RooCategory.h" 
#include "RooRandom.h" 
#include "RooDataSet.h" 
#include "RooMsgService.h" 
#include "RooGaussian.h" 
#include "RooPlot.h" 
#include "RooHist.h" 
#include "RooPullVar.h" 
#include "TMath.h" 
#include "TObjArray.h" 
#include "TCanvas.h" 
#include "TBranch.h" 
#include <algorithm> 

ClassImp(RooHSEventsPDF)



RooHSEventsPDF::RooHSEventsPDF(const RooHSEventsPDF& other, const char* name) :  RooAbsPdf(other,name)
  {
    fIsClone=kTRUE;
    fParent=const_cast<RooHSEventsPDF*>(&other);
    if(other.fEvTree)fEvTree=other.fEvTree->CopyTree("");
    //    if(other.fEvTree)fEvTree=other.fEvTree->CloneTree();
    if(other.fInWeights) fInWeights=other.fInWeights; //probably need to clone this
    fNInt=other.fNInt;
    fGeni=other.fGeni;
    if(other.fEntryList)fEntryList=(TEntryList*)other.fEntryList->Clone();
    fForceConstInt=other.fForceConstInt;
    fConstInt=other.fConstInt;
    fCheckInt=other.fCheckInt;
    fUseWeightsGen=other.fUseWeightsGen;
    fCut=other.fCut;
  }
RooHSEventsPDF::~RooHSEventsPDF(){
  //RooFit clones everything so I need to give the original
    //object the entrylist if I want to use it!
    if(fIsClone&&fParent&&fEntryList){
      if(fMaxValue){//has this clone used generator?
	fParent->SetEntryList(fEntryList);
	fParent->SetGeni(fGeni);
      }
    }
    if(fEntryList) delete fEntryList;
    if(fLast) delete fLast;
    if(fEvTree) delete fEvTree;
    
    if(fWeights){
      fWeights->Save();
      delete fWeights;
    }
    for(UInt_t i=0;i<fVarSet.size();i++)
      delete fVarSet[i];
    
    fVarSet.clear();
  }

void RooHSEventsPDF::InitSets(){
  fNpars=fParSet.size();
  fNvars=fProxSet.size();
  fNcats=fCatSet.size();
  fLast=new Double_t[fNpars+1]; //Number of fit parameters
  for(Int_t i=0;i<fNpars+1;i++)
    fLast[i]=100;
}
RooArgSet RooHSEventsPDF::VarSet(Int_t iset) const{
  RooArgSet aset(Form("VarSet_%d",iset));
  if(iset==0){
    for(UInt_t i=0;i<fProxSet.size();i++){
      aset.add(fProxSet.at(i)->arg());
    }
    for(UInt_t i=0;i<fCatSet.size();i++){
      aset.add(fCatSet.at(i)->arg());
    }
  }
  else{
    for(UInt_t j=0;j<fProxSet.size();j++){//add if not proxy being removed
       if(fProxSet[iset-1]->GetName()!=fProxSet[j]->GetName()) aset.add(fProxSet.at(j)->arg());
    }
    for(UInt_t i=0;i<fCatSet.size();i++){ //Just add all categories
      aset.add(fCatSet.at(i)->arg());
    }
  }
  return aset;
}
  Int_t RooHSEventsPDF::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const
{
  Info("RooHSEventsPDF::getGenerator","Looking for generator");
  if(!fEvTree) return 0; //no MC events to generate from
  //case generate all variables
  if (matchArgs(directVars,generateVars,VarSet(0))) return 1 ;
  return 0;

}
void RooHSEventsPDF::initIntegrator()
{
   //First check if gen branches are available and use them
  fEvTree->ResetBranchAddresses();
  for(UInt_t i=0;i<fProxSet.size();i++){
    if(fEvTree->GetBranch(fProxSet[i]->GetName())){
      fEvTree->SetBranchStatus(fProxSet[i]->GetName(),1);
      fEvTree->SetBranchAddress(fProxSet[i]->GetName(),&fMCVar[i]);
     }
  }
  for(UInt_t i=0;i<fCatSet.size();i++){
    if(fEvTree->GetBranch(fCatSet[i]->GetName())){
      fEvTree->SetBranchStatus(fCatSet[i]->GetName(),1);
      fEvTree->SetBranchAddress(fCatSet[i]->GetName(),&fMCCat[i]);
     }
  }
}
void RooHSEventsPDF::initGenerator(Int_t code)
{
  Info("RooHSEventsPDF::initGenerator","Going to generate starting from %lld with weights %d",fGeni,(Int_t)fUseWeightsGen);
  //First check if gen branches are available and use them
  // if(fEvTree->GetBranch(fgenStr+fProxSet[0]->GetName()))fEvTree->ResetBranchAddresses();
  //fEvTree->ResetBranchAddresses();
  for(UInt_t i=0;i<fProxSet.size();i++){
    //Look for true branches
    if(fEvTree->GetBranch(fgenStr+fProxSet[i]->GetName())){
      fEvTree->ResetBranchAddress(fEvTree->GetBranch(fProxSet[i]->GetName()));
      fEvTree->SetBranchStatus(fgenStr+fProxSet[i]->GetName(),1);
      fEvTree->SetBranchAddress(fgenStr+fProxSet[i]->GetName(),&fMCVar[i]);
    }
    else{
      Warning("RooHSEventsPDF::initGenerator"," No gen branches found will use standard for %s",fProxSet[i]->GetName());
    }
  }
  for(UInt_t i=0;i<fCatSet.size();i++){
    //Look for true branches
    if(fEvTree->GetBranch(fgenStr+fCatSet[i]->GetName())){
      fEvTree->ResetBranchAddress(fEvTree->GetBranch(fCatSet[i]->GetName()));
      fEvTree->SetBranchStatus(fgenStr+fCatSet[i]->GetName(),1);
      fEvTree->SetBranchAddress(fgenStr+fCatSet[i]->GetName(),&fMCCat[i]);
    }
    else{
      Warning("RooHSEventsPDF::initGenerator"," No gen branches found will use standard for %s",fCatSet[i]->GetName());
    }
  }
  Double_t value=0;
  if(code==1){
    //Brute force find maximum value
    //for(Int_t i=0;i<fEvTree->GetEntries();i++){
    for(Int_t i=1;i<fEvTree->GetEntries();i++){
      fEvTree->GetEntry(i);
        value=evaluateMC()*GetIntegralWeight();
      if(value>fMaxValue)fMaxValue=value;
    }
  }
  //construct entry list so can reproduce full tree branches,
  //not jist those loaded as variables
  if(!fEntryList){
    fEntryList=new TEntryList("GenEvents","GenEvents",fEvTree);
  }
  else{
    fEntryList->Reset();
    fEntryList->SetTree(fEvTree);
  }
  if(fUseWeightsGen){
    fWeights=new THSWeights("genWeights");
    fWeights->SetSpecies(GetName());
    fWeights->SetFile(TString(GetName())+"Weights.root");
  }
  Info("RooHSEventsPDF::initGenerator","Max value %lf",fMaxValue);
  cout<<fEvTree->GetEntries()<<endl;

 }
void RooHSEventsPDF::generateEvent(Int_t code){
  // Info("RooHSEventsPDF::generateEvent","Going to generate starting from %ld with ",fGeni);
  
  Double_t value=0;
  if(!fUseWeightsGen){
    while(fGeni<fEvTree->GetEntries()){
      fEvTree->GetEntry(fGeni++);
      value=evaluateMC()*GetIntegralWeight();
      if(value>fMaxValue*RooRandom::uniform()){//accept
	for(Int_t i=0;i<fNvars;i++)
	  (*(fProxSet[i]))=fMCVar[i];
	for(Int_t i=0;i<fNcats;i++)
	  (*(fCatSet[i]))=fMCCat[i];
	
	fEntryList->Enter(fGeni-1);
	return;
      }
    }
  }
  else{
    //using weights
    while(fGeni<fEvTree->GetEntries()){
      fEvTree->GetEntry(fGeni++);
      if(!CheckRange("")) continue;
      value=evaluateMC()*GetIntegralWeight();
      for(Int_t i=0;i<fNvars;i++)
	(*(fProxSet[i]))=fMCVar[i];
     for(Int_t i=0;i<fNcats;i++)
	(*(fCatSet[i]))=fMCCat[i];
      fWeights->FillWeight(fGeni-1,value); 
      fEntryList->Enter(fGeni-1);
      return;
    }
  }
  Error("RooHSEventsPDF::generateEvent","Ran out of events at %lld",fGeni);
  //Used up all the events in the tree!
}
Int_t RooHSEventsPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,const char* rangeName) const
{ 
  if(fForceNumInt) return 0; //might be good to check numerical integral sometimes
  if(!fEvTree&&!fForceConstInt) return 0; //no MC events to integrate over

  if(fProxSet.size()==1){//special case 1 variable
    if (matchArgs(allVars,analVars,VarSet(0))) return 1 ;
  }
  else{//For variables
     //    for(UInt_t i=0;i<1+fProxSet.size();i++){
    for(UInt_t i=0;i<1+fProxSet.size();i++){
      if(!fEvTree&&fForceConstInt&&i==0) return 1; //no tree, but const int
      else if(!fEvTree) return 0; //no const integral for projections
      if (matchArgs(allVars,analVars,VarSet(i))) return i+1 ;
    }
  }
  //Note not implemented for cats
    return 0;
}

Double_t RooHSEventsPDF::analyticalIntegral(Int_t code,const char* rangeName) const
{
  if(code==1&&fForceConstInt&&!fEvTree) {fLast[0]=1;return fLast[0];}
  //sort number of events first in case forced
  Long64_t NEv=0;
  
  if(fNInt>-1) NEv=fNInt;
  else NEv=fEvTree->GetEntries(); 
  if(NEv>fEvTree->GetEntries()) NEv=fEvTree->GetEntries();

  // Info("RooHSEventsPDF::analyticalIntegral","calcing my own integral");
  // return 1;
  //only recalculate if a par changes when all variables included (ie code=1)
  if(code==1)
    if(!CheckChange()) return fLast[0];
  // if(fNInt)fEvTree->SetEntries(fNInt);
  //In case changed for generation
  
  Double_t integral=0;
  if(code==1){
    for(Int_t ie=0;ie<NEv;ie++){
      fEvTree->GetEntry(ie);
      if(!CheckRange(rangeName)) continue;
      integral+=evaluateMC()*GetIntegralWeight();
    }
  }
  else {
    //inegrate over other variables for one variable fixed
    //index given by code -2 (defintion of code in getAnalyticalIntegral
    //This is used for plotting data and PDFS
    //if(dynamic_cast<RooProxyCategory*>fProxSet[code-2]) return 1;
    //cout<<fMCVar.size()<<" "<<fProxSet.size()<<" "<<*(fProxSet[code-2])<<endl;
    Double_t rmax=fProxSet[code-2]->max();
    Double_t rmin=fProxSet[code-2]->min();
    Double_t vval=*(fProxSet[code-2]);
    Double_t vrange=rmax-rmin;
    //  fProxSet[code-2]->Print();
    Int_t nbins=((RooRealVar*)(&(fProxSet[code-2]->arg())))->getBins();
    Double_t delta=vrange/nbins/2;
    //Double_t delta=vrange/nbins;
    for(Int_t ie=0;ie<NEv;ie++){
      fEvTree->GetEntry(ie);
      if(!CheckRange(rangeName)) continue;
     //only inlcude events within same bin as vval in integral
      if(TMath::Abs(fMCVar[code-2]-vval)>delta)continue;
      integral+=evaluateMC()*GetIntegralWeight();
    }
    //correct for delta integration width
    //first 2 case near range limits 
    if((rmax-vval)<delta) delta=delta+rmax-vval;
    else if((vval-rmin)<delta) delta=delta+vval-rmin;
    else delta=delta*2;
    integral=integral/delta;
  }
  //else return 1;

  //Set Last[0] so we can just return that if no parameter changes
  if(fNMCGen) fLast[0]=integral/fNMCGen;
  else fLast[0]=integral/NEv;


  return fLast[0];
}

Bool_t RooHSEventsPDF::CheckRange(const char* rangeName) const{
  for(UInt_t i=0;i<fProxSet.size();i++){
    RooRealVar* var=((RooRealVar*)(&(fProxSet[i]->arg())));
    if(!var->inRange(fMCVar[i],rangeName)) return kFALSE;
  }
  return kTRUE;

}
Bool_t RooHSEventsPDF::CheckChange() const{
  //Note analytical integral is const funtion so can only change data members
  //which are pointed to something, thus need Double_t *fLast
  //and construct a N-D array where we can change elements

  Bool_t hasChanged=0;
  for(Int_t i=1;i<fNpars+1;i++)
    if(fLast[i]!=(*(fParSet[i-1]))) hasChanged=1;
  if(hasChanged){
    for(Int_t i=1;i<fNpars+1;i++)
      fLast[i]=*(fParSet[i-1]);
  }
  return hasChanged;
}
Bool_t RooHSEventsPDF::SetEvTree(TChain* tree,TString cut,Long64_t ngen){
  if(!tree->GetEntries()) return kFALSE;
  return SetEvTree(tree->CloneTree(),cut,ngen);
}
Bool_t RooHSEventsPDF::SetEvTree(TTree* tree,TString cut,Long64_t ngen){
  if(!tree->GetEntries())return kFALSE;
  Info("RooHSEventsPDF::SetEvTree"," with name %s and cut %s",tree->GetName(),cut.Data());
  fCut=cut;
  
  if(fCut.Sizeof()!=1)fEvTree=tree->CopyTree(fCut);
  else fEvTree=tree;
  fCut=""; //remove cut once events filtered as branches may be removed 
  fNMCGen=ngen;
  fConstInt=fEvTree->GetEntries();
   // RooArgList set1=fVarSet.at(0);
  //fMCVar.reserve(fProxSet.size());
  fEvTree->ResetBranchAddresses();
  fEvTree->SetBranchStatus("*",0);
  fBranchStatus=kTRUE;
  fMCCat.clear();
  fMCCat.reserve(fCatSet.size());
  fMCVar.ResizeTo(fProxSet.size());
  //fMCCat.ResizeTo(fCatSet.size());
  for(UInt_t i=0;i<fProxSet.size();i++){
    if(fEvTree->GetBranch(fProxSet[i]->GetName())){
      //  Double_t tempD=0;
      //fMCVar.push_back(tempD); //initalise element
      //fMCVar[i]=0;
      fEvTree->SetBranchStatus(fProxSet[i]->GetName(),1);
      fEvTree->SetBranchAddress(fProxSet[i]->GetName(),&fMCVar[i]);
      if(fEvTree->GetBranch(TString("gen")+fProxSet[i]->GetName())) fEvTree->SetBranchStatus(TString("gen")+fProxSet[i]->GetName(),1);
      }
    else{
      Warning("RooHSEventsPDF::SetEvTree","Branch %s not found",fProxSet[i]->GetName()); //May still get set as prototype data
      fBranchStatus=kFALSE;
    }
  }
  for(UInt_t i=0;i<fCatSet.size();i++){
    cout<<fCatSet[i]->GetName()<<" "<<endl;
    // Int_t tempI=0;
    //fMCCat.push_back(tempI); //initalise element
    if(fEvTree->GetBranch(fCatSet[i]->GetName())){
      fEvTree->SetBranchStatus(fCatSet[i]->GetName(),1);
      fEvTree->SetBranchAddress(fCatSet[i]->GetName(),&fMCCat[i]);
      if(fEvTree->GetBranch(TString("gen")+fCatSet[i]->GetName())) fEvTree->SetBranchStatus(TString("gen")+fCatSet[i]->GetName(),1);
      }
    else{
      Warning("RooHSEventsPDF::SetEvTree","Branch %s not found",fCatSet[i]->GetName()); //May still get set as prototype data
      fBranchStatus=kFALSE;
    }
  }
  // cout<<"CHECK TREE "<<endl;
  fEvTree->GetEntry(0);
  cout<<"SetEvTree "<<fMCVar[0]<<" "<<fBranchStatus<<" "<<fConstInt<<endl;

  if(fCheckInt) CheckIntegralParDep(fCheckInt);
  return fBranchStatus;
}
void  RooHSEventsPDF::CheckIntegralParDep(Int_t Ntests){
  fCheckInt=Ntests;
  if(!fEvTree) return; //will check later when tree is set
  if(!fBranchStatus) return; //missing branches may need to add protodata
  
  Long64_t saveNint=fNInt;
  // fNInt=nint;//REOMVE FOR NOW nint is not used
  fNInt=fEvTree->GetEntries();//just use all entries
  //scale Ntests by Ndimensions
  Ntests=TMath::Power((Double_t)Ntests,(Double_t)fParSet.size());

  Info("RooHSEventsPDF::CheckIntegralParDep","Going to run %d calculations of integral with random parameters",Ntests);
  
  RooRealVar integral("integral","integral",0,0,2);
  integral.setError(sqrt(fNInt)/fNInt); //Error needs to be set before entering in ds
 RooDataSet ds("intds","intds",RooArgSet(integral));
  //want to set random paramter values
  //loop over each paramter and calculate integral
 vector<Double_t> SavedPars;
 for(UInt_t ip=0;ip<fParSet.size();ip++){//first save parameters
   RooRealVar* par=((RooRealVar*)(&(fParSet[ip]->arg())));
   SavedPars.push_back(par->getValV());
 }
 for(Int_t ir=0;ir<Ntests;ir++){ //loop over tests
    for(UInt_t ip=0;ip<fParSet.size();ip++){//loop over parameters
      RooRealVar* par=((RooRealVar*)(&(fParSet[ip]->arg())));
      par->setVal((par->getMax("")-par->getMin(""))*RooRandom::uniform()+par->getMin(""));
    }
    //Now calculate integral
    integral.setVal(analyticalIntegral(1,""));
    ds.add(RooArgSet(integral));
  }

 
  Double_t low=0;
  Double_t high=0;
  ds.getRange(integral,low,high);
  integral.setRange(low,high);
  RooPlot *frame=integral.frame();
  ds.plotOn(frame);

  frame->Draw();

  RooRealVar  mean("mean","mean",ds.mean(integral));
  RooRealVar pvar("IntPull","Integral Pull Dist.",-5,5);
  RooPullVar pull("IntPull","Integral Pull Dist.",integral,mean);
  ds.addColumn(pull,kFALSE);

  ds.Print();
  ds.getRange(pvar,low,high);
  pvar.setRange(low,high);
  RooPlot *framePull=pvar.frame();
  ds.plotOn(framePull);
  RooRealVar mp("mp","mp",0,-5,5) ;
  RooRealVar sp("sp","sp",1,0,100);
  RooGaussian gp("gp","gp",pvar,mp,sp);
  gp.fitTo(ds);
  gp.paramOn(framePull);
  gp.plotOn(framePull);
  
  pvar.Print();
  new TCanvas();
  framePull->Draw();
  
 
  //numerical check gives constant integral, can force const for fit speed
  fConstInt=mean.getVal();
  fNInt=saveNint;
  if(sp.getVal()<2){
    Info("RooHSEventsPDF::CheckIntegralParDep","Numerical integral constant. Will not recalculate during fits to save time, if you want to force calculation do not call this function");
    SetConstInt();
  }
  //reset pars
  for(UInt_t ip=0;ip<fParSet.size();ip++){//first save parameters
    RooRealVar* par=((RooRealVar*)(&(fParSet[ip]->arg())));
    par->setVal(SavedPars[ip]);
  }
  fCheckInt=kFALSE; //only do once
}
Bool_t RooHSEventsPDF::AddProtoData(RooDataSet* data){
  //merge the current tree with data from another dataset
  //Default it will add any branches in data not in fEvTree
  cout<<"RooHSEventsPDF::AddProtoData "<<data<<" "<<fEvTree<<endl;
  if(!fEvTree) return kFALSE;
  if(!fEvTree->GetEntries()) return kFALSE;
  
  //Loop over data branches and check if any missing
  const RooArgSet *dataVars=data->get();
  RooArgSet thisVars =VarSet(0);
  thisVars.Print();
  dataVars->Print();
  
  //Each entry in fEvTree should be given a random value of new variables
  vector<Long64_t> vrandom(data->numEntries());
  for(Long64_t ir=0;ir<data->numEntries();ir++)
    vrandom[ir]=ir;
  cout<<vrandom.size()<<endl;
  std::random_shuffle(vrandom.begin(),vrandom.end());

  TIter iter=dataVars->createIterator();
  vector<Double_t> brD;
  vector<Int_t> brI;
  Int_t ND=0;
  Int_t NI=0;
  vector<Bool_t> isD;
  vector<Int_t> index;
  TObjArray branches;
  //Add new branches to EvTree
  cout<<"Add branches "<<endl;
  while(RooAbsArg* arg=(RooAbsArg*)iter()){
    if(fEvTree->GetBranch(arg->GetName())) continue; //already exists
    if(dynamic_cast<RooAbsCategory*>(arg)){
      index.push_back(NI+ND);
      Int_t tempI=0;
      brI.push_back(tempI);
      branches.Add(fEvTree->Branch(arg->GetName(),&brI[NI++],TString(arg->GetName())+"/I")); 
      isD.push_back(kFALSE);
    }
    else if(dynamic_cast<RooAbsReal*>(arg)){
      index.push_back(NI+ND);
      Int_t tempD=0;
      brD.push_back(tempD);
      branches.Add(fEvTree->Branch(arg->GetName(),&brD[ND++],TString(arg->GetName())+"/D"));
      isD.push_back(kTRUE);
    }
    Info("RooHSEventsPDF::AddProtoData","Added data branch %s",arg->GetName());
  }


  //Loop over data in random order and fill new branches
  Long64_t idata=0;
  if(!(NI+ND)) return kTRUE;
  for(Long64_t id=0;id<fEvTree->GetEntries();id++){
    dataVars=data->get(vrandom[idata]);
    Int_t nd=0;
    Int_t ni=0;
    for(Int_t ib=0;ib<NI+ND;ib++){
      if(isD[ib]){
	brD[nd]=dataVars->getRealValue(branches[ib]->GetName());
	nd++;
      }
      else{
	brI[ni]=dataVars->getCatIndex(branches[ib]->GetName());
	ni++;
      }
      
      ((TBranch*)branches[ib])->Fill();
    }
    if(idata==(Long64_t)vrandom.size()-1){//Need to reuse data until done all MC events
      std::random_shuffle(vrandom.begin(),vrandom.end());
      idata=0;
    }
    idata++;
  }
  
  
  SetEvTree(fEvTree,fCut);
  fEvTree->Print();
  fEvTree->GetEntry(0);
 
  return fBranchStatus;  
}
void RooHSEventsPDF::ResetTree(){

  if(fEvTree) {delete fEvTree;fEvTree=nullptr;}
}
