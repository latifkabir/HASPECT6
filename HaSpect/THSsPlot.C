/**
 * \class THSsPlot
 * 
 * RooFit interface fit class
 * 
 * 
 * @example FitHSSimple.C
 * 
*/

#include "THSsPlot.h"
#include <TH1.h>
#include <TH2.h>
#include <TLegend.h>
#include <algorithm>      // std::sort

ClassImp(THSsPlot);


THSsPlot::THSsPlot() : THSRooFit(),fSPlot(0),fWeights(0){
  fSRange[0]=0;fSRange[1]=0; 
};

THSsPlot::THSsPlot(TString name) : THSRooFit(name){
  fSRange[0]=0;fSRange[1]=0; 
}

THSsPlot::THSsPlot(TString name,RooWorkspace* ws) : THSRooFit(name,ws){
  fSRange[0]=0;fSRange[1]=0; 
}

THSsPlot::THSsPlot(THSsPlot* rf){
   THSRooFit(dynamic_cast<THSRooFit*>(rf));
   //copy constructor, but do not copy the data tree, load that explicitly
   fSaveWeights=rf->fSaveWeights;
}

THSsPlot::~THSsPlot(){
  if(fWeights) delete fWeights;
  //  if(fSPlot) delete fSPlot; //crashes FitSavedBins with weights
}

//////////////////////////////////////////////////////////
/// Function to combine best fit pdfs into single PDF. 
/// This allows sPlot fit to just be dependent on 2 yields
/// where one is defined by fSingleSp and the other is the 
/// sum of the rest.
void THSsPlot::MergeModelSpecies(){
  if(!fWS->pdf(fSingleSp)) {cout<<"THSsPlot::MergeModelSpecies() single Species not found "<<fSingleSp<<endl;}
  RooArgList yields=fYields;
  RooRealVar* sYield=(RooRealVar*)yields.find(fYld+fSingleSp);
  yields.remove(*sYield);
  RooArgSet addpdfs = ((RooAddPdf*)fModel)->pdfList();
  RooAbsPdf* sPdf=(RooAbsPdf*)addpdfs.find(fSingleSp);
  addpdfs.remove(*sPdf);
  RooAbsArg* icoef=0;	
 
  //sPlot Range change must recalculate yields for new range!
   if(fSRange[0]!=fSRange[1]){
    RooRealVar* var=((RooRealVar*)&fVariables[0]);
    var->setRange("sRange",fSRange[0],fSRange[1]);
    RooFIter pdfIter = addpdfs.fwdIterator();
    RooAbsArg* ipdf=0;	
    addpdfs.Print();
    yields.Print();
    Int_t iy=0;
    while((ipdf=pdfIter.next())){ 	     
      Double_t inte1=((RooAbsPdf*)ipdf)->createIntegral(fVariables)->getVal();
      Double_t inte2=((RooAbsPdf*)ipdf)->createIntegral(fVariables,RooFit::Range("sRange"))->getVal();
      Double_t yd=((RooRealVar*)&yields[iy])->getVal();
      yd=yd*inte2/inte1;
      cout<<"THSsPlot::MergeModelSpecies() chenge rangeo of "<<((RooRealVar*)&yields[iy])->GetName()<<" "<<((RooRealVar*)&yields[iy])->getVal()<<" by "<<inte2/inte1<<" to "<<yd<<endl;
      ((RooRealVar*)&yields[iy])->setVal(yd);
      iy++;
    }	
    //change yield of single pdf too
    Double_t inte1=sPdf->createIntegral(fVariables)->getVal();
    Double_t inte2=sPdf->createIntegral(fVariables,RooFit::Range("sRange"))->getVal();
    Double_t yd=sYield->getVal();
    yd=yd*inte2/inte1;
    cout<<"THSsPlot::MergeModelSpecies() chenge rangeo of "<<sYield->GetName()<<" "<<sYield->getVal()<<" by "<<inte2/inte1<<" to "<<yd<<endl;
  
   }
   yields.setAttribAll("Constant");
   RooAddPdf* rap=new RooAddPdf(fSingleSp+"Merged",fSingleSp+"Merged",addpdfs,yields);
   rap->getParameters(*fData)->setAttribAll("Constant");
  //fWS->import(*rap);

  //remove all pdfs and add back single species and summed pdfs
  fPDFs.releaseOwnership();
  fPDFs.removeAll();
  fPDFs.add(*(sPdf));
  fPDFs.add(*rap);

  fYields.releaseOwnership();
  fYields.removeAll();
  fYields.add(*(sYield)); //add back in single species yld
  
  //construct new yield parameter for sum of other contributions
  Double_t Ym=0; //first sum the fitted yields
  for(Int_t iy=0;iy<yields.getSize();iy++){
    Ym+=((RooRealVar*)&yields[iy])->getVal();
    //cout<<iy<<"THSsPlot::MergeModelSpecies() CHECKING COEF "<<((RooRealVar*)&yields[iy])->getVal()<<endl;	
  }
  cout<<"Ym "<< Ym<<endl;
  RooRealVar* ymerge=(RooRealVar*)fWS->factory(fYld+"Merged[0,0,1]");
  (ymerge)->setMin(0,0);
  (ymerge)->setMax(0,2*Ym);
  (ymerge)->setVal(Ym);
  fYields.add(*ymerge);
  ymerge->Print();
  cout<<"THSsPlot::MergeModelSpecies() Print pdfs and yeilds"<<endl;
  fPDFs.Print();
  fYields.Print();
  TotalPDF();
 
}

////////////////////////////////////////
///Perform fit of yields to get sWeights
void THSsPlot::sPlot(){
  if(fTree)
    cout<<"THSsPlot::sPlot() weights will be synchronised to branch "<<fIDBranchName<<" "<<fTree<<" "<<fTree->GetBranch(fIDBranchName)<<endl;
  //Check if there is an ID branch for synchronisation
  //check if ID branch exists
  if(fTree){
    if(fTree->GetBranch(fIDBranchName)){
      fGotID=kTRUE;
      cout<<"THSsPlot::sPlot Loading ID branch data "<<fIDBranchName<<endl;
     }
    else {cout<<"THSsPlot::sPlot() warning No ID branch called will just use tree entry number. This means weights will fail to synchronise with initial tree if you cut events or use bins"<<fIDBranchName<<endl;}
  }
  //Option here to sum background contributions and fix relative yields
  if(fSingleSp!=TString("")) MergeModelSpecies();
  //Set all the model parameters to be constant, so only fitting the species yields
  //calculate the sWeights 
  if(fSRange[0]!=fSRange[1]){//sub range set for sPlot fit, remake dataset
   ((RooRealVar*)&fVariables[0])->setRange(0,fSRange[0],fSRange[1]);
    RooArgSet dataVars(fVariables,fAuxVars);
    if(fID) dataVars.add(*fID);
    RooDataSet* tmpData=0;
    if(fTree) tmpData=new RooDataSet(fTree->GetName(),fTree->GetTitle(),dataVars,RooFit::Import(*fTree));
    else tmpData=new RooDataSet(fData->GetName(),fData->GetTitle(),(RooDataSet*)fData,dataVars);
    fData->SetName("Empty");fData->reset();fData=0;
    fData=tmpData;
    SetDataWeight();//if weighted data reassign the weight
    fFiti++;
    if(fIsPlot){
      PlotDataModel();
      ((TCanvas*)fCanvases->Last())->SetTitle(Form("%s Model on new Range = %lf",GetName(),fChi2));
    }
    
  }
  Info("THSsPlot::sPlot()"," about to start");
  if(fData->sumEntries()<1){cout<<"Warning : void THSsPlot::sPlot() negative data entries no sPlot"<<endl;return ;}//might be -ve with weights, in which case it may crash...also <1
  fParameters.setAttribAll("Constant");
  if(fInWeights) //need to clone datset if want to draw results with in weights
    fSPlot = new RooStats::SPlot("splot_fit", "sPlot Fit", *((RooDataSet*)fData),fModel ,fYields,RooArgSet(),kTRUE,kTRUE);
  //fSPlot = new RooStats::SPlot("splot_fit", "sPlot Fit", *((RooDataSet*)fData),fModel ,fYields);
  else
    fSPlot = new RooStats::SPlot("splot_fit", "sPlot Fit", *((RooDataSet*)fData),fModel ,fYields);
  fParameters.setAttribAll("Constant",kFALSE);
  
  //Check that the fit was succesfull
  Double_t TotalYield=0;
  for(Int_t iy=0;iy<fYields.getSize();iy++)
    TotalYield+=((RooRealVar*)&fYields[iy])->getVal();

  fFiti++;
  if(TotalYield>0){ //got some weights
    fWeights=new THSWeights("HSsWeights");//initialise weights
    fWeights->SetTitle(GetName());
    fWeights->SetFile(fOutDir+TString("Weights")+GetName()+".root");
    ExportWeights();
    if(fIsPlot){
      PlotDataModel();
      ((TCanvas*)fCanvases->Last())->SetTitle(Form("%s sWeights Fit chi2 = %lf",GetName(),fChi2));
    }
  }
  else Warning("THSsPlot::sPlot()"," total weights 0, fit did not converge. Make sure the non-sweight fit to fix parameters was succesful. No weights will be assigned for these events");
  
}

Double_t THSsPlot::GetSWeight(Long64_t event,TString species){
  //Note species should include the prepend fYld
   return fSPlot->GetSWeight(event,species);
}

void THSsPlot::ExportWeights(TString wmname){
  cout<<"THSsPlot:ExportWeights "<<endl;
  //fWeights=new THSWeights(wmname);//initialise map
  Int_t NSpecies=fYields.getSize();
  TVectorD eventW(NSpecies); //initialise weights vector
  for(Int_t iw=0;iw<NSpecies;iw++)//set name for each species, remove leading Yld_	
    fWeights->SetSpecies(TString(fYields.at(iw)->GetName()).Remove(0,4));
  for(Long64_t ev=0;ev<fData->numEntries();ev++){//loop over events
    for(Int_t iw=0;iw<NSpecies;iw++)//loop over species
      eventW[iw]=GetSWeight(ev,fYields.at(iw)->GetName());//get weight for this species
    if(fGotID){//use ID from initial tree
      const RooArgSet* vars=fData->get(ev);
      fWeights->FillWeights((Long64_t)vars->getRealValue(fIDBranchName),eventW);
      } //ID not defined just use entry number on tree
    else fWeights->FillWeights(ev,eventW);
  }
}

void THSsPlot::ExportWeights1(TString wmname){
  cout<<"THSsPlot:ExportWeights1 no sPlot fit done as only 1 species left"<<endl;
  Int_t NSpecies=fYields.getSize();
  TVectorD eventW(NSpecies); //initialise weights vector
  for(Int_t iw=0;iw<NSpecies;iw++)//set name for each species, remove leading Yld_	
    fWeights->SetSpecies(TString(fYields.at(iw)->GetName()).Remove(0,4));
  for(Long64_t ev=0;ev<fData->numEntries();ev++){//loop over events
    for(Int_t iw=0;iw<NSpecies;iw++)//loop over species
      eventW[iw]=1;//get weight for this species
    if(fGotID){//use ID from initial tree
      const RooArgSet* vars=fData->get(ev);
      fWeights->FillWeights((Long64_t)vars->getRealValue(fIDBranchName),eventW);
    } //ID not defined just use entry number on tree
    else fWeights->FillWeights(ev,eventW);
  }
}

///////////////////////////////////////////////////////////////
///This function gives a weight=0 for each entry in the tree
///if event not in fit (outside limits),
///as such the weight tree could be used with TTree::AddFriend.
void THSsPlot::ExportTreeWeights(TString wmname){
  cout<<"THSsPlot::ExportTreeWeights "<<fTree->GetName()<<" "<<fTree->GetEntries()<<" "<<fGotID<<endl;
  if(!fTree) {cout<<"THSsPlot::ExportTreeWeights( No tree found"<<endl;return;}
  Int_t NSpecies=fYields.getSize();
  TVectorD eventW(NSpecies); //initialise weights vector
  for(Int_t iw=0;iw<NSpecies;iw++)//set name for each species,remove "Yld_"
    fWeights->SetSpecies(TString(fYields.at(iw)->GetName()).Remove(0,4));
  Double_t tid=0;
  if(fGotID)  fTree->SetBranchAddress(fIDBranchName,&tid);
  Long64_t Nds=0;
  for(Long64_t ev=0;ev<fTree->GetEntries();ev++){//loop over events
    fTree->GetEntry(ev);
    if(fGotID){//use ID from initial tree
      const RooArgSet* vars=fData->get(Nds);
      Long64_t id=-1;
      if(vars) id=(Long64_t)vars->getRealValue(fIDBranchName);
      if(id==(Long64_t)tid){//event in dataset 
	for(Int_t iw=0;iw<NSpecies;iw++)//loop over species
	  eventW[iw]=GetSWeight(Nds,fYields.at(iw)->GetName());//get weight for this species
	Nds++;//go to next event in dataset
      } 
      else{ //not in dataset no weight calculated for this event set to zero
	for(Int_t iw=0;iw<NSpecies;iw++)//loop over species
	  eventW[iw]=0;//set to zero
	//don't increment Nds
      } //ID not defined just use entry number on tree
      fWeights->FillWeights((Long64_t)tid,eventW); //fill the weight for every tree event
    }
    else {//need to be same events in dataset as tree for this to be meaningful
      for(Int_t iw=0;iw<NSpecies;iw++)//loop over species
		eventW[iw]=GetSWeight(ev,fYields.at(iw)->GetName());//get weight for      
      fWeights->FillWeights(ev,eventW);//no tree id just fill with entry number
    }
  }
}
 
void THSsPlot::ExportWeightsToFile(TString filename){
  Info("THSsPlot::ExportWeightsToFile"," Exporting Weights to %s",filename.Data());
  if(fWeights->GetWeightList())fWeights->SortWeights();//if merged from differen  
   fWeights->Save();
}

/////////////////////////////////////
///Constructs new map if not existing
void THSsPlot::AddWeightMap(THSWeights *Wts){
  Info("THSsPlot::AddWeightMap"," Start add weights %lld",Wts->Size());
  if(!fWeights){//construct new map if not existing
    cout<<"Construct new weightmap"<<endl;
    fWeights=new THSWeights("WeightMap");
    fWeights->SetTitle(GetName());
    fWeights->SetFile(fOutDir+TString("Weights")+GetName()+".root");
  }
  fWeights->Add(Wts);
}

//////////////////////////////////////////////////////
///Look through sub fits and add the maps to the total
///loop over subfits.
void THSsPlot::AddSubWeights(){
  if(!fWeights){//construct new map if not existing
    fWeights=new THSWeights("WeightMap");
    fWeights->SetTitle(GetName());
  }
  for(Int_t i=0;i<fRooFits->GetEntries();i++)
    fWeights->Add(dynamic_cast<THSsPlot*>(fRooFits->At(i))->GetWeights());
}

void THSsPlot::DrawTreeVar(TString VarName,Int_t nbins,Double_t hmin,Double_t hmax){
  Info("THSsPlot::DrawTreeVar"," %s",VarName.Data());
  if(!fWeights->GetTree())fWeights->LoadSaved(fOutDir+TString("Weights")+GetName()+".root","HSsWeights");
  cout<<fTree<<endl;
  if(fTree->GetBranch(fIDBranchName))fGotID=kTRUE;
  if(fWeights->GetTree()->GetEntries()!=fTree->GetEntries()) {cout<<"DrawTreeVar : Tree and Map have different entries!"<<" " <<fWeights->GetTree()->GetEntries()<<" "<<fTree->GetEntries()<<endl;}
  if(!fWeights->IsSorted()) fWeights->SortWeights();
  Double_t dVar=0;
  fTree->SetBranchStatus(VarName,1);
  if(fTree->SetBranchAddress(VarName,&dVar)) {cout<<"No variable in tree "<<VarName<<endl; return;}
  Double_t tid=0;
  if(fGotID)  fTree->SetBranchAddress(fIDBranchName,&tid);
  vector<TH1F* > hists;
  //create histogram for each species
  Int_t Nspecies=fWeights->GetNSpecies();
  for(Int_t isp=0;isp<Nspecies;isp++){
    TH1F* htemp=new TH1F(fWeights->GetSpeciesName(isp)+VarName+GetName(),VarName+" SWeighted "+fWeights->GetSpeciesName(isp),nbins,hmin,hmax);
    htemp->Sumw2();
    hists.push_back(htemp);
    fHists->Add(hists[isp]); 
  }
  Long64_t id=0; //can be diffent from event number, and can be filtered in dataset
  Long64_t nds=0;
  //loop over tree and fill with weight
  fWeights->SetCurrEntry(0);
  for(Long64_t ient=0;ient<fTree->GetEntries();ient++){
    fTree->GetEvent(ient);
   //get the id for the event if not in tree
    if(!fGotID) id=ient;
    else id=(Long64_t)tid; 
    if(!fWeights->GetEntryBinarySearch(id)) continue;//look for this id in weights tree
    for(Int_t isp=0;isp<Nspecies;isp++)//Fill each species
      hists[isp]->Fill(dVar,fWeights->GetWeight(isp));
  }
  //create canvas and draw
  TCanvas *canvas=0;
  TLegend* leg = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg->SetTextSize(0.03);
  leg->SetFillColor(0);
  
  fCanvases->Add(canvas=new TCanvas(VarName+"SW",VarName+"SW"));//create new canvas for drawing on
  hists[0]->Draw(); //First time so need axis
  leg->AddEntry(hists[0], hists[0]->GetTitle(), "l");
  Double_t maxVal=hists[0]->GetMaximum();
  for(Int_t isp=1;isp<Nspecies;isp++){
    hists[isp]->SetLineColor(isp);
    leg->AddEntry(hists[isp], hists[isp]->GetTitle(), "l");
    hists[isp]->Draw("same");
    if(maxVal<hists[isp]->GetMaximum()) maxVal=hists[isp]->GetMaximum();
  }
  hists[0]->SetMaximum(maxVal+0.05*maxVal);
  leg->Draw();
  fTree->ResetBranchAddresses();
}

THSRooFit*  THSsPlot::CreateSubFitBins(TTree* ctree,TString rfname,Bool_t CopyTree){//events already selected

  //create a fit object for a subset of data either by setting cut
  //or by fTree->SetEntryList prior to calling this function 
  //It will be deleted by this object
  cout<<"THSsPlot::CreateSubFitBins with tree "<<ctree->GetName()<<" "<<CopyTree<<endl;
  THSsPlot* RFa=new THSsPlot(rfname);
  if(fBinnedFit)RFa->SetBinnedFit();
  RFa->SetSingleSpecies(fSingleSp);
  RFa->SetBinDir(fBinDir);
  RFa->SetOutDir(fOutDir);
  if(fInWeights){
    RFa->SetInWeights(fInWeights);
    RFa->SetWeightName(fWeightName);
  }
  RFa->SetNStudyTrials(fNStudyTrials);
  RFa->SetStudyPDF(fStudyPDF);
  RFa->SetSPlotRange(fSRange[0],fSRange[1]);//Extra to THSRooFit
   //Done configuring RF
  fRooFits->Add(RFa);
  RFa->LoadWorkSpace(fWS,GetName());
  RFa->GetVariables().Print();
  //  RFa->SetIDBranchName(fIDBranchName);

  for(Int_t ill=0;ill<fFitOptions.GetSize();ill++)
    RFa->AddFitOption(*((RooCmdArg*)fFitOptions.At(ill)));
 
  TDirectory *saveDir=gDirectory;
  ctree->GetDirectory()->cd();
  if(CopyTree)RFa->LoadDataSet(ctree->CopyTree(""));//will use any EntryList
  else RFa->LoadDataSet(ctree);//use whole tree
  saveDir->cd();
  RFa->SetDataWeight();//if defined weights use them for this dataset
  return RFa;
}

void THSsPlot::SaveHists(TString filename){
  TFile* file=new TFile(filename,"recreate");
  fHists->Write();
  file->Close();
  delete file;
}

///////////////////////////////////////////////////
/// Call this function to generate Nfits fits.
/// Fit the model to data with only species yields as free pars,
/// calculate weights and import to WeightMap.
/// Remove any zero yield to be constant to help fit converge
void THSsPlot::FitAndStudy(Int_t Nfits){
  //Create new fit and load the new bin data tree
  if(!fWS->set(TString(GetName())+"PDFs"))DefineSets();
  if(!fModel) TotalPDF();
  FitMany(Nfits);
  //RooStats::sPlot automoatically adds Sumw2Error() option
  for(Int_t iy=0;iy<fYields.getSize();iy++){
    Double_t  thisYield=((RooRealVar*)&fYields[iy])->getVal();
    if(thisYield<1E-2){
      //Need to remove this pdf all weights weights will not be written for this species
      RooRealVar* yield0=(RooRealVar*)&fYields[iy];
      GetPDFsp()->remove(GetPDFs()[iy]);
      GetYieldsp()->remove(*(GetWorkSpace()->var(yield0->GetName())));
      GetWorkSpace()->removeSet("Yields");
      GetWorkSpace()->removeSet("PDFs");
      GetWorkSpace()->defineSet("Yields",GetYields());
      GetWorkSpace()->defineSet("PDFs",GetPDFs());	
      
    }
  }
  
  if(fYields.getSize()==1){//Only 1 species all weights==1
    fWeights=new THSWeights("HSsWeights");//initialise weights
    fWeights->SetTitle(GetName());
    fWeights->SetFile(fOutDir+TString("Weights")+GetName()+".root");
    ExportWeights1();
    GetWeights()->PrintWeight();
    GetWeights()->SortWeights();
    if(fSaveWeights) GetWeights()->Save();
    return;
  }
  sPlot();
  //save any canvases produced
  Info("THSsPlot::FitAndStudy","Save to %s",(fOutDir+TString("Results")+GetName()).Data());
  SavePlots(fOutDir+TString("Results")+GetName()+".root");
  StudyFit();
  //save weights to file
  if(GetWeights()){
    GetWeights()->PrintWeight();
    GetWeights()->SortWeights();
    if(fSaveWeights) GetWeights()->Save();//don't save if single bin so we can draw
  }
  
}
void THSsPlot::FitSavedBins(Int_t Nfits){
  //do standard THSRooFit 
  THSRooFit::FitSavedBins(Nfits);

  //in addition combine the weights into 1 and load them
  THSWeights* wts=new THSWeights("HSsWeights");
  //Note the output file cannot contain the word Weights (because of Merge), hence Tweights!
  wts->Merge(fOutDir+"/Weights",fOutDir+"/"+GetName()+"Tweights.root","HSsWeights");
  fWeights=wts;
  // fWeights->Save();
}
void THSsPlot::DefaultFitOptions(){
  // AddFitOption(RooFit::Extended());
  // if(fData)
  //   if(fData->isNonPoissonWeighted())AddFitOption(RooFit::SumW2Error(kTRUE));
  AddFitOption(RooFit::NumCPU(1));
  AddFitOption(RooFit::Save(kTRUE));
  AddFitOption(RooFit::Warnings(kFALSE));


}
