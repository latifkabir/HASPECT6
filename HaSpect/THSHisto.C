//--Author      DI Glazier 30/06/2014
//--Rev
//--Update
//--Description
//HASPECT Event Reconstruction
//THSHisto
//Facilitates histogram declaration and filling
//Users must implement :
//HistogramList to declare which types of histogram they need
//FillHistogram to fill the histogram with datamembers of their selector class
//If the histogram fHisbins is defined it will be used to split
//the histograms into different (kinematic) bins

#define THSHisto_cxx

#include <TProofOutputFile.h>
#include "THSHisto.h"

void THSHisto::InitialiseHisto(TList *input){ 
  if(input) fBins=dynamic_cast<THSBins*>(input->FindObject("HSBins"));//To get fBins in SlaveBegin
  if(fBins)fBins->InitialiseBins();
  ChangeNames();
  fIsHisto=kTRUE;
}

void THSHisto::LoadCut(TString sCut){
  //Declare the cuts you want to use, just takes a string name
  //the user is responsible for implementing the cut themselves
  //via a call to FillHistograms(sCut,bin) in an apporpriate place
  cout<<"THSHisto::LoadCut "<<sCut<<endl;
  fCurrCut=sCut; //keep the cut name for use in MapHist()
  fMapCuts[sCut]=fVecCuts.size();
  fVecCuts.push_back(sCut); //Enter the cut into the map
}

void THSHisto::LoadHistograms(){
  //Loads histograms for each defined bin and cut
  //The histograms required should be given in the users HistogramList Function
  //Note HistogramList should give the histograms to fOutput, which will take
  //responsibility for deleting them as SetOwner() is called in its contructor
  //GetKinematic bins
  Int_t Nbins=0;
  // if(fHisbins){
  //   if(TString(fHisbins->ClassName()).Contains("1"))Nbins=(fHisbins->GetNbinsX()+2);//+2 for overflow and underflow;
  //   else if(TString(fHisbins->ClassName()).Contains("2"))Nbins=(fHisbins->GetNbinsX()+2)*(fHisbins->GetNbinsY()+2);//+2 for overflow and underflow;
  //   else if(TString(fHisbins->ClassName()).Contains("3"))Nbins=(fHisbins->GetNbinsX()+2)*(fHisbins->GetNbinsY()+2)*(fHisbins->GetNbinsZ()+2);//+2 for overflow and underflow;
  // }
  if(fBins) Nbins=fBins->GetN();
  Info("THSHisto::LoadHistograms()","Nbins = %d",Nbins);
  UInt_t ncut=0;
  if(!Nbins){  //Only 1 bin
    fCurrBin=0;
    //loop over cuts
    for(ncut=0;ncut<fVecCuts.size();ncut++){
      fVecBinNames.resize(1);
      fVecBinNames[0]=TString("_All");
      fCurrCut=fVecCuts[ncut];
    //Call the users Histogram list for each declared cut
      HistogramList(fVecCuts[ncut] + fVecBinNames[0]); //cut + bin gets added to histogram name 
      fVecHistCut.push_back(fHistNameMap); //enter the histograms for this cut into a vector element
      fHistNameMap.clear(); //clear the previous cut histograms
    }
    fVecHistBin.push_back(fVecHistCut);//Add all the histograms for this kin bin
    fVecHistCut.clear(); //clear the previous bin histograms  
    return;
  }
  //keep a list of events in each bin
  // fBinEntryLists=new TObjArray(Nbins);
  // fBinEntryLists->SetName("HSKinBinsList");
  // fBinEntryLists->SetOwner(kTRUE);
  //ifmultiple kinematic bins defined make a histogram for each bin
  //TString sBin;
  for(Int_t ib=0;ib<Nbins;ib++){
    //sBin.Form("_Bin%d",ib); //make bin label
    //Load all the cuts for this bin
    ncut=0;
    //iterate over cuts
    //for(vector< TString >::iterator nextcut=fVecCuts.begin();nextcut!=fVecCuts.end();++nextcut){
    for(ncut=0;ncut<fVecCuts.size();ncut++){
      //Call the users Histogram list for each declared cut
      fCurrCut=fVecCuts[ncut];fCurrBin=ib;
      HistogramList(fVecCuts[ncut]+fVecBinNames[ib]); //construct all histograms for this bin,
      fVecHistCut.push_back(fHistNameMap); //enter the histograms for this cut into a vector element
      fHistNameMap.clear(); //clear the previous cut histograms
      
    }
    fVecHistBin.push_back(fVecHistCut);//Add all the histograms for this kin bin
    fVecHistCut.clear(); //clear the previous bin histograms
    //Add an entry list for each kinematic bin
    // TEntryList* el=new TEntryList(TString("HSBin_")+fVecBinNames[ib],TString("HSBin_")+fVecBinNames[ib]);
    // fBinEntryLists->Add(el);
    // if(dynamic_cast<TSelector*>(this)) dynamic_cast<TSelector*>(this)->GetOutputList()->Add(el);
  }
}
TH1* THSHisto::MapHist(TH1* hist){
  //Place histograms in map for quick look up via their name
  //This should be called in HistogramList and returns the histogram to add to the slector output
  hist->Sumw2(); //Call sumw2 function for all histograms by default; note doesn't make particular sense to do this here but it is most convenient
  //Store histogram type name
  if(std::find(fVecTypes.begin(), fVecTypes.end(),TString(hist->GetName()) ) == fVecTypes.end()) fVecTypes.push_back(TString(hist->GetName()));

  //  if(!fVecTypes.contains(TString(hist->GetName()))) fVecTypes.push_back(TString(hist->GetName()));

  TString sLabel;
  sLabel=fCurrCut+fVecBinNames[fCurrBin];
  hist->SetName(TString(hist->GetName())+"_"+sLabel);
  fHistNameMap[TString(hist->GetName())]=hist;
  return hist;
}
// Int_t THSHisto::GetKinBin(Float_t ev1,Float_t ev2,Float_t ev3){
//   //Function that finds the bin number and converts it to a string
//   Int_t evBin=0;
//   if(fHisbins)evBin=fHisbins->FindBin(ev1,ev2,ev3); //find bin
//   return fKinBin=evBin;
// }
TH1* THSHisto::FindHist(TString hname){
  //function that finds the correct histogram for the user to fill
  //requires the current cut and bin be set beforehand, usually in FillHistograms  
  TString sLabel;
  sLabel=fCurrCut+fVecBinNames[fCurrBin];
  return fVecHistBin[fCurrBin][fMapCuts[fCurrCut]][hname+"_"+sLabel];
}
// TString THSHisto::GetStrBin(Float_t ev1,Float_t ev2){
//   //Function that finds the bin number and converts it to a string
//   Int_t evBin=0;
//   if(fHisbins)evBin=fHisbins->FindBin(ev1,ev2); //find bin
//   TString sBin;
//   sBin.Form("_Bin%d",evBin); //make label
//   return sBin;
// }
void THSHisto::ChangeNames(/*TFile* list*/){
  //make a map from fHisbins containing name tags
  if(!fBins) return;
  //We got a fBins histogram, now make the name map
  if(fBinMap.empty()){
    //UInt_t Naxis=fBins->GetNAxis();
    //Int_t Nbins=fBins->GetN();
    fVecBinNames.resize(fBins->GetN());
    for(Int_t ib=0;ib<fBins->GetN();ib++){
      fVecBinNames[ib]=TString("_")+fBins->GetBinName(ib);}
  }

  return;
      
  //   TString tagname1,tagname2,tagname3,numb,binstr;
  //   UInt_t xyz=0;
  //   //note +2 for overfllow+ underflow
  //   fVecBinNames.resize(Nbins);
  //   for(Int_t xi=0;xi<fHisbins->GetXaxis()->GetNbins()+2;xi++){
  //     tagname1=fHisbins->GetXaxis()->GetTitle();
  //     numb.Form("%1.2f_",fHisbins->GetXaxis()->GetBinCenter(xi));
  //     tagname1+=numb;
  //     if(Naxis==1){
  // 	binstr.Form("_Bin%d",xi);
  // 	//	fBinMap[binstr]=tagname1;
  // 	fVecBinNames[xi]=(tagname1);
  // 	continue;
  //     }
  //     for(Int_t yi=0;yi<fHisbins->GetYaxis()->GetNbins()+2;yi++){
  // 	tagname2=fHisbins->GetYaxis()->GetTitle();
  // 	numb.Form("%1.2f_",fHisbins->GetYaxis()->GetBinCenter(yi));
  // 	tagname2+=numb;
  // 	if(Naxis==2){
  // 	  xyz=fHisbins->GetBin(xi,yi);
  // 	  binstr.Form("Bin%d",xyz);
  // 	  //fBinMap[binstr]=tagname1+tagname2;
  // 	  fVecBinNames[xyz]=(tagname1+tagname2);
  // 	  continue;
  // 	}
  // 	for(Int_t zi=0;zi<fHisbins->GetZaxis()->GetNbins()+2;zi++){
  // 	  tagname3=fHisbins->GetZaxis()->GetTitle();
  // 	  numb.Form("%1.2f_",fHisbins->GetZaxis()->GetBinCenter(zi));
  // 	  tagname3+=numb;
  // 	 if(Naxis==3){
  // 	   xyz=fHisbins->GetBin(xi,yi,zi);
  // 	   binstr.Form("Bin%d",xyz);
  // 	   //fBinMap[binstr]=tagname1+tagname2+tagname3;
  // 	   fVecBinNames[xyz]=(tagname1+tagname2+tagname3);
  // 	   continue;
  // 	 } 
  // 	}
  //     }
  //   }
  // }
  
}
// void THSHisto::SetBinEntryListsTree(TTree* elt){
//   if(!fBinEntryLists)return;
//   for(Int_t ibel=0;ibel<fBinEntryLists->GetEntries();ibel++)
//     ((TEntryList*)(fBinEntryLists->At(ibel)))->SetTree(elt);
// }

/* void THSHisto::DrawAll(){ */
/*   //one canvas per histogram */
/*   //one pad per bin (max 8 per canvas) */
/*   //all cuts drawn on same plot */

/*   //fVecHistBin[fCurrBin][fMapCuts[fCurrCut]][hname]; */

/* } */

void THSHisto::LoadWeights(TString fname,TString wname){
  //Load THSWeights class for weighting histograms with
  //TObjArray* OptList = opt.Tokenize(":");
  //if(OptList->GetEntries()>1){
  // TString fname=((TObjString*)OptList->At(0))->String();
  // TString wname=((TObjString*)OptList->At(1))->String();
    Info("LoadWeights"," Opening weights file %s name %s",fname.Data(),wname.Data());
    fWeights=new THSWeights();
    fWeights->LoadSaved(fname,wname);
    fWeights->PrintWeight();
    // }
  
}
void THSHisto::LoadCutsForWeights(){
  //Optional load a histogram cut for each species in THSWeights
    if(fWeights){
     StrIntMap_t spec=fWeights->GetSpecies();
     for(StrIntMap_t::iterator its=spec.begin();its!=spec.end();++its)
       THSHisto::LoadCut(its->first);
    }
}
void THSHisto::FillCutsForWeights(){
  if(fWeights){//Fill histogram for each species in weights
    StrIntMap_t spec=fWeights->GetSpecies();
    for(StrIntMap_t::iterator itss=spec.begin();itss!=spec.end();++itss){
      SetWeight(itss->first);
      FillHistograms(itss->first);
    }
  }
  fWeight=1;
}


void THSHisto::GetWeightEvent(Double_t id){
  if(fWeights)fWeights->GetEntryBinarySearch((Long64_t)id);
}
void THSHisto::SetWeight(TString species){
  if(fWeights) fWeight=fWeights->GetWeight(species);
}
void THSHisto::OrganiseSavedHists(TFile* saveFile,TList* output){
  //Make directory based on cut and bin to save hist in
  //write hist and reset (set entries to 0)
  if(!fIsHisto) return;
  if(fBins)
    if(fBins->GetN())//if not already called
      fBins->InitialiseBins(); //Get all bin names
  TList* DoneList=new TList();
  for(UInt_t ncut=0;ncut<fVecCuts.size();ncut++){
    saveFile->mkdir(fVecCuts[ncut]);
    saveFile->cd(fVecCuts[ncut]);
    TDirectory* cutDir=gDirectory;
    //Loop over axis(i.e. binned variables) and make dirs
    //	for(Int_t ia=0;ia<fBins->GetNAxis();ia++)
    MakeAxisDir(0,cutDir);
    if(!fBins){//no bins just save hists in cuts dir
      for(UInt_t ity=0;ity<fVecTypes.size();ity++){
	TString hname=fVecTypes[ity]+"_"+fVecCuts[ncut]+fVecBinNames[0];
	TH1* hh=(TH1*)output->FindObject(hname);
	hh->Write(fVecTypes[ity]);
	DoneList->Add(hh);
	hh->Reset("");
      }
      continue; //don't do bins stuff
    }
   //Now make loop over bins, get hist and write
     for(Int_t iB=0;iB<fBins->GetN();iB++){
      cutDir->cd();
      TString dname=fBins->GetBinName(iB);
      dname.ReplaceAll("_","_/");
      //cout<<dname<<endl;
      cutDir->cd(dname);
      TString bname=fBins->GetBinName(iB);
      for(UInt_t ity=0;ity<fVecTypes.size();ity++){
	TString hname=fVecTypes[ity]+"_"+fVecCuts[ncut]+"_"+bname;
	TH1* hh=(TH1*)output->FindObject(hname);
	hh->Write(fVecTypes[ity]);
	DoneList->Add(hh);
	hh->Reset("");
      }//done type loop
    }//done bin loop
  }//done cut loop

  //Look for any extra hists and save them in top level dir
  TIter next(output);
  TObject* ho=0;
  TH1* hf=0;
  saveFile->cd();
  while((ho=next())){
    if((hf=dynamic_cast<TH1*>(ho))){
      if(DoneList->Contains(hf)) continue;
      hf->Write();
      hf->Reset(""); //remove saved histogram entries for next file
    }
    
  }

}

void THSHisto::MakeAxisDir(Int_t iA,TDirectory* curDir){

  if(!fBins) return;
  if(iA>=fBins->GetNAxis()) return;//no more axis
  
  for(Int_t iB=0;iB<fBins->GetAxis(iA).GetNbins();iB++){
    //Make directory, move into it and move to next axis
    curDir->mkdir(fBins->GetPartName(iA,iB));
    curDir->cd(fBins->GetPartName(iA,iB));
    MakeAxisDir(iA+1,gDirectory);
  }
}

void THSHisto::ParticleList(TString name){

  fSelOutput->Add(MapHist(new TH1F(name+" P",name+" P",100,0,10)));
  fSelOutput->Add(MapHist(new TH1F(name+" Theta",name+" Theta",90,0,90)));
  fSelOutput->Add(MapHist(new TH1F(name+" Phi",name+" Phi",90,-180,180)));
  fSelOutput->Add(MapHist(new TH1F(name+" DeltaTime",name+" DeltaTime",100,-10,10)));
  fSelOutput->Add(MapHist(new TH2F(name+" DeltaTimeVP",name+" Theta V P",50,0,10,50,-10,10)));
  fSelOutput->Add(MapHist(new TH2F(name+" ThetaVP",name+" Theta V P",50,0,10,50,0,90)));
  fSelOutput->Add(MapHist(new TH1F(name+" P_RES",name+" P resolution",100,-0.2,0.2)));
  fSelOutput->Add(MapHist(new TH1F(name+" Theta_RES",name+" Theta resolution",100,-10,10)));
  fSelOutput->Add(MapHist(new TH1F(name+" Phi_RES",name+" Phi resolution",100,-10,10)));

}
void THSHisto::FillParticles(TString name,THSParticle* part){
  FindHist(name+" P")->Fill(part->P4p()->Rho(),fWeight);//1D
  FindHist(name+" Theta")->Fill(part->P4p()->Theta()*TMath::RadToDeg(),fWeight);//1D
  FindHist(name+" Phi")->Fill(part->P4p()->Phi()*TMath::RadToDeg(),fWeight);//1D
  FindHist(name+" DeltaTime")->Fill(part->DeltaTime(),fWeight);//1D
  ((TH2D*)FindHist(name+" DeltaTimeVP"))->Fill(part->P4p()->Rho(),part->DeltaTime(),fWeight);
  ((TH2D*)FindHist(name+" ThetaVP"))->Fill(part->P4p()->Rho(),part->P4p()->Theta()*TMath::RadToDeg(),fWeight);
  FindHist(name+" P_RES")->Fill(part->ResRho(),fWeight);//1D
  FindHist(name+" Theta_RES")->Fill(part->ResTheta()*TMath::RadToDeg(),fWeight);//1D
  FindHist(name+" Phi_RES")->Fill(part->ResPhi()*TMath::RadToDeg(),fWeight);//1D

}
