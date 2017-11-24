#ifndef THSROOFIT_h
#define THSROOFIT_h
#include <RooWorkspace.h>
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooLinkedList.h>
#include <RooArgList.h>
#include <RooGaussian.h>
#include <RooCmdArg.h>
#include <RooRealVar.h>
#include <RooAbsArg.h>
#include <RooAbsData.h>
#include <RooPlot.h>
#include <TNamed.h>
#include <TSystem.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TList.h>
#include <TFile.h>
#include "THSBins.h"
#include "THSWeights.h"
#include <vector>
#include <map>
using namespace RooFit;
using namespace std;

class RooHS1StepStudy;

class THSRooFit : public TNamed {
 protected:
  RooWorkspace* fWS=nullptr; //RooFit workspace container
  RooArgList fVariables; //variables in tree, to be fitted
  RooArgList fAuxVars; //variables in tree, not to be fitted
  RooArgList fBinVars; //variables in tree, to be used to make bins for fitting
  RooArgList fYields;//species yields
  RooArgList fPDFs;//species pdfs
  RooArgList fParameters;//model parameters
  RooArgList fConstraints;//constraints on  parameters
  RooRealVar *fID=nullptr;
  RooAbsPdf* fModel=nullptr; //model to be fitted to data
  RooAbsData* fData=nullptr; //dataset to be fitted
  TList* fCanvases=nullptr;  //canvases for plotting fitted variables
  TList* fHists=nullptr;  //histograms for plotting weighted variables
  TList* fRooFits=nullptr;  //collection of sub fits
  TTree* fTree=nullptr; //!     // tree to be imported, don't save as part of class!
  TTree* fMCIntTree=nullptr;//! //tree for calculating MC based NormInt
  TTree* fMCGenTree=nullptr;//! //tree for calculating MC based NormInt
  TString fIDBranchName="DontHaveOneYet";
  TString fYld="Yld_";     //yield variable prepend
  TString fSingleSp;
  TString fOutDir;
  TString fBinDir;
  TString fWeightName; //Input Weight species for this dataset
  TString fStudyPDF; //pdf to be studied
  TString fCut; //Apply cut to data tree and RooHSEventsPDF
  
  RooFitResult* fResult=nullptr;   //RooFit result
  THSWeights* fInWeights=nullptr; //! //input weights for dataset to be fitted
  THSBins* fDataBins=nullptr; //! 
  Bool_t fGotID=kFALSE;
  Bool_t fIsPlot=kTRUE;
  Bool_t fBinnedFit=kFALSE;
  Bool_t ftoWS=kFALSE;
  Bool_t fStudyPlot=kFALSE;
  Bool_t fOwnWorkSpace=kFALSE;
  Double_t fChi2=0;
  Int_t fFiti=0;
  Int_t fNStudyTrials=0;
  
  RooLinkedList fFitOptions;
public:
  THSRooFit() ;//default constructor, must not allocate memory!!!
  THSRooFit(TString name);
  THSRooFit(TString name,RooWorkspace* ws);
  THSRooFit(THSRooFit* rf) ;
  virtual ~THSRooFit();

  void LoadVariable(TString opt);   //Load the variables to be fitted
  void LoadVariable(RooRealVar* var){fVariables.add(*var);fWS->import(*var);}
  void LoadAuxVars(TString opt);   //Load the variables not to be fitted
  void LoadAuxVars(RooRealVar* var){fAuxVars.add(*var);fWS->import(*var);}
  void LoadBinVars(TString opt,Int_t nbins,Double_t min,Double_t max);
  void LoadBinVars(TString opt,Int_t nbins,Double_t* xbins);
  void LoadSpeciesPDF(TString opt,Int_t Scale0=1); //Load a contributing species
  void DefineSets();                  // define sets in the workspace
  void LoadDataSet(TTree* tree,Bool_t toWS=kFALSE);    //Import a tree into a RooDataSet
  void LoadDataSet(RooAbsData* data,Bool_t toWS=kFALSE);    //Import a tree into a RooDataSet
  void LoadMCIntTree(TTree* tree){fMCIntTree=tree;}
  void LoadMCGenTree(TTree* tree){fMCGenTree=tree;}
  void LoadWeights(TString wfile,TString wname);
  void LoadWeights(THSWeights* wts){fInWeights=wts;};
  void SetDataWeight(); //Add a weight to RooFit dataset
  void SetDataWeightFast(); //Add a weight to RooFit dataset without ID 
  void LoadWorkSpace(RooWorkspace* ws,TString rfname=""); //load a workspace without data
  void LoadWorkSpaceData(RooWorkspace* ws,TString rfname=""); //load a workspace with data
  void TotalPDF(); //sum different PDFs if different species
  RooArgList GetVariables(){return fVariables;}
  RooRealVar* GetVariable(TString name){return dynamic_cast<RooRealVar*>(fVariables.find(name));}
  RooRealVar* GetParameter(TString name){return dynamic_cast<RooRealVar*>(fParameters.find(name));}
  RooArgList GetAuxVars(){return fAuxVars;}
  RooArgList GetBinVars(){return fBinVars;}
  RooRealVar* GetAuxVar(TString name){return dynamic_cast<RooRealVar*>(fAuxVars.find(name));}
  RooRealVar* GetBinVar(TString name){return dynamic_cast<RooRealVar*>(fBinVars.find(name));}
  RooRealVar* GetVar(TString name);
  RooArgList GetYields(){return fYields;}
  RooArgList *GetPDFsp(){return &fPDFs;}
  RooArgList *GetYieldsp(){return &fYields;}
  RooArgList GetPDFs(){return fPDFs;}
  RooArgList GetParameters(){return fParameters;}
  RooArgList GetConstraints(){return fConstraints;}
  RooAbsPdf* GetModel(){return fModel;}
  RooAbsData* GetData(){return fData;}
  RooDataSet* GetDataSet(){return dynamic_cast<RooDataSet*>(fData);}
  RooDataHist* GetDataHist(){return dynamic_cast<RooDataHist*>(fData);}
  TList* GetFits(){return fRooFits;}
  RooFitResult* GetResult(){return fResult;}
  TList* GetPlots(){return fCanvases;};
  void AddVariables(RooArgList list){fVariables=list;}
  void AddAuxVars(RooArgList list){fAuxVars=list;}
  void AddParameters(RooArgList list){fParameters=list;}
  void AddYields(RooArgList list){fYields=list;}
  void AddGausConstraint(RooGaussian *pdf){if(!pdf) return; fWS->import(*((RooAbsPdf*)pdf->Clone()));fConstraints.add(*(fWS->pdf(pdf->GetName())));};
  void AddFormulaConstraint(RooFormulaVar *formu){if(!formu) return; fWS->import(*((RooFormulaVar*)formu->Clone()));fConstraints.add(*(fWS->function(formu->GetName())));};
  void SetSingleSpecies(TString ssp){fSingleSp=ssp;};
  void SetIDBranchName(TString str){//Probably just needed for sPlot but keep here
    if(str==fIDBranchName)return;
    fIDBranchName=str;
    fID=dynamic_cast<RooRealVar*>((fWS->factory(str+"[0,9.99999999999999e14]")));
    fWS->defineSet("ID",RooArgSet(*fID));
  }
   void AddCut(TString cut){if(fCut.Sizeof()>1){fCut+="&&";}fCut+=cut;};
  void SetCut(TString cut){fCut=cut;}
  TString GetCut(){return fCut;}
  void CheckRange();
  void SetParVals(RooFitResult *res);
  void SetParVals(RooArgList pars);
  void RemoveDataSet();
  void CleanWSDataSets();
  void RemoveConstraints(){fConstraints.removeAll();};
  TTree* GetTree(){return fTree;}
  RooWorkspace* GetWorkSpace(){return fWS;}
  void SetTree(TTree* tree){fTree=tree;}
  void PlotDataModel();
  void PlotModel(TString varName,TString pdfName);
  void Factory(TString opt){fWS->factory(opt);}
  void Fit(Bool_t randPar=kFALSE);
  virtual Bool_t  InitialiseFit();
  //  THSRooFit*  CreateSubFit(TNamed cut); //allow individual cuts
  //THSRooFit*  CreateSubFitBins(TNamed cut);//from predefined entrylists
  virtual THSRooFit*  CreateSubFitBins(TTree* ctree,TString rfname,Bool_t CopyTree=kTRUE);//from alredy selected tree
  void SavePlots(TString filename);
  virtual void FitMany(Int_t Nfits=1);
  virtual void FitSavedBins(Int_t Nfits,Bool_t cleanup=kTRUE);
  virtual void StudySavedBins(Int_t Nfits,Bool_t cleanup=kTRUE);
  virtual void FitBatchBin(Int_t Nfits);
  virtual void FitAndStudy(Int_t Nfits=1);
  void RandomisePars();
  void MakeBins();
  void MakeBins(TTree* tree,TString name);
  void MakeBinnedTrees(TTree* tree,TString name);
  // virtual void FitWithBins(Int_t Nfits=1);
  // virtual void PrepareForFarm();
  THSBins* GetBins(){return fDataBins;}
  void SetOutDir(TString name){fOutDir=name;fOutDir+="/";gSystem->MakeDirectory(fOutDir);}
  void SetBinDir(TString name){fBinDir=name;fBinDir+="/";}
  void ConfigureSavedBins(TString dirname,TString pdfname="Data");
  TString GetOutDir(){return fOutDir;}
  TString GetBinDir(){return fBinDir;}
  void WriteToFile(TString fname);
  void SetBinnedFit(Bool_t bf=kTRUE){fBinnedFit=bf;}
  void SetInWeights(THSWeights* wts){fInWeights=wts;}
  void SetWeightName(TString WName);
  virtual void DefaultFitOptions();
  void AddFitOption(RooCmdArg cmd){fFitOptions.Add((RooCmdArg*)cmd.Clone());}
  RooLinkedList GetFitOptions(){return fFitOptions;}
  void SetPlot(Bool_t plot=kTRUE){fIsPlot=plot;}
 
  void StudyFit();
  void SetNStudyTrials(Int_t Nt){fNStudyTrials=Nt;}
  void SetStudyPDF(TString pdfn){fStudyPDF=pdfn;}
  void SetStudyPlot(Bool_t plot=kTRUE){fStudyPlot=plot;}
 ClassDef(THSRooFit, 0)  // RooFit interface fit class, 

};


#endif //ifdef THSROOFIT
