#ifndef THSSPLOT_h
#define THSSPLOT_h
#include "THSRooFit.h"
#include <RooStats/SPlot.h>
using namespace RooFit;
using namespace std;

class THSsPlot : public THSRooFit {
 protected:
  RooStats::SPlot* fSPlot=nullptr; //!  //sPlot object
  THSWeights* fWeights=nullptr; //! new weights calcualted in this fit
  Double_t fSRange[2];
  Bool_t fSaveWeights=kTRUE;
public:
  THSsPlot() ;//default constructor, must not allocate memory!!!
  THSsPlot(TString name);
  THSsPlot(TString name,RooWorkspace *ws);
  THSsPlot(THSsPlot* rf) ;
  virtual ~THSsPlot();


  THSWeights* GetWeights(){return fWeights;}
  RooWorkspace* GetWorkSpace(){return fWS;}
  //virtual Bool_t  InitialiseFit();
  //sPlot functions
  void sPlot();
  void SetSPlotRange(Double_t min,Double_t max){fSRange[0]=min,fSRange[1]=max;};
  void MergeModelSpecies();
  void SetSingleSpecies(TString ssp){fSingleSp=ssp;};
  void SetIDBranchName(TString str){
    fIDBranchName=str;
    fID=dynamic_cast<RooRealVar*>((fWS->factory(str+"[0,9.99999999999999e14]")));
    fWS->defineSet("ID",RooArgSet(*fID));
   }
  Double_t GetSWeight(Long64_t event,TString species);
  void ExportWeights(TString wmname="WeightMap");
  void ExportWeights1(TString wmname="WeightMap");
  void ExportTreeWeights(TString wmname="WeightMap");
  void ExportWeightsToFile(TString filename);
  void AddWeightMap(THSWeights* Wts);
  void DrawTreeVar(TString VarName,Int_t nbins,Double_t hmin,Double_t hmax);
// void CorrelationWithVars(TString VarName);
  void AddSubWeights();
  virtual THSRooFit* CreateSubFitBins(TTree* ctree,TString rfname,Bool_t CopyTree);
  void SaveHists(TString filename);
  // virtual void RunWeights(Int_t Nfits=1);
  // virtual void RunSingleWeights(Int_t Nfits);
  void SetWeights(THSWeights* wts){fWeights=wts;}
  void DeleteWeights(){if(fWeights)delete fWeights;fWeights=nullptr;};
  virtual void FitAndStudy(Int_t Nfits);
  virtual void FitSavedBins(Int_t Nfits);
  virtual void DefaultFitOptions();
  void SetSaveWeights(Bool_t save=kTRUE){fSaveWeights=save;}
  
  ClassDef(THSsPlot, 1)  // RooFit interface fit class, 

};


#endif //ifdef THSROOFIT
