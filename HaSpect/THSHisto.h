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

#ifndef THSHisto_h
#define THSHisto_h

#include <TSelector.h>
#include <TH2.h>
#include <TKey.h>
#include <TString.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TEntryList.h>
#include <iostream>
#include <map>
#include "THSWeights.h"
#include "THSBins.h"
#include "THSParticle.h"

using namespace std;

typedef map< TString , TH1* > StrTH1_Map;  //make a type connecting strings to TH1
typedef vector< StrTH1_Map > Vec_StrTH1; //make a type which is a vector of StrTH1_Maps

class THSHisto  {
 public :
  
 THSHisto(TTree * /*tree*/ =0) : fBinEntryLists(0) { }
  virtual ~THSHisto() {if(fWeights)SafeDelete(fWeights);}

  THSWeights* fWeights=nullptr;
  Double_t fWeight=1;
  THSBins* fBins=nullptr;
  virtual void InitialiseHisto(TList *input);
  virtual void LoadCut(TString sCut);
  virtual void LoadHistograms();
  virtual void ChangeNames(/*TFile* list*/);
  virtual TH1* MapHist(TH1* hist);
  virtual TH1* FindHist(TString hname);
  virtual void OrganiseSavedHists(TFile* stdFile,TList* output);
  void MakeAxisDir(Int_t iA,TDirectory* curDir);

  
  virtual void HistogramList(TString sLabel){};// =0; 
  virtual void FillHistograms(TString sCut){};// =0;

  //Functions for THSWeights usage
  void LoadWeights(TString fname,TString wname);
  void LoadCutsForWeights();
  void FillCutsForWeights();
  void GetWeightEvent(Double_t id);
  void SetWeight(TString species);

  //Functions for histogramming THSPArticles
  virtual void ParticleList(TString name);
  virtual void FillParticles(TString sCut,THSParticle* part);

 protected:
  TList* fSelInput=nullptr;  //pointers to the selector input
  TList* fSelOutput=nullptr;  //pointer to the selector ouput
  Bool_t fIsHisto=0;
  StrTH1_Map fHistNameMap; //connect hist names to TH1* objects
  Vec_StrTH1 fVecHistCut; //a vector with elements for each cut
  vector< Vec_StrTH1 > fVecHistBin; //a vector with elements for each bin
  
  TString fCurrCut; //For passing cut names between funcions
  Int_t fCurrBin = 0; //for passing bin number
  vector< TString > fVecTypes; // name of hist without cut or bin
  vector< TString > fVecCuts; // list of defined cuts
  vector< TString > fVecBinNames; // list of histogram bin names
  map< TString, Int_t > fMapCuts; // map of defined cuts to fVecHistCut element
  map< TString,TString > fBinMap; //bin (kinmematic) name to number map

  TObjArray* fBinEntryLists; //arrray of entry lists for each kinematic bin
   
};


#endif 
