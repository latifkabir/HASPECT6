//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Nov 17 16:58:08 2017 by ROOT version 6.08/02
// from TTree HSParticles/data tree
// found on file: /home/dglazier/Dropbox/hsana/Events/HASPECT6/Projects/finalstates/omega/data/hs_4_8_3/out_gemc_32975714.farmpbs14_t-1_s1.root
//////////////////////////////////////////////////////////

#ifndef Histo_h
#define Histo_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include "THSOutput.h"
#include "THSep_omegap.h"

// Header file for the classes stored in the TTree if any.

class Histo : public TSelector, public THSOutput, public THSep_omegap {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

  //data members for new branches
   //you must define how they are processed for each event and set to 0 in constructor
   //e.g.   TLorentzVector  *fp1;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<THSParticle*> *Particles;
   vector<THSParticle*> *Generated;

   // List of branches
   TBranch        *b_Particles;   //!
   TBranch        *b_Generated;   //!

   Histo(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~Histo() {/*//SafeDelete(fp1);//If we add pointer ro classmust delete it*/  }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
  //Additional final state analysis function
  virtual void UserProcess();
   virtual void    Terminate();
//Add THSHisto functions
   virtual void HistogramList(TString sLabel);
   virtual void FillHistograms(TString sCut);

   ClassDef(Histo,0);
};

#endif

#ifdef Histo_cxx
void Histo::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Particles = 0;
   Generated = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Particles", &Particles, &b_Particles);
 THSFinalState::SetDetParts(Particles);
   fChain->SetBranchAddress("Generated", &Generated, &b_Generated);
 if(fChain->GetBranch("Generated"))THSFinalState::SetGenParts(Generated);
}

Bool_t Histo::Notify()
{
   THSOutput::HSNotify(fChain);
  THSOutput::InitOutTree();
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef Histo_cxx
