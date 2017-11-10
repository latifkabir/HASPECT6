//Datareader class for converting Lund format txt files
//into THSParticles
#ifndef THSHIPOREADER_h
#define THSHIPOREADER_h

#include <iostream>
#include <fstream>

#include "THSCLAS12Particle.h"
#include "THSDataManager.h"
#include "THipo.h"

class THSHipoReader: public THSDataManager{

 public :
  THSHipoReader();
  virtual ~THSHipoReader(){if(fHipo) delete fHipo;};

  virtual Bool_t Init(TString filename,TString name="");
  virtual Bool_t ReadEvent(Long64_t entry=0);
  virtual void CloseReadTree();
  void SetUsePID(Bool_t use=kTRUE){fUsePID=use;}
  
 protected :
  
  THipo *fHipo=nullptr;
  THipoBank *fPBank=nullptr;
  THipoBank *fSBank=nullptr;
  THipoBank *fMCBank=nullptr;
  THipoBank *fFTBank=nullptr;
  THipoBank *fEvBank=nullptr;
  //  THipoBank *fCVTBank=nullptr;

  THipoItem* fPid=nullptr;
  THipoItem* fPx=nullptr;
  THipoItem* fPy=nullptr;
  THipoItem* fPz=nullptr;
  THipoItem* fVx=nullptr;
  THipoItem* fVy=nullptr;
  THipoItem* fVz=nullptr;
  THipoItem* fBeta=nullptr;
  THipoItem* fCharge=nullptr;
  
  THipoItem* fSPindex=nullptr;
  THipoItem* fSTime=nullptr;
  THipoItem* fSEnergy=nullptr;
  THipoItem* fSDet=nullptr;
  THipoItem* fSPath=nullptr;

  THipoItem* fFTPindex=nullptr;
  THipoItem* fFTTime=nullptr;
  THipoItem* fFTEnergy=nullptr;
  THipoItem* fFTDet=nullptr;
  THipoItem* fFTPath=nullptr;

  THipoItem* fEvTime=nullptr;
  /* THipoItem* fFTPid=nullptr; */
  /* THipoItem* fFTTime=nullptr; */
  /* THipoItem* fFTPx=nullptr; */
  /* THipoItem* fFTPy=nullptr; */
  /* THipoItem* fFTPz=nullptr; */
  /* THipoItem* fFTCharge=nullptr; */
  /* THipoItem* fFTEnergy=nullptr; */
  
  /* THipoItem* fCVTPx=nullptr; */
  /* THipoItem* fCVTPy=nullptr; */
  /* THipoItem* fCVTPz=nullptr; */
  /* THipoItem* fCVTCharge=nullptr; */
  /* THipoItem* fCVTP=nullptr; */

  THipoItem* fMCPid=nullptr;
  THipoItem* fMCPx=nullptr;
  THipoItem* fMCPy=nullptr;
  THipoItem* fMCPz=nullptr;
  THipoItem* fMCVx=nullptr;
  THipoItem* fMCVy=nullptr;
  THipoItem* fMCVz=nullptr;

  Bool_t fUsePID=kTRUE;
  
 public :
  
 
};

#endif //ifdef THSLUNDREADER
