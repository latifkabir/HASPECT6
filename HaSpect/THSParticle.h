//--Author      DI Glazier 30/06/2014
//--Rev
//--Update
//--Description
//HASPECT Event Reconstruction
//THSParticle
//Persistant Data structure
//Contains reconstructed information required for data analysis

#ifndef __THSParticle_h__
#define __THSParticle_h__

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TDatabasePDG.h"
#include "TMatrixD.h"
#include <vector>
#include <iostream>
using namespace std;

class THSParticle {
 private:
 protected:
  TLorentzVector fP4;  //4-vector 
  TVector3 fVertex;     //particle vertex position
  Int_t fPDGCode;           //PDG number
  /* Double_t fPDGMass; */
  /* Double_t fMeasMass; //Or other PID info */
  /* Double_t fTime; */
  /* Double_t fPath; */
  /* Double_t fDoca; */
  /* Double_t fEdep; */
  Float_t fPDGMass;
  Float_t fMeasMass; //Or other PID info
  Float_t fTime;
  Float_t fPath;
  Float_t fDoca;
  vector<Float_t> fEdep;
  Int_t fDetector=0; //detector code

  vector<Float_t> fCovaMatEntries;
  TMatrixD fCovarianceMatrix;//!

  TLorentzVector fTruthP4;// true generated 4-vector
  TVector3 fTruthV;// true generated vertex
  Int_t fTruthPDG;// true PDG code
  
  Bool_t fTruthOnly=kFALSE;
  
 public:
  THSParticle();  	        //Constructor
  THSParticle(int code);  	        //Constructor
  THSParticle(TString pdgname);  	        //Constructor
  virtual ~THSParticle(){};	     	//Destructor


  //Setting functions
  void SetPDGcode(Int_t code){
    if(TDatabasePDG::Instance()->GetParticle(fPDGCode=code))
      fPDGMass = TDatabasePDG::Instance()->GetParticle(fPDGCode=code)->Mass();
    else fPDGMass=0;
  }
  void SetP4(TLorentzVector v){fP4=v;}
  void SetVectPDG(TLorentzVector v){fP4.SetVectM(v.Vect(),fPDGMass);}
  void SetP4(TLorentzVector *v){fP4=*v;}
  void SetXYZT(Double_t X,Double_t Y,Double_t Z,Double_t T){fP4.SetXYZT(X,Y,Z,T);}
  void SetXYZM(Double_t X,Double_t Y,Double_t Z,Double_t M){fP4.SetXYZM(X,Y,Z,M);}
  void SetVertex(TVector3 v){fVertex=v;}
  void SetVertex(Double_t X,Double_t Y,Double_t Z){fVertex.SetXYZ(X,Y,Z);}
  //void SetPol(TVector3 p){fPol=p;}
  //void SetPol(Double_t X,Double_t Y,Double_t Z){fPol.SetXYZ(X,Y,Z);}
  void SetTime(Double_t time){fTime=time;};
  void SetPath(Double_t path){fPath=path;};
  void SetDoca(Double_t doca){fDoca=doca;};
  void SetEdep(Double_t edep){fEdep.clear();fEdep.push_back(edep);};
  void AddEdep(Double_t edep){fEdep.push_back(edep);};
  void SetDetector(Int_t det){fDetector=det;};
  void SetMeasMass(Double_t mass){fMeasMass=mass;};
  void TakePDGMass(){SetVectPDG(fP4);}; //Preserves momentum
  void TakePDGMassFromE(){Double_t rho=sqrt(fP4.E()*fP4.E()-fPDGMass*fPDGMass);fP4.SetRho(rho);}; //preserves energy
  // void CreateTruth(){fTruth=new THSParticle();};
  void SetTruth(THSParticle* part){fTruthP4=part->P4();fTruthV=part->Vertex();fTruthPDG=part->PDG();};
  void SetTruth(THSParticle part){fTruthP4=part.P4();fTruthV=part.Vertex();fTruthPDG=part.PDG();};
  void SetTruth(TLorentzVector part,TVector3 ver,Int_t pdg){fTruthP4=part;fTruthV=ver;fTruthPDG=pdg;};
  void SetTruthOnly(Bool_t tr=kTRUE){fTruthOnly=tr;}
  //Getting functions
  TLorentzVector P4(){return fP4;}
  TLorentzVector* P4p(){return &fP4;}
  TVector3 Vertex(){return fVertex;}
  //TVector3 Pol(){return fPol;}
  Int_t PDG(){return fPDGCode;}
  Double_t PDGMass(){return fPDGMass;}
  Double_t MeasMass(){return fMeasMass;}
  Double_t Time(){return fTime;}
  Double_t MassDiff(){return fPDGMass-fMeasMass;}
  vector<Float_t> AllEdep(){return fEdep;}
  Double_t Edep(){return fEdep[0];}
  Double_t Edep(Int_t ei){return fEdep[ei];}
  Double_t Doca(){return fDoca;}
  Double_t Path(){return fPath;}
  Double_t Beta(){return fPath/fTime/2.99792e+08*1E9;}//time ns, path m
  Double_t HypBeta(){Double_t pp=fP4.Rho();return pp/sqrt(pp*pp+fPDGMass*fPDGMass);}
  Double_t HypTime(){return fPath/HypBeta()/2.99792e+08*1E9  ;} //in ns
  Double_t DeltaTime(){return HypTime()-fTime;};
  Double_t DeltaTimeVer(){return DeltaTime()+fVertex.Z()/2.99792e+08*1E9;}
  Int_t Charge();
  Int_t Detector(){return fDetector;}
  
  TLorentzVector* TruthP4p(){return &fTruthP4;};
  TLorentzVector TruthP4(){return fTruthP4;};
  TVector3* TruthVer(){return &fTruthV;};
  Int_t TruthPDG(){return fPDGCode;};
  
  Bool_t TruthOnly(){return fTruthOnly;}

  
  void CopyParticle(THSParticle* part,Bool_t andPDG);
  //Utility functions
  virtual void Print(Option_t *option="") const;

  //Add 4-vectors, doca vertices, fix pdg(optional)
  void Add(THSParticle *hsp1, THSParticle *hsp2,Int_t pdg=0);

  //DOCA routines
  Double_t MakeVirtualVertex(THSParticle* p1,THSParticle *p2);
  Double_t Calc_dtfDOCA( TVector3 locVertex1, TVector3 locUnitDir1, TVector3 locVertex2, TVector3 locUnitDir2, TVector3 *result);
  void Calc_dtfInterDOCA( TVector3 locUnitDir1, TVector3 locUnitDir2, TVector3 locVertex1,  TVector3 locVertex2, TVector3 *locInterDOCA1, TVector3 *locInterDOCA2);

 public:
  //need to permutate class
  //order in vector based on particle momentum
  //nb comparitive operator cannot work on pointers so vectors need
  //to be filled with object not pointers for this to work
  friend bool operator< ( const THSParticle& lhs, const THSParticle& rhs ){return lhs.fP4.Rho() < rhs.fP4.Rho(); };

  Double_t p3Distance(TVector3 vec){return (fP4.Vect()-vec).Mag();}

  Double_t ResTheta(){return fP4.Theta()-fTruthP4.Theta();};
  Double_t ResPhi(){return fP4.Phi()-fTruthP4.Phi();};
  Double_t ResRho(){return fP4.Rho()-fTruthP4.Rho();};
  Double_t ResE(){return fP4.E()-fTruthP4.E();};
  ClassDef(THSParticle,2) //class THSParticle
};
inline Int_t THSParticle::Charge(){
  
  if(fPDGCode==1E6) return 1;
  else if(fPDGCode==-1E6) return -1;
  else {
    Int_t charge=TDatabasePDG::Instance()->GetParticle(fPDGCode)->Charge();
    if(charge!=0)charge=(Int_t) charge/TMath::Abs(charge); //just get sign not mag.
    return charge;
  }
}

//inline bool THSParticle::operator<( const THSParticle& rhs ) {cout<<" "<<rhs.fP4.Rho()<<endl;return fP4.Rho() < rhs.fP4.Rho(); }
//inline bool THSParticle::operator<( const THSParticle& rhs ) {cout<<" "<<rhs.fP4.Rho()<<endl;return this->fP4.Rho() < rhs.fP4.Rho(); }
//inline bool THSParticle::IsFaster( const THSParticle& lhs,const THSParticle& rhs ) {cout<<" "<<rhs.fP4.Rho()<<endl;return lhs.fP4.Rho() < rhs.fP4.Rho(); }
#endif
