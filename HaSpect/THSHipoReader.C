#include <sstream>
#include <string>
#include "THSHipoReader.h"

THSHipoReader::THSHipoReader(){
  cout<<"CONFIGURED TO RUN WITh EVENTBUILDER FROM COATJAVA 4a8.3"<<endl;
  cout<<"This specifies the  DetectorType numbering scheme"<<endl;
  cout<<"If the hipo file was created by a different version"<<endl;
  cout<<"The detector information may not be correct"<<endl;
  fHipo = new THipo();
  fHipo->ConfigBank("REC::Particle");
  fHipo->ConfigBank("REC::Scintillator");
  fHipo->ConfigBank("REC::ForwardTagger");
  fHipo->ConfigBank("REC::Event");
  //fHipo->ConfigBank("CVTRec::Tracks");
  //Get the necessary items from Particle Bank
  fPBank=fHipo->GetBank("REC::Particle");
  fPx=fPBank->GetItem("px");
  fPy=fPBank->GetItem("py");
  fPz=fPBank->GetItem("pz");
  fVx=fPBank->GetItem("vx");
  fVy=fPBank->GetItem("vy");
  fVz=fPBank->GetItem("vz");
  fPid=fPBank->GetItem("pid");
  fBeta=fPBank->GetItem("beta");
  fCharge=fPBank->GetItem("charge");
  //Get the necessary items from REC::Scintillator Bank
  fSBank=fHipo->GetBank("REC::Scintillator");
  fSPindex=fSBank->GetItem("pindex");
  fSTime=fSBank->GetItem("time");
  fSSector=fSBank->GetItem("sector");
  fSEnergy=fSBank->GetItem("energy");
  fSPath=fSBank->GetItem("path");
  fSDet=fSBank->GetItem("detector");
  //Get the necessary items from REC::Calorimeter Bank
  fCalBank=fHipo->GetBank("REC::Calorimeter");
  fCalPindex=fCalBank->GetItem("pindex");
  fCalEnergy=fCalBank->GetItem("energy");
 //Get the necessary items from REC::Cherenkov Bank
  fChBank=fHipo->GetBank("REC::Cherenkov");
  fChPindex=fChBank->GetItem("pindex");
  fChEnergy=fChBank->GetItem("nphe");
  //Get the necessary items from REC::FT Bank
  fFTBank=fHipo->GetBank("REC::ForwardTagger");
  fFTPindex=fFTBank->GetItem("pindex");
  fFTTime=fFTBank->GetItem("time");
  fFTEnergy=fFTBank->GetItem("energy");
  fFTPath=fFTBank->GetItem("path");
  fFTDet=fFTBank->GetItem("detector");

  fEvBank=fHipo->GetBank("REC::Event");
  fEvTime=fEvBank->GetItem("STTime");
 

  
  //Make vectors to hold pointers to THSParticles
  //These are just used for speed, new THSParticles
  //will be added when required by that event.
  //In the end this vector will contain max number of particles
  //from all the events
  fReadParticles=new vector<THSParticle*>;
  
 
}
Bool_t THSHipoReader::Init(TString filename,TString name){
  if(fAddGenerated&&!fMCBank){
    fHipo->ConfigBank("MC::Particle");
    fMCBank=fHipo->GetBank("MC::Particle");
    fMCPx=fMCBank->GetItem("px");
    fMCPy=fMCBank->GetItem("py");
    fMCPz=fMCBank->GetItem("pz");
    fMCVx=fMCBank->GetItem("vx");
    fMCVy=fMCBank->GetItem("vy");
    fMCVz=fMCBank->GetItem("vz");
    fMCPid=fMCBank->GetItem("pid");
    fReadGenerated=new vector<THSParticle*>;
  }
  fHipo->InitFile(filename);
  Info("THSHipoReader","Opened file");
  return kTRUE;
}
void THSHipoReader::CloseReadTree(){
  //noting to do for hipo files (I think)
}
Bool_t THSHipoReader::ReadEvent(Long64_t entry){
  //return false at end of file
  // LUND format:
// Header (Event Info):
// # of Particles, # of e, particle id, parent, daughter, p_x, p_y, p_z, E, mass, x vertex, y vertex, z vertex
// type is 1 for particles in the detector
  if(!fHipo->NextEvent()) return kFALSE;
  fEntry++;
  //  cout<<fEntry<<endl;
  fParticles.clear();//reset fParticles
  if(fAddGenerated) fGenerated.clear();//reset fGenerated
  
  if(fAddGenerated&&fMCBank){//First read generated event
    UInt_t NGen=fMCPid->Size();
    while(NGen>fReadGenerated->size()){//make sure enough space in vector
      fReadGenerated->push_back(new THSParticle());
    }
    Int_t ip=0;
    while(fMCBank->NextEntry()){
      THSParticle* particle=fReadGenerated->at(ip++);
      fGenerated.push_back(particle);
      particle->SetXYZM(fMCPx->Val(),fMCPy->Val(),fMCPz->Val(),0);
      particle->SetVertex(fMCVx->Val(),fMCVy->Val(),fMCVz->Val());
      particle->SetPDGcode(fMCPid->Val());
      particle->TakePDGMass();
    }

  }
  //eventbuilder particle bank
  if(fEvBank) fEvBank->NextEntry();
  if(fPBank){ //Filling reconstructed fParticles
    UInt_t Nin=fPid->Size();
    //in case some events have more particles
    while(Nin>fReadParticles->size()){
      fReadParticles->push_back(new THSParticle());
    }
    Int_t ip=0;
    while(fPBank->NextEntry()){
      THSParticle* particle=fReadParticles->at(ip++);
      fParticles.push_back(particle);
      particle->SetXYZM(fPx->Val(),fPy->Val(),fPz->Val(),0);
      particle->SetVertex(fVx->Val(),fVy->Val(),fVz->Val());
      //if(fBeta->Val())particle->SetMeasMass(particle->P4p()->Rho()/fBeta->Val());
      if(!fPid->Val()) particle->SetPDGcode(fCharge->ValI()*1E6); //unknown
      else  particle->SetPDGcode(fPid->ValI());
      particle->SetDetector(100);

      //Now look for the associated detector info
      //we must match the detector pindex to the index of this particle entry
      particle->SetTime(0);
      particle->SetEdep(0);
      particle->SetPath(0);
      particle->SetDetector(0);
      
      
      while(fSPindex->FindEntry(fPBank->GetEntry())){
	//Do something if find a particular detector
	particle->SetTime(fSTime->Val()-fEvTime->Val());
	particle->SetEdep(fSEnergy->Val());
	particle->SetPath(fSPath->Val()/100);
	particle->SetDetector(100+fSSector->Val());
     }
 	
      while(fCalPindex->FindEntry(fPBank->GetEntry())){
	//Do something if find a particular detector
	particle->AddEdep(fCalEnergy->Val());
	particle->SetDetector(particle->Detector()+1000);
     }
      while(fChPindex->FindEntry(fPBank->GetEntry())){
	//Do something if find a particular detector
	particle->AddEdep(fChEnergy->Val());
	particle->SetDetector(particle->Detector()+10000);
     }
 	
      while(fFTPindex->FindEntry(fPBank->GetEntry())){
	//Do something if find a particular detector
	particle->SetTime(fFTTime->Val()-fEvTime->Val());
	particle->SetEdep(fFTEnergy->Val());
	particle->SetPath(fFTPath->Val()/100);
	particle->SetDetector(300);
	particle->SetPDGcode(fCharge->ValI()*1E6);
	
      }
      
      
    }
  }

  return kTRUE;
  
}
