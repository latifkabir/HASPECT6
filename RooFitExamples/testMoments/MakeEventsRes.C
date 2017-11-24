//headers for compilation
#include <TTree.h>
#include <TFile.h>
#include <TRandom.h>
#include <TMath.h>

//Declare variables globally so they can be accessed in all functions
Double_t Z,Phi,UID,E=0;
Double_t genZ,genPhi=0;

//Got to declare acceptance functions here
Bool_t DetectorSim0();
//This is the main function which will be run
//The arguement N= number of generated (not accepeted) events
void MakeEventsRes(Long64_t N=1E6){

  //Make phase space (no acceptance or resolution) tree and file
  TFile *file=new TFile("phasespace_res.root","recreate");
  TTree* tree=new TTree("decayAngles","decayAngles");
  tree->Branch("UID",&UID,"UID/D");
  tree->Branch("Z",&Z,"Z/D");
  tree->Branch("E",&E,"E/D");
  tree->Branch("Phi",&Phi,"Phi/D");

  //Loop over N and generate random flat Z and phi
  for(Int_t i=0;i<N;i++){
    UID=i;
    Z=gRandom->Uniform(-1,1);
    E=gRandom->Uniform(0,10);
    Phi=gRandom->Uniform(-1,1)*TMath::Pi();
    tree->Fill();
  }
  //save tree
  tree->Write();
  file->Close();

  //Make accepted events tree (with resolution and accetpance)
  TFile *filea=new TFile("accepted_res.root","recreate");
  TTree* treea=new TTree("decayAngles","decayAngles");
  
  treea->Branch("UID",&UID,"UID/D");
  treea->Branch("Z",&Z,"Z/D");
  treea->Branch("E",&E,"E/D");
  treea->Branch("Phi",&Phi,"Phi/D");
  treea->Branch("genZ",&genZ,"genZ/D");
  treea->Branch("genPhi",&genPhi,"genPhi/D");

  //Loop over N, generate a random flat genZ and genPhi
  //Pass these to detector simulation (i.e. add resolution and acceptance)
  for(Int_t i=0;i<N;i++){
    UID=i;
    E=gRandom->Uniform(0,10);
    genZ=gRandom->Uniform(-1,1);
    genPhi=gRandom->Uniform(-1,1)*TMath::Pi();
    if(DetectorSim0())treea->Fill(); 
  }
  treea->Write();
  filea->Close();
}

Bool_t DetectorSim0(){
  //Function to mimic detector accetpance and resolutions
  //"Spherical" Acceptance function 
  if(genZ*genZ+genPhi*genPhi/TMath::Pi()/TMath::Pi()<2*gRandom->Uniform(0,1)){
    //detected event, now make reconstructed Z and Phi
    //i.e. Add resolution smearing
    Z=genZ+gRandom->Gaus(0,0.1); //width cos theta 0.1 = 5% or range
    if(TMath::Abs(Z)>1) return kFALSE;//outwith accetpance, there might be more realistic ways than this!!
    Phi=genPhi+gRandom->Gaus(0,2*TMath::DegToRad()); //width 2deg
    //Keep -pi<Phi<pi 
    if(Phi>TMath::Pi()) Phi= -2*TMath::Pi()+Phi;
    if(Phi<-TMath::Pi()) Phi= 2*TMath::Pi()+Phi;
    return kTRUE; 
  }
  return kFALSE;
}
