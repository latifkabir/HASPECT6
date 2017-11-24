{

  TFile *file=new TFile("phasespace.root","recreate");
  TTree* tree=new TTree("decayAngles","decayAngles");
  Double_t Z,Phi,UID=0;
  tree->Branch("UID",&UID,"UID/D");
  tree->Branch("Z",&Z,"Z/D");
  tree->Branch("Phi",&Phi,"Phi/D");

  Long64_t N=1E6;
  for(Int_t i=0;i<N;i++){
    UID=i;
    Z=gRandom->Uniform(-1,1);
    Phi=gRandom->Uniform(-1,1)*TMath::Pi();
    tree->Fill();
  }
  tree->Write();
  file->Close();

  TFile *filea=new TFile("accepted.root","recreate");
  TTree* treea=new TTree("decayAngles","decayAngles");
  
  treea->Branch("UID",&UID,"UID/D");
  treea->Branch("Z",&Z,"Z/D");
  treea->Branch("Phi",&Phi,"Phi/D");

  for(Int_t i=0;i<N;i++){
    UID=i;
    Z=gRandom->Uniform(-1,1);
    Phi=gRandom->Uniform(-1,1)*TMath::Pi();
    if(Z*Z+Phi*Phi/TMath::Pi()/TMath::Pi()<2*gRandom->Uniform(0,1))treea->Fill();
  }
  treea->Write();
  filea->Close();

  
}
