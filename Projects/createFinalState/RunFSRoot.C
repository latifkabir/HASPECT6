//root --hsdata --hsfinal=THSXXX RunFSLund.C
//You need to replace XXX with your final state class name
{
  //Create FinalState
  THSXXX* fs=new THSXXX();
  fs->SetPermutate(); //turn on permuations
  // fs->SetGenerated(); //just analyse generated branch
  //create datamanager
  THSDataManager* dm=new THSDataManager();
  dm->SetReadGenBranch("Generated");
  dm->Init("INPUT.root","HSParticles");
  //connect Project to HSParticles
  fs->SetDetParts(dm->GetParticles());
  fs->SetGenParts(dm->GetGenerated());
  Int_t counter=0;
  
  //create ouput tree
  TFile* outfile=new TFile("OUTPUT.root","recreate");
  TTree* outtree=new TTree("FinalTree","output tree");
  fs->FinalStateOutTree(outtree); //connect ouput tree to project branches
  
  gBenchmark->Start("timer");
  
  while(dm->ReadEvent()){//loop over events
    fs->ProcessEvent();
  }
  gBenchmark->Stop("timer");
  gBenchmark->Print("timer");
  
  outfile->cd();
  outtree->Write();
  delete outfile;
}
