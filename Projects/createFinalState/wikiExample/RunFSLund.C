//root --hsdata --hsfinal=THSep_omegap RunFSLund.C
//You need to replace XXX with your final state class name
{
  //Create FinalState
  THSep_omegap* fs=new THSep_omegap();
  //  fs->SetPermutate(); //turn on permuations
  fs->SetGenerated(); //just analyse generated branch
  //create datamanager
  THSLundReader* dm=new THSLundReader();
  dm->SetWriteGenBranch("Generated");
  dm->Init("Epto_omega_1.ld","HSParticles");
  //connect Project to HSParticles
  // pKK->SetDetParts(dm->GetParticles());//only use generated here
  fs->SetGenParts(dm->GetGenerated());
  Int_t counter=0;
  
  //create ouput tree
  TFile* outfile=new TFile("gentest.root","recreate");
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
