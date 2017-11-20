//root --hsdata --hsfinal=THSep_omegap RunFSRoot.C
//You need to replace ep_omegap with your final state class name
{
  //Create FinalState
  THSep_omegap* fs=new THSep_omegap();
  fs->SetPermutate(); //turn on permuations
  // fs->SetGenerated(); //just analyse generated branch
  //create datamanager
  THSDataManager* dm=new THSDataManager();
  dm->SetReadGenBranch("Generated");
  dm->Init("/home/dglazier/Dropbox/hsana/Events/HASPECT6/Projects/finalstates/omega/data/hs_4_8_3/out_gemc_32975714.farmpbs14_t-1_s1.root","HSParticles");
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
