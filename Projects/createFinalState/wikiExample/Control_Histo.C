{
  //optional Set input file directory, could be done with --hsin or setenv HSIN
  HSin("/home/dglazier/Dropbox/hsana/Events/HASPECT6/Projects/finalstates/omega/data/hs_4_8_3/");
  //optional set ouput directory or file.
  //Directory=>1 out file for every in file
  //file=> 1 total output file
  HSout("testhist.root");
  
  //Set up the input chain to be analysed,could remove HSin() if you want
  TChain* tree=new TChain("HSParticles","datachain");
  tree->Add(HSin()+"out_gemc_32975714.farmpbs14_t-1_s1.root");

  //if using proof
  //Give proof destination file (env variables can't be read directly
  if(gProof) {gProof->AddInput(new TNamed("HSOUT",HSout().Data()));tree->SetProof();}


  tree->Process("Histo.C++");
  
  //Not actually sure if this is needed or not...
  if(gProof)gProof->ClearData(TProof::kUnregistered | TProof::kForceClear);
}
