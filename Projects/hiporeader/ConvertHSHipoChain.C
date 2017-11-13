//Run with root --hsdata ConvertHSHipoChain.C --hsin=/path/to/hipo/files --hsout=/path/to/output/root/files
{
  THSHipoReader* reader=new THSHipoReader();
  reader->SetWriteGenBranch("Generated");

  TChain chain("hipo");
  chain.Add(HSin()+"/*.hipo");

  //Create output directory and make sure it is OK to continue of it
  //already exists
  if(gSystem->MakeDirectory(HSout().Data())==-1){
    cout<<"Warning directory "<<HSout()<<" already exists, is it OK to continue? y/n"<<endl;
    TString ans;
    cin>>ans;
    if(ans!=TString("y"))
      gROOT->ProcessLine(".q"); 
  }
  gBenchmark->Start("time");

  reader->ReadWriteChain(&chain,HSout(),".root");
  
  cout<<"Analysed "<<reader->GetEntry()<<" events"<<endl;
  
  gBenchmark->Stop("time");
  gBenchmark->Print("time");

 }
