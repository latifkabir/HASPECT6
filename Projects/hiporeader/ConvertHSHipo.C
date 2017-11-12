//Run with root --hsdata ConvertHSHipo.C --hsin=out_out.hipo --hsout=file.root
{
  THSHipoReader* reader=new THSHipoReader();
  reader->SetWriteGenBranch("Generated");
  
  gBenchmark->Start("time");

  reader->Init(HSin());
  reader->WriteParticles(HSout());

  cout<<"Analysed "<<reader->GetEntry()<<" events"<<endl;
  
  gBenchmark->Stop("time");
  gBenchmark->Print("time");

  // //loop through tree and add a UID, 
  // reader->SetUID(100); //may want to apply offset
  // reader->AddUID(HSout());
}
