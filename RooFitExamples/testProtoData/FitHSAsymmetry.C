//Run with 
//root --hsfit FitHSAsymmetry.C --SKNewPDF.cxx
{
  //Create THSRooFit manager 
  THSRooFit* RF=new THSRooFit("AFit");
  RF->SetOutDir("out/"); //give output directory
  // RF->SetBinnedFit();
  ///////////////////////////////Load Variables
  RF->LoadVariable("Phi[-180,180]");//going to fit this
  RF->LoadVariable("Pol[0,1]");
  RF->LoadVariable("PolState[Polp=1,Polm=-1]");//Load a category
  
  RF->GetWorkSpace()->var("Phi")->setBins(36);
  RF->GetWorkSpace()->var("Pol")->setBins(5);
  
  /////////////////////////////Make asymmetry model A and B are fit parameters
  RF->Factory("SKNewPDF::SigAsym( Phi,Pol,PolState, A[0.1,-1,1],B[0.6,-1,1] )");
  RooHSEventsPDF* PDF=dynamic_cast<RooHSEventsPDF*>(RF->GetWorkSpace()->pdf("SigAsym"));
  TChain *chainMC=new TChain("MyModel");
  chainMC->Add("SimASymRes10.root");
  cout<<"PDF "<<PDF<<endl;
  PDF->SetEvTree(chainMC,RF->GetCut());
  RF->LoadSpeciesPDF("SigAsym");
  //PDF->SetNInt(1E4);
  //PDF->CheckIntegralParDep(4);
  //PDF->SetConstInt();
  RF->TotalPDF();
  TChain *chainData=new TChain("MyModel");
  chainData->Add("DataASymRes10PS02.root");
  RF->LoadDataSet(chainData);
  gBenchmark->Start("Binned");
  PDF->AddProtoData(RF->GetDataSet());
  PDF->CheckIntegralParDep(10);
  RF->SetStudyPDF("SigAsym"); //study this PDF
  RF->SetNStudyTrials(2);  //Perform 20 trials
  RF->SetStudyPlot();  //Produce plots for each trial fit
  RF->FitAndStudy(1);
  // gBenchmark->Stop("Binned");
  // gBenchmark->Print("Binned");
}
