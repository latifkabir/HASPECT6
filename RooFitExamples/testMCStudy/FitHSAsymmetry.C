//Run with 
//root --hsfit FitHSAsymmetry.C --MyNewPDF.cxx
{
  //Create THSRooFit manager 
  THSRooFit* RF=new THSRooFit("AFit");
  RF->SetOutDir("outRes10/"); //give output directory
  // RF->SetBinnedFit();
  ///////////////////////////////Load Variables
  RF->LoadVariable("Phi[-180,180]");//going to fit this
  RF->LoadVariable("Pol[0,1]");
  RF->GetWorkSpace()->var("Phi")->setBins(36);
  RF->GetWorkSpace()->var("Pol")->setBins(5);
  
  /////////////////////////////Make asymmetry model A and B are fit parameters
  RF->Factory("MyNewPDF::SigAsym( Phi,Pol,A[0.4,-1,1],B[-0.6,-1,1] )");
  RF->LoadSpeciesPDF("SigAsym");
  RooHSEventsPDF* PDF=dynamic_cast<RooHSEventsPDF*>(RF->GetWorkSpace()->pdf("SigAsym"));
  TChain *chainMC=new TChain("MyModel");
  chainMC->Add("SimDataSymRes10.root");
  if(!PDF->SetEvTree(chainMC,RF->GetCut())) exit(0);
  //PDF->SetNInt(5E4);
  PDF->CheckIntegralParDep(4);
  //  PDF->SetConstInt();
  RF->TotalPDF();
  TChain *chainData=new TChain("MyModel");
  chainData->Add("DataSymRes10.root");
  RF->LoadDataSet(chainData);
  gBenchmark->Start("Binned");
  RF->SetStudyPDF("SigAsym");
  RF->SetNStudyTrials(5);
  PDF->SetUseWeightsGen(kFALSE);
  //RF->SetStudyPlot();
  // RF->FitMany(1); //Just fit no study 1=>will only attempt 1 fit, if>1 will randomise paramters and try fitting however many times specified
  RF->FitAndStudy(1);
  gBenchmark->Stop("Binned");
  gBenchmark->Print("Binned");
  RF->GetWorkSpace()->saveSnapshot("fitted",RF->GetParameters(),kTRUE);
}
