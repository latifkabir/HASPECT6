//Run with 
//root --hsfit FitHSFactorisedMmiss.C

{
  THSsPlot* RF=new THSsPlot("SF");
  RF->SetOutDir("outMmiss/");
  ///////////////////////////////Load Variables
  RF->LoadVariable("Mmiss[0.,10]");//should be same name as variable in tree 
  RF->LoadAuxVars("M1[0,10]");//should be same name as variable in tree 
  RF->LoadAuxVars("Eg[3,4]");//Not to be fitted but limits applied to dataset
  RooRealVar * var=RF->GetWorkSpace()->var("Mmiss");
  var->setBins(100); //number of bins used in PDF histogram 

  //////////////////////////////Make signal PDF
  RF->Factory("RooHSEventsHistPDF::Sig(Mmiss,alpha[0,0,20],off[0,-2,2],scale[1,0.8,1.2])");
  RooHSEventsHistPDF* pdf=dynamic_cast<RooHSEventsHistPDF*>(RF->GetWorkSpace()->pdf("Sig"));
  //Attach MC signal data
  TChain *chainmc=new TChain("MyModel","mcsignal");
  chainmc->Add("SigData.root");
  chainmc->Add("BG1Data.root");
  pdf->SetEvTree(chainmc,RF->GetCut());
  RF->LoadSpeciesPDF("Sig",1); 
  //////////////////////////////Make background PDF
  // RF->Factory("RooHSEventsHistPDF::BG1(Mmiss,alphaB[0,0,5],offB[0,0,0],scaleB[1.0,0.8,1.2])");
  // RooHSEventsHistPDF* pdfb=dynamic_cast<RooHSEventsHistPDF*>(RF->GetWorkSpace()->pdf("BG1"));
  // //Attach MC background data
  // TChain *chainmcb=new TChain("MyModel","mcbackground");
  // chainmcb->Add("BG1DataMulti.root");
  // pdfb->SetEvTree(chainmcb);
  // RF->LoadSpeciesPDF("BG1",1);
 //////////////////////////////Make background PDF
  RF->Factory("RooHSEventsHistPDF::BG2(Mmiss,alphaB2[0,0,5],offB2[0,0,0],scaleB2[1.0,0.8,1.2])");
  RooHSEventsHistPDF* pdfb2=dynamic_cast<RooHSEventsHistPDF*>(RF->GetWorkSpace()->pdf("BG2"));
  //Attach MC background data
  TChain *chainmcb2=new TChain("MyModel","mcbackground");
  chainmcb2->Add("BG2Data.root");
  pdfb->SetEvTree(chainmcb2,RF->GetCut());
  RF->LoadSpeciesPDF("BG2",1);

  //Add data to chain
  TChain chain("MyModel");
  chain.Add("Data.root");
  RF->SetIDBranchName("fgID");
  //import to RooFit
  RF->LoadDataSet(&chain);
  gBenchmark->Start("Binned");
  //Add constraints to PDF fudge parameters
  //These are derived from the alpha off and scale parameter initial val and limits
  //i.e. intital val = gaussian mean; range = 10*sigma
  RF->AddGausConstraint(pdf->AlphaConstraint());
  RF->AddGausConstraint(pdf->OffConstraint());
  RF->AddGausConstraint(pdf->ScaleConstraint());
  RF->AddGausConstraint(pdfb->AlphaConstraint());
  RF->AddGausConstraint(pdfb->OffConstraint());
  RF->AddGausConstraint(pdfb->ScaleConstraint());
  RF->AddGausConstraint(pdfb2->AlphaConstraint());
  RF->AddGausConstraint(pdfb2->OffConstraint());
  RF->AddGausConstraint(pdfb2->ScaleConstraint());
  RF->FitAndStudy(1); //argument = number of test fits with random initial pars
  gBenchmark->Stop("Binned");
  gBenchmark->Print("Binned");
  RF->DrawTreeVar("Mmiss",100,0,10);
  RF->DrawTreeVar("M1",100,0,10);
  RF->DrawTreeVar("M2",100,0,10);

}
