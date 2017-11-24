//To run : root --hsfit --RooHSSphHarMoments.cxx TestAccFitBins.C
//Note you need $HSANA set to the HASPECT6 code dir
{

 //Create THSRooFit manager 
  THSRooFit* RF=new THSRooFit("AFit");
  RF->SetOutDir("outAcc3/"); //give output directory

  //Load the variables for fitting, these should be branches in you tree
  //see MakeEventsRes.C
  RF->LoadVariable("Z[-1,1]");//theta
  RF->LoadVariable("Phi[0,-3.14159,3.14159]");//phi
  RF->LoadBinVars("E",5,0,10);//should be same name as variable in tree

  TChain chain("decayAngles");
  chain.Add("accepted_res.root");
  RF->SetIDBranchName("UID");
  RF->MakeBinnedTrees(&chain,"YLM");

  //Here I need to specify the string "YLM" so it oes not look for data by default
  //This just comes from the previous line and was used to tag the split data files
  RF->ConfigureSavedBins(RF->GetOutDir(),"YLM"); 

  // RF->SetCut("TMath::Abs(Z)<0.5");//Make extra cut, variable must have been loaded via LoadAuxVar, if it is not a fit variable
  
  //Set the max L and M for the moments
  Int_t LMAX=4;
  Int_t MMAX=2;
  //Create RooRealVars for all moments
  RooArgSet pars=RooHSSphHarMoments::CreateMoments(LMAX,MMAX);
  RF->GetWorkSpace()->import(pars);
  RF->GetWorkSpace()->defineSet("Moments",pars);

  //Construct Moments PDF for generating/Fitting with
  RooHSSphHarMoments* pdf=new RooHSSphHarMoments("YLM","YLM",*(RF->GetWorkSpace()->var("Z")) ,*(RF->GetWorkSpace()->var("Phi")),LMAX,MMAX,*RF->GetWorkSpace()->set("Moments"));
  //Add data (see MakeEventsRes.C)
  //  pdf->SetEvTree(&chain);
  pdf->SetNInt(10000);//Number of events to use in integration calc.
  pdf->SetUseWeightsGen(kFALSE); //Use accept/reject not weights
  RF->GetWorkSpace()->import(*pdf); //import pdf into workspace

  //Load PDF into RooFit Manager and set to generate 1000 events
  RF->LoadSpeciesPDF("YLM",1000);
  RF->TotalPDF(); //Make extended PDF
  //Choose generated moments
  //Warning some combinations of moments can give negative PDFs
  //It would be good if we could check for this before running?
  RF->GetWorkSpace()->var("MomY_1_1")->setVal(0.3);
  RF->GetWorkSpace()->var("MomY_2_0")->setVal(-0.2);
  RF->GetWorkSpace()->var("MomY_4_2")->setVal(0.1);
  RF->GetWorkSpace()->var("MomY_4_1")->setVal(-0.2);

  //Set to to fit study 
  RF->SetStudyPDF("YLM");
  //Make a plot for every study fit (remove if doing lots of fits)
  RF->SetStudyPlot();
  //Set the number of trials you would like
  RF->SetNStudyTrials(1);

  //Try 1 fit with different starting paramters
  RF->StudySavedBins(1);

 
}
