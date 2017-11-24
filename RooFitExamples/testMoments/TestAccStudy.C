//To run : root --hsfit --RooHSSphHarMoments.cxx TestAccStudy.C
//Note you need $HSANA set to the HASPECT6 code dir
{

 //Create THSRooFit manager 
  THSRooFit* RF=new THSRooFit("AFit");
  RF->SetOutDir("outAcc3/"); //give output directory

  //Load the variables for fitting, these should be branches in you tree
  //see MakeEventsRes.C
  RF->LoadVariable("Z[-1,1]");//theta
  RF->LoadVariable("Phi[0,-3.14159,3.14159]");//phi
  //Set the max L and M for the moments
  Int_t LMAX=4;
  Int_t MMAX=2;
  //Create RooRealVars for all moments
  RooArgSet pars=RooHSSphHarMoments::CreateMoments(LMAX,MMAX);
  RF->GetWorkSpace()->import(pars);
  RF->GetWorkSpace()->defineSet("Moments",pars);

  //Construct Moments PDF for generating/Fitting with
  RooHSSphHarMoments* pdf=new RooHSSphHarMoments("YLM","YLM",*(RF->GetWorkSpace()->var("Z")) ,*(RF->GetWorkSpace()->var("Phi")),LMAX,MMAX,*RF->GetWorkSpace()->set("Moments"));
  //Add data (see MakeEventsRes.C
  TChain chain("decayAngles");
  chain.Add("accepted_res.root");
  pdf->SetEvTree(&chain,RF->GetCut());
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

  //Justr draw the PDF with given moments
  new TCanvas();
  ((TH2F*)RF->GetModel()->createHistogram("Z,Phi",50,50))->Draw("col1");
  //Set to to fit study 
  RF->SetStudyPDF("YLM");
  //Make a plot for every study fit (remove if doing lots of fits)
  RF->SetStudyPlot();
  //Set the number of trials you would like
  RF->SetNStudyTrials(10);
  //Run :
  //Generate events and fit for given number of trials
  RF->StudyFit();

  //Plot the PDF for the last fit result
  new TCanvas();
  ((TH2F*)RF->GetModel()->createHistogram("Z,Phi",50,50))->Draw("col1");

}
