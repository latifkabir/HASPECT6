/**
 \class THSSkeleton
 
 Class to generate selector code
 based on TTree::MakeSelector()
 
 @example MakeHSSelector.C
 
*/

#include "THSSkeleton.h"
#include <RooClassFactory.h>
#include <TSystem.h>
#include <iostream>
using namespace std;


/////////////////////////////////////////////////////////////////////
///Create intial selector code from TTree::MakeSelector()
void THSSkeleton::CreateSelector(TString selname,TString filename,TString treename,TString opt){
  fSelName=selname;
  fFileName=filename;
  fTreeName=treename;
  fOption=opt;
  
  fFile=TFile::Open(fFileName);
  if(!fFile)
    Fatal("THSSkeleton::THSkeleton(TString filename,TString treename)","No file %s found",fFileName.Data());
  fTree=(TTree*)fFile->Get(fTreeName);
  if(!fTree)
    Fatal("THSSkeleton::THSkeleton(TString filename,TString treename)","No tree %s found",fTreeName.Data());

  //check if skelton code exists already
  TString selfile=fSelName+".C";
  TString sOverwrite;
  if(gSystem->FindFile("./",selfile)){Info("THSSkeleton::THSkeleton(TString filename,TString treename)","Selector Code already exists to exit type y");cin>>sOverwrite;};
  if(sOverwrite==TString("y")) exit(0);
  //fTree->MakeSelector(fSelName);
  if(fIsFinalState)fOption+="=legacy";
  fTree->MakeSelector(fSelName,fOption);
  fMadeSelector=kTRUE;


  AddHSOutput();
}

void THSSkeleton::AddHSOutput(){
  fPlace=0;
 
  HSOut_C();
  HSOut_h();
  HSOut_Control();

  if(fIsHisto) HSHisto();
  if(fIsAppendTree) HSAppend();
  if(fIsNewTree) HSNewTree();
  if(fIsWeights) HSWeights();
  if(fNLPS) HSLPS();
  if(fIsFinalState) HSFinalState();
  // AddLineAfter(,,1);
  // AddLineAfter(,,1);
  // AddLineAfter(,,1);
  // AddLineAfter(,,1);
  // AddLineAfter(,,1);
  // AddLineAfter(,,1);

}

void THSSkeleton::HSOut_C(){
  //Open source code as macro
  fCurMacro=TMacro(fSelName+".C");
  //Note bug in TMacro, where it does not read final }character so add it here
  //May cause issues on different systems....
  TList* lines=fCurMacro.GetListOfLines();
  if(!((TObjString*)lines->Last())->String().Contains("}")){
    fPlace=fCurMacro.GetListOfLines()->GetEntries();
    ContinueLineAfter("}");
    ContinueLineAfter(" ");
    fPlace=0;
  }
  //End bug fix
  //second bug fix
  //TTreeReader should call SetEntryLocal for chains
 if(!fOption.Contains("legacy"))  ReplaceMacroText("fReader.SetEntry","fReader.SetLocalEntry");
  //
  //First add shadow Selector function calls
  AddLineAfter("::Begin(TTree * /*tree*/)","   THSOutput::HSBegin(fInput,fOutput);",1);//off=1 to get past { line
  if(fOption.Contains("legacy"))AddLineAfter("::Process(Long64_t entry)","   THSOutput::HSProcessStart(entry);",1);
  else AddLineAfter("fReader.Set","   THSOutput::HSProcessStart(entry);",0);
  AddLineAfter("::SlaveTerminate()","   THSOutput::HSSlaveTerminate();",1);
  AddLineAfter("::Terminate()","   THSOutput::HSTerminate();",1);
  //Case of SlaveBegin, we want to have TString option=GetOPtion line first
  MoveToLine("::SlaveBegin(TTree * /*tree*/)");
  FindNextLineLike("TString option = GetOption();");
  ContinueLineAfter("   THSOutput::HSSlaveBegin(fInput,fOutput);");
  MoveToLine("::Process(Long64_t entry)");
  FindNextLineLike("return kTRUE;");
  ContinueLineAfter("   THSOutput::HSProcessFill();",-1);

  if(fOption.Contains("legacy"))AddLineAfter("HSProcessStart","   GetEntry(entry);",-1);

  fCurMacro.SaveSource(fSelName+".C");
}

void THSSkeleton::HSOut_h(){
  //Open source code as macro
  fCurMacro=TMacro(fSelName+".h");

  AddLineAfter("#include <TSelector.h>","#include \"THSOutput.h\"");
  AddLineAfter("::Notify()","   THSOutput::HSNotify(fChain);",1);
  if(!fOption.Contains("legacy"))AddLineAfter("::Notify()","   fChain=fReader.GetTree();",1);
  
  // AddLineAfter(,,1);
  ReplaceMacroText("public TSelector","public TSelector, public THSOutput");
  
  fCurMacro.SaveSource(fSelName+".h");
}

void THSSkeleton::HSOut_Control(){
  fPlace=0;
  
  TString HSANA=gSystem->Getenv("HSANA");
  //Open a template file 
  fCurMacro=TMacro(HSANA+"/Control_Template.C");
  ReplaceMacroText("INPUT_FILE_DIR",TString(gSystem->DirName(fFileName))+"/");
  ReplaceMacroText("TREENAME",fTreeName);
  ReplaceMacroText("FILENAMES",gSystem->BaseName(fFileName));
  ReplaceMacroText("SELECTOR",fSelName+".C++");

  fCurMacro.SaveSource(TString("Control_")+fSelName+".C");

}

/////////////////////////////////////////////
/// Add examples of how to use the
/// THSHisto class
void THSSkeleton::HSHisto(){
  fPlace=0;
  //////////////////////////////////////////////////////////////////  
  //First deal with .C
  fCurMacro=TMacro(fSelName+".C");
  FindNextLineLike("GetOption();");
  //Define bins if wanted
  ContinueLineAfter("//If you want to split histograms into different kinematic bins, include and configure the lines below");
  ContinueLineAfter("//   fBins=new THSBins(\"HSBins\");");
  ContinueLineAfter("//   fBins->AddAxis(\"AXIS1NAME\",10,-1,1);");
  ContinueLineAfter("//   fBins->AddAxis(\"AXIS2NAME\",10,-1,1);");
  ContinueLineAfter("//   if(!fInput) fInput=new TList();");
  ContinueLineAfter("//   fInput->Add(fBins);");
  
  AddLineAfter("THSOutput::HSSlaveBegin(fInput,fOutput);","  //Initialise and configure histogramming");
  ContinueLineAfter("   THSHisto::InitialiseHisto(fInput);");
  ContinueLineAfter("   THSHisto::LoadCut(\"Cut1\");");
  ContinueLineAfter("   THSHisto::LoadHistograms();");

  //Now fill histograms in process
  AddLineAfter("THSOutput::HSProcessFill();"," //below you can give vars corresponding to fBins axis",-2);
  ContinueLineAfter("   //if(fBins) fCurrBin=fBins->FindBin(var1,var2,...);//if fBins is defined need to give this meaningful arguments");
  ContinueLineAfter("   FillHistograms(\"Cut1\");");

  //Now add new functions at end of file
  fPlace=fCurMacro.GetListOfLines()->GetEntries()+1;
  ContinueLineAfter(TString("void ")+fSelName+"::HistogramList(TString sLabel){");
  ContinueLineAfter("  TDirectory::AddDirectory(kFALSE); //do not add to current directory");
  ContinueLineAfter("  //now define all histograms and add to Output");
  ContinueLineAfter("  //label includes kinematic bin and additional cut name");
  ContinueLineAfter("  // e.g fOutput->Add(MapHist(new TH1F(\"Mp1\",\"M_{p1}\",100,0,2)));");
  ContinueLineAfter("  //end of histogram list");
  ContinueLineAfter("  //ParticleList(\"NAME\"); //Create standard histograms for a THSParticle");
  ContinueLineAfter("  TDirectory::AddDirectory(kTRUE); //back to normal");
  ContinueLineAfter("}");
  
  ContinueLineAfter(TString("void ")+fSelName+"::FillHistograms(TString sCut){");
  ContinueLineAfter("  if(fBins&&fCurrBin==-1) return; //Outwith bin limits so no hists");
  ContinueLineAfter("  fCurrCut=sCut;");
  ContinueLineAfter("  //Get histogram from list");
  ContinueLineAfter("  //Fill histogram");
  ContinueLineAfter("  // e.g. FindHist(\"Mp1\")->Fill(fp1->M());//1D");
  ContinueLineAfter("  // e.g. ((TH2D*)FindHist(\"Mp1\"))->Fill(fp1->M(),fp2->M());//2D");
  ContinueLineAfter("  //Fill THSParticle histogram with values from THSParticle PART");
  ContinueLineAfter("  //FillParticles(\"NAME\",&PART);"); 
  ContinueLineAfter("}");
// ContinueLineAfter();
  fCurMacro.SaveSource(fSelName+".C");
  /////////////////////////////////////////////////////////////////
  //Now with .h
  fPlace=0;
  fCurMacro=TMacro(fSelName+".h");
  AddLineAfter("virtual void    Terminate();","//Add THSHisto functions");
  ContinueLineAfter("   virtual void HistogramList(TString sLabel);");
  ContinueLineAfter("   virtual void FillHistograms(TString sCut);");
  
  fCurMacro.SaveSource(fSelName+".h");
}

//////////////////////////////////////////////////////////////////
/// Adds examples of how to append branches to the input
/// to the input tree.
void THSSkeleton::HSAppend(){
  if(fIsNewTree)Error("THSSkeleton::HSNewTree()","Can't have a new tree and an append tree!");

  fPlace=0;
  //////////////////////////////////////////////////////////////////  
  //First deal with .h
  fCurMacro=TMacro(fSelName+".h");
  AddLineAfter("THSOutput::HSNotify(fChain);","  //If we want to clone the input tree we have to do it here as fChain does not exist in SlaveBegin on PROOF");
  ContinueLineAfter("  if(!fOutTree){//First file, clone the input tree and add branches");
  ContinueLineAfter("   Int_t buff=32000;");
  ContinueLineAfter("   Int_t split=0;//note split is important in the TSelector framework, if increased branches in subsequent selectors will be data members of the THSParticle object rather than the whole object (this can lead to name conflicts)");
  ContinueLineAfter("    fOutTree=fChain->CloneTree(0);");
  ContinueLineAfter("    //e.g. fp1=new TLorentzVector(); //should be declared as class data member");
  ContinueLineAfter("    //e.g. fOutTree->Branch(\"p1\",&fp1,buff,split);");
  ContinueLineAfter("  }");
  ContinueLineAfter("  else {fChain->CopyAddresses(fOutTree);}//reset the branch addresses of the cloned tree to the new file tree");
  ContinueLineAfter("  THSOutput::InitOutTree();");

  AddLineAfter("TTree          *fChain","  //data members for new branches",1);
  ContinueLineAfter("   //you must define how they are processed for each event and set to 0 in constructor");
  ContinueLineAfter("   //e.g.   TLorentzVector  *fp1;");
   TString sline=fCurMacro.GetLineWith("virtual ~")->GetString();
   sline.ReplaceAll("{","{/*SafeDelete(fp1);//If we add pointer ro classmust delete it*/ ");
  fCurMacro.GetLineWith("virtual ~")->SetString(sline);

  fCurMacro.SaveSource(fSelName+".h");
  //////////////////////////////////////////////////////////////////  
  //Now with .C
 
}

////////////////////////////////////////////////////////////
/// Adds example of how to create a new output tree
void THSSkeleton::HSNewTree(){
  if(fIsAppendTree)Error("THSSkeleton::HSNewTree()","Can't have a new tree and an append tree!");
 
  fPlace=0;
  //////////////////////////////////////////////////////////////////  
  //First deal with .C
  fCurMacro=TMacro(fSelName+".C");
  AddLineAfter("THSOutput::HSSlaveBegin(fInput,fOutput);","  //Initialise the ouput tree");
  //ContinueLineAfter("   Int_t buff=32000;");
  // ContinueLineAfter("   Int_t split=0;//note split is important in the TSelector framework, if increased branches in subsequent selectors will be data members of the THSParticle object rather than the whole object (this can lead to name conflicts)");
  ContinueLineAfter("   //Only done once. THSOutput::FinishOutput Resets the tree before new files");
  ContinueLineAfter("   //If you want to save additional branches to a new file, do it here");
  ContinueLineAfter("   //If you are creating a HSParticle tree via fStrParticles then do not create a new tree here, just add the extra branches");
  ContinueLineAfter("   fOutTree=new TTree(\"HSParticles\",\"A tree containing reconstructed particles\"); ");
  ContinueLineAfter("  //e.g. Make some objects to make branch with, these should be declared in the .h file ");
  ContinueLineAfter("  //e.g. fp1=new TLorentzVector();");
  ContinueLineAfter("  //e.g.  fOutTree->Branch(\"p1\",&fp1,buff,split);");
  
  fCurMacro.SaveSource(fSelName+".C");

  
  //////////////////////////////////////////////////////////////////  
  //Now with .h
  fCurMacro=TMacro(fSelName+".h");
  AddLineAfter("TTree          *fChain","  //data members for new branches",1);
  ContinueLineAfter("   //you must define how they are processed for each event and set to 0 in constructor");
  ContinueLineAfter("   //e.g.   TLorentzVector  *fp1;");
  TString sline=fCurMacro.GetLineWith("virtual ~")->GetString();
  sline.ReplaceAll("{","{/*//SafeDelete(fp1);//If we add pointer ro classmust delete it*/ ");
  fCurMacro.GetLineWith("virtual ~")->SetString(sline);
  
  AddLineAfter( "THSOutput::HSNotify(fChain);","  THSOutput::InitOutTree();");
 
  fCurMacro.SaveSource(fSelName+".h");


}

///////////////////////////////////////////////////////////////
/// Adds example of how to use the THSWeights class.
void THSSkeleton::HSWeights(){
  fPlace=0;
  //////////////////////////////////////////////////////////////////  
  //First deal with .C
  fCurMacro=TMacro(fSelName+".C");
  //Load Weights
  AddLineAfter("THSOutput::HSSlaveBegin(fInput,fOutput);","   //Get weights from file");
  ContinueLineAfter("   THSHisto::LoadWeights(\"WEIGHTFILE_DIR+NAME\",\"WEIGHTNAME\");");

  //Use weights for histogramming
  if(fIsHisto){
    AddLineAfter("THSHisto::InitialiseHisto","   //Load histograms for each species that has a weight");
    ContinueLineAfter("   THSHisto::LoadCutsForWeights(); //Loop over weights species and make weighted histograms for each");
    
    
    AddLineAfter("FillHistograms(\"Cut1\");","   fWeight=1; //For unweighted histograms",-2);
    AddLineAfter("THSOutput::HSProcessFill();","   THSHisto::GetWeightEvent(fgID); //get weights for this event, use fgID as works when permutating also",-2);
    ContinueLineAfter("   THSHisto::FillCutsForWeights(); //Fill weighted hists");
    
    AddLineAfter("// e.g. ((TH2D*)FindHist(\"Mp1\"))","  //Find and fill weighted hist below here");
    ContinueLineAfter("  //e.g. FindHist(\"Mp1\")->Fill(fp1->M(),fWeight);");
    ContinueLineAfter("  //e.g. ((TH2F*)FindHist(\"Mp1VMmiss\"))->Fill(fp1->M(),Mmiss,fWeight);");
    ReplaceMacroText("// e.g. ((TH2D*)FindHist(\"Mp1\"))->Fill(fp1->M(),fp2->M());//2D","");
    ReplaceMacroText("// e.g. FindHist(\"Mp1\")->Fill(fp1->M());//1D","");
    
  }
  
  //ContinueLineAfter();
  fCurMacro.SaveSource(fSelName+".C");
  //////////////////////////////////////////////////////////////////  
  //Now with .h
  fCurMacro=TMacro(fSelName+".h");
  AddLineAfter("#include \"THSOutput.h\"","#include \"THSWeights.h\"");
  
  AddLineAfter("TTree          *fChain","  //data members for reading weights",1);
  ContinueLineAfter("   THSWeights* fWeights=0;");
  
  TString sline=fCurMacro.GetLineWith("virtual ~")->GetString();
  sline.ReplaceAll("{","{SafeDelete(fWeights); ");
  fCurMacro.GetLineWith("virtual ~")->SetString(sline);

  fCurMacro.SaveSource(fSelName+".h");
  // //////////////////////////////////////////////////////////////////  
  // //Now with Control
  // fCurMacro=TMacro(TString("Control_")+fSelName+".C");
  // TString cline=fCurMacro.GetLineWith("tree->Process(")->GetString();
  // cline.ReplaceAll(")",",\"WEIGHTSFILEHERE.root:NameOfWeightsObject\")");
  // fCurMacro.GetLineWith("tree->Process(")->SetString(cline);

  // fCurMacro.SaveSource(TString("Control_")+fSelName+".C");
}

///////////////////////////////////////////////////////////////////
/// Adds example of how to use the THSLongPS class.
void THSSkeleton::HSLPS(){
  fPlace=0;
  //////////////////////////////////////////////////////////////////  
  //First deal with .C
  fCurMacro=TMacro(fSelName+".C");
  AddLineAfter("THSOutput::HSSlaveBegin(fInput,fOutput);","   //Construct the LongPS object with number of particles in final state");
  //TString strLPS;
  //strLPS.Form("   fLPS=new THSLongPS(%d);",NLPS);
  ContinueLineAfter(Form("   fLPS=new THSLongPS(%d);",fNLPS));

  AddLineAfter("THSOutput::HSProcessStart(entry)","   fLPS->Reset();");
  ContinueLineAfter("   //below you need to give the final state TLorentzVectors to fLPS. Replace ? by the TLorentzVector object. The order gives the particle indice for the sectors");
  for(Int_t ii=0;ii<fNLPS;ii++)
    ContinueLineAfter("   fLPS->AddParticle(?);");
  
  ContinueLineAfter("   fLPS->Analyse();");
  ContinueLineAfter("   fLPS->GetSector();//Can use the sector number now etc.");
  fCurMacro.SaveSource(fSelName+".C");
  //////////////////////////////////////////////////////////////////  
  //Now with .h
  fCurMacro=TMacro(fSelName+".h");
  AddLineAfter("#include \"THSOutput.h\"","#include \"THSLongPS.h\"");
  AddLineAfter("TTree          *fChain","   THSLongPS* fLPS=0;");

  TString sline=fCurMacro.GetLineWith("virtual ~")->GetString();
  sline.ReplaceAll("{","{SafeDelete(fLPS); ");
  fCurMacro.GetLineWith("virtual ~")->SetString(sline);

  fCurMacro.SaveSource(fSelName+".h");


}

////////////////////////////////////////////////////////
/// Adds examples of how to include an analysis project
/// like e.g. THSFinalTemp.C in HSSelector framework
void THSSkeleton::HSFinalState(){
  fPlace=0;
  //////////////////////////////////////////////////////////////////  
  //First deal with .C
  fCurMacro=TMacro(fSelName+".C");

  AddLineAfter("THSOutput::HSSlaveBegin(fInput,fOutput);","   //SetPermutate();//will permutate like particles");
  
  AddLineAfter("// The return value is currently not used.",TString("   "));
  ContinueLineAfter("  ");
  ContinueLineAfter("  //You must define the UserProcess function if you want to do additional");
  ContinueLineAfter("  //analyssis to making the FinalState tree, e.g. histogramming");
  ContinueLineAfter("  THSFinalState::ProcessEvent();");

  ///Add user process function
  FindNextLineLike("Process(Long64_t entry)");
  fPlace--;
  ContinueLineAfter(Form("void %s::UserProcess(){",fFinalName.Data()));
  ContinueLineAfter("   THSOutput::HSProcessFill();",-1);
  ContinueLineAfter("}");
  ReplaceMacroText("THSOutput::HSProcessFill();","");
  
  //If newtree need to link project tree
  if(fIsNewTree){
    AddLineAfter("fOutTree=new TTree",TString("   ")+fFinalName+"::FinalStateOutTree(fOutTree);");
  }
  
  fCurMacro.SaveSource(fSelName+".C");
  //////////////////////////////////////////////////////////////////  
  //Now with .h
  fCurMacro=TMacro(fSelName+".h");
  AddLineAfter("#include \"THSOutput.h\"",TString("#include \"")+fFinalName+".h\"");
  ReplaceMacroText("public THSOutput",TString("public THSOutput, public ")+fFinalName);
  TString branch=FindNextLineLike("fChain->SetBranchAddress(\"Particles\"");	
  ContinueLineAfter(" THSFinalState::SetDetParts(Particles);");
  branch=FindNextLineLike("fChain->SetBranchAddress(\"Generated\"");
  if(branch.Contains("Generated")){
    ContinueLineAfter(" if(fChain->GetBranch(\"Generated\"))THSFinalState::SetGenParts(Generated);");
  }
  fPlace=0;
  FindNextLineLike("Terminate();");
  ContinueLineAfter("  //Additional final state analysis function");
  ContinueLineAfter(Form("  virtual void UserProcess();"));
 
  fCurMacro.SaveSource(fSelName+".h");
  //////////////////////////////////////////////////////////////////  
  //Now with Control
  fCurMacro=TMacro(TString("Control_")+fSelName+".C");
  AddLineAfter("HSout(","  HSMacPath(\"ADDITIONALMACROPATH_WHEREPROJECTIS\");");
  ContinueLineAfter("  HSfinal(\""+fFinalName+"\");");
  fCurMacro.SaveSource(TString("Control_")+fSelName+".C");
}

// void THSSkeleton::HSHisto(){
//   fPlace=0;
//   //////////////////////////////////////////////////////////////////  
//   //First deal with .C
//   fCurMacro=TMacro(fSelName+".C");
//   fCurMacro.SaveSource(fSelName+".C");
//   //////////////////////////////////////////////////////////////////  
//   //Now with .h
//   fCurMacro=TMacro(fSelName+".h");
//   fCurMacro.SaveSource(fSelName+".h");


// }


void THSSkeleton::CreateRooFitEventsPDF(TString pdfName,TString obsNames,TString parNames){
  TObjArray* obss=obsNames.Tokenize(",");
  TObjArray* pars=parNames.Tokenize(",");
  //Check for categories
  vector<Bool_t> is_cat;
  for(Int_t io=0;io<obss->GetEntries();io++)
    if(TString(obss->At(io)->GetName()).Contains("CAT:")){
      TString newstr=obss->At(io)->GetName();
      is_cat.push_back(kTRUE);
      newstr.ReplaceAll("CAT:","");
      dynamic_cast<TObjString*>(obss->At(io))->SetString(newstr); //Get rid of CAT:
      obsNames.ReplaceAll(TString("CAT:")+newstr,newstr);
    }
    else is_cat.push_back(kFALSE);

  TString varNames=obsNames+","+parNames;
  //Make standard RooFit Pdf
  RooClassFactory::makePdf(pdfName,varNames) ;
  //Open code to add RooHSEventsPDF parts
  fCurMacro=TMacro(pdfName+".cxx");
  FindNextLineLike("RooAbsPdf(name,title)");
  ((TObjString*)fCurMacro.GetListOfLines()->At(fPlace))->SetString("   RooHSEventsPDF(name,title),");
  FindNextLineLike("RooAbsPdf(other,name)");
  ((TObjString*)fCurMacro.GetListOfLines()->At(fPlace))->SetString("   RooHSEventsPDF(other,name),");
  //Now add observables to ProxSet
  AddLineAfter("{","   MakeSets();");

 
  //Edit text
  for(Int_t io=0;io<obss->GetEntries();io++)
    ContinueLineAfter(TString("   ")+obss->At(io)->GetName()+".SetName(_"+obss->At(io)->GetName()+".GetName());");
  for(Int_t ip=0;ip<pars->GetEntries();ip++)
    ContinueLineAfter(TString("   ")+pars->At(ip)->GetName()+".SetName(_"+pars->At(ip)->GetName()+".GetName());");

  FindNextLineLike("{");
  ContinueLineAfter("   MakeSets();");
  for(Int_t io=0;io<obss->GetEntries();io++)
    ContinueLineAfter(TString("   ")+obss->At(io)->GetName()+".SetName(other."+obss->At(io)->GetName()+".GetName());");
  for(Int_t ip=0;ip<pars->GetEntries();ip++)
    ContinueLineAfter(TString("   ")+pars->At(ip)->GetName()+".SetName(other."+pars->At(ip)->GetName()+".GetName());");
  ContinueLineAfter("   if(fEvTree) SetEvTree(fEvTree);//Needs fProxSet filled first");
  //Make make sets function need to iterate over obsNames and parNames
  FindNextLineLike("}");
  ContinueLineAfter("void "+pdfName+"::MakeSets(){");
  for(Int_t io=0;io<obss->GetEntries();io++)
    if(!is_cat[io])ContinueLineAfter(TString("   fProxSet.push_back(&")+obss->At(io)->GetName()+");");
    else ContinueLineAfter(TString("   fCatSet.push_back(&")+obss->At(io)->GetName()+");");
  for(Int_t ip=0;ip<pars->GetEntries();ip++)
    ContinueLineAfter(TString("   fParSet.push_back(&")+pars->At(ip)->GetName()+");");
  // ContinueLineAfter("   fVarSet.push_back(RooArgList(\"AllVars\"));");
  ContinueLineAfter("   InitSets();");
  ContinueLineAfter("}");

  //Now evaluateMC template
  MoveToLine("evaluate()");
  FindNextLineLike("}");
  ContinueLineAfter("Double_t "+pdfName+"::evaluateMC() const {",1);
  ContinueLineAfter("// ENTER IDENTICAL EXPRESSION TO evaluate() IN TERMS OF MC VARIABLE ARGUMENTS HERE");
  Int_t nv=0;
  Int_t nc=0;
  for(Int_t io=0;io<obss->GetEntries();io++)
    if(!is_cat[io]) ContinueLineAfter(TString("  Double_t mc")+obss->At(io)->GetName()+Form("=fMCVar[%d];",nv++));
    else ContinueLineAfter(TString("  Int_t mc")+obss->At(io)->GetName()+Form("=fMCCat[%d];",nc++));

  ContinueLineAfter("   return 1.0;");
  ContinueLineAfter("}");
  // ContinueLineAfter(TString("Int_t ")+pdfName+"::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,const char* rangeName) const");
  // ContinueLineAfter("{ //return 0; //might be good to check numerical integral sometimes");
  // ContinueLineAfter("");
  //  ContinueLineAfter("");
  // ContinueLineAfter("");
  
  //Done .C file

  //fix categories
  for(Int_t io=0;io<obss->GetEntries();io++)
    if(is_cat[io]){
      fPlace=0;
      FindNextLineLike(obss->At(io)->GetName());
      ReplaceInCurrLine("RooAbsReal","RooAbsCategory");

    }
  fPlace=0;

  fCurMacro.SaveSource(pdfName+".cxx");


  ////////////////////////////////////////////////////////////////////  
  //Now with .h
  fCurMacro=TMacro(pdfName+".h");
  AddLineAfter("RooAbsPdf.h","#include \"RooHSEventsPDF.h\"");

  ReplaceMacroText("public RooAbsPdf","public RooHSEventsPDF");

    //fix categories
  for(Int_t io=0;io<obss->GetEntries();io++)
    if(is_cat[io]){
      fPlace=0;
      FindNextLineLike(obss->At(io)->GetName());
      ReplaceInCurrLine("RooAbsReal","RooAbsCategory");
      fPlace++;
      FindNextLineLike(obss->At(io)->GetName());
      ReplaceInCurrLine("RooRealProxy","RooCategoryProxy");

    }
  fPlace=0;

  if(!is_cat[0]) AddLineAfter("protected",TString("  Double_t fMC")+obss->At(0)->GetName()+";");
  else AddLineAfter("protected",TString("  Int_t fMC")+obss->At(0)->GetName()+";");
  for(Int_t io=1;io<obss->GetEntries();io++)
    if(!is_cat[io]) ContinueLineAfter(TString("  Double_t fMC")+obss->At(io)->GetName()+";");
    else  ContinueLineAfter(TString("  Int_t fMC")+obss->At(io)->GetName()+";");
  
  AddLineAfter("Double_t evaluate()","  Double_t evaluateMC() const ;");


  // AddLineAfter("inline virtual ~","  virtual Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars,const char* rangeName) const;");
  ContinueLineAfter("  void MakeSets();");
  fCurMacro.SaveSource(pdfName+".h");

}

////////////////////////////////////////////////////
/// Create a FinalState project. \n
/// Set fFinalTopo with SetFinalStateTopo() \n
/// Set fFinalParts with SetFinalStateParts() \n
/// Set fFinalName with SetFinalState() \n
void THSSkeleton::CreateMyFinalState(){

  TObjArray* topos=0;
  if(fFinalTopo==TString(""))cout<<"Warning No finalstate topologies set, you will have to edit the files yourself"<<endl;
  else{
    topos=fFinalTopo.Tokenize(",");
  }
  TObjArray* finals=0;
  if(fFinalParts==TString(""))cout<<"Warning No finalstate final state set, you will have to edit the files yourself"<<endl;
  else{
    finals=fFinalParts.Tokenize(",");
  }
  TString HSANA=gSystem->Getenv("HSANA");
  gSystem->Exec(Form("cp %s/THSFinalTemp.h THS%s.h",HSANA.Data(),fFinalName.Data()));
  //gSystem->Exec(Form("cp %s/THSFinalTemp.h THS%s.h",".",fFinalName.Data()));
  fCurMacro=TMacro(TString("THS")+fFinalName+".h");
  ReplaceAllMacroText("FinalTemp",fFinalName);
  TString UpperFinalName=fFinalName;
  UpperFinalName.ToUpper();
  ReplaceAllMacroText("PROJTEMP",UpperFinalName);

  FindNextLineLike("//Topology flags for this");
  for(Int_t io=0;io<topos->GetEntries();io++)
    ContinueLineAfter(Form("  Int_t fTID_%d=-1;",io));

  fPlace=0;
  FindNextLineLike("//Init functions");
  for(Int_t io=0;io<topos->GetEntries();io++)
    ContinueLineAfter(Form("  void Init_Iter%d();",io));
  for(Int_t io=0;io<topos->GetEntries();io++)
    ContinueLineAfter(Form("  void Topo_%d();",io));

  fPlace=0;
  FindNextLineLike("//Final Particles");
  for(Int_t io=0;io<finals->GetEntries();io++){
    TString sparticle=finals->At(io)->GetName();
    ContinueLineAfter(Form("  THSParticle f%s=THSParticle(\"%s\");",TString(sparticle(0,sparticle.First(":"))).Data(),TString(sparticle(sparticle.First(":")+1,sparticle.Sizeof())).Data()));
  }
  
  fCurMacro.SaveSource(TString("THS")+fFinalName+".h");
  ///////////////////////////////////////////////////////
  
  gSystem->Exec(Form("cp %s/THSFinalTemp.C THS%s.C",HSANA.Data(),fFinalName.Data()));
  fCurMacro=TMacro(TString("THS")+fFinalName+".C");
  fPlace=0;
  ReplaceAllMacroText("FinalTemp",fFinalName);
  fPlace=0;
  FindNextLineLike("//include topology for analysis");
  for(Int_t io=0;io<topos->GetEntries();io++)
    ContinueLineAfter(Form("  fTID_%d=AddTopology(\"%s\");",io,topos->At(io)->GetName()));

  fPlace=0;
  FindNextLineLike("//Do they correspond to a defined topology?");
  for(Int_t io=0;io<topos->GetEntries();io++)
    ContinueLineAfter(Form("    else if(fCurrTopo==fTID_%d) Topo_%d();",io,io));

  fPlace=0;
  FindNextLineLike("//Define topology Iterator functions");
  for(Int_t io=0;io<topos->GetEntries();io++){
    ContinueLineAfter(Form("void THS%s::Init_Iter%d(){",fFinalName.Data(),io));
    ContinueLineAfter("  //THSParticle iterator initialisation");
    ContinueLineAfter("");
    ContinueLineAfter("//If create iterator must configure it");
    ContinueLineAfter(Form("  //fDetIter[fTID_%d]->ConfigureIters();",io));
    ContinueLineAfter("}");
    
  }
  fPlace=0;
  FindNextLineLike("//Define topology functions");
  for(Int_t io=0;io<topos->GetEntries();io++){
    ContinueLineAfter(Form("void THS%s::Topo_%d(){",fFinalName.Data(),io));
    ContinueLineAfter("  //Reconstruct missing or combined particles");
    ContinueLineAfter("");
    ContinueLineAfter("}");
    
  }
  fPlace=0;
  FindNextLineLike("//Set final state");
  ContinueLineAfter("  //Note for CheckTruth algorithm to work all particles added ");
  ContinueLineAfter("  //to fFinal should have been in Generated ");
  for(Int_t io=0;io<finals->GetEntries();io++){
    TString sparticle=finals->At(io)->GetName();
    ContinueLineAfter(Form("  fFinal.push_back(&f%s);",TString(sparticle(0,sparticle.First(":"))).Data()));
  }
  fPlace=0;
  FindNextLineLike("//Initialise particle iters");
  fPlace++;
  for(Int_t io=0;io<topos->GetEntries();io++)
    ContinueLineAfter(Form("  Init_Iter%d();",io));

  fCurMacro.SaveSource(TString("THS")+fFinalName+".C");

}

//////////////////////////////////////////////////////
void THSSkeleton::AddLineAfter(TString line0,TString line1,Int_t off){
  TList *lines=fCurMacro.GetListOfLines();
  MoveToLine(line0);
  fPlace=fPlace+1+off;
  lines->AddAt(new TObjString(line1),fPlace);
}

void THSSkeleton::ContinueLineAfter(TString line1,Int_t off){
  TList *lines=fCurMacro.GetListOfLines();
  fPlace=fPlace+1+off;
  lines->AddAt(new TObjString(line1),fPlace);
}

void THSSkeleton::MoveToLine(TString line0){
  TList *lines=fCurMacro.GetListOfLines();
  TObject* obj=fCurMacro.GetLineWith(line0);
  if(!obj) Error("THSSkeleton::MoveToLine","Line %s does not exist in %s",line0.Data(),"file");
  fPlace=lines->IndexOf(obj); //get line number
}

TString THSSkeleton::FindNextLineLike(TString linelike){
  TList *lines=fCurMacro.GetListOfLines();
  TObjString* thisline=0;
  Int_t count=0;
  for(count=fPlace;count<lines->GetEntries();count++){
    thisline=dynamic_cast<TObjString*>(lines->At(count));
    if(thisline->String().Contains(linelike))
      break;
  }
  if(count==lines->GetEntries()) {fPlace=count; return "";} //didn't find line go to end
  else fPlace=count; //get line number
  return thisline->String();
}

void THSSkeleton::ReplaceInCurrLine(TString text0,TString text1){
  // TString strline=fCurMacro.GetLineWith(text0)->GetString();
  TList *lines=fCurMacro.GetListOfLines();
  TObjString*  thisline=dynamic_cast<TObjString*>(lines->At(fPlace));
  if(!thisline->String().Contains(text0)) cout<<"Warning : ReplaceInCurrLine text not found "<<text0<<" in line "<<thisline->String()<<endl;
  thisline->String().ReplaceAll(text0,text1);
}

void THSSkeleton::ReplaceMacroText(TString text0,TString text1){
  TString strline=fCurMacro.GetLineWith(text0)->GetString();
  strline.ReplaceAll(text0,text1);
  fCurMacro.GetLineWith(text0)->SetString(strline);
}

void THSSkeleton::ReplaceAllMacroText(TString text0,TString text1){
  TList *lines=fCurMacro.GetListOfLines();
  TObjString* thisline=0;
  Int_t count=0;
  for(count=fPlace;count<lines->GetEntries();count++){
    thisline=dynamic_cast<TObjString*>(lines->At(count));
    if(thisline->String().Contains(text0))
      thisline->SetString(thisline->String().ReplaceAll(text0,text1));
  }

}
