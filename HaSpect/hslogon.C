/*! \file
	Sets up HaSpect environment when included in .rootrc file.
	
	### Usage
	Use flags when starting root to load needed functions \n
	- \--hsfit : Classes for sWeight fits etc
	- \--hssel : Classes for hselector
	- \--hsdata : Classes for hsdata
	- \--hsfinal=MyFinalState : Classes for hsfinalstate where\n
	MyFinalState is the users own derived class
	
	e.g. 
	\code root -l --hsfit FitHSSimple.C \endcode
	to start FitHSSimple.C macro which uses sWeights to fit data.
	\code root --cleanall \endcode
	to remove all .so .d and .pcm files
	\code root RunSomeScript.C --proof=N \endcode
	to initialise a proof session with N workers
	\code root --MyClass.C \endcode
	to load your own class from current directory, or a class from HSANA or HSUSER
	
*/

#include <TProof.h> //Make sure gProof can be seen here

void HSfit();//Load hsfit classes
void HSdata();//Load hsdata classes
void HSselector(); //load hsselector classes
void HSFinal(TString pname); //load hsproject classes
void startproof(Int_t Nw); //intialise proof for hsselector classes
void LoadMacro(TString macro); //Load class via its source code
TString HSin(); //return inout void HSMacPath(TString opt)files directory
TString HSout(); //return inout files directory
TString HSfinal(); //return project classname
void HSin(TString hsin){gSystem->Setenv("HSIN",hsin);} //set in files directory
void HSout(TString hsout){gSystem->Setenv("HSOUT",hsout);} //set out file
void HSfinal(TString hsfinal){gSystem->Setenv("HSFINAL",hsfinal);} //set out file
void CleanAll();
void HSMacPath(TString opt);
void HSUserPath(TString opt);
TString MACPATH; //additional macro path needed for proof
Bool_t gISFARM=kFALSE;
TString THSPARTICLE("THSParticle.C");

///////////////////////////////////////////////////////
///Initialises functions demanded by the user
void hslogon(){
  //Execute system .rootlogon.C script 
  if(gSystem->AccessPathName("~/.rootlogon.C")){};

  //Check if user source code directory defined
  TString HSANA=gSystem->Getenv("HSANA");
  TString HSUSER=gSystem->Getenv("HSUSER");

  if(gSystem->Getenv("HSUSER")){
    gROOT->SetMacroPath(Form("%s:%s",HSUSER.Data(),gROOT->GetMacroPath()));
    gSystem->AddIncludePath(TString("-I")+HSUSER);
  }

  //get command line options first check if makeall
  for(Int_t i=1;i<gApplication->Argc();i++){
    TString opt=gApplication->Argv(i);
    //look for --proof=Nworkers optionif Nworkers not given all cores will be used
    if((opt.Contains("--cleanall"))) CleanAll();
    
  }
  //check if additional macro dir given
   for(Int_t i=1;i<gApplication->Argc();i++){
    TString opt=gApplication->Argv(i);
    //look for --proof=Nworkers optionif Nworkers not given all cores will be used
    if((opt.Contains("--macrodir"))){
      HSMacPath(TString(opt(11,opt.Sizeof())));
    }
   }

  
   //get command line options first check if proof
  for(Int_t i=1;i<gApplication->Argc();i++){
    TString opt=gApplication->Argv(i);
    //look for --proof=Nworkers optionif Nworkers not given all cores will be used
    if((opt.Contains("--proof"))) startproof(TString(opt(8,opt.Sizeof())).Atoi());
    
  }
  
   //get command line options first check if rootbeer
  for(Int_t i=1;i<gApplication->Argc();i++){
    TString opt=gApplication->Argv(i);
    if((opt.Contains("--rootbeer"))) gROOT->ProcessLine(".L $ROOTBEER/RootBeerSetup.cxx");
  }

   //get command line options first check if farm job
  for(Int_t i=1;i<gApplication->Argc();i++){
    TString opt=gApplication->Argv(i);
    //look for --farm and copy all hsroot files to local hsana directory
    if((opt.Contains("--farm"))){
      gISFARM=kTRUE;
      gSystem->Exec(Form("mkdir hsana"));
      gSystem->Exec(Form("cp %s/*.h hsana/.",HSANA.Data()));
      gSystem->Exec(Form("cp %s/*.C hsana/.",HSANA.Data()));
      gSystem->Setenv("HSANA","./");
      if(gSystem->Getenv("RHIPO")){
	TString RHIPO=gSystem->Getenv("RHIPO");
	gSystem->Exec(Form("cp %s/THipo.h hsana/.",RHIPO.Data()));
	gSystem->Exec(Form("cp %s/THipo.C hsana/.",RHIPO.Data()));
	gSystem->Exec(Form("cp %s/Hipo2Root.C hsana/.",RHIPO.Data()));
	gSystem->Setenv("RHIPO","./");	
      }
      if(gSystem->Getenv("CHIPO")){
	TString CHIPO=gSystem->Getenv("CHIPO");
	gSystem->Exec(Form("cp %s/*.h hsana/.",CHIPO.Data()));
	gSystem->Exec(Form("cp %s/*.cpp hsana/.",CHIPO.Data()));
	gSystem->Setenv("CHIPO","./");	
      }
      gROOT->SetMacroPath(TString("hsana/:")+gROOT->GetMacroPath());
    }
  }
  
  //get command line options and call assiciated functions
  for(Int_t i=1;i<gApplication->Argc();i++){
    TString opt=gApplication->Argv(i);
    if(opt.Contains("--hsin")) {HSin(TString(opt(7,opt.Sizeof())));} //Set input tree directory
    if(opt.Contains("--hsout")) {HSout(TString(opt(8,opt.Sizeof())));} //Set input tree directory
    if((opt==TString("--hsfit"))) HSfit(); //Load fit classes
    if((opt==TString("--hsdata"))) HSdata(); //Load data classes
    if((opt==TString("--hssel"))) HSselector(); //Load selector classes
    if((opt.Contains("--hsfinal"))) HSFinal(TString(opt(10,opt.Sizeof()))); //Load finalstate classes
    if(opt.Contains("--")&&opt.Contains(".cxx")){opt.Remove(0,2); cout<<"Loading "<<opt<<endl;LoadMacro(opt);} //Load additional THS classes
    if(opt.Contains("--")&&opt.Contains(".C")){opt.Remove(0,2); cout<<"Loading "<<opt<<endl;LoadMacro(opt);} //Load additional THS classes
    
  }	
  
}


void CleanAll(){

  TString HSANA=gSystem->Getenv("HSANA");
  cout<<"Tidying up "<<HSANA<<endl;
  gSystem->Exec(Form("rm %s/*.so",HSANA.Data()));
  gSystem->Exec(Form("rm %s/*.d",HSANA.Data()));
  gSystem->Exec(Form("rm %s/*.pcm",HSANA.Data()));
  if(gSystem->Getenv("RHIPO")){
    TString RHIPO=gSystem->Getenv("RHIPO");
    gSystem->Exec(Form("rm %s/*.so",RHIPO.Data()));
    gSystem->Exec(Form("rm %s/*.d",RHIPO.Data()));
    gSystem->Exec(Form("rm %s/*.pcm",RHIPO.Data()));
  }
  if(gSystem->Getenv("CHIPO")){
    TString CHIPO=gSystem->Getenv("CHIPO");
    gSystem->Exec(Form("rm %s/*.so",CHIPO.Data()));
    gSystem->Exec(Form("rm %s/*.d",CHIPO.Data()));
    gSystem->Exec(Form("rm %s/*.pcm",CHIPO.Data()));
  }


}

/** Function is called with \--hsfit \n
 * It loads necessary functions to perform RooFit/RooStats fits:\n
 * THSBins.C
 * THSWeights.C
 * RooHSEventsPDF.C
 * RooHSEventsHistPDF.C
 * THSRooFit.C
 * THSsPlot.C
 * RooHS1StepStudy.C
 * RooHSsPlotAndFitStudy.C
 * RooHSStudyManager.C
 */
void HSfit(){
  cout<<"Loading HSFit classes"<<endl;
  //set prompt
  if(!gProof)((TRint*)gROOT->GetApplication())->SetPrompt("hsroot [%d] ");
  //Load extra classes for roofit
  LoadMacro("THSBins.C"); 
  LoadMacro("THSWeights.C");
  LoadMacro("RooHSEventsPDF.C");
  LoadMacro("RooHSEventsHistPDF.C");
  LoadMacro("THSRooFit.C");
  LoadMacro("THSsPlot.C");
  LoadMacro("RooHS1StepStudy.C");
  LoadMacro("RooHSsPlotAndFitStudy.C");
  LoadMacro("RooHSStudyManager.C");
}

/** Function is called with \--hssel \n
 * It loads the following functions:\n
 * THSBins.C
 * THSWeights.C
 * THSHisto.C
 * THSOutput.C
 * THSSkeleton.C
 * THSKinematics.C
 * THSParticle.C (unless explicitly replaced)
 */
void HSselector(){
  cout<<"Loading HSSel classes"<<endl;
 //set prompt
  if(!gProof)((TRint*)gROOT->GetApplication())->SetPrompt("hsroot [%d] ");
   //Load classess for selector analysis
  //particle container class
  LoadMacro("THSBins.C");
  LoadMacro("THSWeights.C");
  LoadMacro(THSPARTICLE);
  LoadMacro("THSHisto.C");
  LoadMacro("THSOutput.C");
  LoadMacro("THSSkeleton.C");
  LoadMacro("THSKinematics.C");
}

/** Function is called with \--hsfinal \n
 * It loads the following functions:\n
 * [pname].C
 * THSWeights.C
 * THSDataManager.C
 * THSFinalState.C
 * THSKinematics.C
 * THSParticle.C
 */
void HSFinal(TString pname){
  HSfinal(pname);
  if(!TClass::GetClass("THSParticle")) LoadMacro("THSParticle.C");
  LoadMacro(THSPARTICLE);
  LoadMacro("THSWeights.C");
  LoadMacro("THSDataManager.C");
  LoadMacro("THSKinematics.C");
  LoadMacro("THSCombitorial.C");
  LoadMacro("THSParticleIter.C");
  LoadMacro("THSFinalState.C");
  LoadMacro(pname+".C");

}

/** Function is called with \--hsdata \n
 * It loads the following functions:\n
 * THSWeights.C
 * THSParticle.C
 * THSDataManager.C
 * THSLundReader.C
 * THSHipoReader.C (if RHIPO set)
 */
void HSdata(){
  if(gSystem->Getenv("RHIPO"))
    gROOT->ProcessLine(TString(".x ")+gSystem->Getenv("RHIPO")+"/Hipo2Root.C");
  cout<<"THSPARTICLE= "<<THSPARTICLE<<endl;
 
  LoadMacro("THSWeights.C");
  //if(!TClass::GetClass("THSParticle"))   LoadMacro("THSParticle.C");
  LoadMacro(THSPARTICLE);
  LoadMacro("THSDataManager.C");
  LoadMacro("THSLundReader.C");
  if(gSystem->Getenv("RHIPO"))
    LoadMacro("THSHipoReader.C");

}

//////////////////////////////////
/// Start proof with Nw workers
void startproof(Int_t Nw){
  cout<<"Loading Proof setup"<<endl;
  //set prompt
  ((TRint*)gROOT->GetApplication())->SetPrompt("hsproof [%d] ");
  TProof::Open("://lite");
  if(Nw) gProof->SetParallel(Nw); //set number of workers if specified on command line
}


void LoadMacro(TString macro){
  //Different methods depending on whether proof or not
  TString HSANA=gSystem->Getenv("HSANA");
  TString HSUSER=gSystem->Getenv("HSUSER");

  cout<<"Loading hsmacro "<<macro<<endl;
  if(gProof) {
    if(gSystem->Which("./",macro)){
      macro=macro;
    }	
    else if(gSystem->Which(HSUSER,macro)){
      macro=HSUSER+"/"+macro;
    }
    else  if(gSystem->Which(MACPATH,macro)){
      macro=MACPATH+"/"+macro;
    }
    else  if(gSystem->Which(HSANA,macro)){
      macro=HSANA+"/"+macro;
    }
    else {cout<<"Can't find macro "<<macro<<endl; return;}
    gProof->Load(macro+"+");
  }
  else gROOT->LoadMacro(macro+"+");//don't use HSANA in case want to overwrite with macro in current directory
  
}

//////////////////////////////////////////
/// Checks for and returns HSPROJ variable
TString HSfinal(){
  if(!gSystem->Getenv("HSFINAL")) cout<<"Warning no HSFINAL env variable defined but hsfinal() called..."<<endl;
  return TString(gSystem->Getenv("HSFINAL"));
}

//////////////////////////////////////////
/// Checks for and returns HSIN variable
TString HSin(){
  if(!gSystem->Getenv("HSIN")) cout<<"Warning no HSIN env variable defined but hsin() called..."<<endl;
  return TString(gSystem->Getenv("HSIN"))+"";
}

//////////////////////////////////////////
/// Checks for and returns HSOUT variable
TString HSout(){
  if(!gSystem->Getenv("HSOUT")) cout<<"Warning no HSOUT env variable defined but hsout() called..."<<endl;
  if(gProof) cout<<"Warning : Using proof so need full path to output"<<endl;
  return gSystem->Getenv("HSOUT");
}

/////////////////////////////////////////
/// Adds opt to MACPATH (macro path)
void HSMacPath(TString opt){
  cout<<"Adding Macro Path "<<opt<<endl;
  MACPATH=opt;
  gROOT->SetMacroPath(Form("%s:%s",gROOT->GetMacroPath(),MACPATH.Data()));
  gSystem->AddIncludePath(TString("-I")+MACPATH);
}
