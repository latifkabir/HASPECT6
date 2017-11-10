/** \file ConfigSel.C
 * 
 * DEPRECATED
 * 
 * For guidance on how to create a THSSelector look at THSSkeleton and its example
 * 
 */
{
  //You need to have the enviroment variable HSANA defined to link to the THS source
  TString HSANA(gSystem->Getenv("HSANA"));
  gROOT->LoadMacro(HSANA+"/MakeHSSelector.C");

  //Set arguments for MakeHSSelector
 IsHisto=kFALSE;  //use THSHisto?
 IsAppendTree=kFALSE; //Append branches to the input tree
 IsNewTree=kFALSE;  //Output a brand new tree
 IsHSTree=kFALSE;   //Use THSOuput to THSParticle interface (probably not)
 IsQval=kFALSE;   //Use Qvalue event weighting algorithm
 IsSW=kFALSE;   //Use sWeights event weighting algorithm
 NLPS=0;   //Use Longitudinal phase space class with NLPS particles

 FileName="";   // The input filename containing the tree
 TreeName="";   // The name of the tree
 OutName="";   // The name of the output directory or file
 SelName="";    // The name of the selector to be produced

 //IsParent=kTRUE; //use connection with parent tree
 //StepN=0;        // If connecting with parent, which step number is its selector in
 //ParSel="";

 //Make the selector skeleton + control files in this directory 
 //This is based on the info given above
 MakeHSSelector();

}
