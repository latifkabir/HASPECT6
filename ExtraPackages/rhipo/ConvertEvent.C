//root $RHIPO/Hipo2Root.C $RHIPO/ConvertSomeJSON.C file1.hipo file2.hipo ...
{
  //Edit this macro to only convert the banks from the selected .json files
  THipo Hipo("EVENT");
  Hipo.AddCommandLineFiles();
  Hipo.ToRoot(); //process files given on command line

  //finished exit here so it doesn't try to process hipo files!
  gROOT->ProcessLine(".q");
}
