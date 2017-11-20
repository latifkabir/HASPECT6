//run with root --THSSkeleton.C CreateMyFinalState.C
{
  THSSkeleton* sk=new THSSkeleton();
  //Give your project a name
  sk->SetFinalState("WillThisWork",kTRUE);  //creating the final state class 
  //Set the detected particle combinations you will analyse
  sk->SetFinalStateTopo("pi+:pi-,pi+:pi-:pi0,proton:pi+:pi-:pi0");
  //Set the actual final state particles of the reaction
  //These will just be used as the data member names
  sk->SetFinalStateParts("Electron:e-,Proton:proton,Pip:pi+,Pim:pi-");
  //Make some code
  sk->CreateMyFinalState();
}
