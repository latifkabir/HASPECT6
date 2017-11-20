/** \file CreateSelector.C
 * 
 * Macro to generate a TSelector
 * using the THSSkeleton class
 * 
 * Usage:
 * root --hssel CreateSelector.C
 * 
 */

void CreateSelector(){
	
	THSSkeleton* sk = new THSSkeleton();
	
	sk->SetHisto();   //Going to make histograms
	sk->SetNewTree();   //Going to make new tree
	
	sk->SetFinalState("THSep_omegap"); //Give final state class name
	
	sk->CreateSelector("Histo","/home/dglazier/Dropbox/hsana/Events/HASPECT6/Projects/finalstates/omega/data/hs_4_8_3/out_gemc_32975714.farmpbs14_t-1_s1.root","HSParticles"); //create selector with all specifications
	
};
