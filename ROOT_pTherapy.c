#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "Riostream.h"
#include <TString.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TPaveText.h>
#include <TLegend.h>
#define pi	3.14159265358979323846

using namespace std;
void ROOT_pTherapy()
{
	Double_t	tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,r;
	TString	name;
	
	Int_t			tree_No,tree_event;
	Double_t	tree_En,tree_Edep,tree_PosX,tree_PosY,tree_Travel,tree_Theta,tree_Phi,tree_MomX,tree_MomY,tree_MomZ;
	TFile* rootFile = new TFile("pTherapy.root","RECREATE");
	TTree* Arbre		=	new TTree("pTherapy", "Primaire");
	Arbre->Branch("Energy",&tree_En);
  Arbre->Branch("X",&tree_PosX);  
  Arbre->Branch("Y",&tree_PosY);
  Arbre->Branch("Travel",&tree_Travel);
  Arbre->Branch("MomX",&tree_MomX);
  Arbre->Branch("MomY",&tree_MomY);
  Arbre->Branch("MomZ",&tree_MomZ);
  Arbre->Branch("Theta",&tree_Theta);
	Arbre->Branch("Phi",&tree_Phi);

  ifstream in;
  name="pTherapy_CPU.fit";
  in.open(Form(name.Data()));
 
  while (1){
    in >> tmp1 >> tmp2 >> tmp3 >> tmp4 >> tmp5 >> tmp6 >> tmp7;
    if (!in.good()) break;
    tree_En			=	tmp1;
    tree_PosX		=	tmp2;
    tree_PosY		=	tmp3;
    tree_Travel	= tmp4;
    tree_MomX		=	tmp5;
    tree_MomY		=	tmp6;
    tree_MomZ		=	tmp7;
    r=sqrt(tmp5*tmp5+tmp6*tmp6+tmp7*tmp7);
    tmp8=acos(tmp7/r)*180./pi;
    r=sqrt(tmp5*tmp5+tmp6*tmp6);
    if(r!=0.)		tmp9=acos(tmp5/r)*180./pi;
    else    		tmp9=0.;
    if(tmp6<0.)	tmp9=-tmp9;		
    tree_Theta	=	tmp8;
    tree_Phi		=	tmp9;
    if(tree_En>0.)
    	Arbre->Fill();
  }
  in.close();
  rootFile->Write();
  rootFile->Close();
}
