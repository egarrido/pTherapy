#include <fstream>
#include <TH1F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TTree.h>
#include <TStyle.h>
#include <TRandom.h>

using namespace std;

void ROOT_pTherapy_dispersion()
{
  //Int_t	bin,minX,maxX;
  Double_t	tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,r;
  //  Double_t	Mean_CPU,Mean_G4,RMS_CPU,RMS_G4;
  TString	name;

  Double_t	tree_En,tree_PosX,tree_PosY,tree_Travel,tree_Theta,tree_Phi,tree_MomX,tree_MomY,tree_MomZ;
  Double_t	tree_EnG4,tree_PosXG4,tree_PosYG4,tree_PosZG4,tree_ThetaG4,tree_PhiG4,tree_MomXG4,tree_MomYG4,tree_MomZG4;
  TFile* rootFile = new TFile("pTherapy.root","RECREATE");
  TTree* CPU =	new TTree("CPU", "Primaire");
  TTree* Geant4	= new TTree("Geant4", "Primaire");
  CPU->Branch("Energy_cinetic",&tree_En);
  CPU->Branch("X",&tree_PosX);
  CPU->Branch("Y",&tree_PosY);
  CPU->Branch("Travel",&tree_Travel);
  CPU->Branch("MomX",&tree_MomX);
  CPU->Branch("MomY",&tree_MomY);
  CPU->Branch("MomZ",&tree_MomZ);
  CPU->Branch("Theta",&tree_Theta);
  CPU->Branch("Phi",&tree_Phi);
  Geant4->Branch("Energy_cinetic",&tree_EnG4);
  Geant4->Branch("X",&tree_PosXG4);
  Geant4->Branch("Y",&tree_PosYG4);
  Geant4->Branch("Z",&tree_PosZG4);
  Geant4->Branch("MomX",&tree_MomXG4);
  Geant4->Branch("MomY",&tree_MomYG4);
  Geant4->Branch("MomZ",&tree_MomZG4);
  Geant4->Branch("Theta",&tree_ThetaG4);
  Geant4->Branch("Phi",&tree_PhiG4);

  ifstream in;
  name="./pTherapy_CPU.fit";
  in.open(Form(name.Data()));
  while (1){
    in >> tmp1 >> tmp2 >> tmp3 >> tmp4 >> tmp5 >> tmp6 >> tmp7;
//    in >> tmp1;
    if (!in.good()) break;
    tree_En	= tmp1;
    tree_PosX	= tmp2;
    tree_PosY	= tmp3;
    tree_Travel	= tmp4;
    tree_MomX	= tmp5;
    tree_MomY	= tmp6;
    tree_MomZ	= tmp7;
    r=sqrt(tmp5*tmp5+tmp6*tmp6+tmp7*tmp7);
    tmp8=acos(tmp7/r)*180./TMath::Pi();
    r=sqrt(tmp5*tmp5+tmp6*tmp6);
    if(r!=0.) tmp9=acos(tmp5/r)*180./TMath::Pi();
    else      tmp9=0.;
    if(tmp6<0.)	tmp9=-tmp9;
    tree_Theta	= tmp8;
    tree_Phi	= tmp9;
    CPU->Fill();
  }
  in.close();

  //name="./pTherapy_CPU_.fit";
  // name="./pTherapy_G4_.fit";
  name="../../../Geant4/pRT_Voxelized_opt3/pTherapy_G4.fit";
  in.open(Form(name.Data()));
  while (1){
    in >> tmp1 >> tmp2 >> tmp3 >> tmp4 >> tmp5 >> tmp6 >> tmp7;
//    in >> tmp1;
    if (!in.good()) break;
    tree_EnG4	= tmp1;
    tree_PosXG4	= tmp2;
    tree_PosYG4	= tmp3;
    tree_PosZG4	= tmp4;
    tree_MomXG4	= tmp5;
    tree_MomYG4	= tmp6;
    tree_MomZG4	= tmp7;
    r=sqrt(tmp5*tmp5+tmp6*tmp6+tmp7*tmp7);
    tmp8=acos(tmp7/r)*180./TMath::Pi();
    r=sqrt(tmp5*tmp5+tmp6*tmp6);
    if(r!=0.) tmp9=acos(tmp5/r)*180./TMath::Pi();
    else      tmp9=0.;
    if(tmp6<0.)	tmp9=-tmp9;
    tree_ThetaG4 = tmp8;
    tree_PhiG4	 = tmp9;
    Geant4->Fill();
  }
  in.close();

  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTickLength(0.01,"XYZ");
{
  TCanvas * c1 = new TCanvas("Comparaison energy","",600,600);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();

  CPU->Draw("Energy_cinetic", "", "P");
  TH1F  *temp_Energy_CPU = (TH1F*)gPad->GetPrimitive("htemp");
  Geant4->Draw("Energy_cinetic", "", "same");
  TH1F  *temp_Energy_G4 = (TH1F*)gPad->GetPrimitive("htemp");
  const double Mean_CPU = temp_Energy_CPU->GetMean();
  const double RMS_CPU  = temp_Energy_CPU->GetRMS();
  double minX=Mean_CPU-5.*RMS_CPU;
  double maxX=Mean_CPU+5.*RMS_CPU;
  const double Mean_G4 = temp_Energy_G4->GetMean();
  const double RMS_G4  = temp_Energy_G4->GetRMS();
  if((Mean_G4-5.*RMS_G4)<minX)
    minX=Mean_G4-5.*RMS_G4;
  if((Mean_G4+5.*RMS_G4)>maxX)
    maxX=Mean_G4+5.*RMS_G4;
  double min_CPU=temp_Energy_CPU->GetXaxis()->GetXmin();
  double min_G4=temp_Energy_G4->GetXaxis()->GetXmin();
  if(min_CPU>min_G4)
  {
    if(minX<min_G4)
      minX=min_G4;
  }
  else
    if(minX<min_CPU)
      minX=min_CPU;  
  double max_CPU=temp_Energy_CPU->GetXaxis()->GetXmax();
  double max_G4=temp_Energy_G4->GetXaxis()->GetXmax();
  if(max_CPU<max_G4)
  {
    if(maxX>max_G4)
      maxX=max_G4;
  }
  else
    if(maxX>max_CPU)
      maxX=max_CPU;  
  // const int bin=(int)(maxX-minX);
  if(RMS_G4==0.||RMS_CPU==0.)
  {
    minX=TMath::Max(min_CPU,min_G4);
    maxX=TMath::Min(max_CPU,max_G4);
  }  
  const int bin=100;

  TH1F * Energy_CPU = new TH1F ("the_Energy_CPU", "", bin, minX, maxX);
  CPU->Draw("Energy_cinetic>>the_Energy_CPU");

  TH1F * Energy_G4 = new TH1F ("the_Energy_G4", "", bin, minX, maxX);
  Geant4->Draw("Energy_cinetic>>the_Energy_G4");

  Energy_CPU->SetLineColor(2);
  Energy_CPU->SetMarkerStyle(4);
  Energy_CPU->SetMarkerColor(2);
  Energy_CPU->GetXaxis()->SetLabelFont(43); //font in pixels
  Energy_CPU->GetXaxis()->SetLabelSize(16); //in pixels
  Energy_CPU->GetYaxis()->SetLabelFont(43); //font in pixels
  Energy_CPU->GetYaxis()->SetLabelSize(16); //in pixels
  Energy_CPU->GetXaxis()->SetTitleFont(43); //font in pixels
  Energy_CPU->GetXaxis()->SetTitleSize(18); //in pixels
  Energy_CPU->GetYaxis()->SetTitleFont(43); //font in pixels
  Energy_CPU->GetYaxis()->SetTitleSize(18); //in pixels
  Energy_CPU->GetYaxis()->CenterTitle(kTRUE);
  Energy_CPU->GetXaxis()->CenterTitle(kTRUE);
  Energy_CPU->GetYaxis()->SetTitleOffset(1.4);
  Energy_CPU->GetXaxis()->SetTitleOffset(3.0);
  Energy_CPU->SetTitle("; Energy [MeV]; Particle [u.a]");

  Energy_CPU->DrawCopy ("P");
  Energy_G4->DrawCopy ("same");

  c1->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  Energy_CPU->Sumw2();
  Energy_CPU->SetStats(0);
  Energy_CPU->Add(Energy_CPU,Energy_G4, 100., -100.);
  Energy_CPU->Divide(Energy_G4);
  //Energy_CPU->GetYaxis()->SetRangeUser (-12, 12);
  Energy_CPU->GetYaxis()->SetTitle ("Dispersion [%]");
  Energy_CPU->SetMarkerStyle(1);
  Energy_CPU->DrawCopy("EP");
}
{
  TCanvas * c1 = new TCanvas("Comparaison theta","",600,600);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();

  CPU->Draw("Theta", "", "P");
  TH1F  *temp_Theta_CPU = (TH1F*)gPad->GetPrimitive("htemp");
  Geant4->Draw("Theta", "", "same");
  TH1F  *temp_Theta_G4 = (TH1F*)gPad->GetPrimitive("htemp");
  const double Mean_CPU = temp_Theta_CPU->GetMean();
  const double RMS_CPU  = temp_Theta_CPU->GetRMS();
  const double Mean_CPU = temp_Theta_CPU->GetMean();
  const double RMS_CPU  = temp_Theta_CPU->GetRMS();
  double minX=Mean_CPU-5.*RMS_CPU;
  double maxX=Mean_CPU+5.*RMS_CPU;
  const double Mean_G4 = temp_Theta_G4->GetMean();
  const double RMS_G4  = temp_Theta_G4->GetRMS();
  if((Mean_G4-5.*RMS_G4)<minX)
    minX=Mean_G4-5.*RMS_G4;
  if((Mean_G4+5.*RMS_G4)>maxX)
    maxX=Mean_G4+5.*RMS_G4;
  double min_CPU=temp_Theta_CPU->GetXaxis()->GetXmin();
  double min_G4=temp_Theta_G4->GetXaxis()->GetXmin();
  if(min_CPU>min_G4)
  {
    if(minX<min_G4)
      minX=min_G4;
  }
  else
    if(minX<min_CPU)
      minX=min_CPU;  
  double max_CPU=temp_Theta_CPU->GetXaxis()->GetXmax();
  double max_G4=temp_Theta_G4->GetXaxis()->GetXmax();
  if(max_CPU<max_G4)
  {
    if(maxX>max_G4)
      maxX=max_G4;
  }
  else
    if(maxX>max_CPU)
      maxX=max_CPU;  
  if(RMS_G4==0.||RMS_CPU==0.)
  {
    minX=TMath::Max(min_CPU,min_G4);
    maxX=TMath::Min(max_CPU,max_G4);
  }  
  const int bin=100;

  TH1F * Theta_CPU = new TH1F ("the_Theta_CPU", "", bin, minX, maxX);
  CPU->Draw("Theta>>the_Theta_CPU");

  TH1F * Theta_G4 = new TH1F ("the_Theta_G4", "", bin, minX, maxX);
  Geant4->Draw("Theta>>the_Theta_G4");

  Theta_CPU->SetLineColor(2);
  Theta_CPU->SetMarkerStyle(4);
  Theta_CPU->SetMarkerColor(2);
  Theta_CPU->GetXaxis()->SetLabelFont(43); //font in pixels
  Theta_CPU->GetXaxis()->SetLabelSize(16); //in pixels
  Theta_CPU->GetYaxis()->SetLabelFont(43); //font in pixels
  Theta_CPU->GetYaxis()->SetLabelSize(16); //in pixels
  Theta_CPU->GetXaxis()->SetTitleFont(43); //font in pixels
  Theta_CPU->GetXaxis()->SetTitleSize(18); //in pixels
  Theta_CPU->GetYaxis()->SetTitleFont(43); //font in pixels
  Theta_CPU->GetYaxis()->SetTitleSize(18); //in pixels
  Theta_CPU->GetYaxis()->CenterTitle(kTRUE);
  Theta_CPU->GetXaxis()->CenterTitle(kTRUE);
  Theta_CPU->GetYaxis()->SetTitleOffset(1.4);
  Theta_CPU->GetXaxis()->SetTitleOffset(3.0);
  Theta_CPU->SetTitle("; Theta [deg]; Particle [u.a]");

  Theta_CPU->DrawCopy ("P");
  Theta_G4->DrawCopy ("same");

  c1->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  Theta_CPU->Sumw2();
  Theta_CPU->SetStats(0);
  Theta_CPU->Add(Theta_CPU,Theta_G4, 100., -100.);
  Theta_CPU->Divide(Theta_G4);
  //Theta_CPU->GetYaxis()->SetRangeUser (-12, 12);
  Theta_CPU->GetYaxis()->SetTitle ("Dispersion [%]");
  Theta_CPU->SetMarkerStyle(1);
  Theta_CPU->DrawCopy("EP");
}
{
  TCanvas * c1 = new TCanvas("Comparaison X","",600,600);
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetBottomMargin(0);
  pad1->Draw();
  pad1->cd();

  CPU->Draw("X", "", "P");
  TH1F  *temp_X_CPU = (TH1F*)gPad->GetPrimitive("htemp");
  Geant4->Draw("X", "", "same");
  TH1F  *temp_X_G4 = (TH1F*)gPad->GetPrimitive("htemp");
  const double Mean_CPU = temp_X_CPU->GetMean();
  const double RMS_CPU  = temp_X_CPU->GetRMS();
  const double Mean_CPU = temp_X_CPU->GetMean();
  const double RMS_CPU  = temp_X_CPU->GetRMS();
  double minX=Mean_CPU-5.*RMS_CPU;
  double maxX=Mean_CPU+5.*RMS_CPU;
  const double Mean_G4 = temp_X_G4->GetMean();
  const double RMS_G4  = temp_X_G4->GetRMS();
  if((Mean_G4-5.*RMS_G4)<minX)
    minX=Mean_G4-5.*RMS_G4;
  if((Mean_G4+5.*RMS_G4)>maxX)
    maxX=Mean_G4+5.*RMS_G4;
  double min_CPU=temp_X_CPU->GetXaxis()->GetXmin();
  double min_G4=temp_X_G4->GetXaxis()->GetXmin();
  if(min_CPU>min_G4)
  {
    if(minX<min_G4)
      minX=min_G4;
  }
  else
    if(minX<min_CPU)
      minX=min_CPU;  
  double max_CPU=temp_X_CPU->GetXaxis()->GetXmax();
  double max_G4=temp_X_G4->GetXaxis()->GetXmax();
  if(max_CPU<max_G4)
  {
    if(maxX>max_G4)
      maxX=max_G4;
  }
  else
    if(maxX>max_CPU)
      maxX=max_CPU;  
  if(RMS_G4==0.||RMS_CPU==0.)
  {
    minX=TMath::Max(min_CPU,min_G4);
    maxX=TMath::Min(max_CPU,max_G4);
  }  
  const int bin=100;

  TH1F * X_CPU = new TH1F ("the_X_CPU", "", bin, minX, maxX);
  CPU->Draw("X>>the_X_CPU");

  TH1F * X_G4 = new TH1F ("the_X_G4", "", bin, minX, maxX);
  Geant4->Draw("X>>the_X_G4");

  X_CPU->SetLineColor(2);
  X_CPU->SetMarkerStyle(4);
  X_CPU->SetMarkerColor(2);
  X_CPU->GetXaxis()->SetLabelFont(43); //font in pixels
  X_CPU->GetXaxis()->SetLabelSize(16); //in pixels
  X_CPU->GetYaxis()->SetLabelFont(43); //font in pixels
  X_CPU->GetYaxis()->SetLabelSize(16); //in pixels
  X_CPU->GetXaxis()->SetTitleFont(43); //font in pixels
  X_CPU->GetXaxis()->SetTitleSize(18); //in pixels
  X_CPU->GetYaxis()->SetTitleFont(43); //font in pixels
  X_CPU->GetYaxis()->SetTitleSize(18); //in pixels
  X_CPU->GetYaxis()->CenterTitle(kTRUE);
  X_CPU->GetXaxis()->CenterTitle(kTRUE);
  X_CPU->GetYaxis()->SetTitleOffset(1.4);
  X_CPU->GetXaxis()->SetTitleOffset(3.0);
  X_CPU->SetTitle("; X [mm]; Particle [u.a]");

  X_CPU->DrawCopy ("P");
  X_G4->DrawCopy ("same");

  c1->cd();
  TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);
  pad2->SetTopMargin(0);
  pad2->SetGridy();
  pad2->SetBottomMargin(0.2);
  pad2->Draw();
  pad2->cd();
  X_CPU->Sumw2();
  X_CPU->SetStats(0);
  X_CPU->Add(X_CPU,X_G4, 100., -100.);
  X_CPU->Divide(X_G4);
  //X_CPU->GetYaxis()->SetRangeUser (-12, 12);
  X_CPU->GetYaxis()->SetTitle ("Dispersion [%]");
  X_CPU->SetMarkerStyle(1);
  X_CPU->DrawCopy("EP");
}

  c1->Modified();
  c1->SetFillColor(0);
  c1->Update();

  rootFile->Write();
  rootFile->Close();
}
