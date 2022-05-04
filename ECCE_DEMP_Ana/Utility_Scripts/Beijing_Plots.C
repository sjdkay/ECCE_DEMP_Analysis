// 26/04/22 - Stephen JD Kay, University of Regina

// A quick script to make some plots for Garth's Beijing talk
// Execute via root -l 'Beijing_Plots.C("INPUT_FILE.root")'
#define Beijing_Plots_cxx
          
// Include relevant stuff
#include <TStyle.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TGaxis.h>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <TSystem.h>
#include <TTree.h>
#include <TArc.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>

void Beijing_Plots(string InFilename = ""){

  gROOT->SetBatch(kTRUE); // Force script to always run without flashing up graphics
  
  TString rootFile;

  if(InFilename == ""){
    cout << "Enter a filename to analyse: ";
    cin >> InFilename;
  }
  
  TString TInFilename = InFilename;
  rootFile = TInFilename;

  if(gSystem->AccessPathName(rootFile) == kTRUE){
    cerr << "!!!!! ERROR !!!!!" << endl << rootFile << " not found" << endl << "!!!!! ERROR !!!!!" << endl;
    exit(0);
  }

  TFile *InFile = new TFile(rootFile);

  TLatex InfoDump;
  InfoDump.SetTextSize(0.1);
  InfoDump.SetTextAlign(12);  //align at centre
  
  TH1F* tHists[8];
  TH1F* ttruthHists[8];
  TH1F* Q2Hists[8];
  TH1F* WHists[8];
  TH2F* Q2tHists[8];
  TH2F* tvst_Q2_Hists[8];
  TH2F* taltvst_Q2_Hists[8];
  TH1F* taltres_Hists[10]; // 22/02/22 - SJDK - t binned t resolution plots  
  
  TH2F* pipThetaTruthHist = (TH2F*)((TH2F*)InFile->Get("Pion_Truth_Info/piTrack_pTheta_Truth"));
  TH2F* epThetaTruthHist = (TH2F*)((TH2F*)InFile->Get("Scattered_Electron_Truth_Info/eTrack_pTheta_Truth"));
  TH2F* npThetaTruthHist = (TH2F*)((TH2F*)InFile->Get("Neutron_Truth_Info/nTrack_pTheta_Truth"));

  TH2F* Q2WHist = (TH2F*)((TH2F*)InFile->Get("Physics_Results_Misc/Q2_W_Result"));
  TH2F* piXYHist = (TH2F*)((TH2F*)InFile->Get("Pion_Info/pi_XY"));
  TH2F* eXYHist = (TH2F*)((TH2F*)InFile->Get("Scattered_Electron_Info/e_XY"));
  TH2F* ZDCHist;
   
  // SJDK - Check which ZDC hist exists, grab the one that does
  if ( InFile->GetListOfKeys()->Contains("ZDC_XY_IP6")){ 
    ZDCHist = (TH2F*)((TH2F*)InFile->Get("ZDC_XY_IP6"));
  }
  else{
    ZDCHist = (TH2F*)((TH2F*)InFile->Get("ZDC_XY_IP8"));
  }

  TH2F* nThetaPhiDiff;
  // 22/02/22 - SJDK - If the t binned t resolution plots exist, grab them
  if (((TDirectory*)InFile->Get("ZDC_Neutron_Info"))->GetListOfKeys()->Contains("n_ThetaPhiDiff")){
    nThetaPhiDiff = (TH2F*)((TH2F*)InFile->Get("ZDC_Neutron_Info/n_ThetaPhiDiff"));
  }
  else{
    nThetaPhiDiff = new TH2F("h2_Empty_Plot", "No Plot", 100, -10, 10, 100, -10, 10);
  }

  // 22/02/22 - SJDK - If the t binned t resolution plots exist, grab them
  if (((TDirectory*)InFile->Get("Physics_Results_Misc"))->GetListOfKeys()->Contains("taltres_result_ttruth_1")){
    for(Int_t A = 0; A <10; A++){
      taltres_Hists[A] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results_Misc/taltres_result_ttruth_%i",(A+1))));
    }
  }
  
  TH1F* Q2EffHist = (TH1F*)((TH1F*)InFile->Get("Detection_Efficiency/Q2_DetEff"));
  TH2F* Q2tEffHist = (TH2F*)((TH2F*)InFile->Get("Detection_Efficiency/Q2_t_DetEff"));
  TH2F* Q2tEffHist_v2 = (TH2F*)((TH2F*)InFile->Get("Detection_Efficiency/Q2_t_DetEff_v2"));

  TH2F* tvstHist = (TH2F*)((TH2F*)InFile->Get("Physics_Results_Misc/t_ttruth_result"));
  TH2F* t_altvstHist = (TH2F*)((TH2F*)InFile->Get("Physics_Results_Misc/t_alt_ttruth_result"));
  TH1F* tresHist =  (TH1F*)((TH1F*)InFile->Get("Physics_Results_Misc/taltres_result"));
  TH1F* MmissHist = (TH1F*)((TH1F*)InFile->Get("Physics_Results_Misc/Mmiss_result"));
  TH1F* MmissSqHist  = (TH1F*)((TH1F*)InFile->Get("Physics_Results_Misc/MmissSq_result"));
  TH1F* Mmiss_truth_Hist = (TH1F*)((TH1F*)InFile->Get("Physics_Results_Misc/Mmiss_truth_result"));
  TH1F* Mmiss_comp_Hist = (TH1F*)((TH1F*)InFile->Get("Physics_Results_Misc/Mmiss_Comp_result"));

  TH1F* piRes_p_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/piRes_p"));
  TH1F* piRes_px_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/piRes_p_{x}"));
  TH1F* piRes_py_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/piRes_p_{y}"));
  TH1F* piRes_pz_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/piRes_p_{z}"));
  TH1F* eRes_p_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/eRes_p"));
  TH1F* eRes_px_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/eRes_p_{x}"));
  TH1F* eRes_py_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/eRes_p_{y}"));
  TH1F* eRes_pz_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/eRes_p_{z}"));
  TH1F* nRes_p_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/nRes_p"));
  TH1F* nRes_px_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/nRes_p_{x}"));
  TH1F* nRes_py_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/nRes_p_{y}"));
  TH1F* nRes_pz_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/nRes_p_{z}"));
  TH1F* pmissDiff_p_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/pmissDiff_p"));
  TH1F* pmissDiff_px_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/pmissDiff_p_{x}"));
  TH1F* pmissDiff_py_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/pmissDiff_p_{y}"));
  TH1F* pmissDiff_pz_Hist = (TH1F*)((TH1F*)InFile->Get("Particle_Momenta_Resolution/pmissDiff_p_{z}"));
  
  for(Int_t A = 1; A <9; A++){
    tHists[A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results/t_cut_Result_Q2_%i",A)));
    ttruthHists[A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results/t_truth_thrown_Result_Q2_%i",A)));
    Q2Hists[A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results/Q2_cut_Result_Q2_%i",A)));
    WHists[A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results/W_cut_Result_Q2_%i",A)));
    Q2tHists[A-1] = (TH2F*)((TH2F*)InFile->Get(Form("Physics_Results/Q2_t_cut_Result_Q2_%i",A)));
    tvst_Q2_Hists[A-1] = (TH2F*)((TH2F*)InFile->Get(Form("Physics_Results/t_ttruth_Result_Q2_%i",A))); 
    taltvst_Q2_Hists[A-1] = (TH2F*)((TH2F*)InFile->Get(Form("Physics_Results/t_alt_ttruth_Result_Q2_%i",A))); 
  }
  
  // Draw reconstructed t distributions with cuts
  TCanvas  *c_Beijing_tRates = new TCanvas("c_Beijing_tRates", "t Binned Rates", 100, 0, 1000, 900);
  c_Beijing_tRates->Divide(3,3);
  for(Int_t A = 0; A < 8; A++){
    c_Beijing_tRates->cd(A+1);
    tHists[A]->SetStats(0); tHists[A]->SetLineWidth(2);
    tHists[A]->Draw("HISTERR");
  }

  c_Beijing_tRates->Print("Beijing_Plots/t_binned_rates_AllCuts.png"); c_Beijing_tRates->Print("Beijing_Plots/t_binned_rates_AllCuts.pdf");
  
  // Section below is for tech note/presentation plots
  
  // TCanvas *c_piXY = new TCanvas("c_piXY", "Pion XY Dist at HCal", 100, 0, 2560, 1920);
  // piXYHist->SetStats(0); piXYHist->SetTitle("#pi X vs Y at z =  375cm (HCal)"); piXYHist->GetXaxis()->SetTitle("x (cm)"); piXYHist->GetYaxis()->SetTitle("y (cm)"); piXYHist->Draw("COLZ"); c_piXY->Print("TechNote_Plots/piXY_5on100.png");
  // TCanvas *c_eXY = new TCanvas("c_eXY", "Electron XY Dist at EMCal", 100, 0, 2560, 1920);
  // eXYHist->SetStats(0); eXYHist->GetXaxis()->SetTitle("x (cm)"); eXYHist->GetYaxis()->SetTitle("y (cm)"); eXYHist->Draw("COLZ"); c_eXY->Print("TechNote_Plots/eXY_5on100.png");
  // TCanvas *c_nXY = new TCanvas("c_nXY", "Neutron XY Dist at ZDC HCal", 100, 0, 2560, 1920);
  // ZDCHist->SetStats(0); ZDCHist->SetTitle("n X vs Y at ZDC"); ZDCHist->GetXaxis()->SetTitle("x (cm)"); ZDCHist->GetYaxis()->SetTitle("y (cm)"); ZDCHist->Draw("COLZ");
  // ZDCHist->GetZaxis()->SetLabelSize(0.03); ZDCHist->GetZaxis()->SetLabelOffset(0); ZDCHist->Draw("COLZ"); c_nXY->Print("ECCE_NIM_Plots/nXY_5on100.png");
  // TCanvas *c_pipTh = new TCanvas("c_pipTh", "Pion Truth p vs Theta", 100, 0, 2560, 1920);
  // pipThetaTruthHist->SetStats(0); pipThetaTruthHist->GetXaxis()->SetTitle("#theta (Deg)"); pipThetaTruthHist->GetYaxis()->SetTitle("P (GeV/c)"); pipThetaTruthHist->Draw("COLZ"); c_pipTh->Print("TechNote_Plots/pipTh_5on100.png");
  // TCanvas *c_epTh = new TCanvas("c_epTh", "Electron Truth p vs Theta", 100, 0, 2560, 1920);
  // epThetaTruthHist->SetStats(0); epThetaTruthHist->GetXaxis()->SetTitle("#theta (Deg)"); epThetaTruthHist->GetYaxis()->SetTitle("P (GeV/c)"); epThetaTruthHist->Draw("COLZ"); c_epTh->Print("TechNote_Plots/epTh_5on100.png");
  // TCanvas *c_npTh = new TCanvas("c_npTh", "Neutron Truth p vs Theta", 100, 0, 2560, 1920);
  // npThetaTruthHist->SetStats(0); npThetaTruthHist->GetXaxis()->SetTitle("#theta (Deg)"); npThetaTruthHist->GetYaxis()->SetTitle("P (GeV/c)"); npThetaTruthHist->Draw("COLZ"); c_npTh->Print("TechNote_Plots/npTh_5on100.png");

  // TCanvas *c_tvstt = new TCanvas("c_tvstt", "t vs t truth", 100, 0, 2560, 1920); 
  // tvstHist->SetStats(0); tvstHist->SetTitle(""); tvstHist->Draw("COLZ"); c_tvstt->Print("ECCE_NIM_Plots/tvstt_5on100.png"); c_tvstt->Print("ECCE_NIM_Plots/tvstt_5on100.pdf");
  // TCanvas *c_taltvstt = new TCanvas("c_taltvstt", "t alt vs t truth", 100, 0, 2560, 1920);
  // t_altvstHist->SetStats(0); t_altvstHist->SetTitle(""); t_altvstHist->Draw("COLZ"); c_taltvstt->Print("ECCE_NIM_Plots/taltvstt_5on100.png"); c_taltvstt->Print("ECCE_NIM_Plots/taltvstt_5on100.pdf");
  // TCanvas *c_tvstt = new TCanvas("c_tvstt", "t vs t truth", 100, 0, 2560, 1920); 
  // tvst_Q2_Hists[4]->SetStats(0); tvst_Q2_Hists[4]->SetTitle(""); tvst_Q2_Hists[4]->Draw("COLZ"); c_tvstt->Print("ECCE_NIM_Plots/tvstt_5on100.png"); c_tvstt->Print("ECCE_NIM_Plots/tvstt_5on100.pdf");
  // TCanvas *c_taltvstt = new TCanvas("c_taltvstt", "t alt vs t truth", 100, 0, 2560, 1920);
  // taltvst_Q2_Hists[4]->SetStats(0); taltvst_Q2_Hists[4]->SetTitle(""); taltvst_Q2_Hists[4]->Draw("COLZ"); c_taltvstt->Print("ECCE_NIM_Plots/taltvstt_5on100.png"); c_taltvstt->Print("ECCE_NIM_Plots/taltvstt_5on100.pdf");
  
  TCanvas *c_Q2tEff = new TCanvas("c_Q2tEff", "Q2 vs t Detection Efficiency", 100, 0, 2560, 1920);
  Q2tEffHist_v2->SetStats(0); Q2tEffHist_v2->GetXaxis()->SetTitle("Q^{2} (GeV/c^{2})"); Q2tEffHist_v2->GetYaxis()->SetTitle("-t (GeV^{2})"); Q2tEffHist_v2->Draw("COLZ");
  c_Q2tEff->Print("Beijing_Plots/Q2tEff_5on100.png"); c_Q2tEff->Print("Beijing_Plots/Q2tEff_5on100.pdf");
  
  // This section is for SIDIS comparison plots.
  // I've chosen the Q2 = 17.5 bin to draw, I scale both histograms to peak at 1 for comparison (actual values are irelevant).
  TFile *SIDISFile = new TFile("LundTest_AllRuns_NoSmear_Analysis_v2.root");
  TH1F* SIDIS_pMiss_Dist = (TH1F*)((TH1F*)SIDISFile->Get("SIDIS_neutron_(Pmiss)_Distributions/Neutron p Dist Q2 = 17.500000")); // Histogram is named very poorly in this file
  TH1F* pn_Dist[8];
  for(Int_t A = 1; A <9; A++){
    pn_Dist[A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results_Cuts/pn_Result_Q2_%i",A)));
  }
  TH1F* CutDemos[3][4];
  for(Int_t A=1; A <5; A++){
    CutDemos[0][A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Cut_Analysis/t_result_cut%i_Low",A)));
    CutDemos[1][A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Cut_Analysis/t_result_cut%i_Mid",A)));
    CutDemos[2][A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Cut_Analysis/t_result_cut%i_High",A)));
  }
  
  TH1F* DEMP_pMiss_Dist = (TH1F*)((TH1F*)InFile->Get("Physics_Results_Cuts/pmiss_Result_Q2_4"));

  // Compairson of DEMP neutrons with DEMP "missing mass"
  //TCanvas *c_DEMPComp = new TCanvas("c_DEMPComp", "DEMP p_{n} vs p_{miss} Comparison", 100, 0, 2560, 1920);
  //pn_Dist[3]->SetStats(0); pn_Dist[3]->SetMinimum(0);
  //pn_Dist[3]->SetTitle("DEMP p_{n} vs p_{miss}, 15 < Q^{2} < 20"); pn_Dist[3]->GetXaxis()->SetTitle("DEMP p_{n}, p_{miss} (GeV/c)"); pn_Dist[3]->GetYaxis()->SetTitle("Rate(Hz)"); // Set titles
  pn_Dist[3]->GetXaxis()->SetRangeUser(40, 100);// Narrow down x range
  pn_Dist[3]->GetYaxis()->SetNdivisions(5,5,0); pn_Dist[3]->GetYaxis()->SetMaxDigits(3); pn_Dist[3]->GetYaxis()->SetLabelSize(0.04); // Tweak ticks and labelling on the Y axis
  //DEMP_pMiss_Dist->SetLineColor(2); pn_Dist[3]->Draw("HIST"); DEMP_pMiss_Dist->Draw("SAMEHIST");
  //TLegend *DEMPCompLegend = new TLegend(0.20, 0.60, 0.50, 0.75); // Define legend
  //DEMPCompLegend->AddEntry(pn_Dist[3], "DEMP Neutron Momentum", "l"); DEMPCompLegend->AddEntry(DEMP_pMiss_Dist, "DEMP Missing Momentum", "l"); // Add entries to legend
  //DEMPCompLegend->Draw("SAME");// Draw legend
  //c_DEMPComp->Print("ECCE_NIM_Plots/DEMP_Comp.png"); c_DEMPComp->Print("ECCE_NIM_Plots/DEMP_Comp.pdf"); // Print to file

  // Write SIDIS plots to root files too
  //TFile* DEMPComp_File = new TFile("ECCE_NIM_Plots/DEMP_Comp.root" , "RECREATE");
  //c_DEMPComp->Write(); DEMPCompLegend->Write(); pn_Dist[3]->Write(); DEMP_pMiss_Dist->Write();
  // DEMPComp_File->Close();
  
  // Scale histograms to 1 - arbitrary scaling to show comparison
  Double_t SIDISMax = SIDIS_pMiss_Dist->GetBinContent(SIDIS_pMiss_Dist->GetMaximumBin()); // Get max value from SIDIS dist
  Double_t DEMPMax = pn_Dist[3]->GetBinContent(pn_Dist[3]->GetMaximumBin()); // Get max value from DEMP dist
  Double_t DEMPMax_pMiss = DEMP_pMiss_Dist->GetBinContent(DEMP_pMiss_Dist->GetMaximumBin()); // Get max value from DEMP pMiss dist
  SIDIS_pMiss_Dist->Scale(1/SIDISMax); // Scale SIDIS hist to 1
  pn_Dist[3]->Scale(1/DEMPMax); // Scale DEMP hist to 1
  DEMP_pMiss_Dist->Scale(1/DEMPMax_pMiss); // Scale DEMP pMiss hist to 1
  
  // Scaled version
  //TCanvas *c_DEMPComp_Scaled = new TCanvas("c_DEMPComp_Scaled", "DEMP p_{n} vs p_{miss} Comparison", 100, 0, 2560, 1920);
  //pn_Dist[3]->SetStats(0); pn_Dist[3]->SetMinimum(0); // Force set minimum to 0 again
  //pn_Dist[3]->GetYaxis()->SetTitle("Arbitrary Scale"); pn_Dist[3]->GetYaxis()->SetTitleOffset(0.75); // Set y title
  //pn_Dist[3]->Draw("HIST"); DEMP_pMiss_Dist->Draw("SAMEHIST"); // Draw histograms
  //DEMPCompLegend->Draw("SAME"); // Draw legend
  //c_DEMPComp_Scaled->Print("ECCE_NIM_Plots/DEMP_Comp_Scaled.png"); c_DEMPComp_Scaled->Print("ECCE_NIM_Plots/DEMP_Comp_Scaled.pdf"); // Print to file

  // Comparison of SIDIS BG with DEMP neutrons
  // TCanvas *c_SIDISComp = new TCanvas("c_SIDISComp", "DEMP vs SIDIS Comparison", 100, 0, 2560, 1920);
  // pn_Dist[3]->SetStats(0); pn_Dist[3]->SetTitle("DEMP p_{n} vs SIDIS p_{miss}, 15 < Q^{2} < 20"); pn_Dist[3]->GetXaxis()->SetTitle("DEMP p_{n}, SIDIS p_{miss} (GeV/c)"); pn_Dist[3]->GetYaxis()->SetTitle("Arbitrary Scale"); // Set titles
  // pn_Dist[3]->SetMinimum(0.); SIDIS_pMiss_Dist->SetLineColor(2);
  // pn_Dist[3]->Draw("HIST"); SIDIS_pMiss_Dist->Draw("SAMEHIST"); // Set y axis to go from 0, set colour of SIDIS to red, draw both
  // TLegend *SIDISCompLegend = new TLegend(0.20, 0.60, 0.50, 0.75); // Define legend
  // SIDISCompLegend->AddEntry(pn_Dist[3], "DEMP Neutron Momentum", "l"); SIDISCompLegend->AddEntry(SIDIS_pMiss_Dist, "SIDIS Missing Momentum", "l"); // Add entries to legend
  // SIDISCompLegend->Draw("SAME");// Draw legend
  //c_SIDISComp->Print("ECCE_NIM_Plots/SIDIS_Comp.png"); c_SIDISComp->Print("ECCE_NIM_Plots/SIDIS_Comp.pdf"); // Print to file
  
  // Write the scaled files to root files
  //TFile* DEMPComp_Scaled_File = new TFile("ECCE_NIM_Plots/DEMP_Comp_Scaled.root", "RECREATE");
  //c_DEMPComp_Scaled->Write(); DEMPCompLegend->Write(); pn_Dist[3]->Write(); DEMP_pMiss_Dist->Write();
  //DEMPComp_Scaled_File->Close();
  // This sucks, but naming is stupid and there's only 8 histograms so just manually do it
  TH1F* SIDIS_Dists[8];
  SIDIS_Dists[0] = (TH1F*)((TH1F*)SIDISFile->Get("SIDIS_neutron_(Pmiss)_Distributions/Neutron p Dist Q2 = 6.250000"));
  SIDIS_Dists[1] = (TH1F*)((TH1F*)SIDISFile->Get("SIDIS_neutron_(Pmiss)_Distributions/Neutron p Dist Q2 = 8.750000"));
  TH1F* SIDIS_tmp1 = (TH1F*)((TH1F*)SIDISFile->Get("SIDIS_neutron_(Pmiss)_Distributions/Neutron p Dist Q2 = 11.250000"));
  TH1F* SIDIS_tmp2 = (TH1F*)((TH1F*)SIDISFile->Get("SIDIS_neutron_(Pmiss)_Distributions/Neutron p Dist Q2 = 13.750000"));
  SIDIS_tmp1->Add(SIDIS_tmp2); // Add these two histograms together to get the 12.5 bin
  SIDIS_Dists[2] = SIDIS_tmp1;
  SIDIS_Dists[3] = (TH1F*)((TH1F*)SIDISFile->Get("SIDIS_neutron_(Pmiss)_Distributions/Neutron p Dist Q2 = 17.500000"));
  SIDIS_Dists[4] = (TH1F*)((TH1F*)SIDISFile->Get("SIDIS_neutron_(Pmiss)_Distributions/Neutron p Dist Q2 = 22.500000"));
  SIDIS_Dists[5] = (TH1F*)((TH1F*)SIDISFile->Get("SIDIS_neutron_(Pmiss)_Distributions/Neutron p Dist Q2 = 27.500000"));
  SIDIS_Dists[6] = (TH1F*)((TH1F*)SIDISFile->Get("SIDIS_neutron_(Pmiss)_Distributions/Neutron p Dist Q2 = 32.500000"));
  SIDIS_Dists[7] = (TH1F*)((TH1F*)SIDISFile->Get("SIDIS_neutron_(Pmiss)_Distributions/Neutron p Dist Q2 = 37.500000"));

  // Plots of Garth's Beijing talk
  // Manually grab all of the inidividual SIDIS plots
  // Draw reconstructed t distributions with cuts
  Double_t Q2Vals[2][8]= {{5, 7.5, 10, 15, 20, 25, 30, 35}, {7.5, 10, 15, 20, 25, 30, 35, 40}};
  Double_t PMissCutVals[8]={96.0, 93.5, 91.0, 87.0, 83.0, 80.0, 77.5, 75.0};
  TLine *PMissCuts[8];
  TCanvas  *c_BeijingDEMPComp = new TCanvas("c_BeijingDEMPComp", "Beijing DEMP vs SIDIS Comparisons", 100, 0, 1000, 900);
  c_BeijingDEMPComp->Divide(3,3);
  for(Int_t A = 0; A < 8; A++){
    c_BeijingDEMPComp->cd(A+1);
    pn_Dist[A]->Rebin(2);
    DEMPMax = pn_Dist[A]->GetBinContent(pn_Dist[A]->GetMaximumBin()); pn_Dist[A]->Scale(1/DEMPMax);
    pn_Dist[A]->SetStats(0); pn_Dist[A]->SetTitle(Form("DEMP p_{n} vs SIDIS p_{miss}, %0.2f < Q^{2} < %0.2f", Q2Vals[0][A], Q2Vals[1][A])); pn_Dist[A]->GetXaxis()->SetTitle("DEMP p_{n}, SIDIS p_{miss} (GeV/c)"); pn_Dist[A]->GetYaxis()->SetTitle("Arbitrary Scale"); // Set titles
    pn_Dist[A]->GetXaxis()->SetRangeUser(40, 100);// Narrow down x range
    pn_Dist[A]->GetYaxis()->SetNdivisions(5,5,0); pn_Dist[A]->GetYaxis()->SetMaxDigits(3); pn_Dist[A]->GetYaxis()->SetLabelSize(0.04); // Tweak ticks and labelling on the Y axis
    pn_Dist[A]->SetLineWidth(2);
    SIDISMax = SIDIS_Dists[A]->GetBinContent(SIDIS_Dists[A]->GetMaximumBin()); SIDIS_Dists[A]->Scale(1/SIDISMax); // Scale SIDIS hist to 1
    SIDIS_Dists[A]->SetLineColor(2); SIDIS_Dists[A]->SetLineWidth(2);
    PMissCuts[A]= new TLine(PMissCutVals[A],gPad->GetUymin(),PMissCutVals[A],gPad->GetUymax()); PMissCuts[A]->SetLineColor(8); PMissCuts[A]->SetLineWidth(3);
    pn_Dist[A]->Draw("HIST"); SIDIS_Dists[A]->Draw("HISTSAME"); PMissCuts[A]->Draw("SAME");
  }
  c_BeijingDEMPComp->cd(9);
  TLegend *SIDISCompLegend2 = new TLegend(0.10, 0.10, 0.90, 0.9); // Define legend
  SIDISCompLegend2->AddEntry(pn_Dist[3], "DEMP p_{n}", "l"); SIDISCompLegend2->AddEntry(SIDIS_Dists[3], "SIDIS p_{Miss}", "l"); SIDISCompLegend2->AddEntry(PMissCuts[0], "p cut", "l"); // Add entries to legend
  SIDISCompLegend2->Draw("SAME");// Draw legend

  c_BeijingDEMPComp->Print("Beijing_Plots/DEMP_SIDIS_Comp.png"); c_BeijingDEMPComp->Print("Beijing_Plots/DEMP_SIDIS_Comp.pdf");

  TCanvas *c_BeijingCutDemoLow = new TCanvas("c_BeijingCutDemoLow", "Beijing DEMP Cut Demonstrations - Low Q2", 100, 0, 1000, 900);
  CutDemos[0][1]->SetStats(0); CutDemos[0][1]->SetMinimum(0); CutDemos[0][1]->SetTitle("-t Dist 7.5 < Q^{2} < 10"); CutDemos[0][1]->GetYaxis()->SetTitle("Rate (Hz)");
  CutDemos[0][1]->SetLineColor(6); CutDemos[0][1]->SetLineWidth(3); CutDemos[0][1]->Draw("HISTERR");
  CutDemos[0][2]->SetLineColor(4); CutDemos[0][2]->SetLineWidth(3); CutDemos[0][2]->Draw("SAMEHISTERR");
  CutDemos[0][3]->SetLineColor(94);CutDemos[0][3]->SetLineWidth(3); CutDemos[0][3]->Draw("SAMEHISTERR");
  TLegend *BeijingCutDemoLegend = new TLegend(0.6, 0.65, 0.95, 0.9); // Define legend
  BeijingCutDemoLegend->AddEntry(CutDemos[0][1], "-t and #theta_{n} cuts", "l"); BeijingCutDemoLegend->AddEntry(CutDemos[0][2], "#splitline{-t, #theta_{n}, #theta_{ZDC}-#theta_{pMiss}}{and #phi_{ZDC}-#phi_{pMiss} cuts}", "l"); // Add entries to legend
  BeijingCutDemoLegend->AddEntry(CutDemos[0][3], "#splitline{-t, #theta_{n}, #theta_{ZDC}-#theta_{pMiss}, #phi_{ZDC}-#phi_{pMiss}}{and p_{miss} cuts}", "l");
  BeijingCutDemoLegend->Draw("SAME");// Draw legend
  c_BeijingCutDemoLow->Print("Beijing_Plots/CutDemo_Low_Q2.png"); c_BeijingCutDemoLow->Print("Beijing_Plots/CutDemo_Low_Q2.pdf");

  TCanvas *c_BeijingCutDemoMid = new TCanvas("c_BeijingCutDemoMid", "Beijing DEMP Cut Demonstrations - Mid Q2", 100, 0, 1000, 900);
  CutDemos[1][1]->SetStats(0); CutDemos[1][1]->SetMinimum(0); CutDemos[1][1]->SetTitle("-t Dist 15 < Q^{2} < 20"); CutDemos[1][1]->GetYaxis()->SetTitle("Rate (Hz)");
  CutDemos[1][1]->SetLineColor(6); CutDemos[1][1]->SetLineWidth(3); CutDemos[1][1]->GetYaxis()->SetTitleOffset(1.15); CutDemos[1][1]->Draw("HISTERR");
  CutDemos[1][2]->SetLineColor(4); CutDemos[1][2]->SetLineWidth(3); CutDemos[1][2]->Draw("SAMEHISTERR");
  CutDemos[1][3]->SetLineColor(94); CutDemos[1][3]->SetLineWidth(3); CutDemos[1][3]->Draw("SAMEHISTERR");
  BeijingCutDemoLegend->Draw("SAME");// Draw legend
  c_BeijingCutDemoMid->Print("Beijing_Plots/CutDemo_Mid_Q2.png"); c_BeijingCutDemoMid->Print("Beijing_Plots/CutDemo_Mid_Q2.pdf");

  TCanvas *c_BeijingCutDemoHigh = new TCanvas("c_BeijingCutDemoHigh", "Beijing DEMP Cut Demonstrations - High Q2", 100, 0, 1000, 900);
  CutDemos[2][1]->SetStats(0); CutDemos[2][1]->SetMinimum(0); CutDemos[2][1]->SetTitle("-t Dist 25 < Q^{2} < 30"); CutDemos[2][1]->GetYaxis()->SetTitle("Rate (Hz)");
  CutDemos[2][1]->SetLineColor(6); CutDemos[2][1]->SetLineWidth(3); CutDemos[2][1]->Draw("HISTERR");
  CutDemos[2][2]->SetLineColor(4); CutDemos[2][2]->SetLineWidth(3); CutDemos[2][2]->Draw("SAMEHISTERR");
  CutDemos[2][3]->SetLineColor(94); CutDemos[2][3]->SetLineWidth(3); CutDemos[2][3]->Draw("SAMEHISTERR");
  BeijingCutDemoLegend->Draw("SAME");// Draw legend
  c_BeijingCutDemoHigh->Print("Beijing_Plots/CutDemo_High_Q2.png"); c_BeijingCutDemoHigh->Print("Beijing_Plots/CutDemo_High_Q2.pdf");

  TH1F* SIDIS_t_Dist_Beijing = (TH1F*)((TH1F*)SIDISFile->Get("SIDIS_Kinematic_Quantity_Distributions/SIDIS t Distribution")); // Histogram is named very poorly in this file
  TH1F* DEMP_t_Dist_Beijing = (TH1F*)((TH1F*)InFile->Get("Kinematics_Info/t_alt_Dist"));
  TH1F* SIDIS_pMiss_Dist_Beijing = (TH1F*)((TH1F*)SIDISFile->Get("SIDIS_neutron_(Pmiss)_Distributions/PMiss Distribution")); // Histogram is named very poorly in this file
  TH1F* DEMP_pMiss_Dist_Beijing = (TH1F*)((TH1F*)InFile->Get("PMiss_Info/pmiss_p"));

  TCanvas *c_BeijingSIDISt = new TCanvas("c_BeijingSIDISt", "Beijing SIDIS t Distribution", 100, 0, 1000, 900);
  SIDISMax = SIDIS_t_Dist_Beijing->GetBinContent(SIDIS_t_Dist_Beijing->GetMaximumBin()); SIDIS_t_Dist_Beijing->Scale(1/SIDISMax); // Scale SIDIS hist to 1
  SIDIS_t_Dist_Beijing->Rebin(4);
  SIDIS_t_Dist_Beijing->SetLineWidth(2);
  SIDIS_t_Dist_Beijing->SetTitle("SIDIS -t Distribution");
  SIDIS_t_Dist_Beijing->GetXaxis()->SetTitle("-t (GeV^{2})");
  SIDIS_t_Dist_Beijing->GetYaxis()->SetTitle("Arbitrary Scale");
  SIDIS_t_Dist_Beijing->Draw("HIST");
  c_BeijingSIDISt->Print("Beijing_Plots/SIDIS_t_Dist.png"); c_BeijingSIDISt->Print("Beijing_Plots/SIDIS_t_Dist.pdf");
  
  TCanvas *c_BeijingDEMPt = new TCanvas("c_BeijingDEMPt", "Beijing DEMP t Distribution", 100, 0, 1000, 900);
  DEMPMax = DEMP_t_Dist_Beijing->GetBinContent(DEMP_t_Dist_Beijing->GetMaximumBin()); DEMP_t_Dist_Beijing->Scale(1/DEMPMax); // Scale DEMP hist to 1
  DEMP_t_Dist_Beijing->SetLineWidth(2);
  DEMP_t_Dist_Beijing->SetTitle("DEMP -t Distribution"); DEMP_t_Dist_Beijing->SetStats(0);
  DEMP_t_Dist_Beijing->GetXaxis()->SetTitle("-t (GeV^{2})");
  DEMP_t_Dist_Beijing->GetYaxis()->SetTitle("Arbitrary Scale");
  DEMP_t_Dist_Beijing->Draw("HIST");
  c_BeijingDEMPt->Print("Beijing_Plots/DEMP_t_Dist.png"); c_BeijingDEMPt->Print("Beijing_Plots/DEMP_t_Dist.pdf");

  TCanvas *c_Beijing_tComp = new TCanvas("c_Beijing_tComp", "Beijing DEMP vs SIDIS t Distribution", 100, 0, 1000, 900);
  DEMP_t_Dist_Beijing->Draw("HIST");
  SIDIS_t_Dist_Beijing->SetLineColor(2); SIDIS_t_Dist_Beijing->Draw("HISTSAME");
  TLegend *Beijing_tDist_Legend = new TLegend(0.6, 0.6, 0.8, 0.8); // Define legend
  Beijing_tDist_Legend->AddEntry(DEMP_t_Dist_Beijing, "DEMP -t", "l"); Beijing_tDist_Legend->AddEntry(SIDIS_t_Dist_Beijing, "SIDIS -t", "l");
  Beijing_tDist_Legend->Draw("SAME");// Draw legend
  c_Beijing_tComp->Print("Beijing_Plots/t_Dist_Comp.png"); c_Beijing_tComp->Print("Beijing_Plots/t_Dist_Comp.pdf");
  
  TCanvas *c_BeijingSIDISpMiss = new TCanvas("c_BeijingSIDISpMiss", "Beijing SIDIS pMiss Distribution", 100, 0, 1000, 900);
  SIDISMax = SIDIS_pMiss_Dist_Beijing->GetBinContent(SIDIS_pMiss_Dist_Beijing->GetMaximumBin()); SIDIS_pMiss_Dist_Beijing->Scale(1/SIDISMax); // Scale SIDIS hist to 1
  SIDIS_pMiss_Dist_Beijing->GetXaxis()->SetRangeUser(40, 110);
  SIDIS_pMiss_Dist_Beijing->SetLineWidth(2);
  SIDIS_pMiss_Dist_Beijing->SetTitle("SIDIS p_{Miss} Distribution");
  SIDIS_pMiss_Dist_Beijing->GetXaxis()->SetTitle("p_{Miss} (GeV)");
  SIDIS_pMiss_Dist_Beijing->GetYaxis()->SetTitle("Arbitrary Scale");
  SIDIS_pMiss_Dist_Beijing->Draw("HIST");
  c_BeijingSIDISpMiss->Print("Beijing_Plots/SIDIS_pMiss_Dist.png"); c_BeijingSIDISpMiss->Print("Beijing_Plots/SIDIS_pMiss_Dist.pdf");
  
  TCanvas *c_BeijingDEMPpMiss = new TCanvas("c_BeijingDEMPpMiss", "Beijing DEMP pMiss Distribution", 100, 0, 1000, 900);
  DEMPMax = DEMP_pMiss_Dist_Beijing->GetBinContent(DEMP_pMiss_Dist_Beijing->GetMaximumBin()); DEMP_pMiss_Dist_Beijing->Scale(1/DEMPMax); // Scale DEMP hist to 1
  DEMP_pMiss_Dist_Beijing->GetXaxis()->SetRangeUser(40, 100);
  DEMP_pMiss_Dist_Beijing->SetLineWidth(2);
  DEMP_pMiss_Dist_Beijing->SetTitle("DEMP p_{Miss} Distribution"); DEMP_pMiss_Dist_Beijing->SetStats(0);
  DEMP_pMiss_Dist_Beijing->GetXaxis()->SetTitle("p_{Miss} (GeV)");
  DEMP_pMiss_Dist_Beijing->GetYaxis()->SetTitle("Arbitrary Scale");
  DEMP_pMiss_Dist_Beijing->Draw("HIST");
  c_BeijingDEMPpMiss->Print("Beijing_Plots/DEMP_pMiss_Dist.png"); c_BeijingDEMPpMiss->Print("Beijing_Plots/DEMP_pMiss_Dist.pdf");

  TCanvas *c_Beijing_pMissComp = new TCanvas("c_Beijing_pMissComp", "Beijing DEMP vs SIDIS t Distribution", 100, 0, 1000, 900);
  DEMP_pMiss_Dist_Beijing->Draw("HIST");
  SIDIS_pMiss_Dist_Beijing->SetLineColor(2); SIDIS_pMiss_Dist_Beijing->Draw("HISTSAME");
  TLegend *Beijing_pMissDist_Legend = new TLegend(0.2, 0.6, 0.4, 0.8); // Define legend
  Beijing_pMissDist_Legend->AddEntry(DEMP_pMiss_Dist_Beijing, "DEMP p_{Miss}", "l"); Beijing_pMissDist_Legend->AddEntry(SIDIS_pMiss_Dist_Beijing, "SIDIS p_{Miss}", "l");
  Beijing_pMissDist_Legend->Draw("SAME");// Draw legend
  c_Beijing_pMissComp->Print("Beijing_Plots/pMiss_Dist_Comp.png"); c_Beijing_pMissComp->Print("Beijing_Plots/pMiss_Dist_Comp.pdf");

  // 20/04/22
  // SIDIS data file had a different crossing angle, so comparison of theta distributions after the -t cut would be misleading at best
  // Regenerating this plot would take a fair bit of time too.
  // Can make a normal theta n plot?

  TH1F *DEMP_nTheta_tCut_Beijing = (TH1F*)((TH1F*)InFile->Get("Cut_Analysis/nTheta_tCut"));
  TCanvas *c_Beijing_nTheta = new TCanvas("c_Beijing_nTheta", "Beijing nTheta Distribution", 100, 0, 1000, 900);
  DEMP_nTheta_tCut_Beijing->SetStats(0);
  DEMP_nTheta_tCut_Beijing->SetLineWidth(2);
  DEMP_nTheta_tCut_Beijing->SetTitle("DEMP #theta_{n} Distribution (-t < 0.4 GeV^{2})");
  DEMP_nTheta_tCut_Beijing->GetXaxis()->SetTitle("#theta_{n} (Deg)"); DEMP_nTheta_tCut_Beijing->GetXaxis()->SetRangeUser(0, 2.4);
  DEMP_nTheta_tCut_Beijing->GetYaxis()->SetTitle("Rate(Hz)");
  DEMP_nTheta_tCut_Beijing->Draw("HIST");
  c_Beijing_nTheta->Print("Beijing_Plots/nTheta_Dist.png"); c_Beijing_nTheta->Print("Beijing_Plots/nTheta_Dist.pdf");

  TH2F* DEMP_eThetaP_Beijing = (TH2F*)((TH2F*)InFile->Get("Scattered_Electron_Info/eTrack_pTheta"));
  TH2F* DEMP_piThetaP_Beijing = (TH2F*)((TH2F*)InFile->Get("Pion_Info/piTrack_pTheta"));
  TH2F* DEMP_nThetaP_Beijing = (TH2F*)((TH2F*)InFile->Get("ZDC_Neutron_Info/nTrack_pTheta"));
  TH2F* DEMP_eThetaP_Truth_Beijing = (TH2F*)((TH2F*)InFile->Get("Scattered_Electron_Truth_Info/eTrack_pTheta_Truth"));
  TH2F* DEMP_piThetaP_Truth_Beijing = (TH2F*)((TH2F*)InFile->Get("Pion_Truth_Info/piTrack_pTheta_Truth"));
  TH2F* DEMP_nThetaP_Truth_Beijing = (TH2F*)((TH2F*)InFile->Get("Neutron_Truth_Info/nTrack_pTheta_Truth"));

  TCanvas *c_Beijing_eThetaP = new TCanvas("c_Beijing_eThetaP", "Beijing e' Theta vs P Distribution", 100, 0, 1000, 900);
  DEMP_eThetaP_Beijing->SetStats(0); DEMP_eThetaP_Beijing->GetXaxis()->SetRangeUser(110, 160);
  DEMP_eThetaP_Beijing->Draw("COLZ");
  c_Beijing_eThetaP->Print("Beijing_Plots/eThetaP_Dist.png"); c_Beijing_eThetaP->Print("Beijing_Plots/eThetaP_Dist.pdf");
  
  TCanvas *c_Beijing_piThetaP = new TCanvas("c_Beijing_piThetaP", "Beijing pi Theta vs P Distribution", 100, 0, 1000, 900);
  DEMP_piThetaP_Beijing->SetStats(0); DEMP_piThetaP_Beijing->GetXaxis()->SetRangeUser(0, 40);
  DEMP_piThetaP_Beijing->Draw("COLZ");
  c_Beijing_piThetaP->Print("Beijing_Plots/piThetaP_Dist.png"); c_Beijing_piThetaP->Print("Beijing_Plots/piThetaP_Dist.pdf");
  
  TCanvas *c_Beijing_nThetaP = new TCanvas("c_Beijing_nThetaP", "Beijing n Theta vs P Distribution", 100, 0, 1000, 900);
  DEMP_nThetaP_Beijing->SetStats(0);
  DEMP_nThetaP_Beijing->Draw("COLZ");
  c_Beijing_nThetaP->Print("Beijing_Plots/nThetaP_Dist.png"); c_Beijing_nThetaP->Print("Beijing_Plots/nThetaP_Dist.pdf");

  TCanvas *c_Beijing_eThetaP_Truth = new TCanvas("c_Beijing_eThetaP_Truth", "Beijing e' Theta vs P Distribution (Truth Info)", 100, 0, 1000, 900);
  DEMP_eThetaP_Truth_Beijing->SetStats(0); DEMP_eThetaP_Truth_Beijing->GetXaxis()->SetRangeUser(110, 160);
  DEMP_eThetaP_Truth_Beijing->Draw("COLZ");
  c_Beijing_eThetaP_Truth->Print("Beijing_Plots/eThetaP_Truth_Dist.png"); c_Beijing_eThetaP_Truth->Print("Beijing_Plots/eThetaP_Truth_Dist.pdf");
  
  TCanvas *c_Beijing_piThetaP_Truth = new TCanvas("c_Beijing_piThetaP_Truth", "Beijing pi Theta vs P Distribution (Truth Info)", 100, 0, 1000, 900);
  DEMP_piThetaP_Truth_Beijing->SetStats(0); DEMP_piThetaP_Truth_Beijing->GetXaxis()->SetRangeUser(0, 40);
  DEMP_piThetaP_Truth_Beijing->Draw("COLZ");
  c_Beijing_piThetaP_Truth->Print("Beijing_Plots/piThetaP_Truth_Dist.png"); c_Beijing_piThetaP_Truth->Print("Beijing_Plots/piThetaP_Truth_Dist.pdf");
  
  TCanvas *c_Beijing_nThetaP_Truth = new TCanvas("c_Beijing_nThetaP_Truth", "Beijing n Theta vs P Distribution (Truth Info)", 100, 0, 1000, 900);
  DEMP_nThetaP_Truth_Beijing->SetStats(0); DEMP_nThetaP_Truth_Beijing->GetXaxis()->SetRangeUser(0, 2.5);
  DEMP_nThetaP_Truth_Beijing->Draw("COLZ");
  c_Beijing_nThetaP_Truth->Print("Beijing_Plots/nThetaP_Truth_Dist.png"); c_Beijing_nThetaP_Truth->Print("Beijing_Plots/nThetaP_Truth_Dist.pdf");

  // t-Resolution plots for the 4th Q2 bin (15-20)
  TH1F* tRes_Dist_Beijing = (TH1F*)((TH1F*)InFile->Get("t_Resolution/t_Resolution_Q2_4"));
  TH1F* taltRes_ZDC_Dist_Beijing = (TH1F*)((TH1F*)InFile->Get("t_Resolution/talt_Resolution_ZDC_Q2_4"));
  TH1F* taltRes_pMiss_Dist_Beijing = (TH1F*)((TH1F*)InFile->Get("t_Resolution/talt_Resolution_pMiss_Q2_4"));
  
  TCanvas *c_Beijing_tRes = new TCanvas("c_Beijing_tRes", "Beijing t Resolution", 100, 0, 1000, 900);
  tRes_Dist_Beijing->GetYaxis()->SetTitleOffset(1.2);
  tRes_Dist_Beijing->SetLineWidth(2);
  tRes_Dist_Beijing->Rebin(3);
  tRes_Dist_Beijing->SetStats(0); tRes_Dist_Beijing->Draw("HIST");
  c_Beijing_tRes->Print("Beijing_Plots/t_Resolution.pdf"); c_Beijing_tRes->Print("Beijing_Plots/t_Resolution.png");
   
  TCanvas *c_Beijing_taltRes_ZDC = new TCanvas("c_Beijing_taltRes_ZDC", "Beijing talt Resolution ZDC", 100, 0, 1000, 900);
  taltRes_ZDC_Dist_Beijing->GetYaxis()->SetTitleOffset(1.5);
  taltRes_ZDC_Dist_Beijing->SetLineWidth(2);
  taltRes_ZDC_Dist_Beijing->SetStats(0); taltRes_ZDC_Dist_Beijing->Draw("HIST");
  c_Beijing_taltRes_ZDC->Print("Beijing_Plots/talt_Resolution_ZDC.pdf"); c_Beijing_taltRes_ZDC->Print("Beijing_Plots/talt_Resolution_ZDC.png");

  TCanvas *c_Beijing_taltRes_pMiss = new TCanvas("c_Beijing_taltRes_pMiss", "Beijing talt Resolution pMiss", 100, 0, 1000, 900);
  taltRes_pMiss_Dist_Beijing->GetYaxis()->SetTitleOffset(1.2);
  taltRes_pMiss_Dist_Beijing->SetLineWidth(2);
  taltRes_pMiss_Dist_Beijing->SetStats(0); taltRes_pMiss_Dist_Beijing->Draw("HIST");
  c_Beijing_taltRes_pMiss->Print("Beijing_Plots/talt_Resolution_pMiss.pdf"); c_Beijing_taltRes_pMiss->Print("Beijing_Plots/talt_Resolution_pMiss.png");

  TCanvas *c_Beijing_taltRes_Comp = new TCanvas("c_Beijing_taltRes_Comp", "Beijing talt Resolution Comparison", 100, 0, 1000, 900);
  taltRes_pMiss_Dist_Beijing->SetTitle("t_{alt}-t_{truth} Dist, 15 < Q^{2} < 20");
  taltRes_pMiss_Dist_Beijing->SetLineColor(4); taltRes_pMiss_Dist_Beijing->Draw("HIST");
  taltRes_ZDC_Dist_Beijing->SetLineColor(2); taltRes_ZDC_Dist_Beijing->Draw("HISTSAME");
  TLegend *Beijing_taltResDist_Legend = new TLegend(0.2, 0.6, 0.4, 0.8); // Define legend
  Beijing_taltResDist_Legend->AddEntry(taltRes_pMiss_Dist_Beijing, "Adjusted p_{Miss}", "l"); Beijing_taltResDist_Legend->AddEntry(taltRes_ZDC_Dist_Beijing, "ZDC Info Only", "l");
  Beijing_taltResDist_Legend->Draw("SAME");// Draw legend
  c_Beijing_taltRes_Comp->Print("Beijing_Plots/talt_Resolution_Comp.pdf"); c_Beijing_taltRes_Comp->Print("Beijing_Plots/talt_Resolution_Comp.png");
  
  SIDISFile->Close(); // Close SIDIS file
  
  InFile->Close();
}
