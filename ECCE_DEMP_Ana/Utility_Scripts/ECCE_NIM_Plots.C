// 26/04/22 - Stephen JD Kay, University of Regina

// A script to extract and plot some stuff for the ECCE NIM paper
// Execute via root -l 'ECCE_NIM_Plots.C("INPUT_FILE.root")'
#define ECCE_NIM_Plots_cxx
          
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

void ECCE_NIM_Plots(string InFilename = ""){

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
    
  // Section below is for tech note/presentation plots
  
  // TCanvas *c_piXY = new TCanvas("c_piXY", "Pion XY Dist at HCal", 100, 0, 2560, 1920);
  // piXYHist->SetStats(0); piXYHist->SetTitle("#pi X vs Y at z =  375cm (HCal)"); piXYHist->GetXaxis()->SetTitle("x (cm)"); piXYHist->GetYaxis()->SetTitle("y (cm)"); piXYHist->Draw("COLZ"); c_piXY->Print("TechNote_Plots/piXY_5on100.png");
  // TCanvas *c_eXY = new TCanvas("c_eXY", "Electron XY Dist at EMCal", 100, 0, 2560, 1920);
  // eXYHist->SetStats(0); eXYHist->GetXaxis()->SetTitle("x (cm)"); eXYHist->GetYaxis()->SetTitle("y (cm)"); eXYHist->Draw("COLZ"); c_eXY->Print("TechNote_Plots/eXY_5on100.png");
  TCanvas *c_nXY = new TCanvas("c_nXY", "Neutron XY Dist at ZDC HCal", 100, 0, 2560, 1920);
  ZDCHist->SetStats(0); ZDCHist->SetTitle(""); ZDCHist->GetXaxis()->SetTitle("x (cm)"); ZDCHist->GetYaxis()->SetTitle("y (cm)"); ZDCHist->Draw("COLZ");
  ZDCHist->GetZaxis()->SetLabelSize(0.03); ZDCHist->GetZaxis()->SetLabelOffset(0); ZDCHist->Draw("COLZ"); c_nXY->Print("ECCE_NIM_Plots/nXY_5on100.png"); c_nXY->Print("ECCE_NIM_Plots/nXY_5on100.pdf");
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
  TCanvas *c_tvstt = new TCanvas("c_tvstt", "t vs t truth", 100, 0, 2560, 1920); 
  tvstHist->SetStats(0); tvstHist->SetTitle(""); tvstHist->Draw("COLZ"); c_tvstt->Print("ECCE_NIM_Plots/tvstt_5on100.png"); c_tvstt->Print("ECCE_NIM_Plots/tvstt_5on100.pdf");
  TCanvas *c_taltvstt = new TCanvas("c_taltvstt", "t alt vs t truth", 100, 0, 2560, 1920);
  t_altvstHist->SetStats(0); t_altvstHist->GetXaxis()->SetRangeUser(0,0.4); t_altvstHist->SetTitle(""); t_altvstHist ->Draw("COLZ"); c_taltvstt->Print("ECCE_NIM_Plots/taltvstt_5on100.png"); c_taltvstt->Print("ECCE_NIM_Plots/taltvstt_5on100.pdf");
  
  TCanvas *c_Q2tEff = new TCanvas("c_Q2tEff", "Q2 vs t Detection Efficiency", 100, 0, 2560, 1920);
  Q2tEffHist_v2->SetStats(0); Q2tEffHist_v2->SetTitle(""); Q2tEffHist_v2->GetXaxis()->SetTitle("Q^{2} (GeV/c^{2})"); Q2tEffHist_v2->GetYaxis()->SetTitle("-t (GeV^{2})"); Q2tEffHist_v2->Draw("COLZ");
  c_Q2tEff->Print("ECCE_NIM_Plots/Q2tEff_5on100.png"); c_Q2tEff->Print("ECCE_NIM_Plots/Q2tEff_5on100.pdf");
  // TCanvas *c_pmissDiffp = new TCanvas("c_pmissDiffp", "Pmiss vs Neutron P Truth", 100, 0, 2560, 1920);
  // pmissDiff_p_Hist->SetStats(0); pmissDiff_p_Hist->Draw("HISTERR"); c_pmissDiffp->Print("TechNote_Plots/pmissDiffp_5on100.png");
  // TCanvas *c_pmissDiffpx = new TCanvas("c_pmissDiffpx", "Pmiss vs Neutron P Truth", 100, 0, 2560, 1920);
  // pmissDiff_px_Hist->SetStats(0); pmissDiff_px_Hist->Draw("HISTERR"); c_pmissDiffpx->Print("TechNote_Plots/pmissDiffpx_5on100.png");
  // TCanvas *c_pmissDiffpy = new TCanvas("c_pmissDiffpy", "Pmiss vs Neutron P Truth", 100, 0, 2560, 1920);
  // pmissDiff_py_Hist->SetStats(0); pmissDiff_py_Hist->Draw("HISTERR"); c_pmissDiffpy->Print("TechNote_Plots/pmissDiffpy_5on100.png");
  // TCanvas *c_pmissDiffpz = new TCanvas("c_pmissDiffpz", "Pmiss vs Neutron P Truth", 100, 0, 2560, 1920);
  // pmissDiff_pz_Hist->SetStats(0); pmissDiff_pz_Hist->Draw("HISTERR"); c_pmissDiffpz->Print("TechNote_Plots/pmissDiffpz_5on100.png");

  TString Q2LabelsInfo[8]={"5 < Q^{2} < 7.5", "7.5 < Q^{2} < 10", "10 < Q^{2} < 15", "15 < Q^{2} < 20", "20 < Q^{2} < 25", "25 < Q^{2} < 30", "30 < Q^{2} < 35", "35 < Q^{2} < 40"}; 
  TLegend *Q2Labels[8];
  
  TCanvas *c_Q2tbinned_rates_v1 = new TCanvas("c_Q2tbinned_rates_v1", "Rates per -t and Q2 Bin", 100, 0, 2560, 1920);
  c_Q2tbinned_rates_v1->Divide(3,3);
  for(Int_t A = 0; A < 8; A++){
    c_Q2tbinned_rates_v1->cd(A+1);
    Q2Labels[A] = new TLegend(0.65, 0.65, 0.85, 0.85);
    Q2Labels[A]->SetHeader(Q2LabelsInfo[A], "C");
    Q2Labels[A]->SetBorderSize(0); Q2Labels[A]->SetFillStyle(0); // Remove the borders and make the legend box transparent
    Q2Labels[A]->SetTextSize(0.06);
    tHists[A]->GetXaxis()->SetLabelSize(0.05); tHists[A]->GetXaxis()->SetTitleOffset(1.20);
    tHists[A]->GetYaxis()->SetNdivisions(5,5,0); tHists[A]->GetYaxis()->SetMaxDigits(2); tHists[A]->GetYaxis()->SetLabelSize(0.05); tHists[A]->GetYaxis()->SetTitleOffset(1.20); // Tweak ticks and labelling on the Y axis
    tHists[A]->SetTitle(""); tHists[A]->SetLineWidth(2);
    //tHists[A]->GetXaxis()->SetTitle(""); tHists[A]->GetYaxis()->SetTitle("");// Optional - Disable titles, have to then be written on manually to the .png
    tHists[A]->SetStats(0); tHists[A]->Draw("HISTERR"); Q2Labels[A]->Draw("SAME");
  }

  c_Q2tbinned_rates_v1->cd(9);  
  InfoDump.DrawLatex(.2,.8,"L = 10^{34} cm^{-2}s^{-1}");
  InfoDump.DrawLatex(.2,.7,"assumed in rate");
  InfoDump.DrawLatex(.2,.6,"calculation");
  c_Q2tbinned_rates_v1->Print("ECCE_NIM_Plots/Q2t_rates_v1_5on100.png"); c_Q2tbinned_rates_v1->Print("ECCE_NIM_Plots/Q2t_rates_v1_5on100.pdf");
   
  TCanvas *c_nResp = new TCanvas("c_nResp", "Neutron p Resolution", 100, 0, 2560, 1920);
  nRes_p_Hist->SetStats(0); nRes_p_Hist->SetTitle(""); nRes_p_Hist->GetXaxis()->SetRangeUser(-3, 3); nRes_p_Hist->GetXaxis()->SetTitle("#frac{#Delta p_{n}}{p_{ntruth}} (%)"); nRes_p_Hist->GetYaxis()->SetTitle("Rate (Hz)");
  nRes_p_Hist->GetXaxis()->SetTitleOffset(1.05); nRes_p_Hist->GetXaxis()->SetLabelSize(0.05); // Tweaks to X axis
  nRes_p_Hist->GetYaxis()->SetTitleOffset(1.0); nRes_p_Hist->GetYaxis()->SetTitleSize(0.05); nRes_p_Hist->GetYaxis()->SetLabelSize(0.05); // Tweaks to y-axis
  nRes_p_Hist->SetLineWidth(2);
  nRes_p_Hist->Draw("HISTERR");
  c_nResp->Print("ECCE_NIM_Plots/nRes_p_5on100.png"); c_nResp->Print("ECCE_NIM_Plots/nRes_p_5on100.pdf");

  // Write our created plots to root files too for easier editing   
  TFile* nXY_File = new TFile("ECCE_NIM_Plots/nXY_5on100.root" , "RECREATE");
  ZDCHist->Write();
  nXY_File->Close();
   
  TFile* tvstHist_File = new TFile("ECCE_NIM_Plots/tvstt_5on100.root" , "RECREATE");
  tvstHist->Write();
  tvstHist_File->Close();

  TFile* t_altvstHist_File = new TFile("ECCE_NIM_Plots/taltvstt_5on100.root" , "RECREATE");
  t_altvstHist->Write();
  t_altvstHist_File->Close();
   
  TFile* Q2tbinned_rates_v1_File = new TFile("ECCE_NIM_Plots/Q2t_rates_v1_5on100.root" , "RECREATE");
  c_Q2tbinned_rates_v1->Write(); // Here the canvas is saved, not the histo
  for(Int_t A = 0; A < 8; A++){ // Save the histos too
    tHists[A]->Write();
  }
  Q2tbinned_rates_v1_File->Close();
   
  TFile* nRes_p_Hist_File = new TFile("ECCE_NIM_Plots/nRes_p_Hist_5on100.root" , "RECREATE");
  nRes_p_Hist->Write();
  nRes_p_Hist_File->Close();

  if (((TDirectory*)InFile->Get("ZDC_Neutron_Info"))->GetListOfKeys()->Contains("n_ThetaPhiDiff")){
    TCanvas *c_nThetaPhiDiff = new TCanvas("c_nThetaPhiDiff", "Neutron ThetaDiff vs PhiDff", 100, 0, 2560, 1920);
    nThetaPhiDiff->Draw("COLZ"); nThetaPhiDiff->SetStats(0); nThetaPhiDiff->SetTitle("");
    c_nThetaPhiDiff->Update();
    TLine *LowerTheta = new TLine(-0.6,gPad->GetUymin(),-0.6,gPad->GetUymax());
    LowerTheta->SetLineColor(kRed); LowerTheta->SetLineWidth(3); LowerTheta->Draw();
    TLine *UpperTheta = new TLine(0.6,gPad->GetUymin(),0.6,gPad->GetUymax());
    UpperTheta->SetLineColor(kRed); UpperTheta->SetLineWidth(3); UpperTheta->Draw();
    TLine *LowerPhi = new TLine(gPad->GetUxmin(), -3, gPad->GetUxmax(), -3);
    LowerPhi->SetLineColor(kRed); LowerPhi->SetLineWidth(3); LowerPhi->Draw();
    TLine *UpperPhi = new TLine(gPad->GetUxmin(), 3, gPad->GetUxmax(), 3);
    UpperPhi->SetLineColor(kRed); UpperPhi->SetLineWidth(3); UpperPhi->Draw();
    c_nThetaPhiDiff->Print("ECCE_NIM_Plots/nThetaPhiDiff_5on100.png"); c_nThetaPhiDiff->Print("ECCE_NIM_Plots/nThetaPhiDiff_5on100.pdf");

    TFile* nThetaPhiDiff_File = new TFile("ECCE_NIM_Plots/nThetaPhiDiff_5on100.root", "RECREATE");
    c_nThetaPhiDiff->Write();
    nThetaPhiDiff->Write();
    nThetaPhiDiff_File->Close();
  }
   
  // Plots for the October EIC talk
  c_tvstt->Print("May22_Talk_Plots/tvstt_5on100.png");
  c_taltvstt->Print("May22_Talk_Plots/taltvstt_5on100.png");
  c_Q2tEff->Print("May22_Talk_Plots/Q2tEff_5on100.png");
  c_tvstt->Print("May22_Talk_Plots/tvstt_5on100.pdf");
  c_taltvstt->Print("May22_Talk_Plots/taltvstt_5on100.pdf");
  c_Q2tEff->Print("May22_Talk_Plots/Q2tEff_5on100.pdf");
  TCanvas *c_Q2Eff = new TCanvas("c_Q2Eff", "Detection Efficiency as a function of Q2", 100, 0, 2560, 1920);
  Q2EffHist->SetStats(0); Q2EffHist->Draw("HISTERR"); c_Q2Eff->Print("May22_Talk_Plots/Q2Eff_5on100.png"); c_Q2Eff->Print("May22_Talk_Plots/Q2Eff_5on100.pdf");

  TCanvas *c_Q2tbinned_rates_v3[3];
  c_Q2tbinned_rates_v3[0] = new TCanvas("c_Q2tbinned_rates_v3_1", "Pres_Results_p1", 100, 0, 5760, 1920);
  c_Q2tbinned_rates_v3[0]->Divide(3,1);
  c_Q2tbinned_rates_v3[0]->cd(1);
  tHists[0]->Draw("HISTERR");
  c_Q2tbinned_rates_v3[0]->cd(2);
  Q2tHists[0]->SetStats(0); Q2tHists[0]->Draw("COLZ");
  c_Q2tbinned_rates_v3[0]->cd(3);
  taltvst_Q2_Hists[0]->GetXaxis()->SetRangeUser(0, 0.4); taltvst_Q2_Hists[0]->SetStats(0);taltvst_Q2_Hists[0]->Draw("COLZ");
  c_Q2tbinned_rates_v3[0]->Print("May22_Talk_Plots/PresResults_5on100_Q2_6p25.png");
  c_Q2tbinned_rates_v3[0]->Print("May22_Talk_Plots/PresResults_5on100_Q2_6p25.pdf");
  
  c_Q2tbinned_rates_v3[1] = new TCanvas("c_Q2tbinned_rates_v3_2", "Pres_Results_p2", 100, 0, 5760, 1920);
  c_Q2tbinned_rates_v3[1]->Divide(3,1);
  c_Q2tbinned_rates_v3[1]->cd(1);
  tHists[3]->Draw("HISTERR");
  c_Q2tbinned_rates_v3[1]->cd(2);
  Q2tHists[3]->SetStats(0); Q2tHists[3]->Draw("COLZ");
  c_Q2tbinned_rates_v3[1]->cd(3);
  taltvst_Q2_Hists[3]->GetXaxis()->SetRangeUser(0, 0.4); taltvst_Q2_Hists[3]->SetStats(0); taltvst_Q2_Hists[3]->Draw("COLZ");
  c_Q2tbinned_rates_v3[1]->Print("May22_Talk_Plots/PresResults_5on100_Q2_17p5.png");
  c_Q2tbinned_rates_v3[1]->Print("May22_Talk_Plots/PresResults_5on100_Q2_17p5.pdf");
 
  c_Q2tbinned_rates_v3[2] = new TCanvas("c_Q2tbinned_rates_v3_3", "Pres_Results_p3", 100, 0, 5760, 1920);
  c_Q2tbinned_rates_v3[2]->Divide(3,1);
  c_Q2tbinned_rates_v3[2]->cd(1);
  tHists[6]->Draw("HISTERR");
  c_Q2tbinned_rates_v3[2]->cd(2);
  Q2tHists[6]->SetStats(0); Q2tHists[6]->Draw("COLZ");
  c_Q2tbinned_rates_v3[2]->cd(3);
  taltvst_Q2_Hists[6]->GetXaxis()->SetRangeUser(0, 0.4); taltvst_Q2_Hists[6]->SetStats(0); taltvst_Q2_Hists[6]->Draw("COLZ");
  c_Q2tbinned_rates_v3[2]->Print("May22_Talk_Plots/PresResults_5on100_Q2_32p5.png");
  c_Q2tbinned_rates_v3[2]->Print("May22_Talk_Plots/PresResults_5on100_Q2_32p5.pdf");

  // This section is for SIDIS comparison plots.
  // I've chosen the Q2 = 17.5 bin to draw, I scale both histograms to peak at 1 for comparison (actual values are irelevant).
  TFile *SIDISFile = new TFile("LundTest_AllRuns_NoSmear_Analysis_v2.root");
  TH1F* SIDIS_pMiss_Dist = (TH1F*)((TH1F*)SIDISFile->Get("SIDIS_neutron_(Pmiss)_Distributions/Neutron p Dist Q2 = 17.500000")); // Histogram is named very poorly in this file
  TH1F* pn_Dist[8];
  for(Int_t A = 1; A <9; A++){
    pn_Dist[A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results_Cuts/pn_Result_Q2_%i",A)));
  }
  
  TH1F* DEMP_pMiss_Dist = (TH1F*)((TH1F*)InFile->Get("Physics_Results_Cuts/pmiss_Result_Q2_4"));

  // Compairson of DEMP neutrons with DEMP "missing mass"
  //TCanvas *c_DEMPComp = new TCanvas("c_DEMPComp", "DEMP p_{n} vs p_{miss} Comparison", 100, 0, 2560, 1920);
  //pn_Dist[3]->SetStats(0); pn_Dist[3]->SetMinimum(0);
  //pn_Dist[3]->SetTitle("DEMP p_{n} vs p_{miss}, 15 < Q^{2} < 20"); pn_Dist[3]->GetXaxis()->SetTitle("DEMP p_{n}, p_{miss} (GeV/c)"); pn_Dist[3]->GetYaxis()->SetTitle("Rate(Hz)"); // Set titles
  pn_Dist[3]->SetLineWidth(2); pn_Dist[3]->Rebin(2);
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
  TCanvas *c_SIDISComp = new TCanvas("c_SIDISComp", "DEMP vs SIDIS Comparison", 100, 0, 2560, 1920);
  pn_Dist[3]->SetStats(0); pn_Dist[3]->SetTitle(""); pn_Dist[3]->GetXaxis()->SetTitle("DEMP p_{n}, SIDIS p_{miss} (GeV/c)"); pn_Dist[3]->GetYaxis()->SetTitle("Arbitrary Scale"); // Set titles
  pn_Dist[3]->SetMinimum(0.); SIDIS_pMiss_Dist->SetLineColor(2); SIDIS_pMiss_Dist->SetLineWidth(2);
  pn_Dist[3]->Draw("HIST"); SIDIS_pMiss_Dist->Draw("SAMEHIST"); // Set y axis to go from 0, set colour of SIDIS to red, draw both
  TLegend *SIDISCompLegend = new TLegend(0.20, 0.60, 0.50, 0.75); // Define legend
  SIDISCompLegend->AddEntry(pn_Dist[3], "DEMP Neutron Momentum", "l"); SIDISCompLegend->AddEntry(SIDIS_pMiss_Dist, "SIDIS Missing Momentum", "l"); // Add entries to legend
  SIDISCompLegend->SetBorderSize(0); SIDISCompLegend->SetFillStyle(0); SIDISCompLegend->SetTextSize(0.04);
  SIDISCompLegend->Draw("SAME");// Draw legend
  c_SIDISComp->Print("ECCE_NIM_Plots/SIDIS_Comp.png"); c_SIDISComp->Print("ECCE_NIM_Plots/SIDIS_Comp.pdf"); // Print to file

  // Write the scaled files to root files
  //TFile* DEMPComp_Scaled_File = new TFile("ECCE_NIM_Plots/DEMP_Comp_Scaled.root", "RECREATE");
  //c_DEMPComp_Scaled->Write(); DEMPCompLegend->Write(); pn_Dist[3]->Write(); DEMP_pMiss_Dist->Write();
  //DEMPComp_Scaled_File->Close();
  TFile* SIDISComp_File = new TFile("ECCE_NIM_Plots/SIDIS_Comp.root", "RECREATE");
  c_SIDISComp->Write(); SIDISCompLegend->Write(); pn_Dist[3]->Write(); SIDIS_pMiss_Dist->Write();
  SIDISComp_File->Close();
  
  SIDISFile->Close(); // Close SIDIS file
  InFile->Close();
}
