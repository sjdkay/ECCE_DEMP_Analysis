// include statements
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

void Plots_5on41(string filename = "") {

	if(filename == "") {
		cout << "Enter a filename to analyze: ";
		cin >> filename;
	}

	TString Tfilename = filename;
	TString rootFilename = Tfilename;
	
	if (gSystem->AccessPathName(rootFilename) == kTRUE) {
		cerr << "error: " << rootFilename << " not found." << endl;
	}

	TFile *file = new TFile(rootFilename);

	// pMiss Histograms
	TH1F* pMiss_Dist[18];
	for (Int_t A = 1; A < 19; A++) {
		pMiss_Dist[A-1] = (TH1F*)((TH1F*)file->Get(Form("Physics_Results_Cuts/pn_Result_Q2_%i", A)));
	}

	Double_t pMiss_CutVals[18] = {38.5, 38.5, 38.5, 37.5, 37.5, 37.5, 37.0, 37.0, 37.0, 36.5, 36.5, 36.5, 35.5, 35.5, 35.5, 34.5, 34.5, 34.5};
	Double_t Q2Vals[2][18] = {{3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5}, {4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0}};
	TLine* pMiss_Lines[18];
	TCanvas* c1[3];
	c1[0] = new TCanvas("c1_Cut_Line", "P_miss Cuts for 5 on 41 with Q2 binning, p1", 100, 0, 1000, 900);
	c1[0]->Divide(2,3);
	
	for (Int_t A = 0; A < 6; A++) {
		c1[0]->cd(A+1);
		pMiss_Dist[A]->Rebin(2);
		pMiss_Dist[A]->SetStats(0); pMiss_Dist[A]->SetTitle(Form("p_{n} with p_{miss} cut, %0.2f < Q^2 < %0.2f", Q2Vals[0][A], Q2Vals[1][A])); pMiss_Dist[A]->GetXaxis()->SetTitle("p_{n} (GeV/c)"); pMiss_Dist[A]->GetYaxis()->SetTitle("Arbitrary Scale");
		pMiss_Dist[A]->GetXaxis()->SetRangeUser(10,45);
		pMiss_Dist[A]->SetLineWidth(2);
		Double_t yMax = pMiss_Dist[A]->GetYaxis()->GetXmax();
		pMiss_Lines[A] = new TLine(pMiss_CutVals[A], gPad->GetUymin(), pMiss_CutVals[A], pMiss_Dist[A]->GetMaximum()); pMiss_Lines[A]->SetLineColor(8); pMiss_Lines[A]->SetLineWidth(3);
		//pMiss_Lines[A] = new TLine(pMiss_CutVals[A], gPad->GetUymin(), pMiss_CutVals[A], gPad->GetUymax()); pMiss_Lines[A]->SetLineColor(8); pMiss_Lines[A]->SetLineWidth(3);
		pMiss_Dist[A]->Scale(1);
		pMiss_Dist[A]->Draw("HIST"); pMiss_Lines[A]->Draw("SAME");
	}
	
	// page 2
	
	c1[1] = new TCanvas("c1_Cut_Line", "P_miss Cuts for 5 on 41 with Q2 binning, p2", 100, 0, 1000, 900);
	c1[1]->Divide(2,3);
	
	for (Int_t A = 6; A < 12; A++) {
		c1[1]->cd(A-5);
		pMiss_Dist[A]->Rebin(2);
		pMiss_Dist[A]->SetStats(0); pMiss_Dist[A]->SetTitle(Form("p_{n} with p_{miss} cut, %0.2f < Q^2 < %0.2f", Q2Vals[0][A], Q2Vals[1][A])); pMiss_Dist[A]->GetXaxis()->SetTitle("p_{n} (GeV/c)"); pMiss_Dist[A]->GetYaxis()->SetTitle("Arbitrary Scale");
		pMiss_Dist[A]->GetXaxis()->SetRangeUser(10,45);
		pMiss_Dist[A]->SetLineWidth(2);
		Double_t yMax = pMiss_Dist[A]->GetYaxis()->GetXmax();
		pMiss_Lines[A] = new TLine(pMiss_CutVals[A], gPad->GetUymin(), pMiss_CutVals[A], pMiss_Dist[A]->GetMaximum()); pMiss_Lines[A]->SetLineColor(8); pMiss_Lines[A]->SetLineWidth(3);
		//pMiss_Lines[A] = new TLine(pMiss_CutVals[A], gPad->GetUymin(), pMiss_CutVals[A], gPad->GetUymax()); pMiss_Lines[A]->SetLineColor(8); pMiss_Lines[A]->SetLineWidth(3);
		pMiss_Dist[A]->Scale(1);
		pMiss_Dist[A]->Draw("HIST"); pMiss_Lines[A]->Draw("SAME");
	}
	
	// page 3
	
	c1[2] = new TCanvas("c1_Cut_Line", "P_miss Cuts for 5 on 41 with Q2 binning, p3", 100, 0, 1000, 900);
	c1[2]->Divide(2,3);
	
	for (Int_t A = 12; A < 18; A++) {
		c1[2]->cd(A-11);
		pMiss_Dist[A]->Rebin(2);
		pMiss_Dist[A]->SetStats(0); pMiss_Dist[A]->SetTitle(Form("p_{n} with p_{miss} cut, %0.2f < Q^2 < %0.2f", Q2Vals[0][A], Q2Vals[1][A])); pMiss_Dist[A]->GetXaxis()->SetTitle("p_{n} (GeV/c)"); pMiss_Dist[A]->GetYaxis()->SetTitle("Arbitrary Scale");
		pMiss_Dist[A]->GetXaxis()->SetRangeUser(10,45);
		pMiss_Dist[A]->SetLineWidth(2);
		Double_t yMax = pMiss_Dist[A]->GetYaxis()->GetXmax();
		pMiss_Lines[A] = new TLine(pMiss_CutVals[A], gPad->GetUymin(), pMiss_CutVals[A], pMiss_Dist[A]->GetMaximum()); pMiss_Lines[A]->SetLineColor(8); pMiss_Lines[A]->SetLineWidth(3);
		//pMiss_Lines[A] = new TLine(pMiss_CutVals[A], gPad->GetUymin(), pMiss_CutVals[A], gPad->GetUymax()); pMiss_Lines[A]->SetLineColor(8); pMiss_Lines[A]->SetLineWidth(3);
		pMiss_Dist[A]->Scale(1);
		pMiss_Dist[A]->Draw("HIST"); pMiss_Lines[A]->Draw("SAME");
	}
	
	// xy E
	TH3F* nTruth_xyE3D = (TH3F*)((TH3F*)file->Get(Form("nTruth_xyE3D")));
	TH3F* ZDC_xyE3D = (TH3F*)((TH3F*)file->Get(Form("ZDC_xyE3D")));
	TH3F* nTruth_Missed_xyE3D = (TH3F*)((TH3F*)file->Get(Form("nTruth_Missed_xyE3D")));
	TCanvas* c2 = new TCanvas("xyE_Dist", "xy vs E ZDC Truth Distributions for 5 on 41", 100, 0, 1000, 900);
	c2->Divide(2,2);
	c2->cd(1);
	nTruth_xyE3D->SetStats(0); nTruth_xyE3D->SetTitle("Truth xy position versus E for thrown events"); nTruth_xyE3D->GetXaxis()->SetTitle("x (cm)"); nTruth_xyE3D->GetYaxis()->SetTitle("y (cm)");
	nTruth_xyE3D->Project3DProfile("yx")->Draw("COLZ");
	c2->cd(2);
	ZDC_xyE3D->SetStats(0); ZDC_xyE3D->SetTitle("Truth xy position versus E for thrown events recorded in ZDC"); ZDC_xyE3D->GetXaxis()->SetTitle("x (cm)"); ZDC_xyE3D->GetYaxis()->SetTitle("y (cm)");
	ZDC_xyE3D->Project3DProfile("yx")->Draw("COLZ");
	c2->cd(3);
	nTruth_Missed_xyE3D->SetStats(0); nTruth_Missed_xyE3D->SetTitle("Truth xy position versus E for missed events"); nTruth_Missed_xyE3D->GetXaxis()->SetTitle("x (cm)"); nTruth_Missed_xyE3D->GetYaxis()->SetTitle("y (cm)");
	nTruth_Missed_xyE3D->Project3DProfile("yx")->Draw("COLZ");
	
	// XY ZDC Position
	
	// neutron Theta.
	
	//TH1F* nTruth_Theta = (TH1F*)((TH1F*)file->Get(Form("nTruth")))
	
	// missing Events
	
	//file->Close();	
}
