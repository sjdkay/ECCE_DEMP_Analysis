// 06/21/22 - Maggie Kerr, modified from prior file made by Stephen Kay

// last updated 06/21/2022

// A quick script to extract values from the bins in a histogram and save them to a .csv file
// Also plots a bunch of histograms too, outputs a pdf and some other plots (comment/uncomment as needed)
// Execute via root -l 'BinExtractor5on41.C("INPUT_FILE.root", "OUTPUT_NAME")'
#define BinExtractor_cxx

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

void BinExtractor5on41(string InFilename = "", string OutFilename = "") {

	gROOT->SetBatch(kTRUE); // Force script to always run without flashing up graphics

	TString rootFile;

	if (InFilename == "") {
		cout << "Enter a filename to analyse: ";
		cin >> InFilename;
	}

	if (OutFilename == "") {
		cout << "Enter a filename to output to: ";
		cin >> OutFilename;
	}

	TString TInFilename = InFilename;
	rootFile = TInFilename;

	if (gSystem->AccessPathName(rootFile) == kTRUE) {
		cerr << "!!! ERROR !!!" << endl << rootFile << " not found" << endl << "!!! ERROR !!!" << endl;
		exit(0);
	}

	TFile *InFile = new TFile(rootFile);
	TString TOutFilename = OutFilename;

	TLatex InfoDump;
	InfoDump.SetTextSize(0.1);
	InfoDump.SetTextAlign(12); // Align at centre

	TH1F* tHists[19];
  	TH1F* ttruthHists[19];
  	TH1F* Q2Hists[19];
  	TH1F* WHists[19];
  	TH2F* Q2tHists[19];
  	TH2F* tvst_Q2_Hists[19];
  	TH2F* taltvst_Q2_Hist[19];
  	TH1F* taltres_Hists[10];

	TH2F* pipThetaTruthHist = (TH2F*)((TH2F*)InFile->Get("Pion_Truth_Info/piTrack_pTheta_Truth"));
  	TH2F* epThetaTruthHist = (TH2F*)((TH2F*)InFile->Get("Scattered_Electron_Truth_Info/eTrack_pTheta_Truth"));
  	TH2F* npThetaTruthHist = (TH2F*)((TH2F*)InFile->Get("Neutron_Truth_Info/nTrack_pTheta_Truth"));

	TH2F* Q2WHist = (TH2F*)((TH2F*)InFile->Get("Physics_Results_Misc/Q2_W_Result"));
  	TH2F* piXYHist = (TH2F*)((TH2F*)InFile->Get("Pion_Info/pi_XY"));
  	TH2F* eXYHist = (TH2F*)((TH2F*)InFile->Get("Scattered_Electron_Info/e_XY"));
  	TH2F* ZDCHist;

	// Check which ZDC histogram exists, and grab the one that does
	if ( InFile->GetListOfKeys()->Contains("ZDC_XY_IP6")){
    		ZDCHist = (TH2F*)((TH2F*)InFile->Get("ZDC_XY_IP6"));
  	}

  	else {
    		ZDCHist = (TH2F*)((TH2F*)InFile->Get("ZDC_XY_IP8"));
  	}

	TH2F* nThetaPhiDiff;
  	// If the t binned t resolution plots exist, grab them
  	if (((TDirectory*)InFile->Get("ZDC_Neutron_Info"))->GetListOfKeys()->Contains("n_ThetaPhiDiff")) {
    		nThetaPhiDiff = (TH2F*)((TH2F*)InFile->Get("ZDC_Neutron_Info/n_ThetaPhiDiff"));
  	}
 	 else {
		nThetaPhiDiff = new TH2F("h2_Empty_Plot", "No Plot", 100, -10, 10, 100, -10, 10);
  	}

	// If the t binned t resolution plots exist, grab them
  	if (((TDirectory*)InFile->Get("Physics_Results_Misc"))->GetListOfKeys()->Contains("taltres_result_ttruth_1")) {
     		for (Int_t A = 0; A <10; A++) {
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

	for (Int_t A = 1; A < 19; A++) {
    		tHists[A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results/t_cut_Result_Q2_%i",A)));
    		ttruthHists[A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results/t_truth_thrown_Result_Q2_%i",A)));
    		Q2Hists[A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results/Q2_cut_Result_Q2_%i",A)));
    		WHists[A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results/W_cut_Result_Q2_%i",A)));
    		Q2tHists[A-1] = (TH2F*)((TH2F*)InFile->Get(Form("Physics_Results/Q2_t_cut_Result_Q2_%i",A)));
    		tvst_Q2_Hists[A-1] = (TH2F*)((TH2F*)InFile->Get(Form("Physics_Results/t_ttruth_Result_Q2_%i",A)));
    		taltvst_Q2_Hists[A-1] = (TH2F*)((TH2F*)InFile->Get(Form("Physics_Results/t_alt_ttruth_Result_Q2_%i",A)));
  	}

	Double_t BinVals[10];
  	Double_t BinErrs[10];
  	TString header = "Q2_Cent,Q2_mean,W_mean,0.02,0.06,0.10,0.14,0.18,0.22,0.26,0.30,0.34,0.38";
  	TString values[19];
  	TString errors[19];
	Double_t Q2BinVal[20] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0};

	for (Int_t A = 0; A < 19; A++) {
		values[A]=Form("%2.2f,", (Q2BinVal[A] + Q2BinVal[A+1])/2);
    		errors[A]=Form("%2.2f,", (Q2BinVal[A] + Q2BinVal[A+1])/2);
   		values[A]+=Form("%2.2f,", (Q2Hists[A]->GetMean()));
    		values[A]+=Form("%2.2f,", (WHists[A]->GetMean()));
    		errors[A]+=Form("%2.2f,", (Q2Hists[A]->GetMean()));
    		errors[A]+=Form("%2.2f,", (WHists[A]->GetMean()));

    		for (Int_t B = 0; B < 10; B++) {
      			BinVals[B]=tHists[A]->GetBinContent(B+1);
      			BinErrs[B]=tHists[A]->GetBinError(B+1);
      			if (B != 9) {
      				values[A]+=Form("%2.10f,",BinVals[B]);
      				errors[A]+=Form("%2.10f,",BinErrs[B]);
     			 }
      			else {
				values[A]+=Form("%2.10f",BinVals[B]);
				errors[A]+=Form("%2.10f",BinErrs[B]);
          		}
    		}
  	}

	ofstream Outfile;
	TString Outpdf = TOutFilename + ".pdf";
	Outfile.open(TOutFilename);
	Outfile << "Rates per -t bin for 19 Q2 settings. -t and Q2 settings are the centres of the bin or Q2 range in each case\n\n";
	Outfile << header << "\n";
	for (Int_t A = 0; A < 19; A++) {
		Outfile << errors[A] << "\n";
	}

	TCanvas *c_Output[35]; 
	// Pages 1 - 19
	for (Int_t A = 0; A < 19; A++) {
		c_Output[A] = new TCanvas(Form("c_Output_%i", (A+1)), Form("Results_p%i", (A+1)), 100, 0, 1000, 900);
		c_Output[A]->Divide(3,2);
		c_Output[A]->cd(1);
		tHists[A]->Draw("HISTERR");
    		c_Output[A]->cd(2);
    		Q2Hists[A]->Draw("HISTERR");
    		c_Output[A]->cd(3);
    		WHists[A]->Draw("HISTERR");
  	 	c_Output[A]->cd(4);
    		Q2tHists[A]->Draw("COLZ");
    		c_Output[A]->cd(5);
    		taltvst_Q2_Hists[A]->Draw("COLZ");
    		c_Output[A]->cd(6);
    		InfoDump.DrawLatex(.2,.8,"L = 10^{34} cm^{-2}s^{-1}");
    		InfoDump.DrawLatex(.2,.7,"assumed in rate");
    		InfoDump.DrawLatex(.2,.6,"calculation");

		if (A == 0) {
			c_Output[A]->Print(Outpdf + '(');
		}
		else {
			c_Output[A]->Print(Outpdf);
		}
	}

	// Page 20
	c_Output[19] = new TCanvas("c_Output_20", "Results_p20", 100, 0, 1000, 900);
  	c_Output[19]->Divide(3,2);
  	c_Output[19]->cd(1);
  	Q2WHist->Draw("COLZ");
  	c_Output[19]->cd(2);
  	ZDCHist->Draw("COLZ");
  	c_Output[19]->cd(3);
  	Q2EffHist->Draw("HISTERR");
  	c_Output[19]->cd(4);
  	Q2tEffHist->Draw("COLZ");
  	c_Output[19]->cd(5);
  	Q2tEffHist_v2->Draw("COLZ");
 	c_Output[19]->cd(6);
  	tresHist->Draw("HISTERR");
  	c_Output[19]->Print(Outpdf);

	// Page 21
	c_Output[20] = new TCanvas("c_Output_21", "Results_p21", 100, 0, 1000, 900);
  	c_Output[20]->Divide(1,2);
  	c_Output[20]->cd(1);
  	tvstHist->Draw("COLZ");
  	c_Output[20]->cd(2);
  	t_altvstHist->GetXaxis()->SetRangeUser(0,0.4); t_altvstHist->Draw("COLZ");
  	c_Output[20]->Print(Outpdf);

	// Page 22
	c_Output[21] = new TCanvas("c_Output_22", "Results_p22", 100, 0, 1000, 900);
   	c_Output[21]->Divide(2,2);
  	c_Output[21]->cd(1);
  	MmissHist->Draw("HISTERR");
  	c_Output[21]->cd(2);
  	MmissSqHist->Draw("HISTERR");
  	c_Output[21]->cd(3);
  	Mmiss_truth_Hist->Draw("HISTERR");
  	c_Output[21]->cd(4);
  	Mmiss_comp_Hist->Draw("HISTERR");
  	c_Output[21]->Print(Outpdf);

	// Page 23
	c_Output[22] = new TCanvas("c_Output_23", "Results_p23", 100, 0, 1000, 900);
        c_Output[22]->Divide(2,2);
  	c_Output[22]->cd(1);
  	piRes_p_Hist->Draw("HISTERR");
  	c_Output[22]->cd(2);
  	piRes_px_Hist->Draw("HISTERR");
  	c_Output[22]->cd(3);
  	piRes_py_Hist->Draw("HISTERR");
  	c_Output[22]->cd(4);
  	piRes_pz_Hist->Draw("HISTERR");
  	c_Output[22]->Print(Outpdf);

	// Page 24
	c_Output[23] = new TCanvas("c_Output_24", "Results_p24", 100, 0, 1000, 900);
  	c_Output[23]->Divide(2,2);
 	c_Output[23]->cd(1);
  	eRes_p_Hist->Draw("HISTERR");
  	c_Output[23]->cd(2);
  	eRes_px_Hist->Draw("HISTERR");
  	c_Output[23]->cd(3);
  	eRes_py_Hist->Draw("HISTERR");
  	c_Output[23]->cd(4);
  	eRes_pz_Hist->Draw("HISTERR");
  	c_Output[23]->Print(Outpdf);

	// Page 25
	c_Output[24] = new TCanvas("c_Output_25", "Results_p25", 100, 0, 1000, 900);
  	c_Output[24]->Divide(2,2);
  	c_Output[24]->cd(1);
  	nRes_p_Hist->Draw("HISTERR");
  	c_Output[24]->cd(2);
  	nRes_px_Hist->Draw("HISTERR");
  	c_Output[24]->cd(3);
  	nRes_py_Hist->Draw("HISTERR");
  	c_Output[24]->cd(4);
  	nRes_pz_Hist->Draw("HISTERR");
  	c_Output[24]->Print(Outpdf);

	// Page 26-30
	// Draw t-truth distributions
	c_Output[25] = new TCanvas("c_Output_26", "Results_p26", 100, 0, 1000, 900);
        c_Output[25]->Divide(2,2);
  	for (Int_t A = 1; A <= 4; A++) {
    		c_Output[25]->cd(A);
                ttruthHists[A-1]->Draw("HISTERR");
  	}
  	c_Output[25]->Print(Outpdf);

        c_Output[26] = new TCanvas("c_Output_27", "Results_p27", 100, 0, 1000, 900);
        c_Output[26]->Divide(2,2);
  	for (Int_t A = 5; A <= 8; A++) {
    		c_Output[26]->cd(A-4);
                ttruthHists[A-1]->Draw("HISTERR");
  	}
  	c_Output[26]->Print(Outpdf);
  	
  	c_Output[27] = new TCanvas("c_Output_28", "Results_p28", 100, 0, 1000, 900);
        c_Output[27]->Divide(2,2);
  	for (Int_t A = 9; A <= 12; A++) {
    		c_Output[27]->cd(A-8);
                ttruthHists[A-1]->Draw("HISTERR");
  	}
  	c_Output[27]->Print(Outpdf);
  	
  	c_Output[28] = new TCanvas("c_Output_29", "Results_p29", 100, 0, 1000, 900);
        c_Output[28]->Divide(2,2);
  	for (Int_t A = 13; A <= 16; A++) {
    		c_Output[28]->cd(A-12);
                ttruthHists[A-1]->Draw("HISTERR");
  	}
  	c_Output[28]->Print(Outpdf);
  	
  	c_Output[29] = new TCanvas("c_Output_30", "Results_p30", 100, 0, 1000, 900);
        c_Output[29]->Divide(2,2);
  	for (Int_t A = 17; A <= 19; A++) {
    		c_Output[29]->cd(A-16);
                ttruthHists[A-1]->Draw("HISTERR");
  	}
  	c_Output[29]->Print(Outpdf);
  	

	// Page 31-35
	// Draw reconstructed t distributions with cuts
	c_Output[30] = new TCanvas("c_Output_31", "Results_p31", 100, 0, 1000, 900);
  	c_Output[30]->Divide(2,2);
  	for (Int_t B = 1; B <= 4; B++) {
    		c_Output[30]->cd(B);
    		tHists[B-1]->SetStats(0);
    		tHists[B-1]->Draw("HISTERR");
 	}
  	c_Output[30]->Print(Outpdf);
 	
	c_Output[31] = new TCanvas("c_Output_32", "Results_p32", 100, 0, 1000, 900);
  	c_Output[31]->Divide(2,2);
 	for (Int_t A = 5; A <= 8; A++) {
    		c_Output[31]->cd(A-4);
    		tHists[A-1]->SetStats(0);
    		tHists[A-1]->Draw("HISTERR");
 	}
 	c_Output[31]->Print(Outpdf);
 	
 	
	c_Output[32] = new TCanvas("c_Output_33", "Results_p33", 100, 0, 1000, 900);
  	c_Output[32]->Divide(2,2);
 	for (Int_t A = 9; A <= 12; A++) {
    		c_Output[32]->cd(A-8);
    		tHists[A-1]->SetStats(0);
    		tHists[A-1]->Draw("HISTERR");
 	}
 	c_Output[32]->Print(Outpdf);
 	
 	c_Output[33] = new TCanvas("c_Output_34", "Results_p34", 100, 0, 1000, 900);
  	c_Output[33]->Divide(2,2);
 	for (Int_t A = 13; A <= 16; A++) {
    		c_Output[33]->cd(A-12);
    		tHists[A-1]->SetStats(0);
    		tHists[A-1]->Draw("HISTERR");
 	}
 	c_Output[33]->Print(Outpdf);
 	
 	
	c_Output[34] = new TCanvas("c_Output_35", "Results_p35", 100, 0, 1000, 900);
  	c_Output[34]->Divide(2,2);
 	for (Int_t A = 17; A <= 19; A++) {
    		c_Output[34]->cd(A+1);
    		tHists[A-1]->SetStats(0);
    		tHists[A-1]->Draw("HISTERR");
 	}
 	

	// Page 36
	// Check if t binned t resolution plots exist, print if they do exist. If they do not exist, p30 is the end of the file.
        if (((TDirectory*)InFile->Get("Physics_Results_Misc"))->GetListOfKeys()->Contains("taltres_result_ttruth_1")) {
        	c_Output[34]->Print(Outpdf);
    		c_Output[35] = new TCanvas("c_Output_36", "Results_p36", 100, 0, 1000, 900);
    		c_Output[35]->Divide(5,2);
	    	for(Int_t A = 0; A < 10; A++) {
	      		c_Output[35]->cd(A+1);
	      		taltres_Hists[A]->Draw("HISTERR");
	    	}
    		c_Output[35]->Print(Outpdf + ')');
  	}
  	else{
    		c_Output[34]->Print(Outpdf + ')');
 	}
  
  InFile->Close();
  Outfile.close();
}
