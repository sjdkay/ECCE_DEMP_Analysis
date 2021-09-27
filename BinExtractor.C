// 22/09/21 - Stephen JD Kay, University of Regina

// A quick script to extract values from the bins in a histogram and save them to a .csv file
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

void BinExtractor(string InFilename = "", string OutFilename = ""){
  
  TString rootFile;

  if(InFilename == ""){
    cout << "Enter a filename to analyse: ";
    cin >> InFilename;
  }
  if(OutFilename == ""){
    cout << "Enter a filename to output to: ";
    cin >> OutFilename;
  }
  
  TString TInFilename = InFilename;
  rootFile = TInFilename;

  if(gSystem->AccessPathName(rootFile) == kTRUE){
    cerr << "!!!!! ERROR !!!!!" << endl << rootFile << " not found" << endl << "!!!!! ERROR !!!!!" << endl;
    exit(0);
  }

  TFile *InFile = new TFile(rootFile);
  TString TOutFilename = OutFilename;

  TH1F* tHists[8];
  TH1F* Q2Hists[8];
  TH1F* WHists[8];
  TH2F* Q2tHists[8];

  TH2F* Q2WHist = (TH2F*)((TH2F*)InFile->Get("Physics_Results/Q2_W_Result"));
  TH2F* ZDCHist = (TH2F*)((TH2F*)InFile->Get("ZDC_XY"));
  TH1F* Q2EffHist = (TH1F*)((TH1F*)InFile->Get("Detection_Efficiency/Q2_DetEff"));
  TH2F* Q2tEffHist = (TH2F*)((TH2F*)InFile->Get("Detection_Efficiency/Q2_t_DetEff"));
  
  for(Int_t A = 1; A <9; A++){
    tHists[A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results/t_cut_Result_Q2_%i",A)));
    Q2Hists[A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results/Q2_cut_Result_Q2_%i",A)));
    WHists[A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results/W_cut_Result_Q2_%i",A)));
    Q2tHists[A-1] = (TH2F*)((TH2F*)InFile->Get(Form("Physics_Results/Q2_t_cut_Result_Q2_%i",A)));
  }

  Double_t BinVals[10];
  Double_t BinErrs[10];
  TString header = "Q2_Cent,Q2_mean,W_mean,0.02,0.06,0.10,0.14,0.18,0.22,0.26,0.30,0.34,0.38";
  TString values[8];
  TString errors[8];
  for(Int_t A = 0; A <8; A++){
    if( A <= 1){
      values[A]=Form("%2.2f,", (6.25+(A*2.5)));
      errors[A]=Form("%2.2f,", (6.25+(A*2.5)));
    }
    else{
      values[A]=Form("%2.2f,", (12.5+((A-2)*5)));
      errors[A]=Form("%2.2f,", (12.5+((A-2)*5)));
    }
    values[A]+=Form("%2.2f,", (Q2Hists[A]->GetMean()));
    values[A]+=Form("%2.2f,", (WHists[A]->GetMean()));
    errors[A]+=Form("%2.2f,", (Q2Hists[A]->GetMean()));
    errors[A]+=Form("%2.2f,", (WHists[A]->GetMean()));

    for(Int_t B = 0; B < 10; B++){
      BinVals[B]=tHists[A]->GetBinContent(B+1);
      BinErrs[B]=tHists[A]->GetBinError(B+1);
      if (B != 9){
      values[A]+=Form("%2.10f,",BinVals[B]);
      errors[A]+=Form("%2.10f,",BinErrs[B]);
      }
      else{
	values[A]+=Form("%2.10f",BinVals[B]);
	errors[A]+=Form("%2.10f",BinErrs[B]);
      }
    }
  }

  ofstream Outfile;
  TString Outpdf = TOutFilename+".pdf";
  Outfile.open(TOutFilename);
  Outfile << "Rates per -t bin for 8 Q2 settings. -t and Q2 settings are the centres of the bin or Q2 range in each case\n\n";
  Outfile << header << "\n";
  for(Int_t A = 0; A <8; A++){
    Outfile << values[A] << "\n";
  }
  Outfile << "\n";
  Outfile << "Errors per -t bin for 8 Q2 settings. -t and Q2 settings are the centres of the bin or Q2 range in each case\n\n";
  Outfile << header << "\n";
  for(Int_t A = 0; A <8; A++){
    Outfile << errors[A] << "\n";
  }
  
  TCanvas *c_Output[9];
  for(Int_t A = 0; A < 8; A++){
    c_Output[A] = new TCanvas(Form("c_Output_%i", (A+1)), Form("Results_p%i", (A+1)), 100, 0, 1000, 900);
    c_Output[A]->Divide(2,2);
    c_Output[A]->cd(1);
    tHists[A]->Draw("HISTERR");
    c_Output[A]->cd(2);
    Q2Hists[A]->Draw("HISTERR");
    c_Output[A]->cd(3);
    WHists[A]->Draw("HISTERR");
    c_Output[A]->cd(4);
    Q2tHists[A]->Draw("COLZ");
    if (A == 0){
      c_Output[A]->Print(Outpdf + '(');
    }
    else if(A < 7){
      c_Output[A]->Print(Outpdf);
    }
    else if(A == 7){
      c_Output[A]->Print(Outpdf);
    }
  }

  c_Output[8] = new TCanvas("c_Output_9", "Results_p9", 100, 0, 1000, 900);
  c_Output[8]->Divide(2,2);
  c_Output[8]->cd(1);
  Q2WHist->Draw("COLZ");
  c_Output[8]->cd(2);
  ZDCHist->Draw("COLZ");
  c_Output[8]->cd(3);
  Q2EffHist->Draw("HISTERR");
  c_Output[8]->cd(4);
  Q2tEffHist->Draw("COLZ");
  
  c_Output[8]->Print(Outpdf + ')');
  
  
  InFile->Close();  
  Outfile.close();

}
