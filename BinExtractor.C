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

  TH1F* Hists[8];

  for(Int_t A = 1; A <9; A++){
    Hists[A-1] = (TH1F*)((TH1F*)InFile->Get(Form("Physics_Results/t_cut_Result_Q2_%i",A)));
  }

  Double_t BinVals[10];
  Double_t BinErrs[10];
  TString header = "Q2,0.02,0.06,0,10,0.14,0.18,0.22,0.26,0.30,0.34,0.38";
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
    for(Int_t B = 0; B < 11; B++){
      BinVals[B]=Hists[A]->GetBinContent(B+1);
      BinErrs[B]=Hists[A]->GetBinError(B+1);
      if (B != 10){
      values[A]+=Form("%1.10f,",BinVals[B]);
      errors[A]+=Form("%1.10f,",BinErrs[B]);
      }
      else{
	values[A]+=Form("%1.10f",BinVals[B]);
	errors[A]+=Form("%1.10f",BinErrs[B]);
      }
    }
  }

  InFile->Close();

  ofstream Outfile;
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

  Outfile.close();

}
