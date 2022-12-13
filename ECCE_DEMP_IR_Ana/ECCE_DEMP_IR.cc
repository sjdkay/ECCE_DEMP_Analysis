#include "ECCE_DEMP_IR.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllHistoManager.h>

#include <phool/PHCompositeNode.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNode.h>  // for PHNode
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>  // for PHObject
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <phool/phool.h>

// G4Hits includes
#include <g4main/PHG4Hit.h>
#include <g4main/PHG4HitContainer.h>

// Track includes
#include <trackbase_historic/SvtxTrackMap.h>

// Jet includes
#include <g4eval/JetEvalStack.h>
#include <g4jets/JetMap.h>
#include <g4eval/SvtxEvalStack.h>

// Cluster includes
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>

/// HEPMC truth includes
#include <HepMC/GenEvent.h>
#include <HepMC/GenVertex.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phhepmc/PHHepMCGenEventMap.h>

/// Fun4All includes
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/EicEventHeader.h>
#include <g4main/PHG4Reco.h>

#include <TFile.h>
#include <TNtuple.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cassert>
#include <sstream>
#include <string>
#include <iostream>
#include <stdexcept>
#include <stdio.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include <phparameter/PHParameters.h>

#include <pdbcalbase/PdbParameterMap.h>
#include <pdbcalbase/PdbParameterMapContainer.h>

using namespace std;

ECCE_DEMP_IR::ECCE_DEMP_IR(const std::string &name, const std::string& filename):
 SubsysReco(name)
 , outfilename(filename)
{
  std::cout << "ECCE_DEMP_IR_example::Diff_Tagg_example(const std::string &name) Calling ctor" << std::endl;

  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(m_RandomGenerator, seed);

}

//____________________________________________________________________________..
ECCE_DEMP_IR::~ECCE_DEMP_IR()
{

  gsl_rng_free(m_RandomGenerator);

  std::cout << "ECCE_DEMP_IR::~ECCE_DEMP_IR() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int ECCE_DEMP_IR::Init(PHCompositeNode *topNode)
{

  static_event_counter = 0;
  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(outfilename.c_str(), "RECREATE");

  std::cout << "ECCE_DEMP_IR::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  event_itt = 0;
  
  // Define beam 4 vectors - Assume IP6 by default, for other IPs, adjust in the Process Event loop (at the top)
  e_beam_energy = 5;
  e_beam_pmag = sqrt(pow(e_beam_energy,2)-pow(mElec,2));
  ion_beam_energy = 100;
  ion_beam_pmag = sqrt((pow(ion_beam_energy,2)-pow(mProt,2)));
  crossing_angle = 0.025; 
  eBeam4Vect.SetPxPyPzE(0,0,-1*e_beam_pmag,e_beam_energy);
  pBeam4Vect.SetPxPyPzE(-ion_beam_pmag*TMath::Sin(crossing_angle),0,ion_beam_pmag*TMath::Cos(crossing_angle),ion_beam_energy);

  // Set cut values for physics analysis
  Thetan_Cent = 1.45; // Cut will be +/- 0.4 from this value
  ThetaDiff_Cut = 0.6;
  PhiDiff_Cut = 3.0;
  // In future, the value below should be determined more accurately - average value of all 100 files?
  if (e_beam_energy == 5){
    if (ion_beam_energy == 100){
      nTried = 17000; // This is the approximate total number of events generated in a single file of this sample (5on100 epi)
    }
    else{
      nTried = 170000; // This is the approximate total number of events generated in a single file of this sample (5on41 epi)
    }
  }
  else if (e_beam_energy == 10){
    nTried = 2000; // This is the approximate total number of events generated in a single file of this sample (10on100 epi)
  }

  // Define some histograms
   gDirectory->mkdir("Electrons");
   gDirectory->cd("Electrons");
   // Full distributions
  h1_eTrack_px[0] = new TH1F("e_px_1", "e' p_{x} Distribution;p_{x} (GeV)", 240, -6, 6);
  h1_eTrack_py[0] = new TH1F("e_py_1", "e' p_{y} Distribution;p_{y} (GeV)", 240, -6, 6);
  h1_eTrack_pz[0] = new TH1F("e_pz_1", "e' p_{z} Distribution;p_{z} (GeV)", 120, -6, 0); 
  h1_eTrack_p[0] = new TH1F("e_p_1", "e' p Distribution;p (GeV)", 160, 0, 8);
  h1_eTrack_E[0] = new TH1F("e_E_1", "e' E Distribution;E (GeV)", 160, 0, 8);
  h1_eTrack_Theta[0] = new TH1F("e_Theta_1", "e' #theta Distribution; #theta (deg)", 200, 110, 160);
  h1_eTrack_Phi[0] = new TH1F("e_Phi_1", "e' #phi Distribution; #phi (deg)", 360, -180, 180);
  h2_eTrack_ThetaPhi[0] = new TH2F("e_ThetaPhi_1", "e' #theta vs #phi; #theta (deg); #phi (deg)", 140, 110, 180, 720, -180, 180);
  h2_eTrack_pTheta[0] = new TH2F("e_pTheta_1", "e' #theta vs P; #theta (deg); P (GeV/c)", 140, 110, 180, 100, 0, 10);
  
  //Q2 binned distributions
  gDirectory->mkdir("Q2_Binned");
  gDirectory->cd("Q2_Binned");
  for(Int_t A = 1; A < 7; A++){
    h1_eTrack_px[A] = new TH1F(Form("e_px_%i", (A+1)), Form("e' p_{x} Distribution, %i < Q^{2} < %i;p_{x} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 240, -6, 6);
    h1_eTrack_py[A] = new TH1F(Form("e_py_%i", (A+1)), Form("e' p_{y} Distribution, %i < Q^{2} < %i;p_{y} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))),  240, -6, 6);
    h1_eTrack_pz[A] = new TH1F(Form("e_pz_%i", (A+1)), Form("e' p_{z} Distribution, %i < Q^{2} < %i;p_{z} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 120, -6, 0); 
    h1_eTrack_p[A] = new TH1F(Form("e_p_%i", (A+1)), Form("e' p Distribution, %i < Q^{2} < %i;p (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 160, 0, 8);
    h1_eTrack_E[A] = new TH1F(Form("e_E_%i", (A+1)), Form("e' E Distribution, %i < Q^{2} < %i;E (GeV)", (5+((A-1)*5)), (10+((A-1)*5))),  160, 0, 8);
    h1_eTrack_Theta[A] = new TH1F(Form("e_Theta_%i", (A+1)), Form("e' #theta Distribution, %i < Q^{2} < %i; #theta (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 200, 110, 160);
    h1_eTrack_Phi[A] = new TH1F(Form("e_Phi_%i", (A+1)), Form("e' #phi Distribution, %i < Q^{2} < %i; #phi (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 360, -180, 180);
    h2_eTrack_ThetaPhi[A] = new TH2F(Form("e_ThetaPhi_%i", (A+1)), Form("e' #theta vs #phi, %i < Q^{2} < %i; #theta (deg); #phi (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 140, 110, 180, 720, -180, 180);
    h2_eTrack_pTheta[A] = new TH2F(Form("e_pTheta_%i", (A+1)), Form("e' #theta vs P, %i < Q^{2} < %i; #theta (deg); P (GeV/c)", (5+((A-1)*5)), (10+((A-1)*5))), 140, 110, 180, 100, 0, 10);
  }
  gDirectory->cd("../");
  gDirectory->cd("../");
    
  gDirectory->mkdir("Pions");
  gDirectory->cd("Pions");
  // Full distributions
  h1_piTrack_px[0] = new TH1F("pi_px_1", "#pi p_{x} Distribution;p_{x} (GeV)", 200, -20, 20);
  h1_piTrack_py[0] = new TH1F("pi_py_1", "#pi p_{y} Distribution;p_{y} (GeV)", 200, -20, 20);
  h1_piTrack_pz[0] = new TH1F("pi_pz_1", "#pi p_{z} Distribution;p_{z} (GeV)", 100, 0, 50); 
  h1_piTrack_p[0] = new TH1F("pi_p_1", "#pi p Distribution;p (GeV)", 200, 0, 50);
  h1_piTrack_E[0] = new TH1F("pi_E_1", "#pi E Distribution;E (GeV)", 200, 0, 50);
  h1_piTrack_Theta[0] = new TH1F("pi_Theta_1", "#pi #theta Distribution; #theta (deg)", 200, 0, 50);
  h1_piTrack_Phi[0] = new TH1F("pi_Phi_1", "#pi #phi Distribution; #phi (deg)", 360, -180, 180);
  h2_piTrack_ThetaPhi[0] = new TH2F("pi_ThetaPhi_1", "#pi #theta vs #phi; #theta (deg); #phi (deg)", 120, 0, 60, 720, -180, 180);
  h2_piTrack_pTheta[0] = new TH2F("pi_pTheta_1", "#pi #theta vs P; #theta (deg); P (GeV/c)", 120, 0, 60, 200, 0, 50);
  
  //Q2 binned distributions
  gDirectory->mkdir("Q2_Binned");
  gDirectory->cd("Q2_Binned");
  for(Int_t A = 1; A < 7; A++){
    h1_piTrack_px[A] = new TH1F(Form("pi_px_%i", (A+1)), Form("#pi p_{x} Distribution, %i < Q^{2} < %i;p_{x} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 200, -20, 20);
    h1_piTrack_py[A] = new TH1F(Form("pi_py_%i", (A+1)), Form("#pi p_{y} Distribution, %i < Q^{2} < %i;p_{y} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))),  200, -20, 20);
    h1_piTrack_pz[A] = new TH1F(Form("pi_pz_%i", (A+1)), Form("#pi p_{z} Distribution, %i < Q^{2} < %i;p_{z} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 100, 0, 50); 
    h1_piTrack_p[A] = new TH1F(Form("pi_p_%i", (A+1)), Form("#pi p Distribution, %i < Q^{2} < %i;p (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 200, 0, 50);
    h1_piTrack_E[A] = new TH1F(Form("pi_E_%i", (A+1)), Form("#pi E Distribution, %i < Q^{2} < %i;E (GeV)", (5+((A-1)*5)), (10+((A-1)*5))),  200, 0, 50);
    h1_piTrack_Theta[A] = new TH1F(Form("pi_Theta_%i", (A+1)), Form("#pi #theta Distribution, %i < Q^{2} < %i; #theta (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 200, 0, 50);
    h1_piTrack_Phi[A] = new TH1F(Form("pi_Phi_%i", (A+1)), Form("#pi #phi Distribution, %i < Q^{2} < %i; #phi (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 360, -180, 180);
    h2_piTrack_ThetaPhi[A] = new TH2F(Form("pi_ThetaPhi_%i", (A+1)), Form("#pi #theta vs #phi, %i < Q^{2} < %i; #theta (deg); #phi (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 120, 0, 60, 720, -180, 180);
    h2_piTrack_pTheta[A] = new TH2F(Form("pi_pTheta_%i", (A+1)), Form("#pi #theta vs P, %i < Q^{2} < %i; #theta (deg); P (GeV/c)", (5+((A-1)*5)), (10+((A-1)*5))), 120, 0, 60, 200, 0, 50);
  }
  gDirectory->cd("../");
  gDirectory->cd("../");

  gDirectory->mkdir("Neutrons");
  gDirectory->cd("Neutrons");
  // Full distributions
  h1_nTrack_px[0] = new TH1F("n_px_1", "n p_{x} Distribution;p_{x} (GeV)", 480, -6, 6);
  h1_nTrack_py[0] = new TH1F("n_py_1", "n p_{y} Distribution;p_{y} (GeV)", 480, -6, 6);
  h1_nTrack_pz[0] = new TH1F("n_pz_1", "n p_{z} Distribution;p_{z} (GeV)", 240, 0, 120); 
  h1_nTrack_p[0] = new TH1F("n_p_1", "n p Distribution;p (GeV)", 240, 0, 120);
  h1_nTrack_E[0] = new TH1F("n_E_1", "n E Distribution;E (GeV)", 240, 0, 120);
  h1_nTrack_Theta[0] = new TH1F("n_Theta_1", "n #theta Distribution; #theta (deg)", 100, 0, 5);
  h1_nTrack_Phi[0] = new TH1F("n_Phi_1", "n #phi Distribution; #phi (deg)", 360, -180, 180);
  h2_nTrack_ThetaPhi[0] = new TH2F("n_ThetaPhi_1", "n #theta vs #phi; #theta (deg); #phi (deg)", 120, 0, 60, 720, -180, 180);
  h2_nTrack_pTheta[0] = new TH2F("n_pTheta_1", "n #theta vs P; #theta (deg); P (GeV/c)", 100, 0, 5, 240, 0, 120);
  
  //Q2 binned distributions
  gDirectory->mkdir("Q2_Binned");
  gDirectory->cd("Q2_Binned");
  for(Int_t A = 1; A < 7; A++){
    h1_nTrack_px[A] = new TH1F(Form("n_px_%i", (A+1)), Form("n p_{x} Distribution, %i < Q^{2} < %i;p_{x} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 480, -6, 6);
    h1_nTrack_py[A] = new TH1F(Form("n_py_%i", (A+1)), Form("n p_{y} Distribution, %i < Q^{2} < %i;p_{y} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))),  480, -6, 6);
    h1_nTrack_pz[A] = new TH1F(Form("n_pz_%i", (A+1)), Form("n p_{z} Distribution, %i < Q^{2} < %i;p_{z} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 240, 0, 120); 
    h1_nTrack_p[A] = new TH1F(Form("n_p_%i", (A+1)), Form("n p Distribution, %i < Q^{2} < %i;p (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 240, 0, 120);
    h1_nTrack_E[A] = new TH1F(Form("n_E_%i", (A+1)), Form("n E Distribution, %i < Q^{2} < %i;E (GeV)", (5+((A-1)*5)), (10+((A-1)*5))),  240, 0, 120);
    h1_nTrack_Theta[A] = new TH1F(Form("n_Theta_%i", (A+1)), Form("n #theta Distribution, %i < Q^{2} < %i; #theta (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 100, 0, 5);
    h1_nTrack_Phi[A] = new TH1F(Form("n_Phi_%i", (A+1)), Form("n #phi Distribution, %i < Q^{2} < %i; #phi (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 360, -180, 180);
    h2_nTrack_ThetaPhi[A] = new TH2F(Form("n_ThetaPhi_%i", (A+1)), Form("n #theta vs #phi, %i < Q^{2} < %i; #theta (deg); #phi (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 500, 0, 5, 720, -180, 180);
    h2_nTrack_pTheta[A] = new TH2F(Form("n_pTheta_%i", (A+1)), Form("n #theta vs P, %i < Q^{2} < %i; #theta (deg); P (GeV/c)", (5+((A-1)*5)), (10+((A-1)*5))), 100, 0, 5, 240, 0, 120);
  }
  gDirectory->cd("../");
  gDirectory->cd("../");
    
  gDirectory->mkdir("en_Coincidences");
  gDirectory->cd("en_Coincidences");
  h1_e_enCoin_px[0] = new TH1F("e_en_px_1", "e' (en Coin) p_{x} Distribution;p_{x} (GeV)", 240, -6, 6);
  h1_e_enCoin_py[0] = new TH1F("e_en_py_1", "e' (en Coin) p_{y} Distribution;p_{y} (GeV)", 240, -6, 6);
  h1_e_enCoin_pz[0] = new TH1F("e_en_pz_1", "e' (en Coin) p_{z} Distribution;p_{z} (GeV)", 120, -6, 0); 
  h1_e_enCoin_p[0] = new TH1F("e_en_p_1", "e' (en Coin) p Distribution;p (GeV)", 160, 0, 8);
  h1_e_enCoin_E[0] = new TH1F("e_en_E_1", "e' (en Coin) E Distribution;E (GeV)", 160, 0, 8);
  h1_e_enCoin_Theta[0] = new TH1F("e_en_Theta_1", "e' (en Coin) #theta Distribution; #theta (deg)", 200, 110, 160);
  h1_e_enCoin_Phi[0] = new TH1F("e_en_Phi_1", "e' (en Coin) #phi Distribution; #phi (deg)", 360, -180, 180);
  h2_e_enCoin_ThetaPhi[0] = new TH2F("e_en_ThetaPhi_1", "e' (en Coin) #theta vs #phi; #theta (deg); #phi (deg)", 140, 110, 180, 720, -180, 180);
  h2_e_enCoin_pTheta[0] = new TH2F("e_en_pTheta_1", "e' (en Coin) #theta vs P; #theta (deg); P (GeV/c)", 140, 110, 180, 100, 0, 10);
  h1_n_enCoin_px[0] = new TH1F("n_en_px_1", "n (en Coin) p_{x} Distribution;p_{x} (GeV)", 480, -6, 6);
  h1_n_enCoin_py[0] = new TH1F("n_en_py_1", "n (en Coin) p_{y} Distribution;p_{y} (GeV)", 480, -6, 6);
  h1_n_enCoin_pz[0] = new TH1F("n_en_pz_1", "n (en Coin) p_{z} Distribution;p_{z} (GeV)", 240, 0, 120); 
  h1_n_enCoin_p[0] = new TH1F("n_en_p_1", "n (en Coin) p Distribution;p (GeV)", 240, 0, 120);
  h1_n_enCoin_E[0] = new TH1F("n_en_E_1", "n (en Coin) E Distribution;E (GeV)", 240, 0, 120);
  h1_n_enCoin_Theta[0] = new TH1F("n_en_Theta_1", "n (en Coin) #theta Distribution; #theta (deg)", 100, 0, 5);
  h1_n_enCoin_Phi[0] = new TH1F("n_en_Phi_1", "n (en Coin) #phi Distribution; #phi (deg)", 360, -180, 180);
  h2_n_enCoin_ThetaPhi[0] = new TH2F("n_en_ThetaPhi_1", "n (en Coin) #theta vs #phi; #theta (deg); #phi (deg)", 120, 0, 60, 720, -180, 180);
  h2_n_enCoin_pTheta[0] = new TH2F("n_en_pTheta_1", "n (en Coin) #theta vs P; #theta (deg); P (GeV/c)", 100, 0, 5, 240, 0, 120);

  // Q2 binned distributions
  gDirectory->mkdir("Q2_Binned");
  gDirectory->cd("Q2_Binned");
  for(Int_t A = 1; A < 7; A++){
    h1_e_enCoin_px[A] = new TH1F(Form("e_en_px_%i", (A+1)), Form("e' (en Coin) p_{x} Distribution, %i < Q^{2} < %i;p_{x} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 240, -6, 6);
    h1_e_enCoin_py[A] = new TH1F(Form("e_en_py_%i", (A+1)), Form("e' (en Coin) p_{y} Distribution, %i < Q^{2} < %i;p_{y} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))),  240, -6, 6);
    h1_e_enCoin_pz[A] = new TH1F(Form("e_en_pz_%i", (A+1)), Form("e' (en Coin) p_{z} Distribution, %i < Q^{2} < %i;p_{z} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 120, -6, 0); 
    h1_e_enCoin_p[A] = new TH1F(Form("e_en_p_%i", (A+1)), Form("e' (en Coin) p Distribution, %i < Q^{2} < %i;p (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 160, 0, 8);
    h1_e_enCoin_E[A] = new TH1F(Form("e_en_E_%i", (A+1)), Form("e' (en Coin) E Distribution, %i < Q^{2} < %i;E (GeV)", (5+((A-1)*5)), (10+((A-1)*5))),  160, 0, 8);
    h1_e_enCoin_Theta[A] = new TH1F(Form("e_en_Theta_%i", (A+1)), Form("e' (en Coin) #theta Distribution, %i < Q^{2} < %i; #theta (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 200, 110, 160);
    h1_e_enCoin_Phi[A] = new TH1F(Form("e_en_Phi_%i", (A+1)), Form("e' (en Coin) #phi Distribution, %i < Q^{2} < %i; #phi (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 360, -180, 180);
    h2_e_enCoin_ThetaPhi[A] = new TH2F(Form("e_en_ThetaPhi_%i", (A+1)), Form("e'(en Coin) #theta vs #phi, %i < Q^{2} < %i; #theta (deg); #phi (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 140, 110, 180, 720, -180, 180);
    h2_e_enCoin_pTheta[A] = new TH2F(Form("e_en_pTheta_%i", (A+1)), Form("e' (en Coin) #theta vs P, %i < Q^{2} < %i; #theta (deg); P (GeV/c)", (5+((A-1)*5)), (10+((A-1)*5))), 140, 110, 180, 100, 0, 10);
    h1_n_enCoin_px[A] = new TH1F(Form("n_en_px_%i", (A+1)), Form("n (en Coin) p_{x} Distribution, %i < Q^{2} < %i;p_{x} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 480, -6, 6);
    h1_n_enCoin_py[A] = new TH1F(Form("n_en_py_%i", (A+1)), Form("n (en Coin) p_{y} Distribution, %i < Q^{2} < %i;p_{y} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))),  480, -6, 6);
    h1_n_enCoin_pz[A] = new TH1F(Form("n_en_pz_%i", (A+1)), Form("n (en Coin) p_{z} Distribution, %i < Q^{2} < %i;p_{z} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 240, 0, 120); 
    h1_n_enCoin_p[A] = new TH1F(Form("n_en_p_%i", (A+1)), Form("n (en Coin) p Distribution, %i < Q^{2} < %i;p (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 240, 0, 120);
    h1_n_enCoin_E[A] = new TH1F(Form("n_en_E_%i", (A+1)), Form("n (en Coin) E Distribution, %i < Q^{2} < %i;E (GeV)", (5+((A-1)*5)), (10+((A-1)*5))),  240, 0, 120);
    h1_n_enCoin_Theta[A] = new TH1F(Form("n_en_Theta_%i", (A+1)), Form("n (en Coin) #theta Distribution, %i < Q^{2} < %i; #theta (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 100, 0, 5);
    h1_n_enCoin_Phi[A] = new TH1F(Form("n_en_Phi_%i", (A+1)), Form("n (en Coin) #phi Distribution, %i < Q^{2} < %i; #phi (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 360, -180, 180);
    h2_n_enCoin_ThetaPhi[A] = new TH2F(Form("n_en_ThetaPhi_%i", (A+1)), Form("n (en Coin) #theta vs #phi, %i < Q^{2} < %i; #theta (deg); #phi (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 500, 0, 5, 720, -180, 180);
    h2_n_enCoin_pTheta[A] = new TH2F(Form("n_en_pTheta_%i", (A+1)), Form("n (en Coin) #theta vs P, %i < Q^{2} < %i; #theta (deg); P (GeV/c)", (5+((A-1)*5)), (10+((A-1)*5))), 100, 0, 5, 240, 0, 120);
  } 
  gDirectory->cd("../");
  gDirectory->cd("../");

  gDirectory->mkdir("en_Coincidences_Pions");
  gDirectory->cd("en_Coincidences_Pions");

  gDirectory->mkdir("1D");
  gDirectory->cd("1D");

  h1_piTrue_enCoin_px[0] = new TH1F("piTrue_en_px_1", "#pi_{True} (en Coin) p_{x} Distribution;p_{x} (GeV)", 200, -20, 20);
  h1_piTrue_enCoin_py[0] = new TH1F("piTrue_en_py_1", "#pi_{True} (en Coin) p_{y} Distribution;p_{y} (GeV)", 200, -20, 20);
  h1_piTrue_enCoin_pz[0] = new TH1F("piTrue_en_pz_1", "#pi_{True} (en Coin) p_{z} Distribution;p_{z} (GeV)", 100, 0, 50); 
  h1_piTrue_enCoin_p[0] = new TH1F("piTrue_en_p_1", "#pi_{True} (en Coin) p Distribution;p (GeV)", 200, 0, 50);
  h1_piTrue_enCoin_E[0] = new TH1F("piTrue_en_E_1", "#pi_{True} (en Coin) E Distribution;E (GeV)", 200, 0, 50);
  h1_piTrue_enCoin_Theta[0] = new TH1F("piTrue_en_Theta_1", "#pi_{True} (en Coin) #theta Distribution; #theta (deg)", 200, 0, 50);
  h1_piTrue_enCoin_Phi[0] = new TH1F("piTrue_en_Phi_1", "#pi_{True} (en Coin) #phi Distribution; #phi (deg)", 360, -180, 180);
  h1_piTrue_Missed_enCoin_px[0] = new TH1F("piTrue_Missed_en_px_1", "#pi_{True} (en Coin, #pi missed) p_{x} Distribution;p_{x} (GeV)", 200, -20, 20);
  h1_piTrue_Missed_enCoin_py[0] = new TH1F("piTrue_Missed_en_py_1", "#pi_{True} (en Coin, #pi missed) p_{y} Distribution;p_{y} (GeV)", 200, -20, 20);
  h1_piTrue_Missed_enCoin_pz[0] = new TH1F("piTrue_Missed_en_pz_1", "#pi_{True} (en Coin, #pi missed) p_{z} Distribution;p_{z} (GeV)", 100, 0, 50); 
  h1_piTrue_Missed_enCoin_p[0] = new TH1F("piTrue_Missed_en_p_1", "#pi_{True} (en Coin, #pi missed) p Distribution;p (GeV)", 200, 0, 50);
  h1_piTrue_Missed_enCoin_E[0] = new TH1F("piTrue_Missed_en_E_1", "#pi_{True} (en Coin, #pi missed) E Distribution;E (GeV)", 200, 0, 50);
  h1_piTrue_Missed_enCoin_Theta[0] = new TH1F("piTrue_Missed_en_Theta_1", "#pi_{True} (en Coin, #pi missed) #theta Distribution; #theta (deg)", 200, 0, 50);
  h1_piTrue_Missed_enCoin_Phi[0] = new TH1F("piTrue_Missed_en_Phi_1", "#pi_{True} (en Coin, #pi missed) #phi Distribution; #phi (deg)", 360, -180, 180);

  gDirectory->cd("../");
  gDirectory->mkdir("2D");
  gDirectory->cd("2D");

  h2_piTrue_enCoin_ThetaPhi[0] = new TH2F("PiTrue_en_ThetaPhi_1", "#pi_{True} (en Coin) #theta vs #phi; #theta (deg); #phi (deg)", 120, 0, 60, 720, -180, 180);
  h2_piTrue_enCoin_pTheta[0] = new TH2F("piTrue_en_pTheta_1", "#pi_{True} (en Coin) #theta vs P; #theta (deg); P (GeV/c)", 120, 0, 60, 200, 0, 50);
  h2_piTrue_Missed_enCoin_ThetaPhi[0] = new TH2F("PiTrue_Missed_en_ThetaPhi_1", "#pi_{True} (en Coin, #pi missed) #theta vs #phi; #theta (deg); #phi (deg)", 120, 0, 60, 720, -180, 180);
  h2_piTrue_Missed_enCoin_pTheta[0] = new TH2F("piTrue_Missed_en_pTheta_1", "#pi_{True} (en Coin, #pi missed) #theta vs P; #theta (deg); P (GeV/c)", 120, 0, 60, 200, 0, 50);

  gDirectory->cd("../");
  gDirectory->mkdir("Q2_Binned");
  gDirectory->cd("Q2_Binned");
  
  for(Int_t A = 1; A < 7; A++){
    h1_piTrue_enCoin_px[A] = new TH1F(Form("piTrue_en_px_%i", (A+1)), Form("#pi_{True} (en Coin) p_{x} Distribution, %i < Q^{2} < %i;p_{x} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 200, -20, 20);
    h1_piTrue_enCoin_py[A] = new TH1F(Form("piTrue_en_py_%i", (A+1)), Form("#pi_{True} (en Coin) p_{y} Distribution, %i < Q^{2} < %i;p_{y} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))),  200, -20, 20);
    h1_piTrue_enCoin_pz[A] = new TH1F(Form("piTrue_en_pz_%i", (A+1)), Form("#pi_{True} (en Coin) p_{z} Distribution, %i < Q^{2} < %i;p_{z} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 100, 0, 50); 
    h1_piTrue_enCoin_p[A] = new TH1F(Form("piTrue_en_p_%i", (A+1)), Form("#pi_{True} (en Coin) p Distribution, %i < Q^{2} < %i;p (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 200, 0, 50);
    h1_piTrue_enCoin_E[A] = new TH1F(Form("piTrue_en_E_%i", (A+1)), Form("#pi_{True} (en Coin) E Distribution, %i < Q^{2} < %i;E (GeV)", (5+((A-1)*5)), (10+((A-1)*5))),  200, 0, 50);
    h1_piTrue_enCoin_Theta[A] = new TH1F(Form("piTrue_en_Theta_%i", (A+1)), Form("#pi_{True} (en Coin) #theta Distribution, %i < Q^{2} < %i; #theta (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 200, 0, 50);
    h1_piTrue_enCoin_Phi[A] = new TH1F(Form("piTrue_en_Phi_%i", (A+1)), Form("#pi_{True} (en Coin) #phi Distribution, %i < Q^{2} < %i; #phi (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 360, -180, 180);
    h2_piTrue_enCoin_ThetaPhi[A] = new TH2F(Form("piTrue_en_ThetaPhi_%i", (A+1)), Form("#pi_{True} (en Coin) #theta vs #phi, %i < Q^{2} < %i; #theta (deg); #phi (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 120, 0, 60, 720, -180, 180);
    h2_piTrue_enCoin_pTheta[A] = new TH2F(Form("piTrue_en_pTheta_%i", (A+1)), Form("#pi_{True} (en Coin) #theta vs P, %i < Q^{2} < %i; #theta (deg); P (GeV/c)", (5+((A-1)*5)), (10+((A-1)*5))), 120, 0, 60, 200, 0, 50);
    h1_piTrue_Missed_enCoin_px[A] = new TH1F(Form("piTrue_Missed_en_px_%i", (A+1)), Form("#pi_{True} (en Coin, #pi missed) p_{x} Distribution, %i < Q^{2} < %i;p_{x} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 200, -20, 20);
    h1_piTrue_Missed_enCoin_py[A] = new TH1F(Form("piTrue_Missed_en_py_%i", (A+1)), Form("#pi_{True} (en Coin, #pi missed) p_{y} Distribution, %i < Q^{2} < %i;p_{y} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))),  200, -20, 20);
    h1_piTrue_Missed_enCoin_pz[A] = new TH1F(Form("piTrue_Missed_en_pz_%i", (A+1)), Form("#pi_{True} (en Coin, #pi missed) p_{z} Distribution, %i < Q^{2} < %i;p_{z} (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 100, 0, 50); 
    h1_piTrue_Missed_enCoin_p[A] = new TH1F(Form("piTrue_Missed_en_p_%i", (A+1)), Form("#pi_{True} (en Coin, #pi missed) p Distribution, %i < Q^{2} < %i;p (GeV)", (5+((A-1)*5)), (10+((A-1)*5))), 200, 0, 50);
    h1_piTrue_Missed_enCoin_E[A] = new TH1F(Form("piTrue_Missed_en_E_%i", (A+1)), Form("#pi_{True} (en Coin, #pi missed) E Distribution, %i < Q^{2} < %i;E (GeV)", (5+((A-1)*5)), (10+((A-1)*5))),  200, 0, 50);
    h1_piTrue_Missed_enCoin_Theta[A] = new TH1F(Form("piTrue_Missed_en_Theta_%i", (A+1)), Form("#pi_{True} (en Coin, #pi missed) #theta Distribution, %i < Q^{2} < %i; #theta (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 200, 0, 50);
    h1_piTrue_Missed_enCoin_Phi[A] = new TH1F(Form("piTrue_Missed_en_Phi_%i", (A+1)), Form("#pi_{True} (en Coin, #pi missed) #phi Distribution, %i < Q^{2} < %i; #phi (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 360, -180, 180);
    h2_piTrue_Missed_enCoin_ThetaPhi[A] = new TH2F(Form("piTrue_Missed_en_ThetaPhi_%i", (A+1)), Form("#pi_{True} (en Coin, #pi missed) #theta vs #phi, %i < Q^{2} < %i; #theta (deg); #phi (deg)", (5+((A-1)*5)), (10+((A-1)*5))), 120, 0, 60, 720, -180, 180);
    h2_piTrue_Missed_enCoin_pTheta[A] = new TH2F(Form("piTrue_Missed_en_pTheta_%i", (A+1)), Form("#pi_{True} (en Coin, #pi missed) #theta vs P, %i < Q^{2} < %i; #theta (deg); P (GeV/c)", (5+((A-1)*5)), (10+((A-1)*5))), 120, 0, 60, 200, 0, 50);  
  }
  gDirectory->cd("../");
  gDirectory->cd("../");
  
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP_IR::InitRun(PHCompositeNode *topNode)
{
  if( static_event_counter == 0){
    encloseure_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_hFarFwdBeamLineEnclosure_0");
    encloseure_nodeparams->Print();
    if (encloseure_nodeparams){
      Enclosure_params.FillFrom(encloseure_nodeparams, 0);
    } 
    else {
      cerr << "There is a issue finding the detector paramter node!" << endl;
    }

    beamlinemagnet_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_BEAMLINEMAGNET");

    if (encloseure_nodeparams){
      BeamLineMagnet_params.FillFrom(beamlinemagnet_nodeparams, 0);
    }
    else{
      cerr << "There is a issue finding the detector paramter node!" << endl;
    }

    zdc_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_ZDCsurrogate");
    // SJDK - 15/02/22 - This analysis doesn't really use any of these detectors, so beyond getting them to determine the IP, we don't reall care
    rp_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_rpTruth");
    rp2_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_rpTruth2");
    b0_nodeparams = findNode::getClass<PdbParameterMapContainer>(topNode, "G4GEOPARAM_b0Truth");

    static_event_counter++;
	
    /// Determining which IP design
    if (zdc_nodeparams) {
      if (rp2_nodeparams) {
	IP_design = "IP8";
      } 
      else {
	IP_design = "IP6";
      }
    } 
    else {
      IP_design = "UNKNOWN";
    }
  }

  cout << " END initialization" << endl;

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP_IR::process_event(PHCompositeNode *topNode)
{
  ZDC_hit = 0;
  EEMC_hit = 0;
  event_itt++; 
  
  if (IP_design == "IP8" ){ // If IP8, need to adjust crossing angle and adjust incoming beams
    crossing_angle = 0.035; // IP8 has a 35 mRad crossing angle
    Thetan_Cent = 2.005; // Adjust the theta cut for IP8 - 2.005 is approx 35 mRad
    eBeam4Vect.SetPxPyPzE(0,0,-1*e_beam_pmag,e_beam_energy);
    pBeam4Vect.SetPxPyPzE(ion_beam_pmag*TMath::Sin(crossing_angle),0,ion_beam_pmag*TMath::Cos(crossing_angle),ion_beam_energy);
  }
  
  if(event_itt%100 == 0)
    std::cout << "Event Processing Counter: " << event_itt << endl;
  // Get event header info for the event (weight)
  EicEventHeader* Evtheader = findNode::getClass<EicEventHeader>(topNode, "EicEventHeader");
  wgt = Evtheader->get_demp_weight();
  // Get MC truth info
  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  // Get the primary particle range
  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
  if (!truthinfo)
    {
      cout << PHWHERE
	   << "PHG4TruthInfoContainer node is missing, can't collect G4 truth particles"
	   << endl;
      return Fun4AllReturnCodes::EVENT_OK;
    }
  // Loop over the G4 truth (stable) particles
  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
       iter != range.second;
       ++iter)
    {
      // Get this truth particle
      const PHG4Particle *truth = iter->second;
      if ( truth->get_pid() == 11){ // PDG 11 -> Scattered electron
	e4VectTruth.SetPxPyPzE(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
      }
      else if (truth->get_pid() == 211){ // PDG 211 -> Pion 
	pi4VectTruth.SetPxPyPzE(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
      }
      else if (truth->get_pid() == 2112){ // PDG 2112 -> Neutron
	n4VectTruth.SetPxPyPzE(truth->get_px(), truth->get_py(), truth->get_pz(), truth->get_e());
      }
    }

  // Truth versions of kinematic quantities
  virtphoton4VectTruth = eBeam4Vect - e4VectTruth;
  //t4VectTruth = virtphoton4VectTruth - pi4VectTruth;
  t4VectTruth = eBeam4Vect - (e4VectTruth + pi4VectTruth);
  t_alt4VectTruth= pBeam4Vect - n4VectTruth;
  pmiss4VectTruth = (eBeam4Vect + pBeam4Vect) - (e4VectTruth+pi4VectTruth);
  pmiss4VectTruth_2 = (eBeam4Vect + pBeam4Vect) - (e4VectTruth+pi4VectTruth+n4VectTruth);
  Q2_truth = -1*(virtphoton4VectTruth.Mag2());
  W_truth = (virtphoton4VectTruth+pBeam4Vect).Mag();
  t_truth = -(t4VectTruth.Mag2());
  t_alt_truth = -(t_alt4VectTruth.Mag2());
  xb_truth =  Q2_truth/(2*(pBeam4Vect.Dot(virtphoton4VectTruth)));
  xi_truth = xb_truth/(2-xb_truth);
  y_inv_truth = (pBeam4Vect.Dot(virtphoton4VectTruth))/(pBeam4Vect.Dot(eBeam4Vect));// Calculation of the fractional energy loss

  // New loop to try and get some B0 hits, B0 tracker and calorimeter
  // B0 info ONLY exists as truth info, apply smearing functions
  PHG4HitContainer* B0_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_b0Truth");
  if(B0_hits){
    PHG4HitContainer::ConstRange B0_hit_range = B0_hits->getHits();
    for (PHG4HitContainer::ConstIterator B0_hit_iter = B0_hit_range.first; B0_hit_iter != B0_hit_range.second; B0_hit_iter++){
      B0_ETrue = B0_hit_iter->second->get_edep();
      B0_PosTrue.SetXYZ(B0_hit_iter->second->get_x(0), B0_hit_iter->second->get_y(0), B0_hit_iter->second->get_z(0));

      B0_ESmear = B0Cal_Energy_Smear(B0_hit_iter->second->get_edep());
      B0_PosSmear.SetXYZ(B0Cal_Position_Smear(B0_hit_iter->second->get_x(0)), B0Cal_Position_Smear(B0_hit_iter->second->get_y(0)), B0Cal_Position_Smear(B0_hit_iter->second->get_z(0)));
    }
  }

  if (Check_e(topNode) == true){
    SvtxTrackMap* e_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

    if (!e_trackmap)
      {
	e_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
	if (!e_trackmap)
	  {
	    cout
	      << "ECCE_DEMP::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
	    exit(-1);
	  }
      }

    for (SvtxTrackMap::Iter e_iter = e_trackmap->begin();
	 e_iter != e_trackmap->end();
	 ++e_iter)
      {
	SvtxTrack* e_track = e_iter->second;

	if (e_track->get_pz() < 0  && e_track->get_charge() == -1 ){ // -ve z direction -> electrons, crappy way of selecting them for now w/o truth info
	  eVect.SetXYZ(e_track->get_px(), e_track->get_py(), e_track->get_pz());
	  e4Vect.SetPxPyPzE(e_track->get_px(), e_track->get_py(), e_track->get_pz(), sqrt(pow(eVect.Mag(), 2)+pow(mElec,2)));
	}

	h1_eTrack_px[0]->Fill(e4Vect.Px(), wgt);
	h1_eTrack_py[0]->Fill(e4Vect.Py(), wgt);
	h1_eTrack_pz[0]->Fill(e4Vect.Pz(), wgt);
	h1_eTrack_p[0]->Fill(e4Vect.P(), wgt);
	h1_eTrack_E[0]->Fill(e4Vect.E(), wgt);
	h1_eTrack_Theta[0]->Fill(e4Vect.Theta()*TMath::RadToDeg(), wgt);
	h1_eTrack_Phi[0]->Fill(e4Vect.Phi()*TMath::RadToDeg(), wgt);
	h2_eTrack_ThetaPhi[0]->Fill((e4Vect.Theta()*TMath::RadToDeg()), (e4Vect.Phi()*TMath::RadToDeg()), wgt);
	h2_eTrack_pTheta[0]->Fill((e4Vect.Theta()*TMath::RadToDeg()), e4Vect.P(), wgt);

	for(Int_t B = 1; B < 7; B++){
	  Q2_low = 5+((B-1)*5.0);
	  Q2_high = 10+((B-1)*5.0);
	  if ( Q2_truth > Q2_low && Q2_truth < Q2_high){
	    h1_eTrack_px[B]->Fill(e4Vect.Px(), wgt);
	    h1_eTrack_py[B]->Fill(e4Vect.Py(), wgt);
	    h1_eTrack_pz[B]->Fill(e4Vect.Pz(), wgt);
	    h1_eTrack_p[B]->Fill(e4Vect.P(), wgt);
	    h1_eTrack_E[B]->Fill(e4Vect.E(), wgt);
	    h1_eTrack_Theta[B]->Fill(e4Vect.Theta()*TMath::RadToDeg(), wgt);
	    h1_eTrack_Phi[B]->Fill(e4Vect.Phi()*TMath::RadToDeg(), wgt);
	    h2_eTrack_ThetaPhi[B]->Fill((e4Vect.Theta()*TMath::RadToDeg()), (e4Vect.Phi()*TMath::RadToDeg()), wgt);
	    h2_eTrack_pTheta[B]->Fill((e4Vect.Theta()*TMath::RadToDeg()), e4Vect.P(), wgt);
	  }
	}
      }
  }

  if (Check_Pi(topNode) == true){
    SvtxTrackMap* pi_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");

    if (!pi_trackmap)
      {
	pi_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
	if (!pi_trackmap)
	  {
	    cout
	      << "ECCE_DEMP::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
	    exit(-1);
	  }
      }

    for (SvtxTrackMap::Iter pi_iter = pi_trackmap->begin();
	 pi_iter != pi_trackmap->end();
	 ++pi_iter)
      {
	SvtxTrack* pi_track = pi_iter->second;

	if ( pi_track->get_pz() > 0 && pi_track->get_charge() == 1){ // +ve z direction -> pions, crappy way of selecting them for now w/o truth info
	  piVect.SetXYZ(pi_track->get_px(), pi_track->get_py(), pi_track->get_pz());
	  pi4Vect.SetPxPyPzE(pi_track->get_px(), pi_track->get_py(), pi_track->get_pz(), sqrt(pow(piVect.Mag(), 2)+pow(mPi,2)));
	}
      }

    h1_piTrack_px[0]->Fill(pi4Vect.Px(), wgt);
    h1_piTrack_py[0]->Fill(pi4Vect.Py(), wgt);
    h1_piTrack_pz[0]->Fill(pi4Vect.Pz(), wgt);
    h1_piTrack_p[0]->Fill(pi4Vect.P(), wgt);
    h1_piTrack_E[0]->Fill(pi4Vect.E(), wgt);
    h1_piTrack_Theta[0]->Fill(pi4Vect.Theta()*TMath::RadToDeg(), wgt);
    h1_piTrack_Phi[0]->Fill(pi4Vect.Phi()*TMath::RadToDeg(), wgt);
    h2_piTrack_ThetaPhi[0]->Fill((pi4Vect.Theta()*TMath::RadToDeg()), (pi4Vect.Phi()*TMath::RadToDeg()), wgt);
    h2_piTrack_pTheta[0]->Fill((pi4Vect.Theta()*TMath::RadToDeg()), pi4Vect.P(), wgt);

    for(Int_t B = 1; B < 7; B++){
      Q2_low = 5+((B-1)*5.0);
      Q2_high = 10+((B-1)*5.0);
      if ( Q2_truth > Q2_low && Q2_truth < Q2_high){
	h1_piTrack_px[B]->Fill(pi4Vect.Px(), wgt);
	h1_piTrack_py[B]->Fill(pi4Vect.Py(), wgt);
	h1_piTrack_pz[B]->Fill(pi4Vect.Pz(), wgt);
	h1_piTrack_p[B]->Fill(pi4Vect.P(), wgt);
	h1_piTrack_E[B]->Fill(pi4Vect.E(), wgt);
	h1_piTrack_Theta[B]->Fill(pi4Vect.Theta()*TMath::RadToDeg(), wgt);
	h1_piTrack_Phi[B]->Fill(pi4Vect.Phi()*TMath::RadToDeg(), wgt);
	h2_piTrack_ThetaPhi[B]->Fill((pi4Vect.Theta()*TMath::RadToDeg()), (pi4Vect.Phi()*TMath::RadToDeg()), wgt);
	h2_piTrack_pTheta[B]->Fill((pi4Vect.Theta()*TMath::RadToDeg()), pi4Vect.P(), wgt);
      }
    }
  }

  // With ONLY ZDC hit, can't correct track
  if (Check_n(topNode) == true){ // Something that looks like a neutron found
    PHG4HitContainer* n_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDCsurrogate");
    if (n_hits) {
      // this returns an iterator to the beginning and the end of our G4Hits
      PHG4HitContainer::ConstRange n_hit_range = n_hits->getHits();
      for (PHG4HitContainer::ConstIterator n_hit_iter = n_hit_range.first; n_hit_iter != n_hit_range.second; n_hit_iter++)
	{	 
	  nZDCPos.SetXYZ(ZDC_Position_Smear(n_hit_iter->second->get_x(0)), ZDC_Position_Smear(n_hit_iter->second->get_y(0)), ZDC_Position_Smear(n_hit_iter->second->get_z(0)));
	  nEDep = ZDC_Energy_Smear_HCAL(n_hit_iter->second->get_edep());
	  nTheta = nZDCPos.Theta();
	  nPhi = nZDCPos.Phi();
	  nPMag = sqrt((pow(nEDep,2)) - (pow(mNeut,2)));
	  n4Vect.SetPxPyPzE(nPMag*sin(nTheta)*cos(nPhi), nPMag*sin(nTheta)*sin(nPhi), nPMag*cos(nTheta), nEDep);
	}
    }	 

    h1_nTrack_px[0]->Fill(n4Vect.Px(), wgt);
    h1_nTrack_py[0]->Fill(n4Vect.Py(), wgt);
    h1_nTrack_pz[0]->Fill(n4Vect.Pz(), wgt);
    h1_nTrack_p[0]->Fill(n4Vect.P(), wgt);
    h1_nTrack_E[0]->Fill(n4Vect.E(), wgt);
    h1_nTrack_Theta[0]->Fill(n4Vect.Theta()*TMath::RadToDeg(), wgt);
    h1_nTrack_Phi[0]->Fill(n4Vect.Phi()*TMath::RadToDeg(), wgt);
    h2_nTrack_ThetaPhi[0]->Fill((n4Vect.Theta()*TMath::RadToDeg()), (n4Vect.Phi()*TMath::RadToDeg()), wgt);
    h2_nTrack_pTheta[0]->Fill((n4Vect.Theta()*TMath::RadToDeg()), n4Vect.P(), wgt);

    for(Int_t B = 1; B < 7; B++){
      Q2_low = 5+((B-1)*5.0);
      Q2_high = 10+((B-1)*5.0);
      if ( Q2_truth > Q2_low && Q2_truth < Q2_high){
	h1_nTrack_px[B]->Fill(n4Vect.Px(), wgt);
	h1_nTrack_py[B]->Fill(n4Vect.Py(), wgt);
	h1_nTrack_pz[B]->Fill(n4Vect.Pz(), wgt);
	h1_nTrack_p[B]->Fill(n4Vect.P(), wgt);
	h1_nTrack_E[B]->Fill(n4Vect.E(), wgt);
	h1_nTrack_Theta[B]->Fill(n4Vect.Theta()*TMath::RadToDeg(), wgt);
	h1_nTrack_Phi[B]->Fill(n4Vect.Phi()*TMath::RadToDeg(), wgt);
	h2_nTrack_ThetaPhi[B]->Fill((n4Vect.Theta()*TMath::RadToDeg()), (n4Vect.Phi()*TMath::RadToDeg()), wgt);
	h2_nTrack_pTheta[B]->Fill((n4Vect.Theta()*TMath::RadToDeg()), n4Vect.P(), wgt);
      }
    }
  }

  if (Check_n(topNode) == true && Check_e(topNode) == true){ // If neutron and electron tracks look ok, plot some quantities
    SvtxTrackMap* enCoin_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    PHG4HitContainer* enCoin_hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDCsurrogate");

    if (!enCoin_trackmap)
      {
	enCoin_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
	if (!enCoin_trackmap)
	  {
	    cout
	      << "ECCE_DEMP_IR::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
	    exit(-1);
	  }
      }

    // Loop over our tracks, assign info to e'
    for (SvtxTrackMap::Iter enCoin_iter = enCoin_trackmap->begin();
	 enCoin_iter != enCoin_trackmap->end();
	 ++enCoin_iter)
      {
	SvtxTrack* enCoin_track = enCoin_iter->second;
	if (enCoin_track->get_pz() < 0  && enCoin_track->get_charge() == -1 ){ // -ve z direction -> electrons, crappy way of selecting them for now w/o truth info
	  eVect.SetXYZ(enCoin_track->get_px(), enCoin_track->get_py(), enCoin_track->get_pz());
	  e4Vect.SetPxPyPzE(enCoin_track->get_px(), enCoin_track->get_py(), enCoin_track->get_pz(), sqrt(pow(eVect.Mag(), 2)+pow(mElec,2)));
	}
      }

    if (enCoin_hits) {
      // this returns an iterator to the beginning and the end of our G4Hits
      PHG4HitContainer::ConstRange enCoin_hit_range = enCoin_hits->getHits();
      for (PHG4HitContainer::ConstIterator enCoin_hit_iter = enCoin_hit_range.first; enCoin_hit_iter != enCoin_hit_range.second; enCoin_hit_iter++)
	{	 
	  nZDCPos.SetXYZ(ZDC_Position_Smear(enCoin_hit_iter->second->get_x(0)), ZDC_Position_Smear(enCoin_hit_iter->second->get_y(0)), ZDC_Position_Smear(enCoin_hit_iter->second->get_z(0)));
	  nEDep = ZDC_Energy_Smear_HCAL(enCoin_hit_iter->second->get_edep());
	  nTheta = nZDCPos.Theta();
	  nPhi = nZDCPos.Phi();
	  nPMag = sqrt((pow(nEDep,2)) - (pow(mNeut,2)));
	  n4Vect.SetPxPyPzE(nPMag*sin(nTheta)*cos(nPhi), nPMag*sin(nTheta)*sin(nPhi), nPMag*cos(nTheta), nEDep);
	}
    }

    // Now fill our histograms
    h1_e_enCoin_px[0]->Fill(e4Vect.Px(), wgt);
    h1_e_enCoin_py[0]->Fill(e4Vect.Py(), wgt);
    h1_e_enCoin_pz[0]->Fill(e4Vect.Pz(), wgt);
    h1_e_enCoin_p[0]->Fill(e4Vect.P(), wgt);
    h1_e_enCoin_E[0]->Fill(e4Vect.E(), wgt);
    h1_e_enCoin_Theta[0]->Fill(e4Vect.Theta()*TMath::RadToDeg(), wgt);
    h1_e_enCoin_Phi[0]->Fill(e4Vect.Phi()*TMath::RadToDeg(), wgt);
    h2_e_enCoin_ThetaPhi[0]->Fill((e4Vect.Theta()*TMath::RadToDeg()), (e4Vect.Phi()*TMath::RadToDeg()), wgt);
    h2_e_enCoin_pTheta[0]->Fill((e4Vect.Theta()*TMath::RadToDeg()), e4Vect.P(), wgt);
    h1_n_enCoin_px[0]->Fill(n4Vect.Px(), wgt);
    h1_n_enCoin_py[0]->Fill(n4Vect.Py(), wgt);
    h1_n_enCoin_pz[0]->Fill(n4Vect.Pz(), wgt);
    h1_n_enCoin_p[0]->Fill(n4Vect.P(), wgt);
    h1_n_enCoin_E[0]->Fill(n4Vect.E(), wgt);
    h1_n_enCoin_Theta[0]->Fill(n4Vect.Theta()*TMath::RadToDeg(), wgt);
    h1_n_enCoin_Phi[0]->Fill(n4Vect.Phi()*TMath::RadToDeg(), wgt);
    h2_n_enCoin_ThetaPhi[0]->Fill((n4Vect.Theta()*TMath::RadToDeg()), (n4Vect.Phi()*TMath::RadToDeg()), wgt);
    h2_n_enCoin_pTheta[0]->Fill((n4Vect.Theta()*TMath::RadToDeg()), n4Vect.P(), wgt);

    h1_piTrue_enCoin_px[0]->Fill(pi4VectTruth.Px(), wgt);
    h1_piTrue_enCoin_py[0]->Fill(pi4VectTruth.Py(), wgt);
    h1_piTrue_enCoin_pz[0]->Fill(pi4VectTruth.Pz(), wgt);
    h1_piTrue_enCoin_p[0]->Fill(pi4VectTruth.P(), wgt);
    h1_piTrue_enCoin_E[0]->Fill(pi4VectTruth.E(), wgt);
    h1_piTrue_enCoin_Theta[0]->Fill(pi4VectTruth.Theta()*TMath::RadToDeg(), wgt);
    h1_piTrue_enCoin_Phi[0]->Fill(pi4VectTruth.Phi()*TMath::RadToDeg(), wgt);
    h2_piTrue_enCoin_ThetaPhi[0]->Fill((pi4VectTruth.Theta()*TMath::RadToDeg()), (pi4VectTruth.Phi()*TMath::RadToDeg()), wgt);
    h2_piTrue_enCoin_pTheta[0]->Fill((pi4VectTruth.Theta()*TMath::RadToDeg()), pi4VectTruth.P(), wgt);

    if(Check_Pi(topNode) == false){
      h1_piTrue_Missed_enCoin_px[0]->Fill(pi4VectTruth.Px(), wgt);
      h1_piTrue_Missed_enCoin_py[0]->Fill(pi4VectTruth.Py(), wgt);
      h1_piTrue_Missed_enCoin_pz[0]->Fill(pi4VectTruth.Pz(), wgt);
      h1_piTrue_Missed_enCoin_p[0]->Fill(pi4VectTruth.P(), wgt);
      h1_piTrue_Missed_enCoin_E[0]->Fill(pi4VectTruth.E(), wgt);
      h1_piTrue_Missed_enCoin_Theta[0]->Fill(pi4VectTruth.Theta()*TMath::RadToDeg(), wgt);
      h1_piTrue_Missed_enCoin_Phi[0]->Fill(pi4VectTruth.Phi()*TMath::RadToDeg(), wgt);
      h2_piTrue_Missed_enCoin_ThetaPhi[0]->Fill((pi4VectTruth.Theta()*TMath::RadToDeg()), (pi4VectTruth.Phi()*TMath::RadToDeg()), wgt);
      h2_piTrue_Missed_enCoin_pTheta[0]->Fill((pi4VectTruth.Theta()*TMath::RadToDeg()), pi4VectTruth.P(), wgt);
    }

    for(Int_t B = 1; B < 7; B++){
      Q2_low = 5+((B-1)*5.0);
      Q2_high = 10+((B-1)*5.0);
      if ( Q2_truth > Q2_low && Q2_truth < Q2_high){
	h1_e_enCoin_px[B]->Fill(e4Vect.Px(), wgt);
	h1_e_enCoin_py[B]->Fill(e4Vect.Py(), wgt);
	h1_e_enCoin_pz[B]->Fill(e4Vect.Pz(), wgt);
	h1_e_enCoin_p[B]->Fill(e4Vect.P(), wgt);
	h1_e_enCoin_E[B]->Fill(e4Vect.E(), wgt);
	h1_e_enCoin_Theta[B]->Fill(e4Vect.Theta()*TMath::RadToDeg(), wgt);
	h1_e_enCoin_Phi[B]->Fill(e4Vect.Phi()*TMath::RadToDeg(), wgt);
	h2_e_enCoin_ThetaPhi[B]->Fill((e4Vect.Theta()*TMath::RadToDeg()), (e4Vect.Phi()*TMath::RadToDeg()), wgt);
	h2_e_enCoin_pTheta[B]->Fill((e4Vect.Theta()*TMath::RadToDeg()), e4Vect.P(), wgt);
	h1_n_enCoin_px[B]->Fill(n4Vect.Px(), wgt);
	h1_n_enCoin_py[B]->Fill(n4Vect.Py(), wgt);
	h1_n_enCoin_pz[B]->Fill(n4Vect.Pz(), wgt);
	h1_n_enCoin_p[B]->Fill(n4Vect.P(), wgt);
	h1_n_enCoin_E[B]->Fill(n4Vect.E(), wgt);
	h1_n_enCoin_Theta[B]->Fill(n4Vect.Theta()*TMath::RadToDeg(), wgt);
	h1_n_enCoin_Phi[B]->Fill(n4Vect.Phi()*TMath::RadToDeg(), wgt);
	h2_n_enCoin_ThetaPhi[B]->Fill((n4Vect.Theta()*TMath::RadToDeg()), (n4Vect.Phi()*TMath::RadToDeg()), wgt);
	h2_n_enCoin_pTheta[B]->Fill((n4Vect.Theta()*TMath::RadToDeg()), n4Vect.P(), wgt);

	h1_piTrue_enCoin_px[B]->Fill(pi4VectTruth.Px(), wgt);
	h1_piTrue_enCoin_py[B]->Fill(pi4VectTruth.Py(), wgt);
	h1_piTrue_enCoin_pz[B]->Fill(pi4VectTruth.Pz(), wgt);
	h1_piTrue_enCoin_p[B]->Fill(pi4VectTruth.P(), wgt);
	h1_piTrue_enCoin_E[B]->Fill(pi4VectTruth.E(), wgt);
	h1_piTrue_enCoin_Theta[B]->Fill(pi4VectTruth.Theta()*TMath::RadToDeg(), wgt);
	h1_piTrue_enCoin_Phi[B]->Fill(pi4VectTruth.Phi()*TMath::RadToDeg(), wgt);
	h2_piTrue_enCoin_ThetaPhi[B]->Fill((pi4VectTruth.Theta()*TMath::RadToDeg()), (pi4VectTruth.Phi()*TMath::RadToDeg()), wgt);
	h2_piTrue_enCoin_pTheta[B]->Fill((pi4VectTruth.Theta()*TMath::RadToDeg()), pi4VectTruth.P(), wgt);

	if(Check_Pi(topNode) == false){
	  h1_piTrue_Missed_enCoin_px[B]->Fill(pi4VectTruth.Px(), wgt);
	  h1_piTrue_Missed_enCoin_py[B]->Fill(pi4VectTruth.Py(), wgt);
	  h1_piTrue_Missed_enCoin_pz[B]->Fill(pi4VectTruth.Pz(), wgt);
	  h1_piTrue_Missed_enCoin_p[B]->Fill(pi4VectTruth.P(), wgt);
	  h1_piTrue_Missed_enCoin_E[B]->Fill(pi4VectTruth.E(), wgt);
	  h1_piTrue_Missed_enCoin_Theta[B]->Fill(pi4VectTruth.Theta()*TMath::RadToDeg(), wgt);
	  h1_piTrue_Missed_enCoin_Phi[B]->Fill(pi4VectTruth.Phi()*TMath::RadToDeg(), wgt);
	  h2_piTrue_Missed_enCoin_ThetaPhi[B]->Fill((pi4VectTruth.Theta()*TMath::RadToDeg()), (pi4VectTruth.Phi()*TMath::RadToDeg()), wgt);
	  h2_piTrue_Missed_enCoin_pTheta[B]->Fill((pi4VectTruth.Theta()*TMath::RadToDeg()), pi4VectTruth.P(), wgt);
	}

      }    
    }
  }

  if (Check_n(topNode) == true && Check_ePi(topNode) == true){ // For event, check if it look like we have an e/pi/n in the event
    // Get track map for e'/pi info
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    // Get ZDC hits for neutron info
    //PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDC");
    // Need to use ZDC surrogate now
    PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDCsurrogate");
    
    if (!trackmap)
      {
	trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
	if (!trackmap)
	  {
	    cout
	      << "ECCE_DEMP_IR::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
	    exit(-1);
	  }
      }
    
    // Loop over our tracks, assign info to e'/pi
    for (SvtxTrackMap::Iter iter = trackmap->begin();
	 iter != trackmap->end();
	 ++iter)
      {
	SvtxTrack* track = iter->second;
	if ( track->get_pz() > 0 && track->get_charge() == 1){ // +ve z direction -> pions, crappy way of selecting them for now w/o truth info
	  piVect.SetXYZ(track->get_px(), track->get_py(), track->get_pz());
	  pi4Vect.SetPxPyPzE(track->get_px(), track->get_py(), track->get_pz(), sqrt(pow(piVect.Mag(), 2)+pow(mPi,2)));
	}

	else if (track->get_pz() < 0  && track->get_charge() == -1 ){ // -ve z direction -> electrons, crappy way of selecting them for now w/o truth info
	  eVect.SetXYZ(track->get_px(), track->get_py(), track->get_pz());
	  e4Vect.SetPxPyPzE(track->get_px(), track->get_py(), track->get_pz(), sqrt(pow(eVect.Mag(), 2)+pow(mElec,2)));
	}
      }
    
    // Loop over the hts in the zdc
    if (hits) {
      // this returns an iterator to the beginning and the end of our G4Hits
      PHG4HitContainer::ConstRange hit_range = hits->getHits();
      for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
	{	 
	  nZDCPos.SetXYZ(ZDC_Position_Smear(hit_iter->second->get_x(0)), ZDC_Position_Smear(hit_iter->second->get_y(0)), ZDC_Position_Smear(hit_iter->second->get_z(0)));
	  nEDep = ZDC_Energy_Smear_HCAL(hit_iter->second->get_edep());
	  nTheta = nZDCPos.Theta();
	  nPhi = nZDCPos.Phi();
	  nPMag = sqrt((pow(nEDep,2)) - (pow(mNeut,2)));
	  n4Vect.SetPxPyPzE(nPMag*sin(nTheta)*cos(nPhi), nPMag*sin(nTheta)*sin(nPhi), nPMag*cos(nTheta), nEDep);
	  
	  // SJDK - 15/02/22 - Implement Bill's conversion to local co-ordinates, fill new local ZDC XY histo
	  PHParameters ZDC_params{"PHG4RP"};
            
	  if (zdc_nodeparams){
	      ZDC_params.FillFrom(zdc_nodeparams, 0);
	  } 
	  else {
	    cerr << "There is a issue finding the detector paramter node!" << endl;
	  }

	  det_x_pos = Enclosure_params.get_double_param("place_x")  + ZDC_params.get_double_param("place_x");
	  det_z_pos = Enclosure_params.get_double_param("place_z")  + ZDC_params.get_double_param("place_z");

	  ZDC_params.set_double_param("place_x", det_x_pos); 
	  ZDC_params.set_double_param("place_z", det_z_pos); 

	  local_x = Get_Local_X(hit_iter->second->get_x(0), hit_iter->second->get_y(0), hit_iter->second->get_z(0), ZDC_params);
	  local_y = hit_iter->second->get_y(0);
	}
    }
    
    // Now have relevant information from this event, fill some histograms and calculate some stuff
    // Calculate kinematic quantities for electroproduction
    virtphoton4Vect = eBeam4Vect - e4Vect;
    //t4Vect = virtphoton4Vect - pi4Vect;
    pmiss4Vect = (eBeam4Vect + pBeam4Vect) - (e4Vect+pi4Vect); 
    // 31/03/22 - Check the difference between the reconstructed and measured neutron angles, veto on this later, for now, just want to plot it
    nTheta_Diff = pmiss4Vect.Theta() - nTheta;
    nPhi_Diff = pmiss4Vect.Phi() - nPhi;

    // Because the ZDC energy resolution is terrible, use the charged tracks to determine the momentum of the neutron
    // 01/04/22 - Construct reconstructed neutron 4-vector from Pmag of Pmiss and the real ZDC angles
    nRec4Vect.SetXYZM(pmiss4Vect.P()*sin(nTheta)*cos(nPhi), pmiss4Vect.P()*sin(nTheta)*sin(nPhi), pmiss4Vect.P()*cos(nTheta), mNeut);
    // 04/04/22 - Set reconstructed vector equal to neutron vector to test how bad getting neutron from ZDC would be
    //nRec4Vect = n4Vect; 

    nRecEDep = nRec4Vect.E();
    nRecPMag = nRec4Vect.P();
    nRecTheta = nRec4Vect.Theta();
    nRecPhi = nRec4Vect.Phi();
    nRecZDCPos.SetXYZ(nRecZDCPos.Z()*tan(nRecTheta)*cos(nRecPhi), nRecZDCPos.Z()*tan(nRecTheta)*sin(nRecPhi), nRecZDCPos.Z()); // Assume Z position from ZDC is true.

    t4Vect = eBeam4Vect - (e4Vect+pi4Vect);
    t_alt4Vect= pBeam4Vect - nRec4Vect;
    t_alt4Vect_ZDC = pBeam4Vect - n4Vect;
    pmiss4Vect_2 = (eBeam4Vect + pBeam4Vect) - (e4Vect+pi4Vect+nRec4Vect);
    Q2 = -1*(virtphoton4Vect.Mag2());
    W = (virtphoton4Vect+pBeam4Vect).Mag();
    t = -(t4Vect.Mag2());
    t_alt = -(t_alt4Vect.Mag2());
    t_alt_ZDC = -(t_alt4Vect_ZDC.Mag2());
    xb =  Q2/(2*(pBeam4Vect.Dot(virtphoton4Vect)));
    xi = xb/(2-xb);
    y_inv = (pBeam4Vect.Dot(virtphoton4Vect))/(pBeam4Vect.Dot(eBeam4Vect));// Calculation of the fractional energy loss
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP_IR::ResetEvent(PHCompositeNode *topNode)
{
//  std::cout << "ECCE_DEMP_IR::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
//
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP_IR::EndRun(const int runnumber)
{
  std::cout << "ECCE_DEMP_IR::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP_IR::End(PHCompositeNode *topNode)
{
  std::cout << "ECCE_DEMP_IR::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  ScalingFact = double(event_itt)/nTried; // This scaling factor is needed to normalise the weighted results
  for(Int_t C = 0; C < 7; C++){
    h1_eTrack_px[C]->Scale((1/ScalingFact));
    h1_eTrack_py[C]->Scale((1/ScalingFact));
    h1_eTrack_pz[C]->Scale((1/ScalingFact));
    h1_eTrack_p[C]->Scale((1/ScalingFact));
    h1_eTrack_E[C]->Scale((1/ScalingFact));
    h1_eTrack_Theta[C]->Scale((1/ScalingFact));
    h1_eTrack_Phi[C]->Scale((1/ScalingFact));
    h2_eTrack_ThetaPhi[C]->Scale((1/ScalingFact));
    h2_eTrack_pTheta[C]->Scale((1/ScalingFact));
    h1_piTrack_px[C]->Scale((1/ScalingFact));
    h1_piTrack_py[C]->Scale((1/ScalingFact));
    h1_piTrack_pz[C]->Scale((1/ScalingFact));
    h1_piTrack_p[C]->Scale((1/ScalingFact));
    h1_piTrack_E[C]->Scale((1/ScalingFact));
    h1_piTrack_Theta[C]->Scale((1/ScalingFact));
    h1_piTrack_Phi[C]->Scale((1/ScalingFact));
    h2_piTrack_ThetaPhi[C]->Scale((1/ScalingFact));
    h2_piTrack_pTheta[C]->Scale((1/ScalingFact));
    h1_nTrack_px[C]->Scale((1/ScalingFact));
    h1_nTrack_py[C]->Scale((1/ScalingFact));
    h1_nTrack_pz[C]->Scale((1/ScalingFact));
    h1_nTrack_p[C]->Scale((1/ScalingFact));
    h1_nTrack_E[C]->Scale((1/ScalingFact));
    h1_nTrack_Theta[C]->Scale((1/ScalingFact));
    h1_nTrack_Phi[C]->Scale((1/ScalingFact));
    h2_nTrack_ThetaPhi[C]->Scale((1/ScalingFact));
    h2_nTrack_pTheta[C]->Scale((1/ScalingFact));
    h1_e_enCoin_px[C]->Scale((1/ScalingFact));
    h1_e_enCoin_py[C]->Scale((1/ScalingFact));
    h1_e_enCoin_pz[C]->Scale((1/ScalingFact));
    h1_e_enCoin_p[C]->Scale((1/ScalingFact));
    h1_e_enCoin_E[C]->Scale((1/ScalingFact));
    h1_e_enCoin_Theta[C]->Scale((1/ScalingFact));
    h1_e_enCoin_Phi[C]->Scale((1/ScalingFact));
    h2_e_enCoin_ThetaPhi[C]->Scale((1/ScalingFact));
    h2_e_enCoin_pTheta[C]->Scale((1/ScalingFact));
    h1_n_enCoin_px[C]->Scale((1/ScalingFact));
    h1_n_enCoin_py[C]->Scale((1/ScalingFact));
    h1_n_enCoin_pz[C]->Scale((1/ScalingFact));
    h1_n_enCoin_p[C]->Scale((1/ScalingFact));
    h1_n_enCoin_E[C]->Scale((1/ScalingFact));
    h1_n_enCoin_Theta[C]->Scale((1/ScalingFact));
    h1_n_enCoin_Phi[C]->Scale((1/ScalingFact));
    h2_n_enCoin_ThetaPhi[C]->Scale((1/ScalingFact));
    h2_n_enCoin_pTheta[C]->Scale((1/ScalingFact));
    
    h1_piTrue_enCoin_px[C]->Scale((1/ScalingFact));
    h1_piTrue_enCoin_py[C]->Scale((1/ScalingFact));
    h1_piTrue_enCoin_pz[C]->Scale((1/ScalingFact));
    h1_piTrue_enCoin_p[C]->Scale((1/ScalingFact));
    h1_piTrue_enCoin_E[C]->Scale((1/ScalingFact));
    h1_piTrue_enCoin_Theta[C]->Scale((1/ScalingFact));
    h1_piTrue_enCoin_Phi[C]->Scale((1/ScalingFact));
    h2_piTrue_enCoin_ThetaPhi[C]->Scale((1/ScalingFact));
    h2_piTrue_enCoin_pTheta[C]->Scale((1/ScalingFact));

    h1_piTrue_Missed_enCoin_px[C]->Scale((1/ScalingFact));
    h1_piTrue_Missed_enCoin_py[C]->Scale((1/ScalingFact));
    h1_piTrue_Missed_enCoin_pz[C]->Scale((1/ScalingFact));
    h1_piTrue_Missed_enCoin_p[C]->Scale((1/ScalingFact));
    h1_piTrue_Missed_enCoin_E[C]->Scale((1/ScalingFact));
    h1_piTrue_Missed_enCoin_Theta[C]->Scale((1/ScalingFact));
    h1_piTrue_Missed_enCoin_Phi[C]->Scale((1/ScalingFact));
    h2_piTrue_Missed_enCoin_ThetaPhi[C]->Scale((1/ScalingFact));
    h2_piTrue_Missed_enCoin_pTheta[C]->Scale((1/ScalingFact));

  }

  outfile->cd();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP_IR::Reset(PHCompositeNode *topNode)
{
 std::cout << "ECCE_DEMP_IR::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void ECCE_DEMP_IR::Print(const std::string &what) const
{
  std::cout << "ECCE_DEMP_IR::Print(const std::string &what) const Printing info for " << what << std::endl;
}


///*****************************************************
/// ZDC Energy and Poisition smearing functions
//
// Energy smearing

float ECCE_DEMP_IR::ZDC_Energy_Smear_EMCAL(float E) {

  float resolution, E_reco;

  resolution = sqrt(.45*.45/E + 0.075*0.075);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}

// Energy smearing

float ECCE_DEMP_IR::ZDC_Energy_Smear_HCAL(float E) {

  float resolution, E_reco;

  //resolution = sqrt(.5*.5/E + 0.1*0.1); // YR Resolution
  resolution = sqrt(.45*.45/E + 0.042*0.042); // Updated Resolution
  //resolution = 0.25/sqrt(E); // Test 25% over sqrt E resolution
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}

// Energy smearing

float ECCE_DEMP_IR::ZDC_Energy_Smear_PbWO4(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

// Position smearing

float ECCE_DEMP_IR::ZDC_Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.15;         /// Position resolution 0.15 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}

///*****************************************************
/// B0 tracker smearing functions

// Energy smearing

float ECCE_DEMP_IR::B0Tracker_Energy_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

// Posision smearing

float ECCE_DEMP_IR::B0Tracker_Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.1;         /// Position resolution 0.1 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}


///*****************************************************
/// B0 Cal smearing functions

// Energy smearing

float ECCE_DEMP_IR::B0Cal_Energy_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

// Posision smearing

float ECCE_DEMP_IR::B0Cal_Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.1;         /// Position resolution 0.1 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}

//***************************************************

bool ECCE_DEMP_IR::Check_ePi(PHCompositeNode* topNode)
{
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap)
    {
      trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
      if (!trackmap)
    	{
    	  cout
    	    << "ECCE_DEMP_IR::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
    	  exit(-1);
    	}
    }
  int nTracks = 0;
  Bool_t ElecTrack = kFALSE;
  Bool_t PionTrack = kFALSE;
  // Iterate over tracks
  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter)
    {
      SvtxTrack* track = iter->second;
      nTracks++;
      if ( track->get_pz() > 0 && track->get_charge() == 1){ // +ve z direction -> pions, crappy way of selecting them for now w/o truth info
	PionTrack = kTRUE;
      }
      else if (track->get_pz() < 0  && track->get_charge() == -1 ){ // -ve z direction -> electrons, crappy way of selecting them for now w/o truth info
	ElecTrack = kTRUE;
      }
    }
  
  if( PionTrack == kTRUE && ElecTrack == kTRUE && nTracks == 2){ // Both a pion and an electron track, only 2 tracks
    return true;
  }
  else{
    return false;
  }
}

//***************************************************
bool ECCE_DEMP_IR::Check_e(PHCompositeNode* topNode)
{
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap)
    {
      trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
      if (!trackmap)
    	{
    	  cout
    	    << "ECCE_DEMP_IR::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
    	  exit(-1);
    	}
    }
  int nTracks_e = 0;
  Bool_t SingElecTrack = kFALSE;
  // Iterate over tracks
  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter)
    {
      SvtxTrack* track = iter->second;
      nTracks_e++;
      if (track->get_pz() < 0  && track->get_charge() == -1){ // -ve z direction -> electrons, crappy way of selecting them for now w/o truth info
	SingElecTrack = kTRUE;
      }
    }
  
  //if(SingElecTrack == kTRUE && nTracks_e == 2){ // An electron track was found
  if(SingElecTrack == kTRUE){ // An electron track was found
    return true;
  }
  else{
    return false;
  }
}
//***************************************************
bool ECCE_DEMP_IR::Check_Pi(PHCompositeNode* topNode)
{
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap)
    {
      trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
      if (!trackmap)
    	{
    	  cout
    	    << "ECCE_DEMP_IR::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
    	  exit(-1);
    	}
    }
  int nTracks_pi = 0;
  Bool_t SingPionTrack = kFALSE;
  // Iterate over tracks
  for (SvtxTrackMap::Iter iter = trackmap->begin();
       iter != trackmap->end();
       ++iter)
    {
      SvtxTrack* track = iter->second;
      nTracks_pi++;
      if ( track->get_pz() > 0 && track->get_charge() == 1){ // +ve z direction -> pions, crappy way of selecting them for now w/o truth info
	SingPionTrack = kTRUE;
      }
    }
  
  //if(SingPionTrack == kTRUE && nTracks_pi == 2){ // A pion track was found
  if(SingPionTrack == kTRUE){ // A pion track was found
    return true;
  }
  else{
    return false;
  }
}

//***************************************************
bool ECCE_DEMP_IR::Check_n(PHCompositeNode* topNode)
{
  // loop over the G4Hits

  //PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDC");
  // Need to use ZDC surrogate now
  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDCsurrogate");

  int ZDCHits = 0;
  Bool_t nZDChit = kFALSE;

  if (hits) {
    // this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
      {
	ZDCHits++;
	//if (hit_iter->second->get_edep() > 40){ // Hit in ZDC with roughly correct energy for neutron, should check the smeared energy instead?
	if (ZDC_Energy_Smear_HCAL(hit_iter->second->get_edep()) > 40){ // Hit in ZDC with roughly correct energy for neutron, use smeared energy
	  nZDChit = kTRUE;
	}
      }
  }

  if( nZDChit == kTRUE ){ // Hit in ZDC with correct energy
    return true;
  }
  else if (nZDChit == kFALSE){
    return false;
  }
  else return false;
}


// SJDK - 15/02/22 - Added in new co-ordinate conversions
// Conversion to local co-ordinate system from global values - from Bill's Diff_Tagg_Ana scripts
//*******************************************

float ECCE_DEMP_IR::Get_Local_X(float global_x, float global_y, float global_z, float det_tilt, float det_rot) {


   TVector3 global_cor(global_x, global_y, global_z);
   float local_x;

   global_cor.RotateY(-det_rot);
   local_x = global_cor.X()/cos(det_tilt - det_rot);
	
   return local_x;

}

//*******************************************

float ECCE_DEMP_IR::Get_Local_Y(float global_x, float global_y, float global_z, float det_tilt, float cross_angle) {

	return global_y;

}

//*******************************************
float ECCE_DEMP_IR::Get_Local_X(float global_x, float global_y, float global_z, PdbParameterMapContainer *det_nodeparams) {

   PHParameters Det_params{"PHDet"};

   if (det_nodeparams)
   {
      Det_params.FillFrom(det_nodeparams, 0);
   } else {
      cerr << "There is a issue finding the detector paramter node!" << endl;
   }

   float det_xCent = Enclosure_params.get_double_param("place_x") + Det_params.get_double_param("place_x");
   float det_zCent = Enclosure_params.get_double_param("place_z") + Det_params.get_double_param("place_z");
   float det_tilt = Det_params.get_double_param("rot_y")/180. * TMath::Pi(); // in Rad

   float det_rot = atan( det_xCent / det_zCent);  // in Rad

   TVector3 global_cor(global_x, global_y, global_z);
   float local_x;

   global_cor.RotateY(-det_rot);
   local_x = global_cor.X()/cos(det_tilt);

   return local_x;

}

//*******************************************

float ECCE_DEMP_IR::Get_Local_X(float global_x, float global_y, float global_z, PHParameters Det_params) {

   float det_xCent = Det_params.get_double_param("place_x");
   float det_zCent = Det_params.get_double_param("place_z");

   float det_tilt = Det_params.get_double_param("rot_y"); // in Rad

   float det_rot = atan( det_xCent / det_zCent);  // in Rad

   TVector3 global_cor(global_x, global_y, global_z);


   float local_x1 = Get_Local_X(global_x, global_y, global_z, det_tilt, det_rot);

   return local_x1;

}

