#include "ECCE_DEMP.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/PHTFileServer.h>

#include <phool/PHCompositeNode.h>

#include <stdio.h>

#include <fun4all/Fun4AllHistoManager.h>

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

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

using namespace std;

ECCE_DEMP::ECCE_DEMP(const std::string &name, const std::string& filename):
 SubsysReco(name)
 , outfilename(filename)
{
  std::cout << "ECCE_DEMP_example::Diff_Tagg_example(const std::string &name) Calling ctor" << std::endl;

  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(m_RandomGenerator, seed);

}

//____________________________________________________________________________..
ECCE_DEMP::~ECCE_DEMP()
{

  gsl_rng_free(m_RandomGenerator);

  std::cout << "ECCE_DEMP::~ECCE_DEMP() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int ECCE_DEMP::Init(PHCompositeNode *topNode)
{
  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(outfilename.c_str(), "RECREATE");

  std::cout << "ECCE_DEMP::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  event_itt = 0;
  // In future, the value below should be determined more accurately - average value of all 100 files?
  nTried = 170000; // This is the approximate total number of events generated in a single file of this sample (5on100 epi)
  
  gDirectory->mkdir("Detection_Efficiency");
  gDirectory->cd("Detection_Efficiency");
  h1_Q2_DetEff_Uncut = new TH1F("Q2_DetEff_Uncut", "Q^{2}_{Truth} for thrown events; Q^{2}", 100, 0, 50); h1_Q2_DetEff_Uncut->Sumw2();
  h1_Q2_DetEff_Cut = new TH1F("Q2_DetEff_Cut", "Q^{2}_{Truth} for detected events; Q^{2}", 100, 0, 50); h1_Q2_DetEff_Cut->Sumw2();
  h1_Q2_DetEff = new TH1F("Q2_DetEff", "Q^{2}_{Truth} detected/thrown ratio; Q^{2}", 100, 0, 50); h1_Q2_DetEff->Sumw2();
  h2_Q2_t_DetEff_Uncut = new TH2F("Q2_t_DetEff_Uncut", "Q^{2}_{Truth} vs -t_{Truth} for thrown events; Q^{2}; -t", 10, 0, 50, 10, 0, 0.5); h2_Q2_t_DetEff_Uncut->Sumw2();
  h2_Q2_t_DetEff_Cut = new TH2F("Q2_t_DetEff_Cut", "Q^{2}_{Truth} vs -t_{Truth} for detected events; Q^{2}; -t", 10, 0, 50, 10, 0, 0.5); h2_Q2_t_DetEff_Cut->Sumw2();
  h2_Q2_t_DetEff = new TH2F("Q2_t_DetEff", "Q^{2}_{Truth} vs -t_{Truth} detected/thrown ratio; Q^{2}; -t", 10, 0, 50, 10, 0, 0.5); h2_Q2_t_DetEff->Sumw2();
  gDirectory->cd("../");  

  gDirectory->mkdir("Particle_Momenta_Resolution");
  gDirectory->cd("Particle_Momenta_Resolution");
  h1_piRes_p = new TH1F("piRes_p", "#pi #frac{#Delta p}{Truth p} Distribution (%); %", 100, -50, 50);
  h1_piRes_px = new TH1F("piRes_px", "#pi #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_piRes_py = new TH1F("piRes_py", "#pi #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_piRes_pz = new TH1F("piRes_pz", "#pi #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_eRes_p = new TH1F("eRes_p", "e' #frac{#Delta p}{Truth p} Distribution (%); %", 100, -50, 50);
  h1_eRes_px = new TH1F("eRes_px", "e' #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_eRes_py = new TH1F("eRes_py", "e' #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_eRes_pz = new TH1F("eRes_pz", "e' #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  h1_nRes_p = new TH1F("nRes_p", "n #frac{#Delta p}{Truth p} Distribution (%); %", 100, -50, 50);
  h1_nRes_px = new TH1F("nRes_px", "n #frac{#Delta px}{Truth px} Distribution (%); %", 100, -50, 50);
  h1_nRes_py = new TH1F("nRes_py", "n #frac{#Delta py}{Truth py} Distribution (%); %", 100, -50, 50);
  h1_nRes_pz = new TH1F("nRes_pz", "n #frac{#Delta pz}{Truth pz} Distribution (%); %", 100, -50, 50);
  gDirectory->cd("../");  

  gDirectory->mkdir("Pion_Info");
  gDirectory->cd("Pion_Info");
  h1_pi_px = new TH1F("pi_px", "#pi p_{x} Distribution;p_{x} [GeV]", 200, -20, 20);
  h1_pi_py = new TH1F("pi_py", "#pi p_{y} Distribution;p_{y} [GeV]", 200, -20, 20);
  h1_pi_pz = new TH1F("pi_pz", "#pi p_{z} Distribution;p_{z} [GeV]", 200, -50, 50); 
  h1_pi_p = new TH1F("pi_p", "#pi p Distribution;p [GeV]", 200, 0, 50);
  h1_pi_E = new TH1F("pi_E", "#pi E Distribution;E [GeV]", 200, 0, 50);
  h1_pi_Theta = new TH1F("pi_Theta", "#pi #theta Distribution; #theta [deg]", 200, 0, 50);
  h1_pi_Phi = new TH1F("pi_Phi", "#pi #phi Distribution; #phi [deg]", 360, -180, 180);
  h2_piTrack_ThetaPhi = new TH2F("piTrack_ThetaPhi", "#pi Track #theta vs #phi; #theta [deg]; #phi [deg]", 120, 0, 60, 720, -180, 180);
  h2_piTrack_pTheta = new TH2F("piTrack_pTheta", "#pi Track #theta vs P; #theta [deg]; P [GeV/c]", 120, 0, 60, 500, 0, 50);
  gDirectory->cd("../");

  gDirectory->mkdir("Pion_Truth_Info");
  gDirectory->cd("Pion_Truth_Info");
  h1_piTruth_px = new TH1F("piTrtuh_px", "#pi Truth p_{x} Distribution;p_{x} [GeV]", 200, -20, 20);
  h1_piTruth_py = new TH1F("piTrtuh_py", "#pi Truth p_{y} Distribution;p_{y} [GeV]", 200, -20, 20);
  h1_piTruth_pz = new TH1F("piTrtuh_pz", "#pi Truth p_{z} Distribution;p_{z} [GeV]", 200, -50, 50); 
  h1_piTruth_p = new TH1F("piTrtuh_p", "#pi Truth p Distribution;p [GeV]", 200, 0, 50);
  h1_piTruth_E = new TH1F("piTrtuh_E", "#pi Truth E Distribution;E [GeV]", 200, 0, 50);
  h1_piTruth_Theta = new TH1F("piTrtuh_Theta", "#pi Truth #theta Distribution; #theta [deg]", 200, 0, 50);

  h2_piTruth_pxpy = new TH2F("piTruth_pxpy", "#pi #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  h2_piTruth_pxpz = new TH2F("piTruth_pxpz", "#pi #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_piTruth_pypz = new TH2F("piTruth_pypz", "#pi #frac{#Delta p_{y}}{Truth p_{y}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  gDirectory->cd("../");
  
  gDirectory->mkdir("Scattered_Electron_Info");
  gDirectory->cd("Scattered_Electron_Info");
  h1_e_px = new TH1F("e_px", "e' p_{x} Distribution;p_{x} [GeV]", 240, -6, 6);
  h1_e_py = new TH1F("e_py", "e' p_{y} Distribution;p_{y} [GeV]", 240, -6, 6);
  h1_e_pz = new TH1F("e_pz", "e' p_{z} Distribution;p_{z} [GeV]", 120, -6, 0); 
  h1_e_p = new TH1F("e_p", "e' p Distribution;p [GeV]", 160, 0, 8);
  h1_e_E = new TH1F("e_E", "e' E Distribution;E [GeV]", 160, 0, 8);
  h1_e_Theta = new TH1F("e_Theta", "e' #theta Distribution; #theta [deg]", 200, 110, 160);
  h1_e_Phi = new TH1F("e_Phi", "e' #phi Distribution; #phi [deg]", 360, -180, 180);
  h2_eTrack_ThetaPhi = new TH2F("eTrack_ThetaPhi", "e' Track #theta vs #phi; #theta [deg]; #phi [deg]", 140, 110, 180, 720, -180, 180);
  h2_eTrack_pTheta = new TH2F("eTrack_pTheta", "e' Track #theta vs P; #theta [deg]; P [GeV/c]", 140, 110, 180, 100, 0, 10);
  h2_pi_XY = new TH2F("pi_XY", "#{pi} X vs Y at z=100cm Dist; x(cm); y(cm)", 200, -1000, 1000, 200, -1000, 1000);
  gDirectory->cd("../");

  gDirectory->mkdir("Scattered_Electron_Truth_Info");
  gDirectory->cd("Scattered_Electron_Truth_Info");
  h1_eTruth_px = new TH1F("eTruth_px", "e' Truth p_{x} Distribution;p_{x} [GeV]", 240, -6, 6);
  h1_eTruth_py = new TH1F("eTruth_py", "e' Truth p_{y} Distribution;p_{y} [GeV]", 240, -6, 6);
  h1_eTruth_pz = new TH1F("eTruth_pz", "e' Truth p_{z} Distribution;p_{z} [GeV]", 120, -6, 0); 
  h1_eTruth_p = new TH1F("eTruth_p", "e' Truth p Distribution;p [GeV]", 160, 0, 8);
  h1_eTruth_E = new TH1F("eTruth_E", "e' Truth E Distribution;E [GeV]", 160, 0, 8);
  h1_eTruth_Theta = new TH1F("eTruth_Theta", "e' Truth #theta Distribution; #theta [deg]", 200, 110, 160);
  h2_eTruth_pxpy = new TH2F("eTruth_pxpy", "e' #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);  
  h2_eTruth_pxpz = new TH2F("eTruth_pxpz", "e' #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);  
  h2_eTruth_pypz = new TH2F("eTruth_pypz", "e' #frac{#Delta p_{y}}{Truth p_{y}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_e_XY = new TH2F("e_XY", "e' X vs Y at z=100cm Dist; x(cm); y(cm)", 200, -1000, 1000, 200, -1000, 1000);
  gDirectory->cd("../");

  gDirectory->mkdir("Neutron_Info");
  gDirectory->cd("Neutron_Info");
  h1_n_px = new TH1F("n_px", "n p_{x} Distribution;p_{x} [GeV]", 480, -6, 6);
  h1_n_py = new TH1F("n_py", "n p_{y} Distribution;p_{y} [GeV]", 200, -2.5, 2.5);
  h1_n_pz = new TH1F("n_pz", "n p_{z} Distribution;p_{z} [GeV]", 240, 0, 120); 
  h1_n_p = new TH1F("n_p", "n p Distribution;p [GeV]", 240, 0, 120);
  h1_n_E = new TH1F("n_E", "n E Distribution;E [GeV]", 240, 0, 120);
  h1_n_Theta = new TH1F("n_Theta", "n #theta Distribution; #theta [deg]", 500, 0, 5);
  h1_n_Phi = new TH1F("n_Phi", "n #phi Distribution; #phi [deg]", 360, -180, 180);
  h2_nTrack_ThetaPhi = new TH2F("nTrack_ThetaPhi", "n Track #theta vs #phi; #theta [deg]; #phi [deg]", 500, 0, 5, 360, -180, 180);
  h2_nTrack_pTheta = new TH2F("nTrack_pTheta", "n Track #theta vs P; #theta [deg]; P [GeV/c]", 500, 0, 5, 1000, 0, 100);
  h2_n_XY = new TH2F("n_XY", "n X vs Y at ZDC Dist; x(cm); y(cm)", 200, -150, -50, 200, -50, 50);
  gDirectory->cd("../");

  gDirectory->mkdir("Neutron_Truth_Info");
  gDirectory->cd("Neutron_Truth_Info");
  h1_nTruth_px = new TH1F("nTruth_px", "nTruth p_{x} Distribution;p_{x} [GeV]", 480, -6, 6);
  h1_nTruth_py = new TH1F("nTruth_py", "nTruth p_{y} Distribution;p_{y} [GeV]", 200, -2.5, 2.5);
  h1_nTruth_pz = new TH1F("nTruth_pz", "nTruth p_{z} Distribution;p_{z} [GeV]", 240, 0, 120); 
  h1_nTruth_p = new TH1F("nTruth_p", "nTruth p Distribution;p [GeV]", 240, 0, 120);
  h1_nTruth_E = new TH1F("nTruth_E", "nTruth E Distribution;E [GeV]", 240, 0, 120);
  h1_nTruth_Theta = new TH1F("nTruth_Theta", "nTruth #theta Distribution; #theta [deg]", 500, 0, 5);
  h2_nTruth_pxpy = new TH2F("nTruth_pxpy", "n #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  h2_nTruth_pxpz = new TH2F("nTruth_pxpz", "n #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_nTruth_pypz = new TH2F("nTruth_pypz", "n #frac{#Delta p_{y}}{Truth p_{y}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  gDirectory->cd("../");

  gDirectory->mkdir("PMiss_Info");
  gDirectory->cd("PMiss_Info");
  h1_pmiss_px = new TH1F("pmiss_px", "p_{miss} p_{x} Distribution", 800, -10, 10);
  h1_pmiss_py = new TH1F("pmiss_py", "p_{miss} p_{y} Distribution", 200, -2.5, 2.5);
  h1_pmiss_pz = new TH1F("pmiss_pz", "p_{miss} p_{z} Distribution", 240, 0, 120); 
  h1_pmiss_p = new TH1F("pmiss_p", "p_{miss} p Distribution", 240, 0, 120);
  h1_pmiss_E = new TH1F("pmiss_E", "p_{miss} E Distribution", 240, 0, 120);
  h1_pmiss_Theta = new TH1F("pmiss_Theta", "p_{miss} #theta Distribution; #theta [deg]", 1000, 0, 10);
  h1_pmiss_Phi = new TH1F("pmiss_Phi", "p_{miss} #phi Distribution; #phi [deg]", 720, -180, 180);
  gDirectory->cd("../");
  
  gDirectory->mkdir("Virtual_Photon_Info");
  gDirectory->cd("Virtual_Photon_Info");
  h1_gamma_px = new TH1F("gamma_px", "#gamma p_{x} Distribution", 200, -10, 10);
  h1_gamma_py = new TH1F("gamma_py", "#gamma p_{y} Distribution", 200, -10, 10);
  h1_gamma_pz = new TH1F("gamma_pz", "#gamma p_{z} Distribution", 200, -10, 0); 
  h1_gamma_p = new TH1F("gamma_p", "#gamma p Distribution", 200, 0, 10);
  h1_gamma_E = new TH1F("gamma_E", "#gamma E Distribution", 200, 0, 10);
  h1_gamma_Theta = new TH1F("gamma_Theta", "#gamma #theta Distribution; #theta [deg]", 360, -180, 180);
  h1_gamma_Phi = new TH1F("gamma_Phi", "#gamma #phi Distribution; #phi [deg]", 360, -180, 180);
  gDirectory->cd("../");
 
  gDirectory->mkdir("Kinematics_Info");
  gDirectory->cd("Kinematics_Info");
  h1_Q2_Dist = new TH1F("Q2_Dist", "Q^{2} Distribution", 200, 0, 50);
  h1_W_Dist = new TH1F("W_Dist", "W Distribution", 500, 0, 50);
  h1_t_Dist = new TH1F("t_Dist", "t Distribution", 100, 0, 10);
  h1_t_alt_Dist = new TH1F("t_alt_Dist", "t (Alternative calculation) Distribution", 100, 0, 1);
  h1_t_comp = new TH1F("t_comp_Dist", "#frac{#Delta t}{t} Distribution; #frac{t_{alt}-t}{t} (%)", 200, -100, 100);
  h1_xb_Dist = new TH1F("xb_Dist", "x_{b} Distribution", 100, 0, 1);
  h1_xi_Dist = new TH1F("xi_Dist", "#xi Distribution", 100, 0, 1);
  gDirectory->cd("../");

  gDirectory->mkdir("Kinematics_Truth_Info");
  gDirectory->cd("Kinematics_Truth_Info");
  h1_Q2Truth_Dist = new TH1F("Q2Truth_Dist", "Q^{2} Truth Distribution", 200, 0, 50);
  h1_WTruth_Dist = new TH1F("WTruth_Dist", "W Truth Distribution", 500, 0, 50);
  h1_tTruth_Dist = new TH1F("tTruth_Dist", "t Truth Distribution", 100, 0, 1);
  h1_t_altTruth_Dist = new TH1F("t_altTruth_Dist", "-t_alt Truth (Alternative Calculation) Distribution", 100, 0, 1);
  h1_xbTruth_Dist = new TH1F("xbTruth_Dist", "x_{b} Truth Distribution", 100, 0, 1);
  h1_xiTruth_Dist = new TH1F("xiTruth_Dist", "#xi Truth Distribution", 100, 0, 1);
  gDirectory->cd("../");

  gDirectory->mkdir("Kinematics_Analysis");
  gDirectory->cd("Kinematics_Analysis");
  h2_t_ep = new TH2F("t_ep", "t vs ScatElec P; t; P_{e'}", 100, 0, 10, 200, 0, 10);
  h2_t_Q2 = new TH2F("t_Q2", "t vs Q^{2}; t; Q^{2}", 100, 0, 10, 200, 0, 50);
  h2_delta_t_t = new TH2F("delta_t_t", "#Delta t vs t; #Delta t (%); t", 200, -100, 100, 100, 0, 1);
  
  for(Int_t A = 0; A < 7; A++){
    h1_t_Q2[A] = new TH1F(Form("t_Q2_%i", (A+1)), Form("t dist, %i < Q^{2} < %i; t", (5 + (A*5)), 10+(A*5)), 100, 0, 10);
    h1_t_alt_Q2[A] = new TH1F(Form("t_alt_Q2_%i", (A+1)), Form("t (Alternative calculation) dist, %i < Q^{2} < %i; t", (5 + (A*5)), 10+(A*5)), 100, 0, 10);
    h2_delta_t_t_Q2[A] = new TH2F(Form("delta_t_t_Q2_%i", (A+1)), Form("#Delta t vs t, %i < Q^{2} < %i; #Delta t (Percent); t", (5 + (A*5)), 10+(A*5)), 200, -100, 100, 100, 0, 1);
  }
  gDirectory->cd("../");

  gDirectory->mkdir("Physics_Results_Cuts");
  gDirectory->cd("Physics_Results_Cuts");
  // Results histograms are binned in Q2, first two bins are special
  for(Int_t A = 0; A < 8; A++){
    if ( A <= 1){
      h1_t_result[A] = new TH1F(Form("t_Result_Q2_%i", (A+1)), Form("-t Dist, %2.1f < Q^{2} < %2.1f; -t", (5+(A*2.5)), (7.5+(A*2.5))), 10, 0, 0.4);
      h1_nTheta_result[A] = new TH1F(Form("nTheta_Result_Q2_%i", (A+1)), Form("#theta_{n} Dist, %2.1f < Q^{2} < %2.1f", (5+(A*2.5)), (7.5+(A*2.5))), 500, 0, 5);
      h1_pmiss_result[A] = new TH1F(Form("pmiss_Result_Q2_%i", (A+1)), Form("p_{miss} Dist, %2.1f < Q^{2} < %2.1f", (5+(A*2.5)), (7.5+(A*2.5))), 240, 0, 120) ;
    }
    else{
      h1_t_result[A] = new TH1F(Form("t_Result_Q2_%i", (A+1)), Form("-t Dist, %i < Q^{2} < %i; -t", (5+((A-1)*5)), (10+((A-1)*5))), 10, 0, 0.4);
      h1_nTheta_result[A] = new TH1F(Form("nTheta_Result_Q2_%i", (A+1)), Form("#theta_{n} Dist, %i < Q^{2} < %i", (5+((A-1)*5)), (10+((A-1)*5))), 500, 0, 5);
      h1_pmiss_result[A] = new TH1F(Form("pmiss_Result_Q2_%i", (A+1)), Form("p_{miss} Dist, %i < Q^{2} < %i", (5+((A-1)*5)), (10+((A-1)*5))), 240, 0, 120) ;
    }
  }

  gDirectory->cd("../");
  gDirectory->mkdir("Physics_Results");
  gDirectory->cd("Physics_Results");
  for(Int_t A = 0; A < 8; A++){
    if ( A <= 1){
      h1_t_cut_result[A] = new TH1F(Form("t_cut_Result_Q2_%i", (A+1)), Form("-t Dist, %2.1f < Q^{2} < %2.1f, with p_{miss}, #theta_{n} cuts; -t (GeV^{2}); Rate(Hz)", (5+(A*2.5)), (7.5+(A*2.5))), 10, 0, 0.4);
      h1_Q2_cut_result[A] = new TH1F(Form("Q2_cut_Result_Q2_%i", (A+1)), Form("Q^{2} Dist, %2.1f < Q^{2} < %2.1f, with p_{miss}, #theta_{n} cuts; Q^{2}(GeV^{2}); Rate(Hz)", (5+(A*2.5)), (7.5+(A*2.5))), 25, (5+(A*2.5)), (7.5+(A*2.5)));
      h1_W_cut_result[A] = new TH1F(Form("W_cut_Result_Q2_%i", (A+1)), Form("W Dist, %2.1f < Q^{2} < %2.1f, with p_{miss}, #theta_{n} cuts; W(GeV); Rate(Hz)", (5+(A*2.5)), (7.5+(A*2.5))), 60, 0, 30);
      h2_Q2_t_result[A] = new TH2F(Form("Q2_t_cut_Result_Q2_%i", (A+1)), Form("Q^{2} vs -t Dist, %2.1f < Q^{2} < %2.1f, with p_{miss}, #theta_{n} cuts; Q^{2}(GeV^{2}); -t (GeV^{2})", (5+(A*2.5)), (7.5+(A*2.5))), 25, (5+(A*2.5)), (7.5+(A*2.5)), 10, 0, 0.4);
    }
    else{
      h1_t_cut_result[A] = new TH1F(Form("t_cut_Result_Q2_%i", (A+1)), Form("-t Dist, %i < Q^{2} < %i, with p_{miss}, #theta_{n} cuts; -t (GeV^{2}); Rate(Hz)", (5+((A-1)*5)), (10+((A-1)*5))), 10, 0, 0.4);;
      h1_Q2_cut_result[A] = new TH1F(Form("Q2_cut_Result_Q2_%i", (A+1)), Form("Q^{2} Dist, %i < Q^{2} < %i, with p_{miss}, #theta_{n} cuts; Q^{2}(GeV^{2}); Rate(Hz)", (5+((A-1)*5)), (10+((A-1)*5))), 50, (5+((A-1)*5)), (10+((A-1)*5)));
      h1_W_cut_result[A] = new TH1F(Form("W_cut_Result_Q2_%i", (A+1)), Form("W Dist, %i < Q^{2} < %i, with p_{miss}, #theta_{n} cuts; W(GeV); Rate(Hz)", (5+((A-1)*5)), (10+((A-1)*5))), 60, 0, 30);
      h2_Q2_t_result[A] = new TH2F(Form("Q2_t_cut_Result_Q2_%i", (A+1)), Form("Q^{2} vs -t Dist, %2.1f < Q^{2} < %2.1f, with p_{miss}, #theta_{n} cuts; Q^{2}(GeV^{2}); -t (GeV^{2})", (5+(A*2.5)), (7.5+(A*2.5))), 50, (5+((A-1)*5)), (10+((A-1)*5)), 10, 0, 0.4);
    }
  }

  gDirectory->cd("../");
  gDirectory->mkdir("Physics_Results_Misc");
  gDirectory->cd("Physics_Results_Misc");
  h2_Q2_W_result = new TH2F("Q2_W_Result", "Q^{2} vs W, with p_{miss}, #theta_{n} cuts; Q^{2}(GeV^{2}); W (GeV)", 200, 0, 40, 60, 0, 30);
  h2_t_t_alt_result = new TH2F("t_t_alt_result", "-t vs -t_{alt}, with p_{miss}, #theta_{n} cuts; -t (GeV^{2}); -t_{alt}(GeV^{2})", 50, 0, 10, 12, 0, 0.48); 
  h1_Mmiss_result = new TH1F("Mmiss_result", "Missing Mass Dist;M_{Miss}(GeV/c^{2})", 100, -5, 5);
  h1_Mmiss_truth_result = new TH1F("Mmiss_truth_result", "Missing Mass (truth) Dist; M_{Miss}(GeV/c^{2})", 40, -2, 2);
  h1_Mmiss_Comp_result = new TH1F("Mmiss_Comp_result", "#Delta M_{Miss} Distribution; #Delta(M_{Miss})(GeV/c^{2})", 100, -5, 5);
  h2_t_ttruth_result = new TH2F("t_ttruth_result", "-t vs -t_{truth} Dist; -t (GeV^{2}); -t_{truth}(GeV^{2})", 50, 0, 10, 50, 0, 0.5);
  h2_t_alt_ttruth_result = new TH2F("t_alt_ttruth_result", "-t_{alt} vs -t_{truth} Dist; -t_{alt} (GeV^{2}); -t_{truth}(GeV^{2})", 50, 0, 0.5, 50, 0, 0.5);

  gDirectory->cd("../");
  h2_ZDC_XY = new TH2F("ZDC_XY", "ZDC XY", 200, -150, -50, 200, -50, 50);

  // Define beam 4 vectors
  e_beam_energy = 5;
  e_beam_pmag = sqrt(pow(e_beam_energy,2)-pow(mElec,2));
  ion_beam_energy = 100;
  ion_beam_pmag = sqrt((pow(ion_beam_energy,2)-pow(mProt,2)));
  crossing_angle = 0.025; 
  //Double_t Pi = TMath::ACos(-1);
  eBeam4Vect.SetPxPyPzE(0,0,-1*e_beam_pmag,e_beam_energy);
  pBeam4Vect.SetPxPyPzE(-ion_beam_pmag*TMath::Sin(crossing_angle),0,ion_beam_pmag*TMath::Cos(crossing_angle),ion_beam_energy);
  // ion_beam_pmag*TMath::Sin(crossing_angle)*TMath::Sin(Pi) was the y component, but this is 0 by a longer route...

  // Set cut values for physics analysis
  Thetan_Cent = 1.45; // Cut will be +/- 0.4 from this value
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP::InitRun(PHCompositeNode *topNode)
{
  std::cout << "ECCE_DEMP::InitRun(PHCompositeNode *topNode) Initializing for Run XXX" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP::process_event(PHCompositeNode *topNode)
{
  ZDC_hit = 0;
  EEMC_hit = 0;
  event_itt++; 

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
  t4VectTruth = virtphoton4VectTruth - pi4VectTruth;
  t_alt4VectTruth= pBeam4Vect - n4VectTruth;
  pmiss4VectTruth = (eBeam4Vect + pBeam4Vect) - (e4VectTruth+pi4VectTruth);
  pmiss4VectTruth_2 = (eBeam4Vect + pBeam4Vect) - (e4VectTruth+pi4VectTruth+n4VectTruth);
  Q2_truth = -1*(virtphoton4VectTruth.Mag2());
  W_truth = (virtphoton4VectTruth+pBeam4Vect).Mag();
  t_truth = -(t4VectTruth.Mag2());
  t_alt_truth = -(t_alt4VectTruth.Mag2());
  xb_truth =  Q2_truth/(2*(pBeam4Vect.Dot(virtphoton4VectTruth)));
  xi_truth = xb_truth/(2-xb_truth);

  h1_Q2_DetEff_Uncut->Fill(Q2_truth, wgt);
  h2_Q2_t_DetEff_Uncut->Fill(Q2_truth, t_truth, wgt);

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
	      << "ECCE_DEMP::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
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
	  nZDCPos.SetXYZ(hit_iter->second->get_x(0), hit_iter->second->get_y(0), hit_iter->second->get_z(0));
	  nEDep = hit_iter->second->get_edep();
	  nTheta = nZDCPos.Theta();
	  nPhi = nZDCPos.Phi();
	  nPMag = sqrt((pow(nEDep,2)) - (pow(mNeut,2)));
	  n4Vect.SetPxPyPzE(nPMag*sin(nTheta)*cos(nPhi), nPMag*sin(nTheta)*sin(nPhi), nPMag*cos(nTheta), nEDep);	  
	}
    }

    // Now have relevant information from this event, fill some histograms and calculate some stuff
    // Calculate kinematic quantities for electroproduction
    virtphoton4Vect = eBeam4Vect - e4Vect;
    t4Vect = virtphoton4Vect - pi4Vect;
    // Alternative calculation of t from different quantities, this seems to work much better
    t_alt4Vect= pBeam4Vect - n4Vect;
    pmiss4Vect = (eBeam4Vect + pBeam4Vect) - (e4Vect+pi4Vect);
    pmiss4Vect_2 = (eBeam4Vect + pBeam4Vect) - (e4Vect+pi4Vect+n4Vect);
    Q2 = -1*(virtphoton4Vect.Mag2());
    W = (virtphoton4Vect+pBeam4Vect).Mag();
    t = -(t4Vect.Mag2());
    t_alt = -(t_alt4Vect.Mag2());
    xb =  Q2/(2*(pBeam4Vect.Dot(virtphoton4Vect)));
    xi = xb/(2-xb);
    //cout << t_alt4Vect[0] << "  " << t_alt4Vect[1] << "  " << t_alt4Vect[2] << "  " << t_alt4Vect[3] << "   " << t_alt << endl;

    // Fill histograms
    // Fill weighted histograms

    h1_Q2_DetEff_Cut->Fill(Q2_truth, wgt);
    h2_Q2_t_DetEff_Cut->Fill(Q2_truth, t_truth, wgt);

    h1_pi_px->Fill(pi4Vect.Px(), wgt);
    h1_pi_py->Fill(pi4Vect.Py(), wgt);
    h1_pi_pz->Fill(pi4Vect.Pz(), wgt);
    h1_pi_p->Fill(pi4Vect.P(), wgt);
    h1_pi_E->Fill(pi4Vect.E(), wgt);
    h1_pi_Theta->Fill(pi4Vect.Theta()*TMath::RadToDeg(), wgt);
    h1_pi_Phi->Fill(pi4Vect.Phi()*TMath::RadToDeg(), wgt);
    h1_e_px->Fill(e4Vect.Px(), wgt);
    h1_e_py->Fill(e4Vect.Py(), wgt);
    h1_e_pz->Fill(e4Vect.Pz(), wgt);
    h1_e_p->Fill(e4Vect.P(), wgt);
    h1_e_E->Fill(e4Vect.E(), wgt);
    h1_e_Theta->Fill(e4Vect.Theta()*TMath::RadToDeg(), wgt);
    h1_e_Phi->Fill(e4Vect.Phi()*TMath::RadToDeg(), wgt);
    h1_n_px->Fill(n4Vect.Px(), wgt);
    h1_n_py->Fill(n4Vect.Py(), wgt);
    h1_n_pz->Fill(n4Vect.Pz(), wgt);
    h1_n_p->Fill(n4Vect.P(), wgt);
    h1_n_E->Fill(n4Vect.E(), wgt);
    h1_n_Theta->Fill(n4Vect.Theta()*TMath::RadToDeg(), wgt);
    h1_n_Phi->Fill(n4Vect.Phi()*TMath::RadToDeg(), wgt);
    h1_pmiss_px->Fill(pmiss4Vect.Px(), wgt);
    h1_pmiss_py->Fill(pmiss4Vect.Py(), wgt);
    h1_pmiss_pz->Fill(pmiss4Vect.Pz(), wgt);
    h1_pmiss_p->Fill(pmiss4Vect.P(), wgt);
    h1_pmiss_E->Fill(pmiss4Vect.E(), wgt);
    h1_pmiss_Theta->Fill(pmiss4Vect.Theta()*TMath::RadToDeg(), wgt);
    h1_pmiss_Phi->Fill(pmiss4Vect.Phi()*TMath::RadToDeg(), wgt);
    h1_gamma_px->Fill(virtphoton4Vect.Px(), wgt);
    h1_gamma_py->Fill(virtphoton4Vect.Py(), wgt);
    h1_gamma_pz->Fill(virtphoton4Vect.Pz(), wgt);
    h1_gamma_p->Fill(virtphoton4Vect.P(), wgt);
    h1_gamma_E->Fill(virtphoton4Vect.E(), wgt);
    h1_gamma_Theta->Fill(virtphoton4Vect.Theta()*TMath::RadToDeg(), wgt);
    h1_gamma_Phi->Fill(virtphoton4Vect.Phi()*TMath::RadToDeg(), wgt);

    h1_Q2_Dist->Fill(Q2, wgt);
    h1_W_Dist->Fill(W, wgt);
    h1_t_Dist->Fill(t, wgt);
    h1_t_alt_Dist->Fill(t_alt, wgt);
    h1_t_comp->Fill(((t_alt-t)/t)*100, wgt);
    h1_xb_Dist->Fill(xb, wgt);
    h1_xi_Dist->Fill(xi, wgt);

    h1_Q2Truth_Dist->Fill(Q2_truth, wgt);
    h1_WTruth_Dist->Fill(W_truth, wgt);
    h1_tTruth_Dist->Fill(t_truth, wgt);
    h1_t_altTruth_Dist->Fill(t_alt_truth, wgt);
    h1_xbTruth_Dist->Fill(xb_truth, wgt);
    h1_xiTruth_Dist->Fill(xi_truth, wgt);

    h2_t_ep->Fill(t, e4Vect.P(), wgt);
    h2_t_Q2->Fill(t,Q2, wgt);
    h2_delta_t_t->Fill(((t - t_truth)/t_truth)*100, t, wgt);
    
    for(Int_t B = 0; B < 7; B++){
      Q2_low = 5+(B*5);
      Q2_high = 10+(B*5);
      if ( Q2_truth > Q2_low && Q2_truth < Q2_high){
	h1_t_Q2[B]->Fill(t, wgt);
	h1_t_alt_Q2[B]->Fill(t_alt, wgt);
	h2_delta_t_t_Q2[B]->Fill(((t - t_truth)/t_truth)*100, t, wgt);
      }
    }

    for(Int_t B = 0; B < 8; B++){
      if ( B <= 1){
	Q2_low = 5+(B*2.5);
	Q2_high = 7.5+(B*2.5);
      }
      else{
	Q2_low = 5+((B-1)*5.0);
	Q2_high = 10+((B-1)*5.0);
      }
      if ( Q2 > Q2_low && Q2 < Q2_high){
	h1_t_result[B]->Fill(t_alt, wgt);
	h1_nTheta_result[B]->Fill(n4Vect.Theta()*TMath::RadToDeg(), wgt);
	h1_pmiss_result[B]->Fill(pmiss4Vect.P(), wgt);
	// Apply other cuts
	if ( (pmiss4Vect.P() < PmissCutVal[B]) && ((n4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (n4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5)))){
	h1_t_cut_result[B]->Fill(t_alt, wgt);
	h1_Q2_cut_result[B]->Fill(Q2, wgt);
	h1_W_cut_result[B]->Fill(W, wgt);
	h2_t_ttruth_result->Fill(t, t_truth, wgt);
	h2_t_alt_ttruth_result->Fill(t_alt, t_truth, wgt);
	h1_Mmiss_result->Fill(pmiss4Vect_2.M(),wgt);
	h1_Mmiss_truth_result->Fill(pmiss4VectTruth_2.M(),wgt);
	h1_Mmiss_Comp_result->Fill((pmiss4VectTruth_2.M()-pmiss4Vect_2.M()), wgt);
	h2_t_t_alt_result->Fill(t, t_alt, wgt);
	h2_Q2_W_result->Fill(Q2, W, wgt);
	h2_Q2_t_result[B]->Fill(Q2, t_alt, wgt);
	}
      }
    }

    h1_piTruth_p->Fill(pi4VectTruth.P(), wgt);
    h1_piTruth_px->Fill(pi4VectTruth.Px(), wgt);
    h1_piTruth_py->Fill(pi4VectTruth.Py(), wgt);
    h1_piTruth_pz->Fill(pi4VectTruth.Pz(), wgt);
    h1_piTruth_E->Fill(pi4VectTruth.E(), wgt);
    h1_piTruth_Theta->Fill(pi4VectTruth.Theta()*TMath::RadToDeg(), wgt);
    h1_eTruth_p->Fill(e4VectTruth.P(), wgt);
    h1_eTruth_px->Fill(e4VectTruth.Px(), wgt);
    h1_eTruth_py->Fill(e4VectTruth.Py(), wgt);
    h1_eTruth_pz->Fill(e4VectTruth.Pz(), wgt);
    h1_eTruth_E->Fill(e4VectTruth.E(), wgt);
    h1_eTruth_Theta->Fill(e4VectTruth.Theta()*TMath::RadToDeg(), wgt);
    h1_nTruth_p->Fill(n4VectTruth.P(), wgt);
    h1_nTruth_px->Fill(n4VectTruth.Px(), wgt);
    h1_nTruth_py->Fill(n4VectTruth.Py(), wgt);
    h1_nTruth_pz->Fill(n4VectTruth.Pz(), wgt);
    h1_nTruth_E->Fill(n4VectTruth.E(), wgt);
    h1_nTruth_Theta->Fill(n4VectTruth.Theta()*TMath::RadToDeg(), wgt);
    
    h1_piRes_p->Fill((pi4Vect.P()-pi4VectTruth.P())/(pi4VectTruth.P())*100, wgt);
    h1_piRes_px->Fill((pi4Vect.Px()-pi4VectTruth.Px())/(pi4VectTruth.Px())*100, wgt);
    h1_piRes_py->Fill((pi4Vect.Py()-pi4VectTruth.Py())/(pi4VectTruth.Py())*100, wgt);
    h1_piRes_pz->Fill((pi4Vect.Pz()-pi4VectTruth.Pz())/(pi4VectTruth.Pz())*100, wgt);
    h1_eRes_p->Fill((e4Vect.P()-e4VectTruth.P())/(e4VectTruth.P())*100, wgt);
    h1_eRes_px->Fill((e4Vect.Px()-e4VectTruth.Px())/(e4VectTruth.Px())*100, wgt);
    h1_eRes_py->Fill((e4Vect.Py()-e4VectTruth.Py())/(e4VectTruth.Py())*100, wgt);
    h1_eRes_pz->Fill((e4Vect.Pz()-e4VectTruth.Pz())/(e4VectTruth.Pz())*100, wgt);
    h1_nRes_p->Fill((n4Vect.P()-n4VectTruth.P())/(n4VectTruth.P())*100, wgt);
    h1_nRes_px->Fill((n4Vect.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100, wgt);
    h1_nRes_py->Fill((n4Vect.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100, wgt);
    h1_nRes_pz->Fill((n4Vect.Pz()-n4VectTruth.Pz())/(n4VectTruth.Pz())*100, wgt);

    h2_ZDC_XY->Fill(nZDCPos.x(), nZDCPos.y(), wgt);

    h2_piTrack_ThetaPhi->Fill((pi4Vect.Theta()*TMath::RadToDeg()), (pi4Vect.Phi()*TMath::RadToDeg()), wgt);
    h2_piTrack_pTheta->Fill((pi4Vect.Theta()*TMath::RadToDeg()), pi4Vect.P(), wgt);
    h2_eTrack_ThetaPhi->Fill((e4Vect.Theta()*TMath::RadToDeg()), (e4Vect.Phi()*TMath::RadToDeg()), wgt);
    h2_eTrack_pTheta->Fill((e4Vect.Theta()*TMath::RadToDeg()), e4Vect.P(), wgt);
    h2_nTrack_ThetaPhi->Fill((n4Vect.Theta()*TMath::RadToDeg()), (n4Vect.Phi()*TMath::RadToDeg()), wgt);
    h2_nTrack_pTheta->Fill((n4Vect.Theta()*TMath::RadToDeg()), n4Vect.P(), wgt);
    
    h2_piTruth_pxpy->Fill((pi4Vect.Px()-pi4VectTruth.Px())/(pi4VectTruth.Px())*100, (pi4Vect.Py()-pi4VectTruth.Py())/(pi4VectTruth.Py())*100, wgt);
    h2_eTruth_pxpy->Fill((e4Vect.Px()-e4VectTruth.Px())/(e4VectTruth.Px())*100, (e4Vect.Py()-e4VectTruth.Py())/(e4VectTruth.Py())*100, wgt);
    h2_nTruth_pxpy->Fill((n4Vect.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100, (n4Vect.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100, wgt);

    h2_piTruth_pxpz->Fill((pi4Vect.Px()-pi4VectTruth.Px())/(pi4VectTruth.Px())*100, (pi4Vect.Pz()-pi4VectTruth.Pz())/(pi4VectTruth.Pz())*100, wgt);
    h2_eTruth_pxpz->Fill((e4Vect.Px()-e4VectTruth.Px())/(e4VectTruth.Px())*100, (e4Vect.Pz()-e4VectTruth.Pz())/(e4VectTruth.Pz())*100, wgt);
    h2_nTruth_pxpz->Fill((n4Vect.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100, (n4Vect.Pz()-n4VectTruth.Pz())/(n4VectTruth.Pz())*100, wgt);

    h2_piTruth_pypz->Fill((pi4Vect.Py()-pi4VectTruth.Py())/(pi4VectTruth.Py())*100, (pi4Vect.Pz()-pi4VectTruth.Pz())/(pi4VectTruth.Pz())*100, wgt);
    h2_eTruth_pypz->Fill((e4Vect.Py()-e4VectTruth.Py())/(e4VectTruth.Py())*100, (e4Vect.Pz()-e4VectTruth.Pz())/(e4VectTruth.Pz())*100, wgt);
    h2_nTruth_pypz->Fill((n4Vect.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100, (n4Vect.Pz()-n4VectTruth.Pz())/(n4VectTruth.Pz())*100, wgt);

  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP::ResetEvent(PHCompositeNode *topNode)
{
//  std::cout << "ECCE_DEMP::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
//
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP::EndRun(const int runnumber)
{
  std::cout << "ECCE_DEMP::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP::End(PHCompositeNode *topNode)
{
  std::cout << "ECCE_DEMP::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  ScalingFact = double(event_itt)/nTried; // This scaling factor is needed to normalise the weighted results
  h1_Q2_DetEff->Divide(h1_Q2_DetEff_Cut, h1_Q2_DetEff_Uncut);
  h2_Q2_t_DetEff->Divide(h2_Q2_t_DetEff_Cut, h2_Q2_t_DetEff_Uncut);
  h2_Q2_W_result->Scale((1/ScalingFact));
  h2_ZDC_XY->Scale((1/ScalingFact));
  h2_t_ttruth_result->Scale((1/ScalingFact));
  h2_t_alt_ttruth_result->Scale((1/ScalingFact));
  h1_Mmiss_result->Scale((1/ScalingFact));
  h1_Mmiss_truth_result->Scale((1/ScalingFact));
  h1_Mmiss_Comp_result->Scale((1/ScalingFact));
  h2_t_t_alt_result->Scale((1/ScalingFact));

  h1_piRes_p->Scale((1/ScalingFact));
  h1_piRes_px->Scale((1/ScalingFact));
  h1_piRes_py->Scale((1/ScalingFact));
  h1_piRes_pz->Scale((1/ScalingFact));
  h1_eRes_p->Scale((1/ScalingFact));
  h1_eRes_px->Scale((1/ScalingFact));
  h1_eRes_py->Scale((1/ScalingFact));
  h1_eRes_pz->Scale((1/ScalingFact));
  h1_nRes_p->Scale((1/ScalingFact));
  h1_nRes_px->Scale((1/ScalingFact));
  h1_nRes_py->Scale((1/ScalingFact));
  h1_nRes_pz->Scale((1/ScalingFact));

  for(Int_t C = 0; C < 8; C++){
    h1_pmiss_result[C]->Scale((1/ScalingFact));
    h1_t_result[C]->Scale((1/ScalingFact));
    h1_t_cut_result[C]->Scale((1/ScalingFact));
    h1_Q2_cut_result[C]->Scale((1/ScalingFact));
    h1_W_cut_result[C]->Scale((1/ScalingFact));
    h2_Q2_t_result[C]->Scale((1/ScalingFact));
  }
  outfile->cd();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP::Reset(PHCompositeNode *topNode)
{
 std::cout << "ECCE_DEMP::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void ECCE_DEMP::Print(const std::string &what) const
{
  std::cout << "ECCE_DEMP::Print(const std::string &what) const Printing info for " << what << std::endl;
}

//***************************************************

float ECCE_DEMP::EMCAL_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.45*.45/E + 0.075*0.075);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}

//*****************************************************

float ECCE_DEMP::HCAL_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.50*.50/E + 0.1*0.1);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}

//*****************************************************

float ECCE_DEMP::PbWO4_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

//*****************************************************

float ECCE_DEMP::Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.1;         /// Position resolution 0.1 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}

//***************************************************

bool ECCE_DEMP::Check_ePi(PHCompositeNode* topNode)
{
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap)
    {
      trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
      if (!trackmap)
    	{
    	  cout
    	    << "ECCE_DEMP::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
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

bool ECCE_DEMP::Check_n(PHCompositeNode* topNode)
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
	if (hit_iter->second->get_edep() > 40){ // Hit in ZDC with roughly correct energy for neutron
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
