#include "ECCE_DEMP10on100.h"

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
#include <TH3.h>
#include <TProfile2D.h>
#include <TString.h>
#include <TTree.h>
#include <TVector3.h>
#include <TCanvas.h>

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

ECCE_DEMP10on100::ECCE_DEMP10on100(const std::string &name, const std::string& filename):
 SubsysReco(name)
 , outfilename(filename)
{
  std::cout << "ECCE_DEMP10on100_example::Diff_Tagg_example(const std::string &name) Calling ctor" << std::endl;

  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(m_RandomGenerator, seed);

}

//____________________________________________________________________________..
ECCE_DEMP10on100::~ECCE_DEMP10on100()
{

  gsl_rng_free(m_RandomGenerator);

  std::cout << "ECCE_DEMP10on100::~ECCE_DEMP10on100() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int ECCE_DEMP10on100::Init(PHCompositeNode *topNode)
{

  static_event_counter = 0;
  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(outfilename.c_str(), "RECREATE");

  std::cout << "ECCE_DEMP10on100::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  event_itt = 0;
  thrownEvents = 0; // mfek 06/22/2022 - new counter variable initialized
  cut1Events = 0; // mfek 06/22/2022 - new counter variable initialized 
  cut2Events = 0; // mfek 06/22/2022 - new counter variable initialized
  cut3Events = 0; // mfek 06/22/2022 - new counter variable initialized
  cut4Events = 0; // mfek 06/22/2022 - new counter variable initialized
  cut5Events = 0; // mfek 06/22/2022 - new counter variable initialized
  count_aftern = 0; // mfek 06/22/2022 - new counter variable initialized
  count_afterePi = 0; // mfek 06/22/2022 - new counter variable initialized
  
  gDirectory->mkdir("Detection_Efficiency");
  gDirectory->cd("Detection_Efficiency");
  h1_Q2_DetEff_Uncut = new TH1F("Q2_DetEff_Uncut", "Q^{2}_{Truth} for thrown events; Q^{2}", 100, 0, 40); h1_Q2_DetEff_Uncut->Sumw2();
  h1_Q2_DetEff_Cut = new TH1F("Q2_DetEff_Cut", "Q^{2}_{Truth} for detected events; Q^{2}", 100, 0, 40); h1_Q2_DetEff_Cut->Sumw2();
  h1_Q2_DetEff = new TH1F("Q2_DetEff", "Q^{2}_{Truth} detected/thrown ratio; Q^{2}", 100, 0, 40); h1_Q2_DetEff->Sumw2();
  h2_Q2_t_DetEff_Uncut = new TH2F("Q2_t_DetEff_Uncut", "Q^{2}_{Truth} vs -t_{Truth} for thrown events; Q^{2}; -t (GeV^{2}); Q^{2} (GeV/c^{2})", 10, 0, 40, 10, 0, 0.5); h2_Q2_t_DetEff_Uncut->Sumw2();
  h2_Q2_t_DetEff_Cut = new TH2F("Q2_t_DetEff_Cut", "Q^{2}_{Truth} vs -t_{Truth} for detected events; Q^{2}; -t (GeV^{2}); Q^{2} (GeV/c^{2})", 10, 0, 40, 10, 0, 0.5); h2_Q2_t_DetEff_Cut->Sumw2();
  h2_Q2_t_DetEff = new TH2F("Q2_t_DetEff", "Q^{2}_{Truth} vs -t_{Truth} detected/thrown ratio; Q^{2}; -t (GeV^{2}); Q^{2} (GeV/c^{2})", 10, 0, 40, 10, 0, 0.5); h2_Q2_t_DetEff->Sumw2();
  h2_Q2_t_DetEff_v2_Uncut = new TH2F("Q2_t_DetEff_v2_Uncut", "Q^{2}_{Truth} vs -t_{Truth} for thrown events; Q^{2}; -t (GeV^{2}); Q^{2} (GeV/c^{2})", 20, 0, 40, 20, 0, 0.5); h2_Q2_t_DetEff_v2_Uncut->Sumw2();
  h2_Q2_t_DetEff_v2_Cut = new TH2F("Q2_t_DetEff_v2_Cut", "Q^{2}_{Truth} vs -t_{Truth} for detected events; Q^{2}; -t (GeV^{2}); Q^{2} (GeV/c^{2})", 20, 0, 40, 20, 0, 0.5); h2_Q2_t_DetEff_v2_Cut->Sumw2();
  h2_Q2_t_DetEff_v2 = new TH2F("Q2_t_DetEff_v2", "Q^{2}_{Truth} vs -t_{Truth} detected/thrown ratio; Q^{2}; -t (GeV^{2}); Q^{2} (GeV/c^{2})", 20, 0, 40, 20, 0, 0.5); h2_Q2_t_DetEff_v2->Sumw2();
  gDirectory->cd("../");  

  gDirectory->mkdir("Particle_Momenta_Resolution");
  gDirectory->cd("Particle_Momenta_Resolution");
  h1_piRes_p = new TH1F("piRes_p", "#pi #frac{#Delta p}{Truth p} Distribution (%); #frac{#Delta p}{Truth p} (%)", 100, -50, 50);
  h1_piRes_px = new TH1F("piRes_p_{x}", "#pi #frac{#Delta p_{x}}{Truth p_{x}} Distribution (%); #frac{#Delta p_{x}}{Truth p_{x}} (%)", 100, -50, 50);
  h1_piRes_py = new TH1F("piRes_p_{y}", "#pi #frac{#Delta p_{y}}{Truth p_{y}} Distribution (%); #frac{#Delta p_{y}}{Truth p_{y}} (%)", 100, -50, 50);
  h1_piRes_pz = new TH1F("piRes_p_{z}", "#pi #frac{#Delta p_{z}}{Truth p_{z}} Distribution (%); #frac{#Delta p_{z}}{Truth p_{z}} (%)", 100, -50, 50);
  h1_eRes_p = new TH1F("eRes_p", "e' #frac{#Delta p}{Truth p} Distribution (%); #frac{#Delta p}{Truth p} (%)", 100, -50, 50);
  h1_eRes_px = new TH1F("eRes_p_{x}", "e' #frac{#Delta p_{x}}{Truth p_{x}} Distribution (%); #frac{#Delta p_{x}}{Truth p_{x}} (%)", 100, -50, 50);
  h1_eRes_py = new TH1F("eRes_p_{y}", "e' #frac{#Delta p_{y}}{Truth p_{y}} Distribution (%); #frac{#Delta p_{y}}{Truth p_{y}} (%)", 100, -50, 50);
  h1_eRes_pz = new TH1F("eRes_p_{z}", "e' #frac{#Delta p_{z}}{Truth p_{z}} Distribution (%); #frac{#Delta p_{z}}{Truth p_{z}} (%)", 100, -50, 50);
  h1_nRes_p = new TH1F("nRes_p", "Neutron #frac{#Delta p_{n}}{p_{ntruth}} Distribution (%); #frac{#Delta p_{n}}{p_{ntruth}} (%)", 100, -5, 5);
  h1_nRes_px = new TH1F("nRes_p_{x}", "n #frac{#Delta p_{x}}{Truth p_{x}} Distribution (%); #frac{#Delta p_{x}}{Truth p_{x}} (%)", 100, -20, 20);
  h1_nRes_py = new TH1F("nRes_p_{y}", "n #frac{#Delta p_{y}}{Truth p_{y}} Distribution (%); #frac{#Delta p_{y}}{Truth p_{y}} (%)", 100, -20, 20);
  h1_nRes_pz = new TH1F("nRes_p_{z}", "n #frac{#Delta p_{z}}{Truth p_{z}} Distribution (%); #frac{#Delta p_{z}}{Truth p_{z}} (%)", 100, -5, 5);
  h1_pmissDiff_p = new TH1F("pmissDiff_p", "#Delta p_{miss} (p_{miss} - p_{nTruth}) Dist; #Delta p_{miss} (GeV/c)", 100, -5, 5);
  h1_pmissDiff_px = new TH1F("pmissDiff_p_{x}", "#Delta p_{xmiss} (p_{xmiss} - p_{xnTruth}) Dist; #Delta p_{xmiss} (GeV/c)", 100, -5, 5);
  h1_pmissDiff_py = new TH1F("pmissDiff_p_{y}", "#Delta p_{ymiss} (p_{ymiss} - p_{ynTruth}) Dist; #Delta p_{ymiss} (GeV/c)", 100, -5, 5);
  h1_pmissDiff_pz = new TH1F("pmissDiff_p_{z}", "#Delta p_{zmiss} (p_{zmiss} - p_{znTruth}) Dist; #Delta p_{ymiss} (GeV/c)", 100, -5, 5);
  gDirectory->cd("../");  

  gDirectory->mkdir("Pion_Info");
  gDirectory->cd("Pion_Info");
  h1_pi_px = new TH1F("pi_px", "#pi p_{x} Distribution;p_{x} (GeV)", 200, -20, 20);
  h1_pi_py = new TH1F("pi_py", "#pi p_{y} Distribution;p_{y} (GeV)", 200, -20, 20);
  h1_pi_pz = new TH1F("pi_pz", "#pi p_{z} Distribution;p_{z} (GeV)", 200, -60, 60); // mfek 06/22/2022 - adjusted range 
  h1_pi_p = new TH1F("pi_p", "#pi p Distribution;p (GeV)", 200, 0, 60); // mfek 06/22/2022 - adjusted range 
  h1_pi_E = new TH1F("pi_E", "#pi E Distribution;E (GeV)", 200, 0, 60); // mfek 06/22/2022 - adjusted range 
  h1_pi_Theta = new TH1F("pi_Theta", "#pi #theta Distribution; #theta (deg)", 200, 0, 50);
  h1_pi_Phi = new TH1F("pi_Phi", "#pi #phi Distribution; #phi (deg)", 360, -180, 180);
  h2_piTrack_ThetaPhi = new TH2F("piTrack_ThetaPhi", "#pi Track #theta vs #phi; #theta (deg); #phi (deg)", 120, 0, 60, 720, -180, 180);
  h2_piTrack_pTheta = new TH2F("piTrack_pTheta", "#pi Track #theta vs P; #theta (deg); P (GeV/c)", 120, 0, 60, 500, 0, 60); // mfek 06/22/2022 - adjusted range 
  h2_pi_XY = new TH2F("pi_XY", "#pi X vs Y at z=375cm (HCal) Dist; x(cm); y(cm)", 120, -300, 300, 120, -300, 300);
  gDirectory->cd("../");

  gDirectory->mkdir("Pion_Unweighted_Info");
  gDirectory->cd("Pion_Unweighted_Info");
  h1_pi_px_Unweighted = new TH1F("pi_px_Unweighted", "#pi Unweighted p_{x} Distribution;p_{x} (GeV)", 200, -20, 20);
  h1_pi_py_Unweighted = new TH1F("pi_py_Unweighted", "#pi Unweighted p_{y} Distribution;p_{y} (GeV)", 200, -20, 20);
  h1_pi_pz_Unweighted = new TH1F("pi_pz_Unweighted", "#pi Unweighted p_{z} Distribution;p_{z} (GeV)", 200, -60, 60); // mfek 06/22/2022 - adjusted range 
  h1_pi_p_Unweighted = new TH1F("pi_p_Unweighted", "#pi Unweighted p Distribution;p (GeV)", 200, 0, 60); // mfek 06/22/2022 - adjusted range 
  h1_pi_E_Unweighted = new TH1F("pi_E_Unweighted", "#pi Unweighted E Distribution;E (GeV)", 200, 0, 60); // mfek 06/22/2022 - adjusted range 
  h1_pi_Theta_Unweighted = new TH1F("pi_Theta_Unweighted", "#pi Unweighted #theta Distribution; #theta (deg)", 200, 0, 50);
  h1_pi_Phi_Unweighted = new TH1F("pi_Phi_Unweighted", "#pi Unweighted #phi Distribution; #phi (deg)", 360, -180, 180);
  gDirectory->cd("../");

  gDirectory->mkdir("Pion_Truth_Info");
  gDirectory->cd("Pion_Truth_Info");
  h1_piTruth_px = new TH1F("piTrtuh_px", "#pi Truth p_{x} Distribution;p_{x} (GeV)", 200, -20, 20);
  h1_piTruth_py = new TH1F("piTrtuh_py", "#pi Truth p_{y} Distribution;p_{y} (GeV)", 200, -20, 20);
  h1_piTruth_pz = new TH1F("piTrtuh_pz", "#pi Truth p_{z} Distribution;p_{z} (GeV)", 200, -60, 60); // mfek 06/22/2022 - adjusted range  
  h1_piTruth_p = new TH1F("piTrtuh_p", "#pi Truth p Distribution;p (GeV)", 200, 0, 60); // mfek 06/22/2022 - adjusted range 
  h1_piTruth_E = new TH1F("piTrtuh_E", "#pi Truth E Distribution;E (GeV)", 200, 0, 60); // mfek 06/22/2022 - adjusted range 
  h1_piTruth_Theta = new TH1F("piTrtuh_Theta", "#pi Truth #theta Distribution; #theta (deg)", 200, 0, 50);

  h2_piTruth_pxpy = new TH2F("piTruth_pxpy", "#pi #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  h2_piTruth_pxpz = new TH2F("piTruth_pxpz", "#pi #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_piTruth_pypz = new TH2F("piTruth_pypz", "#pi #frac{#Delta p_{y}}{Truth p_{y}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_piTrack_pTheta_Truth = new TH2F("piTrack_pTheta_Truth", "#pi Track #theta vs P (Truth); #theta (deg); P (GeV/c)", 120, 0, 60, 500, 0, 60); // mfek 06/22/2022 - adjusted range 
  gDirectory->cd("../");
  
  gDirectory->mkdir("Scattered_Electron_Info");
  gDirectory->cd("Scattered_Electron_Info");
  h1_e_px = new TH1F("e_px", "e' p_{x} Distribution;p_{x} (GeV)", 240, -6, 6);
  h1_e_py = new TH1F("e_py", "e' p_{y} Distribution;p_{y} (GeV)", 240, -6, 6);
  h1_e_pz = new TH1F("e_pz", "e' p_{z} Distribution;p_{z} (GeV)", 120, -15, 0); // mfek 06/22/2022 - adjusted range
  h1_e_p = new TH1F("e_p", "e' p Distribution;p (GeV)", 160, 0, 15); // mfek 06/22/2022 - adjusted range 
  h1_e_E = new TH1F("e_E", "e' E Distribution;E (GeV)", 160, 0, 15); // mfek 06/22/2022 - adjusted range 
  h1_e_Theta = new TH1F("e_Theta", "e' #theta Distribution; #theta (deg)", 200, 110, 180); // mfek 06/22/2022 - adjusted range 
  h1_e_Phi = new TH1F("e_Phi", "e' #phi Distribution; #phi (deg)", 360, -180, 180);
  h2_eTrack_ThetaPhi = new TH2F("eTrack_ThetaPhi", "e' Track #theta vs #phi; #theta (deg); #phi (deg)", 140, 110, 180, 720, -180, 180);
  h2_eTrack_pTheta = new TH2F("eTrack_pTheta", "e' Track #theta vs P; #theta (deg); P (GeV/c)", 140, 110, 180, 100, 0, 15); // mfek 06/22/2022 - adjusted range 
  h2_e_XY = new TH2F("e_XY", "e' X vs Y at z=200cm (EMCal) Dist; x (cm); y (cm)", 240, -600, 600, 240, -600, 600);
  gDirectory->cd("../");

  gDirectory->mkdir("Scattered_Electron_Unweighted_Info");
  gDirectory->cd("Scattered_Electron_Unweighted_Info");
  h1_e_px_Unweighted = new TH1F("e_px_Unweighted", "e' Unweighted p_{x} Distribution;p_{x} (GeV)", 240, -6, 6);
  h1_e_py_Unweighted = new TH1F("e_py_Unweighted", "e' Unweighted p_{y} Distribution;p_{y} (GeV)", 240, -6, 6);
  h1_e_pz_Unweighted = new TH1F("e_pz_Unweighted", "e' Unweighted p_{z} Distribution;p_{z} (GeV)", 120, -15, 0); // mfek 06/22/2022 - adjusted range  
  h1_e_p_Unweighted = new TH1F("e_p_Unweighted", "e' Unweighted p Distribution;p (GeV)", 160, 0, 15); // mfek 06/22/2022 - adjusted range 
  h1_e_E_Unweighted = new TH1F("e_E_Unweighted", "e' Unweighted E Distribution;E (GeV)", 160, 0, 15); // mfek 06/22/2022 - adjusted range 
  h1_e_Theta_Unweighted = new TH1F("e_Theta_Unweighted", "e' Unweighted #theta Distribution; #theta (deg)", 200, 110, 180); // mfek 06/22/2022 - adjusted range 
  h1_e_Phi_Unweighted = new TH1F("e_Phi_Unweighted", "e' Unweighted #phi Distribution; #phi (deg)", 360, -180, 180);
  gDirectory->cd("../");

  gDirectory->mkdir("Scattered_Electron_Truth_Info");
  gDirectory->cd("Scattered_Electron_Truth_Info");
  h1_eTruth_px = new TH1F("eTruth_px", "e' Truth p_{x} Distribution;p_{x} (GeV)", 240, -6, 6);
  h1_eTruth_py = new TH1F("eTruth_py", "e' Truth p_{y} Distribution;p_{y} (GeV)", 240, -6, 6);
  h1_eTruth_pz = new TH1F("eTruth_pz", "e' Truth p_{z} Distribution;p_{z} (GeV)", 120, -15, 0); // mfek 06/22/2022 - adjusted range  
  h1_eTruth_p = new TH1F("eTruth_p", "e' Truth p Distribution;p (GeV)", 160, 0, 15); // mfek 06/22/2022 - adjusted range 
  h1_eTruth_E = new TH1F("eTruth_E", "e' Truth E Distribution;E (GeV)", 160, 0, 15); // mfek 06/22/2022 - adjusted range 
  h1_eTruth_Theta = new TH1F("eTruth_Theta", "e' Truth #theta Distribution; #theta (deg)", 200, 110, 180); // mfek 06/22/2022 - adjusted range 
  h2_eTruth_pxpy = new TH2F("eTruth_pxpy", "e' #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);  
  h2_eTruth_pxpz = new TH2F("eTruth_pxpz", "e' #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);  
  h2_eTruth_pypz = new TH2F("eTruth_pypz", "e' #frac{#Delta p_{y}}{Truth p_{y}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_eTrack_pTheta_Truth = new TH2F("eTrack_pTheta_Truth", "e' Track #theta vs P (Truth); #theta (deg); P (GeV/c)", 140, 110, 180, 100, 0, 15); // mfek 06/22/2022 - adjusted range 
  gDirectory->cd("../");

  gDirectory->mkdir("ZDC_Neutron_Info");
  gDirectory->cd("ZDC_Neutron_Info");
  h1_n_px = new TH1F("n_px", "n p_{x} Distribution;p_{x} (GeV)", 480, -6, 6);
  h1_n_py = new TH1F("n_py", "n p_{y} Distribution;p_{y} (GeV)", 200, -2.5, 2.5);
  h1_n_pz = new TH1F("n_pz", "n p_{z} Distribution;p_{z} (GeV)", 240, 0, 120); 
  h1_n_p = new TH1F("n_p", "n p Distribution;p (GeV)", 240, 0, 120);
  h1_n_E = new TH1F("n_E", "n E Distribution;E (GeV)", 240, 0, 120);
  h1_n_Theta = new TH1F("n_Theta", "n #theta Distribution; #theta (deg)", 500, 0, 5);
  h1_n_Phi = new TH1F("n_Phi", "n #phi Distribution; #phi (deg)", 360, -180, 180);

  h2_nTrack_ThetaPhi = new TH2F("nTrack_ThetaPhi", "n Track #theta vs #phi; #theta (deg); #phi (deg)", 500, 0, 5, 360, -180, 180);
  h2_nTrack_pTheta = new TH2F("nTrack_pTheta", "n Track #theta vs P; #theta (deg); P (GeV/c)", 500, 0, 5, 1000, 0, 120); // mfek 06/22/2022 - adjusted range 
  h2_n_XY = new TH2F("n_XY", "n X vs Y at ZDC Dist; x (cm); y (cm)", 200, -150, -50, 200, -50, 50);
  h1_n_ThetaDiff = new TH1F("n_ThetaDiff", "#theta_{pMiss} - #theta_{ZDC}; #theta_{pMiss}-#theta_{ZDC}(Deg)", 100, -5, 5);
  h1_n_PhiDiff = new TH1F("n_PhiDiff", " #phi_{pMiss} - #phi_{ZDC}; #phi_{pMiss}-#phi_{ZDC}(Deg)", 200, -25, 25);
  h2_n_ThetaPhiDiff = new TH2F("n_ThetaPhiDiff", "#theta_{pMiss} - #theta_{ZDC} vs #phi_{pMiss} - #phi_{ZDC}; #theta_{pMiss}-#theta_{ZDC}; #phi_{pMiss}-#phi_{ZDC}(Deg)",100, -5, 5, 200, -25, 25);
  gDirectory->cd("../");

  gDirectory->mkdir("ZDC_Neutron_Unweighted_Info");
  gDirectory->cd("ZDC_Neutron_Unweighted_Info");
  h1_n_px_Unweighted = new TH1F("n_px_Unweighted", "n p_{x} Distribution (Unweighted);p_{x} (GeV)", 480, -6, 6);
  h1_n_py_Unweighted = new TH1F("n_py_Unweighted", "n p_{y} Distribution (Unweighted);p_{y} (GeV)", 200, -2.5, 2.5);
  h1_n_pz_Unweighted = new TH1F("n_pz_Unweighted", "n p_{z} Distribution (Unweighted);p_{z} (GeV)", 240, 0, 120); 
  h1_n_p_Unweighted = new TH1F("n_p_Unweighted", "n p Distribution (Unweighted);p (GeV)", 240, 0, 120);
  h1_n_E_Unweighted = new TH1F("n_E_Unweighted", "n E Distribution (Unweighted);E (GeV)", 240, 0, 120);
  h1_n_Theta_Unweighted = new TH1F("n_Theta_Unweighted", "n #theta Distribution (Unweighted); #theta (deg)", 500, 0, 5);
  h1_n_Phi_Unweighted = new TH1F("n_Phi_Unweighted", "n #phi Distribution (Unweighted); #phi (deg)", 360, -180, 180);
  h1_n_ThetaDiff_Unweighted = new TH1F("n_ThetaDiff_Unweighted", "#theta_{pMiss} - #theta_{ZDC} (Unweighted); #theta_{pMiss}-#theta_{ZDC}(Deg)", 100, -5, 5);
  h1_n_PhiDiff_Unweighted = new TH1F("n_PhiDiff_Unweighted", "#phi_{pMiss} - #phi_{ZDC} (Unweighted); #phi_{pMiss}-#phi_{ZDC}(Deg)", 200, -25, 25);
  h2_n_ThetaPhiDiff_Unweighted = new TH2F("n_ThetaPhiDiff_Unweighted", "#theta_{pMiss} - #theta_{ZDC} vs #phi_{pMiss} - #phi_{ZDC} (Unweighted); #theta_{pMiss}-#theta_{ZDC}; #phi_{pMiss}-#phi_{ZDC}(Deg)",100, -5, 5, 200, -25, 25);
  gDirectory->cd("../");

  gDirectory->mkdir("Reconstructed_Neutron_Info");
  gDirectory->cd("Reconstructed_Neutron_Info");
  h1_nRec_px = new TH1F("nRec_px", "n_{Rec}p_{x} Distribution;p_{x} (GeV)", 480, -6, 6);
  h1_nRec_py = new TH1F("nRec_py", "n_{Rec}p_{y} Distribution;p_{y} (GeV)", 200, -2.5, 2.5);
  h1_nRec_pz = new TH1F("nRec_pz", "n_{Rec}p_{z} Distribution;p_{z} (GeV)", 240, 0, 120); 
  h1_nRec_p = new TH1F("nRec_p", "n_{Rec}p Distribution;p (GeV)", 240, 0, 120);
  h1_nRec_E = new TH1F("nRec_E", "n_{Rec}E Distribution;E (GeV)", 240, 0, 120);
  h1_nRec_Theta = new TH1F("nRec_Theta", "n_{Rec}#theta Distribution; #theta (deg)", 500, 0, 5);
  h1_nRec_Phi = new TH1F("nRec_Phi", "n_{Rec}#phi Distribution; #phi (deg)", 360, -180, 180);
  gDirectory->cd("../");

  gDirectory->mkdir("Neutron_Truth_Info");
  gDirectory->cd("Neutron_Truth_Info");
  h1_nTruth_px = new TH1F("nTruth_px", "nTruth p_{x} Distribution;p_{x} (GeV)", 480, -6, 6);
  h1_nTruth_py = new TH1F("nTruth_py", "nTruth p_{y} Distribution;p_{y} (GeV)", 200, -2.5, 2.5);
  h1_nTruth_pz = new TH1F("nTruth_pz", "nTruth p_{z} Distribution;p_{z} (GeV)", 240, 0, 120); 
  h1_nTruth_p = new TH1F("nTruth_p", "nTruth p Distribution;p (GeV)", 240, 0, 120);
  h1_nTruth_E = new TH1F("nTruth_E", "nTruth E Distribution;E (GeV)", 240, 0, 120);
  h1_nTruth_Theta = new TH1F("nTruth_Theta", "nTruth #theta Distribution; #theta (deg)", 500, 0, 5);
  h2_nTruth_pxpy = new TH2F("nTruth_pxpy", "n #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{y}}{Truth p_{y}}", 100, -50, 50, 100, -50, 50);
  h2_nTruth_pxpz = new TH2F("nTruth_pxpz", "n #frac{#Delta p_{x}}{Truth p_{x}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{x}}{Truth p_{x}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_nTruth_pypz = new TH2F("nTruth_pypz", "n #frac{#Delta p_{y}}{Truth p_{y}} vs #frac{#Delta p_{z}}{Truth p_{z}}; #frac{#Delta p_{y}}{Truth p_{y}}; #frac{#Delta p_{z}}{Truth p_{z}}", 100, -50, 50, 100, -50, 50);
  h2_nTrack_pTheta_Truth = new TH2F("nTrack_pTheta_Truth", "n Track #theta vs P (Truth); #theta (deg); P (GeV/c)", 500, 0, 5, 1000, 0, 120); // mfek 06/22/2022 - adjusted range 
  gDirectory->cd("../");

  gDirectory->mkdir("PMiss_Info");
  gDirectory->cd("PMiss_Info");
  h1_pmiss_px = new TH1F("pmiss_px", "p_{miss} p_{x} Distribution", 800, -10, 10);
  h1_pmiss_py = new TH1F("pmiss_py", "p_{miss} p_{y} Distribution", 200, -2.5, 2.5);
  h1_pmiss_pz = new TH1F("pmiss_pz", "p_{miss} p_{z} Distribution", 240, 0, 120); 
  h1_pmiss_p = new TH1F("pmiss_p", "p_{miss} p Distribution", 240, 0, 120);
  h1_pmiss_E = new TH1F("pmiss_E", "p_{miss} E Distribution", 240, 0, 120);
  h1_pmiss_Theta = new TH1F("pmiss_Theta", "p_{miss} #theta Distribution; #theta (deg)", 1000, 0, 10);
  h1_pmiss_Phi = new TH1F("pmiss_Phi", "p_{miss} #phi Distribution; #phi (deg)", 720, -180, 180);
  gDirectory->cd("../");
  
  gDirectory->mkdir("Virtual_Photon_Info");
  gDirectory->cd("Virtual_Photon_Info");
  h1_gamma_px = new TH1F("gamma_px", "#gamma p_{x} Distribution", 200, -10, 10);
  h1_gamma_py = new TH1F("gamma_py", "#gamma p_{y} Distribution", 200, -10, 10);
  h1_gamma_pz = new TH1F("gamma_pz", "#gamma p_{z} Distribution", 200, -10, 10); // mfek 06/22/2022 - adjusted range CHECK 
  h1_gamma_p = new TH1F("gamma_p", "#gamma p Distribution", 200, 0, 10);
  h1_gamma_E = new TH1F("gamma_E", "#gamma E Distribution", 200, 0, 10);
  h1_gamma_Theta = new TH1F("gamma_Theta", "#gamma #theta Distribution; #theta (deg)", 360, -180, 180);
  h1_gamma_Phi = new TH1F("gamma_Phi", "#gamma #phi Distribution; #phi (deg)", 360, -180, 180);
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
  h2_t_ep = new TH2F("t_ep", "t vs ScatElec P; t; P_{e'}", 100, 0, 10, 200, 0, 15); // mfek 06/22/2022 - adjusted range 
  h2_t_Q2 = new TH2F("t_Q2", "t vs Q^{2}; t; Q^{2}", 100, 0, 10, 200, 0, 50);
  h2_delta_t_t = new TH2F("delta_t_t", "#Delta t vs t; #Delta t (%); t", 200, -100, 100, 100, 0, 1);

  for(Int_t A = 0; A < 14; A++){ 
    h1_t_Q2[A] = new TH1F(Form("t_Q2_%i", (A+1)), Form("t dist, %2.1f < Q^{2} < %2.1f; t", Q2BinVal[A], Q2BinVal[A+1]), 100, 0, 10);
    h1_t_alt_Q2[A] = new TH1F(Form("t_alt_Q2_%i", (A+1)), Form("t (Alternative calculation) dist, %2.1f < Q^{2} < %2.1f; t", Q2BinVal[A], Q2BinVal[A+1]), 100, 0, 10);
    h2_delta_t_t_Q2[A] = new TH2F(Form("delta_t_t_Q2_%i", (A+1)), Form("#Delta t vs t, %2.1f < Q^{2} < %2.1f; #Delta t (Percent); t", Q2BinVal[A], Q2BinVal[A+1]), 200, -100, 100, 100, 0, 1);
  } // mfek 06/23/2022 - changed Q2 binning
  gDirectory->cd("../");

  gDirectory->mkdir("Physics_Results_Cuts");
  gDirectory->cd("Physics_Results_Cuts");

  for(Int_t A = 0; A < 14; A++){
    h1_t_result[A] = new TH1F(Form("t_Result_Q2_%i", (A+1)), Form("-t Dist, %2.1f < Q^{2} < %2.1f; -t", Q2BinVal[A], Q2BinVal[A+1]), 10, 0, 0.4);
    h1_nTheta_result[A] = new TH1F(Form("nTheta_Result_Q2_%i", (A+1)), Form("#theta_{n} Dist, %2.1f < Q^{2} < %2.1f", Q2BinVal[A], Q2BinVal[A+1]), 500, 0, 5);
    h1_pmiss_result[A] = new TH1F(Form("pmiss_Result_Q2_%i", (A+1)), Form("p_{miss} Dist, %2.1f < Q^{2} < %2.1f", Q2BinVal[A], Q2BinVal[A+1]), 240, 0, 50) ;
    h1_pn_result[A] = new TH1F(Form("pn_Result_Q2_%i", (A+1)), Form("p_{n} Dist, %2.1f < Q^{2} < %2.1f", Q2BinVal[A], Q2BinVal[A+1]), 240, 0, 120) ;
  }// mfek 06/22/2022 - changed binning
  // Results histograms are binned in Q2

  gDirectory->cd("../");
  gDirectory->mkdir("Physics_Results");
  gDirectory->cd("Physics_Results");
  for(Int_t A = 0; A < 14; A++){
    h1_t_cut_result[A] = new TH1F(Form("t_cut_Result_Q2_%i", (A+1)), Form("-t Dist, %2.1f < Q^{2} < %2.1f, with p_{miss}, #theta_{n} cuts; -t (GeV^{2}); Rate(Hz)", Q2BinVal[A], Q2BinVal[A+1]), 10, 0, 0.4);
    h1_t_truth_thrown_result[A] = new TH1F(Form("t_truth_thrown_Result_Q2_%i", (A+1)), Form("-t Dist, %2.1f < Q^{2} < %2.1f; -t", Q2BinVal[A], Q2BinVal[A+1]), 10, 0, 0.4);
    h1_Q2_cut_result[A] = new TH1F(Form("Q2_cut_Result_Q2_%i", (A+1)), Form("Q^{2} Dist, %2.1f < Q^{2} < %2.1f, with p_{miss}, #theta_{n} cuts; Q^{2}(GeV^{2}); Rate(Hz)", Q2BinVal[A], Q2BinVal[A+1]), 25, Q2BinVal[A], Q2BinVal[A+1]);
    h1_W_cut_result[A] = new TH1F(Form("W_cut_Result_Q2_%i", (A+1)), Form("W Dist, %2.1f < Q^{2} < %2.1f, with p_{miss}, #theta_{n} cuts; W(GeV); Rate(Hz)", Q2BinVal[A], Q2BinVal[A+1]), 60, 0, 30);
    h2_Q2_t_result[A] = new TH2F(Form("Q2_t_cut_Result_Q2_%i", (A+1)), Form("Q^{2} vs -t Dist, %2.1f < Q^{2} < %2.1f, with p_{miss}, #theta_{n} cuts; Q^{2}(GeV^{2}); -t (GeV^{2})", Q2BinVal[A], Q2BinVal[A+1]), 25, Q2BinVal[A], Q2BinVal[A+1], 10, 0, 0.4);
    h2_t_ttruth_result_Q2[A] = new TH2F(Form("t_ttruth_Result_Q2_%i",(A+1)), Form("-t vs -t_{truth}, %2.1f < Q^{2} < %2.1f, with p_{miss}, #theta_{n} cuts; -t (GeV^{2}); -t_{truth}(GeV^{2})",Q2BinVal[A], Q2BinVal[A+1]), 50, 0, 10, 50, 0, 0.5); 
    h2_t_alt_ttruth_result_Q2[A] = new TH2F(Form("t_alt_ttruth_Result_Q2_%i",(A+1)), Form("-t_{alt} vs -t_{truth}, %2.1f < Q^{2} < %2.1f, with p_{miss}, #theta_{n} cuts; -t_{alt} (GeV^{2}); -t_{truth}(GeV^{2})", Q2BinVal[A], Q2BinVal[A+1]), 50, 0, 0.5, 50, 0, 0.5);
  } // mfek 06/22/2022 - changed binning

  gDirectory->cd("../");
  gDirectory->mkdir("Physics_Results_Misc");
  gDirectory->cd("Physics_Results_Misc");
  h2_Q2_W_result = new TH2F("Q2_W_Result", "Q^{2} vs W, with p_{miss}, #theta_{n} cuts; Q^{2} (GeV^{2}); W (GeV)", 200, 0, 40, 60, 0, 30);
  h2_t_t_alt_result = new TH2F("t_t_alt_result", "-t vs -t_{alt}, with p_{miss}, #theta_{n} cuts; -t (GeV^{2}); -t_{alt} (GeV^{2})", 50, 0, 10, 12, 0, 0.48); 
  h1_Mmiss_result = new TH1F("Mmiss_result", "Missing Mass Dist;M_{Miss}(GeV/c^{2})", 100, -5, 5);
  h1_MmissSq_result = new TH1F("MmissSq_result", "Missing Mass Squared Dist;M_{Miss}^{2} ((GeV/c^{2})^{2})", 80, 0, 4);
  h1_Mmiss_truth_result = new TH1F("Mmiss_truth_result", "Missing Mass (truth) Dist; M_{Miss} (GeV/c^{2})", 40, -2, 2);
  h1_Mmiss_Comp_result = new TH1F("Mmiss_Comp_result", "#Delta M_{Miss} Distribution; #Delta(M_{Miss}) (GeV/c^{2})", 100, -5, 5);
  h1_taltres_result = new TH1F("taltres_result", "t_{alt} - t_{truth} Dist; t_{alt} - t_{truth} (GeV^{2})", 100, -0.5, 0.5);
  for(Int_t A = 0; A < 10; A++){
    h1_taltres_result_ttruth[A] = new TH1F(Form("taltres_result_ttruth_%i", (A+1)), Form("t_{alt} - t_{truth} Dist, %2.2f < t_{truth} < %2.2f; t_{alt} - t_{truth} (GeV^{2})", (0+(A*0.05)), (0.05+(A*0.05))), 100, -0.3, 0.8);
  }
  h2_t_ttruth_result = new TH2F("t_ttruth_result", "-t vs -t_{truth} Dist; -t (GeV^{2}); -t_{truth} (GeV^{2})", 50, 0, 10, 50, 0, 0.5);
  h2_t_alt_ttruth_result = new TH2F("t_alt_ttruth_result", "-t_{alt} vs -t_{truth} Dist; -t_{alt} (GeV^{2}); -t_{truth} (GeV^{2})", 50, 0, 0.5, 50, 0, 0.5);

  gDirectory->cd("../");
  gDirectory->mkdir("Cut_Analysis");
  gDirectory->cd("Cut_Analysis");
  h1_nTheta_tCut = new TH1F("nTheta_tCut", "n #theta Distribution (-t < 0.4); #theta (deg)", 500, 0, 5);
  h1_t_cut1_Low = new TH1F("t_result_cut1_Low", "t dist 7.5 < Q^{2} < 10, t cut only; -t(GeV^{2})", 10, 0, 0.4);
h1_t_cut2_Low = new TH1F("t_result_cut2_Low", "t dist 7.5 < Q^{2} < 10, t cut, #theta_{n} cut; -t(GeV^{2})", 10, 0, 0.4);
 h1_t_cut3_Low = new TH1F("t_result_cut3_Low", "t dist 7.5 < Q^{2} < 10, t cut, #theta_{n} cut, #theta_{diff} , phi_{diff} cuts; -t(GeV^{2})", 10, 0, 0.4);
h1_t_cut4_Low = new TH1F("t_result_cut4_Low", "t dist 7.5 < Q^{2} < 10,  t cut, #theta_{n} cut, #theta_{diff} , phi_{diff}, p_{miss} cuts; -t(GeV^{2})", 10, 0, 0.4);
  h1_t_cut1_Mid = new TH1F("t_result_cut1_Mid", "t dist 15 < Q^{2} < 20, t cut only; -t(GeV^{2})", 10, 0, 0.4);
h1_t_cut2_Mid = new TH1F("t_result_cut2_Mid", "t dist 15 < Q^{2} < 20, t cut, #theta_{n} cut; -t(GeV^{2})", 10, 0, 0.4);
 h1_t_cut3_Mid = new TH1F("t_result_cut3_Mid", "t dist 15 < Q^{2} < 20, t cut, #theta_{n} cut, #theta_{diff} , phi_{diff} cuts; -t(GeV^{2})", 10, 0, 0.4);
h1_t_cut4_Mid = new TH1F("t_result_cut4_Mid", "t dist 15 < Q^{2} < 20,  t cut, #theta_{n} cut, #theta_{diff} , phi_{diff}, p_{miss} cuts; -t(GeV^{2})", 10, 0, 0.4);
  h1_t_cut1_High = new TH1F("t_result_cut1_High", "t dist 25 < Q^{2} < 30, t cut only; -t(GeV^{2})", 10, 0, 0.4);
h1_t_cut2_High = new TH1F("t_result_cut2_High", "t dist 25 < Q^{2} < 30, t cut, #theta_{n} cut; -t(GeV^{2})", 10, 0, 0.4);
 h1_t_cut3_High = new TH1F("t_result_cut3_High", "t dist 25 < Q^{2} < 30, t cut, #theta_{n} cut, #theta_{diff} , phi_{diff} cuts; -t(GeV^{2})", 10, 0, 0.4);
h1_t_cut4_High = new TH1F("t_result_cut4_High", "t dist 25 < Q^{2} < 30,  t cut, #theta_{n} cut, #theta_{diff} , phi_{diff}, p_{miss} cuts; -t(GeV^{2})", 10, 0, 0.4);

  gDirectory->cd("../");
  gDirectory->mkdir("t_Resolution");
  gDirectory->cd("t_Resolution");
  for(Int_t A = 0; A < 14; A++){  
    h1_t_Resolution[A]=new TH1F(Form("t_Resolution_Q2_%i", (A+1)), Form("t - t_{truth} Dist, %2.1f < Q^{2} < %2.1f;  t - t_{truth} (GeV^{2}); Rate(Hz)", Q2BinVal[A], Q2BinVal[A+1]), 150, -5.5, 0.5);
    h1_talt_Resolution_ZDC[A]=new TH1F(Form("talt_Resolution_ZDC_Q2_%i", (A+1)), Form("t_{alt} - t_{truth} Dist, %2.1f < Q^{2} < %2.1f, ZDC Info Only;  t_{alt} - t_{truth} (GeV^{2}); Rate(Hz)", Q2BinVal[A], Q2BinVal[A+1]), 25, -0.5, 0.5);
    h1_talt_Resolution_pMiss[A]=new TH1F(Form("talt_Resolution_pMiss_Q2_%i", (A+1)), Form("t_{alt} - t_{truth} Dist, %2.1f < Q^{2} < %2.1f, Corrected p_{Miss};  t_{alt} - t_{truth} (GeV^{2}); Rate(Hz)", Q2BinVal[A], Q2BinVal[A+1]), 25, -0.5, 0.5);
  } // mfek 06/22/2022 - changed binning

  gDirectory->cd("../");
  gDirectory->mkdir("pi_e_n_preCut_Truth_Dist"); // mfek 06/23/2022 - new directory for precut truth values
  gDirectory->cd("pi_e_n_preCut_Truth_Dist");

  h1_piTruth_preCut_px = new TH1F("piTruth_preCut_px", "PreCut #pi Truth p_{x} Distribution;p_{x} (GeV)", 200, -20, 20);
  h1_piTruth_preCut_py = new TH1F("piTruth_preCut_py", "PreCut #pi Truth p_{y} Distribution;p_{y} (GeV)", 200, -20, 20);
  h1_piTruth_preCut_pz = new TH1F("piTruth_preCut_pz", "PreCut #pi Truth p_{z} Distribution;p_{z} (GeV)", 200, -60, 60); 
  h1_piTruth_preCut_p = new TH1F("piTruth_preCut_p", "PreCut #pi Truth p Distribution;p (GeV)", 200, 0, 60);
  h1_piTruth_preCut_E = new TH1F("piTruth_preCut_E", "PreCut #pi Truth E Distribution;E (GeV)", 200, 0, 60);
  h1_piTruth_preCut_Theta = new TH1F("piTruth_preCut_Theta", "PreCut #pi Truth #theta Distribution;#theta (deg)", 160, 0, 50);
  h1_piTruth_preCut_Phi = new TH1F("piTruth_preCut_Phi", "PreCut #pi Truth #phi Distribution;#phi (deg)", 360, -180, 180);

  h1_eTruth_preCut_px = new TH1F("eTruth_preCut_px", "PreCut e' Truth p_{x} Distribution;p_{x} (GeV)", 240, -6, 6);
  h1_eTruth_preCut_py = new TH1F("eTruth_preCut_py", "PreCut e' Truth p_{y} Distribution;p_{y} (GeV)", 240, -6, 6);
  h1_eTruth_preCut_pz = new TH1F("eTruth_preCut_pz", "PreCut e' Truth p_{z} Distribution;p_{z} (GeV)", 120, -15, 0); 
  h1_eTruth_preCut_p = new TH1F("eTruth_preCut_p", "PreCut e' Truth p Distribution;p (GeV)", 160, 0, 15);
  h1_eTruth_preCut_E = new TH1F("eTruth_preCut_E", "PreCut e' Truth E Distribution;E (GeV)", 160, 0, 15);
  h1_eTruth_preCut_Theta = new TH1F("eTruth_preCut_Theta", "PreCut e' Truth #theta Distribution;#theta (deg)", 160, 110, 180);
  h1_eTruth_preCut_Phi = new TH1F("eTruth_preCut_Phi", "PreCut e' Truth #phi Distribution;#phi (deg)", 360, -180, 180);

  h1_nTruth_preCut_px = new TH1F("nTruth_preCut_px", "PreCut nTruth p_{x} Distribution;p_{x} (GeV)", 480, -6, 6);
  h1_nTruth_preCut_py = new TH1F("nTruth_preCut_py", "PreCut nTruth p_{y} Distribution;p_{y} (GeV)", 200, -2.5, 2.5);
  h1_nTruth_preCut_pz = new TH1F("nTruth_preCut_pz", "PreCut nTruth p_{z} Distribution;p_{z} (GeV)", 240, 0, 120);
  h1_nTruth_preCut_p = new TH1F("nTruth_preCut_p", "PreCut nTruth p Distribution;p (GeV)", 240, 0, 120);
  h1_nTruth_preCut_E = new TH1F("nTruth_preCut_E", "PreCut nTruth E Distribution;E (GeV)", 240, 0, 120);
  h1_nTruth_preCut_Theta = new TH1F("nTruth_preCut_Theta", "PreCut nTruth #theta Distribution;#theta (deg)", 180, 0, 5);
  h1_nTruth_preCut_Theta_inRange = new TH1F("nTruth_preCut_Theta_inRange", "PreCut nTruth #theta Distribution in Range;#theta (deg)", 50, 1, 2);
  h1_nTruth_preCut_Phi = new TH1F("nTruth_preCut_Phi", "PreCut nTruth #phi Distribution;#phi (deg)", 360, -180, 180);
  h2_nTruth_preCut_XY = new TH2F("nTruth_preCut_XY", "PreCut nTruth XY Distribution; x (cm); y (cm)", 200, -200, 0, 200, -150, 150); 
  //h2_nTruth_preCut_XY_inZDC = new TH2F("nTruth_preCut_XY_inZDC", "PreCut nTruth XY Distribution in ZDC range; x (cm); y (cm)", 200, -200, 0, 200, -150, 150); 
  //h2_nTruth_preCut_XY_outZDC = new TH2F("nTruth_preCut_XY_outZDC", "PreCut nTruth XY Distribution outside ZDC range; x (cm); y (cm)", 200, -200, 0, 200, -150, 150);
  //h2_nTruth_XY_hits = new TH2F("nTruth_preCut_XY_hits", "PreCut nTruth XY Distribution for ZDC hits; x (cm); y (cm)", 200, -200, 0, 200, -150, 150);

  gDirectory->cd("../");
  gDirectory->mkdir("pi_e_n_Missed_Truth_Dist"); // mfek 06/23/2022 - new directory for truth values of missed events
  gDirectory->cd("pi_e_n_Missed_Truth_Dist");

  h1_piTruth_Missed_px = new TH1F("piTruth_Missed_px", "Missed #pi Truth p_{x} Distribution;p_{x} (GeV)", 200, -20, 20);
  h1_piTruth_Missed_py = new TH1F("piTruth_Missed_py", "Missed #pi Truth p_{y} Distribution;p_{y} (GeV)", 200, -20, 20);
  h1_piTruth_Missed_pz = new TH1F("piTruth_Missed_pz", "Missed #pi Truth p_{z} Distribution;p_{z} (GeV)", 200, -60, 60); 
  h1_piTruth_Missed_p = new TH1F("piTruth_Missed_p", "Missed #pi Truth p Distribution;p (GeV)", 200, 0, 60);
  h1_piTruth_Missed_E = new TH1F("piTruth_Missed_E", "Missed #pi Truth E Distribution;E (GeV)", 200, 0, 60);
  h1_piTruth_Missed_Theta = new TH1F("piTruth_Missed_Theta", "Missed #pi Truth #theta Distribution;#theta (deg)", 160, 0, 50);
  h1_piTruth_Missed_Phi = new TH1F("piTruth_Missed_Phi", "Missed #pi Truth #phi Distribution;#phi (deg)", 360, -180, 180);

  h1_eTruth_Missed_px = new TH1F("eTruth_Missed_px", "Missed e' Truth p_{x} Distribution;p_{x} (GeV)", 240, -6, 6);
  h1_eTruth_Missed_py = new TH1F("eTruth_Missed_py", "Missed e' Truth p_{y} Distribution;p_{y} (GeV)", 240, -6, 6);
  h1_eTruth_Missed_pz = new TH1F("eTruth_Missed_pz", "Missed e' Truth p_{z} Distribution;p_{z} (GeV)", 120, -15, 0); 
  h1_eTruth_Missed_p = new TH1F("eTruth_Missed_p", "Missed e' Truth p Distribution;p (GeV)", 160, 0, 15);
  h1_eTruth_Missed_E = new TH1F("eTruth_Missed_E", "Missed e' Truth E Distribution;E (GeV)", 160, 0, 15);
  h1_eTruth_Missed_Theta = new TH1F("eTruth_Missed_Theta", "Missed e' Truth #theta Distribution;#theta (deg)", 160, 110, 180);
  h1_eTruth_Missed_Phi = new TH1F("eTruth_Missed_Phi", "Missed e' Truth #phi Distribution;#phi (deg)", 360, -180, 180);

  h1_nTruth_Missed_px = new TH1F("nTruth_Missed_px", "Missed nTruth p_{x} Distribution;p_{x} (GeV)", 480, -6, 6);
  h1_nTruth_Missed_py = new TH1F("nTruth_Missed_py", "Missed nTruth p_{y} Distribution;p_{y} (GeV)", 200, -2.5, 2.5);
  h1_nTruth_Missed_pz = new TH1F("nTruth_Missed_pz", "Missed nTruth p_{z} Distribution;p_{z} (GeV)", 240, 0, 120);
  h1_nTruth_Missed_p = new TH1F("nTruth_Missed_p", "Missed nTruth p Distribution;p (GeV)", 240, 0, 120);
  h1_nTruth_Missed_E = new TH1F("nTruth_Missed_E", "Missed nTruth E Distribution;E (GeV)", 240, 0, 120);
  h1_nTruth_Missed_Theta = new TH1F("nTruth_Missed_Theta", "Missed nTruth #theta Distribution;#theta (deg)", 160, 0, 5);
  h1_nTruth_Missed_Theta_inRange = new TH1F("nTruth_Missed_Theta_inRange", "Missed nTruth #theta Distribution in Range;#theta (deg)", 50, 1, 2);
  h1_nTruth_Missed_Phi = new TH1F("nTruth_Missed_Phi", "Missed nTruth #phi Distribution;#phi (deg)", 360, -180, 180);
  h2_nTruth_Missed_XY = new TH2F("nTruth_Missed_XY", "Missed nTruth XY Distribution; x (cm); y (cm)", 200, -200, 0, 200, -150, 150);

  h1_Q2Truth_Dist_Missed = new TH1F("Q2Truth_Dist_Missed", "Missed Q^{2} Truth Distribution", 200, 0, 50);
  h1_tTruth_Dist_Missed = new TH1F("tTruth_Dist_Missed", "Missed t Truth Distribution", 100, 0, 1);

  gDirectory->cd("../");
  h2_ZDC_XY_IP6 = new TH2F("ZDC_XY_IP6", "n X vs Y at ZDC; x (cm); y (cm)", 200, -150, -50, 200, -50, 50);
  h2_ZDC_XY_IP8 = new TH2F("ZDC_XY_IP8", "n X vs Y at ZDC; x (cm); y (cm)", 200, 50, 150, 200, -50, 50);
  h2_ZDC_XY_l = new TH2F("ZDC_XY_l", "n X vs Y at ZDC (Local Co-ords); x (cm); y (cm)", 800, -200, 200, 200, -50, 50);
  //nTruth_xyE3D = new TH3F("nTruth_xyE3D", "n X vs Y vs E Truth Distribution; x (cm); y (cm); E (GeV)", 100, -200, 0, 100, -150, 150, 100, 0, 50); // mfek 06/23/2022
  //ZDC_xyE3D = new TH3F("ZDC_xyE3D", "n X vs Y vs E ZDC Distribution; x (cm); y (cm); E (GeV)", 100, -200, 0, 100, -150, 150, 100, 0, 50); // mfek 06/23/2022
  //nTruth_Missed_xyE3D = new TH3F("nTruth_Missed_xyE3D", "n X vs Y vs E Truth Distribution for missed events; x (cm); y (cm); E (GeV)", 100, -200, 0, 100, -150, 150, 100, 0, 50); // mfek 06/23/2022
  //h3nTruth_xyE_pxy = new TProfile2D("h3nTruth_xyE_pxy", "n X vs Y vs E Truth Distribution; y (cm); x (cm)", 100, -150, 150, 100, -200, 0, 0.0, 50); // mfek 06/23/2022
  //h3ZDC_xyE_pxy = new TProfile2D("h3ZDC_xyE_pxy", "n X vs Y vs E ZDC Distribution; y (cm); x (cm)", 100, -150, 150, 100, -200, 0, 0.0, 50); // mfek 06/23/2022
  //h3nTruth_Missed_xyE_pxy = new TProfile2D("h3nTruth_Missed_xyE_pxy", "n X vs Y vs E Truth Distribution for missed events; y (cm); x (cm)", 100, -150, 150, 100, -200, 0, 0.0, 50); // mfek 06/03/2022

  // Define beam 4 vectors - Assume IP6 by default, for other IPs, adjust in the Process Event loop (at the top)
  e_beam_energy = 10;
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

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP10on100::InitRun(PHCompositeNode *topNode)
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
int ECCE_DEMP10on100::process_event(PHCompositeNode *topNode)
{
  ZDC_hit = 0;
  EEMC_hit = 0;
  event_itt++;
  thrownEvents++; // mfek 06/22/2022 
  
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

  // mfek 06/23/2022
  h1_piTruth_preCut_px->Fill(pi4VectTruth.Px(), wgt);
  h1_piTruth_preCut_py->Fill(pi4VectTruth.Py(), wgt);
  h1_piTruth_preCut_pz->Fill(pi4VectTruth.Pz(), wgt);
  h1_piTruth_preCut_p->Fill(pi4VectTruth.P(), wgt);
  h1_piTruth_preCut_E->Fill(pi4VectTruth.E(), wgt);
  h1_piTruth_preCut_Theta->Fill(pi4VectTruth.Theta()*TMath::RadToDeg(), wgt);
  h1_piTruth_preCut_Phi->Fill(pi4VectTruth.Phi()*TMath::RadToDeg(), wgt);

  h1_eTruth_preCut_px->Fill(e4VectTruth.Px(), wgt);
  h1_eTruth_preCut_py->Fill(e4VectTruth.Py(), wgt);
  h1_eTruth_preCut_pz->Fill(e4VectTruth.Pz(), wgt);
  h1_eTruth_preCut_p->Fill(e4VectTruth.P(), wgt);
  h1_eTruth_preCut_E->Fill(e4VectTruth.E(), wgt);
  h1_eTruth_preCut_Theta->Fill(e4VectTruth.Theta()*TMath::RadToDeg(), wgt);
  h1_eTruth_preCut_Phi->Fill(e4VectTruth.Phi()*TMath::RadToDeg(), wgt);

  h1_nTruth_preCut_px->Fill(n4VectTruth.Px(), wgt);
  h1_nTruth_preCut_py->Fill(n4VectTruth.Py(), wgt);
  h1_nTruth_preCut_pz->Fill(n4VectTruth.Pz(), wgt);
  h1_nTruth_preCut_p->Fill(n4VectTruth.P(), wgt);
  h1_nTruth_preCut_E->Fill(n4VectTruth.E(), wgt);
  h1_nTruth_preCut_Theta->Fill(n4VectTruth.Theta()*TMath::RadToDeg(), wgt);
  if (n4VectTruth.Theta()*TMath::RadToDeg() >= 1.2 && n4VectTruth.Theta()*TMath::RadToDeg() <= 1.66) {
    h1_nTruth_preCut_Theta_inRange->Fill(n4VectTruth.Theta()*TMath::RadToDeg(), wgt);
  }
  h1_nTruth_preCut_Phi->Fill(n4VectTruth.Phi()*TMath::RadToDeg(), wgt);
  h2_nTruth_preCut_XY->Fill(3700*TMath::Cos(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()), 3700*TMath::Sin(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()), wgt); // filling neutron truth XY distribution

  /*if (3700*TMath::Cos(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()) > -150 && 3700*TMath::Cos(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()) < -50 && 3700*TMath::Sin(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()) > -50 && 3700*TMath::Sin(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()) < 50) {
    h2_nTruth_preCut_XY_inZDC->Fill(3700*TMath::Cos(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()), 3700*TMath::Sin(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()), wgt);
  }
  else {
    h2_nTruth_preCut_XY_outZDC->Fill(3700*TMath::Cos(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()), 3700*TMath::Sin(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()), wgt);
  }*/

  //nTruth_xyE3D->Fill(3700*TMath::Cos(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()), 3700*TMath::Sin(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()), n4VectTruth.E()/*, wgt*/);

  h1_Q2_DetEff_Uncut->Fill(Q2_truth, wgt);
  h2_Q2_t_DetEff_Uncut->Fill(Q2_truth, t_truth, wgt);
  h2_Q2_t_DetEff_v2_Uncut->Fill(Q2_truth, t_truth, wgt);

  for(Int_t B = 0; B < 14; B++){
    Q2_low = Q2BinVal[B];
    Q2_high = Q2BinVal[B+1];
    if ( Q2_truth > Q2_low && Q2_truth < Q2_high){
      h1_t_truth_thrown_result[B]->Fill(t_truth, wgt);
    }
  } // mfek 06/22/2022 - new binning

  if (Check_n(topNode) == true) {
    count_aftern++;
  } // mfek 06/22/2022 - new counter added

  if (Check_ePi(topNode) == true) {
    count_afterePi++;
  } // mfek 06/22/2022 - new counter added

  if (Check_n(topNode) != true || Check_ePi(topNode) != true) {
    //h1_Q2_DetEff_Missed->Fill(Q2_truth, wgt);
    //h2_Q2_t_DetEff_Missed->Fill(Q2_truth, t_truth, wgt);
    //h2_Q2_t_DetEff_Missed_v2->Fill(Q2_truth, t_truth, wgt);

    h1_piTruth_Missed_p->Fill(pi4VectTruth.P(), wgt);
    h1_piTruth_Missed_px->Fill(pi4VectTruth.Px(), wgt);
    h1_piTruth_Missed_py->Fill(pi4VectTruth.Py(), wgt);
    h1_piTruth_Missed_pz->Fill(pi4VectTruth.Pz(), wgt);
    h1_piTruth_Missed_E->Fill(pi4VectTruth.E(), wgt);
    h1_piTruth_Missed_Theta->Fill(pi4VectTruth.Theta()*TMath::RadToDeg(), wgt);
    h1_piTruth_Missed_Phi->Fill(pi4VectTruth.Phi()*TMath::RadToDeg(), wgt);

    h1_eTruth_Missed_p->Fill(e4VectTruth.P(), wgt);
    h1_eTruth_Missed_px->Fill(e4VectTruth.Px(), wgt);
    h1_eTruth_Missed_py->Fill(e4VectTruth.Py(), wgt);
    h1_eTruth_Missed_pz->Fill(e4VectTruth.Pz(), wgt);
    h1_eTruth_Missed_E->Fill(e4VectTruth.E(), wgt);
    h1_eTruth_Missed_Theta->Fill(e4VectTruth.Theta()*TMath::RadToDeg(), wgt);
    h1_eTruth_Missed_Phi->Fill(e4VectTruth.Phi()*TMath::RadToDeg(), wgt);

    h1_nTruth_Missed_p->Fill(n4VectTruth.P(), wgt);
    h1_nTruth_Missed_px->Fill(n4VectTruth.Px(), wgt);
    h1_nTruth_Missed_py->Fill(n4VectTruth.Py(), wgt);
    h1_nTruth_Missed_pz->Fill(n4VectTruth.Pz(), wgt);
    h1_nTruth_Missed_E->Fill(n4VectTruth.E(), wgt);
    h1_nTruth_Missed_Theta->Fill(n4VectTruth.Theta()*TMath::RadToDeg(), wgt);
    if (n4VectTruth.Theta()*TMath::RadToDeg() >= 1.2 && n4VectTruth.Theta()*TMath::RadToDeg() <= 1.66) {
      h1_nTruth_Missed_Theta_inRange->Fill(n4VectTruth.Theta()*TMath::RadToDeg(), wgt);
    }
    h1_nTruth_Missed_Phi->Fill(n4VectTruth.Phi()*TMath::RadToDeg(), wgt);

    h2_nTruth_Missed_XY->Fill(3700*TMath::Cos(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()), 3700*TMath::Sin(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()), wgt);

    //nTruth_Missed_xyE3D->Fill(3700*TMath::Cos(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()), 3700*TMath::Sin(n4VectTruth.Phi())*TMath::Tan(n4VectTruth.Theta()), n4VectTruth.E()/*, wgt*/);

    h1_Q2Truth_Dist_Missed->Fill(Q2_truth, wgt);
    h1_tTruth_Dist_Missed->Fill(t_truth, wgt);
  } // mfek 06/23/2022


  if (Check_n(topNode) == true && Check_ePi(topNode) == true){ // For event, check if it look like we have an e/pi/n in the event
    // Get track map for e'/pi info
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    // Get ZDC hits for neutron info
    //PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDC");
    // Need to use ZDC surrogate now
    PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDCsurrogate");
    cut1Events++; // mfek 06/22/2022 - new counter added (counts after applying ZDC energy deposit cut)

    if (!trackmap)
      {
	trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
	if (!trackmap)
	  {
	    cout
	      << "ECCE_DEMP10on100::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
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
    
    // Fill weighted histograms
    h1_Q2_DetEff_Cut->Fill(Q2_truth, wgt);
    h2_Q2_t_DetEff_Cut->Fill(Q2_truth, t_truth, wgt);
    h2_Q2_t_DetEff_v2_Cut->Fill(Q2_truth, t_truth, wgt);

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

    h1_nRec_px->Fill(nRec4Vect.Px(), wgt);
    h1_nRec_py->Fill(nRec4Vect.Py(), wgt);
    h1_nRec_pz->Fill(nRec4Vect.Pz(), wgt);
    h1_nRec_p->Fill(nRec4Vect.P(), wgt);
    h1_nRec_E->Fill(nRec4Vect.E(), wgt);
    h1_nRec_Theta->Fill(nRec4Vect.Theta()*TMath::RadToDeg(), wgt);
    h1_nRec_Phi->Fill(nRec4Vect.Phi()*TMath::RadToDeg(), wgt);
    h1_n_ThetaDiff->Fill(nTheta_Diff*TMath::RadToDeg(), wgt);
    h1_n_PhiDiff->Fill(nPhi_Diff*TMath::RadToDeg(), wgt);
    h2_n_ThetaPhiDiff->Fill(nTheta_Diff*TMath::RadToDeg(), nPhi_Diff*TMath::RadToDeg(), wgt);

    h1_pi_px_Unweighted->Fill(pi4Vect.Px());
    h1_pi_py_Unweighted->Fill(pi4Vect.Py());
    h1_pi_pz_Unweighted->Fill(pi4Vect.Pz());
    h1_pi_p_Unweighted->Fill(pi4Vect.P());
    h1_pi_E_Unweighted->Fill(pi4Vect.E());
    h1_pi_Theta_Unweighted->Fill(pi4Vect.Theta()*TMath::RadToDeg());
    h1_pi_Phi_Unweighted->Fill(pi4Vect.Phi()*TMath::RadToDeg());
    h1_e_px_Unweighted->Fill(e4Vect.Px());
    h1_e_py_Unweighted->Fill(e4Vect.Py());
    h1_e_pz_Unweighted->Fill(e4Vect.Pz());
    h1_e_p_Unweighted->Fill(e4Vect.P());
    h1_e_E_Unweighted->Fill(e4Vect.E());
    h1_e_Theta_Unweighted->Fill(e4Vect.Theta()*TMath::RadToDeg());
    h1_e_Phi_Unweighted->Fill(e4Vect.Phi()*TMath::RadToDeg());
    // 31/03/22 - Use unweighted version to compare directly to ATHENA
    h1_n_px_Unweighted->Fill(n4Vect.Px());
    h1_n_py_Unweighted->Fill(n4Vect.Py());
    h1_n_pz_Unweighted->Fill(n4Vect.Pz());
    h1_n_p_Unweighted->Fill(n4Vect.P());
    h1_n_E_Unweighted->Fill(n4Vect.E());
    h1_n_Theta_Unweighted->Fill(n4Vect.Theta()*TMath::RadToDeg());
    h1_n_Phi_Unweighted->Fill(n4Vect.Phi()*TMath::RadToDeg());
    h1_n_ThetaDiff_Unweighted->Fill(nTheta_Diff*TMath::RadToDeg());
    h1_n_PhiDiff_Unweighted->Fill(nPhi_Diff*TMath::RadToDeg());
    h2_n_ThetaPhiDiff_Unweighted->Fill(nTheta_Diff*TMath::RadToDeg(), nPhi_Diff*TMath::RadToDeg());

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

    h1_taltres_result->Fill((t_alt-t_truth),wgt);

    for(Int_t B = 0; B < 10; B++){
      t_low = 0+(B*0.05);
      t_high = 0.05+(B*0.05);
      if ( t_truth > t_low && t_truth < t_high){
	h1_taltres_result_ttruth[B]->Fill((t_alt-t_truth),wgt);
      }
    }

    h2_t_ep->Fill(t, e4Vect.P(), wgt);
    h2_t_Q2->Fill(t,Q2, wgt);
    h2_delta_t_t->Fill(((t - t_truth)/t_truth)*100, t, wgt);

    for(Int_t B = 0; B < 14; B++){
      Q2_low = Q2BinVal[B];
      Q2_high = Q2BinVal[B+1];
      if ( Q2_truth > Q2_low && Q2_truth < Q2_high){
	  h1_t_Q2[B]->Fill(t, wgt);
	  h1_t_alt_Q2[B]->Fill(t_alt, wgt);
	  h2_delta_t_t_Q2[B]->Fill(((t - t_truth)/t_truth)*100, t, wgt);
      }
    } // mfek 06/22/2022 - changed binning for 5 on 41

    for(Int_t B = 0; B < 14; B++){
      Q2_low = Q2BinVal[B];
      Q2_high = Q2BinVal[B+1]; // mfek 06/22/2022 - new binning

      // Fill some histograms under different cut conditions for comparison
      // These are used for some plots Garth wanted
      if (Q2_low == 7.5){
	if ( Q2 > Q2_low && Q2 < Q2_high){
	  if (t_alt < 0.4){ 
	    h1_t_cut1_Low->Fill(t_alt, wgt);
	  }
	  if (t_alt < 0.4 && ((nRec4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (nRec4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5)))){ 
	    h1_t_cut2_Low->Fill(t_alt, wgt);
	  }
	  if (t_alt < 0.4 && ((nRec4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (nRec4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5))) && (abs(nTheta_Diff*TMath::RadToDeg())) < ThetaDiff_Cut && (abs(nPhi_Diff*TMath::RadToDeg())) < PhiDiff_Cut){ 
	    h1_t_cut3_Low->Fill(t_alt, wgt);
	  }
	  if (t_alt < 0.4 && ((nRec4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (nRec4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5))) && (abs(nTheta_Diff*TMath::RadToDeg())) < ThetaDiff_Cut && (abs(nPhi_Diff*TMath::RadToDeg())) < PhiDiff_Cut && (nRec4Vect.P() < PmissCutVal[B])){ 
	    h1_t_cut4_Low->Fill(t_alt, wgt);	  
	  }
	}
      }
      else if (Q2_low == 15.0){
	if ( Q2 > Q2_low && Q2 < Q2_high){
	  if (t_alt < 0.4){ 
	    h1_t_cut1_Mid->Fill(t_alt, wgt);
	  }
	  if (t_alt < 0.4 && ((nRec4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (nRec4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5)))){ 
	    h1_t_cut2_Mid->Fill(t_alt, wgt);
	  }
	  if (t_alt < 0.4 && ((nRec4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (nRec4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5))) && (abs(nTheta_Diff*TMath::RadToDeg())) < ThetaDiff_Cut && (abs(nPhi_Diff*TMath::RadToDeg())) < PhiDiff_Cut){ 
	    h1_t_cut3_Mid->Fill(t_alt, wgt);
	  }
	  if (t_alt < 0.4 && ((nRec4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (nRec4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5))) && (abs(nTheta_Diff*TMath::RadToDeg())) < ThetaDiff_Cut && (abs(nPhi_Diff*TMath::RadToDeg())) < PhiDiff_Cut && (nRec4Vect.P() < PmissCutVal[B])){ 
	    h1_t_cut4_Mid->Fill(t_alt, wgt);	  
	  }
	}
      }
      else if (Q2_low == 25.0){
	if ( Q2 > Q2_low && Q2 < Q2_high){
	  if (t_alt < 0.4){ 
	    h1_t_cut1_High->Fill(t_alt, wgt);
	  }
	  if (t_alt < 0.4 && ((nRec4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (nRec4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5)))){ 
	    h1_t_cut2_High->Fill(t_alt, wgt);
	  }
	  if (t_alt < 0.4 && ((nRec4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (nRec4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5))) && (abs(nTheta_Diff*TMath::RadToDeg())) < ThetaDiff_Cut && (abs(nPhi_Diff*TMath::RadToDeg())) < PhiDiff_Cut){ 
	    h1_t_cut3_High->Fill(t_alt, wgt);
	  }
	  if (t_alt < 0.4 && ((nRec4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (nRec4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5))) && (abs(nTheta_Diff*TMath::RadToDeg())) < ThetaDiff_Cut && (abs(nPhi_Diff*TMath::RadToDeg())) < PhiDiff_Cut && (nRec4Vect.P() < PmissCutVal[B])){ 
	    h1_t_cut4_High->Fill(t_alt, wgt);	  
	  }
	}
      }
      // nTheta distribution with the t cut only
      if(t_alt < 0.4){
	h1_nTheta_tCut->Fill(nRec4Vect.Theta()*TMath::RadToDeg(), wgt);
      }
      
      // t resolution plots under different conditions - we treat each different t value as though it is the "true" value of t
      if ( Q2 > Q2_low && Q2 < Q2_high){
	if ( (pmiss4Vect.P() < PmissCutVal[B]) && ((pmiss4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (pmiss4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5))) && t < 0.4){
	  h1_t_Resolution[B]->Fill(t-t_truth, wgt);
	}
      }
      if ( Q2 > Q2_low && Q2 < Q2_high){
	if ( (n4Vect.P() < PmissCutVal[B]) && ((n4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (n4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5))) && (abs(nTheta_Diff*TMath::RadToDeg())) < ThetaDiff_Cut && (abs(nPhi_Diff*TMath::RadToDeg())) < PhiDiff_Cut && t_alt_ZDC < 0.4){
	  h1_talt_Resolution_ZDC[B]->Fill(t_alt_ZDC-t_truth, wgt);
	}
      }
      if ( Q2 > Q2_low && Q2 < Q2_high){
	if ( (nRec4Vect.P() < PmissCutVal[B]) && ((nRec4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (nRec4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5))) && (abs(nTheta_Diff*TMath::RadToDeg())) < ThetaDiff_Cut && (abs(nPhi_Diff*TMath::RadToDeg())) < PhiDiff_Cut && t_alt < 0.4){
	  h1_talt_Resolution_pMiss[B]->Fill(t_alt-t_truth, wgt);
	}
      }

      // The main analysis is below
      if ( Q2 > Q2_low && Q2 < Q2_high){
	if (t_alt < 0.4){
	  h1_t_result[B]->Fill(t_alt, wgt);
	  h1_nTheta_result[B]->Fill(nRec4Vect.Theta()*TMath::RadToDeg(), wgt);
	  h1_pmiss_result[B]->Fill(pmiss4Vect.P(), wgt);
	  h1_pn_result[B]->Fill(nRec4Vect.P(), wgt);
          cut2Events++; // mfek 06/22/2022 - t cut
	  // Apply other cuts
	  //if ( (pmiss4Vect.P() < PmissCutVal[B]) && ((nSmeared4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (nSmeared4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5)))){
	  // SJDK - 31/03/22
	  // nRec4Vect is EXACTLY the missing momentum vector for now, so comparing it to the pmiss cut value is fine
	  //if ( (nRec4Vect.P() < PmissCutVal[B]) && ((nRec4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (nRec4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5)))){ // Old version without theta/phi diff cuts
	  if ((nRec4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (nRec4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5))) {
	    cut3Events++;
	    if (abs(nTheta_Diff*TMath::RadToDeg()) < ThetaDiff_Cut && (abs(nPhi_Diff*TMath::RadToDeg())) < PhiDiff_Cut) {
	      cut4Events++;
	    } // ThetaDiff/PhiDiff CUt
	  } // Thetan Cut 
	  // mfek 06/22/20222 - added theta diff/phi diff and thetan cut counters
	  if ( (nRec4Vect.P() < PmissCutVal[B]) && ((nRec4Vect.Theta()*TMath::RadToDeg() > (Thetan_Cent-0.5)) && (nRec4Vect.Theta()*TMath::RadToDeg() < (Thetan_Cent+0.5))) && (abs(nTheta_Diff*TMath::RadToDeg())) < ThetaDiff_Cut && (abs(nPhi_Diff*TMath::RadToDeg())) < PhiDiff_Cut){
	    h1_t_cut_result[B]->Fill(t_alt, wgt);
	    h1_Q2_cut_result[B]->Fill(Q2, wgt);
	    h1_W_cut_result[B]->Fill(W, wgt);
	    h2_t_ttruth_result->Fill(t, t_truth, wgt);
	    h2_t_alt_ttruth_result->Fill(t_alt, t_truth, wgt);
	    h1_Mmiss_result->Fill(pmiss4Vect_2.M(),wgt);
	    h1_MmissSq_result->Fill(((pmiss4Vect_2.M()*(pmiss4Vect_2.M()))),wgt);
	    h1_Mmiss_truth_result->Fill(pmiss4VectTruth_2.M(),wgt);
	    h1_Mmiss_Comp_result->Fill((pmiss4VectTruth_2.M()-pmiss4Vect_2.M()), wgt);
	    h2_t_t_alt_result->Fill(t, t_alt, wgt);
	    h2_Q2_W_result->Fill(Q2, W, wgt);
	    h2_Q2_t_result[B]->Fill(Q2, t_alt, wgt);
	    h2_t_ttruth_result_Q2[B]->Fill(t, t_truth, wgt);
	    h2_t_alt_ttruth_result_Q2[B]->Fill(t_alt, t_truth, wgt);
	    cut5Events++; // mfek 06/22/2022 - PMiss cut
	  }
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
    h1_nRes_p->Fill((nRec4Vect.P()-n4VectTruth.P())/(n4VectTruth.P())*100, wgt);
    h1_nRes_px->Fill((nRec4Vect.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100, wgt);
    h1_nRes_py->Fill((nRec4Vect.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100, wgt);
    h1_nRes_pz->Fill((nRec4Vect.Pz()-n4VectTruth.Pz())/(n4VectTruth.Pz())*100, wgt);
    h1_pmissDiff_p->Fill(pmiss4Vect.P()-n4VectTruth.P(), wgt);
    h1_pmissDiff_px->Fill(pmiss4Vect.Px()-n4VectTruth.Px(), wgt);
    h1_pmissDiff_py->Fill(pmiss4Vect.Py()-n4VectTruth.Py(), wgt);
    h1_pmissDiff_pz->Fill(pmiss4Vect.Pz()-n4VectTruth.Pz(), wgt);
    
    // SJDK 30/05/22 The 375 used here is the approximate z position of the central barrel calorimeter in cm. This projects the tracks to that position (since the calorimeter hits didn't work)
    // This isn't exactly coded fantastically, I think this was just here to make a quick plot at one point
    h2_pi_XY->Fill((375*(TMath::Cos(pi4Vect.Phi()))*(TMath::Tan(pi4Vect.Theta()))), (375*(TMath::Sin(pi4Vect.Phi()))*(TMath::Tan(pi4Vect.Theta()))), wgt);
    h2_e_XY->Fill((375*(TMath::Cos(e4Vect.Phi()))*(TMath::Tan(e4Vect.Theta()))), (375*(TMath::Sin(e4Vect.Phi()))*(TMath::Tan(e4Vect.Theta()))), wgt);
    if ( IP_design == "IP6"){
      h2_ZDC_XY_IP6->Fill(nZDCPos.x(), nZDCPos.y(), wgt);
    }
    else if ( IP_design == "IP8"){
      h2_ZDC_XY_IP8->Fill(nZDCPos.x(), nZDCPos.y(), wgt);
    }
    else{
      h2_ZDC_XY_IP6->Fill(nZDCPos.x(), nZDCPos.y(), wgt);
    }

    h2_ZDC_XY_l->Fill(local_x, local_y, wgt);

    //ZDC_xyE3D->Fill(nZDCPos.x(), nZDCPos.y(), n4Vect.E(), wgt); // mfek 06/23/2022
 
    h2_piTrack_ThetaPhi->Fill((pi4Vect.Theta()*TMath::RadToDeg()), (pi4Vect.Phi()*TMath::RadToDeg()), wgt);
    h2_piTrack_pTheta->Fill((pi4Vect.Theta()*TMath::RadToDeg()), pi4Vect.P(), wgt);
    h2_eTrack_ThetaPhi->Fill((e4Vect.Theta()*TMath::RadToDeg()), (e4Vect.Phi()*TMath::RadToDeg()), wgt);
    h2_eTrack_pTheta->Fill((e4Vect.Theta()*TMath::RadToDeg()), e4Vect.P(), wgt);
    h2_nTrack_ThetaPhi->Fill((nRec4Vect.Theta()*TMath::RadToDeg()), (nRec4Vect.Phi()*TMath::RadToDeg()), wgt);
    h2_nTrack_pTheta->Fill((nRec4Vect.Theta()*TMath::RadToDeg()), nRec4Vect.P(), wgt);
    
    h2_piTruth_pxpy->Fill((pi4Vect.Px()-pi4VectTruth.Px())/(pi4VectTruth.Px())*100, (pi4Vect.Py()-pi4VectTruth.Py())/(pi4VectTruth.Py())*100, wgt);
    h2_eTruth_pxpy->Fill((e4Vect.Px()-e4VectTruth.Px())/(e4VectTruth.Px())*100, (e4Vect.Py()-e4VectTruth.Py())/(e4VectTruth.Py())*100, wgt);
    h2_nTruth_pxpy->Fill((nRec4Vect.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100, (nRec4Vect.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100, wgt);

    h2_piTruth_pxpz->Fill((pi4Vect.Px()-pi4VectTruth.Px())/(pi4VectTruth.Px())*100, (pi4Vect.Pz()-pi4VectTruth.Pz())/(pi4VectTruth.Pz())*100, wgt);
    h2_eTruth_pxpz->Fill((e4Vect.Px()-e4VectTruth.Px())/(e4VectTruth.Px())*100, (e4Vect.Pz()-e4VectTruth.Pz())/(e4VectTruth.Pz())*100, wgt);
    h2_nTruth_pxpz->Fill((nRec4Vect.Px()-n4VectTruth.Px())/(n4VectTruth.Px())*100, (nRec4Vect.Pz()-n4VectTruth.Pz())/(n4VectTruth.Pz())*100, wgt);

    h2_piTruth_pypz->Fill((pi4Vect.Py()-pi4VectTruth.Py())/(pi4VectTruth.Py())*100, (pi4Vect.Pz()-pi4VectTruth.Pz())/(pi4VectTruth.Pz())*100, wgt);
    h2_eTruth_pypz->Fill((e4Vect.Py()-e4VectTruth.Py())/(e4VectTruth.Py())*100, (e4Vect.Pz()-e4VectTruth.Pz())/(e4VectTruth.Pz())*100, wgt);
    h2_nTruth_pypz->Fill((nRec4Vect.Py()-n4VectTruth.Py())/(n4VectTruth.Py())*100, (nRec4Vect.Pz()-n4VectTruth.Pz())/(n4VectTruth.Pz())*100, wgt);

    h2_piTrack_pTheta_Truth->Fill((pi4VectTruth.Theta()*TMath::RadToDeg()), pi4VectTruth.P(), wgt);
    h2_eTrack_pTheta_Truth->Fill((e4VectTruth.Theta()*TMath::RadToDeg()), e4VectTruth.P(), wgt);
    h2_nTrack_pTheta_Truth->Fill((n4VectTruth.Theta()*TMath::RadToDeg()), n4VectTruth.P(), wgt); 

  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP10on100::ResetEvent(PHCompositeNode *topNode)
{
//  std::cout << "ECCE_DEMP10on100::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
//
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP10on100::EndRun(const int runnumber)
{
  std::cout << "ECCE_DEMP10on100::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP10on100::End(PHCompositeNode *topNode)
{
  std::cout << "ECCE_DEMP10on100::End(PHCompositeNode *topNode) This is the End..." << std::endl;

  std::cout << "thrown events: " << thrownEvents << std::endl; // mfek 06/22/2022
  std::cout << "events after cut 1: " << cut1Events << std::endl; // mfek 06/22/2022
  std::cout << "events after cut 2: " << cut2Events << std::endl; // mfek 06/22/2022
  std::cout << "events after cut 3: " << cut3Events << std::endl; // mfek 06/22/2022  
  std::cout << "events after cut 4: " << cut4Events << std::endl; // mfek 06/22/2022
  std::cout << "events after cut 5: " << cut5Events << std::endl; // mfek 06/22/2022
  std::cout << "events Check_n(): " << count_aftern << std::endl; // mfek 06/22/2022
  std::cout << "events Check_ePi(): " << count_afterePi << std::endl; // mfek 06/22/2022
  std::cout << "nTried: " << nTried << std::endl; // mfek 05/27/2022

  ScalingFact = double(event_itt)/nTried; // This scaling factor is needed to normalise the weighted results
  h1_Q2_DetEff->Divide(h1_Q2_DetEff_Cut, h1_Q2_DetEff_Uncut);
  h2_Q2_t_DetEff->Divide(h2_Q2_t_DetEff_Cut, h2_Q2_t_DetEff_Uncut);
  h2_Q2_t_DetEff_v2->Divide(h2_Q2_t_DetEff_v2_Cut, h2_Q2_t_DetEff_v2_Uncut);
  h2_Q2_W_result->Scale((1/ScalingFact));
  h2_pi_XY->Scale((1/ScalingFact));  
  h2_e_XY->Scale((1/ScalingFact));
  if ( IP_design == "IP6"){
    delete h2_ZDC_XY_IP8;
    h2_ZDC_XY_IP6->Scale((1/ScalingFact));
  }
  else if ( IP_design == "IP8"){
    delete h2_ZDC_XY_IP6;
    h2_ZDC_XY_IP8->Scale((1/ScalingFact));
  }
  else{
    delete h2_ZDC_XY_IP8;
    h2_ZDC_XY_IP6->Scale((1/ScalingFact));
  }
  h2_ZDC_XY_l->Scale((1/ScalingFact));
  h2_t_ttruth_result->Scale((1/ScalingFact));
  h2_t_alt_ttruth_result->Scale((1/ScalingFact));
  h1_Mmiss_result->Scale((1/ScalingFact));
  h1_MmissSq_result->Scale((1/ScalingFact));
  h1_Mmiss_truth_result->Scale((1/ScalingFact));
  h1_Mmiss_Comp_result->Scale((1/ScalingFact));
  h1_taltres_result->Scale((1/ScalingFact));
  for(Int_t C = 0; C < 10; C++){
    h1_taltres_result_ttruth[C]->Scale((1/ScalingFact));
  }
  h2_t_t_alt_result->Scale((1/ScalingFact));
  h2_piTrack_pTheta_Truth->Scale((1/ScalingFact));
  h2_eTrack_pTheta_Truth->Scale((1/ScalingFact));
  h2_nTrack_pTheta_Truth->Scale((1/ScalingFact));

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
  h1_pmissDiff_p->Scale((1/ScalingFact));
  h1_pmissDiff_px->Scale((1/ScalingFact));
  h1_pmissDiff_py->Scale((1/ScalingFact));
  h1_pmissDiff_pz->Scale((1/ScalingFact));

  for(Int_t C = 0; C < 14; C++){
    h1_pmiss_result[C]->Scale((1/ScalingFact));
    h1_pn_result[C]->Scale((1/ScalingFact));	
    h1_t_result[C]->Scale((1/ScalingFact));
    h1_t_truth_thrown_result[C]->Scale((1/ScalingFact));
    h1_t_cut_result[C]->Scale((1/ScalingFact));
    h1_Q2_cut_result[C]->Scale((1/ScalingFact));
    h1_W_cut_result[C]->Scale((1/ScalingFact));
    h2_Q2_t_result[C]->Scale((1/ScalingFact));
    h2_t_ttruth_result_Q2[C]->Scale((1/ScalingFact));
    h2_t_alt_ttruth_result_Q2[C]->Scale((1/ScalingFact));
    h1_t_Resolution[C]->Scale((1/ScalingFact));
    h1_talt_Resolution_ZDC[C]->Scale((1/ScalingFact));
    h1_talt_Resolution_pMiss[C]->Scale((1/ScalingFact));
  } // mfek 06/22/2022 - changed binning

  h1_nTheta_tCut->Scale((1/ScalingFact));
  h1_t_cut1_Low->Scale((1/ScalingFact));
  h1_t_cut2_Low->Scale((1/ScalingFact));
  h1_t_cut3_Low->Scale((1/ScalingFact));
  h1_t_cut4_Low->Scale((1/ScalingFact));
  h1_t_cut1_Mid->Scale((1/ScalingFact));
  h1_t_cut2_Mid->Scale((1/ScalingFact));
  h1_t_cut3_Mid->Scale((1/ScalingFact));
  h1_t_cut4_Mid->Scale((1/ScalingFact));
  h1_t_cut1_High->Scale((1/ScalingFact));
  h1_t_cut2_High->Scale((1/ScalingFact));
  h1_t_cut3_High->Scale((1/ScalingFact));
  h1_t_cut4_High->Scale((1/ScalingFact));

  //ZDC_xyE3D->Scale((1/ScalingFact));

  /*c = new TCanvas();
  c->Divide(2,2);
  c->cd(1);
  nTruth_xyE3D->Project3DProfile("yx")->Draw("COLZ");
  c->cd(2);
  ZDC_xyE3D->Project3DProfile("yx")->Draw("COLZ");
  c->cd(3);
  nTruth_Missed_xyE3D->Project3DProfile("yx")->Draw("COLZ");*/

  outfile->cd();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP10on100::Reset(PHCompositeNode *topNode)
{
 std::cout << "ECCE_DEMP10on100::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void ECCE_DEMP10on100::Print(const std::string &what) const
{
  std::cout << "ECCE_DEMP10on100::Print(const std::string &what) const Printing info for " << what << std::endl;
}


///*****************************************************
/// ZDC Energy and Poisition smearing functions
//
// Energy smearing

float ECCE_DEMP10on100::ZDC_Energy_Smear_EMCAL(float E) {

  float resolution, E_reco;

  resolution = sqrt(.45*.45/E + 0.075*0.075);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}

// Energy smearing

float ECCE_DEMP10on100::ZDC_Energy_Smear_HCAL(float E) {

  float resolution, E_reco;

  //resolution = sqrt(.5*.5/E + 0.1*0.1); // YR Resolution
  resolution = sqrt(.45*.45/E + 0.042*0.042); // Updated Resolution
  //resolution = 0.25/sqrt(E); // Test 25% over sqrt E resolution
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}

// Energy smearing

float ECCE_DEMP10on100::ZDC_Energy_Smear_PbWO4(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

// Position smearing

float ECCE_DEMP10on100::ZDC_Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.15;         /// Position resolution 0.15 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}

//***************************************************

bool ECCE_DEMP10on100::Check_ePi(PHCompositeNode* topNode)
{
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap)
    {
      trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
      if (!trackmap)
    	{
    	  cout
    	    << "ECCE_DEMP10on100::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
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

bool ECCE_DEMP10on100::Check_n(PHCompositeNode* topNode)
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

float ECCE_DEMP10on100::Get_Local_X(float global_x, float global_y, float global_z, float det_tilt, float det_rot) {


   TVector3 global_cor(global_x, global_y, global_z);
   float local_x;

   global_cor.RotateY(-det_rot);
   local_x = global_cor.X()/cos(det_tilt - det_rot);
	
   return local_x;

}

//*******************************************

float ECCE_DEMP10on100::Get_Local_Y(float global_x, float global_y, float global_z, float det_tilt, float cross_angle) {

	return global_y;

}

//*******************************************
float ECCE_DEMP10on100::Get_Local_X(float global_x, float global_y, float global_z, PdbParameterMapContainer *det_nodeparams) {

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

float ECCE_DEMP10on100::Get_Local_X(float global_x, float global_y, float global_z, PHParameters Det_params) {

   float det_xCent = Det_params.get_double_param("place_x");
   float det_zCent = Det_params.get_double_param("place_z");

   float det_tilt = Det_params.get_double_param("rot_y"); // in Rad

   float det_rot = atan( det_xCent / det_zCent);  // in Rad

   TVector3 global_cor(global_x, global_y, global_z);


   float local_x1 = Get_Local_X(global_x, global_y, global_z, det_tilt, det_rot);

   return local_x1;

}

