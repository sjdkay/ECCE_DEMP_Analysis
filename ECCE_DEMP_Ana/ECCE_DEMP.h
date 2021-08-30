// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ECCE_DEMP_ANA_H
#define ECCE_DEMP_ANA_H

#include <fun4all/SubsysReco.h>

#include <string>
#include <memory>
#include <string>
#include <utility>  // std::pair, std::make_pair                                                                                                                                                            
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"

class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class JetEvalStack;

class ECCE_DEMP : public SubsysReco
{
 public:

  ECCE_DEMP(const std::string &name = "ECCE_DEMP", const std::string &fname = "MyNtuple.root");

  virtual ~ECCE_DEMP();

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  bool Check_ePi(PHCompositeNode *);
  bool Check_n(PHCompositeNode *);

  void use_initial_vertex(const bool b = true) {initial_vertex = b;}

  //private:

 protected:

  //! flag to use initial vertex in track evaluator 
  bool initial_vertex = false;

  std::string detector;
  std::string outfilename;
  Fun4AllHistoManager *hm;

  TFile *outfile;
  unsigned long long int event_itt;
  gsl_rng* m_RandomGenerator;

  //*********************************
  // Energy and Position smearing

  float EMCAL_Smear(float E);
  float HCAL_Smear(float E);
  float PbWO4_Smear(float E);
  float Position_Smear(float E);

  //---------------------
  // From ejana

  double true_q2;
  double true_x;
  double true_s_e;
  double true_xpi;
  double true_ypi;
  double true_tpi;

  double have_true_dis_info = false;
  
  bool  HIT_IN_ZDC;
  bool  HIT_IN_EMCAL;
  bool  HIT_IN_HCAL;
  bool  CLUS_IN_EMCAL;
  bool  CLUS_IN_HCAL;
  bool  HIT_IN_HEC;	

  double e_beam_energy;
  double e_beam_pmag;
  double ion_beam_energy;
  double ion_beam_pmag;

  double crossing_angle;

  double wgt;

  TLorentzVector r_lelectron;
//  TLorentzVector r_lproton;

  TLorentzVector r_lscatelec;
  TLorentzVector r_l_scat_nucleon;

  TLorentzVector lproton;

  // Particle Masses
  Double_t mPi = 0.139570;
  Double_t mElec = 0.000510998950;
  Double_t mNeut = 0.93965420;
  Double_t mProt = 0.93827208816;

  // Quantities we want to determine
  TVector3 eVect;
  TVector3 eVectSmeared;
  TVector3 piVect;
  TVector3 piVectSmeared;
  TVector3 nZDCPos;
  TVector3 nZDCPosSmeared;
  TLorentzVector e4Vect;
  TLorentzVector e4VectSmeared;
  TLorentzVector pi4Vect;
  TLorentzVector pi4VectSmeared;
  TLorentzVector n4Vect;
  TLorentzVector n4VectSmeared;
  TLorentzVector e4VectTruth;
  TLorentzVector pi4VectTruth;
  TLorentzVector n4VectTruth;
  TLorentzVector eBeam4Vect;
  TLorentzVector pBeam4Vect;
  TLorentzVector virtphoton4Vect;
  TLorentzVector t4Vect;
  TLorentzVector t_alt4Vect;
  TLorentzVector pmiss4Vect;
  TLorentzVector virtphoton4VectTruth;
  TLorentzVector t4VectTruth;
  TLorentzVector pmiss4VectTruth;

  Double_t nEDep;
  Double_t nEDepSmeared;
  Double_t nTheta;
  Double_t nThetaSmeared;
  Double_t nPhi;
  Double_t nPhiSmeared;
  Double_t nPMag;
  Double_t nPMagSmeared;

  Double_t Q2;
  Double_t W;
  Double_t t;
  Double_t t_alt;
  Double_t xb;
  Double_t xi;
  Double_t Q2_truth;
  Double_t W_truth;
  Double_t t_truth;
  Double_t xb_truth;
  Double_t xi_truth;

  Int_t ZDC_hit;
  Int_t EEMC_hit;

  // Histogram for coincidence analysis routine

  // 1D distributions for each particle
  TH1F* h1_pi_px;
  TH1F* h1_pi_py;
  TH1F* h1_pi_pz;
  TH1F* h1_pi_p;
  TH1F* h1_pi_E;
  TH1F* h1_pi_Theta;
  TH1F* h1_pi_Phi;
  TH1F* h1_e_px;
  TH1F* h1_e_py;
  TH1F* h1_e_pz;
  TH1F* h1_e_p;
  TH1F* h1_e_E;
  TH1F* h1_e_Theta;
  TH1F* h1_e_Phi;
  TH1F* h1_n_px;
  TH1F* h1_n_py;
  TH1F* h1_n_pz;
  TH1F* h1_n_p;
  TH1F* h1_n_E;
  TH1F* h1_n_Theta;
  TH1F* h1_n_Phi;
  TH1F* h1_pmiss_px;
  TH1F* h1_pmiss_py;
  TH1F* h1_pmiss_pz;
  TH1F* h1_pmiss_p;
  TH1F* h1_pmiss_E;
  TH1F* h1_pmiss_Theta;
  TH1F* h1_pmiss_Phi;
  TH1F* h1_gamma_px;
  TH1F* h1_gamma_py;
  TH1F* h1_gamma_pz;
  TH1F* h1_gamma_p;
  TH1F* h1_gamma_E;
  TH1F* h1_gamma_Theta;
  TH1F* h1_gamma_Phi;
  
  TH1F* h1_Q2_Dist;
  TH1F* h1_W_Dist;
  TH1F* h1_t_Dist;
  TH1F* h1_t_alt_Dist;
  TH1F* h1_t_comp;
  TH1F* h1_xb_Dist;
  TH1F* h1_xi_Dist;

  TH1F* h1_Q2Truth_Dist;
  TH1F* h1_WTruth_Dist;
  TH1F* h1_tTruth_Dist;
  TH1F* h1_xbTruth_Dist;
  TH1F* h1_xiTruth_Dist;

  // Resolution test plots for unsmeared vectors
  TH1F* h1_piTruth_p;
  TH1F* h1_piTruth_px;
  TH1F* h1_piTruth_py;
  TH1F* h1_piTruth_pz;
  TH1F* h1_piTruth_E;
  TH1F* h1_eTruth_p;
  TH1F* h1_eTruth_px;
  TH1F* h1_eTruth_py;
  TH1F* h1_eTruth_pz;
  TH1F* h1_eTruth_E;
  TH1F* h1_nTruth_p;
  TH1F* h1_nTruth_px;
  TH1F* h1_nTruth_py;
  TH1F* h1_nTruth_pz;
  TH1F* h1_nTruth_E;
  // Resolution test plots with smeared vectors
  TH1F* h1_piTruth_p_Smeared;
  TH1F* h1_piTruth_px_Smeared;
  TH1F* h1_piTruth_py_Smeared;
  TH1F* h1_piTruth_pz_Smeared;
  TH1F* h1_piTruth_E_Smeared;
  TH1F* h1_eTruth_p_Smeared;
  TH1F* h1_eTruth_px_Smeared;
  TH1F* h1_eTruth_py_Smeared;
  TH1F* h1_eTruth_pz_Smeared;
  TH1F* h1_eTruth_E_Smeared;
  TH1F* h1_nTruth_p_Smeared;
  TH1F* h1_nTruth_px_Smeared;
  TH1F* h1_nTruth_py_Smeared;
  TH1F* h1_nTruth_pz_Smeared;
  TH1F* h1_nTruth_E_Smeared;

  // 2D distributions 
  TH2F* h2_ZDC_XY;
  TH2F* h2_ZDC_XY_Smeared;
  // Particle Theta/Phi and Theta/p distributions
  TH2F* h2_eTrack_ThetaPhi;
  TH2F* h2_eTrack_pTheta;
  TH2F* h2_piTrack_ThetaPhi;
  TH2F* h2_piTrack_pTheta;
  TH2F* h2_nTrack_ThetaPhi;
  TH2F* h2_nTrack_pTheta;
  TH2F* h2_eTrack_ThetaPhi_Smeared;
  TH2F* h2_eTrack_pTheta_Smeared;
  TH2F* h2_piTrack_ThetaPhi_Smeared;
  TH2F* h2_piTrack_pTheta_Smeared;
  TH2F* h2_nTrack_ThetaPhi_Smeared;
  TH2F* h2_nTrack_pTheta_Smeared;
  // 2D resolution test plots
  TH2F* h2_eTruth_pxpy;
  TH2F* h2_piTruth_pxpy;
  TH2F* h2_nTruth_pxpy;
  TH2F* h2_eTruth_pxpz;
  TH2F* h2_piTruth_pxpz;
  TH2F* h2_nTruth_pxpz;
  TH2F* h2_eTruth_pypz;
  TH2F* h2_piTruth_pypz;
  TH2F* h2_nTruth_pypz;
  TH2F* h2_eTruth_pxpy_Smeared;
  TH2F* h2_piTruth_pxpy_Smeared;
  TH2F* h2_nTruth_pxpy_Smeared;
  
  // 1D Kinematic analysis plots
  TH1F* h1_t_Q2[7];
  TH1F* h1_t_alt_Q2[7];

  // 2D Kinematic analysis plots
  TH2F* h2_t_ep;
  TH2F* h2_t_Q2;
  TH2F* h2_delta_t_t;
  TH2F* h2_delta_t_t_Q2[7];

  // 1D distributions for each particle
  TH1F* h1_pi_px_Weighted;
  TH1F* h1_pi_py_Weighted;
  TH1F* h1_pi_pz_Weighted;
  TH1F* h1_pi_p_Weighted;
  TH1F* h1_pi_E_Weighted;
  TH1F* h1_pi_Theta_Weighted;
  TH1F* h1_pi_Phi_Weighted;
  TH1F* h1_e_px_Weighted;
  TH1F* h1_e_py_Weighted;
  TH1F* h1_e_pz_Weighted;
  TH1F* h1_e_p_Weighted;
  TH1F* h1_e_E_Weighted;
  TH1F* h1_e_Theta_Weighted;
  TH1F* h1_e_Phi_Weighted;
  TH1F* h1_n_px_Weighted;
  TH1F* h1_n_py_Weighted;
  TH1F* h1_n_pz_Weighted;
  TH1F* h1_n_p_Weighted;
  TH1F* h1_n_E_Weighted;
  TH1F* h1_n_Theta_Weighted;
  TH1F* h1_n_Phi_Weighted;
  TH1F* h1_pmiss_px_Weighted;
  TH1F* h1_pmiss_py_Weighted;
  TH1F* h1_pmiss_pz_Weighted;
  TH1F* h1_pmiss_p_Weighted;
  TH1F* h1_pmiss_E_Weighted;
  TH1F* h1_pmiss_Theta_Weighted;
  TH1F* h1_pmiss_Phi_Weighted;
  TH1F* h1_gamma_px_Weighted;
  TH1F* h1_gamma_py_Weighted;
  TH1F* h1_gamma_pz_Weighted;
  TH1F* h1_gamma_p_Weighted;
  TH1F* h1_gamma_E_Weighted;
  TH1F* h1_gamma_Theta_Weighted;
  TH1F* h1_gamma_Phi_Weighted;
  
  TH1F* h1_Q2_Dist_Weighted;
  TH1F* h1_W_Dist_Weighted;
  TH1F* h1_t_Dist_Weighted;
  TH1F* h1_t_alt_Dist_Weighted;
  TH1F* h1_t_comp_Weighted;
  TH1F* h1_xb_Dist_Weighted;
  TH1F* h1_xi_Dist_Weighted;

  TH1F* h1_Q2Truth_Dist_Weighted;
  TH1F* h1_WTruth_Dist_Weighted;
  TH1F* h1_tTruth_Dist_Weighted;
  TH1F* h1_xbTruth_Dist_Weighted;
  TH1F* h1_xiTruth_Dist_Weighted;

  // Resolution test plots for unsmeared vectors
  TH1F* h1_piTruth_p_Weighted;
  TH1F* h1_piTruth_px_Weighted;
  TH1F* h1_piTruth_py_Weighted;
  TH1F* h1_piTruth_pz_Weighted;
  TH1F* h1_piTruth_E_Weighted;
  TH1F* h1_eTruth_p_Weighted;
  TH1F* h1_eTruth_px_Weighted;
  TH1F* h1_eTruth_py_Weighted;
  TH1F* h1_eTruth_pz_Weighted;
  TH1F* h1_eTruth_E_Weighted;
  TH1F* h1_nTruth_p_Weighted;
  TH1F* h1_nTruth_px_Weighted;
  TH1F* h1_nTruth_py_Weighted;
  TH1F* h1_nTruth_pz_Weighted;
  TH1F* h1_nTruth_E_Weighted;
  // Resolution test plots with smeared vectors
  TH1F* h1_piTruth_p_Smeared_Weighted;
  TH1F* h1_piTruth_px_Smeared_Weighted;
  TH1F* h1_piTruth_py_Smeared_Weighted;
  TH1F* h1_piTruth_pz_Smeared_Weighted;
  TH1F* h1_piTruth_E_Smeared_Weighted;
  TH1F* h1_eTruth_p_Smeared_Weighted;
  TH1F* h1_eTruth_px_Smeared_Weighted;
  TH1F* h1_eTruth_py_Smeared_Weighted;
  TH1F* h1_eTruth_pz_Smeared_Weighted;
  TH1F* h1_eTruth_E_Smeared_Weighted;
  TH1F* h1_nTruth_p_Smeared_Weighted;
  TH1F* h1_nTruth_px_Smeared_Weighted;
  TH1F* h1_nTruth_py_Smeared_Weighted;
  TH1F* h1_nTruth_pz_Smeared_Weighted;
  TH1F* h1_nTruth_E_Smeared_Weighted;

  // 2D distributions 
  TH2F* h2_ZDC_XY_Weighted;
  TH2F* h2_ZDC_XY_Smeared_Weighted;
  // Particle Theta/Phi and Theta/p distributions
  TH2F* h2_eTrack_ThetaPhi_Weighted;
  TH2F* h2_eTrack_pTheta_Weighted;
  TH2F* h2_piTrack_ThetaPhi_Weighted;
  TH2F* h2_piTrack_pTheta_Weighted;
  TH2F* h2_nTrack_ThetaPhi_Weighted;
  TH2F* h2_nTrack_pTheta_Weighted;
  TH2F* h2_eTrack_ThetaPhi_Smeared_Weighted;
  TH2F* h2_eTrack_pTheta_Smeared_Weighted;
  TH2F* h2_piTrack_ThetaPhi_Smeared_Weighted;
  TH2F* h2_piTrack_pTheta_Smeared_Weighted;
  TH2F* h2_nTrack_ThetaPhi_Smeared_Weighted;
  TH2F* h2_nTrack_pTheta_Smeared_Weighted;
  // 2D resolution test plots
  TH2F* h2_eTruth_pxpy_Weighted;
  TH2F* h2_piTruth_pxpy_Weighted;
  TH2F* h2_nTruth_pxpy_Weighted;
  TH2F* h2_eTruth_pxpz_Weighted;
  TH2F* h2_piTruth_pxpz_Weighted;
  TH2F* h2_nTruth_pxpz_Weighted;
  TH2F* h2_eTruth_pypz_Weighted;
  TH2F* h2_piTruth_pypz_Weighted;
  TH2F* h2_nTruth_pypz_Weighted;
  TH2F* h2_eTruth_pxpy_Smeared_Weighted;
  TH2F* h2_piTruth_pxpy_Smeared_Weighted;
  TH2F* h2_nTruth_pxpy_Smeared_Weighted;
  
  // 1D Kinematic analysis plots
  TH1F* h1_t_Q2_Weighted[7];
  TH1F* h1_t_alt_Q2_Weighted[7];

  // 2D Kinematic analysis plots
  TH2F* h2_t_ep_Weighted;
  TH2F* h2_t_Q2_Weighted;
  TH2F* h2_delta_t_t_Weighted;
  TH2F* h2_delta_t_t_Q2_Weighted[7];
 
};

#endif // ECCE_DEMP_ANA_H
