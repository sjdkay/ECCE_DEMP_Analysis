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

  double nTried;
  double ScalingFact;
  double wgt;

  TLorentzVector r_lelectron;
  //TLorentzVector r_lproton;

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
  TVector3 piVect;
  TVector3 nZDCPos;
  TLorentzVector e4Vect;
  TLorentzVector pi4Vect;
  TLorentzVector n4Vect;
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
  TLorentzVector t_alt4VectTruth;
  TLorentzVector pmiss4VectTruth;
  TLorentzVector pmiss4Vect_2;
  TLorentzVector pmiss4VectTruth_2;

  Double_t nEDep;
  Double_t nTheta;
  Double_t nPhi;
  Double_t nPMag;

  Double_t Q2;
  Double_t W;
  Double_t t;
  Double_t t_alt;
  Double_t xb;
  Double_t xi;
  Double_t Q2_truth;
  Double_t W_truth;
  Double_t t_truth;
  Double_t t_alt_truth;
  Double_t xb_truth;
  Double_t xi_truth;

  Double_t Q2_low;
  Double_t Q2_high;
  Double_t Thetan_Cent; // Central thetan value to determine cuts from
  // The Pmiss cut values are chosen a little arbitrarily, hard to judge w/o SIDIS to compare with
  // Cut will be anything ABOVE this value for each bin
  Double_t PmissCutVal[8] = {96.0, 93.5, 91.0, 87.0, 83.0, 80.0, 77.5, 75.0}; // Array to store Pmiss cut values in, 5on 100
  Int_t ZDC_hit;
  Int_t EEMC_hit;

  // Histograms for coincidence analysis routine
  TH1F* h1_Q2_DetEff_Uncut;
  TH1F* h1_Q2_DetEff_Cut;
  TH1F* h1_Q2_DetEff;
  TH2F* h2_Q2_t_DetEff_Uncut;
  TH2F* h2_Q2_t_DetEff_Cut;
  TH2F* h2_Q2_t_DetEff;

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

  TH2F* h2_pi_XY;
  TH2F* h2_e_XY;
  TH2F* h2_n_XY;
  
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
  TH1F* h1_t_altTruth_Dist;
  TH1F* h1_xbTruth_Dist;
  TH1F* h1_xiTruth_Dist;

  // Particle Truth info
  TH1F* h1_piTruth_p;
  TH1F* h1_piTruth_px;
  TH1F* h1_piTruth_py;
  TH1F* h1_piTruth_pz;
  TH1F* h1_piTruth_E;
  TH1F* h1_piTruth_Theta;
  TH1F* h1_eTruth_p;
  TH1F* h1_eTruth_px;
  TH1F* h1_eTruth_py;
  TH1F* h1_eTruth_pz;
  TH1F* h1_eTruth_E; 
  TH1F* h1_eTruth_Theta;
  TH1F* h1_nTruth_p;
  TH1F* h1_nTruth_px;
  TH1F* h1_nTruth_py;
  TH1F* h1_nTruth_pz;
  TH1F* h1_nTruth_E;
  TH1F* h1_nTruth_Theta;

  // Particle Resolution histograms
  TH1F* h1_piRes_p;
  TH1F* h1_piRes_px;
  TH1F* h1_piRes_py;
  TH1F* h1_piRes_pz;  
  TH1F* h1_eRes_p;
  TH1F* h1_eRes_px;
  TH1F* h1_eRes_py;
  TH1F* h1_eRes_pz;
  TH1F* h1_nRes_p;
  TH1F* h1_nRes_px;
  TH1F* h1_nRes_py;
  TH1F* h1_nRes_pz;
  
  // 2D distributions 
  TH2F* h2_ZDC_XY;
  // Particle Theta/Phi and Theta/p distributions
  TH2F* h2_eTrack_ThetaPhi;
  TH2F* h2_eTrack_pTheta;
  TH2F* h2_piTrack_ThetaPhi;
  TH2F* h2_piTrack_pTheta;
  TH2F* h2_nTrack_ThetaPhi;
  TH2F* h2_nTrack_pTheta;
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
  
  // 1D Kinematic analysis plots
  TH1F* h1_t_Q2[7];
  TH1F* h1_t_alt_Q2[7];

  // 2D Kinematic analysis plots
  TH2F* h2_t_ep;
  TH2F* h2_t_Q2;
  TH2F* h2_delta_t_t;
  TH2F* h2_delta_t_t_Q2[7];

  // 1D Physics results plots
  TH1F* h1_Mmiss_result;
  TH1F* h1_Mmiss_truth_result;
  TH1F* h1_Mmiss_Comp_result;
  TH1F* h1_t_result[8];
  TH1F* h1_nTheta_result[8];
  TH1F* h1_pmiss_result[8];
  TH1F* h1_t_cut_result[8];
  TH1F* h1_Q2_cut_result[8];
  TH1F* h1_W_cut_result[8];

  // 2D Physics Results Plots
  TH2F* h2_Q2_W_result;
  TH2F* h2_t_ttruth_result;
  TH2F* h2_t_alt_ttruth_result;
  TH2F* h2_t_alt_ttruth_result_Q2[8];
  TH2F* h2_t_t_alt_result;
  TH2F* h2_Q2_t_result[8];

};

#endif // ECCE_DEMP_ANA_H
