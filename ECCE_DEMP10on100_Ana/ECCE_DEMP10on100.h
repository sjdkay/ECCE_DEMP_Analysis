// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ECCE_DEMP10on100_ANA_H
#define ECCE_DEMP10on100_ANA_H

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
#include "TH3.h"
#include "TProfile2D.h"
#include "TCanvas.h"

#include <pdbcalbase/PdbParameterMap.h>
#include <phparameter/PHParameters.h>

class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class JetEvalStack;

class ECCE_DEMP10on100 : public SubsysReco
{
 public:

  ECCE_DEMP10on100(const std::string &name = "ECCE_DEMP10on100", const std::string &fname = "MyNtuple.root");

  virtual ~ECCE_DEMP10on100();

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

  int static_event_counter;

  //*********************************
  // ZDC Energy and Position smearing

  float ZDC_Energy_Smear_EMCAL(float E);
  float ZDC_Energy_Smear_HCAL(float E);
  float ZDC_Energy_Smear_PbWO4(float E);
  float ZDC_Position_Smear(float E);

  //*********************************
  // Coordinate transformation from global to local

  float Get_Local_X(float global_x, float global_y, float global_z, float det_tilt, float det_rot);
  float Get_Local_Y(float global_x, float global_y, float global_z, float det_tilt, float det_rot);
  float Get_Local_X(float global_x, float global_y, float global_z, PdbParameterMapContainer *det_nodeparams);
  float Get_Local_X(float global_x, float global_y, float global_z, PHParameters Det_params);

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

  // Particle Masses in GeV
  Double_t mPi = 0.139570;
  Double_t mElec = 0.000510998950;
  Double_t mNeut = 0.93965420;
  Double_t mProt = 0.93827208816;

  // Quantities we want to determine
  TVector3 eVect;
  TVector3 piVect;
  TVector3 nZDCPos;
  TVector3 nRecZDCPos;
  TLorentzVector e4Vect;
  TLorentzVector pi4Vect;
  TLorentzVector n4Vect;
  TLorentzVector nRec4Vect;
  TLorentzVector e4VectTruth;
  TLorentzVector pi4VectTruth;
  TLorentzVector n4VectTruth;
  TLorentzVector eBeam4Vect;
  TLorentzVector pBeam4Vect;
  TLorentzVector virtphoton4Vect;
  TLorentzVector t4Vect;
  TLorentzVector t_alt4Vect;
  TLorentzVector t_alt4Vect_ZDC;
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

  Double_t nRecEDep;
  Double_t nRecTheta;
  Double_t nRecPhi;
  Double_t nRecPMag;

  Double_t nTheta_Diff;
  Double_t nPhi_Diff;

  float det_x_pos;
  float det_y_pos;
  float det_z_pos;

  float local_x;
  float local_y;

  Double_t Q2;
  Double_t W;
  Double_t t;
  Double_t t_alt;
  Double_t t_alt_ZDC;
  Double_t xb;
  Double_t xi;
  Double_t Q2_truth;
  Double_t W_truth;
  Double_t t_truth;
  Double_t t_alt_truth;
  Double_t xb_truth;
  Double_t xi_truth;

  Double_t t_low;
  Double_t t_high;
  Double_t Q2_low;
  Double_t Q2_high;
  Double_t Thetan_Cent; // Central thetan value to determine cuts from
  Double_t ThetaDiff_Cut;
  Double_t PhiDiff_Cut;
  // The Pmiss cut values are chosen a little arbitrarily, hard to judge w/o SIDIS to compare with
  // Cut will be anything ABOVE this value for each bin
  Double_t PmissCutVal[23] = {96.5, 96.5, 94.5, 92.5, 90.5, 88.5, 86.5, 85.5, 83.5, 81.5, 80.0, 77.5, 76.5, 76.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0, 75.0}; // Array to store Pmiss cut values in
  Double_t Q2BinVal[15] = {2.0, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 32.5, 35.0, 40.0}; // mfek 06/22/2022 - new binning
  // mfek 06/22/2022 - changed binning
  Int_t ZDC_hit;
  Int_t EEMC_hit;

  // Histograms for coincidence analysis routine
  TH1F* h1_Q2_DetEff_Uncut;
  TH1F* h1_Q2_DetEff_Cut;
  TH1F* h1_Q2_DetEff;
  TH2F* h2_Q2_t_DetEff_Uncut;
  TH2F* h2_Q2_t_DetEff_Cut;
  TH2F* h2_Q2_t_DetEff;
  TH2F* h2_Q2_t_DetEff_v2_Uncut;
  TH2F* h2_Q2_t_DetEff_v2_Cut;
  TH2F* h2_Q2_t_DetEff_v2;

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
  TH1F* h1_n_ThetaDiff;
  TH1F* h1_n_PhiDiff;
  TH2F* h2_n_ThetaPhiDiff;

  TH1F* h1_nRec_px;
  TH1F* h1_nRec_py;
  TH1F* h1_nRec_pz;
  TH1F* h1_nRec_p;
  TH1F* h1_nRec_E;
  TH1F* h1_nRec_Theta;
  TH1F* h1_nRec_Phi;

  TH1F* h1_pi_px_Unweighted;
  TH1F* h1_pi_py_Unweighted;
  TH1F* h1_pi_pz_Unweighted;
  TH1F* h1_pi_p_Unweighted;
  TH1F* h1_pi_E_Unweighted;
  TH1F* h1_pi_Theta_Unweighted;
  TH1F* h1_pi_Phi_Unweighted;
  TH1F* h1_e_px_Unweighted;
  TH1F* h1_e_py_Unweighted;
  TH1F* h1_e_pz_Unweighted;
  TH1F* h1_e_p_Unweighted;
  TH1F* h1_e_E_Unweighted;
  TH1F* h1_e_Theta_Unweighted;
  TH1F* h1_e_Phi_Unweighted;
  TH1F* h1_n_px_Unweighted;
  TH1F* h1_n_py_Unweighted;
  TH1F* h1_n_pz_Unweighted;
  TH1F* h1_n_p_Unweighted;
  TH1F* h1_n_E_Unweighted;
  TH1F* h1_n_Theta_Unweighted;
  TH1F* h1_n_Phi_Unweighted;
  TH1F* h1_n_ThetaDiff_Unweighted;
  TH1F* h1_n_PhiDiff_Unweighted;
  TH2F* h2_n_ThetaPhiDiff_Unweighted;

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

  // Particle Truth Info Pre-Cut - mfek 06/23/2022
  TH1F* h1_piTruth_preCut_p;
  TH1F* h1_piTruth_preCut_px;
  TH1F* h1_piTruth_preCut_py;
  TH1F* h1_piTruth_preCut_pz;
  TH1F* h1_piTruth_preCut_E;
  TH1F* h1_piTruth_preCut_Theta;
  TH1F* h1_piTruth_preCut_Phi;

  TH1F* h1_eTruth_preCut_p;
  TH1F* h1_eTruth_preCut_px;
  TH1F* h1_eTruth_preCut_py;
  TH1F* h1_eTruth_preCut_pz;
  TH1F* h1_eTruth_preCut_E;
  TH1F* h1_eTruth_preCut_Theta;
  TH1F* h1_eTruth_preCut_Phi;

  TH1F* h1_nTruth_preCut_p;
  TH1F* h1_nTruth_preCut_px;
  TH1F* h1_nTruth_preCut_py;
  TH1F* h1_nTruth_preCut_pz;
  TH1F* h1_nTruth_preCut_E;
  TH1F* h1_nTruth_preCut_Theta;
  TH1F* h1_nTruth_preCut_Theta_inRange;
  TH1F* h1_nTruth_preCut_Phi;
  TH2F* h2_nTruth_preCut_XY;
  TH2F* h2_nTruth_preCut_XY_inZDC;
  TH2F* h2_nTruth_preCut_XY_outZDC;
  TH2F* h2_nTruth_XY_hits;

  // Particle Truth info for Missed Events - mfek 06/23/2022
  TH1F* h1_piTruth_Missed_p;
  TH1F* h1_piTruth_Missed_px;
  TH1F* h1_piTruth_Missed_py;
  TH1F* h1_piTruth_Missed_pz;
  TH1F* h1_piTruth_Missed_E;
  TH1F* h1_piTruth_Missed_Theta;
  TH1F* h1_piTruth_Missed_Phi;

  TH1F* h1_eTruth_Missed_p;
  TH1F* h1_eTruth_Missed_px;
  TH1F* h1_eTruth_Missed_py;
  TH1F* h1_eTruth_Missed_pz;
  TH1F* h1_eTruth_Missed_E;
  TH1F* h1_eTruth_Missed_Theta;
  TH1F* h1_eTruth_Missed_Phi;

  TH1F* h1_nTruth_Missed_p;
  TH1F* h1_nTruth_Missed_px;
  TH1F* h1_nTruth_Missed_py;
  TH1F* h1_nTruth_Missed_pz;
  TH1F* h1_nTruth_Missed_E;
  TH1F* h1_nTruth_Missed_Theta;
  TH1F* h1_nTruth_Missed_Theta_inRange;
  TH1F* h1_nTruth_Missed_Phi;

  TH2F* h2_nTruth_Missed_XY;

  TH1F* h1_Q2Truth_Dist_Missed;
  TH1F* h1_tTruth_Dist_Missed;

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
  TH1F* h1_pmissDiff_p;
  TH1F* h1_pmissDiff_px;
  TH1F* h1_pmissDiff_py;
  TH1F* h1_pmissDiff_pz;

  // 2D distributions 
  TH2F* h2_ZDC_XY_IP6;
  TH2F* h2_ZDC_XY_IP8;
  TH2F* h2_ZDC_XY_l;
  // Particle Theta/Phi and Theta/p distributions
  TH2F* h2_eTrack_ThetaPhi;
  TH2F* h2_eTrack_pTheta;
  TH2F* h2_eTrack_pTheta_Truth;
  TH2F* h2_piTrack_ThetaPhi;
  TH2F* h2_piTrack_pTheta;
  TH2F* h2_piTrack_pTheta_Truth;
  TH2F* h2_nTrack_ThetaPhi;
  TH2F* h2_nTrack_pTheta;
  TH2F* h2_nTrack_pTheta_Truth;
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
  TH1F* h1_t_Q2[14]; // mfek 06/22/2022 - changed binning
  TH1F* h1_t_alt_Q2[14]; // mfek 06/22/2022 - changed binning

  // 2D Kinematic analysis plots
  TH2F* h2_t_ep;
  TH2F* h2_t_Q2;
  TH2F* h2_delta_t_t;
  TH2F* h2_delta_t_t_Q2[14]; // mfek 06/22/2022 - changed binning

  // 1D Physics results plots
  TH1F* h1_Mmiss_result;
  TH1F* h1_MmissSq_result;
  TH1F* h1_Mmiss_truth_result;
  TH1F* h1_Mmiss_Comp_result;
  TH1F* h1_taltres_result;
  TH1F* h1_taltres_result_ttruth[10]; // Binned in t_truth
  TH1F* h1_t_result[14]; // mfek 06/22/2022 - changed binning
  TH1F* h1_t_truth_thrown_result[14]; // mfek 06/22/2022 - changed binning
  TH1F* h1_nTheta_result[14]; // mfek 06/22/2022 - changed binning
  TH1F* h1_pmiss_result[14]; // mfek 06/22/2022 - changed binning
  TH1F* h1_pn_result[14]; // mfek 06/22/2022 - changed binning
  TH1F* h1_t_cut_result[14]; // mfek 06/22/2022 - changed binning
  TH1F* h1_Q2_cut_result[14]; // mfek 06/22/2022 - changed binning
  TH1F* h1_W_cut_result[14]; // mfek 06/22/2022 - changed binning

  // 2D Physics Results Plots
  TH2F* h2_Q2_W_result;
  TH2F* h2_t_ttruth_result;
  TH2F* h2_t_alt_ttruth_result;
  TH2F* h2_t_ttruth_result_Q2[14]; // mfek 06/22/2022 - changed binning
  TH2F* h2_t_alt_ttruth_result_Q2[14]; // mfek 06/22/2022 - changed binning
  TH2F* h2_t_t_alt_result;
  TH2F* h2_Q2_t_result[14]; // mfek 06/22/2022 - changed binning

  // Cut analysis plots
  TH1F* h1_nTheta_tCut; // nTheta dist with just the -t cut
  TH1F* h1_t_cut1_Low; // Just the -t cut
  TH1F* h1_t_cut2_Low; // -t cut and theta n cut
  TH1F* h1_t_cut3_Low; // -t cut, theta n cut, theta/phi diff cuts
  TH1F* h1_t_cut4_Low; // -t cut, theta n cut, theta/phi diff cuts, pmiss cuts
  TH1F* h1_t_cut1_Mid; // Just the -t cut
  TH1F* h1_t_cut2_Mid; // -t cut and theta n cut
  TH1F* h1_t_cut3_Mid; // -t cut, theta n cut, theta/phi diff cuts
  TH1F* h1_t_cut4_Mid; // -t cut, theta n cut, theta/phi diff cuts, pmiss cuts
  TH1F* h1_t_cut1_High; // Just the -t cut
  TH1F* h1_t_cut2_High; // -t cut and theta n cut
  TH1F* h1_t_cut3_High; // -t cut, theta n cut, theta/phi diff cuts
  TH1F* h1_t_cut4_High; // -t cut, theta n cut, theta/phi diff cuts, pmiss cuts

  TH1F* h1_t_Resolution[14]; // mfek 06/22/2022 - changed binning
  TH1F* h1_talt_Resolution_ZDC[14]; // mfek 06/22/2022 - changed binning
  TH1F* h1_talt_Resolution_pMiss[14]; // mfek 06/22/2022 - changed binning

  int thrownEvents; // mfek 06/22/2022 - new counter variable added
  int cut1Events; // mfek 06/22/2022 - new counter variable added 
  int cut2Events; // mfek 06/22/2022 - new counter variable added
  int cut3Events; // mfek 06/22/2022 - new counter variable added
  int cut4Events; // mfek 06/22/2022 - new counter variable added
  int cut5Events; // mfek 06/22/2022 - new counter variable added

  int count_afterePi; // mfek 06/22/2022 - new counter variable added
  int count_aftern; // mfek 06/22/2022 - new counter variable added

  TH3F* nTruth_xyE3D; // mfek 06/23/2022
  TH3F* ZDC_xyE3D; // mfek 06/23/2022
  TH3F* nTruth_Missed_xyE3D; // mfek 06/23/2022

  TProfile2D* h3nTruth_xyE_pxy; // mfek 06/23/2022
  TProfile2D* h3ZDC_xyE_pxy; // mfek 06/23/2022
  TProfile2D* h3nTruth_Missed_xyE_pxy; // mfek 06/23/2022

  TCanvas* c; // mfek 06/23/2022

  PHParameters Enclosure_params{"PHGEnclosure"};
  PHParameters ZDC_params{"PHG4RP"};
  PHParameters RP_1_params{"PHG4RP"};
  PHParameters RP2_params{"PHG4RP2"};
  PHParameters B0_params{"PHG4B0"};
  PHParameters BeamLineMagnet_params{"PHG4BeamLinMagnet"};

  PdbParameterMapContainer *encloseure_nodeparams; 
  PdbParameterMapContainer *zdc_nodeparams; 
  PdbParameterMapContainer *rp_nodeparams;
  PdbParameterMapContainer *rp2_nodeparams;
  PdbParameterMapContainer *b0_nodeparams;
  PdbParameterMapContainer *beamlinemagnet_nodeparams; 

  TString IP_design;

};

#endif // ECCE_DEMP10on100_ANA_H
