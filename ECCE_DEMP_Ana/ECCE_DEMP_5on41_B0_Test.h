// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef ECCE_DEMP_5on41_B0_Test_ANA_H
#define ECCE_DEMP_5on41_B0_Test_ANA_H

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

#include <pdbcalbase/PdbParameterMap.h>
#include <phparameter/PHParameters.h>

class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;
class JetEvalStack;

class ECCE_DEMP_5on41_B0_Test : public SubsysReco
{
 public:

  ECCE_DEMP_5on41_B0_Test(const std::string &name = "ECCE_DEMP_5on41_B0_Test", const std::string &fname = "MyNtuple.root");

  virtual ~ECCE_DEMP_5on41_B0_Test();

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
  int process_g4hits_B0(PHCompositeNode *);

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
  // B0 Tracker Energy and Position smearing

  float B0Tracker_Energy_Smear(float E);
  float B0Tracker_Position_Smear(float E);

  //*********************************
  // B0 Cal Energy and Position smearing

  float B0Cal_Energy_Smear(float E);
  float B0Cal_Position_Smear(float E);

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
  bool  HIT_IN_B0;
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
  Double_t PmissCutVal[8] = {96.0, 93.5, 91.0, 87.0, 83.0, 80.0, 77.5, 75.0}; // Array to store Pmiss cut values in, 5on 100
  Int_t ZDC_hit;
  Int_t EEMC_hit;

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

#endif // ECCE_DEMP_5on41_B0_Test_ANA_H
