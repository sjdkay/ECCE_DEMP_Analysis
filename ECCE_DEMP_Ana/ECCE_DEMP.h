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

  int process_g4hits_ZDC(PHCompositeNode *);
  int process_g4hits_EEMC(PHCompositeNode *);
  int process_g4hits_CEMC(PHCompositeNode *);
  int process_g4hits_FEMC(PHCompositeNode *);
  int process_g4hits_HCALIN(PHCompositeNode *);
  int process_g4hits_HCALOUT(PHCompositeNode *);
  int process_g4hits_FHCAL(PHCompositeNode *);

  int process_g4hits(PHCompositeNode *, const std::string&);
  int process_g4clusters(PHCompositeNode *, const std::string&);

  int process_g4tracks(PHCompositeNode *);

  void use_initial_vertex(const bool b = true) {initial_vertex = b;}

  //private:

 protected:

  //! flag to use initial vertex in track evaluator 
  bool initial_vertex = false;

  std::string detector;
  std::string outfilename;
  Fun4AllHistoManager *hm;

  TFile *outfile;
  TNtuple *g4hitntuple;
  TNtuple *g4trackntuple;
  TNtuple *ZDChitntuple;
  TNtuple *EEMChitntuple;
  TNtuple *CEMChitntuple;
  TNtuple *FEMChitntuple;
  TNtuple *HCALINhitntuple;
  TNtuple *HCALOUThitntuple;
  TNtuple *FHCALhitntuple;
  TNtuple *EEMCalclusterntuple;
  TNtuple *CEMCalclusterntuple;
  TNtuple *FEMCalclusterntuple;
  TNtuple *HCalInclusterntuple;
  TNtuple *HCalOutclusterntuple;
  TNtuple *FHCalclusterntuple;
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
  double ion_beam_energy;

  double crossing_angle;

  TLorentzVector r_lelectron;
//  TLorentzVector r_lproton;

  TLorentzVector r_lscatelec;
  TLorentzVector r_l_scat_nucleon;

  TLorentzVector lproton;

  Int_t ZDC_hit;
  Int_t EEMC_hit;
  Int_t CEMC_hit;
  Int_t FEMC_hit;
  Int_t HCALIN_hit;
  Int_t HCALOUT_hit;
  Int_t FHCAL_hit;

  TH2F* h2_ZDC_XY; 
  TH2F* h2_ZDC_XY_nEnergy;
  TH2F* h2_ZDC_XY_nEnergy_Smeared;

  TH1F* h1_ZDC_E_dep;
  TH1F* h1_ZDC_E_dep_smeared;

  TH1F* h1_EMCal_E_dep;
  TH1F* h1_EMCal_E_dep_smeared;

  TH1F* h1_HCalIn_E_dep;
  TH1F* h1_HCalIn_E_dep_smeared;

  TH1F* h1_HCalOut_E_dep;
  TH1F* h1_HCalOut_E_dep_smeared;

  TH1F* h1_eTrack_px;
  TH1F* h1_eTrack_py;
  TH1F* h1_eTrack_pz;
  TH1F* h1_eTrack_p;
  TH1F* h1_eTrack_theta;
  TH1F* h1_eTrack_phi;

  TH1F* h1_piTrack_px;
  TH1F* h1_piTrack_py;
  TH1F* h1_piTrack_pz;
  TH1F* h1_piTrack_p;
  TH1F* h1_piTrack_theta;
  TH1F* h1_piTrack_phi;

  TH1F* h1_nTracksDist;

  TH2F* h2_ePiTrackDist;

  TH2F* h2_eTrack_ThetaPhi;
  TH2F* h2_eTrack_pTheta;

  TH2F* h2_piTrack_ThetaPhi;
  TH2F* h2_piTrack_pTheta;

};

#endif // ECCE_DEMP_ANA_H
