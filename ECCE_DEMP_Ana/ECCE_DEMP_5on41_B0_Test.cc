#include "ECCE_DEMP_5on41_B0_Test.h"

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

ECCE_DEMP_5on41_B0_Test::ECCE_DEMP_5on41_B0_Test(const std::string &name, const std::string& filename):
 SubsysReco(name)
 , outfilename(filename)
{
  std::cout << "ECCE_DEMP_5on41_B0_Test_example::Diff_Tagg_example(const std::string &name) Calling ctor" << std::endl;

  unsigned int seed = PHRandomSeed();  // fixed seed is handled in this funtcion
  m_RandomGenerator = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(m_RandomGenerator, seed);

}

//____________________________________________________________________________..
ECCE_DEMP_5on41_B0_Test::~ECCE_DEMP_5on41_B0_Test()
{

  gsl_rng_free(m_RandomGenerator);

  std::cout << "ECCE_DEMP_5on41_B0_Test::~ECCE_DEMP_5on41_B0_Test() Calling dtor" << std::endl;
}

//____________________________________________________________________________..
int ECCE_DEMP_5on41_B0_Test::Init(PHCompositeNode *topNode)
{

  static_event_counter = 0;
  hm = new Fun4AllHistoManager(Name());
  outfile = new TFile(outfilename.c_str(), "RECREATE");

  std::cout << "ECCE_DEMP_5on41_B0_Test::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  event_itt = 0;
  
  // Define beam 4 vectors - Assume IP6 by default, for other IPs, adjust in the Process Event loop (at the top)
  e_beam_energy = 5;
  e_beam_pmag = sqrt(pow(e_beam_energy,2)-pow(mElec,2));
  ion_beam_energy = 41;
  ion_beam_pmag = sqrt((pow(ion_beam_energy,2)-pow(mProt,2)));
  crossing_angle = 0.025; 
  eBeam4Vect.SetPxPyPzE(0,0,-1*e_beam_pmag,e_beam_energy);
  pBeam4Vect.SetPxPyPzE(-ion_beam_pmag*TMath::Sin(crossing_angle),0,ion_beam_pmag*TMath::Cos(crossing_angle),ion_beam_energy);

  // Set cut values for physics analysis
  Thetan_Cent = 1.45; // Cut will be +/- 0.4 from this value
  ThetaDiff_Cut = 0.6;
  PhiDiff_Cut = 3.0;
  nTried = 170000; // This is the approximate total number of events generated in a single file of this sample (5on41 epi) 

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP_5on41_B0_Test::InitRun(PHCompositeNode *topNode)
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
    // SJDK - 15/02/22 - This analysis doesn't really use any of these detectors, so beyond getting them to determine the IP, we don't really care
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
int ECCE_DEMP_5on41_B0_Test::process_event(PHCompositeNode *topNode)
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

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP_5on41_B0_Test::ResetEvent(PHCompositeNode *topNode)
{
//  std::cout << "ECCE_DEMP_5on41_B0_Test::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
//
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP_5on41_B0_Test::EndRun(const int runnumber)
{
  std::cout << "ECCE_DEMP_5on41_B0_Test::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP_5on41_B0_Test::End(PHCompositeNode *topNode)
{
  std::cout << "ECCE_DEMP_5on41_B0_Test::End(PHCompositeNode *topNode) This is the End..." << std::endl;
  ScalingFact = double(event_itt)/nTried; // This scaling factor is needed to normalise the weighted results

  outfile->cd();
  outfile->Write();
  outfile->Close();
  delete outfile;
  hm->dumpHistos(outfilename, "UPDATE");

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int ECCE_DEMP_5on41_B0_Test::Reset(PHCompositeNode *topNode)
{
 std::cout << "ECCE_DEMP_5on41_B0_Test::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void ECCE_DEMP_5on41_B0_Test::Print(const std::string &what) const
{
  std::cout << "ECCE_DEMP_5on41_B0_Test::Print(const std::string &what) const Printing info for " << what << std::endl;
}


///*****************************************************
/// ZDC Energy and Poisition smearing functions
//
// Energy smearing

float ECCE_DEMP_5on41_B0_Test::ZDC_Energy_Smear_EMCAL(float E) {

  float resolution, E_reco;

  resolution = sqrt(.45*.45/E + 0.075*0.075);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}

// Energy smearing

float ECCE_DEMP_5on41_B0_Test::ZDC_Energy_Smear_HCAL(float E) {

  float resolution, E_reco;

  //resolution = sqrt(.5*.5/E + 0.1*0.1); // YR Resolution
  resolution = sqrt(.45*.45/E + 0.042*0.042); // Updated Resolution
  //resolution = 0.25/sqrt(E); // Test 25% over sqrt E resolution
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;
}

// Energy smearing

float ECCE_DEMP_5on41_B0_Test::ZDC_Energy_Smear_PbWO4(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

// Position smearing

float ECCE_DEMP_5on41_B0_Test::ZDC_Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.15;         /// Position resolution 0.15 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}

///*****************************************************
/// B0 tracker smearing functions

// Energy smearing

float diff_tagg_ana::B0Tracker_Energy_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

// Posision smearing

float diff_tagg_ana::B0Tracker_Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.1;         /// Position resolution 0.1 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}


///*****************************************************
/// B0 Cal smearing functions

// Energy smearing

float diff_tagg_ana::B0Cal_Energy_Smear(float E) {

  float resolution, E_reco;

  resolution = sqrt(.25*.25/E + 0.04*0.04);
  E_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * E;

  return E_reco;

}

// Posision smearing

float diff_tagg_ana::B0Cal_Position_Smear(float P) {

  float resolution, P_reco;

  resolution = 0.1;         /// Position resolution 0.1 cm
  P_reco = (1+ gsl_ran_gaussian(m_RandomGenerator, resolution)) * P;

  return P_reco;

}

//***************************************************

bool ECCE_DEMP_5on41_B0_Test::Check_ePi(PHCompositeNode* topNode)
{
  SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
  if (!trackmap)
    {
      trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
      if (!trackmap)
    	{
    	  cout
    	    << "ECCE_DEMP_5on41_B0_Test::process_event - Error can not find DST trackmap node SvtxTrackMap" << endl;
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

bool ECCE_DEMP_5on41_B0_Test::Check_n(PHCompositeNode* topNode)
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
	if (ZDC_Energy_Smear_HCAL(hit_iter->second->get_edep()) > 5){ // Hit in ZDC with roughly correct energy for neutron, use smeared energy
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

float ECCE_DEMP_5on41_B0_Test::Get_Local_X(float global_x, float global_y, float global_z, float det_tilt, float det_rot) {


   TVector3 global_cor(global_x, global_y, global_z);
   float local_x;

   global_cor.RotateY(-det_rot);
   local_x = global_cor.X()/cos(det_tilt - det_rot);
	
   return local_x;

}

//*******************************************

float ECCE_DEMP_5on41_B0_Test::Get_Local_Y(float global_x, float global_y, float global_z, float det_tilt, float cross_angle) {

	return global_y;

}

//*******************************************
float ECCE_DEMP_5on41_B0_Test::Get_Local_X(float global_x, float global_y, float global_z, PdbParameterMapContainer *det_nodeparams) {

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

float ECCE_DEMP_5on41_B0_Test::Get_Local_X(float global_x, float global_y, float global_z, PHParameters Det_params) {

   float det_xCent = Det_params.get_double_param("place_x");
   float det_zCent = Det_params.get_double_param("place_z");

   float det_tilt = Det_params.get_double_param("rot_y"); // in Rad

   float det_rot = atan( det_xCent / det_zCent);  // in Rad

   TVector3 global_cor(global_x, global_y, global_z);


   float local_x1 = Get_Local_X(global_x, global_y, global_z, det_tilt, det_rot);

   return local_x1;

}

