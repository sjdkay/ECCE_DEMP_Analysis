//____________________________________________________________________________..
//
// This is a template for a Fun4All SubsysReco module with all methods from the
// $OFFLINE_MAIN/include/fun4all/SubsysReco.h baseclass
// You do not have to implement all of them, you can just remove unused methods
// here and in ECCE_DEMP.h.
//
// ECCE_DEMP(const std::string &name = "ECCE_DEMP")
// everything is keyed to ECCE_DEMP, duplicate names do work but it makes
// e.g. finding culprits in logs difficult or getting a pointer to the module
// from the command line
//
// ECCE_DEMP::~ECCE_DEMP()
// this is called when the Fun4AllServer is deleted at the end of running. Be
// mindful what you delete - you do loose ownership of object you put on the node tree
//
// int ECCE_DEMP::Init(PHCompositeNode *topNode)
// This method is called when the module is registered with the Fun4AllServer. You
// can create historgrams here or put objects on the node tree but be aware that
// modules which haven't been registered yet did not put antyhing on the node tree
//
// int ECCE_DEMP::InitRun(PHCompositeNode *topNode)
// This method is called when the first event is read (or generated). At
// this point the run number is known (which is mainly interesting for raw data
// processing). Also all objects are on the node tree in case your module's action
// depends on what else is around. Last chance to put nodes under the DST Node
// We mix events during readback if branches are added after the first event
//
// int ECCE_DEMP::process_event(PHCompositeNode *topNode)
// called for every event. Return codes trigger actions, you find them in
// $OFFLINE_MAIN/include/fun4all/Fun4AllReturnCodes.h
//   everything is good:
//     return Fun4AllReturnCodes::EVENT_OK
//   abort event reconstruction, clear everything and process next event:
//     return Fun4AllReturnCodes::ABORT_EVENT; 
//   proceed but do not save this event in output (needs output manager setting):
//     return Fun4AllReturnCodes::DISCARD_EVENT; 
//   abort processing:
//     return Fun4AllReturnCodes::ABORT_RUN
// all other integers will lead to an error and abort of processing
//
// int ECCE_DEMP::ResetEvent(PHCompositeNode *topNode)
// If you have internal data structures (arrays, stl containers) which needs clearing
// after each event, this is the place to do that. The nodes under the DST node are cleared
// by the framework
//
// int ECCE_DEMP::EndRun(const int runnumber)
// This method is called at the end of a run when an event from a new run is
// encountered. Useful when analyzing multiple runs (raw data). Also called at
// the end of processing (before the End() method)
//
// int ECCE_DEMP::End(PHCompositeNode *topNode)
// This is called at the end of processing. It needs to be called by the macro
// by Fun4AllServer::End(), so do not forget this in your macro
//
// int ECCE_DEMP::Reset(PHCompositeNode *topNode)
// not really used - it is called before the dtor is called
//
// void ECCE_DEMP::Print(const std::string &what) const
// Called from the command line - useful to print information when you need it
//
//____________________________________________________________________________..

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
  // create and register your histos (all types) here
  // TH1 *h1 = new TH1F("h1",....)
  // hm->registerHisto(h1);
  outfile = new TFile(outfilename.c_str(), "RECREATE");
  g4hitntuple = new TNtuple("hitntup", "G4Hits", "x0:y0:z0:edep");
  g4trackntuple = new TNtuple("trackntup", "G4Tracks", "px:py:pz:p:theta:phi");
  ZDChitntuple = new TNtuple("ZDChitntup", "G4Hits", "x0:y0:z0:edep");
  EEMChitntuple = new TNtuple("EEMChitntup", "G4Hits", "x0:y0:z0:edep:ESmear");
  EEMCalclusterntuple = new TNtuple ("EEMCalclusterntup", "G4Clusters", "phi:z:edep:nTowers");

  std::cout << "ECCE_DEMP::Init(PHCompositeNode *topNode) Initializing" << std::endl;

  event_itt = 0;

  h2_ZDC_XY = new TH2F("ZDC_XY", "ZDC XY", 200, -50, 50, 200, -50, 50);
  h2_ZDC_XY_nEnergy = new TH2F("ZDC_XY_nEnergy", "ZDC XY - EDep > 70 GeV", 200, -50, 50, 200, -50, 50);
  h2_ZDC_XY_nEnergy_Smeared = new TH2F("ZDC_XY_nEnergy_Smeared", "ZDC XY - EDep > 70 GeV", 200, -50, 50, 200, -50, 50);

  h1_ZDC_E_dep = new TH1F("ZDC_E_dep", "ZDC E Deposit", 200, 0.0, 100.0);
  h1_ZDC_E_dep_smeared = new TH1F("ZDC_E_dep_smeared", "ZDC E Deposit", 240, 0.0, 120.0);

  h1_eTrack_px = new TH1F("eTrack_px", "e' Track P_{x}", 200, -10, 10);
  h1_eTrack_py = new TH1F("eTrack_py", "e' Track P_{y}", 200, -10, 10);
  h1_eTrack_pz = new TH1F("eTrack_pz", "e' Track P_{z}", 100, -10, 0);
  h1_eTrack_p = new TH1F("eTrack_p", "e' Track P", 100, 0, 10);
  h1_eTrack_theta = new TH1F("eTrack_theta", "e' Track #theta", 140, 110, 180);
  h1_eTrack_phi = new TH1F("eTrack_phi", "e' Track #phi", 720, -180, 180);

  h1_piTrack_px = new TH1F("piTrack_px", "#pi Track P_{x}", 200, -10, 10);
  h1_piTrack_py = new TH1F("piTrack_py", "#pi Track P_{y}", 200, -10, 10);
  h1_piTrack_pz = new TH1F("piTrack_pz", "#pi Track P_{z}", 500, 0, 50);
  h1_piTrack_p = new TH1F("piTrack_p", "#pi Track P", 500, 0, 50);
  h1_piTrack_theta = new TH1F("piTrack_theta", "#pi Track #theta", 120, 0, 60);
  h1_piTrack_phi = new TH1F("piTrack_phi", "#pi Track #phi", 720, -180, 180);

  h1_nTracksDist = new TH1F("nTracksDist", "Number of Tracks Per Event", 5, 0, 5);

  h2_ePiTrackDist = new TH2F("ePiTrackDist", " Number of pion tracks vs number of electron tracks; nTracks_{#pi}; nTracks_{e'}", 2, 0, 2, 2, 0, 2);

  h2_eTrack_ThetaPhi = new TH2F("eTrack_ThetaPhi", "e' Track #theta vs #phi; #theta [deg]; #phi [deg]", 140, 110, 180, 720, -180, 180);
  h2_eTrack_pTheta = new TH2F("eTrack_pTheta", "e' Track #theta vs P; #theta [deg]; P [GeV/c]", 140, 110, 180, 100, 0, 10);

  h2_piTrack_ThetaPhi = new TH2F("piTrack_ThetaPhi", "#pi Track #theta vs #phi; #theta [deg]; #phi [deg]", 120, 0, 60, 720, -180, 180);
  h2_piTrack_pTheta = new TH2F("piTrack_pTheta", "#pi Track #theta vs P; #theta [deg]; P [GeV/c]", 120, 0, 60, 500, 0, 50);

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

  process_g4hits_ZDC(topNode);

  //process_g4hits(topNode, "EEMC");
  //process_g4clusters(topNode, "EEMC");
  
  process_g4tracks(topNode);

  if (Check_n(topNode) == true && Check_ePi(topNode) == true){ // For event, check if it look like we have an e/pi/n in the event
    // Get track map for e'/pi info
    SvtxTrackMap* trackmap = findNode::getClass<SvtxTrackMap>(topNode, "SvtxTrackMap");
    // Get ZDC hits for neutron info
    PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDC");
    // Quantities we want to determine
    TVector3 eVect;
    TVector3 piVect;
    Double_t nEDep;

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
	}

	else if (track->get_pz() < 0  && track->get_charge() == -1 ){ // -ve z direction -> electrons, crappy way of selecting them for now w/o truth info
	  eVect.SetXYZ(track->get_px(), track->get_py(), track->get_pz());
	}
      }
    
    // Loop over the hts in the zdc
    if (hits) {
      // this returns an iterator to the beginning and the end of our G4Hits
      PHG4HitContainer::ConstRange hit_range = hits->getHits();
      for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
	{
	  nEDep = hit_iter->second->get_edep();
	}
    }
    cout << piVect.Mag() << "  " << eVect.Mag() << "  " << nEDep << endl;
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

  outfile->cd();
  g4hitntuple->Write();
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

int ECCE_DEMP::process_g4hits_ZDC(PHCompositeNode* topNode)
{
  ostringstream nodename;

  // loop over the G4Hits
  nodename.str("");
  nodename << "G4HIT_" << "ZDC";

  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());

  float smeared_E;

  if (hits) {
//    // this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++)
    {
	ZDC_hit++;
    }

    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {

      // the pointer to the G4Hit is hit_iter->second
      ZDChitntuple->Fill(hit_iter->second->get_x(0),
                        hit_iter->second->get_y(0),
                        hit_iter->second->get_z(0),
                        hit_iter->second->get_edep());

      h2_ZDC_XY->Fill(hit_iter->second->get_x(0)-90, hit_iter->second->get_y(0)); 

      smeared_E = EMCAL_Smear(hit_iter->second->get_edep());
      h1_ZDC_E_dep->Fill(hit_iter->second->get_edep());
      h1_ZDC_E_dep_smeared->Fill(smeared_E);

      // Event has roughly the correct energy for a "real" neutron event
      if( (hit_iter->second->get_edep()) > 70){
	h2_ZDC_XY_nEnergy->Fill(hit_iter->second->get_x(0)-90, hit_iter->second->get_y(0));
      }
    if( smeared_E > 70){
	h2_ZDC_XY_nEnergy_Smeared->Fill(hit_iter->second->get_x(0)-90, hit_iter->second->get_y(0));
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
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
int ECCE_DEMP::process_g4hits(PHCompositeNode* topNode, const string& detector)
{
  ostringstream nodename;

  // loop over the G4Hits
  nodename.str("");
  nodename << "G4HIT_" << detector;
  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, nodename.str().c_str());
  float smeared_E_EEMC;

  if (hits) {
    // this returns an iterator to the beginning and the end of our G4Hits
    PHG4HitContainer::ConstRange hit_range = hits->getHits();
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++){
      if(detector == "EEMC"){
	  EEMC_hit++;
	}
	// else if(detector == "CEMC"){
	//   CEMC_hit++;
	// }
	// else if(detector == "FEMC"){
	//   FEMC_hit++;
	// }
	// else if(detector == "HCALIN"){
	//   HCALIN_hit++;
	// }
	// else if(detector == "HCALOUT"){
	//   HCALOUT_hit++;
	// }
	// else if(detector == "FHCAL"){
	//   FHCAL_hit++;
	// }
      }
    for (PHG4HitContainer::ConstIterator hit_iter = hit_range.first; hit_iter != hit_range.second; hit_iter++) {
      if (detector == "EEMC"){
	smeared_E_EEMC = EMCAL_Smear(hit_iter->second->get_edep());
	// the pointer to the G4Hit is hit_iter->second
	EEMChitntuple->Fill(hit_iter->second->get_x(0),
			    hit_iter->second->get_y(0),
			    hit_iter->second->get_z(0),
			    hit_iter->second->get_edep(),
			    smeared_E_EEMC);
      }
      // else if (detector == "CEMC"){
      // 	smeared_E_CEMC = EMCAL_Smear(hit_iter->second->get_edep());
      // 	// the pointer to the G4Hit is hit_iter->second
      // 	CEMChitntuple->Fill(hit_iter->second->get_x(0),
      // 			    hit_iter->second->get_y(0),
      // 			    hit_iter->second->get_z(0),
      // 			    hit_iter->second->get_edep(),
      // 			    smeared_E_CEMC);    
      // }
      // else if (detector == "FEMC"){
      // 	smeared_E_FEMC = EMCAL_Smear(hit_iter->second->get_edep());
      // 	// the pointer to the G4Hit is hit_iter->second
      // 	FEMChitntuple->Fill(hit_iter->second->get_x(0),
      // 			    hit_iter->second->get_y(0),
      // 			    hit_iter->second->get_z(0),
      // 			    hit_iter->second->get_edep(),
      // 			    smeared_E_FEMC);
      // }
      // else if (detector == "HCALIN"){
      // 	smeared_E_HCALIN = HCAL_Smear(hit_iter->second->get_edep());
      // 	// the pointer to the G4Hit is hit_iter->second
      // 	HCALINhitntuple->Fill(hit_iter->second->get_x(0),
      // 			      hit_iter->second->get_y(0),
      // 			      hit_iter->second->get_z(0),
      // 			      hit_iter->second->get_edep(),
      // 			      smeared_E_HCALIN);
      // }
      // else if (detector == "HCALOUT"){
      // 	smeared_E_HCALOUT = HCAL_Smear(hit_iter->second->get_edep());
      // 	// the pointer to the G4Hit is hit_iter->second
      // 	HCALOUThitntuple->Fill(hit_iter->second->get_x(0),
      // 			      hit_iter->second->get_y(0),
      // 			      hit_iter->second->get_z(0),
      // 			      hit_iter->second->get_edep(),
      // 			      smeared_E_HCALOUT);      
      // }
      // else if (detector == "FHCAL"){
      // 	smeared_E_FHCAL = HCAL_Smear(hit_iter->second->get_edep());
      // 	// the pointer to the G4Hit is hit_iter->second
      // 	FHCALhitntuple->Fill(hit_iter->second->get_x(0),
      // 			      hit_iter->second->get_y(0),
      // 			      hit_iter->second->get_z(0),
      // 			      hit_iter->second->get_edep(),
      // 			      smeared_E_FHCAL);
      // }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int ECCE_DEMP::process_g4clusters(PHCompositeNode* topNode, const string& detector)
{
  ostringstream nodename;

  // loop over the G4Hits
  nodename.str("");
  nodename << "CLUSTER_" << detector;
  RawClusterContainer* clusters = findNode::getClass<RawClusterContainer>(topNode, nodename.str().c_str());
  if (clusters)
  {
    RawClusterContainer::ConstRange cluster_range = clusters->getClusters();
    for (RawClusterContainer::ConstIterator cluster_iter = cluster_range.first; cluster_iter != cluster_range.second; cluster_iter++)
      {// Fill ntuple depending on where clusters are coming from
      if (detector == "EEMC"){
	EEMCalclusterntuple->Fill(cluster_iter->second->get_phi(),
			    cluster_iter->second->get_z(),
			    cluster_iter->second->get_energy(),
			    cluster_iter->second->getNTowers());
      }
      // else if (detector == "CEMC"){
      // 	CEMCalclusterntuple->Fill(cluster_iter->second->get_phi(),
      // 			    cluster_iter->second->get_z(),
      // 			    cluster_iter->second->get_energy(),
      // 			    cluster_iter->second->getNTowers());
      // }
      // else if (detector == "FEMC"){
      // 	FEMCalclusterntuple->Fill(cluster_iter->second->get_phi(),
      // 			    cluster_iter->second->get_z(),
      // 			    cluster_iter->second->get_energy(),
      // 			    cluster_iter->second->getNTowers());
      // }
      // else if (detector == "HCALIN"){
      //   HCalInclusterntuple->Fill(cluster_iter->second->get_phi(),
      // 				 cluster_iter->second->get_z(),
      // 				 cluster_iter->second->get_energy(),
      // 				 cluster_iter->second->getNTowers());
      // }
      // else if (detector == "HCALOUT"){
      //   HCalOutclusterntuple->Fill(cluster_iter->second->get_phi(),
      // 				 cluster_iter->second->get_z(),
      // 				 cluster_iter->second->get_energy(),
      // 				 cluster_iter->second->getNTowers());
      // }
      // else if (detector == "FHCAL"){
      //   FHCalclusterntuple->Fill(cluster_iter->second->get_phi(),
      // 				 cluster_iter->second->get_z(),
      // 				 cluster_iter->second->get_energy(),
      // 				 cluster_iter->second->getNTowers());
      // }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//***************************************************

int ECCE_DEMP::process_g4tracks(PHCompositeNode* topNode)
{
  // Statement to check node actually exists, need to tweak this and include it
  // PHNode *findNode = dynamic_cast<PHNode*>(nodeIter.findFirst("TrackMap"));
  // if (findNode)
  //   {
  //     dst_trackmap = findNode::getClass<SvtxTrackMap>(topNode, "TrackMap");
  //     if (Verbosity() >= VERBOSITY_A_LOT) std::cout << __FILE__ << "::" << __func__ << "::" << __LINE__ << ": Number of tracks: " << dst_trackmap->size() << std::endl; 
  //   }
  // else
  //   {
  //     std::cout << __FILE__ << "::" << __func__ << "::" << __LINE__<< ": TrackMap does not exist" << std::endl;
  //   }

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
	TVector3 piVect(track->get_px(), track->get_py(), track->get_pz());
	g4trackntuple->Fill(track->get_px(),
			    track->get_py(),
			    track->get_pz(),
			    piVect.Mag(),
			    piVect.Theta(),
			    piVect.Phi());

	h1_piTrack_px->Fill(track->get_px());
	h1_piTrack_py->Fill(track->get_py());
	h1_piTrack_pz->Fill(track->get_pz());
	h1_piTrack_p->Fill(piVect.Mag());
	h1_piTrack_theta->Fill(piVect.Theta()*TMath::RadToDeg());
	h1_piTrack_phi->Fill(piVect.Phi()*TMath::RadToDeg());

	h2_piTrack_ThetaPhi->Fill(piVect.Theta()*TMath::RadToDeg(), (piVect.Phi()*TMath::RadToDeg()));
	h2_piTrack_pTheta->Fill((piVect.Theta()*TMath::RadToDeg()), piVect.Mag());
      }

      else if (track->get_pz() < 0  && track->get_charge() == -1 ){ // -ve z direction -> electrons, crappy way of selecting them for now w/o truth info
	ElecTrack = kTRUE;
	TVector3 eVect(track->get_px(), track->get_py(), track->get_pz());
	g4trackntuple->Fill(track->get_px(),
			    track->get_py(),
			    track->get_pz(),
			    eVect.Mag(),
			    eVect.Theta(),
			    eVect.Phi());

	h1_eTrack_px->Fill(track->get_px());
	h1_eTrack_py->Fill(track->get_py());
	h1_eTrack_pz->Fill(track->get_pz());
	h1_eTrack_p->Fill(eVect.Mag());
	h1_eTrack_theta->Fill(eVect.Theta()*TMath::RadToDeg());
	h1_eTrack_phi->Fill(eVect.Phi()*TMath::RadToDeg());
	
	h2_eTrack_ThetaPhi->Fill((eVect.Theta()*TMath::RadToDeg()), (eVect.Phi()*TMath::RadToDeg()));
	h2_eTrack_pTheta->Fill((eVect.Theta()*TMath::RadToDeg()), eVect.Mag());
      }
    }

  h1_nTracksDist->Fill(nTracks);

  if( PionTrack == kFALSE && ElecTrack == kFALSE ){ // No pion or electron track
    h2_ePiTrackDist->Fill(0.,0.);
  }
  else if( PionTrack == kTRUE && ElecTrack == kFALSE ){ // Pion track but no electron track
    h2_ePiTrackDist->Fill(1.,0.);
  }
  else if( PionTrack == kFALSE && ElecTrack == kTRUE ){ // No pion track but an electron track
    h2_ePiTrackDist->Fill(0.,1.);
  }
  else if( PionTrack == kTRUE && ElecTrack == kTRUE ){ // Both a pion and an electron track
    h2_ePiTrackDist->Fill(1.,1.);
  }

  return Fun4AllReturnCodes::EVENT_OK;
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

  PHG4HitContainer* hits = findNode::getClass<PHG4HitContainer>(topNode, "G4HIT_ZDC");

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
