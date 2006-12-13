#ifndef RecoMuon_GlobalMuonProducer_GlobalMuonValidator_H
#define RecoMuon_GlobalMuonProducer_GlobalMuonValidator_H

/** \class GlobalMuonValidator
 *  Analyzer of the StandAlone and Global muon tracks
 *
 *  $Date:  $
 *  $Revision:  $
 *  \author A. Everett     Purdue University
 */

// Base Class Headers
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include <TROOT.h>
#include <TSystem.h>

namespace edm {
  class ParameterSet;
  //  class Event;
  class EventSetup;
  class InputTag;
}

class TFile;
class TH1F;
class TH2F;
//class TrackRef;
//class SimTrackRef;
//class MuonRef;
class MuonServiceProxy;

using namespace std;
using namespace edm;
using namespace reco;

class GlobalMuonValidator : public edm::EDAnalyzer {
public:
  explicit GlobalMuonValidator(const edm::ParameterSet&);
  ~GlobalMuonValidator();
  
  typedef std::pair< TrackRef, SimTrackRef> CandStaSim;
  typedef std::pair< MuonRef,  SimTrackRef> CandMuonSim;
  
private:
  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual float calculateDistance(const math::XYZVector&, const math::XYZVector&);
  virtual TH1F* divideErr(TH1F*, TH1F*, TH1F*);

  // ----------member data ---------------------------

  edm::InputTag TKtrackTags_; 
  edm::InputTag STAtrackTags_; 
  edm::InputTag MuonTags_; 
  edm::InputTag SIMtrackTags_; 

  string out, open;
  double  theMinEta, theMaxEta, theMinPt, thePtCut;
  int theNBins, thePartID;

  Handle<reco::MuonCollection> MuCollection;
  Handle<reco::TrackCollection> TKTrackCollection;
  Handle<reco::TrackCollection> STATrackCollection;
  Handle<edm::SimTrackContainer> SIMTrackCollection;
  
  MuonServiceProxy* theService;

  //ROOT Pointers
  TFile* hFile;
  TStyle* effStyle;

  TH1F* hi_dist_glb_sta;
  TH1F* hi_dist_glb_sim;
  TH1F* hi_dist_sta_sim;

  TH1F* hi_sim_pt  ;
  TH1F* hi_sta_pt  ;
  TH1F* hi_sta2_pt  ;
  TH1F* hi_glb_pt  ;

  TH1F* hi_sim_eta  ;
  TH1F* hi_sta_eta  ;
  TH1F* hi_sta2_eta  ;
  TH1F* hi_glb_eta  ;

  TH1F* hi_glbsim_pt  ;
  TH1F* hi_stasim_pt  ;
  TH1F* hi_glbsta_pt  ;

  TH1F* hi_glbsim_eta  ;
  TH1F* hi_stasim_eta  ;
  TH1F* hi_glbsta_eta  ;

  TH1F* hi_glbsim2_pt  ;
  TH1F* hi_stasim2_pt  ;
  TH1F* hi_glbsta2_pt  ;

  TH1F* hi_glbsim2_eta  ;
  TH1F* hi_stasim2_eta  ;
  TH1F* hi_glbsta2_eta  ;

  TH1F* hi_glbsim_eff_pt  ;
  TH1F* hi_stasim_eff_pt  ;
  TH1F* hi_glbsta_eff_pt  ;

  TH1F* hi_glbsim_eff_eta  ;
  TH1F* hi_stasim_eff_eta  ;
  TH1F* hi_glbsta_eff_eta  ;

  TH1F* hi_glbsim_pur_pt  ;
  TH1F* hi_stasim_pur_pt  ;
  TH1F* hi_glbsta_pur_pt  ;

  TH1F* hi_glbsim_pur_eta  ;
  TH1F* hi_stasim_pur_eta  ;
  TH1F* hi_glbsta_pur_eta  ;

  TH1F* hi_glbsim_ptres;
  TH1F* hi_stasim_ptres;
  TH1F* hi_glbsta_ptres;

  TH1F* hi_glbsim_etares;
  TH1F* hi_stasim_etares;
  TH1F* hi_glbsta_etares;

  TH1F* hi_fails;
};
#endif
