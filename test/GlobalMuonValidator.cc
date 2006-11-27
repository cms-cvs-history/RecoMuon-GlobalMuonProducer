// -*- C++ -*-
//
// Package:    GlobalMuonValidator
// Class:      GlobalMuonValidator
// 
/**\class GlobalMuonValidator GlobalMuonValidator.cc AEverett/GlobalMuonValidator/src/GlobalMuonValidator.cc

 Description: Validator tool to study efficiencies and purities.

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Adam A Everett
//         Created:  Wed Sep 27 14:54:28 EDT 2006
// $Id: GlobalMuonValidator.cc,v 1.1 2006/10/25 14:52:22 aeverett Exp $
//
//


// system include files
#include <memory>
#include <string>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "DataFormats/Common/interface/Ref.h"

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFrame.h>

//
// class decleration
//

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
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GlobalMuonValidator::GlobalMuonValidator(const edm::ParameterSet& iConfig) 
  :
  TKtrackTags_(iConfig.getUntrackedParameter<edm::InputTag>("TKtracks")),
  STAtrackTags_(iConfig.getUntrackedParameter<edm::InputTag>("STAtracks")),
  MuonTags_(iConfig.getUntrackedParameter<edm::InputTag>("Muons")),
  SIMtrackTags_(iConfig.getUntrackedParameter<edm::InputTag>("SIMtracks")),
  out(iConfig.getParameter<string>("out")),
  open(iConfig.getParameter<string>("open")),
  theMinEta(iConfig.getParameter<double>("etaMin")),
  theMaxEta(iConfig.getParameter<double>("etaMax")),
  theMinPt(iConfig.getParameter<double>("simPtMin")),
  thePtCut(iConfig.getParameter<double>("PtCut")),
  theNBins(iConfig.getParameter<int>("nbins")),
  thePartID(iConfig.getParameter<int>("partId"))
{
  //now do what ever initialization is needed

  // service parameters
  ParameterSet serviceParameters = iConfig.getParameter<ParameterSet>("ServiceParameters");
  // the services
  theService = new MuonServiceProxy(serviceParameters);
  
}


GlobalMuonValidator::~GlobalMuonValidator()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
  if (hFile!=0) {
    hFile->Close();
    delete hFile;
  }
  if (theService) delete theService;
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
GlobalMuonValidator::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  //using namespace edm;
  using reco::TrackCollection;
  using reco::MuonCollection;

  // Update the services
  //theService->update(iSetup);

  iEvent.getByLabel(SIMtrackTags_,SIMTrackCollection);
  const SimTrackContainer simTC = *(SIMTrackCollection.product());
  
  iEvent.getByLabel( TKtrackTags_, TKTrackCollection);
  const reco::TrackCollection tkTC = *(TKTrackCollection.product());

  iEvent.getByLabel( STAtrackTags_, STATrackCollection);
  const reco::TrackCollection staTC = *(STATrackCollection.product());
  
  iEvent.getByLabel(MuonTags_,MuCollection);
  const reco::MuonCollection muonC = *(MuCollection.product());

  if(muonC.size() == 0) {
    cout << "*****No GLBMuon in run " << iEvent.id() << " event " << iEvent.id().event() << endl;
    cout << "     GLBMuons " << muonC.size() << " STAMuons " << staTC.size() << " TkTracks  " << tkTC.size() <<endl;
    if(staTC.size() > 0 && tkTC.size() > 0){
      cout << "          IGUANA this event! Run: " << iEvent.id() << " Event: " << iEvent.id().event() << endl;
    }
  }
  
  SimTrackRefVector simMuons1;
  int position = 0;
  SimTrackContainer::const_iterator simTrack;
  for (simTrack = simTC.begin(); simTrack != simTC.end(); ++simTrack){
    position++;
    if (abs((*simTrack).type()) == 13
	&& (*simTrack).momentum().perp() >= theMinPt 
	&& (*simTrack).momentum().eta() <= theMaxEta
	&& (*simTrack).momentum().eta() >= theMinEta
	) {
      SimTrackRef simTrackRef(SIMTrackCollection,position-1);
      simMuons1.push_back(simTrackRef);
      hi_sim_pt->Fill((*simTrack).momentum().perp());
      hi_sim_eta->Fill(((*simTrack).momentum().eta()));
    }    
  }

  reco::TrackRefVector staMuons1,staMuons2;
  TrackCollection::const_iterator staTrack;
  position = 0;
  for (staTrack = staTC.begin(); staTrack != staTC.end(); ++staTrack){
    position++;
    staMuons1.push_back(TrackRef(STATrackCollection,position-1));
    //staMuons2.push_back(TrackRef(STATrackCollection,position-1));
    hi_sta_pt->Fill((*staTrack).pt());
    hi_sta_eta->Fill(((*staTrack).eta()));
    bool staCut = true;
    if ( (*staTrack).pt() < thePtCut || (*staTrack).innerMomentum().Rho() < thePtCut || (*staTrack).innerMomentum().R() < 2.5 ) staCut = false;
    if(staCut && tkTC.size() > 0) {
      staMuons2.push_back(TrackRef(STATrackCollection,position-1));
      hi_sta2_pt->Fill((*staTrack).pt());
      hi_sta2_eta->Fill(((*staTrack).eta()));         
    }
  }

  reco::MuonRefVector glbMuons1, glbMuons2;
  MuonCollection::const_iterator glbMuon;  
  position = 0;
  for (glbMuon = muonC.begin(); glbMuon != muonC.end(); ++glbMuon){
    position++;
    TrackRef glbTrack = glbMuon->combinedMuon();
    glbMuons1.push_back(MuonRef(MuCollection,position-1));
    glbMuons2.push_back(MuonRef(MuCollection,position-1));
    hi_glb_pt->Fill((*glbTrack).pt());
    hi_glb_eta->Fill(((*glbTrack).eta())); 
  }
  
  for (MuonRefVector::const_iterator glbTrack = glbMuons1.begin(); glbTrack != glbMuons1.end(); ++glbTrack){
    for (TrackRefVector::const_iterator staTrack = staMuons1.begin(); staTrack != staMuons1.end(); ++staTrack){
      float D = calculateDistance((*glbTrack)->combinedMuon()->momentum(),(*staTrack)->momentum());
      hi_dist_glb_sta->Fill(D);
    }
  }
  
  for (MuonRefVector::const_iterator glbTrack = glbMuons1.begin(); glbTrack != glbMuons1.end(); ++glbTrack){
    for (SimTrackRefVector::const_iterator simTrack = simMuons1.begin(); simTrack != simMuons1.end(); ++simTrack){
      const math::XYZVector simVect((*simTrack)->momentum().x(),(*simTrack)->momentum().y(),(*simTrack)->momentum().z());
      float D = calculateDistance((*glbTrack)->combinedMuon()->momentum(),simVect);
      hi_dist_glb_sim->Fill(D);
    }
  }
  
  for (TrackRefVector::const_iterator staTrack = staMuons1.begin(); staTrack != staMuons1.end(); ++staTrack){
    for (SimTrackRefVector::const_iterator simTrack = simMuons1.begin(); simTrack != simMuons1.end(); ++simTrack){
      const math::XYZVector simVect((*simTrack)->momentum().x(),(*simTrack)->momentum().y(),(*simTrack)->momentum().z());
      float D = calculateDistance((*staTrack)->momentum(),simVect);
      hi_dist_sta_sim->Fill(D);
    }
  }
  
  //
  //Now make pairs
  //
  //FIXME: add a Distance cut
  
  //Global-Simulated pairs
  vector<CandMuonSim> pairGlbSim;
  for (SimTrackRefVector::const_iterator simTrack = simMuons1.begin(); simTrack != simMuons1.end(); ++simTrack){
    float minDist = 9999.;
    int index = 0;
    int keep = -1;
    MuonRefVector::const_iterator match;
    for (MuonRefVector::const_iterator glbTrack = glbMuons1.begin(); glbTrack != glbMuons1.end(); ++glbTrack){
      //float deta = (*glbTrack)->combinedMuon()->eta() - (*simTrack)->momentum().eta();
      //float dphi = (*glbTrack)->combinedMuon()->phi() - (*simTrack)->momentum().phi();
      //float dR = sqrt(pow(deta,2)+pow(dphi,2));
      const math::XYZVector simVect((*simTrack)->momentum().x(),(*simTrack)->momentum().y(),(*simTrack)->momentum().z());
      float dR = calculateDistance((*glbTrack)->momentum(),simVect); 
      
      if(dR <= minDist) {
	minDist = dR;
	keep = index;
	match = glbTrack;
      }
      index++;
    }
    
    if( keep > -1 ) {
      CandMuonSim tmp = CandMuonSim((*match),(*simTrack));
      pairGlbSim.push_back(tmp);
      glbMuons1.erase(match);
    }
  }
  
  //Standalone-Simulated pairs
  vector<CandStaSim> pairStaSim;
  for (SimTrackRefVector::const_iterator simTrack = simMuons1.begin(); simTrack != simMuons1.end(); ++simTrack){
    float minDist = 9999.;
    int index = 0;
    int keep = -1;
    TrackRefVector::const_iterator match;
    for (TrackRefVector::const_iterator staTrack = staMuons1.begin(); staTrack != staMuons1.end(); ++staTrack){
      //float deta = (*staTrack)->eta() - (*simTrack)->momentum().eta();
      //float dphi = (*staTrack)->phi() - (*simTrack)->momentum().phi();
      //float dR = sqrt(pow(deta,2)+pow(dphi,2));
      const math::XYZVector simVect((*simTrack)->momentum().x(),(*simTrack)->momentum().y(),(*simTrack)->momentum().z());
      float dR = calculateDistance((*staTrack)->momentum(),simVect); 
      
      if(dR <= minDist) {
	minDist = dR;
	keep = index;
	match = staTrack;
      }
      index++;
    }
    
    if( keep > -1 ) {
      CandStaSim tmp = CandStaSim((*match),(*simTrack));
      pairStaSim.push_back(tmp);
      staMuons1.erase(match);
    }
  }
  
  //
  //Fill histograms
  //  
  for (vector<CandMuonSim>::const_iterator pair = pairGlbSim.begin(); pair != pairGlbSim.end(); ++pair){
    hi_glbsim_pt->Fill((*pair).second->momentum().perp());
    hi_glbsim_eta->Fill(((*pair).second->momentum().eta()));       
    hi_glbsim2_pt->Fill((*pair).first->pt());
    hi_glbsim2_eta->Fill(((*pair).first->eta()));       
    
    int qSIM = -1;
    if( (*pair).second->type() == -13 ) qSIM = 1;
    int qGLB = (*pair).first->charge();
    
    float ptSim = (*pair).second->momentum().perp();
    float ptGLB = (*pair).first->pt();
    hi_glbsim_ptres->Fill((qGLB/ptGLB - qSIM/ptSim)/(qSIM/ptSim));
    
    float etaSim = (*pair).second->momentum().eta();
    float etaGLB = (*pair).first->eta();
    hi_glbsim_etares->Fill((qGLB/etaGLB - qSIM/etaSim)/(qSIM/etaSim));
  }

  for (vector<CandStaSim>::const_iterator pair = pairStaSim.begin(); pair != pairStaSim.end(); pair++){
  hi_stasim_pt->Fill((*pair).second->momentum().perp());
  hi_stasim_eta->Fill(((*pair).second->momentum().eta()));
  hi_stasim2_pt->Fill((*pair).first->pt());
  hi_stasim2_eta->Fill(((*pair).first->eta()));
  
  int qSIM = -1;
  if( (*pair).second->type() == -13 ) qSIM = 1;
  int qSTA = (*pair).first->charge();

  float ptSim = (*pair).second->momentum().perp();
  float ptSTA = (*pair).first->pt();
  hi_stasim_ptres->Fill((qSTA/ptSTA-qSIM/ptSim)/(qSIM/ptSim));
  
  float etaSim = (*pair).second->momentum().eta();
  float etaSTA = (*pair).first->eta();
  hi_stasim_etares->Fill((qSTA/etaSTA-qSIM/etaSim)/(qSIM/etaSim));
  }

  for (TrackRefVector::const_iterator staTrack = staMuons2.begin(); staTrack != staMuons2.end(); ++staTrack){
    double chi2 = 9999.;
    int index = 0;
    int keep = -1;
    MuonRefVector::const_iterator match;
    for (MuonRefVector::const_iterator glbTrack = glbMuons2.begin(); glbTrack != glbMuons2.end(); ++glbTrack){
      if((*staTrack) == (*glbTrack)->standAloneMuon()) {
	double tmp_chi2 = (*glbTrack)->combinedMuon()->chi2();
	if(tmp_chi2 < chi2) {
	  chi2 = tmp_chi2;
	  keep = index;
	  match = glbTrack;
	}
	index++;
      }
    }
    if(glbMuons2.size() > 0 && keep > -1) {
      hi_glbsta2_pt->Fill( (*match)->combinedMuon()->pt() );
      hi_glbsta2_eta->Fill( (*match)->combinedMuon()->eta() );
      hi_glbsta_pt->Fill( (*match)->standAloneMuon()->pt() );
      hi_glbsta_eta->Fill( (*match)->standAloneMuon()->eta() );
    }
  }

}


// ------------ method called once each job just before starting event loop  ------------
void 
GlobalMuonValidator::beginJob(const edm::EventSetup&)
{
   hFile = new TFile( out.c_str(), open.c_str() );

   effStyle = new TStyle("effStyle","Efficiency Study Style");   
   effStyle->SetCanvasBorderMode(0);
   effStyle->SetPadBorderMode(1);
   effStyle->SetOptTitle(0);
   effStyle->SetStatFont(42);
   effStyle->SetTitleFont(22);
   effStyle->SetCanvasColor(10);
   effStyle->SetPadColor(0);
   effStyle->SetLabelFont(42,"x");
   effStyle->SetLabelFont(42,"y");
   effStyle->SetHistFillStyle(1001);
   effStyle->SetHistFillColor(0);
   effStyle->SetOptStat(0);
   effStyle->SetOptFit(0111);
   effStyle->SetStatH(0.05);


   hi_dist_glb_sta = new TH1F("hi_dist_glb_sta","Distance GLB-STA",100,0.,100.);
   hi_dist_glb_sim = new TH1F("hi_dist_glb_sim","Distance GLB-SIM",100,0.,100.);
   hi_dist_sta_sim = new TH1F("hi_dist_sta_sim","Distance STA-SIM",100,0.,100.);

   hi_sim_pt  = new TH1F("hi_sim_pt","P_{T}^{sim}",theNBins,0.0,100.);
   hi_sim_pt->Sumw2();
   hi_sta_pt   = new TH1F("hi_sta_pt","P_{T}^{STA}",theNBins,0.0,100.);
   hi_sta_pt->Sumw2();
   hi_sta2_pt   = new TH1F("hi_sta2_pt","P_{T}^{STA}",theNBins,0.0,100.);
   hi_sta2_pt->Sumw2();
   hi_glb_pt   = new TH1F("hi_glb_pt","P_{T}^{GLB}",theNBins,0.0,100.);
   hi_glb_pt->Sumw2();

   hi_sim_eta  = new TH1F("hi_sim_eta","#eta^{sim}",theNBins,theMinEta,theMaxEta);
   hi_sim_eta->Sumw2();
   hi_sta_eta   = new TH1F("hi_sta_eta","#eta^{STA}",theNBins,theMinEta,theMaxEta);
   hi_sta_eta->Sumw2();
   hi_sta2_eta   = new TH1F("hi_sta2_eta","#eta^{STA}",theNBins,theMinEta,theMaxEta);
   hi_sta2_eta->Sumw2();
   hi_glb_eta   = new TH1F("hi_glb_eta","#eta^{GLB}",theNBins,theMinEta,theMaxEta);
   hi_glb_eta->Sumw2();

   hi_glbsim_pt  = new TH1F("hi_glbsim_pt","P_{T}^{GLB,sim}",theNBins,0.0,100.);
   hi_glbsim_pt->Sumw2();
   hi_stasim_pt   = new TH1F("hi_stasim_pt","P_{T}^{STA,sim}",theNBins,0.0,100.);
   hi_stasim_pt->Sumw2();
   hi_glbsta_pt   = new TH1F("hi_glbsta_pt","P_{T}^{GLB,STA}",theNBins,0.0,100.);
   hi_glbsta_pt->Sumw2();

   hi_glbsim_eta  = new TH1F("hi_glbsim_eta","#eta^{GLB,sim}",theNBins,theMinEta,theMaxEta);
   hi_glbsim_eta->Sumw2();
   hi_stasim_eta   = new TH1F("hi_stasim_eta","#eta^{STA,sim}",theNBins,theMinEta,theMaxEta);
   hi_stasim_eta->Sumw2();
   hi_glbsta_eta   = new TH1F("hi_glbsta_eta","#eta^{GLB,STA}",theNBins,theMinEta,theMaxEta);
   hi_glbsta_eta->Sumw2();

   hi_glbsim2_pt  = new TH1F("hi_glbsim2_pt","P_{T}^{GLB,sim}",theNBins,0.0,100.);
   hi_glbsim2_pt->Sumw2();
   hi_stasim2_pt   = new TH1F("hi_stasim2_pt","P_{T}^{STA,sim}",theNBins,0.0,100.);
   hi_stasim2_pt->Sumw2();
   hi_glbsta2_pt   = new TH1F("hi_glbsta2_pt","P_{T}^{GLB,STA}",theNBins,0.0,100.);
   hi_glbsta2_pt->Sumw2();

   hi_glbsim2_eta  = new TH1F("hi_glbsim2_eta","#eta^{GLB,sim}",theNBins,theMinEta,theMaxEta);
   hi_glbsim2_eta->Sumw2();
   hi_stasim2_eta   = new TH1F("hi_stasim2_eta","#eta^{STA,sim}",theNBins,theMinEta,theMaxEta);
   hi_stasim2_eta->Sumw2();
   hi_glbsta2_eta   = new TH1F("hi_glbsta2_eta","#eta^{GLB,STA}",theNBins,theMinEta,theMaxEta);
   hi_glbsta2_eta->Sumw2();

   hi_glbsim_eff_pt  = new TH1F("hi_glbsim_eff_pt","Efficiency GLB,sim",theNBins,0.0,100.);
   hi_stasim_eff_pt   = new TH1F("hi_stasim_eff_pt","Efficiency STA,sim",theNBins,0.0,100.);
   hi_glbsta_eff_pt   = new TH1F("hi_glbsta_eff_pt","Efficiency GLB,STA",theNBins,0.0,100.);

   hi_glbsim_eff_eta  = new TH1F("hi_glbsim_eff_eta","Efficiency GLB,sim",theNBins,theMinEta,theMaxEta);
   hi_stasim_eff_eta   = new TH1F("hi_stasim_eff_eta","Efficiency STA,sim",theNBins,theMinEta,theMaxEta);
   hi_glbsta_eff_eta   = new TH1F("hi_glbsta_eff_eta","Efficiency GLB,STA",theNBins,theMinEta,theMaxEta);

   hi_glbsim_pur_pt  = new TH1F("hi_glbsim_pur_pt","Purity GLB,sim",theNBins,0.0,100.);
   hi_stasim_pur_pt   = new TH1F("hi_stasim_pur_pt","Purity STA,sim",theNBins,0.0,100.);
   hi_glbsta_pur_pt   = new TH1F("hi_glbsta_pur_pt","Purity GLB,STA",theNBins,0.0,100.);

   hi_glbsim_pur_eta  = new TH1F("hi_glbsim_pur_eta","Purity GLB,sim",theNBins,theMinEta,theMaxEta);
   hi_stasim_pur_eta   = new TH1F("hi_stasim_pur_eta","Purity STA,sim",theNBins,theMinEta,theMaxEta);
   hi_glbsta_pur_eta   = new TH1F("hi_glbsta_pur_eta","Purity GLB,STA",theNBins,theMinEta,theMaxEta);


  hi_glbsim_ptres  = new TH1F("hi_glbsim_ptres","P_{T} resolution GLB,sim",theNBins,-0.2,0.2);
  hi_stasim_ptres  = new TH1F("hi_stasim_ptres","P_{T} resolution STA.sim",theNBins,-0.2,0.2);
  hi_glbsta_ptres  = new TH1F("hi_glbsta_ptres","P_{T} resolution GLB,STA",theNBins,-0.2,0.2);

  hi_glbsim_etares  = new TH1F("hi_glbsim_etares","#eta resolution GLB,sim",theNBins,-1.,1.);
  hi_stasim_etares  = new TH1F("hi_stasim_etares","#eta resolution STA.sim",theNBins,-1.,1.);
  hi_glbsta_etares  = new TH1F("hi_glbsta_etares","#eta resolution GLB,STA",theNBins,-0.1,0.1);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GlobalMuonValidator::endJob() {

  hi_glbsim_eff_pt = divideErr(hi_glbsim_pt,hi_sim_pt,hi_glbsim_eff_pt);
  hi_stasim_eff_pt = divideErr(hi_stasim_pt,hi_sim_pt,hi_stasim_eff_pt);
  hi_glbsta_eff_pt = divideErr(hi_glbsta_pt,hi_sta2_pt,hi_glbsta_eff_pt);

  hi_glbsim_eff_eta = divideErr(hi_glbsim_eta,hi_sim_eta,hi_glbsim_eff_eta);
  hi_stasim_eff_eta = divideErr(hi_stasim_eta,hi_sim_eta,hi_stasim_eff_eta);
  hi_glbsta_eff_eta = divideErr(hi_glbsta_eta,hi_sta2_eta,hi_glbsta_eff_eta);

  //FIXME: potential problem with migration
  hi_glbsim_pur_pt = divideErr(hi_glbsim2_pt,hi_glb_pt,hi_glbsim_pur_pt);
  hi_stasim_pur_pt = divideErr(hi_stasim2_pt,hi_sta_pt,hi_stasim_pur_pt);
  hi_glbsta_pur_pt = divideErr(hi_glbsta2_pt,hi_glb_pt,hi_glbsta_pur_pt);

  hi_glbsim_pur_eta = divideErr(hi_glbsim2_eta,hi_glb_eta,hi_glbsim_pur_eta);
  hi_stasim_pur_eta = divideErr(hi_stasim2_eta,hi_sta_eta,hi_stasim_pur_eta);
  hi_glbsta_pur_eta = divideErr(hi_glbsta2_eta,hi_glb_eta,hi_glbsta_pur_eta);

  //...................  
  //effStyle->SetCanvasBorderMode(0);
  //effStyle->SetPadBorderMode(1);
  //effStyle->SetOptTitle(0);
  //effStyle->SetStatFont(42);
  //effStyle->SetTitleFont(22);
  //effStyle->SetCanvasColor(10);
  //effStyle->SetPadColor(0);
  //effStyle->SetLabelFont(42,"x");
  //effStyle->SetLabelFont(42,"y");
  //effStyle->SetHistFillStyle(1001);
  //effStyle->SetHistFillColor(0);
  //effStyle->SetOptStat(0);
  //effStyle->SetOptFit(0111);
  //effStyle->SetStatH(0.05);
  //.................... 

  gROOT->SetStyle("effStyle");

  char* l1string = "Simulated";
  char* l2string = "StandAlone";
  char* l3string = "Global";

  TCanvas* c1 = new TCanvas("eff","Efficiency pt",10,10,700,500);
  c1->SetFillColor(0);
  c1->SetGrid(1);
  c1->SetTicky();
  c1->SetRightMargin(0.03);
  c1->SetTopMargin(0.02);
  c1->cd(); 

  hi_stasim_eff_pt->SetXTitle("p_{T}^{#mu}");
  hi_stasim_eff_pt->SetYTitle("Efficiency");
  hi_stasim_eff_pt->SetTitleOffset(1.1,"x");
  hi_stasim_eff_pt->SetTitleOffset(1.15,"y");
  hi_stasim_eff_pt->SetMaximum(1.02);
  hi_stasim_eff_pt->SetMinimum(0.5);
  hi_stasim_eff_pt->SetStats(false);

  hi_stasim_eff_pt->SetLineWidth(1);
  hi_glbsim_eff_pt->SetLineWidth(1);
  hi_stasim_eff_pt->SetLineColor(2);
  hi_glbsim_eff_pt->SetLineColor(4);
  hi_stasim_eff_pt->SetLineStyle(1);
  hi_glbsim_eff_pt->SetLineStyle(1);

  //hi_stasim_eff_pt->SetFillColor(5);
  //hi_glbsim_eff_pt->SetFillColor(7);

  hi_stasim_eff_pt->SetMarkerStyle(21);
  hi_glbsim_eff_pt->SetMarkerStyle(26);
  hi_stasim_eff_pt->SetMarkerColor(2);
  hi_glbsim_eff_pt->SetMarkerColor(4);

  hi_stasim_eff_pt ->DrawCopy("PE");
  hi_glbsim_eff_pt ->DrawCopy("PEsame");
  //hi_stasim_eff_pt ->DrawCopy("AxisSame");

  TLegend* legend1 = new TLegend(0.6,0.2,0.8,0.4);
  legend1->SetTextAlign(32);
  legend1->SetTextColor(1);
  legend1->SetTextSize(0.04);
  legend1->AddEntry("hi_stasim_eff_pt",l2string,"pl");
  legend1->AddEntry("hi_glbsim_eff_pt",l3string,"pl");
  legend1 ->Draw();
  c1->Write();
  // 

  TCanvas* c1a = new TCanvas("algEff","Algo Efficiency pt",10,10,700,500);
  c1a->SetFillColor(0);
  c1a->SetGrid(1);
  c1a->SetTicky();
  c1a->SetRightMargin(0.03);
  c1a->SetTopMargin(0.02);
  c1a->cd(); 
  
  hi_glbsta_eff_pt->SetXTitle("p_{T}^{#mu}");
  hi_glbsta_eff_pt->SetYTitle("Efficiency");
  hi_glbsta_eff_pt->SetTitleOffset(1.1,"x");
  hi_glbsta_eff_pt->SetStats(false);

  hi_glbsta_eff_pt->SetMaximum(1.01);
  hi_glbsta_eff_pt->SetMinimum(0.8);
  hi_glbsta_eff_pt->SetLineColor(2);
  hi_glbsta_eff_pt->SetMarkerStyle(22);
  hi_glbsta_eff_pt ->DrawCopy("E");
  c1a->Write();


  TCanvas* c2 = new TCanvas("eff_eta","Efficiency eta",10,10,700,500);
  c2->SetFillColor(0);
  c2->SetGrid(1);
  c2->SetTicky();
  c2->SetRightMargin(0.03);
  c2->SetTopMargin(0.02);
  c2->cd(); 

  hi_stasim_eff_eta->SetXTitle("#eta^{#mu}");
  hi_stasim_eff_eta->SetYTitle("Efficiency");
  hi_stasim_eff_eta->SetTitleOffset(1.1,"x");
  hi_stasim_eff_eta->SetTitleOffset(1.15,"y");
  hi_stasim_eff_eta->SetMaximum(1.02);
  hi_stasim_eff_eta->SetMinimum(0.5);
  hi_stasim_eff_eta->SetStats(false);

  hi_stasim_eff_eta->SetLineWidth(1);
  hi_glbsim_eff_eta->SetLineWidth(1);
  hi_stasim_eff_eta->SetLineColor(2);
  hi_glbsim_eff_eta->SetLineColor(4);
  hi_stasim_eff_eta->SetLineStyle(1);
  hi_glbsim_eff_eta->SetLineStyle(1);

  //hi_stasim_eff_eta->SetFillColor(5);
  //hi_glbsim_eff_eta->SetFillColor(7);

  hi_stasim_eff_eta->SetMarkerStyle(21);
  hi_glbsim_eff_eta->SetMarkerStyle(26);
  hi_stasim_eff_eta->SetMarkerColor(2);
  hi_glbsim_eff_eta->SetMarkerColor(4);

  hi_stasim_eff_eta ->DrawCopy("PE");
  hi_glbsim_eff_eta ->DrawCopy("PEsame");
  //hi_stasim_eff_eta ->DrawCopy("AxisSame");

  TLegend* legend2 = new TLegend(0.6,0.2,0.8,0.4);
  legend2->SetTextAlign(32);
  legend2->SetTextColor(1);
  legend2->SetTextSize(0.04);
  legend2->AddEntry("hi_stasim_eff_eta",l2string,"pl");
  legend2->AddEntry("hi_glbsim_eff_eta",l3string,"pl");
  legend2 ->Draw();
  c2->Write();
  // 

  TCanvas* c2a = new TCanvas("algEff_eta","Algo Efficiency eta",10,10,700,500);
  c2a->SetFillColor(0);
  c2a->SetGrid(1);
  c2a->SetTicky();
  c2a->SetRightMargin(0.03);
  c2a->SetTopMargin(0.02);
  c2a->cd(); 
  
  hi_glbsta_eff_eta->SetXTitle("#eta^{#mu}");
  hi_glbsta_eff_eta->SetYTitle("Efficiency");
  hi_glbsta_eff_eta->SetTitleOffset(1.1,"x");
  hi_glbsta_eff_eta->SetStats(false);

  hi_glbsta_eff_eta->SetMaximum(1.01);
  hi_glbsta_eff_eta->SetMinimum(0.8);
  hi_glbsta_eff_eta->SetLineColor(2);
  hi_glbsta_eff_eta->SetMarkerStyle(22);
  hi_glbsta_eff_eta ->DrawCopy("E");
  c2a->Write();


  TCanvas* c3 = new TCanvas("pur","Purity pt",10,10,700,500);
  c3->SetFillColor(0);
  c3->SetGrid(1);
  c3->SetTicky();
  c3->SetRightMargin(0.03);
  c3->SetTopMargin(0.02);
  c3->cd(); 

  hi_stasim_pur_pt->SetXTitle("p_{T}^{#mu}");
  hi_stasim_pur_pt->SetYTitle("Purity");
  hi_stasim_pur_pt->SetTitleOffset(1.1,"x");
  hi_stasim_pur_pt->SetTitleOffset(1.15,"y");
  hi_stasim_pur_pt->SetMaximum(1.02);
  hi_stasim_pur_pt->SetMinimum(0.5);
  hi_stasim_pur_pt->SetStats(false);

  hi_stasim_pur_pt->SetLineWidth(1);
  hi_glbsim_pur_pt->SetLineWidth(1);
  hi_stasim_pur_pt->SetLineColor(2);
  hi_glbsim_pur_pt->SetLineColor(4);
  hi_stasim_pur_pt->SetLineStyle(1);
  hi_glbsim_pur_pt->SetLineStyle(1);

  //hi_stasim_pur_pt->SetFillColor(5);
  //hi_glbsim_pur_pt->SetFillColor(7);

  hi_stasim_pur_pt->SetMarkerStyle(21);
  hi_glbsim_pur_pt->SetMarkerStyle(26);
  hi_stasim_pur_pt->SetMarkerColor(2);
  hi_glbsim_pur_pt->SetMarkerColor(4);

  hi_stasim_pur_pt ->DrawCopy("PE");
  hi_glbsim_pur_pt ->DrawCopy("PEsame");
  //hi_stasim_pur_pt ->DrawCopy("AxisSame");

  TLegend* legend3 = new TLegend(0.6,0.2,0.8,0.4);
  legend3->SetTextAlign(32);
  legend3->SetTextColor(1);
  legend3->SetTextSize(0.04);
  legend3->AddEntry("hi_stasim_pur_pt",l2string,"pl");
  legend3->AddEntry("hi_glbsim_pur_pt",l3string,"pl");
  legend3 ->Draw();
  c3->Write();
  // 

  TCanvas* c3a = new TCanvas("algPur","Algo Purity pt",10,10,700,500);
  c3a->SetFillColor(0);
  c3a->SetGrid(1);
  c3a->SetTicky();
  c3a->SetRightMargin(0.03);
  c3a->SetTopMargin(0.02);
  c3a->cd(); 
  
  hi_glbsta_pur_pt->SetXTitle("p_{T}^{#mu}");
  hi_glbsta_pur_pt->SetYTitle("Purity");
  hi_glbsta_pur_pt->SetTitleOffset(1.1,"x");
  hi_glbsta_pur_pt->SetStats(false);

  hi_glbsta_pur_pt->SetMaximum(1.01);
  hi_glbsta_pur_pt->SetMinimum(0.8);
  hi_glbsta_pur_pt->SetLineColor(2);
  hi_glbsta_pur_pt->SetMarkerStyle(22);
  hi_glbsta_pur_pt ->DrawCopy("E");
  c3a->Write();


  TCanvas* c4 = new TCanvas("pur_eta","Purity eta",10,10,700,500);
  c4->SetFillColor(0);
  c4->SetGrid(1);
  c4->SetTicky();
  c4->SetRightMargin(0.03);
  c4->SetTopMargin(0.02);
  c4->cd(); 

  hi_stasim_pur_eta->SetXTitle("#eta^{#mu}");
  hi_stasim_pur_eta->SetYTitle("Purity");
  hi_stasim_pur_eta->SetTitleOffset(1.1,"x");
  hi_stasim_pur_eta->SetTitleOffset(1.15,"y");
  hi_stasim_pur_eta->SetMaximum(1.02);
  hi_stasim_pur_eta->SetMinimum(0.5);
  hi_stasim_pur_eta->SetStats(false);

  hi_stasim_pur_eta->SetLineWidth(1);
  hi_glbsim_pur_eta->SetLineWidth(1);
  hi_stasim_pur_eta->SetLineColor(2);
  hi_glbsim_pur_eta->SetLineColor(4);
  hi_stasim_pur_eta->SetLineStyle(1);
  hi_glbsim_pur_eta->SetLineStyle(1);

  //hi_stasim_pur_eta->SetFillColor(5);
  //hi_glbsim_pur_eta->SetFillColor(7);

  hi_stasim_pur_eta->SetMarkerStyle(21);
  hi_glbsim_pur_eta->SetMarkerStyle(26);
  hi_stasim_pur_eta->SetMarkerColor(2);
  hi_glbsim_pur_eta->SetMarkerColor(4);

  hi_stasim_pur_eta ->DrawCopy("PE");
  hi_glbsim_pur_eta ->DrawCopy("PEsame");
  //hi_stasim_pur_eta ->DrawCopy("AxisSame");

  TLegend* legend4 = new TLegend(0.6,0.2,0.8,0.4);
  legend4->SetTextAlign(32);
  legend4->SetTextColor(1);
  legend4->SetTextSize(0.04);
  legend4->AddEntry("hi_stasim_pur_eta",l2string,"pl");
  legend4->AddEntry("hi_glbsim_pur_eta",l3string,"pl");
  legend4 ->Draw();
  c4->Write();
  // 

  TCanvas* c4a = new TCanvas("algPur_eta","Algo Purity eta",10,10,700,500);
  c4a->SetFillColor(0);
  c4a->SetGrid(1);
  c4a->SetTicky();
  c4a->SetRightMargin(0.03);
  c4a->SetTopMargin(0.02);
  c4a->cd(); 
  
  hi_glbsta_pur_eta->SetXTitle("#eta^{#mu}");
  hi_glbsta_pur_eta->SetYTitle("Purity");
  hi_glbsta_pur_eta->SetTitleOffset(1.1,"x");
  hi_glbsta_pur_eta->SetStats(false);

  hi_glbsta_pur_eta->SetMaximum(1.01);
  hi_glbsta_pur_eta->SetMinimum(0.8);
  hi_glbsta_pur_eta->SetLineColor(2);
  hi_glbsta_pur_eta->SetMarkerStyle(22);
  hi_glbsta_pur_eta ->DrawCopy("E");
  c4a->Write();

  hFile->Write();
}

float 
GlobalMuonValidator::calculateDistance(const math::XYZVector& vect1, const math::XYZVector& vect2) {
  float dEta = vect1.eta() - vect2.eta();
  float dPhi = fabs(Geom::Phi<float>(vect1.phi()) - Geom::Phi<float>(vect2.phi()));
  float dPt = sqrt(vect1.perp2()) - sqrt(vect2.perp2());
  //float distance = sqrt(pow(dEta,2) + pow(dPhi,2) + pow(dPt,2) );
  float distance = sqrt(pow(dEta,2) + pow(dPhi,2) );

  return distance;
}

//
// return h1/h2 with recalculated errors
//
TH1F* GlobalMuonValidator::divideErr(TH1F* h1, TH1F* h2, TH1F* hout) {

  hout->Reset();
  hout->Divide(h1,h2,1.,1.,"B");

  for (int i = 0; i <= hout->GetNbinsX()+1; i++ ) {
    Float_t tot   = h2->GetBinContent(i) ;
    Float_t tot_e = h2->GetBinError(i);
    Float_t eff = hout->GetBinContent(i) ;
    Float_t Err = 0.;
    if (tot > 0) Err = tot_e / tot * sqrt( eff* (1-eff) );
    if (eff == 1. || isnan(Err) || !isfinite(Err) ) Err=1.e-3;
    hout->SetBinError(i, Err);
  }
  return hout;
}

//define this as a plug-in
DEFINE_FWK_MODULE(GlobalMuonValidator);
