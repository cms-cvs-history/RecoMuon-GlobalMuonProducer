#ifndef RecoMuon_GlobalMuonProducer_GlobalMuonModuleMonitor_H
#define RecoMuon_GlobalMuonProducer_GlobalMuonModuleMonitor_H

/** \class GlobalMuonModuleMonitor
 *
 * Class for Muon inter-module studies.
 *  
 *  $Date: $
 *  $Revision: $
 *
 * \author A. Everett  - Purdue University
 *
 */

#include "RecoMuon/GlobalMuonProducer/src/GlobalMuonMonitorInterface.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"

#include <fstream>
#include <map>
#include <string>
#include <vector>

class TH1;
class TH2;
class TFile;
class ModuleDescription;

using namespace std;

  class GlobalMuonModuleMonitor : public GlobalMuonMonitorInterface {
    
  public:
    
    explicit GlobalMuonModuleMonitor( const edm::ParameterSet& ps,edm::ActivityRegistry& reg);
    
    virtual ~GlobalMuonModuleMonitor();
    
    void postEndJob();
    void preModule(const edm::ModuleDescription& desc);
    
    void book1D(string name, string title, int nchX, 
		double lowX, double highX);
    void book1D(string level, string name, string title, 
		int nchX, double lowX, double highX);
    void book2D(string name, string title, int nchX, 
		double lowX, double highX, int nchY,
		double lowY, double highY);
    void book2D(string level, string name, string title,
		int nchX, double lowX, double highX, int nchY, 
		double lowY, double highY);
    
    void save(string);
    
    void fill1(string, double a, double b=1.);
    void fill1(string, string, double a, double b=1.);
    void fill2(string, double a, double b, double c=1.);
    void fill2(string, string, double a, double b, double c=1.);
    
  private:
    
    TFile* theFile;
    
    bool debug;
    edm::ParameterSet parameters;
    
    // Monitor Elements
    // <moduleLabel, <histoName, histoPointer> >    
    map<string, map<string, TH1*> > m1Histos;
    map<string, map<string, TH2*> > m2Histos;
    
    string outputFileName;
    string curModule;
    
  };
  


#endif
