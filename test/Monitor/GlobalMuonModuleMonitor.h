#ifndef RecoMuon_GlobalMuonProducer_GlobalMuonModuleMonitor_H
#define RecoMuon_GlobalMuonProducer_GlobalMuonModuleMonitor_H

/** \class GlobalMuonModuleMonitor
 *
 * Class for Muon inter-module studies.
 *  
 *  $Date: 2006/12/13 20:22:40 $
 *  $Revision: 1.1 $
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


  class GlobalMuonModuleMonitor : public GlobalMuonMonitorInterface {
    
  public:
    
    explicit GlobalMuonModuleMonitor( const edm::ParameterSet& ps,edm::ActivityRegistry& reg);
    
    virtual ~GlobalMuonModuleMonitor();
    
    void postEndJob();
    void preModule(const edm::ModuleDescription& desc);
    
    void book1D(std::string name, std::string title, int nchX, 
		double lowX, double highX);
    void book1D(std::string level, std::string name, std::string title, 
		int nchX, double lowX, double highX);
    void book2D(std::string name, std::string title, int nchX, 
		double lowX, double highX, int nchY,
		double lowY, double highY);
    void book2D(std::string level, std::string name, std::string title,
		int nchX, double lowX, double highX, int nchY, 
		double lowY, double highY);
    
    void save(std::string);
    
    void fill1(std::string, double a, double b=1.);
    void fill1(std::string, std::string, double a, double b=1.);
    void fill2(std::string, double a, double b, double c=1.);
    void fill2(std::string, std::string, double a, double b, double c=1.);
    
  private:
    
    TFile* theFile;
    
    bool debug;
    edm::ParameterSet parameters;
    
    // Monitor Elements
    // <moduleLabel, <histoName, histoPointer> >    
    std::map<std::string, std::map<std::string, TH1*> > m1Histos;
    std::map<std::string, std::map<std::string, TH2*> > m2Histos;
    
    std::string outputFileName;
    std::string curModule;
    
  };
  


#endif
