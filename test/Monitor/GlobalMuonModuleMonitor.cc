/*
 * \file GlobalMuonModuleMonitor.cc
 * 
 * $Date:  $
 * $Revision:  $
 * \author A. Everett - Purdue University
 *
 */

#include "RecoMuon/GlobalMuonProducer/test/Monitor/GlobalMuonModuleMonitor.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Daemon/interface/MonitorDaemon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/ModuleDescription.h"

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

using namespace std;
using namespace edm;

//namespace service {
  GlobalMuonModuleMonitor::GlobalMuonModuleMonitor(const edm::ParameterSet& ps,edm::ActivityRegistry& reg) {
    
    reg.watchPostEndJob(this,&GlobalMuonModuleMonitor::postEndJob);
    reg.watchPreModule(this,&GlobalMuonModuleMonitor::preModule);
    
    outputFileName =ps.getUntrackedParameter<string>("outputFile","test.root");
    parameters = ps;
    
    theFile = new TFile(outputFileName.c_str(),"RECREATE");
  }
    
  GlobalMuonModuleMonitor::~GlobalMuonModuleMonitor() {
    if(theFile) delete theFile;
  }
  
  void GlobalMuonModuleMonitor::preModule(const edm::ModuleDescription& desc){
    curModule = desc.moduleLabel_;
  }
  
  void GlobalMuonModuleMonitor::postEndJob(){
    save(outputFileName);
  }
  
  void 
  GlobalMuonModuleMonitor::book1D(string name, string title, int nchX, 
  			    double lowX, double highX){
    book1D(curModule,name,title,nchX,lowX,highX);
  }
  
  void 
  GlobalMuonModuleMonitor::book1D(std::string level,
			    std::string name, std::string title, int nchX, 
			    double lowX, double highX){
    theFile->cd();

    if ( ! (m1Histos[level].find(name) != m1Histos[level].end() ) ) {
      string histoName = level + "_" + name;
      (m1Histos[level])[name] = new TH1F(histoName.c_str(),title.c_str(),nchX,lowX,highX);
      //(m1Histos.find(level)->second).find(name)->second->SetDirectory(0);  
    }
  }
  
  void 
  GlobalMuonModuleMonitor::book2D(string name, string title, int nchX, 
  			    double lowX, double highX, int nchY,
			    double lowY, double highY)
  {
    book2D(curModule,name,title,nchX,lowX,highX,nchY,lowY,highY);
  }

  void 
  GlobalMuonModuleMonitor::book2D(string level, 
			    string name, string title, int nchX, 
			    double lowX, double highX, int nchY,
			    double lowY, double highY)
  {  
    theFile->cd();
    if ( ! (m2Histos[level].find(name) != m2Histos[level].end() ) ) {
      string histoName = level + "_" + name;
      (m2Histos[level])[name] = new TH2F(name.c_str(),title.c_str(),
					 nchX, lowX, highX,
					 nchY, lowY, highY);
      //(m2Histos.find(level)->second).find(name)->second->SetDirectory(0);  
    }
  }
  
  void GlobalMuonModuleMonitor::save(string filename) 
  {
    if(theFile->IsZombie())
      {
	cerr << " *** Error! Failed creating filename " << filename << endl;
	return;
      }
    theFile->cd();
    //ADAM
    //map<string, map<string,TH1*> >::iterator iter;
    //map<string,TH1*>::iterator iter2;
    //for(iter = m1Histos.begin(); iter != m1Histos.end(); iter++) {     
    //for(iter2 = (*iter).begin(); iter2 != (*iter).end(); iter2++) {
    //(*iter).second->Write("",TObject::kOverwrite);
    //}
    //}
    //end ADAM
    theFile->Write("",TObject::kOverwrite);
  }
  
  void GlobalMuonModuleMonitor::fill1(string level, string name, double a, double b)
  {
    if (m1Histos[level].find(name) != m1Histos[level].end()) {
      (m1Histos.find(level)->second).find(name)->second->Fill(a,b);  
    }
  }
  
  void GlobalMuonModuleMonitor::fill2(string level, string name, double a, double b, double c)
  {
    if (m2Histos[level].find(name) != m2Histos[level].end()) {
      (m2Histos.find(level)->second).find(name)->second->Fill(a,b,c);  
    }
  }

  
  void GlobalMuonModuleMonitor::fill1(string name, double a, double b)
  {
    if (m1Histos[curModule].find(name) != m1Histos[curModule].end()) {
      (m1Histos.find(curModule)->second).find(name)->second->Fill(a,b);
    }
  }

  void GlobalMuonModuleMonitor::fill2(string name, double a, double b, double c)
  {	
    if (m2Histos[curModule].find(name) != m2Histos[curModule].end()) {
      (m2Histos.find(curModule)->second).find(name)->second->Fill(a,b,c);
    }
  }


  
//}//namespace service

