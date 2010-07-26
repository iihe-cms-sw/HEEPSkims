// -*- C++ -*-
//
// Package:    HEEPSkims
// Class:      HEEPSkim1Ele
// 
/**\class HEEPSkims HEEPSkim1Ele.cc UserCode/HEEPSkims/src/HEEPSkim1Ele.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  local user
//         Created:  Tue May 25 13:12:26 CEST 2010
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCoreFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"

//
// class declaration
//

class HEEPSkim1Ele : public edm::EDFilter {
   public:
      explicit HEEPSkim1Ele(const edm::ParameterSet&);
      ~HEEPSkim1Ele();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      float ptCut_;
      float hoeCut_;
      
      // ----------member data ---------------------------
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
HEEPSkim1Ele::HEEPSkim1Ele(const edm::ParameterSet& iConfig):
  ptCut_(iConfig.getUntrackedParameter<double>("PtCut", 0.)),
  hoeCut_(iConfig.getUntrackedParameter<double>("HoECut", 1.))
{}


HEEPSkim1Ele::~HEEPSkim1Ele()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
HEEPSkim1Ele::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;

  //std::cout<<"run "<<iEvent.id().run()<<" event "<<iEvent.id().event()<<" "<<std::endl;

  edm::Handle<reco::GsfElectronCollection> gsfElectrons;
  iEvent.getByLabel("gsfElectrons", gsfElectrons);
  const reco::GsfElectronCollection* gsfCollection = gsfElectrons.product();
  
//   for(reco::GsfElectronCollection::const_iterator gsfIt = gsfCollection->begin();gsfIt != gsfCollection->end(); gsfIt++)
//     {
//       cout<<"gsfIt->hoe() = "<<gsfIt->hcalOverEcal()<<endl;
//       cout<<"gsfIt->pt()  = "<<gsfIt->pt()<<endl;
//     }
  
  for(reco::GsfElectronCollection::const_iterator gsfIt = gsfCollection->begin();gsfIt != gsfCollection->end(); gsfIt++)
    {
      if ((gsfIt->hcalOverEcal() < hoeCut_) && (gsfIt->caloEnergy()*sin(gsfIt->p4().theta()) > ptCut_)) return true;
    }

  return false;

  

}

// ------------ method called once each job just before starting event loop  ------------
void 
HEEPSkim1Ele::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HEEPSkim1Ele::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HEEPSkim1Ele);
