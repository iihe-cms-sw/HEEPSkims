// -*- C++ -*-
//
// Package:    GsfCheckerTree
// Class:      GsfCheckerTree
// 
/**\class GsfCheckerTree GsfCheckerTree.cc RecoEgamma/GsfCheckerTree/src/GsfCheckerTree.cc

Description: <one line class summary>

Implementation:
<Notes on implementation>
*/
//
// Original Author:  Charaf Otman
//         Created:  Thu Jan 17 14:41:56 CET 2008
// $Id: GsfCheckerTree.cc,v 1.4 2011/06/24 09:00:57 treis Exp $
//
//

//TRIGGER ARNAUD
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"

#include "HLTrigger/HLTanalyzers/interface/HLTrigReport.h"


#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "UserCode/HEEPSkims/interface/GsfCheckerTree.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/EgAmbiguityTools.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronClassification.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronMomentumCorrector.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronEnergyCorrector.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"


//For extrapolation to calo surface
#include "PhysicsTools/IsolationAlgos/interface/PropagateToCal.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#define PI 3.141592654
#define TWOPI 6.283185308

using namespace std;
using namespace reco;


//-----------------------------
//Added from Vincent
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/EDMException.h"


#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "Calibration/EcalCalibAlgos/interface/ElectronRecalibSuperClusterAssociator.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackExtraFwd.h"

#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimG4Core/Generators/interface/HepMCParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCoreFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "Calibration/EcalCalibAlgos/interface/EcalEleCalibLooper.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "Calibration/EcalCalibAlgos/interface/Pi0FixedMassWindowCalibration.h"

#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include "DataFormats/METReco/interface/METCollection.h"

#include "DataFormats/JetReco/interface/CaloJet.h"
//#include "DataFormats/JetReco/interface/CaloJetfwd.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h" 

#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "DataFormats/MuonReco/interface/Muon.h"
// L1 bit
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"

#include "DataFormats/L1Trigger/interface/L1ParticleMap.h"
#include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMap.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
//END L1 bit


// HLT
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/HLTenums.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Provenance/interface/ParameterSetID.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//GEN
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include <iostream>

using namespace edm;
using namespace std;
using namespace reco;

//-----------------------------




HepMC::FourVector momelec,momposi;
HepMC::FourVector unstablemomelec,unstablemomposi;
HepMC::FourVector momboson;

//this is what i added

int fsrposiphotonsize = 0;
int fsrelecphotonsize = 0;

HepMC::FourVector elecaddposigsfgsf = 0;
HepMC::FourVector unstelecaddposigsfgsf = 0;

unsigned ambSortingStrategy_ = 1;
unsigned ambClustersOverlapStrategy_ = 1;

//Method to sort the gsf electrons
bool gsfEtGreater(const reco::GsfElectron &gsf1,const reco::GsfElectron &gsf2)
{
  float et1 = gsf1.caloEnergy()*sin(gsf1.p4().theta());
  float et2 = gsf2.caloEnergy()*sin(gsf2.p4().theta());
  return (et1 > et2);
}

bool scEGreater(const reco::SuperCluster *sc1,const reco::SuperCluster *sc2){return ((sc1->energy()+sc1->preshowerEnergy()) > (sc2->energy()+sc2->preshowerEnergy()));}


bool refScEGreater(reco::SuperClusterRef sc1,reco::SuperClusterRef sc2){return ((sc1->energy()+sc1->preshowerEnergy()) > (sc2->energy()+sc2->preshowerEnergy()));}


// utilities for constructor
float normalized_dphi( float dphi )
{
  if (fabs(dphi)>CLHEP::pi) return (dphi<0?CLHEP::twopi+dphi:dphi-CLHEP::twopi) ;
  else return dphi ;
}


float etacorr(float eta, float pvz, float scz){
  return asinh(sinh(eta)*(1.-pvz/scz));
}


void GsfCheckerTree::setupES(const edm::EventSetup& es) {
}




GsfCheckerTree::GsfCheckerTree(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  eventcounter = 0;

  log_ = iConfig.getParameter<unsigned int>( "logEvents"  );
  src_ = iConfig.getParameter<edm::InputTag>( "src"  );
  hlTriggerResults_ = iConfig.getParameter<edm::InputTag> ("TriggerResultsTag");
  usegendata_ = iConfig.getParameter<bool> ("usegendata");
}


GsfCheckerTree::~GsfCheckerTree()
{

  cout<<"GsfCheckerTree::~GsfCheckerTree"<<endl; 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}





// ------------ method called to for each event  ------------
void
GsfCheckerTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;

  int debugcounter = 0;
  
  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  //setupES(iSetup);

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  eventcounter++;

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  //init all variables
  
  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  //Run and event number
  runnumber = iEvent.id().run();
  eventnumber = iEvent.id().event();
  luminosityBlock = iEvent.id().luminosityBlock(); 

  HLT_Photon20_CaloIdVL_IsoL_v1 = -10;
  HLT_DoublePhoton33_vx = -10;
  HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2 = -10;
  HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2 = -10;
  HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2 = -10;

  HLT_Ele10_SW_EleId_L1R = -10;
  HLT_Ele10_SW_L1R = -10;
  HLT_Ele15_LW_L1R = -10;
  HLT_Ele15_SW_L1R = -10;
  HLT_Ele15_SW_EleId_L1R = -10;
  HLT_Ele15_SW_CaloEleId_L1R  = -10;
  HLT_Ele15_SiStrip_L1R = -10;
  HLT_Ele20_SW_L1R = -10;
  HLT_Ele20_SiStrip_L1R = -10;
  HLT_Ele25_SW_L1R = -10;
  HLT_DoubleEle4_SW_eeRes_L1R = -10;
  HLT_DoubleEle10_SW_L1R = -10;
  HLT_Mu15_Photon20_CaloIdL = -10; //VINCENT

  pthat = -5000.;
  alphaqcd = -5000.;
  alphaqed = -5000.;
  qscale = -5000.;
  processid = -5000;
  weight = -5000.;

  debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;
 
  //Skim Laurent 
  //Final GSF Electron collection
  edm::Handle<reco::GsfElectronCollection> pGsfElectrons;
  bool gsfisvalid = iEvent.getByLabel("gsfElectrons","",pGsfElectrons);
  reco::GsfElectronCollection gsfelectrons(pGsfElectrons->begin(),pGsfElectrons->end());

  //sort all the GSF by transverse energy
  //std::cout<<"sorting the elements by transverse energy"<<std::endl;
  std::sort(gsfelectrons.begin(),gsfelectrons.end(),gsfEtGreater);


  float gsfPtMax=0;
  float gsfPtSecondMax=0; 
  reco::GsfElectronCollection::const_iterator gsfiterbis = gsfelectrons.begin();
  int counter =0; 
  for(;gsfiterbis!=gsfelectrons.end();gsfiterbis++)
    {
      counter ++; 
      if (counter ==1) {
	gsfPtMax = gsfiterbis->caloEnergy()*sin(gsfiterbis->p4().theta());
      }
      if (counter ==2) {
	gsfPtSecondMax = gsfiterbis->caloEnergy()*sin(gsfiterbis->p4().theta());
      } 
    
    }
    if (gsfPtMax>30 && gsfPtSecondMax >30){

  if(usegendata_) {
    edm::Handle<GenEventInfoProduct> GenInfoHandle;
    bool genevtinfovalid = iEvent.getByLabel("generator",GenInfoHandle);

    if(genevtinfovalid) {
      pthat = GenInfoHandle->hasBinningValues() ? (GenInfoHandle->binningValues())[0] : 0.0 ;
      alphaqcd = GenInfoHandle->alphaQCD();
      alphaqed = GenInfoHandle->alphaQED();
      qscale = GenInfoHandle->qScale();
      processid = GenInfoHandle->signalProcessID();
      weight = GenInfoHandle->weight();
    }
  }

  //genelec
  genelechassc = false;
  genelechasgsf = false;

  //genposi
  genposihassc = false;
  genposihasgsf = false;

  //FSR variables
  fsrposiphotonsize = -3;
  fsrelecphotonsize = -3;
  
  for(int k=0;k<10;k++){
    energyfsrelec[k] = -5000.;
    etfsrelec[k] = -5000.;
    etafsrelec[k] = -5000.;
    phifsrelec[k] = -5000.;
  
    energyfsrposi[k] = -5000.;
    etfsrposi[k] = -5000.;
    etafsrposi[k] = -5000.;
    phifsrposi[k] = -5000.;
  }


  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;


  //Supercluster variables
  //fill e,et,eta,phi,charge for every SC in the event
  for(int i=0;i<100;i++) {
    scgsfmatched[i] = -5000.;
    scenergy[i] = -5000.;
    sceta[i] = -5000.;
    scetacorr[i] = -5000.;
    sctheta[i] = -5000.;
    scthetacorr[i] = -5000.;
    scet[i] = -5000.;
    scphi[i] = -5000.;
    scpx[i] = -5000.;
    scpy[i] = -5000.;
    scpz[i] = -5000.;
    scx[i] = -5000.;
    scy[i] = -5000.;
    scz[i] = -5000.;
  }
  

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  scindexforgenelec = -3;
  scindexforgenposi = -3;
  scsize = -5000;
  
  //beam spot variables
  sigmaZ = -5000.;
  sigmaZ0Error = -5000.;
  sq = -5000.; 
  bsposx = -5000.;
  bsposy = -5000.;
  bsposz = -5000.;
   
  // //Primary vertex variables
  //   pvx = -5000.;
  //   pvy = -5000.;  
  //   pvz = -5000.;

  //generated variables for the tree (after FSR)
  genelec_e_var = -5000.;
  genelec_eta_var = -5000.;
  genelec_phi_var = -5000.;
  genelec_et_var = -5000.;
  genposi_e_var = -5000.;
  genposi_eta_var = -5000.;
  genposi_phi_var = -5000.;
  genposi_et_var = -5000.;
  genelec_hassc_var = 0;  
  genposi_hassc_var = 0;  

  //generated variables for the tree (before FSR)
  unstablegenelec_e_var = -5000.;
  unstablegenelec_eta_var = -5000.;
  unstablegenelec_phi_var = -5000.;
  unstablegenelec_et_var = -5000.;
  unstablegenposi_e_var = -5000.;
  unstablegenposi_eta_var = -5000.;
  unstablegenposi_phi_var = -5000.;
  unstablegenposi_et_var = -5000.;
  
  genboson_m_var = -5000.;
  genboson_eta_var = -5000.;
  genboson_phi_var = -5000.;
  genboson_e_var = -5000.;
  genboson_et_var = -5000.;
  genboson_ez_var = -5000.;
  genboson_p_var = -5000.;
  genboson_pt_var = -5000.;
  genboson_pz_var = -5000.;

  x1quark = -5000.;
  x2quark = -5000;;

  //supercluster matching gen variables for the tree
  scelecenergy = -5000.;
  sceleceta = -5000.;
  scelecphi = -5000.;
  scelecgsfmatched = -5000.;
  scposienergy = -5000.;
  scposieta = -5000.;
  scposiphi = -5000.;
  scposigsfmatched = -5000.;


  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  //ELE --- GSF variables
  gsf_size = -3;
  for (int i=0;i<100;i++){
    gsf_p[i]=-1000;
    gsf_e[i]=-1000;
    gsf_pt[i]=-1000;
    gsf_class[i]=-1000;
    gsf_e2x5overe5x5[i]=-1000.;
    gsf_e1x5overe5x5[i]=-1000.;
    gsf_eta[i]=-1000;
    gsf_phi[i]=-1000;
    gsf_px[i]=-1000;
    gsf_py[i]=-1000;
    gsf_pz[i]=-1000;

    gsf_deltaeta[i]=-1000.;
    gsf_deltaphi[i]=-1000.;
    gsf_hovere[i]=-1000.;
	  
    gsf_trackiso[i]=-1000.;
    gsf_ecaliso[i]=-1000.;
    gsf_hcaliso1[i]=-1000.;
    gsf_hcaliso2[i]=-1000.;
	  
    gsf_charge[i]=-1000;
    gsf_sigmaetaeta[i]=-1000.;
    gsf_sigmaIetaIeta[i]=-1000.;
    gsf_isecaldriven[i]=-1000.;
    gsf_istrackerdriven[i]=-1000.;

    gsfsc_e[i]=-1000;
    gsfsc_pt[i]=-1000;
    gsfsc_eta[i]=-1000;
    gsfsc_phi[i]=-1000;
    gsfsc_px[i]=-1000;
    gsfsc_py[i]=-1000;
    gsfsc_pz[i]=-1000;

    gsf_gsfet[i]=-1000;

    gsfpass_ET[i] = 0; 
    gsfpass_PT[i] = 0; 
    gsfpass_DETETA[i] = 0; 
    gsfpass_CRACK[i] = 0; 
    gsfpass_DETAIN[i] = 0; 
    gsfpass_DPHIIN[i] = 0; 
    gsfpass_HADEM[i] = 0; 
    gsfpass_SIGMAIETAIETA [i] = 0;
    gsfpass_E2X5OVER5X5[i] = 0; 
    gsfpass_ISOLEMHADDEPTH1[i] = 0;
    gsfpass_ISOLHADDEPTH2[i] = 0; 
    gsfpass_ISOLPTTRKS[i] = 0; 
    gsfpass_ECALDRIVEN[i] = 0; 
    gsfpass_INVALID[i] = 0;
    gsfpass_NOMISSINGHITS[i] = 0;

    gsfpass_HEEP[i] = 0;

    gsfpass_ID[i] = 0;
    gsfpass_ISO[i] = 0;

    scindexforgsf[i] = -3;

    gsf_theta[i] = -1.;
    gsf_isEB[i] = -1;
    gsf_isEE[i] = -1;
    gsf_deltaEtaATcalo[i] = -1.;
    gsf_deltaPhiATcalo[i] = -1.; 
    gsf_ecalEnergy[i] = -1.;
    gsf_eOVERp[i] = -1.;
    gsf_dxy[i] = -1.;
    gsf_vz[i] = -1.;
    gsf_nHits[i] = -1;
    gsf_nLostInnerHits[i] = -1;
    gsf_nLostOuterHits[i] = -1;
    gsf_fBrem[i] = -1.;
    //gsf_e1OVERe9[i] = -1.;
    gsf_eMax[i] = -1.;
    gsf_SwissCross[i] = -1.;
    gsf_e1x3[i] = -1.;
    gsf_e3x1[i] = -1.;
    gsf_e2x2[i] = -1.;
    gsf_e3x2[i] = -1.;
    gsf_e3x3[i] = -1.;
    gsf_e4x4[i] = -1.;
    gsf_e2x5Right[i] = -1.;
    gsf_e2x5Left[i] = -1.;  
    gsf_e2x5Top[i] = -1.;  
    gsf_e2x5Bottom[i] = -1.;
    gsf_e2x5Max[i] = -1.;
    gsf_eLeft[i] = -1.;
    gsf_eRight[i] = -1.;
    gsf_eTop[i] = -1.;
    gsf_eBottom[i] = -1.;
    gsf_e2nd[i] = -1.;    


    //Charge info
    scpixcharge[i] = -3;
    ctfcharge[i] = -3;
    gsfcharge[i] = -3;
    gsfctfscpixconsistent[i] = false;
    gsfscpixconsistent[i] = false;
    gsfctfconsistent[i] = false;

  }
  gsfindexforgenelec = -3;
  gsfindexforgenposi = -3;


  debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  //------------------------------
  //Added from Vincent
  for (int i=0;i<300;i++) HLTriggers[i]  = -10;

  nJetsAKT_pt15 = -1;
  //  nJetsIC5_pt15 = -1;
  calomet = -1.;
  met = -1.;

  //IC5
//   for (unsigned int i = 0 ; i< 100 ; i++){   
//     jetIC5_pt[i] = -1;
//     jetIC5_eta[i] = -1;
//     jetIC5_phi[i] = -1;
//     jetIC5_em[i] = -1;
//   }

  for (unsigned int i = 0 ; i< 50 ; i++){   
    jetAKT_pt[i] = -1;
    jetAKT_eta[i] = -1;
    jetAKT_phi[i] = -1;
    jetAKT_em[i] = -1;
  }

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;


  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  //MUONS
  muon_size = -3;
  for (unsigned int i = 0 ; i< 100 ; i++){  
    muon_pt[i] = -1.;
    muon_ptError[i] = -1.;
    muon_eta[i] = -1.;
    muon_etaError[i] = -1.;
    muon_phi[i] = -1.;
    muon_phiError[i] = -1.;
    muon_theta[i] = -1.;
    muon_thetaError[i] = -1.;

    muon_outerPt[i] = -1.;
    muon_outerEta[i] = -1.;
    muon_outerPhi[i] = -1.;
    muon_outerTheta[i] = -1.;

    muon_px[i] = -1.;
    muon_py[i] = -1.;
    muon_pz[i] = -1.;

    muon_charge[i] = -1;
    muon_nhitspixel[i] = -1;
    muon_nhitstrack[i] = -1;
    muon_nhitsmuons[i] = -1;
    muon_nhitstotal[i] = -1;
    muon_nSegmentMatch[i] = -1;

    muon_isTrackerMuon[i] = false;

    muon_chi2[i] = -1.;
    muon_ndof[i] = -1.;
    muon_normChi2[i] = -1.;

    muon_d0[i] = -1.;
    muon_d0Error[i] = -1.;
    muon_dz_cmsCenter[i] = -1.;
    muon_dz_beamSpot[i] = -1.;
    muon_dz_firstPVtx[i] = -1.;
    muon_dzError[i] = -1.;
    muon_dxy_cmsCenter[i] = -1.;
    muon_dxy_beamSpot[i] = -1.;
    muon_dxy_firstPVtx[i] = -1.;
    muon_dxyError[i] = -1.; 

    muon_trackIso03[i] = -1.; 
    muon_trackIso05[i] = -1.; 
    muon_trackIso03_ptInVeto[i] = -1.; 
    muon_trackIso05_ptInVeto[i] = -1.; 

    muon_emIso03[i] = -1.; 
    muon_emIso05[i] = -1.; 
    muon_emIso03_ptInVeto[i] = -1.; 
    muon_emIso05_ptInVeto[i] = -1.; 

    muon_hadIso03[i] = -1.; 
    muon_hadIso05[i] = -1.; 
    muon_hadIso03_ptInVeto[i] = -1.; 
    muon_hadIso05_ptInVeto[i] = -1.; 


    muon_innerPosx[i] = -1.; 
    muon_innerPosy[i] = -1.; 
    muon_innerPosz[i] = -1.; 
  }

  //------------------------------


  if(usegendata_){
    datagenerated(iEvent);

    //debugcounter++;
    //cout<<"debug "<<debugcounter<<endl;

    genelec_e_var = momelec.e();
    genelec_eta_var = momelec.eta(); 
    genelec_phi_var = momelec.phi(); 
    genelec_et_var = momelec.e()*sin(momelec.theta()); 
    genposi_e_var = momposi.e();
    genposi_eta_var = momposi.eta(); 
    genposi_phi_var = momposi.phi(); 
    genposi_et_var = momposi.e()*sin(momposi.theta()); 

    unstablegenelec_e_var = unstablemomelec.e();
    unstablegenelec_eta_var = unstablemomelec.eta(); 
    unstablegenelec_phi_var = unstablemomelec.phi(); 
    unstablegenelec_et_var = unstablemomelec.e()*sin(unstablemomelec.theta()); 
    unstablegenposi_e_var = unstablemomposi.e();
    unstablegenposi_eta_var = unstablemomposi.eta(); 
    unstablegenposi_phi_var = unstablemomposi.phi(); 
    unstablegenposi_et_var = unstablemomposi.e()*sin(unstablemomposi.theta()); 

    genboson_m_var = momboson.m();

    genboson_eta_var = momboson.eta(); 
    genboson_phi_var = momboson.phi(); 

    genboson_e_var = momboson.e();
    genboson_et_var = momboson.e()*sin(momboson.theta());
    genboson_ez_var = momboson.e()*cos(momboson.theta()); 

    genboson_p_var = sqrt(momboson.px()*momboson.px()+momboson.py()*momboson.py()+momboson.pz()*momboson.pz());
    genboson_pt_var = sqrt(momboson.px()*momboson.px()+momboson.py()*momboson.py());
    genboson_pz_var = momboson.pz();

    x1quark = (genboson_m_var*genboson_m_var)/(7000.*(genboson_pz_var + sqrt(genboson_pz_var*genboson_pz_var+genboson_m_var*genboson_m_var)));
    x2quark = (genboson_pz_var + sqrt(genboson_pz_var*genboson_pz_var+genboson_m_var*genboson_m_var))/7000.;
  }


  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;




  edm::Handle<TriggerResults> hltResults;
  bool hltriggerresultisvalid = iEvent.getByLabel(hlTriggerResults_, hltResults);


  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  edm::Handle<edm::TriggerResults> hltTriggerResultHandle;
  iEvent.getByLabel(hlTriggerResults_, hltTriggerResultHandle);
 

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  hltCount = 0;   
  if(!hltTriggerResultHandle.isValid()) {
    std::cout << "invalid handle for HLT TriggerResults" << std::endl;
  } else {
    hltCount = hltTriggerResultHandle->size();
    //cout<<"hltCount = "<<hltCount<<endl;
      
    for(int i = 0 ; i < hltCount ; i++) {
      //      aHLTResults[i] = hltTriggerResultHandle->accept(i);
      //cout<<" i , hltTriggerResultHandle->accept(i) = "<<i<<", "<<hltTriggerResultHandle->accept(i)<<endl;
      //prescaling
      //if ( (i!= 1) && (i != 2) && (i != 11) && (i != 18) && (i != 88) && (i != 95) && (i != 96) && (i != 183))
      {
	//if ((i<169) || (i>174))
	{
	  if (hltTriggerResultHandle->accept(i)) HLTriggers[i] = i;
	  //if (hltTriggerResultHandle->accept(i)) HLTriggers[i-1] = i-1;
	}
      }
    }
  } // end HLT
  

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  edm::TriggerResults  triggerResults();
  edm::ParameterSetID psetid_;



  debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  bool gen_bool = false;
  //
//   edm::Handle<EcalRecHitCollection> pEcalRecHitBarrelCollection;
//   //iEvent.getByLabel(ecalHitsProducer_, barrelHits_, pEcalRecHitBarrelCollection);
//   iEvent.getByLabel("ecalRecHit", "EcalRecHitsEB", pEcalRecHitBarrelCollection);
//   const EcalRecHitCollection* ecalRecHitBarrelCollection = pEcalRecHitBarrelCollection.product();


  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  edm::Handle<edm::View<reco::GenParticle> > src; 
  iEvent.getByLabel(src_, src);
  

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  // L1 BITS
  edm::Handle< L1GlobalTriggerReadoutRecord > gtReadoutRecord;
  iEvent.getByLabel( edm::InputTag("gtDigis"), gtReadoutRecord);
  

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = gtReadoutRecord->technicalTriggerWord();
  L1trigger_size = technicalTriggerWordBeforeMask.size();
 

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  for(unsigned int i = 0;i<technicalTriggerWordBeforeMask.size();i++){
    bool bit = technicalTriggerWordBeforeMask.at(i);
    if (bit == 1) L1trigger_bool[i] = 1;   
    if (bit == 0) L1trigger_bool[i] = 0;   
  }

 //  //GEN PARTICLES FOR TTBAR
  
//   int ele_indic1=-1; 
  
//   edm::Handle<HepMCProduct> EvtHandle ;
//   iEvent.getByLabel( "generator", EvtHandle ) ;
//   const HepMC::GenEvent* genEvt = EvtHandle->GetEvent() ;
//   HepMC::GenVertex* GravitonDecVtx = 0 ;
	 
//   for ( HepMC::GenEvent::vertex_const_iterator vit=genEvt->vertices_begin(); vit!=genEvt->vertices_end(); vit++ ) {
    
//     for ( HepMC::GenVertex::particles_out_const_iterator pout=(*vit)->particles_out_const_begin();pout!=(*vit)->particles_out_const_end(); pout++ ){
// 	  //if (  (  ( (*pout)->pdg_id() == 32 ) || ( (*pout)->pdg_id() == 23 ) || ( (*pout)->pdg_id() == 22 ) ) && ( (*pout)->status() == 3))	     {	    
// 	  //if ( ((*pout)->pdg_id() == 23) && ((*pout)->status() == 3) ) {
      
//     //TTbar   
//     //ele1
//     if (fabs( (*pout)->pdg_id() ) == 11 
// 	&& ele_indic1==-1){

//       for(std::set< HepMC::GenParticle*>::const_iterator iter = TopDecVtx->particles_out_const_begin();
// 	 iter != TopDecVtx->particles_out_const_end();iter++) { TopChildren.push_back(*iter); }
      
//       const reco::Candidate * w = (*pout)->particles_out_const_begin(); //faire une boucle
//       //const reco::Candidate * w = (*pout)->mother();
//       if (w != NULL && abs((*w).pdgId()) == 24) { // W
	
// 	const reco::Candidate * t = (*w).mother();
// 	if (t != NULL && abs((*t).pdgId()) == 6) { // top	
// 	  //ele_indic1=i;	
// 	}
//       }
//     }
    
//     // //ele2
// //     if (fabs(genParticles[i].pdgId()) == 11
// // 	&& ele_indic1!=-1 && ele_indic1!=i && ele_indic2==-1){
	
// //       const reco::Candidate * w = genParticles[i].mother();
// //       if (w != NULL && abs((*w).pdgId()) == 24) { // W
	
// // 	const reco::Candidate * t = (*w).mother();
// // 	if (t != NULL && abs((*t).pdgId()) == 6) { // top	    
// // 	  ele_indic2=i;	
// // 	}      
// //       }
// //     }
//     }
//   }
    
  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  //physics declared
  L1GlobalTriggerReadoutRecord const* gtrr = gtReadoutRecord.product();
  L1GtFdlWord fdlWord = gtrr->gtFdlWord();
  //cout << "phys decl. bit=" << fdlWord.physicsDeclared() << endl;
  if (fdlWord.physicsDeclared() == 1) PhysDecl_bool=1;
  else PhysDecl_bool=0;
   
  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  //  //Get MUON
  //   edm::Handle<reco::TrackCollection> muonCollection;
  //   iEvent.getByLabel("globalMuons", muonCollection);
  //   const reco::TrackCollection* muons = muonCollection.product();
  //   //const reco::MuonCollection recoMu = dynamic_cast < const reco::MuonCollection * >(muons);

  // //Get MUON
  edm::Handle<reco::MuonCollection> muonCollection;
  iEvent.getByLabel("muons",muonCollection);
  const reco::MuonCollection* muons = muonCollection.product();

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  bool caloantiktjetisvalid = false;
  edm::Handle<CaloJetCollection> pCaloAntiKtJets;
  caloantiktjetisvalid = iEvent.getByLabel("ak5CaloJets", pCaloAntiKtJets);//Laurent
     const CaloJetCollection *caloAntiKtJets  = pCaloAntiKtJets.product();//Laurent
  

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

// RECO -> AOD
//   bool caloiterconejetisvalid = false;
//   edm::Handle<CaloJetCollection> pCaloIterConeJets;
//   caloiterconejetisvalid = iEvent.getByLabel("iterativeCone5CaloJets", pCaloIterConeJets);
//   const CaloJetCollection *caloIterConeJets  = pCaloIterConeJets.product();
  

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  edm::Handle<CaloMETCollection> pCaloMET;
  bool calometisvalid = iEvent.getByLabel("met", pCaloMET);
  const CaloMETCollection *caloMET  = pCaloMET.product();

  //debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;

  
  edm::Handle<METCollection> pMET;
  bool metisvalid = iEvent.getByLabel("htMetKT4", pMET);
  const METCollection *MET  = pMET.product();

// RECO -> AOD
//   //RECHITS
//   //BARREL
//   const EBRecHitCollection* barrelHitsCollection = 0;
//   edm::Handle<EBRecHitCollection> barrelRecHitsHandle;
//   iEvent.getByLabel("ecalRecHit","EcalRecHitsEB",barrelRecHitsHandle);
//   barrelHitsCollection = barrelRecHitsHandle.product();
//   if(!barrelRecHitsHandle.isValid()) {
//     edm::LogError("reading") << "[CmsShowAnalyzer] barrel rec hits not found";
//     std::cout<<"[CmsShowAnalyzer] barrel rec hits not found"<<std::endl;
//     return;
//   }

//   //debugcounter++;
//   //cout<<"debug "<<debugcounter<<endl;

//   //ENDCAP
//   const EERecHitCollection* endcapHitsCollection = 0;
//   edm::Handle<EERecHitCollection> endcapRecHitsHandle;
//   iEvent.getByLabel("ecalRecHit","EcalRecHitsEE",endcapRecHitsHandle);
//   endcapHitsCollection = endcapRecHitsHandle.product();
//   if (!endcapRecHitsHandle.isValid()) {  
//     edm::LogError("reading") << "[CmsShowAnalyzer] endcap rec hits not found"; 
//     std::cout<<"[CmsShowAnalyzer] endcap rec hits not found"<<std::endl;
//     return;
//   }

   // Triggers  ARNAUD

  // get hold of TriggerResults
  Handle<TriggerResults> HLTR;
  iEvent.getByLabel(hlTriggerResults_,HLTR);
  if (HLTR.isValid()) {
    //cout<<"trigger valid"<<endl;
    if (HLTR->wasrun()) nWasRun_++;
    const bool accept(HLTR->accept());
    LogDebug("HLTrigReport") << "HL TriggerResults decision: " << accept;
    if (accept) ++nAccept_;
    if (HLTR->error() ) nErrors_++;
  } else {
    LogDebug("HLTrigReport") << "HL TriggerResults with label ["+hlTriggerResults_.encode()+"] not found!";
    nErrors_++;
    //cout<<"trigger pas valide"<<endl;
    //return;
  }

  // decision for each HL algorithm
  const unsigned int n(hlNames_.size());
  //cout<<"n = "<<n<<endl;
  for (unsigned int i=0; i!=n; ++i) {
    //cout<<"hlNames(i) = "<<hlNames_.at(i)<<endl;
    //hlNamesTab[i] = hlNames_.at(i)+"\0";
    //hlNamesTab[i] = (hlNames_.at(i)).c_str();
    //cout<<"hlNamesTab[i] = "<<hlNamesTab[i]<<endl;
    if ((hlNames_.at(i)== "HLT_Ele10_SW_EleId_L1R") && (HLTR->accept(i) == 0)) HLT_Ele10_SW_EleId_L1R = 0;
    if ((hlNames_.at(i)== "HLT_Ele10_SW_EleId_L1R") && (HLTR->accept(i) == 1)) HLT_Ele10_SW_EleId_L1R = 1;
    if ((hlNames_.at(i)== "HLT_Ele10_SW_L1R") && (HLTR->accept(i) == 0)) HLT_Ele10_SW_L1R = 0;
    if ((hlNames_.at(i)== "HLT_Ele10_SW_L1R") && (HLTR->accept(i) == 1)) HLT_Ele10_SW_L1R = 1;
    if ((hlNames_.at(i)== "HLT_Ele15_LW_L1R") && (HLTR->accept(i) == 0)) HLT_Ele15_LW_L1R = 0;
    if ((hlNames_.at(i)== "HLT_Ele15_LW_L1R") && (HLTR->accept(i) == 1)) HLT_Ele15_LW_L1R = 1;
    if ((hlNames_.at(i)== "HLT_Ele15_SW_L1R") && (HLTR->accept(i) == 0)) HLT_Ele15_SW_L1R = 0;
    if ((hlNames_.at(i)== "HLT_Ele15_SW_L1R") && (HLTR->accept(i) == 1)) HLT_Ele15_SW_L1R  = 1;
    if ((hlNames_.at(i)== "HLT_Ele15_SW_EleId_L1R") && (HLTR->accept(i) == 0)) HLT_Ele15_SW_EleId_L1R  = 0;
    if ((hlNames_.at(i)== "HLT_Ele15_SW_EleId_L1R") && (HLTR->accept(i) == 1)) HLT_Ele15_SW_EleId_L1R  = 1;
    if ((hlNames_.at(i)== "HLT_Ele15_SW_CaloEleId_L1R") && (HLTR->accept(i) == 0))  HLT_Ele15_SW_CaloEleId_L1R  = 0;
    if ((hlNames_.at(i)== "HLT_Ele15_SW_CaloEleId_L1R") && (HLTR->accept(i) == 1)) HLT_Ele15_SW_CaloEleId_L1R = 1;
    if ((hlNames_.at(i)== "HLT_Ele15_SiStrip_L1R") && (HLTR->accept(i) == 0)) HLT_Ele15_SiStrip_L1R = 0;
    if ((hlNames_.at(i)== "HLT_Ele15_SiStrip_L1R") && (HLTR->accept(i) == 1)) HLT_Ele15_SiStrip_L1R = 1;
    if ((hlNames_.at(i)== "HLT_Ele20_SW_L1R") && (HLTR->accept(i) == 0)) HLT_Ele20_SW_L1R  = 0;
    if ((hlNames_.at(i)== "HLT_Ele20_SW_L1R") && (HLTR->accept(i) == 1)) HLT_Ele20_SW_L1R = 1;
    if ((hlNames_.at(i)== "HLT_Ele20_SiStrip_L1R") && (HLTR->accept(i) == 0)) HLT_Ele20_SiStrip_L1R = 0;
    if ((hlNames_.at(i)== "HLT_Ele20_SiStrip_L1R") && (HLTR->accept(i) == 1)) HLT_Ele20_SiStrip_L1R = 1;
    if ((hlNames_.at(i)== "HLT_Ele25_SW_L1R") && (HLTR->accept(i) == 0)) HLT_Ele25_SW_L1R = 0;
    if ((hlNames_.at(i)== "HLT_Ele25_SW_L1R") && (HLTR->accept(i) == 1)) HLT_Ele25_SW_L1R  = 1;
    if ((hlNames_.at(i)== "HLT_DoubleEle4_SW_eeRes_L1R") && (HLTR->accept(i) == 0)) HLT_DoubleEle4_SW_eeRes_L1R = 0;
    if ((hlNames_.at(i)== "HLT_DoubleEle4_SW_eeRes_L1R") && (HLTR->accept(i) == 1)) HLT_DoubleEle4_SW_eeRes_L1R = 1;
    if ((hlNames_.at(i)== "HLT_DoubleEle10_SW_L1R") && (HLTR->accept(i) == 0)) HLT_DoubleEle10_SW_L1R = 0;
    if ((hlNames_.at(i)== "HLT_DoubleEle10_SW_L1R") && (HLTR->accept(i) == 1)) HLT_DoubleEle10_SW_L1R  = 1;

    if ((hlNames_.at(i)== "HLT_Photon20_CaloIdVL_IsoL_v1") && (HLTR->accept(i) == 1)) HLT_Photon20_CaloIdVL_IsoL_v1  = 1;
    if (((hlNames_.at(i)== "HLT_DoublePhoton33_v1") || (hlNames_.at(i)== "HLT_DoublePhoton33_v2")) && (HLTR->accept(i) == 1))  HLT_DoublePhoton33_vx = 1;
    if ((hlNames_.at(i)== "HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2") && (HLTR->accept(i) == 1)) HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2  = 1;
    if ((hlNames_.at(i)== "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2") && (HLTR->accept(i) == 1)) HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2  = 1;
    if ((hlNames_.at(i)== "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2") && (HLTR->accept(i) == 1)) HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2  = 1;

    //FROM VINCENT
    if (((hlNames_.at(i)== "HLT_Mu15_Photon20_CaloIdL_v1") || (hlNames_.at(i)== "HLT_Mu15_Photon20_CaloIdL_v2") || (hlNames_.at(i)== "HLT_Mu15_Photon20_CaloIdL_v3") || (hlNames_.at(i)== "HLT_Mu15_Photon20_CaloIdL_v4")  
	 || (hlNames_.at(i)== "HLT_Mu15_Photon20_CaloIdL_v5")  || (hlNames_.at(i)== "HLT_Mu15_Photon20_CaloIdL_v6"))&& (HLTR->accept(i) == 0)) HLT_Mu15_Photon20_CaloIdL = 0;
    if (((hlNames_.at(i)== "HLT_Mu15_Photon20_CaloIdL_v1") || (hlNames_.at(i)== "HLT_Mu15_Photon20_CaloIdL_v2") || (hlNames_.at(i)== "HLT_Mu15_Photon20_CaloIdL_v3") || (hlNames_.at(i)== "HLT_Mu15_Photon20_CaloIdL_v4")  
	 || (hlNames_.at(i)== "HLT_Mu15_Photon20_CaloIdL_v5")  || (hlNames_.at(i)== "HLT_Mu15_Photon20_CaloIdL_v6"))&& (HLTR->accept(i) == 1)) HLT_Mu15_Photon20_CaloIdL = 1;


  }
  //hlNamesTab = (hlNames_.at(2)+"\0").c_str();
  //cout<<"hlNamesTab = "<<hlNamesTab<<endl;

  for (unsigned int i=0; i!=n; ++i) {
    //cout<<"i, HLTR->wasrun(i) , HLTR->accept(i) , HLTR->error(i) = "<<i<<", "<<HLTR->wasrun(i)<<", "<<HLTR->accept(i)<<", "<<HLTR->error(i)<<endl;
    if (HLTR->wasrun(i)) hlWasRun_[i]++;
    if (HLTR->wasrun(i)) hlWasRunTab[i] = 1;
    if (!(HLTR->wasrun(i))) cout<<"hlNames(i) = "<<i<<", "<<hlNames_.at(i)<<endl;

    if (HLTR->accept(i)) hlAccept_[i]++;
    if (HLTR->accept(i)) hlAcceptTab[i] = i;
    if (HLTR->error(i) ) hlErrors_[i]++;
    if (HLTR->error(i) ) hlErrorTab[i] = 1;
    const int index(static_cast<int>(HLTR->index(i)));
    if (HLTR->accept(i)) {
      if (index>=posL1s_[i]) hltL1s_[i]++;
      if (index>=posPre_[i]) hltPre_[i]++;
    } else {
      if (index> posL1s_[i]) hltL1s_[i]++;
      if (index> posPre_[i]) hltPre_[i]++;
    }
  }

  // fin triggers ARNAUD


  debugcounter++;
  //cout<<"debug "<<debugcounter<<endl;




  //LOOP ON anti kt jets
  int FnJetsAKT = -1;
  int FnJetsAKT_pt10 = 0;
  int FnJetsAKT_pt15 = 0;
  int FnJetsAKT_pt20 = 0;

  vector <float> VemJetsAKT_pt10;
  vector <float> VemJetsAKT_pt15;
  vector <float> VemJetsAKT_pt20;

  int index_jetAKT = 0;
  if(caloantiktjetisvalid){
      FnJetsAKT = caloAntiKtJets->size();
      for(CaloJetCollection::const_iterator antiktjetiter = caloAntiKtJets->begin();antiktjetiter != caloAntiKtJets->end();antiktjetiter++){
        if(antiktjetiter->et() > 10. && fabs(antiktjetiter->eta()) < 3.) {
  	FnJetsAKT_pt10++;
  	VemJetsAKT_pt10.push_back(antiktjetiter->emEnergyFraction());
        }
        if(antiktjetiter->et() > 15. && fabs(antiktjetiter->eta()) < 3.) {
  	FnJetsAKT_pt15++;
  	VemJetsAKT_pt15.push_back(antiktjetiter->emEnergyFraction());
        }
        if(antiktjetiter->et() > 20. && fabs(antiktjetiter->eta()) < 3.) {
  	FnJetsAKT_pt20++;
  	VemJetsAKT_pt20.push_back(antiktjetiter->emEnergyFraction());
        }

        //FILL TREE
	if(antiktjetiter->et() > 10. && fabs(antiktjetiter->eta()) < 3.) {
        jetAKT_pt[index_jetAKT] = antiktjetiter->et();
        jetAKT_eta[index_jetAKT] = antiktjetiter->eta();
        jetAKT_phi[index_jetAKT] = antiktjetiter->phi();
        jetAKT_em[index_jetAKT] = antiktjetiter->emEnergyFraction();

        index_jetAKT++;
	}
      }
    }

  // LOOP ON reconstructed iterative cone jets

  // RECO -> AOD
  int FnJetsIC5 = -1;
  int FnJetsIC5_pt10 = 0;
  int FnJetsIC5_pt15 = 0;
  int FnJetsIC5_pt20 = 0;

  vector <float> VemJetsIC5_pt10;
  vector <float> VemJetsIC5_pt15;
  vector <float> VemJetsIC5_pt20;

  //  int index_jetIC5 = 0;
//   if(caloiterconejetisvalid){
//     FnJetsIC5 = caloIterConeJets->size();
//     for(CaloJetCollection::const_iterator iterconejetiter = caloIterConeJets->begin();iterconejetiter != caloIterConeJets->end();iterconejetiter++){
//       if(iterconejetiter->et() > 10. && fabs(iterconejetiter->eta()) < 3.) {
// 	FnJetsIC5_pt10++;
// 	VemJetsIC5_pt10.push_back(iterconejetiter->emEnergyFraction());
//       }
//       if(iterconejetiter->et() > 15. && fabs(iterconejetiter->eta()) < 3.) {		
// 	FnJetsIC5_pt15++;
// 	VemJetsIC5_pt15.push_back(iterconejetiter->emEnergyFraction());
//       }
//       if(iterconejetiter->et() > 20. && fabs(iterconejetiter->eta()) < 3.) {		
// 	FnJetsIC5_pt20++;
// 	VemJetsIC5_pt20.push_back(iterconejetiter->emEnergyFraction());
//       }

//        //FILL TREE
//       jetIC5_pt[index_jetIC5] = iterconejetiter->et();
//       jetIC5_eta[index_jetIC5] = iterconejetiter->eta();
//       jetIC5_phi[index_jetIC5] = iterconejetiter->phi();
//       jetIC5_em[index_jetIC5] = iterconejetiter->emEnergyFraction();   

//       index_jetIC5++;
//     }
//   }
  jetAKT_size = index_jetAKT;
    //jetAKT_size = caloAntiKtJets->size();
  //jetIC5_size = caloIterConeJets->size();

   nJetsAKT_pt15 = FnJetsAKT_pt15;
   //JetsIC5_pt15 = FnJetsIC5_pt15;
  
  //CALOMET
  if(calometisvalid){
    for(CaloMETCollection::const_iterator calometiter = caloMET->begin();calometiter != caloMET->end();calometiter++){
      calomet = calometiter->et();
    }
  }  
  //MET
  if(metisvalid){
    for(METCollection::const_iterator metiter = MET->begin();metiter != MET->end();metiter++){
      met = metiter->sumEt();    
    }
  } 
  
  //--------------------------------

  //get calo towers
  edm::Handle<CaloTowerCollection> towersH;
  iEvent.getByLabel("towerMaker", towersH);

  // get the beamspot from the Event:
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByType(theBeamSpot);
  const BeamSpot bs = *theBeamSpot;

  // get all beam spot info
  sigmaZ=theBeamSpot->sigmaZ();
  sigmaZ0Error=theBeamSpot->sigmaZ0Error();
  sq=sqrt(sigmaZ*sigmaZ+sigmaZ0Error*sigmaZ0Error);
  
  bsposx = theBeamSpot->position().x();
  bsposy = theBeamSpot->position().y();
  bsposz = theBeamSpot->position().z();

  //VINCENT
  math::XYZPoint beamspot(theBeamSpot->position().x(),theBeamSpot->position().y(),theBeamSpot->position().z());
  math::XYZPoint firstpvertex(0.,0.,0.);

  //Retrieve primary vertex collection
  Handle<reco::VertexCollection> primaryVertexColl;
  iEvent.getByLabel("offlinePrimaryVertices",primaryVertexColl);
  const reco::VertexCollection* pvcoll = primaryVertexColl.product();

  //We take only the first primary vertex, i.e. the one with the electrons
  if(pvcoll->size() > 0) {
    reco::VertexCollection::const_iterator firstpv = pvcoll->begin();
    firstpvertex.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }
  
  pvsize = pvcoll->size();
  int indexpv = 0;
  for(reco::VertexCollection::const_iterator pvIt = pvcoll->begin();pvIt != pvcoll->end(); pvIt++){
    pvx[indexpv] = pvIt->x();
    pvy[indexpv] = pvIt->y();   
    pvz[indexpv] = pvIt->z();
    
    pv_isValid[indexpv] = pvIt->isValid();
    pv_ndof[indexpv] = pvIt->ndof();
    pv_nTracks[indexpv] = pvIt->nTracks();
    pv_normChi2[indexpv] = pvIt->normalizedChi2();
    pv_totTrackSize[indexpv] = pvIt->tracksSize();

    indexpv++;
  }

  int index_mu = 0;

  for(reco::MuonCollection::const_iterator muIt = muons->begin();muIt != muons->end(); muIt++){
    if (muIt->isGlobalMuon()) index_mu++; 
  }
  
  muon_size = index_mu;

  float muonPtMax=0.;
  
  index_mu = 0;
  //LOOP OVER MUONS
  for(reco::MuonCollection::const_iterator muIt = muons->begin();muIt != muons->end(); muIt++){
    
    if (muIt->isGlobalMuon()){
      //bool isGoodMu = isGoodMuon(const reco::MuonCollection& muIt, muon::SelectionType  GlobalMuonPromptTight);
      
      if (muIt->globalTrack()->pt() > muonPtMax) muonPtMax = muIt->globalTrack()->pt();

      muon_pt[index_mu] = muIt->globalTrack()->pt();
      muon_ptError[index_mu] = muIt->globalTrack()->ptError();
      muon_eta[index_mu] = muIt->globalTrack()->eta();
      muon_etaError[index_mu] = muIt->globalTrack()->etaError();
      muon_phi[index_mu] = muIt->globalTrack()->phi();
      muon_phiError[index_mu] = muIt->globalTrack()->phiError();
      muon_theta[index_mu] = muIt->globalTrack()->theta();
      muon_thetaError[index_mu] = muIt->globalTrack()->thetaError();
      
      muon_outerPt[index_mu] = muIt->globalTrack()->outerPt();
      muon_outerEta[index_mu] = muIt->globalTrack()->outerEta();
      muon_outerPhi[index_mu] = muIt->globalTrack()->outerPhi();
      muon_outerTheta[index_mu] = muIt->globalTrack()->outerTheta();
      
      muon_px[index_mu] = muIt->globalTrack()->px();
      muon_py[index_mu] = muIt->globalTrack()->py();
      muon_pz[index_mu] = muIt->globalTrack()->pz();
      
      muon_charge[index_mu] = muIt->globalTrack()->charge();
      muon_nhitspixel[index_mu] = muIt->globalTrack()->hitPattern().numberOfValidPixelHits();
      muon_nhitstrack[index_mu] = muIt->globalTrack()->hitPattern().numberOfValidTrackerHits();
      muon_nhitsmuons[index_mu] = muIt->globalTrack()->hitPattern().numberOfValidMuonHits();
      muon_nhitstotal[index_mu] = muIt->globalTrack()->numberOfValidHits();
      muon_nlosthits[index_mu] = muIt->globalTrack()->numberOfLostHits();
      muon_nSegmentMatch[index_mu] = muIt->numberOfMatches();

      muon_isTrackerMuon[index_mu] = muIt->isTrackerMuon();

      muon_chi2[index_mu] = muIt->globalTrack()->chi2();
      muon_ndof[index_mu] = muIt->globalTrack()->ndof();
      muon_normChi2[index_mu] = muIt->globalTrack()->normalizedChi2();
      
      muon_d0[index_mu] = muIt->globalTrack()->d0();
      muon_d0Error[index_mu] = muIt->globalTrack()->d0Error();
      muon_dz_cmsCenter[index_mu] = muIt->globalTrack()->dz();
      muon_dz_beamSpot[index_mu] = muIt->globalTrack()->dz(beamspot);
      muon_dz_firstPVtx[index_mu] = muIt->globalTrack()->dz(firstpvertex);
      muon_dzError[index_mu] = muIt->globalTrack()->dzError();
      muon_dxy_cmsCenter[index_mu] = muIt->globalTrack()->dxy();
      muon_dxy_beamSpot[index_mu] = muIt->globalTrack()->dxy(beamspot);
      muon_dxy_firstPVtx[index_mu] = muIt->globalTrack()->dxy(firstpvertex);
      muon_dxyError[index_mu] = muIt->globalTrack()->dxyError();

      muon_innerPosx[index_mu] = muIt->globalTrack()->innerPosition().X();
      muon_innerPosy[index_mu] = muIt->globalTrack()->innerPosition().Y();
      muon_innerPosz[index_mu] = muIt->globalTrack()->innerPosition().Z();

      muon_trackIso03[index_mu] = muIt->isolationR03().sumPt;
      muon_trackIso05[index_mu] = muIt->isolationR05().sumPt;
      muon_trackIso03_ptInVeto[index_mu] = muIt->isolationR03().trackerVetoPt;
      muon_trackIso05_ptInVeto[index_mu] = muIt->isolationR05().trackerVetoPt;
      
      muon_emIso03[index_mu] = muIt->isolationR03().emEt;
      muon_emIso05[index_mu] = muIt->isolationR05().emEt;
      muon_emIso03_ptInVeto[index_mu] = muIt->isolationR03().emVetoEt;
      muon_emIso05_ptInVeto[index_mu] = muIt->isolationR05().emVetoEt;
      
      muon_hadIso03[index_mu] = muIt->isolationR03().hadEt;
      muon_hadIso05[index_mu] = muIt->isolationR05().hadEt;
      muon_hadIso03_ptInVeto[index_mu] = muIt->isolationR03().hadVetoEt;
      muon_hadIso05_ptInVeto[index_mu] = muIt->isolationR05().hadVetoEt;
      
      index_mu++;
    }
  }

  //--------------------------------


  //cout<<"before retrieving SC"<<endl;

  //Get a Handle on different collections

  //Get the superclusters
  edm::Handle<reco::SuperClusterCollection> pHybridSuperClusters;
  edm::Handle<reco::SuperClusterCollection> pIslandSuperClusters;
  try
    {
      iEvent.getByLabel("correctedHybridSuperClusters","",pHybridSuperClusters);
      iEvent.getByLabel("correctedMulti5x5SuperClustersWithPreshower","",pIslandSuperClusters);
    }
  catch(cms::Exception &ex)
    {
      //cout<<" could not get any of the following collections correctedHybridSuperClusters,correctedMulti5x5SuperClustersWithPreshower"<<ex;
    }


  const reco::SuperClusterCollection *hybridSuperClusters = pHybridSuperClusters.product();
  const reco::SuperClusterCollection *islandSuperClusters = pIslandSuperClusters.product();

  //cout<<"size of hybrid collection in the event "<<hybridSuperClusters->size()<<endl;
  //cout<<"size of endcap collection in the event "<<islandSuperClusters->size()<<endl;


  //Merge these two supercluster collections into one (sclusters collection)
  std::vector<const reco::SuperCluster*> sclusters;
  for (reco::SuperClusterCollection::const_iterator hsc = hybridSuperClusters->begin(); 
       hsc != hybridSuperClusters->end(); hsc++ ){sclusters.push_back(&(*hsc));}
  for (reco::SuperClusterCollection::const_iterator isc = islandSuperClusters->begin(); 
       isc != islandSuperClusters->end(); isc++ ){sclusters.push_back(&(*isc));}

  std::vector<reco::SuperClusterRef> refsclusters;
  for(unsigned int i = 0;i<hybridSuperClusters->size();i++)
    {reco::SuperClusterRef hrefsc(reco::SuperClusterRef(pHybridSuperClusters,i));refsclusters.push_back(hrefsc);}
  for(unsigned int i = 0;i<islandSuperClusters->size();i++)
    {reco::SuperClusterRef irefsc(reco::SuperClusterRef(pIslandSuperClusters,i));refsclusters.push_back(irefsc);}


  scsize = sclusters.size();

  //sort all the refSC by energy
  //cout<<"sorting the reference elements by energy"<<endl;
  std::sort(refsclusters.begin(),refsclusters.end(),refScEGreater);
  //cout<<"the elements are now : "<<endl;
  for(std::vector<reco::SuperClusterRef>::const_iterator refsclustersiter = refsclusters.begin();refsclustersiter != refsclusters.end();refsclustersiter++) {
    //cout<<"energy eta phi "<<(*refsclustersiter)->rawEnergy()+(*refsclustersiter)->preshowerEnergy()<<"  "
    //	<<  (*refsclustersiter)->eta()<<"  "
    //	<<  (*refsclustersiter)->phi()<<" "
    //	<<endl;
  }
  


  //sort all the SC by energy
  //cout<<"sorting the elements by energy"<<endl;
  std::sort(sclusters.begin(),sclusters.end(),scEGreater);
  //cout<<"the elements are now : "<<endl;
  for(std::vector<const reco::SuperCluster*>::const_iterator sclustersiter = sclusters.begin();sclustersiter != sclusters.end();sclustersiter++) {
    //cout<<"energy eta phi "<<(*sclustersiter)->rawEnergy()+(*sclustersiter)->preshowerEnergy()<<"  "
    //	<<  (*sclustersiter)->eta()<<"  "
    //	<<  (*sclustersiter)->phi()<<" "
    //	<<endl;
  }
  



  //find the sc associated to the gen electrons

  //   cout<<"before FSR"<<endl;
  //   cout<<"genelec eta  "<<unstablemomelec.eta()<<"  "
  //       <<"genelec phi  "<<unstablemomelec.phi()<<"  "
  //       <<"genelec et  "<<unstablemomelec.e()*sin(unstablemomelec.theta())<<"  "
  //       <<"genelec E  "<<unstablemomelec.e()<<endl;
  //   cout<<"genposi eta  "<<unstablemomposi.eta()<<"  "
  //       <<"genposi phi  "<<unstablemomposi.phi()<<"  "
  //       <<"genposi et  "<<unstablemomposi.e()*sin(unstablemomposi.theta())<<"  "
  //       <<"genposi E  "<<unstablemomposi.e()<<endl;
  
  //   cout<<"after FSR"<<endl;
  //   cout<<"genelec eta  "<<momelec.eta()<<"  "
  //       <<"genelec phi  "<<momelec.phi()<<"  "
  //       <<"genelec et  "<<momelec.e()*sin(momelec.theta())<<"  "
  //       <<"genelec E  "<<momelec.e()<<endl;
  //   cout<<"genposi eta  "<<momposi.eta()<<"  "
  //       <<"genposi phi  "<<momposi.phi()<<"  "
  //       <<"genposi et  "<<momposi.e()*sin(momposi.theta())<<"  "
  //       <<"genposi E  "<<momposi.e()<<endl;

  //   cout<<"Z variables"<<endl;
  //   cout<<"genboson eta  "<<momboson.eta()<<"  "
  //       <<"genboson phi  "<<momboson.phi()<<"  "
  //       <<"genboson et  "<<momboson.e()*sin(momboson.theta())<<"  "
  //       <<"genboson E  "<<momboson.e()<<endl;

  //   cout<<""<<endl;
  //   cout<<""<<endl;
  //   cout<<""<<endl;
  



  int indexelec = -50;
  int indexposi = -50;
  for(unsigned int e=0;e<sclusters.size();e++) {
    double deltaetaelec = sclusters[e]->eta()-momelec.eta();
    double deltaphielec = sclusters[e]->phi()-momelec.phi();
    
    if (deltaphielec>PI) deltaphielec = deltaphielec - 2*PI;
    else if (deltaphielec<-PI) deltaphielec = deltaphielec + 2*PI;
    
    double deltarelec = sqrt( deltaetaelec*deltaetaelec +  deltaphielec*deltaphielec );
    
    double deltaetaposi = sclusters[e]->eta()-momposi.eta();
    double deltaphiposi = sclusters[e]->phi()-momposi.phi();
    
    if (deltaphiposi>PI) deltaphiposi = deltaphiposi - 2*PI;
    else if (deltaphiposi<-PI) deltaphiposi = deltaphiposi + 2*PI;
    
    double deltarposi = sqrt( deltaetaposi*deltaetaposi +  deltaphiposi*deltaphiposi );
    
    if(deltarelec < 0.2 && indexelec == -50) indexelec = e;
    if(deltarposi < 0.2 && indexposi == -50) indexposi = e;
  }

  if(indexelec >= 0) genelec_hassc_var = 1;
  if(indexposi >= 0) genposi_hassc_var = 1;
 
  scindexforgenelec = indexelec;
  scindexforgenposi = indexposi;

  //cout<<" indexelec "<<indexelec<<endl;
  //cout<<" indexposi "<<indexposi<<endl;


  if(usegendata_)  
    {
      if(indexelec >= 0) {
	std::vector<const reco::SuperCluster*>::const_iterator scelec=sclusters.begin()+indexelec;    
	//cout<<" scelec eta "<<(*scelec)->eta()<<"  "
	//    <<" scelec phi "<<(*scelec)->phi()<<"  "
	//    <<" scelec energy "<<(*scelec)->energy()<<endl;
	scelecenergy = (*scelec)->energy();
	sceleceta = (*scelec)->eta();
	scelecphi = (*scelec)->phi();
      }
      
      if(indexposi >= 0) {
	std::vector<const reco::SuperCluster*>::const_iterator scposi=sclusters.begin()+indexposi;    
	//cout<<" scposi eta "<<(*scposi)->eta()<<"  "
	//    <<" scposi phi "<<(*scposi)->phi()<<"  "
	//    <<" scposi energy "<<(*scposi)->energy()<<endl;
	scposienergy = (*scposi)->energy();
	scposieta = (*scposi)->eta();
	scposiphi = (*scposi)->phi();
      }
    }




  //Get the MC Truth (hallelujah)
  //already retrieved in datagenerated
  if(usegendata_)  
    {
      edm::Handle<HepMCProduct> hepMC;
      iEvent.getByLabel("generator","",hepMC);
    } 


  //cout<<"to debug seg violation "<<hybridSuperClusters->size()<<endl;

  int counter = 0;
  std::vector<const reco::SuperCluster*>::const_iterator sciter=sclusters.begin();

  //const reco::SuperCluster* testsc = (*sciter);

  //SuperClusterRef screftest = testsc.castTo<SuperClusterRef>();
  //SuperClusterRef castTo<SuperClusterRef>();

  for(;sciter!=sclusters.end();sciter++)
    {
      sceta[counter] = (*sciter)->eta();
      scetacorr[counter] = etacorr( (*sciter)->eta(), pvz[0], (*sciter)->position().z() );

      sctheta[counter] = 2.*atan(exp(-1.*(*sciter)->eta()));
      scthetacorr[counter] = 2.*atan(exp(-1.*etacorr( (*sciter)->eta(), pvz[0], (*sciter)->position().z() ) ));
 
      scphi[counter] = (*sciter)->phi();
      scenergy[counter] = (*sciter)->rawEnergy()+(*sciter)->preshowerEnergy();
      scet[counter] = scenergy[counter]/cosh(sceta[counter]);

      scpx[counter] = scet[counter]*cos(scphi[counter]);
      scpy[counter] = scet[counter]*sin(scphi[counter]);
      scpz[counter] = scenergy[counter]*tanh(sceta[counter]);
      
      scx[counter] = (*sciter)->position().x();
      scy[counter] = (*sciter)->position().y();
      scz[counter] = (*sciter)->position().z();
      
      counter++;
    }

  std::vector<const reco::SuperCluster*>::const_iterator scelec=sclusters.begin();
  std::vector<const reco::SuperCluster*>::const_iterator scposi=sclusters.begin();

  if(indexelec >= 0) scelec=sclusters.begin()+indexelec;
  if(indexposi >= 0) scposi=sclusters.begin()+indexposi;

 
  //trying to see if the sc is seed associated
  

 
  edm::Handle<GsfTrackCollection> gsfTracksH ;
  iEvent.getByLabel("electronGsfTracks",gsfTracksH) ;
  const GsfTrackCollection *gsftracks = gsfTracksH.product();

  gsftracksize = gsftracks->size();

  int v=0;
  for(GsfTrackCollection::const_iterator gsftrackiter = gsftracks->begin(); 
      gsftrackiter!=gsftracks->end();
      gsftrackiter++)
    {
      //cout<<" track eta "<<gsftrackiter->eta()<<"  "
      //	  <<" track phi "<<gsftrackiter->phi()<<"  "
      //	  <<" track p "<<gsftrackiter->p()
      //	  <<" track pt "<<gsftrackiter->pt()
      //	  <<" track px "<<gsftrackiter->px()
      //	  <<" track py "<<gsftrackiter->py()
      //	  <<" track pz "<<gsftrackiter->pz()
      //	  <<endl;

      gsftracketa[v] = gsftrackiter->eta();
      gsftrackphi[v] = gsftrackiter->phi();  
      gsftrackp[v] = gsftrackiter->p();
      gsftrackpt[v] = gsftrackiter->pt();
      gsftrackpx[v] = gsftrackiter->px();
      gsftrackpy[v] = gsftrackiter->py();
      gsftrackpz[v] = gsftrackiter->pz();
	
      v++;
    }//end of loop on gsf tracks


 
  //Final GSF Electron collection
 //  edm::Handle<reco::GsfElectronCollection> pGsfElectrons;
//   bool gsfisvalid = iEvent.getByLabel("gsfElectrons","",pGsfElectrons);
//   reco::GsfElectronCollection gsfelectrons(pGsfElectrons->begin(),pGsfElectrons->end());

//   //sort all the GSF by transverse energy
//   //std::cout<<"sorting the elements by transverse energy"<<std::endl;
//   std::sort(gsfelectrons.begin(),gsfelectrons.end(),gsfEtGreater);

  //cout<<"size of gsf elec collection "<<gsfelectrons.size()<<endl;
  //cout<<"size of all sc collection "<<sclusters.size()<<endl;

  double gsfsceta = 0.;
  double gsfscphi = 0.;
  double gsfscenergy = 0.;

  //To remove spikes (ECAL CLUSTER LAZY TOOLS)
  edm::Handle<EcalRecHitCollection> EBReducedRecHits;
  iEvent.getByLabel("reducedEcalRecHitsEB",EBReducedRecHits);
  ebRecHits_ = EBReducedRecHits.product();
  edm::Handle<EcalRecHitCollection> EEReducedRecHits;
  iEvent.getByLabel("reducedEcalRecHitsEE",EEReducedRecHits);
  eeRecHits_ = EEReducedRecHits.product();
  EcalClusterLazyTools lazytool(iEvent,iSetup,InputTag("reducedEcalRecHitsEB"),InputTag("reducedEcalRecHitsEE"));
  //EcalSeverityLevelAlgo ecalalgo;

  float gsfPtMax = 0.;

  gsf_size = gsfelectrons.size();
  int e=0;
  reco::GsfElectronCollection::const_iterator gsfiter = gsfelectrons.begin();
  for(;gsfiter!=gsfelectrons.end();gsfiter++)
    {
      gsfsceta = gsfiter->superCluster()->eta();
      gsfscphi = gsfiter->superCluster()->phi();
      gsfscenergy = gsfiter->superCluster()->energy();

      //cout<<" gsf sc eta "<<gsfiter->superCluster()->eta()<<"  "
      //  <<" gsf sc phi "<<gsfiter->superCluster()->phi()<<"  "
      //  <<" gsf sc energy "<<gsfiter->superCluster()->energy()<<endl;

      scindexforgsf[e] = -3;
      //try to get the index for the sc assoc to this gsf
      reco::SuperClusterRef gsfrefsc = gsfiter->superCluster();
      for(unsigned int k=0;k<refsclusters.size();k++){
	if(gsfrefsc == refsclusters[k]){
	  //cout<<"matching ref gsf is done "<<refsclusters[k]->energy()<<"   "<<gsfscenergy<<endl;
	  scindexforgsf[e] = k;
	}
      }
      
      
      //Try to see if the gsf is assoc to the gen elec assoc sc
      if(indexelec >= 0) {
	double elecdeltaeta = sceleceta-gsfsceta;
	double elecdeltaphi = scelecphi-gsfscphi;
	if (elecdeltaphi>PI) elecdeltaphi = elecdeltaphi - 2*PI;
	else if (elecdeltaphi<-PI) elecdeltaphi = elecdeltaphi + 2*PI;
	double elecdeltar = sqrt( elecdeltaeta*elecdeltaeta +  elecdeltaphi*elecdeltaphi );
      
	if(elecdeltar < 0.2) 
	  {
	    //cout<<"matching gsf is done "<<scelecenergy<<"   "<<gsfscenergy<<endl;
	    scgsfmatched[indexelec] = -1.;
	    scelecgsfmatched = 1.;
	  }
      }


      //Try to see if the gsf is assoc to the gen posi assoc sc
      if(indexposi >= 0) {
	double posideltaeta = scposieta-gsfsceta;
	double posideltaphi = scposiphi-gsfscphi;
	if (posideltaphi>PI) posideltaphi = posideltaphi - 2*PI;
	else if (posideltaphi<-PI) posideltaphi = posideltaphi + 2*PI;
	double posideltar = sqrt( posideltaeta*posideltaeta +  posideltaphi*posideltaphi );
      
	if(posideltar < 0.2) 
	  {
	    //cout<<"matching gsf is done "<<scposienergy<<"   "<<gsfscenergy<<endl;
	    scgsfmatched[indexposi] = 1.;
	    scposigsfmatched = 1.;
	  }
      }

      debugcounter++;
      //cout<<"debug "<<debugcounter<<endl;

      double deltaetagenelec = gsfsceta-momelec.eta();
      double deltaphigenelec = gsfscphi-momelec.phi();
      if (deltaphigenelec>PI) deltaphigenelec = deltaphigenelec - 2*PI;
      else if (deltaphigenelec<-PI) deltaphigenelec = deltaphigenelec + 2*PI;
      double deltargenelec = sqrt(deltaetagenelec*deltaetagenelec + deltaphigenelec*deltaphigenelec);


      debugcounter++;
      //cout<<"debug "<<debugcounter<<endl;

      double deltaetagenposi = gsfsceta-momposi.eta();
      double deltaphigenposi = gsfscphi-momposi.phi();
      if (deltaphigenposi>PI) deltaphigenposi = deltaphigenposi - 2*PI;
      else if (deltaphigenposi<-PI) deltaphigenposi = deltaphigenposi + 2*PI;
      double deltargenposi = sqrt(deltaetagenposi*deltaetagenposi + deltaphigenposi*deltaphigenposi);


      debugcounter++;
      //cout<<"debug "<<debugcounter<<endl;

      if(deltargenelec < 0.2 && gsfindexforgenelec == -3) {gsfindexforgenelec = e;}
      if(deltargenposi < 0.2 && gsfindexforgenposi == -3) {gsfindexforgenposi = e;}


      //debugcounter++;
      //cout<<"debug "<<debugcounter<<endl;

      //Fill the gsf related variables
      gsf_e[e] = gsfiter->energy();
      gsf_p[e] = gsfiter->p();
      gsf_pt[e] = gsfiter->pt();
      gsf_class[e] = gsfiter->classification();
      gsf_e2x5overe5x5[e] = gsfiter->scE2x5Max()/gsfiter->scE5x5();
      gsf_e1x5overe5x5[e] = gsfiter->scE1x5()/gsfiter->scE5x5();
      gsf_eta[e] = gsfiter->eta();
      gsf_phi[e] = gsfiter->phi();
      gsf_px[e] = gsfiter->px();
      gsf_py[e] = gsfiter->py();
      gsf_pz[e] = gsfiter->pz();


      //debugcounter++;
      //cout<<"debug "<<debugcounter<<endl;

      //if(gsfiter->phi()==0.) std::cout<<"here funny event"<<std::endl;

      gsf_deltaeta[e] = gsfiter->deltaEtaSuperClusterTrackAtVtx();
      gsf_deltaphi[e] = gsfiter->deltaPhiSuperClusterTrackAtVtx();
      gsf_hovere[e] = gsfiter->hadronicOverEm();

      gsf_trackiso[e] = gsfiter->dr03TkSumPt();
      gsf_ecaliso[e] = gsfiter->dr03EcalRecHitSumEt();
      gsf_hcaliso1[e] = gsfiter->dr03HcalDepth1TowerSumEt();
      gsf_hcaliso2[e] = gsfiter->dr03HcalDepth2TowerSumEt();
  
      gsf_charge[e] = gsfiter->charge();
      gsf_sigmaetaeta[e] = gsfiter->sigmaEtaEta();
      gsf_sigmaIetaIeta[e] = gsfiter->sigmaIetaIeta();
      //gsf_isecaldriven[e] = gsfiter->isEcalDriven();
      //gsf_istrackerdriven[e] = gsfiter->isTrackerDriven();
      gsf_isecaldriven[e] = gsfiter->ecalDrivenSeed();
      gsf_istrackerdriven[e] = gsfiter->trackerDrivenSeed();


      debugcounter++;
      //cout<<"debug "<<debugcounter<<endl;
      //cout<<"ici   "<<endl;

      gsfsc_e[e] = gsfiter->superCluster()->rawEnergy()+gsfiter->superCluster()->preshowerEnergy();
      gsfsc_pt[e] = (gsfiter->superCluster()->rawEnergy()+gsfiter->superCluster()->preshowerEnergy())/cosh(gsfiter->superCluster()->eta());
      gsfsc_eta[e] = gsfiter->superCluster()->eta();
      gsfsc_phi[e] = gsfiter->superCluster()->phi();
      gsfsc_px[e] = gsfsc_pt[e]*cos(gsfsc_phi[e]);
      gsfsc_py[e] = gsfsc_pt[e]*sin(gsfsc_phi[e]);
      gsfsc_pz[e] = (gsfiter->superCluster()->rawEnergy()+gsfiter->superCluster()->preshowerEnergy())*tanh(gsfiter->superCluster()->eta());

      gsf_gsfet[e] = gsfiter->caloEnergy()*sin(gsfiter->p4().theta());

      //HOMEMADE SKIMMING!!!
      if (gsfiter->caloEnergy()*sin(gsfiter->p4().theta()) > gsfPtMax) gsfPtMax = gsfiter->caloEnergy()*sin(gsfiter->p4().theta());

      gsf_theta[e] = gsfiter->theta();
      gsf_isEB[e] = gsfiter->isEB();
      gsf_isEE[e] = gsfiter->isEE();
      gsf_deltaEtaATcalo[e] = gsfiter->deltaEtaSeedClusterTrackAtCalo();
      gsf_deltaPhiATcalo[e] = gsfiter->deltaPhiSeedClusterTrackAtCalo();
      gsf_ecalEnergy[e] = gsfiter->ecalEnergy();
      gsf_eOVERp[e] = gsfiter->eSuperClusterOverP();
      gsf_dxy[e] = gsfiter->gsfTrack()->dxy();
      gsf_vz[e] = gsfiter->gsfTrack()->vz();
      gsf_nHits[e] = gsfiter->gsfTrack()->numberOfValidHits();   
      gsf_nLostInnerHits[e] = gsfiter->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();   
      gsf_nLostOuterHits[e] = gsfiter->gsfTrack()->trackerExpectedHitsOuter().numberOfLostHits();   
      gsf_fBrem[e] = gsfiter->fbrem();
      gsf_e1x5[e] =gsfiter->e1x5() ;
      gsf_e2x5[e] =gsfiter->e2x5Max() ;
      gsf_e5x5[e] =gsfiter->e5x5() ;

      // pbAOD
// Hence, I recommend changing:
//         const reco::CaloClusterPtr seed = elec_iter->superCluster()->seed();
// // seed cluster
//         const DetId seedId = seed->seed();
//         EcalSeverityLevelAlgo severity;
//         double myswissCross = severity.swissCross(seedId, *myRecHits) ;

// to
//          DetId idEB = EcalClusterTools::getMaximum(
// ele->superCluster()->hitsAndFractions(), &(*barrelRecHits) ).first;
//          EcalRecHitCollection::const_iterator thisHitEB = barrelRecHits->find(idEB);
//          double swissCrossNoI85 = EcalSeverityLevelAlgo::swissCross(
// idEB, (*barrelRecHits), 5., true);
//end


      const reco::CaloClusterPtr seed = gsfiter->superCluster()->seed();
      //      DetId id = lazytool.getMaximum(*seed).first;
//       float emax = lazytool.getMaximum(*seed).second;
/////      DetId idEB = EcalClusterTools::getMaximum(gsfiter->superCluster()->hitsAndFractions(), &(*barrelRecHits) ).first;
//       EcalRecHitCollection::const_iterator thisHitEB = barrelRecHits->find(idEB);
//       double swissCrossNoI85 = EcalSeverityLevelAlgo::swissCross(idEB, (*barrelRecHits), 5., true);

//       gsf_eMax[e] = emax;
//       gsf_SwissCross[e] = ecalalgo.swissCross(id, *(getEcalRecHitCollection(*seed)),0.);

      gsf_e1x3[e] = lazytool.e1x3(*seed);
//       gsf_e3x1[e] = lazytool.e3x1(*seed);
//       gsf_e2x2[e] = lazytool.e2x2(*seed);
//       gsf_e3x2[e] = lazytool.e3x2(*seed);
//       gsf_e3x3[e] = lazytool.e3x3(*seed);
//       gsf_e4x4[e] = lazytool.e4x4(*seed);
//       gsf_e2x5Right[e] = lazytool.e2x5Right(*seed);
//       gsf_e2x5Left[e] = lazytool.e2x5Left(*seed);
//       gsf_e2x5Top[e] = lazytool.e2x5Top(*seed);
//       gsf_e2x5Bottom[e] = lazytool.e2x5Bottom(*seed);
//       gsf_e2x5Max[e] = lazytool.e2x5Max(*seed);
//       gsf_eLeft[e] = lazytool.eLeft(*seed);
//       gsf_eRight[e] = lazytool.eRight(*seed);
//       gsf_eTop[e] = lazytool.eTop(*seed);
//       gsf_eBottom[e] = lazytool.eBottom(*seed);
//       gsf_e2nd[e] = lazytool.e2nd(*seed);


      //cout<<"ici 2  "<<endl;

      //debugcounter++;
      //cout<<"debug "<<debugcounter<<endl;

      bool gsfetbarrel = gsf_gsfet[e] > 25.;
      bool gsfetendcap = gsf_gsfet[e] > 25.;

      bool barrelsc = fabs(gsfsc_eta[e]) < 1.442;
      bool endcapsc = (fabs(gsfsc_eta[e]) > 1.56) && (fabs(gsfsc_eta[e]) < 2.5);

      bool deltaetabarrel = fabs(gsf_deltaeta[e]) < 0.005;
      bool deltaetaendcap = fabs(gsf_deltaeta[e]) < 0.007;

      bool deltaphibarrel = fabs(gsf_deltaphi[e]) < 0.09;
      bool deltaphiendcap = fabs(gsf_deltaphi[e]) < 0.09;

      bool hoverebarrel  = gsf_hovere[e] < 0.05;
      bool hovereendcap  = gsf_hovere[e] < 0.05;

      bool sigmaIetaIetabarrel  = true;
      bool sigmaIetaIetaendcap  = gsf_sigmaIetaIeta[e] < 0.03;

      bool e2x5overe5x5barrel  = (gsf_e2x5overe5x5[e] > 0.94) || (gsf_e1x5overe5x5[e] > 0.83);
      bool e2x5overe5x5endcap  = true;

      bool ecalisobarrel = (gsf_ecaliso[e]+gsf_hcaliso1[e]) < (2.+0.03*gsf_gsfet[e]);
      bool ecalisoendcap = true;
      if(gsf_gsfet[e] < 50.) {
	ecalisoendcap = (gsf_ecaliso[e]+gsf_hcaliso1[e]) < 2.5;
      }
      else {
	ecalisoendcap = (gsf_ecaliso[e]+gsf_hcaliso1[e]) < (2.5+0.03*(gsf_gsfet[e]-50.));
      }

      bool hcaliso2barrel  = true;
      bool hcaliso2endcap  = gsf_hcaliso2[e] < 0.5;

      bool trackisobarrel  = gsf_trackiso[e] < 7.5;
      bool trackisoendcap  = gsf_trackiso[e] < 15.;

      bool noMissingHits = (gsf_nLostInnerHits[e] + gsf_nLostOuterHits[e]) == 0;

      //cout<<"ic3  "<<endl;
      //debugcounter++;
      //cout<<"debug "<<debugcounter<<endl;

      //Boolean HEEP cuts
      gsfpass_ET[e] = (gsfetbarrel && barrelsc) || (gsfetendcap && endcapsc); 
      gsfpass_PT[e] = true; 
      gsfpass_DETETA[e] = true; 
      gsfpass_CRACK[e] = true; 
      gsfpass_DETAIN[e] = (deltaetabarrel && barrelsc) || (deltaetaendcap && endcapsc); 
      gsfpass_DPHIIN[e] = (deltaphibarrel && barrelsc) || (deltaphiendcap && endcapsc); 
      gsfpass_HADEM[e] = (hoverebarrel && barrelsc) || (hovereendcap && endcapsc); 
      gsfpass_SIGMAIETAIETA[e] = (sigmaIetaIetabarrel && barrelsc) || (sigmaIetaIetaendcap && endcapsc);
      gsfpass_E2X5OVER5X5[e] = (e2x5overe5x5barrel && barrelsc) || (e2x5overe5x5endcap && endcapsc); 
      gsfpass_ISOLEMHADDEPTH1[e] = (ecalisobarrel && barrelsc) || (ecalisoendcap && endcapsc);
      gsfpass_ISOLHADDEPTH2[e] = (hcaliso2barrel && barrelsc) || (hcaliso2endcap && endcapsc); 
      gsfpass_ISOLPTTRKS[e] = (trackisobarrel && barrelsc) || (trackisoendcap && endcapsc); 
      gsfpass_ECALDRIVEN[e] = gsf_isecaldriven[e]; 
      gsfpass_INVALID[e] = true;
      gsfpass_NOMISSINGHITS[e] = noMissingHits;
      
      gsfpass_ID[e] = (gsfpass_DETAIN[e] && gsfpass_DPHIIN[e] && gsfpass_HADEM[e] && gsfpass_SIGMAIETAIETA[e] && gsfpass_E2X5OVER5X5[e]);
      gsfpass_ISO[e] = (gsfpass_ISOLEMHADDEPTH1[e] && gsfpass_ISOLHADDEPTH2[e] && gsfpass_ISOLPTTRKS[e]);


      //debugcounter++;
      //cout<<"debug "<<debugcounter<<endl;

      gsfpass_HEEP[e] = gsfpass_ET[e] && gsfpass_DETAIN[e] && gsfpass_DPHIIN[e] && gsfpass_HADEM[e] && gsfpass_SIGMAIETAIETA[e] && gsfpass_E2X5OVER5X5[e] && gsfpass_ISOLEMHADDEPTH1[e] && gsfpass_ISOLHADDEPTH2[e] && gsfpass_ISOLPTTRKS[e] && gsfpass_ECALDRIVEN[e];
      

      //debugcounter++;
      //cout<<"debug "<<debugcounter<<endl;

      //charge info
      scpixcharge[e] = gsfiter->scPixCharge();
      if(gsfiter->closestCtfTrackRef().isNonnull()) ctfcharge[e] = gsfiter->closestCtfTrackRef()->charge();
      gsfcharge[e] = gsfiter->gsfTrack()->charge();
      gsfctfscpixconsistent[e] = gsfiter->isGsfCtfScPixChargeConsistent();
      gsfscpixconsistent[e] = gsfiter->isGsfScPixChargeConsistent();
      gsfctfconsistent[e] = gsfiter->isGsfCtfChargeConsistent();

      //cout<<"ici 4  "<<endl;

      //debugcounter++;
      //cout<<"debug "<<debugcounter<<endl;

      //increment index for gsf
      e++;


      //debugcounter++;
      //cout<<"debug "<<debugcounter<<endl;

    }

  debugcounter++;
  //cout<<"outdebug "<<debugcounter<<endl;


  //Have all final info to make effi plots

  //Is the gen elec/posi a SC ??
  if(scindexforgenelec >= 0) genelechassc = true;
  if(scindexforgenposi >= 0) genposihassc = true;

  //cout<<"ici 5  "<<endl;

  //debugcounter++;
  //cout<<"outdebug "<<debugcounter<<endl;

  
  //HOMEMADE SKIMMING!!!
  if (gsfPtMax>25.) mytree->Fill();

  //if (gsfPtMax>25. && muonPtMax>25.) mytree->Fill(); // SKIM 1ele1muon
  //mytree->Fill();

  //cout<<"ici 6  "<<endl;
    }
   
}//end of analyze method



// FROM ARNAUD
void GsfCheckerTree::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  using namespace std;
  using namespace edm;
  
  cout<<"dans beginrun run number"<<iRun.id()<<" and n = "<<hlNames_.size()<<endl;

  bool changed (true);
  if (hltConfig_.init(iRun,iSetup,hlTriggerResults_.process(),changed)) {
    if (changed) {
      // dump previous
      //dumpReport();
      nEvents_=0;
      nWasRun_=0;
      nAccept_=0;
      nErrors_=0;
      // const edm::TriggerNames & triggerNames = iEvent.triggerNames(*HLTR);
      hlNames_=hltConfig_.triggerNames();
      const unsigned int n(hlNames_.size());
      hlWasRun_.resize(n);
      hltL1s_.resize(n);
      hltPre_.resize(n);
      hlAccept_.resize(n);
      hlErrors_.resize(n);
      posL1s_.resize(n);
      posPre_.resize(n);
      for (unsigned int i=0; i<200; i++) {
 	hlWasRunTab[i]=0;
	hlAcceptTab[i]=0;
	hlErrorTab[i]=0;
      }
     for (unsigned int i=0; i!=n; ++i) {
       cout<<"hlNames("<<i<<") = "<<hlNames_.at(i)<<endl;
       hlWasRun_[i]=0;
	hltL1s_[i]=0;
	hltPre_[i]=0;
	hlAccept_[i]=0;
	hlErrors_[i]=0;
	posL1s_[i]=-1;
	posPre_[i]=-1;
	const std::vector<std::string>& moduleLabels(hltConfig_.moduleLabels(i));
	for (unsigned int j=0; j<moduleLabels.size(); ++j) {
	  if (hltConfig_.moduleType(moduleLabels[j])=="HLTLevel1GTSeed") {
	    posL1s_[i]=j;
	  }
	  if (hltConfig_.moduleType(moduleLabels[j])=="HLTPrescaler"   ) {
	    posPre_[i]=j;
	  }
	}
     }
    }
  } else {
    // dump previous
    //dumpReport();
    // clear
    nEvents_=0;
    nWasRun_=0;
    nAccept_=0;
    nErrors_=0;
    hlWasRun_.clear();
    hltL1s_.clear();
    hltPre_.clear();
    hlAccept_.clear();
    hlErrors_.clear();
    posL1s_.clear();
    posPre_.clear();
    hlNames_.clear();
  }
  //return true;
}
      


// ------------ method called once each job just before starting event loop  ------------
void 
GsfCheckerTree::beginJob()
{
  edm::Service<TFileService> fs;
  mytree = fs->make<TTree>("tree","tree");

//   Analysis = new TFile("GsfCheckerTree.root","RECREATE");
//   mytree = new TTree("tree","tree");

  //Added from Vincent
  //------------------------------
  //TRIGGER

  
  mytree->Branch("hltCount",&hltCount,"hltCount/I");
  mytree->Branch("L1trigger_size", &L1trigger_size, "L1trigger_size/I"); 
  mytree->Branch("L1trigger_bool", &L1trigger_bool, "L1trigger_bool[L1trigger_size]/I");
  mytree->Branch("PhysDecl_bool", &PhysDecl_bool, "PhysDecl_bool/I");
  mytree->Branch("HLTriggers", HLTriggers, "HLTriggers[hltCount]/I");

  
  //FROM ARNAUD
  mytree->Branch("nWasRun_",&nWasRun_,"nWasRun_/I");
  mytree->Branch("nAccept_",&nAccept_,"nAccept_/I");
  mytree->Branch("nErrors_",&nErrors_,"nErrors_/I");
  //mytree->Branch("hlWasRun_",std::vector<unsigned int>,"&phlWasRun_");
  mytree->Branch("hlWasRun_",&hlWasRun_,"hlWasRun_/I");
  mytree->Branch("hlWasRunTab",hlWasRunTab,"hlWasRunTab[200]/I");
  mytree->Branch("hlAccept_",&hlAccept_);
  mytree->Branch("hlAcceptTab",hlAcceptTab,"hlAcceptTab[300]/I");
  mytree->Branch("hlErrorTab",hlErrorTab,"hlErrorTab[200]/I");
  //mytree->Branch("hlNamesTab",hlNamesTab,"hlNamesTab[200]/C");
  mytree->Branch("hlNamesTab",&hlNamesTab,"hlNamesTab/C");
  mytree->Branch("hlNames_",&hlNames_);

  mytree->Branch("HLT_Ele10_SW_EleId_L1R",&HLT_Ele10_SW_EleId_L1R,"HLT_Ele10_SW_EleId_L1R/I");
  mytree->Branch("HLT_Ele10_SW_L1R",&HLT_Ele10_SW_L1R,"HLT_Ele10_SW_L1R/I");
  mytree->Branch("HLT_Ele15_LW_L1R",&HLT_Ele15_LW_L1R,"HLT_Ele15_LW_L1R/I");
  mytree->Branch("HLT_Ele15_SW_L1R",&HLT_Ele15_SW_L1R,"HLT_Ele15_SW_L1R/I");
  mytree->Branch("HLT_Ele15_SW_EleId_L1R",&HLT_Ele15_SW_EleId_L1R,"HLT_Ele15_SW_EleId_L1R/I");
  mytree->Branch("HLT_Ele15_SW_CaloEleId_L1R",&HLT_Ele15_SW_CaloEleId_L1R,"HLT_Ele15_SW_CaloEleId_L1R/I");
  mytree->Branch("HLT_Ele15_SiStrip_L1R",&HLT_Ele15_SiStrip_L1R,"HLT_Ele15_SiStrip_L1R/I");
  mytree->Branch("HLT_Ele20_SW_L1R",&HLT_Ele20_SW_L1R,"HLT_Ele20_SW_L1R/I");
  mytree->Branch("HLT_Ele20_SiStrip_L1R",&HLT_Ele20_SiStrip_L1R,"HLT_Ele20_SiStrip_L1R/I");
  mytree->Branch("HLT_Ele25_SW_L1R",&HLT_Ele25_SW_L1R,"HLT_Ele25_SW_L1R/I");
  mytree->Branch("HLT_DoubleEle4_SW_eeRes_L1R",&HLT_DoubleEle4_SW_eeRes_L1R,"HLT_DoubleEle4_SW_eeRes_L1R/I");
  mytree->Branch("HLT_DoubleEle10_SW_L1R",&HLT_DoubleEle10_SW_L1R,"HLT_DoubleEle10_SW_L1R/I");
  mytree->Branch("HLT_Mu15_Photon20_CaloIdL",&HLT_Mu15_Photon20_CaloIdL,"HLT_Mu15_Photon20_CaloIdL/I"); //FROM VINCENT
  //END FROM ARNAUD

  mytree->Branch("HLT_Photon20_CaloIdVL_IsoL_v1", &HLT_Photon20_CaloIdVL_IsoL_v1, "HLT_Photon20_CaloIdVL_IsoL_v1/I");
  mytree->Branch("HLT_DoublePhoton33_vx", & HLT_DoublePhoton33_vx, " HLT_DoublePhoton33_vx/I");
  mytree->Branch("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2", &HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2, "HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2/I");
  mytree->Branch("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2", &HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2, "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2/I");
  mytree->Branch("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2", &HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2, "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2/I");

  //ok
  //GLOBAL 
  mytree->Branch("nJetsAKT_pt15", &nJetsAKT_pt15, "nJetsAKT_pt15/I");
  //mytree->Branch("nJetsIC5_pt15", &nJetsIC5_pt15, "nJetsIC5_pt15/I");
  mytree->Branch("calomet", &calomet, "calomet/F");
  mytree->Branch("met", &met, "met/F");

  // ok
  //JETS

  //  //AKT
    mytree->Branch("jetAKT_size", &jetAKT_size, "jetAKT_size/I");
    mytree->Branch("jetAKT_pt", jetAKT_pt, "jetAKT_pt[jetAKT_size]/F");
    mytree->Branch("jetAKT_eta", jetAKT_eta, "jetAKT_eta[jetAKT_size]/F");
    mytree->Branch("jetAKT_phi", jetAKT_phi, "jetAKT_phi[jetAKT_size]/F");
    mytree->Branch("jetAKT_em", jetAKT_em, "jetAKT_em[jetAKT_size]/F");

  //IC5
 //  mytree->Branch("jetIC5_size", &jetIC5_size, "jetIC5_size/I");
//   mytree->Branch("jetIC5_pt", jetIC5_pt, "jetIC5_pt[jetIC5_size]/F");
//   mytree->Branch("jetIC5_eta", jetIC5_eta, "jetIC5_eta[jetIC5_size]/F");
//   mytree->Branch("jetIC5_phi", jetIC5_phi, "jetIC5_phi[jetIC5_size]/F");
//   mytree->Branch("jetIC5_em", jetIC5_em, "jetIC5_em[jetIC5_size]/F");

  
  
  //MUONS
  mytree->Branch("muon_size", &muon_size, "muon_size/I");
  mytree->Branch("muon_pt", muon_pt, "muon_pt[muon_size]/F");
  mytree->Branch("muon_ptError", muon_ptError, "muon_ptError[muon_size]/F");
  mytree->Branch("muon_eta", muon_eta, "muon_eta[muon_size]/F");
  mytree->Branch("muon_etaError", muon_etaError, "muon_etaError[muon_size]/F");
  mytree->Branch("muon_theta", muon_theta, "muon_theta[muon_size]/F");
  mytree->Branch("muon_thetaError", muon_thetaError, "muon_thetaError[muon_size]/F"); 
  mytree->Branch("muon_phi", muon_phi, "muon_phi[muon_size]/F");
  mytree->Branch("muon_phiError", muon_phiError, "muon_phiError[muon_size]/F");
  mytree->Branch("muon_outerPt", muon_outerPt, "muon_outerPt[muon_size]/F");
  mytree->Branch("muon_outerEta", muon_outerEta, "muon_outerEta[muon_size]/F");
  mytree->Branch("muon_outerPhi", muon_outerPhi, "muon_outerPhi[muon_size]/F");
  mytree->Branch("muon_outerTheta", muon_outerTheta, "muon_outerTheta[muon_size]/F"); 
  mytree->Branch("muon_px", muon_px, "muon_px[muon_size]/F");
  mytree->Branch("muon_py", muon_py, "muon_py[muon_size]/F");
  mytree->Branch("muon_pz", muon_pz, "muon_pz[muon_size]/F");
  mytree->Branch("muon_charge", muon_charge, "muon_charge[muon_size]/I");
  mytree->Branch("muon_nhitspixel", muon_nhitspixel, "muon_nhitspixel[muon_size]/I");
  mytree->Branch("muon_nhitstrack", muon_nhitstrack, "muon_nhitstrack[muon_size]/I");
  mytree->Branch("muon_nhitsmuons", muon_nhitsmuons, "muon_nhitsmuons[muon_size]/I");
  mytree->Branch("muon_nhitstotal", muon_nhitstotal, "muon_nhitstotal[muon_size]/I");
  mytree->Branch("muon_nlosthits", muon_nlosthits, "muon_nlosthits[muon_size]/I");
  mytree->Branch("muon_nSegmentMatch", muon_nSegmentMatch, "muon_nSegmentMatch[muon_size]/I");
  mytree->Branch("muon_isTrackerMuon", muon_isTrackerMuon, "muon_isTrackerMuon[muon_size]/B");
  mytree->Branch("muon_chi2", muon_chi2, "muon_chi2[muon_size]/F");
  mytree->Branch("muon_ndof", muon_ndof, "muon_ndof[muon_size]/I");
  mytree->Branch("muon_normChi2", muon_normChi2, "muon_normChi2[muon_size]/F");
  mytree->Branch("muon_d0", muon_d0, "muon_d0[muon_size]/F");
  mytree->Branch("muon_d0Error", muon_d0Error, "muon_d0Error[muon_size]/F");
  mytree->Branch("muon_dz_cmsCenter", muon_dz_cmsCenter, "muon_dz_cmsCenter[muon_size]/F");
  mytree->Branch("muon_dz_beamSpot", muon_dz_beamSpot, "muon_dz_beamSpot[muon_size]/F");
  mytree->Branch("muon_dz_firstPVtx", muon_dz_firstPVtx, "muon_dz_firstPVtx[muon_size]/F");
  mytree->Branch("muon_dzError", muon_dzError, "muon_dzError[muon_size]/F");
  mytree->Branch("muon_dxy_cmsCenter", muon_dxy_cmsCenter, "muon_dxy_cmsCenter[muon_size]/F");
  mytree->Branch("muon_dxy_beamSpot", muon_dxy_beamSpot, "muon_dxy_beamSpot[muon_size]/F");
  mytree->Branch("muon_dxy_firstPVtx", muon_dxy_firstPVtx, "muon_dxy_firstPVtx[muon_size]/F");
  mytree->Branch("muon_dxyError", muon_dxyError, "muon_dxyError[muon_size]/F");
  mytree->Branch("muon_innerPosx", muon_innerPosx, "muon_innerPosx[muon_size]/F");
  mytree->Branch("muon_innerPosy", muon_innerPosy, "muon_innerPosy[muon_size]/F");
  mytree->Branch("muon_innerPosz", muon_innerPosz, "muon_innerPosz[muon_size]/F");
  mytree->Branch("muon_trackIso03", muon_trackIso03, "muon_trackIso03[muon_size]/F");
  mytree->Branch("muon_trackIso05", muon_trackIso05, "muon_trackIso05[muon_size]/F");
  mytree->Branch("muon_trackIso03_ptInVeto", muon_trackIso03_ptInVeto, "muon_trackIso03_ptInVeto[muon_size]/F");
  mytree->Branch("muon_trackIso05_ptInVeto", muon_trackIso05_ptInVeto, "muon_trackIso05_ptInVeto[muon_size]/F");
  mytree->Branch("muon_emIso03", muon_emIso03, "muon_emIso03[muon_size]/F");
  mytree->Branch("muon_emIso05", muon_emIso05, "muon_emIso05[muon_size]/F");
  mytree->Branch("muon_emIso03_ptInVeto", muon_emIso03_ptInVeto, "muon_emIso03_ptInVeto[muon_size]/F");
  mytree->Branch("muon_emIso05_ptInVeto", muon_emIso05_ptInVeto, "muon_emIso05_ptInVeto[muon_size]/F");
  mytree->Branch("muon_hadIso03", muon_hadIso03, "muon_hadIso03[muon_size]/F");
  mytree->Branch("muon_hadIso05", muon_hadIso05, "muon_hadIso05[muon_size]/F");
  mytree->Branch("muon_hadIso03_ptInVeto", muon_hadIso03_ptInVeto, "muon_hadIso03_ptInVeto[muon_size]/F");
  mytree->Branch("muon_hadIso05_ptInVeto", muon_hadIso05_ptInVeto, "muon_hadIso05_ptInVeto[muon_size]/F");


  //------------------------------
  //pas ok


  //Run and event number
  mytree->Branch("runnumber",&runnumber,"runnumber/I");
  mytree->Branch("eventnumber",&eventnumber,"eventnumber/I");
  mytree->Branch("luminosityBlock",&luminosityBlock,"luminosityBlock/I"); //add Laurent

  mytree->Branch("eventcounter",&eventcounter,"eventcounter/I");

  mytree->Branch("processid",&processid,"processid/I");
  mytree->Branch("pthat",&pthat,"pthat/F");
  mytree->Branch("alphaqcd",&alphaqcd,"alphaqcd/F");
  mytree->Branch("alphaqed",&alphaqed,"alphaqed/F");
  mytree->Branch("qscale",&qscale,"qscale/F");
  mytree->Branch("weight",&weight,"weight/F");

  //Beam spot variables
  mytree->Branch("sigmaZ",&sigmaZ,"sigmaZ/F");
  mytree->Branch("sigmaZ0Error",&sigmaZ0Error,"sigmaZ0Error/F");
  mytree->Branch("sq",&sq,"sq/F");
  mytree->Branch("bsposx",&bsposx,"bsposx/F");
  mytree->Branch("bsposy",&bsposy,"bsposy/F");
  mytree->Branch("bsposz",&bsposz,"bsposz/F");

  //Primary vertex variables
  mytree->Branch("pvsize",&pvsize,"pvsize/I");
  mytree->Branch("pvx",&pvx,"pvx[pvsize]/F");
  mytree->Branch("pvy",&pvy,"pvy[pvsize]/F");
  mytree->Branch("pvz",&pvz,"pvz[pvsize]/F");

  mytree->Branch("pv_isValid",&pv_isValid,"pv_isValid[pvsize]/B");
  mytree->Branch("pv_ndof",&pv_ndof,"pv_ndof[pvsize]/F");
  mytree->Branch("pv_nTracks",&pv_nTracks,"pv_nTracks[pvsize]/I");
  mytree->Branch("pv_normChi2",&pv_normChi2,"pv_normChi2[pvsize]/F");
  mytree->Branch("pv_totTrackSize",&pv_totTrackSize,"pv_totTrackSize[pvsize]/I");

  //Supercluster variables
  //fill e,et,eta,phi,charge for every SC in the event
  mytree->Branch("scsize",&scsize,"scsize/I");

  mytree->Branch("scenergy",scenergy,"scenergy[scsize]/F");
  mytree->Branch("sceta",sceta,"sceta[scsize]/F");
  mytree->Branch("scetacorr",scetacorr,"scetacorr[scsize]/F");
  mytree->Branch("sctheta",sctheta,"sctheta[scsize]/F");
  mytree->Branch("scthetacorr",scthetacorr,"scthetacorr[scsize]/F");
  mytree->Branch("scet",scet,"scet[scsize]/F");
  mytree->Branch("scphi",scphi,"scphi[scsize]/F");
  mytree->Branch("scpx",scpx,"scpx[scsize]/F");
  mytree->Branch("scpy",scpy,"scpy[scsize]/F");
  mytree->Branch("scpz",scpz,"scpz[scsize]/F");
  mytree->Branch("scx",scx,"scx[scsize]/F");
  mytree->Branch("scy",scy,"scy[scsize]/F");
  mytree->Branch("scz",scz,"scz[scsize]/F");

  mytree->Branch("scgsfmatched",scgsfmatched,"scgsfmatched[scsize]/F");


  

  //generated variables for the tree (after FSR)
  mytree->Branch("genelec_e_branch",&genelec_e_var,"genelec_e_branch/F");
  mytree->Branch("genelec_eta_branch",&genelec_eta_var,"genelec_eta_branch/F");
  mytree->Branch("genelec_phi_branch",&genelec_phi_var,"genelec_phi_branch/F");
  mytree->Branch("genelec_et_branch",&genelec_et_var,"genelec_et_branch/F");
  mytree->Branch("genposi_e_branch",&genposi_e_var,"genposi_e_branch/F");
  mytree->Branch("genposi_eta_branch",&genposi_eta_var,"genposi_eta_branch/F");
  mytree->Branch("genposi_phi_branch",&genposi_phi_var,"genposi_phi_branch/F");
  mytree->Branch("genposi_et_branch",&genposi_et_var,"genposi_et_branch/F");
  mytree->Branch("genelec_hassc_branch",&genelec_hassc_var,"genelec_hassc_branch/I");
  mytree->Branch("genposi_hassc_branch",&genposi_hassc_var,"genposi_hassc_branch/I");

  
  //generated variables for the tree (before FSR)
  mytree->Branch("unstablegenelec_e_branch",&unstablegenelec_e_var,"unstablegenelec_e_branch/F");
  mytree->Branch("unstablegenelec_eta_branch",&unstablegenelec_eta_var,"unstablegenelec_eta_branch/F");
  mytree->Branch("unstablegenelec_phi_branch",&unstablegenelec_phi_var,"unstablegenelec_phi_branch/F");
  mytree->Branch("unstablegenelec_et_branch",&unstablegenelec_et_var,"unstablegenelec_et_branch/F");
  mytree->Branch("unstablegenposi_e_branch",&unstablegenposi_e_var,"unstablegenposi_e_branch/F");
  mytree->Branch("unstablegenposi_eta_branch",&unstablegenposi_eta_var,"unstablegenposi_eta_branch/F");
  mytree->Branch("unstablegenposi_phi_branch",&unstablegenposi_phi_var,"unstablegenposi_phi_branch/F");
  mytree->Branch("unstablegenposi_et_branch",&unstablegenposi_et_var,"unstablegenposi_et_branch/F");
  
  //generated variables for the tree (Z variables)
  mytree->Branch("genboson_m_branch",&genboson_m_var,"genboson_m_branch/F");
  mytree->Branch("genboson_eta_branch",&genboson_eta_var,"genboson_eta_branch/F");
  mytree->Branch("genboson_phi_branch",&genboson_phi_var,"genboson_phi_branch/F");
  mytree->Branch("genboson_e_branch",&genboson_e_var,"genboson_e_branch/F");
  mytree->Branch("genboson_et_branch",&genboson_et_var,"genboson_et_branch/F");
  mytree->Branch("genboson_ez_branch",&genboson_ez_var,"genboson_ez_branch/F");
  mytree->Branch("genboson_p_branch",&genboson_p_var,"genboson_p_branch/F");
  mytree->Branch("genboson_pt_branch",&genboson_pt_var,"genboson_pt_branch/F");
  mytree->Branch("genboson_pz_branch",&genboson_pz_var,"genboson_pz_branch/F");

  mytree->Branch("x1quark",&x1quark,"x1quark/F");
  mytree->Branch("x2quark",&x2quark,"x2quark/F");

  //FSR variables
  mytree->Branch("fsrposiphotonsize",&fsrposiphotonsize,"fsrposiphotonsize/I");
  mytree->Branch("fsrelecphotonsize",&fsrelecphotonsize,"fsrelecphotonsize/I");
  
  mytree->Branch("energyfsrelec",&energyfsrelec,"energyfsrelec[fsrelecphotonsize]/F");
  mytree->Branch("etfsrelec",&etfsrelec,"etfsrelec[fsrelecphotonsize]/F");
  mytree->Branch("etafsrelec",&etafsrelec,"etafsrelec[fsrelecphotonsize]/F");
  mytree->Branch("phifsrelec",&phifsrelec,"phifsrelec[fsrelecphotonsize]/F");
  
  mytree->Branch("energyfsrposi",&energyfsrposi,"energyfsrposi[fsrposiphotonsize]/F");
  mytree->Branch("etfsrposi",&etfsrposi,"etfsrposi[fsrposiphotonsize]/F");
  mytree->Branch("etafsrposi",&etafsrposi,"etafsrposi[fsrposiphotonsize]/F");
  mytree->Branch("phifsrposi",&phifsrposi,"phifsrposi[fsrposiphotonsize]/F");
  

  //supercluster matching generated electrons branches for the tree
  mytree->Branch("scelecenergy",&scelecenergy,"scelecenergy/F");
  mytree->Branch("sceleceta",&sceleceta,"sceleceta/F");
  mytree->Branch("scelecphi",&scelecphi,"scelecphi/F");
  mytree->Branch("scelecgsfmatched",&scelecgsfmatched,"scelecgsfmatched/F");
  mytree->Branch("scposienergy",&scposienergy,"scposienergy/F");
  mytree->Branch("scposieta",&scposieta,"scposieta/F");
  mytree->Branch("scposiphi",&scposiphi,"scposiphi/F");
  mytree->Branch("scposigsfmatched",&scposigsfmatched,"scposigsfmatched/F");
  



  mytree->Branch("genelechassc",&genelechassc,"genelechassc/O");
  //mytree->Branch("genelechasgsf",&genelechasgsf,"genelechasgsf/O");

  mytree->Branch("genposihassc",&genposihassc,"genposihassc/O");
  //mytree->Branch("genposihasgsf",&genposihasgsf,"genposihasgsf/O");


  //GSF Electron variables
  //GSF 
  mytree->Branch("gsf_size",&gsf_size, "gsf_size/I");
  mytree->Branch("gsf_theta", gsf_theta, "gsf_theta[gsf_size]/F");
  mytree->Branch("gsf_isEB", gsf_isEB, "gsf_isEB[gsf_size]/I");
  mytree->Branch("gsf_isEE", gsf_isEE, "gsf_isEE[gsf_size]/I");
  mytree->Branch("gsf_deltaEtaATcalo", gsf_deltaEtaATcalo, "gsf_deltaEtaATcalo[gsf_size]/F");
  mytree->Branch("gsf_deltaPhiATcalo", gsf_deltaPhiATcalo, "gsf_deltaPhiATcalo[gsf_size]/F");
  mytree->Branch("gsf_ecalEnergy", gsf_ecalEnergy, "gsf_ecalEnergy[gsf_size]/F");
  mytree->Branch("gsf_eOVERp", gsf_eOVERp, "gsf_eOVERp[gsf_size]/F");
  mytree->Branch("gsf_dxy", gsf_dxy, "gsf_dxy[gsf_size]/F");
  mytree->Branch("gsf_vz", gsf_vz, "gsf_vz[gsf_size]/F");
  mytree->Branch("gsf_nHits", gsf_nHits, "gsf_nHits[gsf_size]/I");
  mytree->Branch("gsf_nLostInnerHits", gsf_nLostInnerHits, "gsf_nLostInnerHits[gsf_size]/I");
  mytree->Branch("gsf_nLostOuterHits", gsf_nLostOuterHits, "gsf_nLostOuterHits[gsf_size]/I");
  mytree->Branch("gsf_fBrem", gsf_fBrem, "gsf_fBrem[gsf_size]/F");
  mytree->Branch("gsf_e1x5", gsf_e1x5, "gsf_e1x5[gsf_size]/F");
  mytree->Branch("gsf_e2x5", gsf_e2x5, "gsf_e2x5[gsf_size]/F");
  mytree->Branch("gsf_e5x5", gsf_e5x5, "gsf_e5x5[gsf_size]/F");
  //mytree->Branch("gsf_e1OVERe9", gsf_e1OVERe9, "gsf_e1OVERe9[gsf_size]/F");

  mytree->Branch("gsf_eMax", gsf_eMax, "gsf_eMax[gsf_size]/F");
  mytree->Branch("gsf_SwissCross", gsf_SwissCross, "gsf_SwissCross[gsf_size]/F");
  
  mytree->Branch("gsf_e1x3",gsf_e1x3,"gsf_e1x3[gsf_size]/F");
  mytree->Branch("gsf_e3x1",gsf_e3x1,"gsf_e3x1[gsf_size]/F");
  mytree->Branch("gsf_e2x2",gsf_e2x2,"gsf_e2x2[gsf_size]/F");
  mytree->Branch("gsf_e3x2",gsf_e3x2,"gsf_e3x2[gsf_size]/F");
  mytree->Branch("gsf_e3x3",gsf_e3x3,"gsf_e3x3[gsf_size]/F");
  mytree->Branch("gsf_e4x4",gsf_e4x4,"gsf_e4x4[gsf_size]/F");
  mytree->Branch("gsf_e2x5Right",gsf_e2x5Right,"gsf_e2x5Right[gsf_size]/F");
  mytree->Branch("gsf_e2x5Left",gsf_e2x5Left,"gsf_e2x5Left[gsf_size]/F");  
  mytree->Branch("gsf_e2x5Top",gsf_e2x5Top,"gsf_e2x5Top[gsf_size]/F");  
  mytree->Branch("gsf_e2x5Bottom",gsf_e2x5Bottom,"gsf_e2x5Bottom[gsf_size]/F");
  mytree->Branch("gsf_e2x5Max",gsf_e2x5Max,"gsf_e2x5Max[gsf_size]/F");
  mytree->Branch("gsf_eLeft",gsf_eLeft,"gsf_eLeft[gsf_size]/F");
  mytree->Branch("gsf_eRight",gsf_eRight,"gsf_eRight[gsf_size]/F");
  mytree->Branch("gsf_eTop",gsf_eTop,"gsf_eTop[gsf_size]/F");
  mytree->Branch("gsf_eBottom",gsf_eBottom,"gsf_eBottom[gsf_size]/F");
  mytree->Branch("gsf_e2nd",gsf_e2nd,"gsf_e2nd[gsf_size]/F");

  mytree->Branch("gsf_p",&gsf_p,"gsf_p[gsf_size]/F");
  mytree->Branch("gsf_e",&gsf_e,"gsf_e[gsf_size]/F");
  mytree->Branch("gsf_pt",&gsf_pt,"gsf_pt[gsf_size]/F");
  mytree->Branch("gsf_class",&gsf_class,"gsf_class[gsf_size]/F");
  mytree->Branch("gsf_e2x5overe5x5",&gsf_e2x5overe5x5,"gsf_e2x5overe5x5[gsf_size]/F");
  mytree->Branch("gsf_e1x5overe5x5",&gsf_e1x5overe5x5,"gsf_e1x5overe5x5[gsf_size]/F");
  mytree->Branch("gsf_eta",&gsf_eta,"gsf_eta[gsf_size]/F");
  mytree->Branch("gsf_phi",&gsf_phi,"gsf_phi[gsf_size]/F");
  mytree->Branch("gsf_px",&gsf_px,"gsf_px[gsf_size]/F");
  mytree->Branch("gsf_py",&gsf_py,"gsf_py[gsf_size]/F");
  mytree->Branch("gsf_pz",&gsf_pz,"gsf_pz[gsf_size]/F");

  mytree->Branch("gsf_deltaeta",&gsf_deltaeta,"gsf_deltaeta[gsf_size]/F");
  mytree->Branch("gsf_deltaphi",&gsf_deltaphi,"gsf_deltaphi[gsf_size]/F");
  mytree->Branch("gsf_hovere",&gsf_hovere,"gsf_hovere[gsf_size]/F");
   
  mytree->Branch("gsf_trackiso",&gsf_trackiso,"gsf_trackiso[gsf_size]/F");
  mytree->Branch("gsf_ecaliso",&gsf_ecaliso,"gsf_ecaliso[gsf_size]/F");
  mytree->Branch("gsf_hcaliso1",&gsf_hcaliso1,"gsf_hcaliso1[gsf_size]/F");
  mytree->Branch("gsf_hcaliso2",&gsf_hcaliso2,"gsf_hcaliso2[gsf_size]/F");

  mytree->Branch("gsf_charge",&gsf_charge,"gsf_charge[gsf_size]/I");
  mytree->Branch("gsf_sigmaetaeta",&gsf_sigmaetaeta,"gsf_sigmaetaeta[gsf_size]/F");
  mytree->Branch("gsf_sigmaIetaIeta",&gsf_sigmaIetaIeta,"gsf_sigmaIetaIeta[gsf_size]/F");
  mytree->Branch("gsf_isecaldriven",&gsf_isecaldriven,"gsf_isecaldriven[gsf_size]/F");
  mytree->Branch("gsf_istrackerdriven",&gsf_istrackerdriven,"gsf_istrackerdriven[gsf_size]/F");
	  
  mytree->Branch("gsfsc_e",&gsfsc_e,"gsfsc_e[gsf_size]/F");
  mytree->Branch("gsfsc_pt",&gsfsc_pt,"gsfsc_pt[gsf_size]/F");
  mytree->Branch("gsfsc_eta",&gsfsc_eta,"gsfsc_eta[gsf_size]/F");
  mytree->Branch("gsfsc_phi",&gsfsc_phi,"gsfsc_phi[gsf_size]/F");
  mytree->Branch("gsfsc_px",&gsfsc_px,"gsfsc_px[gsf_size]/F");
  mytree->Branch("gsfsc_py",&gsfsc_py,"gsfsc_py[gsf_size]/F");
  mytree->Branch("gsfsc_pz",&gsfsc_pz,"gsfsc_pz[gsf_size]/F");
      
  mytree->Branch("gsf_gsfet",&gsf_gsfet,"gsf_gsfet[gsf_size]/F");

  mytree->Branch("scindexforgsf",&scindexforgsf,"scindexforgsf[gsf_size]/I");

  mytree->Branch("gsfindexforgenelec",&gsfindexforgenelec,"gsfindexforgenelec/I"); 
  mytree->Branch("gsfindexforgenposi",&gsfindexforgenposi,"gsfindexforgenposi/I"); 

  mytree->Branch("scindexforgenelec",&scindexforgenelec,"scindexforgenelec/I"); 
  mytree->Branch("scindexforgenposi",&scindexforgenposi,"scindexforgenposi/I"); 

  mytree->Branch("gsfpass_ET",&gsfpass_ET,"gsfpass_ET[gsf_size]/O"); 
  mytree->Branch("gsfpass_PT",&gsfpass_PT,"gsfpass_PT[gsf_size]/O");  
  mytree->Branch("gsfpass_DETETA",&gsfpass_DETETA,"gsfpass_DETETA[gsf_size]/O");  
  mytree->Branch("gsfpass_CRACK",&gsfpass_CRACK,"gsfpass_CRACK[gsf_size]/O"); 
  mytree->Branch("gsfpass_DETAIN",&gsfpass_DETAIN,"gsfpass_DETAIN[gsf_size]/O");  
  mytree->Branch("gsfpass_DPHIIN",&gsfpass_DPHIIN,"gsfpass_DPHIIN[gsf_size]/O");  
  mytree->Branch("gsfpass_HADEM",&gsfpass_HADEM,"gsfpass_HADEM[gsf_size]/O");  
  mytree->Branch("gsfpass_SIGMAIETAIETA",&gsfpass_SIGMAIETAIETA,"gsfpass_SIGMAIETAIETA[gsf_size]/O"); 
  mytree->Branch("gsfpass_E2X5OVER5X5",&gsfpass_E2X5OVER5X5,"gsfpass_E2X5OVER5X5[gsf_size]/O");  
  mytree->Branch("gsfpass_ISOLEMHADDEPTH1",&gsfpass_ISOLEMHADDEPTH1,"gsfpass_ISOLEMHADDEPTH1[gsf_size]/O");  
  mytree->Branch("gsfpass_ISOLHADDEPTH2",&gsfpass_ISOLHADDEPTH2,"gsfpass_ISOLHADDEPTH2[gsf_size]/O");  
  mytree->Branch("gsfpass_ISOLPTTRKS",&gsfpass_ISOLPTTRKS,"gsfpass_ISOLPTTRKS[gsf_size]/O");  
  mytree->Branch("gsfpass_ECALDRIVEN",&gsfpass_ECALDRIVEN,"gsfpass_ECALDRIVEN[gsf_size]/O");  
  mytree->Branch("gsfpass_INVALID",&gsfpass_INVALID,"gsfpass_INVALID[gsf_size]/O");  
  mytree->Branch("gsfpass_NOMISSINGHITS",&gsfpass_NOMISSINGHITS,"gsfpass_NOMISSINGHITS[gsf_size]/O");  

  mytree->Branch("gsfpass_HEEP",&gsfpass_HEEP,"gsfpass_HEEP[gsf_size]/O");  

  mytree->Branch("gsfpass_ID",&gsfpass_ID,"gsfpass_ID[gsf_size]/O");  
  mytree->Branch("gsfpass_ISO",&gsfpass_ISO,"gsfpass_ISO[gsf_size]/O");  

  //charge info
  mytree->Branch("scpixcharge",&scpixcharge,"scpixcharge[gsf_size]/I");
  mytree->Branch("ctfcharge",&ctfcharge,"ctfcharge[gsf_size]/I");
  mytree->Branch("gsfcharge",&gsfcharge,"gsfcharge[gsf_size]/I");
  mytree->Branch("gsfctfscpixconsistent",&gsfctfscpixconsistent,"gsfctfscpixconsistent[gsf_size]/O");
  mytree->Branch("gsfscpixconsistent",&gsfscpixconsistent,"gsfscpixconsistent[gsf_size]/O");
  mytree->Branch("gsfctfconsistent",&gsfctfconsistent,"gsfctfconsistent[gsf_size]/O");


  //Gsf Track information
  mytree->Branch("gsftracksize",&gsftracksize,"gsftracksize/I");

  mytree->Branch("gsftracketa",&gsftracketa,"gsftracketa[gsftracksize]/F");
  mytree->Branch("gsftrackphi",&gsftrackphi,"gsftrackphi[gsftracksize]/F");  
  mytree->Branch("gsftrackp",&gsftrackp,"gsftrackp[gsftracksize]/F");
  mytree->Branch("gsftrackpt",&gsftrackpt,"gsftrackpt[gsftracksize]/F");
  mytree->Branch("gsftrackpx",&gsftrackpx,"gsftrackpx[gsftracksize]/F");
  mytree->Branch("gsftrackpy",&gsftrackpy,"gsftrackpy[gsftracksize]/F");
  mytree->Branch("gsftrackpz",&gsftrackpz,"gsftrackpz[gsftracksize]/F");
  
  

}


// ------------ method called once each job just after ending the event loop  ------------
void 
GsfCheckerTree::endJob() {

  cout<<"GsfCheckerTree::endJob"<<endl;
  //rootfile->cd();
  //rootfile->Write();
  //rootfile->Close();


}









//
// -----------------------------------------------------------
// -----------------------------------------------------------
void GsfCheckerTree::datagenerated(const edm::Event& e) {

  using namespace std;
  using namespace edm;
  

  bool debug = false;

  //Informations about generated data

  //   cout<<"i am in datagenerated "<<endl;
  //   cout<<debug<<endl;

  //   cout<<"before printing "<<endl;

  edm::Handle<HepMCProduct> EvtHandle ;
  e.getByLabel( "generator", EvtHandle ) ;
  const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
  HepMC::GenVertex* GravitonDecVtx = 0 ;
	 
  //   cout<<"before printing "<<endl;
  //   cout<<e.id().event()<<endl;
  //   if  ( e.id().event() <= 10 ) Evt->print();
  //   cout<<"after printing "<<endl;

  //search for e coming from DY or ZP or GR

  for ( HepMC::GenEvent::vertex_const_iterator
	  vit=Evt->vertices_begin(); vit!=Evt->vertices_end(); vit++ )
    {
      for ( HepMC::GenVertex::particles_out_const_iterator
	      pout=(*vit)->particles_out_const_begin();
	    pout!=(*vit)->particles_out_const_end(); pout++ )
	{
	  //if (  (  ( (*pout)->pdg_id() == 32 ) || ( (*pout)->pdg_id() == 23 ) || ( (*pout)->pdg_id() == 22 ) ) && ( (*pout)->status() == 3))	     {	    
	  if ( ((*pout)->pdg_id() == 23) && ((*pout)->status() == 3) ) {
	    if ( (*pout)->end_vertex() != 0 )
	      {
		GravitonDecVtx = (*pout)->end_vertex() ;
		HepMC::FourVector truc = (*pout)->momentum();
		momboson = (*pout)->momentum();
		break ;
	      }
	  }
	}
      if ( GravitonDecVtx != 0 )
	{
	  break ; // break the initial loop over vertices
	}
    }

   
  if ( GravitonDecVtx == 0 ) 
    {
      //cout << " There is NO ZP or Z/g or Graviton in this event ! " << endl ;
      return ;
    }
    
  // Print out informations for the first 10 events 
  int icount=0;
  if  ((debug) &&  e.id().event() <= 10 )
    {
      for( HepMC::GenEvent::particle_const_iterator partIter = Evt->particles_begin() ;
	   partIter != Evt->particles_end() ; ++partIter) //loop over generator particles
	{
	  icount++;
	  if(icount<20) {
	    HepMC::FourVector p = ( *partIter )->momentum() ;
	    //cout<<( *partIter )->momentum()<<endl;
	    // 	    std::cout << " i=  " <<  icount << " status= " << ( *partIter )->status();
	    // 	    std::cout << " id= " << ( *partIter )->pdg_id() ;
	    // 	    std::cout << " m = " <<std::fabs( ( *partIter )->momentum().m() ) ;
	    // 	    std::cout << " e = " <<std::abs( p.e() ) ;
	    // 	    std::cout << " et = " <<std::abs( p.e() * sin(p.theta()) ) ;
	    // 	    std::cout << " pz = " <<std::abs( p.pz() ) ;
	    // 	    std::cout << " eta = " <<std::abs( p.eta() ) ;
	    // 	    std::cout << " phi = " <<std::abs( p.phi() )  << std::endl; 
	  }
	}
    }


  // Print out informations to check out everything is ok  
  if  ((debug) &&  e.id().event() <= 10 )
    {
      //       cout << " ZP or Z/g or Graviton decay found at the vertex " << GravitonDecVtx->barcode() <<" (barcode)" << endl ;
      HepMC::GenVertex::particles_in_const_iterator Gin = GravitonDecVtx->particles_in_const_begin();

      //vector<HepMC::GenParticle*> GravitonChildren = (*Gin)->listChildren() ;      
      vector<HepMC::GenParticle*> GravitonChildren;
      HepMC::GenVertex* outVertex = (*Gin)->end_vertex();


      for(std::vector< HepMC::GenParticle*>::const_iterator iter = outVertex->particles_in_const_begin();
	  //for(std::set< HepMC::GenParticle*>::const_iterator iter = outVertex->particles_in_const_begin();
	  iter != outVertex->particles_in_const_end();iter++)
	GravitonChildren.push_back(*iter);

      //       cout << " Number of ZP or Z/g or Graviton (immediate) children = " << GravitonChildren.size() << endl ;
      //       for (unsigned int ic=0; ic<GravitonChildren.size(); ic++ ) {
      // 	GravitonChildren[ic]->print() ;   }
    }


  // select and store stable descendants of the DY
  //   

  //   vector<GenParticle*> test;

  vector<HepMC::GenParticle*> StableDYDesc ;
  vector<HepMC::GenParticle*> UnstableDYDesc ;
  vector<HepMC::GenParticle*> photonsfsrelec ;
  vector<HepMC::GenParticle*> photonsfsrposi ;
  
  for ( HepMC::GenVertex::particle_iterator
	  des=GravitonDecVtx->particles_begin(HepMC::descendants);
	des!=GravitonDecVtx->particles_end(HepMC::descendants); des++ )
    {
      //       if ( (debug) &&  e.id().event() <=10 ) (*des)->print() ;
      if ( (*des)->status() == 1 ) StableDYDesc.push_back(*des) ;  //means the particle is stable
      if ( (*des)->status() == 3 ) UnstableDYDesc.push_back(*des) ;//means the particle is going to decay
    }
   
     
  // browse the array of STABLE descendants
  // and do 2-e inv.mass
  //
  for ( unsigned int i=0; i<StableDYDesc.size(); i++ )
    {
      // skip other than electron
      if ( abs(StableDYDesc[i]->pdg_id()) != 11 ) continue ; 
	 
      // Searching for unstable electron/positron pairs
      for ( unsigned int j=i+1; j<StableDYDesc.size(); j++ )
	{
	  // skip other than electron
	  if ( abs(StableDYDesc[j]->pdg_id()) != 11 ) continue ;
	  //
	  // skip same charge combo's
	  //
	  if ( (StableDYDesc[i]->pdg_id()*StableDYDesc[j]->pdg_id()) > 0 ) 
	    continue ;
	  //
	  // OK, opposite charges, do the job
	  //skip unstable descendants 
	  if (  (StableDYDesc[i]->status() != 1) || (StableDYDesc[j]->status() != 1)  )
	    continue;
	  //
	  momelec = StableDYDesc[i]->momentum();
	  momposi = StableDYDesc[j]->momentum();

	  // 	  if(debug) cout << " stable counters : " << StableDYDesc[i]->barcode() << " " 
	  // 			 << StableDYDesc[j]->barcode() << endl ;
	}
    }
     
     
  // Searching for UNSTABLE electron/positron pairs
     
  for ( unsigned int i=0; i<UnstableDYDesc.size(); i++ )
    {
      // skip other than electron
      if ( abs(UnstableDYDesc[i]->pdg_id()) != 11 ) continue ; 
	 
      for ( unsigned int j=i+1; j<UnstableDYDesc.size(); j++ )
	{
	  // skip other than electron
	  if ( abs(UnstableDYDesc[j]->pdg_id()) != 11 ) continue ;
	  //
	  // skip same charge combo's
	  if ( (UnstableDYDesc[i]->pdg_id()*UnstableDYDesc[j]->pdg_id()) > 0 ) 
	    continue ;
	  //
	  // OK, opposite charges, do the job
	  // keep unstable descendants 
	  if (  (UnstableDYDesc[i]->status() != 3) || (UnstableDYDesc[j]->status() != 3)  )
	    continue;
	    
	  unstablemomelec = UnstableDYDesc[i]->momentum();
	  unstablemomposi = UnstableDYDesc[j]->momentum();

	  HepMC::GenVertex* ElecDecVtx = UnstableDYDesc[i]->end_vertex();
	  HepMC::GenVertex* PosiDecVtx = UnstableDYDesc[j]->end_vertex();

	  // 	  for(std::set< HepMC::GenParticle*>::const_iterator elecdesc = ElecDecVtx->particles_out_const_begin();
	  // 	     elecdesc  != ElecDecVtx->particles_out_const_end();elecdesc++) {
	  // 	    cout<<"id of desc is "<<(*elecdesc)->pdg_id()
	  // 		<<" status "<<(*elecdesc)->status()
	  // 		<<" energy "<<(*elecdesc)->momentum().e()
	  // 		<<"  "<<(*elecdesc)->momentum().eta()
	  // 		<<"  "<<(*elecdesc)->momentum().phi()
	  // 		<<endl;
	  // // 		<<"  "<<(*elecdesc)->momentum()
	  // // 		<<"  "<<(*elecdesc)->momentum()
	  // 	  }

	  int fsreleccounter = 0;
	  int fsrposicounter = 0;

	  for ( HepMC::GenVertex::particle_iterator
		  fsrelec=ElecDecVtx->particles_begin(HepMC::descendants);
		fsrelec!=ElecDecVtx->particles_end(HepMC::descendants); fsrelec++ )
	    {
	      if((*fsrelec)->pdg_id() == 22) {
		photonsfsrelec.push_back(*fsrelec);
		energyfsrelec[fsreleccounter] = (*fsrelec)->momentum().e();
		etafsrelec[fsreleccounter] = (*fsrelec)->momentum().eta();
		phifsrelec[fsreleccounter] = (*fsrelec)->momentum().phi();
		etfsrelec[fsreleccounter] = (*fsrelec)->momentum().perp();
		// 		energyfsrtest[fsrcounter] = (*fsrelec)->momentum().e();
		fsreleccounter++;
		// 		cout<<"id of desc is "<<(*fsrelec)->pdg_id()
		// 		    <<" status "<<(*fsrelec)->status()
		// 		    <<" energy "<<(*fsrelec)->momentum().e()
		// 		    <<" eta "<<(*fsrelec)->momentum().eta()
		// 		    <<" phi  "<<(*fsrelec)->momentum().phi()
		// 		    <<endl;
		// 		<<"  "<<(*fsrelec)->momentum()
		// 		<<"  "<<(*fsrelec)->momentum()
	      }
	    }

	  for ( HepMC::GenVertex::particle_iterator
		  fsrposi=PosiDecVtx->particles_begin(HepMC::descendants);
		fsrposi!=PosiDecVtx->particles_end(HepMC::descendants); fsrposi++ )
	    {
	      if((*fsrposi)->pdg_id() == 22) {
		photonsfsrposi.push_back(*fsrposi);
		energyfsrposi[fsrposicounter] = (*fsrposi)->momentum().e();
		etafsrposi[fsrposicounter] = (*fsrposi)->momentum().eta();
		phifsrposi[fsrposicounter] = (*fsrposi)->momentum().phi();
		etfsrposi[fsrposicounter] = (*fsrposi)->momentum().perp();
		// 		energyfsrtest[fsrcounter] = (*fsrposi)->momentum().e();
		fsrposicounter++;
		// 		cout<<"id of desc is "<<(*fsrposi)->pdg_id()
		// 		    <<" status "<<(*fsrposi)->status()
		// 		    <<" energy "<<(*fsrposi)->momentum().e()
		// 		    <<" eta "<<(*fsrposi)->momentum().eta()
		// 		    <<" phi  "<<(*fsrposi)->momentum().phi()
		// 		    <<endl;
		// 		<<"  "<<(*fsrposi)->momentum()
		// 		<<"  "<<(*fsrposi)->momentum()
	      }
	    }

	  fsrelecphotonsize = photonsfsrelec.size();
	  fsrposiphotonsize = photonsfsrposi.size();
	  // 	  cout<<" number of fsrelec "<<fsrelecphotonsize<<endl;
	  // 	  cout<<" number of fsrposi "<<fsrposiphotonsize<<endl;


	  // 	  if(debug) cout << " unstable counters : " << UnstableDYDesc[i]->barcode() << " " 
	  // 			 << UnstableDYDesc[j]->barcode() <<endl;
	}
    }




  //   if(  momelec.e()/unstablemomelec.e() < 0.8 ) fsrposi = true;
  //   if(  momposi.e()/unstablemomposi.e() < 0.8 ) fsrposi = true;

  //   HepMC::FourVector elecaddposigsfgsf = 0;
  //   HepMC::FourVector unstelecaddposigsfgsf = 0;

  elecaddposigsfgsf.setE(momposi.e()+momelec.e());
  elecaddposigsfgsf.setPx(momposi.px()+momelec.px());
  elecaddposigsfgsf.setPy(momposi.py()+momelec.py());
  elecaddposigsfgsf.setPz(momposi.pz()+momelec.pz());
    
  unstelecaddposigsfgsf.setE(unstablemomposi.e()+unstablemomelec.e());
  unstelecaddposigsfgsf.setPx(unstablemomposi.px()+unstablemomelec.px());
  unstelecaddposigsfgsf.setPy(unstablemomposi.py()+unstablemomelec.py());
  unstelecaddposigsfgsf.setPz(unstablemomposi.pz()+unstablemomelec.pz());
    
  //   if(debug) {
  //     cout<< "gen electron info: e, et, eta, phi  "<<endl;
  //     cout<< "electron:  " << momelec.e()<<" " << momelec.e() * sin(momelec.theta())<<" "<<momelec.eta()<<" "<<momelec.phi()<<endl;
  //     cout<< "positrons: " << momposi.e()<<" " << momposi.e() * sin(momposi.theta())<<" "<<momposi.eta()<<" "<<momposi.phi()<<endl;
  
  //     cout<<"the ee inv mass is "<<elecaddposigsfgsf.m()<<endl;

  //     cout<< "electron UNSTABLE:  " 
  // 	<< unstablemomelec.e()<<" "
  // 	<<unstablemomelec.e() * sin(unstablemomelec.theta())<<" "
  // 	<<unstablemomelec.eta()<<" "
  // 	<<unstablemomelec.phi()<<endl;
  //     cout<< "positrons UNSTABLE: " 
  // 	<< unstablemomposi.e()<<" "
  // 	<<unstablemomposi.e() * sin(unstablemomposi.theta())<<" "
  // 	<<unstablemomposi.eta()<<" "
  // 	<<unstablemomposi.phi()<<endl;

  //     cout<<"the Z inv mass is "<<unstelecaddposigsfgsf.m()<<endl;
  //   }

}//end of datagenerated





const EcalRecHitCollection * GsfCheckerTree::getEcalRecHitCollection( const reco::BasicCluster &cluster )
{
  if ( cluster.size() == 0 ) {
    throw cms::Exception("InvalidCluster") << "The cluster has no crystals!";
  }
  DetId id = (cluster.hitsAndFractions()[0]).first; // size is by definition > 0 -- FIXME??
  const EcalRecHitCollection *recHits = 0;
  if ( id.subdetId() == EcalBarrel ) {
    recHits = ebRecHits_;
  } else if ( id.subdetId() == EcalEndcap ) {
    recHits = eeRecHits_;
  } else {
    throw cms::Exception("InvalidSubdetector") << "The subdetId() " << id.subdetId() << " does not correspond to EcalBarrel neither EcalEndcap";
  }
  return recHits;
}















