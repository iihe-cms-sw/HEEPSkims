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
//
//Cleaning ladies : Thomas and Laurent
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "UserCode/HEEPSkims/interface/GsfCheckerTree.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
//#include "SHarper/HEEPAnalyzer/interface/HEEPDebug.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
//#include "SHarper/HEEPAnalyzer/interface/HEEPEventHelper.h"
//#include "SHarper/HEEPAnalyzer/interface/HEEPEvent.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
//ECAL 
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/DetId/interface/DetId.h"


#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include <TMath.h>
#define PI 3.141592654
#define TWOPI 6.283185308

using namespace std;
using namespace reco;
using namespace edm;

//Method to sort the gsf electrons
bool 
gsfEtGreater(const reco::GsfElectron &gsf1,const reco::GsfElectron &gsf2)
{
  //float et1 = gsf1.superCluster()->energy() * sin(gsf1.theta());
  //float et2 = gsf1.superCluster()->energy() * sin(gsf2.theta());
  float et1 = gsf1.caloEnergy() * sin(gsf1.p4().theta());
  float et2 = gsf2.caloEnergy() * sin(gsf2.p4().theta());
  return (et1 > et2);
}

bool 
scEtGreater(const reco::SuperCluster *sc1,const reco::SuperCluster *sc2) 
{
  return (((sc1->energy() + sc1->preshowerEnergy()) )/cosh((sc1)->eta()) >((sc2->energy() + sc2->preshowerEnergy()) )/cosh((sc2)->eta()) );
	
}


bool 
refScEtGreater(reco::SuperClusterRef sc1,reco::SuperClusterRef sc2) 
{
  return (((sc1->energy() + sc1->preshowerEnergy()) )/cosh((sc1)->eta()) >((sc2->energy() + sc2->preshowerEnergy()) )/cosh((sc2)->eta()));
  
}

// bool 
// MuonEtGreater(const reco::Muon *sc1,const reco::SuperCluster*sc2) 
// {
//   return (((sc1->energy() + sc1->preshowerEnergy()) )/cosh((*sc1)->eta()) >((sc2->energy() + sc2->preshowerEnergy()) )/cosh((*sc2)->eta()) );
	
// }

float 
etacorr(float eta, float pvz, float scz) 
{
  return asinh(sinh(eta) * (1. - pvz/scz));
}

float
GsfCheckerTree::CalcInvariantMass(const int& iEle1, const int& iEle2)
{
  TLorentzVector ele1;
  TLorentzVector ele2;

  float et1 = gsfsc_e[iEle1] * sin(gsf_theta[iEle1]);
  float et2 = gsfsc_e[iEle2] * sin(gsf_theta[iEle2]);

  ele1.SetPtEtaPhiE(et1, gsf_eta[iEle1], gsf_phi[iEle1], (et1 * cosh(gsf_eta[iEle1])));
  ele2.SetPtEtaPhiE(et2, gsf_eta[iEle2], gsf_phi[iEle2], (et2 * cosh(gsf_eta[iEle2])));

  return (ele1+ele2).Mag();
}

GsfCheckerTree::GsfCheckerTree(const edm::ParameterSet& iConfig):
  //evtHelper_(),heepEvt_(),
  nrPass_(0),nrFail_(0)
{
  //evtHelper_.setup(iConfig);
  //now do what ever initialization is needed
  eventcounter = 0;

  hlTriggerResults_ = iConfig.getParameter<edm::InputTag>("TriggerResultsTag");
  comEnergy_ = iConfig.getParameter<double>("centerOfMassEnergy");
  ScPtMin_ = iConfig.getUntrackedParameter<double>("ScPtMin", 10.);
  bJetPtMin_ = iConfig.getUntrackedParameter<double>("bJetPtMin", 10.);
  bJetEtaMax_ = iConfig.getUntrackedParameter<double>("bJetEtaMax", 3.);
  jetPtMin_ = iConfig.getUntrackedParameter<double>("jetPtMin", 10.);
  jetEtaMax_ = iConfig.getUntrackedParameter<double>("jetEtaMax", 3.);
  GsfPtMin_ = iConfig.getUntrackedParameter<double>("GsfPtMin", 5.);
  GsfTrackPtMin_= iConfig.getUntrackedParameter<double>("GsfTrackPtMin", 5.);
  muPtMin_ = iConfig.getUntrackedParameter<double>("muPtMin", 5.);
  ele1EtMin_ = iConfig.getUntrackedParameter<double>("electron1EtMin", 0.);
  ele1EtMax_ = iConfig.getUntrackedParameter<double>("electron1EtMax", 1.E99);
  ele2EtMin_ = iConfig.getUntrackedParameter<double>("electron2EtMin", 0.);
  ele2EtMax_ = iConfig.getUntrackedParameter<double>("electron2EtMax", 1.E99);
  muonPtMin_ = iConfig.getUntrackedParameter<double>("muonPtMin", 0.);
  muonPtMax_ = iConfig.getUntrackedParameter<double>("muonPtMax", 1.E99);

  // sanity check for cuts
  if (ele1EtMin_ > ele1EtMax_) {
    double helper = ele1EtMax_;
    ele1EtMax_ = ele1EtMin_;
    ele1EtMin_ = helper;
  }
  if (ele2EtMin_ > ele2EtMax_) {
    double helper = ele2EtMax_;
    ele2EtMax_ = ele2EtMin_;
    ele2EtMin_ = helper;
  }
  if (muonPtMin_ > muonPtMax_) {
    double helper = muonPtMax_;
    muonPtMax_ = muonPtMin_;
    muonPtMin_ = helper;
  }
  if (GsfPtMin_ > ele1EtMin_ || GsfPtMin_ > ele2EtMin_) GsfPtMin_ = (ele2EtMin_ > ele1EtMin_) ? ele1EtMin_ : ele2EtMin_;
  if (muPtMin_ > muonPtMin_) muPtMin_ = muonPtMin_;

  hcalCfg.hOverEConeSize = 0.15;
  hcalCfg.useTowers = true;
  hcalCfg.hcalTowers = edm::InputTag("towerMaker");
  hcalCfg.hOverEPtMin = 0;

  hcalHelper = new ElectronHcalHelper(hcalCfg);
  inputTagIsoDepElectrons_ = iConfig.getParameter< std::vector<edm::InputTag> >("IsoDepElectron");
  inputTagIsoValElectronsPFId_   = iConfig.getParameter< std::vector<edm::InputTag> >("IsoValElectronPF");
 
  EcalHcal1EffAreaEndcaps_ = iConfig.getUntrackedParameter<double>("EcalHcal1EffAreaEndcaps", 0.);
  EcalHcal1EffAreaBarrel_ = iConfig.getUntrackedParameter<double>("EcalHcal1EffAreaBarrel", 0.);
}


GsfCheckerTree::~GsfCheckerTree()
{
  //cout<<"GsfCheckerTree::~GsfCheckerTree"<<endl; 
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

// ------------ method called to for each event  ------------
void
GsfCheckerTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{


//   //Test Laurent 
//   evtHelper_.makeHeepEvent(iEvent,iSetup,heepEvt_);

//   // std::cout <<"event pt hat "<<heepEvt_.genEventPtHat()<<std::endl;

//   // std::cout <<"eh"<<std::endl;
 
//   //do what ever you want
//   //lets get the heep electrons and count the number that pass / fail cuts
//   const std::vector<heep::Ele>& eles = heepEvt_.heepEles();
//   for(size_t eleNr=0;eleNr<eles.size();eleNr++){
//     if(eles[eleNr].cutCode()==0x0) nrPass_++;
//     else nrFail_++;
//   }
//   cout << "nb fail "<< nrFail_ << endl; 
//   //End test Laurent 

  hcalHelper->checkSetup(iSetup);
  hcalHelper->readEvent(iEvent);

  bool useGenData_ = !iEvent.isRealData(); 
  
  eventcounter++;

  //Run and event number
  runnumber = iEvent.id().run();
  eventnumber = iEvent.id().event();
  luminosityBlock = iEvent.id().luminosityBlock();
  



  // for skim on pt 
  //Final GSF Electron collection
  edm::Handle<reco::GsfElectronCollection> pGsfElectrons;
  iEvent.getByLabel("gsfElectrons","",pGsfElectrons);
  reco::GsfElectronCollection gsfelectrons(pGsfElectrons->begin(),pGsfElectrons->end());

  //sort all the GSF by transverse energy
  std::sort(gsfelectrons.begin(),gsfelectrons.end(),gsfEtGreater);

  float gsfPtMax = 0.;
  float gsfPtSecondMax = 0.; 
 
  int cntr = 0; 
  for( reco::GsfElectronCollection::const_iterator gsfiterbis = gsfelectrons.begin(); gsfiterbis != gsfelectrons.end(); ++gsfiterbis) {
    cntr++; 
    
    if (cntr == 1) {
      //gsfPtMax = gsfiterbis->superCluster()->energy() * sin(gsfiterbis->theta());
      gsfPtMax = gsfiterbis->caloEnergy()*sin(gsfiterbis->p4().theta());
    }
    if (cntr == 2) {
      //gsfPtSecondMax = gsfiterbis->superCluster()->energy() * sin(gsfiterbis->theta());
      gsfPtSecondMax = gsfiterbis->caloEnergy()*sin(gsfiterbis->p4().theta());
    }
  }
  //Missing hits, Invariant Mass cut

  // Get MUONS
  edm::Handle<reco::MuonCollection> muonCollection;
  iEvent.getByLabel("muons",muonCollection);
  const reco::MuonCollection* muons = muonCollection.product();


  float muonPtMax = 0.;
  // get the two highes pt muons
  for(reco::MuonCollection::const_iterator muIt = muons->begin(); muIt < muons->end(); ++muIt){
    if (muIt->isGlobalMuon()) {
      // get TeV optimized track
      reco::Muon::MuonTrackTypePair tevOptMuTrk = muon::tevOptimized(*muIt, 200, 17., 40., 0.25);
      if (tevOptMuTrk.first->pt() > muonPtMax) {
        muonPtMax = tevOptMuTrk.first->pt();
      }
    }
  }

  // SKIMMING
  if (!(gsfPtMax >= ele1EtMin_ && gsfPtSecondMax >= ele2EtMin_) && !(gsfPtMax >= ele1EtMin_ && muonPtMax >= muonPtMin_)) {return;}
  if (gsfPtMax > ele1EtMax_ || gsfPtSecondMax > ele2EtMax_ || muonPtMax > muonPtMax_) {return;}
  
  //rho variable
  rho = 0;
  //rhoiso = 0;
  edm::Handle<double> rho_;
  //edm::Handle<double> rhoiso_;
  bool isrho; 
  //bool isrhoiso;
  isrho = iEvent.getByLabel(edm::InputTag("kt6PFJets:rho"),rho_);
  //isrhoiso =iEvent.getByLabel(edm::InputTag("kt6PFJetsIso:rho"),rhoiso_);
  if(isrho)   rho =*rho_;
  //if(isrhoiso) rhoiso =*rhoiso_;

  //edm::Handle<double> rho2012_;
  //double rho2012 =0;
  // iEvent.getByLabel(edm::InputTag("kt6PFJets_rho_RECO:double_kt6PFJets_rho_RECO.obj"),rho2012_);
  
  //rho2012 = *rho2012_; 
  //cout << "rho " <<  rho << endl;
  //beam spot variables
  sigmaZ = -5000.;
  sigmaZ0Error = -5000.;
  sq = -5000.; 
  bsposx = -5000.;
  bsposy = -5000.;
  bsposz = -5000.;
   
  trueNVtx = -5000.;
  nVtxBefore = -5000;
  nVtxNow = -5000;
  nVtxAfter = -5000;
 
  // generator information for MC samples
  if (useGenData_) {
    pthat = -5000.;
    alphaqcd = -5000.;
    alphaqed = -5000.;
    qscale = -5000.;
    processid = -5000;
    weight = -5000.;

    edm::Handle<GenEventInfoProduct> GenInfoHandle;
    bool genevtinfovalid = iEvent.getByLabel("generator",GenInfoHandle);
    
    if (genevtinfovalid) {
      pthat = GenInfoHandle->hasBinningValues() ? (GenInfoHandle->binningValues())[0] : 0.0 ;
      alphaqcd = GenInfoHandle->alphaQCD();
      alphaqed = GenInfoHandle->alphaQED();
      qscale = GenInfoHandle->qScale();
      processid = GenInfoHandle->signalProcessID();
      weight = GenInfoHandle->weight();
    }
    
    DataGenPart(iEvent);
    
    // pile up info
    edm::InputTag pileupSrc_("addPileupInfo");
    Handle<std::vector<PileupSummaryInfo> > puInfo;
    iEvent.getByLabel(pileupSrc_, puInfo);

    vector<PileupSummaryInfo>::const_iterator pvi;

    for (pvi = puInfo->begin(); pvi != puInfo->end(); ++pvi) {
      if (pvi->getBunchCrossing() == -1) nVtxBefore = pvi->getPU_NumInteractions();
      else if (pvi->getBunchCrossing() == 0) {
        trueNVtx = pvi->getTrueNumInteractions(); // from Fall11 onwards
        nVtxNow = pvi->getPU_NumInteractions();
      }
      else if (pvi->getBunchCrossing() == 1) nVtxAfter = pvi->getPU_NumInteractions();
      //cout << " Pileup Information: bunchXing, nvtx: " << pvi->getBunchCrossing() << " " << pvi->getTrueNumInteractions() << endl;
    }
  } else {
    genparticles_size = 0;
    hardGenEle_size = 0;
    hardGenMu_size = 0;
    genEle_size = 0;
    genMu_size = 0;
    genquarks_size = 0;
    gengluons_size = 0;

    genPair_mass = 0.; 
  }

  L1TInfo(iEvent);
  HLTInfo(iEvent, iSetup);
  //cout << "fine hltinfo" << endl;
  METData(iEvent);
  //cout << "fine met" << endl; 
  //OLDJetData(iEvent);
  //cout << "fine jet" << endl; 
  JetData(iEvent);
  //cout << "fine bjet" << endl; 

  // get the beamspot from the Event:
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByType(theBeamSpot);
  
  // get all beam spot info
  sigmaZ=theBeamSpot->sigmaZ();
  sigmaZ0Error=theBeamSpot->sigmaZ0Error();
  sq=sqrt(sigmaZ*sigmaZ+sigmaZ0Error*sigmaZ0Error);
  
  bsposx = theBeamSpot->position().x();
  bsposy = theBeamSpot->position().y();
  bsposz = theBeamSpot->position().z();

  math::XYZPoint beamspot(theBeamSpot->position().x(),theBeamSpot->position().y(),theBeamSpot->position().z());
  math::XYZPoint firstpvertex(0.,0.,0.);

  //Retrieve primary vertex collection
  Handle<reco::VertexCollection> primaryVertexColl;
  iEvent.getByLabel("offlinePrimaryVertices",primaryVertexColl);
  const reco::VertexCollection* pvcoll = primaryVertexColl.product();

  pvsize = 0;
  //We take only the first primary vertex, i.e. the one with the electrons
  if(pvcoll->size() > 0) {
    reco::VertexCollection::const_iterator firstpv = pvcoll->begin();
    firstpvertex.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }
  
  ////////////////////////////////////////////////////////////////
  //Trigger Matching infos :https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideHLTAnalysis 
  edm::Handle<edm::TriggerResults> trigResults; //our trigger result object
  edm::InputTag trigResultsTag("TriggerResults","","HLT"); //make sure have correct process on MC
  //data process=HLT, MC depends, Spring11 is REDIGI311X
  iEvent.getByLabel(trigResultsTag,trigResults);

  edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT"); //make sure have correct process on MC
  //data process=HLT, MC depends, Spring11 is REDIGI311X
  edm::Handle<trigger::TriggerEvent> trigEvent; 
  //iEvent.getByLabel(trigResultsTag,trigEvent);
  iEvent.getByLabel(trigEventTag,trigEvent);
  std::string filterName(""); 
  trigger::size_type filterIndex; 
 

  //There are two vtx collections (at least) : offlinePrimaryVertices and offlinePrimaryVerticeswithBS
  //Now storing info about the second one
  math::XYZPoint firstpvertexwithBS(0.,0.,0.);
  Handle<reco::VertexCollection> primaryVertexCollwithBS;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS",primaryVertexCollwithBS);
  const reco::VertexCollection* pvcollwithBS = primaryVertexCollwithBS.product();

  if(pvcollwithBS->size() > 0) {
    reco::VertexCollection::const_iterator firstpv = pvcollwithBS->begin();
    firstpvertexwithBS.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }



  pvsize = pvcoll->size();
  
  pvx = new float [pvsize];
  pvy = new float [pvsize];
  pvz = new float [pvsize];
  pv_isValid = new bool [pvsize];
  pv_ndof = new int [pvsize];
  pv_nTracks = new int [pvsize];
  pv_normChi2 = new float [pvsize];
  pv_totTrackSize = new int [pvsize];
  mytree->GetBranch("pvx")->SetAddress(pvx);
  mytree->GetBranch("pvy")->SetAddress(pvy);
  mytree->GetBranch("pvz")->SetAddress(pvz);
  mytree->GetBranch("pv_isValid")->SetAddress(pv_isValid);
  mytree->GetBranch("pv_ndof")->SetAddress(pv_ndof);
  mytree->GetBranch("pv_nTracks")->SetAddress(pv_nTracks);
  mytree->GetBranch("pv_normChi2")->SetAddress(pv_normChi2);
  mytree->GetBranch("pv_totTrackSize")->SetAddress(pv_totTrackSize);

  int indexpv = 0;
  for(reco::VertexCollection::const_iterator pvIt = pvcoll->begin(); pvIt != pvcoll->end(); ++pvIt){
    pvx[indexpv] = pvIt->x();
    pvy[indexpv] = pvIt->y();   
    pvz[indexpv] = pvIt->z();  
    pv_isValid[indexpv] = pvIt->isValid();
    pv_ndof[indexpv] = pvIt->ndof();
    pv_nTracks[indexpv] = pvIt->nTracks();
    pv_normChi2[indexpv] = pvIt->normalizedChi2();
    pv_totTrackSize[indexpv] = pvIt->tracksSize();

    ++indexpv;
  }

  muon_size = muons->size();
  muon_pt = new float [muon_size];
  muon_ptError = new float [muon_size];
  muon_gTrk_pt = new float [muon_size];
  muon_gTrk_ptError = new float [muon_size];
  muon_eta = new float [muon_size];
  muon_etaError = new float [muon_size];
  muon_phi = new float [muon_size];
  muon_phiError = new float [muon_size];
  muon_theta = new float [muon_size];
  muon_thetaError = new float [muon_size]; 
  muon_outerPt = new float [muon_size];
  muon_outerEta = new float [muon_size];
  muon_outerPhi = new float [muon_size];
  muon_outerTheta = new float [muon_size];
  muon_px = new float [muon_size];
  muon_py = new float [muon_size];
  muon_pz = new float [muon_size];
  muon_charge = new int [muon_size];
  muon_nhitspixel = new int [muon_size];
  muon_nhitstrack = new int [muon_size];
  muon_nhitsmuons = new int [muon_size];
  muon_nhitstotal = new int [muon_size];
  muon_nlayerswithhits = new int [muon_size];
  muon_nlosthits = new int [muon_size];
  muon_nSegmentMatch = new int [muon_size];
  muon_isTrackerMuon = new bool [muon_size];
  muon_isPFMuon = new bool [muon_size];
  muon_isPFIsolationValid = new bool [muon_size];
  muon_chi2 = new float [muon_size];
  muon_ndof = new int [muon_size];
  muon_normChi2 = new float [muon_size];
  muon_d0 = new float [muon_size];
  muon_d0Error = new float [muon_size];
  muon_dz_cmsCenter = new float [muon_size];
  muon_dz_beamSpot = new float [muon_size];
  muon_dz_firstPVtx = new float [muon_size];
  muon_dz_firstPVtxwithBS = new float [muon_size];
  muon_dzError = new float [muon_size];
  muon_dxy_cmsCenter = new float [muon_size];
  muon_dxy_beamSpot = new float [muon_size];
  muon_dxy_firstPVtx = new float [muon_size];
  muon_dxy_firstPVtxwithBS = new float [muon_size];
  muon_dxyError = new float [muon_size]; 
  muon_trackIso03 = new float [muon_size]; 
  muon_trackIso05 = new float [muon_size]; 
  muon_trackIso03_ptInVeto = new float [muon_size]; 
  muon_trackIso05_ptInVeto = new float [muon_size]; 
  muon_emIso03 = new float [muon_size]; 
  muon_emIso05 = new float [muon_size]; 
  muon_emIso03_ptInVeto = new float [muon_size]; 
  muon_emIso05_ptInVeto = new float [muon_size]; 
  muon_hadIso03 = new float [muon_size]; 
  muon_hadIso05 = new float [muon_size]; 
  muon_hadIso03_ptInVeto = new float [muon_size]; 
  muon_hadIso05_ptInVeto = new float [muon_size]; 
  muon_innerPosx = new float [muon_size];
  muon_innerPosy = new float [muon_size];
  muon_innerPosz = new float [muon_size];
  muMatch_hltL1Mu3p5EG12L3Filtered22 = new bool [muon_size]; 
  muMatch_hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q= new bool [muon_size]; 
  mytree->GetBranch("muon_pt")->SetAddress(muon_pt);
  mytree->GetBranch("muon_ptError")->SetAddress(muon_ptError);
  mytree->GetBranch("muon_gTrk_pt")->SetAddress(muon_gTrk_pt);
  mytree->GetBranch("muon_gTrk_ptError")->SetAddress(muon_gTrk_ptError);
  mytree->GetBranch("muon_eta")->SetAddress(muon_eta);
  mytree->GetBranch("muon_etaError")->SetAddress(muon_etaError);
  mytree->GetBranch("muon_phi")->SetAddress(muon_phi);
  mytree->GetBranch("muon_phiError")->SetAddress(muon_phiError);
  mytree->GetBranch("muon_theta")->SetAddress(muon_theta);
  mytree->GetBranch("muon_thetaError")->SetAddress(muon_thetaError); 
  mytree->GetBranch("muon_outerPt")->SetAddress(muon_outerPt);
  mytree->GetBranch("muon_outerEta")->SetAddress(muon_outerEta);
  mytree->GetBranch("muon_outerPhi")->SetAddress(muon_outerPhi);
  mytree->GetBranch("muon_outerTheta")->SetAddress(muon_outerTheta);
  mytree->GetBranch("muon_px")->SetAddress(muon_px);
  mytree->GetBranch("muon_py")->SetAddress(muon_py);
  mytree->GetBranch("muon_pz")->SetAddress(muon_pz);
  mytree->GetBranch("muon_charge")->SetAddress(muon_charge);
  mytree->GetBranch("muon_nhitspixel")->SetAddress(muon_nhitspixel);
  mytree->GetBranch("muon_nhitstrack")->SetAddress(muon_nhitstrack);
  mytree->GetBranch("muon_nhitsmuons")->SetAddress(muon_nhitsmuons);
  mytree->GetBranch("muon_nhitstotal")->SetAddress(muon_nhitstotal);
  mytree->GetBranch("muon_nlayerswithhits")->SetAddress(muon_nlayerswithhits);
  mytree->GetBranch("muon_nlosthits")->SetAddress(muon_nlosthits);
  mytree->GetBranch("muon_nSegmentMatch")->SetAddress(muon_nSegmentMatch);
  mytree->GetBranch("muon_isTrackerMuon")->SetAddress(muon_isTrackerMuon);
  mytree->GetBranch("muon_isPFMuon")->SetAddress(muon_isPFMuon);
  mytree->GetBranch("muon_isPFIsolationValid")->SetAddress(muon_isPFIsolationValid);
  mytree->GetBranch("muon_chi2")->SetAddress(muon_chi2);
  mytree->GetBranch("muon_ndof")->SetAddress(muon_ndof);
  mytree->GetBranch("muon_normChi2")->SetAddress(muon_normChi2);
  mytree->GetBranch("muon_d0")->SetAddress(muon_d0);
  mytree->GetBranch("muon_d0Error")->SetAddress(muon_d0Error);
  mytree->GetBranch("muon_dz_cmsCenter")->SetAddress(muon_dz_cmsCenter);
  mytree->GetBranch("muon_dz_beamSpot")->SetAddress(muon_dz_beamSpot);
  mytree->GetBranch("muon_dz_firstPVtx")->SetAddress(muon_dz_firstPVtx);
  mytree->GetBranch("muon_dz_firstPVtxwithBS")->SetAddress(muon_dz_firstPVtxwithBS);
  mytree->GetBranch("muon_dzError")->SetAddress(muon_dzError);
  mytree->GetBranch("muon_dxy_cmsCenter")->SetAddress(muon_dxy_cmsCenter);
  mytree->GetBranch("muon_dxy_beamSpot")->SetAddress(muon_dxy_beamSpot);
  mytree->GetBranch("muon_dxy_firstPVtx")->SetAddress(muon_dxy_firstPVtx);
  mytree->GetBranch("muon_dxy_firstPVtxwithBS")->SetAddress(muon_dxy_firstPVtxwithBS);
  mytree->GetBranch("muon_dxyError")->SetAddress(muon_dxyError); 
  mytree->GetBranch("muon_trackIso03")->SetAddress(muon_trackIso03); 
  mytree->GetBranch("muon_trackIso05")->SetAddress(muon_trackIso05); 
  mytree->GetBranch("muon_trackIso03_ptInVeto")->SetAddress(muon_trackIso03_ptInVeto); 
  mytree->GetBranch("muon_trackIso05_ptInVeto")->SetAddress(muon_trackIso05_ptInVeto); 
  mytree->GetBranch("muon_emIso03")->SetAddress(muon_emIso03); 
  mytree->GetBranch("muon_emIso05")->SetAddress(muon_emIso05); 
  mytree->GetBranch("muon_emIso03_ptInVeto")->SetAddress(muon_emIso03_ptInVeto); 
  mytree->GetBranch("muon_emIso05_ptInVeto")->SetAddress(muon_emIso05_ptInVeto); 
  mytree->GetBranch("muon_hadIso03")->SetAddress(muon_hadIso03); 
  mytree->GetBranch("muon_hadIso05")->SetAddress(muon_hadIso05); 
  mytree->GetBranch("muon_hadIso03_ptInVeto")->SetAddress(muon_hadIso03_ptInVeto); 
  mytree->GetBranch("muon_hadIso05_ptInVeto")->SetAddress(muon_hadIso05_ptInVeto); 
  mytree->GetBranch("muon_innerPosx")->SetAddress(muon_innerPosx);
  mytree->GetBranch("muon_innerPosy")->SetAddress(muon_innerPosy);
  mytree->GetBranch("muon_innerPosz")->SetAddress(muon_innerPosz);
  mytree->GetBranch("muMatch_hltL1Mu3p5EG12L3Filtered22")->SetAddress(muMatch_hltL1Mu3p5EG12L3Filtered22);
  mytree->GetBranch("muMatch_hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q")->SetAddress(muMatch_hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q);
  //cout << "muons" << endl; 
  int index_mu = 0;
  //LOOP OVER MUONS
  for(reco::MuonCollection::const_iterator muIt = muons->begin(); muIt != muons->end(); ++muIt) {
    
    if (muIt->isGlobalMuon()) {
      // get TeV optimized track
      reco::Muon::MuonTrackTypePair tevOptMuTrk = muon::tevOptimized(*muIt, 200, 17., 40., 0.25);

      if (tevOptMuTrk.first->pt() < muPtMin_) continue;
  
      muon_pt[index_mu] = tevOptMuTrk.first->pt();
      muon_ptError[index_mu] = tevOptMuTrk.first->ptError();
      muon_gTrk_pt[index_mu] = muIt->globalTrack()->pt();
      muon_gTrk_ptError[index_mu] = muIt->globalTrack()->ptError();
      muon_eta[index_mu] = tevOptMuTrk.first->eta();
      muon_etaError[index_mu] = tevOptMuTrk.first->etaError();
      muon_phi[index_mu] = tevOptMuTrk.first->phi();
      muon_phiError[index_mu] = tevOptMuTrk.first->phiError();
      muon_theta[index_mu] = tevOptMuTrk.first->theta();
      muon_thetaError[index_mu] = tevOptMuTrk.first->thetaError();      
      muon_outerPt[index_mu] = muIt->globalTrack()->outerPt();
      muon_outerEta[index_mu] = muIt->globalTrack()->outerEta();
      muon_outerPhi[index_mu] = muIt->globalTrack()->outerPhi();
      muon_outerTheta[index_mu] = muIt->globalTrack()->outerTheta();      
      muon_px[index_mu] = tevOptMuTrk.first->px();
      muon_py[index_mu] = tevOptMuTrk.first->py();
      muon_pz[index_mu] = tevOptMuTrk.first->pz();      
      muon_charge[index_mu] = tevOptMuTrk.first->charge();
      muon_nhitspixel[index_mu] = muIt->innerTrack()->hitPattern().numberOfValidPixelHits();
      muon_nhitstrack[index_mu] = muIt->globalTrack()->hitPattern().numberOfValidTrackerHits();
      muon_nhitsmuons[index_mu] = muIt->globalTrack()->hitPattern().numberOfValidMuonHits();
      muon_nhitstotal[index_mu] = muIt->globalTrack()->numberOfValidHits();
      muon_nlayerswithhits[index_mu] = muIt->track()->hitPattern().trackerLayersWithMeasurement();
      muon_nlosthits[index_mu] = muIt->globalTrack()->numberOfLostHits();
      muon_nSegmentMatch[index_mu] = muIt->numberOfMatchedStations();
      muon_isTrackerMuon[index_mu] = muIt->isTrackerMuon();
      muon_isPFMuon[index_mu] = muIt->isPFMuon();
      muon_isPFIsolationValid[index_mu] = muIt->isPFIsolationValid();
      muon_chi2[index_mu] = muIt->globalTrack()->chi2();
      muon_ndof[index_mu] = muIt->globalTrack()->ndof();
      muon_normChi2[index_mu] = muIt->globalTrack()->normalizedChi2();
      muon_d0[index_mu] = tevOptMuTrk.first->d0();
      muon_d0Error[index_mu] = tevOptMuTrk.first->d0Error();
      muon_dz_cmsCenter[index_mu] = tevOptMuTrk.first->dz();
      muon_dz_beamSpot[index_mu] = tevOptMuTrk.first->dz(beamspot);
      muon_dz_firstPVtx[index_mu] = tevOptMuTrk.first->dz(firstpvertex);
      muon_dz_firstPVtxwithBS[index_mu] = tevOptMuTrk.first->dz(firstpvertexwithBS);
      muon_dzError[index_mu] = tevOptMuTrk.first->dzError();
      muon_dxy_cmsCenter[index_mu] = tevOptMuTrk.first->dxy();
      muon_dxy_beamSpot[index_mu] = tevOptMuTrk.first->dxy(beamspot);
      muon_dxy_firstPVtx[index_mu] = tevOptMuTrk.first->dxy(firstpvertex);
      muon_dxy_firstPVtxwithBS[index_mu] = tevOptMuTrk.first->dxy(firstpvertexwithBS);
      muon_dxyError[index_mu] = tevOptMuTrk.first->dxyError();
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

      // HLT_Mu22_Photon22_CaloIdL muon leg
      filterName ="hltL1Mu3p5EG12L3Filtered22"; 
      muMatch_hltL1Mu3p5EG12L3Filtered22[index_mu] = false;
      //it is important to specify the right HLT process for the filter, not doing this is a common bug
      filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 
      if(filterIndex<trigEvent->sizeFilters()){ 
        const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
        const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
        //now loop of the trigger objects passing filter
        for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt != trigKeys.end(); ++keyIt) { 
          const trigger::TriggerObject& obj = trigObjColl[*keyIt];
          if(deltaR(tevOptMuTrk.first->eta(), tevOptMuTrk.first->phi(), obj.eta(), obj.phi())<0.5){
            muMatch_hltL1Mu3p5EG12L3Filtered22[index_mu] = true;
          }
        }
      }//end filter size check

      // HLT_Mu40_eta2p1
      filterName ="hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q"; 
      muMatch_hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q[index_mu] = false;
      //it is important to specify the right HLT process for the filter, not doing this is a common bug
      filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 
      if(filterIndex<trigEvent->sizeFilters()){ 
        const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
        const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
        //now loop of the trigger objects passing filter
        for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt != trigKeys.end(); ++keyIt) { 
          const trigger::TriggerObject& obj = trigObjColl[*keyIt];
          if(deltaR(tevOptMuTrk.first->eta(), tevOptMuTrk.first->phi(), obj.eta(), obj.phi())<0.5){
            muMatch_hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q[index_mu] = true;
          }
        }
      }//end filter size check
 
      index_mu++;
    }
  }
  muon_size = index_mu;
  //  cout << "muon_size " << muon_size << endl;
  //Get a Handle on different collections
  //Get the superclusters
  edm::Handle<reco::SuperClusterCollection> pHybridSuperClusters;
  edm::Handle<reco::SuperClusterCollection> pIslandSuperClusters;
  
  try
    {
      iEvent.getByLabel("correctedHybridSuperClusters","",pHybridSuperClusters);
      iEvent.getByLabel("correctedMulti5x5SuperClustersWithPreshower","",pIslandSuperClusters);
      }
  catch(cms::Exception &ex){}
  const reco::SuperClusterCollection *hybridSuperClusters = pHybridSuperClusters.product();
  const reco::SuperClusterCollection *islandSuperClusters = pIslandSuperClusters.product();
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
  scsize = 0;
  
  //scsize = sclusters.size();
  
  //sort all the refSC by transverse energy
  std::sort(refsclusters.begin(),refsclusters.end(),refScEtGreater);
 
  
  //sort all the SC by transverse energy
  std::sort(sclusters.begin(),sclusters.end(),scEtGreater);
  
  std::vector<const reco::SuperCluster*>::const_iterator sciterforptcut=sclusters.begin();

  for(; sciterforptcut!=sclusters.end(); ++sciterforptcut)
  {
    if( ((*sciterforptcut)->rawEnergy()+(*sciterforptcut)->preshowerEnergy())/cosh((*sciterforptcut)->eta())> ScPtMin_) scsize ++;
  }
  
  //scgsfmatched = new float [scsize];
  //scseedmatched = new float [scsize];
  scenergy = new float [scsize];
  sceta = new float [scsize];
  scetacorr = new float [scsize];
  sctheta = new float [scsize];
  scthetacorr = new float [scsize];
  scet = new float [scsize];
  scphi = new float [scsize];
  scpx = new float [scsize];
  scpy = new float [scsize];
  scpz = new float [scsize];
  scx = new float [scsize];
  scy = new float [scsize];
  scz = new float [scsize];
  scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter = new bool [scsize];
  //mytree->GetBranch("scgsfmatched")->SetAddress(scgsfmatched);
  //mytree->GetBranch("scseedmatched")->SetAddress(scseedmatched);
  mytree->GetBranch("scenergy")->SetAddress(scenergy);
  mytree->GetBranch("sceta")->SetAddress(sceta);
  mytree->GetBranch("scetacorr")->SetAddress(scetacorr);
  //  mytree->GetBranch("sctheta")->SetAddress(sctheta);
  // mytree->GetBranch("scthetacorr")->SetAddress(scthetacorr);
  mytree->GetBranch("scet")->SetAddress(scet);
  mytree->GetBranch("scphi")->SetAddress(scphi);
  mytree->GetBranch("scpx")->SetAddress(scpx);
  mytree->GetBranch("scpy")->SetAddress(scpy);
  mytree->GetBranch("scpz")->SetAddress(scpz);
  mytree->GetBranch("scx")->SetAddress(scx);
  mytree->GetBranch("scy")->SetAddress(scy);
  mytree->GetBranch("scz")->SetAddress(scz);

  int counter = 0;
  std::vector<const reco::SuperCluster*>::const_iterator sciter=sclusters.begin();

  for(; sciter!=sclusters.end(); ++sciter)
  {
    if( ((*sciter)->rawEnergy()+(*sciter)->preshowerEnergy())/cosh((*sciter)->eta())<= ScPtMin_) continue;
      
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
   
    // Trigger matching 
    // HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50 2nd leg 
    //The second leg is a sc, not a gsf (T&P trigger)
    filterName = "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter";
    
    scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter[counter] =false;
    //it is important to specify the right HLT process for the filter, not doing this is a common bug
    filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 
  
    if(filterIndex<trigEvent->sizeFilters()){ 
      const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
      const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
 
      //now loop of the trigger objects passing filter
      for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
	
	const trigger::TriggerObject& obj = trigObjColl[*keyIt];
	//do what you want with the trigger objects, you have
	//eta,phi,pt,mass,p,px,py,pz,et,energy accessors
      
	if(deltaR((*sciter)->eta(),(*sciter)->phi(),obj.eta(), obj.phi())<0.3){
	   scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter[counter] =true;
	}
      }
    }//end filter size check

    counter++;
  }
  
  //trying to see if the sc is seed associated
  
  edm::Handle<GsfTrackCollection> gsfTracksH ;
  iEvent.getByLabel("electronGsfTracks",gsfTracksH) ;
  const GsfTrackCollection *gsftracks = gsfTracksH.product();

  gsftracksize = 0;
  for(GsfTrackCollection::const_iterator gsftrackiterforptcut = gsftracks->begin(); 
      gsftrackiterforptcut!=gsftracks->end();
      ++gsftrackiterforptcut)
  {
    if(gsftrackiterforptcut->pt()<GsfTrackPtMin_) continue;
    gsftracksize++;
  } 
  //cout << "gsftracksize " <<gsftracksize << endl;
  gsftracketa = new float [gsftracksize];
  gsftrackphi = new float [gsftracksize];
  gsftrackp = new float [gsftracksize];
  gsftrackpt = new float [gsftracksize];
  gsftrackpx = new float [gsftracksize];
  gsftrackpy = new float [gsftracksize];
  gsftrackpz = new float [gsftracksize];
  mytree->GetBranch("gsftracketa")->SetAddress(gsftracketa);
  mytree->GetBranch("gsftrackphi")->SetAddress(gsftrackphi);
  mytree->GetBranch("gsftrackp")->SetAddress(gsftrackp);
  mytree->GetBranch("gsftrackpt")->SetAddress(gsftrackpt);
  mytree->GetBranch("gsftrackpx")->SetAddress(gsftrackpx);
  mytree->GetBranch("gsftrackpy")->SetAddress(gsftrackpy);
  mytree->GetBranch("gsftrackpz")->SetAddress(gsftrackpz);

  int v=0;
  for(GsfTrackCollection::const_iterator gsftrackiter = gsftracks->begin(); 
      gsftrackiter!=gsftracks->end();
      ++gsftrackiter)
  {
    if(gsftrackiter->pt()<GsfTrackPtMin_) continue;
    gsftracketa[v] = gsftrackiter->eta();
    gsftrackphi[v] = gsftrackiter->phi();  
    gsftrackp[v] = gsftrackiter->p();
    gsftrackpt[v] = gsftrackiter->pt();
    gsftrackpx[v] = gsftrackiter->px();
    gsftrackpy[v] = gsftrackiter->py();
    gsftrackpz[v] = gsftrackiter->pz();
    v++;
  }//end of loop on gsf tracks

  //To remove spikes (ECAL CLUSTER LAZY TOOLS)
  EcalClusterLazyTools lazytool(iEvent,iSetup,InputTag("reducedEcalRecHitsEB"),InputTag("reducedEcalRecHitsEE"));

  gsf_size = 0;
  gsf0_crystal_size=0; 
  gsf1_crystal_size=0; 
  pfele_size = 0;
  // rechits test Laurent

  ESHandle<CaloGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry* geo=pG.product();
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  
  Handle<EcalRecHitCollection> EBhits;
  //  event.getByLabel(ebhitcoll_,EBhits);
  iEvent.getByLabel("reducedEcalRecHitsEB",EBhits);
  //const EcalRecHitCollection *ebRecHits=EBhits.product();
  Handle<EcalRecHitCollection> EEhits;
  iEvent.getByLabel("reducedEcalRecHitsEE",EEhits);
  // const EcalRecHitCollection *eeRecHits=EEhits.product();
 
  for( reco::GsfElectronCollection::const_iterator gsfiterforptcut = gsfelectrons.begin(); gsfiterforptcut != gsfelectrons.end(); ++gsfiterforptcut) {
    if( gsfiterforptcut->caloEnergy()*sin(gsfiterforptcut->p4().theta()) <GsfPtMin_ ) continue;
    
    
    if(fabs((*gsfiterforptcut).superCluster()->eta())<1.479){//First : Barrel
      for(reco::CaloCluster_iterator bcIt = (*gsfiterforptcut).superCluster()->clustersBegin();
	  bcIt != (*gsfiterforptcut).superCluster()->clustersEnd();
	  ++bcIt) { 
	for(std::vector< std::pair<DetId, float> >::const_iterator rhIt = (*bcIt)->hitsAndFractions().begin();
	    rhIt != (*bcIt)->hitsAndFractions().end(); 
	    ++rhIt) { //loop over rec hits in basic cluster
    	  for(EcalRecHitCollection::const_iterator it = EBhits->begin();
	      it !=  EBhits->end(); 
	      ++it) { //loop over all rec hits to find the right ones
	    if  (rhIt->first ==  (*it).id() ) { //found the matching rechit
	      if(gsf_size==0) {
		gsf0_crystal_size++; 
	      }
	      if(gsf_size==1) {
		gsf1_crystal_size++; 
	      }
	    }
    	  }
	}
      }
    }
    //Now looking at endcaps rechits
    else{
      for(reco::CaloCluster_iterator bcIt = (*gsfiterforptcut).superCluster()->clustersBegin();
	  bcIt != (*gsfiterforptcut).superCluster()->clustersEnd(); 
	  ++bcIt) {
	for(std::vector< std::pair<DetId, float> >::const_iterator rhIt = (*bcIt)->hitsAndFractions().begin();
	    rhIt != (*bcIt)->hitsAndFractions().end(); ++rhIt) { //loop over rec hits in basic cluster
	  for(EcalRecHitCollection::const_iterator it = EEhits->begin();
	      it !=  EEhits->end(); ++it) { //loop over all rec hits to find the right ones
	    if  (rhIt->first ==  (*it).id() ) { //found the matching rechit
	      if(gsf_size==0) {
		gsf0_crystal_size++; 
	      }
	      if(gsf_size==1) {
		gsf1_crystal_size++; 
	      }
	    }
	  }
	}
      }
    }
    gsf_size++;
  }
  
  conv_size = 0;
  
  edm::Handle<reco::ConversionCollection> hConversions;
  iEvent.getByLabel("allConversions", hConversions);
  for (reco::ConversionCollection::const_iterator conv = hConversions->begin(); conv!= hConversions->end(); ++conv) {
    reco::Vertex vtx = conv->conversionVertex();
    if (vtx.isValid()) {
      for(reco::GsfElectronCollection::const_iterator gsfiterforconv = gsfelectrons.begin(); gsfiterforconv!=gsfelectrons.end(); ++gsfiterforconv) {
        if (ConversionTools::matchesConversion(*gsfiterforconv, *conv)) {
          conv_size++;
          break;
        }
      }
    }
  }

  //cout << "gsf_size " <<  gsf_size << endl;


  //Loop over PF electrons to get the size of the vector 
  //See here : https://savannah.cern.ch/task/?32346
  edm::Handle<reco::PFCandidateCollection> pflowelectrons;
  iEvent.getByLabel(edm::InputTag("particleFlow"),pflowelectrons);
  reco::PFCandidateCollection pfelectrons(pflowelectrons->begin(),pflowelectrons->end());
  for( reco::PFCandidateCollection::const_iterator pfeleiter = pfelectrons.begin(); pfeleiter != pfelectrons.end(); ++pfeleiter) {
    if(pfeleiter->pt()>20 && pfeleiter->particleId()==reco::PFCandidate::e) pfele_size++;
    
  }

  
 
  pfele_pt= new float [pfele_size]; 
  pfele_eta= new float [pfele_size];
  pfele_phi= new float [pfele_size];
  pfele_charge= new int [pfele_size];

  gsf_isEB = new bool [gsf_size];
  gsf_isEE = new bool [gsf_size];
  gsf_px = new float [gsf_size];
  gsf_py = new float [gsf_size];
  gsf_pz = new float [gsf_size];
  gsf_pt = new float [gsf_size];
  //gsf_etSC = new float [gsf_size];
  gsf_eta = new float [gsf_size];
  gsf_phi = new float [gsf_size];
  gsf_theta = new float [gsf_size];
  gsf_charge = new int [gsf_size];
  //gsf_deltaEtaATvtx = new float [gsf_size];
  //gsf_deltaPhiATvtx = new float [gsf_size];
  gsf_deltaEtaATcalo = new float [gsf_size];
  gsf_deltaPhiATcalo = new float [gsf_size];
  gsf_sigmaetaeta = new float [gsf_size];
  gsf_sigmaIetaIeta = new float [gsf_size];
  gsf_ecalEnergy = new float [gsf_size];
  gsf_eOVERp = new float [gsf_size];
  //gsf_ptOVERetsc = new float [gsf_size];
  gsf_dxy = new float [gsf_size];
  gsf_dxy_beamSpot = new float [gsf_size];
  gsf_dxy_firstPVtx = new float [gsf_size];
  gsf_dxy_firstPVtxwithBS = new float [gsf_size];
  gsf_dxyError = new float [gsf_size];
  gsf_dz = new float [gsf_size];
  gsf_dz_beamSpot = new float [gsf_size];
  gsf_dz_firstPVtx = new float [gsf_size];
  gsf_dz_firstPVtxwithBS = new float [gsf_size];
  gsf_dzError = new float [gsf_size];
  gsf_vz = new float [gsf_size];
  gsf_nHits = new int [gsf_size];
  gsf_nLostInnerHits = new int [gsf_size];
  gsf_nLostOuterHits = new int [gsf_size];
  gsf_convFlags = new int [gsf_size];
  gsf_convDist = new float [gsf_size];
  gsf_convDcot = new float [gsf_size];
  gsf_convRadius = new float [gsf_size];
  gsf_fBrem = new float [gsf_size];
  //gsf_e1OVERe9 = new float [gsf_size];
  gsf_e1x5 = new float [gsf_size];
  gsf_e2x5 = new float [gsf_size];
  gsf_e5x5 = new float [gsf_size];
  //gsf_eMax = new float [gsf_size];
  gsf_e1x3 = new float [gsf_size];
  //gsf_e3x1 = new float [gsf_size];
  //gsf_e1x5 = new float [gsf_size];
  //gsf_e2x2 = new float [gsf_size];
  //gsf_e3x2 = new float [gsf_size];
  //gsf_e3x3 = new float [gsf_size];
  //gsf_e4x4 = new float [gsf_size];
  //gsf_e5x5 = new float [gsf_size];
  //gsf_e2x5Right = new float [gsf_size];
  //gsf_e2x5Left = new float [gsf_size];  
  //gsf_e2x5Top = new float [gsf_size];  
  //gsf_e2x5Bottom = new float [gsf_size];
  //gsf_e2x5Max = new float [gsf_size];
  //gsf_eLeft = new float [gsf_size];
  //gsf_eRight = new float [gsf_size];
  //gsf_eTop = new float [gsf_size];
  //gsf_eBottom = new float [gsf_size];
  //gsf_e2nd = new float [gsf_size];
  gsf_p = new float [gsf_size];
  gsf_e = new float [gsf_size];
  gsf_deltaeta = new float [gsf_size];
  gsf_deltaphi = new float [gsf_size];
  gsf_hovere = new float [gsf_size];
  gsf_hdepth1overe = new float [gsf_size];
  gsf_hdepth2overe = new float [gsf_size];
  gsf_hovere2012 = new float [gsf_size];
  gsf_hdepth1overe2012 = new float [gsf_size];
  gsf_hdepth2overe2012 = new float [gsf_size];
  gsf_trackiso = new float [gsf_size];
  gsf_ecaliso = new float [gsf_size];
  gsf_hcaliso1 = new float [gsf_size];
  gsf_hcaliso2 = new float [gsf_size];
  gsf_hcaliso12012 = new float [gsf_size];
  gsf_hcaliso22012 = new float [gsf_size];
  gsf_PFisocharged = new float [gsf_size]; 
  gsf_PFisophoton = new float [gsf_size];
  gsf_PFisoneutral = new float [gsf_size];
  gsf_class = new float [gsf_size];
  gsf_isecaldriven = new bool [gsf_size];
  gsf_istrackerdriven = new bool [gsf_size];
  gsfsc_e = new float [gsf_size];
  gsfsc_pt = new float [gsf_size];
  gsfsc_eta = new float [gsf_size];
  gsfsc_phi = new float [gsf_size];
  gsfsc_px = new float [gsf_size];
  gsfsc_py = new float [gsf_size];
  gsfsc_pz = new float [gsf_size];
  gsf_e2x5overe5x5 = new float [gsf_size];
  gsf_e1x5overe5x5 = new float [gsf_size];
  gsf_gsfet = new float [gsf_size];
  scindexforgsf = new int [gsf_size];
  gsfpass_ET = new bool [gsf_size] ;
  gsfpass_PT = new bool [gsf_size] ;
  gsfpass_DETETA = new bool [gsf_size] ;
  gsfpass_CRACK = new bool [gsf_size] ;
  gsfpass_DETAIN = new bool [gsf_size] ;
  gsfpass_DPHIIN = new bool [gsf_size] ;
  gsfpass_HADEM = new bool [gsf_size] ;
  gsfpass_SIGMAIETAIETA  = new bool [gsf_size];
  gsfpass_E2X5OVER5X5 = new bool [gsf_size] ;
  gsfpass_ISOLEMHADDEPTH1 = new bool [gsf_size];
  gsfpass_ISOLHADDEPTH2 = new bool [gsf_size] ;
  gsfpass_ISOLPTTRKS = new bool [gsf_size] ;
  gsfpass_ECALDRIVEN = new bool [gsf_size] ;
  gsfpass_INVALID = new bool [gsf_size];
  gsfpass_NOMISSINGHITS = new bool [gsf_size];
  gsfpass_NOCONVERSION = new bool [gsf_size];
  gsfpass_DXYFIRSTPV = new bool [gsf_size];
  gsfpass_HEEP = new bool [gsf_size];
  gsfpass_ID = new bool [gsf_size];
  gsfpass_ISO = new bool [gsf_size];
  
  gsfmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter = new bool [gsf_size];
  gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter = new bool [gsf_size];
  gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter = new bool [gsf_size];
  gsfmatch_hltL1sL1SingleEG22 = new bool [gsf_size];
  gsfmatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter = new bool [gsf_size];
  gsfmatch_hltEle33CaloIdLPixelMatchFilter = new bool [gsf_size]; 
  gsfmatch_hltEle27WP80TrackIsoFilter= new bool [gsf_size]; 
  gsfmatch_hltMu22Photon22CaloIdLHEFilter= new bool [gsf_size]; 

  //charge information
  scpixcharge = new int [gsf_size];
  ctfcharge = new int [gsf_size];
  gsfcharge = new int [gsf_size];
  gsfctfscpixconsistent = new bool [gsf_size];
  gsfscpixconsistent = new bool [gsf_size];
  gsfctfconsistent = new bool [gsf_size];

  //Crystal info
  gsf0_crystal_ietaorix = new int [gsf0_crystal_size]; 
  gsf0_crystal_iphioriy = new int [gsf0_crystal_size];   
  gsf0_crystal_energy = new float [gsf0_crystal_size]; 
  gsf0_crystal_eta = new float [gsf0_crystal_size]; 
  gsf1_crystal_ietaorix = new int [gsf1_crystal_size]; 
  gsf1_crystal_iphioriy = new int [gsf1_crystal_size];   
  gsf1_crystal_energy = new float [gsf1_crystal_size]; 
  gsf1_crystal_eta = new float [gsf1_crystal_size];

  //conversion information
  conv_vtxProb = new float [conv_size];
  conv_lxy = new float [conv_size];
  conv_nHitsMax = new int [conv_size];
  conv_eleind = new int [conv_size];

  
  mytree->GetBranch("pfele_pt")->SetAddress(pfele_pt);
  mytree->GetBranch("pfele_eta")->SetAddress(pfele_eta);
  mytree->GetBranch("pfele_phi")->SetAddress(pfele_phi);
  mytree->GetBranch("pfele_charge")->SetAddress(pfele_charge);



  mytree->GetBranch("gsf_py")->SetAddress(gsf_py);
  mytree->GetBranch("gsf_isEB")->SetAddress(gsf_isEB);
  mytree->GetBranch("gsf_isEE")->SetAddress(gsf_isEE);
  mytree->GetBranch("gsf_px")->SetAddress(gsf_px);
  mytree->GetBranch("gsf_py")->SetAddress(gsf_py);
  mytree->GetBranch("gsf_pz")->SetAddress(gsf_pz);
  mytree->GetBranch("gsf_pt")->SetAddress(gsf_pt);
  //mytree->GetBranch("gsf_etSC")->SetAddress(gsf_etSC);
  mytree->GetBranch("gsf_eta")->SetAddress(gsf_eta);
  mytree->GetBranch("gsf_phi")->SetAddress(gsf_phi);
  mytree->GetBranch("gsf_theta")->SetAddress(gsf_theta);
  mytree->GetBranch("gsf_charge")->SetAddress(gsf_charge);
  //mytree->GetBranch("gsf_deltaEtaATvtx")->SetAddress(gsf_deltaEtaATvtx);
  //mytree->GetBranch("gsf_deltaPhiATvtx")->SetAddress(gsf_deltaPhiATvtx);
  mytree->GetBranch("gsf_deltaEtaATcalo")->SetAddress(gsf_deltaEtaATcalo);
  mytree->GetBranch("gsf_deltaPhiATcalo")->SetAddress(gsf_deltaPhiATcalo);
  mytree->GetBranch("gsf_sigmaetaeta")->SetAddress(gsf_sigmaetaeta);
  mytree->GetBranch("gsf_sigmaIetaIeta")->SetAddress(gsf_sigmaIetaIeta);
  mytree->GetBranch("gsf_ecalEnergy")->SetAddress(gsf_ecalEnergy);
  mytree->GetBranch("gsf_eOVERp")->SetAddress(gsf_eOVERp);
  //mytree->GetBranch("gsf_ptOVERetsc")->SetAddress(gsf_ptOVERetsc);
  mytree->GetBranch("gsf_dxy")->SetAddress(gsf_dxy);
  mytree->GetBranch("gsf_dxy_beamSpot")->SetAddress(gsf_dxy_beamSpot);
  mytree->GetBranch("gsf_dxy_firstPVtx")->SetAddress(gsf_dxy_firstPVtx);
  mytree->GetBranch("gsf_dxy_firstPVtxwithBS")->SetAddress(gsf_dxy_firstPVtxwithBS);
  mytree->GetBranch("gsf_dxyError")->SetAddress(gsf_dxyError);
  mytree->GetBranch("gsf_dz")->SetAddress(gsf_dz);
  mytree->GetBranch("gsf_dz_beamSpot")->SetAddress(gsf_dz_beamSpot);
  mytree->GetBranch("gsf_dz_firstPVtx")->SetAddress(gsf_dz_firstPVtx);
  mytree->GetBranch("gsf_dz_firstPVtxwithBS")->SetAddress(gsf_dz_firstPVtxwithBS);
  mytree->GetBranch("gsf_dzError")->SetAddress(gsf_dzError);
  mytree->GetBranch("gsf_vz")->SetAddress(gsf_vz);
  mytree->GetBranch("gsf_nHits")->SetAddress(gsf_nHits);
  mytree->GetBranch("gsf_nLostInnerHits")->SetAddress(gsf_nLostInnerHits);
  mytree->GetBranch("gsf_nLostOuterHits")->SetAddress(gsf_nLostOuterHits);
  mytree->GetBranch("gsf_convFlags")->SetAddress(gsf_convFlags);
  mytree->GetBranch("gsf_convDist")->SetAddress(gsf_convDist);
  mytree->GetBranch("gsf_convDcot")->SetAddress(gsf_convDcot);
  mytree->GetBranch("gsf_convRadius")->SetAddress(gsf_convRadius);
  mytree->GetBranch("gsf_fBrem")->SetAddress(gsf_fBrem);
  //mytree->GetBranch("gsf_e1OVERe9")->SetAddress(gsf_e1OVERe9);
  mytree->GetBranch("gsf_e1x5")->SetAddress(gsf_e1x5);
  mytree->GetBranch("gsf_e2x5")->SetAddress(gsf_e2x5);
  mytree->GetBranch("gsf_e5x5")->SetAddress(gsf_e5x5);
  //mytree->GetBranch("gsf_eMax")->SetAddress(gsf_eMax);
  mytree->GetBranch("gsf_e1x3")->SetAddress(gsf_e1x3);
  //mytree->GetBranch("gsf_e3x1")->SetAddress(gsf_e3x1);
  //mytree->GetBranch("gsf_e1x5")->SetAddress(gsf_e1x5);
  //mytree->GetBranch("gsf_e2x2")->SetAddress(gsf_e2x2);
  //mytree->GetBranch("gsf_e3x2")->SetAddress(gsf_e3x2);
  //mytree->GetBranch("gsf_e3x3")->SetAddress(gsf_e3x3);
  //mytree->GetBranch("gsf_e4x4")->SetAddress(gsf_e4x4);
  //mytree->GetBranch("gsf_e5x5")->SetAddress(gsf_e5x5);
  //mytree->GetBranch("gsf_e2x5Right")->SetAddress(gsf_e2x5Right);
  //mytree->GetBranch("gsf_e2x5Left")->SetAddress(gsf_e2x5Left);  
  //mytree->GetBranch("gsf_e2x5Top")->SetAddress(gsf_e2x5Top);  
  //mytree->GetBranch("gsf_e2x5Bottom")->SetAddress(gsf_e2x5Bottom);
  //mytree->GetBranch("gsf_e2x5Max")->SetAddress(gsf_e2x5Max);
  //mytree->GetBranch("gsf_eLeft")->SetAddress(gsf_eLeft);
  //mytree->GetBranch("gsf_eRight")->SetAddress(gsf_eRight);
  //mytree->GetBranch("gsf_eTop")->SetAddress(gsf_eTop);
  //mytree->GetBranch("gsf_eBottom")->SetAddress(gsf_eBottom);
  //mytree->GetBranch("gsf_e2nd")->SetAddress(gsf_e2nd);
  mytree->GetBranch("gsf_p")->SetAddress(gsf_p);
  mytree->GetBranch("gsf_e")->SetAddress(gsf_e);
  mytree->GetBranch("gsf_deltaeta")->SetAddress(gsf_deltaeta);
  mytree->GetBranch("gsf_deltaphi")->SetAddress(gsf_deltaphi);
  mytree->GetBranch("gsf_hovere")->SetAddress(gsf_hovere);
  mytree->GetBranch("gsf_hdepth1overe")->SetAddress(gsf_hdepth1overe);
  mytree->GetBranch("gsf_hdepth2overe")->SetAddress(gsf_hdepth2overe);
  mytree->GetBranch("gsf_hovere2012")->SetAddress(gsf_hovere2012);
  mytree->GetBranch("gsf_hdepth1overe2012")->SetAddress(gsf_hdepth1overe2012);
  mytree->GetBranch("gsf_hdepth2overe2012")->SetAddress(gsf_hdepth2overe2012);
  mytree->GetBranch("gsf_trackiso")->SetAddress(gsf_trackiso);
  mytree->GetBranch("gsf_ecaliso")->SetAddress(gsf_ecaliso);
  mytree->GetBranch("gsf_hcaliso1")->SetAddress(gsf_hcaliso1);
  mytree->GetBranch("gsf_hcaliso2")->SetAddress(gsf_hcaliso2);
  mytree->GetBranch("gsf_hcaliso12012")->SetAddress(gsf_hcaliso12012);
  mytree->GetBranch("gsf_hcaliso22012")->SetAddress(gsf_hcaliso22012);
  mytree->GetBranch("gsf_PFisocharged")->SetAddress(gsf_PFisocharged); 
  mytree->GetBranch("gsf_PFisophoton")->SetAddress(gsf_PFisophoton) ;
  mytree->GetBranch("gsf_PFisoneutral")->SetAddress(gsf_PFisoneutral);
  mytree->GetBranch("gsf_class")->SetAddress(gsf_class);
  mytree->GetBranch("gsf_isecaldriven")->SetAddress(gsf_isecaldriven);
  mytree->GetBranch("gsf_istrackerdriven")->SetAddress(gsf_istrackerdriven);
  mytree->GetBranch("gsfsc_e")->SetAddress(gsfsc_e);
  mytree->GetBranch("gsfsc_pt")->SetAddress(gsfsc_pt);
  mytree->GetBranch("gsfsc_eta")->SetAddress(gsfsc_eta);
  mytree->GetBranch("gsfsc_phi")->SetAddress(gsfsc_phi);
  mytree->GetBranch("gsfsc_px")->SetAddress(gsfsc_px);
  mytree->GetBranch("gsfsc_py")->SetAddress(gsfsc_py);
  mytree->GetBranch("gsfsc_pz")->SetAddress(gsfsc_pz);
  mytree->GetBranch("gsf_e2x5overe5x5")->SetAddress(gsf_e2x5overe5x5);
  mytree->GetBranch("gsf_e1x5overe5x5")->SetAddress(gsf_e1x5overe5x5);
  mytree->GetBranch("gsf_gsfet")->SetAddress(gsf_gsfet);
  mytree->GetBranch("gsf_hitsinfo")->SetAddress(gsf_hitsinfo);
  mytree->GetBranch("scindexforgsf")->SetAddress(scindexforgsf);
  mytree->GetBranch("gsfpass_ET")->SetAddress(gsfpass_ET); 
  mytree->GetBranch("gsfpass_PT")->SetAddress(gsfpass_PT); 
  mytree->GetBranch("gsfpass_DETETA")->SetAddress(gsfpass_DETETA); 
  mytree->GetBranch("gsfpass_CRACK")->SetAddress(gsfpass_CRACK); 
  mytree->GetBranch("gsfpass_DETAIN")->SetAddress(gsfpass_DETAIN); 
  mytree->GetBranch("gsfpass_DPHIIN")->SetAddress(gsfpass_DPHIIN); 
  mytree->GetBranch("gsfpass_HADEM")->SetAddress(gsfpass_HADEM); 
  mytree->GetBranch("gsfpass_SIGMAIETAIETA")->SetAddress(gsfpass_SIGMAIETAIETA);
  mytree->GetBranch("gsfpass_E2X5OVER5X5")->SetAddress(gsfpass_E2X5OVER5X5); 
  mytree->GetBranch("gsfpass_ISOLEMHADDEPTH1")->SetAddress(gsfpass_ISOLEMHADDEPTH1);
  mytree->GetBranch("gsfpass_ISOLHADDEPTH2")->SetAddress(gsfpass_ISOLHADDEPTH2); 
  mytree->GetBranch("gsfpass_ISOLPTTRKS")->SetAddress(gsfpass_ISOLPTTRKS); 
  mytree->GetBranch("gsfpass_ECALDRIVEN")->SetAddress(gsfpass_ECALDRIVEN); 
  mytree->GetBranch("gsfpass_INVALID")->SetAddress(gsfpass_INVALID);
  mytree->GetBranch("gsfpass_NOMISSINGHITS")->SetAddress(gsfpass_NOMISSINGHITS);
  mytree->GetBranch("gsfpass_NOCONVERSION")->SetAddress(gsfpass_NOCONVERSION);
  mytree->GetBranch("gsfpass_DXYFIRSTPV")->SetAddress(gsfpass_DXYFIRSTPV);
  mytree->GetBranch("gsfpass_HEEP")->SetAddress(gsfpass_HEEP);
  mytree->GetBranch("gsfpass_ID")->SetAddress(gsfpass_ID);
  mytree->GetBranch("gsfpass_ISO")->SetAddress(gsfpass_ISO);
  mytree->GetBranch("gsfmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter")->SetAddress(gsfmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter);
  mytree->GetBranch("gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter")->SetAddress(gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter); 
  mytree->GetBranch("gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter")->SetAddress(gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter); 
  mytree->GetBranch("gsfmatch_hltL1sL1SingleEG22")->SetAddress(gsfmatch_hltL1sL1SingleEG22);
  mytree->GetBranch("gsfmatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter")->SetAddress(gsfmatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter); 
  mytree->GetBranch("gsfmatch_hltEle33CaloIdLPixelMatchFilter")->SetAddress(gsfmatch_hltEle33CaloIdLPixelMatchFilter);
  mytree->GetBranch("scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter")->SetAddress(scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter); 
  mytree->GetBranch("gsfmatch_hltEle27WP80TrackIsoFilter")->SetAddress(gsfmatch_hltEle27WP80TrackIsoFilter);
  mytree->GetBranch("gsfmatch_hltMu22Photon22CaloIdLHEFilter")->SetAddress(gsfmatch_hltMu22Photon22CaloIdLHEFilter);

  //charge information
  mytree->GetBranch("scpixcharge")->SetAddress(scpixcharge);
  mytree->GetBranch("ctfcharge")->SetAddress(ctfcharge);
  mytree->GetBranch("gsfcharge")->SetAddress(gsfcharge);
  mytree->GetBranch("gsfctfscpixconsistent")->SetAddress(gsfctfscpixconsistent);
  mytree->GetBranch("gsfscpixconsistent")->SetAddress(gsfscpixconsistent);
  mytree->GetBranch("gsfctfconsistent")->SetAddress(gsfctfconsistent);

  //Crystal info
 
  mytree->GetBranch("gsf0_crystal_ietaorix")->SetAddress(gsf0_crystal_ietaorix);
  mytree->GetBranch("gsf0_crystal_iphioriy")->SetAddress(gsf0_crystal_iphioriy);
  mytree->GetBranch("gsf0_crystal_energy")->SetAddress(gsf0_crystal_energy);
  mytree->GetBranch("gsf0_crystal_eta")->SetAddress(gsf0_crystal_eta);
  mytree->GetBranch("gsf1_crystal_ietaorix")->SetAddress(gsf1_crystal_ietaorix);
  mytree->GetBranch("gsf1_crystal_iphioriy")->SetAddress(gsf1_crystal_iphioriy);
  mytree->GetBranch("gsf1_crystal_energy")->SetAddress(gsf1_crystal_energy);
  mytree->GetBranch("gsf1_crystal_eta")->SetAddress(gsf1_crystal_eta);


  //Conversion information 
  // mytree->GetBranch("conv_size")->SetAddress(conv_size); 
  mytree->GetBranch("conv_vtxProb")->SetAddress(conv_vtxProb); 
  mytree->GetBranch("conv_lxy")->SetAddress(conv_lxy); 
  mytree->GetBranch("conv_nHitsMax")->SetAddress(conv_nHitsMax); 
  mytree->GetBranch("conv_eleind")->SetAddress(conv_eleind);
 
  int e=0;
  int nHeepEle = 0;
  reco::GsfElectronCollection::const_iterator gsfiter = gsfelectrons.begin();
  for(; gsfiter != gsfelectrons.end(); ++gsfiter) {
    if( gsfiter->caloEnergy()*sin(gsfiter->p4().theta()) < GsfPtMin_) continue;
    scindexforgsf[e] = -3;
    //try to get the index for the sc assoc to this gsf
    reco::SuperClusterRef gsfrefsc = gsfiter->superCluster();
    for(unsigned int k=0;k<refsclusters.size();k++){
      if(gsfrefsc == refsclusters[k]){
        scindexforgsf[e] = k;
      }
    }

    //Adding the new H/E (2012) definition. See here : https://twiki.cern.ch/twiki/bin/viewauth/CMS/HoverE2012 
    std::vector<CaloTowerDetId> hcalTowersBehindClusters = hcalHelper->hcalTowersBehindClusters(*(gsfiter->superCluster()));
 
   
    gsf_hdepth1overe2012[e] = hcalHelper->hcalESumDepth1BehindClusters(hcalTowersBehindClusters)/gsfiter->superCluster()->energy();
    gsf_hdepth2overe2012[e] = hcalHelper->hcalESumDepth2BehindClusters(hcalTowersBehindClusters)/gsfiter->superCluster()->energy();
    gsf_hovere2012[e] = gsf_hdepth1overe2012[e] + gsf_hdepth2overe2012[e] ; 
    // The new H/E definition implies to also change the HCalIso definition 
    gsf_hcaliso12012[e] = gsfiter->dr03HcalDepth1TowerSumEt() + ( gsfiter->hcalDepth1OverEcal()  - gsf_hdepth1overe2012[e] )*gsfiter->superCluster()->energy()/cosh(gsfiter->superCluster()->eta()); 
    gsf_hcaliso22012[e] =  gsfiter->dr03HcalDepth2TowerSumEt() + ( gsfiter->hcalDepth2OverEcal()  - gsf_hdepth2overe2012[e] )*gsfiter->superCluster()->energy()/cosh(gsfiter->superCluster()->eta()); 

    // get the iso deposits. 3 (charged hadrons, photons, neutral hadrons)
    unsigned nTypes=3;
    IsoDepositMaps electronIsoDep(nTypes);

    for (size_t j = 0; j<inputTagIsoDepElectrons_.size(); ++j) {
      iEvent.getByLabel(inputTagIsoDepElectrons_[j], electronIsoDep[j]);
    }

    IsoDepositVals electronIsoValPFId(nTypes);

    // No longer needed. e/g recommendation (04/04/12)
    //  IsoDepositVals electronIsoValNoPFId(nTypes);

    for (size_t j = 0; j<inputTagIsoValElectronsPFId_.size(); ++j) {
      iEvent.getByLabel(inputTagIsoValElectronsPFId_[j], electronIsoValPFId[j]);
    }

    const IsoDepositVals * electronIsoVals =  &electronIsoValPFId  ;
    reco::GsfElectronRef myElectronRef(pGsfElectrons ,e);
    gsf_PFisocharged[e] =  (*(*electronIsoVals)[0])[myElectronRef];
    gsf_PFisophoton[e] =(*(*electronIsoVals)[1])[myElectronRef]  ;
    gsf_PFisoneutral[e] = (*(*electronIsoVals)[2])[myElectronRef];

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
    gsf_deltaeta[e] = gsfiter->deltaEtaSuperClusterTrackAtVtx();
    gsf_deltaphi[e] = gsfiter->deltaPhiSuperClusterTrackAtVtx();
    gsf_hovere[e] = gsfiter->hadronicOverEm();
    gsf_hdepth1overe[e] = gsfiter->hcalDepth1OverEcal();
    gsf_hdepth2overe[e] = gsfiter->hcalDepth2OverEcal();
    gsf_trackiso[e] = gsfiter->dr03TkSumPt();
    gsf_ecaliso[e] = gsfiter->dr03EcalRecHitSumEt();
    gsf_hcaliso1[e] = gsfiter->dr03HcalDepth1TowerSumEt();
    gsf_hcaliso2[e] = gsfiter->dr03HcalDepth2TowerSumEt();

    gsf_charge[e] = gsfiter->charge();
    gsf_sigmaetaeta[e] = gsfiter->sigmaEtaEta();
    gsf_sigmaIetaIeta[e] = gsfiter->sigmaIetaIeta();
    if(gsfiter->ecalDrivenSeed())  gsf_isecaldriven[e] = true; 
    else{gsf_isecaldriven[e] = 0;}
    if(gsfiter->trackerDrivenSeed()) gsf_istrackerdriven[e] = true;
    else{gsf_istrackerdriven[e] = 0;}
    gsfsc_e[e] = gsfiter->superCluster()->energy();//gsfiter->superCluster()->rawEnergy()+gsfiter->superCluster()->preshowerEnergy();
    gsfsc_pt[e] = (gsfiter->superCluster()->rawEnergy()+gsfiter->superCluster()->preshowerEnergy())/cosh(gsfiter->superCluster()->eta());
    gsfsc_eta[e] = gsfiter->superCluster()->eta();
    gsfsc_phi[e] = gsfiter->superCluster()->phi();
    gsfsc_px[e] = gsfsc_pt[e]*cos(gsfsc_phi[e]);
    gsfsc_py[e] = gsfsc_pt[e]*sin(gsfsc_phi[e]);
    gsfsc_pz[e] = (gsfiter->superCluster()->rawEnergy()+gsfiter->superCluster()->preshowerEnergy())*tanh(gsfiter->superCluster()->eta());
    gsf_gsfet[e] = gsfiter->caloEnergy()*sin(gsfiter->p4().theta());
    
    //Laurent SuperStar (if it works)
    //Test Laurent track hits info 
    //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/DataFormats/TrackReco/interface/HitPattern.h?revision=1.32&view=markup

    reco::HitPattern kfHitPattern = gsfiter->gsfTrack()->hitPattern();
    nbtrackhits = kfHitPattern.numberOfHits(); 
        
    for(int hititer =0; hititer<25;hititer++){
      if(e>1) continue;
      if(hititer<nbtrackhits){
	unsigned int myhitbin = kfHitPattern.getHitPattern(hititer);
        gsf_hitsinfo[e][hititer]=myhitbin;

        //for (int j=10; j>=0; j--) {
        //  int bit = (myhitbin >> j) & 0x1;
        //}
        //int NValPixelHit = kfHitPattern.numberOfValidPixelHits();
        //int nhits = gsfiter->gsfTrack()->hitPattern().trackerLayersWithMeasurement();
    
      }
      else{
        gsf_hitsinfo[e][hititer]=0;
      }
    }
    //Try to add info about rechit in the SC 
    //strongly inspired from : http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/DaveC/src/printPhoton.cc

    if(fabs((*gsfiter).superCluster()->eta())<1.479){//First : Barrel
      int iebhit = -1, nclust = 0;
      double amplitot = 0.0;
      double clustot = 0.0;
      
      for(reco::CaloCluster_iterator bcIt = (*gsfiter).superCluster()->clustersBegin();
	  bcIt != (*gsfiter).superCluster()->clustersEnd();
	  ++bcIt) { //loop over basic clusters in SC
	// bcIt seems to be a pointer to a pointer !!!!!!!!!!!!!!!
	 
	double clusterEnergy = (*bcIt)->energy();
	clustot += clusterEnergy;
	nclust +=1;
	for(std::vector< std::pair<DetId, float> >::const_iterator rhIt = (*bcIt)->hitsAndFractions().begin();
	    rhIt != (*bcIt)->hitsAndFractions().end(); 
	    ++rhIt) { //loop over rec hits in basic cluster
    
	  for(EcalRecHitCollection::const_iterator it = EBhits->begin();
	      it !=  EBhits->end(); 
	      ++it) { //loop over all rec hits to find the right ones
	    if  (rhIt->first ==  (*it).id() ) { //found the matching rechit
	 
	      iebhit +=1; 
	      EcalRecHit hit = (*it);
	      EBDetId det    = hit.id(); 
	      //const DetId det2     = hit.id(); 
	      float ampli    = hit.energy();

	      amplitot += ampli;
	      //float time     = hit.time()-toffset; NO OFFSET DECLARED
	      //float time     = hit.time();
	      //float toversig =0;
	      //int   ebflag   = hit.recoFlag();
	      //int sm         = det.ism();
	      //float_t chi2   = hit.chi2();
	      //float_t chi2oot= hit.outOfTimeChi2();
	      
	      GlobalPoint poseb=geo->getPosition(hit.detid());
	      float eta_eb=poseb.eta();
	      //float phi_eb=poseb.phi();
	      //float pf=1.0/cosh(eta_eb);
	      //float eteb=ampli*pf;
	      int ieta=det.ieta();      
	      int iphi=det.iphi();
	      
	      if(e==0) {
		gsf0_crystal_energy[iebhit]=ampli;
		gsf0_crystal_ietaorix[iebhit]=ieta;
		gsf0_crystal_iphioriy[iebhit]=iphi;
		gsf0_crystal_eta[iebhit]=eta_eb;
	      }
	      if(e==1) {
		gsf1_crystal_energy[iebhit]=ampli;
		gsf1_crystal_ietaorix[iebhit]=ieta;
		gsf1_crystal_iphioriy[iebhit]=iphi;
		gsf1_crystal_eta[iebhit]=eta_eb;
	      }
	      //cout << "Barrel electron hit " << iebhit << ", ieta=" << ieta << " iphi= " << iphi << " et=" << eteb << " ampli=" << ampli << " eta = "<< eta_eb<<  endl;
	    }
	  }
	}
      }
      //cout <<  "nb of good hits Barrel : " <<iebhit << endl; 
    }

    //Now looking at endcaps rechits
    else{
      int ieehit = -1, nclustee = 0;
      double amplitotee = 0.0, clustotee = 0.0;
      for(reco::CaloCluster_iterator bcIt = (*gsfiter).superCluster()->clustersBegin();
	  bcIt != (*gsfiter).superCluster()->clustersEnd(); ++bcIt) { //loop over basic clusters in SC
	nclustee +=1;
	//CaloCluster cluster = (*bcIt).superCluster();
	//double clusterEnergyee = cluster.energy();
	
	// bcIt seems to be a pointer to a pointer !!!!!!!!!!!!!!!
	double clusterEnergyee = (*bcIt)->energy();
	clustotee += clusterEnergyee;
	
		
	for(std::vector< std::pair<DetId, float> >::const_iterator rhIt = (*bcIt)->hitsAndFractions().begin();
	    rhIt != (*bcIt)->hitsAndFractions().end(); ++rhIt) { //loop over rec hits in basic cluster
	  
	  for(EcalRecHitCollection::const_iterator it = EEhits->begin();
	      it !=  EEhits->end(); ++it) { //loop over all rec hits to find the right ones
	    
	    if  (rhIt->first ==  (*it).id() ) { //found the matching rechit
	      ieehit += 1;
	      EcalRecHit hit = (*it);
	      EEDetId det = hit.id(); 
	      
	      //int dee=0;
	      float ampli = hit.energy();
	     
	      amplitotee += ampli;
	      //float time     = hit.time()-toffset;
	      //float time     = hit.time();
	      //int   ebflag   = hit.recoFlag();
	      
	      GlobalPoint posee=geo->getPosition(hit.detid());
	      float eta_ee=posee.eta();
	      //float phi_ee=posee.phi();
	      //float pf=1.0/cosh(eta_ee);
	      //float etee=ampli*pf;
	      int ix=det.ix();
	      int iy=det.iy();
	      //int side=det.zside();
	      //int iz=0;

      	      if(e==0) {
		gsf0_crystal_energy[ieehit]=ampli;
		gsf0_crystal_ietaorix[ieehit]=ix;
		gsf0_crystal_iphioriy[ieehit]=iy;
		gsf0_crystal_eta[ieehit]=eta_ee;
	      }
	      if(e==1) {
		gsf1_crystal_energy[ieehit]=ampli;
		gsf1_crystal_ietaorix[ieehit]=ix;
		gsf1_crystal_iphioriy[ieehit]=iy;
		gsf1_crystal_eta[ieehit]=eta_ee;
	      }
	      //cout << "Endcaps electron hit " << ieehit << ", ix=" << ix << " iy= " << iy << " et=" << etee << " ampli=" << ampli << " eta = "<< eta_ee<<  endl;
	    }
	  }
	}
      }
      //      cout <<  "nb of good hits Endcaps : " <<ieehit << endl; 
    }

    gsf_theta[e] = gsfiter->theta();
    gsf_isEB[e] = gsfiter->isEB();
    gsf_isEE[e] = gsfiter->isEE();
    gsf_deltaEtaATcalo[e] = gsfiter->deltaEtaSeedClusterTrackAtCalo();
    gsf_deltaPhiATcalo[e] = gsfiter->deltaPhiSeedClusterTrackAtCalo();
    gsf_ecalEnergy[e] = gsfiter->ecalEnergy();
    gsf_eOVERp[e] = gsfiter->eSuperClusterOverP();
    gsf_dxy[e] = gsfiter->gsfTrack()->dxy();
    gsf_dxy_beamSpot[e] = gsfiter->gsfTrack()->dxy(beamspot);
    gsf_dxy_firstPVtx[e] = gsfiter->gsfTrack()->dxy(firstpvertex);
    gsf_dxy_firstPVtxwithBS[e] = gsfiter->gsfTrack()->dxy(firstpvertexwithBS);
    gsf_dxyError[e] = gsfiter->gsfTrack()->dxyError();
    gsf_dz[e] = gsfiter->gsfTrack()->dz(); 
    gsf_dz_beamSpot[e] = gsfiter->gsfTrack()->dz(beamspot); 
    gsf_dz_firstPVtx[e] = gsfiter->gsfTrack()->dz(firstpvertex); 
    gsf_dz_firstPVtxwithBS[e] = gsfiter->gsfTrack()->dz(firstpvertexwithBS); 
    gsf_dzError[e] = gsfiter->gsfTrack()->dzError(); 
    gsf_vz[e] = gsfiter->gsfTrack()->vz();
    gsf_nHits[e] = gsfiter->gsfTrack()->numberOfValidHits();   
    gsf_nLostInnerHits[e] = gsfiter->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();   
    gsf_nLostOuterHits[e] = gsfiter->gsfTrack()->trackerExpectedHitsOuter().numberOfLostHits();
    gsf_convFlags[e] = gsfiter->convFlags();
    gsf_convDist[e] = gsfiter->convDist();
    gsf_convDcot[e] = gsfiter->convDcot(); 
    gsf_convRadius[e] = gsfiter->convRadius();
    gsf_fBrem[e] = gsfiter->fbrem();
    gsf_e1x5[e] =gsfiter->e1x5() ;
    gsf_e2x5[e] =gsfiter->e2x5Max() ;
    gsf_e5x5[e] =gsfiter->e5x5() ;

    const reco::CaloClusterPtr seed = gsfiter->superCluster()->seed();

    gsf_e1x3[e] = lazytool.e1x3(*seed);
//     gsf_e3x1[e] = lazytool.e3x1(*seed);
//     gsf_e2x2[e] = lazytool.e2x2(*seed);
//     gsf_e3x2[e] = lazytool.e3x2(*seed);
//     gsf_e3x3[e] = lazytool.e3x3(*seed);
//     gsf_e4x4[e] = lazytool.e4x4(*seed);
//     gsf_e2x5Right[e] = lazytool.e2x5Right(*seed);
//     gsf_e2x5Left[e] = lazytool.e2x5Left(*seed);
//     gsf_e2x5Top[e] = lazytool.e2x5Top(*seed);
//     gsf_e2x5Bottom[e] = lazytool.e2x5Bottom(*seed);
//     gsf_e2x5Max[e] = lazytool.e2x5Max(*seed);
//     gsf_eLeft[e] = lazytool.eLeft(*seed);
//     gsf_eRight[e] = lazytool.eRight(*seed);
//     gsf_eTop[e] = lazytool.eTop(*seed);
//     gsf_eBottom[e] = lazytool.eBottom(*seed);
//     gsf_e2nd[e] = lazytool.e2nd(*seed);

    // HEEP selection v4.0  - 07/05/2012 Laurent 
    bool gsfetbarrel = gsf_gsfet[e] > 35.;
    bool gsfetendcap = gsf_gsfet[e] > 35.;
    bool barrelsc = fabs(gsfsc_eta[e]) < 1.442;
    bool endcapsc = (fabs(gsfsc_eta[e]) > 1.56) && (fabs(gsfsc_eta[e]) < 2.5);
    bool deltaetabarrel = fabs(gsf_deltaeta[e]) < 0.005;
    bool deltaetaendcap = fabs(gsf_deltaeta[e]) < 0.007;
    bool deltaphibarrel = fabs(gsf_deltaphi[e]) < 0.06;
    bool deltaphiendcap = fabs(gsf_deltaphi[e]) < 0.06;
    bool hoverebarrel  = gsf_hovere[e] < 0.05;
    bool hovereendcap  = gsf_hovere[e] < 0.05;
    bool sigmaIetaIetabarrel  = true;
    bool sigmaIetaIetaendcap  = gsf_sigmaIetaIeta[e] < 0.03;
    bool e2x5overe5x5barrel  = (gsf_e2x5overe5x5[e] > 0.94) || (gsf_e1x5overe5x5[e] > 0.83);
    bool e2x5overe5x5endcap  = true;
    bool ecalisobarrel = (gsf_ecaliso[e]+gsf_hcaliso1[e]) < (2.+0.03*gsf_gsfet[e] + rho*EcalHcal1EffAreaBarrel_);
    bool ecalisoendcap = true;
    if(gsf_gsfet[e] < 50.) {
      ecalisoendcap = (gsf_ecaliso[e]+gsf_hcaliso1[e]) < 2.5+ rho*EcalHcal1EffAreaEndcaps_;
    }
    else {
      ecalisoendcap = (gsf_ecaliso[e]+gsf_hcaliso1[e]) < (2.5+0.03*(gsf_gsfet[e]-50.)+ rho*EcalHcal1EffAreaEndcaps_ );
    }
    bool hcaliso2barrel  = true;
    bool hcaliso2endcap  = true;
    bool trackisobarrel  = gsf_trackiso[e] < 5.;
    bool trackisoendcap  = gsf_trackiso[e] < 5.;
    bool noMissingHits = gsf_nLostInnerHits[e] <= 1;
    bool noConversion = gsf_convFlags[e] != 3; 
    bool dxyfirstpvbarrel = fabs(gsf_dxy_firstPVtx[e]) <0.02;
    bool dxyfirstpvendcaps =fabs(gsf_dxy_firstPVtx[e]) <0.05;

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
    gsfpass_NOCONVERSION[e] = noConversion;
    gsfpass_DXYFIRSTPV[e] = (dxyfirstpvbarrel && barrelsc) || (dxyfirstpvendcaps && endcapsc); 
    gsfpass_ID[e] = ( gsfpass_DETAIN[e] && gsfpass_DPHIIN[e] && gsfpass_HADEM[e] && gsfpass_SIGMAIETAIETA[e] && gsfpass_E2X5OVER5X5[e]);
    gsfpass_ISO[e] = (gsfpass_ISOLEMHADDEPTH1[e] && gsfpass_ISOLHADDEPTH2[e] && gsfpass_ISOLPTTRKS[e]);

    gsfpass_HEEP[e] = gsfpass_ET[e] && gsfpass_DETAIN[e] && gsfpass_DPHIIN[e] && gsfpass_HADEM[e] && gsfpass_SIGMAIETAIETA[e] && gsfpass_E2X5OVER5X5[e] && gsfpass_ISOLEMHADDEPTH1[e] && gsfpass_ISOLHADDEPTH2[e] && gsfpass_ISOLPTTRKS[e] && gsfpass_NOMISSINGHITS[e]&& gsfpass_DXYFIRSTPV[e];
    if (gsfpass_HEEP[e]) ++nHeepEle;

    //charge info
    scpixcharge[e] = gsfiter->scPixCharge();
    if(gsfiter->closestCtfTrackRef().isNonnull()) ctfcharge[e] = gsfiter->closestCtfTrackRef()->charge();
    gsfcharge[e] = gsfiter->gsfTrack()->charge();
    gsfctfscpixconsistent[e] = gsfiter->isGsfCtfScPixChargeConsistent();
    gsfscpixconsistent[e] = gsfiter->isGsfScPixChargeConsistent();
    gsfctfconsistent[e] = gsfiter->isGsfCtfChargeConsistent();

    ////////////////////////////////////////////////////////////////////////////////////
    // Trigger matching
    ////////////////////////////
    //HLT_DoubleEle33_CaloIdL_GsfTrkIdVL , first leg 
    //This trigger is *not* symmetric because the trigger is built from a L1SingleEG22 
    filterName = "hltEle33CaloIdLPixelMatchFilter";
    gsfmatch_hltEle33CaloIdLPixelMatchFilter[e]=false;

    //std::string filterName("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter");
    //it is important to specify the right HLT process for the filter, not doing this is a common bug
    trigger::size_type filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 
    if(filterIndex<trigEvent->sizeFilters()){ 
      const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
      const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
      //now loop of the trigger objects passing filter
      for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
	const trigger::TriggerObject& obj = trigObjColl[*keyIt];
	//do what you want with the trigger objects, you have
	//eta,phi,pt,mass,p,px,py,pz,et,energy accessors

	if(deltaR(gsfiter->eta(),gsfiter->phi(),obj.eta(), obj.phi())<0.5){
	  gsfmatch_hltEle33CaloIdLPixelMatchFilter[e]=true; 
	}
      }
    }//end filter size check
    
    //HLT_DoubleEle33_CaloIdL_GsfTrkIdVL , second leg 
    filterName = "hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter";
    gsfmatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter[e]=false;
    //it is important to specify the right HLT process for the filter, not doing this is a common bug
    filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 
    if(filterIndex<trigEvent->sizeFilters()){ 
      const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
      const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
      //now loop of the trigger objects passing filter
      for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
	const trigger::TriggerObject& obj = trigObjColl[*keyIt];
	//do what you want with the trigger objects, you have
	//eta,phi,pt,mass,p,px,py,pz,et,energy accessors
	
	if(deltaR(gsfiter->eta(),gsfiter->phi(),obj.eta(), obj.phi())<0.5){
	  gsfmatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter[e]=true; 
	}
      }
    }//end filter size check
    
    //HLT_DoubleEle33_CaloIdL_GsfTrkIdVL , L1 
    //Careful that L1 triggers only have discrete eta phi. Need to be extremely loose. 
    //See here: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/SHarper/SHNtupliser/src/SHTrigInfo.cc?revision=1.5&view=markup&pathrev=HEAD
    filterName ="hltL1sL1SingleEG22";
    gsfmatch_hltL1sL1SingleEG22[e]=false; 
    //it is important to specify the right HLT process for the filter, not doing this is a common bug
    filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 
    if(filterIndex<trigEvent->sizeFilters()){ 
      const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
      const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
      //now loop of the trigger objects passing filter
      for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
	const trigger::TriggerObject& obj = trigObjColl[*keyIt];
	//do what you want with the trigger objects, you have
	//eta,phi,pt,mass,p,px,py,pz,et,energy accessors

	  float objeta = obj.eta(); 
	  float objphi = obj.phi();


	  const double barrelEnd=1.4791;
	  // const double endcapEnd=2.65;
	  const double regionEtaSizeEB=0.522;
	  const double regionEtaSizeEE=1.0;
	  const double regionPhiSize=1.044;

	  double etaBinLow  = 0.;
	  double etaBinHigh = 0.;

	 	  
	  if(fabs(objeta) < barrelEnd){
	    etaBinLow = objeta - regionEtaSizeEB/2.;
	    etaBinHigh = etaBinLow + regionEtaSizeEB;
	  }
	  else{
	    etaBinLow = objeta - regionEtaSizeEE/2.;
	    etaBinHigh = etaBinLow + regionEtaSizeEE;
	  }

	  float deltaPhi=reco::deltaPhi(gsfiter->phi(),objphi);
   
    
	  if(gsfiter->eta() < etaBinHigh && gsfiter->eta() > etaBinLow &&   deltaPhi <regionPhiSize/2. )  {
	    gsfmatch_hltL1sL1SingleEG22[e]=true; 
	  }
      }
    }//end filter size check

    // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, first leg
    filterName ="hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter"; 
    gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter[e]=false;
    //it is important to specify the right HLT process for the filter, not doing this is a common bug
    filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 
    if(filterIndex<trigEvent->sizeFilters()){ 
      const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
      const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
      //now loop of the trigger objects passing filter
      for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
	const trigger::TriggerObject& obj = trigObjColl[*keyIt];
	//do what you want with the trigger objects, you have
	//eta,phi,pt,mass,p,px,py,pz,et,energy accessors
	
	if(deltaR(gsfiter->eta(),gsfiter->phi(),obj.eta(), obj.phi())<0.5){
	  gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter[e]=true; 
	}
      }
    }//end filter size check

    // HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, second leg
    filterName ="hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter"; 
    gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter[e] =false;;
    //it is important to specify the right HLT process for the filter, not doing this is a common bug
    filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 
    if(filterIndex<trigEvent->sizeFilters()){ 
      const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
      const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
      //now loop of the trigger objects passing filter
      for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
	const trigger::TriggerObject& obj = trigObjColl[*keyIt];
	//do what you want with the trigger objects, you have
	//eta,phi,pt,mass,p,px,py,pz,et,energy accessors
      
	if(deltaR(gsfiter->eta(),gsfiter->phi(),obj.eta(), obj.phi())<0.5){
	  gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter[e] =true;
	}
      }
    }//end filter size check

    // HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50 first leg 
    filterName ="hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter"; 
    gsfmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter[e] =false;
    //it is important to specify the right HLT process for the filter, not doing this is a common bug
    filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 
    if(filterIndex<trigEvent->sizeFilters()){ 
      const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
      const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
      //now loop of the trigger objects passing filter
      for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
	const trigger::TriggerObject& obj = trigObjColl[*keyIt];
	//do what you want with the trigger objects, you have
	//eta,phi,pt,mass,p,px,py,pz,et,energy accessors
      
	if(deltaR(gsfiter->eta(),gsfiter->phi(),obj.eta(), obj.phi())<0.5){
	  gsfmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter[e] =true;
	}
      }
    }//end filter size check

    // HLT_Ele27WP80
    filterName ="hltEle27WP80TrackIsoFilter"; 
    gsfmatch_hltEle27WP80TrackIsoFilter[e] =false;
    //it is important to specify the right HLT process for the filter, not doing this is a common bug
    filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 
    if(filterIndex<trigEvent->sizeFilters()){ 
      const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
      const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
      //now loop of the trigger objects passing filter
      for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
	const trigger::TriggerObject& obj = trigObjColl[*keyIt];
	//do what you want with the trigger objects, you have
	//eta,phi,pt,mass,p,px,py,pz,et,energy accessors
      
	if(deltaR(gsfiter->eta(),gsfiter->phi(),obj.eta(), obj.phi())<0.5){
	  gsfmatch_hltEle27WP80TrackIsoFilter[e] =true;
	}
      }
    }//end filter size check

    // HLT_Mu22_Photon22_CaloIdL electron leg
    filterName ="hltMu22Photon22CaloIdLHEFilter"; 
    gsfmatch_hltMu22Photon22CaloIdLHEFilter[e] =false;
    //it is important to specify the right HLT process for the filter, not doing this is a common bug
    filterIndex = trigEvent->filterIndex(edm::InputTag(filterName,"",trigEventTag.process())); 
    if(filterIndex<trigEvent->sizeFilters()){ 
      const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
      const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
      //now loop of the trigger objects passing filter
      for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
	const trigger::TriggerObject& obj = trigObjColl[*keyIt];
	//do what you want with the trigger objects, you have
	//eta,phi,pt,mass,p,px,py,pz,et,energy accessors
      
	if(deltaR(gsfiter->eta(),gsfiter->phi(),obj.eta(), obj.phi())<0.5){
	  gsfmatch_hltMu22Photon22CaloIdLHEFilter[e] =true;
	}
      }
    }//end filter size check

    //increment index for gsf
    e++;
  }
  //Conversion info : https://twiki.cern.ch/twiki/bin/viewauth/CMS/ConversionTools
    
  const reco::BeamSpot &bspot = *theBeamSpot.product();
    
  int iconv=-1;
  for (reco::ConversionCollection::const_iterator conv = hConversions->begin(); conv!= hConversions->end(); ++conv) {
    reco::Vertex vtx = conv->conversionVertex();
    if (vtx.isValid()) {
      int iel=-1;

      for(reco::GsfElectronCollection::const_iterator gsfiterforconv = gsfelectrons.begin(); gsfiterforconv!=gsfelectrons.end(); ++gsfiterforconv) {
        iel++;
        //bool passconversionveto = !ConversionTools::hasMatchedConversion(*gsfiterforconv,hConversions,bspot.position());

        if (ConversionTools::matchesConversion(*gsfiterforconv, *conv)) {
          iconv++;
          conv_eleind[iconv] = iel;
          conv_vtxProb[iconv] = (float)TMath::Prob( vtx.chi2(), vtx.ndof() );
          math::XYZVector mom(conv->refittedPairMomentum());
          double dbsx = vtx.x() - bspot.position().x();   
          double dbsy = vtx.y() - bspot.position().y();
          conv_lxy[iconv] = (float)((mom.x()*dbsx + mom.y()*dbsy)/mom.rho());
          conv_nHitsMax[iconv]=0;
          for (std::vector<uint8_t>::const_iterator it = conv->nHitsBeforeVtx().begin(); it!=conv->nHitsBeforeVtx().end(); ++it) {
            if ((*it)>conv_nHitsMax[iconv]) conv_nHitsMax[iconv] = (int)(*it);
          }
          break;
        }
      }
    }
  }
  //End of conversion info


  //Loop over the PF electron 

  int ctpfele =0; 
  for( reco::PFCandidateCollection::const_iterator pfeleiter = pfelectrons.begin(); pfeleiter != pfelectrons.end(); ++pfeleiter) {
    if(pfeleiter->pt()<=20 || pfeleiter->particleId()!=reco::PFCandidate::e ) continue;
    pfele_pt[ctpfele] = pfeleiter->pt(); 
    pfele_eta[ctpfele] = pfeleiter->eta(); 
    pfele_phi[ctpfele] = pfeleiter->phi(); 
    pfele_charge[ctpfele] = pfeleiter->charge(); 
    ctpfele ++;
  }



  // calculate the invariant mass of two heep electrons if there are any
  heepHeepMass = -100.;
  if (nHeepEle > 1) {
    // first find the two highes et heep electrons
    int iHeep1 = -1;
    int iHeep2 = -1;
    float highestHeepEt = 0;
    for (int n = 0; n < gsf_size; ++n) {
      if (highestHeepEt < gsfsc_e[n] * sin(gsf_theta[n]) && gsfpass_HEEP[n]) {
        iHeep1 = n;
        highestHeepEt = gsfsc_e[n] * sin(gsf_theta[n]);
      }
    }
    highestHeepEt = 0;
    for (int m = 0; m < gsf_size; ++m) {
      if (highestHeepEt < gsfsc_e[m] * sin(gsf_theta[m]) && gsfpass_HEEP[m] && m != iHeep1) {
        iHeep2 = m;
        highestHeepEt = gsfsc_e[m] * sin(gsf_theta[m]);
      }
    }
    // then calculate the invariant mass
    heepHeepMass = CalcInvariantMass(iHeep1, iHeep2); 
  }
 

  mytree->Fill();
 

  delete [] L1trigger_bool;
  delete [] HLTriggers;

  delete [] pvx;
  delete [] pvy;
  delete [] pvz;
  delete [] pv_isValid;
  delete [] pv_ndof;
  delete [] pv_nTracks;
  delete [] pv_normChi2;
  delete [] pv_totTrackSize;

  //delete [] scgsfmatched;
  //delete [] scseedmatched;
  delete [] scenergy;
  delete [] sceta;
  delete [] scetacorr;
  delete [] sctheta;
  delete [] scthetacorr;
  delete [] scet;
  delete [] scphi;
  delete [] scpx;
  delete [] scpy;
  delete [] scpz;
  delete [] scx;
  delete [] scy;
  delete [] scz;


  delete [] pfele_pt;
  delete [] pfele_eta;
  delete [] pfele_phi;
  delete [] pfele_charge;



  delete [] gsf_isEB;
  delete [] gsf_isEE;
  delete [] gsf_px;
  delete [] gsf_py;
  delete [] gsf_pz;
  delete [] gsf_pt;
  //delete [] gsf_etSC;
  delete [] gsf_eta;
  delete [] gsf_phi;
  delete [] gsf_theta;
  delete [] gsf_charge;
  //delete [] gsf_deltaEtaATvtx;
  //delete [] gsf_deltaPhiATvtx;
  delete [] gsf_deltaEtaATcalo;
  delete [] gsf_deltaPhiATcalo;
  delete [] gsf_sigmaetaeta;
  delete [] gsf_sigmaIetaIeta;
  delete [] gsf_ecalEnergy;
  delete [] gsf_eOVERp;
  //delete [] gsf_ptOVERetsc;
  delete [] gsf_dxy;
  delete [] gsf_dxy_beamSpot;
  delete [] gsf_dxy_firstPVtx;
  delete [] gsf_dxy_firstPVtxwithBS;
  delete [] gsf_dxyError;
  delete [] gsf_dz;
  delete [] gsf_dz_beamSpot;
  delete [] gsf_dz_firstPVtx;
  delete [] gsf_dz_firstPVtxwithBS;
  delete [] gsf_dzError;
  delete [] gsf_vz;
  delete [] gsf_nHits;
  delete [] gsf_nLostInnerHits;
  delete [] gsf_nLostOuterHits;
  delete [] gsf_convFlags;
  delete [] gsf_convDist;
  delete [] gsf_convDcot;
  delete [] gsf_convRadius;
  delete [] gsf_fBrem;
  //delete [] gsf_e1OVERe9;
  delete [] gsf_e1x5;
  delete [] gsf_e2x5;
  delete [] gsf_e5x5;
  //delete [] gsf_eMax;
  //delete [] gsf_e1x3;
  //delete [] gsf_e3x1;
  //delete [] gsf_e1x5;
  //delete [] gsf_e2x2;
  //delete [] gsf_e3x2;
  //delete [] gsf_e3x3;
  //delete [] gsf_e4x4;
  //delete [] gsf_e5x5;
  //delete [] gsf_e2x5Right;
  //delete [] gsf_e2x5Left;  
  //delete [] gsf_e2x5Top;  
  //delete [] gsf_e2x5Bottom;
  //delete [] gsf_e2x5Max;
  //delete [] gsf_eLeft;
  //delete [] gsf_eRight;
  //delete [] gsf_eTop;
  //delete [] gsf_eBottom;
  //delete [] gsf_e2nd;
  delete [] gsf_p;
  delete [] gsf_e;
  delete [] gsf_deltaeta;
  delete [] gsf_deltaphi;
  delete [] gsf_hovere;
  delete [] gsf_hdepth1overe;
  delete [] gsf_hdepth2overe;
  delete [] gsf_hovere2012;
  delete [] gsf_hdepth1overe2012;
  delete [] gsf_hdepth2overe2012;
  delete [] gsf_trackiso;
  delete [] gsf_ecaliso;
  delete [] gsf_hcaliso1;
  delete [] gsf_hcaliso2;
  delete [] gsf_hcaliso12012;
  delete [] gsf_hcaliso22012;
  delete [] gsf_PFisocharged; 
  delete [] gsf_PFisophoton;
  delete [] gsf_PFisoneutral;	  
  delete [] gsf_class;
  delete [] gsf_isecaldriven;
  delete [] gsf_istrackerdriven;
  delete [] gsfsc_e;
  delete [] gsfsc_pt;
  delete [] gsfsc_eta;
  delete [] gsfsc_phi;
  delete [] gsfsc_px;
  delete [] gsfsc_py;
  delete [] gsfsc_pz;
  delete [] gsf_e2x5overe5x5;
  delete [] gsf_e1x5overe5x5;
  delete [] gsf_gsfet;
  delete [] scindexforgsf;
  delete [] gsfpass_ET; 
  delete [] gsfpass_PT; 
  delete [] gsfpass_DETETA; 
  delete [] gsfpass_CRACK; 
  delete [] gsfpass_DETAIN; 
  delete [] gsfpass_DPHIIN; 
  delete [] gsfpass_HADEM; 
  delete [] gsfpass_SIGMAIETAIETA;
  delete [] gsfpass_E2X5OVER5X5; 
  delete [] gsfpass_ISOLEMHADDEPTH1;
  delete [] gsfpass_ISOLHADDEPTH2; 
  delete [] gsfpass_ISOLPTTRKS; 
  delete [] gsfpass_ECALDRIVEN; 
  delete [] gsfpass_INVALID;
  delete [] gsfpass_NOMISSINGHITS;
  delete [] gsfpass_NOCONVERSION;
  delete [] gsfpass_DXYFIRSTPV;
  delete [] gsfpass_HEEP;
  delete [] gsfpass_ID;
  delete [] gsfpass_ISO;
  delete [] gsfmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter;
  delete [] gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter;
  delete [] gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter;
  delete [] gsfmatch_hltL1sL1SingleEG22;
  delete [] gsfmatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter;
  delete [] gsfmatch_hltEle33CaloIdLPixelMatchFilter; 
  delete [] scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter;
  delete [] gsfmatch_hltEle27WP80TrackIsoFilter; 
  delete [] gsfmatch_hltMu22Photon22CaloIdLHEFilter; 

  delete [] scpixcharge;
  delete [] ctfcharge;
  delete [] gsfcharge;
  delete [] gsfctfscpixconsistent;
  delete [] gsfscpixconsistent;
  delete [] gsfctfconsistent;
 
  delete [] gsftracketa;
  delete [] gsftrackphi;
  delete [] gsftrackp;
  delete [] gsftrackpt;
  delete [] gsftrackpx;
  delete [] gsftrackpy;
  delete [] gsftrackpz;

  delete [] gsf0_crystal_ietaorix;
  delete [] gsf0_crystal_iphioriy;
  delete [] gsf0_crystal_energy;
  delete [] gsf0_crystal_eta ;
  delete [] gsf1_crystal_ietaorix;
  delete [] gsf1_crystal_iphioriy;
  delete [] gsf1_crystal_energy;
  delete [] gsf1_crystal_eta ;

  delete [] conv_vtxProb;
  delete [] conv_lxy;
  delete [] conv_nHitsMax;
  delete [] conv_eleind;
  
  delete [] muon_pt;
  delete [] muon_ptError;
  delete [] muon_gTrk_pt;
  delete [] muon_gTrk_ptError;
  delete [] muon_eta;
  delete [] muon_etaError;
  delete [] muon_phi;
  delete [] muon_phiError;
  delete [] muon_theta;
  delete [] muon_thetaError; 
  delete [] muon_outerPt;
  delete [] muon_outerEta;
  delete [] muon_outerPhi;
  delete [] muon_outerTheta;
  delete [] muon_px;
  delete [] muon_py;
  delete [] muon_pz;
  delete [] muon_charge;
  delete [] muon_nhitspixel;
  delete [] muon_nhitstrack;
  delete [] muon_nhitsmuons;
  delete [] muon_nhitstotal;
  delete [] muon_nlayerswithhits;
  delete [] muon_nlosthits;
  delete [] muon_nSegmentMatch;
  delete [] muon_isTrackerMuon;
  delete [] muon_isPFMuon;
  delete [] muon_isPFIsolationValid;
  delete [] muon_chi2;
  delete [] muon_ndof;
  delete [] muon_normChi2;
  delete [] muon_d0;
  delete [] muon_d0Error;
  delete [] muon_dz_cmsCenter;
  delete [] muon_dz_beamSpot;
  delete [] muon_dz_firstPVtx;
  delete [] muon_dz_firstPVtxwithBS;
  delete [] muon_dzError;
  delete [] muon_dxy_cmsCenter;
  delete [] muon_dxy_beamSpot;
  delete [] muon_dxy_firstPVtx;
  delete [] muon_dxy_firstPVtxwithBS;
  delete [] muon_dxyError; 
  delete [] muon_trackIso03; 
  delete [] muon_trackIso05; 
  delete [] muon_trackIso03_ptInVeto; 
  delete [] muon_trackIso05_ptInVeto; 
  delete [] muon_emIso03; 
  delete [] muon_emIso05; 
  delete [] muon_emIso03_ptInVeto; 
  delete [] muon_emIso05_ptInVeto; 
  delete [] muon_hadIso03; 
  delete [] muon_hadIso05; 
  delete [] muon_hadIso03_ptInVeto; 
  delete [] muon_hadIso05_ptInVeto; 
  delete [] muon_innerPosx;
  delete [] muon_innerPosy;
  delete [] muon_innerPosz;
  delete [] muMatch_hltL1Mu3p5EG12L3Filtered22; 
  delete [] muMatch_hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q; 

//   delete [] jetAKT_eta;
//   delete [] jetAKT_pt;
//   delete [] jetAKT_phi;
//   delete [] jetAKT_em;

  delete [] pfJet_pt;
  delete [] pfJet_eta;
  delete [] pfJet_phi;
//  delete [] Jet_em;
  delete [] Jet_pt;
  delete [] Jet_eta;
  delete [] Jet_phi;
  delete [] Jet_vx;
  delete [] Jet_vy;
  delete [] Jet_vz;
  delete [] tCHighEffBTags;
  delete [] tCHighPurBTags;
  delete [] jetProbBTags;
  delete [] jetBProbBTags;
  delete [] sSecVertHighEffBTags;
  delete [] sSecVertHighPurBTags;
  delete [] cSecVertBTags;
  delete [] cSecVertMVABTags;
  delete [] ghostTrkBTags;
  delete [] softEleIP3dBTags;
  delete [] softElePtBTags;
  delete [] softMuBTags;
  delete [] softMuIP3dBTags;
  delete [] softMuPtBTags;

  if (useGenData_) {
    delete [] genquark_e;
    delete [] genquark_pt;
    delete [] genquark_px; 
    delete [] genquark_py; 
    delete [] genquark_pz; 
    delete [] genquark_eta; 
    delete [] genquark_phi;
    delete [] genquark_status;
    delete [] genquark_charge;
    delete [] genquark_pdgid;

    delete [] gengluon_e;
    delete [] gengluon_pt;
    delete [] gengluon_px; 
    delete [] gengluon_py; 
    delete [] gengluon_pz; 
    delete [] gengluon_eta; 
    delete [] gengluon_phi;
    delete [] gengluon_status;
    delete [] gengluon_charge;
    delete [] gengluon_pdgid;

    delete [] genele_e;
    delete [] genele_pt;
    delete [] genele_px;
    delete [] genele_py;
    delete [] genele_pz;
    delete [] genele_eta;
    delete [] genele_phi;
    delete [] genele_charge;
    delete [] unstableGenEle_e;
    delete [] unstableGenEle_pt;
    delete [] unstableGenEle_px;
    delete [] unstableGenEle_py;
    delete [] unstableGenEle_pz;
    delete [] unstableGenEle_eta;
    delete [] unstableGenEle_phi;
    delete [] unstableGenEle_charge;
    delete [] hardGenEle_e;
    delete [] hardGenEle_pt;
    delete [] hardGenEle_px;
    delete [] hardGenEle_py;
    delete [] hardGenEle_pz;
    delete [] hardGenEle_eta;
    delete [] hardGenEle_phi;
    delete [] hardGenEle_charge;
    delete [] genelemom_e;
    delete [] genelemom_pt;
    delete [] genelemom_px;
    delete [] genelemom_py;
    delete [] genelemom_pz;
    delete [] genelemom_eta;
    delete [] genelemom_phi;
    delete [] genelemom_charge;
    delete [] genelemom_mass;
    delete [] genelemom_pdgid;

    delete [] genmu_e;
    delete [] genmu_pt;
    delete [] genmu_px;
    delete [] genmu_py;
    delete [] genmu_pz;
    delete [] genmu_eta;
    delete [] genmu_phi;
    delete [] genmu_charge;
    delete [] unstableGenMu_e;
    delete [] unstableGenMu_pt;
    delete [] unstableGenMu_px;
    delete [] unstableGenMu_py;
    delete [] unstableGenMu_pz;
    delete [] unstableGenMu_eta;
    delete [] unstableGenMu_phi;
    delete [] unstableGenMu_charge;
    delete [] hardGenMu_e;
    delete [] hardGenMu_pt;
    delete [] hardGenMu_px;
    delete [] hardGenMu_py;
    delete [] hardGenMu_pz;
    delete [] hardGenMu_eta;
    delete [] hardGenMu_phi;
    delete [] hardGenMu_charge;
    delete [] genmumom_e;
    delete [] genmumom_pt;
    delete [] genmumom_px;
    delete [] genmumom_py;
    delete [] genmumom_pz;
    delete [] genmumom_eta;
    delete [] genmumom_phi;
    delete [] genmumom_charge;
    delete [] genmumom_mass;
    delete [] genmumom_pdgid;

    delete [] x1quark;
    delete [] x2quark;


  }
 }//end of analyze method


void 
GsfCheckerTree::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  edm::LogVerbatim("beginRunHlt")<<"dans beginrun run number"<<iRun.id()<<" and n = "<<hlNames_.size();

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
      edm::LogVerbatim("beginRunHlt") << "HLT menu: " << hltConfig_.tableName();
      for (unsigned int i = 0; i < n; ++i) {
 	hlWasRunTab[i]=0;
	hlAcceptTab[i]=0;
	hlErrorTab[i]=0;
        edm::LogVerbatim("beginRunHlt") << "hlNames(" << i << ") = " << hlNames_.at(i);
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

  //GENERAL TECHNICAL INFOS 
  mytree->Branch("runnumber",&runnumber,"runnumber/i");
  mytree->Branch("eventnumber",&eventnumber,"eventnumber/i");
  mytree->Branch("luminosityBlock",&luminosityBlock,"luminosityBlock/i"); //add Laurent
  mytree->Branch("eventcounter",&eventcounter,"eventcounter/i");
  mytree->Branch("processid",&processid,"processid/I");
  mytree->Branch("pthat",&pthat,"pthat/F");
  mytree->Branch("alphaqcd",&alphaqcd,"alphaqcd/F");
  mytree->Branch("alphaqed",&alphaqed,"alphaqed/F");
  mytree->Branch("qscale",&qscale,"qscale/F");
  mytree->Branch("weight",&weight,"weight/F");

  //TRIGGERS
  mytree->Branch("hltCount",&hltCount,"hltCount/I");
  //mytree->Branch("L1trigger_size", &L1trigger_size, "L1trigger_size/I"); 
  //mytree->Branch("L1trigger_bool", L1trigger_bool, "L1trigger_bool[L1trigger_size]/I");
  mytree->Branch("PhysDecl_bool", &PhysDecl_bool, "PhysDecl_bool/I");

  mytree->Branch("nWasRun_",&nWasRun_,"nWasRun_/I");
  mytree->Branch("nAccept_",&nAccept_,"nAccept_/I");
  mytree->Branch("nErrors_",&nErrors_,"nErrors_/I");
  //mytree->Branch("hlWasRun_",&hlWasRun_,"hlWasRun_/I");
  //mytree->Branch("hlWasRunTab",hlWasRunTab,"hlWasRunTab[400]/I");
  //mytree->Branch("hlAccept_",&hlAccept_);
  //mytree->Branch("hlAcceptTab",hlAcceptTab,"hlAcceptTab[400]/I");
  //mytree->Branch("hlErrorTab",hlErrorTab,"hlErrorTab[200]/I");
  //mytree->Branch("hlNamesTab",&hlNamesTab,"hlNamesTab/C");
  // mytree->Branch("hlNames_",&hlNames_);
  //mytree->Branch("HLTriggers", HLTriggers, "HLTriggers[hltCount]/I");
  mytree->Branch("HLT_Mu15_eta2p1",&HLT_Mu15_eta2p1,"HLT_Mu15_eta2p1/I");
  mytree->Branch("HLT_Mu24_eta2p1",&HLT_Mu24_eta2p1,"HLT_Mu24_eta2p1/I");
  mytree->Branch("HLT_Mu30_eta2p1",&HLT_Mu30_eta2p1,"HLT_Mu30_eta2p1/I");
  mytree->Branch("HLT_Mu40_eta2p1",&HLT_Mu40_eta2p1,"HLT_Mu40_eta2p1/I");
  mytree->Branch("HLT_Mu50_eta2p1",&HLT_Mu50_eta2p1,"HLT_Mu50_eta2p1/I");
  mytree->Branch("HLT_Mu22_TkMu22",&HLT_Mu22_TkMu22,"HLT_Mu22_TkMu22/I");
  mytree->Branch("HLT_Mu22_Photon22_CaloIdL",&HLT_Mu22_Photon22_CaloIdL,"HLT_Mu22_Photon22_CaloIdL/I");
  mytree->Branch("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",&HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL/I");
  mytree->Branch("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",&HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL/I");
  mytree->Branch("HLT_Ele8_CaloIdL_CaloIsoVL", &HLT_Ele8_CaloIdL_CaloIsoVL, "HLT_Ele8_CaloIdL_CaloIsoVL/I");
  mytree->Branch("HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL", &HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL, "HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL/I");
  mytree->Branch("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL", &HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL/I");
  mytree->Branch("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50", &HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50, "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50/I");
  mytree->Branch("HLT_DoubleEle33_CaloIdL",&HLT_DoubleEle33_CaloIdL,"HLT_DoubleEle33_CaloIdL/I");
  mytree->Branch("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL",&HLT_DoubleEle33_CaloIdL_GsfTrkIdVL,"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL/I");
  mytree->Branch("HLT_DoubleEle33_CaloIdT",&HLT_DoubleEle33_CaloIdT,"HLT_DoubleEle33_CaloIdT/I");
  mytree->Branch("HLT_Photon20_CaloIdVL_IsoL", &HLT_Photon20_CaloIdVL_IsoL, "HLT_Photon20_CaloIdVL_IsoL/I");
  mytree->Branch("HLT_Photon30_CaloIdVL", &HLT_Photon30_CaloIdVL, "HLT_Photon30_CaloIdVL/I");
  mytree->Branch("HLT_Photon50_CaloIdVL", &HLT_Photon50_CaloIdVL, "HLT_Photon50_CaloIdVL/I");
  mytree->Branch("HLT_Photon50_CaloIdVL_IsoL", &HLT_Photon50_CaloIdVL_IsoL, "HLT_Photon50_CaloIdVL_IsoL/I");
  mytree->Branch("HLT_Photon75_CaloIdVL", &HLT_Photon75_CaloIdVL, "HLT_Photon75_CaloIdVL/I");
  mytree->Branch("HLT_Photon90_CaloIdVL", &HLT_Photon90_CaloIdVL, "HLT_Photon90_CaloIdVL/I");
  mytree->Branch("HLT_Photon135",&HLT_Photon135,"HLT_Photon135/I");
  mytree->Branch("HLT_Photon150",&HLT_Photon150,"HLT_Photon150/I");
  mytree->Branch("HLT_Photon250_NoHE",&HLT_Photon250_NoHE,"HLT_Photon250_NoHE/I");
  mytree->Branch("HLT_Photon300_NoHE",&HLT_Photon300_NoHE,"HLT_Photon300_NoHE/I");
  mytree->Branch("HLT_Photon26_Photon18",&HLT_Photon26_Photon18,"HLT_Photon26_Photon18/I");
  mytree->Branch("HLT_Photon36_Photon22",&HLT_Photon36_Photon22,"HLT_Photon36_Photon22/I");
  mytree->Branch("HLT_DoublePhoton70",&HLT_DoublePhoton70,"HLT_DoublePhoton70/I");
  mytree->Branch("HLT_DoublePhoton80",&HLT_DoublePhoton80,"HLT_DoublePhoton80/I");
  mytree->Branch("HLT_Ele27_WP80",&HLT_Ele27_WP80,"HLT_Ele27_WP80/I");
  mytree->Branch("prescale_HLT_Mu15_eta2p1",&prescale_HLT_Mu15_eta2p1,"prescale_HLT_Mu15_eta2p1/I");
  mytree->Branch("prescale_HLT_Mu30_eta2p1",&prescale_HLT_Mu30_eta2p1,"prescale_HLT_Mu30_eta2p1/I");
  mytree->Branch("prescale_HLT_Mu40_eta2p1",&prescale_HLT_Mu40_eta2p1,"prescale_HLT_Mu40_eta2p1/I");
  mytree->Branch("prescale_HLT_Mu22_TkMu22",&prescale_HLT_Mu22_TkMu22,"prescale_HLT_Mu22_TkMu22/I");
  mytree->Branch("prescale_HLT_Mu22_Photon22_CaloIdL",&prescale_HLT_Mu22_Photon22_CaloIdL,"prescale_HLT_Mu22_Photon22_CaloIdL/I");
  mytree->Branch("prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",&prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL/I");
  mytree->Branch("prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",&prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL/I");
  mytree->Branch("prescale_HLT_Ele8_CaloIdL_CaloIsoVL",&prescale_HLT_Ele8_CaloIdL_CaloIsoVL,"prescale_HLT_Ele8_CaloIdL_CaloIsoVL/I");
  mytree->Branch("prescale_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL",&prescale_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL,"prescale_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL/I");
  mytree->Branch("prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",&prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL/I");
  mytree->Branch("prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50",&prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50,"prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50/I");
  mytree->Branch("prescale_HLT_DoubleEle33_CaloIdL",&prescale_HLT_DoubleEle33_CaloIdL,"prescale_HLT_DoubleEle33_CaloIdL/I");
  mytree->Branch("prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL",&prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL,"prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL/I");
  mytree->Branch("prescale_HLT_DoubleEle33_CaloIdT",&prescale_HLT_DoubleEle33_CaloIdT,"prescale_HLT_DoubleEle33_CaloIdT/I");
  mytree->Branch("prescale_HLT_Photon20_CaloIdVL_IsoL",&prescale_HLT_Photon20_CaloIdVL_IsoL,"prescale_HLT_Photon20_CaloIdVL_IsoL/I");
  mytree->Branch("prescale_HLT_Photon30_CaloIdVL",&prescale_HLT_Photon30_CaloIdVL,"prescale_HLT_Photon30_CaloIdVL/I");
  mytree->Branch("prescale_HLT_Photon50_CaloIdVL",&prescale_HLT_Photon50_CaloIdVL,"prescale_HLT_Photon50_CaloIdVL/I");
  mytree->Branch("prescale_HLT_Photon50_CaloIdVL_IsoL",&prescale_HLT_Photon50_CaloIdVL_IsoL,"prescale_HLT_Photon50_CaloIdVL_IsoL/I");
  mytree->Branch("prescale_HLT_Photon75_CaloIdVL",&prescale_HLT_Photon75_CaloIdVL,"prescale_HLT_Photon75_CaloIdVL/I");
  mytree->Branch("prescale_HLT_Photon90_CaloIdVL",&prescale_HLT_Photon90_CaloIdVL,"prescale_HLT_Photon90_CaloIdVL/I");
  mytree->Branch("prescale_HLT_Photon135",&prescale_HLT_Photon135,"prescale_HLT_Photon135/I");
  mytree->Branch("prescale_HLT_Photon150",&prescale_HLT_Photon150,"prescale_HLT_Photon150/I");
  mytree->Branch("prescale_HLT_Photon250_NoHE",&prescale_HLT_Photon250_NoHE,"prescale_HLT_Photon250_NoHE/I");
  mytree->Branch("prescale_HLT_Photon300_NoHE",&prescale_HLT_Photon300_NoHE,"prescale_HLT_Photon300_NoHE/I");
  mytree->Branch("prescale_HLT_Photon26_Photon18",&prescale_HLT_Photon26_Photon18,"prescale_HLT_Photon26_Photon18/I");
  mytree->Branch("prescale_HLT_Photon36_Photon22",&prescale_HLT_Photon36_Photon22,"prescale_HLT_Photon36_Photon22/I");
  mytree->Branch("prescale_HLT_DoublePhoton70",&prescale_HLT_DoublePhoton70,"prescale_HLT_DoublePhoton70/I");
  mytree->Branch("prescale_HLT_DoublePhoton80",&prescale_HLT_DoublePhoton80,"prescale_HLT_DoublePhoton80/I");
  mytree->Branch("prescale_HLT_Ele27_WP80",&prescale_HLT_Ele27_WP80,"prescale_HLT_Ele27_WP80/I");

  //GLOBAL PHYSICAL INFO 
  mytree->Branch("rho", &rho, "rho/F");
  //  mytree->Branch("rhoiso", &rhoiso, "rhoiso/F");
  mytree->Branch("calomet", &calomet, "calomet/F");
  mytree->Branch("calomet_phi", &calomet_phi, "calomet_phi/F");
  mytree->Branch("met", &met, "met/F");
  mytree->Branch("pfmet", &pfmet, "pfmet/F");
  mytree->Branch("pfmet_phi", &pfmet_phi, "pfmet_phi/F");
  mytree->Branch("pfmetcor", &pfmetcor, "pfmetcor/F");
  mytree->Branch("pfmetcor_phi", &pfmetcor_phi, "pfmetcor_phi/F");
  //Beam spot variables
  mytree->Branch("sigmaZ",&sigmaZ,"sigmaZ/F");
  mytree->Branch("sigmaZ0Error",&sigmaZ0Error,"sigmaZ0Error/F");
  mytree->Branch("sq",&sq,"sq/F");
  mytree->Branch("bsposx",&bsposx,"bsposx/F");
  mytree->Branch("bsposy",&bsposy,"bsposy/F");
  mytree->Branch("bsposz",&bsposz,"bsposz/F");
  //Primary vertex variables
  mytree->Branch("pvsize", &pvsize, "pvsize/I");
  mytree->Branch("pvx", pvx, "pvx[pvsize]/F");
  mytree->Branch("pvy", pvy, "pvy[pvsize]/F");
  mytree->Branch("pvz", pvz, "pvz[pvsize]/F");
  mytree->Branch("pv_isValid", pv_isValid, "pv_isValid[pvsize]/O");
  mytree->Branch("pv_ndof", pv_ndof, "pv_ndof[pvsize]/I");
  mytree->Branch("pv_nTracks", pv_nTracks, "pv_nTracks[pvsize]/I");
  mytree->Branch("pv_normChi2", pv_normChi2, "pv_normChi2[pvsize]/F");
  mytree->Branch("pv_totTrackSize", pv_totTrackSize, "pv_totTrackSize[pvsize]/I"); 
 
  //AKT JETS 
//   mytree->Branch("jetAKT_size", &jetAKT_size, "jetAKT_size/I");
//   mytree->Branch("jetAKT_pt", jetAKT_pt, "jetAKT_pt[jetAKT_size]/F");
//   mytree->Branch("jetAKT_eta", jetAKT_eta, "jetAKT_eta[jetAKT_size]/F");
//   mytree->Branch("jetAKT_phi", jetAKT_phi, "jetAKT_phi[jetAKT_size]/F");
//   mytree->Branch("jetAKT_em", jetAKT_em, "jetAKT_em[jetAKT_size]/F");
//   mytree->Branch("nJetsAKT_pt15", &nJetsAKT_pt15, "nJetsAKT_pt15/I");

  //IC5
  //  mytree->Branch("jetIC5_size", &jetIC5_size, "jetIC5_size/I");
  //   mytree->Branch("jetIC5_pt", jetIC5_pt, "jetIC5_pt[jetIC5_size]/F");
  //   mytree->Branch("jetIC5_eta", jetIC5_eta, "jetIC5_eta[jetIC5_size]/F");
  //   mytree->Branch("jetIC5_phi", jetIC5_phi, "jetIC5_phi[jetIC5_size]/F");
  //   mytree->Branch("jetIC5_em", jetIC5_em, "jetIC5_em[jetIC5_size]/F");
 
  //PF JET
  mytree->Branch("pfJetColl_size", &pfJetColl_size, "pfJetColl_size/I");
  mytree->Branch("pfJet_pt", pfJet_pt, "pfJet_pt[pfJetColl_size]/F");
  mytree->Branch("pfJet_eta", pfJet_eta, "pfJet_eta[pfJetColl_size]/F");
  mytree->Branch("pfJet_phi", pfJet_phi, "pfJet_phi[pfJetColl_size]/F");
  //BTAG
  mytree->Branch("JetColl_size", &JetColl_size, "JetColl_size/I");
  mytree->Branch("Jet_pt", Jet_pt, "Jet_pt[JetColl_size]/F");
  mytree->Branch("Jet_eta", Jet_eta, "Jet_eta[JetColl_size]/F");
  mytree->Branch("Jet_phi", Jet_phi, "Jet_phi[JetColl_size]/F");
  mytree->Branch("Jet_vx", Jet_vx, "Jet_vx[JetColl_size]/F");
  mytree->Branch("Jet_vy", Jet_vy, "Jet_vy[JetColl_size]/F");
  mytree->Branch("Jet_vz", Jet_vz, "Jet_vz[JetColl_size]/F");
  //  mytree->Branch("Jet_em", Jet_em, "Jet_em[JetColl_size]/F");
  mytree->Branch("tCHighEffBTags", tCHighEffBTags, "tCHighEffBTags[JetColl_size]/F");
  mytree->Branch("tCHighPurBTags", tCHighPurBTags, "tCHighPurBTags[JetColl_size]/F");
  mytree->Branch("jetProbBTags", jetProbBTags, "jetProbBTags[JetColl_size]/F");
  mytree->Branch("jetBProbBTags", jetBProbBTags, "jetBProbBTags[JetColl_size]/F");
  mytree->Branch("sSecVertHighEffBTags", sSecVertHighEffBTags, "sSecVertHighEffBTags[JetColl_size]/F");
  mytree->Branch("sSecVertHighPurBTags", sSecVertHighPurBTags, "sSecVertHighPurBTags[JetColl_size]/F");
  mytree->Branch("cSecVertBTags", cSecVertBTags, "cSecVertBTags[JetColl_size]/F");
  mytree->Branch("cSecVertMVABTags", cSecVertMVABTags, "cSecVertMVABTags[JetColl_size]/F");
  mytree->Branch("ghostTrkBTags", ghostTrkBTags, "ghostTrkBTags[JetColl_size]/F");
  mytree->Branch("softEleIP3dBTags", softEleIP3dBTags, "softEleIP3dBTags[JetColl_size]/F");
  mytree->Branch("softElePtBTags", softElePtBTags, "softElePtBTags[JetColl_size]/F");
  mytree->Branch("softMuBTags", softMuBTags, "softMuBTags[JetColl_size]/F");
  mytree->Branch("softMuIP3dBTags", softMuIP3dBTags, "softMuIP3dBTags[JetColl_size]/F");
  mytree->Branch("softMuPtBTags", softMuPtBTags, "softMuPtBTags[JetColl_size]/F");
  

  //MUONS
  mytree->Branch("muon_size", &muon_size, "muon_size/I");
  mytree->Branch("muon_pt", muon_pt, "muon_pt[muon_size]/F");
  mytree->Branch("muon_ptError", muon_ptError, "muon_ptError[muon_size]/F");
  mytree->Branch("muon_gTrk_pt", muon_gTrk_pt, "muon_gTrk_pt[muon_size]/F");
  mytree->Branch("muon_gTrk_ptError", muon_gTrk_ptError, "muon_gTrk_ptError[muon_size]/F");
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
  mytree->Branch("muon_nlayerswithhits", muon_nlayerswithhits, "muon_nlayerswithhits[muon_size]/I");
  mytree->Branch("muon_nlosthits", muon_nlosthits, "muon_nlosthits[muon_size]/I");
  mytree->Branch("muon_nSegmentMatch", muon_nSegmentMatch, "muon_nSegmentMatch[muon_size]/I");
  mytree->Branch("muon_isTrackerMuon", muon_isTrackerMuon, "muon_isTrackerMuon[muon_size]/O");
  mytree->Branch("muon_isPFMuon", muon_isPFMuon, "muon_isPFMuon[muon_size]/O");
  mytree->Branch("muon_isPFIsolationValid", muon_isPFIsolationValid, "muon_isPFIsolationValid[muon_size]/O");
  mytree->Branch("muon_chi2", muon_chi2, "muon_chi2[muon_size]/F");
  mytree->Branch("muon_ndof", muon_ndof, "muon_ndof[muon_size]/I");
  mytree->Branch("muon_normChi2", muon_normChi2, "muon_normChi2[muon_size]/F");
  mytree->Branch("muon_d0", muon_d0, "muon_d0[muon_size]/F");
  mytree->Branch("muon_d0Error", muon_d0Error, "muon_d0Error[muon_size]/F");
  mytree->Branch("muon_dz_cmsCenter", muon_dz_cmsCenter, "muon_dz_cmsCenter[muon_size]/F");
  mytree->Branch("muon_dz_beamSpot", muon_dz_beamSpot, "muon_dz_beamSpot[muon_size]/F");
  mytree->Branch("muon_dz_firstPVtx", muon_dz_firstPVtx, "muon_dz_firstPVtx[muon_size]/F");
  mytree->Branch("muon_dz_firstPVtxwithBS", muon_dz_firstPVtxwithBS, "muon_dz_firstPVtxwithBS[muon_size]/F");
  mytree->Branch("muon_dzError", muon_dzError, "muon_dzError[muon_size]/F");
  mytree->Branch("muon_dxy_cmsCenter", muon_dxy_cmsCenter, "muon_dxy_cmsCenter[muon_size]/F");
  mytree->Branch("muon_dxy_beamSpot", muon_dxy_beamSpot, "muon_dxy_beamSpot[muon_size]/F");
  mytree->Branch("muon_dxy_firstPVtx", muon_dxy_firstPVtx, "muon_dxy_firstPVtx[muon_size]/F");
  mytree->Branch("muon_dxy_firstPVtxwithBS", muon_dxy_firstPVtxwithBS, "muon_dxy_firstPVtxwithBS[muon_size]/F");
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
  mytree->Branch("muMatch_hltL1Mu3p5EG12L3Filtered22", muMatch_hltL1Mu3p5EG12L3Filtered22, "muMatch_hltL1Mu3p5EG12L3Filtered22[muon_size]/O"); 
  mytree->Branch("muMatch_hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q", muMatch_hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q, "muMatch_hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q[muon_size]/O"); 

  //SC VARIABLES
  mytree->Branch("scsize",&scsize,"scsize/I");
  mytree->Branch("scenergy",scenergy,"scenergy[scsize]/F");
  mytree->Branch("sceta",sceta,"sceta[scsize]/F");
  mytree->Branch("scetacorr",scetacorr,"scetacorr[scsize]/F");
  //  mytree->Branch("sctheta",sctheta,"sctheta[scsize]/F");
  //mytree->Branch("scthetacorr",scthetacorr,"scthetacorr[scsize]/F");
  mytree->Branch("scet",scet,"scet[scsize]/F");
  mytree->Branch("scphi",scphi,"scphi[scsize]/F");
  mytree->Branch("scpx",scpx,"scpx[scsize]/F");
  mytree->Branch("scpy",scpy,"scpy[scsize]/F");
  mytree->Branch("scpz",scpz,"scpz[scsize]/F");
  mytree->Branch("scx",scx,"scx[scsize]/F");
  mytree->Branch("scy",scy,"scy[scsize]/F");
  mytree->Branch("scz",scz,"scz[scsize]/F");
  //mytree->Branch("scgsfmatched",scgsfmatched,"scgsfmatched[scsize]/F");  


  //PF Ele Variables 
  mytree->Branch("pfele_size",&pfele_size, "pfele_size/I");
  mytree->Branch("pfele_pt",pfele_pt, "pfele_pt[pfele_size]/F");
  mytree->Branch("pfele_eta",pfele_eta, "pfele_eta[pfele_size]/F");
  mytree->Branch("pfele_phi",pfele_phi, "pfele_phi[pfele_size]/F");
  mytree->Branch("pfele_charge",pfele_charge, "pfele_charge[pfele_size]/I");

  //GSF VARIABLES
  mytree->Branch("gsf_size",&gsf_size, "gsf_size/I");
  mytree->Branch("gsf_isEB", gsf_isEB, "gsf_isEB[gsf_size]/O");
  mytree->Branch("gsf_isEE", gsf_isEE, "gsf_isEE[gsf_size]/O");
  mytree->Branch("gsf_px", gsf_px, "gsf_px[gsf_size]/F");
  mytree->Branch("gsf_py", gsf_py, "gsf_py[gsf_size]/F");
  mytree->Branch("gsf_pz", gsf_pz, "gsf_pz[gsf_size]/F");
  mytree->Branch("gsf_pt", gsf_pt, "gsf_pt[gsf_size]/F");
  mytree->Branch("gsf_eta", gsf_eta, "gsf_eta[gsf_size]/F");
  mytree->Branch("gsf_phi", gsf_phi, "gsf_phi[gsf_size]/F");
  mytree->Branch("gsf_theta", gsf_theta, "gsf_theta[gsf_size]/F");
  mytree->Branch("gsf_charge", gsf_charge, "gsf_charge[gsf_size]/I");
  mytree->Branch("gsf_deltaEtaATcalo", gsf_deltaEtaATcalo, "gsf_deltaEtaATcalo[gsf_size]/F");
  mytree->Branch("gsf_deltaPhiATcalo", gsf_deltaPhiATcalo, "gsf_deltaPhiATcalo[gsf_size]/F");
  mytree->Branch("gsf_sigmaetaeta", gsf_sigmaetaeta, "gsf_sigmaetaeta[gsf_size]/F");
  mytree->Branch("gsf_sigmaIetaIeta", gsf_sigmaIetaIeta, "gsf_sigmaIetaIeta[gsf_size]/F");
  mytree->Branch("gsf_ecalEnergy", gsf_ecalEnergy, "gsf_ecalEnergy[gsf_size]/F");
  mytree->Branch("gsf_eOVERp", gsf_eOVERp, "gsf_eOVERp[gsf_size]/F");
  mytree->Branch("gsf_dxy", gsf_dxy, "gsf_dxy[gsf_size]/F");
  mytree->Branch("gsf_dxy_beamSpot", gsf_dxy_beamSpot, "gsf_dxy_beamSpot[gsf_size]/F");
  mytree->Branch("gsf_dxy_firstPVtx", gsf_dxy_firstPVtx, "gsf_dxy_firstPVtx[gsf_size]/F");
  mytree->Branch("gsf_dxy_firstPVtxwithBS", gsf_dxy_firstPVtxwithBS, "gsf_dxy_firstPVtxwithBS[gsf_size]/F");
  mytree->Branch("gsf_dxyError", gsf_dxyError, "gsf_dxyError[gsf_size]/F");
  mytree->Branch("gsf_dz", gsf_dz, "gsf_dz[gsf_size]/F");
  mytree->Branch("gsf_dz_beamSpot", gsf_dz_beamSpot, "gsf_dz_beamSpot[gsf_size]/F");
  mytree->Branch("gsf_dz_firstPVtx", gsf_dz_firstPVtx, "gsf_dz_firstPVtx[gsf_size]/F");
  mytree->Branch("gsf_dz_firstPVtxwithBS", gsf_dz_firstPVtxwithBS, "gsf_dz_firstPVtxwithBS[gsf_size]/F");
  mytree->Branch("gsf_dzError", gsf_dzError, "gsf_dzError[gsf_size]/F");
  mytree->Branch("gsf_vz", gsf_vz, "gsf_vz[gsf_size]/F");
  mytree->Branch("gsf_nHits", gsf_nHits, "gsf_nHits[gsf_size]/I");
  mytree->Branch("gsf_nLostInnerHits", gsf_nLostInnerHits, "gsf_nLostInnerHits[gsf_size]/I");
  mytree->Branch("gsf_nLostOuterHits", gsf_nLostOuterHits, "gsf_nLostOuterHits[gsf_size]/I");
  mytree->Branch("gsf_convFlags", gsf_convFlags, "gsf_convFlags[gsf_size]/I");
  mytree->Branch("gsf_convDist", gsf_convDist, "gsf_convDist[gsf_size]/F");
  mytree->Branch("gsf_convDcot", gsf_convDcot, "gsf_convDcot[gsf_size]/F");
  mytree->Branch("gsf_convRadius", gsf_convRadius, "gsf_convRadius[gsf_size]/F");
  mytree->Branch("gsf_fBrem", gsf_fBrem, "gsf_fBrem[gsf_size]/F");

  mytree->Branch("gsf_e1x5", gsf_e1x5, "gsf_e1x5[gsf_size]/F");
  mytree->Branch("gsf_e2x5", gsf_e2x5, "gsf_e2x5[gsf_size]/F");
  mytree->Branch("gsf_e5x5", gsf_e5x5, "gsf_e5x5[gsf_size]/F");
  //mytree->Branch("gsf_eMax", gsf_eMax, "gsf_eMax[gsf_size]/F");
  mytree->Branch("gsf_e1x3", gsf_e1x3, "gsf_e1x3[gsf_size]/F");
  //mytree->Branch("gsf_e3x1", gsf_e3x1, "gsf_e3x1[gsf_size]/F");
  //mytree->Branch("gsf_e2x2", gsf_e2x2, "gsf_e2x2[gsf_size]/F");
  //mytree->Branch("gsf_e3x2", gsf_e3x2, "gsf_e3x2[gsf_size]/F");
  //mytree->Branch("gsf_e3x3", gsf_e3x3, "gsf_e3x3[gsf_size]/F");
  //mytree->Branch("gsf_e4x4", gsf_e4x4, "gsf_e4x4[gsf_size]/F");
  //mytree->Branch("gsf_e2x5Right", gsf_e2x5Right, "gsf_e2x5Right[gsf_size]/F");
  //mytree->Branch("gsf_e2x5Left", gsf_e2x5Left, "gsf_e2x5Left[gsf_size]/F");  
  //mytree->Branch("gsf_e2x5Top", gsf_e2x5Top, "gsf_e2x5Top[gsf_size]/F");  
  //mytree->Branch("gsf_e2x5Bottom", gsf_e2x5Bottom, "gsf_e2x5Bottom[gsf_size]/F");
  //mytree->Branch("gsf_e2x5Max", gsf_e2x5Max, "gsf_e2x5Max[gsf_size]/F");
  //mytree->Branch("gsf_eLeft", gsf_eLeft, "gsf_eLeft[gsf_size]/F");
  //mytree->Branch("gsf_eRight", gsf_eRight, "gsf_eRight[gsf_size]/F");
  //mytree->Branch("gsf_eTop", gsf_eTop, "gsf_eTop[gsf_size]/F");
  //mytree->Branch("gsf_eBottom", gsf_eBottom, "gsf_eBottom[gsf_size]/F");
  //mytree->Branch("gsf_e2nd", gsf_e2nd, "gsf_e2nd[gsf_size]/F");
  mytree->Branch("gsf_p", gsf_p, "gsf_p[gsf_size]/F");
  mytree->Branch("gsf_e", gsf_e, "gsf_e[gsf_size]/F");
  mytree->Branch("gsf_deltaeta", gsf_deltaeta, "gsf_deltaeta[gsf_size]/F");
  mytree->Branch("gsf_deltaphi", gsf_deltaphi, "gsf_deltaphi[gsf_size]/F");
  mytree->Branch("gsf_hovere", gsf_hovere, "gsf_hovere[gsf_size]/F");
  mytree->Branch("gsf_hdepth1overe", gsf_hdepth1overe, "gsf_hdepth1overe[gsf_size]/F");
  mytree->Branch("gsf_hdepth2overe", gsf_hdepth2overe, "gsf_hdepth2overe[gsf_size]/F");
  mytree->Branch("gsf_hovere2012", gsf_hovere2012, "gsf_hovere2012[gsf_size]/F");
  mytree->Branch("gsf_hdepth1overe2012", gsf_hdepth1overe2012, "gsf_hdepth1overe2012[gsf_size]/F");
  mytree->Branch("gsf_hdepth2overe2012", gsf_hdepth2overe2012, "gsf_hdepth2overe2012[gsf_size]/F");
  mytree->Branch("gsf_trackiso", gsf_trackiso, "gsf_trackiso[gsf_size]/F");
  mytree->Branch("gsf_ecaliso", gsf_ecaliso, "gsf_ecaliso[gsf_size]/F");
  mytree->Branch("gsf_hcaliso1", gsf_hcaliso1, "gsf_hcaliso1[gsf_size]/F");
  mytree->Branch("gsf_hcaliso2", gsf_hcaliso2, "gsf_hcaliso2[gsf_size]/F");
  mytree->Branch("gsf_hcaliso12012", gsf_hcaliso12012, "gsf_hcaliso12012[gsf_size]/F");
  mytree->Branch("gsf_hcaliso22012", gsf_hcaliso22012, "gsf_hcaliso22012[gsf_size]/F");
  mytree->Branch("gsf_PFisocharged",gsf_PFisocharged,"gsf_PFisocharged[gsf_size]/F");
  mytree->Branch("gsf_PFisophoton",gsf_PFisophoton,"gsf_PFisophoton[gsf_size]/F");
  mytree->Branch("gsf_PFisoneutral",gsf_PFisoneutral,"gsf_PFisoneutral[gsf_size]/F");
  mytree->Branch("gsf_class", gsf_class, "gsf_class[gsf_size]/F");
  mytree->Branch("gsf_isecaldriven", gsf_isecaldriven, "gsf_isecaldriven[gsf_size]/O");
  mytree->Branch("gsf_istrackerdriven", gsf_istrackerdriven, "gsf_istrackerdriven[gsf_size]/O");
  mytree->Branch("gsfsc_e", gsfsc_e, "gsfsc_e[gsf_size]/F");
  mytree->Branch("gsfsc_pt", gsfsc_pt, "gsfsc_pt[gsf_size]/F");
  mytree->Branch("gsfsc_eta", gsfsc_eta, "gsfsc_eta[gsf_size]/F");
  mytree->Branch("gsfsc_phi", gsfsc_phi, "gsfsc_phi[gsf_size]/F");
  mytree->Branch("gsfsc_px", gsfsc_px, "gsfsc_px[gsf_size]/F");
  mytree->Branch("gsfsc_py", gsfsc_py, "gsfsc_py[gsf_size]/F");
  mytree->Branch("gsfsc_pz", gsfsc_pz, "gsfsc_pz[gsf_size]/F");
  mytree->Branch("gsf_e2x5overe5x5", gsf_e2x5overe5x5, "gsf_e2x5overe5x5[gsf_size]/F");
  mytree->Branch("gsf_e1x5overe5x5", gsf_e1x5overe5x5, "gsf_e1x5overe5x5[gsf_size]/F");
  mytree->Branch("gsf_gsfet", gsf_gsfet, "gsf_gsfet[gsf_size]/F");
  mytree->Branch("gsf_hitsinfo", gsf_hitsinfo, "gsf_hitsinfo[2][25]/i");
  mytree->Branch("scindexforgsf", scindexforgsf, "scindexforgsf[gsf_size]/I");
  mytree->Branch("gsftracksize", &gsftracksize, "gsftracksize/I");
  mytree->Branch("gsftracketa", gsftracketa, "gsftracketa[gsftracksize]/F");
  mytree->Branch("gsftrackphi", gsftrackphi, "gsftrackphi[gsftracksize]/F");  
  mytree->Branch("gsftrackp", gsftrackp, "gsftrackp[gsftracksize]/F");
  mytree->Branch("gsftrackpt", gsftrackpt, "gsftrackpt[gsftracksize]/F");
  mytree->Branch("gsftrackpx", gsftrackpx, "gsftrackpx[gsftracksize]/F");
  mytree->Branch("gsftrackpy", gsftrackpy, "gsftrackpy[gsftracksize]/F");
  mytree->Branch("gsftrackpz", gsftrackpz, "gsftrackpz[gsftracksize]/F");
  mytree->Branch("gsfpass_ET", gsfpass_ET, "gsfpass_ET[gsf_size]/O"); 
  mytree->Branch("gsfpass_PT", gsfpass_PT, "gsfpass_PT[gsf_size]/O");  
  mytree->Branch("gsfpass_DETETA", gsfpass_DETETA, "gsfpass_DETETA[gsf_size]/O");  
  mytree->Branch("gsfpass_CRACK", gsfpass_CRACK, "gsfpass_CRACK[gsf_size]/O"); 
  mytree->Branch("gsfpass_DETAIN", gsfpass_DETAIN, "gsfpass_DETAIN[gsf_size]/O");  
  mytree->Branch("gsfpass_DPHIIN", gsfpass_DPHIIN, "gsfpass_DPHIIN[gsf_size]/O");  
  mytree->Branch("gsfpass_HADEM", gsfpass_HADEM, "gsfpass_HADEM[gsf_size]/O");  
  mytree->Branch("gsfpass_SIGMAIETAIETA", gsfpass_SIGMAIETAIETA, "gsfpass_SIGMAIETAIETA[gsf_size]/O"); 
  mytree->Branch("gsfpass_E2X5OVER5X5", gsfpass_E2X5OVER5X5, "gsfpass_E2X5OVER5X5[gsf_size]/O");  
  mytree->Branch("gsfpass_ISOLEMHADDEPTH1", gsfpass_ISOLEMHADDEPTH1, "gsfpass_ISOLEMHADDEPTH1[gsf_size]/O");  
  mytree->Branch("gsfpass_ISOLHADDEPTH2", gsfpass_ISOLHADDEPTH2, "gsfpass_ISOLHADDEPTH2[gsf_size]/O");  
  mytree->Branch("gsfpass_ISOLPTTRKS", gsfpass_ISOLPTTRKS, "gsfpass_ISOLPTTRKS[gsf_size]/O");  
  mytree->Branch("gsfpass_ECALDRIVEN", gsfpass_ECALDRIVEN, "gsfpass_ECALDRIVEN[gsf_size]/O");  
  mytree->Branch("gsfpass_INVALID", gsfpass_INVALID, "gsfpass_INVALID[gsf_size]/O");  
  mytree->Branch("gsfpass_NOMISSINGHITS", gsfpass_NOMISSINGHITS, "gsfpass_NOMISSINGHITS[gsf_size]/O");  
  mytree->Branch("gsfpass_NOCONVERSION", gsfpass_NOCONVERSION, "gsfpass_NOCONVERSION[gsf_size]/O");  
  mytree->Branch("gsfpass_DXYFIRSTPV", gsfpass_DXYFIRSTPV, "gsfpass_DXYFIRSTPV[gsf_size]/O"); 
  mytree->Branch("gsfpass_HEEP", gsfpass_HEEP, "gsfpass_HEEP[gsf_size]/O");  
  mytree->Branch("gsfpass_ID", gsfpass_ID, "gsfpass_ID[gsf_size]/O");  
  mytree->Branch("gsfpass_ISO", gsfpass_ISO, "gsfpass_ISO[gsf_size]/O");
  mytree->Branch("gsfmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter",gsfmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter,"gsfmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter[gsf_size]/O");
  mytree->Branch("gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter",gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter,"gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter[gsf_size]/O");
  mytree->Branch("gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter,"gsfmatch_hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter[gsf_size]/O");
  mytree->Branch("gsfmatch_hltL1sL1SingleEG22",gsfmatch_hltL1sL1SingleEG22,"gsfmatch_hltL1sL1SingleEG22[gsf_size]/O");
  mytree->Branch("gsfmatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter",gsfmatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter,"gsfmatch_hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter[gsf_size]/O");
  mytree->Branch("gsfmatch_hltEle33CaloIdLPixelMatchFilter",gsfmatch_hltEle33CaloIdLPixelMatchFilter,"gsfmatch_hltEle33CaloIdLPixelMatchFilter[gsf_size]/O"); 
  mytree->Branch("scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter",scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter,"scmatch_hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter[scsize]/O");
  mytree->Branch("gsfmatch_hltEle27WP80TrackIsoFilter",gsfmatch_hltEle27WP80TrackIsoFilter,"gsfmatch_hltEle27WP80TrackIsoFilter[gsf_size]/O"); 
  mytree->Branch("gsfmatch_hltMu22Photon22CaloIdLHEFilter", gsfmatch_hltMu22Photon22CaloIdLHEFilter, "gsfmatch_hltMu22Photon22CaloIdLHEFilter[gsf_size]/O"); 




  mytree->Branch("heepHeepMass", &heepHeepMass, "heepHeepMass/F");

  //CHARGE INFO
  mytree->Branch("scpixcharge", scpixcharge, "scpixcharge[gsf_size]/I");
  mytree->Branch("ctfcharge", ctfcharge, "ctfcharge[gsf_size]/I");
  mytree->Branch("gsfcharge", gsfcharge, "gsfcharge[gsf_size]/I");
  mytree->Branch("gsfctfscpixconsistent", gsfctfscpixconsistent, "gsfctfscpixconsistent[gsf_size]/O");
  mytree->Branch("gsfscpixconsistent", gsfscpixconsistent, "gsfscpixconsistent[gsf_size]/O");
  mytree->Branch("gsfctfconsistent", gsfctfconsistent, "gsfctfconsistent[gsf_size]/O");

  mytree->Branch("genparticles_size", &genparticles_size, "genparticles_size/I");
  mytree->Branch("hardGenEle_size", &hardGenEle_size, "hardGenEle_size/I");
  mytree->Branch("hardGenMu_size", &hardGenMu_size, "hardGenMu_size/I");
  mytree->Branch("genEle_size", &genEle_size, "genEle_size/I");
  mytree->Branch("genMu_size", &genMu_size, "genMu_size/I");
  mytree->Branch("genquarks_size", &genquarks_size, "genquarks_size/I");
  mytree->Branch("gengluons_size", &gengluons_size, "gengluons_size/I");

  // Conversion info 
  mytree->Branch("conv_size",&conv_size,"conv_size/I");
  mytree->Branch("conv_vtxProb",conv_vtxProb,"conv_vtxProb[conv_size]/F");
  mytree->Branch("conv_lxy",conv_lxy,"conv_lxy[conv_size]/F");
  mytree->Branch("conv_nHitsMax",conv_nHitsMax,"conv_nHitsMax[conv_size]/I");
  mytree->Branch("conv_eleind",conv_eleind,"conv_eleind[conv_size]/I");

  // Crystal info 
  mytree->Branch("gsf0_crystal_size",&gsf0_crystal_size,"gsf0_crystal_size/I");
  mytree->Branch("gsf0_crystal_ietaorix",gsf0_crystal_ietaorix,"gsf0_crystal_ietaorix[gsf0_crystal_size]/I");
  mytree->Branch("gsf0_crystal_iphioriy",gsf0_crystal_iphioriy,"gsf0_crystal_iphioriy[gsf0_crystal_size]/I");
  mytree->Branch("gsf0_crystal_energy",gsf0_crystal_energy,"gsf0_crystal_energy[gsf0_crystal_size]/F");
  mytree->Branch("gsf0_crystal_eta",gsf0_crystal_eta,"gsf0_crystal_eta[gsf0_crystal_size]/F");
  mytree->Branch("gsf1_crystal_size",&gsf1_crystal_size,"gsf1_crystal_size/I");
  mytree->Branch("gsf1_crystal_ietaorix",gsf1_crystal_ietaorix,"gsf1_crystal_ietaorix[gsf1_crystal_size]/I");
  mytree->Branch("gsf1_crystal_iphioriy",gsf1_crystal_iphioriy,"gsf1_crystal_iphioriy[gsf1_crystal_size]/I");
  mytree->Branch("gsf1_crystal_energy",gsf1_crystal_energy,"gsf1_crystal_energy[gsf1_crystal_size]/F");
  mytree->Branch("gsf1_crystal_eta",gsf1_crystal_eta,"gsf1_crystal_eta[gsf1_crystal_size]/F");
  //GEN INFO FOR ELE and POSI (after FSR)
  

  mytree->Branch("genquark_e", genquark_e, "genquark_e[genquarks_size]/F");
  mytree->Branch("genquark_eta", genquark_eta, "genquark_eta[genquarks_size]/F");
  mytree->Branch("genquark_phi", genquark_phi, "genquark_phi[genquarks_size]/F");
  mytree->Branch("genquark_pt", genquark_pt, "genquark_pt[genquarks_size]/F");
  mytree->Branch("genquark_px", genquark_px, "genquark_px[genquarks_size]/F");
  mytree->Branch("genquark_py", genquark_py, "genquark_py[genquarks_size]/F");
  mytree->Branch("genquark_pz", genquark_pz, "genquark_pz[genquarks_size]/F");
  mytree->Branch("genquark_status", genquark_status, "genquark_status[genquarks_size]/I");
  mytree->Branch("genquark_charge", genquark_charge, "genquark_charge[genquarks_size]/I");
  mytree->Branch("genquark_pdgid", genquark_pdgid, "genquark_pdgid[genquarks_size]/I");
  mytree->Branch("gengluon_e", gengluon_e, "gengluon_e[gengluons_size]/F");
  mytree->Branch("gengluon_eta", gengluon_eta, "gengluon_eta[gengluons_size]/F");
  mytree->Branch("gengluon_phi", gengluon_phi, "gengluon_phi[gengluons_size]/F");
  mytree->Branch("gengluon_pt", gengluon_pt, "gengluon_pt[gengluons_size]/F");
  mytree->Branch("gengluon_px", gengluon_px, "gengluon_px[gengluons_size]/F");
  mytree->Branch("gengluon_py", gengluon_py, "gengluon_py[gengluons_size]/F");
  mytree->Branch("gengluon_pz", gengluon_pz, "gengluon_pz[gengluons_size]/F");
  mytree->Branch("gengluon_status", gengluon_status, "gengluon_status[gengluons_size]/I");
  mytree->Branch("gengluon_charge", gengluon_charge, "gengluon_charge[gengluons_size]/I");
  mytree->Branch("gengluon_pdgid", gengluon_pdgid, "gengluon_pdgid[gengluons_size]/I");


  mytree->Branch("genele_e", genele_e, "genele_e[genEle_size]/F");
  mytree->Branch("genele_eta", genele_eta, "genele_eta[genEle_size]/F");
  mytree->Branch("genele_phi", genele_phi, "genele_phi[genEle_size]/F");
  mytree->Branch("genele_pt", genele_pt, "genele_pt[genEle_size]/F");
  mytree->Branch("genele_px", genele_px, "genele_px[genEle_size]/F");
  mytree->Branch("genele_py", genele_py, "genele_py[genEle_size]/F");
  mytree->Branch("genele_pz", genele_pz, "genele_pz[genEle_size]/F");
  mytree->Branch("genele_charge", genele_charge, "genele_charge[genEle_size]/I");
  
  //generated variables for the tree (before FSR)   
  mytree->Branch("unstableGenEle_e", unstableGenEle_e, "unstableGenEle_e[genEle_size]/F");
  mytree->Branch("unstableGenEle_eta", unstableGenEle_eta, "unstableGenEle_eta[genEle_size]/F");
  mytree->Branch("unstableGenEle_phi", unstableGenEle_phi, "unstableGenEle_phi[genEle_size]/F");
  mytree->Branch("unstableGenEle_pt", unstableGenEle_pt, "unstableGenEle_pt[genEle_size]/F");
  mytree->Branch("unstableGenEle_px", unstableGenEle_px, "unstableGenEle_px[genEle_size]/F");
  mytree->Branch("unstableGenEle_py", unstableGenEle_py, "unstableGenEle_py[genEle_size]/F");
  mytree->Branch("unstableGenEle_pz", unstableGenEle_pz, "unstableGenEle_pz[genEle_size]/F");
  mytree->Branch("unstableGenEle_charge", unstableGenEle_charge, "unstableGenEle_charge[genEle_size]/I");

  //generated variables for the tree (from hard process)   
  mytree->Branch("hardGenEle_e", hardGenEle_e, "hardGenEle_e[hardGenEle_size]/F");
  mytree->Branch("hardGenEle_eta", hardGenEle_eta, "hardGenEle_eta[hardGenEle_size]/F");
  mytree->Branch("hardGenEle_phi", hardGenEle_phi, "hardGenEle_phi[hardGenEle_size]/F");
  mytree->Branch("hardGenEle_pt", hardGenEle_pt, "hardGenEle_pt[hardGenEle_size]/F");
  mytree->Branch("hardGenEle_px", hardGenEle_px, "hardGenEle_px[hardGenEle_size]/F");
  mytree->Branch("hardGenEle_py", hardGenEle_py, "hardGenEle_py[hardGenEle_size]/F");
  mytree->Branch("hardGenEle_pz", hardGenEle_pz, "hardGenEle_pz[hardGenEle_size]/F");
  mytree->Branch("hardGenEle_charge", hardGenEle_charge, "hardGenEle_charge[hardGenEle_size]/I");
  
  //generated variables for the tree (Z variables)
  mytree->Branch("genelemom_e", genelemom_e, "genelemom_e[genEle_size]/F");
  mytree->Branch("genelemom_eta", genelemom_eta, "genelemom_eta[genEle_size]/F");
  mytree->Branch("genelemom_phi", genelemom_phi, "genelemom_phi[genEle_size]/F");
  mytree->Branch("genelemom_pt", genelemom_pt, "genelemom_pt[genEle_size]/F");
  mytree->Branch("genelemom_px", genelemom_px, "genelemom_px[genEle_size]/F");
  mytree->Branch("genelemom_py", genelemom_py, "genelemom_py[genEle_size]/F");
  mytree->Branch("genelemom_pz", genelemom_pz, "genelemom_pz[genEle_size]/F");
  mytree->Branch("genelemom_charge", genelemom_charge, "genelemom_charge[genEle_size]/I");
  mytree->Branch("genelemom_pdgid", genelemom_pdgid, "genelemom_pdgid[genEle_size]/I");
  mytree->Branch("genelemom_mass", genelemom_mass, "genelemom_mass[genEle_size]/F");

  mytree->Branch("genmu_e", genmu_e, "genmu_e[genMu_size]/F");
  mytree->Branch("genmu_eta", genmu_eta, "genmu_eta[genMu_size]/F");
  mytree->Branch("genmu_phi", genmu_phi, "genmu_phi[genMu_size]/F");
  mytree->Branch("genmu_pt", genmu_pt, "genmu_pt[genMu_size]/F");
  mytree->Branch("genmu_px", genmu_px, "genmu_px[genMu_size]/F");
  mytree->Branch("genmu_py", genmu_py, "genmu_py[genMu_size]/F");
  mytree->Branch("genmu_pz", genmu_pz, "genmu_pz[genMu_size]/F");
  mytree->Branch("genmu_charge", genmu_charge, "genmu_charge[genMu_size]/I");
  
  //generated variables for the tree (before FSR)   
  mytree->Branch("unstableGenMu_e", unstableGenMu_e, "unstableGenMu_e[genMu_size]/F");
  mytree->Branch("unstableGenMu_eta", unstableGenMu_eta, "unstableGenMu_eta[genMu_size]/F");
  mytree->Branch("unstableGenMu_phi", unstableGenMu_phi, "unstableGenMu_phi[genMu_size]/F");
  mytree->Branch("unstableGenMu_pt", unstableGenMu_pt, "unstableGenMu_pt[genMu_size]/F");
  mytree->Branch("unstableGenMu_px", unstableGenMu_px, "unstableGenMu_px[genMu_size]/F");
  mytree->Branch("unstableGenMu_py", unstableGenMu_py, "unstableGenMu_py[genMu_size]/F");
  mytree->Branch("unstableGenMu_pz", unstableGenMu_pz, "unstableGenMu_pz[genMu_size]/F");
  mytree->Branch("unstableGenMu_charge", unstableGenMu_charge, "unstableGenMu_charge[genMu_size]/I");

  //generated variables for the tree (from hard process)   
  mytree->Branch("hardGenMu_e", hardGenMu_e, "hardGenMu_e[hardGenMu_size]/F");
  mytree->Branch("hardGenMu_eta", hardGenMu_eta, "hardGenMu_eta[hardGenMu_size]/F");
  mytree->Branch("hardGenMu_phi", hardGenMu_phi, "hardGenMu_phi[hardGenMu_size]/F");
  mytree->Branch("hardGenMu_pt", hardGenMu_pt, "hardGenMu_pt[hardGenMu_size]/F");
  mytree->Branch("hardGenMu_px", hardGenMu_px, "hardGenMu_px[hardGenMu_size]/F");
  mytree->Branch("hardGenMu_py", hardGenMu_py, "hardGenMu_py[hardGenMu_size]/F");
  mytree->Branch("hardGenMu_pz", hardGenMu_pz, "hardGenMu_pz[hardGenMu_size]/F");
  mytree->Branch("hardGenMu_charge", hardGenMu_charge, "hardGenMu_charge[hardGenMu_size]/I");
  
  //generated variables for the tree (Z variables)
  mytree->Branch("genmumom_e", genmumom_e, "genmumom_e[genMu_size]/F");
  mytree->Branch("genmumom_eta", genmumom_eta, "genmumom_eta[genMu_size]/F");
  mytree->Branch("genmumom_phi", genmumom_phi, "genmumom_phi[genMu_size]/F");
  mytree->Branch("genmumom_pt", genmumom_pt, "genmumom_pt[genMu_size]/F");
  mytree->Branch("genmumom_px", genmumom_px, "genmumom_px[genMu_size]/F");
  mytree->Branch("genmumom_py", genmumom_py, "genmumom_py[genMu_size]/F");
  mytree->Branch("genmumom_pz", genmumom_pz, "genmumom_pz[genMu_size]/F");
  mytree->Branch("genmumom_charge", genmumom_charge, "genmumom_charge[genMu_size]/I");
  mytree->Branch("genmumom_pdgid", genmumom_pdgid, "genmumom_pdgid[genMu_size]/I");
  mytree->Branch("genmumom_mass", genmumom_mass, "genmumom_mass[genMu_size]/F");
  
  //x1 and x2
  mytree->Branch("x1quark", x1quark, "x1quark[genEle_size]/F");
  mytree->Branch("x2quark", x2quark, "x2quark[genEle_size]/F");
  
  mytree->Branch("genPair_mass", &genPair_mass, "genPair_mass/F");
  mytree->Branch("emu_mass", &emu_mass, "emu_mass/F");
  mytree->Branch("res_mass", &res_mass, "res_mass/F");

  mytree->Branch("trueNVtx", &trueNVtx, "trueNVtx/I");
  mytree->Branch("nVtxBefore", &nVtxBefore, "nVtxBefore/I");
  mytree->Branch("nVtxNow", &nVtxNow, "nVtxNow/I");
  mytree->Branch("nVtxAfter", &nVtxAfter, "nVtxAfter/I");
}
 

// ------------ method called once each job just after ending the event loop  ------------
void 
GsfCheckerTree::endJob() {
}

//
void 
GsfCheckerTree::DataGenPart(const edm::Event& e) 
{
  Handle<GenParticleCollection> genParticles;
  e.getByLabel("genParticles", genParticles);
  //cout <<"genParticles->size() " << genParticles->size() << endl;

  genparticles_size = genParticles->size();
  genquarks_size = genParticles->size();
  gengluons_size = genParticles->size();

  genquark_e = new float [genquarks_size];
  genquark_pt = new float [genquarks_size];
  genquark_px = new float [genquarks_size];
  genquark_py = new float [genquarks_size];
  genquark_pz = new float [genquarks_size];
  genquark_eta = new float [genquarks_size];
  genquark_phi = new float [genquarks_size];
  genquark_status = new int [genquarks_size];
  genquark_charge = new int [genquarks_size];
  genquark_pdgid = new int [genquarks_size];
  gengluon_e = new float [gengluons_size];
  gengluon_pt = new float [gengluons_size];
  gengluon_px = new float [gengluons_size];
  gengluon_py = new float [gengluons_size];
  gengluon_pz = new float [gengluons_size];
  gengluon_eta = new float [gengluons_size];
  gengluon_phi = new float [gengluons_size];
  gengluon_status = new int [genquarks_size];
  gengluon_charge = new int [gengluons_size];
  gengluon_pdgid = new int [gengluons_size];

  genele_e = new float [genparticles_size];
  genele_pt = new float [genparticles_size];
  genele_px = new float [genparticles_size];
  genele_py = new float [genparticles_size];
  genele_pz = new float [genparticles_size];
  genele_eta = new float [genparticles_size];
  genele_phi = new float [genparticles_size];
  genele_charge = new int [genparticles_size];
  unstableGenEle_e = new float [genparticles_size];
  unstableGenEle_pt = new float [genparticles_size];
  unstableGenEle_px = new float [genparticles_size];
  unstableGenEle_py = new float [genparticles_size];
  unstableGenEle_pz = new float [genparticles_size];
  unstableGenEle_eta = new float [genparticles_size];
  unstableGenEle_phi = new float [genparticles_size];
  unstableGenEle_charge = new int [genparticles_size];
  hardGenEle_e = new float [genparticles_size];
  hardGenEle_pt = new float [genparticles_size];
  hardGenEle_px = new float [genparticles_size];
  hardGenEle_py = new float [genparticles_size];
  hardGenEle_pz = new float [genparticles_size];
  hardGenEle_eta = new float [genparticles_size];
  hardGenEle_phi = new float [genparticles_size];
  hardGenEle_charge = new int [genparticles_size];
  genelemom_e = new float [genparticles_size];
  genelemom_pt = new float [genparticles_size];
  genelemom_px = new float [genparticles_size];
  genelemom_py = new float [genparticles_size];
  genelemom_pz = new float [genparticles_size];
  genelemom_eta = new float [genparticles_size];
  genelemom_phi = new float [genparticles_size];
  genelemom_charge = new int [genparticles_size];
  genelemom_mass = new float [genparticles_size];
  genelemom_pdgid = new int [genparticles_size];

  genmu_e = new float [genparticles_size];
  genmu_pt = new float [genparticles_size];
  genmu_px = new float [genparticles_size];
  genmu_py = new float [genparticles_size];
  genmu_pz = new float [genparticles_size];
  genmu_eta = new float [genparticles_size];
  genmu_phi = new float [genparticles_size];
  genmu_charge = new int [genparticles_size];
  unstableGenMu_e = new float [genparticles_size];
  unstableGenMu_pt = new float [genparticles_size];
  unstableGenMu_px = new float [genparticles_size];
  unstableGenMu_py = new float [genparticles_size];
  unstableGenMu_pz = new float [genparticles_size];
  unstableGenMu_eta = new float [genparticles_size];
  unstableGenMu_phi = new float [genparticles_size];
  unstableGenMu_charge = new int [genparticles_size];
  hardGenMu_e = new float [genparticles_size];
  hardGenMu_pt = new float [genparticles_size];
  hardGenMu_px = new float [genparticles_size];
  hardGenMu_py = new float [genparticles_size];
  hardGenMu_pz = new float [genparticles_size];
  hardGenMu_eta = new float [genparticles_size];
  hardGenMu_phi = new float [genparticles_size];
  hardGenMu_charge = new int [genparticles_size];
  genmumom_e = new float [genparticles_size];
  genmumom_pt = new float [genparticles_size];
  genmumom_px = new float [genparticles_size];
  genmumom_py = new float [genparticles_size];
  genmumom_pz = new float [genparticles_size];
  genmumom_eta = new float [genparticles_size];
  genmumom_phi = new float [genparticles_size];
  genmumom_charge = new int [genparticles_size];
  genmumom_mass = new float [genparticles_size];
  genmumom_pdgid = new int [genparticles_size];

  x1quark = new float [genparticles_size];
  x2quark = new float [genparticles_size];

  mytree->GetBranch("genquark_e")->SetAddress(genquark_e);
  mytree->GetBranch("genquark_pt")->SetAddress(genquark_pt);
  mytree->GetBranch("genquark_px")->SetAddress(genquark_px);
  mytree->GetBranch("genquark_py")->SetAddress(genquark_py);
  mytree->GetBranch("genquark_pz")->SetAddress(genquark_pz);
  mytree->GetBranch("genquark_eta")->SetAddress(genquark_eta);
  mytree->GetBranch("genquark_phi")->SetAddress(genquark_phi);
  mytree->GetBranch("genquark_status")->SetAddress(genquark_status);
  mytree->GetBranch("genquark_charge")->SetAddress(genquark_charge);
  mytree->GetBranch("genquark_pdgid")->SetAddress(genquark_pdgid);

  mytree->GetBranch("gengluon_e")->SetAddress(gengluon_e);
  mytree->GetBranch("gengluon_pt")->SetAddress(gengluon_pt);
  mytree->GetBranch("gengluon_px")->SetAddress(gengluon_px);
  mytree->GetBranch("gengluon_py")->SetAddress(gengluon_py);
  mytree->GetBranch("gengluon_pz")->SetAddress(gengluon_pz);
  mytree->GetBranch("gengluon_eta")->SetAddress(gengluon_eta);
  mytree->GetBranch("gengluon_phi")->SetAddress(gengluon_phi); 
  mytree->GetBranch("gengluon_status")->SetAddress(gengluon_status);
  mytree->GetBranch("gengluon_charge")->SetAddress(gengluon_charge);
  mytree->GetBranch("gengluon_pdgid")->SetAddress(gengluon_pdgid);

  mytree->GetBranch("genele_e")->SetAddress(genele_e);
  mytree->GetBranch("genele_pt")->SetAddress(genele_pt);
  mytree->GetBranch("genele_px")->SetAddress(genele_px);
  mytree->GetBranch("genele_py")->SetAddress(genele_py);
  mytree->GetBranch("genele_pz")->SetAddress(genele_pz);
  mytree->GetBranch("genele_eta")->SetAddress(genele_eta);
  mytree->GetBranch("genele_phi")->SetAddress(genele_phi);
  mytree->GetBranch("genele_charge")->SetAddress(genele_charge);
  mytree->GetBranch("unstableGenEle_e")->SetAddress(unstableGenEle_e);
  mytree->GetBranch("unstableGenEle_pt")->SetAddress(unstableGenEle_pt);
  mytree->GetBranch("unstableGenEle_px")->SetAddress(unstableGenEle_px);
  mytree->GetBranch("unstableGenEle_py")->SetAddress(unstableGenEle_py);
  mytree->GetBranch("unstableGenEle_pz")->SetAddress(unstableGenEle_pz);
  mytree->GetBranch("unstableGenEle_eta")->SetAddress(unstableGenEle_eta);
  mytree->GetBranch("unstableGenEle_phi")->SetAddress(unstableGenEle_phi);
  mytree->GetBranch("unstableGenEle_charge")->SetAddress(unstableGenEle_charge);
  mytree->GetBranch("hardGenEle_e")->SetAddress(hardGenEle_e);
  mytree->GetBranch("hardGenEle_pt")->SetAddress(hardGenEle_pt);
  mytree->GetBranch("hardGenEle_px")->SetAddress(hardGenEle_px);
  mytree->GetBranch("hardGenEle_py")->SetAddress(hardGenEle_py);
  mytree->GetBranch("hardGenEle_pz")->SetAddress(hardGenEle_pz);
  mytree->GetBranch("hardGenEle_eta")->SetAddress(hardGenEle_eta);
  mytree->GetBranch("hardGenEle_phi")->SetAddress(hardGenEle_phi);
  mytree->GetBranch("hardGenEle_charge")->SetAddress(hardGenEle_charge);
  mytree->GetBranch("genelemom_e")->SetAddress(genelemom_e);
  mytree->GetBranch("genelemom_pt")->SetAddress(genelemom_pt);
  mytree->GetBranch("genelemom_px")->SetAddress(genelemom_px);
  mytree->GetBranch("genelemom_py")->SetAddress(genelemom_py);
  mytree->GetBranch("genelemom_pz")->SetAddress(genelemom_pz);
  mytree->GetBranch("genelemom_eta")->SetAddress(genelemom_eta);
  mytree->GetBranch("genelemom_phi")->SetAddress(genelemom_phi);
  mytree->GetBranch("genelemom_charge")->SetAddress(genelemom_charge);
  mytree->GetBranch("genelemom_mass")->SetAddress(genelemom_mass);
  mytree->GetBranch("genelemom_pdgid")->SetAddress(genelemom_pdgid);

  mytree->GetBranch("genmu_e")->SetAddress(genmu_e);
  mytree->GetBranch("genmu_pt")->SetAddress(genmu_pt);
  mytree->GetBranch("genmu_px")->SetAddress(genmu_px);
  mytree->GetBranch("genmu_py")->SetAddress(genmu_py);
  mytree->GetBranch("genmu_pz")->SetAddress(genmu_pz);
  mytree->GetBranch("genmu_eta")->SetAddress(genmu_eta);
  mytree->GetBranch("genmu_phi")->SetAddress(genmu_phi);
  mytree->GetBranch("genmu_charge")->SetAddress(genmu_charge);
  mytree->GetBranch("unstableGenMu_e")->SetAddress(unstableGenMu_e);
  mytree->GetBranch("unstableGenMu_pt")->SetAddress(unstableGenMu_pt);
  mytree->GetBranch("unstableGenMu_px")->SetAddress(unstableGenMu_px);
  mytree->GetBranch("unstableGenMu_py")->SetAddress(unstableGenMu_py);
  mytree->GetBranch("unstableGenMu_pz")->SetAddress(unstableGenMu_pz);
  mytree->GetBranch("unstableGenMu_eta")->SetAddress(unstableGenMu_eta);
  mytree->GetBranch("unstableGenMu_phi")->SetAddress(unstableGenMu_phi);
  mytree->GetBranch("unstableGenMu_charge")->SetAddress(unstableGenMu_charge);
  mytree->GetBranch("hardGenMu_e")->SetAddress(hardGenMu_e);
  mytree->GetBranch("hardGenMu_pt")->SetAddress(hardGenMu_pt);
  mytree->GetBranch("hardGenMu_px")->SetAddress(hardGenMu_px);
  mytree->GetBranch("hardGenMu_py")->SetAddress(hardGenMu_py);
  mytree->GetBranch("hardGenMu_pz")->SetAddress(hardGenMu_pz);
  mytree->GetBranch("hardGenMu_eta")->SetAddress(hardGenMu_eta);
  mytree->GetBranch("hardGenMu_phi")->SetAddress(hardGenMu_phi);
  mytree->GetBranch("hardGenMu_charge")->SetAddress(hardGenMu_charge);
  mytree->GetBranch("genmumom_e")->SetAddress(genmumom_e);
  mytree->GetBranch("genmumom_pt")->SetAddress(genmumom_pt);
  mytree->GetBranch("genmumom_px")->SetAddress(genmumom_px);
  mytree->GetBranch("genmumom_py")->SetAddress(genmumom_py);
  mytree->GetBranch("genmumom_pz")->SetAddress(genmumom_pz);
  mytree->GetBranch("genmumom_eta")->SetAddress(genmumom_eta);
  mytree->GetBranch("genmumom_phi")->SetAddress(genmumom_phi);
  mytree->GetBranch("genmumom_charge")->SetAddress(genmumom_charge);
  mytree->GetBranch("genmumom_mass")->SetAddress(genmumom_mass);
  mytree->GetBranch("genmumom_pdgid")->SetAddress(genmumom_pdgid);

  mytree->GetBranch("x1quark")->SetAddress(x1quark);
  mytree->GetBranch("x2quark")->SetAddress(x2quark);

  unsigned int counterEle = 0;
  unsigned int counterMu = 0;
  unsigned int counterquark = 0; 
  unsigned int countergluon = 0; 
  unsigned int counterT = 0; 
  unsigned int counterTbar = 0; 
  unsigned int counterHardBoson = 0;
  unsigned int counterHardEle = 0;
  unsigned int counterHardMu = 0;
  int tId = -1;
  int tbarId = -1;
  //int hardBosonId = -1;
  int hardEleId = -1;
  int hardMuId = -1;
  genPair_mass = 0.;
  emu_mass = 0.;
  res_mass = 0.;
  for (size_t i = 0; i < genParticles->size(); ++i) {
    
    const GenParticle & p = (*genParticles)[i];
    int id = p.pdgId();
    int st = p.status(); 

    // find a top and antitop
    if (id == 6) {
      tId = i;
      ++counterT;
    }
    else if (id == -6) {
      tbarId = i;
      ++counterTbar;
    }
  
    // find a hard boson, electron and muon
    if (fabs(id) == 11 && st == 3) {
      hardEleId = i;
      hardGenEle_e[counterHardEle] = p.energy();
      hardGenEle_pt[counterHardEle] = p.pt();
      hardGenEle_px[counterHardEle] = p.px();
      hardGenEle_py[counterHardEle] = p.py();
      hardGenEle_pz[counterHardEle] = p.pz();
      hardGenEle_eta[counterHardEle] = p.eta();
      hardGenEle_phi[counterHardEle] = p.phi();
      hardGenEle_charge[counterHardEle]= p.charge();
      ++counterHardEle;
    }
    else if (fabs(id) == 13 && st == 3) {
      hardMuId = i;
      hardGenMu_e[counterHardMu] = p.energy(); 
      hardGenMu_pt[counterHardMu] = p.pt();
      hardGenMu_px[counterHardMu] = p.px();
      hardGenMu_py[counterHardMu] = p.py();
      hardGenMu_pz[counterHardMu] = p.pz();
      hardGenMu_eta[counterHardMu] = p.eta(); 
      hardGenMu_phi[counterHardMu] = p.phi();
      hardGenMu_charge[counterHardMu]= p.charge();
       ++counterHardMu;
    }
    else if ((fabs(id) == 32 || fabs(id) == 9000006) && st == 3) {
      //hardBosonId = i;
      res_mass = p.mass();
      ++counterHardBoson;
    }
  
    // electrons and their mom
    if (fabs(id) == 11 && st == 1) {
      const Candidate * unstableGenEle = p.clone(); // stable = unstable at the beginning
      const Candidate * mom = p.mother();
 
      while (fabs(mom->pdgId()) == 11) { 
        if(mom->status() ==3 ) unstableGenEle = mom; 
          mom = mom->mother(); 
      }
 
      //cout << "Mass, Pdg Id " <<mom->mass()<< " "<< mom->pdgId()<< endl;
      if(fabs(mom->pdgId()) > 6
         && fabs(mom->pdgId()) != 13 
         && fabs(mom->pdgId()) != 21 
         && fabs(mom->pdgId()) != 22 
         && fabs(mom->pdgId()) != 23 
         && fabs(mom->pdgId()) != 24 
         && fabs(mom->pdgId()) != 32 
         && fabs(mom->pdgId()) != 33 
         && fabs(mom->pdgId()) != 39 
         && fabs(mom->pdgId()) != 5000039 
         && fabs(mom->pdgId()) != 9000006 
         && mom->mass() < 20) continue; 
      if(fabs(mom->pdgId()) == 22 && mom->mass() < 10) continue; //This to remove low mass virtual photons converting to a dielectron pair
      for(int itpart=0; itpart < fabs(mom->numberOfMothers()); ++itpart){
	const Candidate *initpart = mom->mother(itpart); 
	if(fabs(initpart->pdgId()) <9  ) {//Quarks info (Drell-Yan, Z') 
	  genquark_e[counterquark] = initpart->energy(); 
	  genquark_pt[counterquark] = initpart->pt();
	  genquark_px[counterquark] = initpart->px();
	  genquark_py[counterquark] = initpart->py();
	  genquark_pz[counterquark] = initpart->pz();
	  genquark_eta[counterquark] = initpart->eta(); 
	  genquark_phi[counterquark] = initpart->phi();
	  genquark_charge[counterquark]= initpart->charge();
	  genquark_pdgid[counterquark]= initpart->pdgId();
	  genquark_status[counterquark]= initpart->status();
	  counterquark++;
	}
	if(initpart->pdgId() ==21 ) { // For graviton, incoming partons can also be gluons 
	      
	  gengluon_e[countergluon] = initpart->energy(); 
	  gengluon_pt[countergluon] = initpart->pt();
	  gengluon_px[countergluon] = initpart->px();
	  gengluon_py[countergluon] = initpart->py();
	  gengluon_pz[countergluon] = initpart->pz();
	  gengluon_eta[countergluon] = initpart->eta(); 
	  gengluon_phi[countergluon] = initpart->phi();
	  gengluon_charge[countergluon]= initpart->charge();
	  gengluon_pdgid[countergluon]= initpart->pdgId();
	  gengluon_status[countergluon]= initpart->status();
	  countergluon++;
	}
      }
     
      genele_e[counterEle] = p.energy(); 
      genele_pt[counterEle] = p.pt();
      genele_px[counterEle] = p.px();
      genele_py[counterEle] = p.py();
      genele_pz[counterEle] = p.pz();
      genele_eta[counterEle] = p.eta(); 
      genele_phi[counterEle] = p.phi();
      genele_charge[counterEle]= p.charge();
      
      unstableGenEle_e[counterEle] = unstableGenEle->energy(); 
      unstableGenEle_pt[counterEle] = unstableGenEle->pt();
      unstableGenEle_px[counterEle] = unstableGenEle->px();
      unstableGenEle_py[counterEle] = unstableGenEle->py();
      unstableGenEle_pz[counterEle] = unstableGenEle->pz();
      unstableGenEle_eta[counterEle] = unstableGenEle->eta(); 
      unstableGenEle_phi[counterEle] = unstableGenEle->phi();
      unstableGenEle_charge[counterEle]= unstableGenEle->charge();
      
      genelemom_e[counterEle] = mom->energy(); 
      genelemom_pt[counterEle] = mom->pt();
      genelemom_px[counterEle] = mom->px();
      genelemom_py[counterEle] = mom->py();
      genelemom_pz[counterEle] = mom->pz();
      genelemom_eta[counterEle] = mom->eta(); 
      genelemom_phi[counterEle] = mom->phi();
      genelemom_charge[counterEle]= mom->charge();
      genelemom_mass[counterEle]= mom->mass();
      genelemom_pdgid[counterEle]= mom->pdgId();
      
      x1quark[counterEle] = (genelemom_mass[counterEle]*genelemom_mass[counterEle]) / (comEnergy_ * (genelemom_pz[counterEle] + sqrt(genelemom_pz[counterEle]*genelemom_pz[counterEle]+genelemom_mass[counterEle]*genelemom_mass[counterEle] )));
      x2quark[counterEle] = (genelemom_pz[counterEle] + sqrt(genelemom_pz[counterEle]*genelemom_pz[counterEle]+genelemom_mass[counterEle]*genelemom_mass[counterEle] )) / comEnergy_;
      counterEle++;      
    }

    // muons and their mom
    if (fabs(id) == 13 && st == 1) {
      const Candidate * unstableGenMu = p.clone(); // stable = unstable at the beginning
      const Candidate * mom = p.mother();
 
      while (fabs(mom->pdgId()) == 13) { 
        if(mom->status() == 3 ) unstableGenMu = mom; 
          mom = mom->mother(); 
      }
 
      // cut on pdg-id 
      if(fabs(mom->pdgId()) > 6 
         && fabs(mom->pdgId()) != 13 
         && fabs(mom->pdgId()) != 21 
         && fabs(mom->pdgId()) != 22 
         && fabs(mom->pdgId()) != 23 
         && fabs(mom->pdgId()) != 24 
         && fabs(mom->pdgId()) != 32 
         && fabs(mom->pdgId()) != 33 
         && fabs(mom->pdgId()) != 39 
         && fabs(mom->pdgId()) != 5000039
         && fabs(mom->pdgId()) != 9000006
         && mom->mass() < 20) continue; 
      if(fabs(mom->pdgId()) == 22 && mom->mass() < 10) continue;


      genmu_e[counterMu] = p.energy(); 
      genmu_pt[counterMu] = p.pt();
      genmu_px[counterMu] = p.px();
      genmu_py[counterMu] = p.py();
      genmu_pz[counterMu] = p.pz();
      genmu_eta[counterMu] = p.eta(); 
      genmu_phi[counterMu] = p.phi();
      genmu_charge[counterMu]= p.charge();
      
      unstableGenMu_e[counterMu] = unstableGenMu->energy(); 
      unstableGenMu_pt[counterMu] = unstableGenMu->pt();
      unstableGenMu_px[counterMu] = unstableGenMu->px();
      unstableGenMu_py[counterMu] = unstableGenMu->py();
      unstableGenMu_pz[counterMu] = unstableGenMu->pz();
      unstableGenMu_eta[counterMu] = unstableGenMu->eta(); 
      unstableGenMu_phi[counterMu] = unstableGenMu->phi();
      unstableGenMu_charge[counterMu]= unstableGenMu->charge();
      
      genmumom_e[counterMu] = mom->energy(); 
      genmumom_pt[counterMu] = mom->pt();
      genmumom_px[counterMu] = mom->px();
      genmumom_py[counterMu] = mom->py();
      genmumom_pz[counterMu] = mom->pz();
      genmumom_eta[counterMu] = mom->eta(); 
      genmumom_phi[counterMu] = mom->phi();
      genmumom_charge[counterMu]= mom->charge();
      genmumom_mass[counterMu]= mom->mass();
      genmumom_pdgid[counterMu]= mom->pdgId();
      
      counterMu++;      
    }
  }

  genparticles_size=counterEle+counterMu;
  hardGenEle_size=counterHardEle;
  hardGenMu_size=counterHardMu;
  genEle_size=counterEle;
  genMu_size=counterMu;
  genquarks_size = counterquark; 
  gengluons_size = countergluon; 

  // calc invariant mass ot ttbar pair
  //if (counterT > 1 || counterTbar > 1) std::cout << "Found a lot of tops " << counterT << " " << counterTbar << endl;
  if (tId > -1 && tbarId > -1) {
    const GenParticle & top = (*genParticles)[tId];
    const GenParticle & antiTop = (*genParticles)[tbarId];
    TLorentzVector topLv;
    TLorentzVector antiTopLv;

    topLv.SetPxPyPzE(top.px(), top.py(), top.pz(), top.energy());
    antiTopLv.SetPxPyPzE(antiTop.px(), antiTop.py(), antiTop.pz(), antiTop.energy());

    genPair_mass = (float)(topLv + antiTopLv).Mag();
  }

  // calc invariant mass of hard e-mu pair
  if (hardEleId > -1 && hardMuId > -1) {
    const GenParticle & ele = (*genParticles)[hardEleId];
    const GenParticle & mu = (*genParticles)[hardMuId];
    TLorentzVector eleLv;
    TLorentzVector muLv;

    eleLv.SetPxPyPzE(ele.px(), ele.py(), ele.pz(), ele.energy());
    muLv.SetPxPyPzE(mu.px(), mu.py(), mu.pz(), mu.energy());

    emu_mass = (float)(eleLv + muLv).Mag();
    //if (counterHardBoson == 0) emu_mass = 0.;
    //if (counterHardBoson == 0) cout << emu_mass << "    " << (float)(eleLv + muLv).Pt() << "    " << (float)(eleLv + muLv).PseudoRapidity() << "    " << (float)(eleLv + muLv).Phi() << endl;
    //else emu_mass = 0.;
  }
}//end of DataGenPart

//
void
GsfCheckerTree::L1TInfo(const edm::Event &iEvent)
{
  edm::Handle< L1GlobalTriggerReadoutRecord > gtReadoutRecord;
  iEvent.getByLabel( edm::InputTag("gtDigis"), gtReadoutRecord);

  const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = gtReadoutRecord->technicalTriggerWord();
  L1trigger_size = technicalTriggerWordBeforeMask.size();
  //cout << "L1trigger_size " << L1trigger_size << endl;

  L1trigger_bool = new int [L1trigger_size];
  //  mytree->GetBranch("L1trigger_bool")->SetAddress(L1trigger_bool);

  for (unsigned int i = 0; i < technicalTriggerWordBeforeMask.size(); ++i) {
    bool bit = technicalTriggerWordBeforeMask.at(i);
    if (bit == 1) L1trigger_bool[i] = 1;   
    if (bit == 0) L1trigger_bool[i] = 0;   
  }

  //physics declared
  L1GlobalTriggerReadoutRecord const* gtrr = gtReadoutRecord.product();
  L1GtFdlWord fdlWord = gtrr->gtFdlWord();
  if (fdlWord.physicsDeclared() == 1) PhysDecl_bool=1;
  else PhysDecl_bool=0;
} // END of L1TInfo

//
void
GsfCheckerTree::HLTInfo(const edm::Event &iEvent, const edm::EventSetup& iSetup)
{
  HLT_Mu15_eta2p1 = -10;
  HLT_Mu24_eta2p1 = -10;
  HLT_Mu30_eta2p1 = -10;
  HLT_Mu40_eta2p1 = -10;
  HLT_Mu50_eta2p1 = -10;
  HLT_Mu22_TkMu22 = -10;
  HLT_Mu22_Photon22_CaloIdL = -10;
  HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = -10;
  HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = -10;
  HLT_Ele8_CaloIdL_CaloIsoVL = -10;
  HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL = -10;
  HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = -10;
  HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50 = -10;
  HLT_DoubleEle33_CaloIdL = -10;
  HLT_DoubleEle33_CaloIdL_GsfTrkIdVL = -10;
  HLT_DoubleEle33_CaloIdT = -10;
  HLT_Photon20_CaloIdVL_IsoL = -10;
  HLT_Photon30_CaloIdVL = -10;
  HLT_Photon50_CaloIdVL = -10;
  HLT_Photon50_CaloIdVL_IsoL = -10;
  HLT_Photon75_CaloIdVL = -10;
  HLT_Photon90_CaloIdVL = -10;
  HLT_Photon135 = -10;
  HLT_Photon150 = -10;
  HLT_Photon250_NoHE = -10;
  HLT_Photon300_NoHE = -10;
  HLT_Photon26_Photon18 = -10;
  HLT_Photon36_Photon22 = -10;
  HLT_DoublePhoton70 = -10;
  HLT_DoublePhoton80 = -10;
  HLT_Ele27_WP80 = -10;
  
  prescale_HLT_Mu15_eta2p1 = -10;
  prescale_HLT_Mu24_eta2p1 = -10;
  prescale_HLT_Mu30_eta2p1 = -10;
  prescale_HLT_Mu40_eta2p1 = -10;
  prescale_HLT_Mu50_eta2p1 = -10;
  prescale_HLT_Mu22_TkMu22 = -10;
  prescale_HLT_Mu22_Photon22_CaloIdL = -10;
  prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = -10;
  prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = -10;
  prescale_HLT_Ele8_CaloIdL_CaloIsoVL = -10;
  prescale_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL = -10;
  prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = -10;
  prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50 = -10;
  prescale_HLT_DoubleEle33_CaloIdL = -10;
  prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL = -10;
  prescale_HLT_DoubleEle33_CaloIdT = -10;
  prescale_HLT_Photon20_CaloIdVL_IsoL = -10;
  prescale_HLT_Photon30_CaloIdVL = -10;
  prescale_HLT_Photon50_CaloIdVL = -10;
  prescale_HLT_Photon50_CaloIdVL_IsoL = -10;
  prescale_HLT_Photon75_CaloIdVL = -10;
  prescale_HLT_Photon90_CaloIdVL = -10;
  prescale_HLT_Photon135 = -10;
  prescale_HLT_Photon150 = -10;
  prescale_HLT_Photon250_NoHE = -10;
  prescale_HLT_Photon300_NoHE = -10;
  prescale_HLT_Photon26_Photon18 = -10;
  prescale_HLT_Photon36_Photon22 = -10;
  prescale_HLT_DoublePhoton70 = -10;
  prescale_HLT_DoublePhoton80 = -10;
  prescale_HLT_Ele27_WP80 = -10;
 
  hltCount = 0;   
 
  // get hold of TriggerResults
  Handle<TriggerResults> HLTR;
  iEvent.getByLabel(hlTriggerResults_,HLTR);
  if (HLTR.isValid()) {
    hltCount = HLTR->size();
    HLTriggers = new int [hltCount];
    //mytree->GetBranch("HLTriggers")->SetAddress(HLTriggers);
    for (int i = 0; i < hltCount; ++i) HLTriggers[i]  = -10;

    for(int i = 0; i < hltCount; ++i) {
      if (HLTR->accept(i)) HLTriggers[i] = i;
    }

    if (HLTR->wasrun()) nWasRun_++;
    const bool accept(HLTR->accept());
    LogDebug("HLTrigReport") << "HL TriggerResults decision: " << accept;
    if (accept) ++nAccept_;
    if (HLTR->error()) nErrors_++;
  } else {
    LogDebug("HLTrigReport") << "HL TriggerResults with label ["+hlTriggerResults_.encode()+"] not found!";
    nErrors_++;
  }

  // decision for each HL algorithm
  const unsigned int n(hlNames_.size());
  //cout << "hlNames_.size() " << hlNames_.size() << endl; 
  for (unsigned int i=0; i!=n; ++i) {

    if (hlNames_.at(i).find("HLT_Mu15_eta2p1_v") == 0) {
      HLTR->accept(i) ? HLT_Mu15_eta2p1 = 1 : HLT_Mu15_eta2p1 = 0;
      prescale_HLT_Mu15_eta2p1 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Mu24_eta2p1_v") == 0) {
      HLTR->accept(i) ? HLT_Mu24_eta2p1 = 1 : HLT_Mu24_eta2p1 = 0;
      prescale_HLT_Mu24_eta2p1 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Mu30_eta2p1_v") == 0) {
      HLTR->accept(i) ? HLT_Mu30_eta2p1 = 1 : HLT_Mu30_eta2p1 = 0;
      prescale_HLT_Mu30_eta2p1 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Mu40_eta2p1_v") == 0) {
      HLTR->accept(i) ? HLT_Mu40_eta2p1 = 1 : HLT_Mu40_eta2p1 = 0;
      prescale_HLT_Mu40_eta2p1 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Mu50_eta2p1_v") == 0) {
      HLTR->accept(i) ? HLT_Mu50_eta2p1 = 1 : HLT_Mu50_eta2p1 = 0;
      prescale_HLT_Mu50_eta2p1 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Mu22_TkMu22_v") == 0) {
      HLTR->accept(i) ? HLT_Mu22_TkMu22 = 1 : HLT_Mu22_TkMu22 = 0;
      prescale_HLT_Mu22_TkMu22 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Mu22_Photon22_CaloIdL_v") == 0) {
      HLTR->accept(i) ? HLT_Mu22_Photon22_CaloIdL = 1 : HLT_Mu22_Photon22_CaloIdL = 0;
      prescale_HLT_Mu22_Photon22_CaloIdL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") == 0) {
      HLTR->accept(i) ? HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = 1 : HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = 0;
      prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") == 0) {
      HLTR->accept(i) ? HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = 1 : HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = 0;
      prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Ele8_CaloIdL_CaloIsoVL_v") == 0) {
      HLTR->accept(i) ? HLT_Ele8_CaloIdL_CaloIsoVL = 1 : HLT_Ele8_CaloIdL_CaloIsoVL = 0;
      prescale_HLT_Ele8_CaloIdL_CaloIsoVL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL_v") == 0) {
      HLTR->accept(i) ? HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL = 1 : HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL = 0;
      prescale_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") == 0) {
      HLTR->accept(i) ? HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = 1 : HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = 0;
      prescale_HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v") == 0) {
      HLTR->accept(i) ? HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50 = 1 : HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50 = 0;
      prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_DoubleEle33_CaloIdL_v") == 0) {
      HLTR->accept(i) ? HLT_DoubleEle33_CaloIdL = 1 : HLT_DoubleEle33_CaloIdL = 0;
      prescale_HLT_DoubleEle33_CaloIdL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v") == 0) {
      HLTR->accept(i) ? HLT_DoubleEle33_CaloIdL_GsfTrkIdVL = 1 : HLT_DoubleEle33_CaloIdL_GsfTrkIdVL = 0;
      prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_DoubleEle33_CaloIdT_v") == 0) {
      HLTR->accept(i) ? HLT_DoubleEle33_CaloIdT = 1 : HLT_DoubleEle33_CaloIdT = 0;
      prescale_HLT_DoubleEle33_CaloIdT = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon20_CaloIdVL_IsoL_v") == 0) {
      HLTR->accept(i) ? HLT_Photon20_CaloIdVL_IsoL = 1 : HLT_Photon20_CaloIdVL_IsoL = 0;
      prescale_HLT_Photon20_CaloIdVL_IsoL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon30_CaloIdVL_v") == 0) {
      HLTR->accept(i) ? HLT_Photon30_CaloIdVL = 1 : HLT_Photon30_CaloIdVL = 0;
      prescale_HLT_Photon30_CaloIdVL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon50_CaloIdVL_v") == 0) {
      HLTR->accept(i) ? HLT_Photon50_CaloIdVL = 1 : HLT_Photon50_CaloIdVL = 0;
      prescale_HLT_Photon50_CaloIdVL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon50_CaloIdVL_IsoL_v") == 0) {
      HLTR->accept(i) ? HLT_Photon50_CaloIdVL_IsoL = 1 : HLT_Photon50_CaloIdVL_IsoL = 0;
      prescale_HLT_Photon50_CaloIdVL_IsoL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon75_CaloIdVL_v") == 0) {
      HLTR->accept(i) ? HLT_Photon75_CaloIdVL = 1 : HLT_Photon75_CaloIdVL = 0;
      prescale_HLT_Photon75_CaloIdVL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon90_CaloIdVL_v") == 0) {
      HLTR->accept(i) ? HLT_Photon90_CaloIdVL = 1 : HLT_Photon90_CaloIdVL = 0;
      prescale_HLT_Photon90_CaloIdVL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon135_v") == 0) {
      HLTR->accept(i) ? HLT_Photon135 = 1 : HLT_Photon135 = 0;
      prescale_HLT_Photon135 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon150_v") == 0) {
      HLTR->accept(i) ? HLT_Photon150 = 1 : HLT_Photon150 = 0;
      prescale_HLT_Photon150 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon250_NoHE_v") == 0) {
      HLTR->accept(i) ? HLT_Photon250_NoHE = 1 : HLT_Photon250_NoHE = 0;
      prescale_HLT_Photon250_NoHE = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon300_NoHE_v") == 0) {
      HLTR->accept(i) ? HLT_Photon300_NoHE = 1 : HLT_Photon300_NoHE = 0;
      prescale_HLT_Photon300_NoHE = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon26_Photon18_v") == 0) {
      HLTR->accept(i) ? HLT_Photon26_Photon18 = 1 : HLT_Photon26_Photon18 = 0;
      prescale_HLT_Photon26_Photon18 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon36_Photon22_v") == 0) {
      HLTR->accept(i) ? HLT_Photon36_Photon22 = 1 : HLT_Photon36_Photon22 = 0;
      prescale_HLT_Photon36_Photon22 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_DoublePhoton70_v") == 0) {
      HLTR->accept(i) ? HLT_DoublePhoton70 = 1 : HLT_DoublePhoton70 = 0;
      prescale_HLT_DoublePhoton70 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_DoublePhoton80_v") == 0) {
      HLTR->accept(i) ? HLT_DoublePhoton80 = 1 : HLT_DoublePhoton80 = 0;
      prescale_HLT_DoublePhoton80 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Ele27_WP80_v") == 0) {
      HLTR->accept(i) ? HLT_Ele27_WP80 = 1 : HLT_Ele27_WP80 = 0;
      prescale_HLT_Ele27_WP80 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }

  }

  for (unsigned int i = 0; i != n; ++i) {
    if (HLTR->wasrun(i)) {
      hlWasRun_[i]++;
      hlWasRunTab[i] = 1;
    } else {
      edm::LogVerbatim("HltNotRun") << "hlNames(" << i << ") = " << hlNames_.at(i);
    }

    if (HLTR->accept(i)) {
      hlAccept_[i]++;
      hlAcceptTab[i] = i;
    }
    if (HLTR->error(i)) {
      hlErrors_[i]++;
      hlErrorTab[i] = 1;
    }
    const int index(static_cast<int>(HLTR->index(i)));
    if (HLTR->accept(i)) {
      if (index >= posL1s_[i]) hltL1s_[i]++;
      if (index >= posPre_[i]) hltPre_[i]++;
    } else {
      if (index > posL1s_[i]) hltL1s_[i]++;
      if (index > posPre_[i]) hltPre_[i]++;
    }
  }
} // END of HLTInfo

//
void
GsfCheckerTree::METData(const edm::Event &iEvent)
{
  calomet = -1.;
  calomet_phi = -1000.;
  met = -1.;
  pfmet = -1.;
  pfmet_phi = -1000.;
  pfmetcor = -1.;
  pfmetcor_phi = -1000.;

  edm::Handle<CaloMETCollection> pCaloMET;
  bool calometisvalid = iEvent.getByLabel("met", pCaloMET);
  const CaloMETCollection *caloMET  = pCaloMET.product();

  edm::Handle<METCollection> pMET;
  bool metisvalid = iEvent.getByLabel("htMetKT4", pMET);
  const METCollection *MET  = pMET.product();

  edm::Handle<PFMETCollection> pPFMET;
  bool pfmetisvalid = iEvent.getByLabel("pfMet", pPFMET);
  const PFMETCollection *PFMET  = pPFMET.product();

  edm::Handle<PFMETCollection> pPFMETcor;
  bool pfmetcorisvalid = iEvent.getByLabel("pfType1CorrectedMet", pPFMETcor);
  const PFMETCollection *PFMETcor  = pPFMETcor.product();


  //CALOMET
  if (calometisvalid) {
    for (CaloMETCollection::const_iterator calometiter = caloMET->begin(); calometiter != caloMET->end(); ++calometiter) {
      calomet = calometiter->et();
      calomet_phi = calometiter->phi();
    }
  }  
  //MET
  if (metisvalid) {
    for (METCollection::const_iterator metiter = MET->begin(); metiter != MET->end(); ++metiter) {
      met = metiter->sumEt();
    }
  } 
  //PFMET
  if (pfmetisvalid) {
    for(PFMETCollection::const_iterator pfmetiter = PFMET->begin(); pfmetiter != PFMET->end(); ++pfmetiter) {
      pfmet = pfmetiter->et();  
      pfmet_phi = pfmetiter->phi();
    }
  } 
  
 //PFMET Type1 corrected 
  if (pfmetcorisvalid) {
    for(PFMETCollection::const_iterator pfmetcoriter = PFMETcor->begin(); pfmetcoriter != PFMETcor->end(); ++pfmetcoriter) {
      pfmetcor = pfmetcoriter->et();  
      pfmetcor_phi = pfmetcoriter->phi();
    }
  } 

} // END of METData

//
void
GsfCheckerTree::OLDJetData(const edm::Event &iEvent)
{
  //IC5
  //  nJetsIC5_pt15 = -1;
//   for (unsigned int i = 0 ; i< 100 ; i++){   
//     jetIC5_pt[i] = -1;
//     jetIC5_eta[i] = -1;
//     jetIC5_phi[i] = -1;
//     jetIC5_em[i] = -1;
//   }
  
  nJetsAKT_pt15 = -1;

  bool caloantiktjetisvalid = false;
  edm::Handle<CaloJetCollection> pCaloAntiKtJets;
  caloantiktjetisvalid = iEvent.getByLabel("ak5CaloJets", pCaloAntiKtJets);
  const CaloJetCollection *caloAntiKtJets = pCaloAntiKtJets.product();
  
  //LOOP ON anti kt jets
  //int FnJetsAKT = -1;
  int FnJetsAKT_pt10 = 0;
  int FnJetsAKT_pt15 = 0;
  int FnJetsAKT_pt20 = 0;

  vector <float> VemJetsAKT_pt10;
  vector <float> VemJetsAKT_pt15;
  vector <float> VemJetsAKT_pt20;

  int index_jetAKT = 0;
  if(caloantiktjetisvalid){
    //FnJetsAKT = caloAntiKtJets->size();
    //cout << "caloAntiKtJets->size() " << caloAntiKtJets->size()<< endl;

    jetAKT_size = caloAntiKtJets->size();

    jetAKT_eta = new float [jetAKT_size];
    jetAKT_pt = new float [jetAKT_size];
    jetAKT_phi = new float [jetAKT_size];
    jetAKT_em = new float [jetAKT_size];
    mytree->GetBranch("jetAKT_eta")->SetAddress(jetAKT_eta);
    mytree->GetBranch("jetAKT_pt")->SetAddress(jetAKT_pt);
    mytree->GetBranch("jetAKT_phi")->SetAddress(jetAKT_phi);
    mytree->GetBranch("jetAKT_em")->SetAddress(jetAKT_em);

    for(CaloJetCollection::const_iterator antiktjetiter = caloAntiKtJets->begin(); antiktjetiter != caloAntiKtJets->end(); ++antiktjetiter){
      if(antiktjetiter->pt() > 10. && fabs(antiktjetiter->eta()) < 3.) {
      FnJetsAKT_pt10++;
      VemJetsAKT_pt10.push_back(antiktjetiter->emEnergyFraction());
      }
      if(antiktjetiter->pt() > 15. && fabs(antiktjetiter->eta()) < 3.) {
      FnJetsAKT_pt15++;
      VemJetsAKT_pt15.push_back(antiktjetiter->emEnergyFraction());
      }
      if(antiktjetiter->pt() > 20. && fabs(antiktjetiter->eta()) < 3.) {
      FnJetsAKT_pt20++;
      VemJetsAKT_pt20.push_back(antiktjetiter->emEnergyFraction());
      }

      //FILL TREE
      if(antiktjetiter->pt() > 10. && fabs(antiktjetiter->eta()) < 3.) {
      jetAKT_pt[index_jetAKT] = antiktjetiter->pt();
      jetAKT_eta[index_jetAKT] = antiktjetiter->eta();
      jetAKT_phi[index_jetAKT] = antiktjetiter->phi();
      jetAKT_em[index_jetAKT] = antiktjetiter->emEnergyFraction();

      index_jetAKT++;
      }
    }
  }

  jetAKT_size = index_jetAKT;
  nJetsAKT_pt15 = FnJetsAKT_pt15;
}

//
void
GsfCheckerTree::JetData(const edm::Event &event)
{
  edm::Handle<reco::JetTagCollection> tCHighEffBTagHandle;
  event.getByLabel("trackCountingHighEffBJetTags", tCHighEffBTagHandle);
  const reco::JetTagCollection &tCHighEffBTag = *(tCHighEffBTagHandle.product());
  edm::Handle<reco::JetTagCollection> tCHighPurBTagHandle;
  event.getByLabel("trackCountingHighPurBJetTags", tCHighPurBTagHandle);
  const reco::JetTagCollection &tCHighPurBTag = *(tCHighPurBTagHandle.product());
  edm::Handle<reco::JetTagCollection> jetProbBTagHandle;
  event.getByLabel("jetProbabilityBJetTags", jetProbBTagHandle);
  const reco::JetTagCollection &jetProbBTag = *(jetProbBTagHandle.product());
  edm::Handle<reco::JetTagCollection> jetBProbBTagHandle;
  event.getByLabel("jetBProbabilityBJetTags", jetBProbBTagHandle);
  const reco::JetTagCollection &jetBProbBTag = *(jetBProbBTagHandle.product());
  edm::Handle<reco::JetTagCollection> sSecVertHighEffBTagHandle;
  event.getByLabel("simpleSecondaryVertexHighEffBJetTags", sSecVertHighEffBTagHandle);
  const reco::JetTagCollection &sSecVertHighEffBTag = *(sSecVertHighEffBTagHandle.product());
  edm::Handle<reco::JetTagCollection> sSecVertHighPurBTagHandle;
  event.getByLabel("simpleSecondaryVertexHighPurBJetTags", sSecVertHighPurBTagHandle);
  const reco::JetTagCollection &sSecVertHighPurBTag = *(sSecVertHighPurBTagHandle.product());
  edm::Handle<reco::JetTagCollection> cSecVertBTagHandle;
  event.getByLabel("combinedSecondaryVertexBJetTags", cSecVertBTagHandle);
  const reco::JetTagCollection &cSecVertBTag = *(cSecVertBTagHandle.product());
  edm::Handle<reco::JetTagCollection> cSecVertMVABTagHandle;
  event.getByLabel("combinedSecondaryVertexMVABJetTags", cSecVertMVABTagHandle);
  const reco::JetTagCollection &cSecVertMVABTag = *(cSecVertMVABTagHandle.product());
  edm::Handle<reco::JetTagCollection> ghostTrkBTagHandle;
  event.getByLabel("ghostTrackBJetTags", ghostTrkBTagHandle);
  const reco::JetTagCollection &ghostTrkBTag = *(ghostTrkBTagHandle.product());
  edm::Handle<reco::JetTagCollection> softEleIP3dBTagHandle;
  event.getByLabel("softElectronByIP3dBJetTags", softEleIP3dBTagHandle);
  const reco::JetTagCollection &softEleIP3dBTag = *(softEleIP3dBTagHandle.product());
  edm::Handle<reco::JetTagCollection> softElePtBTagHandle;
  event.getByLabel("softElectronByPtBJetTags", softElePtBTagHandle);
  const reco::JetTagCollection &softElePtBTag = *(softElePtBTagHandle.product());
  edm::Handle<reco::JetTagCollection> softMuBTagHandle;
  event.getByLabel("softMuonBJetTags", softMuBTagHandle);
  const reco::JetTagCollection &softMuBTag = *(softMuBTagHandle.product());
  edm::Handle<reco::JetTagCollection> softMuIP3dBTagHandle;
  event.getByLabel("softMuonByIP3dBJetTags", softMuIP3dBTagHandle);
  const reco::JetTagCollection &softMuIP3dBTag = *(softMuIP3dBTagHandle.product());
  edm::Handle<reco::JetTagCollection> softMuPtBTagHandle;
  event.getByLabel("softMuonByPtBJetTags", softMuPtBTagHandle);
  const reco::JetTagCollection &softMuPtBTag = *(softMuPtBTagHandle.product());
  //cout << "tCHighEffBTag.size() " << tCHighEffBTag.size() << endl;

  // PF jet data
  edm::Handle<reco::PFJetCollection> pfJets;
  event.getByLabel("ak5PFJets", pfJets);
  //pfJetColl_size = pfJets->size();
  pfJetColl_size = 0.;
  for (reco::PFJetCollection::const_iterator pfJetIt = pfJets->begin(); pfJetIt != pfJets->end(); ++pfJetIt) {
    if (pfJetIt->pt() > jetPtMin_ && fabs(pfJetIt->eta()) < jetEtaMax_) ++pfJetColl_size;
  }
  pfJet_pt = new float [pfJetColl_size];
  pfJet_eta = new float [pfJetColl_size];
  pfJet_phi = new float [pfJetColl_size];
  mytree->GetBranch("pfJet_pt")->SetAddress(pfJet_pt);
  mytree->GetBranch("pfJet_eta")->SetAddress(pfJet_eta);
  mytree->GetBranch("pfJet_phi")->SetAddress(pfJet_phi);
  int pfjCtr = 0;
  for (reco::PFJetCollection::const_iterator pfJetIt = pfJets->begin(); pfJetIt != pfJets->end(); ++pfJetIt) {
    if (!(pfJetIt->pt() > jetPtMin_ && fabs(pfJetIt->eta()) < jetEtaMax_)) continue;
    pfJet_pt[pfjCtr] = pfJetIt->pt();
    pfJet_eta[pfjCtr] = pfJetIt->eta();
    pfJet_phi[pfjCtr] = pfJetIt->phi();
    ++pfjCtr;
  }  

  // Jet and btag data
  JetColl_size = 0;
  for (unsigned int i = 0; i < tCHighEffBTag.size(); ++i) {
    if (tCHighEffBTag[i].first->pt() > bJetPtMin_ && fabs(tCHighEffBTag[i].first->eta()) < bJetEtaMax_) {
      ++JetColl_size;
    }
  }

  Jet_pt = new float [JetColl_size];
  Jet_eta = new float [JetColl_size];
  Jet_phi = new float [JetColl_size];
  Jet_vx = new float [JetColl_size];
  Jet_vy = new float [JetColl_size];
  Jet_vz = new float [JetColl_size];
  //  Jet_em = new float [JetColl_size];
  tCHighEffBTags = new float [JetColl_size];
  tCHighPurBTags = new float [JetColl_size];
  jetProbBTags = new float [JetColl_size];
  jetBProbBTags = new float [JetColl_size];
  sSecVertHighEffBTags = new float [JetColl_size];
  sSecVertHighPurBTags = new float [JetColl_size];
  cSecVertBTags = new float [JetColl_size];
  cSecVertMVABTags = new float [JetColl_size];
  ghostTrkBTags = new float [JetColl_size];
  softEleIP3dBTags = new float [JetColl_size];
  softElePtBTags = new float [JetColl_size];
  softMuBTags = new float [JetColl_size];
  softMuIP3dBTags = new float [JetColl_size];
  softMuPtBTags = new float [JetColl_size];
  //  mytree->GetBranch("Jet_em")->SetAddress(Jet_em);
  mytree->GetBranch("Jet_pt")->SetAddress(Jet_pt);
  mytree->GetBranch("Jet_eta")->SetAddress(Jet_eta);
  mytree->GetBranch("Jet_phi")->SetAddress(Jet_phi);
  mytree->GetBranch("Jet_vx")->SetAddress(Jet_vx);
  mytree->GetBranch("Jet_vy")->SetAddress(Jet_vy);
  mytree->GetBranch("Jet_vz")->SetAddress(Jet_vz);
  mytree->GetBranch("tCHighEffBTags")->SetAddress(tCHighEffBTags);
  mytree->GetBranch("tCHighPurBTags")->SetAddress(tCHighPurBTags);
  mytree->GetBranch("jetProbBTags")->SetAddress(jetProbBTags);
  mytree->GetBranch("jetBProbBTags")->SetAddress(jetBProbBTags);
  mytree->GetBranch("sSecVertHighEffBTags")->SetAddress(sSecVertHighEffBTags);
  mytree->GetBranch("sSecVertHighPurBTags")->SetAddress(sSecVertHighPurBTags);
  mytree->GetBranch("cSecVertBTags")->SetAddress(cSecVertBTags);
  mytree->GetBranch("cSecVertMVABTags")->SetAddress(cSecVertMVABTags);
  mytree->GetBranch("ghostTrkBTags")->SetAddress(ghostTrkBTags);
  mytree->GetBranch("softEleIP3dBTags")->SetAddress(softEleIP3dBTags);
  mytree->GetBranch("softElePtBTags")->SetAddress(softElePtBTags);
  mytree->GetBranch("softMuBTags")->SetAddress(softMuBTags);
  mytree->GetBranch("softMuIP3dBTags")->SetAddress(softMuIP3dBTags);
  mytree->GetBranch("softMuPtBTags")->SetAddress(softMuPtBTags);

  int btagiter = 0; 
  for (unsigned int i = 0; i < tCHighEffBTag.size(); ++i) {
    if (tCHighEffBTag[i].first->pt() > bJetPtMin_ && fabs(tCHighEffBTag[i].first->eta()) < bJetEtaMax_) {
      
      Jet_pt[btagiter] = tCHighEffBTag[i].first->pt();
      Jet_eta[btagiter] = tCHighEffBTag[i].first->eta();
      Jet_phi[btagiter] = tCHighEffBTag[i].first->phi();
      Jet_vx[btagiter] = tCHighEffBTag[i].first->vx();
      Jet_vy[btagiter] = tCHighEffBTag[i].first->vy();
      Jet_vz[btagiter] = tCHighEffBTag[i].first->vz();
      //      Jet_em[btagiter] = tCHighEffBTag[i].first->emEnergyFraction();
      tCHighEffBTags[btagiter] = tCHighEffBTag[i].second;
      tCHighPurBTags[btagiter] = tCHighPurBTag[i].second;
      jetProbBTags[btagiter] = jetProbBTag[i].second;
      jetBProbBTags[btagiter] = jetBProbBTag[i].second;
      sSecVertHighEffBTags[btagiter] = sSecVertHighEffBTag[i].second;
      sSecVertHighPurBTags[btagiter] = sSecVertHighPurBTag[i].second;
      cSecVertBTags[btagiter] = cSecVertBTag[i].second;
      cSecVertMVABTags[btagiter] = cSecVertMVABTag[i].second;
      ghostTrkBTags[btagiter] = ghostTrkBTag[i].second;
      softEleIP3dBTags[btagiter] = softEleIP3dBTag[i].second;
      softElePtBTags[btagiter] = softElePtBTag[i].second;
      softMuBTags[btagiter] = softMuBTag[i].second;
      softMuIP3dBTags[btagiter] = softMuIP3dBTag[i].second;
      softMuPtBTags[btagiter] = softMuPtBTag[i].second;

      ++btagiter;
    }
  }
} //END of BTagData

DEFINE_FWK_MODULE(GsfCheckerTree);
