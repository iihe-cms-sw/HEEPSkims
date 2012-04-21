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
// $Id: GsfCheckerTree.cc,v 1.24 2012/04/12 10:17:44 treis Exp $
//
//Cleaning ladies : Thomas and Laurent
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "UserCode/HEEPSkims/interface/GsfCheckerTree.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPDebug.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPEventHelper.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPEvent.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#define PI 3.141592654
#define TWOPI 6.283185308

using namespace std;
using namespace reco;
using namespace edm;

//Method to sort the gsf electrons
bool 
gsfEtGreater(const reco::GsfElectron &gsf1,const reco::GsfElectron &gsf2)
{
  float et1 = gsf1.caloEnergy() * sin(gsf1.p4().theta());
  float et2 = gsf2.caloEnergy() * sin(gsf2.p4().theta());
  return (et1 > et2);
}

bool 
scEGreater(const reco::SuperCluster *sc1,const reco::SuperCluster *sc2) 
{
  return ((sc1->energy() + sc1->preshowerEnergy()) > (sc2->energy() + sc2->preshowerEnergy()));
}


bool 
refScEGreater(reco::SuperClusterRef sc1,reco::SuperClusterRef sc2) 
{
  return ((sc1->energy() + sc1->preshowerEnergy()) > (sc2->energy() + sc2->preshowerEnergy()));
}


float 
etacorr(float eta, float pvz, float scz) 
{
  return asinh(sinh(eta) * (1. - pvz/scz));
}


GsfCheckerTree::GsfCheckerTree(const edm::ParameterSet& iConfig):
  evtHelper_(),heepEvt_(),nrPass_(0),nrFail_(0)
{
  evtHelper_.setup(iConfig);
  //now do what ever initialization is needed
  eventcounter = 0;

  hlTriggerResults_ = iConfig.getParameter<edm::InputTag>("TriggerResultsTag");
  comEnergy_ = iConfig.getParameter<double>("centerOfMassEnergy");
  bJetPtMin_ = iConfig.getUntrackedParameter<double>("bJetPtMin", 10.);
  eleEtCut_ = iConfig.getUntrackedParameter<double>("electronEtCut", 0.);
  muPtCut_ = iConfig.getUntrackedParameter<double>("muonPtCut", 0.);

  hcalCfg.hOverEConeSize = 0.15;
  hcalCfg.useTowers = true;
  hcalCfg.hcalTowers = edm::InputTag("towerMaker");
  hcalCfg.hOverEPtMin = 0;

  hcalHelper = new ElectronHcalHelper(hcalCfg);
  inputTagIsoDepElectrons_ = iConfig.getParameter< std::vector<edm::InputTag> >("IsoDepElectron");
  inputTagIsoValElectronsPFId_   = iConfig.getParameter< std::vector<edm::InputTag> >("IsoValElectronPF");
 


}







GsfCheckerTree::~GsfCheckerTree()
{
  cout<<"GsfCheckerTree::~GsfCheckerTree"<<endl; 
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
  reco::GsfElectronCollection::const_iterator gsfiterbis = gsfelectrons.begin();
  int cntr = 0; 
  for(; gsfiterbis != gsfelectrons.end(); ++gsfiterbis) {
    cntr++; 

    if (cntr == 1) {
      gsfPtMax = gsfiterbis->caloEnergy()*sin(gsfiterbis->p4().theta());
    }
    if (cntr == 2) {
      gsfPtSecondMax = gsfiterbis->caloEnergy()*sin(gsfiterbis->p4().theta());
    }
  }
  //Missing hits, Invariant Mass cut

  // Get MUONS
  edm::Handle<reco::MuonCollection> muonCollection;
  iEvent.getByLabel("muons",muonCollection);
  const reco::MuonCollection* muons = muonCollection.product();

  float muonPtMax = 0.;
  float muonPtSecondMax = 0.;
  // get the two highes pt muons
  for(reco::MuonCollection::const_iterator muIt = muons->begin(); muIt < muons->end(); ++muIt){
    if (muIt->isGlobalMuon()) {
      if (muIt->globalTrack()->pt() > muonPtMax) {
        muonPtSecondMax = muonPtMax;
        muonPtMax = muIt->globalTrack()->pt();
      }
      else if (muIt->globalTrack()->pt() > muonPtSecondMax) 
        muonPtSecondMax = muIt->globalTrack()->pt();
    }
  }

  // SKIMMING
  if (!(gsfPtMax > eleEtCut_ && gsfPtSecondMax > eleEtCut_)
      && !(gsfPtMax > eleEtCut_ && muonPtMax > muPtCut_)
      && !(muonPtMax > muPtCut_ && muonPtSecondMax > muPtCut_)
     ) return;



  //rho variable
  rho = 0;
  rhoiso = 0;
  edm::Handle<double> rho_;
  edm::Handle<double> rhoiso_;
  bool isrho; 
  bool isrhoiso;
  isrho = iEvent.getByLabel(edm::InputTag("kt6PFJets:rho"),rho_);
  isrhoiso =iEvent.getByLabel(edm::InputTag("kt6PFJetsIso:rho"),rhoiso_);
  if(isrho)   rho =*rho_;
  if(isrhoiso) rhoiso =*rhoiso_;
  //cout << "rho= " << rho << endl;
  //cout << "rhoiso= " << rhoiso << endl;

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
  }

  L1TInfo(iEvent);
  HLTInfo(iEvent, iSetup);

  METData(iEvent);

  JetData(iEvent);
  BTagData(iEvent);

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
  
  pvsize = pvcoll->size();

  pvx = new float [pvsize];
  pvy = new float [pvsize];
  pvz = new float [pvsize];
  pv_isValid = new bool [pvsize];
  pv_ndof = new float [pvsize];
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
  muon_chi2 = new float [muon_size];
  muon_ndof = new int [muon_size];
  muon_normChi2 = new float [muon_size];
  muon_d0 = new float [muon_size];
  muon_d0Error = new float [muon_size];
  muon_dz_cmsCenter = new float [muon_size];
  muon_dz_beamSpot = new float [muon_size];
  muon_dz_firstPVtx = new float [muon_size];
  muon_dzError = new float [muon_size];
  muon_dxy_cmsCenter = new float [muon_size];
  muon_dxy_beamSpot = new float [muon_size];
  muon_dxy_firstPVtx = new float [muon_size];
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
  mytree->GetBranch("muon_pt")->SetAddress(muon_pt);
  mytree->GetBranch("muon_ptError")->SetAddress(muon_ptError);
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
  mytree->GetBranch("muon_chi2")->SetAddress(muon_chi2);
  mytree->GetBranch("muon_ndof")->SetAddress(muon_ndof);
  mytree->GetBranch("muon_normChi2")->SetAddress(muon_normChi2);
  mytree->GetBranch("muon_d0")->SetAddress(muon_d0);
  mytree->GetBranch("muon_d0Error")->SetAddress(muon_d0Error);
  mytree->GetBranch("muon_dz_cmsCenter")->SetAddress(muon_dz_cmsCenter);
  mytree->GetBranch("muon_dz_beamSpot")->SetAddress(muon_dz_beamSpot);
  mytree->GetBranch("muon_dz_firstPVtx")->SetAddress(muon_dz_firstPVtx);
  mytree->GetBranch("muon_dzError")->SetAddress(muon_dzError);
  mytree->GetBranch("muon_dxy_cmsCenter")->SetAddress(muon_dxy_cmsCenter);
  mytree->GetBranch("muon_dxy_beamSpot")->SetAddress(muon_dxy_beamSpot);
  mytree->GetBranch("muon_dxy_firstPVtx")->SetAddress(muon_dxy_firstPVtx);
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

  int index_mu = 0;
  //LOOP OVER MUONS
  for(reco::MuonCollection::const_iterator muIt = muons->begin(); muIt != muons->end(); ++muIt) {
    
    if (muIt->isGlobalMuon()) {
  
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
      muon_nlayerswithhits[index_mu] = muIt->globalTrack()->hitPattern().trackerLayersWithMeasurement();
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
  muon_size = index_mu;
  //cout << "muon_size " << muon_size << endl;
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
  
  scsize = sclusters.size();
  //cout << "scsize " << scsize << endl;
  //cout << "refscsize " << refsclusters.size() << endl;
  //sort all the refSC by energy
  std::sort(refsclusters.begin(),refsclusters.end(),refScEGreater);

  for(std::vector<reco::SuperClusterRef>::const_iterator refsclustersiter = refsclusters.begin();refsclustersiter != refsclusters.end();refsclustersiter++) {

  }
 
  //sort all the SC by energy
  std::sort(sclusters.begin(),sclusters.end(),scEGreater);

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
  //mytree->GetBranch("scgsfmatched")->SetAddress(scgsfmatched);
  //mytree->GetBranch("scseedmatched")->SetAddress(scseedmatched);
  mytree->GetBranch("scenergy")->SetAddress(scenergy);
  mytree->GetBranch("sceta")->SetAddress(sceta);
  mytree->GetBranch("scetacorr")->SetAddress(scetacorr);
  mytree->GetBranch("sctheta")->SetAddress(sctheta);
  mytree->GetBranch("scthetacorr")->SetAddress(scthetacorr);
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

  //trying to see if the sc is seed associated
  
  edm::Handle<GsfTrackCollection> gsfTracksH ;
  iEvent.getByLabel("electronGsfTracks",gsfTracksH) ;
  const GsfTrackCollection *gsftracks = gsfTracksH.product();

  gsftracksize = gsftracks->size();
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
    gsftracketa[v] = gsftrackiter->eta();
    gsftrackphi[v] = gsftrackiter->phi();  
    gsftrackp[v] = gsftrackiter->p();
    gsftrackpt[v] = gsftrackiter->pt();
    gsftrackpx[v] = gsftrackiter->px();
    gsftrackpy[v] = gsftrackiter->py();
    gsftrackpz[v] = gsftrackiter->pz();
    v++;
  }//end of loop on gsf tracks

  double gsfsceta = 0.;
  double gsfscphi = 0.;
  double gsfscenergy = 0.;

  //To remove spikes (ECAL CLUSTER LAZY TOOLS)
  EcalClusterLazyTools lazytool(iEvent,iSetup,InputTag("reducedEcalRecHitsEB"),InputTag("reducedEcalRecHitsEE"));

  gsf_size = gsfelectrons.size();
  //cout << "gsf_size " <<  gsf_size << endl;

  gsf_isEB = new int [gsf_size];
  gsf_isEE = new int [gsf_size];
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
  gsf_dz = new float [gsf_size];
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
  gsf_isecaldriven = new int [gsf_size];
  gsf_istrackerdriven = new int [gsf_size];
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
  gsfpass_HEEP = new bool [gsf_size];
  gsfpass_ID = new bool [gsf_size];
  gsfpass_ISO = new bool [gsf_size];
  //charge information
  scpixcharge = new int [gsf_size];
  ctfcharge = new int [gsf_size];
  gsfcharge = new int [gsf_size];
  gsfctfscpixconsistent = new bool [gsf_size];
  gsfscpixconsistent = new bool [gsf_size];
  gsfctfconsistent = new bool [gsf_size];
 
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
  mytree->GetBranch("gsf_dz")->SetAddress(gsf_dz);
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
  mytree->GetBranch("gsfpass_HEEP")->SetAddress(gsfpass_HEEP);
  mytree->GetBranch("gsfpass_ID")->SetAddress(gsfpass_ID);
  mytree->GetBranch("gsfpass_ISO")->SetAddress(gsfpass_ISO);
  //charge information
  mytree->GetBranch("scpixcharge")->SetAddress(scpixcharge);
  mytree->GetBranch("ctfcharge")->SetAddress(ctfcharge);
  mytree->GetBranch("gsfcharge")->SetAddress(gsfcharge);
  mytree->GetBranch("gsfctfscpixconsistent")->SetAddress(gsfctfscpixconsistent);
  mytree->GetBranch("gsfscpixconsistent")->SetAddress(gsfscpixconsistent);
  mytree->GetBranch("gsfctfconsistent")->SetAddress(gsfctfconsistent);

  int e=0;
  reco::GsfElectronCollection::const_iterator gsfiter = gsfelectrons.begin();
  for(; gsfiter != gsfelectrons.end(); ++gsfiter) {
    gsfsceta = gsfiter->superCluster()->eta();
    gsfscphi = gsfiter->superCluster()->phi();
    gsfscenergy = gsfiter->superCluster()->energy();

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
    if(gsfiter->ecalDrivenSeed())  gsf_isecaldriven[e] = 1; 
    else{gsf_isecaldriven[e] = 0;}
    if(gsfiter->trackerDrivenSeed()) gsf_istrackerdriven[e] = 1;
    else{gsf_istrackerdriven[e] = 0;}
    gsfsc_e[e] = gsfiter->superCluster()->rawEnergy()+gsfiter->superCluster()->preshowerEnergy();
    gsfsc_pt[e] = (gsfiter->superCluster()->rawEnergy()+gsfiter->superCluster()->preshowerEnergy())/cosh(gsfiter->superCluster()->eta());
    gsfsc_eta[e] = gsfiter->superCluster()->eta();
    gsfsc_phi[e] = gsfiter->superCluster()->phi();
    gsfsc_px[e] = gsfsc_pt[e]*cos(gsfsc_phi[e]);
    gsfsc_py[e] = gsfsc_pt[e]*sin(gsfsc_phi[e]);
    gsfsc_pz[e] = (gsfiter->superCluster()->rawEnergy()+gsfiter->superCluster()->preshowerEnergy())*tanh(gsfiter->superCluster()->eta());
    gsf_gsfet[e] = gsfiter->caloEnergy()*sin(gsfiter->p4().theta());
    gsf_theta[e] = gsfiter->theta();
    gsf_isEB[e] = gsfiter->isEB();
    gsf_isEE[e] = gsfiter->isEE();
    gsf_deltaEtaATcalo[e] = gsfiter->deltaEtaSeedClusterTrackAtCalo();
    gsf_deltaPhiATcalo[e] = gsfiter->deltaPhiSeedClusterTrackAtCalo();
    gsf_ecalEnergy[e] = gsfiter->ecalEnergy();
    gsf_eOVERp[e] = gsfiter->eSuperClusterOverP();
    gsf_dxy[e] = gsfiter->gsfTrack()->dxy();
    gsf_dz[e] = gsfiter->gsfTrack()->dz(); 
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

    // HEEP selection v3.2  - 24/10/2011 Thomas
    bool gsfetbarrel = gsf_gsfet[e] > 35.;
    bool gsfetendcap = gsf_gsfet[e] > 40.;
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
    bool ecalisobarrel = (gsf_ecaliso[e]+gsf_hcaliso1[e]) < (2.+0.03*gsf_gsfet[e]);
    bool ecalisoendcap = true;
    if(gsf_gsfet[e] < 50.) {
      ecalisoendcap = (gsf_ecaliso[e]+gsf_hcaliso1[e]) < 2.5;
    }
    else {
      ecalisoendcap = (gsf_ecaliso[e]+gsf_hcaliso1[e]) < (2.5+0.03*(gsf_gsfet[e]-50.));
    }
    bool hcaliso2barrel  = true;
    bool hcaliso2endcap  = true;
    bool trackisobarrel  = gsf_trackiso[e] < 5.;
    bool trackisoendcap  = gsf_trackiso[e] < 5.;
    bool noMissingHits = gsf_nLostInnerHits[e] == 0;
    bool noConversion = gsf_convFlags[e] != 3; 

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
    gsfpass_ID[e] = (gsfpass_DETAIN[e] && gsfpass_DPHIIN[e] && gsfpass_HADEM[e] && gsfpass_SIGMAIETAIETA[e] && gsfpass_E2X5OVER5X5[e]);
    gsfpass_ISO[e] = (gsfpass_ISOLEMHADDEPTH1[e] && gsfpass_ISOLHADDEPTH2[e] && gsfpass_ISOLPTTRKS[e]);

    gsfpass_HEEP[e] = gsfpass_ET[e] && gsfpass_DETAIN[e] && gsfpass_DPHIIN[e] && gsfpass_HADEM[e] && gsfpass_SIGMAIETAIETA[e] && gsfpass_E2X5OVER5X5[e] && gsfpass_ISOLEMHADDEPTH1[e] && gsfpass_ISOLHADDEPTH2[e] && gsfpass_ISOLPTTRKS[e] && gsfpass_NOMISSINGHITS[e];
    
    //charge info
    scpixcharge[e] = gsfiter->scPixCharge();
    if(gsfiter->closestCtfTrackRef().isNonnull()) ctfcharge[e] = gsfiter->closestCtfTrackRef()->charge();
    gsfcharge[e] = gsfiter->gsfTrack()->charge();
    gsfctfscpixconsistent[e] = gsfiter->isGsfCtfScPixChargeConsistent();
    gsfscpixconsistent[e] = gsfiter->isGsfScPixChargeConsistent();
    gsfctfconsistent[e] = gsfiter->isGsfCtfChargeConsistent();

    //increment index for gsf
    e++;
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
  delete [] gsf_dz;
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
  delete [] gsfpass_HEEP;
  delete [] gsfpass_ID;
  delete [] gsfpass_ISO;
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
  
  delete [] muon_pt;
  delete [] muon_ptError;
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
  delete [] muon_chi2;
  delete [] muon_ndof;
  delete [] muon_normChi2;
  delete [] muon_d0;
  delete [] muon_d0Error;
  delete [] muon_dz_cmsCenter;
  delete [] muon_dz_beamSpot;
  delete [] muon_dz_firstPVtx;
  delete [] muon_dzError;
  delete [] muon_dxy_cmsCenter;
  delete [] muon_dxy_beamSpot;
  delete [] muon_dxy_firstPVtx;
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

  delete [] jetAKT_eta;
  delete [] jetAKT_pt;
  delete [] jetAKT_phi;
  delete [] jetAKT_em;

  delete [] bTagJet_et;
  delete [] bTagJet_pt;
  delete [] bTagJet_eta;
  delete [] bTagJet_phi;
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
    delete [] x1quark;
    delete [] x2quark;
  }
 }//end of analyze method


void 
GsfCheckerTree::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
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
      cout << "HLT menu: " << hltConfig_.tableName() << endl;
      for (unsigned int i = 0; i < n; ++i) {
 	hlWasRunTab[i]=0;
	hlAcceptTab[i]=0;
	hlErrorTab[i]=0;
        cout << "hlNames(" << i << ") = " << hlNames_.at(i) << endl;
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
  mytree->Branch("L1trigger_size", &L1trigger_size, "L1trigger_size/I"); 
  mytree->Branch("L1trigger_bool", L1trigger_bool, "L1trigger_bool[L1trigger_size]/I");
  mytree->Branch("PhysDecl_bool", &PhysDecl_bool, "PhysDecl_bool/I");
  mytree->Branch("HLTriggers", HLTriggers, "HLTriggers[hltCount]/I");
  mytree->Branch("nWasRun_",&nWasRun_,"nWasRun_/I");
  mytree->Branch("nAccept_",&nAccept_,"nAccept_/I");
  mytree->Branch("nErrors_",&nErrors_,"nErrors_/I");
  mytree->Branch("hlWasRun_",&hlWasRun_,"hlWasRun_/I");
  mytree->Branch("hlWasRunTab",hlWasRunTab,"hlWasRunTab[400]/I");
  mytree->Branch("hlAccept_",&hlAccept_);
  mytree->Branch("hlAcceptTab",hlAcceptTab,"hlAcceptTab[400]/I");
  mytree->Branch("hlErrorTab",hlErrorTab,"hlErrorTab[200]/I");
  mytree->Branch("hlNamesTab",&hlNamesTab,"hlNamesTab/C");
  mytree->Branch("hlNames_",&hlNames_);
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
  mytree->Branch("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL", &HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL/I");
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
  mytree->Branch("prescale_HLT_Mu15_eta2p1",&prescale_HLT_Mu15_eta2p1,"prescale_HLT_Mu15_eta2p1/I");
  mytree->Branch("prescale_HLT_Mu30_eta2p1",&prescale_HLT_Mu30_eta2p1,"prescale_HLT_Mu30_eta2p1/I");
  mytree->Branch("prescale_HLT_Mu40_eta2p1",&prescale_HLT_Mu40_eta2p1,"prescale_HLT_Mu40_eta2p1/I");
  mytree->Branch("prescale_HLT_Mu22_TkMu22",&prescale_HLT_Mu22_TkMu22,"prescale_HLT_Mu22_TkMu22/I");
  mytree->Branch("prescale_HLT_Mu22_Photon22_CaloIdL",&prescale_HLT_Mu22_Photon22_CaloIdL,"prescale_HLT_Mu22_Photon22_CaloIdL/I");
  mytree->Branch("prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",&prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL/I");
  mytree->Branch("prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL",&prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL/I");
  mytree->Branch("prescale_HLT_Ele8_CaloIdL_CaloIsoVL",&prescale_HLT_Ele8_CaloIdL_CaloIsoVL,"prescale_HLT_Ele8_CaloIdL_CaloIsoVL/I");
  mytree->Branch("prescale_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL",&prescale_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL,"prescale_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL/I");
  mytree->Branch("prescale_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL",&prescale_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"prescale_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL/I");
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


  //GLOBAL PHYSICAL INFO 
  mytree->Branch("rho", &rho, "rho/F");
  mytree->Branch("rhoiso", &rhoiso, "rhoiso/F");
  mytree->Branch("calomet", &calomet, "calomet/F");
  mytree->Branch("calomet_eta", &calomet_eta, "calomet_eta/F");
  mytree->Branch("calomet_phi", &calomet_phi, "calomet_phi/F");
  mytree->Branch("met", &met, "met/F");
  mytree->Branch("pfmet", &pfmet, "pfmet/F");
  mytree->Branch("pfmet_eta", &pfmet_eta, "pfmet_eta/F");
  mytree->Branch("pfmet_phi", &pfmet_phi, "pfmet_phi/F");
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
  mytree->Branch("pv_ndof", pv_ndof, "pv_ndof[pvsize]/F");
  mytree->Branch("pv_nTracks", pv_nTracks, "pv_nTracks[pvsize]/I");
  mytree->Branch("pv_normChi2", pv_normChi2, "pv_normChi2[pvsize]/F");
  mytree->Branch("pv_totTrackSize", pv_totTrackSize, "pv_totTrackSize[pvsize]/I"); 
 
  //AKT JETS 
  mytree->Branch("jetAKT_size", &jetAKT_size, "jetAKT_size/I");
  mytree->Branch("jetAKT_pt", jetAKT_pt, "jetAKT_pt[jetAKT_size]/F");
  mytree->Branch("jetAKT_eta", jetAKT_eta, "jetAKT_eta[jetAKT_size]/F");
  mytree->Branch("jetAKT_phi", jetAKT_phi, "jetAKT_phi[jetAKT_size]/F");
  mytree->Branch("jetAKT_em", jetAKT_em, "jetAKT_em[jetAKT_size]/F");
  mytree->Branch("nJetsAKT_pt15", &nJetsAKT_pt15, "nJetsAKT_pt15/I");
  //IC5
  //  mytree->Branch("jetIC5_size", &jetIC5_size, "jetIC5_size/I");
  //   mytree->Branch("jetIC5_pt", jetIC5_pt, "jetIC5_pt[jetIC5_size]/F");
  //   mytree->Branch("jetIC5_eta", jetIC5_eta, "jetIC5_eta[jetIC5_size]/F");
  //   mytree->Branch("jetIC5_phi", jetIC5_phi, "jetIC5_phi[jetIC5_size]/F");
  //   mytree->Branch("jetIC5_em", jetIC5_em, "jetIC5_em[jetIC5_size]/F");
 
  //BTAG
  mytree->Branch("bTagJetColl_size", &bTagJetColl_size, "bTagJetColl_size/I");
  mytree->Branch("bTagJet_et", bTagJet_et, "bTagJet_et[bTagJetColl_size]/F");
  mytree->Branch("bTagJet_pt", bTagJet_pt, "bTagJet_pt[bTagJetColl_size]/F");
  mytree->Branch("bTagJet_eta", bTagJet_eta, "bTagJet_eta[bTagJetColl_size]/F");
  mytree->Branch("bTagJet_phi", bTagJet_phi, "bTagJet_phi[bTagJetColl_size]/F");
  mytree->Branch("tCHighEffBTags", tCHighEffBTags, "tCHighEffBTags[bTagJetColl_size]/F");
  mytree->Branch("tCHighPurBTags", tCHighPurBTags, "tCHighPurBTags[bTagJetColl_size]/F");
  mytree->Branch("jetProbBTags", jetProbBTags, "jetProbBTags[bTagJetColl_size]/F");
  mytree->Branch("jetBProbBTags", jetBProbBTags, "jetBProbBTags[bTagJetColl_size]/F");
  mytree->Branch("sSecVertHighEffBTags", sSecVertHighEffBTags, "sSecVertHighEffBTags[bTagJetColl_size]/F");
  mytree->Branch("sSecVertHighPurBTags", sSecVertHighPurBTags, "sSecVertHighPurBTags[bTagJetColl_size]/F");
  mytree->Branch("cSecVertBTags", cSecVertBTags, "cSecVertBTags[bTagJetColl_size]/F");
  mytree->Branch("cSecVertMVABTags", cSecVertMVABTags, "cSecVertMVABTags[bTagJetColl_size]/F");
  mytree->Branch("ghostTrkBTags", ghostTrkBTags, "ghostTrkBTags[bTagJetColl_size]/F");
  mytree->Branch("softEleIP3dBTags", softEleIP3dBTags, "softEleIP3dBTags[bTagJetColl_size]/F");
  mytree->Branch("softElePtBTags", softElePtBTags, "softElePtBTags[bTagJetColl_size]/F");
  mytree->Branch("softMuBTags", softMuBTags, "softMuBTags[bTagJetColl_size]/F");
  mytree->Branch("softMuIP3dBTags", softMuIP3dBTags, "softMuIP3dBTags[bTagJetColl_size]/F");
  mytree->Branch("softMuPtBTags", softMuPtBTags, "softMuPtBTags[bTagJetColl_size]/F");
  

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
  mytree->Branch("muon_nlayerswithhits", muon_nlayerswithhits, "muon_nlayerswithhits[muon_size]/I");
  mytree->Branch("muon_nlosthits", muon_nlosthits, "muon_nlosthits[muon_size]/I");
  mytree->Branch("muon_nSegmentMatch", muon_nSegmentMatch, "muon_nSegmentMatch[muon_size]/I");
  mytree->Branch("muon_isTrackerMuon", muon_isTrackerMuon, "muon_isTrackerMuon[muon_size]/O");
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

  //SC VARIABLES
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
  //mytree->Branch("scgsfmatched",scgsfmatched,"scgsfmatched[scsize]/F");  

  //GSF VARIABLES
  mytree->Branch("gsf_size",&gsf_size, "gsf_size/I");
  mytree->Branch("gsf_isEB", gsf_isEB, "gsf_isEB[gsf_size]/I");
  mytree->Branch("gsf_isEE", gsf_isEE, "gsf_isEE[gsf_size]/I");
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
  mytree->Branch("gsf_dz", gsf_dz, "gsf_dz[gsf_size]/F");
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
  mytree->Branch("gsf_isecaldriven", gsf_isecaldriven, "gsf_isecaldriven[gsf_size]/I");
  mytree->Branch("gsf_istrackerdriven", gsf_istrackerdriven, "gsf_istrackerdriven[gsf_size]/I");
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
  mytree->Branch("gsfpass_HEEP", gsfpass_HEEP, "gsfpass_HEEP[gsf_size]/O");  
  mytree->Branch("gsfpass_ID", gsfpass_ID, "gsfpass_ID[gsf_size]/O");  
  mytree->Branch("gsfpass_ISO", gsfpass_ISO, "gsfpass_ISO[gsf_size]/O");
  

  //CHARGE INFO
  mytree->Branch("scpixcharge", scpixcharge, "scpixcharge[gsf_size]/I");
  mytree->Branch("ctfcharge", ctfcharge, "ctfcharge[gsf_size]/I");
  mytree->Branch("gsfcharge", gsfcharge, "gsfcharge[gsf_size]/I");
  mytree->Branch("gsfctfscpixconsistent", gsfctfscpixconsistent, "gsfctfscpixconsistent[gsf_size]/O");
  mytree->Branch("gsfscpixconsistent", gsfscpixconsistent, "gsfscpixconsistent[gsf_size]/O");
  mytree->Branch("gsfctfconsistent", gsfctfconsistent, "gsfctfconsistent[gsf_size]/O");

  mytree->Branch("genparticles_size", &genparticles_size, "genparticles_size/I");
  //GEN INFO FOR ELE and POSI (after FSR)
  
  mytree->Branch("genele_e", genele_e, "genele_e[genparticles_size]/D");
  mytree->Branch("genele_eta", genele_eta, "genele_eta[genparticles_size]/D");
  mytree->Branch("genele_phi", genele_phi, "genele_phi[genparticles_size]/D");
  mytree->Branch("genele_pt", genele_pt, "genele_pt[genparticles_size]/D");
  mytree->Branch("genele_px", genele_px, "genele_px[genparticles_size]/D");
  mytree->Branch("genele_py", genele_py, "genele_py[genparticles_size]/D");
  mytree->Branch("genele_pz", genele_pz, "genele_pz[genparticles_size]/D");
  mytree->Branch("genele_charge", genele_charge, "genele_charge[genparticles_size]/I");
  
  //generated variables for the tree (before FSR)   
  mytree->Branch("unstableGenEle_e", unstableGenEle_e, "unstableGenEle_e[genparticles_size]/D");
  mytree->Branch("unstableGenEle_eta", unstableGenEle_eta, "unstableGenEle_eta[genparticles_size]/D");
  mytree->Branch("unstableGenEle_phi", unstableGenEle_phi, "unstableGenEle_phi[genparticles_size]/D");
  mytree->Branch("unstableGenEle_pt", unstableGenEle_pt, "unstableGenEle_pt[genparticles_size]/D");
  mytree->Branch("unstableGenEle_px", unstableGenEle_px, "unstableGenEle_px[genparticles_size]/D");
  mytree->Branch("unstableGenEle_py", unstableGenEle_py, "unstableGenEle_py[genparticles_size]/D");
  mytree->Branch("unstableGenEle_pz", unstableGenEle_pz, "unstableGenEle_pz[genparticles_size]/D");
  mytree->Branch("unstableGenEle_charge", unstableGenEle_charge, "unstableGenEle_charge[genparticles_size]/I");
  
  //generated variables for the tree (Z variables)
  mytree->Branch("genelemom_e", genelemom_e, "genelemom_e[genparticles_size]/D");
  mytree->Branch("genelemom_eta", genelemom_eta, "genelemom_eta[genparticles_size]/D");
  mytree->Branch("genelemom_phi", genelemom_phi, "genelemom_phi[genparticles_size]/D");
  mytree->Branch("genelemom_pt", genelemom_pt, "genelemom_pt[genparticles_size]/D");
  mytree->Branch("genelemom_px", genelemom_px, "genelemom_px[genparticles_size]/D");
  mytree->Branch("genelemom_py", genelemom_py, "genelemom_py[genparticles_size]/D");
  mytree->Branch("genelemom_pz", genelemom_pz, "genelemom_pz[genparticles_size]/D");
  mytree->Branch("genelemom_charge", genelemom_charge, "genelemom_charge[genparticles_size]/I");
  mytree->Branch("genelemom_pdgid", genelemom_pdgid, "genelemom_pdgid[genparticles_size]/I");
  mytree->Branch("genelemom_mass", genelemom_mass, "genelemom_mass[genparticles_size]/D");
  
  //x1 and x2
  mytree->Branch("x1quark", x1quark, "x1quark[genparticles_size]/F");
  mytree->Branch("x2quark", x2quark, "x2quark[genparticles_size]/F");

  mytree->Branch("trueNVtx", &trueNVtx, "trueNVtx/F");
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

  genele_e = new double [genparticles_size];
  genele_pt = new double [genparticles_size];
  genele_px = new double [genparticles_size];
  genele_py = new double [genparticles_size];
  genele_pz = new double [genparticles_size];
  genele_eta = new double [genparticles_size];
  genele_phi = new double [genparticles_size];
  genele_charge = new int [genparticles_size];
  unstableGenEle_e = new double [genparticles_size];
  unstableGenEle_pt = new double [genparticles_size];
  unstableGenEle_px = new double [genparticles_size];
  unstableGenEle_py = new double [genparticles_size];
  unstableGenEle_pz = new double [genparticles_size];
  unstableGenEle_eta = new double [genparticles_size];
  unstableGenEle_phi = new double [genparticles_size];
  unstableGenEle_charge = new int [genparticles_size];
  genelemom_e = new double [genparticles_size];
  genelemom_pt = new double [genparticles_size];
  genelemom_px = new double [genparticles_size];
  genelemom_py = new double [genparticles_size];
  genelemom_pz = new double [genparticles_size];
  genelemom_eta = new double [genparticles_size];
  genelemom_phi = new double [genparticles_size];
  genelemom_charge = new int [genparticles_size];
  genelemom_mass = new double [genparticles_size];
  genelemom_pdgid = new int [genparticles_size];
  x1quark = new float [genparticles_size];
  x2quark = new float [genparticles_size];
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
  mytree->GetBranch("x1quark")->SetAddress(x1quark);
  mytree->GetBranch("x2quark")->SetAddress(x2quark);

  unsigned int counter = 0; 
  for (size_t i = 0; i < genParticles->size(); ++i) {
 
    const GenParticle & p = (*genParticles)[i];
    int id = p.pdgId();
    int st = p.status(); 

    if (fabs(id) == 11 && st == 1) {
      const Candidate * unstableGenEle = p.clone(); // stable = unstable at the beginning
      const Candidate * mom = p.mother();
 
      while (fabs(mom->pdgId()) == 11) { 
        if(mom->status() ==3 ) unstableGenEle = mom; 
          mom = mom->mother(); 
      }
 
      // cut on pdg-id 
      if(fabs(mom->pdgId()) != 22 && fabs(mom->pdgId()) != 23 && fabs(mom->pdgId()) != 24 && fabs(mom->pdgId()) != 32 && fabs(mom->pdgId()) != 33 && fabs(mom->pdgId()) != 39 && fabs(mom->pdgId()) !=  13) continue;
      
      genele_e[counter] = p.energy(); 
      genele_pt[counter] = p.pt();
      genele_px[counter] = p.px();
      genele_py[counter] = p.py();
      genele_pz[counter] = p.pz();
      genele_eta[counter] = p.eta(); 
      genele_phi[counter] = p.phi();
      genele_charge[counter]= p.charge();
      
      unstableGenEle_e[counter] = unstableGenEle->energy(); 
      unstableGenEle_pt[counter] = unstableGenEle->pt();
      unstableGenEle_px[counter] = unstableGenEle->px();
      unstableGenEle_py[counter] = unstableGenEle->py();
      unstableGenEle_pz[counter] = unstableGenEle->pz();
      unstableGenEle_eta[counter] = unstableGenEle->eta(); 
      unstableGenEle_phi[counter] = unstableGenEle->phi();
      unstableGenEle_charge[counter]= unstableGenEle->charge();
      
      genelemom_e[counter] = mom->energy(); 
      genelemom_pt[counter] = mom->pt();
      genelemom_px[counter] = mom->px();
      genelemom_py[counter] = mom->py();
      genelemom_pz[counter] = mom->pz();
      genelemom_eta[counter] = mom->eta(); 
      genelemom_phi[counter] = mom->phi();
      genelemom_charge[counter]= mom->charge();
      genelemom_mass[counter]= mom->mass();
      genelemom_pdgid[counter]= mom->pdgId();
      
      x1quark[counter] = (genelemom_mass[counter]*genelemom_mass[counter]) / (comEnergy_ * (genelemom_pz[counter] + sqrt(genelemom_pz[counter]*genelemom_pz[counter]+genelemom_mass[counter]*genelemom_mass[counter] )));
      x2quark[counter] = (genelemom_pz[counter] + sqrt(genelemom_pz[counter]*genelemom_pz[counter]+genelemom_mass[counter]*genelemom_mass[counter] )) / comEnergy_;
      counter++;      
    }
  }  
  genparticles_size=counter;
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
  mytree->GetBranch("L1trigger_bool")->SetAddress(L1trigger_bool);

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
  HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = -10;
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
  prescale_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = -10;
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

  hltCount = 0;   
 
  // get hold of TriggerResults
  Handle<TriggerResults> HLTR;
  iEvent.getByLabel(hlTriggerResults_,HLTR);
  if (HLTR.isValid()) {
    hltCount = HLTR->size();
    HLTriggers = new int [hltCount];
    mytree->GetBranch("HLTriggers")->SetAddress(HLTriggers);
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
    if (hlNames_.at(i).find("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v") == 0) {
      HLTR->accept(i) ? HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = 1 : HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = 0;
      prescale_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
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
  }

  for (unsigned int i = 0; i != n; ++i) {
    if (HLTR->wasrun(i)) {
      hlWasRun_[i]++;
      hlWasRunTab[i] = 1;
    } else {
      cout<<"hlNames(i) = "<<i<<", "<<hlNames_.at(i)<<endl;
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
  calomet_eta = -1000.; 
  calomet_phi = -1000.;
  met = -1.;
  pfmet = -1.;
  pfmet_eta = -1000.; 
  pfmet_phi = -1000.;

  edm::Handle<CaloMETCollection> pCaloMET;
  bool calometisvalid = iEvent.getByLabel("met", pCaloMET);
  const CaloMETCollection *caloMET  = pCaloMET.product();

  edm::Handle<METCollection> pMET;
  bool metisvalid = iEvent.getByLabel("htMetKT4", pMET);
  const METCollection *MET  = pMET.product();

  edm::Handle<PFMETCollection> pPFMET;
  bool pfmetisvalid = iEvent.getByLabel("pfMet", pPFMET);
  const PFMETCollection *PFMET  = pPFMET.product();

  //CALOMET
  if (calometisvalid) {
    for (CaloMETCollection::const_iterator calometiter = caloMET->begin(); calometiter != caloMET->end(); ++calometiter) {
      calomet = calometiter->et();
      calomet_eta = calometiter->eta(); 
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
      pfmet_eta = pfmetiter->eta(); 
      pfmet_phi = pfmetiter->phi();
    }
  } 
} // END of METData

//
void
GsfCheckerTree::JetData(const edm::Event &iEvent)
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
GsfCheckerTree::BTagData(const edm::Event &event)
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

  bTagJet_et = new float [tCHighEffBTag.size()];
  bTagJet_pt = new float [tCHighEffBTag.size()];
  bTagJet_eta = new float [tCHighEffBTag.size()];
  bTagJet_phi = new float [tCHighEffBTag.size()];
  tCHighEffBTags = new float [tCHighEffBTag.size()];
  tCHighPurBTags = new float [tCHighEffBTag.size()];
  jetProbBTags = new float [tCHighEffBTag.size()];
  jetBProbBTags = new float [tCHighEffBTag.size()];
  sSecVertHighEffBTags = new float [tCHighEffBTag.size()];
  sSecVertHighPurBTags = new float [tCHighEffBTag.size()];
  cSecVertBTags = new float [tCHighEffBTag.size()];
  cSecVertMVABTags = new float [tCHighEffBTag.size()];
  ghostTrkBTags = new float [tCHighEffBTag.size()];
  softEleIP3dBTags = new float [tCHighEffBTag.size()];
  softElePtBTags = new float [tCHighEffBTag.size()];
  softMuBTags = new float [tCHighEffBTag.size()];
  softMuIP3dBTags = new float [tCHighEffBTag.size()];
  softMuPtBTags = new float [tCHighEffBTag.size()];
  mytree->GetBranch("bTagJet_et")->SetAddress(bTagJet_et);
  mytree->GetBranch("bTagJet_pt")->SetAddress(bTagJet_pt);
  mytree->GetBranch("bTagJet_eta")->SetAddress(bTagJet_eta);
  mytree->GetBranch("bTagJet_phi")->SetAddress(bTagJet_phi);
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

  bTagJetColl_size = 0;
  for (unsigned int i = 0; i < tCHighEffBTag.size(); ++i) {
    if (tCHighEffBTag[i].first->pt() > bJetPtMin_ && fabs(tCHighEffBTag[i].first->eta()) < 3.) {
      bTagJet_et[bTagJetColl_size] = tCHighEffBTag[i].first->et();
      bTagJet_pt[bTagJetColl_size] = tCHighEffBTag[i].first->pt();
      bTagJet_eta[bTagJetColl_size] = tCHighEffBTag[i].first->eta();
      bTagJet_phi[bTagJetColl_size] = tCHighEffBTag[i].first->phi();

      tCHighEffBTags[bTagJetColl_size] = tCHighEffBTag[i].second;
      tCHighPurBTags[bTagJetColl_size] = tCHighPurBTag[i].second;
      jetProbBTags[bTagJetColl_size] = jetProbBTag[i].second;
      jetBProbBTags[bTagJetColl_size] = jetBProbBTag[i].second;
      sSecVertHighEffBTags[bTagJetColl_size] = sSecVertHighEffBTag[i].second;
      sSecVertHighPurBTags[bTagJetColl_size] = sSecVertHighPurBTag[i].second;
      cSecVertBTags[bTagJetColl_size] = cSecVertBTag[i].second;
      cSecVertMVABTags[bTagJetColl_size] = cSecVertMVABTag[i].second;
      ghostTrkBTags[bTagJetColl_size] = ghostTrkBTag[i].second;
      softEleIP3dBTags[bTagJetColl_size] = softEleIP3dBTag[i].second;
      softElePtBTags[bTagJetColl_size] = softElePtBTag[i].second;
      softMuBTags[bTagJetColl_size] = softMuBTag[i].second;
      softMuIP3dBTags[bTagJetColl_size] = softMuIP3dBTag[i].second;
      softMuPtBTags[bTagJetColl_size] = softMuPtBTag[i].second;

      ++bTagJetColl_size;
    }
  }
} //END of BTagData

DEFINE_FWK_MODULE(GsfCheckerTree);
