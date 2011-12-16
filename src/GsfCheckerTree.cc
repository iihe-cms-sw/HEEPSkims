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
// $Id: GsfCheckerTree.cc,v 1.14 2011/12/12 17:36:00 lathomas Exp $
//
//Cleaning ladies : Thomas and Laurent

#include "UserCode/HEEPSkims/interface/GsfCheckerTree.h"

#define PI 3.141592654
#define TWOPI 6.283185308

using namespace std;
using namespace reco;
using namespace edm;
using namespace std;
using namespace reco;

//Method to sort the gsf electrons
bool gsfEtGreater(const reco::GsfElectron &gsf1,const reco::GsfElectron &gsf2)
{
  float et1 = gsf1.caloEnergy()*sin(gsf1.p4().theta());
  float et2 = gsf2.caloEnergy()*sin(gsf2.p4().theta());
  return (et1 > et2);
}

bool scEGreater(const reco::SuperCluster *sc1,const reco::SuperCluster *sc2) {
  return ((sc1->energy()+sc1->preshowerEnergy()) > (sc2->energy()+sc2->preshowerEnergy()));
}


bool refScEGreater(reco::SuperClusterRef sc1,reco::SuperClusterRef sc2) {
  return ((sc1->energy()+sc1->preshowerEnergy()) > (sc2->energy()+sc2->preshowerEnergy()));
}


float etacorr(float eta, float pvz, float scz){
  return asinh(sinh(eta)*(1.-pvz/scz));
}


GsfCheckerTree::GsfCheckerTree(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  eventcounter = 0;

  hlTriggerResults_ = iConfig.getParameter<edm::InputTag> ("TriggerResultsTag");
  comEnergy_ = iConfig.getParameter<double>("centerOfMassEnergy");
  eleEtCut_ = iConfig.getUntrackedParameter<double>("electronEtCut", 0.);
  muPtCut_ = iConfig.getUntrackedParameter<double>("muonPtCut", 0.);
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
  using namespace edm;
  using namespace std;
  using namespace reco;
  bool useGenData_ = !iEvent.isRealData(); 

  eventcounter++;

  //Run and event number
  runnumber = iEvent.id().run();
  eventnumber = iEvent.id().event();
  luminosityBlock = iEvent.id().luminosityBlock(); 

  HLT_Mu15 = -10;
  HLT_Mu30 = -10;
  HLT_Mu40_eta2p1 = -10;
  HLT_Mu15_Photon20_CaloIdL = -10;
  HLT_Mu8_Ele17_CaloIdT_CaloIsoVL = -10;
  HLT_Mu17_Ele8_CaloIdT_CaloIsoVL = -10;
  HLT_Ele8 = -10;
  HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT = -10;
  HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = -10;
  HLT_Ele32_CaloIdL_CaloIsoVL_SC17 = -10;
  HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17 = -10;
  HLT_DoubleEle33_CaloIdL = -10;
  HLT_DoubleEle33_CaloIdL_CaloIsoT = -10;
  HLT_DoubleEle33_CaloIdT = -10;
  HLT_DoubleEle45_CaloIdL = -10;
  HLT_Photon20_CaloIdVL_IsoL = -10;
  HLT_Photon30_CaloIdVL = -10;
  HLT_Photon50_CaloIdVL = -10;
  HLT_Photon50_CaloIdVL_IsoL = -10;
  HLT_Photon75_CaloIdVL = -10;
  HLT_Photon90_CaloIdVL = -10;
  HLT_Photon125 = -10;
  HLT_Photon135 = -10;
  HLT_Photon200_NoHE = -10;
  HLT_Photon225_NoHE = -10;
  HLT_Photon26_Photon18 = -10;
  HLT_Photon36_Photon22 = -10;
  HLT_DoublePhoton33 = -10;
  HLT_DoublePhoton60 = -10;
  HLT_DoublePhoton70 = -10;
  HLT_DoublePhoton80 = -10;

  prescale_HLT_Mu15 = -10;
  prescale_HLT_Mu30 = -10;
  prescale_HLT_Mu40_eta2p1 = -10;
  prescale_HLT_Mu15_Photon20_CaloIdL = -10;
  prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL = -10;
  prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL = -10;
  prescale_HLT_Ele8 = -10;
  prescale_HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT = -10;
  prescale_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = -10;
  prescale_HLT_Ele32_CaloIdL_CaloIsoVL_SC17 = -10;
  prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17 = -10;
  prescale_HLT_DoubleEle33_CaloIdL = -10;
  prescale_HLT_DoubleEle33_CaloIdL_CaloIsoT = -10;
  prescale_HLT_DoubleEle33_CaloIdT = -10;
  prescale_HLT_DoubleEle45_CaloIdL = -10;
  prescale_HLT_Photon20_CaloIdVL_IsoL = -10;
  prescale_HLT_Photon30_CaloIdVL = -10;
  prescale_HLT_Photon50_CaloIdVL = -10;
  prescale_HLT_Photon50_CaloIdVL_IsoL = -10;
  prescale_HLT_Photon75_CaloIdVL = -10;
  prescale_HLT_Photon90_CaloIdVL = -10;
  prescale_HLT_Photon125 = -10;
  prescale_HLT_Photon135 = -10;
  prescale_HLT_Photon200_NoHE = -10;
  prescale_HLT_Photon225_NoHE = -10;
  prescale_HLT_Photon26_Photon18 = -10;
  prescale_HLT_Photon36_Photon22 = -10;
  prescale_HLT_DoublePhoton33 = -10;
  prescale_HLT_DoublePhoton60 = -10;
  prescale_HLT_DoublePhoton70 = -10;
  prescale_HLT_DoublePhoton80 = -10;


  pthat = -5000.;
  alphaqcd = -5000.;
  alphaqed = -5000.;
  qscale = -5000.;
  processid = -5000;
  weight = -5000.;

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

  //Supercluster variables
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
  scsize = -5000;
  
  //beam spot variables
  sigmaZ = -5000.;
  sigmaZ0Error = -5000.;
  sq = -5000.; 
  bsposx = -5000.;
  bsposy = -5000.;
  bsposz = -5000.;
   
  genparticles_size =-10;
  for (int i=0; i<20; ++i) {
    genele_e[i] =  -5000.;
    genele_pt[i] =  -5000.;
    genele_px[i] = -5000.;
    genele_py[i] = -5000.;
    genele_pz[i] = -5000.;
    genele_eta[i] = -5000.;
    genele_phi[i] = -5000.;
    genele_charge[i]= -5000;
    
    unstableGenEle_e[i] =  -5000.;
    unstableGenEle_pt[i] = -5000.;
    unstableGenEle_px[i] = -5000.;
    unstableGenEle_py[i] = -5000.;
    unstableGenEle_pz[i] = -5000.;
    unstableGenEle_eta[i] = -5000.;
    unstableGenEle_phi[i] =  -5000.;
    unstableGenEle_charge[i]=  -5000;
    
    genelemom_e[i] =  -5000.;
    genelemom_pt[i] = -5000.;
    genelemom_px[i] =  -5000.;
    genelemom_py[i] =  -5000.;
    genelemom_pz[i] =  -5000.;
    genelemom_eta[i] =  -5000.;
    genelemom_phi[i] =  -5000.;
    genelemom_charge[i]=  -5000;
    genelemom_mass[i]= -5000.;
    genelemom_pdgid[i]=  -5000;
  }

  for (int i=0; i<10; ++i) {
  x1quark[i] = -5000.;
  x2quark[i] = -5000.;
  }

  trueNVtx = -5000.;
  nVtxBefore = -5000;
  nVtxNow = -5000;
  nVtxAfter = -5000;
 
  //Primary vertex x,y,z
  pvsize = -5000;
  for (int i = 0; i < 50; ++i) {
    pvx[i] = -5000.;
    pvy[i] = -5000.;
    pvz[i] = -5000.;
    pv_isValid[i] = false;
    pv_ndof[i] = -5000.;
    pv_nTracks[i] = -5000;
    pv_normChi2[i] = -5000.;
    pv_totTrackSize[i] = -5000;
  }

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
    gsf_hdepth1overe[i]=-1000.;
    gsf_hdepth2overe[i]=-1000.;
    gsf_trackiso[i]=-1000.;
    gsf_ecaliso[i]=-1000.;
    gsf_hcaliso1[i]=-1000.;
    gsf_hcaliso2[i]=-1000.;
    gsf_charge[i]=-1000;
    gsf_sigmaetaeta[i]=-1000.;
    gsf_sigmaIetaIeta[i]=-1000.;
    gsf_isecaldriven[i]=-1000;
    gsf_istrackerdriven[i]=-1000;
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
    gsf_dz[i] = -1.;
    gsf_vz[i] = -1.;
    gsf_nHits[i] = -1;
    gsf_convFlags[i]= -1;
    gsf_convDist[i]= -1;
    gsf_convDcot[i]= -1;
    gsf_convRadius[i] = -1;
    gsf_nLostInnerHits[i] = -1;
    gsf_nLostOuterHits[i] = -1;
    gsf_fBrem[i] = -1.;
    gsf_eMax[i] = -1.;
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

  for (int i=0; i<500; ++i) HLTriggers[i]  = -10;

  calomet = -1.;
  calomet_eta = -1000; 
  calomet_phi = -1000;
  met = -1.;
  pfmet = -1;
  pfmet_eta = -1000; 
  pfmet_phi = -1000;

  //IC5
  //  nJetsIC5_pt15 = -1;
//   for (unsigned int i = 0 ; i< 100 ; i++){   
//     jetIC5_pt[i] = -1;
//     jetIC5_eta[i] = -1;
//     jetIC5_phi[i] = -1;
//     jetIC5_em[i] = -1;
//   }
  
  nJetsAKT_pt15 = -1;
  for (unsigned int i = 0 ; i< 50 ; ++i){   
    jetAKT_pt[i] = -1;
    jetAKT_eta[i] = -1;
    jetAKT_phi[i] = -1;
    jetAKT_em[i] = -1;
  }

  //MUONS
  muon_size = -3;
  for (unsigned int i = 0 ; i< 100 ; ++i){  
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
    muon_nlayerswithhits[i] = -1;
    muon_nlosthits[i] = -1;
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

  // generator information for MC samples
  if (useGenData_) {
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
    genparticles_size = 1;
  }


  edm::Handle<TriggerResults> hltResults;
  iEvent.getByLabel(hlTriggerResults_, hltResults);

  edm::Handle<edm::TriggerResults> hltTriggerResultHandle;
  iEvent.getByLabel(hlTriggerResults_, hltTriggerResultHandle);
  
  hltCount = 0;   
  if(!hltTriggerResultHandle.isValid()) {
    std::cout << "invalid handle for HLT TriggerResults" << std::endl;
  } 
  else {
    hltCount = hltTriggerResultHandle->size();
    for(int i = 0 ; i < hltCount ; i++) {
      if (hltTriggerResultHandle->accept(i)) HLTriggers[i] = i;
    }
  } // end HLT
  
  // L1 BITS
  edm::Handle< L1GlobalTriggerReadoutRecord > gtReadoutRecord;
  iEvent.getByLabel( edm::InputTag("gtDigis"), gtReadoutRecord);

  const TechnicalTriggerWord&  technicalTriggerWordBeforeMask = gtReadoutRecord->technicalTriggerWord();
  L1trigger_size = technicalTriggerWordBeforeMask.size();

  for(unsigned int i = 0;i<technicalTriggerWordBeforeMask.size();i++){
    bool bit = technicalTriggerWordBeforeMask.at(i);
    if (bit == 1) L1trigger_bool[i] = 1;   
    if (bit == 0) L1trigger_bool[i] = 0;   
  }

  //physics declared
  L1GlobalTriggerReadoutRecord const* gtrr = gtReadoutRecord.product();
  L1GtFdlWord fdlWord = gtrr->gtFdlWord();
  if (fdlWord.physicsDeclared() == 1) PhysDecl_bool=1;
  else PhysDecl_bool=0;

  bool caloantiktjetisvalid = false;
  edm::Handle<CaloJetCollection> pCaloAntiKtJets;
  caloantiktjetisvalid = iEvent.getByLabel("ak5CaloJets", pCaloAntiKtJets);//Laurent
     const CaloJetCollection *caloAntiKtJets  = pCaloAntiKtJets.product();//Laurent
  
  edm::Handle<CaloMETCollection> pCaloMET;
  bool calometisvalid = iEvent.getByLabel("met", pCaloMET);
  const CaloMETCollection *caloMET  = pCaloMET.product();

  edm::Handle<METCollection> pMET;
  bool metisvalid = iEvent.getByLabel("htMetKT4", pMET);
  const METCollection *MET  = pMET.product();

  edm::Handle<PFMETCollection> pPFMET;
  bool pfmetisvalid = iEvent.getByLabel("pfMet", pPFMET);
  const PFMETCollection *PFMET  = pPFMET.product();

  // Triggers  ARNAUD
  // get hold of TriggerResults
  Handle<TriggerResults> HLTR;
  iEvent.getByLabel(hlTriggerResults_,HLTR);
  if (HLTR.isValid()) {
    if (HLTR->wasrun()) nWasRun_++;
    const bool accept(HLTR->accept());
    LogDebug("HLTrigReport") << "HL TriggerResults decision: " << accept;
    if (accept) ++nAccept_;
    if (HLTR->error() ) nErrors_++;
  } else {
    LogDebug("HLTrigReport") << "HL TriggerResults with label ["+hlTriggerResults_.encode()+"] not found!";
    nErrors_++;
  }

  // decision for each HL algorithm
  const unsigned int n(hlNames_.size());
  for (unsigned int i=0; i!=n; ++i) {

    if (hlNames_.at(i).find("HLT_Mu15_v") == 0) {
      HLTR->accept(i) ? HLT_Mu15 = 1 : HLT_Mu15 = 0;
      prescale_HLT_Mu15 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Mu30_v") == 0) {
      HLTR->accept(i) ? HLT_Mu30 = 1 : HLT_Mu30 = 0;
      prescale_HLT_Mu30 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Mu40_eta2p1_v") == 0) {
      HLTR->accept(i) ? HLT_Mu40_eta2p1 = 1 : HLT_Mu40_eta2p1 = 0;
      prescale_HLT_Mu40_eta2p1 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Mu15_Photon20_CaloIdL_v") == 0) {
      HLTR->accept(i) ? HLT_Mu15_Photon20_CaloIdL = 1 : HLT_Mu15_Photon20_CaloIdL = 0;
      prescale_HLT_Mu15_Photon20_CaloIdL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v") == 0) {
      HLTR->accept(i) ? HLT_Mu8_Ele17_CaloIdT_CaloIsoVL = 1 : HLT_Mu8_Ele17_CaloIdT_CaloIsoVL = 0;
      prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v") == 0) {
      HLTR->accept(i) ? HLT_Mu17_Ele8_CaloIdT_CaloIsoVL = 1 : HLT_Mu17_Ele8_CaloIdT_CaloIsoVL = 0;
      prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Ele8_v") == 0) {
      HLTR->accept(i) ? HLT_Ele8 = 1 : HLT_Ele8 = 0;
      prescale_HLT_Ele8 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v") == 0) {
      HLTR->accept(i) ? HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT = 1 : HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT = 0;
      prescale_HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v") == 0) {
      HLTR->accept(i) ? HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = 1 : HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = 0;
      prescale_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v") == 0) {
      HLTR->accept(i) ? HLT_Ele32_CaloIdL_CaloIsoVL_SC17 = 1 : HLT_Ele32_CaloIdL_CaloIsoVL_SC17 = 0;
      prescale_HLT_Ele32_CaloIdL_CaloIsoVL_SC17 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v") == 0) {
      HLTR->accept(i) ? HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17 = 1 : HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17 = 0;
      prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_DoubleEle33_CaloIdL_v") == 0) {
      HLTR->accept(i) ? HLT_DoubleEle33_CaloIdL = 1 : HLT_DoubleEle33_CaloIdL = 0;
      prescale_HLT_DoubleEle33_CaloIdL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_DoubleEle33_CaloIdL_CaloIsoT_v") == 0) {
      HLTR->accept(i) ? HLT_DoubleEle33_CaloIdL_CaloIsoT = 1 : HLT_DoubleEle33_CaloIdL_CaloIsoT = 0;
      prescale_HLT_DoubleEle33_CaloIdL_CaloIsoT = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_DoubleEle33_CaloIdT_v") == 0) {
      HLTR->accept(i) ? HLT_DoubleEle33_CaloIdT = 1 : HLT_DoubleEle33_CaloIdT = 0;
      prescale_HLT_DoubleEle33_CaloIdT = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_DoubleEle45_CaloIdL_v") == 0) {
      HLTR->accept(i) ? HLT_DoubleEle45_CaloIdL = 1 : HLT_DoubleEle45_CaloIdL = 0;
      prescale_HLT_DoubleEle45_CaloIdL = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
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
    if (hlNames_.at(i).find("HLT_Photon125_v") == 0) {
      HLTR->accept(i) ? HLT_Photon125 = 1 : HLT_Photon125 = 0;
      prescale_HLT_Photon125 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon135_v") == 0) {
      HLTR->accept(i) ? HLT_Photon135 = 1 : HLT_Photon135 = 0;
      prescale_HLT_Photon135 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon200_NoHE_v") == 0) {
      HLTR->accept(i) ? HLT_Photon200_NoHE = 1 : HLT_Photon200_NoHE = 0;
      prescale_HLT_Photon200_NoHE = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon225_NoHE_v") == 0) {
      HLTR->accept(i) ? HLT_Photon225_NoHE = 1 : HLT_Photon225_NoHE = 0;
      prescale_HLT_Photon225_NoHE = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon26_Photon18_v") == 0) {
      HLTR->accept(i) ? HLT_Photon26_Photon18 = 1 : HLT_Photon26_Photon18 = 0;
      prescale_HLT_Photon26_Photon18 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_Photon36_Photon22_v") == 0) {
      HLTR->accept(i) ? HLT_Photon36_Photon22 = 1 : HLT_Photon36_Photon22 = 0;
      prescale_HLT_Photon36_Photon22 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_DoublePhoton33_v") == 0) {
      HLTR->accept(i) ? HLT_DoublePhoton33 = 1 : HLT_DoublePhoton33 = 0;
      prescale_HLT_DoublePhoton33 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
    }
    if (hlNames_.at(i).find("HLT_DoublePhoton60_v") == 0) {
      HLTR->accept(i) ? HLT_DoublePhoton60 = 1 : HLT_DoublePhoton60 = 0;
      prescale_HLT_DoublePhoton60 = hltConfig_.prescaleValue(iEvent, iSetup, hlNames_.at(i));
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

  for (unsigned int i=0; i!=n; ++i) {
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

  jetAKT_size = index_jetAKT;
  nJetsAKT_pt15 = FnJetsAKT_pt15;

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

  //We take only the first primary vertex, i.e. the one with the electrons
  if(pvcoll->size() > 0) {
    reco::VertexCollection::const_iterator firstpv = pvcoll->begin();
    firstpvertex.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }
  
  pvsize = pvcoll->size();
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

  int index_mu = 0;
  //LOOP OVER MUONS
  for(reco::MuonCollection::const_iterator muIt = muons->begin();muIt != muons->end(); muIt++){
    
    if (muIt->isGlobalMuon()){
  
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

  scsize = sclusters.size();

  //sort all the refSC by energy
  std::sort(refsclusters.begin(),refsclusters.end(),refScEGreater);

  for(std::vector<reco::SuperClusterRef>::const_iterator refsclustersiter = refsclusters.begin();refsclustersiter != refsclusters.end();refsclustersiter++) {

  }
 
  //sort all the SC by energy
  std::sort(sclusters.begin(),sclusters.end(),scEGreater);

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

  int v=0;
  for(GsfTrackCollection::const_iterator gsftrackiter = gsftracks->begin(); 
      gsftrackiter!=gsftracks->end();
      gsftrackiter++)
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
  edm::Handle<EcalRecHitCollection> EBReducedRecHits;
  iEvent.getByLabel("reducedEcalRecHitsEB",EBReducedRecHits);
  edm::Handle<EcalRecHitCollection> EEReducedRecHits;
  iEvent.getByLabel("reducedEcalRecHitsEE",EEReducedRecHits);
  EcalClusterLazyTools lazytool(iEvent,iSetup,InputTag("reducedEcalRecHitsEB"),InputTag("reducedEcalRecHitsEE"));

  gsf_size = gsfelectrons.size();
  int e=0;
  reco::GsfElectronCollection::const_iterator gsfiter = gsfelectrons.begin();
  for(;gsfiter!=gsfelectrons.end();gsfiter++)
    {
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
}//end of analyze method


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
  mytree->Branch("L1trigger_bool", &L1trigger_bool, "L1trigger_bool[L1trigger_size]/I");
  mytree->Branch("PhysDecl_bool", &PhysDecl_bool, "PhysDecl_bool/I");
  mytree->Branch("HLTriggers", HLTriggers, "HLTriggers[hltCount]/I");
  mytree->Branch("nWasRun_",&nWasRun_,"nWasRun_/I");
  mytree->Branch("nAccept_",&nAccept_,"nAccept_/I");
  mytree->Branch("nErrors_",&nErrors_,"nErrors_/I");
  mytree->Branch("hlWasRun_",&hlWasRun_,"hlWasRun_/I");
  mytree->Branch("hlWasRunTab",hlWasRunTab,"hlWasRunTab[200]/I");
  mytree->Branch("hlAccept_",&hlAccept_);
  mytree->Branch("hlAcceptTab",hlAcceptTab,"hlAcceptTab[300]/I");
  mytree->Branch("hlErrorTab",hlErrorTab,"hlErrorTab[200]/I");
  mytree->Branch("hlNamesTab",&hlNamesTab,"hlNamesTab/C");
  mytree->Branch("hlNames_",&hlNames_);
  mytree->Branch("HLT_Mu15",&HLT_Mu15,"HLT_Mu15/I");
  mytree->Branch("HLT_Mu30",&HLT_Mu30,"HLT_Mu30/I");
  mytree->Branch("HLT_Mu40_eta2p1",&HLT_Mu40_eta2p1,"HLT_Mu40_eta2p1/I");
  mytree->Branch("HLT_Mu15_Photon20_CaloIdL",&HLT_Mu15_Photon20_CaloIdL,"HLT_Mu15_Photon20_CaloIdL/I");
  mytree->Branch("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL",&HLT_Mu8_Ele17_CaloIdT_CaloIsoVL,"HLT_Mu8_Ele17_CaloIdT_CaloIsoVL/I");
  mytree->Branch("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL",&HLT_Mu17_Ele8_CaloIdT_CaloIsoVL,"HLT_Mu17_Ele8_CaloIdT_CaloIsoVL/I");
  mytree->Branch("HLT_Ele8", &HLT_Ele8, "HLT_Ele8/I");
  mytree->Branch("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT", &HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT/I");
  mytree->Branch("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL", &HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL/I");
  mytree->Branch("HLT_Ele32_CaloIdL_CaloIsoVL_SC17", &HLT_Ele32_CaloIdL_CaloIsoVL_SC17, "HLT_Ele32_CaloIdL_CaloIsoVL_SC17/I");
  mytree->Branch("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17", &HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17, "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17/I");
  mytree->Branch("HLT_DoubleEle33_CaloIdL",&HLT_DoubleEle33_CaloIdL,"HLT_DoubleEle33_CaloIdL/I");
  mytree->Branch("HLT_DoubleEle33_CaloIdL_CaloIsoT",&HLT_DoubleEle33_CaloIdL_CaloIsoT,"HLT_DoubleEle33_CaloIdL_CaloIsoT/I");
  mytree->Branch("HLT_DoubleEle33_CaloIdT",&HLT_DoubleEle33_CaloIdT,"HLT_DoubleEle33_CaloIdT/I");
  mytree->Branch("HLT_DoubleEle45_CaloIdL",&HLT_DoubleEle45_CaloIdL,"HLT_DoubleEle45_CaloIdL/I");
  mytree->Branch("HLT_Photon20_CaloIdVL_IsoL", &HLT_Photon20_CaloIdVL_IsoL, "HLT_Photon20_CaloIdVL_IsoL/I");
  mytree->Branch("HLT_Photon30_CaloIdVL", &HLT_Photon30_CaloIdVL, "HLT_Photon30_CaloIdVL/I");
  mytree->Branch("HLT_Photon50_CaloIdVL", &HLT_Photon50_CaloIdVL, "HLT_Photon50_CaloIdVL/I");
  mytree->Branch("HLT_Photon50_CaloIdVL_IsoL", &HLT_Photon50_CaloIdVL_IsoL, "HLT_Photon50_CaloIdVL_IsoL/I");
  mytree->Branch("HLT_Photon75_CaloIdVL", &HLT_Photon75_CaloIdVL, "HLT_Photon75_CaloIdVL/I");
  mytree->Branch("HLT_Photon90_CaloIdVL", &HLT_Photon90_CaloIdVL, "HLT_Photon90_CaloIdVL/I");
  mytree->Branch("HLT_Photon125",&HLT_Photon125,"HLT_Photon125/I");
  mytree->Branch("HLT_Photon135",&HLT_Photon135,"HLT_Photon135/I");
  mytree->Branch("HLT_Photon200_NoHE",&HLT_Photon200_NoHE,"HLT_Photon200_NoHE/I");
  mytree->Branch("HLT_Photon225_NoHE",&HLT_Photon225_NoHE,"HLT_Photon225_NoHE/I");
  mytree->Branch("HLT_Photon26_Photon18",&HLT_Photon26_Photon18,"HLT_Photon26_Photon18/I");
  mytree->Branch("HLT_Photon36_Photon22",&HLT_Photon36_Photon22,"HLT_Photon36_Photon22/I");
  mytree->Branch("HLT_DoublePhoton33", & HLT_DoublePhoton33, " HLT_DoublePhoton33/I");
  mytree->Branch("HLT_DoublePhoton60",&HLT_DoublePhoton60,"HLT_DoublePhoton60/I");
  mytree->Branch("HLT_DoublePhoton70",&HLT_DoublePhoton70,"HLT_DoublePhoton70/I");
  mytree->Branch("HLT_DoublePhoton80",&HLT_DoublePhoton80,"HLT_DoublePhoton80/I");
  mytree->Branch("prescale_HLT_Mu15",&prescale_HLT_Mu15,"prescale_HLT_Mu15/I");
  mytree->Branch("prescale_HLT_Mu30",&prescale_HLT_Mu30,"prescale_HLT_Mu30/I");
  mytree->Branch("prescale_HLT_Mu40_eta2p1",&prescale_HLT_Mu40_eta2p1,"prescale_HLT_Mu40_eta2p1/I");
  mytree->Branch("prescale_HLT_Mu15_Photon20_CaloIdL",&prescale_HLT_Mu15_Photon20_CaloIdL,"prescale_HLT_Mu15_Photon20_CaloIdL/I");
  mytree->Branch("prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL",&prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL,"prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL/I");
  mytree->Branch("prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL",&prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL,"prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL/I");
  mytree->Branch("prescale_HLT_Ele8",&prescale_HLT_Ele8,"prescale_HLT_Ele8/I");
  mytree->Branch("prescale_HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT",&prescale_HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"prescale_HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT/I");
  mytree->Branch("prescale_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL",&prescale_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"prescale_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL/I");
  mytree->Branch("prescale_HLT_Ele32_CaloIdL_CaloIsoVL_SC17",&prescale_HLT_Ele32_CaloIdL_CaloIsoVL_SC17,"prescale_HLT_Ele32_CaloIdL_CaloIsoVL_SC17/I");
  mytree->Branch("prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17",&prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17,"prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17/I");
  mytree->Branch("prescale_HLT_DoubleEle33_CaloIdL",&prescale_HLT_DoubleEle33_CaloIdL,"prescale_HLT_DoubleEle33_CaloIdL/I");
  mytree->Branch("prescale_HLT_DoubleEle33_CaloIdL_CaloIsoT",&prescale_HLT_DoubleEle33_CaloIdL_CaloIsoT,"prescale_HLT_DoubleEle33_CaloIdL_CaloIsoT/I");
  mytree->Branch("prescale_HLT_DoubleEle33_CaloIdT",&prescale_HLT_DoubleEle33_CaloIdT,"prescale_HLT_DoubleEle33_CaloIdT/I");
  mytree->Branch("prescale_HLT_DoubleEle45_CaloIdL",&prescale_HLT_DoubleEle45_CaloIdL,"prescale_HLT_DoubleEle45_CaloIdL/I");
  mytree->Branch("prescale_HLT_Photon20_CaloIdVL_IsoL",&prescale_HLT_Photon20_CaloIdVL_IsoL,"prescale_HLT_Photon20_CaloIdVL_IsoL/I");
  mytree->Branch("prescale_HLT_Photon30_CaloIdVL",&prescale_HLT_Photon30_CaloIdVL,"prescale_HLT_Photon30_CaloIdVL/I");
  mytree->Branch("prescale_HLT_Photon50_CaloIdVL",&prescale_HLT_Photon50_CaloIdVL,"prescale_HLT_Photon50_CaloIdVL/I");
  mytree->Branch("prescale_HLT_Photon50_CaloIdVL_IsoL",&prescale_HLT_Photon50_CaloIdVL_IsoL,"prescale_HLT_Photon50_CaloIdVL_IsoL/I");
  mytree->Branch("prescale_HLT_Photon75_CaloIdVL",&prescale_HLT_Photon75_CaloIdVL,"prescale_HLT_Photon75_CaloIdVL/I");
  mytree->Branch("prescale_HLT_Photon90_CaloIdVL",&prescale_HLT_Photon90_CaloIdVL,"prescale_HLT_Photon90_CaloIdVL/I");
  mytree->Branch("prescale_HLT_Photon125",&prescale_HLT_Photon125,"prescale_HLT_Photon125/I");
  mytree->Branch("prescale_HLT_Photon135",&prescale_HLT_Photon135,"prescale_HLT_Photon135/I");
  mytree->Branch("prescale_HLT_Photon200_NoHE",&prescale_HLT_Photon200_NoHE,"prescale_HLT_Photon200_NoHE/I");
  mytree->Branch("prescale_HLT_Photon225_NoHE",&prescale_HLT_Photon225_NoHE,"prescale_HLT_Photon225_NoHE/I");
  mytree->Branch("prescale_HLT_Photon26_Photon18",&prescale_HLT_Photon26_Photon18,"prescale_HLT_Photon26_Photon18/I");
  mytree->Branch("prescale_HLT_Photon36_Photon22",&prescale_HLT_Photon36_Photon22,"prescale_HLT_Photon36_Photon22/I");
  mytree->Branch("prescale_HLT_DoublePhoton33",&prescale_HLT_DoublePhoton33,"prescale_HLT_DoublePhoton33/I");
  mytree->Branch("prescale_HLT_DoublePhoton60",&prescale_HLT_DoublePhoton60,"prescale_HLT_DoublePhoton60/I");
  mytree->Branch("prescale_HLT_DoublePhoton70",&prescale_HLT_DoublePhoton70,"prescale_HLT_DoublePhoton70/I");
  mytree->Branch("prescale_HLT_DoublePhoton80",&prescale_HLT_DoublePhoton80,"prescale_HLT_DoublePhoton80/I");


  //GLOBAL PHYSICAL INFO 
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
  mytree->Branch("pvsize",&pvsize,"pvsize/I");
  mytree->Branch("pvx",&pvx,"pvx[pvsize]/F");
  mytree->Branch("pvy",&pvy,"pvy[pvsize]/F");
  mytree->Branch("pvz",&pvz,"pvz[pvsize]/F");
  mytree->Branch("pv_isValid",&pv_isValid,"pv_isValid[pvsize]/O");
  mytree->Branch("pv_ndof",&pv_ndof,"pv_ndof[pvsize]/F");
  mytree->Branch("pv_nTracks",&pv_nTracks,"pv_nTracks[pvsize]/I");
  mytree->Branch("pv_normChi2",&pv_normChi2,"pv_normChi2[pvsize]/F");
  mytree->Branch("pv_totTrackSize",&pv_totTrackSize,"pv_totTrackSize[pvsize]/I"); 

 
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
  mytree->Branch("scgsfmatched",scgsfmatched,"scgsfmatched[scsize]/F");  

 
  //GSF VARIABLES
  mytree->Branch("gsf_size",&gsf_size, "gsf_size/I");
  mytree->Branch("gsf_theta", gsf_theta, "gsf_theta[gsf_size]/F");
  mytree->Branch("gsf_isEB", gsf_isEB, "gsf_isEB[gsf_size]/I");
  mytree->Branch("gsf_isEE", gsf_isEE, "gsf_isEE[gsf_size]/I");
  mytree->Branch("gsf_deltaEtaATcalo", gsf_deltaEtaATcalo, "gsf_deltaEtaATcalo[gsf_size]/F");
  mytree->Branch("gsf_deltaPhiATcalo", gsf_deltaPhiATcalo, "gsf_deltaPhiATcalo[gsf_size]/F");
  mytree->Branch("gsf_ecalEnergy", gsf_ecalEnergy, "gsf_ecalEnergy[gsf_size]/F");
  mytree->Branch("gsf_eOVERp", gsf_eOVERp, "gsf_eOVERp[gsf_size]/F");
  mytree->Branch("gsf_dxy", gsf_dxy, "gsf_dxy[gsf_size]/F");
  mytree->Branch("gsf_dz", gsf_dz, "gsf_dz[gsf_size]/F");
  mytree->Branch("gsf_vz", gsf_vz, "gsf_vz[gsf_size]/F");
  mytree->Branch("gsf_nHits", gsf_nHits, "gsf_nHits[gsf_size]/I");
  mytree->Branch("gsf_nLostInnerHits", gsf_nLostInnerHits, "gsf_nLostInnerHits[gsf_size]/I");
  mytree->Branch("gsf_nLostOuterHits", gsf_nLostOuterHits, "gsf_nLostOuterHits[gsf_size]/I");
  mytree->Branch("gsf_convFlags",gsf_convFlags, "gsf_convFlags[gsf_size]/I");
  mytree->Branch("gsf_convDist",gsf_convDist, "gsf_convDist[gsf_size]/F");
  mytree->Branch("gsf_convDcot",gsf_convDcot, "gsf_convDcot[gsf_size]/F");
  mytree->Branch("gsf_convRadius",gsf_convRadius, "gsf_convRadius[gsf_size]/F");
  mytree->Branch("gsf_fBrem", gsf_fBrem, "gsf_fBrem[gsf_size]/F");
  mytree->Branch("gsf_e1x5", gsf_e1x5, "gsf_e1x5[gsf_size]/F");
  mytree->Branch("gsf_e2x5", gsf_e2x5, "gsf_e2x5[gsf_size]/F");
  mytree->Branch("gsf_e5x5", gsf_e5x5, "gsf_e5x5[gsf_size]/F");
  mytree->Branch("gsf_eMax", gsf_eMax, "gsf_eMax[gsf_size]/F");
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
  mytree->Branch("gsf_hdepth1overe",&gsf_hdepth1overe,"gsf_hdepth1overe[gsf_size]/F");
  mytree->Branch("gsf_hdepth2overe",&gsf_hdepth2overe,"gsf_hdepth2overe[gsf_size]/F");
  mytree->Branch("gsf_trackiso",&gsf_trackiso,"gsf_trackiso[gsf_size]/F");
  mytree->Branch("gsf_ecaliso",&gsf_ecaliso,"gsf_ecaliso[gsf_size]/F");
  mytree->Branch("gsf_hcaliso1",&gsf_hcaliso1,"gsf_hcaliso1[gsf_size]/F");
  mytree->Branch("gsf_hcaliso2",&gsf_hcaliso2,"gsf_hcaliso2[gsf_size]/F");
  mytree->Branch("gsf_charge",&gsf_charge,"gsf_charge[gsf_size]/I");
  mytree->Branch("gsf_sigmaetaeta",&gsf_sigmaetaeta,"gsf_sigmaetaeta[gsf_size]/F");
  mytree->Branch("gsf_sigmaIetaIeta",&gsf_sigmaIetaIeta,"gsf_sigmaIetaIeta[gsf_size]/F");
  mytree->Branch("gsf_isecaldriven",&gsf_isecaldriven,"gsf_isecaldriven[gsf_size]/I");
  mytree->Branch("gsf_istrackerdriven",&gsf_istrackerdriven,"gsf_istrackerdriven[gsf_size]/I");
  mytree->Branch("gsfsc_e",&gsfsc_e,"gsfsc_e[gsf_size]/F");
  mytree->Branch("gsfsc_pt",&gsfsc_pt,"gsfsc_pt[gsf_size]/F");
  mytree->Branch("gsfsc_eta",&gsfsc_eta,"gsfsc_eta[gsf_size]/F");
  mytree->Branch("gsfsc_phi",&gsfsc_phi,"gsfsc_phi[gsf_size]/F");
  mytree->Branch("gsfsc_px",&gsfsc_px,"gsfsc_px[gsf_size]/F");
  mytree->Branch("gsfsc_py",&gsfsc_py,"gsfsc_py[gsf_size]/F");
  mytree->Branch("gsfsc_pz",&gsfsc_pz,"gsfsc_pz[gsf_size]/F");
  mytree->Branch("gsf_gsfet",&gsf_gsfet,"gsf_gsfet[gsf_size]/F");
  mytree->Branch("scindexforgsf",&scindexforgsf,"scindexforgsf[gsf_size]/I");
  mytree->Branch("gsftracksize",&gsftracksize,"gsftracksize/I");
  mytree->Branch("gsftracketa",&gsftracketa,"gsftracketa[gsftracksize]/F");
  mytree->Branch("gsftrackphi",&gsftrackphi,"gsftrackphi[gsftracksize]/F");  
  mytree->Branch("gsftrackp",&gsftrackp,"gsftrackp[gsftracksize]/F");
  mytree->Branch("gsftrackpt",&gsftrackpt,"gsftrackpt[gsftracksize]/F");
  mytree->Branch("gsftrackpx",&gsftrackpx,"gsftrackpx[gsftracksize]/F");
  mytree->Branch("gsftrackpy",&gsftrackpy,"gsftrackpy[gsftracksize]/F");
  mytree->Branch("gsftrackpz",&gsftrackpz,"gsftrackpz[gsftracksize]/F");
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
  mytree->Branch("gsfpass_NOCONVERSION",&gsfpass_NOCONVERSION,"gsfpass_NOCONVERSION[gsf_size]/O");  
  mytree->Branch("gsfpass_HEEP",&gsfpass_HEEP,"gsfpass_HEEP[gsf_size]/O");  
  mytree->Branch("gsfpass_ID",&gsfpass_ID,"gsfpass_ID[gsf_size]/O");  
  mytree->Branch("gsfpass_ISO",&gsfpass_ISO,"gsfpass_ISO[gsf_size]/O");
  

  //CHARGE INFO
  mytree->Branch("scpixcharge",&scpixcharge,"scpixcharge[gsf_size]/I");
  mytree->Branch("ctfcharge",&ctfcharge,"ctfcharge[gsf_size]/I");
  mytree->Branch("gsfcharge",&gsfcharge,"gsfcharge[gsf_size]/I");
  mytree->Branch("gsfctfscpixconsistent",&gsfctfscpixconsistent,"gsfctfscpixconsistent[gsf_size]/O");
  mytree->Branch("gsfscpixconsistent",&gsfscpixconsistent,"gsfscpixconsistent[gsf_size]/O");
  mytree->Branch("gsfctfconsistent",&gsfctfconsistent,"gsfctfconsistent[gsf_size]/O");

  mytree->Branch("genparticles_size",&genparticles_size,"genparticles_size/I");
  //GEN INFO FOR ELE and POSI (after FSR)
  
  mytree->Branch("genele_e",&genele_e,"genele_e[genparticles_size]/D");
  mytree->Branch("genele_eta",&genele_eta,"genele_eta[genparticles_size]/D");
  mytree->Branch("genele_phi",&genele_phi,"genele_phi[genparticles_size]/D");
  mytree->Branch("genele_pt",&genele_pt,"genele_pt[genparticles_size]/D");
  mytree->Branch("genele_px",&genele_px,"genele_px[genparticles_size]/D");
  mytree->Branch("genele_py",&genele_py,"genele_py[genparticles_size]/D");
  mytree->Branch("genele_pz",&genele_pz,"genele_pz[genparticles_size]/D");
  mytree->Branch("genele_charge",&genele_charge,"genele_charge[genparticles_size]/I");
  
  //generated variables for the tree (before FSR)   
  mytree->Branch("unstableGenEle_e",&unstableGenEle_e,"unstableGenEle_e[genparticles_size]/D");
  mytree->Branch("unstableGenEle_eta",&unstableGenEle_eta,"unstableGenEle_eta[genparticles_size]/D");
  mytree->Branch("unstableGenEle_phi",&unstableGenEle_phi,"unstableGenEle_phi[genparticles_size]/D");
  mytree->Branch("unstableGenEle_pt",&unstableGenEle_pt,"unstableGenEle_pt[genparticles_size]/D");
  mytree->Branch("unstableGenEle_px",&unstableGenEle_px,"unstableGenEle_px[genparticles_size]/D");
  mytree->Branch("unstableGenEle_py",&unstableGenEle_py,"unstableGenEle_py[genparticles_size]/D");
  mytree->Branch("unstableGenEle_pz",&unstableGenEle_pz,"unstableGenEle_pz[genparticles_size]/D");
  mytree->Branch("unstableGenEle_charge",&unstableGenEle_charge,"unstableGenEle_charge[genparticles_size]/I");
  
  //generated variables for the tree (Z variables)
  mytree->Branch("genelemom_e",&genelemom_e,"genelemom_e[genparticles_size]/D");
  mytree->Branch("genelemom_eta",&genelemom_eta,"genelemom_eta[genparticles_size]/D");
  mytree->Branch("genelemom_phi",&genelemom_phi,"genelemom_phi[genparticles_size]/D");
  mytree->Branch("genelemom_pt",&genelemom_pt,"genelemom_pt[genparticles_size]/D");
  mytree->Branch("genelemom_px",&genelemom_px,"genelemom_px[genparticles_size]/D");
  mytree->Branch("genelemom_py",&genelemom_py,"genelemom_py[genparticles_size]/D");
  mytree->Branch("genelemom_pz",&genelemom_pz,"genelemom_pz[genparticles_size]/D");
  mytree->Branch("genelemom_charge",&genelemom_charge,"genelemom_charge[genparticles_size]/I");
  mytree->Branch("genelemom_pdgid",&genelemom_pdgid,"genelemom_pdgid[genparticles_size]/I");
  mytree->Branch("genelemom_mass",&genelemom_mass,"genelemom_mass[genparticles_size]/D");
  
  //x1 and x2
  mytree->Branch("x1quark",&x1quark,"x1quark[genparticles_size]/F");
  mytree->Branch("x2quark",&x2quark,"x2quark[genparticles_size]/F");

  mytree->Branch("trueNVtx",&trueNVtx,"trueNVtx/F");
  mytree->Branch("nVtxBefore",&nVtxBefore,"nVtxBefore/I");
  mytree->Branch("nVtxNow",&nVtxNow,"nVtxNow/I");
  mytree->Branch("nVtxAfter",&nVtxAfter,"nVtxAfter/I");
}
 

// ------------ method called once each job just after ending the event loop  ------------
void 
GsfCheckerTree::endJob() {
}

//
void GsfCheckerTree::DataGenPart(const edm::Event& e) {

  using namespace std;
  using namespace edm;
  
  unsigned int counter = 0; 
  Handle<GenParticleCollection> genParticles;
  e.getByLabel("genParticles", genParticles);
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
      
      if(fabs(mom->pdgId())!= 22 && fabs(mom->pdgId())!= 23 && fabs(mom->pdgId())!=24 && fabs(mom->pdgId())!=32 && fabs(mom->pdgId())!=33 && fabs(mom->pdgId())!=39 &&fabs(mom->pdgId())!=  13) continue;
      
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
      counter ++;      
    }
    genparticles_size=counter;
  }  
}//end of DataGenPart

