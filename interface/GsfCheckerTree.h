#ifndef UserCode_HEEPSkims_GsfCheckerTree_h
#define UserCode_HEEPSkims_GsfCheckerTree_h

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
// $Id: GsfCheckerTree.h,v 1.25 2012/04/21 12:37:00 lathomas Exp $
//
//

// system include files
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPEventHelper.h"
#include "SHarper/HEEPAnalyzer/interface/HEEPEvent.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h" 
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TrackingTools/GsfTracking/interface/GsfConstraintAtVertex.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TFile.h"
#include "TTree.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"

class TH1;

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}
const int NbJets = 100;
//
// class decleration
//
class GsfCheckerTree : public edm::EDAnalyzer {

private:
  heep::EventHelper evtHelper_; //this is our magic class where all the nastyness is contained
  heep::Event heepEvt_;

  //the next three variables are simply for the example analysis
  int nrPass_;
  int nrFail_;
  ElectronHcalHelper::Configuration hcalCfg;
  ElectronHcalHelper *hcalHelper;

public:
  explicit GsfCheckerTree(const edm::ParameterSet& iConfig);
  ~GsfCheckerTree();
  typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > > IsoDepositMaps;
  typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);

  void DataGenPart(const edm::Event& e);
  void L1TInfo(const edm::Event& iEvent);
  void HLTInfo(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void METData(const edm::Event& iEvent);
  void JetData(const edm::Event& iEvent);
  void BTagData(const edm::Event& event);



  std::vector<edm::InputTag> inputTagIsoDepElectrons_;
  std::vector<edm::InputTag> inputTagIsoValElectronsNoPFId_;
  std::vector<edm::InputTag> inputTagIsoValElectronsPFId_;   


  // ----------member data ---------------------------

  // config parameters -------------------------------
  edm::InputTag  hlTriggerResults_ ;
  double comEnergy_;
  double bJetPtMin_;
  // parameter for SKIMMING
  double ele1EtMin_;
  double ele1EtMax_;
  double ele2EtMin_;
  double ele2EtMax_;
  double mu1PtMin_;
  double mu1PtMax_;
  double mu2PtMin_;
  double mu2PtMax_;
  // -------------------------------------------------
  
  TTree* mytree;

  //L1TRIGGER
  int L1trigger_size;
  int *L1trigger_bool;

  // HLT
  int hltCount;
  int *HLTriggers;

  int PhysDecl_bool;
  
  //GLOBAL
  float rho;
  float rhoiso;
  float calomet;
  float calomet_eta;
  float calomet_phi;
  float met;
  float pfmet; 
  float pfmet_eta;
  float pfmet_phi;
  float mass;

  //JETS 
  int nJetsAKT_pt15;
  int jetAKT_size;
  float *jetAKT_eta;
  float *jetAKT_pt;
  float *jetAKT_phi;
  float *jetAKT_em;
   //  int nJetsIC5_pt15; 
//   int jetIC5_size;
//   float jetIC5_pt[NbJets];
//   float jetIC5_eta[NbJets];
//   float jetIC5_phi[NbJets];
//   float jetIC5_em[NbJets];

  //BTAG
  int bTagJetColl_size;
  float *bTagJet_et;
  float *bTagJet_pt;
  float *bTagJet_eta;
  float *bTagJet_phi;
  float *tCHighEffBTags;
  float *tCHighPurBTags;
  float *jetProbBTags;
  float *jetBProbBTags;
  float *sSecVertHighEffBTags;
  float *sSecVertHighPurBTags;
  float *cSecVertBTags;
  float *cSecVertMVABTags;
  float *ghostTrkBTags;
  float *softEleIP3dBTags;
  float *softElePtBTags;
  float *softMuBTags;
  float *softMuIP3dBTags;
  float *softMuPtBTags;
  
  // MUON
  int muon_size;
  float *muon_pt;
  float *muon_ptError;
  float *muon_eta;
  float *muon_etaError;
  float *muon_phi;
  float *muon_phiError;
  float *muon_theta;
  float *muon_thetaError; 
  float *muon_outerPt;
  float *muon_outerEta;
  float *muon_outerPhi;
  float *muon_outerTheta;
  float *muon_px;
  float *muon_py;
  float *muon_pz;
  int *muon_charge;
  int *muon_nhitspixel;
  int *muon_nhitstrack;
  int *muon_nhitsmuons;
  int *muon_nhitstotal;
  int *muon_nlayerswithhits;
  int *muon_nlosthits;
  int *muon_nSegmentMatch;
  bool *muon_isTrackerMuon;
  float *muon_chi2;
  int *muon_ndof;
  float *muon_normChi2;
  float *muon_d0;
  float *muon_d0Error;
  float *muon_dz_cmsCenter;
  float *muon_dz_beamSpot;
  float *muon_dz_firstPVtx;
  float *muon_dzError;
  float *muon_dxy_cmsCenter;
  float *muon_dxy_beamSpot;
  float *muon_dxy_firstPVtx;
  float *muon_dxyError; 
  float *muon_trackIso03; 
  float *muon_trackIso05; 
  float *muon_trackIso03_ptInVeto; 
  float *muon_trackIso05_ptInVeto; 
  float *muon_emIso03; 
  float *muon_emIso05; 
  float *muon_emIso03_ptInVeto; 
  float *muon_emIso05_ptInVeto; 
  float *muon_hadIso03; 
  float *muon_hadIso05; 
  float *muon_hadIso03_ptInVeto; 
  float *muon_hadIso05_ptInVeto; 
  float *muon_innerPosx;
  float *muon_innerPosy;
  float *muon_innerPosz;

  unsigned int runnumber;
  unsigned int eventnumber;
  unsigned int eventcounter;
  unsigned int luminosityBlock;

  float pthat;
  float alphaqcd;
  float alphaqed;
  float qscale;
  int processid;
  float weight;

  int genparticles_size;
  //Generated variables (after FSR)
  double *genele_e;
  double *genele_pt;
  double *genele_px; 
  double *genele_py; 
  double *genele_pz; 
  double *genele_eta; 
  double *genele_phi;
  int *genele_charge;
  //Generated variables (before FSR)
  double *unstableGenEle_e;
  double *unstableGenEle_pt;
  double *unstableGenEle_px;
  double *unstableGenEle_py;
  double *unstableGenEle_pz; 
  double *unstableGenEle_eta;
  double *unstableGenEle_phi; 
  int *unstableGenEle_charge;
  //Generated variables (Z variables)
  double *genelemom_e; 
  double *genelemom_pt; 
  double *genelemom_px;
  double *genelemom_py; 
  double *genelemom_pz; 
  double *genelemom_eta;  
  double *genelemom_phi; 
  int *genelemom_charge;
  double *genelemom_mass;
  int *genelemom_pdgid;

  float *x1quark;
  float *x2quark;

  float trueNVtx;
  int nVtxBefore;
  int nVtxNow;
  int nVtxAfter;

  //Beam spot info
  float sigmaZ;
  float sigmaZ0Error;
  float sq;
  float bsposx;
  float bsposy;
  float bsposz;

  //Primary vertex x,y,z
  int pvsize;
  float *pvx;
  float *pvy;
  float *pvz;
  bool *pv_isValid;
  float *pv_ndof;
  int *pv_nTracks;
  float *pv_normChi2;
  int *pv_totTrackSize;

  //Supercluster variables
  int scsize;
  //float *scgsfmatched;
  //float *scseedmatched;
  float *scenergy;
  float *sceta;
  float *scetacorr;
  float *sctheta;
  float *scthetacorr;
  float *scet;
  float *scphi;
  float *scpx;
  float *scpy;
  float *scpz;
  float *scx;
  float *scy;
  float *scz;

  int gsf_size;
  int *gsf_isEB;
  int *gsf_isEE;
  float *gsf_px;
  float *gsf_py;
  float *gsf_pz;
  float *gsf_pt;
  //float *gsf_etSC;
  float *gsf_eta;
  float *gsf_phi;
  float *gsf_theta;
  int *gsf_charge;
  //float *gsf_deltaEtaATvtx;
  //float *gsf_deltaPhiATvtx;
  float *gsf_deltaEtaATcalo;
  float *gsf_deltaPhiATcalo;
  float *gsf_sigmaetaeta;
  float *gsf_sigmaIetaIeta;
  float *gsf_ecalEnergy;
  float *gsf_eOVERp;
  //float *gsf_ptOVERetsc;
  float *gsf_dxy;
  float *gsf_dz;
  float *gsf_vz;
  int *gsf_nHits;
  int *gsf_nLostInnerHits;
  int *gsf_nLostOuterHits;
  int *gsf_convFlags;
  float *gsf_convDist;
  float *gsf_convDcot;
  float *gsf_convRadius;
  float *gsf_fBrem;
  //float *gsf_e1OVERe9;
  float *gsf_e1x5;
  float *gsf_e2x5;
  float *gsf_e5x5;
  //float *gsf_eMax;
  float *gsf_e1x3;
  //float *gsf_e3x1;
  //float *gsf_e1x5;
  //float *gsf_e2x2;
  //float *gsf_e3x2;
  //float *gsf_e3x3;
  //float *gsf_e4x4;
  //float *gsf_e5x5;
  //float *gsf_e2x5Right;
  //float *gsf_e2x5Left;  
  //float *gsf_e2x5Top;  
  //float *gsf_e2x5Bottom;
  //float *gsf_e2x5Max;
  //float *gsf_eLeft;
  //float *gsf_eRight;
  //float *gsf_eTop;
  //float *gsf_eBottom;
  //float *gsf_e2nd;
  float *gsf_p;
  float *gsf_e;
  float *gsf_deltaeta;
  float *gsf_deltaphi;
  float *gsf_hovere;
  float *gsf_hdepth1overe;
  float *gsf_hdepth2overe;
  float *gsf_hovere2012;
  float *gsf_hdepth1overe2012;
  float *gsf_hdepth2overe2012;
  float *gsf_PFisocharged; 
  float *gsf_PFisophoton;
  float *gsf_PFisoneutral;
  float *gsf_trackiso;
  float *gsf_ecaliso;
  float *gsf_hcaliso1;
  float *gsf_hcaliso2;
  float *gsf_hcaliso12012;
  float *gsf_hcaliso22012;
  float *gsf_class;
  int *gsf_isecaldriven;
  int *gsf_istrackerdriven;
  float *gsfsc_e;
  float *gsfsc_pt;
  float *gsfsc_eta;
  float *gsfsc_phi;
  float *gsfsc_px;
  float *gsfsc_py;
  float *gsfsc_pz;
  float *gsf_e2x5overe5x5;
  float *gsf_e1x5overe5x5;
  float *gsf_gsfet;
  int *scindexforgsf;
  bool *gsfpass_ET; 
  bool *gsfpass_PT; 
  bool *gsfpass_DETETA; 
  bool *gsfpass_CRACK; 
  bool *gsfpass_DETAIN; 
  bool *gsfpass_DPHIIN; 
  bool *gsfpass_HADEM; 
  bool *gsfpass_SIGMAIETAIETA ;
  bool *gsfpass_E2X5OVER5X5; 
  bool *gsfpass_ISOLEMHADDEPTH1;
  bool *gsfpass_ISOLHADDEPTH2; 
  bool *gsfpass_ISOLPTTRKS; 
  bool *gsfpass_ECALDRIVEN; 
  bool *gsfpass_INVALID;
  bool *gsfpass_NOMISSINGHITS;
  bool *gsfpass_NOCONVERSION;
  bool *gsfpass_HEEP;
  bool *gsfpass_ID;
  bool *gsfpass_ISO;

  //charge information
  int *scpixcharge;
  int *ctfcharge;
  int *gsfcharge;
  bool *gsfctfscpixconsistent;
  bool *gsfscpixconsistent;
  bool *gsfctfconsistent;
  
  //Gsf Track information
  int gsftracksize;
  float *gsftracketa;
  float *gsftrackphi;
  float *gsftrackp;
  float *gsftrackpt;
  float *gsftrackpx;
  float *gsftrackpy;
  float *gsftrackpz;

  unsigned int  nEvents_;           // number of events processed

  unsigned int  nWasRun_;           // # where at least one HLT was run
  unsigned int  nAccept_;           // # of accepted events
  unsigned int  nErrors_;           // # where at least one HLT had error

  std::vector<unsigned int> hlWasRun_; // # where HLT[i] was run
  int hlWasRunTab[500];
  int hlAcceptTab[500];
  int hlErrorTab[500];
  const char* hlNamesTab;
  std::vector<unsigned int> hltL1s_;   // # of events after L1 seed
  std::vector<unsigned int> hltPre_;   // # of events after HLT prescale
  std::vector<unsigned int> hlAccept_; // # of events accepted by HLT[i]
  std::vector<unsigned int> hlErrors_; // # of events with error in HLT[i]

  std::vector<int> posL1s_;            // pos # of last L1 seed
  std::vector<int> posPre_;            // pos # of last HLT prescale
  std::vector<std::string>  hlNames_;  // name of each HLT algorithm

  HLTConfigProvider hltConfig_;        // to get configuration for L1s/Pre

  //individual triggers
  int HLT_Mu15_eta2p1;
  int HLT_Mu24_eta2p1;
  int HLT_Mu30_eta2p1;
  int HLT_Mu40_eta2p1;
  int HLT_Mu50_eta2p1;
  int HLT_Mu22_TkMu22;
  int HLT_Mu22_Photon22_CaloIdL;
  int HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
  int HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
  int HLT_Ele8_CaloIdL_CaloIsoVL;
  int HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL;
  int HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
  int HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50;
  int HLT_DoubleEle33_CaloIdL;
  int HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
  int HLT_DoubleEle33_CaloIdT;
  int HLT_Photon20_CaloIdVL_IsoL;
  int HLT_Photon30_CaloIdVL;
  int HLT_Photon50_CaloIdVL;
  int HLT_Photon50_CaloIdVL_IsoL;
  int HLT_Photon75_CaloIdVL;
  int HLT_Photon90_CaloIdVL;
  int HLT_Photon135;
  int HLT_Photon150;
  int HLT_Photon250_NoHE;
  int HLT_Photon300_NoHE;
  int HLT_Photon26_Photon18;
  int HLT_Photon36_Photon22;
  int HLT_DoublePhoton70;
  int HLT_DoublePhoton80;

  int prescale_HLT_Mu15_eta2p1;
  int prescale_HLT_Mu24_eta2p1;
  int prescale_HLT_Mu30_eta2p1;
  int prescale_HLT_Mu40_eta2p1;
  int prescale_HLT_Mu50_eta2p1;
  int prescale_HLT_Mu22_TkMu22;
  int prescale_HLT_Mu22_Photon22_CaloIdL;
  int prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
  int prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL;
  int prescale_HLT_Ele8_CaloIdL_CaloIsoVL;
  int prescale_HLT_Ele27_CaloIdL_CaloIsoVL_TrkIdVL_TrkIsoVL;
  int prescale_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
  int prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50;
  int prescale_HLT_DoubleEle33_CaloIdL;
  int prescale_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL;
  int prescale_HLT_DoubleEle33_CaloIdT;
  int prescale_HLT_Photon20_CaloIdVL_IsoL;
  int prescale_HLT_Photon30_CaloIdVL;
  int prescale_HLT_Photon50_CaloIdVL;
  int prescale_HLT_Photon50_CaloIdVL_IsoL;
  int prescale_HLT_Photon75_CaloIdVL;
  int prescale_HLT_Photon90_CaloIdVL;
  int prescale_HLT_Photon135;
  int prescale_HLT_Photon150;
  int prescale_HLT_Photon250_NoHE;
  int prescale_HLT_Photon300_NoHE;
  int prescale_HLT_Photon26_Photon18;
  int prescale_HLT_Photon36_Photon22;
  int prescale_HLT_DoublePhoton70;
  int prescale_HLT_DoublePhoton80;

};
#endif
//define this as a plug-in

