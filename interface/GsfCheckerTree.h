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
// $Id: GsfCheckerTree.h,v 1.17 2012/01/10 17:13:54 treis Exp $
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

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TrackingTools/GsfTracking/interface/GsfConstraintAtVertex.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TFile.h"
#include "TTree.h"

//
// class decleration
//
class GsfCheckerTree : public edm::EDAnalyzer {

public:
  explicit GsfCheckerTree(const edm::ParameterSet&);
  ~GsfCheckerTree();

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

  // ----------member data ---------------------------

  // config parameters -------------------------------
  edm::InputTag  hlTriggerResults_ ;
  double comEnergy_;
  double bJetPtMin_;
  // parameter for SKIMMING
  double eleEtCut_;
  double muPtCut_;
  // -------------------------------------------------

  //L1TRIGGER
  int L1trigger_size;
  int L1trigger_bool[100];

  // HLT
  int hltCount;
  int HLTriggers[500];

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
  float jetAKT_pt[50];
  float jetAKT_eta[50];
  float jetAKT_phi[50];
  float jetAKT_em[50];
   //  int nJetsIC5_pt15; 
//   int jetIC5_size;
//   float jetIC5_pt[100];
//   float jetIC5_eta[100];
//   float jetIC5_phi[100];
//   float jetIC5_em[100];

  //BTAG
  unsigned int bTagJetColl_size;
  float bTagJet_et[50];
  float bTagJet_pt[50];
  float bTagJet_eta[50];
  float bTagJet_phi[50];
  float tCHighEffBTags[50];
  float tCHighPurBTags[50];
  float jetProbBTags[50];
  float jetBProbBTags[50];
  float sSecVertHighEffBTags[50];
  float sSecVertHighPurBTags[50];
  float cSecVertBTags[50];
  float cSecVertMVABTags[50];
  float ghostTrkBTags[50];
  float softEleIP3dBTags[50];
  float softElePtBTags[50];
  float softMuBTags[50];
  float softMuIP3dBTags[50];
  float softMuPtBTags[50];
  
  // MUON
  int muon_size;
  float muon_pt[100];
  float muon_ptError[100];
  float muon_eta[100];
  float muon_etaError[100];
  float muon_phi[100];
  float muon_phiError[100];
  float muon_theta[100];
  float muon_thetaError[100]; 
  float muon_outerPt[100];
  float muon_outerEta[100];
  float muon_outerPhi[100];
  float muon_outerTheta[100];
  float muon_px[100];
  float muon_py[100];
  float muon_pz[100];
  int muon_charge[100];
  int muon_nhitspixel[100];
  int muon_nhitstrack[100];
  int muon_nhitsmuons[100];
  int muon_nhitstotal[100];
  int muon_nlayerswithhits[100];
  int muon_nlosthits[100];
  int muon_nSegmentMatch[100];
  bool muon_isTrackerMuon[100];
  float muon_chi2[100];
  int muon_ndof[100];
  float muon_normChi2[100];
  float muon_d0[100];
  float muon_d0Error[100];
  float muon_dz_cmsCenter[100];
  float muon_dz_beamSpot[100];
  float muon_dz_firstPVtx[100];
  float muon_dzError[100];
  float muon_dxy_cmsCenter[100];
  float muon_dxy_beamSpot[100];
  float muon_dxy_firstPVtx[100];
  float muon_dxyError[100]; 
  float muon_trackIso03[100]; 
  float muon_trackIso05[100]; 
  float muon_trackIso03_ptInVeto[100]; 
  float muon_trackIso05_ptInVeto[100]; 
  float muon_emIso03[100]; 
  float muon_emIso05[100]; 
  float muon_emIso03_ptInVeto[100]; 
  float muon_emIso05_ptInVeto[100]; 
  float muon_hadIso03[100]; 
  float muon_hadIso05[100]; 
  float muon_hadIso03_ptInVeto[100]; 
  float muon_hadIso05_ptInVeto[100]; 
  float muon_innerPosx[100];
  float muon_innerPosy[100];
  float muon_innerPosz[100];

  TTree* mytree;

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
  double genele_e[20];
  double genele_pt[20];
  double genele_px[20]; 
  double genele_py[20]; 
  double genele_pz[20]; 
  double genele_eta[20]; 
  double genele_phi[20];
  int genele_charge[20];
  //Generated variables (before FSR)
  double unstableGenEle_e[20];
  double unstableGenEle_pt[20];
  double unstableGenEle_px[20];
  double unstableGenEle_py[20];
  double unstableGenEle_pz[20]; 
  double unstableGenEle_eta[20];
  double unstableGenEle_phi[20]; 
  int unstableGenEle_charge[20];
  //Generated variables (Z variables)
  double genelemom_e[20]; 
  double genelemom_pt[20]; 
  double genelemom_px[20];
  double genelemom_py[20]; 
  double genelemom_pz[20]; 
  double genelemom_eta[20];  
  double genelemom_phi[20]; 
  int genelemom_charge[20];
  double genelemom_mass[20];
  int genelemom_pdgid[20];

  float x1quark[10];
  float x2quark[10];

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
  float pvx[50];
  float pvy[50];
  float pvz[50];

  bool pv_isValid[50];
  float pv_ndof[50];
  int pv_nTracks[50];
  float pv_normChi2[50];
  int pv_totTrackSize[50];

  //Supercluster variables
  float scgsfmatched[100];
  float scseedmatched[100];
  float scenergy[100];
  float sceta[100];
  float scetacorr[100];
  float sctheta[100];
  float scthetacorr[100];
  float scet[100];
  float scphi[100];
  float scpx[100];
  float scpy[100];
  float scpz[100];
  float scx[100];
  float scy[100];
  float scz[100];
  int scsize;

  int gsf_size;
  int gsf_isEB[100];
  int gsf_isEE[100];
  float gsf_px[100];
  float gsf_py[100];
  float gsf_pz[100];
  float gsf_pt[100];
  float gsf_etSC[100];
  float gsf_eta[100];
  float gsf_phi[100];
  float gsf_theta[100];
  int gsf_charge[100];
  float gsf_deltaEtaATvtx[100];
  float gsf_deltaPhiATvtx[100];
  float gsf_deltaEtaATcalo[100];
  float gsf_deltaPhiATcalo[100];
  float gsf_sigmaetaeta[100];
  float gsf_sigmaIetaIeta[100];
  float gsf_ecalEnergy[100];
  float gsf_eOVERp[100];
  float gsf_ptOVERetsc[100];
  float gsf_dxy[100];
  float gsf_dz[100];
  float gsf_vz[100];
  int gsf_nHits[100];
  int gsf_nLostInnerHits[100];
  int gsf_nLostOuterHits[100];
  int gsf_convFlags[100];
  float gsf_convDist[100];
  float gsf_convDcot[100];
  float gsf_convRadius[100];
  float gsf_fBrem[100];
  //float gsf_e1OVERe9[100];
  float gsf_e1x5[100];
  float gsf_e2x5[100];
  float gsf_e5x5[100];

  float gsf_eMax[100];
  float gsf_e1x3[100];
  float gsf_e3x1[100];
  //float gsf_e1x5[100];
  float gsf_e2x2[100];
  float gsf_e3x2[100];
  float gsf_e3x3[100];
  float gsf_e4x4[100];
  //float gsf_e5x5[100];
  float gsf_e2x5Right[100];
  float gsf_e2x5Left[100];  
  float gsf_e2x5Top[100];  
  float gsf_e2x5Bottom[100];
  float gsf_e2x5Max[100];
  float gsf_eLeft[100];
  float gsf_eRight[100];
  float gsf_eTop[100];
  float gsf_eBottom[100];
  float gsf_e2nd[100];

  int gsf_nb;
  float gsf_p[100];
  float gsf_e[100];

  float gsf_deltaeta[100];
  float gsf_deltaphi[100];
  float gsf_hovere[100];
  float gsf_hdepth1overe[100];
  float gsf_hdepth2overe[100];

  float gsf_trackiso[100];
  float gsf_ecaliso[100];
  float gsf_hcaliso1[100];
  float gsf_hcaliso2[100];

  float gsf_class[100];
  int gsf_isecaldriven[100];
  int gsf_istrackerdriven[100];


  float gsfsc_e[100];
  float gsfsc_pt[100];
  float gsfsc_eta[100];
  float gsfsc_phi[100];
  float gsfsc_px[100];
  float gsfsc_py[100];
  float gsfsc_pz[100];

  float gsf_e2x5overe5x5[100];
  float gsf_e1x5overe5x5[100];

  float gsf_gsfet[100];

  int scindexforgsf[100];

  bool gsfpass_ET[100]; 
  bool gsfpass_PT[100]; 
  bool gsfpass_DETETA[100]; 
  bool gsfpass_CRACK[100]; 
  bool gsfpass_DETAIN[100]; 
  bool gsfpass_DPHIIN[100]; 
  bool gsfpass_HADEM[100]; 
  bool gsfpass_SIGMAIETAIETA [100];
  bool gsfpass_E2X5OVER5X5[100]; 
  bool gsfpass_ISOLEMHADDEPTH1[100];
  bool gsfpass_ISOLHADDEPTH2[100]; 
  bool gsfpass_ISOLPTTRKS[100]; 
  bool gsfpass_ECALDRIVEN[100]; 
  bool gsfpass_INVALID[100];
  bool gsfpass_NOMISSINGHITS[100];
  bool gsfpass_NOCONVERSION[100];
  bool gsfpass_HEEP[100];
  bool gsfpass_ID[100];
  bool gsfpass_ISO[100];

  //charge information
  int scpixcharge[100];
  int ctfcharge[100];
  int gsfcharge[100];
  bool gsfctfscpixconsistent[100];
  bool gsfscpixconsistent[100];
  bool gsfctfconsistent[100];
  
  //Gsf Track information
  int gsftracksize;
  float gsftracketa[100];
  float gsftrackphi[100];
  float gsftrackp[100];
  float gsftrackpt[100];
  float gsftrackpx[100];
  float gsftrackpy[100];
  float gsftrackpz[100];

  unsigned int  nEvents_;           // number of events processed

  unsigned int  nWasRun_;           // # where at least one HLT was run
  unsigned int  nAccept_;           // # of accepted events
  unsigned int  nErrors_;           // # where at least one HLT had error

  std::vector<unsigned int> hlWasRun_; // # where HLT[i] was run
  int hlWasRunTab[450];
  int hlAcceptTab[450];
  int hlErrorTab[450];
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
  int HLT_Mu15;
  int HLT_Mu30;
  int HLT_Mu40_eta2p1;
  int HLT_Mu15_Photon20_CaloIdL;
  int HLT_Mu8_Ele17_CaloIdT_CaloIsoVL;
  int HLT_Mu17_Ele8_CaloIdT_CaloIsoVL;

  int HLT_Ele8;
  int HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT;
  int HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
  int HLT_Ele32_CaloIdL_CaloIsoVL_SC17;
  int HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17;
  int HLT_DoubleEle33_CaloIdL;
  int HLT_DoubleEle33_CaloIdL_CaloIsoT;
  int HLT_DoubleEle33_CaloIdT;
  int HLT_DoubleEle45_CaloIdL;

  int HLT_Photon20_CaloIdVL_IsoL;
  int HLT_Photon30_CaloIdVL;
  int HLT_Photon50_CaloIdVL;
  int HLT_Photon50_CaloIdVL_IsoL;
  int HLT_Photon75_CaloIdVL;
  int HLT_Photon90_CaloIdVL;
  int HLT_Photon125;
  int HLT_Photon135;
  int HLT_Photon200_NoHE;
  int HLT_Photon225_NoHE;
  int HLT_Photon26_Photon18;
  int HLT_Photon36_Photon22;
  int HLT_DoublePhoton33;
  int HLT_DoublePhoton60;
  int HLT_DoublePhoton70;
  int HLT_DoublePhoton80;
 
  int prescale_HLT_Mu15;
  int prescale_HLT_Mu30;
  int prescale_HLT_Mu40_eta2p1;
  int prescale_HLT_Mu15_Photon20_CaloIdL;
  int prescale_HLT_Mu8_Ele17_CaloIdT_CaloIsoVL;
  int prescale_HLT_Mu17_Ele8_CaloIdT_CaloIsoVL;

  int prescale_HLT_Ele8;
  int prescale_HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT;
  int prescale_HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
  int prescale_HLT_Ele32_CaloIdL_CaloIsoVL_SC17;
  int prescale_HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17;
  int prescale_HLT_DoubleEle33_CaloIdL;
  int prescale_HLT_DoubleEle33_CaloIdL_CaloIsoT;
  int prescale_HLT_DoubleEle33_CaloIdT;
  int prescale_HLT_DoubleEle45_CaloIdL;

  int prescale_HLT_Photon20_CaloIdVL_IsoL;
  int prescale_HLT_Photon30_CaloIdVL;
  int prescale_HLT_Photon50_CaloIdVL;
  int prescale_HLT_Photon50_CaloIdVL_IsoL;
  int prescale_HLT_Photon75_CaloIdVL;
  int prescale_HLT_Photon90_CaloIdVL;
  int prescale_HLT_Photon125;
  int prescale_HLT_Photon135;
  int prescale_HLT_Photon200_NoHE;
  int prescale_HLT_Photon225_NoHE;
  int prescale_HLT_Photon26_Photon18;
  int prescale_HLT_Photon36_Photon22;
  int prescale_HLT_DoublePhoton33;
  int prescale_HLT_DoublePhoton60;
  int prescale_HLT_DoublePhoton70;
  int prescale_HLT_DoublePhoton80;
};

//define this as a plug-in
DEFINE_FWK_MODULE(GsfCheckerTree);
