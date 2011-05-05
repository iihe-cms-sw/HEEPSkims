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
// $Id: GsfCheckerTree.h,v 1.2 2011/05/02 09:53:10 agay Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


// #include "DataFormats/CaloRecHit/interface/CaloCluster.h"
// #include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"


// #include "DataFormats/EgammaReco/interface/SuperCluster.h"
// #include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
//#include "DataFormats/EgammaReco/interface/SeedSuperClusterAssociation.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonEnergy.h"
#include "DataFormats/MuonReco/interface/MuonTime.h"
#include "DataFormats/MuonReco/interface/MuonQuality.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
//#include "DataFormats/PixelMatchTrackReco/interface/GsfTrackSeedAssociation.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

//#include "DataFormats/TrackCandidate/interface/TrackCandidateSeedAssociation.h"
//#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include "RecoCaloTools/MetaCollections/interface/CaloRecHitMetaCollections.h"


#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TFile.h"
#include "TTree.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

//--------From Claude Charlot---------------
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"
#include "TrackingTools/GsfTracking/interface/GsfConstraintAtVertex.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronCoreFwd.h"


#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaTowerIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/ElectronTkIsolation.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

#include "DataFormats/EgammaReco/interface/BasicCluster.h"


//-----------------------------------------
//Added from Vincent
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//#include "DQMServices/Core/interface/DQMStore.h"
//#include "DQMServices/Core/interface/MonitorElement.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "TTree.h"
//-----------------------------------------


//
// class decleration
//

class GsfElectronAlgo;
class MultiTrajectoryStateTransform ;
class MultiTrajectoryStateMode ;
class TFile;
class CaloCluster;
//class GsfTrack;

class GsfCheckerTree : public edm::EDAnalyzer {

public:
  explicit GsfCheckerTree(const edm::ParameterSet&);
  ~GsfCheckerTree();
  void datagenerated(const edm::Event& e);


  //typedef edm::RefToBase<CaloCluster> CaloClusterRef;

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  //GlobalVector getTSOS(const GsfTrack &t);
  void setupES(const edm::EventSetup& es);

  typedef std::list<reco::GsfElectron *> GsfElectronPtrCollection ;
  const EcalRecHitCollection * getEcalRecHitCollection( const reco::BasicCluster &cluster );
 


  // ----------member data ---------------------------

  //Added from Vincent
  //------------------------------------
  /// this will be the name of the output file 
  std::string outputFileName_;
  /// number of events for which to print the 
  /// full decay chain to the log output
  unsigned int log_;
  /// generated particle collection src
  edm::InputTag src_;
  edm::InputTag  hlTriggerResults_ ;

  /// event counter for decay chain 
  /// logging
  unsigned int evts_;

  bool ForZee;
  bool ForData;

  bool usegendata_;

  int NbGsf;

  //L1TRIGGER
  int L1trigger_size;
  int L1trigger_bool[100];

  // HLT
  int hltCount;
  int HLTriggers[300];

  int PhysDecl_bool;
  
  //GLOBAL
  int nJetsAKT_pt15;
  int nJetsIC5_pt15;
  float calomet;
  float met;
  float mass;

  //JETS
  /* int jetAKT_size; */

/*   float jetAKT_pt[50]; */
/*   float jetAKT_eta[50]; */
/*   float jetAKT_phi[50]; */
/*   float jetAKT_em[50]; */
  
  int jetIC5_size;

  float jetIC5_pt[100];
  float jetIC5_eta[100];
  float jetIC5_phi[100];
  float jetIC5_em[100];
  

  //MUON
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
  //--------------------------------


  //TFile* rootfile;
  TTree* mytree;

  int runnumber;
  int eventnumber;
  int eventcounter;

  float pthat;
  float alphaqcd;
  float alphaqed;
  float qscale;
  int processid;
  float weight;

  //Generated variables (after FSR)
  float genelec_e_var;
  float genelec_eta_var;
  float genelec_phi_var;
  float genelec_et_var;
  float genposi_e_var;
  float genposi_eta_var;
  float genposi_phi_var;
  float genposi_et_var;
  int genelec_hassc_var;
  int genposi_hassc_var;

  //Generated variables (before FSR)
  float unstablegenelec_e_var;
  float unstablegenelec_eta_var;
  float unstablegenelec_phi_var;
  float unstablegenelec_et_var;
  float unstablegenposi_e_var;
  float unstablegenposi_eta_var;
  float unstablegenposi_phi_var;
  float unstablegenposi_et_var;

  //Generated variables (Z variables)
  float genboson_m_var;
  float genboson_eta_var;
  float genboson_phi_var;
  float genboson_e_var;
  float genboson_ez_var;
  float genboson_et_var;
  float genboson_p_var;
  float genboson_pt_var;
  float genboson_pz_var;

  float x1quark;
  float x2quark;

  //genelec
  bool genelechassc;
  bool genelechasseed;
  bool genelechastrackcand;
  bool genelechasgsftrack;
  bool genelechasgsfcore;
  bool genelechasgsf;

  //genposi
  bool genposihassc;
  bool genposihasseed;
  bool genposihastrackcand;
  bool genposihasgsftrack;
  bool genposihasgsfcore;
  bool genposihasgsf;

  //FSR variables
  int fsrposiphotonsize;
  int fsrelecphotonsize;
  
  int numberfsrelec_var;
  float energyfsrelec_var;
  float energyfsrelec[10];
  float etfsrelec[10]; 
  float etafsrelec[10];
  float phifsrelec[10];
  
  int numberfsrposi_var;
  float energyfsrposi_var;
  float energyfsrposi[10];
  float etfsrposi[10]; 
  float etafsrposi[10];
  float phifsrposi[10];
  

  //Beam spot info
  float sigmaZ;
  float sigmaZ0Error;
  float sq;
  float bsposx;
  float bsposy;
  float bsposz;

  //Primary vertex x,y,z
  int pvsize;
  float pvx[20];
  float pvy[20];
  float pvz[20];

  bool pv_isValid[20];
  float pv_ndof[20];
  int pv_nTracks[20];
  float pv_normChi2[20];
  int pv_totTrackSize[20];

  //Superclusters matching gen electrons
  //In case of several SC matching a gen elec
  //we take the one with higher energy
  //i.e. sorting SC first before matching

  float scelecenergy;
  float sceleceta;
  float scelecphi;
  float scelecgsfmatched;
  float scelecseedmatched;
  float scposienergy;
  float scposieta;
  float scposiphi;
  float scposigsfmatched;
  float scposiseedmatched;



  //Supercluster variables
  //fill e,et,eta,phi,charge for every SC in the event
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
  int ntrackerseedsneg[100];
  int necalseedsneg[100];
  int necaltrackerseedsneg[100];
  int ntrackerseedspos[100];
  int necalseedspos[100];
  int necaltrackerseedspos[100];
  int scsize;

  bool schasseed[100];
  bool schastrackcand[100];
  bool schasgsftrack[100];
  bool schasgsfcore[100];
  bool schasgsf[100];

  //initialize all counters
  int numseedspersc[100];
  int numtrackcandpersc[100];
  int numgsftrackspersc[100];
  int numgsfcorepersc[100];
  int numgsfpersc[100];


  //Matching Supercluster variables
  //fill e,et,eta,phi,charge for every matched SC in the event
  float firstmatchscenergy;
  float firstmatchsceta;
  float firstmatchsctheta;
  float firstmatchscet;
  float firstmatchscphi;
  float firstmatchscpx;
  float firstmatchscpy;
  float firstmatchscpz;
  int firstmatchnseed;

  float secondmatchscenergy;
  float secondmatchsceta;
  float secondmatchsctheta;
  float secondmatchscet;
  float secondmatchscphi;
  float secondmatchscpx;
  float secondmatchscpy;
  float secondmatchscpz;
  int secondmatchnseed;

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
  float gsf_vz[100];
  int gsf_nHits[100];
  float gsf_fBrem[100];
  //float gsf_e1OVERe9[100];
  float gsf_e1x5[100];
  float gsf_e2x5[100];
  float gsf_e5x5[100];

  float gsf_eMax[100];
  float gsf_SwissCross[100];
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

  float gsf_trackiso[100];
  float gsf_ecaliso[100];
  float gsf_hcaliso1[100];
  float gsf_hcaliso2[100];

  float gsf_class[100];
  float gsf_isecaldriven[100];
  float gsf_istrackerdriven[100];


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
  int gsfcoreindexforgsf[100];

  int gsfindexforgenelec;
  int gsfindexforgenposi;

  int scindexforgenelec;
  int scindexforgenposi;

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
  
  //Seed information
  int b_seed_nb;
  int scindexforseed[1000];
  int seedisecaldriven[1000]; 
  int seedistrackerdriven[1000];
  int seedcharge[1000];
  int nseedhits[1000];
  float seedfirsthitx[1000];
  float seedfirsthity[1000];
  float seedfirsthitz[1000];
  float seedsecondhitx[1000];
  float seedsecondhity[1000];
  float seedsecondhitz[1000];
  float seedthirdhitx[1000];
  float seedthirdhity[1000];
  float seedthirdhitz[1000];
  bool seedhastrackcand[1000];
  


  //Track candidate information
  int b_trackcand_nb;
  int scindexfortrcand[1000];
  int seedindexfortrcand[1000];
  int trackcandseedcharge[1000];
  bool trcandhasgsftrack[1000];

  //Gsf Track information
  int gsftracksize;
  float gsftracketa[100];
  float gsftrackphi[100];
  float gsftrackp[100];
  float gsftrackpt[100];
  float gsftrackpx[100];
  float gsftrackpy[100];
  float gsftrackpz[100];

  int scindexforgsftrack[100];
  int trcandindexforgsftrack[100];
  bool gsftrackhasgsfcore[100];

  
  // triggers from Arnaud

  //edm::InputTag hlTriggerResults_;  // Input tag for TriggerResults

  unsigned int  nEvents_;           // number of events processed

  unsigned int  nWasRun_;           // # where at least one HLT was run
  unsigned int  nAccept_;           // # of accepted events
  unsigned int  nErrors_;           // # where at least one HLT had error

  std::vector<unsigned int> hlWasRun_; // # where HLT[i] was run
  int hlWasRunTab[300];
  int hlAcceptTab[300];
  int hlErrorTab[300];
  //TString hlNamesTab[200];
  const char* hlNamesTab;
  //std::vector<unsigned int> *phlWasRun_ = &hlWasRun_; // # where HLT[i] was run
  std::vector<unsigned int> hltL1s_;   // # of events after L1 seed
  std::vector<unsigned int> hltPre_;   // # of events after HLT prescale
  std::vector<unsigned int> hlAccept_; // # of events accepted by HLT[i]
  std::vector<unsigned int> hlErrors_; // # of events with error in HLT[i]

  std::vector<int> posL1s_;            // pos # of last L1 seed
  std::vector<int> posPre_;            // pos # of last HLT prescale
  std::vector<std::string>  hlNames_;  // name of each HLT algorithm

  HLTConfigProvider hltConfig_;        // to get configuration for L1s/Pre

  //individual triggers

  int HLT_Ele10_SW_EleId_L1R;
  int HLT_Ele10_SW_L1R;
  int HLT_Ele15_LW_L1R;
  int HLT_Ele15_SW_L1R;
  int HLT_Ele15_SW_EleId_L1R;
  int HLT_Ele15_SW_CaloEleId_L1R ;
  int HLT_Ele15_SiStrip_L1R;
  int HLT_Ele20_SW_L1R;
  int HLT_Ele20_SiStrip_L1R;
  int HLT_Ele25_SW_L1R;
  int HLT_DoubleEle4_SW_eeRes_L1R;
  int HLT_DoubleEle10_SW_L1R;
  int HLT_Mu15_Photon20_CaloIdL; //VINCENT
 
  int HLT_Photon20_CaloIdVL_IsoL_v1;
  int HLT_DoublePhoton33_vx;
  int HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2;
  int HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2 ;
  int HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2 ;



  //-------------From Claude Charlot--------------------

  edm::ESHandle<MagneticField>                theMagField;
  edm::ESHandle<CaloGeometry>                 theCaloGeom;
  edm::ESHandle<CaloTopology>                 theCaloTopo;
  edm::ESHandle<TrackerGeometry>              trackerHandle_;

  const MultiTrajectoryStateTransform *mtsTransform_;
  const MultiTrajectoryStateMode *mtsMode_;
  GsfConstraintAtVertex *constraintAtVtx_;
  const GsfPropagatorAdapter *geomPropBw_;
  const GsfPropagatorAdapter *geomPropFw_;

  // internal variables
  int subdet_; //subdetector for this cluster
  GlobalPoint sclPos_;
  GlobalVector vtxMom_;
  TrajectoryStateOnSurface innTSOS_;
  TrajectoryStateOnSurface outTSOS_;
  TrajectoryStateOnSurface vtxTSOS_;
  TrajectoryStateOnSurface sclTSOS_;
  TrajectoryStateOnSurface seedTSOS_;
  TrajectoryStateOnSurface eleTSOS_;
  TrajectoryStateOnSurface constrainedVtxTSOS_;
  
  const TrackerGeometry* tracker;

  unsigned long long cacheIDMagField_;
  unsigned long long cacheIDTDGeom_; 
  unsigned long long cacheIDGeom_;
  unsigned long long cacheIDTopo_;

  const EcalRecHitCollection *ebRecHits_;
  const EcalRecHitCollection *eeRecHits_;

};


//define this as a plug-in
DEFINE_FWK_MODULE(GsfCheckerTree);
