#define DataMCTreeAnalysis_cxx
#include "DataMCTreeAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>
#include <stdio.h>
#include <sstream>

#include "TLorentzVector.h"
#include "TPaveLabel.h"
#include "TString.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH1I.h"

//#include "/beo5/charafbis/HomeSafe03062010/CMSSW_2_2_6/src/HEEPRoutine/plostyle.C"

#define PI 3.141592654
#define TWOPI 6.283185308


bool saveplots = true;
//bool saveplots = false;

void DataMCTreeAnalysis::Loop()
{

  std::map< TString, int > histoNBins;
  std::map< TString, Double_t > histoXMin;
  std::map< TString, Double_t > histoXMax;
  std::map< TString, TString> histoName;
  std::map< TString, bool> histoLogScale;
  
  histoNBins["gsf_deltaeta"]=50;histoXMin["gsf_deltaeta"]=-0.02;histoXMax["gsf_deltaeta"]=0.02;
  histoNBins["gsf_deltaphi"]=50;histoXMin["gsf_deltaphi"]=-0.2;histoXMax["gsf_deltaphi"]=0.2;
  histoNBins["gsf_hovere"]=50;histoXMin["gsf_hovere"]=0.;histoXMax["gsf_hovere"]=0.1;
  histoNBins["gsf_trackiso"]=40;histoXMin["gsf_trackiso"]=0.;histoXMax["gsf_trackiso"]=20.;
  histoNBins["gsf_ecaliso"]=40;histoXMin["gsf_ecaliso"]=-10.;histoXMax["gsf_ecaliso"]=10.;
  histoNBins["gsf_hcaliso1"]=20;histoXMin["gsf_hcaliso1"]=0.;histoXMax["gsf_hcaliso1"]=5.;
  histoNBins["gsf_hcaliso2"]=20;histoXMin["gsf_hcaliso2"]=0.;histoXMax["gsf_hcaliso2"]=5.;
  histoNBins["gsf_sigmaetaeta"]=50;histoXMin["gsf_sigmaetaeta"]=0.;histoXMax["gsf_sigmaetaeta"]=0.05;
  histoNBins["gsf_sigmaIetaIeta"]=50;histoXMin["gsf_sigmaIetaIeta"]=0.;histoXMax["gsf_sigmaIetaIeta"]=0.05;
  histoNBins["gsf_e2x5overe5x5"]=50;histoXMin["gsf_e2x5overe5x5"]=0.;histoXMax["gsf_e2x5overe5x5"]=1.;

  histoNBins["gsf_gsfet"]=100;histoXMin["gsf_gsfet"]=0.;histoXMax["gsf_gsfet"]=200.;	
  histoNBins["gsf_pt"]=100;histoXMin["gsf_pt"]=0.;histoXMax["gsf_pt"]=200.;	
  histoNBins["gsf_eta"]=60;histoXMin["gsf_eta"]=-3.;histoXMax["gsf_eta"]=3.;	
  histoNBins["gsf_phi"]=36;histoXMin["gsf_phi"]=-3.1415927;histoXMax["gsf_phi"]=3.1415927;	
  
  histoNBins["gsfgsf_et"]=100;histoXMin["gsfgsf_et"]=0.;histoXMax["gsfgsf_et"]=200.;	
  histoNBins["gsfgsf_pt"]=100;histoXMin["gsfgsf_pt"]=0.;histoXMax["gsfgsf_pt"]=200.;	
  histoNBins["gsfgsf_eta"]=60;histoXMin["gsfgsf_eta"]=-3.;histoXMax["gsfgsf_eta"]=3.;	
  histoNBins["gsfgsf_phi"]=36;histoXMin["gsfgsf_phi"]=-3.1415927;histoXMax["gsfgsf_phi"]=3.1415927;	
  
  //   histoXTitle["gsf_deltaeta"]="#Delta#eta";
  //   histoXTitle["gsf_deltaphi"]="#Delta#phi";
  //   histoXTitle["gsf_hovere"]="H/E";
  //   histoXTitle["gsf_trackiso"]="Track Iso.";
  //   histoXTitle["gsf_ecaliso"]="Ecal Iso.";
  //   histoXTitle["gsf_hcaliso1"]="Hcal Iso. 1";
  //   histoXTitle["gsf_hcaliso2"]="Hcal Iso. 2";
  //   histoXTitle["gsf_sigmaetaeta"]="#sigma_{#eta#eta}";
  //   histoXTitle["gsf_sigmaIetaIeta"]="#sigma_{i#etai#eta}";
  //   histoXTitle["gsf_gsfet"]="E_{t}";
  //   histoXTitle["gsf_pt"]="p_{t}";
  //   histoXTitle["gsf_eta"]="#eta";
  //   histoXTitle["gsf_phi"]="#phi";
  
  //   histoXTitle["gsfgsf_et"]="E_{t}";
  //   histoXTitle["gsfgsf_pt"]="p_{t}";
  //   histoXTitle["gsfgsf_eta"]="#eta";
  //   histoXTitle["gsfgsf_phi"]="#phi";
  
  //float LumiFactor = 35.4; //Lumi in pb-1
  //float LumiFactor = 35.4*(8734./9262.9462); //Lumi in pb-1 (normalize to Z peak 60-120)
  //float LumiFactor = 10.0*(8734./9262.9462); //Lumi in pb-1 (normalize to Z peak 60-120)
  //float LumiFactor = 30.8*(8734./9262.9462); //Lumi in pb-1 (normalize to Z peak 60-120)


  // Mine
  // Electron-Run2010B-PromptReco Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_v3
  //float LumiFactor = 30.782997686*(8734./9262.9462); //Lumi in pb-1 (normalize to Z peak 60-120)
  //float LumiFactor = 35.3*(8734./9262.9462); //Lumi in pb-1 (normalize to Z peak 60-120)
  // Electron-Run2010B-PromptReco EgamGOODNotInOffJson
  //float LumiFactor = 9.64*(8734./9262.9462); //Lumi in pb-1 (normalize to Z peak 60-120)
  // Electron-Run2010B-PromptReco QFLAGSEgamGOOD_DCSNone
  //float LumiFactor = 40.4*(8734./9262.9462); //Lumi in pb-1 (normalize to Z peak 60-120)
  // data 2011
  float LumiFactor = 19.7*(8734./9262.9462); //Lumi in pb-1 (normalize to Z peak 60-120)


  // INPUT FILES
  vector<TFile *> input;
  vector <float> weight;
  vector <float> gridjobeffi;

  vector<TString> legendname;


  // official streamxpress v3
  //input.push_back(new TFile("Electron-Run2010B-PromptReco-v2-RECO-HEEPSkimTwoGSFEleEt20HoE10-Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_v3.root","open"));

  // EgamGOODNotInOffJson
  //input.push_back(new TFile("Electron-Run2010B-PromptReco-v2-RECO-HEEPSkimTwoGSFEleEt20HoE10-Cert_132440-149442_EgamGOODNotInOffJson.txt/res/Electron-Run2010B-PromptReco-v2-RECO-HEEPSkimTwoGSFEleEt20HoE10-Cert_132440-149442_EgamGOODNotInOffJson.root","open"));  

  // QFLAGSEgamGOOD_DCSNone
  //input.push_back(new TFile("Electron-Run2010B-PromptReco-v2-RECO-HEEPSkimTwoGSFEleEt20HoE10-Cert-QFLAGSEgamGOOD_DCSNone_StreamExpress.root","open"));  


  // 2011 data
  //input.push_back(new TFile("Photon-Run2011A-PromptReco-v1-AOD_DCSTRONLY_160433-161312_gsfcheckertree/res/Photon-Run2011A-PromptReco-v1-AOD_DCSTRONLY_160433-161312_gsfcheckertree.root","open"));  
  input.push_back(new TFile("Photon-Run2011A-PromptReco-v1-AOD_DCSTRONLY_160433-161312_gsfcheckertree_sub7/res/Photon-Run2011A-PromptReco-v1-AOD_DCSTRONLY_160433-161312_gsfcheckertree_sub7.root","open"));  

  weight.push_back(1.); //DATA
  gridjobeffi.push_back(1.); //DATA
  legendname.push_back(TString("Data"));
  
  //   input.push_back(new TFile("/user/charaf/CMSSW_3_6_1_patch3/src/UserCode/OCharaf/test/PNFSRootFiles/FullTrees/DATA_0_144114_2773nb_JSON_TwoEle.root","open"));  
  //   weight.push_back(1.); //DATA
  //   gridjobeffi.push_back(1.); //DATA
  //   legendname.push_back(TString("Data"));
  
  input.push_back(new TFile("/user/charaf/CMSSW_3_6_1_patch3/src/UserCode/OCharaf/test/PNFSRootFiles/FullTrees/ZeePowHegNoSkimV6.root","open")); 
  weight.push_back(1522./937224); //ZeeNoSkimPowHeg
  gridjobeffi.push_back(1.09); //ZeeNoSkimPowHeg
  legendname.push_back(TString("Zee PowHeg"));
  
  //   input.push_back(new TFile("/user/charaf/CMSSW_3_6_1_patch3/src/UserCode/OCharaf/test/PNFSRootFiles/FullTrees/ZprimeEE_M750_NoSkim.root","open")); 
  //   weight.push_back(1.63045e-05); //ZeeNoSkim
  //   gridjobeffi.push_back(1.28); //ZeeNoSkim
  //   legendname.push_back(TString("Z' ee M = 750 GeV/c^{2}"));
  
  //   input.push_back(new TFile("/user/charaf/CMSSW_3_6_1_patch3/src/UserCode/OCharaf/test/PNFSRootFiles/FullTrees/ZeeNoSkimV6NEW.root","open")); 
  //   weight.push_back(0.00052811); //ZeeNoSkim
  //   gridjobeffi.push_back(1.28); //ZeeNoSkim
  //   legendname.push_back(TString("Zee"));
  
  input.push_back(new TFile("/user/vdero/ProdTreeSummer2010/CMSSW_3_5_8/src/UserCode/OCharaf/test/Sample13July/FullTrees/TTbarV6_OneEle.root","open"));    
  weight.push_back(0.000064042); //TTbar Madgraph
  gridjobeffi.push_back(1.66*1.1135); //TTbar Madgraph
  legendname.push_back(TString("TTbar+jets"));
  
  //   input.push_back(new TFile("/user/charaf/CMSSW_3_6_1_patch3/src/UserCode/OCharaf/test/PNFSRootFiles/FullTrees/TTbarPythiaV6_OneEle.root","open"));    
  //   weight.push_back(2.6107e-04); //TTbar Pythia
  //   gridjobeffi.push_back(1.66); //TTbar Pythia
  //   legendname.push_back(TString("TTbar"));
  

  //   input.push_back(new TFile("/user/vdero/ProdTreeSummer2010/CMSSW_3_5_8/src/UserCode/OCharaf/test/Sample13July/FullTrees/QCD15V6_OneEle.root","open")); 
  //   weight.push_back(141.53); //QCD15
  //   gridjobeffi.push_back(1.); //QCD15
  //   legendname.push_back(TString("QCD 15-30"));
  
  //   input.push_back(new TFile("/user/vdero/ProdTreeSummer2010/CMSSW_3_5_8/src/UserCode/OCharaf/test/Sample13July/FullTrees/QCD30V6_OneEle.root","open")); 
  //   weight.push_back(11.463);//QCD30
  //   gridjobeffi.push_back(1.);//QCD30
  //   legendname.push_back(TString("QCD 30-80"));
  
  //   input.push_back(new TFile("/user/vdero/ProdTreeSummer2010/CMSSW_3_5_8/src/UserCode/OCharaf/test/Sample13July/FullTrees/QCD80V6_OneEle.root","open")); 
  //   weight.push_back(0.28673); //QCD80
  //   gridjobeffi.push_back(1.0942); //QCD80
  //   legendname.push_back(TString("QCD 80-170"));
  
  //   input.push_back(new TFile("/user/vdero/ProdTreeSummer2010/CMSSW_3_5_8/src/UserCode/OCharaf/test/Sample13July/FullTrees/QCD170V6_OneEle.root","open")); 
  //   weight.push_back(0.0080297);//QCD170
  //   gridjobeffi.push_back(1.3312);//QCD170
  //   legendname.push_back(TString("QCD 170"));


  input.push_back(new TFile("/user/charaf/CMSSW_3_6_1_patch3/src/UserCode/OCharaf/test/PNFSRootFiles/FullTrees/QCDEmEnrichedPt20to30V6_TwoEle.root","open")); 
  weight.push_back(6.924*0.0073);//QCDEm20to30
  gridjobeffi.push_back(1.0302);//QCDEm20to30
  legendname.push_back(TString("QCD Em 20-30"));
  
  input.push_back(new TFile("/user/charaf/CMSSW_3_6_1_patch3/src/UserCode/OCharaf/test/PNFSRootFiles/FullTrees/QCDEmEnrichedPt30to80V6_TwoEle.root","open"));   
  weight.push_back(1.399*0.059);//QCDEm30to80
  gridjobeffi.push_back(42377278./41529754);//QCDEm30to80
  legendname.push_back(TString("QCD Em 30-80"));
  
  input.push_back(new TFile("/user/charaf/CMSSW_3_6_1_patch3/src/UserCode/OCharaf/test/PNFSRootFiles/FullTrees/QCDEmEnrichedPt80to170V6_TwoEle.root","open")); 
  weight.push_back(0.1648*0.148);//QCDEm80to170
  gridjobeffi.push_back(1.);//QCDEm80to170
  legendname.push_back(TString("QCD Em 80-170"));
  


  input.push_back(new TFile("/user/charaf/CMSSW_3_6_1_patch3/src/UserCode/OCharaf/test/PNFSRootFiles/FullTrees/QCDBCtoEPt20to30V6_OneEle.root","open")); 
  weight.push_back(85.294*0.00046);//QCDBCtoE20to30
  gridjobeffi.push_back(1.5249);//QCDBCtoE20to30
  legendname.push_back(TString("QCD b/c 20-30"));

  input.push_back(new TFile("/user/charaf/CMSSW_3_6_1_patch3/src/UserCode/OCharaf/test/PNFSRootFiles/FullTrees/QCDBCtoEPt30to80V6_TwoEle.root","open"));   
  weight.push_back(23.953*0.00234);//QCDBCtoE30to80
  gridjobeffi.push_back(1.3792);//QCDBCtoE30to80
  legendname.push_back(TString("QCD b/c 30-80"));
  
  input.push_back(new TFile("/user/charaf/CMSSW_3_6_1_patch3/src/UserCode/OCharaf/test/PNFSRootFiles/FullTrees/QCDBCtoEPt80to170V6_TwoEle.root","open")); 
  weight.push_back(0.7495*0.0104);//QCDBCtoE80to170
  gridjobeffi.push_back(1.0031);//QCDBCtoE80to170
  legendname.push_back(TString("QCD b/c 80-170"));
  

  input.push_back(new TFile("/user/vdero/ProdTreeSummer2010/CMSSW_3_5_8/src/UserCode/OCharaf/test/Sample13July/FullTrees/ZmumuV6_OneEle.root","open")); 
  weight.push_back(6.1574e-04);//Zmumu
  gridjobeffi.push_back(1.28);//Zmumu
  legendname.push_back(TString("Z#mu#mu"));

  input.push_back(new TFile("/user/vdero/ProdTreeSummer2010/CMSSW_3_5_8/src/UserCode/OCharaf/test/Sample13July/FullTrees/ZtautauV6_OneEle.root","open")); 
  weight.push_back(5.9218e-04);//Ztautau
  gridjobeffi.push_back(1.28);//Ztautau
  legendname.push_back(TString("Z#tau#tau"));

  input.push_back(new TFile("/user/vdero/ProdTreeSummer2010/CMSSW_3_5_8/src/UserCode/OCharaf/test/SamplesJuly/FullTrees/WJetV5_OneEle.root","open"));
  weight.push_back(2.4004e-03);//WJets
  gridjobeffi.push_back(1.30*2.2593);//WJets
  legendname.push_back(TString("W+jets"));

  input.push_back(new TFile("/user/vdero/ProdTreeSummer2010/CMSSW_3_5_8/src/UserCode/OCharaf/test/Sample13July/FullTrees/WWV6_OneEle.root","open")); 
  weight.push_back(0.000027787);//WW
  gridjobeffi.push_back(1.);//WW
  legendname.push_back(TString("WW"));

  input.push_back(new TFile("/user/charaf/CMSSW_3_6_1_patch3/src/UserCode/OCharaf/test/PNFSRootFiles/FullTrees/tWV6_TwoEle.root","open")); 
  weight.push_back(0.0000000213);//tW
  gridjobeffi.push_back(1.);//tW
  legendname.push_back(TString("tW"));


  int nbFile = input.size();

  //   //WEIGHTS FOR 1pb-1


  vector<TH1F *> histopvsize;
  vector<TH1F *> histopvx;
  vector<TH1F *> histopvy;
  vector<TH1F *> histopvz;


  vector<TH1F *> Zbosonpt;
  vector<TH1F *> ZbosonptBBandEB;
  vector<TH1F *> ZbosonptBB;
  vector<TH1F *> ZbosonptEB;
  vector<TH1F *> ZbosonptEE;
  vector<TH1F *> ZbosonptSameSign;
  vector<TH1F *> ZbosonptOppSign;

  vector<TH1F *> ZbosonptZoom;
  vector<TH1F *> ZbosonptZoomBBandEB;
  vector<TH1F *> ZbosonptZoomBB;
  vector<TH1F *> ZbosonptZoomEB;
  vector<TH1F *> ZbosonptZoomEE;
  vector<TH1F *> ZbosonptZoomSameSign;
  vector<TH1F *> ZbosonptZoomOppSign;

  vector<TH1F *> Zbosonpz;
  vector<TH1F *> ZbosonpzBBandEB;
  vector<TH1F *> ZbosonpzBB;
  vector<TH1F *> ZbosonpzEB;
  vector<TH1F *> ZbosonpzEE;
  vector<TH1F *> ZbosonpzSameSign;
  vector<TH1F *> ZbosonpzOppSign;


  vector<TH1F *> ee_mass;
  vector<TH1F *> ee_massBBandEB;
  vector<TH1F *> ee_massBB;
  vector<TH1F *> ee_massEB;
  vector<TH1F *> ee_massEE;
  vector<TH1F *> ee_massSameSign;
  vector<TH1F *> ee_massOppSign;
  vector<TH1F *> gsfgsf_mass;

  vector<TH1F *> ee_massZoom;
  vector<TH1F *> ee_massZoomBBandEB;
  vector<TH1F *> ee_massZoomBB;
  vector<TH1F *> ee_massZoomEB;
  vector<TH1F *> ee_massZoomEE;
  vector<TH1F *> ee_massZoomSameSign;
  vector<TH1F *> ee_massZoomOppSign;

  vector<TH1F *> histogsf_deltaeta;
  vector<TH1F *> histogsf_deltaphi;
  vector<TH1F *> histogsf_hovere;
  vector<TH1F *> histogsf_trackiso;
  vector<TH1F *> histogsf_ecaliso;
  vector<TH1F *> histogsf_hcaliso1;
  vector<TH1F *> histogsf_hcaliso2;
  vector<TH1F *> histogsf_sigmaetaeta;
  vector<TH1F *> histogsf_sigmaIetaIeta;
  vector<TH1F *> histogsf_e2x5overe5x5;

  vector<TH1F *> histogsfbarrel_deltaeta;
  vector<TH1F *> histogsfbarrel_deltaphi;
  vector<TH1F *> histogsfbarrel_hovere;
  vector<TH1F *> histogsfbarrel_trackiso;
  vector<TH1F *> histogsfbarrel_ecaliso;
  vector<TH1F *> histogsfbarrel_hcaliso1;
  vector<TH1F *> histogsfbarrel_hcaliso2;
  vector<TH1F *> histogsfbarrel_sigmaetaeta;
  vector<TH1F *> histogsfbarrel_sigmaIetaIeta;
  vector<TH1F *> histogsfbarrel_e2x5overe5x5;

  vector<TH1F *> histogsfendcap_deltaeta;
  vector<TH1F *> histogsfendcap_deltaphi;
  vector<TH1F *> histogsfendcap_hovere;
  vector<TH1F *> histogsfendcap_trackiso;
  vector<TH1F *> histogsfendcap_ecaliso;
  vector<TH1F *> histogsfendcap_hcaliso1;
  vector<TH1F *> histogsfendcap_hcaliso2;
  vector<TH1F *> histogsfendcap_sigmaetaeta;
  vector<TH1F *> histogsfendcap_sigmaIetaIeta;
  vector<TH1F *> histogsfendcap_e2x5overe5x5;

  vector<TH1F *> histogsf_gsfet;
  vector<TH1F *> histogsf_pt;
  vector<TH1F *> histogsf_eta;
  vector<TH1F *> histogsf_phi;
  vector<TH1F *> histogsfsc_eta;
  vector<TH1F *> histogsfsc_phi;

  vector<TH1F *> histogsfbarrel_gsfet;
  vector<TH1F *> histogsfbarrel_pt;
  vector<TH1F *> histogsfbarrel_eta;
  vector<TH1F *> histogsfbarrel_phi;
  vector<TH1F *> histogsfscbarrel_eta;
  vector<TH1F *> histogsfscbarrel_phi;

  vector<TH1F *> histogsfendcap_gsfet;
  vector<TH1F *> histogsfendcap_pt;
  vector<TH1F *> histogsfendcap_eta;
  vector<TH1F *> histogsfendcap_phi;
  vector<TH1F *> histogsfscendcap_eta;
  vector<TH1F *> histogsfscendcap_phi;

  vector<TH1F *> histogsf_gsfet_M120;
  vector<TH1F *> histogsfsc_eta_M120;
  vector<TH1F *> histogsfsc_phi_M120;

  vector<TH1F *> histogsfgsf_et;
  vector<TH1F *> histogsfgsf_pt;
  vector<TH1F *> histogsfgsf_eta;
  vector<TH1F *> histogsfgsf_phi;

  vector<TH1F *> histoprobept;
  vector<TH1F *> histoprobebarrelpt;
  vector<TH1F *> histoprobeendcappt;
  vector<TH1F *> histoprobeeta;
  vector<TH1F *> histoprobesceta;
  vector<TH1F *> histoprobephi;
  vector<TH1F *> histoprobescphi;
  vector<TH1F *> histoprobeenergy;
  vector<TH1F *> histoprobescenergy;
  vector<TH1F *> histoprobeEoverP;

  vector<TH1F *> histoprobe_deltaeta;
  vector<TH1F *> histoprobe_deltaphi;
  vector<TH1F *> histoprobe_hovere;
  vector<TH1F *> histoprobe_trackiso;
  vector<TH1F *> histoprobe_ecaliso;
  vector<TH1F *> histoprobe_hcaliso1;
  vector<TH1F *> histoprobe_hcaliso2;
  vector<TH1F *> histoprobe_sigmaetaeta;
  vector<TH1F *> histoprobe_sigmaIetaIeta;

  vector<TH1F *> histotagpt;
  vector<TH1F *> histotageta;
  vector<TH1F *> histotagsceta;
  vector<TH1F *> histotagphi;
  vector<TH1F *> histotagscphi;
  vector<TH1F *> histotagenergy;
  vector<TH1F *> histotagscenergy;
  vector<TH1F *> histotagEoverP;

  vector<TH1F *> histotag_deltaeta;
  vector<TH1F *> histotag_deltaphi;
  vector<TH1F *> histotag_hovere;
  vector<TH1F *> histotag_trackiso;
  vector<TH1F *> histotag_ecaliso;
  vector<TH1F *> histotag_hcaliso1;
  vector<TH1F *> histotag_hcaliso2;
  vector<TH1F *> histotag_sigmaetaeta;
  vector<TH1F *> histotag_sigmaIetaIeta;
  vector<TH1F *> histotag_gsfet;

  vector<TH1F *> histoprobeHEEPpt;
  vector<TH1F *> histoprobeHEEPbarrelpt;
  vector<TH1F *> histoprobeHEEPendcappt;
  vector<TH1F *> histoprobeHEEPeta;
  vector<TH1F *> histoprobeHEEPsceta;
  vector<TH1F *> histoprobeHEEPphi;
  vector<TH1F *> histoprobeHEEPscphi;
  vector<TH1F *> histoprobeHEEPenergy;
  vector<TH1F *> histoprobeHEEPscenergy;

  vector<TH1F *> histoEFFIpt;
  vector<TH1F *> histoEFFIeta;
  vector<TH1F *> histoEFFIphi;
  vector<TH1F *> histoEFFIenergy;

  vector<TH1F *> tagprobe_mass;
  vector<TH1F *> tagprobe_mass_80_100;
  vector<TH1F *> tagprobe_num;
  vector<TH1F *> tagprobe_num_80_100;

  vector<TH1F *> CaloMet;
  vector<TH1F *> Met;

  vector<TH1F *> tagprobedeltaphi;

  vector<TH1F *> NJetsAK;
  vector<TH1F *> NJetsIC;

  for(unsigned int i=0;i<nbFile;i++){
    TString wholename((input[i])->GetName());
    int index1 = wholename.Index("FullTrees/")+10;
    int index2 = wholename.Index(".root");
    TString name(wholename(index1,index2-index1));


    histopvx.push_back(new TH1F(TString("pvx_indic").Append(name.Data()),"",40,-2.,2.));
    histopvy.push_back(new TH1F(TString("pvy_indic").Append(name.Data()),"",40,-2.,2.));
    histopvz.push_back(new TH1F(TString("pvz_indic").Append(name.Data()),"",500,-50.,50.));
    histopvsize.push_back(new TH1F(TString("pvsize_indic").Append(name.Data()),"",20,0.,20.));


    Zbosonpt.push_back(new TH1F(TString("Zbosonpt_indic").Append(name.Data()),"",100,0.,500.));
    ZbosonptBBandEB.push_back(new TH1F(TString("ZbosonptBBandEB_indic").Append(name.Data()),"",100,0.,500.));
    ZbosonptBB.push_back(new TH1F(TString("ZbosonptBB_indic").Append(name.Data()),"",100,0.,500.));
    ZbosonptEB.push_back(new TH1F(TString("ZbosonptEB_indic").Append(name.Data()),"",100,0.,500.));
    ZbosonptEE.push_back(new TH1F(TString("ZbosonptEE_indic").Append(name.Data()),"",100,0.,500.));
    ZbosonptSameSign.push_back(new TH1F(TString("ZbosonptSameSign_indic").Append(name.Data()),"",100,0.,500.));
    ZbosonptOppSign.push_back(new TH1F(TString("ZbosonptOppSign_indic").Append(name.Data()),"",100,0.,500.));

    ZbosonptZoom.push_back(new TH1F(TString("ZbosonptZoom_indic").Append(name.Data()),"",100,0.,100.));
    ZbosonptZoomBBandEB.push_back(new TH1F(TString("ZbosonptZoomBBandEB_indic").Append(name.Data()),"",100,0.,100.));
    ZbosonptZoomBB.push_back(new TH1F(TString("ZbosonptZoomBB_indic").Append(name.Data()),"",100,0.,100.));
    ZbosonptZoomEB.push_back(new TH1F(TString("ZbosonptZoomEB_indic").Append(name.Data()),"",100,0.,100.));
    ZbosonptZoomEE.push_back(new TH1F(TString("ZbosonptZoomEE_indic").Append(name.Data()),"",100,0.,100.));
    ZbosonptZoomSameSign.push_back(new TH1F(TString("ZbosonptZoomSameSign_indic").Append(name.Data()),"",100,0.,100.));
    ZbosonptZoomOppSign.push_back(new TH1F(TString("ZbosonptZoomOppSign_indic").Append(name.Data()),"",100,0.,100.));

    Zbosonpz.push_back(new TH1F(TString("Zbosonpz_indic").Append(name.Data()),"",100,-1000.,1000.));
    ZbosonpzBBandEB.push_back(new TH1F(TString("ZbosonpzBBandEB_indic").Append(name.Data()),"",100,-1000.,1000.));
    ZbosonpzBB.push_back(new TH1F(TString("ZbosonpzBB_indic").Append(name.Data()),"",100,-1000.,1000.));
    ZbosonpzEB.push_back(new TH1F(TString("ZbosonpzEB_indic").Append(name.Data()),"",100,-1000.,1000.));
    ZbosonpzEE.push_back(new TH1F(TString("ZbosonpzEE_indic").Append(name.Data()),"",100,-1000.,1000.));
    ZbosonpzSameSign.push_back(new TH1F(TString("ZbosonpzSameSign_indic").Append(name.Data()),"",100,-1000.,1000.));
    ZbosonpzOppSign.push_back(new TH1F(TString("ZbosonpzOppSign_indic").Append(name.Data()),"",100,-1000.,1000.));

    ee_mass.push_back(new TH1F(TString("ee_mass_indic").Append(name.Data()),"",150,0.,1000.));
    ee_massBBandEB.push_back(new TH1F(TString("ee_massBBandEB_indic").Append(name.Data()),"",150,0.,1000.));
    ee_massBB.push_back(new TH1F(TString("ee_massBB_indic").Append(name.Data()),"",150,0.,1000.));
    ee_massEB.push_back(new TH1F(TString("ee_massEB_indic").Append(name.Data()),"",150,0.,1000.));
    ee_massEE.push_back(new TH1F(TString("ee_massEE_indic").Append(name.Data()),"",150,0.,1000.));
    ee_massSameSign.push_back(new TH1F(TString("ee_massSameSign_indic").Append(name.Data()),"",150,0.,1000.));
    ee_massOppSign.push_back(new TH1F(TString("ee_massOppSign_indic").Append(name.Data()),"",150,0.,1000.));
    gsfgsf_mass.push_back(new TH1F(TString("gsfgsf_mass_indic").Append(name.Data()),"",150,0.,1000.));

    ee_massZoom.push_back(new TH1F(TString("ee_massZoom_indic").Append(name.Data()),"",100,40.,140.));
    ee_massZoomBBandEB.push_back(new TH1F(TString("ee_massZoomBBandEB_indic").Append(name.Data()),"",100,40.,140.));
    ee_massZoomBB.push_back(new TH1F(TString("ee_massZoomBB_indic").Append(name.Data()),"",100,40.,140.));
    ee_massZoomEB.push_back(new TH1F(TString("ee_massZoomEB_indic").Append(name.Data()),"",100,40.,140.));
    ee_massZoomEE.push_back(new TH1F(TString("ee_massZoomEE_indic").Append(name.Data()),"",100,40.,140.));
    ee_massZoomSameSign.push_back(new TH1F(TString("ee_massZoomSameSign_indic").Append(name.Data()),"",100,40.,140.));
    ee_massZoomOppSign.push_back(new TH1F(TString("ee_massZoomOppSign_indic").Append(name.Data()),"",100,40.,140.));

    histogsf_deltaeta.push_back(new TH1F(TString("histogsf_deltaeta_indic").Append(name.Data()),"",histoNBins["gsf_deltaeta"],histoXMin["gsf_deltaeta"],histoXMax["gsf_deltaeta"]));
    histogsf_deltaphi.push_back(new TH1F(TString("histogsf_deltaphi_indic").Append(name.Data()),"",histoNBins["gsf_deltaphi"],histoXMin["gsf_deltaphi"],histoXMax["gsf_deltaphi"]));
    histogsf_hovere.push_back(new TH1F(TString("histogsf_hovere_indic").Append(name.Data()),"",histoNBins["gsf_hovere"],histoXMin["gsf_hovere"],histoXMax["gsf_hovere"]));
    histogsf_trackiso.push_back(new TH1F(TString("histogsf_trackiso_indic").Append(name.Data()),"",histoNBins["gsf_trackiso"],histoXMin["gsf_trackiso"],histoXMax["gsf_trackiso"]));
    histogsf_ecaliso.push_back(new TH1F(TString("histogsf_ecaliso_indic").Append(name.Data()),"",histoNBins["gsf_ecaliso"],histoXMin["gsf_ecaliso"],histoXMax["gsf_ecaliso"]));
    histogsf_hcaliso1.push_back(new TH1F(TString("histogsf_hcaliso1_indic").Append(name.Data()),"",histoNBins["gsf_hcaliso1"],histoXMin["gsf_hcaliso1"],histoXMax["gsf_hcaliso1"]));
    histogsf_hcaliso2.push_back(new TH1F(TString("histogsf_hcaliso2_indic").Append(name.Data()),"",histoNBins["gsf_hcaliso2"],histoXMin["gsf_hcaliso2"],histoXMax["gsf_hcaliso2"]));
    histogsf_sigmaetaeta.push_back(new TH1F(TString("histogsf_sigmaetaeta_indic").Append(name.Data()),"",histoNBins["gsf_sigmaetaeta"],histoXMin["gsf_sigmaetaeta"],histoXMax["gsf_sigmaetaeta"]));
    histogsf_sigmaIetaIeta.push_back(new TH1F(TString("histogsf_sigmaIetaIeta_indic").Append(name.Data()),"",histoNBins["gsf_sigmaIetaIeta"],histoXMin["gsf_sigmaIetaIeta"],histoXMax["gsf_sigmaIetaIeta"]));
    histogsf_e2x5overe5x5.push_back(new TH1F(TString("histogsf_e2x5overe5x5_indic").Append(name.Data()),"",histoNBins["gsf_e2x5overe5x5"],histoXMin["gsf_e2x5overe5x5"],histoXMax["gsf_e2x5overe5x5"]));

    histogsfbarrel_deltaeta.push_back(new TH1F(TString("histogsfbarrel_deltaeta_indic").Append(name.Data()),"",histoNBins["gsf_deltaeta"],histoXMin["gsf_deltaeta"],histoXMax["gsf_deltaeta"]));
    histogsfbarrel_deltaphi.push_back(new TH1F(TString("histogsfbarrel_deltaphi_indic").Append(name.Data()),"",histoNBins["gsf_deltaphi"],histoXMin["gsf_deltaphi"],histoXMax["gsf_deltaphi"]));
    histogsfbarrel_hovere.push_back(new TH1F(TString("histogsfbarrel_hovere_indic").Append(name.Data()),"",histoNBins["gsf_hovere"],histoXMin["gsf_hovere"],histoXMax["gsf_hovere"]));
    histogsfbarrel_trackiso.push_back(new TH1F(TString("histogsfbarrel_trackiso_indic").Append(name.Data()),"",histoNBins["gsf_trackiso"],histoXMin["gsf_trackiso"],histoXMax["gsf_trackiso"]));
    histogsfbarrel_ecaliso.push_back(new TH1F(TString("histogsfbarrel_ecaliso_indic").Append(name.Data()),"",histoNBins["gsf_ecaliso"],histoXMin["gsf_ecaliso"],histoXMax["gsf_ecaliso"]));
    histogsfbarrel_hcaliso1.push_back(new TH1F(TString("histogsfbarrel_hcaliso1_indic").Append(name.Data()),"",histoNBins["gsf_hcaliso1"],histoXMin["gsf_hcaliso1"],histoXMax["gsf_hcaliso1"]));
    histogsfbarrel_hcaliso2.push_back(new TH1F(TString("histogsfbarrel_hcaliso2_indic").Append(name.Data()),"",histoNBins["gsf_hcaliso2"],histoXMin["gsf_hcaliso2"],histoXMax["gsf_hcaliso2"]));
    histogsfbarrel_sigmaetaeta.push_back(new TH1F(TString("histogsfbarrel_sigmaetaeta_indic").Append(name.Data()),"",histoNBins["gsf_sigmaetaeta"],histoXMin["gsf_sigmaetaeta"],histoXMax["gsf_sigmaetaeta"]));
    histogsfbarrel_sigmaIetaIeta.push_back(new TH1F(TString("histogsfbarrel_sigmaIetaIeta_indic").Append(name.Data()),"",histoNBins["gsf_sigmaIetaIeta"],histoXMin["gsf_sigmaIetaIeta"],histoXMax["gsf_sigmaIetaIeta"]));
    histogsfbarrel_e2x5overe5x5.push_back(new TH1F(TString("histogsfbarrel_e2x5overe5x5_indic").Append(name.Data()),"",histoNBins["gsf_e2x5overe5x5"],histoXMin["gsf_e2x5overe5x5"],histoXMax["gsf_e2x5overe5x5"]));

    histogsfendcap_deltaeta.push_back(new TH1F(TString("histogsfendcap_deltaeta_indic").Append(name.Data()),"",histoNBins["gsf_deltaeta"],histoXMin["gsf_deltaeta"],histoXMax["gsf_deltaeta"]));
    histogsfendcap_deltaphi.push_back(new TH1F(TString("histogsfendcap_deltaphi_indic").Append(name.Data()),"",histoNBins["gsf_deltaphi"],histoXMin["gsf_deltaphi"],histoXMax["gsf_deltaphi"]));
    histogsfendcap_hovere.push_back(new TH1F(TString("histogsfendcap_hovere_indic").Append(name.Data()),"",histoNBins["gsf_hovere"],histoXMin["gsf_hovere"],histoXMax["gsf_hovere"]));
    histogsfendcap_trackiso.push_back(new TH1F(TString("histogsfendcap_trackiso_indic").Append(name.Data()),"",histoNBins["gsf_trackiso"],histoXMin["gsf_trackiso"],histoXMax["gsf_trackiso"]));
    histogsfendcap_ecaliso.push_back(new TH1F(TString("histogsfendcap_ecaliso_indic").Append(name.Data()),"",histoNBins["gsf_ecaliso"],histoXMin["gsf_ecaliso"],histoXMax["gsf_ecaliso"]));
    histogsfendcap_hcaliso1.push_back(new TH1F(TString("histogsfendcap_hcaliso1_indic").Append(name.Data()),"",histoNBins["gsf_hcaliso1"],histoXMin["gsf_hcaliso1"],histoXMax["gsf_hcaliso1"]));
    histogsfendcap_hcaliso2.push_back(new TH1F(TString("histogsfendcap_hcaliso2_indic").Append(name.Data()),"",histoNBins["gsf_hcaliso2"],histoXMin["gsf_hcaliso2"],histoXMax["gsf_hcaliso2"]));
    histogsfendcap_sigmaetaeta.push_back(new TH1F(TString("histogsfendcap_sigmaetaeta_indic").Append(name.Data()),"",histoNBins["gsf_sigmaetaeta"],histoXMin["gsf_sigmaetaeta"],histoXMax["gsf_sigmaetaeta"]));
    histogsfendcap_sigmaIetaIeta.push_back(new TH1F(TString("histogsfendcap_sigmaIetaIeta_indic").Append(name.Data()),"",histoNBins["gsf_sigmaIetaIeta"],histoXMin["gsf_sigmaIetaIeta"],histoXMax["gsf_sigmaIetaIeta"]));
    histogsfendcap_e2x5overe5x5.push_back(new TH1F(TString("histogsfendcap_e2x5overe5x5_indic").Append(name.Data()),"",histoNBins["gsf_e2x5overe5x5"],histoXMin["gsf_e2x5overe5x5"],histoXMax["gsf_e2x5overe5x5"]));

    histogsf_gsfet.push_back(new TH1F(TString("histogsf_gsfet_indic").Append(name.Data()),"",histoNBins["gsf_gsfet"],histoXMin["gsf_gsfet"],histoXMax["gsf_gsfet"]));	
    histogsf_pt.push_back(new TH1F(TString("histogsf_pt_indic").Append(name.Data()),"",histoNBins["gsf_pt"],histoXMin["gsf_pt"],histoXMax["gsf_pt"]));	
    histogsf_eta.push_back(new TH1F(TString("histogsf_eta_indic").Append(name.Data()),"",histoNBins["gsf_eta"],histoXMin["gsf_eta"],histoXMax["gsf_eta"]));	
    histogsf_phi.push_back(new TH1F(TString("histogsf_phi_indic").Append(name.Data()),"",histoNBins["gsf_phi"],histoXMin["gsf_phi"],histoXMax["gsf_phi"]));	
    histogsfsc_eta.push_back(new TH1F(TString("histogsfsc_eta_indic").Append(name.Data()),"",histoNBins["gsf_eta"],histoXMin["gsf_eta"],histoXMax["gsf_eta"]));	
    histogsfsc_phi.push_back(new TH1F(TString("histogsfsc_phi_indic").Append(name.Data()),"",histoNBins["gsf_phi"],histoXMin["gsf_phi"],histoXMax["gsf_phi"]));	

    histogsfbarrel_gsfet.push_back(new TH1F(TString("histogsfbarrel_gsfet_indic").Append(name.Data()),"",histoNBins["gsf_gsfet"],histoXMin["gsf_gsfet"],histoXMax["gsf_gsfet"]));	
    histogsfbarrel_pt.push_back(new TH1F(TString("histogsfbarrel_pt_indic").Append(name.Data()),"",histoNBins["gsf_pt"],histoXMin["gsf_pt"],histoXMax["gsf_pt"]));	
    histogsfbarrel_eta.push_back(new TH1F(TString("histogsfbarrel_eta_indic").Append(name.Data()),"",histoNBins["gsf_eta"],histoXMin["gsf_eta"],histoXMax["gsf_eta"]));	
    histogsfbarrel_phi.push_back(new TH1F(TString("histogsfbarrel_phi_indic").Append(name.Data()),"",histoNBins["gsf_phi"],histoXMin["gsf_phi"],histoXMax["gsf_phi"]));	
    histogsfscbarrel_eta.push_back(new TH1F(TString("histogsfscbarrel_eta_indic").Append(name.Data()),"",histoNBins["gsf_eta"],histoXMin["gsf_eta"],histoXMax["gsf_eta"]));	
    histogsfscbarrel_phi.push_back(new TH1F(TString("histogsfscbarrel_phi_indic").Append(name.Data()),"",histoNBins["gsf_phi"],histoXMin["gsf_phi"],histoXMax["gsf_phi"]));	

    histogsfendcap_gsfet.push_back(new TH1F(TString("histogsfendcap_gsfet_indic").Append(name.Data()),"",histoNBins["gsf_gsfet"],histoXMin["gsf_gsfet"],histoXMax["gsf_gsfet"]));	
    histogsfendcap_pt.push_back(new TH1F(TString("histogsfendcap_pt_indic").Append(name.Data()),"",histoNBins["gsf_pt"],histoXMin["gsf_pt"],histoXMax["gsf_pt"]));	
    histogsfendcap_eta.push_back(new TH1F(TString("histogsfendcap_eta_indic").Append(name.Data()),"",histoNBins["gsf_eta"],histoXMin["gsf_eta"],histoXMax["gsf_eta"]));	
    histogsfendcap_phi.push_back(new TH1F(TString("histogsfendcap_phi_indic").Append(name.Data()),"",histoNBins["gsf_phi"],histoXMin["gsf_phi"],histoXMax["gsf_phi"]));	
    histogsfscendcap_eta.push_back(new TH1F(TString("histogsfscendcap_eta_indic").Append(name.Data()),"",histoNBins["gsf_eta"],histoXMin["gsf_eta"],histoXMax["gsf_eta"]));	
    histogsfscendcap_phi.push_back(new TH1F(TString("histogsfscendcap_phi_indic").Append(name.Data()),"",histoNBins["gsf_phi"],histoXMin["gsf_phi"],histoXMax["gsf_phi"]));	

    histogsf_gsfet_M120.push_back(new TH1F(TString("histogsf_gsfet_M120_indic").Append(name.Data()),"",histoNBins["gsf_gsfet"]/5,histoXMin["gsf_gsfet"],histoXMax["gsf_gsfet"]));	
    histogsfsc_eta_M120.push_back(new TH1F(TString("histogsfsc_eta_M120_indic").Append(name.Data()),"",histoNBins["gsf_eta"]/5,histoXMin["gsf_eta"],histoXMax["gsf_eta"]));	
    histogsfsc_phi_M120.push_back(new TH1F(TString("histogsfsc_phi_M120_indic").Append(name.Data()),"",histoNBins["gsf_phi"]/2,histoXMin["gsf_phi"],histoXMax["gsf_phi"]));	

    histogsfgsf_et.push_back(new TH1F(TString("histogsfgsf_et_indic").Append(name.Data()),"",histoNBins["gsfgsf_et"],histoXMin["gsfgsf_et"],histoXMax["gsfgsf_et"]));	
    histogsfgsf_pt.push_back(new TH1F(TString("histogsfgsf_pt_indic").Append(name.Data()),"",histoNBins["gsfgsf_pt"],histoXMin["gsfgsf_pt"],histoXMax["gsfgsf_pt"]));	
    histogsfgsf_eta.push_back(new TH1F(TString("histogsfgsf_eta_indic").Append(name.Data()),"",histoNBins["gsfgsf_eta"],histoXMin["gsfgsf_eta"],histoXMax["gsfgsf_eta"]));	
    histogsfgsf_phi.push_back(new TH1F(TString("histogsfgsf_phi_indic").Append(name.Data()),"",histoNBins["gsfgsf_phi"],histoXMin["gsfgsf_phi"],histoXMax["gsfgsf_phi"]));	


    //For tag and probe studies
    histotag_gsfet.push_back(new TH1F(TString("histotag_gsfet_indic").Append(name.Data()),"",50,0.,500.)); 
    histotagpt.push_back(new TH1F(TString("histotagpt_indic").Append(name.Data()),"",50,0.,500.)); 
    histotageta.push_back(new TH1F(TString("histotageta_indic").Append(name.Data()),"",60,-3.,3.));
    histotagsceta.push_back(new TH1F(TString("histotagsceta_indic").Append(name.Data()),"",60,-3.,3.));
    histotagphi.push_back(new TH1F(TString("histotagphi_indic").Append(name.Data()),"",36,-PI,PI));
    histotagscphi.push_back(new TH1F(TString("histotagscphi_indic").Append(name.Data()),"",36,-PI,PI));
    histotagenergy.push_back(new TH1F(TString("histotagenergy_indic").Append(name.Data()),"",50,0.,500.));
    histotagscenergy.push_back(new TH1F(TString("histotagscenergy_indic").Append(name.Data()),"",50,0.,500.));
    histotagEoverP.push_back(new TH1F(TString("histotagEoverP_indic").Append(name.Data()),"",500,0.,50.));
    
    histotag_deltaeta.push_back(new TH1F(TString("histotag_deltaeta_indic").Append(name.Data()),"",histoNBins["gsf_deltaeta"],histoXMin["gsf_deltaeta"],histoXMax["gsf_deltaeta"]));
    histotag_deltaphi.push_back(new TH1F(TString("histotag_deltaphi_indic").Append(name.Data()),"",histoNBins["gsf_deltaphi"],histoXMin["gsf_deltaphi"],histoXMax["gsf_deltaphi"]));
    histotag_hovere.push_back(new TH1F(TString("histotag_hovere_indic").Append(name.Data()),"",histoNBins["gsf_hovere"],histoXMin["gsf_hovere"],histoXMax["gsf_hovere"]));
    histotag_trackiso.push_back(new TH1F(TString("histotag_trackiso_indic").Append(name.Data()),"",histoNBins["gsf_trackiso"],histoXMin["gsf_trackiso"],histoXMax["gsf_trackiso"]));
    histotag_ecaliso.push_back(new TH1F(TString("histotag_ecaliso_indic").Append(name.Data()),"",histoNBins["gsf_ecaliso"],histoXMin["gsf_ecaliso"],histoXMax["gsf_ecaliso"]));
    histotag_hcaliso1.push_back(new TH1F(TString("histotag_hcaliso1_indic").Append(name.Data()),"",histoNBins["gsf_hcaliso1"],histoXMin["gsf_hcaliso1"],histoXMax["gsf_hcaliso1"]));
    histotag_hcaliso2.push_back(new TH1F(TString("histotag_hcaliso2_indic").Append(name.Data()),"",histoNBins["gsf_hcaliso2"],histoXMin["gsf_hcaliso2"],histoXMax["gsf_hcaliso2"]));
    histotag_sigmaetaeta.push_back(new TH1F(TString("histotag_sigmaetaeta_indic").Append(name.Data()),"",histoNBins["gsf_sigmaetaeta"],histoXMin["gsf_sigmaetaeta"],histoXMax["gsf_sigmaetaeta"]));
    histotag_sigmaIetaIeta.push_back(new TH1F(TString("histotag_sigmaIetaIeta_indic").Append(name.Data()),"",histoNBins["gsf_sigmaIetaIeta"],histoXMin["gsf_sigmaIetaIeta"],histoXMax["gsf_sigmaIetaIeta"]));



    histoprobept.push_back(new TH1F(TString("histoprobept_indic").Append(name.Data()),"",50,0.,500.)); 
    histoprobebarrelpt.push_back(new TH1F(TString("histoprobebarrelpt_indic").Append(name.Data()),"",50,0.,500.)); 
    histoprobeendcappt.push_back(new TH1F(TString("histoprobeendcappt_indic").Append(name.Data()),"",50,0.,500.)); 
    histoprobeeta.push_back(new TH1F(TString("histoprobeeta_indic").Append(name.Data()),"",60,-3.,3.));
    histoprobesceta.push_back(new TH1F(TString("histoprobesceta_indic").Append(name.Data()),"",60,-3.,3.));
    histoprobephi.push_back(new TH1F(TString("histoprobephi_indic").Append(name.Data()),"",36,-PI,PI));
    histoprobescphi.push_back(new TH1F(TString("histoprobescphi_indic").Append(name.Data()),"",36,-PI,PI));
    histoprobeenergy.push_back(new TH1F(TString("histoprobeenergy_indic").Append(name.Data()),"",50,0.,500.));
    histoprobescenergy.push_back(new TH1F(TString("histoprobescenergy_indic").Append(name.Data()),"",50,0.,500.));
    histoprobeEoverP.push_back(new TH1F(TString("histoprobeEoverP_indic").Append(name.Data()),"",500,0.,50.));
    
    histoprobe_deltaeta.push_back(new TH1F(TString("histoprobe_deltaeta_indic").Append(name.Data()),"",histoNBins["gsf_deltaeta"],histoXMin["gsf_deltaeta"],histoXMax["gsf_deltaeta"]));
    histoprobe_deltaphi.push_back(new TH1F(TString("histoprobe_deltaphi_indic").Append(name.Data()),"",histoNBins["gsf_deltaphi"],histoXMin["gsf_deltaphi"],histoXMax["gsf_deltaphi"]));
    histoprobe_hovere.push_back(new TH1F(TString("histoprobe_hovere_indic").Append(name.Data()),"",histoNBins["gsf_hovere"],histoXMin["gsf_hovere"],histoXMax["gsf_hovere"]));
    histoprobe_trackiso.push_back(new TH1F(TString("histoprobe_trackiso_indic").Append(name.Data()),"",histoNBins["gsf_trackiso"],histoXMin["gsf_trackiso"],histoXMax["gsf_trackiso"]));
    histoprobe_ecaliso.push_back(new TH1F(TString("histoprobe_ecaliso_indic").Append(name.Data()),"",histoNBins["gsf_ecaliso"],histoXMin["gsf_ecaliso"],histoXMax["gsf_ecaliso"]));
    histoprobe_hcaliso1.push_back(new TH1F(TString("histoprobe_hcaliso1_indic").Append(name.Data()),"",histoNBins["gsf_hcaliso1"],histoXMin["gsf_hcaliso1"],histoXMax["gsf_hcaliso1"]));
    histoprobe_hcaliso2.push_back(new TH1F(TString("histoprobe_hcaliso2_indic").Append(name.Data()),"",histoNBins["gsf_hcaliso2"],histoXMin["gsf_hcaliso2"],histoXMax["gsf_hcaliso2"]));
    histoprobe_sigmaetaeta.push_back(new TH1F(TString("histoprobe_sigmaetaeta_indic").Append(name.Data()),"",histoNBins["gsf_sigmaetaeta"],histoXMin["gsf_sigmaetaeta"],histoXMax["gsf_sigmaetaeta"]));
    histoprobe_sigmaIetaIeta.push_back(new TH1F(TString("histoprobe_sigmaIetaIeta_indic").Append(name.Data()),"",histoNBins["gsf_sigmaIetaIeta"],histoXMin["gsf_sigmaIetaIeta"],histoXMax["gsf_sigmaIetaIeta"]));


    histoprobeHEEPpt.push_back(new TH1F(TString("histoprobeHEEPpt_indic").Append(name.Data()),"",50,0.,500.)); 
    histoprobeHEEPbarrelpt.push_back(new TH1F(TString("histoprobeHEEPbarrelpt_indic").Append(name.Data()),"",50,0.,500.)); 
    histoprobeHEEPendcappt.push_back(new TH1F(TString("histoprobeHEEPendcappt_indic").Append(name.Data()),"",50,0.,500.)); 
    histoprobeHEEPeta.push_back(new TH1F(TString("histoprobeHEEPeta_indic").Append(name.Data()),"",60,-3.,3.));
    histoprobeHEEPsceta.push_back(new TH1F(TString("histoprobeHEEPsceta_indic").Append(name.Data()),"",60,-3.,3.));
    histoprobeHEEPphi.push_back(new TH1F(TString("histoprobeHEEPphi_indic").Append(name.Data()),"",36,-PI,PI));
    histoprobeHEEPscphi.push_back(new TH1F(TString("histoprobeHEEPscphi_indic").Append(name.Data()),"",36,-PI,PI));
    histoprobeHEEPenergy.push_back(new TH1F(TString("histoprobeHEEPenergy_indic").Append(name.Data()),"",50,0.,500.));
    histoprobeHEEPscenergy.push_back(new TH1F(TString("histoprobeHEEPscenergy_indic").Append(name.Data()),"",50,0.,500.));
    
    histoEFFIpt.push_back(new TH1F(TString("histoEFFIpt_indic").Append(name.Data()),"",50,0.,500.)); 
    histoEFFIeta.push_back(new TH1F(TString("histoEFFIeta_indic").Append(name.Data()),"",60,-3.,3.));
    histoEFFIphi.push_back(new TH1F(TString("histoEFFIphi_indic").Append(name.Data()),"",36,-PI,PI));
    histoEFFIenergy.push_back(new TH1F(TString("histoEFFIenergy_indic").Append(name.Data()),"",50,0.,500.));
    
    tagprobe_mass.push_back(new TH1F(TString("tagprobe_mass_indic").Append(name.Data()),"",50,0.,500.)); 
    tagprobe_mass_80_100.push_back(new TH1F(TString("tagprobe_mass_80_100_indic").Append(name.Data()),"",50,0.,500.)); 
    tagprobe_num.push_back(new TH1F(TString("tagprobe_num_indic").Append(name.Data()),"",100,0.,10.)); 
    tagprobe_num_80_100.push_back(new TH1F(TString("tagprobe_num_80_100_indic").Append(name.Data()),"",100,0.,10.)); 

    CaloMet.push_back(new TH1F(TString("CaloMet_indic").Append(name.Data()),"",50,0.,200.)); 
    Met.push_back(new TH1F(TString("Met_indic").Append(name.Data()),"",50,0.,200.)); 
    tagprobedeltaphi.push_back(new TH1F(TString("tagprobedeltaphi_indic").Append(name.Data()),"",36,-PI,PI)); 
    NJetsAK.push_back(new TH1F(TString("NJetsAK_indic").Append(name.Data()),"",200,0.,20.)); 
    NJetsIC.push_back(new TH1F(TString("NJetsIC_indic").Append(name.Data()),"",200,0.,20.)); 
  }
  

  for (unsigned int p = 0;p<nbFile;p++){

    (histopvx[p])->Sumw2();
    (histopvy[p])->Sumw2();
    (histopvz[p])->Sumw2();
    (histopvsize[p])->Sumw2();

    (Zbosonpt[p])->Sumw2();
    (ZbosonptBBandEB[p])->Sumw2();
    (ZbosonptBB[p])->Sumw2();
    (ZbosonptEB[p])->Sumw2();
    (ZbosonptEE[p])->Sumw2();
    (ZbosonptSameSign[p])->Sumw2();
    (ZbosonptOppSign[p])->Sumw2();

    (ZbosonptZoom[p])->Sumw2();
    (ZbosonptZoomBBandEB[p])->Sumw2();
    (ZbosonptZoomBB[p])->Sumw2();
    (ZbosonptZoomEB[p])->Sumw2();
    (ZbosonptZoomEE[p])->Sumw2();
    (ZbosonptZoomSameSign[p])->Sumw2();
    (ZbosonptZoomOppSign[p])->Sumw2();

    (Zbosonpz[p])->Sumw2();
    (ZbosonpzBBandEB[p])->Sumw2();
    (ZbosonpzBB[p])->Sumw2();
    (ZbosonpzEB[p])->Sumw2();
    (ZbosonpzEE[p])->Sumw2();
    (ZbosonpzSameSign[p])->Sumw2();
    (ZbosonpzOppSign[p])->Sumw2();

    (ee_mass[p])->Sumw2();
    (ee_massBBandEB[p])->Sumw2();
    (ee_massBB[p])->Sumw2();
    (ee_massEB[p])->Sumw2();
    (ee_massEE[p])->Sumw2();
    (ee_massSameSign[p])->Sumw2();
    (ee_massOppSign[p])->Sumw2();
    (gsfgsf_mass[p])->Sumw2();

    (ee_massZoom[p])->Sumw2();
    (ee_massZoomBBandEB[p])->Sumw2();
    (ee_massZoomBB[p])->Sumw2();
    (ee_massZoomEB[p])->Sumw2();
    (ee_massZoomEE[p])->Sumw2();
    (ee_massZoomSameSign[p])->Sumw2();
    (ee_massZoomOppSign[p])->Sumw2();

    (histogsf_deltaeta[p])->Sumw2(); 
    (histogsf_deltaphi[p])->Sumw2(); 
    (histogsf_hovere[p])->Sumw2(); 
    (histogsf_trackiso[p])->Sumw2(); 
    (histogsf_ecaliso[p])->Sumw2(); 
    (histogsf_hcaliso1[p])->Sumw2(); 
    (histogsf_hcaliso2[p])->Sumw2(); 
    (histogsf_sigmaetaeta[p])->Sumw2(); 
    (histogsf_sigmaIetaIeta[p])->Sumw2(); 
    (histogsf_e2x5overe5x5[p])->Sumw2(); 

    (histogsfbarrel_deltaeta[p])->Sumw2(); 
    (histogsfbarrel_deltaphi[p])->Sumw2(); 
    (histogsfbarrel_hovere[p])->Sumw2(); 
    (histogsfbarrel_trackiso[p])->Sumw2(); 
    (histogsfbarrel_ecaliso[p])->Sumw2(); 
    (histogsfbarrel_hcaliso1[p])->Sumw2(); 
    (histogsfbarrel_hcaliso2[p])->Sumw2(); 
    (histogsfbarrel_sigmaetaeta[p])->Sumw2(); 
    (histogsfbarrel_sigmaIetaIeta[p])->Sumw2(); 
    (histogsfbarrel_e2x5overe5x5[p])->Sumw2(); 

    (histogsfendcap_deltaeta[p])->Sumw2(); 
    (histogsfendcap_deltaphi[p])->Sumw2(); 
    (histogsfendcap_hovere[p])->Sumw2(); 
    (histogsfendcap_trackiso[p])->Sumw2(); 
    (histogsfendcap_ecaliso[p])->Sumw2(); 
    (histogsfendcap_hcaliso1[p])->Sumw2(); 
    (histogsfendcap_hcaliso2[p])->Sumw2(); 
    (histogsfendcap_sigmaetaeta[p])->Sumw2(); 
    (histogsfendcap_sigmaIetaIeta[p])->Sumw2(); 
    (histogsfendcap_e2x5overe5x5[p])->Sumw2(); 

    (histogsf_gsfet[p])->Sumw2(); 
    (histogsf_pt[p])->Sumw2(); 
    (histogsf_eta[p])->Sumw2(); 
    (histogsf_phi[p])->Sumw2(); 
    (histogsfsc_eta[p])->Sumw2(); 
    (histogsfsc_phi[p])->Sumw2(); 

    (histogsfbarrel_gsfet[p])->Sumw2(); 
    (histogsfbarrel_pt[p])->Sumw2(); 
    (histogsfbarrel_eta[p])->Sumw2(); 
    (histogsfbarrel_phi[p])->Sumw2(); 
    (histogsfscbarrel_eta[p])->Sumw2(); 
    (histogsfscbarrel_phi[p])->Sumw2(); 

    (histogsfendcap_gsfet[p])->Sumw2(); 
    (histogsfendcap_pt[p])->Sumw2(); 
    (histogsfendcap_eta[p])->Sumw2(); 
    (histogsfendcap_phi[p])->Sumw2(); 
    (histogsfscendcap_eta[p])->Sumw2(); 
    (histogsfscendcap_phi[p])->Sumw2(); 

    (histogsf_gsfet_M120[p])->Sumw2(); 
    (histogsfsc_eta_M120[p])->Sumw2(); 
    (histogsfsc_phi_M120[p])->Sumw2(); 

    (histogsfgsf_et[p])->Sumw2(); 
    (histogsfgsf_pt[p])->Sumw2(); 
    (histogsfgsf_eta[p])->Sumw2(); 
    (histogsfgsf_phi[p])->Sumw2(); 
    
    (histoprobept[p])->Sumw2();
    (histoprobebarrelpt[p])->Sumw2();
    (histoprobeendcappt[p])->Sumw2();
    (histoprobeeta[p])->Sumw2();
    (histoprobesceta[p])->Sumw2();
    (histoprobephi[p])->Sumw2();
    (histoprobescphi[p])->Sumw2();
    (histoprobeenergy[p])->Sumw2();
    (histoprobescenergy[p])->Sumw2();
    (histoprobeEoverP[p])->Sumw2();
        
    (histoprobe_deltaeta[p])->Sumw2(); 
    (histoprobe_deltaphi[p])->Sumw2(); 
    (histoprobe_hovere[p])->Sumw2(); 
    (histoprobe_trackiso[p])->Sumw2(); 
    (histoprobe_ecaliso[p])->Sumw2(); 
    (histoprobe_hcaliso1[p])->Sumw2(); 
    (histoprobe_hcaliso2[p])->Sumw2(); 
    (histoprobe_sigmaetaeta[p])->Sumw2(); 
    (histoprobe_sigmaIetaIeta[p])->Sumw2(); 

    (histotagpt[p])->Sumw2();
    (histotageta[p])->Sumw2();
    (histotagsceta[p])->Sumw2();
    (histotagphi[p])->Sumw2();
    (histotagscphi[p])->Sumw2();
    (histotagenergy[p])->Sumw2();
    (histotagscenergy[p])->Sumw2();
    (histotagEoverP[p])->Sumw2();
    
    (histotag_deltaeta[p])->Sumw2(); 
    (histotag_deltaphi[p])->Sumw2(); 
    (histotag_hovere[p])->Sumw2(); 
    (histotag_trackiso[p])->Sumw2(); 
    (histotag_ecaliso[p])->Sumw2(); 
    (histotag_hcaliso1[p])->Sumw2(); 
    (histotag_hcaliso2[p])->Sumw2(); 
    (histotag_sigmaetaeta[p])->Sumw2(); 
    (histotag_sigmaIetaIeta[p])->Sumw2(); 
    (histotag_gsfet[p])->Sumw2(); 

    (histoprobeHEEPpt[p])->Sumw2();
    (histoprobeHEEPbarrelpt[p])->Sumw2();
    (histoprobeHEEPendcappt[p])->Sumw2();
    (histoprobeHEEPeta[p])->Sumw2();
    (histoprobeHEEPsceta[p])->Sumw2();
    (histoprobeHEEPphi[p])->Sumw2();
    (histoprobeHEEPscphi[p])->Sumw2();
    (histoprobeHEEPenergy[p])->Sumw2();
    (histoprobeHEEPscenergy[p])->Sumw2();
    
    (histoEFFIpt[p])->Sumw2();
    (histoEFFIeta[p])->Sumw2();
    (histoEFFIphi[p])->Sumw2();
    (histoEFFIenergy[p])->Sumw2();

    (tagprobe_mass[p])->Sumw2();
    (tagprobe_mass_80_100[p])->Sumw2();
    (tagprobe_num[p])->Sumw2();
    (tagprobe_num_80_100[p])->Sumw2();

    (CaloMet[p])->Sumw2();
    (Met[p])->Sumw2();
    (tagprobedeltaphi[p])->Sumw2();
    (NJetsAK[p])->Sumw2();
    (NJetsIC[p])->Sumw2();
  }

  for (int p = 0; p<nbFile; p++) {
    
    cout<<"file "<<(input[p])->GetName()<<endl;
    (input[p])->cd();
    TTree *thetree;
    //if(p == 0) thetree = (TTree*)(input[p])->Get("tree");
    if(p == 0) thetree = (TTree*)(input[p])->Get("gsfcheckerjob/tree");
    if(p != 0) thetree = (TTree*)(input[p])->Get("gsfcheckerjob/tree");
    Init(thetree);
    Long64_t nentries = (*thetree).GetEntries();
    cout<<nentries<<" entries"<<endl;

    float BarrelEnergyCorrectionFactor = 1.00;
    float EndcapEnergyCorrectionFactor = 1.00;
    if(p == 0) BarrelEnergyCorrectionFactor = 1.00;
    if(p == 0) EndcapEnergyCorrectionFactor = 1.04;

    int nbheep = 0;
    int indexele1 = -1;
    int indexele2 = -1;
    bool isbarrel = false;
    bool isendcap = false;

    bool tagisbarrel = false;
    bool tagisendcap = false;
    bool probeisbarrel = false;
    bool probeisendcap = false;

    int nbevents120toInf = 0;
    int nbevents200toInf = 0;
    int nbevents20toInf = 0;
    int nbevents40toInf = 0;
    int nbevents80to100 = 0;
    int nbevents60to120 = 0;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;

      thetree->GetEntry(jentry);
      if( (jentry%50000) == 0 ) cout<<"entry "<<jentry<<endl;
      //if(jentry > 1000) break;

      nbheep = 0;

      indexele1 = -1;
      indexele2 = -1;

      isbarrel = false;
      isendcap = false;

      if(gsf_size < 2) continue;
//       if(p == 0 && runnumber > 149442) continue;
//       if(p == 0 && runnumber < 141950) continue;

      TString samplename((input[p])->GetName());

      if(samplename.Contains("QCD15V6") && pthat > 30.) continue;
      if(samplename.Contains("QCD30V6") && pthat > 80.) continue;
      if(samplename.Contains("QCD80V6") && pthat > 170.) continue;
      

      TLorentzVector gsfele1;
      TLorentzVector gsfele2;
      
      //Put back deltaeta cut in the endcaps

      bool gsfele1barrel = fabs(gsfsc_eta[0]) < 1.442 ;
      bool gsfele1endcap = fabs(gsfsc_eta[0]) > 1.56 && fabs(gsfsc_eta[0]) < 2.5;

      bool gsfele2barrel = fabs(gsfsc_eta[1]) < 1.442 ;
      bool gsfele2endcap = fabs(gsfsc_eta[1]) > 1.56 && fabs(gsfsc_eta[1]) < 2.5;


      bool gsf1preselbarrel = (gsfele1barrel && gsf_gsfet[0] > 25. && gsfsc_pt[0] > 4. && gsf_isecaldriven[0] && fabs(gsf_deltaeta[0]) < 0.02 && fabs(gsf_deltaphi[0]) < 0.15 && gsf_hovere[0] < 0.1 && gsf_SwissCross[0] < 0.95);
      bool gsf1preselendcap = (gsfele1endcap && gsf_gsfet[0] > 25. && gsfsc_pt[0] > 4. && gsf_isecaldriven[0] && fabs(gsf_deltaeta[0]) < 0.02 && fabs(gsf_deltaphi[0]) < 0.15 && gsf_hovere[0] < 0.1 && gsf_SwissCross[0] < 0.95);

      bool gsf2preselbarrel = (gsfele2barrel && gsf_gsfet[1] > 25. && gsfsc_pt[1] > 4. && gsf_isecaldriven[1] && fabs(gsf_deltaeta[1]) < 0.02 && fabs(gsf_deltaphi[1]) < 0.15 && gsf_hovere[1] < 0.1 && gsf_SwissCross[1] < 0.95);
      bool gsf2preselendcap = (gsfele2endcap && gsf_gsfet[1] > 25. && gsfsc_pt[1] > 4. && gsf_isecaldriven[1] && fabs(gsf_deltaeta[1]) < 0.02 && fabs(gsf_deltaphi[1]) < 0.15 && gsf_hovere[1] < 0.1 && gsf_SwissCross[1] < 0.95);


      bool gsfpreselection = ((gsf1preselbarrel || gsf1preselendcap) && (gsf2preselbarrel || gsf2preselendcap));
      
      if(gsfpreselection && fabs(pvz[0]) < 24.) {

	if(gsfele1endcap) gsfele1.SetPtEtaPhiE(gsf_gsfet[0]*EndcapEnergyCorrectionFactor,gsf_eta[0],gsf_phi[0],(gsf_gsfet[0]*EndcapEnergyCorrectionFactor*cosh(gsf_eta[0])));
	if(gsfele2endcap) gsfele2.SetPtEtaPhiE(gsf_gsfet[1]*EndcapEnergyCorrectionFactor,gsf_eta[1],gsf_phi[1],(gsf_gsfet[1]*EndcapEnergyCorrectionFactor*cosh(gsf_eta[1])));
	if(gsfele1barrel) gsfele1.SetPtEtaPhiE(gsf_gsfet[0]*BarrelEnergyCorrectionFactor,gsf_eta[0],gsf_phi[0],(gsf_gsfet[0]*BarrelEnergyCorrectionFactor*cosh(gsf_eta[0])));
	if(gsfele2barrel) gsfele2.SetPtEtaPhiE(gsf_gsfet[1]*BarrelEnergyCorrectionFactor,gsf_eta[1],gsf_phi[1],(gsf_gsfet[1]*BarrelEnergyCorrectionFactor*cosh(gsf_eta[1])));
	
	float gsfmass = (gsfele1+gsfele2).Mag();
	gsfgsf_mass[p]->Fill(gsfmass);
	
	if(gsfele1endcap) (histogsfgsf_et[p])->Fill(gsf_gsfet[0]*EndcapEnergyCorrectionFactor); 
	if(gsfele1barrel) (histogsfgsf_et[p])->Fill(gsf_gsfet[0]*BarrelEnergyCorrectionFactor); 
	(histogsfgsf_pt[p])->Fill(gsf_pt[0]); 
	(histogsfgsf_eta[p])->Fill(gsf_eta[0]); 
	(histogsfgsf_phi[p])->Fill(gsf_phi[0]); 
	
	if(gsfele2endcap) (histogsfgsf_et[p])->Fill(gsf_gsfet[1]*EndcapEnergyCorrectionFactor); 
	if(gsfele2barrel) (histogsfgsf_et[p])->Fill(gsf_gsfet[1]*BarrelEnergyCorrectionFactor); 
	(histogsfgsf_pt[p])->Fill(gsf_pt[1]); 
	(histogsfgsf_eta[p])->Fill(gsf_eta[1]); 
	(histogsfgsf_phi[p])->Fill(gsf_phi[1]); 
	
	//if(p == 0) cout<<"gsfgsf mass"<<gsfmass<<endl;
      }

      //FOR INFO, IN THE CODE
      //gsfpass_ID[e] = (gsfpass_DETAIN[e] && gsfpass_DPHIIN[e] && gsfpass_HADEM[e] && gsfpass_SIGMAIETAIETA[e] && gsfpass_E2X5OVER5X5[e]);
      //gsfpass_ISO[e] = (gsfpass_ISOLEMHADDEPTH1[e] && gsfpass_ISOLHADDEPTH2[e] && gsfpass_ISOLPTTRKS[e]);

      for(int i=0;i<gsf_size;i++) {

	if(fabs(pvz[0]) > 24.) continue;
	if(gsf_SwissCross[i] > 0.95) continue;
	if(gsf_gsfet[i] < 25.) continue;
	isbarrel = fabs(gsfsc_eta[i]) < 1.442;
	isendcap = (fabs(gsfsc_eta[i]) > 1.56) && (fabs(gsfsc_eta[i]) < 2.5);
	if(!(isbarrel || isendcap)) continue;

	if(!gsfpass_ECALDRIVEN[i]) continue;

	bool heepbarrel = (isbarrel && gsfpass_ISO[i] && gsfpass_HADEM[i] && gsfpass_DETAIN[i] && gsfpass_DPHIIN[i] && gsfpass_SIGMAIETAIETA[i] && gsfpass_E2X5OVER5X5[i]);
	bool heependcap = (isendcap && gsfpass_ISO[i] && gsfpass_HADEM[i] && gsfpass_DETAIN[i] && gsfpass_DPHIIN[i] && gsfpass_SIGMAIETAIETA[i] && gsfpass_E2X5OVER5X5[i]);

	if(heepbarrel || heependcap){
	  if(indexele1 >= 0 && indexele2 == -1) indexele2 = i;
	  if(indexele1 == -1 && indexele2 == -1) indexele1 = i;
	}
	
      }
      
      //cout<<eventnumber<<" "<<indexele1<<" "<<indexele2<<endl;
      
      if(indexele1 >= 0 && indexele2 >= 0) {
	TLorentzVector ele1;
	TLorentzVector ele2;
	
	bool ele1barrel = fabs(gsfsc_eta[indexele1]) < 1.442;
	bool ele1endcap = fabs(gsfsc_eta[indexele1]) > 1.56 && fabs(gsfsc_eta[indexele1]) < 2.5;

	bool ele2barrel = fabs(gsfsc_eta[indexele2]) < 1.442;
	bool ele2endcap = fabs(gsfsc_eta[indexele2]) > 1.56 && fabs(gsfsc_eta[indexele2]) < 2.5;

	if(ele1barrel) ele1.SetPtEtaPhiE(gsf_gsfet[indexele1]*BarrelEnergyCorrectionFactor,gsf_eta[indexele1],gsf_phi[indexele1],(gsf_gsfet[indexele1]*BarrelEnergyCorrectionFactor*cosh(gsf_eta[indexele1])));
	if(ele1endcap) ele1.SetPtEtaPhiE(gsf_gsfet[indexele1]*EndcapEnergyCorrectionFactor,gsf_eta[indexele1],gsf_phi[indexele1],(gsf_gsfet[indexele1]*EndcapEnergyCorrectionFactor*cosh(gsf_eta[indexele1])));
	if(ele2barrel) ele2.SetPtEtaPhiE(gsf_gsfet[indexele2]*BarrelEnergyCorrectionFactor,gsf_eta[indexele2],gsf_phi[indexele2],(gsf_gsfet[indexele2]*BarrelEnergyCorrectionFactor*cosh(gsf_eta[indexele2])));
	if(ele2endcap) ele2.SetPtEtaPhiE(gsf_gsfet[indexele2]*EndcapEnergyCorrectionFactor,gsf_eta[indexele2],gsf_phi[indexele2],(gsf_gsfet[indexele2]*EndcapEnergyCorrectionFactor*cosh(gsf_eta[indexele2])));

	float mass = (ele1+ele2).Mag();
	float pz = (ele1+ele2).Pz();
	float pt = (ele1+ele2).Perp();

	if(ele1barrel) (histogsf_gsfet[p])->Fill(gsf_gsfet[indexele1]*BarrelEnergyCorrectionFactor); 
	if(ele1endcap) (histogsf_gsfet[p])->Fill(gsf_gsfet[indexele1]*EndcapEnergyCorrectionFactor); 
	(histogsf_pt[p])->Fill(gsf_pt[indexele1]); 
	(histogsf_eta[p])->Fill(gsf_eta[indexele1]); 
	(histogsf_phi[p])->Fill(gsf_phi[indexele1]); 
	(histogsfsc_eta[p])->Fill(gsfsc_eta[indexele1]); 
	(histogsfsc_phi[p])->Fill(gsfsc_phi[indexele1]); 

	if(ele1barrel){ 
	  (histogsfbarrel_gsfet[p])->Fill(gsf_gsfet[indexele1]*BarrelEnergyCorrectionFactor); 
	  (histogsfbarrel_pt[p])->Fill(gsf_pt[indexele1]); 
	  (histogsfbarrel_eta[p])->Fill(gsf_eta[indexele1]); 
	  (histogsfbarrel_phi[p])->Fill(gsf_phi[indexele1]); 
	  (histogsfscbarrel_eta[p])->Fill(gsfsc_eta[indexele1]); 
	  (histogsfscbarrel_phi[p])->Fill(gsfsc_phi[indexele1]); 
	}

	if(ele1endcap){ 
	  (histogsfendcap_gsfet[p])->Fill(gsf_gsfet[indexele1]*EndcapEnergyCorrectionFactor); 
	  (histogsfendcap_pt[p])->Fill(gsf_pt[indexele1]); 
	  (histogsfendcap_eta[p])->Fill(gsf_eta[indexele1]); 
	  (histogsfendcap_phi[p])->Fill(gsf_phi[indexele1]); 
	  (histogsfscendcap_eta[p])->Fill(gsfsc_eta[indexele1]); 
	  (histogsfscendcap_phi[p])->Fill(gsfsc_phi[indexele1]); 
	}

	if(ele2barrel) (histogsf_gsfet[p])->Fill(gsf_gsfet[indexele2]*BarrelEnergyCorrectionFactor); 
	if(ele2endcap) (histogsf_gsfet[p])->Fill(gsf_gsfet[indexele2]*EndcapEnergyCorrectionFactor); 
	(histogsf_pt[p])->Fill(gsf_pt[indexele2]); 
	(histogsf_eta[p])->Fill(gsf_eta[indexele2]); 
	(histogsf_phi[p])->Fill(gsf_phi[indexele2]); 
	(histogsfsc_eta[p])->Fill(gsfsc_eta[indexele2]); 
	(histogsfsc_phi[p])->Fill(gsfsc_phi[indexele2]); 

	if(ele2barrel){ 
	  (histogsfbarrel_gsfet[p])->Fill(gsf_gsfet[indexele2]*BarrelEnergyCorrectionFactor); 
	  (histogsfbarrel_pt[p])->Fill(gsf_pt[indexele2]); 
	  (histogsfbarrel_eta[p])->Fill(gsf_eta[indexele2]); 
	  (histogsfbarrel_phi[p])->Fill(gsf_phi[indexele2]); 
	  (histogsfscbarrel_eta[p])->Fill(gsfsc_eta[indexele2]); 
	  (histogsfscbarrel_phi[p])->Fill(gsfsc_phi[indexele2]); 
	}

	if(ele2endcap){ 
	  (histogsfendcap_gsfet[p])->Fill(gsf_gsfet[indexele2]*EndcapEnergyCorrectionFactor); 
	  (histogsfendcap_pt[p])->Fill(gsf_pt[indexele2]); 
	  (histogsfendcap_eta[p])->Fill(gsf_eta[indexele2]); 
	  (histogsfendcap_phi[p])->Fill(gsf_phi[indexele2]); 
	  (histogsfscendcap_eta[p])->Fill(gsfsc_eta[indexele2]); 
	  (histogsfscendcap_phi[p])->Fill(gsfsc_phi[indexele2]); 
	}

	//if(p == 0) cout<<" mass"<<mass<<endl;
	(histopvx[p])->Fill(pvx[0]);
	(histopvy[p])->Fill(pvy[0]);
	(histopvz[p])->Fill(pvz[0]);
	(histopvsize[p])->Fill(pvsize*1.);

	ee_mass[p]->Fill(mass);
	Zbosonpz[p]->Fill(pz);
	Zbosonpt[p]->Fill(pt);

	if(!(ele1endcap && ele2endcap)){
	  if (mass>600) {
	    cout<<" runnumber = "<<runnumber<<" eventnumber = "<<eventnumber<<endl;
	    cout<<"mass = "<<mass<<endl;
	    cout<<"gsfsc_eta[indexele1] = "<<gsfsc_eta[indexele1]<<endl;
	    cout<<"gsfsc_eta[indexele2] = "<<gsfsc_eta[indexele2]<<endl;
	    cout<<"gsfsc_phi[indexele1] = "<<gsfsc_phi[indexele1]<<endl;
	    cout<<"gsfsc_phi[indexele2] = "<<gsfsc_phi[indexele2]<<endl;
	    cout<<"gsf_et[indexele1] = "<<gsf_gsfet[indexele1]<<endl;
	    cout<<"gsf_et[indexele2] = "<<gsf_gsfet[indexele2]<<endl;
	    cout<<"gsf_pt[indexele1] = "<<gsf_pt[indexele1]<<endl;
	    cout<<"gsf_pt[indexele2] = "<<gsf_pt[indexele2]<<endl;
	    //	    cout<<" = "<<<<endl;
	  }

	  ee_massBBandEB[p]->Fill(mass);
	  ZbosonpzBBandEB[p]->Fill(pz);
	  ZbosonptBBandEB[p]->Fill(pt);
	}

	if(ele1barrel && ele2barrel){ 
	  ee_massBB[p]->Fill(mass);
	  ZbosonpzBB[p]->Fill(pz);
	  ZbosonptBB[p]->Fill(pt);
	}

	if( (ele1barrel && ele2endcap) || (ele1endcap && ele2barrel) ){
	  ee_massEB[p]->Fill(mass);
	  ZbosonpzEB[p]->Fill(pz);
	  ZbosonptEB[p]->Fill(pt);
	}
	
	if(ele1endcap && ele2endcap){
	  ee_massEE[p]->Fill(mass);
	  ZbosonpzEE[p]->Fill(pz);
	  ZbosonptEE[p]->Fill(pt);
	}

	if( (gsf_charge[indexele2]*gsf_charge[indexele1]) > 0 ){ 
	  ee_massSameSign[p]->Fill(mass);
	  ZbosonpzSameSign[p]->Fill(pz);
	  ZbosonptSameSign[p]->Fill(pt);
	}

	if( (gsf_charge[indexele2]*gsf_charge[indexele1]) < 0 ){ 
	  ee_massOppSign[p]->Fill(mass);
	  ZbosonpzOppSign[p]->Fill(pz);
	  ZbosonptOppSign[p]->Fill(pt);
	}


	//Zoom on mass 0-200
	if(mass > 40 && mass < 140.){
	  ee_massZoom[p]->Fill(mass);
	  
	  if(!(ele1endcap && ele2endcap)){
	    ee_massZoomBBandEB[p]->Fill(mass);
	  }
	  
	  if(ele1barrel && ele2barrel){ 
	    ee_massZoomBB[p]->Fill(mass);
	  }
	  
	  if( (ele1barrel && ele2endcap) || (ele1endcap && ele2barrel) ){
	    ee_massZoomEB[p]->Fill(mass);
	  }
	  
	  if(ele1endcap && ele2endcap){
	    ee_massZoomEE[p]->Fill(mass);
	  }
	  
	  if( (gsf_charge[indexele2]*gsf_charge[indexele1]) > 0 ){ 
	    ee_massZoomSameSign[p]->Fill(mass);
	  }
	  
	  if( (gsf_charge[indexele2]*gsf_charge[indexele1]) < 0 ){ 
	    ee_massZoomOppSign[p]->Fill(mass);
	  }
	}


	//Zoom on pt 0-100
	if(pt < 100){
	  ZbosonptZoom[p]->Fill(pt);
	  
	  if(!(ele1endcap && ele2endcap)){
	    ZbosonptZoomBBandEB[p]->Fill(pt);
	  }
	  
	  if(ele1barrel && ele2barrel){ 
	    ZbosonptZoomBB[p]->Fill(pt);
	  }
	  
	  if( (ele1barrel && ele2endcap) || (ele1endcap && ele2barrel) ){
	    ZbosonptZoomEB[p]->Fill(pt);
	  }
	  
	  if(ele1endcap && ele2endcap){
	    ZbosonptZoomEE[p]->Fill(pt);
	  }
	  
	  if( (gsf_charge[indexele2]*gsf_charge[indexele1]) > 0 ){ 
	    ZbosonptZoomSameSign[p]->Fill(pt);
	  }
	  
	  if( (gsf_charge[indexele2]*gsf_charge[indexele1]) < 0 ){ 
	    ZbosonptZoomOppSign[p]->Fill(pt);
	  }
	}




	// 	if(mass > 120.) {
	// 	  (histogsf_gsfet_M120[p])->Fill(gsf_gsfet[indexele1]); 
	// 	  (histogsfsc_eta_M120[p])->Fill(gsfsc_eta[indexele1]); 
	// 	  (histogsfsc_phi_M120[p])->Fill(gsfsc_phi[indexele1]); 
	  
	// 	  (histogsf_gsfet_M120[p])->Fill(gsf_gsfet[indexele2]); 
	// 	  (histogsfsc_eta_M120[p])->Fill(gsfsc_eta[indexele2]); 
	// 	  (histogsfsc_phi_M120[p])->Fill(gsfsc_phi[indexele2]); 
	// 	}

	//	if(!(TString((input[p])->GetName()).Contains("Zee"))){
	if(p==0){
	  std::ostringstream eventinfo;
	  eventinfo<<"sample "<<(input[p])->GetName()<<" eventcounter "<<eventcounter<<" runnumber "<<runnumber<<" eventnumber "<<eventnumber<<" mass "<<mass<<" \n ";
	  eventinfo<<"********************General Event Info********************"<<" \n ";
	  eventinfo<<" ele1 pt "<<gsf_gsfet[indexele1]
		   <<" ele1 eta "<<gsf_eta[indexele1]
		   <<" ele1 phi "<<gsf_phi[indexele1]
		   <<" ele1 e "<<gsf_gsfet[indexele1]*cosh(gsf_eta[indexele1])
		   <<" ele1 charge "<<gsf_charge[indexele1]
		   <<" \n ";
	  
	  eventinfo<<" ele2 pt "<<gsf_gsfet[indexele2]
		   <<" ele2 eta "<<gsf_eta[indexele2]
		   <<" ele2 phi "<<gsf_phi[indexele2]
		   <<" ele2 e "<<gsf_gsfet[indexele2]*cosh(gsf_eta[indexele2])
		   <<" ele2 charge "<<gsf_charge[indexele2]
		   <<" \n ";

	  eventinfo<<" elesc1 pt "<<gsfsc_pt[indexele1]
		   <<" elesc1 eta "<<gsfsc_eta[indexele1]
		   <<" elesc1 phi "<<gsfsc_phi[indexele1]
		   <<" elesc1 e "<<gsfsc_e[indexele1]
		   <<" \n ";

	  eventinfo<<" elesc2 pt "<<gsfsc_pt[indexele2]
		   <<" elesc2 eta "<<gsfsc_eta[indexele2]
		   <<" elesc2 phi "<<gsfsc_phi[indexele2]
		   <<" elesc2 e "<<gsfsc_e[indexele2]
		   <<" \n ";
	  eventinfo<<"****************************************"<<" \n ";

	  //cout<<eventinfo.str()<<endl;
	}

	if(mass > 60. && mass < 120.) nbevents60to120++;
	if(mass > 20.) nbevents20toInf++;
	if(mass > 120.) nbevents120toInf++;
	if(mass > 200.) nbevents200toInf++;
	if(mass > 40.) nbevents40toInf++;
	if(mass > 80. && mass < 100.) nbevents80to100++;

      }
      


      //----------------N-1 distributions----------------------
      std::vector<bool> heepcutsbarrel;
      std::vector<bool> heepcutsendcap;
      //Be careful of the number of cuts considered
      //Update this when necessary
      //heepcuts.reserve(8);
      //for(int cutnumber = 0;cutnumber < heepcuts.size();cutnumber++) {

      for(int cutnumber = 0;cutnumber < 8;cutnumber++) {
	
	indexele1 = -1;
	indexele2 = -1;
      
	for(int i=0;i<gsf_size;i++) {
	  
	  if(fabs(pvz[0]) > 24.) continue;
	  if(gsf_SwissCross[i] > 0.95) continue;
	  if(gsf_gsfet[i] < 25.) continue;
	  isbarrel = fabs(gsfsc_eta[i]) < 1.442;
	  isendcap = (fabs(gsfsc_eta[i]) > 1.56) && (fabs(gsfsc_eta[i]) < 2.5);
	  if(!(isbarrel || isendcap)) continue;
	  
	  //cout<<"before ecaldriven cut"<<endl;

	  if(!gsfpass_ECALDRIVEN[i]) continue;
	  
	  //Init all the HEEP booleans for this GSF
	  heepcutsbarrel.erase(heepcutsbarrel.begin(),heepcutsbarrel.end());
	  heepcutsendcap.erase(heepcutsendcap.begin(),heepcutsendcap.end());

	  heepcutsbarrel.push_back(gsfpass_ISOLEMHADDEPTH1[i]); 
	  heepcutsbarrel.push_back(gsfpass_ISOLHADDEPTH2[i]); 
	  heepcutsbarrel.push_back(gsfpass_ISOLPTTRKS[i]); 
	  heepcutsbarrel.push_back(gsfpass_HADEM[i]); 
	  heepcutsbarrel.push_back(gsfpass_DETAIN[i]); 
	  heepcutsbarrel.push_back(gsfpass_DPHIIN[i]); 
	  heepcutsbarrel.push_back(gsfpass_SIGMAIETAIETA[i]); 
	  heepcutsbarrel.push_back(gsfpass_E2X5OVER5X5[i]);

	  heepcutsendcap.push_back(gsfpass_ISOLEMHADDEPTH1[i]); 
	  heepcutsendcap.push_back(gsfpass_ISOLHADDEPTH2[i]); 
	  heepcutsendcap.push_back(gsfpass_ISOLPTTRKS[i]); 
	  heepcutsendcap.push_back(gsfpass_HADEM[i]); 
	  //heepcutsendcap.push_back(true); 
	  heepcutsendcap.push_back(gsfpass_DETAIN[i]); 
	  heepcutsendcap.push_back(gsfpass_DPHIIN[i]); 
	  heepcutsendcap.push_back(gsfpass_SIGMAIETAIETA[i]); 
	  heepcutsendcap.push_back(gsfpass_E2X5OVER5X5[i]);

	  //cout<<"after heepcuts push back"<<endl;

	  bool heepN1barrel = isbarrel;
	  bool heepN1endcap = isendcap;
	  
	  for(int cutcounter = 0;cutcounter < heepcutsbarrel.size();cutcounter++) {
	    if(cutnumber == cutcounter) continue;
	    //cout<<"cutcounter "<<cutcounter<<endl;
	    heepN1barrel = heepN1barrel && (heepcutsbarrel[cutcounter]);
	    heepN1endcap = heepN1endcap && (heepcutsendcap[cutcounter]);
	  }
	  
	  //cout<<"after HEEPN1barrel"<<endl;

	  if(heepN1barrel || heepN1endcap){
	    if(indexele1 >= 0 && indexele2 == -1) indexele2 = i;
	    if(indexele1 == -1 && indexele2 == -1) indexele1 = i;
	  }
	
	}
	
	if(indexele1 >= 0 && indexele2 >= 0) {
	  //Be careful with the ordering
	  
	  if(cutnumber == 0){(histogsf_ecaliso[p])->Fill(gsf_ecaliso[indexele1]+gsf_hcaliso1[indexele1]);}//ISOLEMHADDEPTH1 
	  if(cutnumber == 1){(histogsf_hcaliso2[p])->Fill(gsf_hcaliso2[indexele1]);}//ISOLHADDEPTH2 
	  if(cutnumber == 2){(histogsf_trackiso[p])->Fill(gsf_trackiso[indexele1]);}//ISOLPTTRKS 
	  if(cutnumber == 3){(histogsf_hovere[p])->Fill(gsf_hovere[indexele1]);}//HADEM 
	  if(cutnumber == 4){(histogsf_deltaeta[p])->Fill(gsf_deltaeta[indexele1]);}//DETAIN 
	  if(cutnumber == 5){(histogsf_deltaphi[p])->Fill(gsf_deltaphi[indexele1]);}//DPHIIN 
	  if(cutnumber == 6){(histogsf_sigmaIetaIeta[p])->Fill(gsf_sigmaIetaIeta[indexele1]);}//SIGMAIETAIETA 
	  if(cutnumber == 7){(histogsf_e2x5overe5x5[p])->Fill(gsf_e2x5overe5x5[indexele1]);}//E2X5OVER5X5
	  
	  if(cutnumber == 0){(histogsf_ecaliso[p])->Fill(gsf_ecaliso[indexele2]+gsf_hcaliso1[indexele2]);}//ISOLEMHADDEPTH1 
	  if(cutnumber == 1){(histogsf_hcaliso2[p])->Fill(gsf_hcaliso2[indexele2]);}//ISOLHADDEPTH2 
	  if(cutnumber == 2){(histogsf_trackiso[p])->Fill(gsf_trackiso[indexele2]);}//ISOLPTTRKS 
	  if(cutnumber == 3){(histogsf_hovere[p])->Fill(gsf_hovere[indexele2]);}//HADEM 
	  if(cutnumber == 4){(histogsf_deltaeta[p])->Fill(gsf_deltaeta[indexele2]);}//DETAIN 
	  if(cutnumber == 5){(histogsf_deltaphi[p])->Fill(gsf_deltaphi[indexele2]);}//DPHIIN 
	  if(cutnumber == 6){(histogsf_sigmaIetaIeta[p])->Fill(gsf_sigmaIetaIeta[indexele2]);}//SIGMAIETAIETA 
	  if(cutnumber == 7){(histogsf_e2x5overe5x5[p])->Fill(gsf_e2x5overe5x5[indexele2]);}//E2X5OVER5X5


	  if(fabs(gsfsc_eta[indexele1]) < 1.442){
	    if(cutnumber == 0){(histogsfbarrel_ecaliso[p])->Fill(gsf_ecaliso[indexele1]+gsf_hcaliso1[indexele1]);}//ISOLEMHADDEPTH1 
	    if(cutnumber == 1){(histogsfbarrel_hcaliso2[p])->Fill(gsf_hcaliso2[indexele1]);}//ISOLHADDEPTH2 
	    if(cutnumber == 2){(histogsfbarrel_trackiso[p])->Fill(gsf_trackiso[indexele1]);}//ISOLPTTRKS 
	    if(cutnumber == 3){(histogsfbarrel_hovere[p])->Fill(gsf_hovere[indexele1]);}//HADEM 
	    if(cutnumber == 4){(histogsfbarrel_deltaeta[p])->Fill(gsf_deltaeta[indexele1]);}//DETAIN 
	    if(cutnumber == 5){(histogsfbarrel_deltaphi[p])->Fill(gsf_deltaphi[indexele1]);}//DPHIIN 
	    if(cutnumber == 6){(histogsfbarrel_sigmaIetaIeta[p])->Fill(gsf_sigmaIetaIeta[indexele1]);}//SIGMAIETAIETA 
	    if(cutnumber == 7){(histogsfbarrel_e2x5overe5x5[p])->Fill(gsf_e2x5overe5x5[indexele1]);}//E2X5OVER5X5
	  }

	  if(fabs(gsfsc_eta[indexele2]) < 1.442){
	    if(cutnumber == 0){(histogsfbarrel_ecaliso[p])->Fill(gsf_ecaliso[indexele2]+gsf_hcaliso1[indexele2]);}//ISOLEMHADDEPTH1 
	    if(cutnumber == 1){(histogsfbarrel_hcaliso2[p])->Fill(gsf_hcaliso2[indexele2]);}//ISOLHADDEPTH2 
	    if(cutnumber == 2){(histogsfbarrel_trackiso[p])->Fill(gsf_trackiso[indexele2]);}//ISOLPTTRKS 
	    if(cutnumber == 3){(histogsfbarrel_hovere[p])->Fill(gsf_hovere[indexele2]);}//HADEM 
	    if(cutnumber == 4){(histogsfbarrel_deltaeta[p])->Fill(gsf_deltaeta[indexele2]);}//DETAIN 
	    if(cutnumber == 5){(histogsfbarrel_deltaphi[p])->Fill(gsf_deltaphi[indexele2]);}//DPHIIN 
	    if(cutnumber == 6){(histogsfbarrel_sigmaIetaIeta[p])->Fill(gsf_sigmaIetaIeta[indexele2]);}//SIGMAIETAIETA 
	    if(cutnumber == 7){(histogsfbarrel_e2x5overe5x5[p])->Fill(gsf_e2x5overe5x5[indexele2]);}//E2X5OVER5X5
	  }

	  if(fabs(gsfsc_eta[indexele1]) > 1.56 && fabs(gsfsc_eta[indexele1]) < 2.5){
	    if(cutnumber == 0){(histogsfendcap_ecaliso[p])->Fill(gsf_ecaliso[indexele1]+gsf_hcaliso1[indexele1]);}//ISOLEMHADDEPTH1 
	    if(cutnumber == 1){(histogsfendcap_hcaliso2[p])->Fill(gsf_hcaliso2[indexele1]);}//ISOLHADDEPTH2 
	    if(cutnumber == 2){(histogsfendcap_trackiso[p])->Fill(gsf_trackiso[indexele1]);}//ISOLPTTRKS 
	    if(cutnumber == 3){(histogsfendcap_hovere[p])->Fill(gsf_hovere[indexele1]);}//HADEM 
	    if(cutnumber == 4){(histogsfendcap_deltaeta[p])->Fill(gsf_deltaeta[indexele1]);}//DETAIN 
	    if(cutnumber == 5){(histogsfendcap_deltaphi[p])->Fill(gsf_deltaphi[indexele1]);}//DPHIIN 
	    if(cutnumber == 6){(histogsfendcap_sigmaIetaIeta[p])->Fill(gsf_sigmaIetaIeta[indexele1]);}//SIGMAIETAIETA 
	    if(cutnumber == 7){(histogsfendcap_e2x5overe5x5[p])->Fill(gsf_e2x5overe5x5[indexele1]);}//E2X5OVER5X5
	  }
	  
	  if(fabs(gsfsc_eta[indexele2]) > 1.56 && fabs(gsfsc_eta[indexele2]) < 2.5){
	    if(cutnumber == 0){(histogsfendcap_ecaliso[p])->Fill(gsf_ecaliso[indexele2]+gsf_hcaliso1[indexele2]);}//ISOLEMHADDEPTH1 
	    if(cutnumber == 1){(histogsfendcap_hcaliso2[p])->Fill(gsf_hcaliso2[indexele2]);}//ISOLHADDEPTH2 
	    if(cutnumber == 2){(histogsfendcap_trackiso[p])->Fill(gsf_trackiso[indexele2]);}//ISOLPTTRKS 
	    if(cutnumber == 3){(histogsfendcap_hovere[p])->Fill(gsf_hovere[indexele2]);}//HADEM 
	    if(cutnumber == 4){(histogsfendcap_deltaeta[p])->Fill(gsf_deltaeta[indexele2]);}//DETAIN 
	    if(cutnumber == 5){(histogsfendcap_deltaphi[p])->Fill(gsf_deltaphi[indexele2]);}//DPHIIN 
	    if(cutnumber == 6){(histogsfendcap_sigmaIetaIeta[p])->Fill(gsf_sigmaIetaIeta[indexele2]);}//SIGMAIETAIETA 
	    if(cutnumber == 7){(histogsfendcap_e2x5overe5x5[p])->Fill(gsf_e2x5overe5x5[indexele2]);}//E2X5OVER5X5
	  }
	  
	}

      }
	
      //----------------N-1 distributions----------------------






      int numtagandprobepairs = 0;
      int numtagandprobepairs80_100 = 0;

      //----------------TAG AND PROBE-------------------
      for(int n=0;n<gsf_size;n++)
	{
	  if(gsf_SwissCross[n] > 0.95) continue;
	  if(gsf_gsfet[n] < 25.) continue;
	  probeisbarrel = fabs(gsfsc_eta[n]) < 1.442;
	  probeisendcap = (fabs(gsfsc_eta[n]) > 1.56) && (fabs(gsfsc_eta[n]) < 2.5);
	  if(!(probeisbarrel || probeisendcap)) continue;

	  if(!gsfpass_ECALDRIVEN[n]) continue;

	  bool probepreselbarrel = probeisbarrel && (gsfsc_pt[n] > 4.&& fabs(gsf_deltaeta[n]) < 0.02 &&  fabs(gsf_deltaphi[n]) < 0.15 && gsf_hovere[n] < 0.1);
	  bool probepreselendcap = probeisendcap && (gsfsc_pt[n] > 4.&& fabs(gsf_deltaphi[n]) < 0.15 && gsf_hovere[n] < 0.1);
      
	  bool gsfpreselection  = (probepreselbarrel || probepreselendcap);

	  if(!gsfpreselection) continue;

	  //Here normally we have a probe
	  //Look for a tag
	  for(int l=0;l<gsf_size;l++)
	    {
	      //Take care not to take the same object
	      if(n == l) continue;

	      if(gsf_gsfet[l] < 25.) continue;

	      tagisbarrel = fabs(gsfsc_eta[l]) < 1.442;
	      tagisendcap = (fabs(gsfsc_eta[l]) > 1.56) && (fabs(gsfsc_eta[l]) < 2.5);
	      if(!(tagisbarrel || tagisendcap)) continue;

	      if(gsf_SwissCross[l] > 0.95) continue;

	      if(!gsfpass_ECALDRIVEN[l]) continue;

	      bool tagselectbarrel = tagisbarrel && (gsfpass_ISO[l] && gsfpass_HADEM[l] && gsfpass_DETAIN[l] && gsfpass_DPHIIN[l] && gsfpass_SIGMAIETAIETA[l] && gsfpass_E2X5OVER5X5[l]);
	      bool tagselectendcap = tagisendcap && (gsfpass_ISO[l] && gsfpass_HADEM[l] && gsfpass_DPHIIN[l] && gsfpass_SIGMAIETAIETA[l] && gsfpass_E2X5OVER5X5[l]);

	      bool tagselect = (tagselectbarrel || tagselectendcap);

	      if(!tagselect) continue;

	      //--------------maybe additional cut on the tag
	      if( (gsfsc_e[l]/gsftrackp[l]) > 1.5) continue;
	      //--------------

	      //Here normally we have a tag and probe pair
	      //Constraint its invariant mass to 80-100

	      numtagandprobepairs++;

	      TLorentzVector probe;
	      TLorentzVector tag;
	      
	      probe.SetPtEtaPhiE(gsf_gsfet[n],gsf_eta[n],gsf_phi[n],(gsf_gsfet[n]*cosh(gsf_eta[n])));
	      tag.SetPtEtaPhiE(gsf_gsfet[l],gsf_eta[l],gsf_phi[l],(gsf_gsfet[l]*cosh(gsf_eta[l])));
	
	      float tagprobemass = (tag+probe).Mag();

	      (tagprobe_mass[p])->Fill(tagprobemass);

	      if(tagprobemass < 80. || tagprobemass > 100.) continue;

	      numtagandprobepairs80_100++;

	      if(TString((input[p])->GetName()).Contains("QCDEmEnriched")) cout<<"tagprobe info "<<(input[p])->GetName()<<" eventnumber "<<eventnumber<<endl;

	      (tagprobe_mass_80_100[p])->Fill(tagprobemass);

	      (CaloMet[p])->Fill(calomet);
	      (Met[p])->Fill(met);

	      float dphitagprobe = gsf_phi[n]-gsf_phi[l];
	      if(dphitagprobe > PI) dphitagprobe = dphitagprobe - 2.*PI;
	      if(dphitagprobe < -PI) dphitagprobe = dphitagprobe + 2.*PI;

	      (tagprobedeltaphi[p])->Fill(dphitagprobe);

	      (NJetsAK[p])->Fill(nJetsAKT_pt15*1.);
	      (NJetsIC[p])->Fill(nJetsIC5_pt15*1.);

	      //Now we have a tag and probe pair with mass in the correct range
	      //!!!! Go for efficiencies

	      (histoprobept[p])->Fill(gsf_gsfet[n]);
	      if(fabs(gsfsc_eta[n]) < 1.442) (histoprobebarrelpt[p])->Fill(gsf_gsfet[n]);
	      if(fabs(gsfsc_eta[n]) < 2.5 && fabs(gsfsc_eta[n]) > 1.56) (histoprobeendcappt[p])->Fill(gsf_gsfet[n]);
	      (histoprobeeta[p])->Fill(gsf_eta[n]);
	      (histoprobesceta[p])->Fill(gsfsc_eta[n]);
	      (histoprobephi[p])->Fill(gsf_phi[n]);
	      (histoprobescphi[p])->Fill(gsfsc_phi[n]);
	      (histoprobeenergy[p])->Fill(gsf_gsfet[n]*cosh(gsf_eta[n]));
	      (histoprobescenergy[p])->Fill(gsfsc_e[n]);
	      (histoprobeEoverP[p])->Fill(gsfsc_e[n]/gsftrackp[n]);

	      (histoprobe_deltaeta[p])->Fill(gsf_deltaeta[n]); 
	      (histoprobe_deltaphi[p])->Fill(gsf_deltaphi[n]); 
	      (histoprobe_hovere[p])->Fill(gsf_hovere[n]); 
	      (histoprobe_trackiso[p])->Fill(gsf_trackiso[n]); 
	      (histoprobe_ecaliso[p])->Fill(gsf_ecaliso[n]); 
	      (histoprobe_hcaliso1[p])->Fill(gsf_hcaliso1[n]); 
	      (histoprobe_hcaliso2[p])->Fill(gsf_hcaliso2[n]); 
	      (histoprobe_sigmaetaeta[p])->Fill(gsf_sigmaetaeta[n]); 
	      (histoprobe_sigmaIetaIeta[p])->Fill(gsf_sigmaIetaIeta[n]); 
	      
	      (histotagpt[p])->Fill(gsf_gsfet[l]);
	      (histotageta[p])->Fill(gsf_eta[l]);
	      (histotagsceta[p])->Fill(gsfsc_eta[l]);
	      (histotagphi[p])->Fill(gsf_phi[l]);
	      (histotagscphi[p])->Fill(gsfsc_phi[l]);
	      (histotagenergy[p])->Fill(gsf_gsfet[l]*cosh(gsf_eta[l]));
	      (histotagscenergy[p])->Fill(gsfsc_e[l]);
	      (histotagEoverP[p])->Fill(gsfsc_e[l]/gsftrackp[l]);

	      (histotag_deltaeta[p])->Fill(gsf_deltaeta[l]); 
	      (histotag_deltaphi[p])->Fill(gsf_deltaphi[l]); 
	      (histotag_hovere[p])->Fill(gsf_hovere[l]); 
	      (histotag_trackiso[p])->Fill(gsf_trackiso[l]); 
	      (histotag_ecaliso[p])->Fill(gsf_ecaliso[l]); 
	      (histotag_hcaliso1[p])->Fill(gsf_hcaliso1[l]); 
	      (histotag_hcaliso2[p])->Fill(gsf_hcaliso2[l]); 
	      (histotag_sigmaetaeta[p])->Fill(gsf_sigmaetaeta[l]); 
	      (histotag_sigmaIetaIeta[p])->Fill(gsf_sigmaIetaIeta[l]); 
	      (histotag_gsfet[p])->Fill(gsf_gsfet[l]); 
	      
	      bool probebarrelisheep = probeisbarrel && (gsfpass_ISO[n] && gsfpass_HADEM[n] && gsfpass_DETAIN[n] && gsfpass_DPHIIN[n] && gsfpass_SIGMAIETAIETA[n] && gsfpass_E2X5OVER5X5[n]);
	      bool probeendapisheep = probeisendcap && (gsfpass_ISO[n] && gsfpass_HADEM[n] && gsfpass_DPHIIN[n] && gsfpass_SIGMAIETAIETA[n] && gsfpass_E2X5OVER5X5[n]);

	      bool probeisheep = (probebarrelisheep || probeendapisheep);

	      //if(gsfpass_ISO[n] && gsfpass_HADEM[n] && gsfpass_SIGMAIETAIETA[n] && gsfpass_E2X5OVER5X5[n] && gsfpass_ECALDRIVEN[n]){
	      if(probeisheep){
		(histoprobeHEEPpt[p])->Fill(gsf_gsfet[n]);
		if(fabs(gsfsc_eta[n]) < 1.442) (histoprobeHEEPbarrelpt[p])->Fill(gsf_gsfet[n]);
		if(fabs(gsfsc_eta[n]) < 2.5 && fabs(gsfsc_eta[n]) > 1.56) (histoprobeHEEPendcappt[p])->Fill(gsf_gsfet[n]);
		(histoprobeHEEPeta[p])->Fill(gsf_eta[n]);
		(histoprobeHEEPsceta[p])->Fill(gsfsc_eta[n]);
		(histoprobeHEEPphi[p])->Fill(gsf_phi[n]);
		(histoprobeHEEPscphi[p])->Fill(gsfsc_phi[n]);
		(histoprobeHEEPenergy[p])->Fill(gsf_gsfet[n]*cosh(gsf_eta[n]));		
		(histoprobeHEEPscenergy[p])->Fill(gsfsc_e[n]);
	      }

	    }

	}
      //----------------END OF TAG AND PROBE-------------------


      (tagprobe_num[p])->Fill(numtagandprobepairs*1.);
      (tagprobe_num_80_100[p])->Fill(numtagandprobepairs80_100*1.);
      
    }//end of loop over entries
    
    std::ostringstream nbeventinfo;

    nbeventinfo<<" ------------------------------nbevents with M > 20 "<<nbevents20toInf<<" \n";
    if(p!=0) nbeventinfo<<" and integral "<<(nbevents20toInf*weight[p]*gridjobeffi[p]*LumiFactor)<<" \n ";

    nbeventinfo<<" ------------------------------nbevents with M > 120 "<<nbevents120toInf<<" \n";
    if(p!=0) nbeventinfo<<" and integral "<<(nbevents120toInf*weight[p]*gridjobeffi[p]*LumiFactor)<<" \n ";

    nbeventinfo<<" ------------------------------nbevents with M > 200 "<<nbevents200toInf<<" \n";
    if(p!=0) nbeventinfo<<" and integral "<<(nbevents200toInf*weight[p]*gridjobeffi[p]*LumiFactor)<<" \n ";

    nbeventinfo<<" ------------------------------nbevents with M > 40 "<<nbevents40toInf<<" \n";
    if(p!=0) nbeventinfo<<" and integral "<<(nbevents40toInf*weight[p]*gridjobeffi[p]*LumiFactor)<<" \n ";

    nbeventinfo<<" ------------------------------nbevents with 80 < M < 100 "<<nbevents80to100<<" \n";
    if(p!=0) nbeventinfo<<" and integral "<<(nbevents80to100*weight[p]*gridjobeffi[p]*LumiFactor)<<" \n ";

    nbeventinfo<<" ------------------------------nbevents with 60 < M < 120 "<<nbevents60to120<<" \n";
    if(p!=0) nbeventinfo<<" and integral "<<(nbevents60to120*weight[p]*gridjobeffi[p]*LumiFactor)<<" \n ";

    cout<<nbeventinfo.str()<<endl;

    //histos tag and probe

    (histoEFFIpt[p])->Divide((histoprobeHEEPpt[p]),(histoprobept[p]),1.,1.,"B");
    (histoEFFIeta[p])->Divide((histoprobeHEEPeta[p]),(histoprobeeta[p]),1.,1.,"B");
    (histoEFFIphi[p])->Divide((histoprobeHEEPphi[p]),(histoprobephi[p]),1.,1.,"B");
    (histoEFFIenergy[p])->Divide((histoprobeHEEPenergy[p]),(histoprobeenergy[p]),1.,1.,"B");

  }//end of loop over pfile


  //   cout<<"draw histos"<<endl;
  // //   DrawHistosBis(ee_mass,weight, gridjobeffi, LumiFactor, "mass",TString("M (GeV/c^{2})"));
  // //   DrawHistosBis(gsfgsf_mass,weight, gridjobeffi, LumiFactor, "gsfgsfmass",TString("M (GeV/c^{2})"));
  // //   DrawHistosBis(tagprobe_mass, weight, gridjobeffi, LumiFactor,"TPMass",TString("Tag and Probe Mass (GeV/c^{2})"));
  // //   DrawHistosBis(tagprobe_mass_80_100, weight, gridjobeffi, LumiFactor,"TPMass80100",TString("Tag and Probe Mass 80 100 (GeV/c^{2})"));


  //   DrawHistos(ee_mass,weight, gridjobeffi, LumiFactor, "mass",TString("M (GeV/c^{2})"),legendname);
  //   DrawHistos(ee_massSameSign,weight, gridjobeffi, LumiFactor, "mass_SameSign_",TString("M (GeV/c^{2}) (Same Sign)"),legendname);
  //   DrawHistos(ee_massOppSign,weight, gridjobeffi, LumiFactor, "mass_OppSign_",TString("M (GeV/c^{2}) (Opp Sign)"),legendname);
  //   //DrawHistos(gsfgsf_mass,weight, gridjobeffi, LumiFactor, "gsfgsfmass",TString("M (GeV/c^{2})"),legendname);

  //   //N-1 distributions
  //   DrawHistos(histogsf_deltaeta,weight, gridjobeffi, LumiFactor, "deltaeta",TString("#Delta#eta"),legendname);
  //   DrawHistos(histogsf_deltaphi,weight, gridjobeffi, LumiFactor,"deltaphi",TString("#Delta#phi (rad)"),legendname);
  //   DrawHistos(histogsf_hovere,weight, gridjobeffi, LumiFactor,"hovere",TString("H/E"),legendname);
  //   DrawHistos(histogsf_trackiso,weight, gridjobeffi, LumiFactor,"trackiso",TString("Track Iso. (GeV/c)"),legendname);
  //   DrawHistos(histogsf_ecaliso,weight, gridjobeffi, LumiFactor,"ecaliso",TString("Ecal+Hcal1 Iso. (GeV)"),legendname);
  //   //DrawHistos(histogsf_hcaliso1,weight, gridjobeffi, LumiFactor,"hcaliso1",TString("Hcal Iso. 1 (GeV)"),legendname);
  //   DrawHistos(histogsf_hcaliso2,weight, gridjobeffi, LumiFactor,"hcaliso2",TString("Hcal Iso. 2 (GeV)"),legendname);
  //   //DrawHistos(histogsf_sigmaetaeta,weight, gridjobeffi, LumiFactor,"sigmaetaeta",TString("#sigma_{#eta#eta}"),legendname);
  //   DrawHistos(histogsf_sigmaIetaIeta,weight, gridjobeffi, LumiFactor,"sigmaIetaIeta",TString("#sigma_{i#etai#eta}"),legendname);
  //   DrawHistos(histogsf_gsfet,weight, gridjobeffi, LumiFactor,"gsfet",TString("E_{t} (GeV)"),legendname);
  //   DrawHistos(histogsf_pt,weight, gridjobeffi, LumiFactor,"gsfpt",TString("p_{t} (GeV/c)"),legendname);
  //   DrawHistos(histogsf_eta,weight, gridjobeffi, LumiFactor,"gsfeta",TString("#eta"),legendname);
  //   DrawHistos(histogsf_phi,weight, gridjobeffi, LumiFactor,"gsfphi",TString("#phi (rad)"),legendname);
  //   DrawHistos(histogsfsc_eta,weight, gridjobeffi, LumiFactor,"gsfsceta",TString("#eta_{SC}"),legendname);
  //   DrawHistos(histogsfsc_phi,weight, gridjobeffi, LumiFactor,"gsfscphi",TString("#phi_{SC} (rad)"),legendname);
  //   //N-1 distributions

  // //   DrawHistos(histogsfgsf_et,weight, gridjobeffi, LumiFactor,"gsfgsfet",TString("E_{t}"),legendname);
  // //   DrawHistos(histogsfgsf_pt,weight, gridjobeffi, LumiFactor,"gsfgsfpt",TString("p_{t}"),legendname);
  // //   DrawHistos(histogsfgsf_eta,weight, gridjobeffi, LumiFactor,"gsfgsfeta",TString("#eta"),legendname);
  // //   DrawHistos(histogsfgsf_phi,weight, gridjobeffi, LumiFactor,"gsfgsfphi",TString("#phi"),legendname);

  //   DrawHistos(tagprobe_mass, weight, gridjobeffi, LumiFactor,"TPMass",TString("Tag and Probe Mass (GeV/c^{2})"),legendname);
  //   DrawHistos(tagprobe_mass_80_100, weight, gridjobeffi, LumiFactor,"TPMass80100",TString("Tag and Probe Mass 80 100 (GeV/c^{2})"),legendname);

  // //   DrawHistos(CaloMet, weight, gridjobeffi, LumiFactor,"CaloMet",TString("CaloMet"),legendname);
  // //   DrawHistos(Met, weight, gridjobeffi, LumiFactor,"Met",TString("Met"),legendname);
  // //   DrawHistos(tagprobedeltaphi, weight, gridjobeffi, LumiFactor,"tagprobedeltaphi",TString("tagprobedeltaphi"),legendname);
  // //   DrawHistos(NJetsAK, weight, gridjobeffi, LumiFactor,"NJetsAK",TString("NJetsAK"),legendname);
  // //   DrawHistos(NJetsIC, weight, gridjobeffi, LumiFactor,"NJetsIC",TString("NJetsIC"),legendname);

  // //   DrawHistos(histoprobept, weight, gridjobeffi, LumiFactor,"ProbePt",TString("ProbePt"),legendname);
  // //   DrawHistos(histoprobebarrelpt, weight, gridjobeffi, LumiFactor,"Probebarrelpt",TString("Probebarrelpt"),legendname);
  // //   DrawHistos(histoprobeendcappt, weight, gridjobeffi, LumiFactor,"Probeendcappt",TString("Probeendcappt"),legendname);
  // //   DrawHistos(histoprobeeta, weight, gridjobeffi, LumiFactor,"ProbeEta",TString("ProbeEta"),legendname);
  // //   DrawHistos(histoprobesceta, weight, gridjobeffi, LumiFactor,"ProbeScEta",TString("Probe Sc Eta"),legendname);
  // //   DrawHistos(histoprobephi, weight, gridjobeffi, LumiFactor,"ProbePhi",TString("ProbePhi"),legendname);
  // //   DrawHistos(histoprobescphi, weight, gridjobeffi, LumiFactor,"ProbeScPhi",TString("Probe Sc Phi"),legendname);
  // //   DrawHistos(histoprobeenergy, weight, gridjobeffi, LumiFactor,"ProbeEnergy",TString("ProbeEnergy"),legendname);
  // //   DrawHistos(histoprobescenergy, weight, gridjobeffi, LumiFactor,"ProbeScEnergy",TString("Probe Sc Energy"),legendname);
  // //   DrawHistos(histoprobeEoverP, weight, gridjobeffi, LumiFactor,"ProbeEoverP",TString("ProbeEoverP"),legendname);

  // //   DrawHistos(histoprobe_deltaeta,weight, gridjobeffi, LumiFactor, "probedeltaeta",TString("#Delta#eta"),legendname);
  // //   DrawHistos(histoprobe_deltaphi,weight, gridjobeffi, LumiFactor,"probedeltaphi",TString("#Delta#phi"),legendname);
  // //   DrawHistos(histoprobe_hovere,weight, gridjobeffi, LumiFactor,"probehovere",TString("H/E"),legendname);
  // //   DrawHistos(histoprobe_trackiso,weight, gridjobeffi, LumiFactor,"probetrackiso",TString("Track Iso."),legendname);
  // //   DrawHistos(histoprobe_ecaliso,weight, gridjobeffi, LumiFactor,"probeecaliso",TString("Ecal Iso."),legendname);
  // //   DrawHistos(histoprobe_hcaliso1,weight, gridjobeffi, LumiFactor,"probehcaliso1",TString("Hcal Iso. 1"),legendname);
  // //   DrawHistos(histoprobe_hcaliso2,weight, gridjobeffi, LumiFactor,"probehcaliso2",TString("Hcal Iso. 2"),legendname);
  // //   DrawHistos(histoprobe_sigmaetaeta,weight, gridjobeffi, LumiFactor,"probesigmaetaeta",TString("#sigma_{#eta#eta}"),legendname);
  // //   DrawHistos(histoprobe_sigmaIetaIeta,weight, gridjobeffi, LumiFactor,"probesigmaIetaIeta",TString("#sigma_{i#etai#eta}"),legendname);

  // //   DrawHistos(histotagpt, weight, gridjobeffi, LumiFactor,"TagPt",TString("TagPt"),legendname);
  // //   DrawHistos(histotageta, weight, gridjobeffi, LumiFactor,"TagEta",TString("TagEta"),legendname);
  // //   DrawHistos(histotagsceta, weight, gridjobeffi, LumiFactor,"TagScEta",TString("Tag Sc Eta"),legendname);
  // //   DrawHistos(histotagphi, weight, gridjobeffi, LumiFactor,"TagPhi",TString("TagPhi"),legendname);
  // //   DrawHistos(histotagscphi, weight, gridjobeffi, LumiFactor,"TagScPhi",TString("Tag Sc Phi"),legendname);
  // //   DrawHistos(histotagenergy, weight, gridjobeffi, LumiFactor,"TagEnergy",TString("TagEnergy"),legendname);
  // //   DrawHistos(histotagscenergy, weight, gridjobeffi, LumiFactor,"TagScEnergy",TString("Tag Sc Energy"),legendname);
  // //   DrawHistos(histotagEoverP, weight, gridjobeffi, LumiFactor,"TagEoverP",TString("TagEoverP"),legendname);

  // //   DrawHistos(histotag_deltaeta,weight, gridjobeffi, LumiFactor, "tagdeltaeta",TString("#Delta#eta"),legendname);
  // //   DrawHistos(histotag_deltaphi,weight, gridjobeffi, LumiFactor,"tagdeltaphi",TString("#Delta#phi"),legendname);
  // //   DrawHistos(histotag_hovere,weight, gridjobeffi, LumiFactor,"taghovere",TString("H/E"),legendname);
  // //   DrawHistos(histotag_trackiso,weight, gridjobeffi, LumiFactor,"tagtrackiso",TString("Track Iso."),legendname);
  // //   DrawHistos(histotag_ecaliso,weight, gridjobeffi, LumiFactor,"tagecaliso",TString("Ecal Iso."),legendname);
  // //   DrawHistos(histotag_hcaliso1,weight, gridjobeffi, LumiFactor,"taghcaliso1",TString("Hcal Iso. 1"),legendname);
  // //   DrawHistos(histotag_hcaliso2,weight, gridjobeffi, LumiFactor,"taghcaliso2",TString("Hcal Iso. 2"),legendname);
  // //   DrawHistos(histotag_sigmaetaeta,weight, gridjobeffi, LumiFactor,"tagsigmaetaeta",TString("#sigma_{#eta#eta}"),legendname);
  // //   DrawHistos(histotag_sigmaIetaIeta,weight, gridjobeffi, LumiFactor,"tagsigmaIetaIeta",TString("#sigma_{i#etai#eta}"),legendname);
  // //   DrawHistos(histotag_gsfet,weight, gridjobeffi, LumiFactor,"taggsfet",TString("E_{t}"),legendname);


  // //   DrawHistos(histoprobeHEEPpt, weight, gridjobeffi, LumiFactor,"ProbeHEEPPt",TString("ProbeHEEPPt"),legendname);
  // //   DrawHistos(histoprobeHEEPbarrelpt, weight, gridjobeffi, LumiFactor,"ProbeHEEPBarrelPt",TString("ProbeHEEPBarrelPt"),legendname);
  // //   DrawHistos(histoprobeHEEPendcappt, weight, gridjobeffi, LumiFactor,"ProbeHEEPEndcapPt",TString("ProbeHEEPEndcapPt"),legendname);
  // //   DrawHistos(histoprobeHEEPeta, weight, gridjobeffi, LumiFactor,"ProbeHEEPEta",TString("ProbeHEEPEta"),legendname);
  // //   DrawHistos(histoprobeHEEPsceta, weight, gridjobeffi, LumiFactor,"ProbeHEEPScEta",TString("ProbeHEEP Sc Eta"),legendname);
  // //   DrawHistos(histoprobeHEEPphi, weight, gridjobeffi, LumiFactor,"ProbeHEEPPhi",TString("ProbeHEEPPhi"),legendname);
  // //   DrawHistos(histoprobeHEEPscphi, weight, gridjobeffi, LumiFactor,"ProbeHEEPScPhi",TString("ProbeHEEP Sc Phi"),legendname);
  // //   DrawHistos(histoprobeHEEPenergy, weight, gridjobeffi, LumiFactor,"ProbeHEEPEnergy",TString("ProbeHEEPEnergy"),legendname);
  // //   DrawHistos(histoprobeHEEPscenergy, weight, gridjobeffi, LumiFactor,"ProbeHEEPScEnergy",TString("ProbeHEEP Sc Energy"),legendname);


  // //   //   DrawTagAndProbeHistos(histoEFFIenergy,weight, gridjobeffi, LumiFactor,"gsfgsfet",TString("E"));
  // //   //   DrawTagAndProbeHistos(histoEFFIpt,weight, gridjobeffi, LumiFactor,"gsfgsfpt",TString("p_{t}"));
  // //   //   DrawTagAndProbeHistos(histoEFFIeta,weight, gridjobeffi, LumiFactor,"gsfgsfeta",TString("#eta"));
  // //   //   DrawTagAndProbeHistos(histoEFFIphi,weight, gridjobeffi, LumiFactor,"gsfgsfphi",TString("#phi"));




  cout<<"draw histos in category"<<endl;

  DrawHistosCategory(histopvx,weight, gridjobeffi, LumiFactor, "pvxCategory",TString("x_{PV} (cm)"),TString(""),legendname);
  DrawHistosCategory(histopvy,weight, gridjobeffi, LumiFactor, "pvyCategory",TString("y_{PV} (cm)"),TString(""),legendname);
  DrawHistosCategory(histopvz,weight, gridjobeffi, LumiFactor, "pvzCategory",TString("z_{PV} (cm)"),TString(""),legendname);
  DrawHistosCategory(histopvsize,weight, gridjobeffi, LumiFactor, "pvsizeCategory",TString("#PV"),TString(""),legendname);

  DrawHistosCategory(Zbosonpt,weight, gridjobeffi, LumiFactor, "ZbosonptCategory",TString("p_{t,Z} (GeV/c)"),TString("# Events / 5 GeV/c"),legendname);
  DrawHistosCategory(ZbosonptBBandEB,weight, gridjobeffi, LumiFactor, "ZbosonptCategoryBBandEB",TString("p_{t,Z} (GeV/c)"),TString("# Events / 5 GeV/c"),legendname);
  DrawHistosCategory(ZbosonptBB,weight, gridjobeffi, LumiFactor, "ZbosonptCategoryBB",TString("p_{t,Z} (GeV/c)"),TString("# Events / 5 GeV/c"),legendname);
  DrawHistosCategory(ZbosonptEB,weight, gridjobeffi, LumiFactor, "ZbosonptCategoryEB",TString("p_{t,Z} (GeV/c)"),TString("# Events / 5 GeV/c"),legendname);
  DrawHistosCategory(ZbosonptEE,weight, gridjobeffi, LumiFactor, "ZbosonptCategoryEE",TString("p_{t,Z} (GeV/c)"),TString("# Events / 5 GeV/c"),legendname);
  DrawHistosCategory(ZbosonptSameSign,weight, gridjobeffi, LumiFactor, "Zbosonpt_SameSign_Category",TString("p_{t,Z} (GeV/c) (Same Sign)"),TString("# Events / 5 GeV/c"),legendname);
  DrawHistosCategory(ZbosonptOppSign,weight, gridjobeffi, LumiFactor, "Zbosonpt_OppSign_Category",TString("p_{t,Z} (GeV/c) (Opp Sign)"),TString("# Events / 5 GeV/c"),legendname);

  DrawHistosCategory(ZbosonptZoom,weight, gridjobeffi, LumiFactor, "ZbosonptZoomCategory",TString("p_{t,Z} (GeV/c)"),TString("# Events / 1 GeV/c"),legendname);
  DrawHistosCategory(ZbosonptZoomBBandEB,weight, gridjobeffi, LumiFactor, "ZbosonptZoomCategoryBBandEB",TString("p_{t,Z} (GeV/c)"),TString("# Events / 1 GeV/c"),legendname);
  DrawHistosCategory(ZbosonptZoomBB,weight, gridjobeffi, LumiFactor, "ZbosonptZoomCategoryBB",TString("p_{t,Z} (GeV/c)"),TString("# Events / 1 GeV/c"),legendname);
  DrawHistosCategory(ZbosonptZoomEB,weight, gridjobeffi, LumiFactor, "ZbosonptZoomCategoryEB",TString("p_{t,Z} (GeV/c)"),TString("# Events / 1 GeV/c"),legendname);
  DrawHistosCategory(ZbosonptZoomEE,weight, gridjobeffi, LumiFactor, "ZbosonptZoomCategoryEE",TString("p_{t,Z} (GeV/c)"),TString("# Events / 1 GeV/c"),legendname);
  DrawHistosCategory(ZbosonptZoomSameSign,weight, gridjobeffi, LumiFactor, "ZbosonptZoom_SameSign_Category",TString("p_{t,Z} (GeV/c) (Same Sign)"),TString("# Events / 1 GeV/c"),legendname);
  DrawHistosCategory(ZbosonptZoomOppSign,weight, gridjobeffi, LumiFactor, "ZbosonptZoom_OppSign_Category",TString("p_{t,Z} (GeV/c) (Opp Sign)"),TString("# Events / 1 GeV/c"),legendname);

  DrawHistosCategory(Zbosonpz,weight, gridjobeffi, LumiFactor, "ZbosonpzCategory",TString("p_{z,Z} (GeV/c)"),TString("# Events / 20 GeV/c"),legendname);
  DrawHistosCategory(ZbosonpzBBandEB,weight, gridjobeffi, LumiFactor, "ZbosonpzCategoryBBandEB",TString("p_{z,Z} (GeV/c)"),TString("# Events / 20 GeV/c"),legendname);
  DrawHistosCategory(ZbosonpzBB,weight, gridjobeffi, LumiFactor, "ZbosonpzCategoryBB",TString("p_{z,Z} (GeV/c)"),TString("# Events / 20 GeV/c"),legendname);
  DrawHistosCategory(ZbosonpzEB,weight, gridjobeffi, LumiFactor, "ZbosonpzCategoryEB",TString("p_{z,Z} (GeV/c)"),TString("# Events / 20 GeV/c"),legendname);
  DrawHistosCategory(ZbosonpzEE,weight, gridjobeffi, LumiFactor, "ZbosonpzCategoryEE",TString("p_{z,Z} (GeV/c)"),TString("# Events / 20 GeV/c"),legendname);
  DrawHistosCategory(ZbosonpzSameSign,weight, gridjobeffi, LumiFactor, "Zbosonpz_SameSign_Category",TString("p_{z,Z} (GeV/c) (Same Sign)"),TString("# Events / 20 GeV/c"),legendname);
  DrawHistosCategory(ZbosonpzOppSign,weight, gridjobeffi, LumiFactor, "Zbosonpz_OppSign_Category",TString("p_{z,Z} (GeV/c) (Opp Sign)"),TString("# Events / 20 GeV/c"),legendname);

  DrawHistosCategory(ee_mass,weight, gridjobeffi, LumiFactor, "massCategory",TString("M (GeV/c^{2})"),TString("# Events / 5 GeV/c^{2}"),legendname);
  DrawHistosCategory(ee_massBBandEB,weight, gridjobeffi, LumiFactor, "massCategoryBBandEB",TString("M (GeV/c^{2})"),TString("# Events / 5 GeV/c^{2}"),legendname);
  DrawHistosCategory(ee_massBB,weight, gridjobeffi, LumiFactor, "massCategoryBB",TString("M (GeV/c^{2})"),TString("# Events / 5 GeV/c^{2}"),legendname);
  DrawHistosCategory(ee_massEB,weight, gridjobeffi, LumiFactor, "massCategoryEB",TString("M (GeV/c^{2})"),TString("# Events / 5 GeV/c^{2}"),legendname);
  DrawHistosCategory(ee_massEE,weight, gridjobeffi, LumiFactor, "massCategoryEE",TString("M (GeV/c^{2})"),TString("# Events / 5 GeV/c^{2}"),legendname);
  DrawHistosCategory(ee_massSameSign,weight, gridjobeffi, LumiFactor, "mass_SameSign_Category",TString("M (GeV/c^{2}) (Same Sign)"),TString("# Events / 5 GeV/c^{2}"),legendname);
  DrawHistosCategory(ee_massOppSign,weight, gridjobeffi, LumiFactor, "mass_OppSign_Category",TString("M (GeV/c^{2}) (Opp Sign)"),TString("# Events / 5 GeV/c^{2}"),legendname);
  //DrawHistosCategory(gsfgsf_mass,weight, gridjobeffi, LumiFactor, "gsfgsfmassCategory",TString("M (GeV/c^{2})"),TString("# Events / 5 GeV/c^{2}"),legendname);

  DrawHistosCategory(ee_massZoom,weight, gridjobeffi, LumiFactor, "massZoomCategory",TString("M (GeV/c^{2})"),TString("# Events / 1 GeV/c^{2}"),legendname);
  DrawHistosCategory(ee_massZoomBBandEB,weight, gridjobeffi, LumiFactor, "massZoomCategoryBBandEB",TString("M (GeV/c^{2})"),TString("# Events / 1 GeV/c^{2}"),legendname);
  DrawHistosCategory(ee_massZoomBB,weight, gridjobeffi, LumiFactor, "massZoomCategoryBB",TString("M (GeV/c^{2})"),TString("# Events / 1 GeV/c^{2}"),legendname);
  DrawHistosCategory(ee_massZoomEB,weight, gridjobeffi, LumiFactor, "massZoomCategoryEB",TString("M (GeV/c^{2})"),TString("# Events / 1 GeV/c^{2}"),legendname);
  DrawHistosCategory(ee_massZoomEE,weight, gridjobeffi, LumiFactor, "massZoomCategoryEE",TString("M (GeV/c^{2})"),TString("# Events / 1 GeV/c^{2}"),legendname);
  DrawHistosCategory(ee_massZoomSameSign,weight, gridjobeffi, LumiFactor, "massZoom_SameSign_Category",TString("M (GeV/c^{2}) (Same Sign)"),TString("# Events / 1 GeV/c^{2}"),legendname);
  DrawHistosCategory(ee_massZoomOppSign,weight, gridjobeffi, LumiFactor, "massZoom_OppSign_Category",TString("M (GeV/c^{2}) (Opp Sign)"),TString("# Events / 1 GeV/c^{2}"),legendname);

  //N-1 distributions
  DrawHistosCategory(histogsf_deltaeta,weight, gridjobeffi, LumiFactor, "deltaetaCategory",TString("#Delta#eta"),TString(""),legendname);
  DrawHistosCategory(histogsf_deltaphi,weight, gridjobeffi, LumiFactor,"deltaphiCategory",TString("#Delta#phi (rad)"),TString(""),legendname);
  DrawHistosCategory(histogsf_hovere,weight, gridjobeffi, LumiFactor,"hovereCategory",TString("H/E"),TString(""),legendname);
  DrawHistosCategory(histogsf_trackiso,weight, gridjobeffi, LumiFactor,"trackisoCategory",TString("Track Iso. (GeV/c)"),TString(""),legendname);
  DrawHistosCategory(histogsf_ecaliso,weight, gridjobeffi, LumiFactor,"ecalisoCategory",TString("Ecal+Hcal1 Iso. (GeV)"),TString(""),legendname);
  //DrawHistosCategory(histogsf_hcaliso1,weight, gridjobeffi, LumiFactor,"hcaliso1Category",TString("Hcal Iso. 1 (GeV)"),TString(""),legendname);
  DrawHistosCategory(histogsf_hcaliso2,weight, gridjobeffi, LumiFactor,"hcaliso2Category",TString("Hcal Iso. 2 (GeV)"),TString(""),legendname);
  //DrawHistosCategory(histogsf_sigmaetaeta,weight, gridjobeffi, LumiFactor,"sigmaetaetaCategory",TString("#sigma_{#eta#eta}"),TString(""),legendname);
  DrawHistosCategory(histogsf_sigmaIetaIeta,weight, gridjobeffi, LumiFactor,"sigmaIetaIetaCategory",TString("#sigma_{i#etai#eta}"),TString(""),legendname);
  DrawHistosCategory(histogsf_e2x5overe5x5,weight, gridjobeffi, LumiFactor,"e2x5overe5x5Category",TString("E2X5/E5X5"),TString(""),legendname);

  DrawHistosCategory(histogsfbarrel_deltaeta,weight, gridjobeffi, LumiFactor, "deltaetaBarrelCategory",TString("#Delta#eta"),TString(""),legendname);
  DrawHistosCategory(histogsfbarrel_deltaphi,weight, gridjobeffi, LumiFactor,"deltaphiBarrelCategory",TString("#Delta#phi (rad)"),TString(""),legendname);
  DrawHistosCategory(histogsfbarrel_hovere,weight, gridjobeffi, LumiFactor,"hovereBarrelCategory",TString("H/E"),TString(""),legendname);
  DrawHistosCategory(histogsfbarrel_trackiso,weight, gridjobeffi, LumiFactor,"trackisoBarrelCategory",TString("Track Iso. (GeV/c)"),TString(""),legendname);
  DrawHistosCategory(histogsfbarrel_ecaliso,weight, gridjobeffi, LumiFactor,"ecalisoBarrelCategory",TString("Ecal+Hcal1 Iso. (GeV)"),TString(""),legendname);
  //DrawHistosCategory(histogsfbarrel_hcaliso1,weight, gridjobeffi, LumiFactor,"hcaliso1BarrelCategory",TString("Hcal Iso. 1 (GeV)"),TString(""),legendname);
  DrawHistosCategory(histogsfbarrel_hcaliso2,weight, gridjobeffi, LumiFactor,"hcaliso2BarrelCategory",TString("Hcal Iso. 2 (GeV)"),TString(""),legendname);
  //DrawHistosCategory(histogsfbarrel_sigmaetaeta,weight, gridjobeffi, LumiFactor,"sigmaetaetaBarrelCategory",TString("#sigma_{#eta#eta}"),TString(""),legendname);
  DrawHistosCategory(histogsfbarrel_sigmaIetaIeta,weight, gridjobeffi, LumiFactor,"sigmaIetaIetaBarrelCategory",TString("#sigma_{i#etai#eta}"),TString(""),legendname);
  DrawHistosCategory(histogsfbarrel_e2x5overe5x5,weight, gridjobeffi, LumiFactor,"e2x5overe5x5BarrelCategory",TString("E2X5/E5X5"),TString(""),legendname);

  DrawHistosCategory(histogsfendcap_deltaeta,weight, gridjobeffi, LumiFactor, "deltaetaEndcapCategory",TString("#Delta#eta"),TString(""),legendname);
  DrawHistosCategory(histogsfendcap_deltaphi,weight, gridjobeffi, LumiFactor,"deltaphiEndcapCategory",TString("#Delta#phi (rad)"),TString(""),legendname);
  DrawHistosCategory(histogsfendcap_hovere,weight, gridjobeffi, LumiFactor,"hovereEndcapCategory",TString("H/E"),TString(""),legendname);
  DrawHistosCategory(histogsfendcap_trackiso,weight, gridjobeffi, LumiFactor,"trackisoEndcapCategory",TString("Track Iso. (GeV/c)"),TString(""),legendname);
  DrawHistosCategory(histogsfendcap_ecaliso,weight, gridjobeffi, LumiFactor,"ecalisoEndcapCategory",TString("Ecal+Hcal1 Iso. (GeV)"),TString(""),legendname);
  //DrawHistosCategory(histogsfendcap_hcaliso1,weight, gridjobeffi, LumiFactor,"hcaliso1EndcapCategory",TString("Hcal Iso. 1 (GeV)"),TString(""),legendname);
  DrawHistosCategory(histogsfendcap_hcaliso2,weight, gridjobeffi, LumiFactor,"hcaliso2EndcapCategory",TString("Hcal Iso. 2 (GeV)"),TString(""),legendname);
  //DrawHistosCategory(histogsfendcap_sigmaetaeta,weight, gridjobeffi, LumiFactor,"sigmaetaetaEndcapCategory",TString("#sigma_{#eta#eta}"),TString(""),legendname);
  DrawHistosCategory(histogsfendcap_sigmaIetaIeta,weight, gridjobeffi, LumiFactor,"sigmaIetaIetaEndcapCategory",TString("#sigma_{i#etai#eta}"),TString(""),legendname);
  DrawHistosCategory(histogsfendcap_e2x5overe5x5,weight, gridjobeffi, LumiFactor,"e2x5overe5x5EndcapCategory",TString("E2X5/E5X5"),TString(""),legendname);

  DrawHistosCategory(histogsf_gsfet,weight, gridjobeffi, LumiFactor,"gsfetCategory",TString("E_{t} (GeV)"),TString(""),legendname);
  DrawHistosCategory(histogsf_pt,weight, gridjobeffi, LumiFactor,"gsfptCategory",TString("p_{t} (GeV/c)"),TString(""),legendname);
  DrawHistosCategory(histogsf_eta,weight, gridjobeffi, LumiFactor,"gsfetaCategory",TString("#eta"),TString(""),legendname);
  DrawHistosCategory(histogsf_phi,weight, gridjobeffi, LumiFactor,"gsfphiCategory",TString("#phi (rad)"),TString(""),legendname);
  DrawHistosCategory(histogsfsc_eta,weight, gridjobeffi, LumiFactor,"gsfscetaCategory",TString("#eta_{SC}"),TString(""),legendname);
  DrawHistosCategory(histogsfsc_phi,weight, gridjobeffi, LumiFactor,"gsfscphiCategory",TString("#phi_{SC} (rad)"),TString(""),legendname);

  DrawHistosCategory(histogsfbarrel_gsfet,weight, gridjobeffi, LumiFactor,"gsfetBarrelCategory",TString("E_{t} (GeV)"),TString(""),legendname);
  DrawHistosCategory(histogsfbarrel_pt,weight, gridjobeffi, LumiFactor,"gsfptBarrelCategory",TString("p_{t} (GeV/c)"),TString(""),legendname);
  DrawHistosCategory(histogsfbarrel_eta,weight, gridjobeffi, LumiFactor,"gsfetaBarrelCategory",TString("#eta"),TString(""),legendname);
  DrawHistosCategory(histogsfbarrel_phi,weight, gridjobeffi, LumiFactor,"gsfphiBarrelCategory",TString("#phi (rad)"),TString(""),legendname);
  DrawHistosCategory(histogsfscbarrel_eta,weight, gridjobeffi, LumiFactor,"gsfscetaBarrelCategory",TString("#eta_{SC}"),TString(""),legendname);
  DrawHistosCategory(histogsfscbarrel_phi,weight, gridjobeffi, LumiFactor,"gsfscphiBarrelCategory",TString("#phi_{SC} (rad)"),TString(""),legendname);

  DrawHistosCategory(histogsfendcap_gsfet,weight, gridjobeffi, LumiFactor,"gsfetEndcapCategory",TString("E_{t} (GeV)"),TString(""),legendname);
  DrawHistosCategory(histogsfendcap_pt,weight, gridjobeffi, LumiFactor,"gsfptEndcapCategory",TString("p_{t} (GeV/c)"),TString(""),legendname);
  DrawHistosCategory(histogsfendcap_eta,weight, gridjobeffi, LumiFactor,"gsfetaEndcapCategory",TString("#eta"),TString(""),legendname);
  DrawHistosCategory(histogsfendcap_phi,weight, gridjobeffi, LumiFactor,"gsfphiEndcapCategory",TString("#phi (rad)"),TString(""),legendname);
  DrawHistosCategory(histogsfscendcap_eta,weight, gridjobeffi, LumiFactor,"gsfscetaEndcapCategory",TString("#eta_{SC}"),TString(""),legendname);
  DrawHistosCategory(histogsfscendcap_phi,weight, gridjobeffi, LumiFactor,"gsfscphiEndcapCategory",TString("#phi_{SC} (rad)"),TString(""),legendname);

  //   DrawHistosCategory(histogsf_gsfet_M120,weight, gridjobeffi, LumiFactor,"gsfetCategory_M120",TString("E_{t} (GeV)"),TString(""),legendname);
  //   DrawHistosCategory(histogsfsc_eta_M120,weight, gridjobeffi, LumiFactor,"gsfscetaCategory_M120",TString("#eta_{SC}"),TString(""),legendname);
  //   DrawHistosCategory(histogsfsc_phi_M120,weight, gridjobeffi, LumiFactor,"gsfscphiCategory_M120",TString("#phi_{SC} (rad)"),TString(""),legendname);
  //N-1 distributions





  // //   DrawHistosCategory(histogsfgsf_et,weight, gridjobeffi, LumiFactor,"gsfgsfetCategory",TString("E_{t}"),TString(""),legendname);
  // //   DrawHistosCategory(histogsfgsf_pt,weight, gridjobeffi, LumiFactor,"gsfgsfptCategory",TString("p_{t}"),TString(""),legendname);
  // //   DrawHistosCategory(histogsfgsf_eta,weight, gridjobeffi, LumiFactor,"gsfgsfetaCategory",TString("#eta"),TString(""),legendname);
  // //   DrawHistosCategory(histogsfgsf_phi,weight, gridjobeffi, LumiFactor,"gsfgsfphiCategory",TString("#phi"),TString(""),legendname);

  //   DrawHistosCategory(tagprobe_mass, weight, gridjobeffi, LumiFactor,"TPMassCategory",TString("Tag and Probe Mass (GeV/c^{2})"),TString(""),legendname);
  //   DrawHistosCategory(tagprobe_mass_80_100, weight, gridjobeffi, LumiFactor,"TPMass80100Category",TString("Tag and Probe Mass 80 100 (GeV/c^{2})"),TString(""),legendname);

  //   DrawHistosCategory(tagprobe_num, weight, gridjobeffi, LumiFactor,"TPNumCategory",TString("Number of TP pairs per event"),TString(""),legendname);
  //   DrawHistosCategory(tagprobe_num_80_100, weight, gridjobeffi, LumiFactor,"TPNum80100Category",TString("Number of TP pairs per event [80-100]"),TString(""),legendname);

  // //   DrawHistosCategory(CaloMet, weight, gridjobeffi, LumiFactor,"CaloMetCategory",TString("CaloMet"),TString(""),legendname);
  // //   DrawHistosCategory(Met, weight, gridjobeffi, LumiFactor,"MetCategory",TString("Met"),TString(""),legendname);
  // //   DrawHistosCategory(tagprobedeltaphi, weight, gridjobeffi, LumiFactor,"tagprobedeltaphiCategory",TString("tagprobedeltaphi"),TString(""),legendname);
  // //   DrawHistosCategory(NJetsAK, weight, gridjobeffi, LumiFactor,"NJetsAKCategory",TString("NJetsAK"),TString(""),legendname);
  // //   DrawHistosCategory(NJetsIC, weight, gridjobeffi, LumiFactor,"NJetsICCategory",TString("NJetsIC"),TString(""),legendname);

  //   DrawHistosCategory(histoprobept, weight, gridjobeffi, LumiFactor,"ProbePtCategory",TString("ProbePt"),TString(""),legendname);
  //   DrawHistosCategory(histoprobebarrelpt, weight, gridjobeffi, LumiFactor,"ProbebarrelptCategory",TString("Probebarrelpt"),TString(""),legendname);
  //   DrawHistosCategory(histoprobeendcappt, weight, gridjobeffi, LumiFactor,"ProbeendcapptCategory",TString("Probeendcappt"),TString(""),legendname);
  //   DrawHistosCategory(histoprobeeta, weight, gridjobeffi, LumiFactor,"ProbeEtaCategory",TString("ProbeEta"),TString(""),legendname);
  //   DrawHistosCategory(histoprobesceta, weight, gridjobeffi, LumiFactor,"ProbeScEtaCategory",TString("Probe Sc Eta"),TString(""),legendname);
  //   DrawHistosCategory(histoprobephi, weight, gridjobeffi, LumiFactor,"ProbePhiCategory",TString("ProbePhi"),TString(""),legendname);
  //   DrawHistosCategory(histoprobescphi, weight, gridjobeffi, LumiFactor,"ProbeScPhiCategory",TString("Probe Sc Phi"),TString(""),legendname);
  //   DrawHistosCategory(histoprobeenergy, weight, gridjobeffi, LumiFactor,"ProbeEnergyCategory",TString("ProbeEnergy"),TString(""),legendname);
  //   DrawHistosCategory(histoprobescenergy, weight, gridjobeffi, LumiFactor,"ProbeScEnergyCategory",TString("Probe Sc Energy"),TString(""),legendname);
  //   DrawHistosCategory(histoprobeEoverP, weight, gridjobeffi, LumiFactor,"ProbeEoverPCategory",TString("ProbeEoverP"),TString(""),legendname);

  //   DrawHistosCategory(histoprobe_deltaeta,weight, gridjobeffi, LumiFactor, "probedeltaetaCategory",TString("#Delta#eta"),TString(""),legendname);
  //   DrawHistosCategory(histoprobe_deltaphi,weight, gridjobeffi, LumiFactor,"probedeltaphiCategory",TString("#Delta#phi"),TString(""),legendname);
  //   DrawHistosCategory(histoprobe_hovere,weight, gridjobeffi, LumiFactor,"probehovereCategory",TString("H/E"),TString(""),legendname);
  //   DrawHistosCategory(histoprobe_trackiso,weight, gridjobeffi, LumiFactor,"probetrackisoCategory",TString("Track Iso."),TString(""),legendname);
  //   DrawHistosCategory(histoprobe_ecaliso,weight, gridjobeffi, LumiFactor,"probeecalisoCategory",TString("Ecal Iso."),TString(""),legendname);
  //   DrawHistosCategory(histoprobe_hcaliso1,weight, gridjobeffi, LumiFactor,"probehcaliso1Category",TString("Hcal Iso. 1"),TString(""),legendname);
  //   DrawHistosCategory(histoprobe_hcaliso2,weight, gridjobeffi, LumiFactor,"probehcaliso2Category",TString("Hcal Iso. 2"),TString(""),legendname);
  //   DrawHistosCategory(histoprobe_sigmaetaeta,weight, gridjobeffi, LumiFactor,"probesigmaetaetaCategory",TString("#sigma_{#eta#eta}"),TString(""),legendname);
  //   DrawHistosCategory(histoprobe_sigmaIetaIeta,weight, gridjobeffi, LumiFactor,"probesigmaIetaIetaCategory",TString("#sigma_{i#etai#eta}"),TString(""),legendname);

  //   DrawHistosCategory(histotagpt, weight, gridjobeffi, LumiFactor,"TagPtCategory",TString("TagPt"),TString(""),legendname);
  //   DrawHistosCategory(histotageta, weight, gridjobeffi, LumiFactor,"TagEtaCategory",TString("TagEta"),TString(""),legendname);
  //   DrawHistosCategory(histotagsceta, weight, gridjobeffi, LumiFactor,"TagScEtaCategory",TString("Tag Sc Eta"),TString(""),legendname);
  //   DrawHistosCategory(histotagphi, weight, gridjobeffi, LumiFactor,"TagPhiCategory",TString("TagPhi"),TString(""),legendname);
  //   DrawHistosCategory(histotagscphi, weight, gridjobeffi, LumiFactor,"TagScPhiCategory",TString("Tag Sc Phi"),TString(""),legendname);
  //   DrawHistosCategory(histotagenergy, weight, gridjobeffi, LumiFactor,"TagEnergyCategory",TString("TagEnergy"),TString(""),legendname);
  //   DrawHistosCategory(histotagscenergy, weight, gridjobeffi, LumiFactor,"TagScEnergyCategory",TString("Tag Sc Energy"),TString(""),legendname);
  //   DrawHistosCategory(histotagEoverP, weight, gridjobeffi, LumiFactor,"TagEoverPCategory",TString("TagEoverP"),TString(""),legendname);

  //   DrawHistosCategory(histotag_deltaeta,weight, gridjobeffi, LumiFactor, "tagdeltaetaCategory",TString("#Delta#eta"),TString(""),legendname);
  //   DrawHistosCategory(histotag_deltaphi,weight, gridjobeffi, LumiFactor,"tagdeltaphiCategory",TString("#Delta#phi"),TString(""),legendname);
  //   DrawHistosCategory(histotag_hovere,weight, gridjobeffi, LumiFactor,"taghovereCategory",TString("H/E"),TString(""),legendname);
  //   DrawHistosCategory(histotag_trackiso,weight, gridjobeffi, LumiFactor,"tagtrackisoCategory",TString("Track Iso."),TString(""),legendname);
  //   DrawHistosCategory(histotag_ecaliso,weight, gridjobeffi, LumiFactor,"tagecalisoCategory",TString("Ecal Iso."),TString(""),legendname);
  //   DrawHistosCategory(histotag_hcaliso1,weight, gridjobeffi, LumiFactor,"taghcaliso1Category",TString("Hcal Iso. 1"),TString(""),legendname);
  //   DrawHistosCategory(histotag_hcaliso2,weight, gridjobeffi, LumiFactor,"taghcaliso2Category",TString("Hcal Iso. 2"),TString(""),legendname);
  //   DrawHistosCategory(histotag_sigmaetaeta,weight, gridjobeffi, LumiFactor,"tagsigmaetaetaCategory",TString("#sigma_{#eta#eta}"),TString(""),legendname);
  //   DrawHistosCategory(histotag_sigmaIetaIeta,weight, gridjobeffi, LumiFactor,"tagsigmaIetaIetaCategory",TString("#sigma_{i#etai#eta}"),TString(""),legendname);
  //   DrawHistosCategory(histotag_gsfet,weight, gridjobeffi, LumiFactor,"taggsfetCategory",TString("E_{t}"),TString(""),legendname);


  //   DrawHistosCategory(histoprobeHEEPpt, weight, gridjobeffi, LumiFactor,"ProbeHEEPPtCategory",TString("ProbeHEEPPt"),TString(""),legendname);
  //   DrawHistosCategory(histoprobeHEEPbarrelpt, weight, gridjobeffi, LumiFactor,"ProbeHEEPBarrelPtCategory",TString("ProbeHEEPBarrelPt"),TString(""),legendname);
  //   DrawHistosCategory(histoprobeHEEPendcappt, weight, gridjobeffi, LumiFactor,"ProbeHEEPEndcapPtCategory",TString("ProbeHEEPEndcapPt"),TString(""),legendname);
  //   DrawHistosCategory(histoprobeHEEPeta, weight, gridjobeffi, LumiFactor,"ProbeHEEPEtaCategory",TString("ProbeHEEPEta"),TString(""),legendname);
  //   DrawHistosCategory(histoprobeHEEPsceta, weight, gridjobeffi, LumiFactor,"ProbeHEEPScEtaCategory",TString("ProbeHEEP Sc Eta"),TString(""),legendname);
  //   DrawHistosCategory(histoprobeHEEPphi, weight, gridjobeffi, LumiFactor,"ProbeHEEPPhiCategory",TString("ProbeHEEPPhi"),TString(""),legendname);
  //   DrawHistosCategory(histoprobeHEEPscphi, weight, gridjobeffi, LumiFactor,"ProbeHEEPScPhiCategory",TString("ProbeHEEP Sc Phi"),TString(""),legendname);
  //   DrawHistosCategory(histoprobeHEEPenergy, weight, gridjobeffi, LumiFactor,"ProbeHEEPEnergyCategory",TString("ProbeHEEPEnergy"),TString(""),legendname);
  //   DrawHistosCategory(histoprobeHEEPscenergy, weight, gridjobeffi, LumiFactor,"ProbeHEEPScEnergyCategory",TString("ProbeHEEP Sc Energy"),TString(""),legendname);


}//end of method






void DataMCTreeAnalysis::DrawHistos(vector<TH1F *> histos, vector <float> weight, vector <float> gridjobeffi, float LumiFactor, TString canvasname, TString XTitle, vector<TString> legendname)
{

  //setPlotStyle();
  //E-E SPECTRUM
  TCanvas *c0 = new TCanvas(canvasname.Data(),canvasname.Data(),100,100,1200,725);
  //gStyle->SetOptStat("emruo");

  int nbhistos = histos.size();
  float maxy = 0.;

  TLegend *legend = new TLegend(0.38,0.5,0.6,0.85);
  legend->SetTextSize(0.03);
  legend->SetLineColor(kWhite);
  legend->SetFillColor(kWhite);

  for (int p = 1; p<nbhistos; p++) {
    (histos[p])->Scale((weight[p]*gridjobeffi[p])*LumiFactor);
    (histos[p])->SetFillColor(p+1);
    if(p >= 9) (histos[p])->SetFillColor(p+2);

    //     TString histoname((histos[p])->GetName());
    //     int index1 = histoname.Index("indic")+5;
    //     int index2 = histoname.Index("V");
    //     TString legendname(histoname(index1,index2-index1));

    //legend->AddEntry((histos[p]),(histos[p])->GetName(),"f");
    legend->AddEntry((histos[p]), (legendname[p]).Data(),"f");
  }

  legend->AddEntry((histos[0]),"Data","l");

  cout<<"--------------------------------"<<endl;
  cout<<"----------HISTO INFO------------"<<endl;
  for(unsigned int f=0;f<histos.size();f++)
    {
      cout<<(histos[f])->GetName()<<"  integral "<<(histos[f])->Integral()<<"  entries "<<(histos[f])->GetEntries()<<endl;
    }
  cout<<"--------------------------------"<<endl;
  cout<<"----------HISTO INFO------------"<<endl;


  for (int p = 1; p<nbhistos; p++) { 
    for (int q = p+1; q<nbhistos; q++) {(histos[p])->Add(histos[q]);}
  }

  for (int p = 0; p<nbhistos; p++) {
    if( (histos[p])->GetMaximum() > maxy) maxy = (histos[p])->GetMaximum()+(histos[p])->GetBinError((histos[p])->GetMaximumBin());
  }

  //cout<<"maxy "<<maxy<<endl;

  for (int p = 0; p<nbhistos; p++) {
    (histos[p])->SetMaximum(maxy*3.);
    (histos[p])->SetMinimum(0.001);
    (histos[p])->GetXaxis()->SetTitle(XTitle.Data());
  }

  (histos[0])->SetLineWidth(2);
  (histos[0])->SetMarkerStyle(20);

  c0->cd();
  c0->SetLogy();
  for(unsigned int f=0;f<histos.size();f++)
    {
      if (f==1) (histos[f])->Draw("HIST");
      if (f >1)(histos[f])->Draw("HISTsames");
      //histos[f]->Draw("sames");
    }

  gPad->RedrawAxis(); 

  (histos[0])->SetFillColor(0);
  (histos[0])->Draw("sames");

  legend->Draw("sames");

  TPaveLabel *label0 = new TPaveLabel(0.59,0.64,0.85,0.74,"#sqrt{s} = 7 TeV, #int L dt = 19.7 pb^{-1}","brNDC");
  label0->SetFillColor(0);
  label0->SetFillStyle(0);
  label0->SetBorderSize(0);
  label0->SetTextSize(0.35);
  label0->Draw("sames");
  

  TPaveLabel *label1 = new TPaveLabel(0.592,0.75,0.732,0.82,"CMS preliminary","brNDC");
  label1->SetFillColor(0);
  label1->SetFillStyle(0);
  label1->SetBorderSize(0);
  label1->SetTextSize(0.45);
  label1->Draw("sames");


  TString jpgcanvasfilename(TString(canvasname.Data()).Append(".jpg"));
  TString gifcanvasfilename(TString(canvasname.Data()).Append(".gif"));
  TString epscanvasfilename(TString(canvasname.Data()).Append(".eps"));
  TString cxxcanvasfilename(TString(canvasname.Data()).Append(".C"));

  if(saveplots){
    c0->Print(jpgcanvasfilename.Data(),"JPG");
    c0->Print(gifcanvasfilename.Data(),"GIF");
    c0->Print(epscanvasfilename.Data(),"EPS");
    c0->Print(cxxcanvasfilename.Data(),"cxx");
  }


}
















void DataMCTreeAnalysis::DrawHistosCategory(vector<TH1F *> histos, vector <float> weight, vector <float> gridjobeffi, float LumiFactor, TString canvasname, TString XTitle, TString YTitle, vector<TString> legendname)
{

  //setPlotStyle();
  int log = 1;
  //E-E SPECTRUM
  TCanvas *c0 = new TCanvas(canvasname.Data(),canvasname.Data(),100,100,1200,725);
  //gStyle->SetOptStat("emruo");


  int nbhistos = histos.size();
  float maxy = 0.;


  TLegend *legend = new TLegend(0.38,0.65,0.6,0.85);
  legend->SetTextSize(0.03);
  legend->SetLineColor(kWhite);
  legend->SetFillColor(kWhite);


  for (int p = 1; p<nbhistos; p++) {
    (histos[p])->Scale((weight[p]*gridjobeffi[p])*LumiFactor);
  }


  int nbinshisto = (histos[1])->GetNbinsX();
  float xminhisto = (histos[1])->GetXaxis()->GetXmin();
  float xmaxhisto =  (histos[1])->GetXaxis()->GetXmax();

  cout<<"nbinshisto "<<nbinshisto<<endl;
  cout<<"xminhisto "<<xminhisto<<endl;
  cout<<"xmaxhisto "<<xmaxhisto<<endl;

  vector<TH1F *> newhistos;
  vector<TString> newlegendname;

  TString histoname((histos[0])->GetName());
  int index1 = histoname.Index("indic");
  TString histotitle(histoname(0,index1));


  newhistos.push_back(new TH1F(TString(histotitle).Append("DATA"),"",nbinshisto,xminhisto,xmaxhisto));
  //newhistos.push_back(new TH1F(TString(histotitle).Append("Zprime"),"",nbinshisto,xminhisto,xmaxhisto));
  newhistos.push_back(new TH1F(TString(histotitle).Append("Zee"),"",nbinshisto,xminhisto,xmaxhisto));
  newhistos.push_back(new TH1F(TString(histotitle).Append("TTbar"),"",nbinshisto,xminhisto,xmaxhisto));
  newhistos.push_back(new TH1F(TString(histotitle).Append("QCD"),"",nbinshisto,xminhisto,xmaxhisto));
  newhistos.push_back(new TH1F(TString(histotitle).Append("EWK"),"",nbinshisto,xminhisto,xmaxhisto));

  for(int k=0;k<newhistos.size();k++) {(histos[k])->Sumw2();}

  (newhistos[0])->Add((histos[0]));
  for (int p = 1; p<nbhistos; p++) { 
    cout<<"name histo "<<(histos[p])->GetName()<<" integral "<<(histos[p])->Integral()<<endl;
    //     if( TString((histos[p])->GetName()).Contains("Zprime") ) {(newhistos[1])->Add((histos[p]));cout<<"name "<<(newhistos[1])->GetName()<<endl;}
    //     if( TString((histos[p])->GetName()).Contains("Zee") ) {(newhistos[2])->Add((histos[p]));cout<<"name "<<(newhistos[2])->GetName()<<endl;}
    //     if( TString((histos[p])->GetName()).Contains("TTbar") ) {(newhistos[3])->Add((histos[p]));cout<<"name "<<(newhistos[3])->GetName()<<endl;}
    //     if( TString((histos[p])->GetName()).Contains("QCD") ) {(newhistos[4])->Add((histos[p]));cout<<"name "<<(newhistos[4])->GetName()<<endl;}
    //     if( TString((histos[p])->GetName()).Contains("Zmumu") || 
    // 	TString((histos[p])->GetName()).Contains("Ztautau") || 
    // 	TString((histos[p])->GetName()).Contains("WJet") || 
    // 	TString((histos[p])->GetName()).Contains("WW") || 
    // 	TString((histos[p])->GetName()).Contains("tW") ) {(newhistos[5])->Add((histos[p]));cout<<"name "<<(newhistos[5])->GetName()<<endl;}
    //if( TString((histos[p])->GetName()).Contains("Zprime") ) {(newhistos[1])->Add((histos[p]));cout<<"name "<<(newhistos[1])->GetName()<<endl;}
    if( TString((histos[p])->GetName()).Contains("Zee") ) {(newhistos[1])->Add((histos[p]));cout<<"name "<<(newhistos[1])->GetName()<<endl;}
    if( TString((histos[p])->GetName()).Contains("TTbar") ) {(newhistos[2])->Add((histos[p]));cout<<"name "<<(newhistos[2])->GetName()<<endl;}
    if( TString((histos[p])->GetName()).Contains("QCD") ) {(newhistos[3])->Add((histos[p]));cout<<"name "<<(newhistos[3])->GetName()<<endl;}
    if( TString((histos[p])->GetName()).Contains("Zmumu") || 
	TString((histos[p])->GetName()).Contains("Ztautau") || 
	TString((histos[p])->GetName()).Contains("WJet") || 
	TString((histos[p])->GetName()).Contains("WW") || 
	TString((histos[p])->GetName()).Contains("tW") ) {(newhistos[4])->Add((histos[p]));cout<<"name "<<(newhistos[4])->GetName()<<endl;}
  }
  

  newlegendname.push_back(TString("Data"));
  //newlegendname.push_back(TString("Z' ee M = 750 GeV/c^{2}"));
  newlegendname.push_back(TString("Zee PowHeg"));
  newlegendname.push_back(TString("TTbar+jets"));
  newlegendname.push_back(TString("QCD (Em + b/c)"));
  newlegendname.push_back(TString("EWK"));



  for (int p = 1; p<newhistos.size(); p++) {
    cout<<"debug "<<p<<endl;
    (newhistos[p])->SetFillColor(p+1);
    cout<<"debug "<<p<<endl;
    legend->AddEntry((newhistos[p]),(newlegendname[p]).Data(),"f");
    cout<<"debug "<<p<<endl;
  }
  legend->AddEntry((newhistos[0]), "Data","lp");


  cout<<"--------------------------------"<<endl;
  cout<<"----------HISTO INFO------------"<<endl;
  for(unsigned int f=0;f<newhistos.size();f++)
    {
      cout<<(newhistos[f])->GetName()<<"  integral "<<(newhistos[f])->Integral()<<"  entries "<<(newhistos[f])->GetEntries()<<endl;
    }
  cout<<"--------------------------------"<<endl;
  cout<<"----------HISTO INFO------------"<<endl;



  for (int p = 1; p<newhistos.size(); p++) { 
    for (int q = p+1; q<newhistos.size(); q++) {(newhistos[p])->Add(newhistos[q]);}
  }

  for (int p = 0; p<newhistos.size(); p++) {
    if( (newhistos[p])->GetMaximum() > maxy) maxy = (newhistos[p])->GetMaximum()+(newhistos[p])->GetBinError((newhistos[p])->GetMaximumBin());
  }

  cout<<"maxy "<<maxy<<endl;

  for (int p = 0; p<newhistos.size(); p++) {
    //(newhistos[p])->SetMaximum(maxy*15.);
    (newhistos[p])->SetMaximum(maxy*1.2);
    if (histoname.Contains("Zbosonpz") && (log==0)) {cout<<"HERE"<<endl;(newhistos[p])->SetMaximum(maxy*1.5);}
    (newhistos[p])->SetMinimum(0.01);
    (newhistos[p])->GetXaxis()->SetTitle(XTitle.Data());
    (newhistos[p])->GetYaxis()->SetTitle(YTitle.Data());
  }

  (newhistos[0])->SetLineWidth(2);
  (newhistos[0])->SetMarkerStyle(20);

  c0->cd();
  if (log==1) c0->SetLogy();
  for(unsigned int f=0;f<newhistos.size();f++)
    {
      if (f==1) (newhistos[f])->Draw("HIST");
      if (f>1)(newhistos[f])->Draw("HISTsames");
      //histos[f]->Draw("sames");
    }

  gPad->RedrawAxis(); 

  (newhistos[0])->SetFillColor(0);
  (newhistos[0])->Draw("sames");

  legend->Draw("sames");

  TPaveLabel *label0 = new TPaveLabel(0.59,0.64,0.85,0.74,"#sqrt{s} = 7 TeV, #int L dt = 19.7 pb^{-1}","brNDC");
  label0->SetFillColor(0);
  label0->SetFillStyle(0);
  label0->SetBorderSize(0);
  label0->SetTextSize(0.35);
  label0->Draw("sames");
  

  TPaveLabel *label1 = new TPaveLabel(0.592,0.75,0.732,0.82,"CMS preliminary","brNDC");
  label1->SetFillColor(0);
  label1->SetFillStyle(0);
  label1->SetBorderSize(0);
  label1->SetTextSize(0.45);
  label1->Draw("sames");


  TString jpgcanvasfilenameLOG(TString(canvasname.Data()).Append("LOG.jpg"));
  TString gifcanvasfilenameLOG(TString(canvasname.Data()).Append("LOG.gif"));
  TString epscanvasfilenameLOG(TString(canvasname.Data()).Append("LOG.eps"));
  TString cxxcanvasfilenameLOG(TString(canvasname.Data()).Append("LOG.C"));

  TString jpgcanvasfilename(TString(canvasname.Data()).Append(".jpg"));
  TString gifcanvasfilename(TString(canvasname.Data()).Append(".gif"));
  TString epscanvasfilename(TString(canvasname.Data()).Append(".eps"));
  TString cxxcanvasfilename(TString(canvasname.Data()).Append(".C"));

  if(saveplots){
    //     c0->cd();
    //     c0->SetLogy();

    //     c0->Print(jpgcanvasfilenameLOG.Data(),"JPG");
         c0->Print(gifcanvasfilenameLOG.Data(),"GIF");
    //     c0->Print(epscanvasfilenameLOG.Data(),"EPS");
    //     c0->Print(cxxcanvasfilenameLOG.Data(),"cxx");

    //     c0->cd();
    //     c0->SetLogy(0);

    //     //c0->Print(jpgcanvasfilename.Data(),"JPG");
    //     //c0->Print(gifcanvasfilename.Data(),"GIF");
    //     //c0->Print(epscanvasfilename.Data(),"EPS");
    //     c0->Print(cxxcanvasfilename.Data(),"cxx");
  }


}





































void DataMCTreeAnalysis::DrawTagAndProbeHistos(vector<TH1F *> histos, vector <float> weight, vector <float> gridjobeffi, float LumiFactor, TString canvasname, TString XTitle)
{

  //setPlotStyle();
  //E-E SPECTRUM
  TCanvas *c0 = new TCanvas(canvasname.Data(),canvasname.Data(),100,100,1200,725);
  gStyle->SetOptStat("emruo");

  int nbhistos = histos.size();
  float maxy = 0.;

  TLegend *legend = new TLegend(0.38,0.5,0.53,0.8);
  legend->SetTextSize(0.03);
  legend->SetLineColor(kWhite);
  legend->SetFillColor(kWhite);

  for (int p = 1; p<nbhistos; p++) {
    (histos[p])->SetLineColor(p+1);
    if(p >= 9) (histos[p])->SetLineColor(p+2);

    TString histoname((histos[p])->GetName());
    int index1 = histoname.Index("indic")+5;
    int index2 = histoname.Index("V");
    TString legendname(histoname(index1,index2-index1));

    //legend->AddEntry((histos[p]),(histos[p])->GetName(),"f");
    legend->AddEntry((histos[p]), legendname,"lf");
  }

  TString histoname((histos[0])->GetName());
  int index1 = histoname.Index("indic")+5;
  int index2 = histoname.Index("V");
  TString legendname(histoname(index1,index2-index1));
  //legend->AddEntry((histos[0]),(histos[0])->GetName(),"l");
  legend->AddEntry((histos[0]),legendname,"l");

  cout<<"--------------------------------"<<endl;
  cout<<"----------HISTO INFO------------"<<endl;
  for(unsigned int f=0;f<histos.size();f++)
    {
      cout<<(histos[f])->GetName()<<"  integral "<<(histos[f])->Integral()<<"  entries "<<(histos[f])->GetEntries()<<endl;
    }
  cout<<"--------------------------------"<<endl;
  cout<<"----------HISTO INFO------------"<<endl;

  for (int p = 0; p<nbhistos; p++) {
    if( (histos[p])->GetMaximum() > maxy) maxy = (histos[p])->GetMaximum()+(histos[p])->GetBinError((histos[p])->GetMaximumBin());
  }

  //cout<<"maxy "<<maxy<<endl;

  for (int p = 0; p<nbhistos; p++) {
    (histos[p])->SetMaximum(maxy*1.5);
    (histos[p])->SetMinimum(0.);
    (histos[p])->GetXaxis()->SetTitle(XTitle.Data());
  }

  c0->cd();
  //c0->SetLogy();
  for(unsigned int f=0;f<histos.size();f++)
    {
      //       if (f==1) (histos[f])->Draw("HIST");
      //       if (f >1)(histos[f])->Draw("HISTsames");
      if (f==1) (histos[f])->Draw("");
      if (f >1)(histos[f])->Draw("sames");
      //histos[f]->Draw("sames");
    }

  (histos[0])->SetFillColor(0);
  (histos[0])->Draw("sames");

  legend->Draw("sames");

  gPad->RedrawAxis(); 

  TString jpgcanvasfilename(TString(canvasname.Data()).Append(".jpg"));
  TString gifcanvasfilename(TString(canvasname.Data()).Append(".gif"));
  TString epscanvasfilename(TString(canvasname.Data()).Append(".eps"));
  TString cxxcanvasfilename(TString(canvasname.Data()).Append(".C"));

  if(saveplots){
    c0->Print(jpgcanvasfilename.Data(),"JPG");
    c0->Print(gifcanvasfilename.Data(),"GIF");
    c0->Print(epscanvasfilename.Data(),"EPS");
    c0->Print(cxxcanvasfilename.Data(),"cxx");
  }



}





























































void DataMCTreeAnalysis::DrawHistosBis(vector<TH1F *> histos, vector <float> weight, vector <float> gridjobeffi, float LumiFactor, TString canvasname, TString XTitle)
{

  //setPlotStyle();
  //E-E SPECTRUM
  TCanvas *c0 = new TCanvas(canvasname.Data(),canvasname.Data(),100,100,1200,725);
  gStyle->SetOptStat("emruoi");

  int nbhistos = histos.size();
  float maxy = 0.;

  TLegend *legend = new TLegend(0.38,0.5,0.53,0.8);
  legend->SetTextSize(0.03);
  legend->SetLineColor(kWhite);
  legend->SetFillColor(kWhite);

  for (int p = 0; p<nbhistos; p++) {
    (histos[p])->Scale((weight[p]*gridjobeffi[p])*LumiFactor);
    (histos[p])->SetFillColor(p+2);
    if(p >= 9) (histos[p])->SetFillColor(p+3);
    (histos[p])->SetLineColor(p+2);
    if(p >= 9) (histos[p])->SetLineColor(p+3);

    TString histoname((histos[p])->GetName());
    int index1 = histoname.Index("indic")+5;
    int index2 = histoname.Index("V");
    TString legendname(histoname(index1,index2-index1));

    legend->AddEntry((histos[p]), legendname,"f");
  }

  cout<<"--------------------------------"<<endl;
  cout<<"----------HISTO INFO------------"<<endl;
  for(unsigned int f=0;f<histos.size();f++)
    {
      cout<<(histos[f])->GetName()<<"  integral "<<(histos[f])->Integral()<<"  entries "<<(histos[f])->GetEntries()<<endl;
    }
  cout<<"--------------------------------"<<endl;
  cout<<"----------HISTO INFO------------"<<endl;

  for (int p = 0; p<nbhistos; p++) {
    if( (histos[p])->GetMaximum() > maxy) maxy = (histos[p])->GetMaximum()+(histos[p])->GetBinError((histos[p])->GetMaximumBin());
  }

  for (int p = 0; p<nbhistos; p++) {
    (histos[p])->SetMaximum(maxy*1.5);
    (histos[p])->SetMinimum(0.001);
    (histos[p])->GetXaxis()->SetTitle(XTitle.Data());
    (histos[p])->SetLineWidth(2);
  }

  c0->cd();
  c0->SetLogy();
  for(unsigned int f=0;f<histos.size();f++)
    {
      //       if (f==0) (histos[f])->Draw("HIST");
      //       if (f >0)(histos[f])->Draw("HISTsames");
      if (f==0) (histos[f])->Draw("");
      if (f >0)(histos[f])->Draw("sames");
      //histos[f]->Draw("sames");
    }

  gPad->RedrawAxis(); 

  legend->Draw("sames");


  TString jpgcanvasfilename(TString(canvasname.Data()).Append(".jpg"));
  TString gifcanvasfilename(TString(canvasname.Data()).Append(".gif"));
  TString epscanvasfilename(TString(canvasname.Data()).Append(".eps"));

  //   c0->Print(jpgcanvasfilename.Data(),"JPG");
  //   c0->Print(gifcanvasfilename.Data(),"GIF");
  //   c0->Print(epscanvasfilename.Data(),"EPS");



}








