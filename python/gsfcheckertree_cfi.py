import FWCore.ParameterSet.Config as cms

from SHarper.HEEPAnalyzer.HEEPEventParameters_cfi import *

gsfcheckerjob = cms.EDAnalyzer("GsfCheckerTree",
  heepEventPara,
  TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
  #TriggerResultsTag = cms.InputTag("TriggerResults","","RECO"),
  centerOfMassEnergy = cms.double(7000.),
  bJetPtMin = cms.untracked.double(10.),

  # for skimming on 2ele, 1ele-1mu and 2mu
  electronEtCut = cms.untracked.double(35),
  muonPtCut = cms.untracked.double(35),
               
)


