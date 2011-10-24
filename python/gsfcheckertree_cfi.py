import FWCore.ParameterSet.Config as cms

gsfcheckerjob = cms.EDAnalyzer("GsfCheckerTree",
  src = cms.InputTag('genParticles'),
  logEvents = cms.uint32(10),                                        
  TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
#  TriggerResultsTag = cms.InputTag("TriggerResults","","RECO"),
  usegendata = cms.bool(False),
##  usegendata = cms.bool(True),

  # for skimming on 2ele, 1ele-1mu and 2mu
  electronEtCut = cms.untracked.double(35.),
  muonPtCut = cms.untracked.double(35.),
)


