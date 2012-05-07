import FWCore.ParameterSet.Config as cms

from SHarper.HEEPAnalyzer.HEEPEventParameters_cfi import *

gsfcheckerjob = cms.EDAnalyzer("GsfCheckerTree",
  heepEventPara,
  TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
  #TriggerResultsTag = cms.InputTag("TriggerResults","","RECO"),
  centerOfMassEnergy = cms.double(8000.),

  #Pt cut for stored infos. This is *NOT* the skimming part ! 
  bJetPtMin = cms.untracked.double(10.),
  ScPtMin =  cms.untracked.double(10.),
  GsfPtMin=   cms.untracked.double(10.),
  GsfTrackPtMin=   cms.untracked.double(5.),
                               
  # for skimming on 2ele, 1ele-1mu
  electron1EtMin = cms.untracked.double(35),
  electron1EtMax = cms.untracked.double(1.E99),
  electron2EtMin = cms.untracked.double(35),
  electron2EtMax = cms.untracked.double(1.E99),
  muonPtMin = cms.untracked.double(35),
  muonPtMax = cms.untracked.double(1.E99),                        

  # Effective areas
  EcalHcal1EffAreaEndcaps=cms.untracked.double(0.28),
  EcalHcal1EffAreaBarrel=cms.untracked.double(0.28),

                               
  IsoDepElectron = cms.VInputTag(cms.InputTag('elPFIsoDepositChargedPFIso'),
                   cms.InputTag('elPFIsoDepositGammaPFIso'),
                   cms.InputTag('elPFIsoDepositNeutralPFIso')),
  IsoValElectronPF = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03PFIdPFIso'),
                   cms.InputTag('elPFIsoValueGamma03PFIdPFIso'),
                   cms.InputTag('elPFIsoValueNeutral03PFIdPFIso')),
  IsoValElectronNoPF = cms.VInputTag(cms.InputTag('elPFIsoValueCharged03NoPFIdPFIso'),
                       cms.InputTag('elPFIsoValueGamma03NoPFIdPFIso'),
                       cms.InputTag('elPFIsoValueNeutral03NoPFIdPFIso'))                 
)


