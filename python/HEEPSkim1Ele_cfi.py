import FWCore.ParameterSet.Config as cms

HEEPSkim1Ele = cms.EDFilter("HEEPSkim1Ele",
                   PtCut = cms.untracked.double(15.0),                 
                   HoECut = cms.untracked.double(0.1),                 
                   )

HEEPSkim1ElePath = cms.Path(HEEPSkim1Ele)
