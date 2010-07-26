import FWCore.ParameterSet.Config as cms

HEEPSkim2Ele = cms.EDFilter("HEEPSkim2Ele",
                   PtCut = cms.untracked.double(15.0),                 
                   HoECut = cms.untracked.double(0.1),                 
                   )

HEEPSkim2ElePath = cms.Path(HEEPSkim2Ele)
