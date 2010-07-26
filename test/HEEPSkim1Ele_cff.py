import FWCore.ParameterSet.Config as cms

# Process, how many events, inout files, ...
process = cms.Process("USER")
process.maxEvents = cms.untracked.PSet(
      input = cms.untracked.int32(-1)
      #input = cms.untracked.int32(20)
)


process.source = cms.Source("PoolSource",
##       debugVerbosity = cms.untracked.uint32(0),
##       debugFlag = cms.untracked.bool(False),
      fileNames = cms.untracked.vstring(

       '/store/data/Commissioning10/MinimumBias/RECO/v9/000/135/802/FA0CE9B9-4B65-DF11-BC5F-003048D47A40.root'

)
)


process.load("UserCode.HEEPSkims.HEEPSkim1Ele_cfi")


process.myEventContent = cms.PSet(
      outputCommands = cms.untracked.vstring(
            'keep *'
      )
)

process.SkimOutput = cms.OutputModule("PoolOutputModule",
      #process.AODSIMEventContent,
      process.myEventContent,
      SelectEvents = cms.untracked.PSet(
            SelectEvents = cms.vstring('HEEPSkim1ElePath')
      ),
      fileName = cms.untracked.string('output.root')
)

process.end = cms.EndPath(process.SkimOutput)


