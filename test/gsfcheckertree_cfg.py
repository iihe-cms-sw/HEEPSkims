import FWCore.ParameterSet.Config as cms

process = cms.Process("gsfcheckertree")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("RecoTracker.Configuration.RecoTracker_cff")


process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
##process.GlobalTag.globaltag = 'START42_V13::All'
process.GlobalTag.globaltag = 'GR_R_42_V18::All'

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

process.source = cms.Source("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

readFiles.extend( [
##        '/store/data/Run2011A/DoubleMu/AOD/May10ReReco-v1/0000/06F51A88-807C-E011-AC86-001A92971ACE.root',
##        '/store/data/Run2011A/DoubleMu/AOD/May10ReReco-v1/0000/06CBB6AF-D28B-E011-B5EE-00248C0BE013.root',
##        '/store/data/Run2011A/DoubleMu/AOD/May10ReReco-v1/0000/0662A4B8-687C-E011-B2C4-003048678C3A.root',
##        '/store/data/Run2011A/DoubleMu/AOD/May10ReReco-v1/0000/063975F5-B17B-E011-AAB2-003048678BAA.root',
##        '/store/data/Run2011A/DoubleMu/AOD/May10ReReco-v1/0000/060CAB29-B17B-E011-BCE6-002618943877.root',
##        '/store/data/Run2011A/DoubleMu/AOD/May10ReReco-v1/0000/046174B6-C17B-E011-A76D-00248C55CC3C.root',

##        '/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/875/24DC5B25-06DC-E011-9F25-003048D2C020.root',
##        '/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/874/582F2395-1FDC-E011-964A-003048F0258C.root',
##        '/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/874/40616F63-F9DC-E011-89B5-001D09F23F2A.root',
##        '/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/873/80730A1C-85DD-E011-A7B2-BCAEC53296F9.root',
##        '/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/872/6E90F149-37DE-E011-886B-003048D2C0F2.root',
##        '/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/866/FCC1BADC-D7DB-E011-B711-001D09F244BB.root',
##        '/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/866/D6DDD0A9-C2DB-E011-A408-001D09F2512C.root',
##        '/store/data/Run2011B/DoubleMu/AOD/PromptReco-v1/000/175/866/501F8B60-C8DB-E011-A823-003048D2C0F4.root',

##        '/store/data/Run2011B/SingleMu/AOD/PromptReco-v1/000/175/860/1EECB4B6-6BDC-E011-AA4D-E0CB4E4408C4.root',
##        '/store/data/Run2011B/SingleMu/AOD/PromptReco-v1/000/175/858/F2547E16-6FDC-E011-A361-BCAEC518FF7A.root',
##        '/store/data/Run2011B/SingleMu/AOD/PromptReco-v1/000/175/857/4676663E-44DC-E011-83E8-003048D2C020.root',
##        '/store/data/Run2011B/SingleMu/AOD/PromptReco-v1/000/175/837/C825E183-A8DB-E011-A133-BCAEC53296FF.root',
##        '/store/data/Run2011B/SingleMu/AOD/PromptReco-v1/000/175/837/C4408E75-CADB-E011-BF0E-001D09F24303.root',
##        '/store/data/Run2011B/SingleMu/AOD/PromptReco-v1/000/175/837/A2003FE6-B0DB-E011-AEFD-BCAEC518FF8E.root',

##        '/store/mc/Summer11/DYToEE_M-200_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S3_START42_V11-v2/0000/F41C0BF0-C988-E011-9D24-0017A477102C.root',
##        '/store/mc/Summer11/DYToEE_M-200_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S3_START42_V11-v2/0000/BE0BA381-B088-E011-919A-001E0B5FE542.root',
##        '/store/mc/Summer11/DYToEE_M-200_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S3_START42_V11-v2/0000/66FE9A12-AB88-E011-906E-D8D3855BBDC4.root',
##        '/store/mc/Summer11/DYToEE_M-200_TuneZ2_7TeV-pythia6-tauola/AODSIM/PU_S3_START42_V11-v2/0000/12228F90-9B88-E011-B2F5-0017A477102C.root',

##        '/store/mc/Summer11/DYToEE_M-200_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S4_START42_V11-v2/0000/EE1827DC-FC9D-E011-A5C4-0025B3E05CDE.root',
##        '/store/mc/Summer11/DYToEE_M-200_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S4_START42_V11-v2/0000/DE7B3F74-719E-E011-9C77-002590200934.root',
##        '/store/mc/Summer11/DYToEE_M-200_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S4_START42_V11-v2/0000/C81695E4-B29E-E011-8277-002481E14E2C.root',
##        '/store/mc/Summer11/DYToEE_M-200_TuneZ2_7TeV-powheg-pythia/GEN-SIM-RECO/PU_S4_START42_V11-v2/0000/B870F4BB-009E-E011-BE41-00E08178C0CD.root',

        '/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/380EDCD0-CFFA-E011-8B63-002618943834.root',
])

##process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.out = cms.OutputModule("PoolOutputModule",
    ##process.FEVTSIMEventContent,
    fileName = cms.untracked.string('gsfchecker_TEST.root')
)

process.options = cms.untracked.PSet(
    #fileMode = cms.untracked.string('NOMERGE')
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.maxEvents = cms.untracked.PSet(
    #input = cms.untracked.int32(1000)
    input = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('gsfcheckertree_Test.root')
)

## # Primary vertex filter
## # https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCollisionsDataAnalysis
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32 (4),
                                           maxAbsZ = cms.double (15),
                                           maxd0 = cms.double (2)
                                           )
process.primaryVertexPath = cms.Path(process.primaryVertexFilter)


## The next three lines are for rho computation (energy density, highly correlated to PU), see here :
## https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRecipesFor2011#FastJet_based_pile_up_isolation

process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)


process.load("UserCode.HEEPSkims.gsfcheckertree_cfi")
process.otherStuff = cms.Sequence( process.kt6PFJets ) 
process.p1 = cms.Path(process.otherStuff*process.primaryVertexFilter*process.gsfcheckerjob) 
##process.outpath = cms.EndPath(process.out)
