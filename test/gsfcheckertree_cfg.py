import FWCore.ParameterSet.Config as cms

process = cms.Process("gsfcheckertree")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("RecoTracker.Configuration.RecoTracker_cff")


process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'START52_V8::All'
process.GlobalTag.globaltag = 'GR_R_52_V7::All'


readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.source = cms.Source("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

readFiles.extend( [
    '/store/mc/Summer12/DYToEE_M_20_TuneZ2star_8TeV_pythia6/AODSIM/PU_S7_START50_V15-v1/0000/FCCBFAC4-847E-E111-A077-002618943829.root'
##    'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/lathomas/Photon/LaurentPhoton-Run2011BSkim2ElePt35/319d9d50ddc1c21c2a4623a85e06b6f6/output_77_2_49z.root'
       ##   '/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/380EDCD0-CFFA-E011-8B63-002618943834.root',
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

process.hltPhysicsDeclared = cms.EDFilter('HLTPhysicsDeclared',
                                  invert = cms.bool(False),
                                  L1GtReadoutRecordTag = cms.InputTag('gtDigis')
                                  )

## # Primary vertex filter and no scraping events
## # https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCollisionsDataAnalysis
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32 (4),
                                           maxAbsZ = cms.double (15),
                                           maxd0 = cms.double (2)
                                           )
process.primaryVertexPath = cms.Path(process.primaryVertexFilter)

process.noscraping = cms.EDFilter("FilterOutScraping",
                                applyfilter = cms.untracked.bool(True),
                                debugOn = cms.untracked.bool(False),
                                numtrack = cms.untracked.uint32(10),
                                thresh = cms.untracked.double(0.25)
                                )
## from RecoJets.Configuration.CaloTowersRec_cff import *

## from RecoJets.JetProducers.CaloJetParameters_cfi import *
## from RecoJets.JetProducers.AnomalousCellParameters_cfi import *
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
  MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
)

## The next three lines are for rho computation (energy density, highly correlated to PU), see here :
## https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRecipesFor2011#FastJet_based_pile_up_isolation

process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)


from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
process.muIsoSequence = setupPFMuonIso(process, 'muons')

process.load("UserCode.HEEPSkims.gsfcheckertree_cfi")
process.otherStuff = cms.Sequence( process.kt6PFJets )


process.p1 = cms.Path(process.otherStuff * process.hltPhysicsDeclared * process.noscraping * process.primaryVertexFilter * process.pfParticleSelectionSequence * process.eleIsoSequence * process.muIsoSequence * process.gsfcheckerjob) 
##process.outpath = cms.EndPath(process.out)
