import FWCore.ParameterSet.Config as cms

process = cms.Process("gsfcheckertree")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("RecoTracker.Configuration.RecoTracker_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.INFO.limit = 100000
process.source = cms.Source("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

readFiles.extend( [
#    '/store/user/treis/McGenSim3_DY/McReco3_DY/ba742de6463696a09655542c9e24f59b/DyToEE_madgraph_pu_7_1_DOp.root'
    '/store/mc/Summer12_DR53X/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6/AODSIM/PU_S10_START53_V7A-v1/0000/682F90A9-5FF0-E111-91D7-485B39800BBE.root'
#    '/store/mc/Summer12_DR53X/DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/0002/EE14FCF2-15F6-E111-8261-0025B3E05DDA.root'
#    '/store/mc/Summer12_DR53X/TTWJets_8TeV-madgraph/AODSIM/PU_S10_START53_V7A-v1/0000/0C570F59-B7DA-E111-A245-003048D437BA.root'
##   'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/lathomas/Photon/LaurentPhoton-Run2011BSkim2ElePt35/319d9d50ddc1c21c2a4623a85e06b6f6/output_77_2_49z.root'
##   '/store/mc/Fall11/TTJets_TuneZ2_7TeV-madgraph-tauola/AODSIM/PU_S6_START42_V14B-v1/0000/380EDCD0-CFFA-E011-8B63-002618943834.root',
])

# PFMET Type 1 (JEC) correction
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")

# select global tag and pfJetMET correction automatically from datasetpath when using multicrab
noDataset = True
import sys
for arg in sys.argv:
    if arg.startswith("-CMSSW.datasetpath"):
        dataset = arg[19:]
        noDataset =  False
        break
if noDataset:
    dataset = 'none'

if dataset.endswith('/AOD'):
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual") #this for data

    if dataset.find("13Jul2012") != -1:
        process.GlobalTag.globaltag = 'FT_53_V6C_AN4::All'  # this one for july13 rereco
    elif dataset.find("06Aug2012") != -1:
        process.GlobalTag.globaltag = 'FT_53_V6C_AN4::All' # this one for aug06 rereco
    elif dataset.find("24Aug2012") != -1:
        process.GlobalTag.globaltag = 'FT_53_V10A_AN4::All' # this one for aug24 rereco
    elif dataset.find("11Dec2012") != -1:
        process.GlobalTag.globaltag = 'FT_P_V42C_AN4::All' # this one for 11dec rereco
    elif dataset.find("Run2012C-PromptReco-v1") != -1:
        process.GlobalTag.globaltag = 'GR_P_V40_AN3::All'  # this one for run2012C v1 with cmssw version < 533
    elif dataset.find("Run2012C-PromptReco-v2") != -1:
        process.GlobalTag.globaltag = 'GR_P_V42_AN4::All'  # this one for run2012C v2 with cmssw version >= 533
    elif dataset.find("Run2012D") != -1:
        process.GlobalTag.globaltag = 'GR_P_V42_AN4::All'  # this one for run2012D
    else:
        process.GlobalTag.globaltag = 'GR_P_V39_AN2::All'  # this one for run2012A and B
else:
    process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3")  #this for MC
    process.GlobalTag.globaltag = 'START53_V7G::All'

print "Global Tag is ", process.GlobalTag.globaltag

##process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.out = cms.OutputModule("PoolOutputModule",
    ##process.FEVTSIMEventContent,
    fileName = cms.untracked.string('gsfchecker_TEST.root')
)

process.options = cms.untracked.PSet(
    #fileMode = cms.untracked.string('NOMERGE')
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    #wantSummary = cms.untracked.bool(True)
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
                                           maxAbsZ = cms.double (24),
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
#process.gsfcheckerjob.electron1EtMin = 0.
#process.gsfcheckerjob.electron2EtMin = 0.
process.otherStuff = cms.Sequence( process.kt6PFJets )

process.load("RecoMET.METFilters.ecalLaserCorrFilter_cfi")
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.MessageLogger.suppressError = cms.untracked.vstring ('ecalLaserCorrFilter') 

process.p1 = cms.Path(process.otherStuff * process.hltPhysicsDeclared * process.eeBadScFilter * process.ecalLaserCorrFilter * process.noscraping * process.primaryVertexFilter * process.pfParticleSelectionSequence * process.eleIsoSequence * process.muIsoSequence * process.producePFMETCorrections * process.gsfcheckerjob)
