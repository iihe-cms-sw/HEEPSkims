import FWCore.ParameterSet.Config as cms

process = cms.Process("gsfcheckertree")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("RecoTracker.Configuration.RecoTracker_cff")


process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START3X_V26::All'

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

process.source = cms.Source("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

readFiles.extend( [
##     'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/agay/Photon/Photon-Run2010B-PromptReco-v2-RECO-HEEPSkimTwoGSFEleEt20HoE10-NoCert/75e3d59fbd0f207e610c6308251433be/skim2eleRECOPAT_207_1_ARd.root',
##     'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/agay/Photon/Photon-Run2010B-PromptReco-v2-RECO-HEEPSkimTwoGSFEleEt20HoE10-NoCert/75e3d59fbd0f207e610c6308251433be/skim2eleRECOPAT_206_1_KTJ.root',
##     'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/agay/Photon/Photon-Run2010B-PromptReco-v2-RECO-HEEPSkimTwoGSFEleEt20HoE10-NoCert/75e3d59fbd0f207e610c6308251433be/skim2eleRECOPAT_205_1_Cwe.root',
##     'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/agay/Photon/Photon-Run2010B-PromptReco-v2-RECO-HEEPSkimTwoGSFEleEt20HoE10-NoCert/75e3d59fbd0f207e610c6308251433be/skim2eleRECOPAT_204_1_mBe.root',
##     'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/agay/Photon/Photon-Run2010B-PromptReco-v2-RECO-HEEPSkimTwoGSFEleEt20HoE10-NoCert/75e3d59fbd0f207e610c6308251433be/skim2eleRECOPAT_203_1_nPe.root',

## 'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/agay/Electron/Electron-Run2010B-PromptReco-v2-RECO-HEEPSkimTwoGSFEleEt20HoE10-NoCert/75e3d59fbd0f207e610c6308251433be/skim2eleRECOPAT_209_1_gia.root',
## 'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/agay/Electron/Electron-Run2010B-PromptReco-v2-RECO-HEEPSkimTwoGSFEleEt20HoE10-NoCert/75e3d59fbd0f207e610c6308251433be/skim2eleRECOPAT_208_1_7hA.root',
## 'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/agay/Electron/Electron-Run2010B-PromptReco-v2-RECO-HEEPSkimTwoGSFEleEt20HoE10-NoCert/75e3d59fbd0f207e610c6308251433be/skim2eleRECOPAT_207_1_CFX.root',
## 'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/agay/Electron/Electron-Run2010B-PromptReco-v2-RECO-HEEPSkimTwoGSFEleEt20HoE10-NoCert/75e3d59fbd0f207e610c6308251433be/skim2eleRECOPAT_206_1_q2u.root',
## 'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/user/agay/Electron/Electron-Run2010B-PromptReco-v2-RECO-HEEPSkimTwoGSFEleEt20HoE10-NoCert/75e3d59fbd0f207e610c6308251433be/skim2eleRECOPAT_205_1_Nfo.root',

#Photon2011 AOD

###'dcap://maite.iihe.ac.be/pnfs/iihe/cms/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/994/386ED4FB-6F55-E011-8517-003048D2C01A.root'

        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/994/386ED4FB-6F55-E011-8517-003048D2C01A.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/956/7C4706C9-7155-E011-AF5B-003048D37456.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/955/FC17319C-DC55-E011-BA49-003048F1C424.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/955/848E2E62-3E55-E011-B1D7-001617C3B76E.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/954/A4EBB80B-2455-E011-BA11-0019B9F709A4.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/943/02257F10-FE54-E011-A8E8-0016177CA7A0.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/942/188031CD-F054-E011-A2C7-003048F118C6.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/940/185B0B39-0C55-E011-A949-001617C3B5D8.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/939/54C71398-1555-E011-AF89-000423D98B6C.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/938/1E35CB0E-1455-E011-B359-0030487D1BCC.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/937/BA8D1BAA-1A55-E011-8698-0030487C8CB6.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/936/206466F0-ED54-E011-A935-001617E30CC2.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/935/50691794-0955-E011-B5B1-000423D987E0.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/916/D033E56D-A954-E011-8133-001D09F24934.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/915/7052855B-7B54-E011-8200-003048F024E0.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/914/FADBFB7B-9454-E011-A9EC-001617C3B706.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/913/62F5C5EA-9854-E011-BB00-001617C3B77C.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/911/7AD46E20-C354-E011-B77E-003048F1C832.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/911/52EC45A4-7354-E011-A09C-003048F024E0.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/907/BC213A56-8D54-E011-B72B-0030487C90D4.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/898/4E78CEA5-F853-E011-A940-001D09F24353.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/894/A4FF7DC6-5E54-E011-AFE2-000423D94908.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/890/888230D4-3254-E011-A55B-001D09F28D54.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/890/504A6CAF-6854-E011-A898-003048F024F6.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/888/B433199F-2E54-E011-860D-0016177CA778.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/888/7AB9DB70-6754-E011-A099-003048D37538.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/877/A616636B-E453-E011-8488-00304879BAB2.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/876/DAF405AD-D453-E011-B905-003048F117B4.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/875/98507C6F-1854-E011-ACD0-003048D2BB90.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/874/525BC22F-0C54-E011-B0FB-00304879EE3E.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/873/0ED09777-0454-E011-8247-0030487CD6B4.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/872/98EE27F3-EB53-E011-9B9C-0030487CD76A.root',
        '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/871/1EB150FD-1554-E011-A7FF-0030487C90EE.root'
])

## replace PoolSource.fileNames = {
##         '/store/data/Run2011A/Photon/AOD/PromptReco-v1/000/160/994/386ED4FB-6F55-E011-8517-003048D2C01A.root'

##         }


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
    input = cms.untracked.int32(20000)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('gsfcheckertree_Test.root')
)


process.load("UserCode.HEEPSkims.gsfcheckertree_cfi")
process.p1 = cms.Path(process.gsfcheckerjob)
##process.outpath = cms.EndPath(process.out)
