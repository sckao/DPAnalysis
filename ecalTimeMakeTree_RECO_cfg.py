import FWCore.ParameterSet.Config as cms

process = cms.Process("TIMECALIBANALYSIS")

# gfworks: to get clustering 
process.load('Configuration/StandardSequences/GeometryExtended_cff')

# Geometry
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi") # gfwork: need this?


# Global Tag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_noesprefer_cff")
process.GlobalTag.globaltag = 'GR_R_42_V19::All'


# Trigger
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerScalesConfig_cff")
process.load("L1TriggerConfig.L1ScalesProducers.L1MuTriggerPtScaleConfig_cff")
process.load("L1TriggerConfig.L1GtConfigProducers.L1GtBoardMapsConfig_cff")
process.load("L1TriggerConfig.L1GtConfigProducers.L1GtConfig_cff")
process.load("L1TriggerConfig.L1GtConfigProducers.Luminosity.startup.L1Menu_startup2_v2_Unprescaled_cff")
import FWCore.Modules.printContent_cfi
process.dumpEv = FWCore.Modules.printContent_cfi.printContent.clone()

import EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi
process.gtDigis = EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi.l1GtUnpack.clone()


# this is the ntuple producer
process.load("CalibCalorimetry.EcalTiming.ecalTimePhyTree_cfi")
process.ecalTimePhyTree.fileName = 'EcalTimePhyTree'
process.ecalTimePhyTree.barrelEcalRecHitCollection = cms.InputTag("reducedEcalRecHitsEB","")
process.ecalTimePhyTree.endcapEcalRecHitCollection = cms.InputTag("reducedEcalRecHitsEE","")
process.ecalTimePhyTree.barrelBasicClusterCollection = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters")
# if you want cleaned EB BC use this: process.ecalTimePhyTree.barrelBasicClusterCollection = cms.InputTag("hybridSuperClusters","uncleanOnlyHybridBarrelBasicClusters")
process.ecalTimePhyTree.endcapBasicClusterCollection = cms.InputTag("multi5x5BasicClusters","multi5x5EndcapBasicClusters")
process.ecalTimePhyTree.barrelSuperClusterCollection = cms.InputTag("correctedHybridSuperClusters","")
#  if you want cleaned EB SC use this: process.ecalTimePhyTree.barrelSuperClusterCollection = cms.InputTag("hybridSuperClusters","uncleanOnlyHybridSuperClusters")
process.ecalTimePhyTree.endcapSuperClusterCollection = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower","")
process.ecalTimePhyTree.muonCollection = cms.InputTag("muons")
process.ecalTimePhyTree.runNum = 0
# Set up cuts for physics objects
# jet cuts                                           pt    eta  nJets
process.ecalTimePhyTree.jetCuts       = cms.vdouble( 25. , 2.4, 3 )
process.ecalTimePhyTree.metCuts       = cms.vdouble( 20  )
# photon cuts                                        pt  eta  dR   nPhoton
process.ecalTimePhyTree.photonCuts    = cms.vdouble( 30, 2.4, 0.3, 1 )
process.ecalTimePhyTree.electronCuts  = cms.vdouble( 25, 2.4, 0.15, 0.3 )
process.ecalTimePhyTree.muonCuts      = cms.vdouble( 25, 2.1, 0.2, 0.3 )

# Set up cuts for physics objects
process.ecalTimePhyTree.patMuonSource = cms.InputTag("muons")
process.ecalTimePhyTree.patJetSource = cms.InputTag("ak5PFJets")
process.ecalTimePhyTree.patMETSource = cms.InputTag("pfMet")
process.ecalTimePhyTree.patElectronSource = cms.InputTag("gsfElectrons")
process.ecalTimePhyTree.patPhotonSource = cms.InputTag("photons")

process.dumpEvContent = cms.EDAnalyzer("EventContentAnalyzer")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.p = cms.Path(
    # process.dumpEvContent  *
    process.ecalTimePhyTree
    )

process.options   = cms.untracked.PSet(
                    wantSummary = cms.untracked.bool(True),
                    SkipEvent = cms.untracked.vstring('ProductNotFound')
)


process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    ),
    categories = cms.untracked.vstring('ecalTimePhyTree'),
    destinations = cms.untracked.vstring('cout')
)
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)




# GF: some legacy reco files to test; replace w/ collision data
# dbs search --query "find file where dataset=/ExpressPhysics/BeamCommissioning09-Express-v2/FEVT and run=124020" | grep store | awk '{printf "\"%s\",\n", $1}'
process.source = cms.Source(
    "PoolSource",
    skipEvents = cms.untracked.uint32(0),
    
    # a few files from:    /MinimumBias/Commissioning10-GR_R_35X_V7A_SD_EG-v2/RECO
    fileNames = (cms.untracked.vstring(
    #'/store/data/Commissioning10/MinimumBias/RAW-RECO/v9/000/135/494/A4C5C9FA-C462-DF11-BC35-003048D45F7A.root',
    #'/store/data/Run2010A/EG/RECO/v4/000/144/114/EEC21BFA-25B4-DF11-840A-001617DBD5AC.root'
    #'file:AOD.root'
    #'file:/data/franzoni/data/Run2011A-MinimumBias-RECO-PromptReco-v1-run160406-0C132C90-434F-E011-8FF6-003048D2BF1C.root'
    #'file:/data/franzoni/data/Run2011A_DoubleElectron_AOD_PromptReco-v4_000_166_946_CE9FBCFF-4B98-E011-A6C3-003048F11C58.root'
    #'dcache:/pnfs/cms/WAX/11/store/mc/Summer11/G_Pt-50to80_TuneZ2_7TeV_pythia6/AODSIM/PU_S3_START42_V11-v2/0000/FC76E112-417F-E011-BAD8-0022649F01AA.root'
    'file:testPF2PAT.root'
    )
                 )
    )
