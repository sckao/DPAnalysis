import FWCore.ParameterSet.Config as cms

process = cms.Process("test")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",

    fileNames = cms.untracked.vstring(
#'file:/mnt/hadoop/store/user/sckao/SinglePhoton/EXO_DisplacedPhoton2ndSKIM/176d86b18a84f276e44dc86f253a35ed/skim_c_247_1_KQl.root'
#'dcache:/pnfs/cms/WAX/11/store/data/Run2012D/Cosmics25ns/RECO/PromptReco-v1/000/208/980/8CEAA5A0-0347-E211-B91D-5404A63886B4.root'
'dcache:/pnfs/cms/WAX/11/store/data/Run2012D/Cosmics/RECO/22Jan2013-v1/10000/FC86048C-2C8B-E211-BD18-003048D47768.root'
#'dcache:/pnfs/cms/WAX/11/store/data/Run2012D/SinglePhoton/AOD/PromptReco-v1/000/208/888/0CB9AF2E-8D45-E211-9DD9-001D09F29114.root'

#'file:/local/cms/phedex/store/data/Run2012C/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v3/000/200/190/00000/18BB0794-8CDF-E111-B9B0-0025B31E3D3C.root'
    ),

    # explicitly drop photons resident in AOD/RECO, to make sure only those locally re-made (uncleaned photons) are used
    inputCommands = cms.untracked.vstring('keep *'
                                          #,'drop  *_photonCore_*_RECO' # drop hfRecoEcalCandidate as remade in this process
                                          #, 'drop *_photons_*_RECO' # drop photons as remade in this process
                                          )

)


#import EXO.DPAnalysis.skim2012c as fileList
#process.source.fileNames = fileList.fileNames

process.options   = cms.untracked.PSet(
                    wantSummary = cms.untracked.bool(True),  
                    SkipEvent = cms.untracked.vstring('ProductNotFound')
)   

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.ana = cms.EDAnalyzer('CosmicAnalysis',
    rootFileName     = cms.untracked.string('cosmic_test.root'),
    triggerName      = cms.vstring('HLT_Photon50_CaloIdVL_IsoL','HLT_DisplacedPhoton65_CaloIdVL_IsoL_PFMET25'),
    #L1GTSource       = cms.string('L1_SingleEG22'),
    L1GTSource       = cms.string('L1_SingleMu12'),
    CSCSegmentCollection = cms.InputTag("cscSegments"),
    DTSegmentCollection = cms.InputTag("dt4DSegmentsT0Seg"),
    #DTSegmentCollection = cms.InputTag("dt4DSegments"),
    #DTSegmentCollection = cms.InputTag("dt4DCosmicSegments"),
    EBRecHitCollection = cms.InputTag("ecalRecHit", "EcalRecHitsEB"),
    EERecHitCollection = cms.InputTag("ecalRecHit", "EcalRecHitsEE"),
    bcSource           = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters"),
    scSource           = cms.InputTag("correctedHybridSuperClusters"),
    trigSource          = cms.InputTag("TriggerResults","","HLT"),
    metSource           = cms.InputTag("tcMet"),
    #photonSource     = cms.InputTag("myphotons"),
    photonSource     = cms.InputTag("photons"),
    jetSource        = cms.InputTag("ak5CaloJets"),
    muonSource       = cms.InputTag("muons"),
    staMuons         = cms.InputTag("standAloneMuons"),

    L1Select         = cms.bool( True ),
    # Set up cuts for physics objects
    # photon cuts                pt   eta  sMajMax,  sMinMin, sMinMax,   dR,  Num  leadingPt  
    photonCuts    = cms.vdouble( 10,  2.4,     99.,      -1.,     99.,   0.0,  1,    10  ),
    # jet cuts                   pt    eta    dR,  nJets
    jetCuts       = cms.vdouble( 10. , 2.4,  0.3,    0 ),
    metCuts       = cms.vdouble( 0. ),
    # muon cuts                  pt  eta  Iso  dR   
    muonCuts      = cms.vdouble( 15, 2.1, 0.2, 0.3 ),

)

###########  USE UNCLEANED SUPERCLUSTERS  ######################### 
# Global Tag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_noesprefer_cff")
#process.GlobalTag.globaltag = 'GR_R_53_V18::All'
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag( process.GlobalTag, 'GR_R_53_V18::All' )


# to get clustering 
process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration/StandardSequences/GeometryExtended_cff')

# Geometry
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi") # gfwork: need this?


#process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder")

#process.load("RecoEcal.EgammaClusterProducers.uncleanSCRecovery_cfi")
#process.uncleanSCRecovered.cleanScCollection=cms.InputTag ("correctedHybridSuperClusters")

process.out = cms.OutputModule("PoolOutputModule" ,
                fileName = cms.untracked.string( 'patTuple_data.root' ) ,
		outputCommands = cms.untracked.vstring(
			'keep *'
			#               'keep *_cscSegments_*_*'
			#               *patEventContentNoCleaning
			)
																                 )


# this function will modify the PAT sequences.
#from PhysicsTools.PatAlgos.tools.pfTools import *

#postfix = "PFlow"

#usePF2PAT( process
#		, runPF2PAT = True
#		, jetAlgo   = 'AK5'
#		, runOnMC   = False
#		, postfix   = postfix
		# for MC
		#, jetCorrections=('AK5PFchs', ['L1FastJet','L2Relative','L3Absolute'])
		# for data
#		, jetCorrections=('AK5PFchs', ['L2L3Residual'])
#	 )



process.p = cms.Path(
#                     process.uncleanPhotons*
#		     getattr(process,"patPF2PATSequence"+postfix)*
#                     process.producePFMETCorrections *
                     process.ana
                    )

# top projections in PF2PAT:
#getattr(process,"pfNoPileUp"+postfix).enable = True
#getattr(process,"pfNoMuon"+postfix).enable = True
#getattr(process,"pfNoElectron"+postfix).enable = True
#getattr(process,"pfNoTau"+postfix).enable = False
#getattr(process,"pfNoJet"+postfix).enable = True

# verbose flags for the PF2PAT modules
#getattr(process,"pfNoMuon"+postfix).verbose = False

# enable delta beta correction for muon selection in PF2PAT?
#getattr(process,"pfIsolatedMuons"+postfix).doDeltaBetaCorrection = False

process.out.outputCommands.extend( [ 'drop *_*_*_*' ] )

