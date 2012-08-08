import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( -1 ) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       #'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/PromptReco-v1/000/195/530/0CDCB0BB-A8B0-E111-9059-BCAEC518FF8E.root',
       #'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/PromptReco-v1/000/195/530/14CFC608-B4B0-E111-83C6-0025B32036D2.root'
       #'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/MET/RECO/PromptReco-v1/000/195/552/04EA349F-5DB1-E111-A002-E0CB4E553651.root' 
       #'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SingleMu/RECO/PromptReco-v1/000/194/151/020D9700-FE9F-E111-9535-002481E0E56C.root'

        #'dcache:/pnfs/cms/WAX/11/store/data/Run2011B/Photon/AOD/PromptReco-v1/000/179/558/5886DB5E-A5FF-E011-92AE-BCAEC5364CFB.root'
       'dcache:/pnfs/cms/WAX/11/store/mc/Summer12/GMSB_Lambda-140_CTau-2000_TuneZ2star_8TeV-pythia6/GEN-SIM-RECO/PU_S9_START52_V9-v2/0000/003B2D4F-1C94-E111-B785-002618943976.root',
       'dcache:/pnfs/cms/WAX/11/store/mc/Summer12/GMSB_Lambda-140_CTau-2000_TuneZ2star_8TeV-pythia6/GEN-SIM-RECO/PU_S9_START52_V9-v2/0000/049FA3DE-1C94-E111-A57E-0026189438A0.root',
       'dcache:/pnfs/cms/WAX/11/store/mc/Summer12/GMSB_Lambda-140_CTau-2000_TuneZ2star_8TeV-pythia6/GEN-SIM-RECO/PU_S9_START52_V9-v2/0000/069BFC4C-2094-E111-809A-00248C55CC97.root'

 
    ),
    # explicitly drop photons resident in AOD/RECO, to make sure only those locally re-made (uncleaned photons) are used
    inputCommands = cms.untracked.vstring('keep *'
                                          #,'drop  *_photonCore_*_RECO' # drop hfRecoEcalCandidate as remade in this process
                                          #, 'drop *_photons_*_RECO' # drop photons as remade in this process
                                          )

)

process.options   = cms.untracked.PSet(
                    wantSummary = cms.untracked.bool(True),  
                    SkipEvent = cms.untracked.vstring('ProductNotFound')
)   


process.ana = cms.EDAnalyzer('DPAnalysis',
    rootFileName     = cms.untracked.string('Photon_2012B_test.root'),
    #triggerName      = cms.vstring('HLT_IsoMu30_v', 'HLT_DisplacedPhoton65_CaloIdVL_IsoL_PFMET25'),
    triggerName      = cms.vstring('HLT_Photon50_CaloIdVL_IsoL', 'HLT_DisplacedPhoton65_CaloIdVL_IsoL_PFMET25'),
    L1GTSource       = cms.string('L1_SingleEG22'),
    L1Select         = cms.bool( True ),
    isData           = cms.bool(True),
    trigSource = cms.InputTag("TriggerResults","","HLT"),
    jetSource   = cms.InputTag("ak5PFJets"),
    metSource   = cms.InputTag("pfMet"),
    muonSource  = cms.InputTag("muons"),
    trackSource = cms.InputTag("generalTracks"),
    electronSource   = cms.InputTag("gsfElectrons"),
    photonSource     = cms.InputTag("myphotons"),
    pvSource         = cms.InputTag("offlinePrimaryVerticesWithBS"),
    beamSpotSource   = cms.InputTag("offlineBeamSpot"),
    EBRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    EERecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    tau                = cms.double( 1000 ), 
    genParticles = cms.InputTag("genParticles"),

    # Set up cuts for physics objects
    # vertex cuts                z   ndof   d0 
    vtxCuts       = cms.vdouble( 99,    0,  99 ),
    # photon cuts                pt   eta  sMajMax,  sMinMin, sMinMax, trkVeto  Num  leadingPt 
    photonCuts    = cms.vdouble( 45,  2.4,    99.,      -1,     99.,    0.,   1 ,     80   ),
    # photon isolation           trk,  ecalSumEt, ecalR, hcalSumEt, hcalR 
    photonIso     = cms.vdouble(  0.2,       4.5,   0.1,       4.0,   0.1 ),
    # jet cuts                   pt    eta    dR,  nJets
    jetCuts       = cms.vdouble( 30. , 2.4,  0.3,    0 ),
    metCuts       = cms.vdouble(  0. ),
    # electron cuts              pt  eta  EBIso  EEIso nLostHit  
    electronCuts  = cms.vdouble( 25, 2.4,  0.15,   0.1,      2 ),
    # muon cuts                  pt  eta  Iso  dR   
    muonCuts      = cms.vdouble( 25, 2.1, 0.2, 0.3 ),

)

###########  USE UNCLEANED SUPERCLUSTERS  ######################### 
# Global Tag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_noesprefer_cff")
process.GlobalTag.globaltag = 'GR_P_V35::All'
#process.GlobalTag.globaltag = 'START52_V9::All'

# to get clustering 
process.load('Configuration/StandardSequences/GeometryExtended_cff')

# Geometry
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi") # gfwork: need this?


process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder")

process.load("RecoEcal.EgammaClusterProducers.uncleanSCRecovery_cfi")
process.uncleanSCRecovered.cleanScCollection=cms.InputTag ("correctedHybridSuperClusters")

# myPhoton sequence
process.load("RecoEgamma.PhotonIdentification.photonId_cff")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
#process.ecalSeverityLevel.timeThresh = cms.double(2.0)

import RecoEgamma.EgammaPhotonProducers.photonCore_cfi
import RecoEgamma.EgammaPhotonProducers.photons_cfi

process.myphotonCores=RecoEgamma.EgammaPhotonProducers.photonCore_cfi.photonCore.clone()
process.myphotonCores.scHybridBarrelProducer=cms.InputTag ("uncleanSCRecovered:uncleanHybridSuperClusters")

from RecoEgamma.PhotonIdentification.isolationCalculator_cfi import*
newisolationSumsCalculator = isolationSumsCalculator.clone()
newisolationSumsCalculator.barrelEcalRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
newisolationSumsCalculator.endcapEcalRecHitCollection = cms.InputTag('reducedEcalRecHitsEE')

process.myphotons=RecoEgamma.EgammaPhotonProducers.photons_cfi.photons.clone()
process.myphotons.barrelEcalHits=cms.InputTag("reducedEcalRecHitsEB")
process.myphotons.endcapEcalHits=cms.InputTag("reducedEcalRecHitsEE")
process.myphotons.isolationSumsCalculatorSet=newisolationSumsCalculator
process.myphotons.photonCoreProducer=cms.InputTag("myphotonCores")

process.myPhotonSequence = cms.Sequence(process.myphotonCores+
                                        process.myphotons)
# photonID sequence
from RecoEgamma.PhotonIdentification.photonId_cfi import *
process.myPhotonIDSequence = cms.Sequence(PhotonIDProd)
process.PhotonIDProd.photonProducer=cms.string("myphotons")

process.uncleanPhotons = cms.Sequence(
               process.uncleanSCRecovered *
               process.myPhotonSequence *
               process.myPhotonIDSequence
               )

process.p = cms.Path(
                     process.uncleanPhotons*
                     process.ana
                    )


