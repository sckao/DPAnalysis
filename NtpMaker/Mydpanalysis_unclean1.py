import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

##Technique for Extending Files
#readFiles = cms.untracked.vstring()
#secFiles = cms.untracked.vstring() 
#source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
#readFiles.extend( [
#       '/store/data/Run2012B/SinglePhoton/RAW-RECO/20Nov2012-v2/00000/C039A863-8734-E211-90FD-002354EF3BE0.root' ] );
#       secFiles.extend( [
#                      ] )


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )
process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'dcache:/pnfs/cms/WAX/11/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/191/247/002CDB97-4588-E111-8BF4-003048D37694.root',
      # 'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SinglePhoton/RECO/PromptReco-v1/000/193/998/10950429-439D-E111-B454-BCAEC5329705.root'


#        'file:/local/cms/phedex/store/data/Run2012C/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v3/000/198/969/00000/AA013C6D-9BD1-E111-AD92-001E67398412.root',
#	'file:/local/cms/phedex/store/data/Run2012C/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v3/000/198/969/00000/AC283CE3-9CD1-E111-962F-001E6739722E.root',
#   'file:/local/cms/phedex/store/data/Run2012B/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v1/000/195/163/0000/D0A14BCB-AFAC-E111-8605-001E67396FD1.root',
 #      'file:/local/cms/phedex/store/data/Run2012B/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v1/000/195/163/0000/E0E59635-B0AC-E111-A426-0025B3E06438.root',
 #          'file:/local/cms/phedex/store/data/Run2012B/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v1/000/195/163/0000/F6BD244C-AFAC-E111-B7A4-9C8E991A143E.root',
#	        'file:/local/cms/phedex/store/data/Run2012B/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v1/000/195/165/0000/102D08C2-72AC-E111-BF88-003048D47794.root',
#		     'file:/local/cms/phedex/store/data/Run2012B/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v1/000/195/165/0000/10BE5BE7-72AC-E111-967E-00304866C4AA.root',
#		          'file:/local/cms/phedex/store/data/Run2012B/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v1/000/195/165/0000/20240AE1-72AC-E111-86D7-003048D47A44.root',
#			       'file:/local/cms/phedex/store/data/Run2012B/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v1/000/195/165/0000/36AABE72-72AC-E111-8318-003048D479F0.root',
#			            'file:/local/cms/phedex/store/data/Run2012B/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v1/000/195/165/0000/B06AB000-72AC-E111-BE7C-D8D385FF4A7C.root',
#				         'file:/local/cms/phedex/store/data/Run2012B/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v1/000/195/165/0000/F285BEB7-72AC-E111-9738-003048D4772A.root',
#'file:/local/cms/phedex/store/data/Run2012B/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v1/000/195/247/0000/2AE05816-F1AC-E111-A456-002590200AB8.root',
#'file:/local/cms/phedex/store/data/Run2012B/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v1/000/195/249/0000/FC133B1E-25AC-E111-814E-001E6739811F.root',
#'file:/local/cms/phedex/store/data/Run2012B/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v1/000/195/251/0000/060A9527-DFAC-E111-A456-001E67398412.root',
#'file:/local/cms/phedex/store/data/Run2012B/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v1/000/195/251/0000/08BA2387-DFAC-E111-99AF-0025902009B8.root'
#'file:/local/cms/phedex/store/data/Run2012C/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v3/000/202/178/00000/06C5A54E-1DF9-E111-8830-0025902008D8.root',
#'file:/local/cms/phedex/store/data/Run2012C/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v3/000/202/178/00000/0C876365-1CF9-E111-BFA0-003048D476AA.root',
#'file:/local/cms/phedex/store/data/Run2012C/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v3/000/202/178/00000/246EEDFF-1BF9-E111-8128-0025B3E05C2C.root'
 
 'file:/local/cms/user/norbert/data/SinglePhoton-Run2012D-EXODisplacedPhoton-19Dec2012-v1-RECO/96A9B5BE-1566-E211-9523-003048D47A80.root',
       #'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/MET/RECO/PromptReco-v1/000/195/552/04EA349F-5DB1-E111-A002-E0CB4E553651.root' 
       #'dcache:/pnfs/cms/WAX/11/store/data/Run2012B/SingleMu/RECO/PromptReco-v1/000/194/151/020D9700-FE9F-E111-9535-002481E0E56C.root'

        #'dcache:/pnfs/cms/WAX/11/store/data/Run2011B/Photon/AOD/PromptReco-v1/000/179/558/5886DB5E-A5FF-E011-92AE-BCAEC5364CFB.root'
        #'dcache:/pnfs/cms/WAX/resilient/sckao/MC2012/GMSB_L100_CTau1000_8TeV_RECO_10_1_qX9.root'
        #'dcache:/pnfs/cms/WAX/11/store/data/Run2011A/Photon/AOD/PromptReco-v6/000/173/389/B8BBAB75-A4CA-E011-84BA-BCAEC53296F8.root'
 
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
    #rootFileName     = cms.untracked.string('Photon_2012B_test.root'),
    rootFileName     = cms.untracked.string('Test.root'),
    #triggerName      = cms.vstring('HLT_IsoMu30_v', 'HLT_DisplacedPhoton65_CaloIdVL_IsoL_PFMET25'),
    triggerName      = cms.vstring('HLT_Photon50_CaloIdVL_IsoL', 'HLT_DisplacedPhoton65_CaloIdVL_IsoL_PFMET25'),
    L1GTSource       = cms.string('L1_SingleEG22'),
    L1Select         = cms.bool( True ),
    isData           = cms.bool(True),
    cscHaloData      = cms.InputTag("CSCHaloData"),
    staMuons = cms.InputTag("standAloneMuons"),
    CSCSegmentCollection = cms.InputTag("cscSegments"),
    DTSegmentCollection = cms.InputTag("dt4DCosmicSegments"),
    muonSource = cms.InputTag("muonsFromCosmics"),
    trigSource = cms.InputTag("TriggerResults","","HLT"),
    jetSource   = cms.InputTag("ak5PFJets"),
    patJetSource = cms.InputTag("selectedPatJetsPFlow"),
    metSource   = cms.InputTag("pfMet"),
    #muonSource  = cms.InputTag("muons"),
    trackSource = cms.InputTag("generalTracks"),
    electronSource   = cms.InputTag("gsfElectrons"),
    photonSource     = cms.InputTag("myphotons"),
    pvSource         = cms.InputTag("offlinePrimaryVerticesWithBS"),
    beamSpotSource   = cms.InputTag("offlineBeamSpot"),
    EBRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    EERecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    tau                = cms.double( 1000 ), 
    genParticles = cms.InputTag("genParticles"),
    hbhereco     =cms.InputTag("hbhereco"),
    # Set up cuts for physics objects
    # vertex cuts                z   ndof   d0 
    vtxCuts       = cms.vdouble( 99,    0,  99 ),
    # photon cuts                pt   eta  sMajMax,  sMinMin, sMinMax, trkVeto  Num  leadingPt 
    photonCuts    = cms.vdouble( 45,  2.5,    999.,      -1,     99.,    0.05,   1 ,     60   ),
    # photon isolation           trk,  ecalSumEt, ecalR, hcalSumEt, hcalR 
    photonIso     = cms.vdouble(  0.2,       4.5,   0.1,       4.0,   0.1 ),
    # jet cuts                   pt    eta    dR,  nJets
    jetCuts       = cms.vdouble( 35. , 2.5,  0.3,    0 ),
    metCuts       = cms.vdouble(  0. ),
    # electron cuts              pt  eta  EBIso  EEIso nLostHit  
    electronCuts  = cms.vdouble( 25, 2.5,  0.15,   0.1,      2 ),
    # muon cuts                  pt  eta  Iso  dR   
    muonCuts      = cms.vdouble( 25, 2.1, 0.2, 0.3 ),

)

###########  USE UNCLEANED SUPERCLUSTERS  ######################### 
# Global Tag
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_noesprefer_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'GR_P_V35::All'
process.GlobalTag.globaltag = 'GR_R_70_V2::All'
#process.GlobalTag.globaltag = 'START52_V9::All'
## Choose Auto GlobalTag
#from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')


# to get clustering 
#process.load("Configuration.StandardSequences.Geometry_cff")
process.load('Configuration/StandardSequences/GeometryExtended_cff')

# Geometry
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cff")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.EcalMapping.EcalMapping_cfi")
process.load("Geometry.EcalMapping.EcalMappingRecord_cfi")
process.load("Geometry.MuonNumbering.muonNumberingInitialization_cfi") # gfwork: need this?

 # Specify IdealMagneticField ESSource
process.load("MagneticField.Engine.uniformMagneticField_cfi")
process.load("Geometry.CSCGeometryBuilder.cscGeometry_cfi")

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


