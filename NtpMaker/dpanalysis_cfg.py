import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'dcache:/pnfs/cms/WAX/11/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/191/247/002CDB97-4588-E111-8BF4-003048D37694.root',
        'dcache:/pnfs/cms/WAX/11/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/191/247/069259D0-2C88-E111-9729-5404A63886C3.root'
        #'dcache:/pnfs/cms/WAX/11/store/data/Run2012A/Photon/AOD/PromptReco-v1/000/191/247/'
    )
)

process.demo = cms.EDAnalyzer('DPAnalysis',
    rootFileName     = cms.untracked.string('Photon_2012A_test.root'),
    triggerName      = cms.untracked.string('HLT_Photon75_CaloIdVL_IsoL_v16'),
    isData           = cms.untracked.bool(True),
    trigSource = cms.InputTag("TriggerResults","","HLT"),
    jetSource   = cms.InputTag("ak5PFJets"),
    metSource   = cms.InputTag("pfMet"),
    muonSource  = cms.InputTag("muons"),
    electronSource   = cms.InputTag("gsfElectrons"),
    photonSource     = cms.InputTag("photons"),
    pvSource         = cms.InputTag("offlinePrimaryVertices"),
    beamSpotSource   = cms.InputTag("offlineBeamSpot"),
    EBRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    EERecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    tau                = cms.double( 1000 ), 
    genParticles = cms.InputTag("genParticles"),

    # Set up cuts for physics objects
    # vertex cuts                z   ndof   d0 
    vtxCuts       = cms.vdouble( 99,    0,  99 ),
    # jet cuts                   pt    eta    dR,  nJets
    jetCuts       = cms.vdouble( 30. , 2.4,  0.3,    0 ),
    # photon cuts                pt   eta  sMajMax,  sMinMin, sMinMax,  Num  
    photonCuts    = cms.vdouble( 50,  2.4,    999.,      0.0,     99.,    1  ),
    # photon isolation           trk,  ecalSumEt, ecalR, hcalSumEt, hcalR 
    photonIso     = cms.vdouble(  0.2,       4.5,   0.1,       4.0,   0.1 ),
    # electron cuts              pt  eta  EBIso  EEIso nLostHit  
    electronCuts  = cms.vdouble( 25, 2.4,  0.15,   0.1,      2 ),
    metCuts       = cms.vdouble(  0. ),
    muonCuts      = cms.vdouble( 25, 2.1, 0.2, 0.3 ),

)


process.p = cms.Path(process.demo)
