// -*- C++ -*-
// Package:    DPAnalysis
// Class:      DPAnalysis
// 
/**\class DPAnalysis DPAnalysis.cc EXO/DPAnalysis/src/DPAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Sat Oct  8 06:50:16 CDT 2011
// $Id$
//
//


// system include files
#include "DPAnalysis.h"
#include "Ntuple.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"

// For DT Segment
#include "Geometry/DTGeometry/interface/DTChamber.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"

// For PFIsolation
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"

// global tracking geometry
//#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
//#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

using namespace cms ;
using namespace edm ;
using namespace std ;

//static bool HtDecreasing( VtxInfo s1, VtxInfo s2) { return ( s1.ht > s2.ht ); }
//static bool ZDecreasing( VtxInfo s1, VtxInfo s2) { return ( s1.z > s2.z ); }
//static bool Z0Decreasing( TrkInfo s1, TrkInfo s2) { return ( s1.dz > s2.dz ); }

// constants, enums and typedefs
// static data member definitions

// constructors and destructor
DPAnalysis::DPAnalysis(const edm::ParameterSet& iConfig){

   //now do what ever initialization is needed
   rootFileName         = iConfig.getUntrackedParameter<string> ("rootFileName");
   trigSource           = iConfig.getParameter<edm::InputTag> ("trigSource");
   l1GTSource           = iConfig.getParameter<string> ("L1GTSource");
   pvSource             = iConfig.getParameter<edm::InputTag> ("pvSource");
   beamSpotSource       = iConfig.getParameter<edm::InputTag> ("beamSpotSource");
   muonSource           = iConfig.getParameter<edm::InputTag> ("muonSource");
   electronSource       = iConfig.getParameter<edm::InputTag> ("electronSource");
   photonSource         = iConfig.getParameter<edm::InputTag> ("photonSource");
   metSource            = iConfig.getParameter<edm::InputTag> ("metSource");
   type1metSource       = iConfig.getParameter<edm::InputTag> ("type1metSource");
   jetSource            = iConfig.getParameter<edm::InputTag> ("jetSource");
   patJetSource         = iConfig.getParameter<edm::InputTag> ("patJetSource");
   trackSource          = iConfig.getParameter<edm::InputTag> ("trackSource");

   EBRecHitCollection   = iConfig.getParameter<edm::InputTag> ("EBRecHitCollection") ;
   EERecHitCollection   = iConfig.getParameter<edm::InputTag> ("EERecHitCollection") ;
   DTSegmentTag         = iConfig.getParameter<edm::InputTag> ("DTSegmentCollection") ;
   CSCSegmentTag        = iConfig.getParameter<edm::InputTag> ("CSCSegmentCollection") ;
   cscHaloTag           = iConfig.getParameter<edm::InputTag> ("cscHaloData");
   staMuons             = iConfig.getParameter<edm::InputTag> ("staMuons");
   EBSuperClusterCollection = iConfig.getParameter<edm::InputTag> ("EBSuperClusterCollection") ;
   EESuperClusterCollection = iConfig.getParameter<edm::InputTag> ("EESuperClusterCollection") ;

   //pileupSource         = iConfig.getParameter<edm::InputTag>("addPileupInfo");
   vtxCuts              = iConfig.getParameter<std::vector<double> >("vtxCuts");
   jetCuts              = iConfig.getParameter<std::vector<double> >("jetCuts");
   metCuts              = iConfig.getParameter<std::vector<double> >("metCuts");
   photonCuts           = iConfig.getParameter<std::vector<double> >("photonCuts");
   photonIso            = iConfig.getParameter<std::vector<double> >("photonIso");
   electronCuts         = iConfig.getParameter<std::vector<double> >("electronCuts");
   muonCuts             = iConfig.getParameter<std::vector<double> >("muonCuts");  
   //triggerPatent        = iConfig.getUntrackedParameter<string> ("triggerName");
   triggerPatent        = iConfig.getParameter< std::vector<string> >("triggerName");
   L1Select             = iConfig.getParameter<bool> ("L1Select");
   isData               = iConfig.getParameter<bool> ("isData");
   useRECO              = iConfig.getParameter<bool> ("useRECO");
   tau                  = iConfig.getParameter<double> ("tau");

   const InputTag TrigEvtTag("hltTriggerSummaryAOD","","HLT");
   trigEvent            = iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag", TrigEvtTag);

   gen = new GenStudy( iConfig );

   theFile  = new TFile( rootFileName.c_str(), "RECREATE") ;
   theFile->cd () ;
   theTree  = new TTree ( "DPAnalysis","DPAnalysis" ) ;
   setBranches( theTree, leaves ) ;

   CutFlowTree = new TTree( "CutFlow", "CutFlow") ;
   CutFlowTree->Branch("counter",   counter,    "counter[12]/I");

   h_z0 = new TH1D("h_z0", " z0 for tracks", 121, -181.5, 181.5 ) ;

   targetTrig = 0 ;
   firedTrig.clear() ;
   for ( size_t i=0; i< triggerPatent.size(); i++ ) firedTrig.push_back(-1) ;  

   // reset the counter
   for ( int i=0; i< 12 ; i++) counter[i] = 0 ;

   runID_ = 0 ;
   debugT = false ;
   rhoIso = 0 ;
   beamspot = 0 ;

   // PF isolation for photon
   isolator.initializePhotonIsolation(kTRUE);
   isolator.setConeSize(0.3);

}


DPAnalysis::~DPAnalysis()
{
   // do anything here that needs to be done at desctruction time

   delete gen ;
   cout<<"All:"<< counter[0]<<" Trigger:"<<counter[1]<<" Vertex:"<< counter[2] ;
   cout<<" Fiducial: "<< counter[3]<<" Conversion: "<< counter[4]<<" sMaj_sMin: "<< counter[5]<<" dR(g,trk):"<< counter[6] ;
   cout<<" LeadingPt:"<<counter[7] ;
   cout<<" beamHalo:"<< counter[8] <<" Jet:"<< counter[9] <<" MET:"<<counter[10] <<" Pre-Selection:"<<counter[11] <<endl ;
   CutFlowTree->Fill() ;

   theFile->cd () ;
   theTree->Write() ; 
   CutFlowTree->Write() ; 
   h_z0->Write() ;

   theFile->Close() ;

}

//
// member functions
//

// ------------ method called for each event  ------------
void DPAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

   // get calibration service
   // IC's
   iSetup.get<EcalIntercalibConstantsRcd>().get(ical);
   // ADCtoGeV
   iSetup.get<EcalADCToGeVConstantRcd>().get(agc);
   // transp corrections
   iSetup.get<EcalLaserDbRecord>().get(laser);
   // Geometry
   //iSetup.get<CaloGeometryRecord> ().get (pGeometry) ;
   //theGeometry = pGeometry.product() ;

   // event time
   eventTime = iEvent.time() ;
   // Initialize ntuple branches
   initializeBranches( theTree, leaves );

   leaves.bx          = iEvent.bunchCrossing();
   leaves.lumiSection = iEvent.id().luminosityBlock();
   leaves.orbit       = iEvent.orbitNumber();
   leaves.runId       = iEvent.id().run() ;
   leaves.eventId     = iEvent.id().event() ;

   /* 
   Handle<std::vector< PileupSummaryInfo > >  PileUpInfo;
   iEvent.getByLabel(pileupSource, PupInfo);
   for( std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
       std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
   }
   */

   // For conversion veto
   iEvent.getByLabel( beamSpotSource, bsHandle);
   beamspot = bsHandle.product();

   iEvent.getByLabel("allConversions", hConversions);
   iEvent.getByLabel( electronSource, electrons);

   // For PFIso
   iEvent.getByLabel( "kt6PFJets", "rho", rho_ );
   rhoIso = *(rho_.product());

   // Get the JES Uncertainty
   edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
   iSetup.get<JetCorrectionsRecord>().get("AK5PFchs",JetCorParColl); 
   JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
   jecUnc = new JetCorrectionUncertainty(JetCorPar);

   if (counter[0] == 0 )  PrintTriggers( iEvent ) ;

   int run_id    = iEvent.id().run()  ;

   // Global Tracking Geometry
   //ESHandle<GlobalTrackingGeometry> trackingGeometry;
   //iSetup.get<GlobalTrackingGeometryRecord>().get(trackingGeometry);

   counter[0]++ ;  // All events

   // L1 Trigger Selection
   passL1 = L1TriggerSelection( iEvent, iSetup ) ;

   // HLT trigger analysis
   Handle<edm::TriggerResults> triggers;
   iEvent.getByLabel( trigSource, triggers );
   const edm::TriggerNames& trgNameList = iEvent.triggerNames( *triggers ) ;

   TriggerTagging( triggers, trgNameList, run_id, firedTrig ) ;
   passHLT = TriggerSelection( triggers, firedTrig ) ;

   // Using L1 or HLT to select events ?!
   bool passTrigger = ( L1Select ) ? passL1 : passHLT  ;

   if ( passTrigger ) counter[1]++ ;   // Pass trigger cut

   // get the generator information
   if ( !isData && tau > -0.1 ) { 
      gen->GetGenEvent( iEvent, leaves, true );
      //gen->GetGen( iEvent, leaves );
   }
   //if ( !isData ) gen->PrintGenEvent( iEvent );

   bool pass = EventSelection( iEvent, iSetup ) ;
   if ( pass && passTrigger ) counter[11]++ ;   

   // fill the ntuple
   if ( !isData && pass ) theTree->Fill();
   if ( pass && isData && passTrigger ) theTree->Fill();
   delete jecUnc ;  
}

bool DPAnalysis::EventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

   Handle<reco::BeamHaloSummary>       beamHaloSummary ;
   Handle<edm::TriggerResults>         triggers;
   Handle<reco::VertexCollection>      recVtxs;
   Handle<reco::PhotonCollection>      photons; 
   //Handle<reco::GsfElectronCollection> electrons; 
   Handle<reco::MuonCollection>        muons; 
   Handle<reco::PFJetCollection>       jets; 
   //Handle<reco::GenJetCollection>      genjets; 
   Handle<std::vector<pat::Jet> >      patjets;
   Handle<reco::PFMETCollection>       met0; 
   Handle<reco::PFMETCollection>       met; 
   Handle<EcalRecHitCollection>        recHitsEB ;
   Handle<EcalRecHitCollection>        recHitsEE ;
   Handle<reco::TrackCollection>       tracks; 
   Handle<reco::PFCandidateCollection>           pfCand ;
   Handle<reco::SuperClusterCollection>  scEB ;
   Handle<reco::SuperClusterCollection>  scEE ;

   iEvent.getByLabel( trigSource,     triggers );
   iEvent.getByLabel( pvSource,       recVtxs  );
   iEvent.getByLabel( photonSource,   photons  );
   //iEvent.getByLabel( electronSource, electrons);
   iEvent.getByLabel( muonSource,     muons );
   iEvent.getByLabel( jetSource,      jets  );
   iEvent.getByLabel( patJetSource,   patjets);
   //iEvent.getByLabel( "ak5GenJets",   genjets  );
   iEvent.getByLabel( metSource,      met0  );
   iEvent.getByLabel( type1metSource, met  );
   iEvent.getByLabel( EBRecHitCollection,     recHitsEB );
   iEvent.getByLabel( EERecHitCollection,     recHitsEE );
   iEvent.getByLabel( EBSuperClusterCollection , scEB) ;
   iEvent.getByLabel( EESuperClusterCollection , scEE) ;
   iEvent.getByLabel("BeamHaloSummary", beamHaloSummary) ;
   iEvent.getByLabel( trackSource,    tracks  );
   iEvent.getByLabel( "particleFlow", pfCand ) ;
   
   bool passEvent = true ;

   // find trigger matched objects
   const reco::Photon rPho ;
   GetTrgMatchObject(  rPho , iEvent,  photonSource ) ;
   const reco::PFMET pfMet_ = (*met)[0] ;
   GetTrgMatchObject( pfMet_, iEvent,  metSource ) ;

   Track_Z0( tracks ) ;

   bool hasGoodVtx = VertexSelection( recVtxs );
   if ( !hasGoodVtx ) passEvent = false ;
   if ( passEvent )   counter[2]++ ;  // pass vertex cuts

   selectedPhotons.clear() ;
   PhotonSelection( photons, recHitsEB, recHitsEE, tracks, selectedPhotons ) ;

   // Check event flow with photon cuts 
   if ( passEvent ) {
      if ( gcounter[2] > 0 ) counter[3]++  ;  // pt & fiducial 
      if ( gcounter[3] > 0 ) counter[4]++  ;  // conversion 
      if ( gcounter[5] > 0 ) counter[5]++ ;   // sMaj & sMin
      if ( gcounter[6] > 0 ) counter[6]++ ;   // dR(phot, track)
   }

   if ( selectedPhotons.size() < (size_t)photonCuts[6] )  passEvent = false ;
   if ( passEvent )   counter[7]++ ;  // pass photon cuts

   // Stupid PFIso 
   if ( passEvent ) { 
      reco::VertexRef vtxRef(recVtxs, 0);
      PhotonPFIso( selectedPhotons, &(*pfCand), vtxRef, recVtxs ) ;
   }

   if( beamHaloSummary.isValid() ) {
     const reco::BeamHaloSummary TheSummary = (*beamHaloSummary.product() ); 
     if( !TheSummary.CSCTightHaloId() && passEvent ) { 
       counter[8]++ ;  
     } else {
       passEvent = false ;
     }
   } else {
       counter[8]++ ;
   }
   CSCHaloCleaning( iEvent, selectedPhotons ) ;

   // for cscsegments and halo muon/photon studies 
   if ( useRECO ) { 
      Handle<CSCSegmentCollection>       cscSegments ;
      iEvent.getByLabel( CSCSegmentTag,  cscSegments );
      BeamHaloMatch( cscSegments, selectedPhotons, iSetup ) ;

      Handle<DTRecSegment4DCollection>   dtSegments ;
      iEvent.getByLabel( DTSegmentTag,   dtSegments );
      CosmicRayMatch( dtSegments, selectedPhotons, iSetup ) ;
   }
   //CosmicRayMatch( muons, selectedPhotons ) ;
   
   //IsoPhotonSelection( selectedPhotons ) ;
   //if ( selectedPhotons.size() < photonCuts[5] )  passEvent = false ;

   selectedJets.clear() ;
   //JetSelection( jets, selectedPhotons, selectedJets );
   JetSelection( patjets, selectedPhotons, scEB, scEE, recHitsEB, recHitsEE, selectedJets );
   //bool isGammaJets = GammaJetVeto( selectedPhotons, selectedJets ) ;
   //if ( isGammaJets ) passEvent = false ;
   if ( selectedJets.size() < jetCuts[3] )   passEvent = false ;
   if ( passEvent )   counter[9]++ ;   
   
   //JERUncertainty( patjets ) ;
   selectedElectrons.clear() ;
   ElectronSelection( electrons, selectedElectrons ) ;

   selectedMuons.clear() ;
   MuonSelection( muons, selectedMuons );

   //HLTMET( jets, selectedMuons );   
   const reco::PFMET pfMet0 = (*met0)[0] ;
   leaves.met0   = pfMet0.et() ;
   leaves.met0Px = pfMet0.px() ;
   leaves.met0Py = pfMet0.py() ;

   const reco::PFMET pfMet = (*met)[0] ;
   leaves.met   = pfMet.et() ;
   leaves.metPx = pfMet.px() ;
   leaves.metPy = pfMet.py() ;

   if ( pfMet.pt() < metCuts[0]  ) passEvent = false;
   if ( passEvent )  counter[10]++ ;  

   return passEvent ;
}


bool DPAnalysis::L1TriggerSelection( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

    // Get L1 Trigger menu
    ESHandle<L1GtTriggerMenu> menuRcd;
    iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
    const L1GtTriggerMenu* menu = menuRcd.product();
    // Get L1 Trigger record  
    Handle< L1GlobalTriggerReadoutRecord > gtRecord;
    iEvent.getByLabel( edm::InputTag("gtDigis"), gtRecord);
    // Get dWord after masking disabled bits
    const DecisionWord dWord = gtRecord->decisionWord();
 
   bool l1_accepted = menu->gtAlgorithmResult( l1GTSource , dWord);
   //int passL1 =  ( l1SingleEG22 ) ? 1 : 0 ; 

   //cout<<" pass L1 EG22 ? "<<  passL1 <<endl ;
   if ( l1_accepted ) leaves.L1a = 1 ;

   return l1_accepted ;
}

void DPAnalysis::CSCHaloCleaning( const edm::Event& iEvent, vector<const reco::Photon*>& selectedPhotons ) { 

   Handle<reco::CSCHaloData> cschalo;
   iEvent.getByLabel( cscHaloTag, cschalo ); 

   if ( cschalo.isValid() ) {
      const reco::CSCHaloData cscData = *(cschalo.product()) ;
      int nOutTimeHits       = cscData.NumberOfOutTimeHits() ;
      int nMinusHTrks        = cscData.NHaloTracks( reco::HaloData::minus ) ; 
      int nPlusHTrks         = cscData.NHaloTracks( reco::HaloData::plus )  ; 
      //int nTrksSmallBeta     = cscData.NTracksSmallBeta()    ; 
      //int nHaloSegs          = cscData.NFlatHaloSegments() ;
       
      //cout<<" nOutT: "<< nOutTimeHits  <<" nTrkBeta:"<< nTrksSmallBeta <<" nHaloSeg: "<< nHaloSegs ;
      //cout<<" N Minus tracks : "<< nMinusHTrks <<" N Plus tracks : "<< nPlusHTrks << endl ;
      RefVector< reco::TrackCollection > trkRef  = cscData.GetTracks() ;
      //cout<<" NTrkRefs = "<<  trkRef.size()  << endl ;
      /*
      // the first track is closed to the ImpactPosition
      for( RefVector< reco::TrackCollection >::const_iterator it =  trkRef.begin();  it != trkRef.end() ; ++it ) {
         const vector<reco::Track>* trks = it->product() ;
         for ( size_t j=0; j< trks->size(); j++ ) {
             if ( ! (*trks)[j].innerOk() ) continue ;
    	     math::XYZPoint tXYZ = (*trks)[j].innerPosition() ;
    	     cout<<"    --> ( "<< tXYZ.Phi() <<", " << tXYZ.Rho() <<", " << tXYZ.Z() <<") "<<endl ;
         }
      }
      */

      leaves.nOutTimeHits = nOutTimeHits ;
      leaves.nHaloTrack   = nMinusHTrks + nPlusHTrks ;

      //std::vector<GlobalPoint> gp = cscData.GetCSCTrackImpactPositions() ;
      //cout<<" impact gp sz : "<< gp.size() <<" photon sz:"<< selectedPhotons.size() << endl ;

      //for (vector<GlobalPoint>::const_iterator it = gp.begin(); it != gp.end() ; ++it ) {
      //    double rho = sqrt(  (it->x()*it->x()) + (it->y()*it->y()) );
          //cout<<"    ==>  phi:"<<  it->phi() <<" rho: "<< rho << " z: "<< it->z() <<endl ;
      //    leaves.haloPhi = it->phi() ;
      //    leaves.haloRho = rho ;
 
          /*
          for ( size_t i=0; i< selectedPhotons.size() ; i++ ) { 
              double dPhi = fabs( selectedPhotons[i]->phi() - it->phi() ) ;
              double rhoG = sqrt(  ( selectedPhotons[i]->p4().x()*selectedPhotons[i]->p4().x() ) 
                                 + ( selectedPhotons[i]->p4().y()*selectedPhotons[i]->p4().y() )  );
              double dRho = fabs( rho - rhoG ) ;
              cout<<"    ("<< i << ") ==> dPhi : "<< dPhi <<" dRho : "<< dRho << endl ;
          }
          */
      //}
   }

}

void DPAnalysis::TriggerTagging( Handle<edm::TriggerResults> triggers, const edm::TriggerNames& trgNameList, int RunId, vector<int>& firedTrig ) {

   if ( runID_ != RunId )  {
      for (size_t j=0; j< triggerPatent.size(); j++ ) firedTrig[j] = -1;

      // loop through trigger menu
      for ( size_t i =0 ; i < trgNameList.size(); i++ ) {
          string tName  = trgNameList.triggerName( i );
          // loop through desired triggers
          for ( size_t j=0; j< triggerPatent.size(); j++ )  {
              if ( strncmp( tName.c_str(), triggerPatent[j].c_str(), triggerPatent[j].size() ) ==0 ) {
                  firedTrig[j] = i;
                  //cout<<" Trigger Found ("<<j <<"):  "<<tName ;
                  //cout<<" Idx: "<< i <<" triggers "<<endl;
              }
          }
      }
      runID_ = RunId ;
   }

}

// current used method 
bool DPAnalysis::TriggerSelection( Handle<edm::TriggerResults> triggers, vector<int> firedTrigID ) {

   bool pass =false ;
   uint32_t trgbits = 0 ;
   for ( size_t i=0; i< firedTrigID.size(); i++ ) {
       if ( firedTrigID[i] == -1 ) continue ; 
       if ( triggers->accept( firedTrigID[i] ) == 1  ) trgbits |= ( 1 << i ) ;
       //`cout<<" ("<< i <<") Trigger Found : "<< firedTrigID[i] <<" pass ? "<< triggers->accept( firedTrigID[i] ) <<" trigbit = "<< trgbits << endl; 
   }

   if ( trgbits != 0 ) {
      pass = true ;
   }
   leaves.triggered = (int)(trgbits) ;
   targetTrig = (int)(trgbits) ;

   return pass ;
}


template<class object >
bool DPAnalysis::GetTrgMatchObject( object, const edm::Event& iEvent,  InputTag inputProducer_ ) {

    bool findMatch = false ;
    // Get the input collection
    Handle<edm::View<object> > candHandle;
    iEvent.getByLabel( inputProducer_ , candHandle);

    Handle<trigger::TriggerEvent> trgEvent;
    iEvent.getByLabel( trigEvent, trgEvent);

    // get trigger object
    const trigger::TriggerObjectCollection& TOC( trgEvent->getObjects() );

    // find how many objects there are
    unsigned int nobj=0;
    double mindR = 999. ;
    double minP4[5] = { 0, 0, 0, 0, 0 } ;
    bool testPho = false ;
    bool testMET = false ;
    for( typename edm::View< object>::const_iterator j = candHandle->begin(); j != candHandle->end(); ++j, ++nobj) {

       // find out what filter is  
       //cout<<" trg : "<< inputProducer_.label() <<" filter sz: "<< trgEvent->sizeFilters() <<" nObj:"<< nobj << endl;
       for ( size_t ia = 0; ia < trgEvent->sizeFilters(); ++ia ) {

           // get the hlt filter
	   string fullname = trgEvent->filterTag(ia).encode();
	   size_t p = fullname.find_first_of(':');
	   string filterName = ( p != std::string::npos) ? fullname.substr(0, p) : fullname ; 
           //cout<<"    ... filterName : " << fullname <<" ==> "<< filterName << endl ;

	   bool trigPhoton = false ;
	   bool trigPfMet  = false ;
	   if ( strncmp( filterName.c_str(), "hltPhoton65CaloIdVLIsoLTrackIsoFilter", filterName.size() ) ==0 ) trigPhoton = true ;
	   if ( strncmp( filterName.c_str(), "hltPFMET25", filterName.size() ) ==0 ) trigPfMet = true ;
	   if ( strncmp( filterName.c_str(), "hltPFMET25Filter", filterName.size() ) ==0 ) trigPfMet = true ;

	   testPho = ( strncmp( inputProducer_.label().c_str(), "myphotons", 9 ) == 0 ) ;
	   testMET = ( strncmp( inputProducer_.label().c_str(), "pfMet", 5 ) == 0 ) ;
	   if ( !trigPhoton && !trigPfMet ) continue ;
	   if (  trigPhoton && !testPho   ) continue ;
	   if (  trigPfMet  && !testMET   ) continue ;

	   //if ( trigPhoton ) cout<<" Photon filter name: "<< filterName ;
	   //if ( trigPfMet  ) cout<<" PfMet  filter name: "<< filterName ;

	   // Get cut decision for each candidate
           const trigger::Keys& KEYS( trgEvent->filterKeys( ia ) );
	   int nKey = (int) KEYS.size() ;
	   //cout<<" nKey : "<< nKey << endl ;
           
	   for (int ipart = 0; ipart != nKey; ++ipart) { 
               const trigger::TriggerObject& TO = TOC[KEYS[ipart]];       
	       double dRval = deltaR( j->eta(), j->phi(), TO.eta(), TO.phi() );
	       //cout<<" Reco pt:"<< j->pt() <<"   TO pt: "<< TO.pt() <<endl ;
	       //cout<<"   -- > dR = " << dRval <<endl ;
               if ( dRval < mindR  ) { 
                  mindR = dRval ;
                  minP4[0] = TO.px() ;
                  minP4[1] = TO.py() ;
                  minP4[2] = TO.pz() ;
                  minP4[3] = TO.energy() ;
                  minP4[4] = TO.pt() ;
               }
           }
       }
    }
   
    if ( mindR < 0.5 ) findMatch = true ;
    if ( mindR < 9. && testMET ) {
       leaves.t_metPx = minP4[0] ;
       leaves.t_metPy = minP4[1] ;
       leaves.t_met   = minP4[4] ;
       leaves.t_metdR = mindR ;
       //if ( mindR < 0.5 ) cout<<" ** GOT MET => "<< minP4[4] << endl ;
    }   
    if ( mindR < 9. && testPho ) {
       leaves.t_phoPx = minP4[0] ;
       leaves.t_phoPy = minP4[1] ;
       leaves.t_phoPz = minP4[2] ;
       leaves.t_phoE  = minP4[3] ;
       leaves.t_phodR = mindR    ;
       //if ( mindR < 0.5 ) cout<<" ** GOT Photon => "<< minP4[4] << endl ;
    }   

    return findMatch ;
}

void DPAnalysis::PrintTriggers( const edm::Event& iEvent ) {

   Handle<edm::TriggerResults> triggers;
   iEvent.getByLabel( trigSource, triggers );

   cout<<" ** Trigger size = "<< triggers->size() <<endl;
   const edm::TriggerNames& trgNames = iEvent.triggerNames( *triggers );

   for ( size_t i =0 ; i < trgNames.size(); i++ ) {
       string tName  = trgNames.triggerName( i );
       int trgIndex  = trgNames.triggerIndex(tName);
       int trgResult = triggers->accept(trgIndex);
       cout<<" name: "<< tName <<" ("<< i <<")  idx:"<< trgIndex <<"  accept:"<< trgResult <<endl;
       for ( size_t j=0; j< triggerPatent.size(); j++) {
           if ( strncmp( tName.c_str(), triggerPatent[j].c_str(), triggerPatent[j].size() ) ==0 ) {
              //TriggerName = tName ;
              cout<<" Trigger Found : "<< tName <<" accepted ? "<< triggers->accept(i) <<endl;
           }
       }
       //string triggered = triggers->accept(i) ? "Yes" : "No" ;
       //cout<<" path("<<i<<") accepted ? "<< triggered ;
   }
}


void DPAnalysis::Track_Z0( Handle<reco::TrackCollection> tracks ) {

    //vector<TrkInfo> trkColl ;
    double in_trk = 0 ;
    double out_trk = 0 ;
    for (reco::TrackCollection::const_iterator it = tracks->begin(); it != tracks->end(); it++ )  {
        if ( fabs(it->d0()) >= vtxCuts[2] ) continue ;
         /*
         TrkInfo tf ;
         tf.dz  = it->dz() ;
         tf.dsz = it->dsz() ;
         tf.d0  = it->d0() ;
         tf.pt  = it->pt() ;
         tf.vz  = it->vz() ;
         tf.vr  = sqrt( (it->vx()*it->vx()) + (it->vy()*it->vy()) ) ;
         trkColl.push_back(tf) ;
         */
         h_z0->Fill( it->dz() ) ;
         if ( fabs( it->dz() ) < 20. )  in_trk += 1. ;
         else                          out_trk += 1. ;
         /*
         int ibin = (int)( (it->dz() + 155. ) / 10.) + 1;
         if ( it->dz() < -155. ) ibin = 0 ;
         if ( it->dz() >  155. ) ibin = 32 ;
         leaves.nTrkZ0[ ibin ] += 1 ;
         */
    }
    leaves.z0Ratio = ( in_trk > 0. ) ? (out_trk / in_trk) : -1 ;
    /*
    sort( trkColl.begin(), trkColl.end(), Z0Decreasing ) ;
    
    //printf("=============================== \n") ;
    vector<TrkInfo> selTrk ;
    double ddz = 99. ;
    for (size_t i=0 ; i < trkColl.size() ; i++) {
        //printf( "= dZ = %f , dsz: %f, d0: %f, pt: %f ", trkColl[i].dz , trkColl[i].dsz, trkColl[i].d0, trkColl[i].pt ) ;
        //printf( " v_z: %f , v_r: %f ", trkColl[i].vz , trkColl[i].vr ) ;
        ddz = ( i > 0 ) ? fabs( trkColl[i].dz - selTrk[ selTrk.size() -1 ].dz ) : 99. ;
        if ( ddz > 1. ) {
           selTrk.push_back( trkColl[i] ) ;
           //printf("  -- add \n") ;
        } else  {
           if (  trkColl[i].pt >  selTrk[ selTrk.size() -1 ].pt ) {
              //printf("  \n") ;
              //printf("  =>  delete track : dZ: %f, pt: %f \n", selTrk[ selTrk.size()-1].dz ,  selTrk[ selTrk.size()-1].pt ) ;
              selTrk.erase( selTrk.end() - 1 );
              selTrk.push_back( trkColl[i] ) ;
              //printf("  =>  add    track : dZ: %f, pt: %f \n", selTrk[ selTrk.size()-1].dz ,  selTrk[ selTrk.size()-1].pt ) ;
           } else {
              //printf("  -- skip \n") ;
           }
        } 
    }
    printf("  ----------------- \n") ;
    for (size_t i=0 ; i < selTrk.size() ; i++) {
        printf( "* dZ = %f , dsz: %f, d0: %f, pt: %f \n", selTrk[i].dz , selTrk[i].dsz, selTrk[i].d0, selTrk[i].pt ) ;
    }
    */
}

bool DPAnalysis::VertexSelection( Handle<reco::VertexCollection> vtx ) {

    bool hasGoodVertex = true ;
    int totalN_vtx = 0 ;
    int i = 0 ;
    for (reco::VertexCollection::const_iterator v=vtx->begin();  v!=vtx->end() ; v++){

        if ( ! v->isValid() ||  v->isFake() ) continue ;

        int ntrk = 0;
        for (reco::Vertex::trackRef_iterator itrk = v->tracks_begin(); itrk != v->tracks_end(); ++itrk) {
            double d0 = (*itrk)->d0() ;
            if ( d0 >= vtxCuts[2] ) continue ;
            ntrk++ ;
        }
        //if ( fabs(v->z()) >= vtxCuts[0] ) continue ; 
        if (   v->ndof()   < vtxCuts[1] ) continue ;

        // counting real number of vertices
        totalN_vtx++ ;

        if ( i >= MAXVTX ) continue ;
	leaves.vtxNTracks[i] = ntrk ;
	leaves.vtxChi2[i]    = v->normalizedChi2()  ;
	leaves.vtxNdof[i]    = v->ndof()   ;
	leaves.vtxRho[i]     = sqrt( ( v->y()*v->y() ) + ( v->x()*v->x() ) );
	leaves.vtxZ[i]       = v->z()  ;
        i++ ;

     }

     leaves.nVertices = i ;
     leaves.totalNVtx = totalN_vtx ;
 
     if ( i < 1 )   hasGoodVertex = false ;
     return hasGoodVertex ;
}

bool DPAnalysis::PhotonSelection( Handle<reco::PhotonCollection> photons, Handle<EcalRecHitCollection> recHitsEB, Handle<EcalRecHitCollection> recHitsEE, Handle<reco::TrackCollection> tracks, vector<const reco::Photon*>& selectedPhotons ) {

   for ( int i=0; i<7; i++) gcounter[i] = 0 ;
   int k= 0 ;
   double maxPt = 0 ;
   double met_dx(0), met_dy(0) ; 
   for(reco::PhotonCollection::const_iterator it = photons->begin(); it != photons->end(); it++) {

       // calculate met uncertainty
       if ( it->pt() < 10. ) continue ;
       met_dx += it->px() ;
       met_dy += it->py() ;
       //double ptscale = ( it->isEB() ) ? 0.994 : 0.985;
       double ptscale = ( it->isEB() ) ? 1.006 : 1.015 ;
       met_dx -= ( it->px() * ptscale ) ;
       met_dy -= ( it->py() * ptscale ) ;

       gcounter[0]++ ;
       // fiducial cuts
       if ( k >= MAXPHO ) continue ;
       gcounter[1]++ ;
       if ( it->pt() < photonCuts[0] || fabs( it->eta() ) > photonCuts[1] ) continue ;
       gcounter[2]++ ;
       //float hcalIsoRatio = it->hcalTowerSumEtConeDR04() / it->pt() ;
       //if  ( ( hcalIsoRatio + it->hadronicOverEm() )*it->energy() >= 6.0 ) continue ;
      
       // pixel veto - replaced by Conversion Veto
       //if ( it->hasPixelSeed() ) continue ;
       //if ( ConversionVeto( &(*it) ) ) cout<<" Got Conversion Case !! " << endl ; 
       if ( ConversionVeto( &(*it) ) ) continue; 
       gcounter[3]++ ;

       // S_Minor Cuts from the seed cluster
       reco::CaloClusterPtr SCseed = it->superCluster()->seed() ;
       const EcalRecHitCollection* rechits = ( it->isEB()) ? recHitsEB.product() : recHitsEE.product() ;

       Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*SCseed, *rechits);
       float sMin =  moments.sMin  ;
       float sMaj =  moments.sMaj  ;

       // seed Time 
       pair<DetId, float> maxRH = EcalClusterTools::getMaximum( *SCseed, rechits );
       DetId seedCrystalId = maxRH.first;
       EcalRecHitCollection::const_iterator seedRH = rechits->find(seedCrystalId);
       float seedTime    = (float)seedRH->time();
       float seedTimeErr = (float)seedRH->timeError();
       float swissX = EcalTools::swissCross( seedCrystalId, *rechits , 0., true ) ;
       float seedE       = seedRH->energy() ;
       // sMin and sMaj cuts
       if ( sMaj  > photonCuts[2] ) continue ;
       gcounter[4]++ ;
       if ( sMin <= photonCuts[3] || sMin >= photonCuts[4] ) continue ;
       gcounter[5]++ ;

       // Isolation Cuts 
       float ecalSumEt = it->ecalRecHitSumEtConeDR04();
       float hcalSumEt = it->hcalTowerSumEtConeDR04();
       float trkSumPt  = it->trkSumPtSolidConeDR04();  

       //bool trkIso  = ( ( trkSumPt / it->pt())     < photonIso[0] ) ; 
       //bool ecalIso = ( (ecalSumEt / it->energy()) < photonIso[2] && ecalSumEt < photonIso[1] ) ; 
       //bool hcalIso = ( (hcalSumEt / it->energy()) < photonIso[4] && hcalSumEt < photonIso[3] ) ; 
       //if ( !trkIso || !ecalIso || !hcalIso ) continue ;

       // PF Isolation - still not working , useless
       //float cHadIso = max( it->chargedHadronIso() - RhoCorrection( 1, it->eta() ), 0. );
       //float nHadIso = max( it->neutralHadronIso() - RhoCorrection( 2, it->eta() ), 0. );
       //float photIso = max( it->photonIso()        - RhoCorrection( 3, it->eta() ), 0. );
       //printf(" cHad: %.3f, nHad: %.3f, phot: %.3f \n", it->chargedHadronIso(),  it->neutralHadronIso(),  it->photonIso() ) ;

       // Track Veto 
       int nTrk = 0 ;
       double minDR = 99. ;
       double trkPt = 0 ;
       for (reco::TrackCollection::const_iterator itrk = tracks->begin(); itrk != tracks->end(); itrk++ )  {
           if ( itrk->pt() < 3. ) continue ;
	   LorentzVector trkP4( itrk->px(), itrk->py(), itrk->pz(), itrk->p() ) ;
	   double dR =  ROOT::Math::VectorUtil::DeltaR( trkP4 , it->p4()  ) ;
           if ( dR < minDR ) {
              minDR = dR ;
              trkPt = itrk->pt() ;
           }
	   if ( dR < photonCuts[5] )  nTrk++ ;
       }
       if ( nTrk > 0 ) continue ;
       gcounter[6]++ ;

       // check leading photon pt  
       maxPt = ( it->pt() > maxPt ) ? it->pt() : maxPt ;

       // Timing Calculation
       pair<double,double> AveXtalTE =  ClusterTime( it->superCluster(), recHitsEB , recHitsEE );

       PhoInfo phoTmp ;
       phoTmp.t      = AveXtalTE.first ;
       phoTmp.dt     = AveXtalTE.second ;
       phoTmp.nchi2  = 0 ;
       phoTmp.nxtals = 0 ;
       phoTmp.nBC    = 0 ;
       //cout<<" 1st xT : "<< aveXtalTime <<"  xTE : "<< aveXtalTimeErr << endl;
       // Only use the seed cluster
       //if ( seedTime > 5. ) debugT = true ;
       if ( debugT ) printf("===== seedT: %.2f, 1st Ave.T: %.2f =====\n", seedTime,  AveXtalTE.first ) ;
       ClusterTime( it->superCluster(), recHitsEB , recHitsEE, phoTmp );
       //cout<<" 2nd xT : "<< aveXtalTime <<"  xTE : "<< aveXtalTimeErr << endl;
       leaves.aveTime1[k]     = phoTmp.t ;    // weighted ave. time of seed cluster
       leaves.aveTimeErr1[k]  = phoTmp.dt ;
       leaves.timeChi2[k]     = phoTmp.nchi2 ;
       leaves.nXtals[k]       = phoTmp.nxtals ;
       leaves.nBC[k]          = phoTmp.nBC ;
       leaves.seedSwissX[k]   = swissX ;
 
       debugT = false ;
       // refine the timing to exclude out-of-time xtals
       ClusterTime( it->superCluster(), recHitsEB , recHitsEE, phoTmp );
       //cout<<" 3rd xT : "<< aveXtalTime <<"  xTE : "<< aveXtalTimeErr << endl;

       leaves.phoPx[k] = it->p4().Px() ;
       leaves.phoPy[k] = it->p4().Py() ;
       leaves.phoPz[k] = it->p4().Pz() ;
       leaves.phoE[k]  = it->p4().E() ;
       //leaves.phoHoverE[k]  = it->hadronicOverEm() ;
       leaves.phoHoverE[k]  = it->hadTowOverEm() ;

       leaves.phoEcalIso[k] = ecalSumEt ;
       leaves.phoHcalIso[k] = hcalSumEt ;
       leaves.phoTrkIso[k]  = trkSumPt ;

       leaves.dR_TrkPho[k]  = minDR ;
       leaves.pt_TrkPho[k]  = trkPt ;
       leaves.sMinPho[k]    = sMin ;
       leaves.sMajPho[k]    = sMaj ;

       leaves.seedE[k]        = seedE ;
       leaves.seedTime[k]     = seedTime ;
       leaves.seedTimeErr[k]  = seedTimeErr ;
       leaves.aveTime[k]      = phoTmp.t ;       // weighted ave. time of seed cluster 
       leaves.aveTimeErr[k]   = phoTmp.dt ;
       leaves.sigmaEta[k]     = it->sigmaEtaEta() ;
       leaves.sigmaIeta[k]    = it->sigmaIetaIeta() ;
       selectedPhotons.push_back( &(*it) ) ;
       k++ ;
   }
   leaves.nPhotons = k ;
   leaves.met_dx3  += met_dx ;
   leaves.met_dy3  += met_dy ;
   //leaves.nPhotons = (int)( selectedPhotons.size() ) ;

   if ( selectedPhotons.size() > 0 && maxPt >= photonCuts[7] )  return true ; 
   else                               return false ;    

}

// !!! Only work for RECO - dtSegment is not available for AOD 
bool DPAnalysis::CosmicRayMatch( Handle<DTRecSegment4DCollection> dtSeg, vector<const reco::Photon*>& selectedPhotons, const EventSetup& iSetup ) {

   ESHandle<DTGeometry> dtGeom;
   iSetup.get<MuonGeometryRecord>().get(dtGeom);

   //bool cosmicmatch = false ;
   for ( size_t i=0; i< selectedPhotons.size() ; i++ ) {

       double dPhi = 99. ;
       double dEta = 99. ;
       double dR   = 99. ;
       for (DTRecSegment4DCollection::const_iterator it = dtSeg->begin(); it != dtSeg->end(); it++) {
           if ( !it->isValid()) continue ;

           // Get the corresponding DTChamber
           DetId id = it->geographicalId();
           DTChamberId chamberId(id.rawId());
           if ( chamberId.station() > 1 ) continue ;  // only look at segment from inner most chambers
           const DTChamber* dtchamber = dtGeom->chamber( chamberId ) ;

           // Get segment position in DT
	   GlobalPoint  gp = dtchamber->toGlobal( it->localPosition()  );
           //double gpMag = sqrt( (gp.x()*gp.x()) + (gp.y()*gp.y()) + (gp.z()*gp.z()) ) ;
	   //LorentzVector segP4( gp.x(), gp.y(), gp.z(), gpMag ) ;

           // Get segment direction in DT
           GlobalVector gv = dtchamber->toGlobal( it->localDirection() ) ;
           double dx = gp.x() - selectedPhotons[i]->caloPosition().x() ;
           double dy = gp.y() - selectedPhotons[i]->caloPosition().y() ;
           double dz = gp.z() - selectedPhotons[i]->caloPosition().z() ;
           double dr = sqrt( (dx*dx) + (dy*dy) + (dz*dz) ) ;
           //printf(" dir( %f, %f, %f, %f ) \n ", gv.x() , gv.y(), gv.z() , gv.mag() ) ;

           // Project segment position to Ecal 
           double projX =  gp.x() - ( dr* gv.x() / gv.mag()) ;
           double projY =  gp.y() - ( dr* gv.y() / gv.mag()) ;
           double projZ =  gp.z() - ( dr* gv.z() / gv.mag()) ;
           double projR = sqrt( (projX*projX) + (projY*projY) + (projZ*projZ) ) ;

           LorentzVector segEBP4( projX, projY, projZ, projR ) ;

           double dR_   = ROOT::Math::VectorUtil::DeltaR( segEBP4, selectedPhotons[i]->p4() ) ;
           double dPhi_ = ROOT::Math::VectorUtil::DeltaPhi( segEBP4, selectedPhotons[i]->p4() ) ;
           double dEta_ = fabs( segEBP4.Eta() - selectedPhotons[i]->p4().Eta() ) ;
           
           if ( dR_ < dR ) {
               dR = dR_ ;
               dPhi = fabs( dPhi_ ) ;
               dEta = dEta_ ;
           }
       }
       leaves.dtdEta[i]  = dEta ;
       leaves.dtdPhi[i]  = dPhi ;
  }
  return true ;

}

// !!! Only work for RECO - cscSegment is not available for AOD 
bool DPAnalysis::BeamHaloMatch( Handle<CSCSegmentCollection> cscSeg, vector<const reco::Photon*>& selectedPhotons, const EventSetup& iSetup ) {

   ESHandle<CSCGeometry> cscGeom;
   iSetup.get<MuonGeometryRecord>().get(cscGeom);

   bool halomatch = false ;
   for ( size_t i=0; i< selectedPhotons.size() ; i++ ) {

       double dPhi    = 99. ;
       double cscRho  = -1 ;
       double cscTime = 99. ;
       for (CSCSegmentCollection::const_iterator it = cscSeg->begin(); it != cscSeg->end(); it++) {
           if ( !it->isValid() ) continue ;
           CSCDetId DetId = it->cscDetId();
	   const CSCChamber* cscchamber = cscGeom->chamber( DetId );
	   GlobalPoint  gp = cscchamber->toGlobal( it->localPosition()  );
           double gpMag = sqrt( (gp.x()*gp.x()) + (gp.y()*gp.y()) + (gp.z()*gp.z()) ) ;
	   LorentzVector segP4( gp.x(), gp.y(), gp.z(), gpMag ) ;

           if ( fabs( gp.eta() ) < 1.6  ) continue ;
           //double dEta_ = fabs( gp.eta() - selectedPhotons[i]->eta() )  ;
           //dEta = ( dEta_ < dEta ) ? dEta_ : dEta ;
           //double dPhi_ = fabs( gp.phi() - selectedPhotons[i]->phi() )  ;
           double dPhi_ = ROOT::Math::VectorUtil::DeltaPhi( segP4, selectedPhotons[i]->p4() ) ;
	   double rho = sqrt( (gp.x()*gp.x()) + (gp.y()*gp.y()) ) ;
           if ( fabs(dPhi_) < dPhi ) {
              dPhi = fabs( dPhi_ ) ;
              cscRho  = rho ;
              cscTime = it->time() ;
           }
       }
       if ( dPhi < 0.1 )  halomatch = true ;
       leaves.cscRho[i]   = cscRho ;
       leaves.cscdPhi[i]  = dPhi ;
       leaves.cscTime[i]  = cscTime ;
   }

   return halomatch ;
}


// return time, timeError
pair<double,double> DPAnalysis::ClusterTime( reco::SuperClusterRef scRef, Handle<EcalRecHitCollection> recHitsEB, Handle<EcalRecHitCollection> recHitsEE, float ECut, bool useAllClusters ) {

  const EcalIntercalibConstantMap& icalMap = ical->getMap();
  float adcToGeV_EB = float(agc->getEBValue());
  float adcToGeV_EE = float(agc->getEEValue());

  double xtime = 0 ;
  double xtimeErr = 0 ;

  // 1. loop all the basic clusters 
  for ( reco::CaloCluster_iterator  clus = scRef->clustersBegin() ;  clus != scRef->clustersEnd();  ++clus) {

      // only use seed basic cluster  
      if ( *clus != scRef->seed() && !useAllClusters ) continue ;
      // GFdoc clusterDetIds holds crystals that participate to this basic cluster 
      // 2. loop on xtals in cluster
      std::vector<std::pair<DetId, float> > clusterDetIds = (*clus)->hitsAndFractions() ; //get these from the cluster
      //cout<<" --------------- "<<endl ;
      int nXtl = 0 ;
      for (std::vector<std::pair<DetId, float> >::const_iterator detitr = clusterDetIds.begin () ; 
           detitr != clusterDetIds.end () ; ++detitr) { 

             // Here I use the "find" on a recHit collection... I have been warned...   (GFdoc: ??)
   	     // GFdoc: check if DetId belongs to ECAL; if so, find it among those if this basic cluster
    	     if ( (detitr -> first).det () != DetId::Ecal)  { 
   	          cout << " det is " << (detitr -> first).det () << " (and not DetId::Ecal)" << endl ;
	          continue ;
	     }
             bool isEB = ( (detitr -> first).subdetId () == EcalBarrel)  ? true : false ;
	   
	     // GFdoc now find it!
	     EcalRecHitCollection::const_iterator thishit = (isEB) ? recHitsEB->find( (detitr->first) ) : recHitsEE->find( (detitr->first) );
	     if (thishit == recHitsEB->end () &&  isEB )  continue ;
	     if (thishit == recHitsEE->end () && !isEB )  continue ;

	     // GFdoc this is one crystal in the basic cluster
	     EcalRecHit myhit = (*thishit) ;
	   
             double thisamp  = myhit.energy () ;
             EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(detitr->first);
             EcalIntercalibConstant icalconst = 1;
             if( icalit!=icalMap.end() ) {
               icalconst = (*icalit);
             } 
             // get laser coefficient
             float lasercalib = laser->getLaserCorrection( detitr->first, eventTime );
             float adcToGeV = ( isEB ) ? adcToGeV_EB : adcToGeV_EE ;
             // discard rechits with A/sigma < 12
             if ( thisamp/(icalconst*lasercalib*adcToGeV) < (1.1*12) ) continue;
             // SIC Feb 14 2011 -- Add check on RecHit flags (takes care of spike cleaning in 42X)
             //if ( !( myhit.checkFlag(EcalRecHit::kGood) || myhit.checkFlag(EcalRecHit::kOutOfTime) || 
             //       myhit.checkFlag(EcalRecHit::kPoorCalib)  ) )  continue;
             if ( !( myhit.checkFlag(EcalRecHit::kGood) || myhit.checkFlag(EcalRecHit::kOutOfTime) ) )  continue;
             //GlobalPoint pos = theGeometry->getPosition( (myhit).detid() );
             if ( myhit.energy() < ECut ) continue ; 
             nXtl++ ;

             // time and time correction
	     double thistime = myhit.time();

             // get time error 
             double xtimeErr_ = ( myhit.isTimeErrorValid() ) ?  myhit.timeError() : 999 ;
 
             xtime     += thistime / pow( xtimeErr_ , 2 ) ;
             xtimeErr  += 1/ pow( xtimeErr_ , 2 ) ;
      }
      //cout<<" total Xtl = " << nXtl << endl ;
  }
  double wAveTime = ( xtimeErr < 0.000001 ) ? 0. : xtime / xtimeErr ;
  double wAveTimeErr =  ( xtimeErr < 0.000001 ) ? 0. : 1. / sqrt( xtimeErr) ;
  pair<double, double> wAveTE( wAveTime, wAveTimeErr ) ;
  return wAveTE ;  

}

// re-calculate time and timeError as well as normalized chi2
//void DPAnalysis::ClusterTime( reco::SuperClusterRef scRef, Handle<EcalRecHitCollection> recHitsEB, Handle<EcalRecHitCollection> recHitsEE, double& aveTime, double& aveTimeErr, double& nChi2, bool useAllClusters ) {

void DPAnalysis::ClusterTime( reco::SuperClusterRef scRef, Handle<EcalRecHitCollection> recHitsEB, Handle<EcalRecHitCollection> recHitsEE, PhoInfo& phoTmp, float ECut, bool useAllClusters ) {

  const EcalIntercalibConstantMap& icalMap = ical->getMap();
  float adcToGeV_EB = float(agc->getEBValue());
  float adcToGeV_EE = float(agc->getEEValue());

  double xtime    = 0 ;
  double xtimeErr = 0 ;
  double chi2_bc  = 0 ;
  double ndof     = 0 ;
  int    nBC      = 0 ;
  int    nXtl     = 0 ;
  for ( reco::CaloCluster_iterator  clus = scRef->clustersBegin() ;  clus != scRef->clustersEnd();  ++clus) {

      nBC++ ;
      // only use seed basic cluster  
      //bool isSeed = ( *clus == scRef->seed() ) ;
      if ( *clus != scRef->seed() && !useAllClusters ) continue ;

      // GFdoc clusterDetIds holds crystals that participate to this basic cluster 
      //loop on xtals in cluster
      std::vector<std::pair<DetId, float> > clusterDetIds = (*clus)->hitsAndFractions() ; //get these from the cluster
      for (std::vector<std::pair<DetId, float> >::const_iterator detitr = clusterDetIds.begin () ; 
           detitr != clusterDetIds.end () ; ++detitr) { 
	      // Here I use the "find" on a recHit collection... I have been warned...   (GFdoc: ??)
   	      // GFdoc: check if DetId belongs to ECAL; if so, find it among those if this basic cluster
    	     if ( (detitr -> first).det () != DetId::Ecal)  { 
   	          cout << " det is " << (detitr -> first).det () << " (and not DetId::Ecal)" << endl ;
	          continue ;
	     }
             bool isEB = ( (detitr -> first).subdetId () == EcalBarrel)  ? true : false ;
	   
	     // GFdoc now find it!
	     EcalRecHitCollection::const_iterator thishit = (isEB) ? recHitsEB->find( (detitr->first) ) : recHitsEE->find( (detitr->first) );
	     if (thishit == recHitsEB->end () &&  isEB )  continue ;
	     if (thishit == recHitsEE->end () && !isEB )  continue ;

	     // GFdoc this is one crystal in the basic cluster
	     EcalRecHit myhit = (*thishit) ;
	   
             // SIC Feb 14 2011 -- Add check on RecHit flags (takes care of spike cleaning in 42X)
             if ( !( myhit.checkFlag(EcalRecHit::kGood) || myhit.checkFlag(EcalRecHit::kOutOfTime) || 
                    myhit.checkFlag(EcalRecHit::kPoorCalib)  ) )  continue;

             //if ( myhit.checkFlag(EcalRecHit::kWeird) || myhit.checkFlag(EcalRecHit::kDiWeird) ) continue ;
             //bool gotSpike = ( myhit.checkFlag(EcalRecHit::kWeird) || myhit.checkFlag(EcalRecHit::kDiWeird) )  ;
             // swiss cross cleaning 
             //float swissX = (isEB) ? EcalTools::swissCross(detitr->first, *recHitsEB , 0., true ) : 
             //                        EcalTools::swissCross(detitr->first, *recHitsEE , 0., true ) ;

             // thisamp is the EB amplitude of the current rechit
             double thisamp  = myhit.energy () ;
             EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(detitr->first);
             EcalIntercalibConstant icalconst = 1;
             if( icalit!=icalMap.end() ) {
               icalconst = (*icalit);
             } else {
               edm::LogError("EcalTimePhyTreeMaker") << "No intercalib const found for xtal " << (detitr->first).rawId();
             }

             // get laser coefficient
             float lasercalib = laser->getLaserCorrection( detitr->first, eventTime );
             float adcToGeV = ( isEB ) ? adcToGeV_EB : adcToGeV_EE ;
             // discard rechits with A/sigma < 12
             if ( thisamp/(icalconst*lasercalib*adcToGeV) < (1.1*12) ) continue;
             if ( myhit.energy() < ECut ) continue ; 
             //GlobalPoint pos = theGeometry->getPosition((myhit).detid());

             // time and time correction
	     double thistime = myhit.time();
	     //thistime += theTimeCorrector_.getCorrection((float) thisamp/(icalconst*lasercalib*adcToGeV), pos.eta()  );

             // get time error 
             double xtimeErr_ = ( myhit.isTimeErrorValid() ) ?  myhit.timeError() : 999 ;

             // calculate chi2 for the BC of the seed
             double chi2_x = pow( ((thistime - phoTmp.t) / xtimeErr_ ) , 2 ) ; 

             chi2_bc += chi2_x ;
             ndof += 1 ;
             nXtl++ ;
             // remove un-qualified hits 
             if ( fabs ( thistime - phoTmp.t ) > 3.*phoTmp.dt ) continue ;
 
             xtime     += thistime / pow( xtimeErr_ , 2 ) ;
             xtimeErr  += 1/ pow( xtimeErr_ , 2 ) ;
      }
  }
  if ( debugT ) printf("--- sum_chi2: %.2f, ndof: %.1f norm_chi2: %.2f ---\n", chi2_bc, ndof, chi2_bc/ndof );
  // update ave. time and error

  phoTmp.t     = ( xtimeErr < 0.000001 ) ? 0. : xtime / xtimeErr ;
  phoTmp.dt    = ( xtimeErr < 0.000001 ) ? 0. : 1. / sqrt( xtimeErr) ;
  phoTmp.nchi2 = ( ndof != 0 ) ? chi2_bc / ndof : 9999999 ;     
  phoTmp.nxtals = nXtl ;
  phoTmp.nBC    = nBC ;

}

bool DPAnalysis::JetSelection( Handle<reco::PFJetCollection> jets, 
                      vector<const reco::Photon*>& selectedPhotons, vector<const reco::PFJet*>& selectedJets) {

   int k = 0 ;
   for(reco::PFJetCollection::const_iterator it = jets->begin(); it != jets->end(); it++) {
       // fiducial cuts
       if ( it->pt() < jetCuts[0] || fabs( it->eta() ) > jetCuts[1] ) continue ;

       // Jet ID cuts
       /*
       if ( it->numberOfDaughters() < 2 )               continue ;
       if ( it->chargedEmEnergyFraction() >= 0.99 )     continue ;
       if ( it->neutralHadronEnergyFraction() >= 0.99 ) continue ;
       if ( it->neutralEmEnergyFraction() >= 0.99 )     continue ;
       if ( fabs( it->eta() ) < 2.4 && it->chargedHadronEnergyFraction() <=0 ) continue ;
       if ( fabs( it->eta() ) < 2.4 && it->chargedMultiplicity() <=0 ) continue ;
       */
       // dR cuts 
       double dR = 999 ;
       for (size_t j=0; j < selectedPhotons.size(); j++ ) {
           double dR_ =  ROOT::Math::VectorUtil::DeltaR( it->p4(), selectedPhotons[j]->p4() ) ;
           if ( dR_ < dR ) dR = dR_ ;
       }
       if ( dR <= jetCuts[2] ) continue ;

       vector<double> uncV = JECUncertainty( it->pt(), it->eta(), jecUnc ) ;

       if ( k >= MAXJET ) break ;

       selectedJets.push_back( &(*it) ) ;
       leaves.jetPx[k] = it->p4().Px() ;
       leaves.jetPy[k] = it->p4().Py() ;
       leaves.jetPz[k] = it->p4().Pz() ;
       leaves.jetE[k]  = it->p4().E()  ;
       leaves.jetNDau[k] = it->numberOfDaughters() ;
       leaves.jetCM[k]   = it->chargedMultiplicity() ;
       leaves.jetCEF[k]  = it->chargedEmEnergyFraction() ;
       leaves.jetNHF[k]  = it->neutralHadronEnergyFraction() ;  
       leaves.jetNEF[k]  = it->neutralEmEnergyFraction() ;
       //leaves.jecUncU[k]  = uncV[0] ;
       //leaves.jecUncD[k]  = uncV[1] ;
       leaves.jecUnc[k]  = uncV[2] ;
       k++ ;
   }
   leaves.nJets = (int)( selectedJets.size() ) ;

   if ( selectedJets.size() > 0 )  return true ; 
   else                            return false ;    

}

bool DPAnalysis::JetSelection( Handle< vector<pat::Jet> > patjets, vector<const reco::Photon*>& selectedPhotons, 
                                          vector< pat_Jet* >& selectedJets_ ) {

   int k = 0 ;
   double met_dx(0), met_dy(0) ;
   double met_sx(0), met_sy(0) ;
   double SF = 1. ;
   for (std::vector<pat::Jet>::const_iterator it = patjets->begin(); it != patjets->end(); it++) {
 
       bool passPtCut = it->pt() > jetCuts[0] ;
       // calculate JER uncertainty 
       double ptscale = 1 ;
       double dPt = 0 ;
       if ( !isData ) {
          const reco::GenJet* matchedGenJet = it->genJet() ;
	  if ( it->pt() < 10. ) continue ;
	  if ( matchedGenJet == NULL ) continue ;

	  met_dx += it->px() ;
	  met_dy += it->py() ;

	  // This is data/MC ratio
	  if ( fabs(it->eta()) < 0.5 ) SF = 1.052 ;
	  if ( fabs(it->eta()) >= 0.5 && fabs(it->eta()) < 1.1 ) SF = 1.057 ;
	  if ( fabs(it->eta()) >= 1.1 && fabs(it->eta()) < 1.7 ) SF = 1.096 ;
	  if ( fabs(it->eta()) >= 1.7 && fabs(it->eta()) < 2.3 ) SF = 1.134 ;
	  if ( fabs(it->eta()) >= 2.3 && fabs(it->eta()) < 2.5 ) SF = 1.288 ;
	  dPt = ( it->pt() - matchedGenJet->pt() )*(1-SF) ;
	  ptscale = max( 0.0, ( it->pt() + dPt)/it->pt() ) ;

	  met_dx -= ( it->px() * ptscale ) ;
	  met_dy -= ( it->py() * ptscale ) ;
       }
       if ( (it->pt() + dPt) > jetCuts[0] || (it->pt() - dPt) > jetCuts[0] ) passPtCut = true ;

       // Calculate JES uncertainty
       vector<double> uncV = JECUncertainty( it->pt(), it->eta(), jecUnc ) ;
       if ( (it->pt() + (uncV[2]*it->pt()) ) > jetCuts[0] || (it->pt() - (uncV[2]*it->pt()) ) > jetCuts[0] ) passPtCut = true ;
       met_sx += it->px() ;
       met_sy += it->py() ;
       double jes_sc = ( it->pt() + uncV[2] ) / it->pt() ;
       met_sx -= it->px() * jes_sc ;
       met_sy -= it->py() * jes_sc ;

       // Pt and Fiducial cuts - include those within JES and JER range
       if ( !passPtCut || fabs( it->eta() ) > jetCuts[1] ) continue ;

       // Jet ID cuts
       //if ( it->numberOfDaughters() < 2 )               continue ;
       //if ( it->chargedEmEnergyFraction() >= 0.99 )     continue ;
       //if ( it->neutralHadronEnergyFraction() >= 0.99 ) continue ;
       //if ( it->neutralEmEnergyFraction() >= 0.99 )     continue ;
       //if ( fabs( it->eta() ) < 2.4 && it->chargedHadronEnergyFraction() <=0 ) continue ;
       //if ( fabs( it->eta() ) < 2.4 && it->chargedMultiplicity() <=0 ) continue ;
       
       // dR cuts 
       double dR = 999 ;
       for (size_t j=0; j < selectedPhotons.size(); j++ ) {
           double dR_ =  ROOT::Math::VectorUtil::DeltaR( it->p4(), selectedPhotons[j]->p4() ) ;
           if ( dR_ < dR ) dR = dR_ ;
       }
       if ( dR <= jetCuts[2] ) continue ;


       if ( k >= MAXJET ) continue ;
       selectedJets_.push_back( &(*it) ) ;
       leaves.jetPx[k] = it->p4().Px() ;
       leaves.jetPy[k] = it->p4().Py() ;
       leaves.jetPz[k] = it->p4().Pz() ;
       leaves.jetE[k]  = it->p4().E()  ;
       leaves.jetNDau[k] = it->numberOfDaughters() ;
       leaves.jetCM[k]   = it->chargedMultiplicity() ;
       leaves.jetCEF[k]  = it->chargedEmEnergyFraction() ;
       leaves.jetNHF[k]  = it->neutralHadronEnergyFraction() ;  
       leaves.jetNEF[k]  = it->neutralEmEnergyFraction() ;
       //leaves.jecUncU[k]  = uncV[0] ;
       //leaves.jecUncD[k]  = uncV[1] ;
       leaves.jecUnc[k]  = uncV[2] ;
       leaves.jerUnc[k]  = dPt ;
       k++ ;
   }
   leaves.nJets = (int)( selectedJets.size() ) ;
   leaves.met_dx1    = met_dx  ;
   leaves.met_dy1    = met_dy  ;
   leaves.met_dx2    = met_sx  ;
   leaves.met_dy2    = met_sy  ;

   if ( selectedJets.size() > 0 )  return true ; 
   else                            return false ;    

}

bool DPAnalysis::JetSelection( Handle< vector<pat::Jet> > patjets, vector<const reco::Photon*>& selectedPhotons, 
                               Handle<reco::SuperClusterCollection> scEB, Handle<reco::SuperClusterCollection> scEE, 
                               Handle<EcalRecHitCollection> recHitsEB, Handle<EcalRecHitCollection> recHitsEE,
                               vector< pat_Jet* >& selectedJets_ ) {

   int k = 0 ;
   double met_dx(0), met_dy(0) ;
   double met_sx(0), met_sy(0) ;
   double SF = 1. ;
   for ( std::vector<pat::Jet>::const_iterator it = patjets->begin(); it != patjets->end(); it++ ) {
 
       bool passPtCut = it->pt() > jetCuts[0] ;
       // calculate JER uncertainty 
       double ptscale = 1 ;
       double dPt = 0 ;
       double genJetP4[5] = { 0,0,0,0, -1 };
       if ( !isData ) {
          const reco::GenJet* matchedGenJet = it->genJet() ;
	  if ( it->pt() < 10. ) continue ;
	  if ( matchedGenJet == NULL ) continue ;

	  met_dx += it->px() ;
	  met_dy += it->py() ;

	  // This is data/MC ratio
	  if ( fabs(it->eta()) < 0.5 ) SF = 1.052 ;
	  if ( fabs(it->eta()) >= 0.5 && fabs(it->eta()) < 1.1 ) SF = 1.057 ;
	  if ( fabs(it->eta()) >= 1.1 && fabs(it->eta()) < 1.7 ) SF = 1.096 ;
	  if ( fabs(it->eta()) >= 1.7 && fabs(it->eta()) < 2.3 ) SF = 1.134 ;
	  if ( fabs(it->eta()) >= 2.3 && fabs(it->eta()) < 2.5 ) SF = 1.288 ;
	  dPt = ( it->pt() - matchedGenJet->pt() )*(1-SF) ;
	  ptscale = max( 0.0, ( it->pt() + dPt)/it->pt() ) ;

	  met_dx -= ( it->px() * ptscale ) ;
	  met_dy -= ( it->py() * ptscale ) ;

          // find the corresponding genJets
	  genJetP4[0] = matchedGenJet->energy()  ;
	  genJetP4[1] = matchedGenJet->px() ;
	  genJetP4[2] = matchedGenJet->py() ;
	  genJetP4[3] = matchedGenJet->pz();
       }
       if ( (it->pt() + dPt) > jetCuts[0] || (it->pt() - dPt) > jetCuts[0] ) passPtCut = true ;

       // Calculate JES uncertainty
       vector<double> uncV = JECUncertainty( it->pt(), it->eta(), jecUnc ) ;
       if ( (it->pt() + (uncV[2]*it->pt()) ) > jetCuts[0] || (it->pt() - (uncV[2]*it->pt()) ) > jetCuts[0] ) passPtCut = true ;
       met_sx += it->px() ;
       met_sy += it->py() ;
       double jes_sc = ( it->pt() + uncV[2] ) / it->pt() ;
       met_sx -= it->px() * jes_sc ;
       met_sy -= it->py() * jes_sc ;

       // Pt and Fiducial cuts - include those within JES and JER range
       if ( !passPtCut || fabs( it->eta() ) > jetCuts[1] ) continue ;

       // Jet ID cuts
       //if ( it->numberOfDaughters() < 2 )               continue ;
       //if ( it->chargedEmEnergyFraction() >= 0.99 )     continue ;
       //if ( it->neutralHadronEnergyFraction() >= 0.99 ) continue ;
       //if ( it->neutralEmEnergyFraction() >= 0.99 )     continue ;
       //if ( fabs( it->eta() ) < 2.4 && it->chargedHadronEnergyFraction() <=0 ) continue ;
       //if ( fabs( it->eta() ) < 2.4 && it->chargedMultiplicity() <=0 ) continue ;
       
       // dR cuts - exclude fake jet which looks like photon
       double dR = 999 ;
       for (size_t j=0; j < selectedPhotons.size(); j++ ) {
           double dR_ =  ROOT::Math::VectorUtil::DeltaR( it->p4(), selectedPhotons[j]->p4() ) ;
           if ( dR_ < dR ) dR = dR_ ;
       }
       if ( dR <= jetCuts[2] ) continue ;


       int iSC = 0;
       float maxEnergy = 0 ;
       pair<double,double> jetTime ;
       for (reco::SuperClusterCollection::const_iterator scIt = scEB->begin() ; scIt != scEB->end(); scIt++) { 

           if ( fabs( it->p4().Eta() ) > 1.54 ) break ;
           double dEta = it->p4().Eta() - scIt->position().eta() ;
           double dPhi = it->p4().Phi() - scIt->position().phi() ;
           double dR = sqrt( dEta*dEta + dPhi*dPhi ) ;

           if ( dR < 0.3  && scIt->energy() > maxEnergy ) {
              jetTime = ClusterTime( reco::SuperClusterRef(scEB, iSC), recHitsEB, recHitsEE, 10, true ) ;
              maxEnergy = scIt->energy() ;
           } else {
              continue ;
           }
           iSC++ ;
       }

       iSC = 0 ;
       for (reco::SuperClusterCollection::const_iterator scIt = scEE->begin() ; scIt != scEE->end(); scIt++) { 

           if ( fabs( it->p4().Eta() ) < 1.4 ) break ;
           double dEta = it->p4().Eta() - scIt->position().eta() ;
           double dPhi = it->p4().Phi() - scIt->position().phi() ;
           double dR = sqrt( dEta*dEta + dPhi*dPhi ) ;

           pair<double,double> jetTime ;
           if ( dR < 0.3  && scIt->energy() > maxEnergy ) {
              jetTime = ClusterTime( reco::SuperClusterRef(scEE, iSC), recHitsEB, recHitsEE, 10, true ) ;
              maxEnergy = scIt->energy() ;
           } else {
              continue ;
           }
           iSC++ ;
       }

       if ( k >= MAXJET ) continue ;
       selectedJets_.push_back( &(*it) ) ;
       leaves.jetPx[k] = it->p4().Px() ;
       leaves.jetPy[k] = it->p4().Py() ;
       leaves.jetPz[k] = it->p4().Pz() ;
       leaves.jetE[k]  = it->p4().E()  ;
       leaves.jetNDau[k] = it->numberOfDaughters() ;
       leaves.jetCM[k]   = it->chargedMultiplicity() ;
       leaves.jetCEF[k]  = it->chargedEmEnergyFraction() ;
       leaves.jetNHF[k]  = it->neutralHadronEnergyFraction() ;  
       leaves.jetNEF[k]  = it->neutralEmEnergyFraction() ;
       leaves.jetTime[k]    = jetTime.first ;
       leaves.jetTimeErr[k] = jetTime.second ;
       leaves.jetSCE[k]     = maxEnergy ;
       //leaves.jecUncU[k]  = uncV[0] ;
       //leaves.jecUncD[k]  = uncV[1] ;
       leaves.jecUnc[k]  = uncV[2] ;
       leaves.jerUnc[k]  = dPt ;
       leaves.genJetE[k]  = genJetP4[0] ;
       leaves.genJetPx[k] = genJetP4[1] ;
       leaves.genJetPy[k] = genJetP4[2] ;
       leaves.genJetPz[k] = genJetP4[3] ;
       k++ ;
   }
   leaves.nJets = (int)( selectedJets.size() ) ;
   leaves.met_dx1    = met_dx  ;
   leaves.met_dy1    = met_dy  ;
   leaves.met_dx2    = met_sx  ;
   leaves.met_dy2    = met_sy  ;

   if ( selectedJets.size() > 0 )  return true ; 
   else                            return false ;    

}

vector<double> DPAnalysis::JECUncertainty( double jetpt, double jeteta, JetCorrectionUncertainty* unc ) {

      unc->setJetPt( jetpt);
      unc->setJetEta( jeteta);
      double sup = unc->getUncertainty(true); // up variation
      unc->setJetPt(jetpt);
      unc->setJetEta(jeteta);
      double sdw = unc->getUncertainty(false); // down variation

      double Unc_up  =  sup  ;
      double Unc_dw  =  sdw  ;
      double Unc_max =  max(sup,sdw)  ;

      vector<double> UncV ;
      UncV.push_back( Unc_up ) ;
      UncV.push_back( Unc_dw ) ;
      UncV.push_back( Unc_max ) ;

      //cout<<" method 1 ->  max unc : "<< Unc_max <<" up = "<< Unc_up << " dw = "<< Unc_dw << endl ;   
 
      return UncV ;
}

vector<double> DPAnalysis::JECUncertainty( double jetpt, double jeteta ) {

  const int nsrc = 19;
  const char* srcnames[nsrc] =
  {"Absolute", "HighPtExtra", "SinglePionECAL", "SinglePionHCAL", "Flavor", "Time", 
   "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtEC1", "RelativePtEC2", "RelativePtHF",
   "RelativeStatEC2", "RelativeStatHF","PileUpDataMC", "PileUpBias", "PileUpPtBB", "PileUpPtEC", "PileUpPtHF"};

  std::vector<JetCorrectionUncertainty*> vsrc(nsrc);

  double sum2_max(0), sum2_up(0), sum2_dw(0);
  for (int isrc = 0; isrc < nsrc; isrc++) {
      const char *name = srcnames[isrc];
      JetCorrectorParameters *p = new JetCorrectorParameters("/uscms/home/sckao/work/Exotica/CMSSW_5_3_7_patch4/src/EXO/DPAnalysis/test/GMSB/Fall12_V6_DATA_UncertaintySources_AK5PFchs.txt", name);
      JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
      unc->setJetPt( jetpt);
      unc->setJetEta( jeteta);
      double sup = unc->getUncertainty(true); // up variation
      unc->setJetPt(jetpt);
      unc->setJetEta(jeteta);
      double sdw = unc->getUncertainty(false); // down variation

      sum2_up += pow(sup,2);
      sum2_dw += pow(sdw,2);
      sum2_max += pow(max(sup,sdw),2);

      //delete p ;
      //delete unc ;
  } 
  double Unc_up  = sqrt( sum2_up ) ;
  double Unc_dw  = sqrt( sum2_dw ) ;
  double Unc_max = sqrt( sum2_max ) ;

  //cout<<" method 2 ->  max unc : "<< Unc_max <<" up = "<< Unc_up << " dw = "<< Unc_dw << endl ;   
  //cout<<" "<<endl ;

  vector<double> UncV ;
  UncV.push_back( Unc_up ) ;
  UncV.push_back( Unc_dw ) ;
  UncV.push_back( Unc_max ) ;
    
  return UncV ;

}


//double DPAnalysis::JERUncertainty( double jetpt, double jeteta ) {
void DPAnalysis::JERUncertainty( Handle< std::vector<pat::Jet> > patjets ) {

   //cout<<" ============= "<<endl ;
   double met_x = 0 ;
   double met_y = 0 ;
   double SF = 1. ;
   for (std::vector<pat::Jet>::const_iterator it = patjets->begin(); it != patjets->end(); it++) {
 
       const reco::GenJet* matchedGenJet = it->genJet() ;
       if ( it->pt() < 10 ) continue ;
       if ( matchedGenJet == NULL ) continue ;

       met_x += it->px() ;
       met_y += it->py() ;
       
       // This is data/MC ratio
       if ( fabs(it->eta()) < 0.5 ) SF = 1.052 ;
       if ( fabs(it->eta()) >= 0.5 && fabs(it->eta()) < 1.1 ) SF = 1.057 ;
       if ( fabs(it->eta()) >= 1.1 && fabs(it->eta()) < 1.7 ) SF = 1.096 ;
       if ( fabs(it->eta()) >= 1.7 && fabs(it->eta()) < 2.3 ) SF = 1.134 ;
       if ( fabs(it->eta()) >= 2.3 && fabs(it->eta()) < 2.5 ) SF = 1.288 ;
       double dPt = ( it->pt() - matchedGenJet->pt() )*(1-SF) ;
       double ptscale = max( 0.0, ( it->pt() + dPt)/it->pt() ) ;

       double newPt = it->pt() * ptscale ;
       met_x -= ( it->px() * ptscale ) ;
       met_y -= ( it->py() * ptscale ) ;
       printf(" reco pt: %.1f, gen pt: %.1f  newPt: %.1f \n ", it->pt() , matchedGenJet->pt(), newPt ) ;
       /*
       if ( matchedGenJet != NULL ) printf(" reco pt: %.1f, gen pt: %.1f \n ", it->pt() , matchedGenJet->pt() ) ;

       if ( it->isPFJet() && matchedGenJet == NULL ) { 
          printf(" reco pt: %.1f but no gen jets \n", it->pt() ) ;
	  if ( it->isMuon() ) cout<<" it's muon "<<endl ; 
	  if ( it->isElectron() ) cout<<" it's electron "<<endl ; 
	  if ( it->isPhoton() ) cout<<" it's photon "<<endl ; 
       }
       */
   }
	  
}


bool DPAnalysis::ElectronSelection( Handle<reco::GsfElectronCollection> electrons, 
                                    vector<const reco::GsfElectron*>& selectedElectrons ) {

   // Electron Identification Based on Simple Cuts
   // https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID#Selections_and_How_to_use_them
   int k = 0 ;
   double met_dx(0), met_dy(0) ;
   for(reco::GsfElectronCollection::const_iterator it = electrons->begin(); it != electrons->end(); it++) {

       // calculate met uncertainty
       if ( it->pt() < 10. ) continue ;
       met_dx += it->px() ;
       met_dy += it->py() ;
       double ptscale = ( it->isEB() ) ? 1.006 : 1.015;
       met_dx -= ( it->px() * ptscale ) ;
       met_dy -= ( it->py() * ptscale ) ;


       if ( it->pt() < electronCuts[0] || fabs( it->eta() ) > electronCuts[1] ) continue ;

       // Isolation Cuts
       float ecalSumEt = ( it->isEB() ) ? max(0., it->dr03EcalRecHitSumEt() - 1. ) : it->dr03EcalRecHitSumEt();
       float hcalSumEt = it->dr03HcalTowerSumEt();
       float trkSumPt  = it->dr03TkSumPt();  
       //double relIso   = (ecalSumEt + hcalSumEt + trkSumPt) / it->pt() ;

       // obsoleted
       //if ( relIso > electronCuts[2] &&  it->isEB() ) continue ;
       //if ( relIso > electronCuts[3] && !it->isEB() ) continue ;

       double nLost = it->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() ;
       if ( nLost >= electronCuts[4]  ) continue ;
       if ( k >= MAXELE ) break ;
       selectedElectrons.push_back( &(*it) ) ;
       leaves.elePx[k] = it->p4().Px() ;
       leaves.elePy[k] = it->p4().Py() ;
       leaves.elePz[k] = it->p4().Pz() ;
       leaves.eleE[k]  = it->p4().E() ;
       leaves.eleEcalIso[k] = ecalSumEt ;
       leaves.eleHcalIso[k] = hcalSumEt ;
       leaves.eleTrkIso[k]  = trkSumPt ;
       leaves.eleNLostHits[k]  = nLost ;
       leaves.e_cHadIso[k]  = it->pfIsolationVariables().chargedHadronIso ;
       leaves.e_nHadIso[k]  = it->pfIsolationVariables().neutralHadronIso ;
       leaves.e_photIso[k]  = it->pfIsolationVariables().photonIso ;
       k++;
   }
   leaves.nElectrons = (int)( selectedElectrons.size() ) ;
   leaves.met_dx3  += met_dx ;
   leaves.met_dy3  += met_dy ;

   if ( selectedElectrons.size() > 0 )  return true ; 
   else                                 return false ;    

}

bool DPAnalysis::MuonSelection( Handle<reco::MuonCollection> muons, vector<const reco::Muon*>& selectedMuons ) {

   int k = 0;
   for(reco::MuonCollection::const_iterator it = muons->begin(); it != muons->end(); it++) {
       if ( it->pt() < muonCuts[0] || fabs( it->eta() ) > muonCuts[1] ) continue ;
       // Isolation for PAT muon
       //double relIso =  ( it->chargedHadronIso()+ it->neutralHadronIso() + it->photonIso () ) / it->pt();
       // Isolation for RECO muon
       double relIso = 99. ;
       if ( it->isIsolationValid() ) {
	 relIso = ( it->isolationR05().emEt + it->isolationR05().hadEt + it->isolationR05().sumPt ) / it->pt();
       }
       if ( relIso > muonCuts[2] ) continue ;
       /*
       double dR = 999. ;
       for (size_t j=0; j < selectedJets.size(); j++ ) {
           double dR_ =  ROOT::Math::VectorUtil::DeltaR( it->p4(), selectedJets[j]->p4() ) ; 
           if ( dR_ < dR ) dR = dR_ ;
       }
       if ( dR <= muonCuts[3] ) continue ;
       */
       if ( k >= MAXMU ) break ;
       selectedMuons.push_back( &(*it) ) ;
       leaves.muPx[k] = it->p4().Px() ;
       leaves.muPy[k] = it->p4().Py() ;
       leaves.muPz[k] = it->p4().Pz() ;
       leaves.muE[k]  = it->p4().E() ;
       leaves.muIso[k] = relIso ;
       k++ ;
   }
   leaves.nMuons = (int)( selectedMuons.size() ) ;

   if ( selectedMuons.size() > 0 )  return true ; 
   else                             return false ;    

}


bool DPAnalysis::sMinorSelection( vector<const reco::Photon*>& selectedPhotons,  Handle<EcalRecHitCollection> recHitsEB, 
                                  Handle<EcalRecHitCollection> recHitsEE ) {

    // sMinor and sMajor are from 
    // CMSSW/JetMETCorrections/GammaJet/src/GammaJetAnalyzer.cc
    vector<float> sMinV ;
 
    size_t sz = selectedPhotons.size() ;
    for ( size_t i=0; i < selectedPhotons.size(); i++ ) {

        // S_Minor Cuts from the seed cluster
        reco::CaloClusterPtr SCseed = selectedPhotons[i]->superCluster()->seed() ;
        const EcalRecHitCollection* rechits = ( selectedPhotons[i]->isEB()) ? recHitsEB.product() : recHitsEE.product() ;
        Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*SCseed, *rechits);
        float sMin =  moments.sMin  ;
        //float sMaj =  moments.sMaj  ;

        // seed Time 
        /* 
        pair<DetId, float> maxRH = EcalClusterTools::getMaximum( *SCseed, rechits );
        DetId seedCrystalId = maxRH.first;
        EcalRecHitCollection::const_iterator seedRH = rechits->find(seedCrystalId);
        float seedTime = (float)seedRH->time();
        */

        //if ( sMin < 0.  ) selectedPhotons.erase( selectedPhotons.begin() + i ) ;
        if ( sMin <= photonCuts[3] || sMin >= photonCuts[4] ) selectedPhotons.erase( selectedPhotons.begin() + i ) ;
        sMinV.push_back( sMin );
    }
    if ( sMinV.size() > 0 ) sMin_ = sMinV[0] ; 
    
    if ( sz != selectedPhotons.size() ) return true ;
    else                                return false ;
}



bool DPAnalysis::IsoPhotonSelection( vector<const reco::Photon*>& selectedPhotons ) {

    // Another photon Isolation also can be done by using the EgammaIsolationAlogs
    // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoEgamma/EgammaIsolationAlgos/interface/

    size_t sz = selectedPhotons.size() ;
    for ( size_t i=0; i < selectedPhotons.size(); i++ ) {

        if ( fabs( selectedPhotons[i]->eta() ) > 1.3 ) {
           selectedPhotons.erase( selectedPhotons.begin() + i ) ;
           continue ;
        }
        // Isolation Cuts 
        float ecalSumEt = selectedPhotons[i]->ecalRecHitSumEtConeDR04();
	float hcalSumEt = selectedPhotons[i]->hcalTowerSumEtConeDR04();
	float trkSumPt  = selectedPhotons[i]->trkSumPtSolidConeDR04();  
	bool trkIso  = ( ( trkSumPt / selectedPhotons[i]->pt())     < photonIso[0] ) ; 
	bool ecalIso = ( (ecalSumEt / selectedPhotons[i]->energy()) < photonIso[2] && ecalSumEt < photonIso[1] ) ; 
	bool hcalIso = ( (hcalSumEt / selectedPhotons[i]->energy()) < photonIso[4] && hcalSumEt < photonIso[3] ) ; 
	if ( !trkIso || !ecalIso || !hcalIso ) selectedPhotons.erase( selectedPhotons.begin() + i ) ;
    }

    if ( sz != selectedPhotons.size() ) return true ;
    else                                return false ;

}


bool DPAnalysis::GammaJetVeto( vector<const reco::Photon*>& selectedPhotons, vector<const reco::PFJet*>& selectedJets) {

     bool isGammaJets = false ;
     
     if (  selectedJets.size() > 0 && selectedPhotons.size() > 0  ) {
       double dR      = ROOT::Math::VectorUtil::DeltaR( selectedJets[0]->p4(), selectedPhotons[0]->p4() ) ;
       double PtRatio = selectedJets[0]->pt() / selectedPhotons[0]->pt() ;
       if ( dR > (2.*3.1416/3.) && PtRatio > 0.7 && PtRatio < 1.3 )  isGammaJets = true ;
     }
     return isGammaJets ;
}

// return true if it's a conversion case
bool DPAnalysis::ConversionVeto( const reco::Photon* thePhoton  ) {

   bool passVeto = ConversionTools::hasMatchedPromptElectron( thePhoton->superCluster(), electrons, hConversions, beamspot->position() );
   return passVeto ;
}

// type : 1: chargedHadron 2:neutralHadron 3:Photon
double DPAnalysis::RhoCorrection( int type , double eta ) {

     double EA = 0 ;
     if ( type == 1 ) {
        if (                      fabs(eta) < 1.0   ) EA = 0.012 ;
        if ( fabs(eta) > 1.0   && fabs(eta) < 1.479 ) EA = 0.010 ;
        if ( fabs(eta) > 1.479 && fabs(eta) < 2.0   ) EA = 0.014 ;
        if ( fabs(eta) > 2.0   && fabs(eta) < 2.2   ) EA = 0.012 ;
        if ( fabs(eta) > 2.2   && fabs(eta) < 2.3   ) EA = 0.016 ;
        if ( fabs(eta) > 2.3   && fabs(eta) < 2.4   ) EA = 0.020 ;
        if ( fabs(eta) > 2.4                        ) EA = 0.012 ;
     }
     if ( type == 2 ) {
        if (                      fabs(eta) < 1.0   ) EA = 0.030 ;
        if ( fabs(eta) > 1.0   && fabs(eta) < 1.479 ) EA = 0.057 ;
        if ( fabs(eta) > 1.479 && fabs(eta) < 2.0   ) EA = 0.039 ;
        if ( fabs(eta) > 2.0   && fabs(eta) < 2.2   ) EA = 0.015 ;
        if ( fabs(eta) > 2.2   && fabs(eta) < 2.3   ) EA = 0.024 ;
        if ( fabs(eta) > 2.3   && fabs(eta) < 2.4   ) EA = 0.039 ;
        if ( fabs(eta) > 2.4                        ) EA = 0.072 ;
     }
     if ( type == 3 ) {
        if (                      fabs(eta) < 1.0   ) EA = 0.148 ;
        if ( fabs(eta) > 1.0   && fabs(eta) < 1.479 ) EA = 0.130 ;
        if ( fabs(eta) > 1.479 && fabs(eta) < 2.0   ) EA = 0.112 ;
        if ( fabs(eta) > 2.0   && fabs(eta) < 2.2   ) EA = 0.216 ;
        if ( fabs(eta) > 2.2   && fabs(eta) < 2.3   ) EA = 0.262 ;
        if ( fabs(eta) > 2.3   && fabs(eta) < 2.4   ) EA = 0.260 ;
        if ( fabs(eta) > 2.4                        ) EA = 0.266 ;
     }
     return EA*rhoIso ;

}

void DPAnalysis::PhotonPFIso( std::vector<const reco::Photon*> thePhotons, const reco::PFCandidateCollection* pfParticlesColl, reco::VertexRef vtxRef, Handle< reco::VertexCollection > vtxColl ) {

    for ( size_t k=0; k< thePhotons.size(); k++ ) { 
        isolator.fGetIsolation( thePhotons[k], pfParticlesColl, vtxRef, vtxColl );
        //cout<<"PF  :  "<<isolator.getIsolationCharged()<<" : "<<isolator.getIsolationPhoton()<<" : "<<isolator.getIsolationNeutral()<<endl;
        leaves.cHadIso[k] = max( isolator.getIsolationCharged() - RhoCorrection( 1, thePhotons[k]->eta() ), 0.) ;
	leaves.nHadIso[k] = max( isolator.getIsolationNeutral() - RhoCorrection( 2, thePhotons[k]->eta() ), 0.) ;
	leaves.photIso[k] = max( isolator.getIsolationPhoton()  - RhoCorrection( 3, thePhotons[k]->eta() ), 0.) ;
        //printf(" cHad: %.3f, nHad: %.3f, phot: %.3f \n", leaves.cHadIso[k],  leaves.nHadIso[k],  leaves.photIso[k] ) ;
    }
}

/*
void DPAnalysis::VertexFinder( Handle<reco::TrackCollection> tracks ) {

    // prepare a histogram vector with 101 elements, bin size =  3cm -> -151.5 ~ 151.5 cm
    int hTrk[101] = { 0 } ;
    for (reco::TrackCollection::const_iterator it = tracks->begin(); it != tracks->end(); it++ )  {
        if ( fabs(it->d0()) > 2. ) continue ;
         TrkInfo tf ;
         tf.dz  = it->dz() ;
         tf.dsz = it->dsz() ;
         tf.d0  = it->d0() ;
         tf.pt  = it->pt() ;
         tf.vz  = it->vz() ;
         tf.vr  = sqrt( (it->vx()*it->vx()) + (it->vy()*it->vy()) ) ;

         int ibin = (int)( ( it->dz() + 151.5 ) / 3.)  ; 
         hTrk[ ibin ]++ ;
    }
}
*/


//define this as a plug-in
DEFINE_FWK_MODULE(DPAnalysis);
