 // -*- C++ -*-
//
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
// $Id: DPAnalysis.cc,v 1.18 2012/09/24 22:43:11 sckao Exp $
//
//


// system include files
#include "DPAnalysis.h"
#include "Ntuple.h"

using namespace cms ;
using namespace edm ;
using namespace std ;

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
   jetSource            = iConfig.getParameter<edm::InputTag> ("jetSource");
   trackSource          = iConfig.getParameter<edm::InputTag> ("trackSource");

   EBRecHitCollection   = iConfig.getParameter<edm::InputTag> ("EBRecHitCollection") ;
   EERecHitCollection   = iConfig.getParameter<edm::InputTag> ("EERecHitCollection") ;
   CSCSegmentTag        = iConfig.getParameter<edm::InputTag> ("CSCSegmentCollection") ;
   cscHaloTag           = iConfig.getParameter<edm::InputTag> ("cscHaloData");
   //pileupSource         = iConfig.getParameter<edm::InputTag>("addPileupInfo");

   hbhereco             = iConfig.getParameter<edm::InputTag> ("hbhereco") ; // HE Hallo Tag
   vtxCuts              = iConfig.getParameter<std::vector<double> >("vtxCuts");
   jetCuts              = iConfig.getParameter<std::vector<double> >("jetCuts");
   metCuts              = iConfig.getParameter<std::vector<double> >("metCuts");
   photonCuts           = iConfig.getParameter<std::vector<double> >("photonCuts");
   photonIso            = iConfig.getParameter<std::vector<double> >("photonIso");
   electronCuts         = iConfig.getParameter<std::vector<double> >("electronCuts");
   muonCuts             = iConfig.getParameter<std::vector<double> >("muonCuts");  
   //triggerPatent        = iConfig.getUntrackedParameter<string> ("triggerName");
   triggerPatent        = iConfig.getParameter< std::vector<string> >("triggerName");
   isData               = iConfig.getParameter<bool> ("isData");
   L1Select             = iConfig.getParameter<bool> ("L1Select");

   const InputTag TrigEvtTag("hltTriggerSummaryAOD","","HLT");
   trigEvent            = iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag", TrigEvtTag);

   gen = new GenStudy( iConfig );

   theFile  = new TFile( rootFileName.c_str(), "RECREATE") ;
   theFile->cd () ;
   theTree  = new TTree ( "DPAnalysis","DPAnalysis" ) ;
   setBranches( theTree, leaves ) ;

   targetTrig = 0 ;
   firedTrig.clear() ;
   for ( size_t i=0; i< triggerPatent.size(); i++ ) firedTrig.push_back(-1) ;  

   // reset the counter
   for ( int i=0; i< 10 ; i++) counter[i] = 0 ;

   // initialize the time corrector
   theTimeCorrector_.initEB("EB");
   theTimeCorrector_.initEE("EElow");
   runID_ = 0 ;

   debugT = false ;
}


DPAnalysis::~DPAnalysis()
{
   // do anything here that needs to be done at desctruction time

   delete gen ;
   cout<<"All:"<< counter[0]<<" Trigger:"<<counter[1]<<" Vertex:"<< counter[2] <<" Photon:"<<counter[3] ;
   cout<<" beamHalo:"<< counter[4] <<" Jet:"<< counter[5] <<" MET:"<<counter[6] <<" Pre-Selection:"<<counter[7] <<endl ;

   theFile->cd () ;
   theTree->Write() ; 
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
   iSetup.get<CaloGeometryRecord> ().get (pGeometry) ;
   theGeometry = pGeometry.product() ;
   // event time
   eventTime = iEvent.time() ;

   initializeBranches( theTree, leaves );

   leaves.bx          = iEvent.bunchCrossing();
   leaves.lumiSection = iEvent.id().luminosityBlock();
   leaves.orbit       = iEvent.orbitNumber();
   leaves.runId       = iEvent.id ().run () ;
   leaves.eventId     = iEvent.id ().event () ;

   /* 
   Handle<std::vector< PileupSummaryInfo > >  PupInfo;
   iEvent.getByLabel(pileupSource, PupInfo);
   for( std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
       std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
   }
   */

   if (counter[0] == 0 )  PrintTriggers( iEvent ) ;

   int run_id    = iEvent.id().run()  ;

   counter[0]++ ;  
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

   if ( passTrigger ) counter[1]++ ;  

   // get the generator information
   if ( !isData ) { 
      gen->GetGenEvent( iEvent, leaves );
      //gen->GetGen( iEvent, leaves );
   }
   //if ( !isData ) gen->PrintGenEvent( iEvent );

   bool pass = EventSelection( iEvent, iSetup ) ;
   if ( pass && passTrigger ) counter[7]++ ;   

  // // for cscsegments and halo muon/photon studies

 //  Handle<CSCSegmentCollection>       cscSegments ;
  //   iEvent.getByLabel( CSCSegmentTag,  cscSegments );
 //  bool haloMuon = BeamHaloMatch( cscSegments, selectedPhotons, iSetup ) ;
   
 //  if ( haloMuon ) cout<<" haloMuon ! "<<endl ;
 

  
   // fill the ntuple
   if ( pass && !isData ) theTree->Fill();
   if ( pass && isData && passTrigger ) theTree->Fill();
}

bool DPAnalysis::EventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

   Handle<reco::BeamHaloSummary>       beamHaloSummary ;
   Handle<edm::TriggerResults>         triggers;
   Handle<reco::VertexCollection>      recVtxs;
   Handle<reco::PhotonCollection>      photons; 
   Handle<reco::GsfElectronCollection> electrons; 
   Handle<reco::MuonCollection>        muons; 
   Handle<reco::PFJetCollection>       jets; 
   Handle<reco::PFMETCollection>       met; 
   Handle<EcalRecHitCollection>        recHitsEB ;
   Handle<EcalRecHitCollection>        recHitsEE ;
   Handle<reco::TrackCollection>       tracks; 

   iEvent.getByLabel( trigSource,     triggers );
   iEvent.getByLabel( pvSource,       recVtxs  );
   iEvent.getByLabel( photonSource,   photons  );
   iEvent.getByLabel( electronSource, electrons);
   iEvent.getByLabel( muonSource,     muons );
   iEvent.getByLabel( jetSource,      jets  );
   iEvent.getByLabel( metSource,      met  );
   iEvent.getByLabel( EBRecHitCollection,     recHitsEB );
   iEvent.getByLabel( EERecHitCollection,     recHitsEE );
   iEvent.getByLabel("BeamHaloSummary", beamHaloSummary) ;
   iEvent.getByLabel( trackSource,    tracks  );

   bool passEvent = true ;

   // find trigger matched objects
   //cout<<" ~~~~~~~~~~~~~~~~~ "<<endl ;
   const reco::Photon rPho ;
   GetTrgMatchObject(  rPho , iEvent,  photonSource ) ;
   //cout<<" ----------------- "<<endl ;
   const reco::PFMET pfMet_ = (*met)[0] ;
   GetTrgMatchObject( pfMet_, iEvent,  metSource ) ;
   //cout<<" ================= "<<endl ;
   //cout<<" "<<endl ;
  

   bool hasGoodVtx = VertexSelection( recVtxs );
   if ( !hasGoodVtx ) passEvent = false ;
   if ( passEvent )   counter[2]++ ;  

   selectedPhotons.clear() ;
   PhotonSelection( photons, recHitsEB, recHitsEE, tracks, selectedPhotons ) ;
   if ( selectedPhotons.size() < (size_t)photonCuts[6] )  passEvent = false ;
   if ( passEvent )   counter[3]++ ;  

   if( beamHaloSummary.isValid() ) {
     const reco::BeamHaloSummary TheSummary = (*beamHaloSummary.product() );
    // Tag CSCTightHaloEvents
 if(TheSummary.EcalTightHaloId() && TheSummary.HcalTightHaloId() && TheSummary.CSCTightHaloId() && TheSummary.GlobalTightHaloId()) leaves.IsBeamHaloIDTightTag= true;

//   CSCTighthaloId Tagging
 if(TheSummary.CSCTightHaloId()) leaves.IscscHaloTight_Tag = true;
 
     if( !TheSummary.CSCTightHaloId() && passEvent ) { 
       counter[4]++ ;  
     } else {
       passEvent = false ;
     }
   } else {
       counter[4]++ ;
   }
   CSCHaloCleaning( iEvent, selectedPhotons ) ;

    // for cscsegments and halo muon/photon studies
    Handle<CSCSegmentCollection>       cscSegments ;
    iEvent.getByLabel( CSCSegmentTag,  cscSegments );
    bool haloMuon = BeamHaloMatch( cscSegments, selectedPhotons, iSetup ) ;
   
   if ( haloMuon ) cout<<" haloMuon ! "<<endl ;

   // HE Halo Matching studies!

    // Grab HBHE Rechits
   //  Handle<HBHERecHitCollection> hcalRecHitHandle;
   //  const HBHERecHitCollection *hbhe =  (hcalRecHitHandle.product()); 
   //
     edm::Handle<edm::SortedCollection<HBHERecHit> > hcalRecHitHandle;
     iEvent.getByLabel("hbhereco", hcalRecHitHandle);                       
   //  const edm::SortedCollection<HBHERecHit> *hbhe =  (hcalRecHitHandle.product()); 
     bool  hehalo = HEBeamHaloMatch( hcalRecHitHandle , selectedPhotons, iSetup );
      if(hehalo) cout <<" HE Halo Muon  ! " << endl;
   
   //IsoPhotonSelection( selectedPhotons ) ;
   //if ( selectedPhotons.size() < photonCuts[5] )  passEvent = false ;

   selectedJets.clear() ;
   JetSelection( jets, selectedPhotons, selectedJets );
   //bool isGammaJets = GammaJetVeto( selectedPhotons, selectedJets ) ;
   //if ( isGammaJets ) passEvent = false ;
   if ( selectedJets.size() < jetCuts[3] )   passEvent = false ;
   if ( passEvent )   counter[5]++ ;   
   
   selectedElectrons.clear() ;
   ElectronSelection( electrons, selectedElectrons ) ;

   selectedMuons.clear() ;
   MuonSelection( muons, selectedMuons );

   //HLTMET( jets, selectedMuons );   

   const reco::PFMET pfMet = (*met)[0] ;
   leaves.met   = pfMet.et() ;
   leaves.metPx = pfMet.px() ;
   leaves.metPy = pfMet.py() ;
   if ( pfMet.pt() < metCuts[0]  ) passEvent = false;
   if ( passEvent )   counter[6]++ ;  

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

// HE Halo Cleaning right here!
//
//
bool DPAnalysis::HEBeamHaloMatch(Handle<HBHERecHitCollection> hhit, vector<const reco::Photon*> & selectedPhotons, const EventSetup &iSetup){

  // Grab rechit geometry(position)
   ESHandle<CaloGeometry> geoHandle;
   iSetup.get<CaloGeometryRecord>().get(geoHandle);
   const CaloGeometry& geometry = *geoHandle;
  
 // //Get & use rechit collection for removing  hot crystals
 // Handle<EcalRecHitCollection> EBrecHitColl;
 // iEvent.getByLabel("ecalRecHit","EcalRecHitsEB",EBrecHitColl);
 // const EcalRecHitCollection *rechit = EBrecHitColl.product();
 

    bool Hehalomatch = false ;
    // Loop over selected Photons
    for ( size_t i=0; i< selectedPhotons.size() ; i++ ) {
           
      // Important Variable to keep!
      float heGdphi = 99.;
      float heRho  = -1 ;
      float  hetime = 99999.;
      float heEnergy = 999999.;
      float heXYradius = 9999.;
      float   hedphi  = 99.;
                                                                              
     // Get Selected Photons energy and Phi
  // double scPhoEnergy = selectedPhotons[i]->superCluster()->energy();
   double scPhoPhi = selectedPhotons[i]->superCluster()->phi();

    /////////////////////////////////////////////
    ///// Loop HE Rechits to grab info! and do Matching ////////
   ////////////////////////////////////////////
   //Loop over HBHE rechit collection, extract he rechit information
   for (HBHERecHitCollection::const_iterator hhIt = hhit->begin(); hhIt != hhit->end(); hhIt++){
      HcalDetId id(hhIt->detid());
      //returns the subdetector, either HB(==1) or HE(==2)
      int hbheDetId = id.subdet();

      //get geometry info of HBHE rechits
      const CaloCellGeometry *hbhe_cell = (geometry).getGeometry(hhIt->id());
      GlobalPoint hbhePosition = hbhe_cell->getPosition();

      //only looking in HE members
      if (hbheDetId == 2){
       float herhEnergy = hhIt->energy();
       float herhTime = hhIt->time();
       float herhPhi = hbhePosition.phi();
       float radius = hbhePosition.perp(); // sqrt(x*x + y*y)
       float mag  = hbhePosition.mag();   // sqrt(x*x + y*y + z*z)

       LorentzVector HEp4(hbhePosition.x(), hbhePosition.y(), hbhePosition.z(), hbhePosition.mag() ) ; 

       //compare phi of highest Et supercluster in event with phi's of all he rechits
      float  dPhi = ROOT::Math::VectorUtil::DeltaPhi( HEp4,   selectedPhotons[i]->p4() );
      float  AbsdPhi = fabs(scPhoPhi-herhPhi);   // Delta Phi between Selected Photon and HE rechit Position
      //phi wrapping
      if (AbsdPhi > TMath::Pi()){AbsdPhi = TMath::Pi()*2 - AbsdPhi;}
      if (dPhi < 0.05 ) Hehalomatch = true;

           heGdphi = dPhi ;
//        if( fabs(dPhi) < heGdphi )  {
           heXYradius = radius;
           hetime = herhTime;
           heRho  =  mag ;
           heEnergy = herhEnergy;
           hedphi  = AbsdPhi;

//              }
            }//end of HE loop
           }//end of HBHE rechit loop
   
        // Fill leaves
         leaves.HERho[i] = heRho;
         leaves.HETime[i] = hetime;
         leaves.HEGPhi[i] = heGdphi;
         leaves.HERadius[i] = heXYradius;
         leaves.HEEnergy[i] = heEnergy;
         leaves.HEdphi[i]   = hedphi;
 
    } // End loop photons

    return Hehalomatch;
}

void DPAnalysis::CSCHaloCleaning( const edm::Event& iEvent, vector<const reco::Photon*>& selectedPhotons ) { 

   Handle<reco::CSCHaloData> cschalo;
   iEvent.getByLabel( cscHaloTag, cschalo ); 

   if ( cschalo.isValid() ) {
      const reco::CSCHaloData cscData = *(cschalo.product()) ;
      int nOutTimeHits       = cscData.NumberOfOutTimeHits() ;
      int nMinusHTrks        = cscData.NHaloTracks( reco::HaloData::minus ) ; 
      int nPlusHTrks         = cscData.NHaloTracks( reco::HaloData::plus )  ; 
      //int nTrksSmallBeta     = cscData.NTracksSmallBeta(); 
     // int nHaloSegs          = cscData.NFlatHaloSegments();
       
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

      std::vector<GlobalPoint> gp = cscData.GetCSCTrackImpactPositions() ;
      //cout<<" impact gp sz : "<< gp.size() <<" photon sz:"<< selectedPhotons.size() << endl ;

      for (vector<GlobalPoint>::const_iterator it = gp.begin(); it != gp.end() ; ++it ) {
          double rho = sqrt(  (it->x()*it->x()) + (it->y()*it->y()) );
          //cout<<"    ==>  phi:"<<  it->phi() <<" rho: "<< rho << " z: "<< it->z() <<endl ;
          leaves.haloPhi = it->phi() ;
          leaves.haloRho = rho;
           
          
          for ( size_t i=0; i< selectedPhotons.size() ; i++ ) { 
              double dPhi = fabs( selectedPhotons[i]->phi() - it->phi() ) ;
              double rhoG = sqrt(  ( selectedPhotons[i]->p4().x()*selectedPhotons[i]->p4().x() ) 
                                 + ( selectedPhotons[i]->p4().y()*selectedPhotons[i]->p4().y() )  );
              double dRho = fabs( rho - rhoG ) ;
              cout<<"    ("<< i << ") ==> dPhi : "<< dPhi <<" dRho : "<< dRho << endl ;

           // My  CSC Trk Matching condition
           if( dPhi < 0.18 ) leaves.IscscHaloTrk_Tag = true;
          }
         
      }
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
                  cout<<" Trigger Found ("<<j <<"):  "<<tName ;
                  cout<<" Idx: "<< i <<" triggers "<<endl;
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

bool DPAnalysis::VertexSelection( Handle<reco::VertexCollection> vtx ) {

    int thisVertex=0;
    bool hasGoodVertex = true ;
    int totalN_vtx = 0 ;

    for(reco::VertexCollection::const_iterator v=vtx->begin();  v!=vtx->end() ; v++){

       if ( ! v->isValid() ||  v->isFake() ) continue ;
       if ( fabs(v->z()) >= vtxCuts[0] ) continue ; 
       if (   v->ndof()   < vtxCuts[1] ) continue ;
       double d0 = sqrt( ( v->x()*v->x() ) + ( v->y()*v->y() ) );
       if ( d0 >= vtxCuts[2] ) continue ;
       // counting real number of vertices
       totalN_vtx++ ;

       if ( thisVertex >= MAXVTX ) continue ;
       leaves.vtxNTracks[thisVertex]= v->tracksSize();
       leaves.vtxChi2[thisVertex] =   v->chi2();
       leaves.vtxNdof[thisVertex] =   v->ndof();
       leaves.vtxX[thisVertex] =      v->x();
       leaves.vtxY[thisVertex] =      v->y();
       leaves.vtxZ[thisVertex] =      v->z();
       leaves.vtxDx[thisVertex] =     v->xError();
       leaves.vtxDy[thisVertex] =     v->yError();
       leaves.vtxDz[thisVertex] =     v->zError();
       
       thisVertex++ ;
     }
     leaves.nVertices = thisVertex ;
     leaves.totalNVtx = totalN_vtx ;
 
     if ( thisVertex < 1 )   hasGoodVertex = false ;
     return hasGoodVertex ;
}

bool DPAnalysis::PhotonSelection( Handle<reco::PhotonCollection> photons, Handle<EcalRecHitCollection> recHitsEB, Handle<EcalRecHitCollection> recHitsEE, Handle<reco::TrackCollection> tracks, vector<const reco::Photon*>& selectedPhotons ) {

   int k= 0 ;
   double maxPt = 0 ;
   for(reco::PhotonCollection::const_iterator it = photons->begin(); it != photons->end(); it++) {

       // fiducial cuts
       if ( k >= MAXPHO ) break ;
       if ( it->pt() < photonCuts[0] || fabs( it->eta() ) > photonCuts[1] ) continue ;
       //float hcalIsoRatio = it->hcalTowerSumEtConeDR04() / it->pt() ;
       //if  ( ( hcalIsoRatio + it->hadronicOverEm() )*it->energy() >= 6.0 ) continue ;
       // pixel veto
       if ( it->hasPixelSeed() ) continue ;
 
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

       // sMin and sMaj cuts
       if ( sMaj  > photonCuts[2] ) continue ;
       if ( sMin <= photonCuts[3] || sMin >= photonCuts[4] ) continue ;

       // Isolation Cuts 
       float ecalSumEt = it->ecalRecHitSumEtConeDR04();
       float hcalSumEt = it->hcalTowerSumEtConeDR04();
       float trkSumPt  = it->trkSumPtSolidConeDR04();  

       bool trkIso  = ( ( trkSumPt / it->pt())     < photonIso[0] ) ; 
       bool ecalIso = ( (ecalSumEt / it->energy()) < photonIso[2] && ecalSumEt < photonIso[1] ) ; 
       bool hcalIso = ( (hcalSumEt / it->energy()) < photonIso[4] && hcalSumEt < photonIso[3] ) ; 
       if ( !trkIso || !ecalIso || !hcalIso ) continue ;

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
       phoTmp.fSpike = -1 ;
       phoTmp.maxSX  = -1 ;
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
       leaves.fSpike[k]       = phoTmp.fSpike ;
       leaves.maxSwissX[k]    = phoTmp.maxSX ; 
       leaves.seedSwissX[k]   = swissX ;
 
       debugT = false ;
       // refine the timing to exclude out-of-time xtals
       ClusterTime( it->superCluster(), recHitsEB , recHitsEE, phoTmp );
       //cout<<" 3rd xT : "<< aveXtalTime <<"  xTE : "<< aveXtalTimeErr << endl;

       leaves.phoPx[k] = it->p4().Px() ;
       leaves.phoPy[k] = it->p4().Py() ;
       leaves.phoPz[k] = it->p4().Pz() ;
       leaves.phoE[k]  = it->p4().E() ;
       leaves.phoHoverE[k]  = it->hadronicOverEm() ;
       leaves.phoEcalIso[k] = ecalSumEt ;
       leaves.phoHcalIso[k] = hcalSumEt ;
       leaves.phoTrkIso[k]  = trkSumPt ;
       leaves.dR_TrkPho[k]  = minDR ;
       leaves.pt_TrkPho[k]  = trkPt ;
       leaves.sMinPho[k]    = sMin ;
       leaves.sMajPho[k]    = sMaj ;

       leaves.seedTime[k]     = seedTime ;
       leaves.seedTimeErr[k]  = seedTimeErr ;
       leaves.aveTime[k]      = phoTmp.t ;       // weighted ave. time of seed cluster 
       leaves.aveTimeErr[k]   = phoTmp.dt ;
       leaves.sigmaEta[k]     = it->sigmaEtaEta() ;
       leaves.sigmaIeta[k]    = it->sigmaIetaIeta() ;
       leaves.SigmaIetaIeta[k]    = it->sigmaIetaIeta() ;  // Careless Repetition
       leaves.SigmaEtaEta[k]     = it->sigmaEtaEta() ;
   
       selectedPhotons.push_back( &(*it) ) ;
       k++ ;
   }
   leaves.nPhotons = k ;
   //leaves.nPhotons = (int)( selectedPhotons.size() ) ;

   if ( selectedPhotons.size() > 0 && maxPt >= photonCuts[7] )  return true ; 
   else                               return false ;    

}

// Trial BeamHalo Matching.
bool DPAnalysis::BeamHaloMatch( Handle<CSCSegmentCollection> cscSeg, vector<const reco::Photon*>& selectedPhotons, const EventSetup& iSetup ) {

//Get Geometry from Event
   edm::ESHandle<CSCGeometry> cscGeom;
   iSetup.get<MuonGeometryRecord>().get(cscGeom);
// iSetup.get<CSCRecoGeometryRcd>().get(cscGeom);
if(!cscGeom.isValid()){ std::cout <<"Unable to find CSCMuonGeometryRecord in event!"<< std::endl;}

   bool halomatch = false ;
   for ( size_t i=0; i< selectedPhotons.size() ; i++ ) {

 //      if ( leaves.seedTime[i] > -3. ) continue ;   // Look @ All Photons
       double dPhi = 99. ;
       double cscT = 99. ;
       double dR   = 99. ;
       double cscEta = 99. ;
       double dphi = 99. ;
       for (CSCSegmentCollection::const_iterator it = cscSeg->begin(); it != cscSeg->end(); it++) {
           if ( !it->isValid() ) continue ;
           CSCDetId DetId = it->cscDetId();
	   const CSCChamber* cscchamber = cscGeom->chamber( DetId );
//	   const CSCChamber* cscchamber = TheCSCGeometry->chamber( DetId );
           
	   GlobalPoint  gp = cscchamber->toGlobal( it->localPosition()  );
           double gpMag = sqrt( (gp.x()*gp.x()) + (gp.y()*gp.y()) + (gp.z()*gp.z()) ) ;
	   LorentzVector segP4( gp.x(), gp.y(), gp.z(), gpMag ) ;
           
           dR = sqrt( (gp.x()*gp.x()) + (gp.y()*gp.y()) );

           if ( fabs( gp.eta() ) < 1.6  ) continue ;  // Extend to EB to include CSC Segments from EB
           
           cscEta = gp.eta();
           //double dEta_ = fabs( gp.eta() - selectedPhotons[i]->eta() )  ;
           //dEta = ( dEta_ < dEta ) ? dEta_ : dEta ;
           //double dPhi_ = fabs( gp.phi() - selectedPhotons[i]->phi() )  ;
           double dPhi_ = ROOT::Math::VectorUtil::DeltaPhi( segP4, selectedPhotons[i]->p4() ) ; // Matching CSCSegment to Selected Photon Here!
           dphi = dPhi_;
           if ( fabs(dPhi_) < dPhi ) {
              dPhi = fabs( dPhi_ ) ;
              cscT = it->time() ;
           }
       }
       if ( dPhi < 0.0717 )  { halomatch = true; leaves.IscscHaloSeg_Tag = true; }
          leaves.cscTime[i]  = cscT ;
          leaves.cscdPhi[i]  = dphi ;   // Keep the real value not just the absolute value. 
          leaves.cscR[i]          = dR ;
          leaves.cscEta[i]        = cscEta ;      
    
   }
   return halomatch ;
}


// return time, timeError
pair<double,double> DPAnalysis::ClusterTime( reco::SuperClusterRef scRef, Handle<EcalRecHitCollection> recHitsEB, Handle<EcalRecHitCollection> recHitsEE ) {

  const EcalIntercalibConstantMap& icalMap = ical->getMap();
  float adcToGeV = float(agc->getEBValue());

  double xtime = 0 ;
  double xtimeErr = 0 ;

  // 1. loop all the basic clusters 
  for ( reco::CaloCluster_iterator  clus = scRef->clustersBegin() ;  clus != scRef->clustersEnd();  ++clus) {

      // only use seed basic cluster  
      if ( *clus != scRef->seed() ) continue ;
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
	   
             // SIC Feb 14 2011 -- Add check on RecHit flags (takes care of spike cleaning in 42X)
             if ( !( myhit.checkFlag(EcalRecHit::kGood) || myhit.checkFlag(EcalRecHit::kOutOfTime) || 
                    myhit.checkFlag(EcalRecHit::kPoorCalib)  ) )  continue;

             // swiss cross cleaning 
             //float swissX = (isEB) ? EcalTools::swissCross(detitr->first, *recHitsEB , 0., true ) : 
             //                        EcalTools::swissCross(detitr->first, *recHitsEE , 0., true ) ;

             //if ( swissX > 0.95 ) { 
             //if ( myhit.checkFlag(EcalRecHit::kWeird) || myhit.checkFlag(EcalRecHit::kDiWeird) ) {
                //cout<<" swissX = "<< swissX <<" @ "<< nXtl <<endl ;
                //continue ;
             //}
             nXtl++ ;

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

	     // discard rechits with A/sigma < 12
	     if ( thisamp/(icalconst*lasercalib*adcToGeV) < (1.1*12) ) continue;

	     GlobalPoint pos = theGeometry->getPosition((myhit).detid());

             // time and time correction
	     double thistime = myhit.time();
	     thistime += theTimeCorrector_.getCorrection((float) thisamp/(icalconst*lasercalib*adcToGeV), pos.eta()  );

             // get time error 
             double xtimeErr_ = ( myhit.isTimeErrorValid() ) ?  myhit.timeError() : 999999 ;
 
             xtime     += thistime / pow( xtimeErr_ , 2 ) ;
             xtimeErr  += 1/ pow( xtimeErr_ , 2 ) ;
      }
      //cout<<" total Xtl = " << nXtl << endl ;
  }
  double wAveTime = xtime / xtimeErr ;
  double wAveTimeErr = 1. / sqrt( xtimeErr) ;
  pair<double, double> wAveTE( wAveTime, wAveTimeErr ) ;
  return wAveTE ;  

}

// re-calculate time and timeError as well as normalized chi2
//void DPAnalysis::ClusterTime( reco::SuperClusterRef scRef, Handle<EcalRecHitCollection> recHitsEB, Handle<EcalRecHitCollection> recHitsEE, double& aveTime, double& aveTimeErr, double& nChi2, bool useAllClusters ) {

void DPAnalysis::ClusterTime( reco::SuperClusterRef scRef, Handle<EcalRecHitCollection> recHitsEB, Handle<EcalRecHitCollection> recHitsEE, PhoInfo& phoTmp, bool useAllClusters ) {

  const EcalIntercalibConstantMap& icalMap = ical->getMap();
  float adcToGeV = float(agc->getEBValue());

  double xtime    = 0 ;
  double xtimeErr = 0 ;
  double chi2_bc  = 0 ;
  double ndof     = 0 ;
  double maxSwissX = 0 ;
  int    nBC      = 0 ;
  int    nXtl     = 0 ;
  int    nSpike   = 0 ; 
  int    nSeedXtl = 0 ;
  for ( reco::CaloCluster_iterator  clus = scRef->clustersBegin() ;  clus != scRef->clustersEnd();  ++clus) {

      nBC++ ;
      // only use seed basic cluster  
      bool isSeed = ( *clus == scRef->seed() ) ;
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
             bool gotSpike = ( myhit.checkFlag(EcalRecHit::kWeird) || myhit.checkFlag(EcalRecHit::kDiWeird) )  ;

             // swiss cross cleaning 
             float swissX = (isEB) ? EcalTools::swissCross(detitr->first, *recHitsEB , 0., true ) : 
                                     EcalTools::swissCross(detitr->first, *recHitsEE , 0., true ) ;
             maxSwissX = ( isSeed && swissX  > maxSwissX ) ? swissX : maxSwissX ;
             if ( gotSpike && isSeed ) nSpike++  ;
             if ( isSeed             ) nSeedXtl++  ;
             //if ( gotSpike ) continue ;

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

	     // discard rechits with A/sigma < 12
	     if ( thisamp/(icalconst*lasercalib*adcToGeV) < (1.1*12) ) continue;

	     GlobalPoint pos = theGeometry->getPosition((myhit).detid());

             // time and time correction
	     double thistime = myhit.time();
	     thistime += theTimeCorrector_.getCorrection((float) thisamp/(icalconst*lasercalib*adcToGeV), pos.eta()  );

             // get time error 
             double xtimeErr_ = ( myhit.isTimeErrorValid() ) ?  myhit.timeError() : 999999 ;

             // calculate chi2 for the BC of the seed
             double chi2_x = pow( ((thistime - phoTmp.t) / xtimeErr_ ) , 2 ) ; 
             string EBorEE = ( isEB ) ? "EB" : "EE" ;
             if ( debugT ) printf(" %s xtal(%d)  t: %.2f, dt: %f,  chi2: %.2f , amp: %.2f  \n", 
                                 EBorEE.c_str(),  (int)ndof, thistime, xtimeErr_,  chi2_x, thisamp );
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
  //cout<<" nSpike = "<<  nSpike <<" nXtl = "<< nSeedXtl <<"  maxSwissX = "<< maxSwissX  << endl ;
  // update ave. time and error
  phoTmp.t     = xtime / xtimeErr ;
  phoTmp.dt    = 1. / sqrt( xtimeErr) ;
  phoTmp.nchi2 = ( ndof != 0 ) ? chi2_bc / ndof : 9999999 ;     
  phoTmp.fSpike = ( nSeedXtl > 0 ) ? (nSpike*1.) / (nSeedXtl*1.) : -1 ;
  phoTmp.nxtals = nXtl ;
  phoTmp.nBC    = nBC ;
  phoTmp.maxSX  = maxSwissX ;

}
/*
void DPAnalysis::EventTime( const edm::Event& iEvent ) {

   Handle<EcalRecHitCollection>        recHitsEB ;
   Handle<EcalRecHitCollection>        recHitsEE ;
   Handle<reco::BasicClusterCollection> EBClusters ;
   Handle<reco::BasicClusterCollection> EEBClusters ;

   iEvent.getByLabel( EBRecHitCollection,     recHitsEB );
   iEvent.getByLabel( EERecHitCollection,     recHitsEE );
   iEvent.getByLabel( EBBasicClusterCollection, EBClusters) ;
   iEvent.getByLabel( EEBasicClusterCollection, EEClusters) ;

   // Barrel BasicClusters
   const reco::BasicClusterCollection* theEBClusters = EBClusters.product () ;
   
   // Endcap BasicClusters
   const reco::BasicClusterCollection* theEEClusters = EEClusters.product () ;

   for (reco::BasicClusterCollection::const_iterator it = theEBClusters->begin () ;
        it != theEBClusters->end ()  ;  ++it)  {

       for (std::vector<std::pair<DetId, float> >::const_iterator detitr = clusterDetIds.begin () ; 
	      detitr != clusterDetIds.end () ; ++detitr) {// loop on rechics of barrel basic clusters

       }
   }

}
*/
/*
double DPAnalysis::HLTMET( Handle<reco::PFJetCollection> jets, vector<const reco::Muon*>& selectedMuons , bool addMuon ) {

   double hltHT  = 0. ;
   double hltMEx = 0. ;
   double hltMEy = 0. ;
   int nj_mht = 0 ;
   for(reco::PFJetCollection::const_iterator it = jets->begin(); it != jets->end(); it++) {
       // fiducial cuts
       if ( it->pt() < 40. || fabs( it->eta() ) > 999.  ) continue ;
       hltHT += it->pt() ;
       hltMEx -= it->p4().Px() ;
       hltMEy -= it->p4().Py() ;
       nj_mht++ ;
   }
 
   if ( addMuon ) {  
      for( size_t i =0 ; i < selectedMuons.size() ; i++) {
         if ( selectedMuons[i]->pt() < 40 || fabs( selectedMuons[i]->eta() ) > 999. ) continue ;
         hltMEx -= selectedMuons[i]->p4().Px() ;
         hltMEy -= selectedMuons[i]->p4().Py() ;
         nj_mht++ ;
      }
   }

   double hltMET = ( nj_mht < 1 ) ? 0 : sqrt( (hltMEx*hltMEx) + (hltMEy*hltMEy) ) ;
   leaves.hltMet   = hltMET ;
   leaves.hltMetPx = hltMEx ;
   leaves.hltMetPy = hltMEy ;

   return hltMET ;
}
*/

bool DPAnalysis::JetSelection( Handle<reco::PFJetCollection> jets, vector<const reco::Photon*>& selectedPhotons, 
                               vector<const reco::PFJet*>& selectedJets) {

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

       if ( k >= MAXJET ) break ;
       selectedJets.push_back( &(*it) ) ;
       leaves.jetPx[k] = it->p4().Px() ;
       leaves.jetPy[k] = it->p4().Py() ;
       leaves.jetPz[k] = it->p4().Pz() ;
       leaves.jetE[k]  = it->p4().E()  ;
       leaves.jetNDau[k] = it->numberOfDaughters() ;
       leaves.jetCM[k]   = it->chargedMultiplicity() ;
       leaves.jetCEF[k]  = it->chargedEmEnergyFraction() ;
       leaves.jetCHEF[k] = it->chargedHadronEnergyFraction();  // New Variable for Jet ID Offline.
       leaves.jetNHF[k]  = it->neutralHadronEnergyFraction() ;  
       leaves.jetNEF[k]  = it->neutralEmEnergyFraction() ;
       leaves.jetEta[k]  = it->eta();                         // New Variable for Jet ID offline.
       k++ ;
   }
   leaves.nJets = (int)( selectedJets.size() ) ;

   if ( selectedJets.size() > 0 )  return true ; 
   else                            return false ;    

}

bool DPAnalysis::ElectronSelection( Handle<reco::GsfElectronCollection> electrons, 
                                    vector<const reco::GsfElectron*>& selectedElectrons ) {

   // Electron Identification Based on Simple Cuts
   // https://twiki.cern.ch/twiki/bin/view/CMS/SimpleCutBasedEleID#Selections_and_How_to_use_them
   int k = 0 ;
   for(reco::GsfElectronCollection::const_iterator it = electrons->begin(); it != electrons->end(); it++) {
       if ( it->pt() < electronCuts[0] || fabs( it->eta() ) > electronCuts[1] ) continue ;
       // Isolation Cuts
       float ecalSumEt = ( it->isEB() ) ? max(0., it->dr03EcalRecHitSumEt() - 1. ) : it->dr03EcalRecHitSumEt();
       float hcalSumEt = it->dr03HcalTowerSumEt();
       float trkSumPt  = it->dr03TkSumPt();  
       double relIso   = (ecalSumEt + hcalSumEt + trkSumPt) / it->pt() ;
       if ( relIso > electronCuts[2] &&  it->isEB() ) continue ;
       if ( relIso > electronCuts[3] && !it->isEB() ) continue ;

       // code for CMSSW < 71X
       // double nLost = it->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() ;
       //if ( nLost >= electronCuts[4]  ) continue ;
       // Code for CMSSW > 71X
       double  nLost =  it->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
       if( nLost >= electronCuts[4] ) continue ;

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
       k++;
   }
   leaves.nElectrons = (int)( selectedElectrons.size() ) ;

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
       double relIso =0. ;
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
     bool IsBack2Back = false ; 
     
     if (  selectedJets.size() > 0 && selectedPhotons.size() > 0  ) {
       double d_R      = ROOT::Math::VectorUtil::DeltaR( selectedJets[0]->p4(), selectedPhotons[0]->p4() ) ;
           
       double dphi_ = ROOT::Math::VectorUtil::DeltaPhi( selectedJets[0]->p4(), selectedPhotons[0]->p4() ) ;
       double PtRatio = selectedJets[0]->pt() / selectedPhotons[0]->pt() ;
       if ( d_R > (2.*3.1416/3.) && PtRatio > 0.7 && PtRatio < 1.3 )  isGammaJets = true ;
       
       if( (fabs(dphi_) >= 7*3.1416/8) && (d_R >= 0.5)) IsBack2Back = true;  // Loose Tagging of Gamma Jet Events;
        
       leaves.IsMyGammaJet_Tag = IsBack2Back ;
       leaves.IsSCGammaJet_Tag = isGammaJets ;
       leaves.phi_Gamma_Jet = dphi_;
       leaves.dR_Gamma_Jet  = d_R ; 
 
     }
     return isGammaJets ;
}


//define this as a plug-in
DEFINE_FWK_MODULE(DPAnalysis);
