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
// $Id$
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
   isData               = iConfig.getParameter<bool> ("isData");
   L1Select             = iConfig.getParameter<bool> ("L1Select");

   gen = new GenStudy( iConfig );

   theFile  = new TFile( rootFileName.c_str(), "RECREATE") ;
   theFile->cd () ;
   theTree  = new TTree ( "DPAnalysis","DPAnalysis" ) ;
   setBranches( theTree, leaves ) ;

   firedTrig.clear() ;
   for ( size_t i=0; i< triggerPatent.size(); i++ ) firedTrig.push_back(-1) ;  

   // reset the counter
   for ( int i=0; i< 10 ; i++) counter[i] = 0 ;

   // initialize the time corrector
   theTimeCorrector_.initEB("EB");
   theTimeCorrector_.initEE("EElow");
   runID_ = 0 ;
}


DPAnalysis::~DPAnalysis()
{
   // do anything here that needs to be done at desctruction time

   delete gen ;
   cout<<"All:"<< counter[0]<<" Trigger:"<<counter[1]<<" Vertex:"<< counter[2] <<" photon:"<<counter[3] ;
   cout<<" beamHalo:"<< counter[4] <<" BeamHalo:"<< counter[5] <<" Jet:"<<counter[6] <<" MET:"<<counter[7] <<endl ;

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
   bool passL1 = L1TriggerSelection( iEvent, iSetup ) ;

   // HLT trigger analysis
   Handle<edm::TriggerResults> triggers;
   iEvent.getByLabel( trigSource, triggers );
   const edm::TriggerNames& trgNameList = iEvent.triggerNames( *triggers ) ;

   TriggerTagging( triggers, trgNameList, run_id, firedTrig ) ;
   bool passHLT = TriggerSelection( triggers, firedTrig ) ;

   // Using L1 or HLT to select events ?!
   bool passTrigger = ( L1Select ) ? passL1 : passHLT  ;

   if ( passTrigger ) counter[1]++ ;  

   if ( !isData ) { 
      gen->GetGenEvent( iEvent, leaves );
      //gen->GetGen( iEvent, leaves );
   }
   //if ( !isData ) gen->PrintGenEvent( iEvent );

   bool pass = EventSelection( iEvent ) ;
   
   if ( pass && !isData ) theTree->Fill();
   if ( pass && isData && passTrigger ) theTree->Fill();
}

bool DPAnalysis::EventSelection(const edm::Event& iEvent ) {

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

   bool hasGoodVtx = VertexSelection( recVtxs );
   if ( !hasGoodVtx ) passEvent = false ;
   if ( passEvent )   counter[2]++ ;  

   selectedPhotons.clear() ;
   PhotonSelection( photons, recHitsEB, recHitsEE, tracks, selectedPhotons ) ;
   if ( selectedPhotons.size() < (size_t)photonCuts[6] )  passEvent = false ;
   if ( passEvent )   counter[3]++ ;  

   //sMinorSelection( selectedPhotons, recHitsEB, recHitsEE ) ;
   //if ( selectedPhotons.size() < (size_t)photonCuts[5] )  passEvent = false ;
   /*
   if ( passEvent ) {  counter[4]++ ;  
                       printf("id:%d, nPho:%d, Pt:%.3f sMin:%.6f \n", 
                iEvent.id().event(), (int)selectedPhotons.size(), selectedPhotons[0]->pt(), sMin_ ) ;
   }   
   */
   if( beamHaloSummary.isValid() ) {
     const reco::BeamHaloSummary TheSummary = (*beamHaloSummary.product() ); 
     if( !TheSummary.CSCTightHaloId() && passEvent ) { 
       counter[4]++ ;  
     } else {
       passEvent = false ;
     }
   } else {
       counter[4]++ ;
   }
 
   //IsoPhotonSelection( selectedPhotons ) ;
   //if ( selectedPhotons.size() < photonCuts[5] )  passEvent = false ;
   if ( passEvent )   counter[5]++ ;  

   selectedJets.clear() ;
   JetSelection( jets, selectedPhotons, selectedJets );
   //bool isGammaJets = GammaJetVeto( selectedPhotons, selectedJets ) ;
   //if ( isGammaJets ) passEvent = false ;
   if ( selectedJets.size() < jetCuts[3] )   passEvent = false ;
   if ( passEvent )   counter[6]++ ;   

   
   selectedElectrons.clear() ;
   ElectronSelection( electrons, selectedElectrons ) ;

   selectedMuons.clear() ;
   MuonSelection( muons, selectedMuons );
   
   const reco::PFMET pfMet = (*met)[0] ;
   leaves.met   = pfMet.et() ;
   leaves.metPx = pfMet.px() ;
   leaves.metPy = pfMet.py() ;
   if ( pfMet.pt() < metCuts[0]  ) passEvent = false;
   if ( passEvent )   counter[7]++ ;  

   return passEvent ;
}

bool DPAnalysis::L1TriggerSelection( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

    // Get L1 Trigger menu
    edm::ESHandle<L1GtTriggerMenu> menuRcd;
    iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
    const L1GtTriggerMenu* menu = menuRcd.product();
    // Get L1 Trigger record  
    edm::Handle< L1GlobalTriggerReadoutRecord > gtRecord;
    iEvent.getByLabel( edm::InputTag("gtDigis"), gtRecord);
    // Get dWord after masking disabled bits
    const DecisionWord dWord = gtRecord->decisionWord();
 
   bool l1_accepted = menu->gtAlgorithmResult( l1GTSource , dWord);
   //int passL1 =  ( l1SingleEG22 ) ? 1 : 0 ; 

   //cout<<" pass L1 EG22 ? "<<  passL1 <<endl ;
   if ( l1_accepted ) leaves.L1a = 1 ;

   return l1_accepted ;
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
      leaves.triggered = (int)(trgbits) ;
      pass = true ;
   }

   return pass ;
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
       cout<<" name: "<< tName <<"  idx:"<< trgIndex <<"  accept:"<< trgResult <<endl;
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
       ClusterTime( it->superCluster(), recHitsEB , recHitsEE, phoTmp );
       //cout<<" 2nd xT : "<< aveXtalTime <<"  xTE : "<< aveXtalTimeErr << endl;
       leaves.aveTime1[k]     = phoTmp.t ;    // weighted ave. time of seed cluster
       leaves.aveTimeErr1[k]  = phoTmp.dt ;
       leaves.timeChi2[k]     = phoTmp.nchi2 ;
       leaves.nXtals[k]       = phoTmp.nxtals ;
       leaves.nBC[k]          = phoTmp.nBC ;
       leaves.fSpike[k]       = phoTmp.fSpike ;
       leaves.maxSwissX[k]    = phoTmp.maxSX ; 

       // refine the timing 
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
       leaves.aveTime[k]      = phoTmp.t ;       // weighted ave. time of all clusters
       leaves.aveTimeErr[k]   = phoTmp.dt ;

       selectedPhotons.push_back( &(*it) ) ;
       k++ ;
   }
   leaves.nPhotons = k ;
   //leaves.nPhotons = (int)( selectedPhotons.size() ) ;

   if ( selectedPhotons.size() > 0 && maxPt >= photonCuts[7] )  return true ; 
   else                               return false ;    

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
             double chi2_x = pow( (thistime - phoTmp.t / xtimeErr_ ) , 2 ) ; 
             chi2_bc += chi2_x ;
             ndof += 1 ;
             // remove un-qualified hits 
             if ( fabs ( thistime - phoTmp.t ) > 3.*phoTmp.dt ) continue ;
 
             xtime     += thistime / pow( xtimeErr_ , 2 ) ;
             xtimeErr  += 1/ pow( xtimeErr_ , 2 ) ;
             nXtl++ ;
      }
  }

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
       leaves.jetNHF[k]  = it->neutralHadronEnergyFraction() ;  
       leaves.jetNEF[k]  = it->neutralEmEnergyFraction() ;
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
     
     if (  selectedJets.size() > 0 && selectedPhotons.size() > 0  ) {
       double dR      = ROOT::Math::VectorUtil::DeltaR( selectedJets[0]->p4(), selectedPhotons[0]->p4() ) ;
       double PtRatio = selectedJets[0]->pt() / selectedPhotons[0]->pt() ;
       if ( dR > (2.*3.1416/3.) && PtRatio > 0.7 && PtRatio < 1.3 )  isGammaJets = true ;
     }
     return isGammaJets ;
}


//define this as a plug-in
DEFINE_FWK_MODULE(DPAnalysis);
