// -*- C++ -*-
// Package:    CosmicAnalysis
// Class:      CosmicAnalysis
// 
/**\class CosmicAnalysis CosmicAnalysis.cc EXO/CosmicAnalysis/src/CosmicAnalysis.cc

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
#include "CosmicAnalysis.h"
#include "Ntuple.h"

#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"

// For DT Segment
#include "Geometry/DTGeometry/interface/DTChamber.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"

// global tracking geometry
//#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
//#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

using namespace cms ;
using namespace edm ;
using namespace std ;

// constants, enums and typedefs
// static data member definitions
static bool EDecreasing( iCluster s1, iCluster s2) { return ( s1.E > s2.E ); }
static bool dRIncreasing( iDT s1, iDT s2) { return ( s1.dR < s2.dR ); }
static bool dPhiIncreasing( iDT s1, iDT s2) { return ( s1.dPhi < s2.dPhi ); }

// constructors and destructor
CosmicAnalysis::CosmicAnalysis(const edm::ParameterSet& iConfig){

   //now do what ever initialization is needed
   rootFileName         = iConfig.getUntrackedParameter<string> ("rootFileName");
   trigSource           = iConfig.getParameter<edm::InputTag> ("trigSource");
   l1GTSource           = iConfig.getParameter<string> ("L1GTSource");
   muonSource           = iConfig.getParameter<edm::InputTag> ("muonSource");
   metSource            = iConfig.getParameter<edm::InputTag> ("metSource");
   //photonSource         = iConfig.getParameter<edm::InputTag> ("photonSource");
   //jetSource            = iConfig.getParameter<edm::InputTag> ("jetSource");
   scSource             = iConfig.getParameter<edm::InputTag> ("scSource");
   bcSource             = iConfig.getParameter<edm::InputTag> ("bcSource");
   EBRecHitCollection   = iConfig.getParameter<edm::InputTag> ("EBRecHitCollection") ;
   EERecHitCollection   = iConfig.getParameter<edm::InputTag> ("EERecHitCollection") ;
   DTSegmentTag         = iConfig.getParameter<edm::InputTag> ("DTSegmentCollection") ;
   CSCSegmentTag        = iConfig.getParameter<edm::InputTag> ("CSCSegmentCollection") ;
   muonCuts             = iConfig.getParameter<std::vector<double> >("muonCuts");

   triggerPatent        = iConfig.getParameter< std::vector<string> >("triggerName");
   L1Select             = iConfig.getParameter<bool> ("L1Select");

   const InputTag TrigEvtTag("hltTriggerSummaryAOD","","HLT");
   trigEvent            = iConfig.getUntrackedParameter<edm::InputTag>("triggerEventTag", TrigEvtTag);

   theFile  = new TFile( rootFileName.c_str(), "RECREATE") ;
   theFile->cd() ;
   theTree  = new TTree ( "CosmicAnalysis","CosmicAnalysis" ) ;

   // SetBranches
   setCosmicBranches( theTree, leaves ) ;

   targetTrig = 0 ;
   firedTrig.clear() ;
   for ( size_t i=0; i< triggerPatent.size(); i++ ) firedTrig.push_back(-1) ;  

   runID_ = 0 ;
   debugT = false ;

}


CosmicAnalysis::~CosmicAnalysis()
{
   // do anything here that needs to be done at desctruction time
   theTree->Print() ;
   theFile->cd () ;
   theTree->Write() ; 
   theFile->Close() ;

}

//
// member functions
//

// ------------ method called for each event  ------------
void CosmicAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

   // get calibration service
   // IC's
   iSetup.get<EcalIntercalibConstantsRcd>().get(ical);
   // ADCtoGeV
   iSetup.get<EcalADCToGeVConstantRcd>().get(agc);
   // transp corrections
   iSetup.get<EcalLaserDbRecord>().get(laser);
   // Geometry
   iSetup.get<CaloGeometryRecord> ().get(pGeometry) ;
   theGeometry = pGeometry.product() ;

   // event time
   eventTime = iEvent.time() ;
   // Initialize ntuple branches
   InitializeCosmicBranches();
   
   leaves.runId       = iEvent.id().run() ;
   leaves.eventId     = iEvent.id().event() ;

   if ( runID_ == 0 )  PrintTriggers( iEvent ) ;
   int run_id    = iEvent.id().run()  ;

   // L1 Trigger Selection
   passL1 = L1TriggerSelection( iEvent, iSetup ) ;
   //if ( passL1 ) cout<<" Pass L1 "<<endl ;

   // HLT trigger analysis
   Handle<edm::TriggerResults> triggers;
   iEvent.getByLabel( trigSource, triggers );
   const edm::TriggerNames& trgNameList = iEvent.triggerNames( *triggers ) ;

   TriggerTagging( triggers, trgNameList, run_id, firedTrig ) ;
   passHLT = TriggerSelection( triggers, firedTrig ) ;

   // Using L1 or HLT to select events ?!
   bool passTrigger = ( L1Select ) ? passL1 : passHLT  ;
   bool pass        = EventSelection( iEvent, iSetup ) ;

   //if ( pass ) cout<<"  5.0 " ;
   //else        cout<<"  5.1 " ;
   // fill the ntuple
 
   if ( pass ) {
      int err_code = theTree->Fill();
      if ( err_code < 1) cout<<" err: " << err_code <<endl; 
   }
   //cout<<"  6 "<<endl ;
}

bool CosmicAnalysis::EventSelection(const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

   Handle<edm::TriggerResults>         triggers;
   //Handle<reco::PhotonCollection>      photons; 
   Handle<reco::MuonCollection>        muons; 
   //Handle<reco::CaloJetCollection>     jets; 
   Handle<reco::METCollection>         met; 
   Handle<EcalRecHitCollection>        recHitsEB ;
   Handle<EcalRecHitCollection>        recHitsEE ;
   Handle<reco::SuperClusterCollection>  pSuperClusters;
   Handle<reco::BasicClusterCollection>  pBarrelBasicClusters ;

   iEvent.getByLabel( trigSource,     triggers );
   //iEvent.getByLabel( photonSource,   photons  );
   iEvent.getByLabel( muonSource,     muons );
   //iEvent.getByLabel( jetSource,      jets  );
   iEvent.getByLabel( metSource,      met  );
   iEvent.getByLabel( EBRecHitCollection,     recHitsEB );
   iEvent.getByLabel( EERecHitCollection,     recHitsEE );
   iEvent.getByLabel( scSource,   pSuperClusters ) ;
   iEvent.getByLabel (bcSource,   pBarrelBasicClusters) ;
   
   bool passEvent = true ;

   //cout<<" ----- new event ------ "<<endl ;
   // find trigger matched objects
   const reco::MET MET_ = (*met)[0] ;

   // for cscsegments and halo muon/photon studies 
   Handle<CSCSegmentCollection>       cscSegments ;
   iEvent.getByLabel( CSCSegmentTag,  cscSegments );
   //BeamHaloMatch( cscSegments, selectedPhotons, iSetup ) ;
   Handle<DTRecSegment4DCollection>   dtSegments ;
   iEvent.getByLabel( DTSegmentTag,   dtSegments );

   //cout<<" =================================== "<<endl ;
   //cout<<"  Muon size: "<< distance( muons->begin() , muons->end() )  ;
   //cout<<"  CSC size: "<< distance( cscSegments->begin() , cscSegments->end() ) ;
   //cout<<"  DT  size: "<< distance( dtSegments->begin() ,  dtSegments->end() ) <<endl ;
   

   // SuperCluster Info
   selectSC.clear() ;
   ClusterInfo( pSuperClusters, recHitsEB, recHitsEE, selectSC ) ;
   // DT and CSC Segment Info
   selectDT.clear() ;
   DTSegmentInfo( dtSegments, iSetup, selectDT ) ;
   selectCSC.clear() ;
   CSCSegmentInfo( cscSegments, iSetup, selectCSC ) ;
   
   //printf(" Selected => nSC: %d, nDT: %d , nCSC %d \n ", (int)selectSC.size() , (int)selectDT.size(), (int)selectCSC.size() ) ;
   CosmicRayMatch( selectDT, selectCSC, selectSC ) ;

   int ii = 0 ; 
   for ( size_t i=0; i < selectDT.size(); i++ ) {
       if ( (i+1) > MAXDT ) break ;
       leaves.dtX[i] = selectDT[i].x ;
       leaves.dtY[i] = selectDT[i].y ;
       leaves.dtZ[i] = selectDT[i].z ;
       leaves.dtdX[i] = selectDT[i].dx ;
       leaves.dtdY[i] = selectDT[i].dy ;
       leaves.dtdZ[i] = selectDT[i].dz ;
       leaves.dt_dR[i] = selectDT[i].dR ;
       ii++ ; 
   }
   leaves.nDT   = ii ;
   //cout<<"  Save "<< ii <<" DT " ;
  
   int jj = 0 ;
   for ( size_t i=0; i < selectCSC.size(); i++ ) {
       if ( (i+1) > MAXCSC ) break ;
       leaves.cscX[i] = selectCSC[i].x ;
       leaves.cscY[i] = selectCSC[i].y ;
       leaves.cscZ[i] = selectCSC[i].z ;
       leaves.cscdX[i] = selectCSC[i].dx ;
       leaves.cscdY[i] = selectCSC[i].dy ;
       leaves.cscdZ[i] = selectCSC[i].dz ;
       leaves.csc_dPhi[i] = selectCSC[i].dPhi ;
       jj++ ;
   } 
   leaves.nCSC  = jj ;
   //cout<<" and "<< jj <<" CSC "<<endl ;

   MuonSelection( muons ) ;

   // MET Info
   leaves.met   = MET_.et() ;
   leaves.metPx = MET_.px() ;
   leaves.metPy = MET_.py() ;
   leaves.nSC   = ( selectSC.size() > 15 ) ? 15 : (int)selectSC.size()  ;
   if ( selectSC.size() < 1 ) passEvent = false ;

   /*
   printf(" ---- Leaves Info ----  \n" ) ;
   printf("   nSC: %d, nCSC: %d , nDT: %d , nMu: %d \n", leaves.nSC, leaves.nCSC, leaves.nDT, leaves.nMuons ) ;
   int ik = ( ii < 1 ) ? ii : ii -1 ;
   int jk = ( jj < 1 ) ? jj : jj -1 ;
   int sk = ( leaves.nSC < 1 ) ? leaves.nSC : leaves.nSC - 1 ;
   int mk = ( leaves.nMuons < 1 ) ? leaves.nMuons : leaves.nMuons - 1 ;
   printf("   xtalE: %.2f, csc_dPhi: %.2f, dt_dR: %.2f \n", leaves.xtalE[sk], leaves.csc_dPhi[jk], leaves.dt_dR[ik] ) ;
   printf("   muE: %.2f, mu_nCSC: %d, mu_nDT: %d \n", leaves.muE[mk], leaves.mu_nCSC[mk], leaves.mu_nDT[mk] ) ;
   */
   return passEvent ;
}


bool CosmicAnalysis::L1TriggerSelection( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

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

void CosmicAnalysis::TriggerTagging( Handle<edm::TriggerResults> triggers, const edm::TriggerNames& trgNameList, int RunId, vector<int>& firedTrig ) {

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
bool CosmicAnalysis::TriggerSelection( Handle<edm::TriggerResults> triggers, vector<int> firedTrigID ) {

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


void CosmicAnalysis::PrintTriggers( const edm::Event& iEvent ) {

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

bool CosmicAnalysis::CosmicRayMatch( vector<iDT>& dtsegV, vector<iDT>& cscsegV, vector<iCluster>& scV ) {

   //bool cosmicmatch = false ;
   int k = 0 ;
   for ( size_t i=0; i< scV.size() ; i++ ) {

       //cout<<" fk - "<<i <<" dt size : "<<  distance( dtSeg->begin(), dtSeg->end() ) <<endl ;
       //cout<<"  ecalR: " << sqrt( (scV[i].x*scV[i].x) + (scV[i].y*scV[i].y) ) << endl ;
       float dPhi = 99. ;
       float dEta = 99. ;
       float dR   = 99. ;
       if ( i > 14 ) break ;
       
       float ecalRho = sqrt( (scV[i].x*scV[i].x) + (scV[i].y*scV[i].y) );
       if  ( ecalRho < 128. || fabs(scV[i].z) > 390. ) continue ;

       for ( size_t j=0 ; j < dtsegV.size(); j++ ) {

           //double dx = dtsegV[j].x - scV[i].x ;
           //double dy = dtsegV[j].y - scV[i].y ;
           //double dz = dtsegV[j].z - scV[i].z ;
           //printf(" dir( %f, %f, %f, %f ) \n ", gv.x() , gv.y(), gv.z() , gv.mag() ) ;

           // Propagate to ECAL Surface  
           float pjX = dtsegV[j].x ;
           float pjY = dtsegV[j].y ;
           float pjZ = dtsegV[j].z ;
           float pjRho = sqrt( (pjX*pjX) + (pjY*pjY) ) ;
           float step  = 10. ;
           bool hitECAL = false ;
           float gvMag = sqrt( (dtsegV[j].dx*dtsegV[j].dx) + (dtsegV[j].dy*dtsegV[j].dy) + (dtsegV[j].dz*dtsegV[j].dz) ) ;
           //printf("== Dir ( %.1f, %.1f, %.1f) \n", dtsegV[j].dx , dtsegV[j].dy , dtsegV[j].dz ) ;
           float rho_ = pjRho ;
           while ( !hitECAL )  {
                pjX =  pjX - ( step* dtsegV[j].dx / gvMag ) ;
                pjY =  pjY - ( step* dtsegV[j].dy / gvMag ) ;
                pjZ =  pjZ - ( step* dtsegV[j].dz / gvMag ) ;
                pjRho = sqrt( (pjX*pjX) + (pjY*pjY) ) ;
                if ( pjRho < rho_ ) rho_    = pjRho ;
                if ( pjRho > rho_ ) hitECAL = true ; // actually not hit Ecal but out of ECAL range
                //printf("     ( %.1f, %.1f, %.1f ) \n", pjX, pjY, pjZ ) ;
                if ( pjRho < ecalRho || fabs(pjZ) > 341. ) hitECAL = true ;
           } ;
           float pjR = sqrt( (pjX*pjX) + (pjY*pjY) + (pjZ*pjZ) ) ;
           //printf("*** ECALR: %.1f , ProjR :%.1f, z: %.1f \n", ecalRho, pjR, scV[i].z ) ;

           LorentzVector segEBP4( pjX, pjY, pjZ, pjR ) ;
           float calo_n =  scV[i].E / sqrt( (scV[i].x*scV[i].x) + (scV[i].y*scV[i].y) + (scV[i].z*scV[i].z) ) ; 
           LorentzVector caloP4( scV[i].x*calo_n, scV[i].y*calo_n, scV[i].z*calo_n, scV[i].E );  

           float dR_   = ROOT::Math::VectorUtil::DeltaR( segEBP4, caloP4 ) ;
           float dPhi_ = ROOT::Math::VectorUtil::DeltaPhi( segEBP4, caloP4 ) ;
           float dEta_ = fabs( segEBP4.Eta() - scV[i].eta ) ;
           dtsegV[j].dR = dR_ ;
           if ( dR_ < dR ) {
               dR   = dR_ ;
               dPhi = fabs( dPhi_ ) ;
               dEta = dEta_ ;
           }
           
       }
       sort( dtsegV.begin() , dtsegV.end(), dRIncreasing ) ;
       leaves.dtdEta[k]  = dEta ;
       leaves.dtdPhi[k]  = dPhi ;

       dPhi = 99. ;
       float cscRho  = -1 ;
       for ( size_t j=0 ; j < cscsegV.size(); j++ ) {
           LorentzVector segP4( cscsegV[j].x, cscsegV[j].y, cscsegV[j].z, cscsegV[j].L ) ;
           float calo_n =  scV[i].E / sqrt( (scV[i].x*scV[i].x) + (scV[i].y*scV[i].y) + (scV[i].z*scV[i].z) ) ; 
           LorentzVector caloP4( scV[i].x*calo_n, scV[i].y*calo_n, scV[i].z*calo_n, scV[i].E );  
           
           float dPhi_ = ROOT::Math::VectorUtil::DeltaPhi( segP4, caloP4 ) ;
           float rho = sqrt( (cscsegV[j].x*cscsegV[j].x) + (cscsegV[j].y*cscsegV[j].y) ) ;
           cscsegV[j].dPhi = fabs( dPhi_ ) ;

           if ( fabs(dPhi_) < dPhi ) {
              dPhi = fabs( dPhi_ ) ;
              cscRho  = rho ;
           }
       }
       sort( cscsegV.begin() , cscsegV.end(), dPhiIncreasing ) ;
       leaves.cscRho[k]   = cscRho ;
       leaves.cscdPhi[k]  = dPhi ;

       leaves.xtalEta[k]    = scV[i].eta ;
       leaves.xtalPhi[k]    = scV[i].phi ;
       leaves.xtalSigma[k]  = scV[i].calib ;
       leaves.xtalSwissX[k] = scV[i].swissX ;
       leaves.xtalChi2[k]   = scV[i].chi2 ;
       leaves.xtal_x[k]     = scV[i].x ;
       leaves.xtal_y[k]     = scV[i].y ;
       leaves.xtal_z[k]     = scV[i].z ;
       leaves.xtal_t[k]     = scV[i].t ;
       leaves.xtalE[k]      = scV[i].E ;
       leaves.bc_nXtals[k]  = scV[i].nXtals ;
       leaves.bc_nBC[k]     = scV[i].nBC ;
       leaves.bc_sMaj[k]    = scV[i].sMaj ;
       leaves.bc_sMin[k]    = scV[i].sMin ;
       leaves.bc_E[k]       = scV[i].scE ;
       k++ ;
  }
  return true ;

}

void CosmicAnalysis::CSCSegmentInfo( Handle<CSCSegmentCollection> cscSeg, const EventSetup& iSetup, vector<iDT>& selectCSC ) {

       ESHandle<CSCGeometry> cscGeom;
       iSetup.get<MuonGeometryRecord>().get(cscGeom);

       //int k = 0 ;
       for ( CSCSegmentCollection::const_iterator it = cscSeg->begin(); it != cscSeg->end(); it++) {

           if ( !it->isValid() ) continue ;
           CSCDetId DetId = it->cscDetId();
           const CSCChamber* cscchamber = cscGeom->chamber( DetId );
           GlobalPoint  gp = cscchamber->toGlobal( it->localPosition()  );
           float gpMag = sqrt( (gp.x()*gp.x()) + (gp.y()*gp.y()) + (gp.z()*gp.z()) ) ;

           //LorentzVector segP4( gp.x(), gp.y(), gp.z(), gpMag ) ;
           // Get segment direction in DT
           GlobalVector gv = cscchamber->toGlobal( it->localDirection() ) ;

           iDT theSegment ;
           theSegment.x = gp.x() ;
           theSegment.y = gp.y() ;
           theSegment.z = gp.z() ;
           theSegment.L = gpMag ;
           theSegment.dx = gv.x() ;
           theSegment.dy = gv.y() ;
           theSegment.dz = gv.z() ;
           theSegment.dPhi = 3.16 ;
           selectCSC.push_back( theSegment ) ;

           /*
           if ( k >= MAXDT ) continue ;
           leaves.cscX[k] = gp.x() ;
           leaves.cscY[k] = gp.y() ;
           leaves.cscZ[k] = gp.z() ;
           leaves.cscdX[k] = gv.x() ;
           leaves.cscdY[k] = gv.y() ;
           leaves.cscdZ[k] = gv.z() ;
           */
           //double dEta_ = fabs( gp.eta() - selectedPhotons[i]->eta() )  ;
           //dEta = ( dEta_ < dEta ) ? dEta_ : dEta ;
           //double dPhi_ = fabs( gp.phi() - selectedPhotons[i]->phi() )  ;
           //k++ ;
       }

}

// !!! Only work for RECO - dtSegment is not available for AOD 
void CosmicAnalysis::DTSegmentInfo( Handle<DTRecSegment4DCollection> dtSeg, const EventSetup& iSetup, vector<iDT>& selectDT ) {

       ESHandle<DTGeometry> dtGeom;
       iSetup.get<MuonGeometryRecord>().get(dtGeom);
       
       //int k = 0 ;
       //cout<<" dt size: "<< distance( dtSeg->begin(), dtSeg->end() ) ;
       for (DTRecSegment4DCollection::const_iterator it = dtSeg->begin(); it != dtSeg->end(); it++) {
           if ( !it->isValid()) continue ;
           if ( !it->hasPhi() ) continue ;
           if ( !it->hasZed() ) continue ;
           
           // Get the corresponding DTChamber
           DetId id = it->geographicalId();
           DTChamberId chamberId(id.rawId());
           if ( chamberId.station() > 1 ) continue ;  // only look at segment from inner most chambers
           const DTChamber* dtchamber = dtGeom->chamber( chamberId ) ;

           // Get segment position in DT
	   GlobalPoint  gp = dtchamber->toGlobal( it->localPosition()  );
           float gpMag = sqrt( (gp.x()*gp.x()) + (gp.y()*gp.y()) + (gp.z()*gp.z()) ) ;
	   //LorentzVector segP4( gp.x(), gp.y(), gp.z(), gpMag ) ;

           // Get segment direction in DT
           GlobalVector gv = dtchamber->toGlobal( it->localDirection() ) ;

           iDT theSegment ;
           theSegment.x = gp.x() ;
           theSegment.y = gp.y() ;
           theSegment.z = gp.z() ;
           theSegment.L = gpMag ;
           theSegment.dx = gv.x() ;
           theSegment.dy = gv.y() ;
           theSegment.dz = gv.z() ;
           theSegment.dR = 99. ;
           selectDT.push_back( theSegment ) ;

           /*
           if ( k >= MAXDT ) continue ;
           leaves.dtX[k] = gp.x() ;
           leaves.dtY[k] = gp.y() ;
           leaves.dtZ[k] = gp.z() ;
           leaves.dtdX[k] = gv.x() ;
           leaves.dtdY[k] = gv.y() ;
           leaves.dtdZ[k] = gv.z() ;
           k++ ;
           */
       }

}

void CosmicAnalysis::ClusterInfo( Handle<reco::SuperClusterCollection> pSuperClusters, Handle<EcalRecHitCollection> recHitsEB, Handle<EcalRecHitCollection> recHitsEE, vector<iCluster>& sc_V ) {

  const reco::SuperClusterCollection* scCollection = pSuperClusters.product();
  const EcalIntercalibConstantMap& icalMap = ical->getMap();
  float adcToGeV_EB = float(agc->getEBValue());
  float adcToGeV_EE = float(agc->getEEValue());

  iCluster theCluster ; 
  sc_V.clear() ;
  for (reco::SuperClusterCollection::const_iterator scIt = scCollection->begin();   scIt != scCollection->end(); scIt++) {
      //printf(" == SC E: %.1f \n", scIt->energy() ) ;

      // Pick-up seed basic cluster
      // for ( reco::CaloCluster_iterator bcIt = scIt->clustersBegin() ; bcIt != scIt->clustersEnd() ;  ++bcIt) 
      int nBC = distance( scIt->clustersBegin(),  scIt->clustersEnd() ) ;
      theCluster.nBC = nBC ;
      theCluster.scE = scIt->energy() ;

      reco::CaloClusterPtr bcPtr = scIt->seed() ;
      //printf(" === BC E: %.1f \n", (*bcPtr)->energy() ) ;
      //printf(" === BC E: %.1f \n", bcPtr->energy() ) ;

      // Get DetID
      vector< std::pair<DetId, float> > clusterDetIds = bcPtr->hitsAndFractions() ; 
      int NGoodXtals = 0 ;
      //loop on xtals in cluster 
      for ( size_t i=0; i < clusterDetIds.size() ; i++ ) {
          // Find the seed crystal
          if ( bcPtr->seed() != clusterDetIds[i].first ) continue ;
          if ( clusterDetIds[i].first.det()   != DetId::Ecal ) continue ;
          //if ( clusterDetIds[i].first.subdetId() != EcalBarrel  ) continue ;             

	  bool isEB = ( clusterDetIds[i].first.subdetId () == EcalBarrel)  ? true : false ;

          // Get Rechits
	  EcalRecHit hit = ( isEB) ?  *(recHitsEB->find( clusterDetIds[i].first ) ) : 
                                      *(recHitsEE->find( clusterDetIds[i].first ) ) ;

	  if( !( hit.checkFlag(EcalRecHit::kGood) || hit.checkFlag(EcalRecHit::kOutOfTime) ) )   continue;

          // Get sMaj and sMin
          const EcalRecHitCollection* rechitCollection = ( isEB ) ? recHitsEB.product() : recHitsEE.product() ;
          Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*bcPtr, *rechitCollection );
          float sMin =  moments.sMin  ;
          float sMaj =  moments.sMaj  ;
          theCluster.sMin = sMin ;
          theCluster.sMaj = sMaj ;

	  float amp  = hit.energy() ;
	  float chi2 = hit.chi2() ;
	  float time = hit.time() ;
          theCluster.E    = amp ;
          theCluster.chi2 = chi2 ;
          theCluster.t    = time ;
	  // Swiss cross cleaning 
	  float swissX = (isEB) ? EcalTools::swissCross( clusterDetIds[i].first, *recHitsEB , 0., true ) : 
		                  EcalTools::swissCross( clusterDetIds[i].first, *recHitsEE , 0., true ) ;
          theCluster.swissX = swissX ;

	  // Get inter-calibration
	  EcalIntercalibConstantMap::const_iterator icalit = icalMap.find( clusterDetIds[i].first);
	  EcalIntercalibConstant icalconst = ( icalit!=icalMap.end() )  ?  (*icalit) : 1 ;
	  // get laser coefficient
	  float lasercalib = laser->getLaserCorrection( clusterDetIds[i].first , eventTime );
	  float adcToGeV = ( isEB ) ? adcToGeV_EB : adcToGeV_EE ;
          float calib = icalconst*lasercalib*adcToGeV ;
          theCluster.calib = calib ;

	  // Get Xtal position
	  GlobalPoint pos = theGeometry->getPosition( hit.detid() );
          theCluster.eta =  pos.eta() ;
          theCluster.phi =  pos.phi() ;
          theCluster.x =  pos.x() ;
          theCluster.y =  pos.y() ;
          theCluster.z =  pos.z() ;

          if ( amp/calib >= (12*1.1) ) NGoodXtals++ ;
	  // discard rechits with A/sigma < 12
	  //printf(" ==== Hit E: %.1f , chi2: %.2f , t: %.2f , sX: %.2f , eta: %.2f \n", amp, chi2, time, swissX, pos.eta() ) ;

      } // end of Xtal loop
      theCluster.nXtals =  NGoodXtals ;
      sc_V.push_back( theCluster ) ;
  

      /*
      if ( k >= MAXPHO ) continue ;
      leaves.xtalEta[k] = theCluster.eta ;
      leaves.xtalPhi[k] = theCluster.phi ;
      leaves.xtalSigma[k]  = theCluster.calib ;
      leaves.xtalSwissX[k] = theCluster.swissX ;
      leaves.xtalChi2[k]   = theCluster.chi2 ;
      leaves.xtal_x[k] = theCluster.x ;
      leaves.xtal_y[k] = theCluster.y ;
      leaves.xtal_z[k] = theCluster.z ;
      leaves.xtal_t[k] = theCluster.t ;
      leaves.xtalE[k]  = theCluster.E ;
      leaves.bc_nXtals[k] = theCluster.nXtals ;
      leaves.bc_nBC[k]    = theCluster.nBC ;
      leaves.bc_sMaj[k]   = theCluster.sMaj ;
      leaves.bc_sMin[k]   = theCluster.sMin ;
      leaves.bc_E[k]      = theCluster.scE ;
      */

  } // end of SC loop
  //cout<<" N of SC : "<< sc_V.size() << endl ; 
  sort( sc_V.begin(), sc_V.end() , EDecreasing ) ;
  
}

void CosmicAnalysis::MuonSelection( Handle<reco::MuonCollection> muons ) {

   int k = 0;
   for(reco::MuonCollection::const_iterator it = muons->begin(); it != muons->end(); it++) {
       //if ( ! it->isStandAloneMuon() ) continue ;
       if ( k >= MAXCMU ) break ;
       if ( it->pt() < muonCuts[0] || fabs( it->eta() ) > muonCuts[1] ) continue ;
       int nDTseg = 0 ;
       int nCSCseg = 0 ;
       for ( int station = 1; station <= 4; station++ ) { 
           nDTseg  += it->numberOfSegments(station, MuonSubdetId::DT, reco::Muon::SegmentArbitration) ;
           if ( station > 3 ) continue ;
           nCSCseg += it->numberOfSegments(station, MuonSubdetId::CSC, reco::Muon::SegmentArbitration) ;
       }
       //printf(" muon(%d) nDT_Seg : %d, nCSC_Seg : %d \n", k+1, nDTseg, nCSCseg );
       leaves.muPx[k] = it->p4().Px() ;
       leaves.muPy[k] = it->p4().Py() ;
       leaves.muPz[k] = it->p4().Pz() ;
       leaves.muE[k]  = it->p4().E() ;
       leaves.mu_nDT[k] = nDTseg ;
       leaves.mu_nCSC[k] = nCSCseg ;
       k++ ;
   }
   leaves.nMuons = k ;
}


void CosmicAnalysis::InitializeCosmicBranches() {

     leaves.runId   = 0 ;
     leaves.eventId = 0 ;

     leaves.triggered = 0 ;
     leaves.L1a = 0 ;

     leaves.nSC = 0 ;
     leaves.nDT = 0 ;
     leaves.nCSC = 0 ;
     leaves.nMuons = 0 ;

     leaves.met = -1. ;
     leaves.metPx = -1. ;
     leaves.metPy = -1. ;

     for ( int i=0; i< MAXSC ;  i++ ) {
         leaves.cscRho[i]  = -1. ;
         leaves.cscdPhi[i] = 3.16 ;
         leaves.dtdEta[i]  = 9. ;
         leaves.dtdPhi[i]  = 9. ;
         leaves.xtalEta[i] = 9. ;
         leaves.xtalPhi[i] = 9. ;
         leaves.xtalE[i]   = 0. ;
         leaves.xtalSigma[i]  = 0. ;
         leaves.xtalSwissX[i] = -9. ;
         leaves.xtalChi2[i]   = -9. ;
         leaves.xtal_x[i] = 199.;
         leaves.xtal_y[i] = 199.;
         leaves.xtal_z[i] = 399.;
         leaves.xtal_t[i] = -199.;
         leaves.bc_E[i]    = 0. ;
         leaves.bc_sMaj[i] = 0. ;
         leaves.bc_sMin[i] = 0. ;
         leaves.bc_nXtals[i] = 0;
         leaves.bc_nBC[i]    = 0;
     }
 
     for ( int i=0; i< MAXCSC ; i++ ) {
         leaves.cscX[i] = 0. ;
         leaves.cscY[i] = 0. ;
         leaves.cscZ[i] = 0. ;
         leaves.cscdX[i] = 0. ;
         leaves.cscdY[i] = 0. ;
         leaves.cscdZ[i] = 0. ;
         leaves.csc_dPhi[i] = 3.16 ;
     }
     for ( int i=0; i< MAXDT ; i++ ) {
         leaves.dtX[i] = 0. ;
         leaves.dtY[i] = 0. ;
         leaves.dtZ[i] = 0. ;
         leaves.dtdX[i] = 0. ;
         leaves.dtdY[i] = 0. ;
         leaves.dtdZ[i] = 0. ;
         leaves.dt_dR[i] = 99. ;
     }
     for ( int i=0; i< MAXCMU ; i++ ) {
         leaves.mu_nDT[i] = 0 ; 
         leaves.mu_nCSC[i] = 0 ; 
         leaves.muPx[i] = 0. ;
         leaves.muPy[i] = 0. ;
         leaves.muPz[i] = 0. ;
         leaves.muE[i]  = 0. ;
     }
}

//define this as a plug-in
DEFINE_FWK_MODULE(CosmicAnalysis);
