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

using namespace edm ;
using namespace std ;

// constants, enums and typedefs
// static data member definitions

// constructors and destructor
DPAnalysis::DPAnalysis(const edm::ParameterSet& iConfig){

   //now do what ever initialization is needed
   rootFileName         = iConfig.getUntrackedParameter<string> ("rootFileName");
   trigSource           = iConfig.getParameter<edm::InputTag> ("trigSource");
   pvSource             = iConfig.getParameter<edm::InputTag> ("pvSource");
   beamSpotSource       = iConfig.getParameter<edm::InputTag> ("beamSpotSource");
   muonSource           = iConfig.getParameter<edm::InputTag> ("muonSource");
   electronSource       = iConfig.getParameter<edm::InputTag> ("electronSource");
   photonSource         = iConfig.getParameter<edm::InputTag> ("photonSource");
   metSource            = iConfig.getParameter<edm::InputTag> ("metSource");
   jetSource            = iConfig.getParameter<edm::InputTag> ("jetSource");

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
   triggerPatent        = iConfig.getUntrackedParameter<string> ("triggerName");
   isData               = iConfig.getUntrackedParameter<bool> ("isData");

   gen = new GenStudy( iConfig );

   theFile  = new TFile( rootFileName.c_str(), "RECREATE") ;
   theFile->cd () ;
   theTree  = new TTree ( "DPAnalysis","DPAnalysis" ) ;
   setBranches( theTree, leaves ) ;

   TriggerName = "" ;
   // reset the counter
   for ( int i=0; i< 10 ; i++) counter[i] = 0 ;

}


DPAnalysis::~DPAnalysis()
{
   // do anything here that needs to be done at desctruction time

   delete gen ;
   cout<<"All:"<< counter[0]<<" dumper:"<<counter[1]<<" Vertex:"<< counter[2] <<" photon:"<<counter[3] ;
   cout<<" sMinor:"<< counter[4] <<" BeamHalo:"<< counter[5] <<" Iso:"<<counter[6] <<" GJet:"<<counter[7] ;
   cout<<" Jets:"<< counter[8] << " G100:"<< counter[9] <<endl;
   theFile->cd () ;
   theTree->Write() ; 
   theFile->Close() ;

}

//
// member functions
//

// ------------ method called for each event  ------------
void DPAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

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

   counter[0]++ ;  
   bool passTrigger = TriggerSelection( iEvent ) ;

   if (  passTrigger ) leaves.triggered    = 1 ;
   if ( !passTrigger ) leaves.triggered    = 0 ;

   if ( !isData ) { 
      gen->GetGenEvent( iEvent, leaves );
      //gen->GetGen( iEvent, leaves );
   }
   //if ( !isData ) gen->PrintGenEvent( iEvent );

   bool pass = EventSelection( iEvent ) ;
   
   if ( pass ) theTree->Fill();
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

   bool passEvent = true ;

   if ( passEvent )   counter[1]++ ;  
   bool hasGoodVtx = VertexSelection( recVtxs );
   if ( !hasGoodVtx ) passEvent = false ;
   if ( passEvent )   counter[2]++ ;  

   selectedPhotons.clear() ;
   PhotonSelection( photons, recHitsEB, recHitsEE, selectedPhotons ) ;
   if ( selectedPhotons.size() < (size_t)photonCuts[5] )  passEvent = false ;
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
       counter[5]++ ;  
     } else {
       passEvent = false ;
     }
   } else {
       counter[5]++ ;
   }
 
   //IsoPhotonSelection( selectedPhotons ) ;
   //if ( selectedPhotons.size() < photonCuts[5] )  passEvent = false ;
   if ( passEvent )   counter[6]++ ;  

   selectedJets.clear() ;
   JetSelection( jets, selectedPhotons, selectedJets );
   //bool isGammaJets = GammaJetVeto( selectedPhotons, selectedJets ) ;
   //if ( isGammaJets ) passEvent = false ;
   if ( passEvent )   counter[7]++ ;   
   if ( selectedJets.size() < jetCuts[3] )   passEvent = false ;
   if ( passEvent )   counter[8]++ ;  

   
   selectedElectrons.clear() ;
   ElectronSelection( electrons, selectedElectrons ) ;

   selectedMuons.clear() ;
   MuonSelection( muons, selectedMuons );
   
   const reco::PFMET pfMet = (*met)[0] ;
   leaves.met   = pfMet.et() ;
   leaves.metPx = pfMet.px() ;
   leaves.metPy = pfMet.py() ;
   if ( pfMet.pt() < metCuts[0]  ) passEvent = false;

   return passEvent ;
}

bool DPAnalysis::TriggerSelection( const edm::Event& iEvent ) {

   bool pass =false ;
   Handle<edm::TriggerResults> triggers;
   iEvent.getByLabel( trigSource, triggers );
   const edm::TriggerNames& trgNameList = iEvent.triggerNames( *triggers );

   if ( TriggerName.size() == 0 ) { 
      for ( size_t i =0 ; i < trgNameList.size(); i++ ) {
          string tName  = trgNameList.triggerName( i );
          if ( strncmp( tName.c_str(), triggerPatent.c_str(), triggerPatent.size() ) ==0 ) {
             TriggerName = tName ; 
             cout<<" Trigger Found : "<< tName <<" accepted ? "<< triggers->accept(i) <<endl;
             pass = ( triggers->accept(i) == 1 ) ? true : false ;
          }
      }
   } else {

          int trgIndex  = trgNameList.triggerIndex( TriggerName );
          pass = ( triggers->accept(trgIndex) == 1 ) ? true : false ;
          //if ( pass ) cout<<" Trigger Found : "<< TriggerName  << endl; 
   } 

   return pass ;
}
/*
bool DPAnalysis::TriggerSelection( const edm::Event& iEvent ) {

   Handle<edm::TriggerResults> triggers;
   iEvent.getByLabel( trigSource, triggers );

   const edm::TriggerNames& trgNames = iEvent.triggerNames( *triggers );

   int trgIndex  = trgNames.triggerIndex(triggerPatent);
   int trgResult = 0;
   if ( trgIndex == (int)(trgNames.size()) ) {
          cout<<" NO Matched Trigger -- Turn TriggerSelection Off "<<endl;
          cout<<"" <<endl;
          trgResult = 1 ;
   } else {
          trgResult = triggers->accept(trgIndex);
          leaves.triggered = trgResult ;
   }

   bool pass =  ( trgResult == 1 ) ? true : false ;
   return pass ;
}
*/
int DPAnalysis::TriggerSelection( const edm::Event& iEvent, int cutVal, string str_head, string str_body ) {

   Handle<edm::TriggerResults> triggers;
   iEvent.getByLabel( trigSource, triggers );

   cout<<" ** Trigger size = "<< triggers->size() <<endl;
   const edm::TriggerNames& trgNames = iEvent.triggerNames( *triggers );

   int cutSize = ( cutVal > 99 ) ? 3 : 2 ;
   int trgWordSize = str_head.size() + str_body.size() + cutSize + 3 ; //  3 -> _vX

   int trgResult = -1 ;
   int trgResultSum = 0 ;
   for ( size_t i =0 ; i < trgNames.size(); i++ ) {
       string tName  = trgNames.triggerName( i );
       if  ( (int) tName.size() < trgWordSize ) continue ;
       string strhead = tName.substr(0, str_head.size() ) ;
       string strPt   = tName.substr( str_head.size() , cutSize ) ;
       int cutPt      =  atoi ( strPt.c_str() ) ;
       string strbody = tName.substr( str_head.size() + cutSize , str_body.size() ) ;
       string strend = tName.substr( trgWordSize-3 ,3) ;
       if ( strhead == str_head && strbody == str_body && cutPt >= cutVal) {
          int trgIndex  = trgNames.triggerIndex(tName);
          int accept    = triggers->accept(trgIndex);
          if ( accept == 1 ) trgResult = cutPt ; 
          if ( accept == 0 ) trgResult = 0 ; 
          
          if ( accept == 1 ) trgResultSum += cutPt ; 
          cout<<" -> "<< strhead << cutPt << strbody << strend <<" = "<< trgResult <<"  sum = "<<trgResultSum <<endl;
       }
   } 
   return trgResult ;
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
       if ( strncmp( tName.c_str(), triggerPatent.c_str(), triggerPatent.size() ) ==0 ) {
          TriggerName = tName ;
       }
       //string triggered = triggers->accept(i) ? "Yes" : "No" ;
       //cout<<" path("<<i<<") accepted ? "<< triggered ;
   }
}

bool DPAnalysis::VertexSelection( Handle<reco::VertexCollection> vtx ) {

    int thisVertex=0;
    bool hasGoodVertex = true ;

    for(reco::VertexCollection::const_iterator v=vtx->begin();  v!=vtx->end() ; v++){

       if ( fabs(v->z()) >= vtxCuts[0] ) continue ; 
       if (   v->ndof()   < vtxCuts[1] ) continue ;
       double d0 = sqrt( ( v->x()*v->x() ) + ( v->y()*v->y() ) );
       if ( d0 >= vtxCuts[2] ) continue ;

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
 
     if ( thisVertex < 1 )   hasGoodVertex = false ;
     return hasGoodVertex ;
}

bool DPAnalysis::PhotonSelection( Handle<reco::PhotonCollection> photons, Handle<EcalRecHitCollection> recHitsEB, Handle<EcalRecHitCollection> recHitsEE, vector<const reco::Photon*>& selectedPhotons ) {

   int k= 0 ;
   for(reco::PhotonCollection::const_iterator it = photons->begin(); it != photons->end(); it++) {
       // fiducial cuts
       if ( k >= MAXPHO ) break ;
       if ( it->pt() < photonCuts[0] || fabs( it->eta() ) > photonCuts[1] ) continue ;
       float hcalIsoRatio = it->hcalTowerSumEtConeDR04() / it->pt() ;
       if  ( ( hcalIsoRatio + it->hadronicOverEm() )*it->energy() >= 6.0 ) continue ;
 
       // S_Minor Cuts from the seed cluster
       reco::CaloClusterPtr SCseed = it->superCluster()->seed() ;
       const EcalRecHitCollection* rechits = ( it->isEB()) ? recHitsEB.product() : recHitsEE.product() ;
       //const EBRecHitCollection* rechits = ( it->isEB()) ? recHitsEB.product() : recHitsEE.product() ;

       Cluster2ndMoments moments = EcalClusterTools::cluster2ndMoments(*SCseed, *rechits);
       float sMin =  moments.sMin  ;
       float sMaj =  moments.sMaj  ;

       // seed Time 
       pair<DetId, float> maxRH = EcalClusterTools::getMaximum( *SCseed, rechits );
       DetId seedCrystalId = maxRH.first;
       EcalRecHitCollection::const_iterator seedRH = rechits->find(seedCrystalId);
       float seedTime = (float)seedRH->time();
      
       if ( sMaj > photonCuts[2] ) continue ;
       if ( sMin <= photonCuts[3] || sMin >= photonCuts[4] ) continue ;

       // Isolation Cuts 
       float ecalSumEt = it->ecalRecHitSumEtConeDR04();
       float hcalSumEt = it->hcalTowerSumEtConeDR04();
       float trkSumPt  = it->trkSumPtSolidConeDR04();  

       bool trkIso  = ( ( trkSumPt / it->pt())     < photonIso[0] ) ; 
       bool ecalIso = ( (ecalSumEt / it->energy()) < photonIso[2] && ecalSumEt < photonIso[1] ) ; 
       bool hcalIso = ( (hcalSumEt / it->energy()) < photonIso[4] && hcalSumEt < photonIso[3] ) ; 
       if ( !trkIso || !ecalIso || !hcalIso ) continue ;

       leaves.phoPx[k] = it->p4().Px() ;
       leaves.phoPy[k] = it->p4().Py() ;
       leaves.phoPz[k] = it->p4().Pz() ;
       leaves.phoE[k]  = it->p4().E() ;
       leaves.phoHoverE[k]  = it->hadronicOverEm() ;
       leaves.phoEcalIso[k] = ecalSumEt ;
       leaves.phoHcalIso[k] = hcalSumEt ;
       leaves.phoTrkIso[k]  = trkSumPt ;

       leaves.sMinPho[k] = sMin ;
       leaves.sMajPho[k] = sMaj ;
       leaves.phoTime[k] = seedTime ;

       selectedPhotons.push_back( &(*it) ) ;
       k++ ;
   }
   leaves.nPhotons = k ;
   //leaves.nPhotons = (int)( selectedPhotons.size() ) ;

   if ( selectedPhotons.size() > 0 )  return true ; 
   else                               return false ;    

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
   float eidx = 11. ;
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
       eidx += 0.1 ;
       //if ( !it->superCluster().isNull() ) sclist.push_back( make_pair(it->superCluster(), eidx ) );
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
