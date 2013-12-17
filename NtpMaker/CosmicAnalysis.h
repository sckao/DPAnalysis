#ifndef CosmicAnalysis_H
#define CosmicAnalysis_H
// -*- C++ -*-
//
// Package:    CosmicAnalysis
// Class:      CosmicAnalysis
// 
/**\class CosmicAnalysis CosmicAnalysis.cc Exotica/CosmicAnalysis/src/CosmicAnalysis.h

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Shih-Chuan Kao
//         Created:  Thu Sep 29 05:26:22 CDT 2011
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"
#include <DataFormats/JetReco/interface/CaloJet.h>
#include "DataFormats/EgammaCandidates/interface/Photon.h"

// AOD Objects
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

// L1 Trigger 
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

// for ECAL cluster
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"

// for CSC Segment
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/METReco/interface/CSCHaloData.h"

// For DT Segment
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"

// Calibration services
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include <algorithm>


#include <TMath.h>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "Ntuple.h"

#include <Math/VectorUtil.h>

using namespace std ;
typedef math::XYZTLorentzVector LorentzVector;

//
// class declaration
//

struct iCluster {
          // from crystals 
          float eta ;
          float phi ;
          float x ;
          float y ;
          float z ;
          float chi2 ;
          float calib ;
          float t ;
          float E ;
          float swissX ;
          // from cluster 
          float sMaj ;
          float sMin ;
          float scE ;
          int nXtals ;
          int nBC ;
} ;

struct iDT {
        float x ;
        float y ;
        float z ;
        float L ;
        float dx ;
        float dy ;
        float dz ;
} ;

class CosmicAnalysis : public edm::EDAnalyzer {

   public:
      explicit CosmicAnalysis(const edm::ParameterSet&);
      ~CosmicAnalysis();

      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      bool EventSelection( const edm::Event& iEvent, const edm::EventSetup& iSetup );

      bool L1TriggerSelection( const edm::Event& iEvent, const edm::EventSetup& iSetup ) ;

      void TriggerTagging( edm::Handle<edm::TriggerResults> triggers, const edm::TriggerNames& trgNameList, int RunID, vector<int>& firedTrig ) ;
      bool TriggerSelection( edm::Handle<edm::TriggerResults> triggers, vector<int> firedTrig ) ;
 
      void PrintTriggers( const edm::Event& iEvent ) ;

      bool CosmicRayMatch( vector<iDT>& dtsegV, vector<iDT>& cscsegV, vector<iCluster>& scV );

      void DTSegmentInfo( edm::Handle<DTRecSegment4DCollection> dtSeg, const edm::EventSetup& iSetup, vector<iDT>& segV );

      void CSCSegmentInfo( edm::Handle<CSCSegmentCollection> cscSeg, const edm::EventSetup& iSetup, vector<iDT>& selectCSC ) ;

      void ClusterInfo( edm::Handle<reco::SuperClusterCollection> scCollection, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE, vector<iCluster>& scV ) ;

      void InitializeCosmicBranches() ;

   private:

      Ctuple leaves ;
      TTree *theTree;
      TFile *theFile;

      // ----------member data ---------------------------
      string rootFileName;
      std::vector<string> triggerPatent ;
      bool L1Select ;
      string l1GTSource ;

      edm::InputTag trigSource;
      edm::InputTag trigEvent;
      edm::InputTag muonSource;
      edm::InputTag photonSource;
      edm::InputTag metSource;
      edm::InputTag jetSource;

      edm::InputTag bcSource;
      edm::InputTag scSource;
      edm::InputTag EBRecHitCollection;
      edm::InputTag EERecHitCollection;
      edm::InputTag DTSegmentTag ;
      edm::InputTag CSCSegmentTag ;
      edm::InputTag staMuons ;

      //edm::InputTag pileupSource ;
      edm::ESHandle<EcalIntercalibConstants> ical;
      edm::ESHandle<EcalADCToGeVConstant> agc;
      edm::ESHandle<EcalLaserDbService> laser;
      edm::ESHandle<CaloGeometry> pGeometry ;
      const CaloGeometry * theGeometry ;
      //edm::ESHandle<GlobalTrackingGeometry> trackingGeometry;

      std::vector<double> muonCuts ;
      std::vector<double> photonCuts ;
      std::vector<double> metCuts ;
      std::vector<double> jetCuts ; 

      //std::vector<const reco::Muon*> selectedMuons ;

      bool passEvent ;
      int runID_ ;

//      timeCorrector theTimeCorrector_;
      edm::Timestamp eventTime ;

      std::vector<int> firedTrig ;
      int targetTrig ;
      //std::vector<int> firedTrigID ;
      ///string TriggerName ;
      bool passL1 ;
      bool passHLT ;

      bool debugT ; 

      vector<iCluster> selectSC ;
      vector<iDT>      selectDT ;
      vector<iDT>      selectCSC ;

};

#endif
