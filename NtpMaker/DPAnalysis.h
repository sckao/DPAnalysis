#ifndef DPAnalysis_H
#define DPAnalysis_H
// -*- C++ -*-
//
// Package:    DPAnalysis
// Class:      DPAnalysis
// 
/**\class DPAnalysis DPAnalysis.cc Exotica/DPAnalysis/src/DPAnalysis.h

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

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

// AOD Objects
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"

#include "DataFormats/METReco/interface/BeamHaloSummary.h"

// for ECAL cluster
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include <algorithm>

// Calibration services
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "CalibCalorimetry/EcalTiming/interface/timeVsAmpliCorrector.h"

// PU SummeryInfo
//#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include <TMath.h>
#include "TFile.h"
#include "TTree.h"

#include "Ntuple.h"
#include "GenStudy.h"

#include <Math/VectorUtil.h>

using namespace std ;

//
// class declaration
//
typedef std::pair<reco::SuperClusterRef, float> ParticleSC  ;


class DPAnalysis : public edm::EDAnalyzer {
   public:
      explicit DPAnalysis(const edm::ParameterSet&);
      ~DPAnalysis();

      virtual void analyze(const edm::Event&, const edm::EventSetup&);

      bool EventSelection( const edm::Event& iEvent );

      int  TriggerSelection( const edm::Event& iEvent, int cutVal, string str_head = "HLT_Photon", string str_body = "_CaloIdVL_IsoL" ) ;
      bool TriggerSelection( const edm::Event& iEvent ) ;
 
      bool VertexSelection( edm::Handle<reco::VertexCollection> vtx ) ;

      bool PhotonSelection(  edm::Handle<reco::PhotonCollection> photons, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE, vector<const reco::Photon*>& selectedPhotons ) ;

      pair<double,double> ClusterTime( reco::SuperClusterRef scRef, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE ) ;
      void ClusterTime( reco::SuperClusterRef scRef, edm::Handle<EcalRecHitCollection> recHitsEB, edm::Handle<EcalRecHitCollection> recHitsEE, double& aveTime, double& aveTimeErr ) ;

      bool JetSelection( edm::Handle<reco::PFJetCollection> jets, vector<const reco::Photon*>& selectedPhotons,
                                                                     vector<const reco::PFJet*>& selectedJets ) ;
      bool ElectronSelection( edm::Handle<reco::GsfElectronCollection> electrons, 
                              vector<const reco::GsfElectron*>& selectedElectrons ) ;
      bool MuonSelection( edm::Handle<reco::MuonCollection> muons, vector<const reco::Muon*>& selectedMuons ) ;
      void PrintTriggers( const edm::Event& iEvent ) ;

      bool sMinorSelection( vector<const reco::Photon*>& selectedPhotons,  edm::Handle<EcalRecHitCollection> recHitsEB,     
                              edm:: Handle<EcalRecHitCollection> recHitsEE ) ;

      bool IsoPhotonSelection( vector<const reco::Photon*>& selectedPhotons ) ; 

      bool GammaJetVeto( vector<const reco::Photon*>& selectedPhotons, vector<const reco::PFJet*>& selectedJets) ;

   private:

      Ntuple leaves ;

      TTree *theTree;

      TFile *theFile;
      GenStudy *gen ; 

      // ----------member data ---------------------------
      string rootFileName;
      std::vector<string> triggerPatent ;
      bool isData ;

      edm::InputTag trigSource;
      edm::InputTag pvSource;
      edm::InputTag beamSpotSource;
      edm::InputTag muonSource;
      edm::InputTag electronSource;
      edm::InputTag photonSource;
      edm::InputTag metSource;
      edm::InputTag jetSource;

      edm::InputTag EBRecHitCollection;
      edm::InputTag EERecHitCollection;
      //edm::InputTag pileupSource ;

      edm::ESHandle<EcalIntercalibConstants> ical;
      edm::ESHandle<EcalADCToGeVConstant> agc;
      edm::ESHandle<EcalLaserDbService> laser;
      edm::ESHandle<CaloGeometry> pGeometry ;
      const CaloGeometry * theGeometry ;

      std::vector<double> muonCuts ;
      std::vector<double> electronCuts ;
      std::vector<double> photonCuts ;
      std::vector<double> photonIso ;
      std::vector<double> metCuts ;
      std::vector<double> jetCuts ; 
      std::vector<double> vtxCuts ; 

      std::vector<const reco::PFJet*> selectedJets ;
      std::vector<const reco::GsfElectron*> selectedElectrons ;
      std::vector<const reco::Muon*> selectedMuons ;
      std::vector<const reco::Photon*> selectedPhotons ;

      bool passEvent ;
      int counter[10] ; 
      float sMin_ ;

      timeCorrector theTimeCorrector_;
      edm::Timestamp eventTime ;

      std::vector<int> firedTrig ;
      ///string TriggerName ;

};

#endif
