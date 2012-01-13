#ifndef PhotonAna_H
#define PhotonAna_H

#include "TObject.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TString.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include "TLorentzVector.h"

#include "AnaInput.h"
#include "Selection.h"
#include "Histograms.h"
#include "MathFunctions.h"
#include "timeVsAmpliCorrector.h"

#define MAXSC 50
#define MAXC 200
#define MAXXTALINC 25 
#define MAXOBJ 10
#define MAXVTX 40

class PhotonAna : public TObject {

public:

   PhotonAna( string datacardfile = "DataCard.txt");     
   ~PhotonAna();     
   
   void ReadTree();
   
   vector<double> BinningFit( TH2D* h2, string hName, int xbinMin, int xbinMax, bool debugPlots = false );
   void BinningFitScan( TH2D* h2, vector<double>& xV, vector<double>& yV, vector<double>& yWidthV, vector<double>& yErrV, 
                        bool debugPlots = false,  int rbin = 1, int startBin = 1, int finalBin = -1) ;

   // drawOpt = 0 : default , 1:SetLogx , 2:SetLogy 
   void ScalarPlotter( TH2D* h2, TString hname, double yMin, double yMax, int rbin = 1, int drawOpt = 0, bool debugPlots = false );
   void ScalarPlotList() ;

   double DeltaR( float eta1, float phi1, float eta2, float phi2 ) ;

   bool GammaJetsBackground( TLorentzVector gP4, vector<TLorentzVector> jP4s ) ; 

private:

   AnaInput*     Input;
   Selection*    select;

   string hfolder ;
   string plotType ;
   int ProcessEvents ;

   hJetTime hJets ;
   string debugStr ;

   vector<objID> jetV ;
   vector<objID> eleV ;
   vector<objID> phoV ;
   vector<float> eIdV ;

   vector<double> vtxCuts ;
   vector<double> photonCuts ;
   vector<double> electronCuts ;
   vector<double> jetCuts ;
   vector<double> BasicClusterCuts ;
   vector<double> XtalCuts ;
   int selectBackground ;
   int split ;
   int doTimeCorrection ;

   timeCorrector theTimeCorrector_;

   float jetPx[MAXOBJ], jetPy[MAXOBJ], jetPz[MAXOBJ], jetE[MAXOBJ] ;
   float phoPx[MAXOBJ], phoPy[MAXOBJ], phoPz[MAXOBJ], phoE[MAXOBJ] ;
   float elePx[MAXOBJ], elePy[MAXOBJ], elePz[MAXOBJ], eleE[MAXOBJ] ;
   float CPIdx[MAXC], clusterTime[MAXC], clusterEnergy[MAXC];
   float xtalInBCEnergy[MAXC][MAXXTALINC], xtalInBCTime[MAXC][MAXXTALINC], xtalInBCTimeErr[MAXC][MAXXTALINC];
   float xtalInBCEta[MAXC][MAXXTALINC],    xtalInBCPhi[MAXC][MAXXTALINC],  xtalADC[MAXC][MAXXTALINC] ;
   float xtalChi2[MAXC][MAXXTALINC], xtalOutTimeChi2[MAXC][MAXXTALINC];
   float clusterPhi[MAXC], clusterEta[MAXC] ;
   float vtxX[MAXVTX], vtxY[MAXVTX], vtxZ[MAXVTX], vtxChi2[MAXVTX], vtxNdof[MAXVTX];
   int   nXtalInBC[MAXC], clusterXtals[MAXC] ;
   float metPx, metPy, metE ;
   int   nJets, nPhotons, nElectrons, nVertices, eventId, nClusters ;


   //ClassDef(PhotonAna, 1);
};

//#if !defined(__CINT__)
//    ClassImp(PhotonAna);
#endif

