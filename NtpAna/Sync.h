#ifndef Sync_H
#define Sync_H

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
#include "MathFunctions.h"

#define MAXSC 50
#define MAXC 200
#define MAXXTALINC 25 
#define MAXOBJ 10
#define MAXPHO 40
#define MAXVTX 40

//class Sync : public TObject {
class Sync {

public:

   Sync( string datacardfile = "DataCard.txt");     
   ~Sync();     
   
   void Init() ;
   void ReadTree() ;
   bool GammaJetsBackground( TLorentzVector gP4, vector<TLorentzVector> jP4s ) ; 

private:

   AnaInput*       Input ;
   Selection*      select ;

   string hfolder ;
   string plotType ;
   int ProcessEvents ;

   string debugStr ;

   vector<TLorentzVector> jetV ;
   vector<TLorentzVector> eleV ;
   vector<float> eIdV ;

   vector<double> vtxCuts ;
   vector<double> photonCuts ;
   vector<double> photonIso ;
   vector<double> electronCuts ;
   vector<double> jetCuts ;
   int selectBackground ;
   int split ;

   TTree* tr ;
   float vtxX[MAXVTX],    vtxY[MAXVTX],  vtxZ[MAXVTX],   vtxChi2[MAXVTX], vtxNdof[MAXVTX];
   float jetPx[MAXOBJ],   jetPy[MAXOBJ], jetPz[MAXOBJ],  jetE[MAXOBJ] ;
   float jetNDau[MAXOBJ], jetCM[MAXOBJ], jetCEF[MAXOBJ], jetNHF[MAXOBJ], jetNEF[MAXOBJ];
   float phoPx[MAXPHO],      phoPy[MAXPHO],      phoPz[MAXPHO],     phoE[MAXPHO], phoTime[MAXPHO] ;
   float phoEcalIso[MAXPHO], phoHcalIso[MAXPHO], phoTrkIso[MAXPHO], phoHovE[MAXPHO], phoSmin[MAXPHO], phoSmaj[MAXPHO] ;
   float elePx[MAXOBJ], elePy[MAXOBJ], elePz[MAXOBJ], eleE[MAXOBJ] ;
   float eleEcalIso[MAXOBJ], eleHcalIso[MAXOBJ], eleTrkIso[MAXOBJ], eleNLostHits[MAXOBJ] ;
   float metPx, metPy, metE ;
   int   nJets, nPhotons, nElectrons, nVertices, eventId, trgCut ;

   int counter[10] ;

   //ClassDef(Sync, 1);
};

//#if !defined(__CINT__)
//    ClassImp(Sync);
#endif

