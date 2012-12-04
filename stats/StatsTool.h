#ifndef StatsTool_H
#define StatsTool_H

#include "TObject.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <TMath.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <THStack.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TRandom.h>
#include <TLeaf.h>
#include <TString.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFrame.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TGraphAsymmErrors.h>
#include "TLorentzVector.h"

#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooRandom.h"
#include "RooPlot.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

#include <algorithm>

#include "AnaInput.h"
#include "hDraw.h"
#include "StatsLib.h"

class StatsTool : public TObject {

public:

   StatsTool( string datacardfile = "DataCard.txt");     
   ~StatsTool();     
   
   void Test() ;
   void Demo() ;
   void EffCalculator( string mcTag, int norm_id ) ;
   void DPCLsTest( string mcTag, int norm_id, double pValue_B = -1 ) ;
   void MyCLsTest() ;
   void Toys() ;
   void DrawLimit() ;
   void SpecialPlot() ;
   
private:

   AnaInput*     Input;
   hDraw*        h_draw ;
   StatsLib*     statslib ;

   string rfolder  ;
   string hfolder  ;
   string plotType ;
   string hFileName ;
   string logFileName ;
   double searchLimit ;

   vector<double> timeScope ;
   vector<string> observables ;
   double lumi ;
   vector<double> sgXsec ;
   vector<double> nGen ;
   //int ProcessEvents ;
   string dataFileNames;
   string mcFileNames;

  
   //Input histogram file name

   //ClassDef(StatsTool, 1);
};

//#if !defined(__CINT__)
//    ClassImp(StatsTool);
#endif

