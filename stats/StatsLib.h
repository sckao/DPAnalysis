#ifndef StatsLib_H
#define StatsLib_H

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
#include <TRandom3.h>
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

//#include <stdlib.h>
#include <algorithm>
#include <memory>

#include "AnaInput.h"
#include "hDraw.h"


class StatsLib : public TObject {

public:

   StatsLib( string datacardfile = "DataCard.txt");     
   ~StatsLib();     
   
   double Q1( vector<double>& xV , double s, double b, TH1D* hPs, TH1D* hPb )  ;
   double Q1( TH1D* hX , double s, double b, TH1D* hPs, TH1D* hPb )  ;

   void   pdfQ( TH1D* hPdfQ_s, TH1D* hPdfQ_b, int m,  double s, double b, TH1D* hPs, TH1D* hPb, int n = 2000 ) ;
   void   pdfQ( string type, TH1D* hPdfQ, int m,  double s, double b, TH1D* hPs, TH1D* hPb, int n = 2000 ) ;
   void   hPdfQ( TH1D* hPdfQ_s, TH1D* hPdfQ_b, int m,  double s, double b, TH1D* hPs, TH1D* hPb, int n = 10000 ) ;

   double ExpectedLimit( int nObs, double s, double b, TH1D* hPs, TH1D* hPb, double alpha = 0.05, double pVal_b = 0.5, int nToys = 1500, string subfolder = "" ) ;
   double ExpectedLimit1( int nObs, double s, double b, TH1D* hPs, TH1D* hPb, double alpha = 0.05, double pVal_b = 0.5, int nToys = 1500, string subfolder = "" ) ;
   vector<double> ExpectedLimit2( int nObs, double s, double b, TH1D* hPs, TH1D* hPb, double alpha = 0.05, int nToys = 1500, string subfolder = "" ) ;

   bool   CLPropagator( double expS, double expB, double expCL, double CL0, double& CL, double& step, vector<double> CLsCollV, vector<double> sNormV ) ; 
   double Interpolate( vector<double> xV, vector<double> yV, double interX ) ;

   double CLx( TH1D* hPdfQ,  vector<double>& xV,  double s, double b, TH1D* hPs, TH1D* hPb ) ; 
   double SBModelPValue( TH1D* hPdfQ_sb, TH1D* hPdfQ_b, double pVal_b ) ;
   vector<double> SBModelPValue( TH1D* hPdfQ_sb, TH1D* hPdfQ_b ) ;

   void Sample( TH1D* hSample, string pdfName, int n, double mean, double sigam = 0. ) ;
   vector<double> toyMC( int n, TH1D* hPDF, int si = 0 ) ;
   void toyMC( TH1D* hOut, int n, TH1D* hPDF, bool debug =false ) ;
   double hProb( double x , TH1D* hPDF ) ;

   TH1D* SubHistogram( TH1D* h1, double minX, double maxX ) ;
   TH1D* SmearHistogram( TH1D* h1, double resol ) ;
   vector<double> GausSmear( double mean, double sigma, double bW ) ;

   double KSTest( TH1D* h1, TH1D* h2, int minBin, int maxBin, double hminX = 0 , double binWidth = 1) ;

   TH1D* RebinHistogram( TH1D* h1, string newHistoName, double minBC ) ;

private:

   AnaInput*     Input;
   hDraw*        h_draw ;

   string hfolder  ;
   string plotType ;

   int randomSeed ;
 
   //ClassDef(StatsLib, 1);
};

//#if !defined(__CINT__)
//    ClassImp(StatsLib);
#endif

