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
#include "Histograms.h"

#define MAXSC 50
#define MAXC 200
#define MAXXTALINC 25 
#define MAXVTX 40
#define MAXHCALRECHITS 100
#define MAXCALOTOWERS 100
#define MAXMU 20
#define MAXTOWERSINTPGSUMMARY 100
#define MAXL1OBJS 100
#define MAXJET 10
#define MAXELE 10
#define MAXPHO 10
#define MAXOBJ 50


class PhotonAna : public TObject {

private:

   AnaInput*       Input;

   string hfolder ;
   string plotType ;

public:

   PhotonAna();     
   ~PhotonAna();     
   
   void test();
   void ReadTree();

   //ClassDef(PhotonAna, 1);

};

//#if !defined(__CINT__)
//    ClassImp(PhotonAna);
#endif

