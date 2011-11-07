#ifndef Histograms_H
#define Histograms_H

#include <TObject.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TString.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>

using namespace std;

class hObj {

  public:

  hObj(){}
  hObj( string fsuffix ){}
  virtual ~hObj(){}
  virtual void Draw( string hfolder, string plotType ) = 0 ;
  //virtual void Fill1D( TString hName, double weight, double scale = 1. )  = 0;
  //virtual void Fill2D( TString hName, double weight, double scale = 1. )  = 0;
  //virtual vector<TH2D*> Output2D() = 0;

};

class hJetTime : public hObj {

  public:

  hJetTime(){

   nJets    = new TH1D("nJets",    " N of jets ", 11,  0, 11);
   nPhotons = new TH1D("nPhotons", " N of photons ", 11,  0, 11);
   pho0Pt   = new TH1D("pho0Pt",   " leading photon Pt spectrum ", 50,   0., 1000.);
   pho0T    = new TH1D("pho0T",    " leading photon Time spectrum ", 200,   -5., 5.);
   jet0Pt   = new TH1D("jet0Pt",   " leading jet Pt spectrum ", 50,   0., 1000.);
   jet1Pt   = new TH1D("jet1Pt",   " 2nd jet Pt spectrum ", 50,   0., 1000.);
   jet0T    = new TH1D("jet0T",    " leading jet Time spectrum ", 200,   -5., 5.);
   jet1T    = new TH1D("jet1T",    " 2nd jet Time spectrum ", 200,   -5., 5.);
   jet0xT   = new TH1D("jet0xT",   " leading jet Time spectrum from xtals", 200,   -5., 5.);
   jet1xT   = new TH1D("jet1xT",   " 2nd jet Time spectrum from xtals", 200,   -5., 5.);

   jet0Pt_T = new TH2D("jet0Pt_T", " leading jet Pt vs Time", 40,   0., 800., 80,  -4., 4.) ;

  }

  virtual ~hJetTime(){
    delete nJets ;
    delete nPhotons ;
    delete jet0Pt ;
    delete jet1Pt ;
    delete jet0T ;
    delete jet1T ;
    delete pho0Pt ;
    delete pho0T ;
   
    delete jet0xT ;
    delete jet1xT ;
    delete jet0Pt_T ;    

    delete c1 ;
    delete c2 ;
  }

  void Draw( string hfolder, string plotType ) {

      c1 = new TCanvas("c1","", 800, 600);
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->Divide(2,2);

      c1->cd(1);
      nJets->Draw() ;
      c1->Update();

      c1->cd(2);
      nPhotons->Draw() ;
      c1->Update();

      c1->cd(3);
      jet0Pt->SetLineColor(2) ;
      jet0Pt->Draw() ;
      c1->Update();
      jet1Pt->SetLineColor(4) ;
      jet1Pt->Draw("same") ;
      c1->Update();

      c1->cd(4);
      pho0Pt->Draw() ;
      c1->Update();

      TString plotname1 = hfolder + "test1" + "."+plotType ;
      c1->Print( plotname1 );

      gStyle->SetOptStat("erm");
      c2 = new TCanvas("c2","", 800, 600);
      c2->SetFillColor(10);
      c2->SetFillColor(10);
      c2->Divide(2,2);

      c2->cd(1);
      gStyle->SetStatY(0.9);
      jet0T->Draw() ;
      c2->Update();

      gStyle->SetStatY(0.7);
      jet0xT->SetLineColor(2) ;
      jet0xT->Draw("same") ;
      c2->Update();

      c2->cd(2);
      gStyle->SetStatY(0.9);
      jet1T->Draw() ;
      c2->Update();

      gStyle->SetStatY(0.7);
      jet1xT->SetLineColor(2) ;
      jet1xT->Draw("same") ;
      c2->Update();

      c2->cd(3);
      gStyle->SetStatY(0.9);
      pho0T->Draw() ;
      c2->Update();

      c2->cd(4);
      gStyle->SetPalette(1);
      gStyle->SetNumberContours(5);
      gStyle->SetStatX(0.85);
      jet0Pt_T->Draw("COLZ") ;
      c2->Update();

      TString plotname2 = hfolder + "test2" + "."+plotType ;
      c2->Print( plotname2 );
  }

  TH1D* nJets ;
  TH1D* nPhotons ;
  TH1D* jet0Pt ;
  TH1D* jet1Pt ;
  TH1D* pho0Pt ;
  TH1D* jet0T ;
  TH1D* jet1T ;
  TH1D* jet0xT ;
  TH1D* jet1xT ;
  TH1D* pho0T ;

  TH2D* jet0Pt_T ;

  TCanvas* c1 ;
  TCanvas* c2 ;
  //ClassDef(hJetTime, 1);
};
#endif
