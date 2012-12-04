#include "hDraw.h"
//#include "MathFunctions.h"
hDraw::hDraw( string datacardfile, string canvasName ){

  Input  = new AnaInput( datacardfile );

  Input->GetParameters("PlotType",      &plotType ) ;
  Input->GetParameters("Path",          &hfolder ) ;
  defaultPath = hfolder ;

  gSystem->mkdir( hfolder.c_str() );

  string cname[4] ;
  for ( int i=1; i<= 4; i++) {
      char cnameBuffer[15] ; 
      sprintf(cnameBuffer, "c%d_%s", i, canvasName.c_str() ) ;
      cname[i-1] = cnameBuffer ;
  }

  c1  = new TCanvas( cname[0].c_str(),"", 800, 600);
  c2  = new TCanvas( cname[1].c_str(),"", 800, 600);
  c3  = new TCanvas( cname[2].c_str(),"", 800, 600);
  c4  = new TCanvas( cname[3].c_str(),"", 800, 700);

  func1 = NULL ;

}

hDraw::~hDraw(){

  delete c1 ;
  delete c2 ;
  delete c3 ;
  delete Input ;
  cout<<" Draw ! "<<endl ;

}


void hDraw::SetNewPath( string newpath ) {

  if ( newpath.size() > 0 ) {
     hfolder = newpath ;
     //gSystem->cd( "../" ) ;
     printf(" Set New Path to %s \n ", hfolder.c_str() ) ;
  } else {
     hfolder = defaultPath ;
     printf(" Set to Default Path => %s \n ", hfolder.c_str() ) ;
  }

}

void hDraw::Draw( TH1D* h1_, string plotName, string xTitle, string yTitle, string logY, float statY, int color, double scale_, TLegend* leg ) {

      TH1D* h1 = (TH1D*)h1_->Clone() ;
      h1->Scale( scale_ ); 

      c1->Clear();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->SetLogy(0);

      gStyle->SetOptStat("eoumir");
      if ( strncasecmp( "logY", logY.c_str(), logY.size() ) ==0 && logY.size() > 0 ) c1->SetLogy() ;

      gStyle->SetStatY( statY  );
      gStyle->SetStatTextColor( color );
      gStyle->SetStatTextColor( color );

      h1->GetXaxis()->SetTitle( xTitle.c_str() );
      h1->GetYaxis()->SetTitle( yTitle.c_str() );
      h1->SetLineColor( color ) ;
      h1_->SetLineColor( color ) ;

      c1->cd();
      h1->Draw() ;
      c1->Update();

      if ( leg != NULL ) leg->Draw("sames") ;

      if ( plotName.size() > 0 ) {
         TString plotname1 = hfolder + plotName + "."+plotType ;
         c1->Print( plotname1 );
         //if ( func1 != NULL ) func1 = NULL ;
      }

}

void hDraw::DrawAppend( TH1D* h1_, string plotName, float statY, int color, double scale_, TLegend* leg ) {

      TH1D* h1 = (TH1D*)h1_->Clone() ;
      h1->Scale( scale_ ); 

      gStyle->SetStatY( statY  );
      gStyle->SetStatTextColor( color );

      h1->SetLineColor( color ) ;
      h1_->SetLineColor( color ) ;

      c1->cd();
      h1->DrawCopy("sames") ;
      c1->Update();

      if ( leg != NULL ) leg->Draw("sames") ;

      if ( plotName.size() > 0 ) {
         TString plotname1 = hfolder + plotName + "."+plotType ;
         c1->Print( plotname1 );
         //if ( func1 != NULL ) func1 = NULL ;
      }
}

void hDraw::Draw2D( TH2D* h2, string plotName, string xTitle, string yTitle, string logZ, int nColor, float titleFontSize, float statFontSize ) {

      TCanvas* c4  = new TCanvas("c4","", 800, 600 );
      c4->SetFillColor(10);
      c4->SetFillColor(10);
      //c4->SetLogy(0);

      if ( strncasecmp( "logZ", logZ.c_str(), logZ.size() ) ==0  && logZ.size() > 0 ) gPad->SetLogz(1) ;

      SetHistoAtt( h2 ) ;
      gStyle->SetPalette(1);
      gStyle->SetStatFontSize( statFontSize ) ;
      gStyle->SetTitleFontSize( titleFontSize ) ;
      gStyle->SetNumberContours( nColor );
      gStyle->SetStatX(0.9);

      h2->GetXaxis()->SetTitle( xTitle.c_str() );
      h2->GetYaxis()->SetTitle( yTitle.c_str() );

      c4->cd() ;
      h2->Draw("COLZ") ;
      c4->Update();

      if ( plotName.size() > 0 ) {
         TString plotname1 = hfolder + plotName + "."+plotType ;
         c4->Print( plotname1 );
      }
      delete c4 ;
}

void hDraw::CreateNxM( string plotName , int N, int M ) {

      c2->Clear();
      c2->SetFillColor(10);
      c2->SetFillColor(10);
      c2->Divide(N, M) ;

      gStyle->SetOptStat("eroumi");

      plotname2 = hfolder + plotName + "."+plotType ;
}

void hDraw::DrawNxM( int id, TH1D* h1, string xTitle, string yTitle, string logY, int color, bool close  ) {

      c2->cd( id );

      gPad->SetLogy(0) ;
      if ( strncasecmp( "logY", logY.c_str(), logY.size() ) ==0  && logY.size() > 0 ) gPad->SetLogy(1) ;

      //SetHistoAtt( h1 ) ;
      gStyle->SetStatTextColor( color );
      h1->SetLineColor( color ) ;

      h1->GetXaxis()->SetTitle( xTitle.c_str() );
      h1->GetYaxis()->SetTitle( yTitle.c_str() );

      h1->Draw() ;
      c2->Update();

      if ( close ) {
         c2->Print( plotname2 );
      }
}

void hDraw::DrawNxM( int id, TH2D* h2, string xTitle, string yTitle, string logZ, int nColor, float titleFontSize, float statFontSize, bool close  ) {

      c2->cd( id );

      gPad->SetLogy(0) ;
      if ( strncasecmp( "logZ", logZ.c_str(), logZ.size() ) ==0  && logZ.size() > 0 ) gPad->SetLogz(1) ;

      SetHistoAtt( h2 ) ;
      gStyle->SetStatFontSize( statFontSize ) ;
      gStyle->SetTitleFontSize( titleFontSize ) ;
      gStyle->SetNumberContours( nColor );

      h2->GetXaxis()->SetTitle( xTitle.c_str() );
      h2->GetYaxis()->SetTitle( yTitle.c_str() );

      h2->Draw("COLZ") ;
      c2->Update();

      if ( close ) {
         c2->Print( plotname2 );
      }
}

void hDraw::DrawGraph( TGraph* gr , string graphName,  string xTitle, string yTitle, double yMax, double yMin ) {

    // Draw the result
     c4->SetFillColor(10);
     c4->SetFillColor(10);
     gPad->SetGridx();
     gPad->SetGridy();

     c4->cd();
     gStyle->SetOptStat(kTRUE);
     gStyle->SetOptFit(111);

     if ( yMax <  9999999 ) gr->SetMaximum( yMax );
     if ( yMin > -9999999 ) gr->SetMinimum( yMin);
     gr->SetMarkerStyle(20);
     gr->SetMarkerSize(1);
     gr->SetLineWidth(2);
     gStyle->SetTitleFontSize(0.04) ;
     gr->GetXaxis()->SetTitleOffset(1.34);
     gr->GetYaxis()->SetTitleOffset(1.39);
     gr->GetXaxis()->SetTitleFont(42);
     gr->GetYaxis()->SetTitleFont(42);
     gr->GetXaxis()->SetTitleSize(0.04);
     gr->GetXaxis()->SetTitle( xTitle.c_str() ) ;
     gr->GetYaxis()->SetTitle( yTitle.c_str() ) ;
     gr->Draw("ALP") ;

     c4->Update();
     TString plotname = hfolder + graphName + "." + plotType ;
     c4->Print( plotname );

}

void hDraw::FitNDraw2( TH1D* h1, TH1D* h2, string plotName, string xTitle, string yTitle, string logY, float statY, TLegend* leg ) {

      // Plot the histogram
      c1->Clear();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->SetLogy(0);

      //gStyle->SetOptStat("emio");
      if ( strncasecmp( "logY", logY.c_str(), logY.size() ) ==0 && logY.size() > 0 ) c1->SetLogy() ;

      gStyle->SetStatY( statY  );
      gStyle->SetStatTextColor( 2 );

      h1->GetXaxis()->SetTitle( xTitle.c_str() );
      h1->GetYaxis()->SetTitle( yTitle.c_str() );
      h1->SetLineColor( 2 ) ;

      c1->cd();
      h1->Draw() ;
      c1->Update();

      // Define the fitting function and fit
      h1->Fit( func1, "R" );
      func1->SetLineColor(2) ;
      func1->Draw("same") ;
      c1->Update();

      gStyle->SetStatY( statY - 0.2 );
      gStyle->SetStatTextColor( 4 );
      h2->SetLineColor( 4 ) ;

      h2->DrawCopy("sames") ;
      c1->Update();
      //func1->SetLineColor(4) ;
      //h2->Fit( func1, "R" ) ;
      //func1->Draw("same") ;
      //c1->Update();

      if ( leg != NULL ) leg->Draw("sames") ;

      if ( plotName.size() > 0 ) {
         TString plotname1 = hfolder + plotName + "."+plotType ;
         c1->Print( plotname1 );
      }
}

void hDraw::FitNDraw( TH1D* h1, string plotName, string xTitle, string yTitle, string logY, float statY, int color, TLegend* leg ) {

      // Plot the histogram
      c1->Clear();
      c1->SetFillColor(10);
      c1->SetFillColor(10);
      c1->SetLogy(0);

      //gStyle->SetOptStat("emio");
      if ( strncasecmp( "logY", logY.c_str(), logY.size() ) ==0 && logY.size() > 0 ) c1->SetLogy() ;

      gStyle->SetStatY( statY  );
      gStyle->SetStatTextColor( color );

      h1->GetXaxis()->SetTitle( xTitle.c_str() );
      h1->GetYaxis()->SetTitle( yTitle.c_str() );
      h1->SetLineColor( color ) ;

      c1->cd();
      h1->Draw() ;
      c1->Update();

      // Define the fitting function and fit
      h1->Fit( func1, "R" );
      func1->Draw("same") ;
      c1->Update();

      if ( leg != NULL ) leg->Draw("sames") ;

      if ( plotName.size() > 0 ) {
         TString plotname1 = hfolder + plotName + "."+plotType ;
         c1->Print( plotname1 );
      }
}

void hDraw::FitNDrawAppend( TH1D* h1, string plotName, float statY, int color, TLegend* leg ) {

      gStyle->SetStatY( statY  );
      gStyle->SetStatTextColor( color );

      h1->SetLineColor( color ) ;

      c1->cd();
      h1->DrawCopy("sames") ;
      c1->Update();

      h1->Fit( func1, "R" );
      func1->Draw("same") ;
      c1->Update();

      if ( leg != NULL ) leg->Draw("sames") ;

      if ( plotName.size() > 0 ) {
         TString plotname1 = hfolder + plotName + "."+plotType ;
         c1->Print( plotname1 );
         if ( func1 != NULL ) func1 = NULL ;
      }
}

void hDraw::SetFitParameter( string fitFunc_, double fitMin_, double fitMax_ , int nFitPara_,  double* initVals, int color ) {

   fitMin = fitMin_ ;
   fitMax = fitMax_ ;
   nFitPara = nFitPara_ ;

   gStyle->SetOptFit(111);

   if ( strncasecmp( "Gaus", fitFunc_.c_str(), fitFunc_.size() ) == 0  )  {
      func1 = new TF1("func1", &hDraw::fitGS, fitMin, fitMax, nFitPara );  
      func1->SetParameter( 0, initVals[0] );
      func1->SetParameter( 1, initVals[1] );
      func1->SetParameter( 2, initVals[2] );
      func1->SetLineColor( color );
   }
   if ( strncasecmp( "Exp", fitFunc_.c_str(), fitFunc_.size() ) == 0  )  {
      func1 = new TF1("func1", &hDraw::fExp , fitMin, fitMax, nFitPara );
      func1->SetParameter( 0, initVals[0] );
      func1->SetParameter( 1, initVals[1] );
      func1->SetLineColor( color );
   }

}

void hDraw::SetFitParameter( string fitFunc_, TH1D* h1, double fitMin_, double fitMax_ , int nFitPara_, int color ) {

   if ( fitMin_ != fitMax_ ) fitMin = fitMin_ ;
   if ( fitMin_ != fitMax_ ) fitMax = fitMax_ ;
   nFitPara = nFitPara_ ;

   gStyle->SetOptFit(111);

   if ( strncasecmp( "Gaus", fitFunc_.c_str(), fitFunc_.size() ) == 0  )  {
      double maxv = h1->GetMaximum() ;
      double mean = h1->GetMean() ;
      double rms  = h1->GetRMS() ;
      if ( fitMin_ == fitMax_ ) fitMin = mean - (2.*rms) ;
      if ( fitMin_ == fitMax_ ) fitMax = mean + (2.*rms) ;

      func1 = new TF1("func1", &hDraw::fitGS, fitMin, fitMax, nFitPara );  
      func1->SetParameter( 0, maxv );
      func1->SetParameter( 1, mean );
      func1->SetParameter( 2, rms  );
      func1->SetLineColor( color );
   }

}

// Drawing Efficiency Plot from two histograms    
void hDraw::EffPlot( TH1D* hCut, TH1D* hAll, string xlable, double minBinContent, int beginBin, int endBin, string graphName ) {

   endBin = ( endBin == -1 ) ? hAll->GetNbinsX() : endBin ;

   double ba = 0 ;
   double bc = 0 ;
   double x  = 0 ;
   double rbin  = 0 ;
   double rbin_ = 0 ;
   vector<double> baV ;
   vector<double> bcV ;
   vector<double> xV ;
   vector<double> xW ;
   vector<double> errL;
   vector<double> errH;
   double bW = hCut->GetBinWidth(1) ;
   cout<<" bin width = " << bW <<" end bin = "<< endBin << endl ;
   // Accumulate bin contain
   for ( int i= beginBin ; i<= endBin; i++ ) {
       double bc_ = hCut->GetBinContent(i) ;
       double ba_ = hAll->GetBinContent(i) ;
       double x_  = hAll->GetBinCenter(i) ;
       // rebin the histogram in order to have consistence statistic for each bin
       if ( ba < minBinContent || ba < bc || rbin < rbin_ ) {
          ba += ba_ ;
          bc += bc_ ;
          x  +=  x_ ;
          rbin++ ;
          if ( ba >= minBinContent && ba >= bc && rbin >= rbin_ ) {
             bcV.push_back( bc ) ;
	     baV.push_back( ba ) ;
	     xV.push_back( x / rbin )  ;
	     xW.push_back( rbin * bW / 2.) ;
	     // sc's method to calculate error
	     pair<double,double> errs = EffError( ba, bc ) ;
	     errH.push_back( errs.first ) ;
	     errL.push_back( errs.second ) ;
	     rbin_ = rbin ;
	     //cout<<" x: "<< x/rbin <<" rb: "<< rbin <<" bW:"<< (rbin*bW) / 2. << "  bc: "<< bc <<"  ba: "<< ba ;
	     //cout<<" eff:"<< bc/ba <<" + "<< errs.first <<" - "<< errs.second << endl ;
	     bc   = 0 ;
	     ba   = 0 ;
	     x    = 0 ;
	     rbin = 0 ;
          }
       }
       else {

          bcV.push_back( bc ) ;
          baV.push_back( ba ) ;
          xV.push_back( x / rbin )  ;
          xW.push_back( rbin * bW / 2.) ;
          // sc's method to calculate error
          pair<double,double> errs = EffError( ba, bc ) ;
          errH.push_back( errs.first ) ;
          errL.push_back( errs.second ) ;
          rbin_ = rbin ;
          //cout<<" x: "<< x/rbin <<" rb: "<< rbin <<" bW:"<< (rbin*bW) / 2. << "  bc: "<< bc <<"  ba: "<< ba ;
          //cout<<" eff:"<< bc/ba <<" + "<< errs.first <<" - "<< errs.second << endl ;
          bc   = 0 ;
          ba   = 0 ;
          x    = 0 ;
          rbin = 0 ;
       }
   }

   cout<<" nbins = "<< xV.size() <<endl ;

   // calculate proper error
   int nbins = (int) xV.size() ;
   TH1D* hTop      = new TH1D("hTop", " ",  nbins,  0, nbins );
   TH1D* hBot      = new TH1D("hBot", " ",  nbins,  0, nbins );
   for ( int i=0; i< nbins ; i++ ) {
       hTop->SetBinContent( i+1 , bcV[i] ) ;
       hBot->SetBinContent( i+1 , baV[i] ) ;
   }
   TGraphAsymmErrors* gRatio = new TGraphAsymmErrors();
   gRatio->BayesDivide( hTop, hBot );

   const int sz = (int) xV.size() ;
   Double_t eYL[sz], eYH[sz], fX[sz], eX[sz];
   Double_t* fY = gRatio->GetY() ;
   Double_t yA[sz], eY_L[sz], eY_H[sz] ;  // for SC's method
   for ( Int_t i=0 ; i< sz ; i++) {
       fX[i]  = xV[i] ;
       eX[i]  = xW[i] ;
       eYL[i] = gRatio->GetErrorYlow(i) ;
       eYH[i] = gRatio->GetErrorYhigh(i) ;
       // for SC's method
       yA[i]   = bcV[i] / baV[i] ;
       eY_L[i] = errL[i] ;
       eY_H[i] = errH[i] ;
   }

   TGraphAsymmErrors* gr  = new TGraphAsymmErrors( sz, fX, fY, eX, eX, eYL, eYH );
   TGraphAsymmErrors* gr1 = new TGraphAsymmErrors( sz, fX, yA, eX, eX, eY_L, eY_H ); // SC's result !

   TCanvas* c0 = new TCanvas("c0","", 800, 700);
   c0->SetFillColor(10);
   c0->SetFillColor(10);
   gPad->SetGridx();
   gPad->SetGridy();
   c0->SetLeftMargin(0.15);
   c0->SetRightMargin(0.12);
   c0->SetTopMargin(0.1);
   c0->SetBottomMargin(0.12);

   c0->cd();
   gStyle->SetOptStat(kTRUE);
   gStyle->SetOptFit(111);

   gr->SetMaximum( 1.1 );
   gr->SetMinimum( 0.0 );
   gr->SetMarkerStyle(20);
   gr->SetMarkerSize(1);
   gr->SetLineWidth(2);
   gStyle->SetTitleFontSize(0.04) ;
   gr->SetTitle(" Efficiency") ;
   gr->GetXaxis()->SetTitleOffset(1.34);
   gr->GetYaxis()->SetTitleOffset(1.39);
   gr->GetXaxis()->SetTitleFont(42);
   gr->GetYaxis()->SetTitleFont(42);
   gr->GetXaxis()->SetTitleSize(0.04);
   gr->GetYaxis()->SetTitleSize(0.04);
   gr->GetXaxis()->SetTitle( xlable.c_str() ) ;
   gr->GetYaxis()->SetTitle(" Efficiency ") ;
   gr->Draw("AP");

   c0->Update();
   TString plotname = hfolder + graphName + "Asym."+plotType ;
   c0->Print( plotname );
   
   c0->cd();
   gr1->SetTitle(" Efficiency ") ;
   gr1->SetMaximum( 1.1 );
   gr1->SetMinimum( 0.0 );
   gr1->SetMarkerColor(4);
   gr1->SetMarkerStyle(22);
   gr1->SetMarkerSize(1);
   gr1->SetLineWidth(2);
   gr1->GetXaxis()->SetTitle( xlable.c_str() ) ;
   gr1->GetYaxis()->SetTitle(" Efficiency ") ;
   gr1->Draw("AP");

   c0->Update();
   plotname = hfolder + graphName + "Asym1."+plotType ;
   c0->Print( plotname );
   
   delete hTop ;
   delete hBot ;
   delete c0 ;
   delete gr ;
   delete gr1 ;
}

// return asymmetry errors <H,L>
pair<double, double> hDraw::EffError( double N_all, double N_pass ) {

    double eff0 = N_pass / N_all ;
    pair<double,double> theErr ;
    if ( eff0 > 1 ) {
       theErr = make_pair( 0 , 0 );
       return theErr ;
    }
    //cout<<" N_All: "<<N_all <<"  N_pass: "<< N_pass << endl;
    int nStep = 1000 ;
    double step = 1. / nStep ;
    //cout<<" step = "<< step <<endl;
    Double_t par[3] = { 1, N_all, N_pass } ;
    Double_t xL[1] = { eff0 } ;
    Double_t xR[1] = { eff0 } ;
    double IntEff = (N_all+ 1) * step * BinomialErr( xR, par ) ;
    //cout<<" Eff 0th : "<< BinomialErr( xR, par ) << endl ;  

    bool skipR = false ;
    bool skipL = false ;
    Double_t pR = 0 ;
    Double_t pL = 0 ;
    while ( IntEff < 0.683)  {
        if ( !skipR && xR[0] < 1 ) {
           xR[0] +=  step ;
           pR = BinomialErr( xR, par ) ;
           IntEff += (pR*step*(N_all+1) ) ;
           //cout<<" ("<< xR[0] <<") --> R : "<< IntEff <<"  pR = "<< pR <<endl ;
        }
        if ( !skipL && xL[0] > 0 ) {
           xL[0] -=  step ;
           pL = BinomialErr( xL, par ) ;
           IntEff += (pL*step*(N_all+1) ) ;
           //cout<<" ("<< xL[0] <<") <-- L : "<< IntEff <<"  pL = "<< pL <<endl;
        }
        //cout<<" ------ "<<endl; 
        skipR = ( pL > pR ) ? true : false ;
        skipL = ( pL < pR ) ? true : false ;
        if ( pL == pR ) {
           skipR = false ;
           skipL = false ;
        }
    }
    //cout<<"  ["<< N_pass/N_all <<"] Prob = "<< IntEff <<endl ; 
    //cout<<"                 - "<< (N_pass/N_all) - xL[0] <<endl ;
    //cout<<"                 + "<< xR[0] - (N_pass/N_all) <<endl ;
    theErr      = make_pair( xR[0] - eff0 , eff0 - xL[0] );

    return theErr ;
}

Double_t hDraw::BinomialErr( Double_t* x, Double_t* par ) {

  Double_t N_all  = par[1] ;
  Double_t N_pass = par[2] ;

  //Double_t Bxy = TMath::Beta( ( N_pass + 1 ), ( N_all - N_pass + 1 ) ) ;
  //cout<< " Beta(x,y) = "<< Bxy <<endl ;
  //Double_t Cnk = pow(x[0], N_pass ) * pow( (1-x[0]) , (N_all - N_pass) ) ;
  //Double_t prob = par[0]*Cnk / ( Bxy * (N_all + 1. ) );
  
  double betaPDF = TMath::BetaDist( x[0],  (N_pass + 1) , (N_all - N_pass + 1) ) ;
  Double_t prob = (par[0] / (N_all + 1.) ) * betaPDF ;
  //cout<<" x = "<< x[0] <<" betaPDF:"<<  betaPDF << " p = "<< prob <<endl ;

  if ( x[0] < 0 || x[0] > 1 ) prob = 0 ;

  return prob ;

}


Double_t hDraw::fExp(Double_t *v, Double_t *par) {
      Double_t arg = v[0] /par[1];

      Double_t fitval = par[0]*TMath::Exp( -1*arg);
      return fitval;
}

Double_t hDraw::fitGS(Double_t *x, Double_t *par) {

     Double_t gs_Value = TMath::Gaus(x[0],par[1],par[2]) ;
     Double_t fitV = par[0]*gs_Value ;
     return fitV;
}

void hDraw::EffProbPlot( double N_all, string graphName ){

  cout<<" N_All = "<< N_all <<endl ;
  TF1* fn1 = new TF1("fn1", hDraw::BinomialErr, 0., 1., 3);

  fn1->SetParameter( 0, 1. ) ;

   TCanvas* c0 = new TCanvas("c0","", 800, 600);
   c0->SetFillColor(10);
   c0->SetFillColor(10);
   gPad->SetGridx();

   for ( int i = 0 ; i <= N_all; i++ ) {

      if ( N_all >= 50 && i%5 != 0 ) continue ;
      c0->cd();
      fn1->SetParameter( 1, N_all ) ;
      fn1->SetParameter( 2, i ) ;
      if ( i == 0 ) fn1->Draw() ;
      fn1->DrawCopy("sames") ;
      c0->Update() ;
      //double intP = fn1->Integral( 0, 1 ) ; 
      //cout<<" Int = "<< intP <<endl ; 
   }
   TString plotname = hfolder + graphName + "."+plotType ;
   c0->Print( plotname );

   delete fn1 ;
   delete c0 ;
}

void hDraw::SetHistoAtt( string axis, float labelSize, float tickLength, float titleSize, float titleOffset ){

     if ( axis == "x" || axis == "X" ) {
        labelSize_x   = labelSize ;
        tickLength_x  = tickLength ;
        titleSize_x   = titleSize ;
        titleOffset_x = titleOffset ;
        if ( labelSize   == (float)0 ) labelSize_x   = 0.05 ;
        if ( tickLength  == (float)0 ) tickLength_x  = 0.03 ;
        if ( titleSize   == (float)0 ) titleSize_x   = 0.04 ;
        if ( titleOffset == (float)0 ) titleOffset_x = 1.0 ;
     }
     if ( axis == "y" || axis == "Y" ) {
        labelSize_y   = labelSize ;
        tickLength_y  = tickLength ;
        titleSize_y   = titleSize ;
        titleOffset_y = titleOffset ;
        if ( labelSize   == (float)0 ) labelSize_y   = 0.05 ;
        if ( tickLength  == (float)0 ) tickLength_y  = 0.03 ;
        if ( titleSize   == (float)0 ) titleSize_y   = 0.04 ;
        if ( titleOffset == (float)0 ) titleOffset_y = 1.0 ;
     }
}

void hDraw::SetHistoInfo( TH1D* h1, TString newTitle, int lineWidth ) {

     if ( newTitle != "" ) h1->SetTitle( newTitle ) ;
     h1->SetLineWidth( lineWidth );

}

void hDraw::SetHistoAtt( TH1D* h1 ){

        h1->GetXaxis()->SetLabelSize( labelSize_x );
        h1->GetXaxis()->SetTickLength( tickLength_x );
        h1->GetXaxis()->SetTitleSize( titleSize_x );
        h1->GetXaxis()->SetTitleOffset( titleOffset_x );
        //h1->GetYaxis()->SetTitleFont(42);

        h1->GetYaxis()->SetLabelSize( labelSize_y );
        h1->GetYaxis()->SetTickLength( tickLength_y );
        h1->GetYaxis()->SetTitleSize( titleSize_y );
        h1->GetYaxis()->SetTitleOffset( titleOffset_y );
}

void hDraw::SetHistoAtt( TH2D* h2 ){

        h2->GetXaxis()->SetLabelSize( labelSize_x );
        h2->GetXaxis()->SetTickLength( tickLength_x );
        h2->GetXaxis()->SetTitleSize( titleSize_x );
        h2->GetXaxis()->SetTitleOffset( titleOffset_x );
        //h1->GetYaxis()->SetTitleFont(42);

        h2->GetYaxis()->SetLabelSize( labelSize_y );
        h2->GetYaxis()->SetTickLength( tickLength_y );
        h2->GetYaxis()->SetTitleSize( titleSize_y );
        h2->GetYaxis()->SetTitleOffset( titleOffset_y );
}

