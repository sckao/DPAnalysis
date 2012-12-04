#include "StatsTool.h"

using namespace RooFit ;
using namespace RooStats ;

StatsTool::StatsTool( string datacardfile ) {

  Input  = new AnaInput( datacardfile );
  h_draw = new hDraw( datacardfile ) ; 
  statslib = new StatsLib( datacardfile ) ;

  //SkipEvents = 0 ;
  Input->GetParameters("PlotType",      &plotType ) ; 
  Input->GetParameters("Path",          &hfolder ) ; 
  Input->GetParameters("RootFiles",     &rfolder ) ; 
  Input->GetParameters("HFileName",     &hFileName ) ;
  Input->GetParameters("LogFile",       &logFileName ) ;

  Input->GetParameters("Lumi",          &lumi ) ;
  Input->GetParameters("N_GenSignal",   &nGen ) ;
  Input->GetParameters("SignalXsec",    &sgXsec ) ;

  Input->GetParameters("SearchLimit",   &searchLimit ) ;
  Input->GetParameters("TimeScope",     &timeScope );
  Input->GetParameters("Observables",   &observables );
 
  Input->GetParameters("TheData",       &dataFileNames );
  Input->GetParameters("TheMC",         &mcFileNames );
  //Input->GetParameters("ProcessEvents", &ProcessEvents ) ; 
}

StatsTool::~StatsTool(){

  delete statslib ;
  delete h_draw ;
  delete Input ;
  cout<<" done ! "<<endl ;

}

// analysis template
void StatsTool::Demo() { 

   RooWorkspace ws("ws");
   ws.Print() ;

   // mystery fucking FactoryWSTool
   ws.factory("Gaussian:normal(x[-10,10],mu[-1,1],sigmab[1])") ;
   ws.defineSet("poi","mu") ;
   ws.defineSet("obs","x") ;

   // set up the strange configuration object
   ModelConfig *modelCfg = new ModelConfig (" Example G(x|mu ,1)") ;
   modelCfg->SetWorkspace ( ws );
   modelCfg->SetPdf ( *(ws.pdf("normal")) );
   modelCfg->SetParametersOfInterest ( *ws.set("poi") );
   modelCfg->SetObservables ( *ws.set("obs") ) ;

   // create a fake dataset 
   RooDataSet *data = ws.pdf("normal")->generate( *ws.set("obs") , 100) ;
   data->Print() ;

   // for convenience later on
   RooRealVar *x  = ws.var("x");

   // set confidence level
   double confidenceLevel = 0.95;

   // example use profile likelihood calculator
   ProfileLikelihoodCalculator plc( *data , *modelCfg );
   plc.SetConfidenceLevel( confidenceLevel );
   LikelihoodInterval* plInt = plc.GetInterval() ;

   // some plots
   TCanvas *c1 = new TCanvas("c1", " ", 600 ,400) ;
   c1->Divide (2) ;

   // plot the data
   c1->cd(1) ;
   RooPlot *frame = x->frame () ;
   data->plotOn( frame );
   data->statOn( frame );
   frame->Draw() ;

   c1->cd(2) ;
   LikelihoodIntervalPlot plot( plInt );
   plot.Draw() ;

   c1->Print( "atest.png" ) ;
   delete c1 ;
}  


// not working yet
void StatsTool::Test() { 

   // get the tree from ntuple
   TTree* tr_data = Input->TreeMap( dataFileNames );
   TTree* tr_mc = Input->TreeMap( mcFileNames );

   // create a dataset with following variables
   RooRealVar bgTime("bgTime","bgTime", -5., 20.) ;
   RooRealVar g1Time("g1Time","g1Time", -5., 20.) ;
   RooRealVar g1Pt("g1Pt","g1Pt", 10, 500.) ;
   RooRealVar met("met","met", 0., 500.) ;
   RooRealVar nJets("nJets","nJets", 0., 500.) ;
   RooArgSet dataVars( bgTime, g1Time, g1Pt, met, nJets, "run2012C" );

   cout<<" print ArgSet "<<endl ;
   dataVars.Print() ;

   RooDataSet *data = new RooDataSet( "data", "data", tr_data, dataVars ) ;
   data->Print() ;
   cout<<" print data "<<endl ;

   // create a MC dataset with following variables
   RooRealVar g1Time_s("g1Time", "g1Time", -5., 20.) ;
   RooRealVar met_s("met", "met", 0., 500.) ;
   RooRealVar nJets_s("nJets","nJets", 0., 500.) ;
   RooArgSet mcVars( g1Time_s, met_s, nJets_s, "gmsb_2000" );
   mcVars.Print() ;
   RooDataSet *sgmc = new RooDataSet( "sgmc", "sgmc", tr_mc, mcVars ) ;
   sgmc->Print() ;
   cout<<" print sgmc "<<endl ;

   // create a background only PDF , RooDataHist is a garbage 
   RooDataHist h_bgTime("h_bgTime", "h_bgTime", bgTime, *data ) ;
   RooHistPdf  bgPDF( "bgPDF", "bgPDF", bgTime, h_bgTime ) ;

   //TH1D* h1_bg = (TH1D*) bgTime.createHistogram("h1_bg") ;

   // create a signal+background PDF
   RooDataHist h_sbTime("h_sbTime", "h_sbTime", g1Time_s, *sgmc ) ;
   RooHistPdf  sbPDF( "sbPDF", "sbPDF", g1Time_s, h_sbTime ) ;

   // create the  mystery fucking Workspace, 
   RooWorkspace ws("ws");
   ws.import( *data ) ;
   ws.Print() ;
   cout<<" print workspace "<<endl ;

   // set up the strange configuration object
   /*
   ModelConfig *modelCfg = new ModelConfig ("test model") ;
   modelCfg->SetWorkspace( ws );
   modelCfg->SetPdf( timePDF );
   modelCfg->SetParametersOfInterest ( *ws.set("g1Time") );
   modelCfg->SetObservables ( *ws.set("g1Time") ) ;

   // set confidence level
   double confidenceLevel = 0.95;

   // example use profile likelihood calculator
   ProfileLikelihoodCalculator plc( *data , *modelCfg );
   plc.SetConfidenceLevel( confidenceLevel );
   LikelihoodInterval* plInt = plc.GetInterval() ;
   */

   // some plots
   TCanvas *c0 = new TCanvas("c1", " ", 800 ,700) ;
   c0->Divide (2,2) ;

   // plot the data
   c0->cd(1) ;
   RooPlot* frame1 = g1Time.frame() ;
   data->plotOn( frame1, Binning(50) ) ;
   data->statOn( frame1 );
   frame1->Draw() ;

   //RooPlot *frame = x->frame () ;
   //data->plotOn( frame );
   //frame->Draw() ;

   c0->cd(2) ;
   RooPlot* frame2 = g1Time_s.frame() ;
   sgmc->plotOn( frame2, Binning(50) ) ;
   sgmc->statOn( frame2 );
   frame2->Draw() ;
   
   c0->cd(3) ;
   RooPlot* frame3 = bgTime.frame() ;
   h_bgTime.plotOn( frame3 ) ;
   h_bgTime.statOn( frame3 );
   frame3->Draw() ;

   c0->Print( "test.png" ) ;
   //delete c0 ;
}

void StatsTool::EffCalculator( string mcTag, int norm_id ) {

   TString Path_fName = rfolder + hFileName + ".root" ;
   TFile* theFile = new TFile( Path_fName, "READ" );
   theFile->cd() ;

   char mchName[20] ;
   sprintf( mchName, "h_sgTimeA_%s",  mcTag.c_str() ) ;
   TH1D* hS = (TH1D*) theFile->Get( mchName ) ;

   double expS0 = hS->Integral() ;
   double Eff   = expS0 / nGen[norm_id] ;
   printf(" %4s Eff : %.5f   from %.1f / %.1f \n", mcTag.c_str() , Eff, expS0, nGen[norm_id] ) ;

   theFile->Close() ;
}

void StatsTool::SpecialPlot() {

   TString Path_fName = rfolder + hFileName + ".root" ;
   TFile* theFile = new TFile( Path_fName, "READ" );
   theFile->cd() ;

   // Open histogram to get signal and background pdf of observable x 
   TH1D* hS_ = (TH1D*) theFile->Get( "h_sgTime_20"  ) ;
   TH1D* hB_ = (TH1D*) theFile->Get( "h_bgTime" ) ;
   //TH1D* hData_ = (TH1D*) theFile->Get( "h_dataTime" ) ;

   // turn the observation scope at certain timing range
   TH1D* hS = statslib->SubHistogram( hS_ ,   -4., 4 ) ;
   TH1D* hB = statslib->SubHistogram( hB_ ,   -4., 4 ) ;
   //TH1D* hData = statslib->SubHistogram( hData_ ,  -4., 4 ) ;

   // smear the signal timing distribution
   TH1D* hS_1 = statslib->SmearHistogram( hS , 0.4 ) ;

   TCanvas* c_a  = new TCanvas( "c_a","", 800, 600);
   c_a->SetFillColor(10);
   c_a->SetFillColor(10);
   //c_a->SetLogy();

   // Define the fitting function and fit
   c_a->cd();

   gStyle->SetOptStat("");
   gStyle->SetOptFit(0011);
   int mode = 2 ;

   if ( mode == 1 ) {
      gStyle->SetStatY( 0.90  );
      gStyle->SetStatTextColor( 2 );
      gStyle->SetStatFontSize( 0.015 ) ;
      c_a->Update();

      hS->SetLineColor(2) ;
      hS->SetLineWidth(2) ;
      hS->Draw() ;
      c_a->Update();

      // fitting results
      TF1*  fc1 = new TF1("fc1", &hDraw::fitGS, -2, 2, 3 );
      fc1->SetParameter( 0, 800. );
      fc1->SetParameter( 1,  0.  );
      fc1->SetParameter( 2,  0.5 );
      fc1->SetLineColor( 2 );
      hS->Fit( "fc1", "R", "sames" );
      c_a->Update();
   }

   // background plot
   if ( mode == 2 ) {
      gStyle->SetStatY( 0.9 );
      gStyle->SetStatTextColor( 8 );
      gStyle->SetStatFontSize( 0.015 ) ;
      c_a->Update();
      hB->SetLineColor( 8 ) ;
      hB->SetLineWidth( 2 ) ;
      hB->Draw() ;
      c_a->Update();


      TF1*  fc3 = new TF1("fc3", &hDraw::fitGS, -2, 2, 3 );
      fc3->SetParameter( 0, 800. );
      fc3->SetParameter( 1,  0.  );
      fc3->SetParameter( 2,  0.5 );
      fc3->SetLineColor( 8 );
      hB->Fit( fc3, "R", "sames" );
      c_a->Update();
   }

   // 2nd plot
   gStyle->SetStatY( 0.7 );
   gStyle->SetStatTextColor( 4 );
   gStyle->SetStatFontSize( 0.015 ) ;
   c_a->Update();
   hS_1->SetLineColor( 4 ) ;
   hS_1->SetLineWidth( 2 ) ;
   hS_1->DrawCopy("SAME") ;
   c_a->Update();


   TF1*  fc2 = new TF1("fc2", &hDraw::fitGS, -2, 2, 3 );
   fc2->SetParameter( 0, 800. );
   fc2->SetParameter( 1,  0.  );
   fc2->SetParameter( 2,  0.5 );
   fc2->SetLineColor( 4 );
   hS_1->Fit( fc2, "R", "sames" );
   c_a->Update();


   TString plotname_a = hfolder +"SmearMC_bg."+plotType ;
   c_a->Print( plotname_a );


}

void StatsTool::DPCLsTest( string mcTag, int norm_id, double pValue_B ) {

   TString Path_fName = rfolder + hFileName + ".root" ;
   TFile* theFile = new TFile( Path_fName, "READ" );
   theFile->cd() ;

   // set the folder to record all the result
   string currentPath = h_draw->CurrentPath() ;
   gSystem->cd( currentPath.c_str() );
   gSystem->mkdir( mcTag.c_str() );
   gSystem->cd( mcTag.c_str() );
   gSystem->mkdir( "center" );
   gSystem->mkdir( "u1sig" );
   gSystem->mkdir( "u2sig" );
   gSystem->mkdir( "d1sig" );
   gSystem->mkdir( "d2sig" );
   gSystem->mkdir( "data" );
   gSystem->mkdir( "test" );
   gSystem->cd( "../../" );

   string newPath = currentPath ;
   newPath +=  mcTag  ;
   newPath.append( "/" ) ;
   h_draw->SetNewPath( newPath ) ;  
   cout<<"set sub path "<< endl ;
   char subPath_c[55] , subPath_u1[55] , subPath_u2[55] , subPath_d1[55] , subPath_d2[55], subPath_0[55], subPath_t[55] ; 
   sprintf( subPath_c,  "%scenter/", newPath.c_str() ) ;
   sprintf( subPath_u1, "%su1sig/",  newPath.c_str() ) ;
   sprintf( subPath_u2, "%su2sig/",  newPath.c_str() ) ;
   sprintf( subPath_d1, "%sd1sig/",  newPath.c_str() ) ;
   sprintf( subPath_d2, "%sd2sig/",  newPath.c_str() ) ;
   sprintf( subPath_0,  "%sdata/",   newPath.c_str() ) ;
   sprintf( subPath_t,  "%stest/",   newPath.c_str() ) ;

   // text file to record final number
   string logfilePath = currentPath ;
   logfilePath += logFileName ;
   FILE* logfile = fopen( logfilePath.c_str() ,"a"); 

   char mchName[20] ;
   sprintf( mchName, "%s_%s",  observables[0].c_str(), mcTag.c_str() ) ;

   // Open histogram to get signal and background pdf of observable x 
   TH1D* hS_ = (TH1D*) theFile->Get( mchName ) ;
   TH1D* hB_ = (TH1D*) theFile->Get( observables[1].c_str() ) ;
   TH1D* hData_ = (TH1D*) theFile->Get( observables[2].c_str() ) ;

   // Smear MC samples
   /*
   string mchName2 ;
   mchName2 += "_s" ;
   TH1D* hS_1 = statslib->SmearHistogram( hS_ , 0.4 ) ;
   TH1D* hS_2 = statslib->RebinHistogram( hS_1, mchName2, 1 ) ;
   h_draw->Draw(         hS_,        "", " EcalTime ", "", "logY", 0.95, 2 ) ;
   h_draw->DrawAppend(  hS_1, "SmearMC",  0.73, 4 ) ;
   */
   string mchName2 = mchName ;
   mchName2 += "_s" ;
   TH1D* hS_2 = statslib->RebinHistogram( hS_, mchName2, 1 ) ;

 
   // normalize background for the center peak 
   double nLT3J = hB_->Integral(10, 26) ;
   double nGE3J = hData_->Integral(10, 26) ;
   double bgNorm = nGE3J/ nLT3J ;
   cout<<" background normalization "<< bgNorm <<" from "<< hB_->Integral() <<" => " ;
   hB_->Scale( bgNorm ) ;
   cout<< hB_->Integral() << endl ;

   // turn the observation scope at certain timing range
   TH1D* hS = statslib->SubHistogram( hS_2 ,      timeScope[0] , timeScope[1] ) ;
   TH1D* hB = statslib->SubHistogram( hB_ ,       timeScope[0] , timeScope[1] ) ;
   TH1D* hData = statslib->SubHistogram( hData_ , timeScope[0] , timeScope[1] ) ;

   // normalization : 0.0223 = 0.1 * 11180(/pb) / 50112(generated events)
   // assuming 0.1 pb for this signal model 
   double norm  = lumi*sgXsec[norm_id] / nGen[norm_id] ;
   double expS0 = hS->Integral() ;
   double expS  = norm*expS0 ;
   double expB  = hB->Integral() ;
   int    nExp  = (int)(expB) ;
   int    nObs  = (int) hData->Integral() ;
   double Eff   = expS0 / nGen[norm_id] ;
   cout<<" normal: "<< norm <<" Expected Signal: "<< expS <<" Background: "<< expB ;
   cout<<" nOBs: "<< nObs <<" Eff: "<< Eff <<endl ;

   // normalized signal and background model
   hS->Scale( 1./ hS->Integral() ) ;
   hB->Scale( 1./ hB->Integral() ) ;

   //double ks = statslib->KSTest( hS, hB, 12, 80 , -1.5, 0.25 ) ;
   h_draw->Draw(       hS, "", " EcalTime ", "", "logY", 0.95, 2 ) ;
   h_draw->DrawAppend( hB, "PDF_SB",  0.75, 1 ) ;
   h_draw->Draw(   hData, "hData", " EcalTime ", "", "logY", 0.95, 4 ) ;

   // calculate results
   double alpha = 1 - searchLimit  ;
   int nToys ;
   Input->GetParameters("NToys",    &nToys );
   double limit0, limit1u, limit1d, limit2u, limit2d, limit ;


   if ( pValue_B > 0. ) {
         cout<<" Calculating for background-only Expectation "<<endl;
	 double limit_ ; 
         if ( nObs < 100 ) limit_ = statslib->ExpectedLimit( nExp, expS, expB, hS, hB, alpha, pValue_B, nToys, subPath_t ) ;
         else              limit_ = statslib->ExpectedLimit1( nExp, expS, expB, hS, hB, alpha, pValue_B, nToys, subPath_t ) ;
	 fprintf(logfile," %.4f  %.3f  %.3f  %.5f \n", pValue_B, limit_, expS, Eff );

   } else if ( pValue_B < 0.00001 && pValue_B > -1 ) {
         cout<<" Calculating for Data "<<endl ;
	 double limit_ ; 
         if ( nObs < 100 ) limit_ = statslib->ExpectedLimit( nObs, expS, expB, hS, hB, alpha, 0.4999, nToys, subPath_0 ) ;
         else              limit_ = statslib->ExpectedLimit1( nObs, expS, expB, hS, hB, alpha, 0.4999, nToys, subPath_0 ) ;
	 fprintf(logfile," %.3f  %.3f  %.5f \n", limit_, expS, Eff );

   } else {
         cout<<" Scanning All - data and background-only model ! "<<endl ;

      if ( nObs < 100 ) {
          
	 limit0  = statslib->ExpectedLimit( nExp, expS, expB, hS, hB, alpha, 0.5 ,    nToys, subPath_c ) ;
	 limit1u = statslib->ExpectedLimit( nExp, expS, expB, hS, hB, alpha, 0.1587 , nToys, subPath_u1 ) ;
	 limit1d = statslib->ExpectedLimit( nExp, expS, expB, hS, hB, alpha, 0.8413 , nToys, subPath_d1 ) ;
	 limit2u = statslib->ExpectedLimit( nExp, expS, expB, hS, hB, alpha, 0.0228 , nToys, subPath_u2 ) ;
	 limit2d = statslib->ExpectedLimit( nExp, expS, expB, hS, hB, alpha, 0.9772 , nToys, subPath_d2 ) ;
	 limit   = statslib->ExpectedLimit( nObs, expS, expB, hS, hB, alpha, 0.4999,  nToys, subPath_0  ) ;
         fprintf(logfile," %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.5f \n", 
                 limit0 - limit2d, limit0 - limit1d, limit0 , limit1u - limit0, limit2u - limit0, limit, expS, Eff );
 
      } else {
         time_t t1_ini, t1_end, t2_ini, t2_end ;

         time( &t1_ini) ;
         // data
	 limit   = statslib->ExpectedLimit1( nObs, expS, expB, hS, hB, alpha, 0.4999,  nToys, subPath_0  ) ;
         time( &t1_end) ;

         time( &t2_ini) ;
         /*
         limit0  = statslib->ExpectedLimit1( nExp, expS, expB, hS, hB, alpha, 0.5 ,    nToys, subPath_c ) ;
	 limit1u = statslib->ExpectedLimit1( nExp, expS, expB, hS, hB, alpha, 0.1587 , nToys, subPath_u1 ) ;
	 limit1d = statslib->ExpectedLimit1( nExp, expS, expB, hS, hB, alpha, 0.8413 , nToys, subPath_d1 ) ;
	 limit2u = statslib->ExpectedLimit1( nExp, expS, expB, hS, hB, alpha, 0.0228 , nToys, subPath_u2 ) ;
	 limit2d = statslib->ExpectedLimit1( nExp, expS, expB, hS, hB, alpha, 0.9772 , nToys, subPath_d2 ) ;
	 fprintf(logfile," %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.5f \n", 
                 limit0 - limit2d, limit0 - limit1d, limit0 , limit1u - limit0, limit2u - limit0, limit, expS, Eff );
         */
         // mc background-only 
         vector<double> limitV = statslib->ExpectedLimit2( nExp, expS, expB, hS, hB, alpha,  nToys, subPath_c  ) ;
         time( &t2_end) ;
	 fprintf(logfile," %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.5f \n", 
                limitV[2] - limitV[0], limitV[2] - limitV[1], limitV[2], limitV[3] - limitV[2], limitV[4] - limitV[2], limit, expS, Eff );

         double dT1 = difftime( t1_end, t1_ini) ;
         double dT2 = difftime( t2_end, t2_ini) ;
         printf(" **** dT1 = %.3f", dT1 ) ;
         printf(" **** dT2 = %.3f", dT2 ) ;
      }
   }
  

   //TH1D* h95 = new TH1D("h95", "EcalTime for 95 CL Uplimit ", 80, -4.5, 15.5 );
   //statslib->toyMC( h95, limit, hS, true ) ;
   //h_draw->Draw( h95, "h95Data", " EcalTime ", "", "logY", 0.95, 2 ) ;

   h_draw->SetNewPath( "" ) ;
   theFile->Close() ;
   fclose( logfile ) ;
} 

// Example
void StatsTool::MyCLsTest() {

   // generate samples for signal and background pdf of observable x 
   TH1D* hS = new TH1D("hS", "Signal", 50, 0, 10. );
   statslib->Sample( hS, "Gaus", 10000, 7, 1. ) ;
   TH1D* hB = new TH1D("hB", "Background", 50, 0, 10. );
   //statslib->Sample( hB, "Gaus", 10000, 5, 2. ) ;
   statslib->Sample( hB, "Exp", 10000, 3 ) ;

   hS->Scale( 1./ hS->Integral() ) ;
   hB->Scale( 1./ hB->Integral() ) ;

   h_draw->Draw(       hS, "", " Signal & Background ", "", "", 0.95, 2 ) ;
   h_draw->DrawAppend( hB, "PDF_SB",  0.75, 1 ) ;

   // generated fake data
   int nObs = 55  ;
   vector<double> xV = statslib->toyMC( nObs, hB, 0 );
   TH1D* hData = new TH1D("hData", " Data ", 50, 0, 10. );
   for ( size_t i=0; i< xV.size() ; i++) {
       hData->Fill( xV[i] ) ;
   }
   h_draw->Draw(   hData, "hData", " Data ", "", "", 0.95, 4 ) ;

   // calculate results
   // Find the expected limit
   //TH1D* hPdfQ_sb = new TH1D("hPdfQ_sb", " PDF(Q|H) ", 200,  -10, 10 );
   //TH1D* hPdfQ_b  = new TH1D("hPdfQ_b", " PDF(Q|Hb) ", 200,  -10, 10 );

   // Loop through different signal normalizations
   double expS  = 3.5 ;
   double expB  = 50 ;
   double alpha = 0.05 ;
   double limit0 = statslib->ExpectedLimit( nObs, expS, expB, hS, hB, alpha, 0.5 ) ;
   double limit1a = statslib->ExpectedLimit( nObs, expS, expB, hS, hB, alpha, 0.1587 ) ;
   double limit1b = statslib->ExpectedLimit( nObs, expS, expB, hS, hB, alpha, 0.8413 ) ;
   double limit2a = statslib->ExpectedLimit( nObs, expS, expB, hS, hB, alpha, 0.0228 ) ;
   double limit2b = statslib->ExpectedLimit( nObs, expS, expB, hS, hB, alpha, 0.9772 ) ;
   cout<<" ====== Report Limit2a : "<< limit2a << endl ;
   cout<<" ====== Report Limit1a : "<< limit1a << endl ;
   cout<<" ====== Report Limit0  : "<< limit0 << endl ;
   cout<<" ====== Report Limit1b : "<< limit1b << endl ;
   cout<<" ====== Report Limit2b : "<< limit2b << endl ;

} 

void StatsTool::Toys() {

   TH1D* hS = new TH1D("hS", "Signal", 50, 0, 10. );
   statslib->Sample( hS, "Gaus", 10000, 7, 1. ) ;

   hS->Scale( 1./ hS->Integral() ) ;

   h_draw->Draw( hS, "PDF", " PDF", "", "", 0.95, 2 ) ;
  
   TH1D* hToy = (TH1D*) hS->Clone() ;
   for (int i=0; i< 10 ; i++) {
       char hname[15] ;
       sprintf(hname, "Toy_S%02d", i ) ;
       hToy->Reset() ;
       statslib->toyMC( hToy, 10000, hS, i ) ;
       h_draw->Draw( hToy, hname, " toy ", "", "", 0.95, 4 ) ;
   }

}

void StatsTool::DrawLimit( ) {

   // open the text file for reading the result
   string currentPath = h_draw->CurrentPath() ;
   string logfilePath = currentPath ;
   logfilePath += "limitReport.log" ;
   
   FILE* fitlog = fopen( logfilePath.c_str() ,"r");

   const int sz = (int) nGen.size() ;

   cout<<" Testing "<< sz <<" models "<<endl ;
   vector<double> testModel ;
   Input->GetParameters("TestModel",   &testModel ) ; 
 
   int r1 ;
   float c0, u1, u2, d1, d2, ob, ep, ef;
   double ac0[sz],  au1[sz], au2[sz], ad1[sz], ad2[sz], aob[sz], aep[sz] ;
   double ax[sz], aex[sz] ;

   for ( int i=0; i< sz ; i++ ) {
       ax[i] = testModel[i] ;
       aex[i] = 0. ;
       r1 = fscanf(fitlog, "%f", &d2 );
       r1 = fscanf(fitlog, "%f", &d1 );
       r1 = fscanf(fitlog, "%f", &c0 );
       r1 = fscanf(fitlog, "%f", &u1 );
       r1 = fscanf(fitlog, "%f", &u2 );
       r1 = fscanf(fitlog, "%f", &ob );
       r1 = fscanf(fitlog, "%f", &ep );
       r1 = fscanf(fitlog, "%f", &ef );
       au2[i] = u2 / (ef*lumi) ;
       au1[i] = u1 / (ef*lumi) ;
       ac0[i] = c0 / (ef*lumi) ;
       ad1[i] = d1 / (ef*lumi) ;
       ad2[i] = d2 / (ef*lumi) ;
       aob[i] = ob / (ef*lumi) ;
       //aep[i] = ep / (ef*lumi) ;
       aep[i] = sgXsec[i] ;
       //printf(" c0: %f, ef: %f , ac0: %f \n", c0, ef, ac0[i] ) ;
       //printf(" %f, %f, %f, %f, %f \n", d2, d1, c0, u1, u2 ) ;
       if ( r1 != 1 ) cout<<" reading error "<<endl ;
   }
   

   //TCanvas *c1a = new TCanvas("c1a","c1a",200,10,600,400);
   TCanvas *c1a = new TCanvas("c1a","c1a" ,800,600);

   c1a->SetLogy() ;
   TMultiGraph *mg = new TMultiGraph();
   mg->SetTitle( " 95% CLs Limit " );

   TGraphAsymmErrors* gbe = new TGraphAsymmErrors(sz, ax, ac0, aex, aex, ad2, au2);
   gbe->SetFillStyle(3001);
   gbe->SetFillColor(kYellow);
   //gbe->GetXaxis()->SetLimits(1.,5.) ;
   
   TGraphAsymmErrors* gae = new TGraphAsymmErrors(sz, ax, ac0, aex, aex, ad1, au1);
   gae->SetFillStyle(3001);
   gae->SetFillColor(kSpring);
   //gae->GetXaxis()->SetRangeUser(1.,5.) ;
   //gae->GetXaxis()->SetLimits(1.,5.) ;
   
   mg->Add( gbe ) ;
   mg->Add( gae ) ;
   mg->Draw("a3") ;
 
   // this modify the x-axis range, must set after Draw() ...stupid ROOT 
   gPad->Modified();
   mg->GetXaxis()->SetLimits( ax[0], ax[sz-1] ) ;
   mg->SetMaximum(1.0);
   mg->SetMinimum(0.001);
   mg->GetXaxis()->SetTitleOffset(1.);
   mg->GetYaxis()->SetTitleOffset(1.);
   mg->GetXaxis()->SetTitleFont(42);
   mg->GetYaxis()->SetTitleFont(42);
   mg->GetXaxis()->SetTitleSize(0.04);
   mg->GetYaxis()->SetTitleSize(0.04);
   mg->GetXaxis()->SetTitle( "#chi^{0} Lifetime" ) ;
   mg->GetYaxis()->SetTitle(" #sigma_{Upper Limit} (pb) at 95% CL_{s}") ;

   // expected
   TGraph*  g1  = new TGraph(sz, ax, ac0 ) ;
   g1->SetLineColor(1) ;
   g1->SetLineWidth(3) ;
   g1->SetLineStyle(7) ;
   //g1->GetXaxis()->SetLimits(1.,5.) ;
   g1->Draw("LP") ;

   // observation
   TGraph*  g0  = new TGraph(sz, ax, aob ) ;
   g0->SetLineColor(2) ;
   g0->SetLineWidth(2) ;
   g0->SetLineStyle(1) ;
   //g0->GetXaxis()->SetLimits(1.,5.) ;
   //g0->SetMarkerStyle(20) ;
   g0->Draw("LP") ;
   c1a->Update() ;

   // theoratical suggestion
   TGraph*  gth = new TGraph(sz, ax, aep ) ;
   gth->SetLineColor(4) ;
   gth->SetLineWidth(3) ;
   gth->SetLineStyle(1) ;
   gth->GetXaxis()->SetLimits( ax[0], ax[sz-1] ) ;
   gth->Draw("L") ;

   TLegend* leg1  = new TLegend(.65, .73, .90, .90 );
   leg1->SetFillColor(10) ;
   leg1->SetTextSize(0.030) ;
   leg1->AddEntry( gbe, "#pm 2#sigma" ,  "F");
   leg1->AddEntry( gae, "#pm 1#sigma" ,  "F");
   leg1->AddEntry( g1,  "Expected" , "L");
   leg1->AddEntry( gth, "Theory" , "L");
   leg1->AddEntry( g0,  "Observed" , "L");
   leg1->Draw("same") ;

   TString gPlotname = hfolder +  "Uplimit."+plotType ;
   c1a->Print( gPlotname ) ;
    
}

