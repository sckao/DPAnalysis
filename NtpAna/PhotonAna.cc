#include "PhotonAna.h"

PhotonAna::PhotonAna(){

  Input = new AnaInput();
  
  Input->GetParameters("PlotType", &plotType ) ; 
  Input->GetParameters("Path",     &hfolder ) ; 
  Input->GetParameters("ProcessEvents",     &ProcessEvents ) ; 

  Input->GetParameters( "PhotonCuts", &photonCuts );

  //Input->GetParameters("Debug", &debugStr ) ; 
  //if ( debugStr == "True" ) debug = true ;
}

PhotonAna::~PhotonAna(){

  cout<<" finish !" <<endl;
  delete Input;
  cout<<" done ! "<<endl ;

}

// analysis template
void PhotonAna::test() { 

  double lumi = 0 ;
  Input->GetParameters( "Lumi", &lumi ) ;
  cout <<" lumi = "<< lumi <<endl;

}
 
void PhotonAna::ReadTree() { 

   //vector<string> rFiles ;
   //Input->GetParameters( "TheData", &rFiles );
   //TTree* tr = Input->GetTree( rFiles[0],"EcalTimeAnalysis" );
   TTree* tr = Input->TreeMap( "data+" );

   float jetPx[10], jetPy[10], jetPz[10], jetE[10] ;
   float phoPx[10], phoPy[10], phoPz[10], phoE[10] ;
   float elePx[10], elePy[10], elePz[10], eleE[10] ;
   float CPIdx[MAXC], clusterTime[MAXC], clusterEnergy[MAXC]; 
   float xtalInBCEnergy[MAXC][MAXXTALINC],  xtalInBCTime[MAXC][MAXXTALINC] , xtalInBCTimeErr[MAXC][MAXXTALINC];
   float xtalInBCEta[MAXC][MAXXTALINC], xtalADC[MAXC][MAXXTALINC] ;
   int clusterXtals[MAXC] ;
   int nJets, nPhotons, nElectrons, eventId ;

   tr->SetBranchAddress("eventId",    &eventId);
   tr->SetBranchAddress("nJets",      &nJets);
   tr->SetBranchAddress("nPhotons",   &nPhotons);
   tr->SetBranchAddress("nElectrons", &nElectrons);
   tr->SetBranchAddress("jetPx",       jetPx );
   tr->SetBranchAddress("jetPy",       jetPy );
   tr->SetBranchAddress("jetPz",       jetPz );
   tr->SetBranchAddress("jetE",        jetE );
   tr->SetBranchAddress("phoPx",       phoPx );
   tr->SetBranchAddress("phoPy",       phoPy );
   tr->SetBranchAddress("phoPz",       phoPz );
   tr->SetBranchAddress("phoE",        phoE );
   tr->SetBranchAddress("elePx",       elePx );
   tr->SetBranchAddress("elePy",       elePy );
   tr->SetBranchAddress("elePz",       elePz );
   tr->SetBranchAddress("eleE",        eleE );

   tr->SetBranchAddress("CPIdx",        CPIdx );
   tr->SetBranchAddress("clusterTime",  clusterTime );
   tr->SetBranchAddress("clusterEnergy",  clusterEnergy );
   tr->SetBranchAddress("clusterXtals",   clusterXtals );
  
   tr->SetBranchAddress("xtalInBCEnergy",  xtalInBCEnergy );
   tr->SetBranchAddress("xtalInBCTime",    xtalInBCTime );
   tr->SetBranchAddress("xtalInBCTimeErr", xtalInBCTimeErr );
   tr->SetBranchAddress("xtalInBCEta",     xtalInBCEta );
   tr->SetBranchAddress("xtalInBCAmplitudeADC", xtalADC );

   int totalN = tr->GetEntries();
   cout<<" total entries = "<< totalN <<" Process "<< ProcessEvents <<endl;
    
   for ( int i=0; i< tr->GetEntries(); i++ ) {
       if ( ProcessEvents > 0 && i > ( ProcessEvents - 1 ) ) break;
       tr->GetEntry( i );
     
       //cout<<" Event "<< eventId <<endl ;
       hJets.nJets->Fill( nJets );
       hJets.nPhotons->Fill( nPhotons );
       hJets.nElectrons->Fill( nElectrons );

       if ( nPhotons < 1 || nJets < 3 ) continue ;
      
           //// 1. get the associated bc for photon
           TLorentzVector gp4( phoPx[0], phoPy[0], phoPz[0], phoE[0] ) ;
           if ( gp4.Pt() < photonCuts[0] ) continue ;
           hJets.pho0Pt->Fill( gp4.Pt() );

           //// 1.1 get the leading electron info, just check
           if ( nElectrons >  0 ) {
              TLorentzVector ep4( elePx[0], elePy[0], elePz[0], eleE[0] ) ;
              hJets.ele0Pt->Fill( ep4.Pt() );
           }

           //cout<<" photon 0  pt : "<< gp4.Pt() <<endl;
           float gTime  = 0. ;
           float gxTime = 0. ;
	   float gleadEnergy = 0 ;
	   float gleadTime   = 0 ;
	   float ngX = 0. ;
	   float ngC = 0. ;
           for (int x =0; x< 200; x++) {
               if ( CPIdx[x] < 22. || CPIdx[x] > 23. ) continue ;
               float gidx = 22.1 ;
               if ( CPIdx[x] != gidx || fabs(clusterTime[x]) > 10.  ) continue ;
               //cout<<"  cluster "<< CPIdx[x] <<" time : "<< clusterTime[x] ;
               //cout<<" E: "<< clusterEnergy[x]<<" N: "<< clusterXtals[x]  << endl;

	       float seedEnergy = 0 ;
	       float seedTime   = 0 ;
	       float seedTErr   = 999999 ;
               for (int y=0; y< 25; y++) {
                   if ( xtalInBCEnergy[x][y] < 0.01 ) continue ;
		   if ( xtalInBCTimeErr[x][y] < 0.2 || xtalInBCTimeErr[x][y] > 10. ) continue ;
		   hJets.phoXtalTime->Fill( xtalInBCTime[x][y] ); 
		   hJets.phoXtalTErr->Fill( xtalInBCTimeErr[x][y] ); 
		   if ( xtalInBCEnergy[x][y] > gleadEnergy ) {                         
                      gleadEnergy = xtalInBCEnergy[x][y] ;
                      gleadTime   = xtalInBCTime[x][y] ;
                   }
                   if ( xtalInBCEnergy[x][y] > seedEnergy ) {
                      seedEnergy = clusterEnergy[x] ; 
                      seedTime   = xtalInBCTime[x][y] ;
                      seedTErr   = xtalInBCTimeErr[x][y] ;
                   }
		   gxTime += xtalInBCTime[x][y] / pow( xtalInBCTimeErr[x][y], 2) ;
		   ngX += 1./ pow( xtalInBCTimeErr[x][y],2) ;
               }
               gTime +=  seedTime / pow( seedTErr, 2) ;
               ngC   +=       1. / pow( seedTErr, 2) ;
           }
           if ( ngX > 0. && ngC > 0. ) {
              gTime  = gTime  / ngC ;
              gxTime = gxTime / ngX ;
              hJets.pho0T->Fill( gTime );
              hJets.pho0xT->Fill( gxTime ); 
              hJets.pho0lT->Fill( gleadTime ); 
           }

           //// 2. get the associated bc for jet
           for ( int j=0 ; j< nJets; j++) {
               if ( j > 1 ) break ; 
               TLorentzVector jp4( jetPx[j], jetPy[j], jetPz[j], jetE[j] ) ;
	       if ( j == 0 ) hJets.jet0Pt->Fill( jp4.Pt() );
	       if ( j == 1 ) hJets.jet1Pt->Fill( jp4.Pt() );

	       //cout<<" jet"<<j<<"  pt : "<< jp4.Pt() <<endl;
               float emF = 0 ;
	       float nC = 0 ;
               float nX = 0 ;
               float xE = 0 ;
               float cE = 0 ;
	       float jleadEnergy = 0 ;
	       float jleadTime   = 0 ;
	       float jTime  = 0 ;
	       float jxTime = 0 ;
               float xADC = 0 ;  
               int   nXtal = 0 ; 
               // loop the clusters
	       for (int x =0; x< 200; x++) {
                  if ( CPIdx[x] < 100. || CPIdx[x] > 101. ) continue ;

		  float jidx =  100. + (j+1.)*0.1 ;
		  if ( CPIdx[x] != jidx  || fabs(clusterTime[x]) > 10. ) continue ;

                  float seedEnergy = 0 ;
	          float seedTime   = 0 ;
	          float seedTErr   = 999999 ;
		  //cout<<"  cluster "<<CPIdx[x] <<" time : "<< clusterTime[x] ;
		  //cout<<" E: "<< clusterEnergy[x]<<" N: "<< clusterXtals[x]  << endl;
		  float xADC_BC = 0. ;
		  // loop the xtals
                  for (int y=0; y< 25; y++) {
                      if ( xtalInBCEnergy[x][y] < 0.01 ) continue ;
		      if ( xtalInBCTimeErr[x][y] < 0.2 || xtalInBCTimeErr[x][y] > 10. ) continue ;
		      hJets.jetXtalTime->Fill( xtalInBCTime[x][y] ); 
		      hJets.jetXtalTErr->Fill( xtalInBCTimeErr[x][y] ); 

		      if ( xtalInBCEnergy[x][y] > jleadEnergy ) {                         
                         jleadEnergy = xtalInBCEnergy[x][y] ;
                         jleadTime   = xtalInBCTime[x][y] ;
                      }
		      if ( xtalInBCEnergy[x][y] > seedEnergy ) {                         
                         seedEnergy = xtalInBCEnergy[x][y] ;
                         seedTime   = xtalInBCTime[x][y] ;
                         seedTErr   = xtalInBCTimeErr[x][y] ;
                      }
                      jxTime += xtalInBCTime[x][y] / pow( xtalInBCTimeErr[x][y], 2) ;
		      nX += 1./ pow( xtalInBCTimeErr[x][y],2) ;
		      //double sinTheta   = fabs( sin( 2. *atan( exp(-1*xtalInBCEta[x][y]  ) ) ) );
		      emF += xtalInBCEnergy[x][y] ;
		      xE  += xtalInBCEnergy[x][y] ;
		      xADC_BC = ( xtalADC[x][y] > xADC_BC ) ? xtalADC[x][y] : xADC_BC ;
                      nXtal++ ;
                  }
		  jTime += seedTime/ pow( seedTErr, 2 ) ;
		  nC    += 1./ pow( seedTErr, 2) ;
                  xADC = ( xADC_BC > xADC ) ? xADC_BC : xADC ;
                  if ( nXtal > 0 ) cE += clusterEnergy[x] ;
                      
               }
               if ( nC > 0. && nXtal > 0 ) {
                  jTime  = jTime / nC ;
                  jxTime = jxTime / nX ;
		  emF    = emF / jp4.E() ;
                  hJets.jet_emF_T->Fill( emF , jTime ) ;
                  hJets.xE_cE->Fill( xE , cE ) ;
                  hJets.jet_Pt_nXtal->Fill( jp4.Pt() , nXtal ) ;
		  if ( j == 0 ) { 
                     hJets.jet0T->Fill( jTime ); 
                     hJets.jet0xT->Fill( jxTime ); 
                     hJets.jet0lT->Fill( jleadTime ); 
                     hJets.jet0Pt_T->Fill( jp4.Pt() , jTime ) ;
                     hJets.jet0Eta_T->Fill( jp4.Eta() , jTime ) ;
                     hJets.jet0Pt_ADC->Fill( jp4.Pt() , xADC ) ;
                     hJets.jet0Eta_ADC->Fill( jp4.Eta() , xADC ) ;
                     //cout<<"  jxTime: "<< jxTime <<" nX: "<< nX <<" nC: "<<nC<<" jTime:"<< jTime << endl ;
                  }
		  if ( j == 1 ) {
                     hJets.jet1T->Fill( jTime );
                     hJets.jet1xT->Fill( jxTime );
                     hJets.jet1lT->Fill( jleadTime ); 
                  }
	       }
           }
   }

   hJets.Draw( hfolder, plotType ) ;

}

void PhotonAna::ScalarPlotList() {

   bool debugPlots = false ;
   ScalarPlotter( hJets.jet0Eta_ADC,  "jet_Eta_ADC",   0., 1000., 2 ) ;
   ScalarPlotter( hJets.jet0Pt_ADC,   "jet_Pt_ADC",    0., 1000., 2, debugPlots ) ;
   ScalarPlotter( hJets.jet_Pt_nXtal, "jet_Pt_nXtal",  0., 250., 2 ) ;
   ScalarPlotter( hJets.jet0Eta_T,    "jet_Eta_Time",  -2,   2 ) ;
   ScalarPlotter( hJets.jet0Pt_T,     "jet_Pt_Time",   -2,   2 ) ;

}

void PhotonAna::ScalarPlotter( TH2D* h2, TString hname, double yMin, double yMax, int rbin, bool debugPlots ) {

    TCanvas* c1 = new TCanvas("c1","", 800, 600);
    c1->SetFillColor(10);
    c1->SetFillColor(10);
    gPad->SetGridx();
    gPad->SetGridy();

    vector<double> xV ;
    vector<double> yV ;
    vector<double> yErrV ;
    BinningFitScan( h2, xV, yV, yErrV, debugPlots, rbin ) ; 

    const int sz = static_cast<const int>( xV.size() );
    float  x[sz] , xErr[sz], y[sz], yErr[sz] ;
    for ( int i=0; i< sz; i++) {
        x[i] = xV[i] ;
        y[i] = yV[i] ;
        yErr[i] = yErrV[i] ;
        xErr[i] = 0.  ;
    }

    c1->cd();
    gStyle->SetOptStat(kTRUE);
    gStyle->SetOptFit(111);

    TGraphErrors* g1 = new TGraphErrors( sz, x, y,  xErr, yErr );
    g1->SetMaximum( yMax );
    g1->SetMinimum( yMin );
    g1->SetMarkerColor(4);
    g1->SetMarkerStyle(20);
    g1->SetLineWidth(2);
    g1->SetLineColor(4);
    g1->SetTitle( hname );
    g1->Draw("AP");
 
    c1->Update();
    TString plotname = hfolder + hname +  "_pj."+plotType ;
    c1->Print( plotname );


    delete g1 ;
    delete c1 ;
}


void PhotonAna::BinningFitScan( TH2D* h2, vector<double>& xV, vector<double>& yV, vector<double>& yErrV, bool debugPlots, int rbin, int startBin, int finalBin ) {

    int      sz = h2->GetNbinsX();
    for ( int i=0; i< sz; i++) {
        int b1 = startBin + (i*rbin) ;
        int b2 = b1 + rbin - 1 ;
        if ( b1 == finalBin ) break;
        vector<double> results = BinningFit( h2, "h2Py", b1, b2, debugPlots ) ;

        if ( results[2] == 0. ) continue ;
        xV.push_back( h2->GetBinCenter(b1) ) ;
        yV.push_back( results[1] ) ;
        yErrV.push_back( results[2] ) ;
    }

}


vector<double> PhotonAna::BinningFit( TH2D* h2, string hName, int xbinMin, int xbinMax, bool debugPlots ) { 

   TH1D* h1_pjy = h2->ProjectionY( hName.c_str() , xbinMin, xbinMax ) ;

   double p1_ = h1_pjy->GetMean(1) ;
   double p2_ = h1_pjy->GetRMS(1) ;
   double yMin = p1_ - (3*p2_) ;
   double yMax = p1_ + (3*p2_) ;
   //cout<<" ***  m = "<<p1_<<" s = "<< p2_ <<" h2 m:"<< h2->GetMean(2) << endl ;
   TF1 *func = new TF1("func",MathFunctions::fitGS, yMin, yMax, 3 );
   func->SetParLimits(1, p1_-p2_ , p1_+p2_ );
   func->SetParLimits(2, p2_*0.8 , p2_*1.2 );

   h1_pjy->Fit(func, "R", "", yMin, yMax);
 
   double amp   = func->GetParameter(0) ;
   double mean  = func->GetParameter(1) ;
   double sigma = func->GetParameter(2) ;

   vector<double> result ;
   result.push_back( amp );
   result.push_back( mean );
   result.push_back( sigma );

   TString plotname0 ;
   if ( debugPlots ) {
      cout <<" debugging slice&fit"<< endl; 
      gSystem->cd( hfolder.c_str() ) ;
      gSystem->mkdir( "debug" );
      gSystem->cd( "../" );
     
      TCanvas* c0 = new TCanvas("c0","", 800, 600);
      c0->SetFillColor(10);
      c0->SetFillColor(10);
      c0->cd();
      h1_pjy->Draw() ;
      c0->Update();
      ostringstream sfxStr ;
      sfxStr << "_py" ;
      sfxStr << xbinMin ;

      plotname0= hfolder + "debug/" + hName +  sfxStr.str() +"."+plotType ;
      c0->Print( plotname0 );
      delete c0 ;
   }

   delete h1_pjy ;
   delete func ;
   return result ;

}

