#include "PhotonAna.h"

PhotonAna::PhotonAna(){

  Input = new AnaInput();
  
  Input->GetParameters("PlotType", &plotType ) ; 
  Input->GetParameters("Path",     &hfolder ) ; 
  Input->GetParameters("ProcessEvents",     &ProcessEvents ) ; 

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
   TTree* tr = Input->TreeMap( "phyRh-Photon-Run2011A-PromptReco-v6-AOD-Wjson-CORR-nov2.HADDED" );

   float jetPx[10], jetPy[10], jetPz[10], jetE[10] ;
   float phoPx[10], phoPy[10], phoPz[10], phoE[10] ;
   float CPIdx[MAXC], clusterTime[MAXC], clusterEnergy[MAXC]; 
   float xtalInBCEnergy[MAXC][MAXXTALINC],  xtalInBCTime[MAXC][MAXXTALINC] ;
   int clusterXtals[MAXC] ;
   int nJets, nPhotons, eventId ;

   tr->SetBranchAddress("eventId",    &eventId);
   tr->SetBranchAddress("nJets",      &nJets);
   tr->SetBranchAddress("nPhotons",   &nPhotons);
   tr->SetBranchAddress("jetPx",       jetPx );
   tr->SetBranchAddress("jetPy",       jetPy );
   tr->SetBranchAddress("jetPz",       jetPz );
   tr->SetBranchAddress("jetE",        jetE );
   tr->SetBranchAddress("phoPx",       phoPx );
   tr->SetBranchAddress("phoPy",       phoPy );
   tr->SetBranchAddress("phoPa",       phoPz );
   tr->SetBranchAddress("phoE",        phoE );

   tr->SetBranchAddress("CPIdx",        CPIdx );
   tr->SetBranchAddress("clusterTime",  clusterTime );
   tr->SetBranchAddress("clusterEnergy",  clusterEnergy );
   tr->SetBranchAddress("clusterXtals",   clusterXtals );
  
   tr->SetBranchAddress("xtalInBCEnergy", xtalInBCEnergy );
   tr->SetBranchAddress("xtalInBCTime",   xtalInBCTime );

   hJetTime hJets ;

   int totalN = tr->GetEntries();
   cout<<" total entries = "<< totalN <<" Process "<< ProcessEvents <<endl;
    
   for ( int i=0; i< tr->GetEntries(); i++ ) {
       if ( ProcessEvents > 0 && i > ( ProcessEvents - 1 ) ) break;
       tr->GetEntry( i );
     
       hJets.nJets->Fill( nJets );
       hJets.nPhotons->Fill( nPhotons );
       //cout<<" Event "<< eventId <<endl ;

       if ( nPhotons < 1 || nJets < 3 ) continue ;
      
           // get the associated bc for photon
           TLorentzVector gp4( phoPx[0], phoPy[0], phoPz[0], phoE[0] ) ;
           hJets.pho0Pt->Fill( gp4.Pt() );
           //cout<<" photon 0  pt : "<< gp4.Pt() <<endl;
           float gTime = -999. ;
           float cE = 0. ;
           for (int x =0; x< 200; x++) {
               if ( CPIdx[x] < 22. || CPIdx[x] > 23. ) continue ;
               float gidx = 22.1 ;
               if ( CPIdx[x] == gidx && fabs(clusterTime[x]) < 10.  ) {
                  //cout<<"  cluster "<< CPIdx[x] <<" time : "<< clusterTime[x] ;
                  //cout<<" E: "<< clusterEnergy[x]<<" N: "<< clusterXtals[x]  << endl;
                  if ( clusterEnergy[x] > cE ) {
                     cE    = clusterEnergy[x] ; 
                     gTime = clusterTime[x] ;
                  }
               }
           }
           if ( gTime != -999 ) hJets.pho0T->Fill( gTime );

           // get the associated bc for jet
           for ( int j=0 ; j< nJets; j++) {
               if ( j > 2 ) break ; 
               TLorentzVector jp4( jetPx[j], jetPy[j], jetPz[j], jetE[j] ) ;
	       if ( j == 0 ) hJets.jet0Pt->Fill( jp4.Pt() );
	       if ( j == 1 ) hJets.jet1Pt->Fill( jp4.Pt() );

	       //cout<<" jet"<<j<<"  pt : "<< jp4.Pt() <<endl;
	       float jTime = 0 ;
	       float nC = 0 ;
	       float jxTime = 0 ;
               float nX = 0 ;
	       for (int x =0; x< 200; x++) {
                  if ( CPIdx[x] < 100. || CPIdx[x] > 101. ) continue ;
		  float jidx =  100. + (j+1.)*0.1 ;
		  if ( CPIdx[x] == jidx  && fabs(clusterTime[x]) < 10. ) {
		      //cout<<"  cluster "<<CPIdx[x] <<" time : "<< clusterTime[x] ;
		      //cout<<" E: "<< clusterEnergy[x]<<" N: "<< clusterXtals[x]  << endl;
		      jTime += clusterTime[x] ;
		      nC += 1. ;
                      for (int y=0; y< 25; y++) {
                          if ( xtalInBCEnergy[x][y] == 0 ) continue ;
                          //cout<<"      Xtal"<< y <<"  T: "<< xtalInBCTime[x][y]<<" E:"<<xtalInBCEnergy[x][y]<<endl;
                          jxTime += xtalInBCTime[x][y] ;
                          nX += 1. ;
                      }
                  }
               }
               if ( nC != 0 || nX != 0) {
                  jTime = jTime/nC ;
                  jxTime = jxTime/nX ;
		  if ( j == 0 ) { 
                     hJets.jet0T->Fill( jTime ); 
                     hJets.jet0xT->Fill( jxTime ); 
                     hJets.jet0Pt_T->Fill( jp4.Pt() , jTime ) ;
                  }
		  if ( j == 1 ) {
                     hJets.jet1T->Fill( jTime );
                     hJets.jet1xT->Fill( jxTime );
                  }
	       }
           }
   }

   hJets.Draw( hfolder, plotType ) ;
   cout<<" fk3 "<<endl ; 

}
