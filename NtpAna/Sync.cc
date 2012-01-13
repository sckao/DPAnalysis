#include "Sync.h"

Sync::Sync( string datacardfile ){

  Input = new AnaInput( datacardfile );

  select = new Selection( datacardfile ) ;
  
  Input->GetParameters("PlotType",      &plotType ) ; 
  Input->GetParameters("Path",          &hfolder ) ; 
  Input->GetParameters("ProcessEvents", &ProcessEvents ) ; 

  Input->GetParameters( "SelectBackground",  &selectBackground );
  Input->GetParameters( "SplitEvent",        &split );

  //Input->GetParameters("Debug", &debugStr ) ; 
  //if ( debugStr == "True" ) debug = true ;
  for (int i=0; i<10; i++)  counter[i] = 0 ;

  Init();

}

Sync::~Sync(){


  delete Input ;
  delete select ;
  cout<<" done ! "<<endl ;

}

void Sync::Init() { 

   //vector<string> rFiles ;
   //Input->GetParameters( "TheData", &rFiles );
   //TTree* tr = Input->GetTree( rFiles[0],"EcalTimeAnalysis" );
   tr = Input->TreeMap( "data+" );

   tr->SetBranchAddress("eventId",    &eventId);
   tr->SetBranchAddress("trgCut",     &trgCut);
   tr->SetBranchAddress("nJets",      &nJets);
   tr->SetBranchAddress("nPhotons",   &nPhotons);
   tr->SetBranchAddress("nElectrons", &nElectrons);
   tr->SetBranchAddress("nVertices",  &nVertices );
   tr->SetBranchAddress("metPx",      &metPx );
   tr->SetBranchAddress("metPy",      &metPy );
   tr->SetBranchAddress("met",        &metE );

   tr->SetBranchAddress("jetPx",       jetPx );
   tr->SetBranchAddress("jetPy",       jetPy );
   tr->SetBranchAddress("jetPz",       jetPz );
   tr->SetBranchAddress("jetE",        jetE );
   tr->SetBranchAddress("jetNDau",     jetNDau );
   tr->SetBranchAddress("jetCM",       jetCM );
   tr->SetBranchAddress("jetCEF",      jetCEF );
   tr->SetBranchAddress("jetNHF",      jetNHF );
   tr->SetBranchAddress("jetNEF",      jetNEF );

   tr->SetBranchAddress("phoPx",       phoPx );
   tr->SetBranchAddress("phoPy",       phoPy );
   tr->SetBranchAddress("phoPz",       phoPz );
   tr->SetBranchAddress("phoE",        phoE );
   tr->SetBranchAddress("phoEcalIso",  phoEcalIso );
   tr->SetBranchAddress("phoHcalIso",  phoHcalIso );
   tr->SetBranchAddress("phoTrkIso",   phoTrkIso );
   tr->SetBranchAddress("phoHovE",     phoHovE );
   tr->SetBranchAddress("phoSmin",     phoSmin );
   tr->SetBranchAddress("phoSmaj",     phoSmaj );
   tr->SetBranchAddress("phoTime",     phoTime );

   tr->SetBranchAddress("elePx",        elePx );
   tr->SetBranchAddress("elePy",        elePy );
   tr->SetBranchAddress("elePz",        elePz );
   tr->SetBranchAddress("eleE",         eleE );
   tr->SetBranchAddress("eleEcalIso",   eleEcalIso );
   tr->SetBranchAddress("eleHcalIso",   eleHcalIso );
   tr->SetBranchAddress("eleTrkIso",    eleTrkIso );
   tr->SetBranchAddress("eleNLostHits", eleNLostHits );

   tr->SetBranchAddress("vtxX",       vtxX );
   tr->SetBranchAddress("vtxY",       vtxY );
   tr->SetBranchAddress("vtxZ",       vtxZ );
   tr->SetBranchAddress("vtxChi2",    vtxChi2 );
   tr->SetBranchAddress("vtxNdof",    vtxNdof );

}

// analysis template
void Sync::ReadTree() { 

   // link the variables for selection module
   select->Init( tr ) ;

   int totalN = tr->GetEntries();
   cout<<" total entries = "<< totalN <<" Process "<< ProcessEvents <<endl;
 
   for ( int i=0; i< tr->GetEntries(); i++ ) {
       if ( ProcessEvents > 0 && i > ( ProcessEvents - 1 ) ) break;
       tr->GetEntry( i );

       bool passEvent = true ;

       if ( i%2 == split ) continue ;    

       counter[0]++ ;   
       bool passPhoton = select->PhotonFilter( false );
       passEvent = ( passPhoton && passEvent ) ? true : false ;
       if ( passEvent ) counter[1]++ ;   

       bool passVertex = select->VertexFilter();
       passEvent = ( passVertex && passEvent ) ? true : false ;
       if ( passEvent ) counter[2]++ ;   

       select->ResetCuts( "PhotonCuts", 1, 1.4 ) ;
       bool passIso = select->PhotonFilter( true );
       passEvent = ( passIso && passEvent ) ? true : false ;
       if ( passEvent ) counter[3]++ ;   

       bool isGJets = select->GammaJetsBackground() ;
       passEvent = ( !isGJets && passEvent ) ? true : false ;
       if ( passEvent ) counter[4]++ ;   

       select->ResetCuts( "PhotonCuts", 0, 100. ) ;
       bool passTight = select->PhotonFilter( true );
       passEvent = ( passTight && passEvent ) ? true : false ;
       if ( passEvent ) counter[5]++ ;   

       bool passJet    = select->JetMETFilter();
       passEvent = ( passJet && passEvent ) ? true : false ;
       if ( passEvent ) counter[6]++ ;   
      
       bool passHLT    = select->HLTFilter();
       passEvent = ( passHLT && passEvent ) ? true : false ;
       if ( passEvent ) counter[7]++ ;   
       
       select->ResetCuts() ;
       select->ResetCollection() ;

   } // end of event looping

   cout<<" all:"<< counter[0] <<" pho:"<<counter[1]<<" vtx:"<< counter[2]<<" Iso:"<<counter[3]<<" noGjets:"<<counter[4] ;
   cout<<" tight:"<< counter[5] <<" jet:"<<counter[6]<<" hlt:"<<counter[7]<<endl; 

}

bool Sync::GammaJetsBackground( TLorentzVector gP4, vector<TLorentzVector> jP4s ) {

    bool pass = true ;

    if ( jP4s.size() < 1 ) return false ;

    double dPhi   = gP4.DeltaPhi( jP4s[0] ) ;
    double ratio1 = jP4s[0].Pt() / gP4.Pt() ;
    double ratio2 = ( jP4s.size() > 1 ) ? jP4s[1].Pt() / gP4.Pt() : 0. ;

    if ( dPhi <= (2*3.141592/3) )         pass = false ;
    if ( ratio1 >= 1.3 || ratio1 <= 0.7 ) pass = false ;
    if ( ratio2 >= 0.1 )                  pass = false ;

    return pass ;

}

