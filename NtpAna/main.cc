#include <iostream> 
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <TString.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TMinuit.h>

#include "MathFunctions.h"
#include "AnaInput.h"
#include "PhotonAna.h" 
#include "Sync.h"
#include "timeVsAmpliCorrector.h"

using namespace std; 

int main( int argc, const char* argv[] ) { 

  string datacardfile = ( argc > 1 ) ? argv[1] : "DataCard.txt";
  AnaInput        *Input = new AnaInput( datacardfile );

  vector<int> module;
  Input->GetParameters( "Module", &module );

  Input->GetForest("TheData", "EcalTimeAnalysis");

  if ( module[0] == 1 ) {
     PhotonAna   *photonAna = new PhotonAna( datacardfile );
     photonAna->ReadTree() ;
     photonAna->ScalarPlotList() ;
  }
  if ( module[1] == 1 ) {
     Sync        *sync      = new Sync( datacardfile ) ;
     sync->ReadTree() ;
  }
  return 0;

}

