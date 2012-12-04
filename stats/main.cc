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

#include "AnaInput.h"
#include "StatsTool.h"

using namespace std; 

int main( int argc, const char* argv[] ) { 

  string datacardfile = ( argc > 1 ) ? argv[1] : "DataCard.txt";
  AnaInput  *Input = new AnaInput( datacardfile );

  // method to read root files
  Input->LinkForests("DP");

  int module = -1 ;
  Input->GetParameters( "Module", & module ) ;

  StatsTool  *statool  = new StatsTool( datacardfile ) ;
  if ( module == 0 ) {
     statool->Test();
  }
  if ( module == 1 ) {
     std::vector<string> testModel ;
     Input->GetParameters( "TestModel",  &testModel  ) ;

     for ( size_t i=0 ; i < testModel.size() ; i++ ) {
         statool->DPCLsTest( testModel[i] , (int)i );
     } 
  }
  if ( module == 2 ) {
     statool->SpecialPlot();
  }
  if ( module == 3 ) {
     statool->MyCLsTest();
  }
  if ( module == 4 ) {
     statool->DrawLimit();
  }
  if ( module == 5 ) {
     std::vector<string> testModel ;
     Input->GetParameters( "TestModel",  &testModel  ) ;

     for ( size_t i=0 ; i < testModel.size() ; i++ ) {
         statool->EffCalculator( testModel[i] , (int)i ) ;
     }
  }
  if ( module == 6 ) {
     std::vector<string> testModel ;
     Input->GetParameters( "TestModel",  &testModel  ) ;
     double pVal_b ;
     Input->GetParameters( "P_Value_B",  &pVal_b ) ;
     statool->DPCLsTest( testModel[0] , 0, pVal_b ) ;
  }

  delete statool ;

  cout<<" finished "<<endl ;

  delete Input ;

  return 0;

}


