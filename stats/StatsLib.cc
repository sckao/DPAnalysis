#include "StatsLib.h"
#include <memory>

StatsLib::StatsLib( string datacardfile ) {

  Input  = new AnaInput( datacardfile );
  h_draw = new hDraw( datacardfile, "a" ) ;

  Input->GetParameters("PlotType",      &plotType ) ; 
  Input->GetParameters("Path",          &hfolder ) ; 
  Input->GetParameters("RandomSeed",    &randomSeed );

}

StatsLib::~StatsLib(){

  delete h_draw ;
  delete Input ;
  cout<<" done ! "<<endl ;

}

// ******************************************************************************
// Test Statistic Q ,  This is a one channel case
// xV : the set of measurements , contain n measurements 
// s : N. of expected signal , b: N. of expected background 
// hPs : signal PDF of the observable, hPb : background PDF of the observable
// hPs and hPb should be normalized to 1 !
// ******************************************************************************
double StatsLib::Q1( vector<double>& xV , double s, double b, TH1D* hPs, TH1D* hPb ) {

    double R = 0. ;
    for(size_t i=0; i < xV.size(); i++ ) {  
       
       double sS = s*hProb( xV[i], hPs ) ;
       double bB = b*hProb( xV[i], hPb ) ;
       double ri =  1. + ( sS / bB ) ;
       //cout<<" R:"<<R <<" ri:"<< ri <<" pS:"<< hProb( xV[i], hPs ) <<" pB:"<< hProb( xV[i], hPb ) <<endl ;
       R += log( ri ) ;
    }

    double Q = -2.0*( R - s )  ; 

    return Q ;
}

// ******************************************************************************
//  Histogram method ; an approximation 
// ******************************************************************************
double StatsLib::Q1( TH1D* hX , double s, double b, TH1D* hPs, TH1D* hPb ) {

    double R = 0. ;
    double w = hX->GetBinWidth(1) ;
    for( int i=0; i < hX->GetNbinsX(); i++ ) {  

       int n_X = hX->GetBinContent(i+1) ;
       double xV  = hX->GetBinCenter( i+1 ) ;
       double xVp = hX->GetBinCenter( i+1 ) + 0.25*w ;
       double xVn = hX->GetBinCenter( i+1 ) - 0.25*w ;

       double sS  = s*hProb( xV, hPs ) ;
       double bB  = b*hProb( xV, hPb ) ;
       double ri  =  1. + ( sS / bB ) ;

       int n_X1 = n_X / 4 ;  
       int n_X0 = n_X - (2*n_X1) ;  
       double sSp = s*hProb( xVp, hPs ) ;
       double bBp = b*hProb( xVp, hPb ) ;
       double rip =  1. + ( sSp / bBp ) ;

       double sSn = s*hProb( xVn, hPs ) ;
       double bBn = b*hProb( xVn, hPb ) ;
       double rin =  1. + ( sSn / bBn ) ;

       double R_i = ( n_X0*log( ri ) ) + ( n_X1*log( rip ) ) + ( n_X1*log( rin ) ) ;
       //cout<<" R:"<<R <<" ri:"<< ri <<" pS:"<< hProb( xV[i], hPs ) <<" pB:"<< hProb( xV[i], hPb ) <<endl ;
       //R += ( n_X*log( ri ) ) ;
       R += R_i  ;
    }

    double Q = -2.0*( R - s )  ; 

    return Q ;
}

// ******************************************************************************
// m : number of toyMC events for computing Q
// s : expected signal , b : expected background
// n : number of sample generated for PDF(Q) 
// ******************************************************************************
void StatsLib::pdfQ( TH1D* hPdfQ_sb, TH1D* hPdfQ_b, int m,  double s, double b, TH1D* hPs, TH1D* hPb, int n ) {

     // Reset the hPdfQ
     hPdfQ_sb->Reset() ;
     hPdfQ_b->Reset() ;
     double minQ_sb = hPdfQ_sb->GetBinCenter(1) ;
     double maxQ_sb = hPdfQ_sb->GetBinCenter( hPdfQ_sb->GetNbinsX() ) ;
     double minQ_b  = hPdfQ_b->GetBinCenter(1) ;
     double maxQ_b  = hPdfQ_b->GetBinCenter( hPdfQ_b->GetNbinsX() ) ;

     // build s+b model 
     TH1D* hPsb = (TH1D*) hPb->Clone() ;
     TH1D* hPs_ = (TH1D*) hPs->Clone() ; 
     hPs_->Scale( s/b ) ;
     hPsb->Add( hPs_ ) ;
     hPsb->Scale( 1./ hPsb->Integral() ) ;

    TRandom3* tRan = new TRandom3();
    tRan->SetSeed( randomSeed );

     vector<double>  bV ;
     vector<double> sbV ;
     for ( int j=0 ; j < n ; j++ ) {
           int rndmSeed1 = j*7 +  tRan->Rndm() ;
	   int rndmSeed2 = j*13 + tRan->Rndm() ;
           bV.clear() ;
           sbV.clear() ;

           // generate some toyMC in order to get Q distribution under the hypothesis b or s+b
	   bV  = toyMC( m,  hPb, rndmSeed1 );
	   sbV = toyMC( m, hPsb, rndmSeed2 );
	   double  qb = Q1(  bV, s, b, hPs, hPb ) ;
	   double qsb = Q1( sbV, s, b, hPs, hPb ) ;
	   //cout<<" sbV: "<< sbV.size() <<" bV: "<< bV.size() ;
	   //cout<<" qsb = "<< qsb <<" qb = "<< qb << endl ;
           if ( qb  < minQ_b  ) qb  = minQ_b ;
           if ( qb  > maxQ_b  ) qb  = maxQ_b ;
           if ( qsb < minQ_sb ) qsb = minQ_sb ;
           if ( qsb > maxQ_sb ) qsb = maxQ_sb ;
  
	   hPdfQ_b->Fill(qb) ;
	   hPdfQ_sb->Fill(qsb) ;
     }
     hPdfQ_b->Scale( 1./ hPdfQ_b->Integral() ) ;
     hPdfQ_sb->Scale( 1./ hPdfQ_sb->Integral() ) ;

     delete tRan ;
}

// ******************************************************************************
//  Histogram version of PDF_Q 
// ******************************************************************************
void StatsLib::hPdfQ( TH1D* hPdfQ_sb, TH1D* hPdfQ_b, int m,  double s, double b, TH1D* hPs, TH1D* hPb, int n ) {

     // Reset the hPdfQ
     hPdfQ_sb->Reset() ;
     hPdfQ_b->Reset() ;
     double minQ_sb = hPdfQ_sb->GetBinCenter(1) ;
     double maxQ_sb = hPdfQ_sb->GetBinCenter( hPdfQ_sb->GetNbinsX() ) ;
     double minQ_b  = hPdfQ_b->GetBinCenter(1) ;
     double maxQ_b  = hPdfQ_b->GetBinCenter( hPdfQ_b->GetNbinsX() ) ;

     // build s+b model 
     TH1D* hPsb = (TH1D*) hPb->Clone() ;
     TH1D* hPs_ = (TH1D*) hPs->Clone(); 
     hPs_->Scale( s/b ) ;
     hPsb->Add( hPs_ ) ;
     hPsb->Scale( 1./ hPsb->Integral() ) ;

     TH1D* hbV  = (TH1D*) hPb->Clone() ;
     TH1D* hsbV = (TH1D*) hPs->Clone(); 
     for ( int j=0 ; j < n ; j++ ) {
           hbV->Reset() ;
           hsbV->Reset() ;

           // generate some toyMC in order to get Q distribution under the hypothesis b or s+b
	   toyMC( hbV,  m,  hPb );
	   toyMC( hsbV, m, hPsb );
           /*
           if ( j < 10 ) {
              char hname_sb[15], hname_b[15] ;
              sprintf(hname_sb, "h_SB_toy_%02d", j ) ;
              sprintf(hname_b,  "h_B_toy_%02d", j ) ;
              h_draw->Draw(   hsbV, hname_sb, " x ", "", "logY", 0.95, 4 ) ;
              h_draw->Draw(    hbV, hname_b,  " x ", "", "logY", 0.95, 2 ) ;
              printf("(%d) sb: %f, b: %f \n" , j , hsbV->Integral( 29, 80), hbV->Integral(29, 80)  ) ;
           }
           */
	   double  qb = Q1(  hbV, s, b, hPs, hPb ) ;
	   double qsb = Q1( hsbV, s, b, hPs, hPb ) ;
	   //cout<<" sbV: "<< sbV.size() <<" bV: "<< bV.size() ;
	   //cout<<" qsb = "<< qsb <<" qb = "<< qb << endl ;
           if ( qb  < minQ_b  ) qb  = minQ_b ;
           if ( qb  > maxQ_b  ) qb  = maxQ_b ;
           if ( qsb < minQ_sb ) qsb = minQ_sb ;
           if ( qsb > maxQ_sb ) qsb = maxQ_sb ;
  
	   hPdfQ_b->Fill(qb) ;
	   hPdfQ_sb->Fill(qsb) ;
     }
     hPdfQ_b->Scale( 1./ hPdfQ_b->Integral() ) ;
     hPdfQ_sb->Scale( 1./ hPdfQ_sb->Integral() ) ;

}

// ******************************************************************************
//  Construct the desired pdfQ 
// ******************************************************************************
void StatsLib::pdfQ( string type, TH1D* hPdfQ, int m,  double s, double b, TH1D* hPs, TH1D* hPb, int n ) {

     // build s+b model 
     TH1D* hPx ; 
     if (  strncasecmp( "sb", type.c_str(), type.size() ) ==0 ) {
        hPx  = (TH1D*) hPb->Clone() ;
        TH1D* hPs_ = (TH1D*) hPs->Clone(); 
	hPs_->Scale( s/b ) ;
	hPx->Add( hPs_ ) ;
	hPx->Scale( 1./ hPx->Integral() ) ;
     }
     if (  strncasecmp( "b", type.c_str(), type.size() ) ==0 ) {
	hPx  = (TH1D*) hPb->Clone() ;
     }
     if (  strncasecmp( "s", type.c_str(), type.size() ) ==0 ) {
	hPx  = (TH1D*) hPs->Clone() ;
     }

     for ( int j=0 ; j < n ; j++ ) {
           int rndmSeed1 = j*7 +  j ;

           // generate some toyMC in order to get Q distribution under the hypothesis b or s+b
	   vector<double>  xV = toyMC( m,  hPb, rndmSeed1 );
	   double  qx = Q1(  xV, s, b, hPs, hPb ) ;

	   //cout<<" sbV: "<< sbV.size() <<" bV: "<< bV.size() ;
	   //cout<<" qsb = "<< qsb <<" qb = "<< qb << endl ;
	   hPdfQ->Fill(qx) ;
     }
     hPdfQ->Scale( 1./ hPdfQ->Integral() ) ;

}

// ******************************************************************************
// Calculate CL_s+b or CL_s ... according to hPdfQ
//  -> Integral( -inf  ~ q_obs ) pdf(Q|H) dq 
// ******************************************************************************
double StatsLib::CLx( TH1D* hPdfQ,  vector<double>& xV,  double s, double b, TH1D* hPs, TH1D* hPb ) {

     double q = Q1( xV, s, b, hPs, hPb ) ;

     double IntP = 0 ;
     int    ibin = 1 ;
     int    qbin = hPdfQ->FindBin( q ) ; 
     //cout<<" Q = "<< q <<" @ bin: "<< qbin ;
     do {

        double pi = hPdfQ->GetBinContent( ibin ) ;
        if ( ibin == qbin ) pi = pi / 2. ; 
        IntP += pi ;
        ibin++ ;
     } while ( ibin < qbin ) ;  
     //cout<<" stop @ "<< ibin <<endl ;

     return (1. - IntP) ;
}

// ******************************************************************************
//  Find the P-Value for S+B Model while P-Value for B model = pVal_b
// ******************************************************************************
double StatsLib::SBModelPValue( TH1D* hPdfQ_sb, TH1D* hPdfQ_b, double pVal_b ) {

     double IntP_b  = 0 ;
     double IntP_sb = 0 ;
     double last_pVal_b = 0 ;
     double last_pVal_sb = 0 ;
     int    ibin = 1 ;
     do {
        last_pVal_b  = IntP_b ;
        last_pVal_sb = IntP_sb ;
        double pi_b  = hPdfQ_b->GetBinContent( ibin ) ;
        double pi_sb = hPdfQ_sb->GetBinContent( ibin ) ;
        IntP_b  += pi_b ;
        IntP_sb += pi_sb ;
        ibin++ ;
     } while ( IntP_b < pVal_b ) ;

     double ratio = ( pVal_b - last_pVal_b) / ( IntP_b - last_pVal_b);
     IntP_sb = last_pVal_sb + ( ratio*(IntP_sb - last_pVal_sb) ) ;

     printf(" IntP_sb: %.3f / IntP_b : %.3f  --> CLs = %.3f \n", IntP_sb, IntP_b, (1 - IntP_sb)/(1- pVal_b) ) ;
     return  (1 - IntP_sb)/ ( 1- pVal_b );
}

// Calculate the p-value for S+B model as well as its +/- 1 sigam and 2 sigma
vector<double> StatsLib::SBModelPValue( TH1D* hPdfQ_sb, TH1D* hPdfQ_b ) {

     double flags[5] = { 0.0228, 0.1587, 0.5, 0.8413, 0.9772 } ;
     vector<double> plist ;

     double IntP_b  = 0 ;
     double IntP_sb = 0 ;
     double last_pVal_b = 0 ;
     double last_pVal_sb = 0 ;
     int    nbin = hPdfQ_b->GetNbinsX()  ;
     for ( int i=0; i < nbin; i++) {
         last_pVal_b  = IntP_b ;
         last_pVal_sb = IntP_sb ;
	 double pi_b  = hPdfQ_b->GetBinContent( i ) ;
	 double pi_sb = hPdfQ_sb->GetBinContent( i ) ;
	 IntP_b  += pi_b  ;
	 IntP_sb += pi_sb ;
         if ( IntP_b > flags[ plist.size() ] ) {
            double ratio = (IntP_sb - last_pVal_sb)/ ( IntP_b - last_pVal_b ) ;
            double thePVal = last_pVal_sb + ( ratio* ( flags[ plist.size() ] - last_pVal_b) ) ;
            double theCLs  = (1 - thePVal )/ ( 1- flags[ plist.size() ]  ) ;
            printf("@ p-Val-B: %.4f ,  IntP_sb: %3f / IntP_b : %3f -->  %3f \n", flags[plist.size()],  IntP_sb, IntP_b, theCLs ) ;
            plist.push_back( theCLs ) ;
         } 
         if ( plist.size() > 4 ) break ;
     } 

     return  plist ;
}


// ******************************************************************************
//  Find the expected S for desired Confidence Level 
//  alpha : significance for the limit 
//  pVal_b : 50% for nominal limit, 84.13% and 15.865% for +/- 1 sigma 
// ******************************************************************************
double StatsLib::ExpectedLimit( int nObs, double s, double b, TH1D* hPs, TH1D* hPb, double alpha, double pVal_b, int nToys, string subfolder ) {

     // create a folder to store histograms for your reference
     string origPath = h_draw->CurrentPath() ;
     if ( subfolder.size() > 1 )   {
        h_draw->SetNewPath( subfolder ) ;
     }

     // determine the scale of the PDF_Q's range
     //double ks = KSTest( hPs, hPb, 12, 80 , -1.5, 0.25 ) ;
     //double qs = ( ks > 0.99999 ) ? 0.8 : 1.8 ;
     double qs =  2. ;
     double qs1 = 2. ;

     char SBQbuffer[20] ;
     char BQbuffer[20] ;
     sprintf( SBQbuffer, "hPdfQ_sb_%.4f", pVal_b ) ;
     sprintf( BQbuffer,  "hPdfQ_b_%.4f",  pVal_b ) ;
     TH1D* hPdfQ_sb = new TH1D( SBQbuffer, " PDF(Q|H) ", 200,  -10, 10 ) ;
     TH1D* hPdfQ_b  = new TH1D( BQbuffer, " PDF(Q|Hb) ", 200,  -10, 10 ) ;

     double step = 0. ;
     bool   next = true ;
     double CL = 0 ;
     double expCL = 1. - alpha ;
     double expS = s ;
     double expB = b ;
     vector<double> CLsCollV ;
     vector<double> sNormV   ;

     bool badPDF = false ;
     double sS = expS * hPs->GetBinContent(hPs->GetMaximumBin()) / hPs->Integral() ;
     double bB = expB * hPb->GetBinContent(hPs->GetMaximumBin()) / hPb->Integral() ;
     double minPQ = (qs*expS) -  ( qs* nObs * log(  1.0 + (sS/bB) ) )  ;
     double maxPQ =  qs1*expS ;

     for ( int k=0; (k < 50) && next ; k++ ) {

        double pn = ( step < 0 ) ? -1. : 1. ;
        expS  = ( badPDF ) ?  expS : expS*( 1. + step )  ;
        minPQ = ( badPDF || k < 1 ) ? minPQ : minPQ - (pn*fabs( minPQ/sqrt(1.*nObs) )); 
        maxPQ = ( badPDF || k < 1 ) ? maxPQ : maxPQ + (pn*fabs( maxPQ/sqrt(1.*nObs) )) ;
        //minPQ = ( badPDF ) ? minPQ : (qs*expS) -  ( qs* nObs * log(  1.0 + (sS/bB) ) )  ; 
        //maxPQ = ( badPDF ) ? maxPQ : qs1*expS ;
        if ( maxPQ < minPQ ) {
           badPDF = true ;
           maxPQ = maxPQ + 0.5*fabs( maxPQ ) ; 
           continue ;
        } 
        cout<<"  *** ExpS = "<<  expS <<" minX: "<< minPQ <<" ~  maxX: "<< maxPQ << endl ;

        // Reset bin range and number
        hPdfQ_sb->SetBins( 500, minPQ, maxPQ ) ;
	hPdfQ_b->SetBins( 500, minPQ, maxPQ ) ;

	// Calculate PDF(Q|H) , H -> H_b or H_sb
        // re-do if the PDF_Q has too many over/under-flow 
        badPDF = false ;         
	pdfQ( hPdfQ_sb, hPdfQ_b, nObs , expS, expB , hPs, hPb, nToys );

	if ( hPdfQ_sb->GetBinContent(1) > pVal_b*hPdfQ_sb->Integral() ) badPDF = true ;
	if ( hPdfQ_b->GetBinContent(1) > pVal_b*hPdfQ_b->Integral() )   badPDF = true ;
        if ( badPDF ) {
           minPQ = minPQ - ( qs* nObs * log(  1.0 + (sS/bB) ) ) ;
           hPdfQ_sb->Reset() ;
           hPdfQ_b->Reset() ;
           continue ;
        }
	int nbin = hPdfQ_sb->GetNbinsX() ;
	if ( hPdfQ_sb->GetBinContent( nbin ) > (1.-pVal_b)*hPdfQ_sb->Integral() ) badPDF = true  ;
	if ( hPdfQ_b->GetBinContent( nbin ) > (1.-pVal_b)*hPdfQ_b->Integral() )   badPDF = true  ;
        if ( badPDF ) {
           maxPQ = maxPQ + 0.5*fabs( maxPQ ) ; 
           hPdfQ_sb->Reset() ;
           hPdfQ_b->Reset() ;
           continue ;
        }

        double pValue_SB = SBModelPValue( hPdfQ_sb, hPdfQ_b, pVal_b ) ;
        double CL_0 = 1 - pValue_SB ;
        cout<<" CL_0 = "<< CL_0 <<" step = "<< step <<endl ;

	// Draw the PdfQ distribution
	char hnameBuffer[15] ;
	sprintf(hnameBuffer, "PDF_Q_SB%02d", k ) ;
	h_draw->Draw(       hPdfQ_b, "", " Q ", "", "", 0.98, 1 ) ;
	h_draw->DrawAppend( hPdfQ_sb, hnameBuffer,  0.70, 2 ) ;

        if ( CLsCollV.size() < 1 ) {
           CL = CL_0 ;
           step = ( CL > expCL ) ? -0.3 : 0.3 ; 
        }
        if ( step > 0. && CL_0 > CL ) {
           CL = CL_0 ;
           if ( CL > expCL ) next = false ;
        } 
        if ( step < 0. && CL_0 < CL ) {
           CL = CL_0 ;
           if ( CL < expCL ) next = false ;
        }

        if ( CLsCollV.size() > 0 && step > 0 && CL_0 < expCL ) {
           double dCLs  = CL_0 - CLsCollV[ CLsCollV.size() -1 ] ;
           double dSg   = expS - sNormV[ sNormV.size() - 1 ] ;
           if ( fabs(dCLs) > 0.0001 ) { 
              double nextS = ( CL_0 < 0.8*expCL  ) ? ((0.9*expCL) - CL_0)*dSg / dCLs : (expCL - CL_0)*dSg / dCLs ;
              step = ( (nextS/expS) < 0.05 ) ? 0.05 : (nextS/expS) ;
           }
           if ( step > 1.5 || step < -0.3 ) step = 1.5 ;
        }
        if ( CLsCollV.size() > 0 && step < 0 && CL_0 > expCL ) {
           double dCLs  = CL_0 - CLsCollV[ CLsCollV.size() -1 ] ;
           double dSg   = expS - sNormV[ sNormV.size() - 1 ] ;
           if ( fabs(dCLs) > 0.0001 ) { 
               double nextS = ( expCL - CL_0 )*dSg / dCLs ;
               step = ( fabs(nextS/expS) < 0.1 ) ? -0.1 : (0.8*nextS/expS) ;
           } else if ( fabs(dCLs) < 0.0001 ) { 
              step = -0.5 ;
           }
           if ( step < -0.5 || step > 0.3 ) step = -0.4 ;
        }
        if ( expS > 10.*expB ) next = false ;
        cout <<" "<< hnameBuffer <<" expS : "<< expS << " CL : "<< CL <<" step : "<< step << endl ;
        cout <<" "<<endl ;
        CLsCollV.push_back( CL_0 );
        sNormV.push_back( expS ) ;
        // Draw the S+B model 
        if ( !next ) {
           TH1D* hPsb = (TH1D*) hPb->Clone() ;
	   TH1D* hPs_ = (TH1D*) hPs->Clone();
           TH1D* hPb_ = (TH1D*) hPb->Clone() ;
	   hPsb->SetTitle(" Signal+Background model");
	   hPs_->Scale( expS/expB ) ;
	   hPsb->Add( hPs_ ) ;
	   hPsb->Scale( nObs*1.0 / hPsb->Integral() ) ;
	   hPb_->Scale( nObs*1.0 ) ;
	   h_draw->Draw(  hPsb,     "", " x ", "", "logY", 0.95, 4 ) ;
	   h_draw->DrawAppend( hPb_, "hP_SB",  0.75, 2 ) ;
        }
     } 

     // Calculate the exact (1 - alpha) CL 
     const int sz = CLsCollV.size() ;
     double slope   = ( sNormV[sz-1] - sNormV[sz-2] ) / ( CLsCollV[sz-1] - CLsCollV[sz-2] ) ;
     double upLimit = slope*( 1-alpha - CLsCollV[sz-2] ) +  sNormV[sz-2]  ;
     cout<<" Up limit signal model: "<< upLimit << " @ " <<1 - alpha<< " CL "<< endl ;

     double CLsColl[sz]  ;
     double sNorm[sz]    ;
     for ( size_t i=0; i< CLsCollV.size() ; i++ ) {
         CLsColl[i] = CLsCollV[i] ;
         sNorm[i]   = sNormV[i] ;
     }

     // Draw the result
     TGraph* gCLs  = new TGraph( sz, sNorm, CLsColl );
     h_draw->DrawGraph( gCLs , "CLs",  " Expected Signal ", "CLs" ) ;

     hPdfQ_sb->Reset() ;
     hPdfQ_b->Reset() ;
     if ( subfolder.size() > 1 )  {
        h_draw->SetNewPath( origPath ) ;
     }

     hPdfQ_sb->Reset() ;
     hPdfQ_b->Reset() ;

     delete hPdfQ_sb ;
     delete hPdfQ_b ;

     return upLimit ;
}

// ******************************************************************************
//  Histogram version of limit setting
// ******************************************************************************
double StatsLib::ExpectedLimit1( int nObs, double s, double b, TH1D* hPs, TH1D* hPb, double alpha, double pVal_b, int nToys, string subfolder ) {

     // create a folder to store histograms for your reference
     string origPath = h_draw->CurrentPath() ;
     if ( subfolder.size() > 1 )   {
        h_draw->SetNewPath( subfolder ) ;
     }

     // determine the scale of the PDF_Q's range
     //double ks = KSTest( hPs, hPb, 12, 80 , -1.5, 0.25 ) ;
     //double qs = ( ks > 0.99999 ) ? 0.8 : 1.8 ;
     double qs  = 2. ;
     double qs1 = 1.8 ;

     char SBQbuffer[20] ;
     char BQbuffer[20] ;
     sprintf( SBQbuffer, "hPdfQ_sb_%.4f", pVal_b ) ;
     sprintf( BQbuffer,  "hPdfQ_b_%.4f",  pVal_b ) ;
     TH1D* hPdfQ_sb = new TH1D( SBQbuffer, " PDF(Q|H) ", 200,  -10, 10 ) ;
     TH1D* hPdfQ_b  = new TH1D( BQbuffer, " PDF(Q|Hb) ", 200,  -10, 10 ) ;

     double step = 0 ;
     bool   next = true ;
     double CL = 0 ;
     double expCL = 1. - alpha ;
     double expS = s ;
     double expB = b ;
     vector<double> CLsCollV ;
     vector<double> sNormV   ;

     bool badPDF = false ;
     double sS = expS * hPs->GetBinContent(hPs->GetMaximumBin()) / hPs->Integral() ;
     double bB = expB * hPb->GetBinContent(hPs->GetMaximumBin()) / hPb->Integral() ;
     double minPQ = (qs*expS) -  ( qs* nObs* log(  1.0 + (sS/bB) ) )  ;
     double maxPQ =  qs1*expS ;
     double newMinPQ = minPQ ;
     double newMaxPQ = maxPQ ;

     for ( int k=0; (k < 50) && next ; k++ ) {
 
        // PDF_Q's range
        //double pn = ( step < 0. ) ? -1. : 1. ;
        expS  = ( badPDF ) ?  expS : expS + step   ;
        minPQ = ( badPDF || k < 1 ) ? minPQ : newMinPQ ; 
        maxPQ = ( badPDF || k < 1 ) ? maxPQ : newMaxPQ ;
        cout<<"  *** ExpS = "<<  expS <<" minX: "<< minPQ <<" ~  maxX: "<< maxPQ <<" step: "<< step <<endl ;
        if ( maxPQ < minPQ ) {
           badPDF = true ;
           maxPQ = maxPQ + 0.5*fabs( maxPQ ) ; 
           continue ;
        } 

        // Reset bin range and number
        hPdfQ_sb->SetBins( 200, minPQ, maxPQ ) ;
	hPdfQ_b->SetBins( 200, minPQ, maxPQ ) ;

	// Calculate PDF(Q|H) , H -> H_b or H_sb
        // re-do if the PDF_Q has too many over/under-flow 
        badPDF = false ;
	hPdfQ( hPdfQ_sb, hPdfQ_b, nObs , expS, expB , hPs, hPb, nToys );
	if ( hPdfQ_sb->GetBinContent(1) > pVal_b*hPdfQ_sb->Integral() ) badPDF = true ;
	if ( hPdfQ_b->GetBinContent(1)  > pVal_b*hPdfQ_b->Integral() ) badPDF = true ;
        if ( badPDF ) {
           minPQ = minPQ - ( nObs * log(  1.0 + (sS/bB) ) ) ;
           hPdfQ_sb->Reset() ;
           hPdfQ_b->Reset() ;
           continue ;
        }
	int nbin = hPdfQ_sb->GetNbinsX() ;
	if ( hPdfQ_sb->GetBinContent( nbin ) > (1-pVal_b)*hPdfQ_sb->Integral() ) badPDF = true  ;
	if ( hPdfQ_b->GetBinContent( nbin )  > (1-pVal_b)*hPdfQ_b->Integral() )   badPDF = true  ;
        if ( badPDF ) {
	   maxPQ = maxPQ + 0.5*fabs( maxPQ ) ; 
           hPdfQ_sb->Reset() ;
           hPdfQ_b->Reset() ;
           continue ;
        }
        // Adjust the range of PDFQ for next step
	cout<<" Valid PdfQ Range ...  the new range is ... " <<endl ;
	newMinPQ = hPdfQ_sb->GetMean() - (5.*hPdfQ_sb->GetRMS()) ;
	newMaxPQ =  hPdfQ_b->GetMean() + (5.*hPdfQ_b->GetRMS() ) ;
	if ( hPdfQ_sb->GetMaximumBin() < 2 || hPdfQ_b->GetMaximumBin() < 2 ) newMinPQ -= 3.*hPdfQ_b->GetRMS() ;
	if ( hPdfQ_sb->GetMaximumBin() > nbin-1 || hPdfQ_b->GetMaximumBin() > nbin-1 ) newMaxPQ += 3.*hPdfQ_sb->GetRMS() ;
	printf("   min = %.2f - 5x %.2f \n", hPdfQ_sb->GetMean() , hPdfQ_sb->GetRMS() ) ;
	printf("   max = %.2f + 5x %.2f \n", hPdfQ_b->GetMean()  , hPdfQ_b->GetRMS()  ) ;

        // Calculating the CLs value
        double pValue_SB = SBModelPValue( hPdfQ_sb, hPdfQ_b, pVal_b ) ;
        double CL_0 = 1 - pValue_SB ;

	// Draw the PdfQ distribution
	char hnameBuffer[15] ;
	sprintf(hnameBuffer, "PDF_Q_SB%02d", k ) ;
	h_draw->Draw(       hPdfQ_b, "", " Q ", "", "", 0.98, 1 ) ;
	h_draw->DrawAppend( hPdfQ_sb, hnameBuffer,  0.72, 2 ) ;

        next = CLPropagator( expS, expB, expCL, CL_0, CL, step, CLsCollV, sNormV ) ;
        cout<<" CL_0 = "<< CL_0 <<" next step = "<< step <<endl ;

        CLsCollV.push_back( CL_0 );
        sNormV.push_back( expS ) ;
        // Draw the S+B model 
        if ( !next ) {
           TH1D* hPsb = (TH1D*) hPb->Clone() ;
	   TH1D* hPs_ = (TH1D*) hPs->Clone();
           TH1D* hPb_ = (TH1D*) hPb->Clone() ;
	   hPsb->SetTitle(" Signal+Background model");
	   hPs_->Scale( expS/expB ) ;
	   hPsb->Add( hPs_ ) ;
	   hPsb->Scale( nObs*1.0 / hPsb->Integral() ) ;
	   hPb_->Scale( nObs*1.0 ) ;
	   h_draw->Draw(  hPsb,     "", " x ", "", "logY", 0.95, 4 ) ;
	   h_draw->DrawAppend( hPb_, "hP_SB",  0.75, 2 ) ;
        }

     }

     // Calculate the exact (1 - alpha) CL 
     double upLimit = Interpolate( CLsCollV, sNormV, 1. - alpha ) ;
     cout<<" Up limit signal model: "<< upLimit << " @ " <<1 - alpha<< " CL "<< endl ;

     // Reset PDFQ and path
     hPdfQ_sb->Reset() ;
     hPdfQ_b->Reset() ;
     if ( subfolder.size() > 1 )  {
        h_draw->SetNewPath( origPath ) ;
     }

     delete hPdfQ_sb ;
     delete hPdfQ_b ;

     return upLimit ;
}

// ******************************************************************************
//  Histogram version of limit setting - this is used for scanning expected value
// ******************************************************************************
vector<double> StatsLib::ExpectedLimit2( int nObs, double s, double b, TH1D* hPs, TH1D* hPb, double alpha, int nToys, string subfolder ) {

     // create a folder to store histograms for your reference
     string origPath = h_draw->CurrentPath() ;
     if ( subfolder.size() > 1 ) {
        h_draw->SetNewPath( subfolder ) ;
     }

     // determine the scale of the PDF_Q's range
     //double ks = KSTest( hPs, hPb, 12, 80 , -1.5, 0.25 ) ;
     //double qs = ( ks > 0.99999 ) ? 0.8 : 1.8 ;

     TH1D* hPdfQ_sb = new TH1D( "hPdfQ_sb_test", " PDF(Q|H) ", 200,  -10, 10 ) ;
     TH1D* hPdfQ_b  = new TH1D( "hPdfQ_b_test", " PDF(Q|Hb) ", 200,  -10, 10 ) ;

     bool   next[5]  = {   true,   true, true,   true, true } ;
     double flags[5] = { 0.0228, 0.1587,  0.5, 0.8413, 0.9772 } ;
     double CL    = 0. ;
     double expCL = 1. - alpha ;
     double expS = 1  ;
     double expB = b ;
     double step = sqrt( s ) ;
     vector<double> upLimitV ;
     vector<double> CLsV, CLsU1V, CLsU2V, CLsD1V, CLsD2V ;
     vector<double> esV  ;

     bool badPDF = false ;
     double sS = expS * hPs->GetBinContent(hPs->GetMaximumBin()) / hPs->Integral() ;
     double bB = expB * hPb->GetBinContent(hPs->GetMaximumBin()) / hPb->Integral() ;
     double minPQ = (2*expS) -  ( 3 * nObs * log( 1.0 + (sS/bB) ) )  ;
     double maxPQ =  2*expS ;
     double newMinPQ = minPQ ;
     double newMaxPQ = maxPQ ;

     for ( int k=0; (k < 50) && next[0] ; k++ ) {
 
        // PDF_Q's range
        expS  = ( badPDF || k < 1 ) ?  expS : expS + step   ; 
        minPQ = ( badPDF || k < 1 ) ? minPQ : newMinPQ ; 
        maxPQ = ( badPDF || k < 1 ) ? maxPQ : newMaxPQ ;
        cout<<"  *** ExpS = "<<  expS <<" minX: "<< minPQ <<" ~  maxX: "<< maxPQ <<" step :"<< step<< endl ;
        if ( maxPQ <= minPQ ) {
           badPDF = true ;
           maxPQ = maxPQ + 0.5*fabs( maxPQ ) ; 
           continue ;
        } 

        // Reset bin range and number
        hPdfQ_sb->SetBins( 200, minPQ, maxPQ ) ;
	hPdfQ_b->SetBins( 200, minPQ, maxPQ ) ;

        // Creating PDFQ
	hPdfQ( hPdfQ_sb, hPdfQ_b, nObs , expS, expB , hPs, hPb, nToys );

	// Calculate PDF(Q|H) , H -> H_b or H_sb
        // re-do if the PDF_Q has too many over/under-flow
        double bound = 0.5 ;
        for ( int j=4; j>= 0; j-- ) {
            if ( next[j] ) {
               bound = flags[j] ; 
               break ;
            }
        }
        badPDF = false ;
	if ( hPdfQ_sb->GetBinContent(1) > bound*hPdfQ_sb->Integral() ) badPDF = true ;
	if ( hPdfQ_b->GetBinContent(1) >  bound*hPdfQ_b->Integral() ) badPDF = true ;
        if ( badPDF ) {
           minPQ = minPQ - ( nObs * log(  1.0 + (sS/bB) ) ) ;
           hPdfQ_sb->Reset() ;
           hPdfQ_b->Reset() ;
           continue ;
        }
	int nbin = hPdfQ_sb->GetNbinsX() ;
	if ( hPdfQ_sb->GetBinContent( nbin ) > (1 - bound)*hPdfQ_sb->Integral() ) badPDF = true  ;
	if ( hPdfQ_b->GetBinContent( nbin ) >  (1 - bound)*hPdfQ_b->Integral() )   badPDF = true  ;
        if ( badPDF ) {
	   maxPQ = maxPQ + 0.5*fabs( maxPQ ) ; 
           hPdfQ_sb->Reset() ;
           hPdfQ_b->Reset() ;
           continue ;
        }

        // Adjust the range of PDFQ , skip the first step 
        if ( k > 0 ) { 
           cout<<" Valid PdfQ Range ..." << bound <<"  the new range is ... " <<endl ;
           newMinPQ = hPdfQ_sb->GetMean() - (5.*hPdfQ_sb->GetRMS()) ;
	   newMaxPQ =  hPdfQ_b->GetMean() + (5.*hPdfQ_b->GetRMS() ) ;
           if ( hPdfQ_sb->GetMaximumBin() < 2 || hPdfQ_b->GetMaximumBin() < 2 ) newMinPQ -= 3.*hPdfQ_b->GetRMS() ;
           if ( hPdfQ_sb->GetMaximumBin() > nbin-1 || hPdfQ_b->GetMaximumBin() > nbin-1 ) newMaxPQ += 3.*hPdfQ_sb->GetRMS() ;
	   printf("   min = %.2f - 5x %.2f \n", hPdfQ_sb->GetMean() , hPdfQ_sb->GetRMS() ) ;
	   printf("   max = %.2f + 5x %.2f \n", hPdfQ_b->GetMean()  , hPdfQ_b->GetRMS()  ) ;
        } else {
           newMinPQ -= 3.*hPdfQ_b->GetRMS() ;
           newMaxPQ += 3.*hPdfQ_sb->GetRMS() ;
        }

        // Calculate the CLs values
        vector<double> pValue_SB = SBModelPValue( hPdfQ_sb, hPdfQ_b ) ;

	// Draw the PdfQ distribution
	char hnameBuffer[15] ;
	sprintf(hnameBuffer, "PDF_Q_SB%02d", k ) ;
	h_draw->Draw(       hPdfQ_b, "", " Q ", "", "", 0.98, 1 ) ;
	h_draw->DrawAppend( hPdfQ_sb, hnameBuffer,  0.72, 2 ) ;

        // Propogate CL_0
        if ( next[4] ) {
           next[4] = CLPropagator( expS, expB, expCL, 1 - pValue_SB[4], CL, step, CLsD2V, esV ) ;
	   CLsD2V.push_back( 1-pValue_SB[4] );
	   CLsD1V.push_back( 1-pValue_SB[3] );
	   CLsV.push_back(   1-pValue_SB[2] );
	   CLsU1V.push_back( 1-pValue_SB[1] );
	   CLsU2V.push_back( 1-pValue_SB[0] );
        } else if ( next[3] ) {
           cout<<" > go to D1 " << endl ;
           CL = CLsD1V[ CLsD1V.size() -1 ] ; 
           next[3] = CLPropagator( expS, expB, expCL, 1 - pValue_SB[3], CL, step, CLsD1V, esV ) ;
	   CLsD1V.push_back( 1-pValue_SB[3] );
	   CLsV.push_back(   1-pValue_SB[2] );
	   CLsU1V.push_back( 1-pValue_SB[1] );
	   CLsU2V.push_back( 1-pValue_SB[0] );
        } else if ( next[2] ) {
           cout<<" > go to center " << endl ;
           CL = CLsV[ CLsV.size() -1 ] ; 
           next[2] = CLPropagator( expS, expB, expCL, 1 - pValue_SB[2], CL, step, CLsV,   esV ) ;
	   CLsV.push_back(   1-pValue_SB[2] );
	   CLsU1V.push_back( 1-pValue_SB[1] );
	   CLsU2V.push_back( 1-pValue_SB[0] );
        } else if ( next[1] ) {
           cout<<" > go to U1 " << endl ;
           CL = CLsU1V[ CLsU1V.size() -1 ] ; 
           next[1] = CLPropagator( expS, expB, expCL, 1 - pValue_SB[1], CL, step, CLsU1V, esV ) ;
	   CLsU1V.push_back( 1-pValue_SB[1] );
	   CLsU2V.push_back( 1-pValue_SB[0] );
        } else if ( next[0] ) {
           cout<<" > go to U2 " << endl ;
           CL = CLsU2V[ CLsU2V.size() -1 ] ; 
           next[0] = CLPropagator( expS, expB, expCL, 1 - pValue_SB[0], CL, step, CLsU2V, esV ) ;
	   CLsU2V.push_back( 1-pValue_SB[0] );
        }

        esV.push_back( expS ) ;
        cout <<" "<< hnameBuffer <<" expS : "<< expS << " CL : "<< CL <<" step : "<< step << endl ;
        cout <<" "<< endl ;

        if ( !next[0] ) {
           // Calculate the exact (1 - alpha) CL 
           double upLimit_d2 = Interpolate( CLsD2V, esV, 1-alpha );
           double upLimit_d1 = Interpolate( CLsD1V, esV, 1-alpha );
           double upLimit_0  = Interpolate( CLsV,   esV, 1-alpha );
           double upLimit_u1 = Interpolate( CLsU1V, esV, 1-alpha );
           double upLimit_u2 = Interpolate( CLsU2V, esV, 1-alpha );
           printf(" ========= %.2f Up Limit Report ================ \n", 1- alpha  ) ;
           printf("  - %.2f  - %.2f  %.2f  +%.2f  +%.2f \n", upLimit_d2, upLimit_d1, upLimit_0, upLimit_u1, upLimit_u2 ) ;
           printf(" =============================================== \n ");
           upLimitV.push_back( upLimit_d2 ) ;
           upLimitV.push_back( upLimit_d1 ) ;
           upLimitV.push_back( upLimit_0 ) ;
           upLimitV.push_back( upLimit_u1 ) ;
           upLimitV.push_back( upLimit_u2 ) ;

           // Draw the S+B model 
           TH1D* hPsb = (TH1D*) hPb->Clone() ;
	   TH1D* hPs_ = (TH1D*) hPs->Clone();
           TH1D* hPb_ = (TH1D*) hPb->Clone() ;
	   hPsb->SetTitle(" Signal+Background model");
	   hPs_->Scale( expS/expB ) ;
	   hPsb->Add( hPs_ ) ;
	   hPsb->Scale( nObs*1.0 / hPsb->Integral() ) ;
	   hPb_->Scale( nObs*1.0 ) ;
	   h_draw->Draw(  hPsb,     "", " x ", "", "logY", 0.95, 4 ) ;
	   h_draw->DrawAppend( hPb_, "hP_SB",  0.75, 2 ) ;
        }
     }


     hPdfQ_sb->Reset() ;
     hPdfQ_b->Reset() ;
     if ( subfolder.size() > 1 )  {
        h_draw->SetNewPath( origPath ) ;
     }

     delete hPdfQ_sb ;
     delete hPdfQ_b ;

     return upLimitV ;
}

bool StatsLib::CLPropagator( double expS, double expB, double expCL, double CL_0, double& CL, double& step, vector<double> CLsCollV, vector<double> sNormV ) {

        bool next = true ;
        if ( CLsCollV.size() < 1 ) {
           CL = CL_0 ;
           step = ( CL_0 > expCL ) ? -1*sqrt( expS ) : sqrt( expS ) ; 
           if ( CL_0 < 0.5 ) step = 2.*sqrt( expS ) ;
        }
        // across target CL value from small value
        if ( step > 0. && CL_0 > CL ) {
           CL = CL_0 ;
           if ( CL > expCL ) next = false ;
           if ( CL > expCL ) step = sqrt( expS ) ;
        } 
        // across target CL value from large value
        if ( step < 0. && CL_0 < CL ) {
           CL = CL_0 ;
           if ( CL < expCL ) next = false ;
           if ( CL < expCL ) step = -1*sqrt( expS ) ;
        }
        if ( CLsCollV.size() > 0 && step > 0 && CL_0 < expCL ) {
           double dCLs  = CL_0 - CLsCollV[ CLsCollV.size() -1 ] ;
           double dSg   = expS - sNormV[ sNormV.size() - 1 ] ;
           if ( fabs(dCLs) >= 0.001 ) { 
              double nextS = ( CL_0 < 0.8*expCL  ) ? ((0.9*expCL) - CL_0)*dSg / dCLs : (expCL - CL_0)*dSg / dCLs;
              step = ( nextS < 1 ) ? 1 : nextS ;
           } else if ( fabs(dCLs) < 0.001 && CL_0 >= 0.9*expCL ) {
              step = sqrt( expS ) ;
           } else if ( fabs(dCLs) < 0.001 && CL_0 < 0.9*expCL ) {
              step = 2.*sqrt( expS ) ;
           }

           if ( step > 5*sqrt( expS ) || step < 0 ) step = 5.*sqrt( expS ) ;
        }
        if ( CLsCollV.size() > 0 && step < 0 && CL_0 > expCL ) {
           double dCLs  = CL_0 - CLsCollV[ CLsCollV.size() -1 ] ;
           double dSg   = expS - sNormV[ sNormV.size() - 1 ] ;
           if ( fabs(dCLs) > 0.001 ) { 
              double nextS = ( expCL - CL_0 )*dSg / dCLs ;
              step = ( fabs(nextS) < 1 ) ? -1. : nextS ;
           } else if ( fabs(dCLs) < 0.001 && CL_0 >= 0.999 ) { 
              step = -3*sqrt( expS ) ;
           } else if ( fabs(dCLs) < 0.001 && CL_0 < 0.999 ) { 
              step = -1*sqrt( expS ) ;
           }
           if ( step < -3*sqrt(expS) || step > 0 ) step = -3*sqrt(expS) ;
        }

        if ( expS > 10.*sqrt(expB) + expB ) next = false ;
        if ( expS < 0 ) next = false ;

        return next  ;
}

// *******************************************************************
// Calculate the interpolate value for given interX. 
// If X is not in the range of xV, a projective value will be given 
// *******************************************************************
double StatsLib::Interpolate( vector<double> xV, vector<double> yV, double interX ) {

     const int sz   = xV.size() ;
     bool increasing = (  xV[sz -1] > interX && interX > xV[sz-2] ) ;
     bool decreasing = (  xV[sz -1] < interX && interX < xV[sz-2] ) ;

     double slope  = 0 ;
     double interY = 0 ;
     if ( increasing || decreasing ) {
        cout<<" using last interval ... "<<endl ;
        slope = ( yV[sz-1] - yV[sz-2] ) / ( xV[sz-1] - xV[sz-2] ) ;
        interY = slope*( interX - xV[sz -2 ] ) + yV[sz-2] ;
     } else {
        cout<<"searching the interval ... "<<endl ;
	double dxL = -9999 ;
	double dxH =  9999 ;
	int il = 0 ;
	int ih = sz-1 ;
        for ( int i =sz-1; i >=0  ; i-- ) {
            double dx_ = xV[i] - interX ;
            if ( dx_ < 0 && fabs(dx_) < fabs(dxL) ) {
               dxL = dx_ ;
               il = i ;
            }
            if ( dx_ > 0 && fabs(dx_) < fabs(dxH) ) {
               dxH = dx_ ;
               ih = i ;
            }
        }
        slope = ( yV[ ih ] - yV[ il ] ) / ( xV[ ih ] - xV[ il ] ) ;
        interY = slope*( interX - xV[il] ) + yV[il] ;
     }
     
      // Draw the result
      double xA[sz], yA[sz]  ;
      for ( int i=0; i< sz ; i++ ) {
          xA[i] = xV[i] ;
          yA[i] = yV[i] ;
      }
      char gnameBuffer[20] ;
      sprintf( gnameBuffer, "CLs_%05d", (int)interY ) ;
      TGraph* gCLs  = new TGraph( sz, yA, xA );
      h_draw->DrawGraph( gCLs , gnameBuffer ,  " Expected Signal ", "CLs" ) ;
 
      return interY ;
}

// ******************************************************************************
//  Generate fake sample, within range 0 ~ 10
// ******************************************************************************
void StatsLib::Sample( TH1D* hSample, string pdfName, int n , double mean, double sigma  ) {

    TRandom3* tRan = new TRandom3();
    tRan->SetSeed( randomSeed );

    for ( int i=0; i< n; i++) {

        double val = 0. ;
        if ( strncasecmp( "Gaus", pdfName.c_str(), pdfName.size() ) ==0 ) val = tRan->Gaus( mean , sigma ) ;
        if ( strncasecmp( "Pois", pdfName.c_str(), pdfName.size() ) ==0 ) val = tRan->Poisson( mean ) ;
        if ( strncasecmp( "Exp", pdfName.c_str(), pdfName.size() ) ==0 ) val = tRan->Exp( mean ) ;
        if ( strncasecmp( "Rndm", pdfName.c_str(), pdfName.size() ) ==0 ) val = tRan->Rndm() ;

        if ( val <  0 ) val =  0 ;
        if ( val > 10 ) val = 10 ;
        hSample->Fill( val ) ;
    } 

}

// ******************************************************************************
// Generate toy MC
// hPDF : PDF of the measurement , 
// n : desired number of events generated 
// ******************************************************************************
vector<double> StatsLib::toyMC( int n, TH1D* hPDF, int si ) {

  // set up the random number function
  TRandom3* tRan = new TRandom3();
  tRan->SetSeed( randomSeed + si );

  vector<double> xV ;
  for ( int i= 1; i<= n ; i++) {
      double xi = hPDF->GetRandom() ;
      xV.push_back( xi ) ;
  }

  delete tRan ;
  return xV ;

}

void StatsLib::toyMC( TH1D* hOut, int n, TH1D* hPDF, bool debug ) {

  // set up the random number function
  //TRandom3* tRan = new TRandom3();

  hOut->Reset() ;
  hOut->FillRandom( hPDF, n ) ;
  if ( debug ) {
     for ( int i=1 ; i <= hOut->GetNbinsX() ; i++ ) { 
         int ni_ = hOut->GetBinContent( i ) ;
	 float pi_ = hPDF->GetBinContent( i ) ;
	 printf(" (%d) ni: %d * pi: %f \n ", i, ni_ , pi_ ) ;
     }
     cout<<" h_Integral : "<< hOut->Integral() << endl ;
  }

}

// ******************************************************************************
// Take the histogram as PDF, return the probability of input x 
// ******************************************************************************
double StatsLib::hProb( double x, TH1D* hPDF ) {

    TH1D* hPDF_ = (TH1D*) hPDF->Clone() ;    
    if ( hPDF_->Integral() > 1.00001 || hPDF_->Integral() < 0.99999 ) {
       double norm = hPDF_->Integral() ;
       hPDF_->Scale( 1./ norm ) ;
    }

    int xbin = hPDF_->FindBin( x ) ;
    double prob = hPDF_->GetBinContent( xbin ) ;
    //cout<<" x : "<< x <<" bin : "<< xbin <<" cont:"<< prob <<endl ;

    delete hPDF_ ;
    return prob ;
}

TH1D* StatsLib::SubHistogram( TH1D* h1, double minX, double maxX ) {

     char hname1[20] ;
     sprintf( hname1,  "%s_s",  h1->GetName() ) ;
     int minBin = h1->FindBin( minX ) ;
     int maxBin = h1->FindBin( maxX ) ;
     double minBinX = h1->GetBinCenter( minBin ) - (0.5*h1->GetBinWidth( minBin ) ) ;
     double maxBinX = h1->GetBinCenter( maxBin ) + (0.5*h1->GetBinWidth( maxBin ) ) ;
     int nbin = maxBin - minBin + 1 ;
     //printf("  h1 %d -> %d :  x( %.3f ~ %.3f ) ", minBin, maxBin , minBinX, maxBinX ) ;

     TH1D* hsub = new TH1D( hname1, "", nbin , minBinX, maxBinX ) ;

     for ( int i = 1 ; i <= nbin ; i++ ) {
         hsub->SetBinContent( i, h1->GetBinContent( i + minBin - 1 ) ) ;
     }
     //printf(" x( %.3f ~ %.3f ) with %d bins => %.3f \n ",  minBinX, maxBinX, nbin, hsub->GetBinWidth(1) ) ;
 
     return hsub ;
}


// Convolute histogram with gaussian distribution
TH1D* StatsLib::SmearHistogram( TH1D* h1, double resol ) {

      int nBins = h1->GetNbinsX() ;
      vector<double> bV ;

      TH1D* h1S = (TH1D*) h1->Clone() ; 
      h1S->Reset() ;

      for ( int i=1; i<= nBins ; i++ ) {
          
          bV.clear() ;
          double xx = h1->GetBinCenter(i) ;
          double bw = h1->GetBinWidth(i) ;
          double bc = h1->GetBinContent(i) ;
          if ( bc < 0.0000001 ) {
             h1S->SetBinContent( i , bc ) ;
             continue ;
          }

          bV = GausSmear( xx, resol, bw ) ;
          int sz = bV.size() ;
          for ( int j = 0 ; j < sz; j++ ) {
              if ( j ==0 ) { 
                 double bc_ = h1S->GetBinContent(i) + ( bV[j]*2.*bc ) ;
                 h1S->SetBinContent( i, bc_ ) ;
              } else {
                 int ui = ( i + j > nBins ) ? nBins : i + j ;
                 int di = ( i - j <   1   ) ?   1   : i - j ;
                 double bc_u = h1S->GetBinContent(ui) + ( bV[j]*bc ) ;
                 double bc_d = h1S->GetBinContent(di) + ( bV[j]*bc ) ;
                 h1S->SetBinContent( ui , bc_u ) ;
                 h1S->SetBinContent( di , bc_d ) ;
              }
          }
      }

      return h1S ;
}

// Get the integrated gaussian probability for each bin with given bin width and sigam
vector<double> StatsLib::GausSmear( double mean, double sigma, double bW ) {

     vector<double> bV ;
     double step = 3.*sigma / 1000. ;
     double intP = 0. ;
     double tailP = 0.5 ;
     double x_   = mean + (step/2.) ; 
     double bound = mean + (bW/2.) ;
     int    nbins = 1 ;
     for ( int j=0 ; j<1000; j++) {
         double pl = step * TMath::Gaus( x_ , mean, sigma, true );
         x_ = x_ + step ;
         if ( x_ > bound ) {
            bV.push_back( intP ) ;
            tailP -= intP ;
            bound += bW ;
            nbins++ ;
            intP = 0 ;
         }
         intP += pl ;
     }
     bV.push_back( tailP ) ;

     return bV ;

}

double StatsLib::KSTest( TH1D* h1, TH1D* h2, int minBin, int maxBin, double hminX, double binWidth ) {

     int nBin = maxBin - minBin + 1 ;
     double minX = hminX + (minBin*binWidth) ;
     double maxX = minX + (binWidth*nBin) ;

     char hname1[20], hname2[20] ;
     sprintf( hname1,  "%s_ks",  h1->GetName() ) ;
     sprintf( hname2,  "%s_ks",  h2->GetName() ) ;

     TH1D* h1t = new TH1D( hname1, "KS Test 1" , nBin,  minX, maxX );
     TH1D* h2t = new TH1D( hname2, "KS Test 2" , nBin,  minX, maxX );
     for ( int i = 0 ; i < nBin ; i++) {
         h1t->SetBinContent( i+1 , h1->GetBinContent( minBin + i  ) ) ;
         h2t->SetBinContent( i+1 , h2->GetBinContent( minBin + i  ) ) ;
     }
     //h2t->Scale( h1t->Integral() / h2t->Integral() ) ;

     double ksVal = h1t->KolmogorovTest( h2t ) ;
     //TLegend* leg3  = new TLegend(.5, .7, .77, .9 );
     //leg3->AddEntry( h1t, " 2012B ", "L" ) ;
     //leg3->AddEntry( h2t, " 2012C ", "L" ) ;

     //h_draw->Draw(       h1t, "", "Ecal Time (ns)", "", "logY", 0.95, 2  ) ;
     //h_draw->DrawAppend( h2t, "KsTest", 0.7, 4 ) ;

     printf(" KS Val : %f ", ksVal ) ;
 
     delete h1t ; 
     delete h2t ; 
     return ksVal ;

}

TH1D* StatsLib::RebinHistogram( TH1D* h1, string newHistoName, double minBC ) {

     // accumuate bin information
     vector<double> xV ;
     vector<double> yV ;
     double aveBC = 0 ;
     double sumBC = 0 ;
     double rbin_ = 0 ;
     for ( int i= 1 ; i<= h1->GetNbinsX() ; i++ ) {

         double bc_   = h1->GetBinContent(i) ;
         double x_    = h1->GetBinCenter(i) ;
         double bc1_  = ( i == h1->GetNbinsX() ) ? minBC+1 : h1->GetBinContent(i+1) ;
         //cout<<" x: "<< x_ <<" y:"<< bc_   ;
         if ( aveBC < minBC ) {
            sumBC += bc_ ;
            rbin_ += 1. ;
            aveBC = sumBC / rbin_ ;
            if ( (aveBC >  minBC && bc1_ > minBC) || i == h1->GetNbinsX() ) {
               xV.push_back( x_ ) ;
               yV.push_back( aveBC );
               //cout<<" aveBC: "<< aveBC<<"  from sum:"<< sumBC<<"/"<<rbin_   ;
               sumBC = 0 ;
               rbin_ = 0 ;
            }
            aveBC = 0 ;
         }
         //cout<<" "<<endl ; 
     }
     //cout<<" ====================== "<<endl ; 
     // refill the histogram 
     TH1D* h2 = (TH1D*) h1->Clone() ;
     h2->SetName( newHistoName.c_str() ) ;
     int ii = 0 ;
     for ( int i= 1 ; i<= h1->GetNbinsX() ; i++ ) {
         h2->SetBinContent(i, 0 ) ;
         for ( size_t j= ii ; j< xV.size(); j++ ) {
              double x_    = h2->GetBinCenter(i) ;
              if ( x_ <= xV[j] ) {
                 h2->SetBinContent(i, yV[j] ) ;
                 ii = j ;
                 //cout<<" x: "<< x_ <<" xB:"<< xV[j]<<" y: "<< yV[j] <<endl ;
                 break ;
              }
         }
     }
     return h2 ;
}

