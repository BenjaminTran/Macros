#include <TLatex.h>
#include <TStyle.h>
#include "TF1.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TVirtualFitter.h"

#define piMass 0.13957018
#define lambdaMass 1.115683
//Define background fit functions
Double_t bgfunc1( Double_t *x, Double_t *par )
{
    Double_t xx = x[0];
    //Double_t q = xx - ( piMass + lambdaMass );
    Double_t bg = par[0]*TMath::Power(xx - ( piMass + lambdaMass ),par[1]);
    return bg;
}


/*
Double_t bgfunc2( Double_t *x, Double_t *par )
{
    Double_t xx = x[0];
    //Double_t q = xx - ( piMass + lambdaMass );
    Double_t bg = par[0]*TMath::Power(xx - ( piMass + lambdaMass ),par[1]);
    return bg;
}
*/


//Define double Gaussian fit function
Double_t sgfunc( Double_t *x, Double_t *par )
{
    Double_t xx1 = 0;
    Double_t xx2 = 0;
    if( par[2]!=0 ) xx1 = ( x[0] - par[1] )/par[2];
    if( par[4]!=0 ) xx2 = ( x[0] - par[1] )/par[4];
    Double_t sgaus1 = par[0]*TMath::Exp( -0.5*xx1*xx1 );
    Double_t sgaus2 = par[3]*TMath::Exp( -0.5*xx2*xx2 );
    Double_t bgsig  = par[5]*TMath::Power( x[0] - ( piMass + lambdaMass ),par[6] );
    return sgaus1 + sgaus2 + bgsig;
}

//Define Total fit
/*
Double_t total( Double_t *x, Double_t *par )
{
    return bgfunc1( x,par ) + bgfunc2( x,&par[2] ) + sgfunc( x,&par[4] );
}
*/

Double_t total( Double_t *x, Double_t *par )
{
    return bgfunc1( x,par ) + sgfunc( x,&par[2] );
}



void InvMassFit(  )
{
    
    gStyle->SetOptFit( 1111 );
    TFile *f = new TFile(
            "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/XiAnalysisCorrelation.root " );
    /*
    //Pull 2D Histograms
    TH2F *hbackground = (TH2F*)   f->Get( "xiCorrelation/Background" );
    TH2F *hsignal     = ( TH2F* ) f->Get( "xiCorrelation/Signal" );

    //Project azimuth
    TH1F* hbPhi = hbackground->ProjectionY(  );
    TH1F* hsPhi = hsignal->ProjectionY(  );
    */

    //Pull InvMass
    TH1F* invMass = ( TH1F* ) f->Get( "xiCorrelation/InvMassXi" );


    TF1 *bgfit = new TF1( "bgfit", bgfunc1, 1.26,1.305,2 );
    bgfit->SetNpx( 60 );
    bgfit->SetParameter( 0,3.17496e1 );
    bgfit->SetParameter( 1,5.97019e-1 );
    bgfit->SetParNames( "C","D" );

    TF1 *sigfit = new TF1( "sigfit", sgfunc,1.26,1.4,7 );
    //sigfit->SetNpx( 50 );
    sigfit->SetNpx( 250 );
    sigfit->SetParNames( "gaus1Norm","gausMean","gaus1std","gaus2Norm","gaus2std","sbkgScale","sbkgPow" );
    sigfit->SetParameter( 0,6.6688e2 );
    sigfit->SetParameter( 1,1.32219 );
    sigfit->SetParameter( 2,3.15655e-3 );
    sigfit->SetParameter( 3,3.67793e3 );
    sigfit->SetParameter( 4,3.69916e-3 );
    sigfit->SetParameter( 5,3.21847e3 );
    sigfit->SetParameter( 6,5.99562-1 );
    

    //TF1 *totalfit = new TF1( "totalfit", total , 1.26,1.4,11 );
    //totalfit->SetNpx( 250 );

    TF1 *testfunc = new TF1( "testfunc", total, 1.26, 1.4, 9 );
    testfunc->SetNpx( 250 );
    testfunc->SetParNames( "bkgScale","bkgPow","gaus1Norm","gausMean","gaus1std","gaus2Norm","gaus2std","sbkgScale","sbkgPow" );

    /*
    TF1 *bgfit2 = new TF1( "bgfit2", bgfunc2, 1.34,1.4,2 );
    bgfit2->SetNpx( 60 );
    bgfit2->SetParameter( 0,6.05272e2 );
    bgfit2->SetParameter( 1,4.88774e-2 );
    bgfit2->SetParNames( "E","F" );
    */

    TVirtualFitter::SetMaxIterations( 300000 );

    //Double_t par[11];
    //for test function
    Double_t partest[9];

    //invMass->Fit( "bgfit2","L ","",1.34,1.4 );
    //invMass->Fit( "bgfit","L","",1.26,1.31 );
    invMass->Fit( "sigfit", "L","",1.26,1.4 ); //1.335

    /*
    bgfit->GetParameters( &par[0] );
    bgfit2->GetParameters( &par[2] );
    sigfit->GetParameters( &par[4] );
    */
    
    bgfit->GetParameters( &partest[0] );
    sigfit->GetParameters( &partest[2] );

    //totalfit->SetParameters( par );
    testfunc->SetParameters( partest );
    

    //invMass->Fit( "totalfit","L ","",1.26,1.4);
    //invMass->Fit( "testfunc", "L","", 1.26, 1.4 );
    


}
