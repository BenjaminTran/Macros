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


Double_t bgfunc2( Double_t *x, Double_t *par )
{
    Double_t xx = x[0];
    //Double_t q = xx - ( piMass + lambdaMass );
    Double_t bg = par[0]*TMath::Power(xx - ( piMass + lambdaMass ),par[1]);
    return bg;
}


//Define double Gaussian fit function
Double_t sgfunc( Double_t *x, Double_t *par )
{
    Double_t xx1 = 0;
    Double_t xx2 = 0;
    if( par[2]!=0 ) xx1 = ( x[0] - par[1] )/par[2];
    if( par[4]!=0 ) xx2 = ( x[0] - par[1] )/par[4];
    Double_t sgaus1 = par[0]*TMath::Exp( -0.5*xx1*xx1 );
    Double_t sgaus2 = par[3]*TMath::Exp( -0.5*xx2*xx2 );
    return sgaus1 + sgaus2;
}

//Define Total fit
Double_t total( Double_t *x, Double_t *par )
{
    return bgfunc2( x,par ) + sgfunc( x,&par[2] );
}

void InvMassCascadeFit(  )
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
    bgfit->SetParameter( 0,3.17496e3 );
    bgfit->SetParameter( 1,5.97019e-1 );
    bgfit->SetParNames( "C","D" );

    TF1 *sigfit = new TF1( "sigfit", sgfunc,1.26,1.4,5 );
    sigfit->SetNpx( 50 );
    sigfit->SetParameter( 0,6.6688e2 );
    sigfit->SetParameter( 1,1.32219 );
    sigfit->SetParameter( 2,6.15655e8 );
    sigfit->SetParameter( 3,3.67793e3 );
    sigfit->SetParameter( 4,3.69916e-3 );

    TF1 *totalfit = new TF1( "totalfit", total , 1.26,1.4,7 );
    totalfit->SetNpx( 250 );
    totalfit->SetParameter( 0,8.61548e3 );
    totalfit->SetParameter( 1,1.71519e-2 );
    totalfit->SetParameter( 2,-7.73080e3 );
    totalfit->SetParameter( 3,1.32217 );
    totalfit->SetParameter( 4,6.15655e8 );
    totalfit->SetParameter( 5,3.64537e3 );
    totalfit->SetParameter( 6,4.16903e-3 );

    /* 2
    totalfit->SetParameter( 0,8.61548e3 );
    totalfit->SetParameter( 1,1.71519e-2 );
    totalfit->SetParameter( 2,-7.73080e3 );
    totalfit->SetParameter( 3,1.32217 );
    totalfit->SetParameter( 4,6.15655e8 );
    totalfit->SetParameter( 5,3.64537e3 );
    totalfit->SetParameter( 6,4.16903e-3 );
    */
   /* 1
    totalfit->SetParameter( 0,3.17496e3 );
    totalfit->SetParameter( 1,5.97019e-1 );
    totalfit->SetParameter( 2,6.6688e2 );
    totalfit->SetParameter( 3,1.32219 );
    totalfit->SetParameter( 4,6.15655e8 );
    totalfit->SetParameter( 5,3.67793e3 );
    totalfit->SetParameter( 6,3.69916e-3 );
    */

    TF1 *bgfit2 = new TF1( "bgfit2", bgfunc2, 1.34,1.4,2 );
    bgfit2->SetNpx( 60 );
    bgfit2->SetParameter( 0,1 );
    bgfit2->SetParameter( 1,1 );
    bgfit2->SetParNames( "E","F" );
    
    //TMath::MinimizerOptions::SetDefaultMaxFunctionCalls( 3000);
    //TVirtualFitter::SetMaxIterations( 20000 );

    Double_t par[7];

    invMass->Fit( "bgfit2","R","",1.34,1.4 );
    invMass->Fit( "bgfit","R+","",1.26,1.31 );
    invMass->Fit( "sigfit", "R+","",1.31,1.34 ); //1.335

    bgfit->GetParameters( &par[0] );
    


    //invMass->Fit( "totalfit","L","",1.26,1.4);
    


}
