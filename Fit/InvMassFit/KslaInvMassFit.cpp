#include <TLatex.h>
#include <TStyle.h>
#include "TF1.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TLegend.h"
#include "TMathText.h"
#include "TImage.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TPad.h"
#include "TFile.h"
#include <TString.h>
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGaxis.h"

#include <vector>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdio>

#define piMass 0.13957018
#define lambdaMass 1.115683
#define PI 3.1416


//Define background fit functions

Double_t V0bgfunc( Double_t *x, Double_t *par )
{
    Double_t xx = x[0];
    Double_t bg = par[0] + par[1]*xx + par[2]*xx*xx + par[3]*TMath::Power( xx, 3 ) + par[4]*TMath::Power( xx,4 ); 
    return bg;
}

//Define double Gaussian fit function
Double_t V0sgfunc( Double_t *x, Double_t *par )
{
    Double_t xx1 = 0;
    Double_t xx2 = 0;
    if( par[2]!=0 ) xx1 = ( x[0] - par[1] )/par[2];
    if( par[4]!=0 ) xx2 = ( x[0] - par[1] )/par[4];
    Double_t sgaus1 = par[0]*TMath::Exp( -0.5*xx1*xx1 );
    Double_t sgaus2 = par[3]*TMath::Exp( -0.5*xx2*xx2 );
    return sgaus1 + sgaus2;
}

Double_t V0total( Double_t *x, Double_t *par )
{
    return V0bgfunc( x,par ) + V0sgfunc( x,&par[5] );
}

void KslaInvMassFit(  )
{
    using namespace std;

    int numPtBins_ks = 12;
    int numPtBins_la = 11;
    double pks[] = {0.6,0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,9.0,12.0}; //if number of bins changes make sure you change numPtBins
    double pla[] = {0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,9.0,12.0}; //if number of bins changes make sure you change numPtBins

    TFile* f = new TFile( "/Volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/Ksla/kslaMassPtJL1.root" );

    std::ostringstream os;
    TH1::SetDefaultSumw2(  );

    std::vector<double> PtBin_ks( pks, pks+numPtBins_ks );
    std::vector<double> PtBin_la( pla, pla+numPtBins_la );

    bool publish = false; // Adds Latex labels
    bool Correlation = false;

    int PtBinSize_ks = PtBin_ks.size(  ) - 1;
    int PtBinSize_la = PtBin_la.size(  ) - 1;

    std::vector<TH1D*> InvMassPtBinned_ks;
    std::vector<TF1*> BgFit_ks;
    std::vector<TF1*> SigFit_ks;
    std::vector<TF1*> FitTot_ks;
    std::vector<TF1*> FitFcn_pt_ks;
    std::vector<double> Fsig_ks;
    std::vector<double> sig_ks;
    std::vector<double> MassMean_ks;
    std::vector<double> MassStd_ks;

    std::vector<TH1D*> InvMassPtBinned_la;
    std::vector<TF1*> BgFit_la;
    std::vector<TF1*> SigFit_la;
    std::vector<TF1*> FitTot_la;
    std::vector<TF1*> FitFcn_pt_la;
    std::vector<double> Fsig_la;
    std::vector<double> sig_la;
    std::vector<double> MassMean_la;
    std::vector<double> MassStd_la;

    TH2D* MassPt_ks;
    TH2D* MassPt_la;

    if( Correlation )
    {
        MassPt_ks = ( TH2D* )f->Get( "kslaCorrelation/ksMassPt" );
        MassPt_la = ( TH2D* )f->Get( "kslaCorrelation/laMassPt" );
    }
    else
    {
        MassPt_ks = ( TH2D* )f->Get( "KslaMassPt/KsMassPt" );
        MassPt_la = ( TH2D* )f->Get( "KslaMassPt/LaMassPt" );
    }


    TVirtualFitter::SetMaxIterations( 300000 );

    //Fit Total Inv Mass

    TH1D* InvMass_ks = ( TH1D* )MassPt_ks->ProjectionX(  Form( "InvMass_ks_%d",0 ), 0, 100 );
    InvMass_ks->SetMarkerStyle( 21 );
    InvMass_ks->SetMarkerSize( 0.5 );
    TH1D* InvMass_la = ( TH1D* )MassPt_la->ProjectionX(  Form( "InvMass_la_%d",0 ), 0, 100 );
    InvMass_la->SetMarkerStyle( 21 );
    InvMass_la->SetMarkerSize( 0.5 );

    TF1* bgfit_ks = new TF1( Form( "bgfit_ks_%d",0 ), V0bgfunc, 0.432,0.54,5 );
    bgfit_ks->SetNpx( 60 );
    bgfit_ks->SetParameter( 0,-3e2 );
    bgfit_ks->SetParameter( 1,0 );
    bgfit_ks->SetParameter( 2,0 );
    bgfit_ks->SetParameter( 3,0 );
    bgfit_ks->SetParameter( 4,0 );
    bgfit_ks->SetParNames( "const_ks","pow1_ks","pow2_ks","pow3_ks","pow4_ks" );
    bgfit_ks->SetLineColor( kBlue );
    TF1* bgfit_la = new TF1( Form( "bgfit_la_%d",0 ), V0bgfunc, 1.08,1.155,5 );
    bgfit_la->SetNpx( 60 );
    bgfit_la->SetParameter( 0,1e6 );
    bgfit_la->SetParameter( 1,100 );
    bgfit_la->SetParameter( 2,100 );
    bgfit_la->SetParameter( 3,100 );
    bgfit_la->SetParameter( 4,100 );
    bgfit_la->SetParNames( "const_la","pow1_la","pow2_la","pow3_la","pow4_la" );

    TF1* sigfit_ks = new TF1( Form( "sigfit_ks_%d",0 ), V0sgfunc,0.43,0.565,5 );
    sigfit_ks->SetNpx( 250 );
    sigfit_ks->SetParNames( "gaus1Norm_ks","gausMean_ks","gaus1std_ks","gaus2Norm_ks","gaus2std_ks");
    sigfit_ks->SetParameter( 0,6.6688e3 );
    sigfit_ks->SetParameter( 1,1.32219 );
    sigfit_ks->SetParameter( 2,3.15655e-3 );
    sigfit_ks->SetParameter( 3,3.67793e3 );
    sigfit_ks->SetParameter( 4,3.69916e-3 );
    TF1* sigfit_la = new TF1( Form( "sigfit_la_%d",0 ), V0sgfunc,0.43,0.565,5 );
    sigfit_la->SetNpx( 250 );
    sigfit_la->SetParNames( "gaus1Norm_la","gausMean_la","gaus1std_la","gaus2Norm_la","gaus2std_la");
    sigfit_la->SetParameter( 0,6.6688e3 );
    sigfit_la->SetParameter( 1,1.32219 );
    sigfit_la->SetParameter( 2,3.15655e-3 );
    sigfit_la->SetParameter( 3,3.67793e3 );
    sigfit_la->SetParameter( 4,3.69916e-3 );

    TF1* fitTot_ks = new TF1( Form( "fitTot_ks_%d",0 ), V0total, 0.43,0.565, 10 );
    fitTot_ks->SetNpx( 250 );
    fitTot_ks->SetParNames( "const_ks","pow1_ks","pow2_ks","pow3_ks","pow4_ks","gaus1Norm_ks","gausMean_ks","gaus1std_ks","gaus2Norm_ks","gaus2std_ks");
    TF1* fitTot_la = new TF1( Form( "fitTot_la_%d",0 ), V0total, 1.08,1.155, 10 );
    fitTot_la->SetNpx( 250 );
    fitTot_la->SetParNames( "const_la","pow1_la","pow2_la","pow3_la","pow4_la","gaus1Norm_la","gausMean_la","gaus1std_la","gaus2Norm_la","gaus2std_la");

    //for test function
    Double_t partestTot_ks[10];
    Double_t partestTot_la[10];

    InvMass_ks->Fit( Form( "bgfit_ks_%d",0 ),"L","",0.43,0.565 );
    InvMass_ks->Fit( Form( "sigfit_ks_%d",0 ), "L R","",0.43,0.565 );
    InvMass_la->Fit( Form( "bgfit_la_%d",0 ),"L","",1.08,1.155 );
    InvMass_la->Fit( Form( "sigfit_la_%d",0 ), "L R","",1.08,1.155 );

    bgfit_ks->GetParameters( &partestTot_ks[0] );
    sigfit_ks->GetParameters( &partestTot_ks[5] );
    bgfit_la->GetParameters( &partestTot_la[0] );
    sigfit_la->GetParameters( &partestTot_la[5] );

    fitTot_ks->SetParameters( partestTot_ks );
    fitTot_la->SetParameters( partestTot_la );


    InvMass_ks->SetTitle( "" );
    InvMass_ks->GetYaxis(  )->SetTitleOffset( 1.5 );
    InvMass_ks->GetYaxis(  )->SetTitleSize( 0.035 );
    InvMass_ks->GetYaxis(  )->CenterTitle( true );
    InvMass_ks->GetYaxis(  )->SetTitle( "Candidates/0.5 MeV" );
    InvMass_ks->GetYaxis(  )->SetRangeUser( 0,3e6 );
    InvMass_ks->SetTitleOffset( 1.5, "X" );
    InvMass_ks->GetXaxis(  )->SetTitleSize( 0.035 );
    InvMass_ks->GetXaxis(  )->CenterTitle( true );
    InvMass_ks->GetXaxis(  )->SetTitle( "#pi^{#plus}#pi^{#minus} invariant mass (GeV)" );
    InvMass_ks->SetStats( !publish );
    //Drawing the fit function again for the sake of getting it onto the legend
    InvMass_ks->Fit( Form( "fitTot_ks_%d",0 ), "L0","", 0.43,0.565 );
    Double_t parfordrawTot_ks[10];
    fitTot_ks->GetParameters( &parfordrawTot_ks[0] );
    TF1* fitFcn_ks = new TF1( Form( "fitFcn_ks_%d",0 ), V0total, 0.43,0.565,10 );
    fitFcn_ks->SetNpx( 250 );
    fitFcn_ks->SetParameters( parfordrawTot_ks );
    fitFcn_ks->SetLineColor( kRed );
    /* Only if signal fit needs a background piece as well
       TF1* backFcn_ks = new TF1( Form( "backFcn_ks_%d",0 ), xi_bkgfuncdisplay, 1.26,1.4,4 );
       backFcn_->SetParameter( 0,fitTot_->GetParameter( 0 ) );
       backFcn_->SetParameter( 1,fitTot_->GetParameter( 1 ) );
       backFcn_->SetParameter( 2,fitTot_->GetParameter( 7 ) );
       backFcn_->SetParameter( 3,fitTot_->GetParameter( 8 ) );
       backFcn_->SetLineColor( kBlue );
       */

    InvMass_la->SetTitle( "" );
    InvMass_la->GetYaxis(  )->SetTitleOffset( 1.5 );
    InvMass_la->GetYaxis(  )->SetTitleSize( 0.035 );
    InvMass_la->GetYaxis(  )->CenterTitle( true );
    InvMass_la->GetYaxis(  )->SetTitle( "Candidates/0.5 MeV" );
    InvMass_la->GetYaxis(  )->SetRangeUser( 0,2e6 );
    InvMass_la->SetTitleOffset( 1.5, "X" );
    InvMass_la->GetXaxis(  )->SetTitleSize( 0.035 );
    InvMass_la->GetXaxis(  )->CenterTitle( true );
    InvMass_la->GetXaxis(  )->SetTitle( "p#pi^{#minus} + charge conjugate invariant mass (GeV)" );
    InvMass_la->SetStats( !publish );
    //Drawing the fit function again for the sake of getting it onto the legend
    InvMass_la->Fit( Form( "fitTot_la_%d",0 ), "L0","", 1.08,1.155 );
    Double_t parfordrawTot_la[10];
    fitTot_la->GetParameters( &parfordrawTot_la[0] );
    TF1* fitFcn_la = new TF1( Form( "fitFcn_la_%d",0 ), V0total, 1.08,1.155,10 );
    fitFcn_la->SetNpx( 250 );
    fitFcn_la->SetParameters( parfordrawTot_la );
    fitFcn_la->SetLineColor( kRed );

    double rmsSigmaTot_ks = TMath::Sqrt( 0.5*fitTot_ks->GetParameter( "gaus1std_ks" )*fitTot_ks->GetParameter( "gaus1std_ks" ) + 0.5*fitTot_ks->GetParameter( "gaus2std_ks" )*fitTot_ks->GetParameter( "gaus2std_ks" ) );
    double rmsSigmaTot_la = TMath::Sqrt( 0.5*fitTot_la->GetParameter( "gaus1std_la" )*fitTot_la->GetParameter( "gaus1std_la" ) + 0.5*fitTot_la->GetParameter( "gaus2std_la" )*fitTot_la->GetParameter( "gaus2std_la" ) );
    double meanTot_ks = fitTot_ks->GetParameter( "gausMean_ks" );
    double meanTot_la = fitTot_la->GetParameter( "gausMean_la" );

    /*
    for( int i=0; i<PtBinSize_ks; i++ ){
        TH1D* InvMassPtBinned_k = ( TH1D* )MassPt_ks->ProjectionX( Form( "InvMass_pT_ks_%d",i ), PtBin_ks[i]*10+1, PtBin_ks[i+1]*10 );
        os << PtBin_ks[i] << "_Pt_" << PtBin_ks[i+1];
        InvMassPtBinned_k->SetTitle( os.str(  ).c_str(  ) );
        os.str( std::string(  ) );


        //Perform Fits Kshort
        TF1* bgFit_ks = new TF1(  Form( "bgfit_ks_%d",i ), V0bgfunc, 0.43,0.565,5  );
        bgFit_ks->SetNpx( 60 );
        bgFit_ks->SetParameter( 0,1 );
        bgFit_ks->SetParameter( 1,1 );
        bgFit_ks->SetParameter( 2,1 );
        bgFit_ks->SetParameter( 3,1 );
        bgFit_ks->SetParameter( 4,1 );
        bgFit_ks->SetParNames( "const_ks","pow1_ks","pow2_ks","pow3_ks","pow4_ks" );

        TF1* sigFit_ks = new TF1( Form( "sigfit_ks_%d",i ), V0sgfunc, 0.43,0.565,5 ) ;
        sigFit_ks->SetParameter( 0,1 );
        sigFit_ks->SetParameter( 1,1 );
        sigFit_ks->SetParameter( 2,1 );
        sigFit_ks->SetParameter( 3,1 );
        sigFit_ks->SetParameter( 4,1 );
        sigFit_ks->SetParNames( "gaus1Norm_ks","gausMean_ks","gaus1std_ks","gaus2Norm_ks","gaus2std_ks");

        TF1* FitTot_pt_ks = new TF1( Form( "fitTot_ks_%d",i ), V0total, 0.43,0.565, 10 ) ;
        FitTot_pt_ks->SetNpx( 250 );
        FitTot_pt_ks->SetParNames( "const_ks","pow1_ks","pow2_ks","pow3_ks","pow4_ks","gaus1Norm_ks","gausMean_ks","gaus1std_ks","gaus2Norm_ks","gaus2std_ks");

        Double_t partest_ks[10];

        InvMassPtBinned_k->Fit( Form( "bgfit_ks_%d",i ), "L","",0.43,0.565 );
        InvMassPtBinned_k->Fit( Form( "sigfit_ks_%d",i ), "L R", "",0.43,0.565 );

        bgFit_ks ->GetParameters( &partest_ks[0] );
        sigFit_ks->GetParameters( &partest_ks[5] );

        FitTot_pt_ks->SetParameters( partest_ks );

        InvMassPtBinned_k->Fit( Form( "FitTot_pt_ks_%d",i ), "L","",0.43,0.565 );

        //Draw signal and background fit also need bkg fit for fsig calc

        Double_t parfordraw_ks[10];
        FitTot_pt_ks->GetParameters( &parfordraw_ks[0] );
        TF1* fitFcn_pt_ks = new TF1( Form( "fitFcn_pt_ks_%d",i ), V0total, 0.43,0.565,10 );
        fitFcn_pt_ks->SetNpx( 250 );
        fitFcn_pt_ks->SetParameters( parfordraw_ks );
        fitFcn_pt_ks->SetLineColor( kRed );
        fitFcn_pt_ks->Draw( "same" );

        //Calculate values

        //Sigma
        double rmsSigma_ks = TMath::Sqrt( 0.5*FitTot_pt_ks->GetParameter( "gaus1std_ks" )*FitTot_pt_ks->GetParameter( "gaus1std_ks" ) + 0.5*FitTot_pt_ks->GetParameter( "gaus2std_ks" )*FitTot_pt_ks->GetParameter( "gaus2std_ks" ) );      

        //Signal fraction yield
        double totInt_ks = -999.0;
        double bkgInt_ks = -999.0;

        totInt_ks = FitTot_pt_ks->Integral( FitTot_pt_ks->GetParameter( "gausMean_ks" ) - 2*rmsSigma_ks, FitTot_pt_ks->GetParameter( "gausMean_ks" ) + 2*rmsSigma_ks );

        bkgInt_ks = bgFit_ks->Integral( FitTot_pt_ks->GetParameter( "gausMean_ks" ) - 2*rmsSigma_ks, FitTot_pt_ks->GetParameter( "gausMean_ks" ) + 2*rmsSigma_ks );

        Fsig_ks.push_back( ( totInt_ks - bkgInt_ks )/totInt_ks );
        sig_ks.push_back( totInt_ks - bkgInt_ks );
        MassMean_ks.push_back( FitTot_pt_ks->GetParameter( "gausMean_ks" ) );
        MassStd_ks.push_back( rmsSigma_ks );

        BgFit_ks.push_back( bgFit_ks );
        SigFit_ks.push_back( sigFit_ks );
        FitTot_ks.push_back( FitTot_pt_ks );
        FitFcn_pt_ks.push_back( fitFcn_pt_ks );
        InvMassPtBinned_ks.push_back( InvMassPtBinned_k );
    }

    
    for( int i=0; i<PtBinSize_la; i++ ){
        TH1D* InvMassPtBinned_l = ( TH1D* )MassPt_la->ProjectionX( Form( "InvMass_pT_la_%d",i ), PtBin_la[i]*10+1, PtBin_la[i+1]*10 );
        os << PtBin_la[i] << "_Pt_" << PtBin_la[i+1];
        InvMassPtBinned_l->SetTitle( os.str(  ).c_str(  ) );
        os.str( std::string(  ) );

        //Perform Fits Lambda
        TF1* bgFit_la = new TF1(  Form( "bgfit_la_%d",i ), V0bgfunc, 0.43,0.565,5  );
        bgFit_la->SetNpx( 60 );
        bgFit_la->SetParameter( 0,1 );
        bgFit_la->SetParameter( 1,1 );
        bgFit_la->SetParameter( 2,1 );
        bgFit_la->SetParameter( 3,1 );
        bgFit_la->SetParameter( 4,1 );
        bgFit_la->SetParNames( "const_la","pow1_la","pow2_la","pow3_la","pow4_la" );

        TF1* sigFit_la = new TF1( Form( "sigfit_la_%d",i ), V0sgfunc, 0.43,0.565,5 ) ;
        sigFit_la->SetParameter( 0,1 );
        sigFit_la->SetParameter( 1,1 );
        sigFit_la->SetParameter( 2,1 );
        sigFit_la->SetParameter( 3,1 );
        sigFit_la->SetParameter( 4,1 );
        sigFit_la->SetParNames( "gaus1Norm_la","gausMean_la","gaus1std_la","gaus2Norm_la","gaus2std_la");

        TF1* FitTot_pt_la = new TF1( Form( "fitTot_la_%d",i ), V0total, 0.43,0.565, 10 ) ;
        FitTot_pt_la->SetNpx( 250 );
        FitTot_pt_la->SetParNames( "const_la","pow1_la","pow2_la","pow3_la","pow4_la","gaus1Norm_la","gausMean_la","gaus1std_la","gaus2Norm_la","gaus2std_la");

        Double_t partest_la[10];

        InvMassPtBinned_l->Fit( Form( "bgfit_la_%d",i ), "L","",0.43,0.565 );
        InvMassPtBinned_l->Fit( Form( "sigfit_la_%d",i ), "L R", "",0.43,0.565 );

        bgFit_la ->GetParameters( &partest_la[0] );
        sigFit_la->GetParameters( &partest_la[5] );

        FitTot_pt_la->SetParameters( partest_la );

        InvMassPtBinned_l->Fit( Form( "FitTot_pt_la_%d",i ), "L","",0.43,0.565 );

        //Draw signal and background fit also need bkg fit for fsig calc

        Double_t parfordraw_la[10];
        FitTot_pt_la->GetParameters( &parfordraw_la[0] );
        TF1* fitFcn_pt_la = new TF1( Form( "fitFcn_pt_la_%d",i ), V0total, 0.43,0.565,10 );
        fitFcn_pt_la->SetNpx( 250 );
        fitFcn_pt_la->SetParameters( parfordraw_la );
        fitFcn_pt_la->SetLineColor( kRed );
        fitFcn_pt_la->Draw( "same" );
        //Calculate values

        //Sigma
        double rmsSigma_la = TMath::Sqrt( 0.5*FitTot_pt_la->GetParameter( "gaus1std_la" )*FitTot_pt_la->GetParameter( "gaus1std_la" ) + 0.5*FitTot_pt_la->GetParameter( "gaus2std_la" )*FitTot_pt_la->GetParameter( "gaus2std_la" ) );            

        //Signal fraction yield
        double totInt_la = -999.0;
        double bkgInt_la = -999.0;

        totInt_la = FitTot_pt_la->Integral( FitTot_pt_la->GetParameter( "gausMean_la" ) - 2*rmsSigma_la, FitTot_pt_la->GetParameter( "gausMean_la" ) + 2*rmsSigma_la );

        bkgInt_la = bgFit_la->Integral( FitTot_pt_la->GetParameter( "gausMean_la" ) - 2*rmsSigma_la, FitTot_pt_la->GetParameter( "gausMean_la" ) + 2*rmsSigma_la );

        Fsig_la.push_back( ( totInt_la - bkgInt_la )/totInt_la );
        sig_la.push_back( totInt_la - bkgInt_la );
        MassMean_la.push_back( FitTot_pt_la->GetParameter( "gausMean_la" ) );
        MassStd_la.push_back( rmsSigma_la );

        BgFit_la.push_back( bgFit_la );
        SigFit_la.push_back( sigFit_la );
        FitTot_la.push_back( FitTot_pt_la );
        FitFcn_pt_la.push_back( fitFcn_pt_la );
        InvMassPtBinned_la.push_back( InvMassPtBinned_l );


    }
        
    

    int PtBinCounter=0;

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "K SHORT K SHORT K SHORT" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    // Sig yield values
    cout << "==========================================================" << endl;
    cout << "Signal Yield" << endl;
    cout << "==========================================================" << endl;

    for( std::vector<double>::iterator it = sig_ks.begin(  ); it != sig_ks.end(  ); ++it ){ 
        cout << PtBin_ks[PtBinCounter] << " < Pt =< " << PtBin_ks[PtBinCounter+1] << ": " << *it << endl;
        PtBinCounter++;
    }
    PtBinCounter = 0;

    // Fsig values
    cout << "==========================================================" << endl;
    cout << "Signal Fraction Yield" << endl;
    cout << "==========================================================" << endl;
    for( std::vector<double>::iterator it = Fsig_ks.begin(  ); it != Fsig_ks.end(  ); ++it ){ 
        cout << PtBin_ks[PtBinCounter] << " < Pt =< " << PtBin_ks[PtBinCounter+1] << ": " << *it << endl;
        PtBinCounter++;
    }
    PtBinCounter = 0;

    // MassMean
    cout << "==========================================================" << endl;
    cout << "Mean Mass" << endl;
    cout << "==========================================================" << endl;

    for( std::vector<double>::iterator it = MassMean_ks.begin(  ); it != MassMean_ks.end(  ); ++it ){ 
        cout << PtBin_ks[PtBinCounter] << " < Pt =< " << PtBin_ks[PtBinCounter+1] << ": " << *it << endl;
        PtBinCounter++;
    }
    PtBinCounter = 0;

    // RMS Sigma
    cout << "==========================================================" << endl;
    cout << "RMS Sigma" << endl;
    cout << "==========================================================" << endl;

    for( std::vector<double>::iterator it = MassStd_ks.begin(  ); it != MassStd_ks.end(  ); ++it ){ 
        cout << PtBin_ks[PtBinCounter] << " < Pt =< " << PtBin_ks[PtBinCounter+1] << ": " << *it << endl;
        PtBinCounter++;
    }
    PtBinCounter = 0;


    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
    cout << "LAMBDA LAMBDA LAMBDA" << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

    // Sig yield values
    cout << "==========================================================" << endl;
    cout << "Signal Yield" << endl;
    cout << "==========================================================" << endl;

    for( std::vector<double>::iterator it = sig_la.begin(  ); it != sig_la.end(  ); ++it ){ 
        cout << PtBin_la[PtBinCounter] << " < Pt =< " << PtBin_la[PtBinCounter+1] << ": " << *it << endl;
        PtBinCounter++;
    }
    PtBinCounter = 0;

    // Fsig values
    cout << "==========================================================" << endl;
    cout << "Signal Fraction Yield" << endl;
    cout << "==========================================================" << endl;
    for( std::vector<double>::iterator it = Fsig_la.begin(  ); it != Fsig_la.end(  ); ++it ){ 
        cout << PtBin_la[PtBinCounter] << " < Pt =< " << PtBin_la[PtBinCounter+1] << ": " << *it << endl;
        PtBinCounter++;
    }
    PtBinCounter = 0;

    // MassMean
    cout << "==========================================================" << endl;
    cout << "Mean Mass" << endl;
    cout << "==========================================================" << endl;

    for( std::vector<double>::iterator it = MassMean_la.begin(  ); it != MassMean_la.end(  ); ++it ){ 
        cout << PtBin_la[PtBinCounter] << " < Pt =< " << PtBin_la[PtBinCounter+1] << ": " << *it << endl;
        PtBinCounter++;
    }
    PtBinCounter = 0;

    // RMS Sigma
    cout << "==========================================================" << endl;
    cout << "RMS Sigma" << endl;
    cout << "==========================================================" << endl;

    for( std::vector<double>::iterator it = MassStd_la.begin(  ); it != MassStd_la.end(  ); ++it ){ 
        cout << PtBin_la[PtBinCounter] << " < Pt =< " << PtBin_la[PtBinCounter+1] << ": " << *it << endl;
        PtBinCounter++;
    }
    PtBinCounter = 0;
    

    // Draw Histgorams
    TCanvas* cHist_ks = new TCanvas( "Canvas_ks","KS",1600,800 );
    cHist_ks->Divide( 4,3 );
    for( int i=0; i<numPtBins_ks; i++ ){
        cHist_ks->cd( i+1 );
        InvMassPtBinned_ks[i]->Draw( "E1" );
    }

    TCanvas* cHist_la = new TCanvas( "Canvas_la","KS",1600,800 );
    cHist_la->Divide( 4,3 );
    for( int i=0; i<numPtBins_la; i++ ){
        cHist_la->cd( i+1 );
        InvMassPtBinned_la[i]->Draw( "E1" );
    }
    */

    TH1D* PtSpecra_ks = ( TH1D* )MassPt_ks->ProjectionY( "Pt Spectrum",0,150 );
    TCanvas* Pt_ks = new TCanvas( "Canvas_ptks","", 800, 800 );
    Pt_ks->cd(  );
    PtSpecra_ks->Draw( "E1" );

    TH1D* PtSpecra_la = ( TH1D* )MassPt_la->ProjectionY( "Pt Spectrum",0,150 );
    TCanvas* Pt_la = new TCanvas( "Canvas_ptla","", 800, 800 );
    Pt_la->cd(  );
    PtSpecra_la->Draw( "E1" );

    TCanvas *cMass_ks = new TCanvas( "cMass_ks", "", 800,800 );
    cMass_ks->SetFillColor      (0);
    cMass_ks->SetBorderMode     (0);
    cMass_ks->SetBorderSize     (10);
    cMass_ks->SetLeftMargin     (0.20);
    cMass_ks->SetRightMargin    (0.06);
    cMass_ks->SetTopMargin      (0.11);
    cMass_ks->SetBottomMargin   (0.15);
    cMass_ks->SetFrameFillStyle (0);
    cMass_ks->SetFrameLineStyle (0);
    cMass_ks->SetFrameBorderMode(0);
    cMass_ks->SetFrameBorderSize(10);
    cMass_ks->SetFrameFillStyle (0);
    cMass_ks->SetFrameLineStyle (0);
    cMass_ks->SetFrameBorderMode(0);
    cMass_ks->SetFrameBorderSize(10);

    TCanvas *cMass_la = new TCanvas( "cMass_la", "", 800,800 );
    cMass_la->SetFillColor      (0);
    cMass_la->SetBorderMode     (0);
    cMass_la->SetBorderSize     (10);
    cMass_la->SetLeftMargin     (0.20);
    cMass_la->SetRightMargin    (0.06);
    cMass_la->SetTopMargin      (0.11);
    cMass_la->SetBottomMargin   (0.15);
    cMass_la->SetFrameFillStyle (0);
    cMass_la->SetFrameLineStyle (0);
    cMass_la->SetFrameBorderMode(0);
    cMass_la->SetFrameBorderSize(10);
    cMass_la->SetFrameFillStyle (0);
    cMass_la->SetFrameLineStyle (0);
    cMass_la->SetFrameBorderMode(0);
    cMass_la->SetFrameBorderSize(10);

    TGaxis::SetMaxDigits( 3 );

    cMass_la->cd(  );
    gPad->SetTickx(  );
    gPad->SetTicky(  );
    InvMass_la->Draw( "E1 9" );
    fitFcn_la->Draw( "same" );
    bgfit_la->Draw( "same" );

    cMass_ks->cd(  );
    gPad->SetTickx(  );
    gPad->SetTicky(  );
    InvMass_ks->Draw( "E1 9" );
    //fitFcn_ks->Draw( "same" );
    bgfit_ks->Draw( "same" );

    if( publish ) 
    {
        cMass_ks->cd(  );
        TLatex* ltx1 = new TLatex(  );
        TLatex* ltx2 = new TLatex(  );
        ltx1->SetTextSize( 0.032 );
        ltx2->SetTextSize( 0.045 );
        ltx1->SetNDC( kTRUE );
        ltx2->SetNDC( kTRUE );

        os << "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = " << " 8.16 TeV";
        ltx1->DrawLatex( 0.23, 0.82, os.str(  ).c_str(  ) );
        os.str( std::string(  ) ); 
        os << "L_{#lower[-0.25]{int}} = " << 62 << " nb^{#font[122]{\55}1}";
        ltx1->DrawLatex( 0.23, 0.74, os.str(  ).c_str(  ) );
        os.str( std::string(  ) );
        os << "#splitline{Mean: " << fixed << std::setprecision( 4 ) <<  meanTot_ks << " GeV}{#sigma : " << rmsSigmaTot_ks << " GeV}";
        ltx1->DrawLatex( 0.65, 0.75, os.str(  ).c_str(  ) );
        os.str( std::string(  ) ); // clears ostringstream
        os << "185 #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220";
        ltx1->DrawLatex( 0.23, 0.67, os.str(  ).c_str(  ) );
        os.str( std::string(  ) );
        ltx2->DrawLatex( 0.23, 0.60, "#K^{0}_{S} " );

        /*
        TLegend* leg = new TLegend( 0.6,0.5,0.85,0.6);
        leg->AddEntry( "fitFcn", "Complete fit function", "1L" );
        leg->AddEntry( "", "Background fit", "1L" );
        leg->SetBorderSize( 0 );
        leg->SetFillColor( 0 );
        leg->Draw(  );
        */
         
        cMass_la->cd(  );

        os << "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = " << " 8.16 TeV";
        ltx1->DrawLatex( 0.23, 0.82, os.str(  ).c_str(  ) );
        os.str( std::string(  ) ); 
        os << "L_{#lower[-0.25]{int}} = " << 62 << " nb^{#font[122]{\55}1}";
        ltx1->DrawLatex( 0.23, 0.74, os.str(  ).c_str(  ) );
        os.str( std::string(  ) );
        os << "#splitline{Mean: " << fixed << std::setprecision( 4 ) <<  meanTot_la << " GeV}{#sigma : " << rmsSigmaTot_la << " GeV}";
        ltx1->DrawLatex( 0.65, 0.75, os.str(  ).c_str(  ) );
        os.str( std::string(  ) ); // clears ostringstream
        os << "185 #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220";
        ltx1->DrawLatex( 0.23, 0.67, os.str(  ).c_str(  ) );
        os.str( std::string(  ) );
        ltx2->DrawLatex( 0.23, 0.60, "#Lambda/ #overline{#Lambda} " );

        TLegend* leg = new TLegend( 0.6,0.5,0.85,0.6);
        leg->AddEntry( "fitFcn", "Complete fit function", "1L" );
        leg->AddEntry( "backFcn", "Background fit", "1L" );
        leg->SetBorderSize( 0 );
        leg->SetFillColor( 0 );
        leg->Draw(  );
    }
}


