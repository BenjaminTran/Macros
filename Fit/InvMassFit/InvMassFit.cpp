#include <TLatex.h>
#include <TStyle.h>
#include "TF1.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TLegend.h"
#include "TMathText.h"

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
Double_t bgfunc1( Double_t *x, Double_t *par )
{
    Double_t xx = x[0];
    Double_t bg = par[0]*TMath::Power(xx - ( piMass + lambdaMass ),par[1]);

    return bg;
}

Double_t bkgfuncdisplay( Double_t *x, Double_t *par )
{
    Double_t bg = par[0]*TMath::Power( x[0] - ( piMass + lambdaMass ), par[1] );
    Double_t bgsig = par[2]*TMath::Power( x[0] - ( piMass + lambdaMass ), par[3] );
    return bg + bgsig;
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
    Double_t bgsig  = par[5]*TMath::Power( x[0] - ( piMass + lambdaMass ),par[6] );
    return sgaus1 + sgaus2 + bgsig;
}

Double_t total( Double_t *x, Double_t *par )
{
    return bgfunc1( x,par ) + sgfunc( x,&par[2] );
}


void InvMassFit(  )
{
    std::ostringstream os;
    TH1::SetDefaultSumw2(  );

    TFile* f = new TFile( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/8TeV/MassPt8TeVPD1.root" );
    //TFile* f = new TFile( "/volumes/MacHD/Users/blt1/research/TestRootFiles/MassPt8TeV_1.root" );


    int numPtBins = 7;
    TH1D* InvMassPtBinned[numPtBins];
    //TH1D* PtSpecra;
    double p[] = {1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0}; //if number of bins changes make sure you change numPtBins
    std::vector<double> PtBin( p, p+8 );
    int PtBinSize = PtBin.size(  ) - 1;
    //double MassMean[numPtBins];
    //double MassStd[numPtBins];
    std::vector<double> MassMean;
    std::vector<double> MassStd;
    std::vector<double> Fsig;

    TF1* BgFit[numPtBins];
    TF1* SigFit[numPtBins];
    TF1* FitTot[numPtBins];
    TF1* FitFcn[numPtBins]; //For drawing and fsig calc
    TF1* BackFcn[numPtBins];

    TH2D* MassPt = ( TH2D* )f->Get( "xiMassPt/MassPt" );

    TVirtualFitter::SetMaxIterations( 300000 );

    for( int i=0; i<PtBinSize; i++ ){
        cout << i << endl;
            InvMassPtBinned[i] = MassPt->ProjectionX( Form( "InvMass_pT_%d",i ), PtBin[i]*10+1, PtBin[i+1]*10 );
            os << PtBin[i] << "_Pt_" << PtBin[i+1];
            InvMassPtBinned[i]->SetTitle( os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            
            //Perform Fits
            BgFit[i] = new TF1( Form( "bgfit_%d",i ), bgfunc1, 1.26,1.305,2 );
            BgFit[i]->SetNpx( 60 );
            BgFit[i]->SetParameter( 0,3.17496e1 );
            BgFit[i]->SetParameter( 1,5.97019e-1 );
            BgFit[i]->SetParNames( "bkgScale","bkgPow" );

            SigFit[i] = new TF1( Form( "sigfit_%d",i ), sgfunc, 1.26,1.4,7 );
            SigFit[i]->SetNpx( 250 );
            SigFit[i]->SetParNames( "gaus1Norm","gausMean","gaus1std","gaus2Norm","gaus2std","sbkgScale","sbkgPow" );
            SigFit[i]->SetParameter( 0,6.6688e3 );
            SigFit[i]->SetParameter( 1,1.32219 );
            SigFit[i]->SetParameter( 2,3.15655e-3 );
            SigFit[i]->SetParameter( 2,6.15655e-3 ); //for pT range 1-2
            SigFit[i]->SetParameter( 3,3.67793e3 );
            SigFit[i]->SetParameter( 4,3.69916e-3 );
            SigFit[i]->SetParameter( 5,3.21847e3 );
            SigFit[i]->SetParameter( 6,5.99562-1 );

            FitTot[i] = new TF1( Form( "fitTot_%d",i ), total, 1.26, 1.4, 9 );
            FitTot[i]->SetNpx( 250 );
            FitTot[i]->SetParNames( "bkgScale","bkgPow","gaus1Norm","gausMean","gaus1std","gaus2Norm","gaus2std","sbkgScale","sbkgPow" );

            
            Double_t partest[9];
            InvMassPtBinned[i]->Fit( Form( "bgfit_%d",i ), "L","",1.26,1.31 );
            InvMassPtBinned[i]->Fit( Form( "sigfit_%d",i ), "L R", "",1.26,1.4 );

            BgFit[i]->GetParameters( &partest[0] );
            SigFit[i]->GetParameters( &partest[2] );
            
            FitTot[i]->SetParameters( partest );

            InvMassPtBinned[i]->Fit( Form( "fitTot_%d",i ), "L","",1.26,1.4 );

            //Draw signal and background fit also need bkg fit for fsig calc
            Double_t parfordraw[9];
            FitTot[i]->GetParameters( &parfordraw[0] );
            FitFcn[i] = new TF1( Form( "FitFcn_%d",i ), total, 1.26,1.4,9 );
            FitFcn[i]->SetNpx( 250 );
            FitFcn[i]->SetParameters( parfordraw );
            FitFcn[i]->SetLineColor( kRed );
            FitFcn[i]->Draw( "same" );
            BackFcn[i] = new TF1( Form( "backFcn_%d",i ), bkgfuncdisplay, 1.26,1.4,4 );
            BackFcn[i]->SetParameter( 0,FitTot[i]->GetParameter( 0 ) );
            BackFcn[i]->SetParameter( 1,FitTot[i]->GetParameter( 1 ) );
            BackFcn[i]->SetParameter( 2,FitTot[i]->GetParameter( 7 ) );
            BackFcn[i]->SetParameter( 3,FitTot[i]->GetParameter( 8 ) );
            BackFcn[i]->SetLineColor( kBlue );
            BackFcn[i]->Draw( "same" );

            //Calculate values

            //Sigma
            double rmsSigma = TMath::Sqrt( 0.5*FitTot[i]->GetParameter( "gaus1std" )*FitTot[i]->GetParameter( "gaus1std" ) + 0.5*FitTot[i]->GetParameter( "gaus2std" )*FitTot[i]->GetParameter( "gaus2std" ) );            

            //Signal fraction yield
            double totInt = -999.0;
            double bkgInt = -999.0;

            totInt = FitTot[i]->Integral( FitTot[i]->GetParameter( "gausMean" ) - 2*rmsSigma, FitTot[i]->GetParameter( "gausMean" ) + 2*rmsSigma );

            bkgInt = BackFcn[i]->Integral( FitTot[i]->GetParameter( "gausMean" ) - 2*rmsSigma, FitTot[i]->GetParameter( "gausMean" ) + 2*rmsSigma );

            //Fsig[i]     = ( totInt - bkgInt )/totInt;
            Fsig.insert(Fsig.begin(  ), ( totInt - bkgInt )/totInt );
            MassMean.insert( MassMean.begin(  ), FitTot[i]->GetParameter( "gausMean" ) );
            MassStd.insert( MassStd.begin(  ), rmsSigma );
    }

    int PtBinCounter=0;

    // Fsig values
    cout << "==========================================================" << endl;
    cout << "Signal Fraction Yield" << endl;
    cout << "==========================================================" << endl;

    for( std::vector<double>::iterator it = Fsig.begin(  ); it != Fsig.end(  ); ++it ){ 
        cout << PtBin[PtBinCounter] << " < Pt =< " << PtBin[PtBinCounter+1] << ": " << *it << endl;
        PtBinCounter++;
    }

    // MassMean
    cout << "==========================================================" << endl;
    cout << "Mean Mass" << endl;
    cout << "==========================================================" << endl;

    PtBinCounter = 0;
    for( std::vector<double>::iterator it = MassMean.begin(  ); it != MassMean.end(  ); ++it ){ 
        cout << PtBin[PtBinCounter] << " < Pt =< " << PtBin[PtBinCounter+1] << ": " << *it << endl;
        PtBinCounter++;
    }

    cout << "==========================================================" << endl;
    cout << "RMS Sigma" << endl;
    cout << "==========================================================" << endl;

    PtBinCounter = 0;
    for( std::vector<double>::iterator it = MassStd.begin(  ); it != MassStd.end(  ); ++it ){ 
        cout << PtBin[PtBinCounter] << " < Pt =< " << PtBin[PtBinCounter+1] << ": " << *it << endl;
        PtBinCounter++;
    }

    // Draw Histgorams
    TCanvas* c1 = new TCanvas( "Canvas","",1600,800 );
    c1->Divide( 4,2 );
    for( int i=0; i<numPtBins; i++ ){
        c1->cd( i+1 );
        InvMassPtBinned[i]->Draw( "E1" );
    }

    TH1D* PtSpecra = ( TH1D* )MassPt->ProjectionY( "Pt Spectrum",0,150 );
    TCanvas* c2 = new TCanvas( "Canvas2","", 800, 800 );
    c2->cd(  );
    PtSpecra->Draw( "E1" );

    /*
    for( int i=0; i<7; i++ ){
    
    }
    */

    /*
    TCanvas* c1 = new TCanvas( "Canvas","",1600,800 );
    c1->Divide( 4,2 );
    c1->cd( 1 );
    InvMassPtBinned[0]->Draw(  );
    c1->cd( 2 );
    InvMassPtBinned[1]->Draw(  );
    c1->cd( 3 );
    InvMassPtBinned[2]->Draw(  );
    c1->cd( 4 );
    InvMassPtBinned[3]->Draw(  );
    c1->cd( 5 );
    InvMassPtBinned[4]->Draw(  );
    c1->cd( 6 );
    InvMassPtBinned[5]->Draw(  );
    c1->cd( 7 );
    InvMassPtBinned[6]->Draw(  );
    */
}

