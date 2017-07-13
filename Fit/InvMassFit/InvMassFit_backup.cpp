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
Double_t xi_bgfunc( Double_t *x, Double_t *par )
{
    Double_t xx = x[0];
    Double_t bg = par[0]*TMath::Power(xx - ( piMass + lambdaMass ),par[1]);

    return bg;
}

Double_t V0bgfunc( Double_t *x, Double_t *par )
{
    Double_t xx = x[0];
    Double_t bg = par[0] + par[1]*xx*xx + par[2]*TMath::Power( xx, 3 ) + par[3]*TMath::Power( xx,4 ); 
    return bg;
}


Double_t xi_bkgfuncdisplay( Double_t *x, Double_t *par )
{
    Double_t bg = par[0]*TMath::Power( x[0] - ( piMass + lambdaMass ), par[1] );
    Double_t bgsig = par[2]*TMath::Power( x[0] - ( piMass + lambdaMass ), par[3] );
    return bg + bgsig;
}

//Define double Gaussian fit function
Double_t xi_sgfunc( Double_t *x, Double_t *par )
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

Double_t V0_sgfunc( Double_t *x, Double_t *par )
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

Double_t xi_total( Double_t *x, Double_t *par )
{
    return xi_bgfunc( x,par ) + xi_sgfunc( x,&par[2] );
}

Double_t V0_total( Double_t *x, Double_t *par )
{
    return V0bgfunc( x,par ) + V0_sgfunc( x,&par[4] );
}

void InvMassFit(  )
{
    using namespace std;
    /*
    int numFiles = 2;
    // Make sure to change numFiles
    const char* FN[] = {
        "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL1.root",
        "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL2.root"
    };
    */
    int numPtBins = 8;
    double p[] = {1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0}; //if number of bins changes make sure you change numPtBins

    std::ostringstream os;
    TH1::SetDefaultSumw2(  );

    //Containers for Big for loop
    std::vector< std::vector<double> > Cont_MassMean;
    std::vector< std::vector<double> > Cont_MassStd;
    std::vector< std::vector<double> > Cont_Fsig;
    std::vector< std::vector<double> > Cont_sig;
    std::vector< std::vector<TH1D*> > Cont_InvMassPtBinned;
    std::vector<TH1D*> Cont_PtSpectra;
    std::vector<TF1*> BgFit;
    std::vector<TF1*> SigFit;
    std::vector<TF1*> FitTot;
    std::vector<TF1*> FitFcn;
    std::vector<TF1*> BackFcn;
    std::vector<TCanvas*> cHist;
    std::vector<TCanvas*> cPt;
    std::vector<TCanvas*> cMass;
    //

    //Fit for total invariant mass peak
    std::vector<TH2D*> MassPt;
    std::vector<TF1*> bgfit;
    std::vector<TF1*> sigfit; 
    std::vector<TH1D*> InvMass;
    std::vector<TF1*> fitTot;
    std::vector<TF1*> fitFcn;
    std::vector<TF1*> backFcn;
    std::vector<double> rmsSigmaTot;
    std::vector<double> meanTot;
    //

    std::vector<TFile*> FILES;
    //std::vector<const char*> Fn( FN, FN+numFiles );
    std::vector<const char*> Fn;
    std::vector<double> PtBin( p, p+numPtBins );


    //Fn.push_back( "/Volumes/MacHD/Users/blt1/research/TestRootFiles/XiAnalysisCorrelationPtCut8TeVPD1_4_ForFinal.root" );
    
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL4.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL5.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL6.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL7.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL8.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL9.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL11.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL12.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL13.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL14.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL15.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL16.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL17.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL18.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL19.root" );
    Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL20.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL21.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL22.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL23.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL24.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL25.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL26.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL27.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL28.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL29.root" );
    //Fn.push_back( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL20.root" );
    
        

    int PtBinSize = PtBin.size(  ) - 1;

    for( std::vector<const char*>::iterator it=Fn.begin(  ); it!=Fn.end(  ); ++it ){
        TFile* f = new TFile( *it );
        FILES.push_back( f );
    }

    cout << PtBin.size(  ) << endl;
    
    //TFile* f = new TFile( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL5.root" );

    //TH1D* InvMassPtBinned[PtBinSize];
    bool publish = false; // Adds Latex labels

    /*
    TF1* BgFit[PtBinSize];
    TF1* SigFit[PtBinSize];
    TF1* FitTot[PtBinSize];
    TF1* FitFcn[PtBinSize]; //For drawing and fsig calc
    TF1* BackFcn[PtBinSize];
    */

    //Big For Loop HERE
    //
    for( unsigned j=0; j<FILES.size(  ); j++ )
    {
        std::vector<TH1D*> InvMassPtBinned;
        /*
           std::vector<TF1*> BgFit;
           std::vector<TF1*> SigFit;
           std::vector<TF1*> FitTot;
           std::vector<TF1*> FitFcn;
           std::vector<TF1*> BackFcn;
           */
        std::vector<double> Fsig;
        std::vector<double> sig;
        std::vector<double> MassMean;
        std::vector<double> MassStd;

        MassPt.push_back( ( TH2D* )FILES[j]->Get( "xiMassPt/MassPt" ) );
        //MassPt.push_back( ( TH2D* )FILES[j]->Get( "xiCorrelation/MassPt" ) );

        TVirtualFitter::SetMaxIterations( 300000 );

        //Fit Total Inv Mass

        InvMass.push_back( MassPt[j]->ProjectionX( Form( "InvMass_%d",j ), 0, 100 ) );
        InvMass[j]->SetMarkerStyle( 21 );
        InvMass[j]->SetMarkerSize( 0.5 );

        TF1* bgfit_ = new TF1( Form( "bgfit_%d",j ), xi_bgfunc, 1.26,1.305,2 );
        bgfit_->SetNpx( 60 );
        bgfit_->SetParameter( 0,3.17496e1 );
        bgfit_->SetParameter( 1,5.97019e-1 );
        bgfit_->SetParNames( "bkgScale","bkgPow" );

        TF1* sigfit_ = new TF1( Form( "sigfit_%d",j ), xi_sgfunc,1.26,1.4,5 );
        //sigfit_->SetNpx( 50 );
        sigfit_->SetNpx( 250 );
        sigfit_->SetParNames( "gaus1Norm","gausMean","gaus1std","gaus2Norm","gaus2std","sbkgScale","sbkgPow" );
        sigfit_->SetParameter( 0,6.6688e3 );
        sigfit_->SetParameter( 1,1.32219 );
        sigfit_->SetParameter( 2,3.15655e-3 );
        //sigfit_->SetParameter( 2,6.15655e-3 ); //for pT range 1-2
        sigfit_->SetParameter( 3,3.67793e3 );
        sigfit_->SetParameter( 4,3.69916e-3 );
        sigfit_->SetParameter( 5,3.21847e3 );
        sigfit_->SetParameter( 6,5.99562-1 );

        TF1* fitTot_ = new TF1( Form( "fitTot_%d",j ), xi_total, 1.26, 1.4, 9 );
        fitTot_->SetNpx( 250 );
        fitTot_->SetParNames( "bkgScale","bkgPow","gaus1Norm","gausMean","gaus1std","gaus2Norm","gaus2std","sbkgScale","sbkgPow" );

        //for test function
        Double_t partestTot[9];

        InvMass[j]->Fit( Form( "bgfit_%d",j ),"L","",1.26,1.31 );
        InvMass[j]->Fit( Form( "sigfit_%d",j ), "L R","",1.26,1.4 ); //1.335

        bgfit_->GetParameters( &partestTot[0] );
        sigfit_->GetParameters( &partestTot[2] );

        fitTot_->SetParameters( partestTot );


        InvMass[j]->SetTitle( "" );
        InvMass[j]->GetYaxis(  )->SetTitleOffset( 1.5 );
        InvMass[j]->GetYaxis(  )->SetTitleSize( 0.035 );
        InvMass[j]->GetYaxis(  )->CenterTitle( true );
        InvMass[j]->GetYaxis(  )->SetTitle( "Candidates/0.001 GeV" );
        InvMass[j]->GetYaxis(  )->SetRangeUser( 0,34000 );
        InvMass[j]->SetTitleOffset( 1.5, "X" );
        InvMass[j]->GetXaxis(  )->SetTitleSize( 0.035 );
        InvMass[j]->GetXaxis(  )->CenterTitle( true );
        InvMass[j]->GetXaxis(  )->SetTitle( "#Lambda #pi^{#minus} invariant mass (GeV)" );
        InvMass[j]->SetStats( !publish );
        //Drawing the fit function again for the sake of getting it onto the legend
        InvMass[j]->Fit( Form( "fitTot_%d",j ), "L0","", 1.26, 1.4 );
        Double_t parfordrawTot[9];
        fitTot_->GetParameters( &parfordrawTot[0] );
        TF1* fitFcn_ = new TF1( Form( "fitFcn_%d",j ), xi_total, 1.26,1.4,9 );
        fitFcn_->SetNpx( 250 );
        fitFcn_->SetParameters( parfordrawTot );
        fitFcn_->SetLineColor( kRed );
        TF1* backFcn_ = new TF1( Form( "backFcn_%d",j ), xi_bkgfuncdisplay, 1.26,1.4,4 );
        backFcn_->SetParameter( 0,fitTot_->GetParameter( 0 ) );
        backFcn_->SetParameter( 1,fitTot_->GetParameter( 1 ) );
        backFcn_->SetParameter( 2,fitTot_->GetParameter( 7 ) );
        backFcn_->SetParameter( 3,fitTot_->GetParameter( 8 ) );
        backFcn_->SetLineColor( kBlue );

        rmsSigmaTot.push_back( TMath::Sqrt( 0.5*fitTot_->GetParameter( "gaus1std" )*fitTot_->GetParameter( "gaus1std" ) + 0.5*fitTot_->GetParameter( "gaus2std" )*fitTot_->GetParameter( "gaus2std" ) ) );
        meanTot.push_back( fitTot_->GetParameter( "gausMean" ) );


        bgfit.push_back( bgfit_ );
        sigfit.push_back( sigfit_ ); 
        fitTot.push_back( fitTot_ );
        fitFcn.push_back( fitFcn_ );
        backFcn.push_back( backFcn_ );



        for( int i=0; i<PtBinSize; i++ ){
            InvMassPtBinned.push_back( MassPt[j]->ProjectionX( Form( "InvMass_pT_%d",j*PtBinSize + i ), PtBin[i]*10+1, PtBin[i+1]*10 ) );
            os << PtBin[i] << "_Pt_" << PtBin[i+1];
            InvMassPtBinned[i]->SetTitle( os.str(  ).c_str(  ) );
            os.str( std::string(  ) );

            //Perform Fits
            BgFit.push_back( new TF1( Form( "bgfit_%d",j*PtBinSize + i ), xi_bgfunc, 1.26,1.305,2 ) );
            BgFit[j*PtBinSize + i]->SetNpx( 60 );
            BgFit[j*PtBinSize + i]->SetParameter( 0,3.17496e1 );
            BgFit[j*PtBinSize + i]->SetParameter( 1,5.97019e-1 );
            BgFit[j*PtBinSize + i]->SetParNames( "bkgScale","bkgPow" );

            SigFit.push_back( new TF1( Form( "sigfit_%d",j*PtBinSize + i ), xi_sgfunc, 1.26,1.4,7 ) );
            SigFit[j*PtBinSize + i]->SetNpx( 250 );
            SigFit[j*PtBinSize + i]->SetParNames( "gaus1Norm","gausMean","gaus1std","gaus2Norm","gaus2std","sbkgScale","sbkgPow" );
            SigFit[j*PtBinSize + i]->SetParameter( 0,6.6688e3 );
            SigFit[j*PtBinSize + i]->SetParameter( 1,1.32219 );
            SigFit[j*PtBinSize + i]->SetParameter( 2,3.15655e-3 );
            SigFit[j*PtBinSize + i]->SetParameter( 2,6.15655e-3 ); //for pT range 1-2
            SigFit[j*PtBinSize + i]->SetParameter( 3,3.67793e3 );
            SigFit[j*PtBinSize + i]->SetParameter( 4,3.69916e-3 );
            SigFit[j*PtBinSize + i]->SetParameter( 5,3.21847e3 );
            SigFit[j*PtBinSize + i]->SetParameter( 6,5.99562-1 );

            FitTot.push_back( new TF1( Form( "fitTot_%d",j*PtBinSize + i ), xi_total, 1.26, 1.4, 9 ) );
            FitTot[j*PtBinSize + i]->SetNpx( 250 );
            FitTot[j*PtBinSize + i]->SetParNames( "bkgScale","bkgPow","gaus1Norm","gausMean","gaus1std","gaus2Norm","gaus2std","sbkgScale","sbkgPow" );

            Double_t partest[9];
            InvMassPtBinned[i]->Fit( Form( "bgfit_%d",j*PtBinSize + i ), "L","",1.26,1.31 );
            InvMassPtBinned[i]->Fit( Form( "sigfit_%d",j*PtBinSize + i ), "L R", "",1.26,1.4 );

            BgFit[j*PtBinSize + i]->GetParameters( &partest[0] );
            SigFit[j*PtBinSize + i]->GetParameters( &partest[2] );

            FitTot[j*PtBinSize + i]->SetParameters( partest );

            InvMassPtBinned[i]->Fit( Form( "fitTot_%d",j*PtBinSize + i ), "L","",1.26,1.4 );

            //Draw signal and background fit also need bkg fit for fsig calc
            Double_t parfordraw[9];
            FitTot[j*PtBinSize + i]->GetParameters( &parfordraw[0] );
            FitFcn.push_back( new TF1( Form( "FitFcn_%d",j*PtBinSize + i ), xi_total, 1.26,1.4,9 ) );
            FitFcn[j*PtBinSize + i]->SetNpx( 250 );
            FitFcn[j*PtBinSize + i]->SetParameters( parfordraw );
            FitFcn[j*PtBinSize + i]->SetLineColor( kRed );
            FitFcn[j*PtBinSize + i]->Draw( "same" );
            BackFcn.push_back( new TF1( Form( "backFcn_%d",j*PtBinSize + i ), xi_bkgfuncdisplay, 1.26,1.4,4 ) );
            BackFcn[j*PtBinSize + i]->SetParameter( 0,FitTot[j*PtBinSize + i]->GetParameter( 0 ) );
            BackFcn[j*PtBinSize + i]->SetParameter( 1,FitTot[j*PtBinSize + i]->GetParameter( 1 ) );
            BackFcn[j*PtBinSize + i]->SetParameter( 2,FitTot[j*PtBinSize + i]->GetParameter( 7 ) );
            BackFcn[j*PtBinSize + i]->SetParameter( 3,FitTot[j*PtBinSize + i]->GetParameter( 8 ) );
            BackFcn[j*PtBinSize + i]->SetLineColor( kBlue );
            BackFcn[j*PtBinSize + i]->Draw( "same" );

            //Calculate values

            //Sigma
            double rmsSigma = TMath::Sqrt( 0.5*FitTot[j*PtBinSize + i]->GetParameter( "gaus1std" )*FitTot[j*PtBinSize + i]->GetParameter( "gaus1std" ) + 0.5*FitTot[j*PtBinSize + i]->GetParameter( "gaus2std" )*FitTot[j*PtBinSize + i]->GetParameter( "gaus2std" ) );            

            //Signal fraction yield
            double totInt = -999.0;
            double bkgInt = -999.0;

            totInt = FitTot[j*PtBinSize + i]->Integral( FitTot[j*PtBinSize + i]->GetParameter( "gausMean" ) - 2*rmsSigma, FitTot[j*PtBinSize + i]->GetParameter( "gausMean" ) + 2*rmsSigma );

            bkgInt = BackFcn[j*PtBinSize + i]->Integral( FitTot[j*PtBinSize + i]->GetParameter( "gausMean" ) - 2*rmsSigma, FitTot[j*PtBinSize + i]->GetParameter( "gausMean" ) + 2*rmsSigma );

            Fsig.push_back( ( totInt - bkgInt )/totInt );
            sig.push_back( totInt - bkgInt );
            MassMean.push_back( FitTot[j*PtBinSize + i]->GetParameter( "gausMean" ) );
            MassStd.push_back( rmsSigma );
        }
        Cont_Fsig.push_back( Fsig );
        Cont_sig.push_back( sig );
        Cont_MassMean.push_back( MassMean );
        Cont_MassStd.push_back( MassStd );
        Cont_InvMassPtBinned.push_back( InvMassPtBinned );

        // Dont need these, just keep all in one big vector hence j*PtBinSize + i
        // These don't need a vector of vector because these functions aren't
        // actually called for later
        //Cont_BgFit[j].push_back( BgFit ); 
        //Cont_SigFit[j].push_back( SigFit );
        //Cont_FitTot[j].push_back( FitTot );

    }
    //End for loop HERE

    // Another for loop to print out everything at the end of the output
    for( unsigned j=0; j<FILES.size(  ); j++ )
    {
        int PtBinCounter=0;

        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << "PRINTING FOR FILE: " << Fn[j] << endl;
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

        // Sig yield values
        cout << "==========================================================" << endl;
        cout << "Signal Yield" << endl;
        cout << "==========================================================" << endl;

        for( std::vector<double>::iterator it = Cont_sig[j].begin(  ); it != Cont_sig[j].end(  ); ++it ){ 
            cout << PtBin[PtBinCounter] << " < Pt =< " << PtBin[PtBinCounter+1] << ": " << *it << endl;
            PtBinCounter++;
        }

        // Fsig values
        cout << "==========================================================" << endl;
        cout << "Signal Fraction Yield" << endl;
        cout << "==========================================================" << endl;
        PtBinCounter = 0;
        for( std::vector<double>::iterator it = Cont_Fsig[j].begin(  ); it != Cont_Fsig[j].end(  ); ++it ){ 
            cout << PtBin[PtBinCounter] << " < Pt =< " << PtBin[PtBinCounter+1] << ": " << *it << endl;
            PtBinCounter++;
        }

        // MassMean
        cout << "==========================================================" << endl;
        cout << "Mean Mass" << endl;
        cout << "==========================================================" << endl;

        PtBinCounter = 0;
        for( std::vector<double>::iterator it = Cont_MassMean[j].begin(  ); it != Cont_MassMean[j].end(  ); ++it ){ 
            cout << PtBin[PtBinCounter] << " < Pt =< " << PtBin[PtBinCounter+1] << ": " << *it << endl;
            PtBinCounter++;
        }

        // RMS Sigma
        cout << "==========================================================" << endl;
        cout << "RMS Sigma" << endl;
        cout << "==========================================================" << endl;

        PtBinCounter = 0;
        for( std::vector<double>::iterator it = Cont_MassStd[j].begin(  ); it != Cont_MassStd[j].end(  ); ++it ){ 
            cout << PtBin[PtBinCounter] << " < Pt =< " << PtBin[PtBinCounter+1] << ": " << *it << endl;
            PtBinCounter++;
        }
    }

    // Draw Histgorams
    int JLcounter=4;
    for( unsigned j=0; j<FILES.size(  ); j++ )
    {
        if( JLcounter==10 ) {JLcounter++; continue; }
        TCanvas* Hist = new TCanvas( Form( "Canvas_%d",j ), Form( "JL%d",JLcounter ),1600,800 );
        Hist->Divide( 4,2 );
        for( int i=0; i<PtBinSize; i++ )
        {
            Hist->cd( i+1 );
            Cont_InvMassPtBinned[j][i]->Draw( "E1" );
        }
        cHist.push_back( Hist );
        /*
        TImage *img = TImage::Create(  );
        img->FromPad( Hist );
        img->WriteImage( Form( "PtBinned_JL%d.png",JLcounter ) );
        */

        //Hist->Print( Form( "PtBinned_JL%d.png",JLcounter ) );

        TH1D* PtSpectra = ( TH1D* )MassPt[j]->ProjectionY( Form( "Pt Spectrum_%d",j ),0,150 );
        TCanvas* Pt = new TCanvas( Form("Pt_%d",j), Form( "JL%d",JLcounter ), 800, 800 );
        Pt->cd(  );
        PtSpectra->GetXaxis(  )->SetRangeUser( 0,2 );
        PtSpectra->Draw( "E1" );
        Cont_PtSpectra.push_back( PtSpectra );
        cPt.push_back( Pt );

        PtSpectra->Print( Form( "PtSpectra_JL%d",JLcounter ) );

        TCanvas *Mass = new TCanvas( Form( "cMass_%d",j ), Form( "JL%d",JLcounter ), 800,800 );
        Mass->SetFillColor      (0);
        Mass->SetBorderMode     (0);
        Mass->SetBorderSize     (10);
        Mass->SetLeftMargin     (0.20);
        Mass->SetRightMargin    (0.06);
        Mass->SetTopMargin      (0.11);
        Mass->SetBottomMargin   (0.15);
        Mass->SetFrameFillStyle (0);
        Mass->SetFrameLineStyle (0);
        Mass->SetFrameBorderMode(0);
        Mass->SetFrameBorderSize(10);
        Mass->SetFrameFillStyle (0);
        Mass->SetFrameLineStyle (0);
        Mass->SetFrameBorderMode(0);
        Mass->SetFrameBorderSize(10);

        TGaxis::SetMaxDigits( 3 );
        Mass->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        InvMass[j]->Draw( "E1 9" );
        fitFcn[j]->Draw( "same" );
        backFcn[j]->Draw( "same" );




        if( publish ) 
        {
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
            os << "#splitline{Mean: " << fixed << std::setprecision( 4 ) <<  meanTot[j] << " GeV}{#sigma : " << rmsSigmaTot[j] << " GeV}";
            ltx1->DrawLatex( 0.65, 0.75, os.str(  ).c_str(  ) );
            os.str( std::string(  ) ); // clears ostringstream
            os << "185 #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220";
            ltx1->DrawLatex( 0.23, 0.67, os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            ltx2->DrawLatex( 0.23, 0.60, "#Xi^{#plus}/ #Xi^{#minus} " );
            //ltx1->DrawLatex( 0.6, 0.75, "#splitline{Mean: 1.3222 GeV}{#sigma : 0.0053 GeV}" );
            //ltx1->DrawLatex( 0.6, 0.75, tLatex( "#splitline{Mean:", "GeV}{#sigma :", testint, 0.0053  ));

            //ltx1->DrawLatex( 0.2, 0.82, "CMS pPb  #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
            //ltx1->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
            //ltx1->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            //ltx2->DrawLatex( 0.2, 0.60, "#Xi^{#plus}/ #Xi^{#minus} " );

            TLegend* leg = new TLegend( 0.6,0.5,0.85,0.6);
            leg->AddEntry( "fitFcn", "Complete fit function", "1L" );
            leg->AddEntry( "backFcn", "Background fit", "1L" );
            leg->SetBorderSize( 0 );
            leg->SetFillColor( 0 );
            leg->Draw(  );
        }
        cMass.push_back( Mass );
        JLcounter++;
    }

}

