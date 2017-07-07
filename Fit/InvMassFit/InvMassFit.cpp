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

    //TFile* f = new TFile( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/8TeV/MassPt8TeVPD1.root" );
    //TFile* f = new TFile( "/volumes/MacHD/Users/blt1/research/TestRootFiles/XiAnalysisCorrelationPtCut8TeVPD1_4_Inc.root" );
    //TFile* f = new TFile( "/volumes/MacHD/Users/blt1/research/TestRootFiles/XiAnalysisCorrelationPtCut8TeVPD1_4_ForFinal.root" );
    //TFile* f = new TFile( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL1.root" );
    //TFile* f = new TFile( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL2.root" );
    //TFile* f = new TFile( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL3.root" );
    //TFile* f = new TFile( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL4.root" );
    TFile* f = new TFile( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/Flow/CasCutLoose/CasCutLooseJL5.root" );

    int numPtBins = 7;
    int numFiles = 5;

    TFile* Files[numFiles];

    TH1D* InvMassPtBinned[numPtBins];
    double p[] = {1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0}; //if number of bins changes make sure you change numPtBins
    std::vector<double> PtBin( p, p+8 );
    std::vector<double> MassMean;
    std::vector<double> MassStd;
    std::vector<double> Fsig;
    std::vector<double> sig;
    int PtBinSize = PtBin.size(  ) - 1;
    bool publish = true; // Adds Latex labels

    TF1* BgFit[numPtBins];
    TF1* SigFit[numPtBins];
    TF1* FitTot[numPtBins];
    TF1* FitFcn[numPtBins]; //For drawing and fsig calc
    TF1* BackFcn[numPtBins];

    TH2D* MassPt = ( TH2D* )f->Get( "xiMassPt/MassPt" );
    //TH2D* MassPt = ( TH2D* )f->Get( "xiCorrelation/MassPt" );

    TVirtualFitter::SetMaxIterations( 300000 );

    //Fit Total Inv Mass

    TH1D* InvMass = MassPt->ProjectionX( "InvMass", 0, 100 );
    InvMass->SetMarkerStyle( 21 );
    InvMass->SetMarkerSize( 0.5 );

    TF1 *bgfit = new TF1( "bgfit", bgfunc1, 1.26,1.305,2 );
    //TF1 *bgfit = new TF1( "bgfit", bgfunc1, 1.26,1.4,2 );
    bgfit->SetNpx( 60 );
    bgfit->SetParameter( 0,3.17496e1 );
    bgfit->SetParameter( 1,5.97019e-1 );
    bgfit->SetParNames( "bkgScale","bkgPow" );

    TF1 *sigfit = new TF1( "sigfit", sgfunc,1.26,1.4,7 ); 
    //TF1 *sigfit = new TF1( "sigfit", sgfunc,1.26,1.4,5 );
    //sigfit->SetNpx( 50 );
    sigfit->SetNpx( 250 );
    sigfit->SetParNames( "gaus1Norm","gausMean","gaus1std","gaus2Norm","gaus2std","sbkgScale","sbkgPow" );
    sigfit->SetParameter( 0,6.6688e3 );
    sigfit->SetParameter( 1,1.32219 );
    sigfit->SetParameter( 2,3.15655e-3 );
    //sigfit->SetParameter( 2,6.15655e-3 ); //for pT range 1-2
    sigfit->SetParameter( 3,3.67793e3 );
    sigfit->SetParameter( 4,3.69916e-3 );
    sigfit->SetParameter( 5,3.21847e3 );
    sigfit->SetParameter( 6,5.99562-1 );
    
    TF1 *fitTot = new TF1( "fitTot", total, 1.26, 1.4, 9 );
    fitTot->SetNpx( 250 );
    fitTot->SetParNames( "bkgScale","bkgPow","gaus1Norm","gausMean","gaus1std","gaus2Norm","gaus2std","sbkgScale","sbkgPow" );

    //for test function
    Double_t partestTot[9];

    InvMass->Fit( "bgfit","L","",1.26,1.31 );
    InvMass->Fit( "sigfit", "L R","",1.26,1.4 ); //1.335

    bgfit->GetParameters( &partestTot[0] );
    sigfit->GetParameters( &partestTot[2] );


    fitTot->SetParameters( partestTot );

    InvMass->SetTitle( "" );
    InvMass->GetYaxis(  )->SetTitleOffset( 1.5 );
    InvMass->GetYaxis(  )->SetTitleSize( 0.035 );
    InvMass->GetYaxis(  )->CenterTitle( true );
    InvMass->GetYaxis(  )->SetTitle( "Candidates/0.001 GeV" );
    InvMass->GetYaxis(  )->SetRangeUser( 0,34000 );
    InvMass->SetTitleOffset( 1.5, "X" );
    InvMass->GetXaxis(  )->SetTitleSize( 0.035 );
    InvMass->GetXaxis(  )->CenterTitle( true );
    InvMass->GetXaxis(  )->SetTitle( "#Lambda #pi^{#minus} invariant mass (GeV)" );
    InvMass->SetStats( !publish );
    //Drawing the fit function again for the sake of getting it onto the legend
    InvMass->Fit( "fitTot", "L0","", 1.26, 1.4 );
    Double_t parfordrawTot[9];
    fitTot->GetParameters( &parfordrawTot[0] );
    TF1 *fitFcn = new TF1( "fitFcn", total, 1.26,1.4,9 );
    fitFcn->SetNpx( 250 );
    fitFcn->SetParameters( parfordrawTot );
    fitFcn->SetLineColor( kRed );
    TF1 *backFcn = new TF1( "backFcn", bkgfuncdisplay, 1.26,1.4,4 );
    backFcn->SetParameter( 0,fitTot->GetParameter( 0 ) );
    backFcn->SetParameter( 1,fitTot->GetParameter( 1 ) );
    backFcn->SetParameter( 2,fitTot->GetParameter( 7 ) );
    backFcn->SetParameter( 3,fitTot->GetParameter( 8 ) );
    backFcn->SetLineColor( kBlue );

    double rmsSigmaTot = TMath::Sqrt( 0.5*fitTot->GetParameter( "gaus1std" )*fitTot->GetParameter( "gaus1std" ) + 0.5*fitTot->GetParameter( "gaus2std" )*fitTot->GetParameter( "gaus2std" ) );
    double meanTot = fitTot->GetParameter( "gausMean" );



    for( int i=0; i<PtBinSize; i++ ){
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
 ;
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

            Fsig.push_back( ( totInt - bkgInt )/totInt );
            sig.push_back( totInt - bkgInt );
            MassMean.push_back( FitTot[i]->GetParameter( "gausMean" ) );
            MassStd.push_back( rmsSigma );
    }

    int PtBinCounter=0;

    // Sig yield values
    cout << "==========================================================" << endl;
    cout << "Signal Yield" << endl;
    cout << "==========================================================" << endl;

    for( std::vector<double>::iterator it = sig.begin(  ); it != sig.end(  ); ++it ){ 
        cout << PtBin[PtBinCounter] << " < Pt =< " << PtBin[PtBinCounter+1] << ": " << *it << endl;
        PtBinCounter++;
    }

    // Fsig values
    cout << "==========================================================" << endl;
    cout << "Signal Fraction Yield" << endl;
    cout << "==========================================================" << endl;
    PtBinCounter = 0;
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
    TCanvas* cHist = new TCanvas( "Canvas","",1600,800 );
    cHist->Divide( 4,2 );
    for( int i=0; i<numPtBins; i++ ){
        cHist->cd( i+1 );
        InvMassPtBinned[i]->Draw( "E1" );
    }

    TH1D* PtSpecra = ( TH1D* )MassPt->ProjectionY( "Pt Spectrum",0,150 );
    TCanvas* c2 = new TCanvas( "Canvas2","", 800, 800 );
    c2->cd(  );
    PtSpecra->Draw( "E1" );

    TCanvas *cMass = new TCanvas( "cMass", "", 800,800 );
    cMass->SetFillColor      (0);
    cMass->SetBorderMode     (0);
    cMass->SetBorderSize     (10);
    cMass->SetLeftMargin     (0.20);
    cMass->SetRightMargin    (0.06);
    cMass->SetTopMargin      (0.11);
    cMass->SetBottomMargin   (0.15);
    cMass->SetFrameFillStyle (0);
    cMass->SetFrameLineStyle (0);
    cMass->SetFrameBorderMode(0);
    cMass->SetFrameBorderSize(10);
    cMass->SetFrameFillStyle (0);
    cMass->SetFrameLineStyle (0);
    cMass->SetFrameBorderMode(0);
    cMass->SetFrameBorderSize(10);

    TGaxis::SetMaxDigits( 3 );
    cMass->cd(  );
    gPad->SetTickx(  );
    gPad->SetTicky(  );
    InvMass->Draw( "E1 9" );
    fitFcn->Draw( "same" );
    backFcn->Draw( "same" );


    TLatex* ltx1 = new TLatex(  );
    TLatex* ltx2 = new TLatex(  );
    ltx1->SetTextSize( 0.032 );
    ltx2->SetTextSize( 0.045 );
    ltx1->SetNDC( kTRUE );
    ltx2->SetNDC( kTRUE );
    
    
    if( publish ) 
    {
        os << "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = " << " 8.16 TeV";
        ltx1->DrawLatex( 0.23, 0.82, os.str(  ).c_str(  ) );
        os.str( std::string(  ) ); 
        os << "L_{#lower[-0.25]{int}} = " << 62 << " nb^{#font[122]{\55}1}";
        ltx1->DrawLatex( 0.23, 0.74, os.str(  ).c_str(  ) );
        os.str( std::string(  ) );
        os << "#splitline{Mean: " << fixed << std::setprecision( 4 ) <<  meanTot << " GeV}{#sigma : " << rmsSigmaTot << " GeV}";
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

}

