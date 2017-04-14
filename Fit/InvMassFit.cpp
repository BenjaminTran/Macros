#include <TLatex.h>
#include <TStyle.h>
#include "TF1.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TMathText.h"
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdio>

#define piMass 0.13957018
#define lambdaMass 1.115683
#define PI 3.1416

//Convert number to const char* for TLatex
const char* DoubleToConstChar( double value )
{
    std::stringstream ss;
    ss << value;
    const char* str = ss.str(  ).c_str(  );
    return str;
}

//Define background fit functions
Double_t bgfunc1( Double_t *x, Double_t *par )
{
    Double_t xx = x[0];
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
    Double_t bgsig  = par[5]*TMath::Power( x[0] - ( piMass + lambdaMass ),par[6] );
    return sgaus1 + sgaus2 + bgsig;
}

Double_t total( Double_t *x, Double_t *par )
{
    return bgfunc1( x,par ) + sgfunc( x,&par[2] );
}

Double_t FourierHad( Double_t *x, Double_t *par )
{
    //Double_t xx1 = par[0]/(2*PI);
    Double_t xx1 = par[0];
    Double_t xx2 = 1 + 2*(par[1]*TMath::Cos( x[0] ) + par[2]*TMath::Cos( 2*x[0] )  + par[3]*TMath::Cos( 3*x[0] ) );
    return xx1*xx2;
}

Double_t bkgfuncdisplay( Double_t *x, Double_t *par )
{
    Double_t bg = par[0]*TMath::Power( x[0] - ( piMass + lambdaMass ), par[1] );
    Double_t bgsig = par[2]*TMath::Power( x[0] - ( piMass + lambdaMass ), par[3] );
    return bg + bgsig;
}

void InvMassFit(  )
{
    //TLatex
    std::ostringstream os; // stringstream for making dynamic TLatex labels
    double SNN = 5.02;
    int Lint = 35;
    int Nmin = 185;
    int Nmax = 220;
    int pTassMin = 1;
    int pTassMax = 3;
    int longRange = 2;


    //For Enabling TLatex labels
    //Bool_t publish = kTRUE;
    Bool_t publish = kFALSE; 
    
    gStyle->SetOptFit( 1111 );
    gStyle->SetErrorX( 0 ); //removes horizontal error bars

    //HIST APPEARANCE
    gStyle->SetPadLeftMargin( 0.15 );
    gStyle->SetPadBottomMargin( 0.15 );
    TGaxis::SetMaxDigits( 3 ); // Forces exponents after n number of digits

    //Font and Size
    gStyle->SetTextSize( 20 );
    gStyle->SetTextFont( 42 ); //2=times-bold-r-normal, 2=precision for TLatex to work


    //TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/XiAnalysisCorrelation.root " );
    //TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/NoPtCut/XiAnalysisCorrelationNoPtCutTotal.root " );
    TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/NoPtCutPeakAndSide/XiAnalysisCorrelationNoPtCutPeakAndSideTotal.root " );

    //Pull InvMass
    TH1F* invMass = ( TH1F* ) f->Get( "xiCorrelation/InvMassXi" );
    invMass->SetMarkerStyle( 21 );
    invMass->SetMarkerSize( 0.5 );

    //TH1F* invMass13 = ( TH1F* ) f->Get( "xiCorrelation/pT_13" );
    //TH1F* invMass12 = ( TH1F* ) f->Get( "xiCorrelation/pT_12" );


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


    TVirtualFitter::SetMaxIterations( 300000 );

    //for test function
    Double_t partest[9];

    invMass->Fit( "bgfit","L","",1.26,1.31 );
    invMass->Fit( "sigfit", "L R","",1.26,1.4 ); //1.335

    bgfit->GetParameters( &partest[0] );
    sigfit->GetParameters( &partest[2] );


    fitTot->SetParameters( partest );

    TCanvas *c1 = new TCanvas( "c1", "", 800,800 );
    c1->cd(  );
    gPad->SetTickx(  );
    gPad->SetTicky(  );
    invMass->SetTitle( "" );
    invMass->SetTitleOffset( 2, "Y" );
    invMass->GetYaxis(  )->SetTitleSize( 0.03 );
    invMass->GetYaxis(  )->CenterTitle( true );
    invMass->GetYaxis(  )->SetTitle( "Candidates/0.001 GeV" );
    invMass->GetYaxis(  )->SetRangeUser( 0, 18000 );
    invMass->SetTitleOffset( 1.5, "X" );
    invMass->GetXaxis(  )->SetTitleSize( 0.035 );
    invMass->GetXaxis(  )->CenterTitle( true );
    invMass->GetXaxis(  )->SetTitle( "#Lambda #pi^{#minus} invariant mass (GeV)" );
    invMass->SetStats( !publish );
    invMass->Draw( "E1 9" );
    //Drawing the fit function again for the sake of getting it onto the legend
    invMass->Fit( "fitTot", "L0","", 1.26, 1.4 );
    Double_t parfordraw[9];
    fitTot->GetParameters( &parfordraw[0] );
    TF1 *fitFcn = new TF1( "fitFcn", total, 1.26,1.4,9 );
    fitFcn->SetNpx( 250 );
    fitFcn->SetParameters( parfordraw );
    fitFcn->SetLineColor( kRed );
    fitFcn->Draw( "same" );
    TF1 *backFcn = new TF1( "backFcn", bkgfuncdisplay, 1.26,1.4,4 );
    backFcn->SetParameter( 0,fitTot->GetParameter( 0 ) );
    backFcn->SetParameter( 1,fitTot->GetParameter( 1 ) );
    backFcn->SetParameter( 2,fitTot->GetParameter( 7 ) );
    backFcn->SetParameter( 3,fitTot->GetParameter( 8 ) );
    backFcn->SetLineColor( kBlue );
    backFcn->Draw( "same" );

    
    TLatex* ltx1 = new TLatex(  );
    TLatex* ltx2 = new TLatex(  );
    ltx1->SetTextSize( 0.035 );
    ltx2->SetTextSize( 0.045 );
    ltx1->SetNDC( kTRUE );
    ltx2->SetNDC( kTRUE );

    double rmsSigma = TMath::Sqrt( 0.5*fitTot->GetParameter( "gaus1std" )*fitTot->GetParameter( "gaus1std" ) + 0.5*fitTot->GetParameter( "gaus2std" )*fitTot->GetParameter( "gaus2std" ) );
    double mean = fitTot->GetParameter( "gausMean" );

    if( publish ) 
    {
        os << "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = " << SNN << " TeV";
        ltx1->DrawLatex( 0.2, 0.82, os.str(  ).c_str(  ) );
        os.str( std::string(  ) ); 
        os << "L_{#lower[-0.25]{int}} = " << Lint << " nb^{-1}";
        ltx1->DrawLatex( 0.2, 0.74, os.str(  ).c_str(  ) );
        os.str( std::string(  ) );
        os << "#splitline{Mean: " << fixed << std::setprecision( 4 ) << mean << " GeV}{#sigma : " << fixed << std::setprecision( 4 ) << rmsSigma << " GeV}";
        ltx1->DrawLatex( 0.6, 0.75, os.str(  ).c_str(  ) );
        os.str( std::string(  ) ); // clears ostringstream
        os << Nmin << "  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< " << Nmax;
        ltx1->DrawLatex( 0.2, 0.67, os.str(  ).c_str(  ) );
        os.str( std::string(  ) );
        ltx2->DrawLatex( 0.2, 0.60, "#Xi^{#plus}/ #Xi^{#minus} " );
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


    //Calculating signal yield in peak region
    double totInt = -999.0;
    double bkgInt = -999.0;

    cout << "rms sigma " << rmsSigma << endl; 
    //totInt = fitTot->Integral( 1.31138364, 1.33291636 );
    totInt = fitTot->Integral( fitTot->GetParameter( "gausMean" ) - 2*rmsSigma, fitTot->GetParameter( "gausMean" ) + 2*rmsSigma );
    cout << "Total Fit integral " << totInt << endl;  

    bkgInt = backFcn->Integral( fitTot->GetParameter( "gausMean" ) - 2*rmsSigma, fitTot->GetParameter( "gausMean" ) + 2*rmsSigma );
    cout << "back fcn integral " << bkgInt << endl;

    double sigYieldFrac = -999.9;
    sigYieldFrac = ( totInt - bkgInt )/totInt;
    cout << "The signal yield fraction in the peak region is " << sigYieldFrac << endl;

    //Calculating number of cascades in peak region
    int bin1 = 70;
    int bin2 = 90;
    double XiCounter = 0;

    for( int i = bin1; i<=bin2; i++ ){
        XiCounter += invMass->GetBinContent( i );   
    }

    cout << "Number of cascades " << XiCounter << endl;
    cout << "Number of signal cascades " << XiCounter*sigYieldFrac << endl;

    

    //FITTING FOR V2
    //
    //Define divided hist
    bool Peak = true;
    //bool Peak = false;
    if( Peak ){
        TH1D *dPhiPeak = new TH1D( "dPhiPeak", "#Xi - h^{#pm} ", 31, -( 0.5 -
                    1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI  );
        TH1D *dPhiHad = new TH1D( "dPhiHad", "h^{#pm}- h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
        //Pull 2D Histograms
        TH2D *hbackgroundPeak = (TH2D*)   f->Get( "xiCorrelation/BackgroundPeak" );
        TH2D *hsignalPeak     = (TH2D*)   f->Get( "xiCorrelation/SignalPeak" );
        TH2D *hBackgroundHad = ( TH2D* ) f->Get( "xiCorrelation/BackgroundHad" );
        TH2D *hSignalHad = ( TH2D* ) f->Get( "xiCorrelation/SignalHad" );

        TH1::SetDefaultSumw2(  );
        //Project Phi
        /* For doing only one shoulder
        TH1D* hbPhiTotPeak = hbackgroundPeak->ProjectionY( "PhiBkgTotPeak", 1, 10 );
        TH1D* hsPhiTotPeak = hsignalPeak->ProjectionY( "PhiSigTotPeak", 1, 10 );
        TH1D* hbHadPhiTot = hBackgroundHad->ProjectionY( "PhiBkgHadTot", 1, 10 );
        TH1D* hsHadPhiTot = hSignalHad->ProjectionY( "PhiSigHadTot", 1, 10 );
        */


        // For projecting both shoulders
        TH1D* hbPhiTotPeak = hbackgroundPeak->ProjectionY( "PhiBkgTotPeak", 1, 10 );
        TH1D* hbPhiOthPeak = hbackgroundPeak->ProjectionY( "PhiBkgOthPeak", 23, -1 );
        TH1D* hsPhiTotPeak = hsignalPeak->ProjectionY( "PhiSigTotPeak", 1, 10 );
        TH1D* hsPhiOthPeak = hsignalPeak->ProjectionY( "PhiSigOthPeak", 23, -1 );
        TH1D* hbHadPhiTot = hBackgroundHad->ProjectionY( "PhiBkgHadTot", 1, 10 );
        TH1D* hbHadPhiOth = hBackgroundHad->ProjectionY( "PhiBkgHadOth", 23, -1 );
        TH1D* hsHadPhiTot = hSignalHad->ProjectionY( "PhiSigHadTot", 1, 10 );
        TH1D* hsHadPhiOth = hSignalHad->ProjectionY( "PhiSigHadOth", 23, -1 );

        hbPhiTotPeak->Add( hbPhiOthPeak );
        hsPhiTotPeak->Add( hsPhiOthPeak );

        hbHadPhiTot->Add( hbHadPhiOth );
        hsHadPhiTot->Add( hsHadPhiOth );

        
        //Divide
        dPhiPeak->Divide( hsPhiTotPeak, hbPhiTotPeak );
        dPhiHad->Divide( hsHadPhiTot, hbHadPhiTot );

        //Clone histograms for display without fit functions
        TH1D* dPhiFourierPeak = ( TH1D* )dPhiPeak->Clone(  );
        TH1D* dPhiHadFourier = ( TH1D* )dPhiHad->Clone(  );

        TF1 *FourierFitXi = new TF1( "FourierFitXi", FourierHad, -1.5, 5, 4 );
        FourierFitXi->SetNpx( 250 );
        FourierFitXi->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

        TCanvas *c4 = new TCanvas( "c4", "", 800,800 );
        c4->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        //dPhiFourierPeak->Fit( "FourierFitXi","","",0,PI );
        dPhiFourierPeak->Fit( "FourierFitXi");
        dPhiFourierPeak->SetStats( kFALSE );
        cout << "---------------------------------" << endl;
        cout << "Peak V2 for xi-h is " << FourierFitXi->GetParameter( 2 ) << endl;
        cout << "---------------------------------" << endl;

        TF1 *FourierFitHad = new TF1( "FourierFitHad", FourierHad, -1.5, 5, 4 );
        FourierFitHad->SetNpx( 250 );
        FourierFitHad->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

        TCanvas *c5 = new TCanvas( "c5", "", 800,800 );
        c5->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        //dPhiHadFourier->Fit( "FourierFitHad","","",0,PI );
        dPhiHadFourier->Fit( "FourierFitHad");
        dPhiHadFourier->SetStats( kFALSE );
        cout << "---------------------------------" << endl;
        cout << "The V2 for h-h is " << FourierFitHad->GetParameter( 2 ) << endl;
        cout << "---------------------------------" << endl;

        double maxBinContent = dPhiFourierPeak->GetBinContent( dPhiFourierPeak->GetMaximumBin(  ) );
        double minBinContent = dPhiFourierPeak->GetBinContent( dPhiFourierPeak->GetMinimumBin(  ) );
        double minRange = minBinContent - 0.005*minBinContent;
        double maxRange = minRange + 2*( maxBinContent - minBinContent );

        //ZYAM FITS
        TCanvas *c2 = new TCanvas( "c2", "", 800,800 );
        c2->cd(  );

        gPad->SetTickx(  );
        gPad->SetTicky(  );
        dPhiPeak->SetMarkerStyle( 21 );
        dPhiPeak->SetMarkerColor( 4 );
        dPhiPeak->SetTitleOffset( 2, "Y" );
        dPhiPeak->SetTitle( "Peak" );
        dPhiPeak->GetYaxis(  )->SetRangeUser( minRange , maxRange );
        dPhiPeak->GetYaxis(  )->SetTitleSize( 0.03 );
        dPhiPeak->GetYaxis(  )->CenterTitle( true );
        dPhiPeak->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} " );
        dPhiPeak->SetTitleOffset( 1.5, "X" );
        dPhiPeak->GetXaxis(  )->SetTitleSize( 0.035 );
        dPhiPeak->GetXaxis(  )->CenterTitle( true );
        dPhiPeak->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
        dPhiPeak->Draw( "E1" );
        dPhiPeak->Fit( "pol2","","", 0.4,2.4 );
        dPhiPeak->SetStats( !publish );

        TLatex *ltx3 = new TLatex(  );
        ltx3->SetTextSize( 0.035 );
        ltx3->SetNDC( kTRUE );
        ltx3->SetTextFont( 42 );

        if( publish )
        {
            os << "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = " << fixed << std::setprecision( 2 ) << SNN << " TeV";
            ltx3->DrawLatex( 0.2, 0.82, os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            os << "L_{#lower[-0.25]{int}} = " << Lint << " nb^{-1}";
            ltx3->DrawLatex( 0.2, 0.74, os.str( ).c_str(  ) );
            os.str( std::string(  ) );
            os << Nmin << "  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< " << Nmax;
            ltx3->DrawLatex( 0.2, 0.67, os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            os << pTassMin << " < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < " << pTassMax << " GeV";
            ltx3->DrawLatex( 0.2, 0.60, os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            os << "Long range (|#Delta#eta| > " <<  longRange << ")";
            ltx3->DrawLatex( 0.2, 0.53, os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );

            //ltx3->DrawLatex( 0.2, 0.82, "CMS pPb  #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
            //ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
            //ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            //ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
            //ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
            //ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
        }

        maxBinContent = dPhiHadFourier->GetBinContent( dPhiHadFourier->GetMaximumBin(  ) );
        minBinContent = dPhiHadFourier->GetBinContent( dPhiHadFourier->GetMinimumBin(  ) );
        minRange = minBinContent - 0.005*minBinContent;
        maxRange = minRange + 2*( maxBinContent - minBinContent );

        TCanvas *c3 = new TCanvas( "c3", "", 800,800 );
        c3->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        dPhiHad->SetMarkerStyle( 34 );
        dPhiHad->SetMarkerSize( 1.5 );
        dPhiHad->SetTitleOffset( 2, "Y" );
        dPhiHad->SetTitle( "" );
        dPhiHad->GetYaxis(  )->SetRangeUser( minRange , maxRange );
        dPhiHad->GetYaxis(  )->SetTitleSize( 0.03 );
        dPhiHad->GetYaxis(  )->CenterTitle( true );
        dPhiHad->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi}" );
        dPhiHad->SetTitleOffset( 1.5, "X" );
        dPhiHad->GetXaxis(  )->SetTitleSize( 0.035 );
        dPhiHad->GetXaxis(  )->CenterTitle( true );
        dPhiHad->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
        dPhiHad->Draw( "hist E1" );
        dPhiHad->Fit( "pol2","","", 0.4,2 );
        dPhiHad->SetStats( !publish );

        if( publish )
        {
            ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
            ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
            ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
            ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
            ltx2->DrawLatex( 0.7, 0.74, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
        }

        TF1 *dPhiFitPeak = dPhiPeak->GetFunction( "pol2" );
        TF1 *dPhiHadFit = dPhiHad->GetFunction( "pol2" );

        double dPhiFitMinPeak = dPhiFitPeak->GetMinimum(  );
        double dPhiHadFitMin = dPhiHadFit->GetMinimum(  );

        for( int i = 1; i < 32; i++ )
        {
            dPhiFourierPeak->AddBinContent( i, -dPhiFitMinPeak );
        }

        for( int i = 1; i < 32; i++ )
        {
            dPhiHadFourier->AddBinContent( i, -dPhiHadFitMin );
        }

        c4->cd(  );
        dPhiFourierPeak->SetMarkerStyle( 21 );
        dPhiFourierPeak->SetMarkerColor( 4 );
        dPhiFourierPeak->Draw( "E1" );
        dPhiFourierPeak->SetStats( kFALSE );
        dPhiFourierPeak->Fit( "FourierFitXi" );
        dPhiFourierPeak->SetTitle( "Peak" );
        dPhiFourierPeak->SetTitleOffset( 2, "Y" );
        dPhiFourierPeak->GetYaxis(  )->CenterTitle( true );
        dPhiFourierPeak->GetYaxis(  )->SetTitleSize( 0.03 );
        dPhiFourierPeak->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
        dPhiFourierPeak->GetYaxis(  )->SetRangeUser( -0.0004, 0.008 );
        dPhiFourierPeak->SetTitleOffset( 1.5, "X" );
        dPhiFourierPeak->GetXaxis(  )->SetTitleSize( 0.035 );
        dPhiFourierPeak->GetXaxis(  )->CenterTitle( true );
        dPhiFourierPeak->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );


        if( publish )
        {
            ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
            ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
            ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
            ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
            ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
        }

        c5->cd(  );
        dPhiHadFourier->SetMarkerStyle( 34 );
        dPhiHadFourier->SetMarkerSize( 1.5 );
        dPhiHadFourier->Draw( "E1" );
        dPhiHadFourier->SetStats( kFALSE );
        dPhiHadFourier->Fit( "FourierFitHad" );
        dPhiHadFourier->SetTitle( "" );
        dPhiHadFourier->SetTitleOffset( 2, "Y" );
        dPhiHadFourier->GetYaxis(  )->CenterTitle( true );
        dPhiHadFourier->GetYaxis(  )->SetTitleSize( 0.03 );
        dPhiHadFourier->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
        dPhiHadFourier->GetYaxis(  )->SetRangeUser( -0.0004 , 0.008);
        dPhiHadFourier->SetTitleOffset( 1.5, "X" );
        dPhiHadFourier->GetXaxis(  )->SetTitleSize( 0.035 );
        dPhiHadFourier->GetXaxis(  )->CenterTitle( true );
        dPhiHadFourier->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );

        if( publish )
        {
            ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
            ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
            ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
            ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
            ltx2->DrawLatex( 0.7, 0.74, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
    }
    }
    else{
        TH1D *dPhiSide = new TH1D( "dPhiSide", "#Xi - h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI  );
        TH1D *dPhiHad = new TH1D( "dPhiHad", "h^{#pm}- h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
        TH2D *hbackgroundSide = (TH2D*)   f->Get( "xiCorrelation/BackgroundSide" );
        TH2D *hsignalSide     = (TH2D*)   f->Get( "xiCorrelation/SignalSide" );
        TH2D *hBackgroundHad = ( TH2D* ) f->Get( "xiCorrelation/BackgroundHad" );
        TH2D *hSignalHad = ( TH2D* ) f->Get( "xiCorrelation/SignalHad" );
        //Project Phi
        TH1D* hbPhiTotSide = hbackgroundSide->ProjectionY( "PhiBkgTot", 0, 10 );
        TH1D* hsPhiTotSide = hsignalSide->ProjectionY( "PhiSigTot", 0, 10 );
        TH1D* hbHadPhiTot = hBackgroundHad->ProjectionY( "PhiBkgHadTot", 0, 10 );
        TH1D* hsHadPhiTot = hSignalHad->ProjectionY( "PhiSigHadTot", 0, 10 );

        TH1::SetDefaultSumw2(  );

        //Divide
        dPhiSide->Divide( hsPhiTotSide, hbPhiTotSide );
        dPhiHad->Divide( hsHadPhiTot, hbHadPhiTot );

        //Clone histograms for display without fit functions
        TH1D* dPhiFourierSide = ( TH1D* )dPhiSide->Clone(  );
        TH1D* dPhiHadFourier = ( TH1D* )dPhiHad->Clone(  );

        TF1 *FourierFitXi = new TF1( "FourierFitXi", FourierHad, -1.5, 5, 4 );
        FourierFitXi->SetNpx( 250 );
        FourierFitXi->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

        TCanvas *c4 = new TCanvas( "c4", "", 800,800 );
        c4->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        //dPhiFourierSide->Fit( "FourierFitXi","","",0,PI );
        dPhiFourierSide->Fit( "FourierFitXi");
        dPhiFourierSide->SetStats( kFALSE );
        cout << "---------------------------------" << endl;
        cout << "Side V2 for xi-h is " << FourierFitXi->GetParameter( 2 ) << endl;
        cout << "---------------------------------" << endl;

        TF1 *FourierFitHad = new TF1( "FourierFitHad", FourierHad, -1.5, 5, 4 );
        FourierFitHad->SetNpx( 250 );
        FourierFitHad->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

        TCanvas *c5 = new TCanvas( "c5", "", 800,800 );
        c5->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        //dPhiHadFourier->Fit( "FourierFitHad","","",0,PI );
        dPhiHadFourier->Fit( "FourierFitHad");
        dPhiHadFourier->SetStats( kFALSE );
        cout << "---------------------------------" << endl;
        cout << "The V2 for h-h is " << FourierFitHad->GetParameter( 2 ) << endl;
        cout << "---------------------------------" << endl;

        double maxBinContent = dPhiFourierSide->GetBinContent( dPhiFourierSide->GetMaximumBin(  ) );
        double minBinContent = dPhiFourierSide->GetBinContent( dPhiFourierSide->GetMinimumBin(  ) );
        double minRange = minBinContent - 0.005*minBinContent;
        double maxRange = minRange + 2*( maxBinContent - minBinContent );


        //ZYAM FITS
        TCanvas *c2 = new TCanvas( "c2", "", 800,800 );
        c2->cd(  );

        gPad->SetTickx(  );
        gPad->SetTicky(  );
        dPhiSide->SetMarkerStyle( 21 );
        dPhiSide->SetMarkerColor( 4 );
        dPhiSide->SetTitleOffset( 2, "Y" );
        dPhiSide->SetTitle( "Sideband" );
        dPhiSide->GetYaxis(  )->SetRangeUser( minRange , maxRange );
        dPhiSide->GetYaxis(  )->SetTitleSize( 0.03 );
        dPhiSide->GetYaxis(  )->CenterTitle( true );
        dPhiSide->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} " );
        dPhiSide->SetTitleOffset( 1.5, "X" );
        dPhiSide->GetXaxis(  )->SetTitleSize( 0.035 );
        dPhiSide->GetXaxis(  )->CenterTitle( true );
        dPhiSide->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
        dPhiSide->Draw( "E1" );
        dPhiSide->Fit( "pol2","","", 0.4,2.4 );
        dPhiSide->SetStats( !publish );

        TLatex *ltx3 = new TLatex(  );
        ltx3->SetTextSize( 0.035 );
        ltx3->SetNDC( kTRUE );
        ltx3->SetTextFont( 42 );

        if( publish )
        {
            os << "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = " << fixed << std::setprecision( 2 ) << SNN << " TeV";
            ltx3->DrawLatex( 0.2, 0.82, os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            os << "L_{#lower[-0.25]{int}} = " << Lint << " nb^{-1}";
            ltx3->DrawLatex( 0.2, 0.74, os.str( ).c_str(  ) );
            os.str( std::string(  ) );
            os << Nmin << "  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< " << Nmax;
            ltx3->DrawLatex( 0.2, 0.67, os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            os << pTassMin << " < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < " << pTassMax << " GeV";
            ltx3->DrawLatex( 0.2, 0.60, os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            os << "Long range (|#Delta#eta| > " <<  longRange << ")";
            ltx3->DrawLatex( 0.2, 0.53, os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );

            //ltx3->DrawLatex( 0.2, 0.82, "CMS pPb  #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
            //ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
            //ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            //ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
            //ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
            //ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
        }


        maxBinContent = dPhiHadFourier->GetBinContent( dPhiHadFourier->GetMaximumBin(  ) );
        minBinContent = dPhiHadFourier->GetBinContent( dPhiHadFourier->GetMinimumBin(  ) );
        minRange = minBinContent - 0.005*minBinContent;
        maxRange = minRange + 2*( maxBinContent - minBinContent );

        TCanvas *c3 = new TCanvas( "c3", "", 800,800 );
        c3->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        dPhiHad->SetMarkerStyle( 34 );
        dPhiHad->SetMarkerSize( 1.5 );
        dPhiHad->SetTitleOffset( 2, "Y" );
        dPhiHad->SetTitle( "" );
        dPhiHad->GetYaxis(  )->SetRangeUser( minRange , maxRange );
        dPhiHad->GetYaxis(  )->SetTitleSize( 0.03 );
        dPhiHad->GetYaxis(  )->CenterTitle( true );
        dPhiHad->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi}" );
        dPhiHad->SetTitleOffset( 1.5, "X" );
        dPhiHad->GetXaxis(  )->SetTitleSize( 0.035 );
        dPhiHad->GetXaxis(  )->CenterTitle( true );
        dPhiHad->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
        dPhiHad->Draw( "hist E1" );
        dPhiHad->Fit( "pol2","","", 0.4,2 );
        dPhiHad->SetStats( !publish );

        if( publish )
        {
            ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
            ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
            ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
            ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
            ltx2->DrawLatex( 0.7, 0.74, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
        }

        TF1 *dPhiFitSide = dPhiSide->GetFunction( "pol2" );
        TF1 *dPhiHadFit = dPhiHad->GetFunction( "pol2" );

        double dPhiFitMinSide = dPhiFitSide->GetMinimum(  );
        double dPhiHadFitMin = dPhiHadFit->GetMinimum(  );

        for( int i = 1; i < 32; i++ )
        {
            dPhiFourierSide->AddBinContent( i, -dPhiFitMinSide );
        }

        for( int i = 1; i < 32; i++ )
        {
            dPhiHadFourier->AddBinContent( i, -dPhiHadFitMin );
        }

        c4->cd(  );
        dPhiFourierSide->SetMarkerStyle( 21 );
        dPhiFourierSide->SetMarkerColor( 4 );
        dPhiFourierSide->Draw( "E1" );
        dPhiFourierSide->SetStats( kFALSE );
        dPhiFourierSide->Fit( "FourierFitXi" );
        dPhiFourierSide->SetTitle( "Sideband" );
        dPhiFourierSide->SetTitleOffset( 2, "Y" );
        dPhiFourierSide->GetYaxis(  )->CenterTitle( true );
        dPhiFourierSide->GetYaxis(  )->SetTitleSize( 0.03 );
        dPhiFourierSide->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
        dPhiFourierSide->GetYaxis(  )->SetRangeUser( -0.0004, 0.008 );
        dPhiFourierSide->SetTitleOffset( 1.5, "X" );
        dPhiFourierSide->GetXaxis(  )->SetTitleSize( 0.035 );
        dPhiFourierSide->GetXaxis(  )->CenterTitle( true );
        dPhiFourierSide->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );

        if( publish )
        {
            ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
            ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
            ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
            ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
            ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
        }

        c5->cd(  );
        dPhiHadFourier->SetMarkerStyle( 34 );
        dPhiHadFourier->SetMarkerSize( 1.5 );
        dPhiHadFourier->Draw( "E1" );
        dPhiHadFourier->SetStats( kFALSE );
        dPhiHadFourier->Fit( "FourierFitHad" );
        dPhiHadFourier->SetTitle( "" );
        dPhiHadFourier->SetTitleOffset( 2, "Y" );
        dPhiHadFourier->GetYaxis(  )->CenterTitle( true );
        dPhiHadFourier->GetYaxis(  )->SetTitleSize( 0.03 );
        dPhiHadFourier->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
        dPhiHadFourier->GetYaxis(  )->SetRangeUser( -0.0004 , 0.008);
        dPhiHadFourier->SetTitleOffset( 1.5, "X" );
        dPhiHadFourier->GetXaxis(  )->SetTitleSize( 0.035 );
        dPhiHadFourier->GetXaxis(  )->CenterTitle( true );
        dPhiHadFourier->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );

        if( publish )
        {
            ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
            ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
            ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
            ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
            ltx2->DrawLatex( 0.7, 0.74, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
        }
    }

    
    //Pull 2D Histograms
    //TH2D *hbackground = (TH2D*)   f->Get( "xiCorrelation/Background" );
    //TH2D *hsignal     = (TH2D*)   f->Get( "xiCorrelation/Signal" );








    
    
    


    


    //Determine Range User values as twice the difference between max and min
    //double maxBinContent = dPhiFourier->GetBinContent( dPhiFourier->GetMaximumBin(  ) );
    //double minBinContent = dPhiFourier->GetBinContent( dPhiFourier->GetMinimumBin(  ) );
    //double minRange = minBinContent - 0.005*minBinContent;
    //double maxRange = minRange + 2*( maxBinContent - minBinContent );






    


    



    /*
    //Calculate min and max for RangeUser
    maxBinContent = dPhiFourier->GetBinContent( dPhiFourier->GetMaximumBin(  ) );
    minBinContent = dPhiFourier->GetBinContent( dPhiFourier->GetMinimumBin(  ) );
    minRange = -0.0004;
    maxRange = minRange + 3*( maxBinContent - minBinContent );
    */




    /*
    maxBinContent = dPhiHadFourier->GetBinContent( dPhiHadFourier->GetMaximumBin(  ) );
    minBinContent = dPhiHadFourier->GetBinContent( dPhiHadFourier->GetMinimumBin(  ) );
    minRange = -0.0004;
    maxRange = minRange + 3*( maxBinContent - minBinContent );
    */



}
