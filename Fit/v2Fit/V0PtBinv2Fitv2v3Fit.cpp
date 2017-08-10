//Includes
#include <TLatex.h>
#include <TStyle.h>
#include "TColor.h"
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
#include "TPad.h"
#include <TString.h>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <fstream>

#define PI 3.1416

Double_t FourierHad( Double_t *x, Double_t *par )
{
    //Double_t xx1 = par[0]/(2*PI);
    Double_t xx1 = par[0];
    Double_t xx2 = 1 + 2*(par[1]*TMath::Cos( x[0] ) + par[2]*TMath::Cos( 2*x[0] )  + par[3]*TMath::Cos( 3*x[0] ) );
    return xx1*xx2;
}


void V0PtBinv2Fit(  )
{
    bool Peak = true;
	//bool Peak = false;
	//Aesthetics

    //TLatex labels
    std::ostringstream os; // stringstream for making dynamic TLatex labels
    double SNN = 8.16;
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

	//Files
	//TFiles
    TFile *f = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationJL7_8.root" );
    TFile *fhad = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Thesis/XiAnalysisCorrelationPtCut8TeVPD1_4_ForFinal.root" ); //For v2 of hadron

	//Txt files
	ofstream myfile;
	if(Peak) myfile.open("v2_v3Peak.txt");
	else myfile.open("v2_v3Sideband.txt");

    TVirtualFitter::SetMaxIterations( 300000 );
    TH1::SetDefaultSumw2(  );
    std::vector<double> PtBin_ks = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0}; //if number of bins changes make sure you change numPtBins
    std::vector<double> PtBin_la = {0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0}; //if number of bins changes make sure you change numPtBins
    int numPtBins_ks = PtBin_ks.size() - 1;
    int numPtBins_la = PtBin_la.size() - 1;
    TH1D* dPhiFourierPeak_ks[numPtBins_ks];
    TH1D* dPhiFourierSide_ks[numPtBins_ks];
    TH1D* dPhiFourierPeak_la[numPtBins_la];
    TH1D* dPhiFourierSide_la[numPtBins_la];

    TH1D* dPhiPeak_ks[numPtBins_ks];
    TH1D* dPhiSide_ks[numPtBins_ks];
    TH1D* dPhiPeak_la[numPtBins_la];
    TH1D* dPhiSide_la[numPtBins_la];

    TF1* FourierFit_ks[numPtBins_ks];
    TF1* FourierFit_la[numPtBins_la];
    std::vector<double> v2values_ks;
    std::vector<double> v2error_ks;
	std::vector<double> v3values_ks;
	std::vector<double> v3error_ks;
    double v2value_h = 0;
    double v3value_h = 0;

    std::vector<double> v2values_la;
    std::vector<double> v2error_la;
	std::vector<double> v3values_la;
	std::vector<double> v3error_la;
    double v2error_h = 0;
    double v3error_h = 0;

    TLatex* ltx2 = new TLatex(  );
    ltx2->SetTextSize( 0.045 );
    ltx2->SetNDC( kTRUE );

    //FITTING FOR V2
	//KSHORT
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
	cout << "KSHORT KSHORT KSHORT"                                                             << endl;
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
    for( int i=0; i<numPtBins_ks; i++ )
    {
        if( Peak ){
            dPhiPeak_ks[i] = new TH1D( Form( "dPhiPeak_ks%d",i ), "K_{S}^{0} - h^{#pm} ", 31, -( 0.5 -
                        1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI  );
            TH1D *dPhiHad = new TH1D( "dPhiHad", "h^{#pm}- h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
            //Pull 2D Histograms
            TH2D *hbackgroundPeak = (TH2D*) f->Get( Form( "v0Correlation/backgroundkshort_pt%d",i ) );
            TH2D *hsignalPeak     = (TH2D*) f->Get( Form( "v0Correlation/signalkshort_pt%d",i ) );
            TH2D *hBackgroundHad  = (TH2D*) fhad->Get( "xiCorrelation/BackgroundHad" );
            TH2D *hSignalHad      = (TH2D*) fhad->Get( "xiCorrelation/SignalHad" );

            TH1::SetDefaultSumw2(  );
            //Project Phi

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
            dPhiPeak_ks[i]->Divide( hsPhiTotPeak, hbPhiTotPeak );
            dPhiHad->Divide( hsHadPhiTot, hbHadPhiTot );

            //Clone histograms for display without fit functions
            dPhiFourierPeak_ks[i] = ( TH1D* )dPhiPeak_ks[i]->Clone(  );
            TH1D* dPhiHadFourier = ( TH1D* )dPhiHad->Clone(  );

            FourierFit_ks[i] = new TF1( Form( "FourierFit_ks%d",i ), FourierHad, -1.5, 5, 4 );
            FourierFit_ks[i]->SetNpx( 250 );
            FourierFit_ks[i]->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *c4_ks = new TCanvas( "c4_ks", "", 800,800 );
            c4_ks->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            //dPhiFourierPeak_ks->Fit( "FourierFit_ks","","",0,PI );
            dPhiFourierPeak_ks[i]->Fit( Form( "FourierFit_ks%d",i ) );
            dPhiFourierPeak_ks[i]->SetStats( kFALSE );
            v2values_ks.push_back( FourierFit_ks[i]->GetParameter( 2 ) );
            v2error_ks.push_back( FourierFit_ks[i]->GetParError( 2 ) );
			v3values_ks.push_back(FourierFit_ks[i]->GetParameter(3));
			v3error_ks.push_back(FourierFit_ks[i]->GetParError(3));
            cout << "---------------------------------" << endl;
            cout << "Peak V2 for xi-h is " << FourierFit_ks[i]->GetParameter( 2 ) << endl;
            cout << "---------------------------------" << endl;

            TF1 *FourierFitHad = new TF1( "FourierFitHad", FourierHad, -1.5, 5, 4 );
            FourierFitHad->SetNpx( 250 );
            FourierFitHad->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *c5_ks = new TCanvas( "c5_ks", "", 800,800 );
            c5_ks->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            //dPhiHadFourier->Fit( "FourierFitHad","","",0,PI );
            dPhiHadFourier->Fit( "FourierFitHad");
            dPhiHadFourier->SetStats( kFALSE );
            v2value_h = FourierFitHad->GetParameter(2);
            v2error_h = FourierFitHad->GetParError(2);
			v3value_h = FourierFitHad->GetParameter(3);
			v3error_h = FourierFitHad->GetParError(3);
            cout << "---------------------------------" << endl;
            cout << "The V2 for h-h is " << FourierFitHad->GetParameter( 2 ) << endl;
            cout << "---------------------------------" << endl;

            double maxBinContent = dPhiFourierPeak_ks[i]->GetBinContent( dPhiFourierPeak_ks[i]->GetMaximumBin(  ) );
            double minBinContent = dPhiFourierPeak_ks[i]->GetBinContent( dPhiFourierPeak_ks[i]->GetMinimumBin(  ) );
            double minRange = minBinContent - 0.005*minBinContent;
            double maxRange = minRange + 2*( maxBinContent - minBinContent );

            //ZYAM FITS
            TCanvas *c2_ks = new TCanvas( "c2_ks", "", 800,800 );
            c2_ks->cd(  );

            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiPeak_ks[i]->SetMarkerStyle( 21 );
            dPhiPeak_ks[i]->SetMarkerColor( 4 );
            dPhiPeak_ks[i]->SetTitleOffset( 2, "Y" );
            dPhiPeak_ks[i]->SetTitle( "Peak" );
            dPhiPeak_ks[i]->GetYaxis(  )->SetRangeUser( minRange , maxRange );
            dPhiPeak_ks[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiPeak_ks[i]->GetYaxis(  )->CenterTitle( true );
            dPhiPeak_ks[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} " );
            dPhiPeak_ks[i]->SetTitleOffset( 1.5, "X" );
            dPhiPeak_ks[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiPeak_ks[i]->GetXaxis(  )->CenterTitle( true );
            dPhiPeak_ks[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
            //dPhiPeak_ks[i]->Draw( "E1" );
            dPhiPeak_ks[i]->Fit( "pol2","","", 0.4,2.4 );
            dPhiPeak_ks[i]->SetStats( !publish );

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

            TCanvas *c3_ks = new TCanvas( "c3_ks", "", 800,800 );
            c3_ks->cd(  );
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
            //dPhiHad->Draw( "hist E1" );
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

            TF1 *dPhiFitPeak = dPhiPeak_ks[i]->GetFunction( "pol2" );
            TF1 *dPhiHadFit = dPhiHad->GetFunction( "pol2" );

            double dPhiFitMinPeak = dPhiFitPeak->GetMinimum(  );
            double dPhiHadFitMin = dPhiHadFit->GetMinimum(  );

            for( int j = 1; j < 32; j++ )
            {
                dPhiFourierPeak_ks[i]->AddBinContent( j, -dPhiFitMinPeak );
            }

            for( int j = 1; j < 32; j++ )
            {
                dPhiHadFourier->AddBinContent( j, -dPhiHadFitMin );
            }

            c4_ks->cd(  );
            dPhiFourierPeak_ks[i]->SetMarkerStyle( 21 );
            dPhiFourierPeak_ks[i]->SetMarkerColor( 4 );
            //dPhiFourierPeak_ks[i]->Draw( "E1" );
            dPhiFourierPeak_ks[i]->SetStats( kFALSE );
            dPhiFourierPeak_ks[i]->Fit( Form( "FourierFit_ks%d",i ) );
            os << "Peak " << PtBin_ks[i] << "_Pt_" << PtBin_ks[i+1];
            dPhiFourierPeak_ks[i]->SetTitle( os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            dPhiFourierPeak_ks[i]->SetTitleOffset( 2, "Y" );
            dPhiFourierPeak_ks[i]->GetYaxis(  )->CenterTitle( true );
            dPhiFourierPeak_ks[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiFourierPeak_ks[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
            dPhiFourierPeak_ks[i]->GetYaxis(  )->SetRangeUser( -0.0004, 0.008 );
            dPhiFourierPeak_ks[i]->SetTitleOffset( 1.5, "X" );
            dPhiFourierPeak_ks[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiFourierPeak_ks[i]->GetXaxis(  )->CenterTitle( true );
            dPhiFourierPeak_ks[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );


            if( publish )
            {
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }

            c5_ks->cd(  );
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
            dPhiSide_ks[i] = new TH1D( Form( "dPhiSide_ks%d",i ), "K_{S}^{0} - h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI  );
            TH1D *dPhiHad = new TH1D( "dPhiHad", "h^{#pm}- h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
            TH2D *hbackgroundSide = (TH2D*) f->Get( Form( "v0Correlation/backgroundkshort_bkg_pt%d",i ) );
            TH2D *hsignalSide     = (TH2D*) f->Get( Form( "v0Correlation/signalkshort_bkg_pt%d",i ) );
            TH2D *hBackgroundHad  = (TH2D*) fhad->Get( "xiCorrelation/BackgroundHad" );
            TH2D *hSignalHad      = (TH2D*) fhad->Get( "xiCorrelation/SignalHad" );

            TH1::SetDefaultSumw2(  );

            //Project Phi
            TH1D* hbPhiTotSide = hbackgroundSide->ProjectionY( "PhiBkgTot", 0, 10 );
            TH1D* hbPhiOthSide = hbackgroundSide->ProjectionY( "PhiBkgOthPeak", 23, -1 );
            TH1D* hsPhiTotSide = hsignalSide->ProjectionY( "PhiSigTot", 0, 10 );
            TH1D* hsPhiOthSide = hsignalSide->ProjectionY( "PhiSigOthPeak", 23, -1 );
            TH1D* hbHadPhiTot = hBackgroundHad->ProjectionY( "PhiBkgHadTot", 0, 10 );
            TH1D* hbHadPhiOth = hBackgroundHad->ProjectionY( "PhiBkgHadOth", 23, -1 );
            TH1D* hsHadPhiTot = hSignalHad->ProjectionY( "PhiSigHadTot", 0, 10 );
            TH1D* hsHadPhiOth = hSignalHad->ProjectionY( "PhiSigHadOth", 23, -1 );


            hbPhiTotSide->Add( hbPhiOthSide );
            hsPhiTotSide->Add( hsPhiOthSide );

            hbHadPhiTot->Add( hbHadPhiOth );
            hsHadPhiTot->Add( hsHadPhiOth );

            //Divide
            dPhiSide_ks[i]->Divide( hsPhiTotSide, hbPhiTotSide );
            dPhiHad->Divide( hsHadPhiTot, hbHadPhiTot );

            //Clone histograms for display without fit functions
            dPhiFourierSide_ks[i] = ( TH1D* )dPhiSide_ks[i]->Clone(  );
            TH1D* dPhiHadFourier = ( TH1D* )dPhiHad->Clone(  );

            FourierFit_ks[i] = new TF1( Form( "FourierFit_ks%d",i ) , FourierHad, -1.5, 5, 4 );
            FourierFit_ks[i]->SetNpx( 250 );
            FourierFit_ks[i]->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *c4_ks = new TCanvas( "c4_ks", "", 800,800 );
            c4_ks->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            //dPhiFourierSide_ks->Fit( "FourierFit_ks","","",0,PI );
            dPhiFourierSide_ks[i]->Fit( Form( "FourierFit_ks%d",i ) );
            dPhiFourierSide_ks[i]->SetStats( kFALSE );
            v2values_ks.push_back( FourierFit_ks[i]->GetParameter( 2 ) );
            v2error_ks.push_back( FourierFit_ks[i]->GetParError( 2 ) );
            v3values_ks.push_back( FourierFit_ks[i]->GetParameter( 3 ) );
            v3error_ks.push_back( FourierFit_ks[i]->GetParError( 3 ) );
            cout << "---------------------------------" << endl;
            cout << "Side V2 for xi-h is " << FourierFit_ks[i]->GetParameter( 2 ) << endl;
            cout << "---------------------------------" << endl;

            TF1 *FourierFitHad = new TF1( "FourierFitHad", FourierHad, -1.5, 5, 4 );
            FourierFitHad->SetNpx( 250 );
            FourierFitHad->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *c5_ks = new TCanvas( "c5_ks", "", 800,800 );
            c5_ks->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            //dPhiHadFourier->Fit( "FourierFitHad","","",0,PI );
            dPhiHadFourier->Fit( "FourierFitHad");
            dPhiHadFourier->SetStats( kFALSE );
            cout << "---------------------------------" << endl;
            cout << "The V2 for h-h is " << FourierFitHad->GetParameter( 2 ) << endl;
            cout << "---------------------------------" << endl;

            double maxBinContent = dPhiFourierSide_ks[i]->GetBinContent( dPhiFourierSide_ks[i]->GetMaximumBin(  ) );
            double minBinContent = dPhiFourierSide_ks[i]->GetBinContent( dPhiFourierSide_ks[i]->GetMinimumBin(  ) );
            double minRange = minBinContent - 0.005*minBinContent;
            double maxRange = minRange + 2*( maxBinContent - minBinContent );


            //ZYAM FITS
            TCanvas *c2_ks = new TCanvas( "c2_ks", "", 800,800 );
            c2_ks->cd(  );

            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiSide_ks[i]->SetMarkerStyle( 21 );
            dPhiSide_ks[i]->SetMarkerColor( 4 );
            dPhiSide_ks[i]->SetTitleOffset( 2, "Y" );
            dPhiSide_ks[i]->SetTitle( "Sideband" );
            dPhiSide_ks[i]->GetYaxis(  )->SetRangeUser( minRange , maxRange );
            dPhiSide_ks[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiSide_ks[i]->GetYaxis(  )->CenterTitle( true );
            dPhiSide_ks[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} " );
            dPhiSide_ks[i]->SetTitleOffset( 1.5, "X" );
            dPhiSide_ks[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiSide_ks[i]->GetXaxis(  )->CenterTitle( true );
            dPhiSide_ks[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
            //dPhiSide_ks[i]->Draw( "E1" );
            dPhiSide_ks[i]->Fit( "pol2","","", 0.4,2.4 );
            dPhiSide_ks[i]->SetStats( !publish );

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

            TCanvas *c3_ks = new TCanvas( "c3_ks", "", 800,800 );
            c3_ks->cd(  );
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

            TF1 *dPhiFitSide = dPhiSide_ks[i]->GetFunction( "pol2" );
            TF1 *dPhiHadFit = dPhiHad->GetFunction( "pol2" );

            double dPhiFitMinSide = dPhiFitSide->GetMinimum(  );
            double dPhiHadFitMin = dPhiHadFit->GetMinimum(  );

            for( int j = 1; j < 32; j++ )
            {
                dPhiFourierSide_ks[i]->AddBinContent( j, -dPhiFitMinSide );
            }

            for( int j = 1; j < 32; j++ )
            {
                dPhiHadFourier->AddBinContent( j, -dPhiHadFitMin );
            }

            c4_ks->cd(  );
            dPhiFourierSide_ks[i]->SetMarkerStyle( 21 );
            dPhiFourierSide_ks[i]->SetMarkerColor( 4 );
            //dPhiFourierSide_ks[i]->Draw( "E1" );
            dPhiFourierSide_ks[i]->SetStats( kFALSE );
            dPhiFourierSide_ks[i]->Fit( Form( "FourierFit_ks%d",i ) );
            os << "SideBand " << PtBin_ks[i] << "_Pt_" << PtBin_ks[i+1];
            dPhiFourierSide_ks[i]->SetTitle( os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            dPhiFourierSide_ks[i]->SetTitleOffset( 2, "Y" );
            dPhiFourierSide_ks[i]->GetYaxis(  )->CenterTitle( true );
            dPhiFourierSide_ks[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiFourierSide_ks[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
            dPhiFourierSide_ks[i]->GetYaxis(  )->SetRangeUser( -0.0004, 0.008 );
            dPhiFourierSide_ks[i]->SetTitleOffset( 1.5, "X" );
            dPhiFourierSide_ks[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiFourierSide_ks[i]->GetXaxis(  )->CenterTitle( true );
            dPhiFourierSide_ks[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );

            if( publish )
            {
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }

            c5_ks->cd(  );
            dPhiHadFourier->SetMarkerStyle( 34 );
            dPhiHadFourier->SetMarkerSize( 1.5 );
            dPhiHadFourier->Draw( "E1" );
            dPhiHadFourier->SetStats( kFALSE );
            dPhiHadFourier->Fit( "FourierFitHad" );
            dPhiHadFourier->SetTitle( "h^{#pm} Pairing" );
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
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.74,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.67, "0.3 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{h}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.60, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.53, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }
        }
    }

	//LAMBDA
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
	cout << "LAMBDA LAMBDA LAMBDA"                                                             << endl;
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
    for( int i=0; i<numPtBins_la; i++ )
    {
        if( Peak ){
            dPhiPeak_la[i] = new TH1D( Form( "dPhiPeak_la%d",i ), "#Lambda - h^{#pm} ", 31, -( 0.5 -
                        1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI  );
            TH1D *dPhiHad = new TH1D( "dPhiHad", "h^{#pm}- h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
            //Pull 2D Histograms
            TH2D *hbackgroundPeak = (TH2D*) f->Get( Form( "v0Correlation/backgroundlambda_pt%d",i ) );
            TH2D *hsignalPeak     = (TH2D*) f->Get( Form( "v0Correlation/signallambda_pt%d",i ) );
            TH2D *hBackgroundHad  = (TH2D*) fhad->Get( "xiCorrelation/BackgroundHad" );
            TH2D *hSignalHad      = (TH2D*) fhad->Get( "xiCorrelation/SignalHad" );

            TH1::SetDefaultSumw2(  );
            //Project Phi

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
            dPhiPeak_la[i]->Divide( hsPhiTotPeak, hbPhiTotPeak );
            dPhiHad->Divide( hsHadPhiTot, hbHadPhiTot );

            //Clone histograms for display without fit functions
            dPhiFourierPeak_la[i] = ( TH1D* )dPhiPeak_la[i]->Clone(  );
            TH1D* dPhiHadFourier = ( TH1D* )dPhiHad->Clone(  );

            FourierFit_la[i] = new TF1( Form( "FourierFit_la%d",i ), FourierHad, -1.5, 5, 4 );
            FourierFit_la[i]->SetNpx( 250 );
            FourierFit_la[i]->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *c4_la = new TCanvas( "c4_la", "", 800,800 );
            c4_la->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            //dPhiFourierPeak_la->Fit( "FourierFit_la","","",0,PI );
            dPhiFourierPeak_la[i]->Fit( Form( "FourierFit_la%d",i ) );
            dPhiFourierPeak_la[i]->SetStats( kFALSE );
            v2values_la.push_back(FourierFit_la[i]->GetParameter(2));
            v2error_la.push_back(FourierFit_la[i]->GetParError(2));
			v3values_la.push_back(FourierFit_la[i]->GetParameter(3));
			v3error_la.push_back(FourierFit_la[i]->GetParError(3));
            cout << "---------------------------------" << endl;
            cout << "Peak V2 for xi-h is " << FourierFit_la[i]->GetParameter( 2 ) << endl;
            cout << "---------------------------------" << endl;

            TF1 *FourierFitHad = new TF1( "FourierFitHad", FourierHad, -1.5, 5, 4 );
            FourierFitHad->SetNpx( 250 );
            FourierFitHad->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *c5_la = new TCanvas( "c5_la", "", 800,800 );
            c5_la->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            //dPhiHadFourier->Fit( "FourierFitHad","","",0,PI );
            dPhiHadFourier->Fit( "FourierFitHad");
            dPhiHadFourier->SetStats( kFALSE );
            cout << "---------------------------------" << endl;
            cout << "The V2 for h-h is " << FourierFitHad->GetParameter( 2 ) << endl;
            cout << "---------------------------------" << endl;

            double maxBinContent = dPhiFourierPeak_la[i]->GetBinContent( dPhiFourierPeak_la[i]->GetMaximumBin(  ) );
            double minBinContent = dPhiFourierPeak_la[i]->GetBinContent( dPhiFourierPeak_la[i]->GetMinimumBin(  ) );
            double minRange = minBinContent - 0.005*minBinContent;
            double maxRange = minRange + 2*( maxBinContent - minBinContent );

            //ZYAM FITS
            TCanvas *c2_la = new TCanvas( "c2_la", "", 800,800 );
            c2_la->cd(  );

            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiPeak_la[i]->SetMarkerStyle( 21 );
            dPhiPeak_la[i]->SetMarkerColor( 4 );
            dPhiPeak_la[i]->SetTitleOffset( 2, "Y" );
            dPhiPeak_la[i]->SetTitle( "Peak" );
            dPhiPeak_la[i]->GetYaxis(  )->SetRangeUser( minRange , maxRange );
            dPhiPeak_la[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiPeak_la[i]->GetYaxis(  )->CenterTitle( true );
            dPhiPeak_la[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} " );
            dPhiPeak_la[i]->SetTitleOffset( 1.5, "X" );
            dPhiPeak_la[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiPeak_la[i]->GetXaxis(  )->CenterTitle( true );
            dPhiPeak_la[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
            //dPhiPeak_la[i]->Draw( "E1" );
            dPhiPeak_la[i]->Fit( "pol2","","", 0.4,2.4 );
            dPhiPeak_la[i]->SetStats( !publish );

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

            TCanvas *c3_la = new TCanvas( "c3_la", "", 800,800 );
            c3_la->cd(  );
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
            //dPhiHad->Draw( "hist E1" );
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

            TF1 *dPhiFitPeak = dPhiPeak_la[i]->GetFunction( "pol2" );
            TF1 *dPhiHadFit = dPhiHad->GetFunction( "pol2" );

            double dPhiFitMinPeak = dPhiFitPeak->GetMinimum(  );
            double dPhiHadFitMin = dPhiHadFit->GetMinimum(  );

            for( int j = 1; j < 32; j++ )
            {
                dPhiFourierPeak_la[i]->AddBinContent( j, -dPhiFitMinPeak );
            }

            for( int j = 1; j < 32; j++ )
            {
                dPhiHadFourier->AddBinContent( j, -dPhiHadFitMin );
            }

            c4_la->cd(  );
            dPhiFourierPeak_la[i]->SetMarkerStyle( 21 );
            dPhiFourierPeak_la[i]->SetMarkerColor( 4 );
            //dPhiFourierPeak_la[i]->Draw( "E1" );
            dPhiFourierPeak_la[i]->SetStats( kFALSE );
            dPhiFourierPeak_la[i]->Fit( Form( "FourierFit_la%d",i ) );
            os << "Peak " << PtBin_la[i] << "_Pt_" << PtBin_la[i+1];
            dPhiFourierPeak_la[i]->SetTitle( os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            dPhiFourierPeak_la[i]->SetTitleOffset( 2, "Y" );
            dPhiFourierPeak_la[i]->GetYaxis(  )->CenterTitle( true );
            dPhiFourierPeak_la[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiFourierPeak_la[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
            dPhiFourierPeak_la[i]->GetYaxis(  )->SetRangeUser( -0.0004, 0.008 );
            dPhiFourierPeak_la[i]->SetTitleOffset( 1.5, "X" );
            dPhiFourierPeak_la[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiFourierPeak_la[i]->GetXaxis(  )->CenterTitle( true );
            dPhiFourierPeak_la[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );


            if( publish )
            {
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }

            c5_la->cd(  );
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
            dPhiSide_la[i] = new TH1D( Form( "dPhiSide_la%d",i ), "#Xi - h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI  );
            TH1D *dPhiHad = new TH1D( "dPhiHad", "h^{#pm}- h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
            TH2D *hbackgroundSide = (TH2D*) f->Get( Form( "v0Correlation/backgroundlambda_bkg_pt%d",i ) );
            TH2D *hsignalSide     = (TH2D*) f->Get( Form( "v0Correlation/signallambda_bkg_pt%d",i ) );
            TH2D *hBackgroundHad  = (TH2D*) fhad->Get( "xiCorrelation/BackgroundHad" );
            TH2D *hSignalHad      = (TH2D*) fhad->Get( "xiCorrelation/SignalHad" );

            TH1::SetDefaultSumw2(  );

            //Project Phi
            TH1D* hbPhiTotSide = hbackgroundSide->ProjectionY( "PhiBkgTot", 0, 10 );
            TH1D* hbPhiOthSide = hbackgroundSide->ProjectionY( "PhiBkgOthPeak", 23, -1 );
            TH1D* hsPhiTotSide = hsignalSide->ProjectionY( "PhiSigTot", 0, 10 );
            TH1D* hsPhiOthSide = hsignalSide->ProjectionY( "PhiSigOthPeak", 23, -1 );
            TH1D* hbHadPhiTot = hBackgroundHad->ProjectionY( "PhiBkgHadTot", 0, 10 );
            TH1D* hbHadPhiOth = hBackgroundHad->ProjectionY( "PhiBkgHadOth", 23, -1 );
            TH1D* hsHadPhiTot = hSignalHad->ProjectionY( "PhiSigHadTot", 0, 10 );
            TH1D* hsHadPhiOth = hSignalHad->ProjectionY( "PhiSigHadOth", 23, -1 );


            hbPhiTotSide->Add( hbPhiOthSide );
            hsPhiTotSide->Add( hsPhiOthSide );

            hbHadPhiTot->Add( hbHadPhiOth );
            hsHadPhiTot->Add( hsHadPhiOth );

            //Divide
            dPhiSide_la[i]->Divide( hsPhiTotSide, hbPhiTotSide );
            dPhiHad->Divide( hsHadPhiTot, hbHadPhiTot );

            //Clone histograms for display without fit functions
            dPhiFourierSide_la[i] = ( TH1D* )dPhiSide_la[i]->Clone(  );
            TH1D* dPhiHadFourier = ( TH1D* )dPhiHad->Clone(  );

            FourierFit_la[i] = new TF1( Form( "FourierFit_la%d",i ) , FourierHad, -1.5, 5, 4 );
            FourierFit_la[i]->SetNpx( 250 );
            FourierFit_la[i]->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *c4_la = new TCanvas( "c4_la", "", 800,800 );
            c4_la->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            //dPhiFourierSide_la->Fit( "FourierFit_la","","",0,PI );
            dPhiFourierSide_la[i]->Fit( Form( "FourierFit_la%d",i ) );
            dPhiFourierSide_la[i]->SetStats( kFALSE );
            v2values_la.push_back( FourierFit_la[i]->GetParameter( 2 ) );
            v2error_la.push_back( FourierFit_la[i]->GetParError( 2 ) );
            v3values_la.push_back( FourierFit_la[i]->GetParameter( 3 ) );
            v3error_la.push_back( FourierFit_la[i]->GetParError( 3 ) );
            cout << "---------------------------------" << endl;
            cout << "Side V2 for xi-h is " << FourierFit_la[i]->GetParameter( 2 ) << endl;
            cout << "---------------------------------" << endl;

            TF1 *FourierFitHad = new TF1( "FourierFitHad", FourierHad, -1.5, 5, 4 );
            FourierFitHad->SetNpx( 250 );
            FourierFitHad->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *c5_la = new TCanvas( "c5_la", "", 800,800 );
            c5_la->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            //dPhiHadFourier->Fit( "FourierFitHad","","",0,PI );
            dPhiHadFourier->Fit( "FourierFitHad");
            dPhiHadFourier->SetStats( kFALSE );
            cout << "---------------------------------" << endl;
            cout << "The V2 for h-h is " << FourierFitHad->GetParameter( 2 ) << endl;
            cout << "---------------------------------" << endl;

            double maxBinContent = dPhiFourierSide_la[i]->GetBinContent( dPhiFourierSide_la[i]->GetMaximumBin(  ) );
            double minBinContent = dPhiFourierSide_la[i]->GetBinContent( dPhiFourierSide_la[i]->GetMinimumBin(  ) );
            double minRange = minBinContent - 0.005*minBinContent;
            double maxRange = minRange + 2*( maxBinContent - minBinContent );


            //ZYAM FITS
            TCanvas *c2_la = new TCanvas( "c2_la", "", 800,800 );
            c2_la->cd(  );

            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiSide_la[i]->SetMarkerStyle( 21 );
            dPhiSide_la[i]->SetMarkerColor( 4 );
            dPhiSide_la[i]->SetTitleOffset( 2, "Y" );
            dPhiSide_la[i]->SetTitle( "Sideband" );
            dPhiSide_la[i]->GetYaxis(  )->SetRangeUser( minRange , maxRange );
            dPhiSide_la[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiSide_la[i]->GetYaxis(  )->CenterTitle( true );
            dPhiSide_la[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} " );
            dPhiSide_la[i]->SetTitleOffset( 1.5, "X" );
            dPhiSide_la[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiSide_la[i]->GetXaxis(  )->CenterTitle( true );
            dPhiSide_la[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
            //dPhiSide_la[i]->Draw( "E1" );
            dPhiSide_la[i]->Fit( "pol2","","", 0.4,2.4 );
            dPhiSide_la[i]->SetStats( !publish );

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

            TCanvas *c3_la = new TCanvas( "c3_la", "", 800,800 );
            c3_la->cd(  );
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

            TF1 *dPhiFitSide = dPhiSide_la[i]->GetFunction( "pol2" );
            TF1 *dPhiHadFit = dPhiHad->GetFunction( "pol2" );

            double dPhiFitMinSide = dPhiFitSide->GetMinimum(  );
            double dPhiHadFitMin = dPhiHadFit->GetMinimum(  );

            for( int j = 1; j < 32; j++ )
            {
                dPhiFourierSide_la[i]->AddBinContent( j, -dPhiFitMinSide );
            }

            for( int j = 1; j < 32; j++ )
            {
                dPhiHadFourier->AddBinContent( j, -dPhiHadFitMin );
            }

            c4_la->cd(  );
            dPhiFourierSide_la[i]->SetMarkerStyle( 21 );
            dPhiFourierSide_la[i]->SetMarkerColor( 4 );
            //dPhiFourierSide_la[i]->Draw( "E1" );
            dPhiFourierSide_la[i]->SetStats( kFALSE );
            dPhiFourierSide_la[i]->Fit( Form( "FourierFit_la%d",i ) );
            os << "SideBand " << PtBin_la[i] << "_Pt_" << PtBin_la[i+1];
            dPhiFourierSide_la[i]->SetTitle( os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            dPhiFourierSide_la[i]->SetTitleOffset( 2, "Y" );
            dPhiFourierSide_la[i]->GetYaxis(  )->CenterTitle( true );
            dPhiFourierSide_la[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiFourierSide_la[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
            dPhiFourierSide_la[i]->GetYaxis(  )->SetRangeUser( -0.0004, 0.008 );
            dPhiFourierSide_la[i]->SetTitleOffset( 1.5, "X" );
            dPhiFourierSide_la[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiFourierSide_la[i]->GetXaxis(  )->CenterTitle( true );
            dPhiFourierSide_la[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );

            if( publish )
            {
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }

            c5_la->cd(  );
            dPhiHadFourier->SetMarkerStyle( 34 );
            dPhiHadFourier->SetMarkerSize( 1.5 );
            dPhiHadFourier->Draw( "E1" );
            dPhiHadFourier->SetStats( kFALSE );
            dPhiHadFourier->Fit( "FourierFitHad" );
            dPhiHadFourier->SetTitle( "h^{#pm} Pairing" );
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
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.74,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.67, "0.3 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{h}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.60, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.53, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }
        }
    }

    // V2 values Kshort
    cout << "==========================================================" << endl;
    cout << "Kshort V2 values" << endl;
    cout << "==========================================================" << endl;

	myfile << "Kshort v2 values\n";

    int PtBinCounter=0;

    for( std::vector<double>::iterator it = v2values_ks.begin(  ); it != v2values_ks.end(  ); ++it ){
        cout << PtBin_ks[PtBinCounter] << " < Pt =< " << PtBin_ks[PtBinCounter + 1] << ": " << *it << endl;
		myfile << *it << "\n";
        PtBinCounter++;
    }

    PtBinCounter=0;

    // V2 errors Kshort
    cout << "==========================================================" << endl;
    cout << "Kshort V2 errors" << endl;
    cout << "==========================================================" << endl;

	myfile << "Kshort v2 errors\n";

    for( std::vector<double>::iterator it = v2error_ks.begin(  ); it != v2error_ks.end(  ); ++it ){
        cout << PtBin_ks[PtBinCounter] << " < Pt =< " << PtBin_ks[PtBinCounter + 1] << ": " << *it << endl;
		myfile << *it << "\n";
        PtBinCounter++;
    }

    PtBinCounter=0;

	// V3 values Kshort
    cout << "==========================================================" << endl;
    cout << "Kshort V3 values" << endl;
    cout << "==========================================================" << endl;

	myfile << "Kshort v3 values\n";

    for( std::vector<double>::iterator it = v3values_ks.begin(  ); it != v3values_ks.end(  ); ++it ){
        cout << PtBin_ks[PtBinCounter] << " < Pt =< " << PtBin_ks[PtBinCounter + 1] << ": " << *it << endl;
		myfile << *it << "\n";
        PtBinCounter++;
    }

    PtBinCounter=0;

    // V3 errors Kshort
    cout << "==========================================================" << endl;
    cout << "Kshort V3 errors" << endl;
    cout << "==========================================================" << endl;

	myfile << "Kshort v3 errors\n";

    for( std::vector<double>::iterator it = v3error_ks.begin(  ); it != v3error_ks.end(  ); ++it ){
        cout << PtBin_ks[PtBinCounter] << " < Pt =< " << PtBin_ks[PtBinCounter + 1] << ": " << *it << endl;
		myfile << *it << "\n";
        PtBinCounter++;
    }

    PtBinCounter=0;

	// V2 values Lambda
    cout << "==========================================================" << endl;
    cout << "Lambda V2 values" << endl;
    cout << "==========================================================" << endl;

	myfile << "Lambda v2 values\n";

    for( std::vector<double>::iterator it = v2values_la.begin(  ); it != v2values_la.end(  ); ++it ){
        cout << PtBin_la[PtBinCounter] << " < Pt =< " << PtBin_la[PtBinCounter + 1] << ": " << *it << endl;
		myfile << *it << "\n";
        PtBinCounter++;
    }

    PtBinCounter=0;

    // V2 errors Lambda
    cout << "==========================================================" << endl;
    cout << "Lambda V2 errors" << endl;
    cout << "==========================================================" << endl;

	myfile << "Lambda v2 errors\n";

    for( std::vector<double>::iterator it = v2error_la.begin(  ); it != v2error_la.end(  ); ++it ){
        cout << PtBin_la[PtBinCounter] << " < Pt =< " << PtBin_la[PtBinCounter + 1] << ": " << *it << endl;
		myfile << *it << "\n";
        PtBinCounter++;
    }

    PtBinCounter=0;

	// V3 values Lambda
    cout << "==========================================================" << endl;
    cout << "Lambda V3 values" << endl;
    cout << "==========================================================" << endl;

	myfile << "Lambda v3 values\n";

    for( std::vector<double>::iterator it = v3values_la.begin(  ); it != v3values_la.end(  ); ++it ){
        cout << PtBin_la[PtBinCounter] << " < Pt =< " << PtBin_la[PtBinCounter + 1] << ": " << *it << endl;
		myfile << *it << "\n";
        PtBinCounter++;
    }

    PtBinCounter=0;

    // V3 errors Lambda
    cout << "==========================================================" << endl;
    cout << "Lambda V3 errors" << endl;
    cout << "==========================================================" << endl;

	myfile << "Lambda v3 errors\n";

    for( std::vector<double>::iterator it = v3error_la.begin(  ); it != v3error_la.end(  ); ++it ){
        cout << PtBin_la[PtBinCounter] << " < Pt =< " << PtBin_la[PtBinCounter + 1] << ": " << *it << endl;
		myfile << *it << "\n";
        PtBinCounter++;
    }

    PtBinCounter=0;

	// V2 values Hadron
    cout << "==========================================================" << endl;
    cout << "Hadron V2 value" << endl;
    cout << "==========================================================" << endl;

    cout << v2value_h << endl;
	myfile << "Hadron v2 value\n";
    myfile << v2value_h << "\n";

    // V2 errors Hadron
    cout << "==========================================================" << endl;
    cout << "Hadron V2 error" << endl;
    cout << "==========================================================" << endl;

    cout << v2error_h << endl;
	myfile << "Hadron v2 errors\n";
    myfile << v2error_h << "\n";

	// V3 values Hadron
    cout << "==========================================================" << endl;
    cout << "Hadron V3 value" << endl;
    cout << "==========================================================" << endl;

    cout << v3value_h << endl;
	myfile << "Hadron v3 value\n";
    myfile << v3value_h << "\n";

    // V3 errors Hadron
    cout << "==========================================================" << endl;
    cout << "Hadron V3 error" << endl;
    cout << "==========================================================" << endl;

    cout << v3error_h << endl;
	myfile << "Hadron v3 errors\n";
    myfile << v3error_h << "\n";


    TCanvas* Fourier_ks = new TCanvas( "Fourier_ks", "Fourier_ks", 1600,800 );
    Fourier_ks->Divide( 5,2 );
    if( Peak ){
        for( int i=0; i<numPtBins_ks; i++ ){
            Fourier_ks->cd( i+1 );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiFourierPeak_ks[i]->Draw( "E1" );
        }
    }
    else{
        for( int i=0; i<numPtBins_ks; i++ ){
            Fourier_ks->cd( i+1 );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiFourierSide_ks[i]->Draw( "E1" );
        }
    }

				cout << "HI" << endl;
    TCanvas* Fourier_la = new TCanvas( "Fourier_la", "Fourier_la", 1600,800 );
    Fourier_la->Divide( 5,2 );
    if( Peak ){
        for( int i=0; i<numPtBins_la; i++ ){
            Fourier_la->cd( i+1 );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiFourierPeak_la[i]->Draw( "E1" );
        }
    }
    else{
        for( int i=0; i<numPtBins_la; i++ ){
            Fourier_la->cd( i+1 );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiFourierSide_la[i]->Draw( "E1" );
        }
    }

    //Output Publication plots
	if( publish )
	{
    //1D correlation functions
	//Kshort
    if( Peak ){
        TCanvas* PubFourier_ks = new TCanvas( "PubFourier_ks", "Pub_ks", 800,800 );
        PubFourier_ks->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        TH1D* dPhiFourierPeakCopy_ks = ( TH1D* )dPhiFourierPeak_ks[4]->Clone(  );
        dPhiFourierPeakCopy_ks->SetTitle( "Peak" );
        dPhiFourierPeakCopy_ks->Draw( "E1" );

        TLatex *ltx3 = new TLatex(  );
        ltx3->SetTextSize( 0.035 );
        ltx3->SetNDC( kTRUE );
        ltx3->SetTextFont( 42 );
        if( publish )
        {
            ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}" );
            ltx3->DrawLatex( 0.2, 0.74,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            ltx3->DrawLatex( 0.2, 0.67, "0.3 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
            ltx3->DrawLatex( 0.2, 0.60, "2.8 < p_{T}^{#Xi} < 3.6 GeV" );
            ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
            ltx3->SetTextSize( 0.045 );
            ltx3->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
        }
    }
    else{
        TCanvas* PubFourier_ks = new TCanvas( "PubFourier_ks", "Pub_ks", 800,800 );
        PubFourier_ks->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        TH1D* dPhiFourierSideCopy_ks = ( TH1D* )dPhiFourierSide_ks[4]->Clone(  );
        dPhiFourierSideCopy_ks->SetTitle( "SideBand" );
        dPhiFourierSideCopy_ks->Draw( "E1" );


        TLatex *ltx3 = new TLatex(  );
        ltx3->SetTextSize( 0.035 );
        ltx3->SetNDC( kTRUE );
        ltx3->SetTextFont( 42 );
        if( publish )
        {
            ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}" );
            ltx3->DrawLatex( 0.2, 0.74,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            ltx3->DrawLatex( 0.2, 0.67, "0.3 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
            ltx3->DrawLatex( 0.2, 0.60, "2.8 < p_{T}^{#Xi} < 3.6 GeV" );
            ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
            ltx3->SetTextSize( 0.045 );
            ltx3->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
        }
    }


	// Lambda
    if( Peak ){
        TCanvas* PubFourier_la = new TCanvas( "PubFourier_la", "Pub_la", 800,800 );
        PubFourier_la->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        TH1D* dPhiFourierPeakCopy_la = ( TH1D* )dPhiFourierPeak_la[4]->Clone(  );
        dPhiFourierPeakCopy_la->SetTitle( "Peak" );
        dPhiFourierPeakCopy_la->Draw( "E1" );

        TLatex *ltx3 = new TLatex(  );
        ltx3->SetTextSize( 0.035 );
        ltx3->SetNDC( kTRUE );
        ltx3->SetTextFont( 42 );
        if( publish )
        {
            ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}" );
            ltx3->DrawLatex( 0.2, 0.74,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            ltx3->DrawLatex( 0.2, 0.67, "0.3 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
            ltx3->DrawLatex( 0.2, 0.60, "2.8 < p_{T}^{#Xi} < 3.6 GeV" );
            ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
            ltx3->SetTextSize( 0.045 );
            ltx3->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
        }
    }
    else{
        TCanvas* PubFourier_la = new TCanvas( "PubFourier_la", "Pub_la", 800,800 );
        PubFourier_la->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        TH1D* dPhiFourierSideCopy_la = ( TH1D* )dPhiFourierSide_la[4]->Clone(  );
        dPhiFourierSideCopy_la->SetTitle( "SideBand" );
        dPhiFourierSideCopy_la->Draw( "E1" );


        TLatex *ltx3 = new TLatex(  );
        ltx3->SetTextSize( 0.035 );
        ltx3->SetNDC( kTRUE );
        ltx3->SetTextFont( 42 );
        if( publish )
        {
            ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}" );
            ltx3->DrawLatex( 0.2, 0.74,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            ltx3->DrawLatex( 0.2, 0.67, "0.3 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
            ltx3->DrawLatex( 0.2, 0.60, "2.8 < p_{T}^{#Xi} < 3.6 GeV" );
            ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
            ltx3->SetTextSize( 0.045 );
            ltx3->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
        }
    }
	}

    //2D Correlation function 1-3 GeV associated
    TLatex *ltx0 = new TLatex(  );
    ltx0->SetTextSize( 0.031 );
    ltx0->SetNDC( kTRUE );
    ltx0->SetTextFont( 42 );

	//Kshort
    TCanvas* TwoDCorrelation_ks = new TCanvas( "TwoDCorrelation_ks", "", 1000, 1000 );
    TwoDCorrelation_ks->SetLeftMargin( 0.2 );

    TH2D* Signal_ks = ( TH2D* )f->Get( "v0Correlation/signalkshort_pt2" );
    TH2D* Background_ks = ( TH2D* )f->Get( "v0Correlation/backgroundkshort_pt2" );

    TGaxis::SetMaxDigits( 1 );

    TH2D* Correlation_ks = ( TH2D* )Signal_ks->Clone(  );
    Correlation_ks->Divide( Background_ks );
    Correlation_ks->GetXaxis(  )->SetRangeUser( -4.0,4.0 );
    Correlation_ks->GetYaxis(  )->SetRangeUser( -PI/2.0,4.5 );
    Correlation_ks->GetXaxis(  )->SetTitle( "#Delta#eta" );
    Correlation_ks->GetXaxis(  )->SetTitleOffset( 1.4 );
    Correlation_ks->GetXaxis(  )->CenterTitle( true );
    Correlation_ks->GetYaxis(  )->SetTitle( "#Delta#phi (radians)" );
    Correlation_ks->GetYaxis(  )->SetTitleOffset( 1.4 );
    Correlation_ks->GetYaxis(  )->CenterTitle( true );
    Correlation_ks->GetZaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{d^{2}N^{pair}}{d#Delta#eta d#Delta#phi}" );
    Correlation_ks->GetZaxis(  )->SetTitleOffset( 2.3 );
    Correlation_ks->GetZaxis(  )->CenterTitle( true );
    Correlation_ks->GetXaxis(  )->SetNdivisions( 405 );
    Correlation_ks->GetYaxis(  )->SetNdivisions( 405 );
    Correlation_ks->GetZaxis(  )->SetNdivisions( 4 );
    Correlation_ks->SetTitle( "" );
    Correlation_ks->SetStats( kFALSE );
    Correlation_ks->Scale( 10 );

    const Int_t NRGBs = 5;
    const Int_t NCont = 20;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    //gStyle->SetNumberContours( 90 );

    TwoDCorrelation_ks->cd(  );
    Correlation_ks->Draw( "SURF1 FB " );


    ltx0->DrawLatex( 0.05, 0.95, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}" );
    ltx0->DrawLatex( 0.05, 0.88, "185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
    ltx0->DrawLatex( 0.05, 0.81, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
    ltx0->DrawLatex( 0.05, 0.75, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{trig}}} < 3 GeV" );
    ltx0->SetTextSize( 0.04 );
    ltx0->DrawLatex( 0.85, 0.88, "K_{S}^{0}#kern[-0.3]{#lower[0.02]{{}^{#pm}}}- h^{#pm}" );


	//Lambda
    TCanvas* TwoDCorrelation_la = new TCanvas( "TwoDCorrelation_la", "", 1000, 1000 );
    TwoDCorrelation_la->SetLeftMargin( 0.2 );

    TH2D* Signal_la = ( TH2D* )f->Get( "v0Correlation/signallambda_pt2" );
    TH2D* Background_la = ( TH2D* )f->Get( "v0Correlation/backgroundlambda_pt2" );

    TGaxis::SetMaxDigits( 1 );

    TH2D* Correlation_la = ( TH2D* )Signal_la->Clone(  );
    Correlation_la->Divide( Background_la );
    Correlation_la->GetXaxis(  )->SetRangeUser( -4.0,4.0 );
    Correlation_la->GetYaxis(  )->SetRangeUser( -PI/2.0,4.5 );
    Correlation_la->GetXaxis(  )->SetTitle( "#Delta#eta" );
    Correlation_la->GetXaxis(  )->SetTitleOffset( 1.4 );
    Correlation_la->GetXaxis(  )->CenterTitle( true );
    Correlation_la->GetYaxis(  )->SetTitle( "#Delta#phi (radians)" );
    Correlation_la->GetYaxis(  )->SetTitleOffset( 1.4 );
    Correlation_la->GetYaxis(  )->CenterTitle( true );
    Correlation_la->GetZaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{d^{2}N^{pair}}{d#Delta#eta d#Delta#phi}" );
    Correlation_la->GetZaxis(  )->SetTitleOffset( 2.3 );
    Correlation_la->GetZaxis(  )->CenterTitle( true );
    Correlation_la->GetXaxis(  )->SetNdivisions( 405 );
    Correlation_la->GetYaxis(  )->SetNdivisions( 405 );
    Correlation_la->GetZaxis(  )->SetNdivisions( 4 );
    Correlation_la->SetTitle( "" );
    Correlation_la->SetStats( kFALSE );
    Correlation_la->Scale( 10 );

    TwoDCorrelation_la->cd(  );
    Correlation_la->Draw( "SURF1 FB " );


    ltx0->DrawLatex( 0.05, 0.95, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}" );
    ltx0->DrawLatex( 0.05, 0.88, "185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
    ltx0->DrawLatex( 0.05, 0.81, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
    ltx0->DrawLatex( 0.05, 0.75, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{trig}}} < 3 GeV" );
    ltx0->SetTextSize( 0.04 );
    ltx0->DrawLatex( 0.85, 0.88, "#Lambda#kern[-0.3]{#lower[0.02]{{}^{#pm}}}- h^{#pm}" );


	/*
    TCanvas* SigAndBkg = new TCanvas( "SigAndBkg", "", 1600, 800 );
    SigAndBkg->Divide( 2,1 );
    SigAndBkg->SetLeftMargin( 0.2 );

    SigAndBkg->cd( 1 );
    Signal->GetXaxis(  )->SetRangeUser( -4.0,4.0 );
    Signal->GetYaxis(  )->SetRangeUser( -PI/2.0,4.5 );
    Signal->GetXaxis(  )->SetTitle( "#Delta#eta" );
    Signal->GetXaxis(  )->SetTitleOffset( 1.4 );
    Signal->GetXaxis(  )->CenterTitle( true );
    Signal->GetYaxis(  )->SetTitle( "#Delta#phi (radians)" );
    Signal->GetYaxis(  )->SetTitleOffset( 1.4 );
    Signal->GetYaxis(  )->CenterTitle( true );
    Signal->GetZaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{d^{2}N^{pair}}{d#Delta#eta d#Delta#phi}" );
    Signal->GetZaxis(  )->SetTitleOffset( 2.5 );
    Signal->GetZaxis(  )->CenterTitle( true );
    Signal->GetXaxis(  )->SetNdivisions( 405 );
    Signal->GetYaxis(  )->SetNdivisions( 405 );
    Signal->GetZaxis(  )->SetNdivisions( 4 );
    Signal->SetTitle( "" );
    Signal->SetStats( kFALSE );
    Signal->Draw( "Surf1 FB" );

    SigAndBkg->cd( 2 );
    Background->GetXaxis(  )->SetRangeUser( -4.0,4.0 );
    Background->GetYaxis(  )->SetRangeUser( -PI/2.0,4.5 );
    Background->GetXaxis(  )->SetTitle( "#Delta#eta" );
    Background->GetXaxis(  )->SetTitleOffset( 1.4 );
    Background->GetXaxis(  )->CenterTitle( true );
    Background->GetYaxis(  )->SetTitle( "#Delta#phi (radians)" );
    Background->GetYaxis(  )->SetTitleOffset( 1.4 );
    Background->GetYaxis(  )->CenterTitle( true );
    Background->GetZaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{d^{2}N^{pair}}{d#Delta#eta d#Delta#phi}" );
    Background->GetZaxis(  )->SetTitleOffset( 2.5 );
    Background->GetZaxis(  )->CenterTitle( true );
    Background->GetXaxis(  )->SetNdivisions( 405 );
    Background->GetYaxis(  )->SetNdivisions( 405 );
    Background->GetZaxis(  )->SetNdivisions( 4 );
    Background->SetTitle( "" );
    Background->SetStats( kFALSE );
    Background->Draw( "Surf1 FB" );
	*/


}
