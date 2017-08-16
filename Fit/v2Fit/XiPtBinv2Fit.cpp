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
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>

#define PI 3.1416

Double_t FourierHad( Double_t *x, Double_t *par )
{
    //Double_t xx1 = par[0]/(2*PI);
    Double_t xx1 = par[0];
    Double_t xx2 = 1 + 2*(par[1]*TMath::Cos( x[0] ) + par[2]*TMath::Cos( 2*x[0] )  + par[3]*TMath::Cos( 3*x[0] ) );
    return xx1*xx2;
}

void OutputVnValues(int degree, std::string region, std::string V0IDname, std::vector<double> vnvalues, std::vector<double> ptbin, std::string file)
{
    std::ostringstream output;
    std::ofstream myfile;
    myfile.open(file.c_str(),std::ios_base::app);

    output << V0IDname << " V" << degree << " values " << region << "\n";
    cout << "==========================================================" << endl;
    cout << output.str();
    cout << "==========================================================" << endl;

    myfile << output.str();

    int PtBinCounter=0;
    for(std::vector<double>::iterator it = vnvalues.begin(); it != vnvalues.end(); ++it)
    {
        cout << ptbin[PtBinCounter] << " < Pt =< " << ptbin[PtBinCounter + 1] << ": " << *it << endl;
        myfile << *it << "\n";
        PtBinCounter++;
    }
    //myfile.close();
}

void OutputVnErrors(int degree, std::string region, std::string V0IDname, std::vector<double> vnerrors, std::vector<double> ptbin, std::string file)
{
    std::ostringstream output;
    std::ofstream myfile;
    myfile.open(file.c_str(),std::ios_base::app);

    output << V0IDname << " V" << degree << " errors " << region << "\n";
    cout << "==========================================================" << endl;
    cout << output.str();
    cout << "==========================================================" << endl;

    myfile << output.str();

    int PtBinCounter=0;
    for(std::vector<double>::iterator it = vnerrors.begin(); it != vnerrors.end(); ++it)
    {
        cout << ptbin[PtBinCounter] << " < Pt =< " << ptbin[PtBinCounter + 1] << ": " << *it << endl;
        myfile << *it << "\n";
        PtBinCounter++;
    }
}

void vnCalculate(int degree, std::string V0IDname, std::vector<double> vnvalues_peak, std::vector<double> vnerrors_peak, std::vector<double> vnvalues_side, std::vector<double> vnerrors_side, std::vector<double> vnvalues_h, std::vector<double> vnerrors_h, std::vector<double> fsig, std::string file)
{
    std::ostringstream output;
    std::ofstream myfile;
    myfile.open(file.c_str(),std::ios_base::app);

    std::vector<double> sigvalues;
    std::vector<double> sigerrors;

    output << V0IDname << " signal v" << degree << " values\n";
    cout << output.str();

    for(unsigned i=0; i<fsig.size(); i++)
    {
        cout << "Pt Bin " << i+1 << endl;
        double vnObs = vnvalues_peak[i]/TMath::Sqrt(vnvalues_h[degree]);
        double vnBkg = vnvalues_side[i]/TMath::Sqrt(vnvalues_h[degree]);
        double sig = (vnObs - (1 - fsig[i])*vnBkg)/fsig[i];

        double vnObsError = vnObs*TMath::Sqrt(TMath::Power(vnerrors_peak[i]/vnvalues_peak[i],2) + TMath::Power(0.5*vnerrors_h[degree]/vnvalues_h[degree],2));
        double vnBkgError = vnBkg*TMath::Sqrt(TMath::Power(vnerrors_side[i]/vnvalues_side[i],2) + TMath::Power(0.5*vnerrors_h[degree]/vnvalues_h[degree],2));
        double sigError = TMath::Sqrt(vnObsError*vnObsError + TMath::Power(vnBkgError*(1-fsig[i]),2))/fsig[i];

        cout << "Sig: " << sig << endl;
        cout << "Error:" << sigError << endl;

        sigvalues.push_back(sig);
        sigerrors.push_back(sigError);
    }

    myfile << output.str();
    output.str(std::string());
    for(unsigned i=0; i<sigvalues.size(); i++) myfile << sigvalues[i] << "\n";

    output << V0IDname << " signal v" << degree << " errors\n";
    myfile << output.str();
    output.str(std::string());
    for(unsigned i=0; i<sigerrors.size(); i++) myfile << sigerrors[i] << "\n";
}

void Xiv2Fit(  )
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

    std::ofstream Xiv2Peak;
    std::ofstream Xiv2Side;
    std::ofstream Xiv2Calculator;
    Xiv2Peak.open("Xiv2Peak.txt");
    Xiv2Side.open("Xiv2Side.txt");
    Xiv2Calculator.open("Xiv2Signal.txt");

    //TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/XiAnalysisCorrelation.root " );
    //TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/NoPtCut/XiAnalysisCorrelationNoPtCutTotal.root " );
    //TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/NoPtCutPeakAndSide/XiAnalysisCorrelationNoPtCutPeakAndSideTotal.root " );
    //TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/NoPtCutPeakAndSide/XiAnalysisSeparated.root " );
    TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationPD1-6reverseJL10-15_08_15_2017.root" );

    TVirtualFitter::SetMaxIterations( 300000 );
    TH1::SetDefaultSumw2(  );
    std::vector<double> PtBin = {1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 10.0, 20.0};
    std::vector<double> fsig_xi = {0.954019 ,0.973881 ,0.976705 ,0.97829 ,0.978074 ,0.978057 ,0.978603 ,0.974693 ,0.976012};
    int numPtBins = PtBin.size()-1;
    TH1D* dPhiFourierPeak[numPtBins];
    TH1D* dPhiFourierSide[numPtBins];

    TH1D* dPhiPeak[numPtBins];
    TH1D* dPhiSide[numPtBins];

    TF1* FourierFitXi[numPtBins];
    std::vector<double> v2values_peak;
    std::vector<double> v2values_side;
    std::vector<double> v2error_peak;
    std::vector<double> v2error_side;
    std::vector<double> v2value_h; //Need to use vector for vnCalculate function
    std::vector<double> v2error_h;


    TLatex* ltx2 = new TLatex(  );
    ltx2->SetTextSize( 0.045 );
    ltx2->SetNDC( kTRUE );

    //FITTING FOR V2
    //
    //Define divided hist
    bool Peak = true;
    //bool Peak = false;
    for( int i=0; i<numPtBins; i++ )
    {
            dPhiPeak[i] = new TH1D( Form( "dPhiPeak%d",i ), "#Xi - h^{#pm} ", 31, -( 0.5 -
                        1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI  );
            TH1D *dPhiHad = new TH1D( "dPhiHad", "h^{#pm}- h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
            //Pull 2D Histograms
            TH2D *hbackgroundPeak = (TH2D*) f->Get( Form( "xiCorrelationRapidity/BackgroundPeak_pt%d",i ) );
            TH2D *hsignalPeak     = (TH2D*) f->Get( Form( "xiCorrelationRapidity/SignalPeak_pt%d",i ) );
            TH2D *hBackgroundHad  = (TH2D*) f->Get( "xiCorrelationRapidity/BackgroundHad" );
            TH2D *hSignalHad      = (TH2D*) f->Get( "xiCorrelationRapidity/SignalHad" );

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
            dPhiPeak[i]->Divide( hsPhiTotPeak, hbPhiTotPeak );
            dPhiHad->Divide( hsHadPhiTot, hbHadPhiTot );

            //Clone histograms for display without fit functions
            dPhiFourierPeak[i] = ( TH1D* )dPhiPeak[i]->Clone(  );
            TH1D* dPhiHadFourier = ( TH1D* )dPhiHad->Clone(  );

            FourierFitXi[i] = new TF1( Form( "FourierFitXi%d",i ), FourierHad, -1.5, 5, 4 );
            FourierFitXi[i]->SetNpx( 250 );
            FourierFitXi[i]->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *FourierPeak = new TCanvas( "FourierPeak", "Fourier Peak", 800,800 );
            FourierPeak->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            //dPhiFourierPeak->Fit( "FourierFitXi","","",0,PI );
            dPhiFourierPeak[i]->Fit( Form( "FourierFitXi%d",i ) );
            dPhiFourierPeak[i]->SetStats( kFALSE );
            v2values_peak.push_back( FourierFitXi[i]->GetParameter( 2 ) );
            v2error_peak.push_back( FourierFitXi[i]->GetParError( 2 ) );
            cout << "---------------------------------" << endl;
            cout << "Peak V2 for xi-h is " << FourierFitXi[i]->GetParameter( 2 ) << endl;
            cout << "---------------------------------" << endl;

            TF1 *FourierFitHad = new TF1( "FourierFitHad", FourierHad, -1.5, 5, 4 );
            FourierFitHad->SetNpx( 250 );
            FourierFitHad->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *FourierHadron = new TCanvas( "FourierHadron", "Fourier Hadron", 800,800 );
            FourierHadron->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            //dPhiHadFourier->Fit( "FourierFitHad","","",0,PI );
            dPhiHadFourier->Fit( "FourierFitHad");
            dPhiHadFourier->SetStats( kFALSE );
            v2value_h.push_back(FourierFitHad->GetParameter(2));
            v2error_h.push_back(FourierFitHad->GetParError(2));

            double maxBinContent = dPhiFourierPeak[i]->GetBinContent( dPhiFourierPeak[i]->GetMaximumBin(  ) );
            double minBinContent = dPhiFourierPeak[i]->GetBinContent( dPhiFourierPeak[i]->GetMinimumBin(  ) );
            double minRange = minBinContent - 0.005*minBinContent;
            double maxRange = minRange + 2*( maxBinContent - minBinContent );

            //ZYAM FITS
            TCanvas *ZYAMFitPeak = new TCanvas( "ZYAMFitPeak", "ZYAM Fit Peak", 800,800 );
            ZYAMFitPeak->cd(  );

            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiPeak[i]->SetMarkerStyle( 21 );
            dPhiPeak[i]->SetMarkerColor( 4 );
            dPhiPeak[i]->SetTitleOffset( 2, "Y" );
            dPhiPeak[i]->SetTitle( "Peak" );
            dPhiPeak[i]->GetYaxis(  )->SetRangeUser( minRange , maxRange );
            dPhiPeak[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiPeak[i]->GetYaxis(  )->CenterTitle( true );
            dPhiPeak[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} " );
            dPhiPeak[i]->SetTitleOffset( 1.5, "X" );
            dPhiPeak[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiPeak[i]->GetXaxis(  )->CenterTitle( true );
            dPhiPeak[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
            //dPhiPeak[i]->Draw( "E1" );
            dPhiPeak[i]->Fit( "pol2","","", 0.4,2.4 );
            dPhiPeak[i]->SetStats( !publish );

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

            TF1 *dPhiFitPeak = dPhiPeak[i]->GetFunction( "pol2" );
            TF1 *dPhiHadFit = dPhiHad->GetFunction( "pol2" );

            double dPhiFitMinPeak = dPhiFitPeak->GetMinimum(  );
            double dPhiHadFitMin = dPhiHadFit->GetMinimum(  );

            for( int j = 1; j < 32; j++ )
            {
                dPhiFourierPeak[i]->AddBinContent( j, -dPhiFitMinPeak );
            }

            for( int j = 1; j < 32; j++ )
            {
                dPhiHadFourier->AddBinContent( j, -dPhiHadFitMin );
            }

            FourierPeak->cd(  );
            dPhiFourierPeak[i]->SetMarkerStyle( 21 );
            dPhiFourierPeak[i]->SetMarkerColor( 4 );
            //dPhiFourierPeak[i]->Draw( "E1" );
            dPhiFourierPeak[i]->SetStats( kFALSE );
            dPhiFourierPeak[i]->Fit( Form( "FourierFitXi%d",i ) );
            os << "Peak " << PtBin[i] << "_Pt_" << PtBin[i+1];
            dPhiFourierPeak[i]->SetTitle( os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            dPhiFourierPeak[i]->SetTitleOffset( 2, "Y" );
            dPhiFourierPeak[i]->GetYaxis(  )->CenterTitle( true );
            dPhiFourierPeak[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiFourierPeak[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
            dPhiFourierPeak[i]->GetYaxis(  )->SetRangeUser( -0.0004, 0.008 );
            dPhiFourierPeak[i]->SetTitleOffset( 1.5, "X" );
            dPhiFourierPeak[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiFourierPeak[i]->GetXaxis(  )->CenterTitle( true );
            dPhiFourierPeak[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );


            if( publish )
            {
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }

            FourierHadron->cd(  );
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
            //side
            dPhiSide[i] = new TH1D( Form( "dPhiSide%d",i ), "#Xi - h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI  );
            TH2D *hbackgroundSide = (TH2D*) f->Get( Form( "xiCorrelationRapidity/BackgroundSide_pt%d",i ) );
            TH2D *hsignalSide     = (TH2D*) f->Get( Form( "xiCorrelationRapidity/SignalSide_pt%d",i ) );

            TH1::SetDefaultSumw2(  );

            //Project Phi
            TH1D* hbPhiTotSide = hbackgroundSide->ProjectionY( "PhiBkgTot", 0, 10 );
            TH1D* hbPhiOthSide = hbackgroundSide->ProjectionY( "PhiBkgOthPeak", 23, -1 );
            TH1D* hsPhiTotSide = hsignalSide->ProjectionY( "PhiSigTot", 0, 10 );
            TH1D* hsPhiOthSide = hsignalSide->ProjectionY( "PhiSigOthPeak", 23, -1 );

            hbPhiTotSide->Add( hbPhiOthSide );
            hsPhiTotSide->Add( hsPhiOthSide );

            hbHadPhiTot->Add( hbHadPhiOth );
            hsHadPhiTot->Add( hsHadPhiOth );

            //Divide
            dPhiSide[i]->Divide( hsPhiTotSide, hbPhiTotSide );

            //Clone histograms for display without fit functions
            dPhiFourierSide[i] = ( TH1D* )dPhiSide[i]->Clone(  );

            FourierFitXi[i] = new TF1( Form( "FourierFitXi%d",i ) , FourierHad, -1.5, 5, 4 );
            FourierFitXi[i]->SetNpx( 250 );
            FourierFitXi[i]->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *FourierSide = new TCanvas( "FourierSide", "Fourier Side", 800,800 );
            FourierSide->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiFourierSide[i]->Fit( Form( "FourierFitXi%d",i ) );
            dPhiFourierSide[i]->SetStats( kFALSE );
            v2values_side.push_back( FourierFitXi[i]->GetParameter( 2 ) );
            v2error_side.push_back( FourierFitXi[i]->GetParError( 2 ) );
            cout << "---------------------------------" << endl;
            cout << "Side V2 for xi-h is " << FourierFitXi[i]->GetParameter( 2 ) << endl;
            cout << "---------------------------------" << endl;

            maxBinContent = dPhiFourierSide[i]->GetBinContent( dPhiFourierSide[i]->GetMaximumBin(  ) );
            minBinContent = dPhiFourierSide[i]->GetBinContent( dPhiFourierSide[i]->GetMinimumBin(  ) );
            minRange = minBinContent - 0.005*minBinContent;
            maxRange = minRange + 2*( maxBinContent - minBinContent );


            //ZYAM FITS
            TCanvas *ZYAMFitSide = new TCanvas( "ZYAMFitSide", "ZYAM Fit Side", 800,800 );
            ZYAMFitSide->cd(  );

            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiSide[i]->SetMarkerStyle( 21 );
            dPhiSide[i]->SetMarkerColor( 4 );
            dPhiSide[i]->SetTitleOffset( 2, "Y" );
            dPhiSide[i]->SetTitle( "Sideband" );
            dPhiSide[i]->GetYaxis(  )->SetRangeUser( minRange , maxRange );
            dPhiSide[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiSide[i]->GetYaxis(  )->CenterTitle( true );
            dPhiSide[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} " );
            dPhiSide[i]->SetTitleOffset( 1.5, "X" );
            dPhiSide[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiSide[i]->GetXaxis(  )->CenterTitle( true );
            dPhiSide[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
            //dPhiSide[i]->Draw( "E1" );
            dPhiSide[i]->Fit( "pol2","","", 0.4,2.4 );
            dPhiSide[i]->SetStats( !publish );

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

            TF1 *dPhiFitSide = dPhiSide[i]->GetFunction( "pol2" );

            double dPhiFitMinSide = dPhiFitSide->GetMinimum(  );

            for( int j = 1; j < 32; j++ )
            {
                dPhiFourierSide[i]->AddBinContent( j, -dPhiFitMinSide );
            }


            FourierSide->cd(  );
            dPhiFourierSide[i]->SetMarkerStyle( 21 );
            dPhiFourierSide[i]->SetMarkerColor( 4 );
            //dPhiFourierSide[i]->Draw( "E1" );
            dPhiFourierSide[i]->SetStats( kFALSE );
            dPhiFourierSide[i]->Fit( Form( "FourierFitXi%d",i ) );
            os << "SideBand " << PtBin[i] << "_Pt_" << PtBin[i+1];
            dPhiFourierSide[i]->SetTitle( os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            dPhiFourierSide[i]->SetTitleOffset( 2, "Y" );
            dPhiFourierSide[i]->GetYaxis(  )->CenterTitle( true );
            dPhiFourierSide[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiFourierSide[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
            dPhiFourierSide[i]->GetYaxis(  )->SetRangeUser( -0.0004, 0.008 );
            dPhiFourierSide[i]->SetTitleOffset( 1.5, "X" );
            dPhiFourierSide[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiFourierSide[i]->GetXaxis(  )->CenterTitle( true );
            dPhiFourierSide[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );

            if( publish )
            {
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }

    }

    // V2 values Peak
    cout << "==========================================================" << endl;
    cout << "V2 values peak" << endl;
    cout << "==========================================================" << endl;

    int PtBinCounter=0;

    for( std::vector<double>::iterator it = v2values_peak.begin(  ); it != v2values_peak.end(  ); ++it ){
        cout << PtBin[PtBinCounter] << " < Pt =< " << PtBin[PtBinCounter + 1] << ": " << *it << endl;
        PtBinCounter++;
    }


    // V2 errors
    cout << "==========================================================" << endl;
    cout << "V2 errors peak" << endl;
    cout << "==========================================================" << endl;

    PtBinCounter=0;

    for( std::vector<double>::iterator it = v2error_peak.begin(  ); it != v2error_peak.end(  ); ++it ){
        cout << PtBin[PtBinCounter] << " < Pt =< " << PtBin[PtBinCounter + 1] << ": " << *it << endl;
        PtBinCounter++;
    }

    // V2 values Side
    cout << "==========================================================" << endl;
    cout << "V2 values side" << endl;
    cout << "==========================================================" << endl;

    PtBinCounter = 0;

    for( std::vector<double>::iterator it = v2values_side.begin(  ); it != v2values_side.end(  ); ++it ){
        cout << PtBin[PtBinCounter] << " < Pt =< " << PtBin[PtBinCounter + 1] << ": " << *it << endl;
        PtBinCounter++;
    }


    // V2 errors
    cout << "==========================================================" << endl;
    cout << "V2 errors side" << endl;
    cout << "==========================================================" << endl;

    PtBinCounter=0;

    for( std::vector<double>::iterator it = v2error_side.begin(  ); it != v2error_side.end(  ); ++it ){
        cout << PtBin[PtBinCounter] << " < Pt =< " << PtBin[PtBinCounter + 1] << ": " << *it << endl;
        PtBinCounter++;
    }

    //Calculate Flow
    vnCalculate(2,"Xi",v2values_peak,v2error_peak,v2values_side,v2error_side,v2value_h,v2error_h,fsig_xi,"Xiv2Signal.txt");

    TCanvas* FourierPeakComp = new TCanvas( "FourierPeakComp", "Fourier Peak Composite", 1200,1000 );
    FourierPeakComp->Divide( 3,3 );
    for( int i=0; i<numPtBins; i++ ){
        FourierPeakComp->cd( i+1 );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        dPhiFourierPeak[i]->Draw( "E1" );
    }
    TCanvas* FourierSideComp = new TCanvas( "FourierSideComp", "Fourier Side Composite", 1200,1000 );
    FourierSideComp->Divide( 3,3 );
    for( int i=0; i<numPtBins; i++ ){
        FourierSideComp->cd( i+1 );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        dPhiFourierSide[i]->Draw( "E1" );
    }

    //Output Publication plots
    //1D correlation functions
    if( Peak ){
        TCanvas* PubFourier = new TCanvas( "PubFourier", "Pub", 800,800 );
        PubFourier->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        TH1D* dPhiFourierPeakCopy = ( TH1D* )dPhiFourierPeak[4]->Clone(  );
        dPhiFourierPeakCopy->SetTitle( "Peak" );
        dPhiFourierPeakCopy->Draw( "E1" );

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
        TCanvas* PubFourier = new TCanvas( "PubFourier", "Pub", 800,800 );
        PubFourier->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        TH1D* dPhiFourierSideCopy = ( TH1D* )dPhiFourierSide[4]->Clone(  );
        dPhiFourierSideCopy->SetTitle( "SideBand" );
        dPhiFourierSideCopy->Draw( "E1" );


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

    //2D Correlation function 1-3 GeV associated
    TCanvas* TwoDCorrelation = new TCanvas( "TwoDCorrelation", "", 1000, 1000 );
    TwoDCorrelation->SetLeftMargin( 0.2 );

    TH2D* Signal = ( TH2D* )f->Get( "xiCorrelationRapidity/SignalXiHad" );
    TH2D* Background = ( TH2D* )f->Get( "xiCorrelationRapidity/BackgroundXiHad" );

    TGaxis::SetMaxDigits( 1 );
    
    TH2D* Correlation = ( TH2D* )Signal->Clone(  );
    Correlation->Divide( Background );
    Correlation->GetXaxis(  )->SetRangeUser( -4.0,4.0 );
    Correlation->GetYaxis(  )->SetRangeUser( -PI/2.0,4.5 );
    Correlation->GetXaxis(  )->SetTitle( "#Delta#eta" );
    Correlation->GetXaxis(  )->SetTitleOffset( 1.4 );
    Correlation->GetXaxis(  )->CenterTitle( true );
    Correlation->GetYaxis(  )->SetTitle( "#Delta#phi (radians)" );
    Correlation->GetYaxis(  )->SetTitleOffset( 1.4 );
    Correlation->GetYaxis(  )->CenterTitle( true );
    Correlation->GetZaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{d^{2}N^{pair}}{d#Delta#eta d#Delta#phi}" );
    Correlation->GetZaxis(  )->SetTitleOffset( 2.3 );
    Correlation->GetZaxis(  )->CenterTitle( true );
    Correlation->GetXaxis(  )->SetNdivisions( 405 );
    Correlation->GetYaxis(  )->SetNdivisions( 405 );
    Correlation->GetZaxis(  )->SetNdivisions( 4 );
    Correlation->SetTitle( "" );
    Correlation->SetStats( kFALSE );
    Correlation->Scale( 10 );

    const Int_t NRGBs = 5;
    const Int_t NCont = 20;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    //gStyle->SetNumberContours( 90 );
    
    TwoDCorrelation->cd(  );
    Correlation->Draw( "SURF1 FB " );

    TLatex *ltx0 = new TLatex(  );
    ltx0->SetTextSize( 0.031 );
    ltx0->SetNDC( kTRUE );
    ltx0->SetTextFont( 42 );

    ltx0->DrawLatex( 0.05, 0.95, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}" );
    ltx0->DrawLatex( 0.05, 0.88, "185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
    ltx0->DrawLatex( 0.05, 0.81, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
    ltx0->DrawLatex( 0.05, 0.75, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{trig}}} < 3 GeV" );
    ltx0->SetTextSize( 0.04 );
    ltx0->DrawLatex( 0.85, 0.88, "#Xi#kern[-0.3]{#lower[0.02]{{}^{#pm}}}- h^{#pm}" );

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

    

}

