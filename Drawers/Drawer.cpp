#include <iostream>
#include "TH1D.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGaxis.h"

#define PI 3.1416

void HM2DKs()
{
    TH1::SetDefaultSumw2();
    gStyle->SetTitleFontSize(0.04);

    TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityCorrectMultB_09_19_17.root");
    //
    //================================================================================
    //KET Calculations
    //================================================================================
    int i = 13;
    TH1D* hKetKs = (TH1D*)f->Get("v0CorrelationRapidity/KETkshort_pt10");
    TH1D* hKetKs_bkg = (TH1D*)f->Get("v0CorrelationRapidity/KETkshort_bkg_pt10");

    cout << "First" <<  hKetKs->FindFirstBinAbove( 0 , 1) << endl;
    cout << "Last" <<  hKetKs->FindLastBinAbove( 0, 1 ) << endl;

    TLatex *ltx0 = new TLatex(  );
    ltx0->SetNDC( kTRUE );
    ltx0->SetTextFont( 62 );
    //int nEntries = 0;
    //double KetTotal = 0;
    //for(int j=hKetKs; j<20000; j++)
    //{
        //double nKet = hKetKs->GetBinContent(i);
        //double nKet_bkg = hKetKs_bkg->GetBinContent(i);
        //double Ket = nKet*(hKetKs->GetBinCenter(i));
        //double Ket_bkg = nKet_bkg*(hKetKs_bkg->GetBinCenter(i));
        //double TotalKet = Ket + Ket_bkg;
        //if(TotalKet == 0) continue;
        //KetTotal += TotalKet;
        //nEntries++;
    //}
    //AvgKetKs.push_back(KetTotal/nEntries);
    //
    //
    const Int_t NRGBs = 5;
    const Int_t NCont = 20;
    Double_t stops[NRGBs] = { 0.00, 0.15, 0.51, 0.75, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.77, 1.00, 0.91 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.91, 1.00, 0.06, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

    TCanvas* c1 = new TCanvas("c1","c1",800,800);
    gPad->SetTheta(50);
    gPad->SetPhi(45);
    gPad->Update();
    c1->SetLeftMargin(0.15);
    TH2D* hsignal_ks = (TH2D*)f->Get("v0CorrelationRapidity/signalkshort_pt8");
    TH2D* hbackground_ks = (TH2D*)f->Get("v0CorrelationRapidity/backgroundkshort_pt8");
    TH1D* hmult = (TH1D*)f->Get("v0CorrelationRapidity/mult_ks_pt8");

    double BW2D = (9.9/33)*((2-1.0/16)*3.1416/31);
    double nEvents = hmult->Integral(2,10000);
    double b0 = hbackground_ks->GetBinContent(hbackground_ks->FindBin(0,0));
    hsignal_ks->Divide(hbackground_ks);
    hsignal_ks->Scale(b0/nEvents/BW2D);

    hsignal_ks->GetXaxis()->SetRangeUser(-PI,PI);
    hsignal_ks->GetYaxis()->SetRangeUser(-PI/2,1.5*PI);
    hsignal_ks->GetXaxis()->CenterTitle(1);
    hsignal_ks->GetYaxis()->CenterTitle(1);
    hsignal_ks->GetZaxis()->CenterTitle(1);
    hsignal_ks->SetMaximum(8.3);
    hsignal_ks->GetXaxis()->SetTitleOffset(1.6);
    hsignal_ks->GetYaxis()->SetTitleOffset(1.6);
    hsignal_ks->GetZaxis()->SetTitleOffset(1.6);
    hsignal_ks->GetXaxis()->SetTitle("#Delta#eta");
    hsignal_ks->GetYaxis()->SetTitle("#Delta#phi (radians)");
    hsignal_ks->GetZaxis()->SetTitle("#frac{1}{N_{trig}} #frac{d^{2} N^{pair}}{d#Delta#eta d#Delta#phi}");
    hsignal_ks->GetXaxis()->SetNdivisions(407);
    hsignal_ks->GetYaxis()->SetNdivisions(407);
    hsignal_ks->GetZaxis()->SetNdivisions(8);
    hsignal_ks->SetStats(kFALSE);

    double xpos = 0.02;
    double ypos = 1.02;
    double increment = 0.06;
    hsignal_ks->Draw("SURF1 FB");

    ltx0->SetTextSize(0.040);
    ltx0->DrawLatex(xpos,ypos-=increment,"CMS pPb #sqrt{S_{NN}} = 8.16 TeV");
    ltx0->SetTextSize( 0.030 );
    ltx0->DrawLatex(xpos,ypos-=increment,"185 #leq N_{trk}^{Offline} < 250");
    ltx0->DrawLatex(xpos,ypos-=increment,"2.8 < p_{T} < 3.6 GeV/c");
    ltx0->DrawLatex(xpos,ypos-=increment,"0.3 < p_{T}^{assoc} < 3.0 GeV/c");
    ltx0->DrawLatex(xpos,ypos-=increment,"|y| < 1");

    ltx0->SetTextSize(0.045);
    ltx0->DrawLatex(0.78,0.89,"K_{S}^{0} - h^{#pm}");

    c1->Print("KsHM2DForZhenyu.pdf");
}

void MB2DKs()
{
    TH1::SetDefaultSumw2();
    gStyle->SetTitleFontSize(0.04);

    TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/MBCorr/PeripheralSubtractionMB_0_n_20_V0Only.root");

    TLatex *ltx0 = new TLatex(  );
    ltx0->SetNDC( kTRUE );
    ltx0->SetTextFont( 62 );

    const Int_t NRGBs = 5;
    const Int_t NCont = 20;
    Double_t stops[NRGBs] = { 0.00, 0.15, 0.51, 0.75, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.77, 1.00, 0.91 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.91, 1.00, 0.06, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

    TCanvas* c1 = new TCanvas("c1","c1",800,800);
    gPad->SetTheta(50);
    gPad->SetPhi(45);
    gPad->Update();
    c1->SetLeftMargin(0.15);
    TH2D* hsignal_ks = (TH2D*)f->Get("v0CasCorrelationRapidityPeriSub/signalkshort_pt8");
    TH2D* hbackground_ks = (TH2D*)f->Get("v0CasCorrelationRapidityPeriSub/backgroundkshort_pt8");
    TH1D* hmult = (TH1D*)f->Get("v0CasCorrelationRapidityPeriSub/mult_ks_pt8");

    double BW2D = (9.9/33)*((2-1.0/16)*3.1416/31);
    double nEvents = hmult->Integral(2,10000);
    double b0 = hbackground_ks->GetBinContent(hbackground_ks->FindBin(0,0));
    hsignal_ks->Divide(hbackground_ks);
    hsignal_ks->Scale(b0/nEvents/BW2D);

    hsignal_ks->GetXaxis()->SetRangeUser(-PI,PI);
    hsignal_ks->GetYaxis()->SetRangeUser(-PI/2,1.5*PI);
    hsignal_ks->GetXaxis()->CenterTitle(1);
    hsignal_ks->GetYaxis()->CenterTitle(1);
    hsignal_ks->GetZaxis()->CenterTitle(1);
    hsignal_ks->SetMaximum(0.95);
    hsignal_ks->GetXaxis()->SetTitleOffset(1.6);
    hsignal_ks->GetYaxis()->SetTitleOffset(1.6);
    hsignal_ks->GetZaxis()->SetTitleOffset(1.6);
    hsignal_ks->GetXaxis()->SetTitle("#Delta#eta");
    hsignal_ks->GetYaxis()->SetTitle("#Delta#phi (radians)");
    hsignal_ks->GetZaxis()->SetTitle("#frac{1}{N_{trig}} #frac{d^{2} N^{pair}}{d#Delta#eta d#Delta#phi}");
    hsignal_ks->GetXaxis()->SetNdivisions(407);
    hsignal_ks->GetYaxis()->SetNdivisions(407);
    hsignal_ks->GetZaxis()->SetNdivisions(8);
    hsignal_ks->SetStats(kFALSE);

    double xpos = 0.02;
    double ypos = 1.02;
    double increment = 0.06;
    hsignal_ks->Draw("SURF1 FB");

    ltx0->SetTextSize(0.040);
    ltx0->DrawLatex(xpos,ypos-=increment,"CMS pPb #sqrt{S_{NN}} = 8.16 TeV");
    ltx0->SetTextSize( 0.030 );
    ltx0->DrawLatex(xpos,ypos-=increment,"0 #leq N_{trk}^{Offline} < 20");
    ltx0->DrawLatex(xpos,ypos-=increment,"2.8 < p_{T} < 3.6 GeV/c");
    ltx0->DrawLatex(xpos,ypos-=increment,"0.3 < p_{T}^{assoc} < 3.0 GeV/c");
    ltx0->DrawLatex(xpos,ypos-=increment,"|y| < 1");

    ltx0->SetTextSize(0.045);
    ltx0->DrawLatex(0.78,0.89,"K_{S}^{0} - h^{#pm}");

    c1->Print("KsMB2DForZhenyu.pdf");
}
