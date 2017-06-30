#include "TH1D.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGaxis.h"

void Drawer()
{
    TH1::SetDefaultSumw2(  );
    gStyle->SetTitleFontSize(0.04);

    //TFile* f1 = new TFile("/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/XiPtBins.root");
    //TFile* f1 = new TFile("/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/HLT185_220_XiMass_Pt.root");
    //TFile* f1 = new TFile("/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/XiAnalysisCorrelation.root");
    //TFile* f1 = new TFile("/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/NoPtCut/XiAnalysisCorrelationNoPtCutTotal.root");
    //TFile* f1 = new TFile("/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/NoPtCut/XiAnalysisCorrelationNoPtCutTotal_Prelim.root");
    //TFile* f1 = new TFile( "/volumes/MacHD/Users/blt1/research/TestRootFiles/XiAnalysisCorrelationNoPtCut8TeV_1.root" );
    TFile* f1 = new TFile("/volumes/MacHD/Users/blt1/research/TestRootFiles/XiAnalysisCorrelationPtCut8TeVPD1_4_ForFinal.root");

    TH1D* Inv_Mass = (TH1D*) f1->Get("xiCorrelation/InvMassXi");

    Inv_Mass->SetTitle( "Invariant Mass of Charged #Xi" );
    Inv_Mass->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Inv_Mass->GetYaxis(  )->SetTitleOffset( 1.4 );
    Inv_Mass->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Inv_Mass->SetMarkerStyle(20);
    Inv_Mass->SetMarkerSize( 0.7 );

    TH1D* Pt_12 = (TH1D*) f1->Get("xiCorrelation/pT_12");
    Pt_12->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_12->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_12->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Pt_12->SetMarkerStyle(20);
    Pt_12->SetMarkerSize( 0.7 );

    TH1D* Pt_13 = (TH1D*) f1->Get("xiCorrelation/pT_13");
    Pt_13->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_13->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_13->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Pt_13->SetMarkerStyle(20);
    Pt_13->SetMarkerSize( 0.7 );

    TH1D* Pt_14 = (TH1D*) f1->Get("xiCorrelation/pT_14");
    Pt_14->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_14->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_14->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Pt_14->SetMarkerStyle(20);
    Pt_14->SetMarkerSize( 0.7 );

    TH1D* Pt_15 = (TH1D*) f1->Get("xiCorrelation/pT_15");
    Pt_15->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_15->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_15->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Pt_15->SetMarkerStyle(20);
    Pt_15->SetMarkerSize( 0.7 );
    
    TH1D* Pt_16 = (TH1D*) f1->Get("xiCorrelation/pT_16");
    Pt_16->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_16->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_16->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Pt_16->SetMarkerStyle(20);
    Pt_16->SetMarkerSize( 0.7 );

    TH1D* Pt_23 = (TH1D*) f1->Get("xiCorrelation/pT_23");
    Pt_23->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_23->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_23->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Pt_23->SetMarkerStyle(20);
    Pt_23->SetMarkerSize( 0.7 );

    TH1D* Pt_24 = (TH1D*) f1->Get("xiCorrelation/pT_24");
    Pt_24->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_24->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_24->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Pt_24->SetMarkerStyle(20);
    Pt_24->SetMarkerSize( 0.7 );

    TH1D* Pt_25 = (TH1D*) f1->Get("xiCorrelation/pT_25");
    Pt_25->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_25->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_25->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Pt_25->SetMarkerStyle(20);
    Pt_25->SetMarkerSize( 0.7 );

    TH1D* Pt_26 = (TH1D*) f1->Get("xiCorrelation/pT_26");
    Pt_26->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_26->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_26->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Pt_26->SetMarkerStyle(20);
    Pt_26->SetMarkerSize( 0.7 );

    TH1D* Pt_34 = (TH1D*) f1->Get("xiCorrelation/pT_34");
    Pt_34->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_34->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_34->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Pt_34->SetMarkerStyle(20);
    Pt_34->SetMarkerSize( 0.7 );

    TH1D* Pt_35 = (TH1D*) f1->Get("xiCorrelation/pT_35");
    Pt_35->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_35->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_35->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Pt_35->SetMarkerStyle(20);
    Pt_35->SetMarkerSize( 0.7 );

    TH1D* Pt_36 = (TH1D*) f1->Get("xiCorrelation/pT_36");
    Pt_36->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_36->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_36->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Pt_36->SetMarkerStyle(20);
    Pt_36->SetMarkerSize( 0.7 );

    TH1D* Pt_45 = (TH1D*) f1->Get("xiCorrelation/pT_45");
    Pt_45->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_45->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_45->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Pt_45->SetMarkerStyle(20);
    Pt_45->SetMarkerSize( 0.7 );

    TH1D* Pt_46 = (TH1D*) f1->Get("xiCorrelation/pT_46");
    Pt_46->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_46->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_46->GetXaxis(  )->SetTitle("Invariant Mass (GeV)");
    Pt_46->SetMarkerStyle(20);
    Pt_46->SetMarkerSize( 0.7 );

    TH1D* Pt_Spec = (TH1D*) f1->Get("xiCorrelation/pT_Xi");
    Pt_Spec->GetYaxis(  )->SetTitle("Candidates/0.001 (GeV)");
    Pt_Spec->GetYaxis(  )->SetTitleOffset( 1.4 );
    Pt_Spec->GetXaxis(  )->SetTitle("pT (GeV)");
    Pt_Spec->SetMarkerStyle(20);
    Pt_Spec->SetMarkerSize( 0.7 );

    /*
    TH2D* Background = ( TH2D* ) f1->Get( "xiCorrelation/CorrelationSide" );
    Background->SetTitle( "Background" );
    Background->GetXaxis(  )->SetRangeUser( -4,4 );
    Background->GetYaxis(  )->SetRangeUser( -4,4.5 );
    Background->GetXaxis(  )->SetTitle( "#Delta#eta" );
    Background->GetXaxis(  )->SetTitleOffset( 1.4 );
    Background->GetYaxis(  )->SetTitle( "#Delta#phi" );
    Background->GetYaxis(  )->SetTitleOffset( 1.4 );

    TH2D* Signal = ( TH2D* ) f1->Get( "xiCorrelation/CorrelationPeak" );
    Signal->SetTitle( "Signal" );
    Signal->GetXaxis(  )->SetRangeUser( -4,4 );
    Signal->GetYaxis(  )->SetRangeUser( -4,4.5 );
    Signal->GetXaxis(  )->SetTitle( "#Delta#eta" );
    Signal->GetXaxis(  )->SetTitleOffset( 1.4 );
    Signal->GetYaxis(  )->SetTitle( "#Delta#phi" );
    Signal->GetYaxis(  )->SetTitleOffset( 1.4 );

    TH2D* Correlation = ( TH2D* ) f1->Get( "xiCorrelation/CorrelationHad" );
    Correlation->GetXaxis(  )->SetRangeUser( -4,4 );
    Correlation->GetYaxis(  )->SetRangeUser( -4,4.5 );
    Correlation->GetXaxis(  )->SetTitle( "#Delta#eta" );
    Correlation->GetXaxis(  )->SetTitleOffset( 1.4 );
    Correlation->GetYaxis(  )->SetTitle( "#Delta#phi" );
    Correlation->GetYaxis(  )->SetTitleOffset( 1.4 );
    */

    TGaxis::SetMaxDigits( 1 );
    TH2D* Correlation = ( TH2D* ) f1->Get( "xiCorrelation/Correlation" );
    Correlation->GetXaxis(  )->SetRangeUser( -4,4 );
    Correlation->GetYaxis(  )->SetRangeUser( -4,4.5 );
    Correlation->GetXaxis(  )->SetTitle( "#Delta#eta" );
    Correlation->GetXaxis(  )->SetTitleOffset( 1.4 );
    Correlation->GetYaxis(  )->SetTitle( "#Delta#phi" );
    Correlation->GetYaxis(  )->SetTitleOffset( 1.4 );
    Correlation->GetZaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{d^{2}N^{pair}}{d#Delta#phi d#Delta#eta}" );
    Correlation->GetZaxis(  )->SetTitleOffset( 2 );
    Correlation->GetZaxis(  )->SetNdivisions( 4, kFALSE );
    Correlation->SetTitle( "" );
    Correlation->SetStats( kFALSE );


    //InvMass by itself
    TCanvas* c1 = new TCanvas("c1", "", 800, 600);
    c1->cd();
    gPad->SetTickx();
    gPad->SetTicky();
    Inv_Mass->Draw("hist p E1");

    //InvMass with pT spectrum
    TCanvas* c2 = new TCanvas("c2", "", 1600, 600);
    c2->Divide( 2,1 );
    c2->cd( 1 );
    gPad->SetTickx();
    gPad->SetTicky();
    Inv_Mass->Draw("hist p E1");
    c2->cd( 2 );
    gPad->SetTickx();
    gPad->SetTicky();
    Pt_Spec->Draw( "hist p E1" );

    //pT interval 1's
    TCanvas* c3 = new TCanvas("c3", "", 1600, 800);
    c3->Divide( 3,2 );

    c3->cd( 1 );
    gPad->SetTickx();
    gPad->SetTicky();
    Pt_12->Draw( "hist p E1" );

    c3->cd( 2 );
    gPad->SetTickx();
    gPad->SetTicky();
    Pt_13->Draw( "hist p E1" );

    c3->cd( 3 );
    gPad->SetTickx();
    gPad->SetTicky();
    Pt_14->Draw( "hist p E1" );

    c3->cd( 4 );
    gPad->SetTickx();
    gPad->SetTicky();
    Pt_15->Draw( "hist p E1" );

    c3->cd( 5 );
    gPad->SetTickx();
    gPad->SetTicky();
    Pt_16->Draw( "hist p E1" );

    //pT interval 2's
    TCanvas* c4 = new TCanvas("c4", "", 1600, 800);
    c4->Divide( 2,2 );

    c4->cd( 1 );
    gPad->SetTickx();
    gPad->SetTicky();
    Pt_23->Draw( "hist p E1" );

    c4->cd( 2 );
    gPad->SetTickx();
    gPad->SetTicky();
    Pt_24->Draw( "hist p E1" );

    c4->cd( 3 );
    gPad->SetTickx();
    gPad->SetTicky();
    Pt_25->Draw( "hist p E1" );

    c4->cd( 4 );
    gPad->SetTickx();
    gPad->SetTicky();
    Pt_26->Draw( "hist p E1" );

    //pT interval 3's
    TCanvas* c5 = new TCanvas("c5", "", 1600, 800);
    c5->Divide( 2,2 );

    c5->cd( 1 );
    gPad->SetTickx();
    gPad->SetTicky();
    Pt_34->Draw( "hist p E1" );

    c5->cd( 2 );
    gPad->SetTickx();
    gPad->SetTicky();
    Pt_35->Draw( "hist p E1" );

    c5->cd( 3 );
    gPad->SetTickx();
    gPad->SetTicky();
    Pt_36->Draw( "hist p E1" );

    //Background and Signal
    TCanvas* c6 = new TCanvas( "c6", "", 1600,800 );
    c6->Divide( 2,1 );

    /*
    c6->cd( 1 );
    //gPad->SetTickx();
    //gPad->SetTicky();
    Background->Draw( "SURF1" );

    c6->cd( 2 );
    //gPad->SetTickx();
    //gPad->SetTicky();
    Signal->Draw( "SURF1" );
    */

    //Correlation
    TCanvas* c7 = new TCanvas( "c7", "", 1000,1000 );
    c7->cd(  );
    c7->SetLeftMargin( 0.2 );
    Correlation->Draw( "SURF1" );

    TLatex *ltx1 = new TLatex(  );
    ltx1->SetTextSize( 0.035 );
    ltx1->SetNDC( kTRUE );
    ltx1->SetTextFont( 42 );

    ltx1->DrawLatex( 0.05, 0.95, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV, L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
    ltx1->DrawLatex( 0.05, 0.88, "185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
    ltx1->DrawLatex( 0.05, 0.81, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
    ltx1->DrawLatex( 0.05, 0.75, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{trig}}} < 3 GeV" );
    ltx1->DrawLatex( 0.85, 0.88, "#Xi^{#pm} - h^{#pm}" );

    //TCanvas* c6 = new TCanvas("c5", "", 1600, 800);
    //c5->Divide( 2,2 );
    //TCanvas* c3 = new TCanvas("c3", "", 800, 600);
    //c3->cd(  );
    //gPad->SetTickx();
    //gPad->SetTicky();
    //Pt_14->Draw( "hist p E1" );

    //TCanvas* c4 = new TCanvas("c4", "", 800, 600);
    //c4->cd(  );
    //gPad->SetTickx();
    //gPad->SetTicky();
    //Pt_24->Draw( "hist p E1" );

    //TCanvas* c5 = new TCanvas("c5", "", 800, 600);
    //c5->cd(  );
    //gPad->SetTickx();
    //gPad->SetTicky();
    //Pt_25->Draw( "hist p E1" );
}
