#include "GetGraphFromFile.C"
#include "MITStyle.C"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"

void FinalResult(  )
{
    //MITStyle(  );
    //double v2Cas[7] = {0.056368, 0.064876, 0.109999, 0.126307, 0.153165, 0.182298, 0.217702};
    //double pTCas[7] = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055};
    //double v2CasE[7] = {0.017074, 0.008531, 0.007160, 0.005156, 0.004894, 0.006175, 0.010446};
    double v2Cas[7] = {0.051751, 0.071787, 0.107110, 0.126828, 0.157772, 0.181781, 0.217510};
    double pTCas[7] = {1.267, 1.62, 2.007, 2.5, 3.173, 4.029, 5.054};
    double v2CasE[7] = {0.015362, 0.007650, 0.006419, 0.004632, 0.004404, 0.005550, 0.009393};

    // Pull TGraph for Kshort and lambda

    TFile* file_pPbv2 = TFile::Open( "lrgraphv2_v3_pPb_185-220.root" );
    TFile* file_hadv2 = TFile::Open( "lrgraphv2_v3_pPb_hadron_185-above.root" );

    TGraphErrors* ks_v2 = ( TGraphErrors* )file_pPbv2->Get( "kshortv2true" );
    TGraphErrors* la_v2 = ( TGraphErrors* )file_pPbv2->Get( "lambdav2true" );
    //TGraphErrors* ha_v2 = ( TGraphErrors* )file_hadv2->Get( "hadronv2" );
    TGraphErrors* ha_v2 = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033pPb_v2.txt"),1,28,1.2);
    
    TGraphErrors* xi_v2  = new TGraphErrors(7,pTCas,v2Cas,0,v2CasE);

    ks_v2->SetMarkerColor( kRed-3 );
    ks_v2->SetLineColor( kRed-3 );
    //ks_v2->SetMarkerStyle( 33 );
    ks_v2->SetMarkerStyle( 25 );
    ks_v2->SetMarkerSize( 1.4 );
    la_v2->SetMarkerColor( kBlue-4 );
    //la_v2->SetMarkerStyle( 22 );
    la_v2->SetMarkerStyle( 26 );
    la_v2->SetMarkerSize( 1.3 );
    la_v2->SetLineColor( kBlue );
    ha_v2->SetMarkerStyle( 28 );
    ha_v2->SetMarkerSize( 1.3 );
    xi_v2->SetMarkerColor( kRed );
    //xi_v2->SetMarkerColor( kTeal-2 );
    //xi_v2->SetMarkerColor( kMagenta-4 );
    //xi_v2->SetMarkerColor( kGreen+2 );
    xi_v2->SetMarkerStyle( 20 );
    xi_v2->SetMarkerSize( 1.5 );
    xi_v2->SetLineColor( kRed );

    MakeCanvas( "c1", "Plot" );
    c1->cd(  );
    c1->SetLeftMargin( 0.12 );

    // draw the frame using a histogram frame

    TH1F* frame = c1->DrawFrame( 0,-0.01,6,0.31 );
    gPad->SetTickx(  );
    gPad->SetTicky(  );
    frame->GetXaxis(  )->CenterTitle( 1 );
    frame->GetYaxis(  )->CenterTitle( 1 );
    frame->GetXaxis(  )->SetTitleSize( 0.05 );
    frame->GetXaxis(  )->SetTitle( "p_{T} (GeV)" );
    frame->GetYaxis(  )->SetTitle( "v#kern[-0.3]{_{2}}" );
    frame->GetYaxis(  )->SetTitleSize( 0.05 );
    frame->SetTitleOffset( 1.2,"Y" );
    frame->SetTitleOffset( 1.2,"X" );



    // Plot systematic errors for ks and la
    TBox* box;
    double syst_pPb = 0.069;
    double syst_pPbH = 0.039;
    double xOffset = 0.08;

    int n = ks_v2->GetN(  );
    for( int j=0; j<n; j++ )
    {
        double xk,yk,xl,yl,xh,yh = 0;
        if(la_v2->GetN()-1 >= j){
            la_v2->GetPoint(j,xl,yl);
            double percent1L = syst_pPb;
        }
        if(ha_v2->GetN()-1 >= j)
        {
            ha_v2->GetPoint(j,xh,yh);
            double percentH = syst_pPbH;
        }
        ks_v2->GetPoint(j,xk,yk);
        double percent1K = syst_pPb;
        double yLerr, yKerr, yHerr;

        yKerr=fabs(yk)*percent1K;
        yLerr=fabs(yl)*percent1L;
        yHerr=fabs(yh)*percentH;

        if(la_v2->GetN()-1 >= j)
        {
            box = new TBox(xl-xOffset,yl-yLerr,xl+xOffset,yl+yLerr);
            box->SetFillColor(17);
            box->SetFillStyle(1001);
            box->SetLineWidth(0);
            box->Draw("SAME");
        }
        if(ha_v2->GetN()-1 >= j){
            box = new TBox(xh-xOffset,yh-yHerr,xh+xOffset,yh+yHerr);
            box->SetFillColor(17);
            box->SetFillStyle(1001);
            box->SetLineWidth(0);
            box->Draw("SAME");
        }
        box = new TBox(xk-xOffset,yk-yKerr,xk+xOffset,yk+yKerr);
        box->SetFillColor(17);
        box->SetFillStyle(1001);
        box->SetLineWidth(0);
        box->Draw("SAME");

    }
    //
    // draw the graph with axis, continuous line, and put
    ks_v2->Draw( "PESAME" );
    la_v2->Draw( "PESAME" );
    ha_v2->Draw( "PESAME" );
    xi_v2->Draw("P");
    Int_t oldColor = 38;
    Int_t newColor = 46;

    TLegend* leg = new TLegend( 0.2,0.45,0.51,0.68 );
    leg->SetFillColor( 10 );
    leg->SetFillStyle( 0 );
    leg->SetBorderSize( 0.035 );
    leg->SetTextFont( 42 );
    leg->SetTextSize( 0.05 );
    leg->AddEntry( ks_v2,"#color[38]{K_{#lower[-0.3]{S}}#kern[-1.05]{#lower[0.1]{{}^{0}}}}", "P" );
    leg->AddEntry( la_v2, "#color[38]{#Lambda / #lower[0.1]{#bar{  }}#kern[-1.2]{#Lambda}}","P" );
    leg->AddEntry( ha_v2, "#color[38]{h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}}", "P" );
    leg->AddEntry( xi_v2, "#color[46]{#Xi#kern[-0.3]{#lower[0.1]{{}^{+}}}/ #Xi#kern[-0.3]{#lower[0.1]{{}^{-}}}}", "P" );
    leg->Draw(  );

    TLatex *tex = new TLatex(  );
    tex->SetNDC(  );
    tex->SetTextFont( 42 );
    tex->SetTextSize( 0.05 );
    tex->DrawLatex( 0.23,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = #color[38]{5.02} TeV, #color[46]{8.16} TeV" );
    tex->SetTextSize( 0.045 );
    tex->DrawLatex( 0.23,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}" );
   // tex->DrawLatex( 0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}" );
    tex->DrawLatex( 0.55,0.23,"185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
    //tex->DrawLatex( 0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}" );
    TLine* line = new TLine( 0,0,6,0 );
    line->SetLineStyle( 8 );
    line->Draw(  );

}
