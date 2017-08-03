//Includes
#include "TFile.h"
#include "TGraph.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TH2.h"
#include "GetGraphFromFile.C"
#include "MITStyle.C"

void KslaCrossCheckv3(  )
{
    MITStyle(  );
	const int ks_npoints = 10;
    double v2Ks8[ks_npoints] = {0.005070706 ,0.013541808 ,0.017894052 ,0.026917717 ,0.036227128 ,0.045597329 ,0.051020426 ,0.045881397 ,0.037791439 ,0.014758437};
    double pTKs8[ks_npoints] = {0.5213, 0.7084, 0.9036, 1.201, 1.591, 1.986, 2.465, 3.136, 4.008, 5.071};
    double v2Ks8E[ks_npoints] = {0.003412278 ,0.002043 ,0.001651186 ,0.001085994 ,0.001197198 ,0.001463774 ,0.001608019 ,0.002147082 ,0.003280849 ,0.005174682};

	const int la_npoints = 8;
    double v2La8[la_npoints] = {0.011261892 ,0.020795145 ,0.033555289 ,0.044550178 ,0.058043265 ,0.060781616 ,0.073580978 ,0.046895594};
    double pTLa8[la_npoints] = {0.9145, 1.221, 1.605, 1.997, 2.481, 3.149, 4.004, 5.041};
    double v2La8E[la_npoints] = {0.007693167 ,0.003225583 ,0.00266863 ,0.002604339 ,0.002402621 ,0.002821647 ,0.004176484 ,0.007135051};

    // Pull TGraph for Kshort and lambda

    TFile* file_pPbv2 = TFile::Open( "lrgraphv2_v3_pPb_185-220.root" );
    TFile* file_hadv2 = TFile::Open( "lrgraphv2_v3_pPb_hadron_185-above.root" );

    TGraphErrors* ks_v2 = ( TGraphErrors* )file_pPbv2->Get( "kshortv3true" );
    TGraphErrors* la_v2 = ( TGraphErrors* )file_pPbv2->Get( "lambdav3true" );
    //TGraphErrors* ha_v2 = ( TGraphErrors* )file_hadv2->Get( "hadronv2" );
    /*TGraphErrors* ha_v2 = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033pPb_v2.txt"),1,28,1.2);*/

    TGraphErrors* ks8_v2  = new TGraphErrors(ks_npoints,pTKs8,v2Ks8,0,v2Ks8E);
	TGraphErrors* la8_v2  = new TGraphErrors(la_npoints,pTLa8,v2La8,0,v2La8E);

		ks_v2->SetMarkerColor( kBlue-4 );
		ks_v2->SetLineColor( kBlue-4 );
		ks_v2->SetMarkerStyle( 24 );
		ks_v2->SetMarkerSize( 1.4 );
		la_v2->SetMarkerColor( kBlue-4 );
		//la_v2->SetMarkerStyle( 22 );
		la_v2->SetMarkerStyle( 26 );
		la_v2->SetMarkerSize( 1.3 );
		la_v2->SetLineColor( kBlue );
		ks8_v2->SetMarkerColor( kRed );
		//xi_v2->SetMarkerColor( kTeal-2 );
		//xi_v2->SetMarkerColor( kMagenta-4 );
		//xi_v2->SetMarkerColor( kGreen+2 );
		ks8_v2->SetMarkerStyle( 20 );
		ks8_v2->SetMarkerSize( 1.5 );
		ks8_v2->SetLineColor( kRed );

		la8_v2->SetMarkerColor( kGreen+2 );
		//xi_v2->SetMarkerColor( kTeal-2 );
		//xi_v2->SetMarkerColor( kMagenta-4 );
		//xi_v2->SetMarkerColor( kGreen+2 );
		la8_v2->SetMarkerStyle( 22 );
		la8_v2->SetMarkerSize( 1.5 );
		la8_v2->SetLineColor( kBlue-4 );

		TCanvas* c1 = MakeCanvas( "c1", "Plot" );
		c1->cd(  );
		c1->SetLeftMargin( 0.12 );

		// draw the frame using a histogram frame

		TH1F* frame = c1->DrawFrame( 0,-0.01,6,0.15 );
		gPad->SetTickx(  );
		gPad->SetTicky(  );
		frame->GetXaxis(  )->CenterTitle( 1 );
		frame->GetYaxis(  )->CenterTitle( 1 );
		frame->GetXaxis(  )->SetTitleSize( 0.05 );
		frame->GetXaxis(  )->SetTitle( "p_{T} (GeV)" );
		frame->GetYaxis(  )->SetTitle( "v#kern[-0.3]{_{3}}" );
		frame->GetYaxis(  )->SetTitleSize( 0.05 );
		frame->SetTitleOffset( 1.2,"Y" );
		frame->SetTitleOffset( 1.2,"X" );



		// Plot systematic errors for ks and la
		/*TBox* box;*/
		/*double syst_pPb = 0.069;*/
		/*double syst_pPbH = 0.039;*/
		/*double xOffset = 0.08;*/

		/*int n = ks_v2->GetN(  );*/
		/*for( int j=0; j<n; j++ )*/
		/*{*/
				/*double xk,yk,xl,yl,xh,yh = 0;*/
				/*if(la_v2->GetN()-1 >= j){*/
						/*la_v2->GetPoint(j,xl,yl);*/
						/*double percent1L = syst_pPb;*/
				/*}*/
				/*if(ha_v2->GetN()-1 >= j)*/
				/*{*/
						/*ha_v2->GetPoint(j,xh,yh);*/
						/*double percentH = syst_pPbH;*/
				/*}*/
				/*ks_v2->GetPoint(j,xk,yk);*/
				/*double percent1K = syst_pPb;*/
				/*double yLerr, yKerr, yHerr;*/

				/*yKerr=fabs(yk)*percent1K;*/
				/*yLerr=fabs(yl)*percent1L;*/
				/*yHerr=fabs(yh)*percentH;*/

				/*if(la_v2->GetN()-1 >= j)*/
				/*{*/
						/*box = new TBox(xl-xOffset,yl-yLerr,xl+xOffset,yl+yLerr);*/
						/*box->SetFillColor(17);*/
						/*box->SetFillStyle(1001);*/
						/*box->SetLineWidth(0);*/
						/*box->Draw("SAME");*/
				/*}*/
				/*if(ha_v2->GetN()-1 >= j){*/
						/*box = new TBox(xh-xOffset,yh-yHerr,xh+xOffset,yh+yHerr);*/
						/*box->SetFillColor(17);*/
						/*box->SetFillStyle(1001);*/
						/*box->SetLineWidth(0);*/
                        /*box->Draw("SAME");*/
				/*}*/
				/*box = new TBox(xk-xOffset,yk-yKerr,xk+xOffset,yk+yKerr);*/
				/*box->SetFillColor(17);*/
				/*box->SetFillStyle(1001);*/
				/*box->SetLineWidth(0);*/
				/*box->Draw("SAME");*/

		/*}*/
		//
		// draw the graph with axis, continuous line, and put
		ks_v2->Draw( "PESAME" );
		la_v2->Draw( "PESAME" );
		//ha_v2->Draw( "PESAME" );
		ks8_v2->Draw("P");
		la8_v2->Draw( "P" );
		Int_t oldColor = 38;
		Int_t newColor = 46;

		TLegend* leg = new TLegend( 0.15,0.45,0.51,0.68 );
		leg->SetFillColor( 10 );
		leg->SetFillStyle( 0 );
		leg->SetBorderSize( 0.035 );
		leg->SetTextFont( 42 );
		leg->SetTextSize( 0.05 );
		leg->AddEntry( ks_v2,"#color[38]{K_{#lower[-0.3]{S}}#kern[-1.05]{#lower[0.1]{{}^{0}}}}", "P" );
		leg->AddEntry( la_v2, "#color[38]{#Lambda / #lower[0.1]{#bar{  }}#kern[-1.2]{#Lambda}}","P" );
		//leg->AddEntry( ha_v2, "#color[38]{h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}}", "P" );
		leg->AddEntry( ks8_v2, "#color[46]{K_{#lower[-0.3]{S}}#kern[-1.05]{#lower[0.1]{{}^{0}}}}", "P" );
		leg->AddEntry( la8_v2, "#color[46]{#Lambda / #lower[0.1]{#bar{  }}#kern[-1.2]{#Lambda}}", "P" );
		leg->Draw(  );

		TLatex *tex = new TLatex(  );
		tex->SetNDC(  );
		tex->SetTextFont( 42 );
		tex->SetTextSize( 0.05 );
		tex->DrawLatex( 0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = #color[38]{5.02} TeV, #color[46]{8.16} TeV" );
		tex->SetTextSize( 0.045 );
		//tex->DrawLatex( 0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}" );
		// tex->DrawLatex( 0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}" );
		tex->DrawLatex( 0.45,0.73,"185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
		//tex->DrawLatex( 0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}" );
		TLine* line = new TLine( 0,0,6,0 );
		line->SetLineStyle( 8 );
		line->Draw(  );
	}

