#include "TFile.h"
#include "TGraph.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TH2.h"
#include "interface/GetGraphFromFile.C"
#include "interface/MITStyle.C"

void V0CrossCheck_v2()
{
    MITStyle();
	bool Newdata = true;

	const int xi_npoints = 7;
    double v2Xi8[xi_npoints]  = {0.056368, 0.064876, 0.109999, 0.126307, 0.153165, 0.182298, 0.217702};
    double pTXi8[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055};
    double v2Xi8E[xi_npoints] = {0.017074, 0.008531, 0.007160, 0.005156, 0.004894, 0.006175, 0.010446};

	const int ks_npoints = 11;
    double v2Ks8[ks_npoints]  = {0.012813549, 0.028133917,0.043212, 0.05876162, 0.081288146 ,0.105507591 ,0.123968514 ,0.135179572 ,0.144205576 ,0.142730876 ,0.129697566};
    double pTKs8[ks_npoints]  = {0.3401, 0.5213, 0.7084, 0.9036, 1.201, 1.591, 1.986, 2.465, 3.136, 4.008, 5.071};
    double v2Ks8E[ks_npoints] = {0.00355244, 0.00110708,0.000619229, 0.000500387 ,0.000329161 ,0.000362995 ,0.000443819 ,0.000487576 ,0.000650721 ,0.00099496 ,0.00156744};

	const int la_npoints = 8;
    double v2La8[la_npoints]  = {0.031657 ,0.055453976 ,0.080209188 ,0.105853753 ,0.137719142 ,0.170879559 ,0.190002697 ,0.190002697};
    double pTLa8[la_npoints]  = {0.9145, 1.221, 1.605, 1.997, 2.481, 3.149, 4.004, 5.041};
    double v2La8E[la_npoints] = {0.002330326 ,0.000977656 ,0.000808941 ,0.000789458 ,0.000724636 ,0.000850324 ,0.001268046 ,0.002144964};

    // Pull TGraph for Kshort and lambda

    TFile* file_pPbv2 = TFile::Open("lrgraphv2_v3_pPb_185-220.root");
    TFile* file_hadv2 = TFile::Open("lrgraphv2_v3_pPb_hadron_185-above.root");

    TGraphErrors* ks_v2 = (TGraphErrors*)file_pPbv2->Get("kshortv2true");
    TGraphErrors* la_v2 = (TGraphErrors*)file_pPbv2->Get("lambdav2true");
    //TGraphErrors* ha_v2 = (TGraphErrors*)file_hadv2->Get("hadronv2");
    TGraphErrors* ha_v2 = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033pPb_v2.txt"),1,28,1.2);

    TGraphErrors* xi8_v2  = new TGraphErrors(xi_npoints,pTXi8,v2Xi8,0,v2Xi8E);
    TGraphErrors* ks8_v2  = new TGraphErrors(ks_npoints,pTKs8,v2Ks8,0,v2Ks8E);
	TGraphErrors* la8_v2  = new TGraphErrors(la_npoints,pTLa8,v2La8,0,v2La8E);

    ks_v2  ->SetMarkerColor(kBlue-4);
    ks_v2  ->SetLineColor(kBlue-4);
    ks_v2  ->SetMarkerStyle(25);
    ks_v2  ->SetMarkerSize(1.4);
    la_v2  ->SetMarkerColor(kBlue-4);
    la_v2  ->SetMarkerStyle(26);
    la_v2  ->SetMarkerSize(1.3);
    la_v2  ->SetLineColor(kBlue);
    ha_v2  ->SetMarkerStyle(28);
    ha_v2  ->SetMarkerSize(1.3);
    ks8_v2 ->SetMarkerColor(kRed);
    ks8_v2 ->SetMarkerStyle(20);
    ks8_v2 ->SetMarkerSize(1.5);
    ks8_v2 ->SetLineColor(kRed);

    la8_v2->SetMarkerColor(kGreen+2);
    la8_v2->SetMarkerStyle(22);
    la8_v2->SetMarkerSize(1.5);
    la8_v2->SetLineColor(kRed);

    TCanvas* c1 = MakeCanvas("c1","Plot");
    c1->cd();
    c1->SetLeftMargin(0.12);

    // draw the frame using a histogram frame

    TH1F* frame = c1->DrawFrame(0,-0.01,6,0.31);
    gPad->SetTickx();
    gPad->SetTicky();
    frame->GetXaxis()->CenterTitle(1);
    frame->GetYaxis()->CenterTitle(1);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame->GetYaxis()->SetTitle("v#kern[-0.3]{_{2}}");
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->SetTitleOffset(1.2,"Y");
    frame->SetTitleOffset(1.2,"X");



    // Plot systematic errors for ks and la
    TBox* box;
    double syst_pPb  = 0.069;
    double syst_pPbH = 0.039;
    double xOffset   = 0.08;
    double percent1L = 999;
    double percentH  = 999;

    int n = ks_v2->GetN();
    for(int j=0; j<n; j++)
    {
        double xk=0,yk=0,xl=0,yl=0,xh=0,yh=0;
        if(la_v2->GetN()-1 >= j){
            la_v2->GetPoint(j,xl,yl);
            percent1L = syst_pPb;
        }
        if(ha_v2->GetN()-1 >= j)
        {
            ha_v2->GetPoint(j,xh,yh);
            percentH = syst_pPbH;
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
            //box->Draw("SAME");
        }
        box = new TBox(xk-xOffset,yk-yKerr,xk+xOffset,yk+yKerr);
        box->SetFillColor(17);
        box->SetFillStyle(1001);
        box->SetLineWidth(0);
        box->Draw("SAME");

    }
    //
    // draw the graph with axis, continuous line, and put
    ks_v2->Draw("PESAME");
    la_v2->Draw("PESAME");
    //ha_v2->Draw("PESAME");
    ks8_v2->Draw("P");
    la8_v2->Draw("P");
    Int_t oldColor = 38;
    Int_t newColor = 46;

    TLegend* leg = new TLegend(0.15,0.45,0.51,0.68);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->AddEntry(ks_v2,"#color[38]{K_{#lower[-0.3]{S}}#kern[-1.05]{#lower[0.1]{{}^{0}}}}", "P");
    leg->AddEntry(la_v2, "#color[38]{#Lambda / #lower[0.1]{#bar{  }}#kern[-1.2]{#Lambda}}","P");
    //leg->AddEntry(ha_v2, "#color[38]{h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}}", "P");
    leg->AddEntry(ks8_v2, "#color[46]{K_{#lower[-0.3]{S}}#kern[-1.05]{#lower[0.1]{{}^{0}}}}", "P");
    leg->AddEntry(la8_v2, "#color[46]{#Lambda / #lower[0.1]{#bar{  }}#kern[-1.2]{#Lambda}}", "P");
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.05);
    tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = #color[38]{5.02} TeV, #color[46]{8.16} TeV");
    tex->SetTextSize(0.045);
    //tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
    // tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
    tex->DrawLatex(0.55,0.23,"185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");
    TLine* line = new TLine(0,0,6,0);
    line->SetLineStyle(8);
    line->Draw();
}

void V0CrossCheck_v3()
{
    MITStyle();
	const int ks_npoints = 10;
    double v2Ks8[ks_npoints]  = {0.005070706 ,0.013541808 ,0.017894052 ,0.026917717 ,0.036227128 ,0.045597329 ,0.051020426 ,0.045881397 ,0.037791439 ,0.014758437};
    double pTKs8[ks_npoints]  = {0.5213, 0.7084, 0.9036, 1.201, 1.591, 1.986, 2.465, 3.136, 4.008, 5.071};
    double v2Ks8E[ks_npoints] = {0.003412278 ,0.002043 ,0.001651186 ,0.001085994 ,0.001197198 ,0.001463774 ,0.001608019 ,0.002147082 ,0.003280849 ,0.005174682};

	const int la_npoints = 8;
    double v2La8[la_npoints]  = {0.011261892 ,0.020795145 ,0.033555289 ,0.044550178 ,0.058043265 ,0.060781616 ,0.073580978 ,0.046895594};
    double pTLa8[la_npoints]  = {0.9145, 1.221, 1.605, 1.997, 2.481, 3.149, 4.004, 5.041};
    double v2La8E[la_npoints] = {0.007693167 ,0.003225583 ,0.00266863 ,0.002604339 ,0.002402621 ,0.002821647 ,0.004176484 ,0.007135051};

    // Pull TGraph for Kshort and lambda

    TFile* file_pPbv2 = TFile::Open("lrgraphv2_v3_pPb_185-220.root");
    TFile* file_hadv2 = TFile::Open("lrgraphv2_v3_pPb_hadron_185-above.root");

    TGraphErrors* ks_v2 = (TGraphErrors*)file_pPbv2->Get("kshortv3true");
    TGraphErrors* la_v2 = (TGraphErrors*)file_pPbv2->Get("lambdav3true");
    //TGraphErrors* ha_v2 = (TGraphErrors*)file_hadv2->Get("hadronv2");
    /*TGraphErrors* ha_v2 = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033pPb_v2.txt"),1,28,1.2);*/

    TGraphErrors* ks8_v2 = new TGraphErrors(ks_npoints,pTKs8,v2Ks8,0,v2Ks8E);
	TGraphErrors* la8_v2 = new TGraphErrors(la_npoints,pTLa8,v2La8,0,v2La8E);

		ks_v2   ->SetMarkerColor(kBlue-4);
		ks_v2   ->SetLineColor(kBlue-4);
		ks_v2   ->SetMarkerStyle(24);
		ks_v2   ->SetMarkerSize(1.4);
		la_v2   ->SetMarkerColor(kBlue-4);
		//la_v2 ->SetMarkerStyle(22);
		la_v2   ->SetMarkerStyle(26);
		la_v2   ->SetMarkerSize(1.3);
		la_v2   ->SetLineColor(kBlue);
		ks8_v2  ->SetMarkerColor(kRed);
		//xi_v2 ->SetMarkerColor(kTeal-2);
		//xi_v2 ->SetMarkerColor(kMagenta-4);
		//xi_v2 ->SetMarkerColor(kGreen+2);
		ks8_v2  ->SetMarkerStyle(20);
		ks8_v2  ->SetMarkerSize(1.5);
		ks8_v2  ->SetLineColor(kRed);

		la8_v2  ->SetMarkerColor(kGreen+2);
		//xi_v2 ->SetMarkerColor(kTeal-2);
		//xi_v2 ->SetMarkerColor(kMagenta-4);
		//xi_v2 ->SetMarkerColor(kGreen+2);
		la8_v2  ->SetMarkerStyle(22);
		la8_v2  ->SetMarkerSize(1.5);
		la8_v2  ->SetLineColor(kBlue-4);

		TCanvas* c1 = MakeCanvas("c1", "Plot");
		c1->cd();
		c1->SetLeftMargin(0.12);

		// draw the frame using a histogram frame

		TH1F* frame = c1->DrawFrame(0,-0.01,6,0.15);
		gPad->SetTickx();
		gPad->SetTicky();
		frame->GetXaxis()->CenterTitle(1);
		frame->GetYaxis()->CenterTitle(1);
		frame->GetXaxis()->SetTitleSize(0.05);
		frame->GetXaxis()->SetTitle("p_{T} (GeV)");
		frame->GetYaxis()->SetTitle("v#kern[-0.3]{_{3}}");
		frame->GetYaxis()->SetTitleSize(0.05);
		frame->SetTitleOffset(1.2,"Y");
		frame->SetTitleOffset(1.2,"X");



		// Plot systematic errors for ks and la
		/*TBox* box;*/
		/*double syst_pPb = 0.069;*/
		/*double syst_pPbH = 0.039;*/
		/*double xOffset = 0.08;*/

		/*int n = ks_v2->GetN();*/
		/*for(int j=0; j<n; j++)*/
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
		ks_v2->Draw("PESAME");
		la_v2->Draw("PESAME");
		//ha_v2->Draw("PESAME");
		ks8_v2->Draw("P");
		la8_v2->Draw("P");
		Int_t oldColor = 38;
		Int_t newColor = 46;

		TLegend* leg = new TLegend(0.15,0.45,0.51,0.68);
		leg->SetFillColor(10);
		leg->SetFillStyle(0);
		leg->SetBorderSize(0);
		leg->SetTextFont(42);
		leg->SetTextSize(0.05);
		leg->AddEntry(ks_v2,"#color[38]{K_{#lower[-0.3]{S}}#kern[-1.05]{#lower[0.1]{{}^{0}}}}", "P");
		leg->AddEntry(la_v2, "#color[38]{#Lambda / #lower[0.1]{#bar{  }}#kern[-1.2]{#Lambda}}","P");
		//leg->AddEntry(ha_v2, "#color[38]{h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}}", "P");
		leg->AddEntry(ks8_v2, "#color[46]{K_{#lower[-0.3]{S}}#kern[-1.05]{#lower[0.1]{{}^{0}}}}", "P");
		leg->AddEntry(la8_v2, "#color[46]{#Lambda / #lower[0.1]{#bar{  }}#kern[-1.2]{#Lambda}}", "P");
		leg->Draw();

		TLatex *tex = new TLatex();
		tex->SetNDC();
		tex->SetTextFont(42);
		tex->SetTextSize(0.05);
		tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = #color[38]{5.02} TeV, #color[46]{8.16} TeV");
		tex->SetTextSize(0.045);
		//tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
		// tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
		tex->DrawLatex(0.45,0.73,"185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
		//tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");
		TLine* line = new TLine(0,0,6,0);
		line->SetLineStyle(8);
		line->Draw();
}

void Cascade_v2()
{
    MITStyle();

	const int xi_npoints = 7;
    double v2Xi8[xi_npoints]  = {0.056368, 0.064876, 0.109999, 0.126307, 0.153165, 0.182298, 0.217702};
    double pTXi8[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055};
    double v2Xi8E[xi_npoints] = {0.017074, 0.008531, 0.007160, 0.005156, 0.004894, 0.006175, 0.010446};

	const int ks_npoints = 11;
    double v2Ks8[ks_npoints]  = {0.012813549, 0.028133917,0.043212, 0.05876162, 0.081288146 ,0.105507591 ,0.123968514 ,0.135179572 ,0.144205576 ,0.142730876 ,0.129697566};
    double pTKs8[ks_npoints]  = {0.3401, 0.5213, 0.7084, 0.9036, 1.201, 1.591, 1.986, 2.465, 3.136, 4.008, 5.071};
    double v2Ks8E[ks_npoints] = {0.00355244, 0.00110708,0.000619229, 0.000500387 ,0.000329161 ,0.000362995 ,0.000443819 ,0.000487576 ,0.000650721 ,0.00099496 ,0.00156744};

	const int la_npoints = 8;
    double v2La8[la_npoints]  = {0.031657 ,0.055453976 ,0.080209188 ,0.105853753 ,0.137719142 ,0.170879559 ,0.190002697 ,0.190002697};
    double pTLa8[la_npoints]  = {0.9145, 1.221, 1.605, 1.997, 2.481, 3.149, 4.004, 5.041};
    double v2La8E[la_npoints] = {0.002330326 ,0.000977656 ,0.000808941 ,0.000789458 ,0.000724636 ,0.000850324 ,0.001268046 ,0.002144964};

    // Pull TGraph for Kshort and lambda

    TFile* file_pPbv2 = TFile::Open("lrgraphv2_v3_pPb_185-220.root");
    TFile* file_hadv2 = TFile::Open("lrgraphv2_v3_pPb_hadron_185-above.root");

    TGraphErrors* ks_v2 = (TGraphErrors*)file_pPbv2->Get("kshortv2true");
    TGraphErrors* la_v2 = (TGraphErrors*)file_pPbv2->Get("lambdav2true");
    //TGraphErrors* ha_v2 = (TGraphErrors*)file_hadv2->Get("hadronv2");
    TGraphErrors* ha_v2 = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033pPb_v2.txt"),1,28,1.2);

    TGraphErrors* xi8_v2 = new TGraphErrors(xi_npoints,pTXi8,v2Xi8,0,v2Xi8E);
    TGraphErrors* ks8_v2 = new TGraphErrors(ks_npoints,pTKs8,v2Ks8,0,v2Ks8E);
	TGraphErrors* la8_v2 = new TGraphErrors(la_npoints,pTLa8,v2La8,0,v2La8E);

    ha_v2->SetMarkerStyle(28);
    ha_v2->SetMarkerSize(1.3);

    ks8_v2->SetMarkerColor(kRed);
    ks8_v2->SetMarkerStyle(20);
    ks8_v2->SetMarkerSize(1.5);
    ks8_v2->SetLineColor(kRed);

    xi8_v2->SetMarkerColor(kMagenta-4);
    xi8_v2->SetMarkerStyle(21);
    xi8_v2->SetMarkerSize(1.5);
    xi8_v2->SetLineColor(kMagenta-4);

    la8_v2->SetMarkerColor(kGreen+2);
    la8_v2->SetMarkerStyle(22);
    la8_v2->SetMarkerSize(1.5);
    la8_v2->SetLineColor(kGreen+2);

    TCanvas* c1 = MakeCanvas("c1", "Plot");
    c1->cd();
    c1->SetLeftMargin(0.12);

    // draw the frame using a histogram frame

    TH1F* frame = c1->DrawFrame(0,-0.01,6,0.31);
    gPad->SetTickx();
    gPad->SetTicky();
    frame->GetXaxis()->CenterTitle(1);
    frame->GetYaxis()->CenterTitle(1);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame->GetYaxis()->SetTitle("v#kern[-0.3]{_{2}}");
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->SetTitleOffset(1.2,"Y");
    frame->SetTitleOffset(1.2,"X");

    //ha_v2->Draw("PESAME");
    ks8_v2->Draw("P");
    la8_v2->Draw("P");
    xi8_v2->Draw("P");

    // Write points into rootfile
    TFile out("8TeVgraphv2_pPb185-220.root","RECREATE");
    ks8_v2->Write("kshortv2");
    la8_v2->Write("lambdav2");
    xi8_v2->Write("cascadev2");

    TLegend* leg = new TLegend(0.15,0.45,0.51,0.68);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    //leg->AddEntry(ha_v2, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}", "P");
    leg->AddEntry(ks8_v2, "K_{#lower[-0.3]{S}}#kern[-1.05]{#lower[0.1]{{}^{0}}}", "P");
    leg->AddEntry(la8_v2, "#Lambda / #lower[0.1]{#bar{  }}#kern[-1.2]{#Lambda}", "P");
    leg->AddEntry(xi8_v2, "#Xi#kern[-0.3]{#lower[0.1]{{}^{+}}}/ #Xi#kern[-0.3]{#lower[0.1]{{}^{-}}}", "P");
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.05);
    tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->SetTextSize(0.045);
    //tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
    // tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
    tex->DrawLatex(0.55,0.23,"185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");
    TLine* line = new TLine(0,0,6,0);
    line->SetLineStyle(8);
    line->Draw();
}

void Rap_v2sig()
{
    MITStyle();

    const int xi_npoints = 8;
    double v2Xi8[xi_npoints]  = {0.0371521 ,0.0669904 ,0.0877339 ,0.123959 ,0.163467 ,0.183918 ,0.209974 ,0.17582};// ,0.178101 ,0.244044};
    /*double v2Xi8[xi_npoints]  = {0.012384 ,0.0223301 ,0.0292446 ,0.0413197 ,0.0544889 ,0.0613061 ,0.0699914 ,0.0586068};// v2/nq*/

    double pTXi8[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.474, 8.047, 11.52};
    /*double pTXi8[xi_npoints]  = {0.515319 ,0.768532 ,1.07944 ,1.50413 ,2.12034 ,2.92705 ,3.98122 ,5.28277};//KET*/
    /*double pTXi8[xi_npoints]  = {0.171773 ,0.256177 ,0.359813 ,0.501378 ,0.70678 ,0.975684 ,1.32707 ,1.76092};//KET/nq*/

    double v2Xi8E[xi_npoints] = {0.00932224 ,0.00484125 ,0.0041527 ,0.00326934 ,0.00309064 ,0.00361334 ,0.00517473 ,0.0113186};// ,0.0179293 ,0.0515298};
    /*double v2Xi8E[xi_npoints] = {0.00310741 ,0.00161375 ,0.00138423 ,0.00108978 ,0.00103021 ,0.00120445 ,0.00172491 ,0.00377288};// v2/nq erorrs*/

	/*const int xi_npoints = 11;*/
    /*double v2Xi8[xi_npoints]  = {0.0272316 ,0.0735904 ,0.108707 ,0.119905 ,0.166593 ,0.186789 ,0.213248 ,0.190285 ,0.165164 ,0.288952 ,0.25353};*/
    /*double pTXi8[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.275, 7.726, 9.098, 11.59};*/
    /*double v2Xi8E[xi_npoints] = {0.0152834 ,0.00784469 ,0.00667958 ,0.00524615 ,0.00493004 ,0.00577176 ,0.00822074 ,0.018064 ,0.0323931 ,0.0626886 ,0.0868749};*/


	const int ks_npoints = 13;
    double v2Ks8[ks_npoints]  = {0.0135543 ,0.0293041 ,0.0433881 ,0.0594498 ,0.082389 ,0.106184 ,0.124363 ,0.137749 ,0.146475 ,0.145182 ,0.135266 ,0.125505 ,0.124708};// ,0.11911 ,0.133803 ,0.0390893};
    /*double v2Ks8[ks_npoints]  = {0.00677717 ,0.0146521 ,0.021694 ,0.0297249 ,0.0411945 ,0.0530921 ,0.0621816 ,0.0688744 ,0.0732373 ,0.0725908 ,0.0676328 ,0.0627523 ,0.062354};// ,0.11911 ,0.133803 ,0.0390893}; //v2/nq*/
    double pTKs8[ks_npoints]  = {0.3666, 0.5309, 0.711, 0.9046, 1.202, 1.591, 1.986, 2.465, 3.136, 4.008, 5.142, 6.431, 7.619, 9.142, 11.64, 16.86};
    /*double pTKs8[ks_npoints]  = {0.120493 ,0.230933 ,0.370907 ,0.535241 ,0.804172 ,1.17001 ,1.55031 ,2.01776 ,2.67887 ,3.54188 ,4.67169 ,5.95357 ,7.13835};//KET*/
    /*double pTKs8[ks_npoints]  = {0.0602465 ,0.115467 ,0.185454 ,0.26762 ,0.402086 ,0.585004 ,0.775155 ,1.00888 ,1.33943 ,1.77094 ,2.33584 ,2.97678 ,3.56917};//KET/nqhhhhhhhhhhhhhhhh*/
    /*double v2Ks8E[ks_npoints] = {0.00217719 ,0.000522063 ,0.000305164 ,0.000248705 ,0.000161989 ,0.000173887 ,0.000208535 ,0.000225897 ,0.000298544 ,0.000448324 ,0.000705022 ,0.00146462 ,0.00192261};// ,0.00315165 ,0.0042099 ,0.0308448};*/
    double v2Ks8E[ks_npoints] = {0.0010886 ,0.000261031 ,0.000152582 ,0.000124353 ,8.09944e-05 ,8.69433e-05 ,0.000104268 ,0.000112948 ,0.000149272 ,0.000224162 ,0.000352511 ,0.00073231 ,0.000961305};// ,0.00315165 ,0.0042099 ,0.0308448}; ///nq

	const int la_npoints = 10;
    /*double v2La8[la_npoints]  = {0.0345669 ,0.0510091 ,0.0771498 ,0.106513 ,0.138504 ,0.171975 ,0.19385 ,0.202724 ,0.19716 ,0.178541};// ,0.168714 ,0.148544 ,0.122694};*/
    double v2La8[la_npoints]  = {0.0115223 ,0.017003 ,0.0257166 ,0.0355043 ,0.0461679 ,0.0573251 ,0.0646166 ,0.0675745 ,0.0657199 ,0.0595136};// ,0.168714 ,0.148544 ,0.122694}; 
    /*double pTLa8[la_npoints]  = {0.9252, 1.224, 1.603, 1.995, 2.485, 3.156, 4.007, 5.116, 6.414, 7.588, 9.117, 11.56, 16.89};*/
    /*double pTLa8[la_npoints]  = {0.326395 ,0.537952 ,0.836813 ,1.16991 ,1.60772 ,2.23092 ,3.0437 ,4.1234 ,5.39811 ,6.55991};//KET*/
    double pTLa8[la_npoints]  = {0.108798 ,0.179317 ,0.278938 ,0.389971 ,0.535908 ,0.74364 ,1.01457 ,1.37447 ,1.79937 ,2.18664};//KET/nq

    /*double v2La8E[la_npoints] = {0.00106898 ,0.000437383 ,0.000387893 ,0.000396151 ,0.000360204 ,0.000398977 ,0.000561998 ,0.000940498 ,0.00223511 ,0.0032202};// ,0.0063431 ,0.00742138 ,0.0622323};*/
    double v2La8E[la_npoints] = {0.000356326 ,0.000145794 ,0.000129298 ,0.00013205 ,0.000120068 ,0.000132992 ,0.000187333 ,0.000313499 ,0.000745036 ,0.0010734};// ,0.0063431 ,0.00742138 ,0.0622323};



    // Pull TGraph for Kshort and lambda

    TFile* file_pPbv2 = TFile::Open("lrgraphv2_v3_pPb_185-220.root");
    TFile* file_hadv2 = TFile::Open("lrgraphv2_v3_pPb_hadron_185-above.root");

    TGraphErrors* ks_v2 = (TGraphErrors*)file_pPbv2->Get("kshortv2true");
    TGraphErrors* la_v2 = (TGraphErrors*)file_pPbv2->Get("lambdav2true");
    //TGraphErrors* ha_v2 = (TGraphErrors*)file_hadv2->Get("hadronv2");
    TGraphErrors* ha_v2 = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033pPb_v2.txt"),1,28,1.2);

    TGraphErrors* xi8_v2 = new TGraphErrors(xi_npoints,pTXi8,v2Xi8,0,v2Xi8E);
    TGraphErrors* ks8_v2 = new TGraphErrors(ks_npoints,pTKs8,v2Ks8,0,v2Ks8E);
	TGraphErrors* la8_v2 = new TGraphErrors(la_npoints,pTLa8,v2La8,0,v2La8E);

    ha_v2->SetMarkerStyle(28);
    ha_v2->SetMarkerSize(1.3);

    ks8_v2->SetMarkerColor(kRed);
    ks8_v2->SetMarkerStyle(20);
    ks8_v2->SetMarkerSize(1.5);
    ks8_v2->SetLineColor(kRed);

    xi8_v2->SetMarkerColor(kGreen+2);
    xi8_v2->SetMarkerStyle(21);
    xi8_v2->SetMarkerSize(1.5);
    xi8_v2->SetLineColor(kGreen+2);

    la8_v2->SetMarkerColor(kBlue-4);
    la8_v2->SetMarkerStyle(22);
    la8_v2->SetMarkerSize(1.5);
    la8_v2->SetLineColor(kBlue-4);

    TCanvas* c1 = MakeCanvas("c1", "Plot");
    c1->cd();
    /*c1->SetLogy();*/
    c1->SetLeftMargin(0.12);

    // draw the frame using a histogram frame

    TH1F* frame = c1->DrawFrame(0,-0.01,6,0.2);
    /*TH1F* frame = c1->DrawFrame(0,0.01,20,1);*/
    gPad->SetTickx();
    gPad->SetTicky();
    frame->GetXaxis()->CenterTitle(1);
    frame->GetYaxis()->CenterTitle(1);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitle("KET/n_{q} (GeV)");
    frame->GetYaxis()->SetTitle("v_{2}^{sig}/n_{q}");
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->SetTitleOffset(1.1,"Y");
    frame->SetTitleOffset(1.2,"X");

    //ha_v2->Draw("PESAME");
    ks8_v2->Draw("P");
    la8_v2->Draw("P");
    xi8_v2->Draw("P");

    // Write points into rootfile
    TFile out("8TeVv2sigGraphv2_pPb185-250FullStats.root","RECREATE");
    ks8_v2->Write("kshortv2");
    la8_v2->Write("lambdav2");
    xi8_v2->Write("cascadev2");

    TLegend* leg = new TLegend(0.15,0.55,0.27,0.75);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    //leg->AddEntry(ha_v2, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}", "P");
    leg->AddEntry(ks8_v2, "K_{S}^{0}", "P");
    leg->AddEntry(la8_v2, "#Lambda / #bar{#Lambda}", "P");
    /*leg->AddEntry(xi8_v2, "#Xi#kern[-0.3]{#lower[0.1]{{}^{+}}}/ #Xi#kern[-0.3]{#lower[0.1]{{}^{-}}}", "P");*/
    leg->AddEntry(xi8_v2, "#Xi^{+}/ #Xi^{-}", "P");
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(62);
    tex->SetTextSize(0.05);
    tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->SetTextSize(0.045);
    //tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
    // tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
    tex->SetTextFont(42);
    tex->DrawLatex(0.40,0.24,"185 #leq N_{trk}^{offline} < 250");
    /*tex->DrawLatex(0.15,0.74,"|y| < 1");*/
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");
    TLine* line = new TLine(0,0,6,0);
    line->SetLineStyle(2);
    line->Draw("same");

    c1->Print("v2SigRapidityDividednq.pdf");

}

void Rap_v2obs()
{
    MITStyle();

    const int xi_npoints = 9;
    double v2Xi8[xi_npoints]  = {0.028288 ,0.0734096 ,0.109459 ,0.119556 ,0.166368 ,0.187297 ,0.21356 ,0.192063 ,0.16625};
    double pTXi8[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.938, 11.31};
    double v2Xi8E[xi_npoints] = {0.0145154 ,0.00762373 ,0.00651275 ,0.00512493 ,0.00481561 ,0.00563771 ,0.00803405 ,0.0175756 ,0.0315233};

	/*const int xi_npoints = 11;*/
    /*double v2Xi8[xi_npoints]  = {0.0272316 ,0.0735904 ,0.108707 ,0.119905 ,0.166593 ,0.186789 ,0.213248 ,0.190285 ,0.165164 ,0.288952 ,0.25353};*/
    /*double pTXi8[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.275, 7.726, 9.098, 11.59};*/
    /*double v2Xi8E[xi_npoints] = {0.0152834 ,0.00784469 ,0.00667958 ,0.00524615 ,0.00493004 ,0.00577176 ,0.00822074 ,0.018064 ,0.0323931 ,0.0626886 ,0.0868749};*/


	const int ks_npoints = 16;
    double v2Ks8[ks_npoints]  = {0.0183679 ,0.0289959 ,0.0428532 ,0.0598327 ,0.082664 ,0.10635 ,0.124404 ,0.138096 ,0.146914 ,0.143748 ,0.135454 ,0.120565 ,0.12986 ,0.11892 ,0.134655 ,0.129064};
    double pTKs8[ks_npoints]  = {0.3666, 0.5309, 0.711, 0.9046, 1.202, 1.591, 1.986, 2.465, 3.136, 4.008, 5.142, 6.431, 7.619, 9.142, 11.64, 16.86};
    double v2Ks8E[ks_npoints] = {0.00375874 ,0.000902538 ,0.000527948 ,0.000430599 ,0.000281917 ,0.000303551 ,0.000363822 ,0.00039436 ,0.000519449 ,0.000777174 ,0.00122037 ,0.00253245 ,0.00329789 ,0.00543151 ,0.00622973 ,0.017103};

	const int la_npoints = 13;
    double v2La8[la_npoints]  = {0.0354139 ,0.0516292 ,0.077228 ,0.107637 ,0.138136 ,0.171592 ,0.193789 ,0.200912 ,0.200429 ,0.176394 ,0.180825 ,0.169005 ,0.202348};
    double pTLa8[la_npoints]  = {0.9252, 1.224, 1.603, 1.995, 2.485, 3.156, 4.007, 5.116, 6.414, 7.588, 9.117, 11.56, 16.89};
    double v2La8E[la_npoints] = {0.00184566 ,0.000756422 ,0.000671435 ,0.000686081 ,0.000624802 ,0.000692926 ,0.000974103 ,0.00162467 ,0.00386568 ,0.00554866 ,0.0101038 ,0.0126748 ,0.0355234};



    // Pull TGraph for Kshort and lambda

    //TGraphErrors* ha_v2 = (TGraphErrors*)file_hadv2->Get("hadronv2");
    TGraphErrors* ha_v2 = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033pPb_v2.txt"),1,28,1.2);

    TGraphErrors* xi8_v2 = new TGraphErrors(xi_npoints,pTXi8,v2Xi8,0,v2Xi8E);
    TGraphErrors* ks8_v2 = new TGraphErrors(ks_npoints,pTKs8,v2Ks8,0,v2Ks8E);
	TGraphErrors* la8_v2 = new TGraphErrors(la_npoints,pTLa8,v2La8,0,v2La8E);

    ha_v2->SetMarkerStyle(28);
    ha_v2->SetMarkerSize(1.3);

    ks8_v2->SetMarkerColor(kRed);
    ks8_v2->SetMarkerStyle(20);
    ks8_v2->SetMarkerSize(1.5);
    ks8_v2->SetLineColor(kRed);

    xi8_v2->SetMarkerColor(kGreen+2);
    xi8_v2->SetMarkerStyle(21);
    xi8_v2->SetMarkerSize(1.5);
    xi8_v2->SetLineColor(kGreen+2);

    la8_v2->SetMarkerColor(kBlue-4);
    la8_v2->SetMarkerStyle(22);
    la8_v2->SetMarkerSize(1.5);
    la8_v2->SetLineColor(kBlue-4);

    TCanvas* c1 = MakeCanvas("c1", "Plot");
    c1->cd();
    /*c1->SetLogy();*/
    c1->SetLeftMargin(0.12);

    // draw the frame using a histogram frame

    TH1F* frame = c1->DrawFrame(0,-0.01,20,0.45);
    /*TH1F* frame = c1->DrawFrame(0,0.01,20,1);*/
    gPad->SetTickx();
    gPad->SetTicky();
    frame->GetXaxis()->CenterTitle(1);
    frame->GetYaxis()->CenterTitle(1);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame->GetYaxis()->SetTitle("v_{2}^{obs}");
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->SetTitleOffset(1.1,"Y");
    frame->SetTitleOffset(1.2,"X");

    //ha_v2->Draw("PESAME");
    ks8_v2->Draw("P");
    la8_v2->Draw("P");
    xi8_v2->Draw("P");

    // Write points into rootfile
    TFile out("8TeVv2obsGraphv2_pPb185-250.root","RECREATE");
    ks8_v2->Write("kshortv2");
    la8_v2->Write("lambdav2");
    xi8_v2->Write("cascadev2");

    TLegend* leg = new TLegend(0.15,0.55,0.27,0.75);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    //leg->AddEntry(ha_v2, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}", "P");
    leg->AddEntry(ks8_v2, "K_{S}^{0}", "P");
    leg->AddEntry(la8_v2, "#Lambda / #bar{#Lambda}", "P");
    /*leg->AddEntry(xi8_v2, "#Xi#kern[-0.3]{#lower[0.1]{{}^{+}}}/ #Xi#kern[-0.3]{#lower[0.1]{{}^{-}}}", "P");*/
    leg->AddEntry(xi8_v2, "#Xi^{+}/ #Xi^{-}", "P");
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(62);
    tex->SetTextSize(0.05);
    tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->SetTextSize(0.045);
    //tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
    // tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
    tex->SetTextFont(42);
    tex->DrawLatex(0.40,0.24,"185 #leq N_{trk}^{offline} < 250");
    /*tex->DrawLatex(0.15,0.74,"|y| < 1");*/
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");
    TLine* line = new TLine(0,0,6,0);
    line->SetLineStyle(2);
    line->Draw("same");

    c1->Print("v2ObsRapidity.pdf");

}

void Rap_v2bkg()
{
    MITStyle();

    const int xi_npoints = 9;
    double v2Xi8[xi_npoints]  = {0.0502074 ,0.0666694 ,0.140997 ,0.103829 ,0.156365 ,0.209947 ,0.227812 ,0.260397 ,0.2087};
    double pTXi8[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.938, 11.31};
    double v2Xi8E[xi_npoints] = {0.0168844 ,0.0209328 ,0.0128893 ,0.014907 ,0.0126887 ,0.0123909 ,0.0188829 ,0.030211 ,0.0632523};

	/*const int xi_npoints = 11;*/
    /*double v2Xi8[xi_npoints]  = {0.0272316 ,0.0735904 ,0.108707 ,0.119905 ,0.166593 ,0.186789 ,0.213248 ,0.190285 ,0.165164 ,0.288952 ,0.25353};*/
    /*double pTXi8[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.275, 7.726, 9.098, 11.59};*/
    /*double v2Xi8E[xi_npoints] = {0.0152834 ,0.00784469 ,0.00667958 ,0.00524615 ,0.00493004 ,0.00577176 ,0.00822074 ,0.018064 ,0.0323931 ,0.0626886 ,0.0868749};*/


	const int ks_npoints = 16;
    double v2Ks8[ks_npoints]  = {-0.0276159 ,0.0505831 ,0.0633616 ,0.0659727 ,0.0935834 ,0.116099 ,0.132347 ,0.148838 ,0.163849 ,0.169637 ,0.176803 ,0.158785 ,0.167225 ,0.129859 ,0.155915 ,0.21973};
    double pTKs8[ks_npoints]  = {0.3666, 0.5309, 0.711, 0.9046, 1.202, 1.591, 1.986, 2.465, 3.136, 4.008, 5.142, 6.431, 7.619, 9.142, 11.64, 16.86};
    double v2Ks8E[ks_npoints] = {0.0148551 ,0.00407726 ,0.00286001 ,0.00296167 ,0.00170357 ,0.00164754 ,0.00185762 ,0.00176369 ,0.00190469 ,0.00258989 ,0.0033231 ,0.00672046 ,0.00846984 ,0.0149719 ,0.0140114 ,0.021461};

	const int la_npoints = 13;
    double v2La8[la_npoints]  = {0.0837146 ,0.0974558 ,0.118958 ,0.13642 ,0.1587 ,0.176193 ,0.189751 ,0.192593 ,0.17456 ,0.15699 ,0.188704 ,0.207833 ,0.217665};
    double pTLa8[la_npoints]  = {0.9252, 1.224, 1.603, 1.995, 2.485, 3.156, 4.007, 5.116, 6.414, 7.588, 9.117, 11.56, 16.89};
    double v2La8E[la_npoints] = {0.000934595 ,0.000809648 ,0.00102125 ,0.00135576 ,0.00144641 ,0.00189797 ,0.00280107 ,0.00447446 ,0.010036 ,0.0129962 ,0.0179161 ,0.0155427 ,0.0443472};



    // Pull TGraph for Kshort and lambda

    TFile* file_pPbv2 = TFile::Open("lrgraphv2_v3_pPb_185-220.root");
    TFile* file_hadv2 = TFile::Open("lrgraphv2_v3_pPb_hadron_185-above.root");

    TGraphErrors* ks_v2 = (TGraphErrors*)file_pPbv2->Get("kshortv2true");
    TGraphErrors* la_v2 = (TGraphErrors*)file_pPbv2->Get("lambdav2true");
    //TGraphErrors* ha_v2 = (TGraphErrors*)file_hadv2->Get("hadronv2");
    TGraphErrors* ha_v2 = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033pPb_v2.txt"),1,28,1.2);

    TGraphErrors* xi8_v2 = new TGraphErrors(xi_npoints,pTXi8,v2Xi8,0,v2Xi8E);
    TGraphErrors* ks8_v2 = new TGraphErrors(ks_npoints,pTKs8,v2Ks8,0,v2Ks8E);
	TGraphErrors* la8_v2 = new TGraphErrors(la_npoints,pTLa8,v2La8,0,v2La8E);

    ha_v2->SetMarkerStyle(28);
    ha_v2->SetMarkerSize(1.3);

    ks8_v2->SetMarkerColor(kRed);
    ks8_v2->SetMarkerStyle(20);
    ks8_v2->SetMarkerSize(1.5);
    ks8_v2->SetLineColor(kRed);

    xi8_v2->SetMarkerColor(kGreen+2);
    xi8_v2->SetMarkerStyle(21);
    xi8_v2->SetMarkerSize(1.5);
    xi8_v2->SetLineColor(kGreen+2);

    la8_v2->SetMarkerColor(kBlue-4);
    la8_v2->SetMarkerStyle(22);
    la8_v2->SetMarkerSize(1.5);
    la8_v2->SetLineColor(kBlue-4);

    TCanvas* c1 = MakeCanvas("c1", "Plot");
    c1->cd();
    /*c1->SetLogy();*/
    c1->SetLeftMargin(0.12);

    // draw the frame using a histogram frame

    TH1F* frame = c1->DrawFrame(0,-0.05,20,0.45);
    /*TH1F* frame = c1->DrawFrame(0,0.01,20,1);*/
    gPad->SetTickx();
    gPad->SetTicky();
    frame->GetXaxis()->CenterTitle(1);
    frame->GetYaxis()->CenterTitle(1);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame->GetYaxis()->SetTitle("v_{2}^{bkg}");
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->SetTitleOffset(1.1,"Y");
    frame->SetTitleOffset(1.2,"X");

    //ha_v2->Draw("PESAME");
    ks8_v2->Draw("P");
    la8_v2->Draw("P");
    xi8_v2->Draw("P");

    // Write points into rootfile
    TFile out("8TeVv2bkgGraphv2_pPb185-250.root","RECREATE");
    ks8_v2->Write("kshortv2bkg");
    la8_v2->Write("lambdav2bkg");
    xi8_v2->Write("cascadev2bkg");

    TLegend* leg = new TLegend(0.15,0.55,0.27,0.75);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    //leg->AddEntry(ha_v2, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}", "P");
    leg->AddEntry(ks8_v2, "K_{S}^{0}", "P");
    leg->AddEntry(la8_v2, "#Lambda / #bar{#Lambda}", "P");
    /*leg->AddEntry(xi8_v2, "#Xi#kern[-0.3]{#lower[0.1]{{}^{+}}}/ #Xi#kern[-0.3]{#lower[0.1]{{}^{-}}}", "P");*/
    leg->AddEntry(xi8_v2, "#Xi^{+}/ #Xi^{-}", "P");
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(62);
    tex->SetTextSize(0.05);
    tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->SetTextSize(0.045);
    //tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
    // tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
    tex->SetTextFont(42);
    tex->DrawLatex(0.40,0.24,"185 #leq N_{trk}^{offline} < 250");
    /*tex->DrawLatex(0.15,0.74,"|y| < 1");*/
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");
    TLine* line = new TLine(0,0,6,0);
    line->SetLineStyle(2);
    line->Draw("same");

    c1->Print("v2BkgRapidity.pdf");

}

void Rap_fsig()
{
    MITStyle();

    const int xi_npoints = 9;
    double fsig_Xi8[xi_npoints]  = {0.954019 ,0.973881 ,0.976705 ,0.97829 ,0.978074 ,0.978057 ,0.978603 ,0.974693 ,0.976012};
    double pTXi8[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.938, 11.31};

	const int ks_npoints = 16;
    double fsig_Ks8[ks_npoints]  = {0.999476 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,0.999999 ,0.999999 ,0.999988 ,0.999997 ,0.992327 ,0.99836 ,0.890106,0.481433};
    double pTKs8[ks_npoints]  = {0.3666, 0.5309, 0.711, 0.9046, 1.202, 1.591, 1.986, 2.465, 3.136, 4.008, 5.142, 6.431, 7.619, 9.142, 11.64, 16.86};

	const int la_npoints = 13;
    double fsig_La8[la_npoints]  = {0.99882 ,0.999987 ,1 ,0.999524 ,0.999632 ,0.999855 ,0.999698 ,0.998783 ,0.999771 ,0.997088 ,0.92718 ,0.990913,0.421011};
    double pTLa8[la_npoints]  = {0.9252, 1.224, 1.603, 1.995, 2.485, 3.156, 4.007, 5.116, 6.414, 7.588, 9.117, 11.56, 16.89};

    // Pull TGraph for Kshort and lambda

    TGraphErrors* xi8_fsig_ = new TGraphErrors(xi_npoints,pTXi8,fsig_Xi8,0,0);
    TGraphErrors* ks8_fsig_ = new TGraphErrors(ks_npoints,pTKs8,fsig_Ks8,0,0);
	TGraphErrors* la8_fsig_ = new TGraphErrors(la_npoints,pTLa8,fsig_La8,0,0);

    ks8_fsig_->SetMarkerColor(kRed);
    ks8_fsig_->SetMarkerStyle(20);
    ks8_fsig_->SetMarkerSize(1.5);
    ks8_fsig_->SetLineColor(kRed);

    xi8_fsig_->SetMarkerColor(kGreen+2);
    xi8_fsig_->SetMarkerStyle(21);
    xi8_fsig_->SetMarkerSize(1.5);
    xi8_fsig_->SetLineColor(kGreen+2);

    la8_fsig_->SetMarkerColor(kBlue-4);
    la8_fsig_->SetMarkerStyle(22);
    la8_fsig_->SetMarkerSize(1.5);
    la8_fsig_->SetLineColor(kBlue-4);

    TCanvas* c1 = MakeCanvas("c1", "Plot");
    c1->cd();
    /*c1->SetLogy();*/
    c1->SetLeftMargin(0.12);

    // draw the frame using a histogram frame

    TH1F* frame = c1->DrawFrame(0,0.2,20,1.5);
    /*TH1F* frame = c1->DrawFrame(0,0.01,20,1);*/
    gPad->SetTickx();
    gPad->SetTicky();
    frame->GetXaxis()->CenterTitle(1);
    frame->GetYaxis()->CenterTitle(1);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame->GetYaxis()->SetTitle("f_{sig}");
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->SetTitleOffset(1.1,"Y");
    frame->SetTitleOffset(1.2,"X");

    //ha_fsig_->Draw("PESAME");
    ks8_fsig_->Draw("P");
    la8_fsig_->Draw("P");
    xi8_fsig_->Draw("P");

    // Write points into rootfile
    TFile out("8TeVfsig_GraphpPb185-250.root","RECREATE");
    ks8_fsig_->Write("kshortfsig");
    la8_fsig_->Write("lambdafsig");
    xi8_fsig_->Write("cascadefsig");

    TLegend* leg = new TLegend(0.70,0.65,0.90,0.85);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    //leg->AddEntry(ha_v2, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}", "P");
    leg->AddEntry(ks8_fsig_, "K_{S}^{0}", "P");
    leg->AddEntry(la8_fsig_, "#Lambda / #bar{#Lambda}", "P");
    /*leg->AddEntry(xi8_v2, "#Xi#kern[-0.3]{#lower[0.1]{{}^{+}}}/ #Xi#kern[-0.3]{#lower[0.1]{{}^{-}}}", "P");*/
    leg->AddEntry(xi8_fsig_, "#Xi^{+}/ #Xi^{-}", "P");
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(62);
    tex->SetTextSize(0.045);
    tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->SetTextSize(0.04);
    //tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
    // tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
    tex->SetTextFont(42);
    tex->DrawLatex(0.40,0.24,"185 #leq N_{trk}^{offline} < 250");
    /*tex->DrawLatex(0.15,0.74,"|y| < 1");*/
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");
    TLine* line = new TLine(0,1,20,1);
    line->SetLineStyle(2);
    line->Draw("same");

    c1->Print("fsig_Rapidity.pdf");

}

void RapSys_RecoCuts()
{
    MITStyle();
    gStyle->SetTitleAlign(33);
    const int ks_npoints = 16;
    double v2loose_ks[ks_npoints]  = {0.0103546 ,0.0302139 ,0.0433535 ,0.0596093 ,0.082594 ,0.106931 ,0.124934 ,0.138023 ,0.147092 ,0.146583 ,0.135002 ,0.122919 ,0.127094 ,0.110311 ,0.134414 ,-0.0317351};
    double pTloose_ks[ks_npoints]  = {0.3666, 0.5309, 0.711, 0.9046, 1.202, 1.591, 1.986, 2.465, 3.136, 4.008, 5.142, 6.431, 7.619, 9.142, 11.64, 16.86};
    double v2loose_ksE[ks_npoints] = {0.00362787 ,0.000865776 ,0.000504115 ,0.000409552 ,0.000265851 ,0.000284811 ,0.00034122 ,0.000369488 ,0.000488002 ,0.000731881 ,0.00114944 ,0.00238476 ,0.00312644 ,0.00512804 ,0.00682114 ,0.0492771};

    double v2standard_ks[ks_npoints]  = {0.0135543 ,0.0293041 ,0.0433881 ,0.0594498 ,0.082389 ,0.106184 ,0.124363 ,0.137749 ,0.146475 ,0.145182 ,0.135266 ,0.125505 ,0.124708 ,0.11911 ,0.133803 ,0.0390893};
    double pTstandard_ks[ks_npoints]  = {0.3666, 0.5309, 0.711, 0.9046, 1.202, 1.591, 1.986, 2.465, 3.136, 4.008, 5.142, 6.431, 7.619, 9.142, 11.64, 16.86};
    double v2standard_ksE[ks_npoints] = {0.00217719 ,0.000522063 ,0.000305164 ,0.000248705 ,0.000161989 ,0.000173887 ,0.000208535 ,0.000225897 ,0.000298544 ,0.000448324 ,0.000705022 ,0.00146462 ,0.00192261 ,0.00315165 ,0.0042099 ,0.0308448};

    double v2tight_ks[ks_npoints]  = {0.00784946 ,0.0303067 ,0.0438613 ,0.0596843 ,0.0831109 ,0.107091 ,0.124913 ,0.137904 ,0.147068 ,0.146607 ,0.134644 ,0.122733 ,0.126291 ,0.109094 ,0.131269 ,0.077138};
    double pTtight_ks[ks_npoints]  = {0.3666, 0.5309, 0.711, 0.9046, 1.202, 1.591, 1.986, 2.465, 3.136, 4.008, 5.142, 6.431, 7.619, 9.142, 11.64, 16.86};
    double v2tight_ksE[ks_npoints] = {0.00646683 ,0.00128322 ,0.000643845 ,0.000482569 ,0.000295253 ,0.000304724 ,0.000358735 ,0.000384387 ,0.000504543 ,0.000755624 ,0.00118735 ,0.00246939 ,0.00324637 ,0.00534846 ,0.00762112 ,0.07372};


    // Pull TGraph for Kshort and lambda

    TGraphErrors* loose_v2_ks = new TGraphErrors(ks_npoints,pTloose_ks,v2loose_ks,0,v2loose_ksE);
    TGraphErrors* standard_v2_ks = new TGraphErrors(ks_npoints,pTstandard_ks,v2standard_ks,0,v2standard_ksE);
	TGraphErrors* tight_v2_ks = new TGraphErrors(ks_npoints,pTtight_ks,v2tight_ks,0,v2tight_ksE);

    standard_v2_ks->SetMarkerColor(kRed);
    standard_v2_ks->SetMarkerStyle(20);
    standard_v2_ks->SetMarkerSize(1.5);
    standard_v2_ks->SetLineColor(kRed);

    loose_v2_ks->SetMarkerColor(kRed);
    loose_v2_ks->SetMarkerStyle(25);
    loose_v2_ks->SetMarkerSize(1.5);
    loose_v2_ks->SetLineColor(kRed);

    tight_v2_ks->SetMarkerColor(kRed);
    tight_v2_ks->SetMarkerStyle(26);
    tight_v2_ks->SetMarkerSize(1.5);
    tight_v2_ks->SetLineColor(kRed);

    TCanvas* c1_ks = MakeCanvas("c1_ks", "Plot_ks");
    c1_ks->cd();
    /*c1_ks->SetLogy();*/
    c1_ks->SetLeftMargin(0.12);

    // draw the frame_ks using a histogram frame_ks

    TH1F* frame_ks = c1_ks->DrawFrame(0,-0.08,20,0.45);
    /*TH1F* frame_ks = c1_ks->DrawFrame(0,0.01,20,1);*/
    gPad->SetTickx();
    gPad->SetTicky();
    frame_ks->SetTitle("K_{S}^{0} Reconstruction Cuts");
    frame_ks->SetTitleSize(0.055,"t");
    frame_ks->GetXaxis()->CenterTitle(1);
    frame_ks->GetYaxis()->CenterTitle(1);
    frame_ks->GetXaxis()->SetTitleSize(0.05);
    frame_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_ks->GetYaxis()->SetTitle("v_{2}^{sig}");
    frame_ks->GetYaxis()->SetTitleSize(0.05);
    frame_ks->SetTitleOffset(1.1,"Y");
    frame_ks->SetTitleOffset(1.2,"X");

    //ha_v2->Draw("PESAME");
    standard_v2_ks->Draw("P");
    tight_v2_ks->Draw("P");
    loose_v2_ks->Draw("P");

    // Write points into rootfile

    TLegend* leg_ks = new TLegend(0.15,0.55,0.27,0.75);
    leg_ks->SetFillColor(10);
    leg_ks->SetFillStyle(0);
    leg_ks->SetBorderSize(0);
    leg_ks->SetTextFont(42);
    leg_ks->SetTextSize(0.03);
    leg_ks->AddEntry(standard_v2_ks, "Standard reconstruction", "P");
    leg_ks->AddEntry(tight_v2_ks, "Tight reconstruction", "P");
    leg_ks->AddEntry(loose_v2_ks, "Loose reconstruction", "P");
    leg_ks->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(62);
    tex->SetTextSize(0.045);
    tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->SetTextSize(0.04);
    tex->SetTextFont(42);
    tex->DrawLatex(0.40,0.32,"185 #leq N_{trk}^{offline} < 250");
    TLine* line = new TLine(0,0,6,0);
    line->SetLineStyle(2);
    line->Draw("same");

    c1_ks->Print("v2KshortRecoSystematics.pdf");

    //Do lambda
    const int la_npoints = 13;
    double v2loose_la[la_npoints]  = {0.0364876 ,0.0503007 ,0.0766704 ,0.106931 ,0.138954 ,0.172568 ,0.19213 ,0.201382 ,0.193664 ,0.174807 ,0.178977 ,0.14746 ,-0.0336306};
    double pTloose_la[la_npoints]  = {0.9252, 1.224, 1.603, 1.995, 2.485, 3.156, 4.007, 5.116, 6.414, 7.588, 9.117, 11.56, 16.89};
    double v2loose_laE[la_npoints] = {0.00177145 ,0.000720517 ,0.000636639 ,0.000648237 ,0.000588564 ,0.0006516 ,0.000916771 ,0.00153283 ,0.00362961 ,0.00522624 ,0.0102896 ,0.0119938 ,0.0987539};

    double v2standard_la[la_npoints]  = {0.0345669 ,0.0510091 ,0.0771498 ,0.106513 ,0.138504 ,0.171975 ,0.19385 ,0.202724 ,0.19716 ,0.178541 ,0.168714 ,0.148544 ,0.122694};
    double pTstandard_la[la_npoints]  = {0.9252, 1.224, 1.603, 1.995, 2.485, 3.156, 4.007, 5.116, 6.414, 7.588, 9.117, 11.56, 16.89};
    double v2standard_laE[la_npoints] = {0.00106898 ,0.000437383 ,0.000387893 ,0.000396151 ,0.000360204 ,0.000398977 ,0.000561998 ,0.000940498 ,0.00223511 ,0.0032202 ,0.0063431 ,0.00742138 ,0.0622323};

    double v2tight_la[la_npoints]  = {0.0312339 ,0.0484565 ,0.0762154 ,0.106016 ,0.138673 ,0.17217 ,0.192093 ,0.201094 ,0.193343 ,0.172166 ,0.178001 ,0.150506 ,-0.120121};
    double pTtight_la[la_npoints]  = {0.9252, 1.224, 1.603, 1.995, 2.485, 3.156, 4.007, 5.116, 6.414, 7.588, 9.117, 11.56, 16.89};
    double v2tight_laE[la_npoints] = {0.00224476 ,0.000827416 ,0.000706016 ,0.000709656 ,0.000637924 ,0.000699588 ,0.000977282 ,0.00162535 ,0.003845 ,0.00555049 ,0.0111427 ,0.0131359 ,0.135443};


    // Pull TGraph for Kshort and lambda

    TGraphErrors* loose_v2_la = new TGraphErrors(la_npoints,pTloose_la,v2loose_la,0,v2loose_laE);
    TGraphErrors* standard_v2_la = new TGraphErrors(la_npoints,pTstandard_la,v2standard_la,0,v2standard_laE);
	TGraphErrors* tight_v2_la = new TGraphErrors(la_npoints,pTtight_la,v2tight_la,0,v2tight_laE);

    standard_v2_la->SetMarkerColor(kBlue-4);
    standard_v2_la->SetMarkerStyle(20);
    standard_v2_la->SetMarkerSize(1.5);
    standard_v2_la->SetLineColor(kBlue-4);

    loose_v2_la->SetMarkerColor(kBlue-4);
    loose_v2_la->SetMarkerStyle(25);
    loose_v2_la->SetMarkerSize(1.5);
    loose_v2_la->SetLineColor(kBlue-4);

    tight_v2_la->SetMarkerColor(kBlue-4);
    tight_v2_la->SetMarkerStyle(26);
    tight_v2_la->SetMarkerSize(1.5);
    tight_v2_la->SetLineColor(kBlue-4);

    TCanvas* c1_la = MakeCanvas("c1_la", "Plot_la");
    c1_la->cd();
    /*c1_la->SetLogy();*/
    c1_la->SetLeftMargin(0.12);

    // draw the frame_la using a histogram frame_la

    TH1F* frame_la = c1_la->DrawFrame(0,-0.15,20,0.60);
    /*TH1F* frame_la = c1_la->DrawFrame(0,0.01,20,1);*/
    gPad->SetTickx();
    gPad->SetTicky();
    /*frame_la->SetTitle("#Lambda Reconstruction Cuts");*/
    frame_la->GetXaxis()->CenterTitle(1);
    frame_la->GetYaxis()->CenterTitle(1);
    frame_la->GetXaxis()->SetTitleSize(0.05);
    frame_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_la->GetYaxis()->SetTitle("v_{2}^{sig}");
    frame_la->GetYaxis()->SetTitleSize(0.05);
    frame_la->SetTitleOffset(1.1,"Y");
    frame_la->SetTitleOffset(1.2,"X");

    //ha_v2->Draw("PESAME");
    standard_v2_la->Draw("P");
    tight_v2_la->Draw("P");
    loose_v2_la->Draw("P");

    TLegend* leg_la = new TLegend(0.15,0.55,0.27,0.75);
    leg_la->SetFillColor(10);
    leg_la->SetFillStyle(0);
    leg_la->SetBorderSize(0);
    leg_la->SetTextFont(42);
    leg_la->SetTextSize(0.03);
    leg_la->AddEntry(standard_v2_la, "Standard reconstruction", "P");
    leg_la->AddEntry(tight_v2_la, "Tight reconstruction", "P");
    leg_la->AddEntry(loose_v2_la, "Loose reconstruction", "P");
    leg_la->Draw();

    tex->SetTextFont(62);
    tex->SetTextSize(0.04);
    tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->SetTextSize(0.04);
    tex->SetTextFont(42);
    tex->DrawLatex(0.40,0.32,"185 #leq N_{trk}^{offline} < 250");
    line->Draw("same");

    c1_la->Print("v2LambdaRecoSystematics.pdf");

    //Combined
    TCanvas* c1_co = MakeCanvas("c1_co", "Plot_co");
    c1_co->cd();
    /*c1_co->SetLogy();*/
    c1_co->SetLeftMargin(0.12);

    // draw the frame_co using a histogram frame_co

    TH1F* frame_co = c1_co->DrawFrame(0,-0.15,20,0.60);
    /*TH1F* frame_co = c1_co->DrawFrame(0,0.01,20,1);*/
    gPad->SetTickx();
    gPad->SetTicky();
    frame_co->GetXaxis()->CenterTitle(1);
    frame_co->GetYaxis()->CenterTitle(1);
    frame_co->GetXaxis()->SetTitleSize(0.05);
    frame_co->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_co->GetYaxis()->SetTitle("v_{2}^{sig}");
    frame_co->GetYaxis()->SetTitleSize(0.05);
    frame_co->SetTitleOffset(1.1,"Y");
    frame_co->SetTitleOffset(1.2,"X");

    standard_v2_ks->Draw("P");
    loose_v2_ks->Draw("P");
    tight_v2_ks->Draw("P");

    standard_v2_la->Draw("P");
    loose_v2_la->Draw("P");
    tight_v2_la->Draw("P");

    TLegend* leg_co1 = new TLegend(0.15,0.55,0.27,0.75);
    leg_co1->SetFillColor(10);
    leg_co1->SetFillStyle(0);
    leg_co1->SetBorderSize(0);
    leg_co1->SetTextFont(42);
    leg_co1->SetTextSize(0.04);
    leg_co1->AddEntry(standard_v2_ks, "K_{S}^{0}", "P");
    leg_co1->AddEntry(standard_v2_la, "#Lambda / #bar{#Lambda}", "P");
    leg_co1->Draw();

    TLegend* leg_co2 = new TLegend(0.6,0.65,0.85,0.85);
    leg_co2->SetFillColor(10);
    leg_co2->SetFillStyle(0);
    leg_co2->SetBorderSize(0);
    leg_co2->SetTextFont(42);
    leg_co2->SetTextSize(0.03);
    leg_co2->AddEntry(standard_v2_la, "Standard reconstruction", "P");
    leg_co2->AddEntry(tight_v2_la, "Tight reconstruction", "P");
    leg_co2->AddEntry(loose_v2_la, "Loose reconstruction", "P");
    leg_co2->Draw();

    tex->SetTextFont(62);
    tex->SetTextSize(0.045);
    tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->SetTextSize(0.04);
    tex->SetTextFont(42);
    tex->DrawLatex(0.40,0.32,"185 #leq N_{trk}^{offline} < 250");
    line->Draw("same");

    c1_co->Print("v2CombinedRecoSystematics.pdf");

    //Calculate Ratios
    std::vector<double> v2RatioLooseKs;
    std::vector<double> v2RatioTightKs;

    std::vector<double> v2RatioLooseLa;
    std::vector<double> v2RatioTightLa;

    std::vector<double> v2ErrorRatioLooseKs;
    std::vector<double> v2ErrorRatioTightKs;

    std::vector<double> v2ErrorRatioLooseLa;
    std::vector<double> v2ErrorRatioTightLa;

    //Kshort
    for(int i=0; i<ks_npoints; i++)
    {
        v2RatioLooseKs.push_back(v2loose_ks[i]/v2standard_ks[i]);
        v2ErrorRatioLooseKs.push_back(v2loose_ksE[i]/v2standard_ks[i]);

        v2RatioTightKs.push_back(v2tight_ks[i]/v2standard_ks[i]);
        v2ErrorRatioTightKs.push_back(v2tight_ksE[i]/v2standard_ks[i]);
    }

    //Lambda
    for(int i=0; i<la_npoints; i++)
    {
        v2RatioLooseLa.push_back(v2loose_la[i]/v2standard_la[i]);
        v2ErrorRatioLooseLa.push_back(v2loose_laE[i]/v2standard_la[i]);

        v2RatioTightLa.push_back(v2tight_la[i]/v2standard_la[i]);
        v2ErrorRatioTightLa.push_back(v2tight_laE[i]/v2standard_la[i]);
    }

    double* av2RatioLooseKs = &v2RatioLooseKs[0];
    double* av2RatioLooseLa = &v2RatioLooseLa[0];

    double* av2RatioTightKs = &v2RatioTightKs[0];
    double* av2RatioTightLa = &v2RatioTightLa[0];

    double* av2ErrorRatioLooseKs = &v2ErrorRatioLooseKs[0];
    double* av2ErrorRatioLooseLa = &v2ErrorRatioLooseLa[0];

    double* av2ErrorRatioTightKs = &v2ErrorRatioTightKs[0];
    double* av2ErrorRatioTightLa = &v2ErrorRatioTightLa[0];

    TGraphErrors* Ratioloose_v2_la = new TGraphErrors(la_npoints,pTloose_la,av2RatioLooseLa,0,av2ErrorRatioLooseLa);
	TGraphErrors* Ratiotight_v2_la = new TGraphErrors(la_npoints,pTtight_la,av2RatioTightLa,0,av2ErrorRatioTightLa);

    TGraphErrors* Ratioloose_v2_ks = new TGraphErrors(ks_npoints,pTloose_ks,av2RatioLooseKs,0,av2ErrorRatioLooseKs);
	TGraphErrors* Ratiotight_v2_ks = new TGraphErrors(ks_npoints,pTtight_ks,av2RatioTightKs,0,av2ErrorRatioTightKs);

    Ratioloose_v2_ks->SetMarkerColor(kRed);
    Ratioloose_v2_ks->SetMarkerStyle(20);
    Ratioloose_v2_ks->SetMarkerSize(1.5);
    Ratioloose_v2_ks->SetLineColor(kRed);

    Ratiotight_v2_ks->SetMarkerColor(kBlue-4);
    Ratiotight_v2_ks->SetMarkerStyle(21);
    Ratiotight_v2_ks->SetMarkerSize(1.5);
    Ratiotight_v2_ks->SetLineColor(kBlue-4);

    Ratioloose_v2_la->SetMarkerColor(kRed);
    Ratioloose_v2_la->SetMarkerStyle(20);
    Ratioloose_v2_la->SetMarkerSize(1.5);
    Ratioloose_v2_la->SetLineColor(kRed);

    Ratiotight_v2_la->SetMarkerColor(kBlue-4);
    Ratiotight_v2_la->SetMarkerStyle(21);
    Ratiotight_v2_la->SetMarkerSize(1.5);
    Ratiotight_v2_la->SetLineColor(kBlue-4);

    TCanvas* c1_ratio_ks = MakeCanvas("c1_ratio_ks", "Plot_ratio_ks");
    TH1F* frame_Ratio_ks = c1_ratio_ks->DrawFrame(0,0.9,8,1.1);
    /*TH1F* frame_co = c1_co->DrawFrame(0,0.01,20,1);*/
    gPad->SetTickx();
    gPad->SetTicky();
    /*frame_Ratio_ks->SetTitle("K_{S}^{-1} Reconstruction Cut Ratio");*/
    frame_Ratio_ks->GetXaxis()->CenterTitle(1);
    frame_Ratio_ks->GetYaxis()->CenterTitle(1);
    frame_Ratio_ks->GetXaxis()->SetTitleSize(0.05);
    frame_Ratio_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Ratio_ks->GetYaxis()->SetTitle("v_{2}^{sig}");
    frame_Ratio_ks->GetYaxis()->SetTitleSize(0.05);
    frame_Ratio_ks->SetTitleOffset(1.2,"Y");
    frame_Ratio_ks->SetTitleOffset(1.2,"X");

    TLine* LineRatio_ks = new TLine(0,1,8,1);
    /*LineRatio_ks->SetLineStyle(2);*/
    LineRatio_ks->Draw();

    TLine* LineRatio_min_ks = new TLine(0,0.98,8,0.98);
    LineRatio_min_ks->SetLineStyle(2);
    LineRatio_min_ks->Draw();

    TLine* LineRatio_max_ks = new TLine(0,1.02,8,1.02);
    LineRatio_max_ks->SetLineStyle(2);
    LineRatio_max_ks->Draw();

    Ratioloose_v2_ks->Draw("P");
    Ratiotight_v2_ks->Draw("P");

    TLegend* leg_ratio_ks = new TLegend(0.45,0.75,0.57,0.85);
    leg_ratio_ks->SetFillColor(10);
    leg_ratio_ks->SetFillStyle(0);
    leg_ratio_ks->SetBorderSize(0);
    leg_ratio_ks->SetTextFont(42);
    leg_ratio_ks->SetTextSize(0.03);
    leg_ratio_ks->AddEntry(Ratioloose_v2_ks, "Loose K_{S}^{0} Reconstruction", "P");
    leg_ratio_ks->AddEntry(Ratiotight_v2_ks, "Tight K_{S}^{0} Reconstruction", "P");
    leg_ratio_ks->Draw();

    TCanvas* c1_ratio_la = MakeCanvas("c1_ratio_la", "Plot_ratio_la");
    TH1F* frame_Ratio_la = c1_ratio_la->DrawFrame(0,0.85,8,1.15);
    /*TH1F* frame_co = c1_co->DrawFrame(0,0.01,20,1);*/
    gPad->SetTickx();
    gPad->SetTicky();
    /*frame_Ratio_la->SetTitle("#Lambda Reconstruction Cut Ratio");*/
    frame_Ratio_la->GetXaxis()->CenterTitle(1);
    frame_Ratio_la->GetYaxis()->CenterTitle(1);
    frame_Ratio_la->GetXaxis()->SetTitleSize(0.05);
    frame_Ratio_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Ratio_la->GetYaxis()->SetTitle("v_{2}^{sig}");
    frame_Ratio_la->GetYaxis()->SetTitleSize(0.05);
    frame_Ratio_la->SetTitleOffset(1.1,"Y");
    frame_Ratio_la->SetTitleOffset(1.2,"X");

    TLine* LineRatio_la = new TLine(0,1,8,1);
    //LineRatio_la->SetLineStyle(2);
    LineRatio_la->Draw();

    TLine* LineRatio_min_la = new TLine(0,0.99,8,0.99);
    LineRatio_min_la->SetLineStyle(2);
    LineRatio_min_la->Draw();

    TLine* LineRatio_max_la = new TLine(0,1.01,8,1.01);
    LineRatio_max_la->SetLineStyle(2);
    LineRatio_max_la->Draw();

    Ratioloose_v2_la->Draw("P");
    Ratiotight_v2_la->Draw("P");

    TLegend* leg_ratio_la = new TLegend(0.45,0.75,0.57,0.85);
    leg_ratio_la->SetFillColor(10);
    leg_ratio_la->SetFillStyle(0);
    leg_ratio_la->SetBorderSize(0);
    leg_ratio_la->SetTextFont(42);
    leg_ratio_la->SetTextSize(0.03);
    leg_ratio_la->AddEntry(Ratioloose_v2_la, "Loose #Lambda / #bar{#Lambda} Reconstruction", "P");
    leg_ratio_la->AddEntry(Ratiotight_v2_la, "Tight #Lambda / #bar{#Lambda} Reconstruction", "P");
    leg_ratio_la->Draw();

    c1_ratio_ks->Print("v2RecoRatioSystematicsKs.pdf");
    c1_ratio_la->Print("v2RecoRatioSystematicsLa.pdf");

    // Write points into rootfile
    TFile out("V0ReconstructionSys.root","RECREATE");
    standard_v2_la->Write("StandardLa");
    tight_v2_la->Write("TightLa");
    loose_v2_la->Write("LooseLa");
    standard_v2_ks->Write("StandardKs");
    tight_v2_ks->Write("TightKs");
    loose_v2_ks->Write("LooseKs");
    Ratioloose_v2_ks->Write("RatioLooseKs");
    Ratiotight_v2_ks->Write("RatioTightKs");
    Ratioloose_v2_la->Write("RatioLooseLa");
    Ratiotight_v2_la->Write("RatioTightLa");
}

void V0_v4()
{
    MITStyle();
	const int ks_npoints = 11;
    double v2Ks8[ks_npoints]  = {0.00699898 ,0.00109379 ,0.00927869 ,0.00427604 ,0.0139798 ,0.0270377 ,0.0247767 ,0.033899 ,0.0343803 ,0.0364174 ,0.0397623};
    double pTKs8[ks_npoints]  = {0.3401, 0.5213, 0.7084, 0.9036, 1.201, 1.591, 1.986, 2.465, 3.136, 4.008, 5.071};
    double v2Ks8E[ks_npoints] = {0.0240717 ,0.00749964 ,0.00449023 ,0.00362903 ,0.00238665 ,0.00263102 ,0.00321623 ,0.00353347 ,0.00471759 ,0.00720811 ,0.0113695};

	const int la_npoints = 8;
    double v2La8[la_npoints]  = {-0.0199048 ,0.00945889 ,0.0146571 ,0.00929574 ,0.0102335 ,0.0408502 ,0.0596036 ,0.0702546};
    double pTLa8[la_npoints]  = {0.9145, 1.221, 1.605, 1.997, 2.481, 3.149, 4.004, 5.041};
    double v2La8E[la_npoints] = {0.016915 ,0.00708958 ,0.00586441 ,0.005723 ,0.00527798 ,0.00619893 ,0.00917289 ,0.015673};

    // Pull TGraph for Kshort and lambda


    TGraphErrors* ks8_v2 = new TGraphErrors(ks_npoints,pTKs8,v2Ks8,0,v2Ks8E);
	TGraphErrors* la8_v2 = new TGraphErrors(la_npoints,pTLa8,v2La8,0,v2La8E);

    ks8_v2->SetMarkerColor(kRed);
    ks8_v2->SetMarkerStyle(20);
    ks8_v2->SetMarkerSize(1.5);
    ks8_v2->SetLineColor(kRed);

    /*la8_v2->SetMarkerColor(kGreen+2);*/
    la8_v2->SetMarkerColor(kBlue-4);
    la8_v2->SetMarkerStyle(22);
    la8_v2->SetMarkerSize(1.5);
    la8_v2->SetLineColor(kBlue-4);

    TCanvas* c1 = MakeCanvas("c1", "Plot");
    c1->cd();
    //c1->SetLogy();
    c1->SetLeftMargin(0.12);

    // draw the frame using a histogram frame

    TH1F* frame = c1->DrawFrame(0,-0.01,6,0.1);
    gPad->SetTickx();
    gPad->SetTicky();
    frame->GetXaxis()->CenterTitle(1);
    frame->GetYaxis()->CenterTitle(1);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame->GetYaxis()->SetTitle("v#kern[-0.3]{_{4}}");
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->SetTitleOffset(1.2,"Y");
    frame->SetTitleOffset(1.2,"X");

    //
    // draw the graph with axis, continuous line, and put
    ks8_v2->Draw("P");
    la8_v2->Draw("P");
    Int_t oldColor = 38;
    Int_t newColor = 46;

    TLegend* leg = new TLegend(0.15,0.65,0.61,0.78);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->AddEntry(ks8_v2, "#color[46]{K_{#lower[-0.3]{S}}#kern[-1.05]{#lower[0.1]{{}^{0}}}}", "P");
    leg->AddEntry(la8_v2, "#color[46]{#Lambda / #bar{#Lambda}}", "P");
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.05);
    tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = #color[46]{8.16} TeV");
    tex->SetTextSize(0.045);
    //tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
    // tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
    tex->DrawLatex(0.45,0.73,"185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");
    TLine* line = new TLine(0,0,6,0);
    line->SetLineStyle(8);
    line->Draw();

    c1->Print("V0_v4.pdf");
}

void V0_v5()
{
    MITStyle();
	const int ks_npoints = 11;
    double v2Ks8[ks_npoints]  = {0.0158506 ,-0.00552973 ,0.00581856 ,0.0114181 ,0.0174554 ,-0.00281132 ,0.0152802 ,0.00859869 ,0.00261752 ,0.0231879 ,0.00377595};
    double pTKs8[ks_npoints]  = {0.3401, 0.5213, 0.7084, 0.9036, 1.201, 1.591, 1.986, 2.465, 3.136, 4.008, 5.071};
    double v2Ks8E[ks_npoints] = {0.0653367 ,0.0203568 ,0.0121898 ,0.00985446 ,0.00648768 ,0.00714284 ,0.00873752 ,0.00959518 ,0.012811 ,0.0195802 ,0.0308744};


	const int la_npoints = 8;
    double v2La8[la_npoints]  = {-0.0139506 ,-0.0214405 ,0.0102338 ,-0.012811 ,0.0264369 ,0.0301995 ,0.0309084 ,0.0818709};
    double pTLa8[la_npoints]  = {0.9145, 1.221, 1.605, 1.997, 2.481, 3.149, 4.004, 5.041};
    double v2La8E[la_npoints] = {0.0459104 ,0.0192491 ,0.0159232 ,0.0155425 ,0.0143434 ,0.0168456 ,0.0249281 ,0.0425993};

    // Pull TGraph for Kshort and lambda


    TGraphErrors* ks8_v2 = new TGraphErrors(ks_npoints,pTKs8,v2Ks8,0,v2Ks8E);
	TGraphErrors* la8_v2 = new TGraphErrors(la_npoints,pTLa8,v2La8,0,v2La8E);

    ks8_v2->SetMarkerColor(kRed);
    ks8_v2->SetMarkerStyle(20);
    ks8_v2->SetMarkerSize(1.5);
    ks8_v2->SetLineColor(kRed);

    /*la8_v2->SetMarkerColor(kGreen+2);*/
    la8_v2->SetMarkerColor(kBlue-4);
    la8_v2->SetMarkerStyle(22);
    la8_v2->SetMarkerSize(1.5);
    la8_v2->SetLineColor(kBlue-4);

    TCanvas* c1 = MakeCanvas("c1", "Plot");
    c1->cd();
    //c1->SetLogy();
    c1->SetLeftMargin(0.12);

    // draw the frame using a histogram frame

    TH1F* frame = c1->DrawFrame(0,-0.05,6,0.15);
    gPad->SetTickx();
    gPad->SetTicky();
    frame->GetXaxis()->CenterTitle(1);
    frame->GetYaxis()->CenterTitle(1);
    frame->GetXaxis()->SetTitleSize(0.05);
    frame->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame->GetYaxis()->SetTitle("v#kern[-0.3]{_{5}}");
    frame->GetYaxis()->SetTitleSize(0.05);
    frame->SetTitleOffset(1.2,"Y");
    frame->SetTitleOffset(1.2,"X");

    //
    // draw the graph with axis, continuous line, and put
    ks8_v2->Draw("P");
    la8_v2->Draw("P");
    Int_t oldColor = 38;
    Int_t newColor = 46;

    TLegend* leg = new TLegend(0.15,0.65,0.61,0.78);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->AddEntry(ks8_v2, "#color[46]{K_{#lower[-0.3]{S}}#kern[-1.05]{#lower[0.1]{{}^{0}}}}", "P");
    leg->AddEntry(la8_v2, "#color[46]{#Lambda / #bar{#Lambda}}", "P");
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.05);
    tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = #color[46]{8.16} TeV");
    tex->SetTextSize(0.045);
    //tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
    // tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
    tex->DrawLatex(0.45,0.73,"185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");
    TLine* line = new TLine(0,0,6,0);
    line->SetLineStyle(8);
    line->Draw();

    c1->Print("V0_v5.pdf");
}
