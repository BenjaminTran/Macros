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
    leg->SetBorderSize(0.035);
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
		leg->SetBorderSize(0.035);
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
    leg->SetBorderSize(0.035);
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

void Rap_Cascade_v2()
{
    MITStyle();

	const int xi_npoints = 9;
    double v2Xi8[xi_npoints]  = {0.0485332 ,0.0602203 ,0.0815388 ,0.126116 ,0.161428 ,0.184352 ,0.208616 ,0.168331 ,0.281171};
    double pTXi8[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.938, 11.31};
    double v2Xi8E[xi_npoints] = {0.0117314 ,0.00614764 ,0.00530833 ,0.00418641 ,0.0039623 ,0.00462749 ,0.00664224 ,0.0123018 ,0.0657875};

	const int ks_npoints = 16;
    double v2Ks8[ks_npoints]  = {0.0165183 ,0.0256269 ,0.0427305 ,0.0592655 ,0.0817776 ,0.105675 ,0.124798 ,0.138062 ,0.146527 ,0.14193 ,0.138978 ,0.124995 ,0.131851 ,0.12689 ,0.149737 ,0.1126};
    double pTKs8[ks_npoints]  = {0.3666, 0.5309, 0.711, 0.9046, 1.202, 1.591, 1.986, 2.465, 3.136, 4.008, 5.142, 6.431, 7.619, 9.142, 11.64, 16.86};
    double v2Ks8E[ks_npoints] = {0.00664835 ,0.00160315 ,0.000940095 ,0.000767227 ,0.000500488 ,0.000537786 ,0.000645228 ,0.000699308 ,0.000924693 ,0.00138942 ,0.00218948 ,0.00454549 ,0.00597703 ,0.00975685 ,0.0131384 ,0.0957745};

	const int la_npoints = 13;
    double v2La8[la_npoints]  = {0.031126 ,0.0516416 ,0.0770812 ,0.107901 ,0.136777 ,0.170426 ,0.194357 ,0.197845 ,0.205666 ,0.184392 ,0.150626 ,0.176849 ,0.0713852};
    double pTLa8[la_npoints]  = {0.9252, 1.224, 1.603, 1.995, 2.485, 3.156, 4.007, 5.116, 6.414, 7.588, 9.117, 11.56, 16.89};
    double v2La8E[la_npoints] = {0.00329462 ,0.00134879 ,0.00119965 ,0.00122695 ,0.0011159 ,0.00123679 ,0.00174121 ,0.0029136 ,0.00695908 ,0.00996631 ,0.0198574 ,0.0229987 ,0.19794};


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
    /*c1->SetLogy();*/
    c1->SetLeftMargin(0.12);

    // draw the frame using a histogram frame

    TH1F* frame = c1->DrawFrame(0,-0.01,20,0.41);
    /*TH1F* frame = c1->DrawFrame(0,0.01,20,1);*/
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

    TLegend* leg = new TLegend(0.6532847,0.6084211,0.8540146,0.8589474);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0.035);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    //leg->AddEntry(ha_v2, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}", "P");
    leg->AddEntry(ks8_v2, "K_{S}^{0}", "P");
    leg->AddEntry(la8_v2, "#Lambda / #bar{#Lambda}", "P");
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
    tex->DrawLatex(0.1441606,0.6673684,"185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 250");
    tex->DrawLatex(0.15,0.74,"|y| < 1");
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");
    TLine* line = new TLine(0,0,6,0);
    line->SetLineStyle(8);
    line->Draw();

    c1->Print("v2Rapidity.pdf");

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
    leg->SetBorderSize(0.035);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->AddEntry(ks8_v2, "#color[46]{K_{#lower[-0.3]{S}}#kern[-1.05]{#lower[0.1]{{}^{0}}}}", "P");
    leg->AddEntry(la8_v2, "#color[46]{#Lambda / #lower[0.1]{#bar{  }}#kern[-1.2]{#Lambda}}", "P");
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
    leg->SetBorderSize(0.035);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->AddEntry(ks8_v2, "#color[46]{K_{#lower[-0.3]{S}}#kern[-1.05]{#lower[0.1]{{}^{0}}}}", "P");
    leg->AddEntry(la8_v2, "#color[46]{#Lambda / #lower[0.1]{#bar{  }}#kern[-1.2]{#Lambda}}", "P");
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
}
