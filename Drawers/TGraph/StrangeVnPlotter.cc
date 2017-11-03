#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TH2.h"
#include "TVirtualFitter.h"
#include "interface/GetGraphFromFile.C"
#include "interface/MITStyle.C"
#include "RooChebychev.h"

std::string v2RootFileName = "v2V0sGraphs185_250D0Ana.root";
using namespace RooFit;

TGraphErrors* TGDivideSameX(TGraphErrors* TG1, TGraphErrors* TG2) // For dividing tgraph 1 by tgraph2. Need to have the same X axis.
{
    TGraphErrors* TGDiv;
    double* X = TG1->GetX();
    double* TG1_Y = TG1->GetY();
    double* TG1_EY = TG1->GetEY();
    double* TG2_Y = TG2->GetY();

    std::vector<double> TGDiv_Y;
    std::vector<double> TGDiv_EY;

    for(int i=0; i<TG1->GetN(); i++)
    {
        TGDiv_Y.push_back(TG1_Y[i]/TG2_Y[i]);
        TGDiv_EY.push_back(TG1_EY[i]/TG2_Y[i]);
    }
    TGDiv = new TGraphErrors(TG1->GetN(),X,&TGDiv_Y[0],0,&TGDiv_EY[0]);
    return TGDiv;
}

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

	const int om_npoints = 9;
    double v2Om8[om_npoints]  = {0.0191671 ,0.0770592 ,0.0777969 ,0.101627 ,0.146736 ,0.201549 ,0.251998 ,0.261003 ,0.189401};
    double pTOm8[om_npoints]  = {1.319, 1.658, 2.002, 2.492, 3.164, 4.026, 5.138, 6.476, 8.05};
    double v2Om8E[om_npoints] = {0.0365897 ,0.0213553 ,0.0133332 ,0.00922245 ,0.00840292 ,0.00953781 ,0.0132644 ,0.0286688 ,0.043903};

	const int ks_npoints = 11;
    double v2Ks8[ks_npoints]  = {0.012813549, 0.028133917,0.043212, 0.05876162, 0.081288146 ,0.105507591 ,0.123968514 ,0.135179572 ,0.144205576 ,0.142730876 ,0.129697566};
    double pTKs8[ks_npoints]  = {0.3401, 0.5213, 0.7084, 0.9036, 1.201, 1.591, 1.986, 2.465, 3.136, 4.008, 5.071};
    double v2Ks8E[ks_npoints] = {0.00355244, 0.00110708,0.000619229, 0.000500387 ,0.000329161 ,0.000362995 ,0.000443819 ,0.000487576 ,0.000650721 ,0.00099496 ,0.00156744};

	const int la_npoints = 8;
    double v2La8[la_npoints]  = {0.031657 ,0.055453976 ,0.080209188 ,0.105853753 ,0.137719142 ,0.170879559 ,0.190002697 ,0.190002697};
    double pTLa8[la_npoints]  = {0.9145, 1.221, 1.605, 1.997, 2.481, 3.149, 4.004, 5.041};
    double v2La8E[la_npoints] = {0.002330326 ,0.000977656 ,0.000808941 ,0.000789458 ,0.000724636 ,0.000850324 ,0.001268046 ,0.002144964};

    // Pull TGraph for Kshort and lambda

    //TGraphErrors* ha_v2 = (TGraphErrors*)file_hadv2->Get("hadronv2");
    TGraphErrors* ha_v2 = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033pPb_v2.txt"),1,28,1.2);

    TGraphErrors* xi8_v2 = new TGraphErrors(xi_npoints,pTXi8,v2Xi8,0,v2Xi8E);
    TGraphErrors* om8_v2 = new TGraphErrors(om_npoints,pTOm8,v2Om8,0,v2Om8E);
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

    om8_v2->SetMarkerColor(kMagenta-4);
    om8_v2->SetMarkerStyle(29);
    om8_v2->SetMarkerSize(1.5);
    om8_v2->SetLineColor(kMagenta-4);

    la8_v2->SetMarkerColor(kBlue-4);
    la8_v2->SetMarkerStyle(22);
    la8_v2->SetMarkerSize(1.5);
    la8_v2->SetLineColor(kBlue-4);

    TCanvas* c1 = MakeCanvas("c1", "Plot");
    c1->cd();
    c1->SetLeftMargin(0.12);

    // draw the frame using a histogram frame

    TH1F* frame = c1->DrawFrame(0,-0.01,10,0.61);
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
    om8_v2->Draw("P");

    // Write points into rootfile
    TFile out("8TeVgraphv2_pPb185-220.root","RECREATE");
    ks8_v2->Write("kshortv2");
    la8_v2->Write("lambdav2");
    xi8_v2->Write("xiv2");
    om8_v2->Write("omegav2");

    TLegend* leg = new TLegend(0.15,0.45,0.51,0.68);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    //leg->AddEntry(ha_v2, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}", "P");
    leg->AddEntry(ks8_v2, "K_{S}^{0}", "P");
    leg->AddEntry(la8_v2, "#Lambda/#bar{#Lambda}", "P");
    leg->AddEntry(xi8_v2, "#Xi^{#pm}", "P");
    leg->AddEntry(om8_v2, "#Omega^{#pm}", "P");
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
    TLine* line = new TLine(0,0,10,0);
    line->SetLineStyle(8);
    line->Draw();

    c1->Print("V2all.pdf");
    c1->Print("V2all.png");
}

void Rap_v2sig_pPb()
{
    MITStyle();
    TCanvas* c1 = MakeCanvas("c1", "Plot");
    /*c1->SetLogy();*/
    c1->SetLeftMargin(0.12);

    TCanvas* c2 = MakeCanvas("c2", "Plot");
    c2->SetLeftMargin(0.12);

    c1->cd();

    // draw the frame using a histogram frame
    TH1F* frame1;
    TH1F* frame2;

    frame1 = c1->DrawFrame(0,-0.01,9,0.45);
    gPad->SetTickx();
    gPad->SetTicky();
    frame1->GetXaxis()->CenterTitle(1);
    frame1->GetYaxis()->CenterTitle(1);
    frame1->GetXaxis()->SetTitleSize(0.05);
    frame1->GetYaxis()->SetTitleSize(0.05);
    frame1->SetTitleOffset(1.1,"Y");
    frame1->SetTitleOffset(1.2,"X");
    frame1->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame1->GetYaxis()->SetTitle("v_{2}^{sig}");
    const int om_npoints = 8;
    const int xi_npoints = 9;
    const int ks_npoints = 13;
    const int la_npoints = 10;

    //frame = c1->DrawFrame(0,-0.01,8,0.5);
    //gPad->SetTickx();
    //gPad->SetTicky();
    //frame->GetXaxis()->CenterTitle(1);
    //frame->GetYaxis()->CenterTitle(1);
    //frame->GetXaxis()->SetTitleSize(0.05);
    //frame->GetYaxis()->SetTitleSize(0.05);
    //frame->SetTitleOffset(1.1,"Y");
    //frame->SetTitleOffset(1.2,"X");
    //frame->GetXaxis()->SetTitle("KE_{T} (GeV)");
    //frame->GetYaxis()->SetTitle("v_{2}^{sig}");

    // Pull TGraph for Kshort and lambda
    TFile* file_pidv2 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/v2valuesRapidity_pPb_185_250_10_31_17.root");
    //TFile* file_pidv2_Omega = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/v2valuesRapidityOmegaCheck.root");


    //TGraphErrors* om8_v2_Fixed = (TGraphErrors*)file_pidv2_Omega->Get("v2omega");
    TGraphErrors* om8_v2 = (TGraphErrors*)file_pidv2->Get("v2omega");
    TGraphErrors* xi8_v2 = (TGraphErrors*)file_pidv2->Get("v2xi");
    TGraphErrors* ks8_v2 = (TGraphErrors*)file_pidv2->Get("v2kshort");
    TGraphErrors* la8_v2 = (TGraphErrors*)file_pidv2->Get("v2lambda");

    TGraphErrors* om8_v2kn = (TGraphErrors*)file_pidv2->Get("v2omega_ket_nq");
    TGraphErrors* xi8_v2kn = (TGraphErrors*)file_pidv2->Get("v2xi_ket_nq");
    TGraphErrors* ks8_v2kn = (TGraphErrors*)file_pidv2->Get("v2kshort_ket_nq");
    TGraphErrors* la8_v2kn = (TGraphErrors*)file_pidv2->Get("v2lambda_ket_nq");

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

    om8_v2->SetMarkerColor(kMagenta);
    om8_v2->SetMarkerStyle(29);
    om8_v2->SetMarkerSize(1.5);
    om8_v2->SetLineColor(kMagenta);

    //om8_v2_Fixed->SetMarkerColor(kBlack);
    //om8_v2_Fixed->SetMarkerStyle(29);
    //om8_v2_Fixed->SetMarkerSize(1.5);
    //om8_v2_Fixed->SetLineColor(kBlack);

    ks8_v2kn->SetMarkerColor(kRed);
    ks8_v2kn->SetMarkerStyle(20);
    ks8_v2kn->SetMarkerSize(1.5);
    ks8_v2kn->SetLineColor(kRed);

    xi8_v2kn->SetMarkerColor(kGreen+2);
    xi8_v2kn->SetMarkerStyle(21);
    xi8_v2kn->SetMarkerSize(1.5);
    xi8_v2kn->SetLineColor(kGreen+2);

    la8_v2kn->SetMarkerColor(kBlue-4);
    la8_v2kn->SetMarkerStyle(22);
    la8_v2kn->SetMarkerSize(1.5);
    la8_v2kn->SetLineColor(kBlue-4);

    om8_v2kn->SetMarkerColor(kMagenta);
    om8_v2kn->SetMarkerStyle(29);
    om8_v2kn->SetMarkerSize(1.5);
    om8_v2kn->SetLineColor(kMagenta);

    c1->cd();

    TLegend* leg = new TLegend(0.15,0.55,0.27,0.75);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    //leg->AddEntry(ha_v2, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}", "P");
    leg->AddEntry(ks8_v2, "K_{S}^{0}", "P");
    leg->AddEntry(la8_v2, "#Lambda / #bar{#Lambda}", "P");
    leg->AddEntry(xi8_v2, "#Xi^{#pm}", "P");
    leg->AddEntry(om8_v2, "#Omega^{#pm}", "P");
    leg->Draw();


    //ha_v2->Draw("PESAME");
    ks8_v2->Draw("P");
    la8_v2->Draw("P");
    xi8_v2->Draw("P");
    om8_v2->Draw("P");
    //om8_v2_Fixed->Draw("P");

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(62);
    tex->SetTextSize(0.05);
    tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->SetTextSize(0.045);
    tex->SetTextFont(42);
    tex->DrawLatex(0.34,0.72,"185 #leq N_{trk}^{offline} < 250");
    /*tex->DrawLatex(0.15,0.74,"|y| < 1");*/
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");

    c2->cd();

    frame2 = c2->DrawFrame(0,-0.01,4.5,0.3);
    gPad->SetTickx();
    gPad->SetTicky();
    frame2->GetXaxis()->CenterTitle(1);
    frame2->GetYaxis()->CenterTitle(1);
    frame2->GetXaxis()->SetTitleSize(0.05);
    frame2->GetYaxis()->SetTitleSize(0.05);
    frame2->SetTitleOffset(1.1,"Y");
    frame2->SetTitleOffset(1.2,"X");
    frame2->GetXaxis()->SetTitle("KE_{T}/n_{q} (GeV)");
    frame2->GetYaxis()->SetTitle("v_{2}^{sig}/n_{q}");

    TLegend* legkn = new TLegend(0.15,0.55,0.27,0.75);
    legkn->SetFillColor(10);
    legkn->SetFillStyle(0);
    legkn->SetBorderSize(0);
    legkn->SetTextFont(42);
    legkn->SetTextSize(0.05);
    legkn->AddEntry(ks8_v2kn, "K_{S}^{0}", "P");
    legkn->AddEntry(la8_v2kn, "#Lambda / #bar{#Lambda}", "P");
    legkn->AddEntry(xi8_v2kn, "#Xi^{#pm}", "P");
    legkn->AddEntry(om8_v2kn, "#Omega^{#pm}", "P");
    legkn->Draw();

    ks8_v2kn->Draw("P");
    la8_v2kn->Draw("P");
    xi8_v2kn->Draw("P");
    om8_v2kn->Draw("P");

    // Draw Legend and write points into rootfile

    tex->SetTextFont(62);
    tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->SetTextSize(0.045);
    tex->SetTextFont(42);
    tex->DrawLatex(0.40,0.24,"185 #leq N_{trk}^{offline} < 250");
    /*tex->DrawLatex(0.15,0.74,"|y| < 1");*/
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");

    c1->Print("v2SigRapiditypPb.pdf");
    c1->Print("v2SigRapiditypPb.png");
    //c1->Print("v2SigRapidityKET.pdf");
    c2->Print("v2SigRapidityDividednqpPb.pdf");
    c2->Print("v2SigRapidityDividednqpPb.png");
}
void Rap_v2sig_PbPb()
{
    MITStyle();
    TCanvas* c1 = MakeCanvas("c1", "Plot");
    /*c1->SetLogy();*/
    c1->SetLeftMargin(0.12);

    TCanvas* c2 = MakeCanvas("c2", "Plot");
    c2->SetLeftMargin(0.12);

    c1->cd();

    // draw the frame using a histogram frame
    TH1F* frame1;
    TH1F* frame2;

    frame1 = c1->DrawFrame(0,-0.01,9,0.6);
    gPad->SetTickx();
    gPad->SetTicky();
    frame1->GetXaxis()->CenterTitle(1);
    frame1->GetYaxis()->CenterTitle(1);
    frame1->GetXaxis()->SetTitleSize(0.05);
    frame1->GetYaxis()->SetTitleSize(0.05);
    frame1->SetTitleOffset(1.1,"Y");
    frame1->SetTitleOffset(1.2,"X");
    frame1->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame1->GetYaxis()->SetTitle("v_{2}^{sig}");
    const int om_npoints = 8;
    const int xi_npoints = 9;
    const int ks_npoints = 13;
    const int la_npoints = 10;

    //frame = c1->DrawFrame(0,-0.01,8,0.5);
    //gPad->SetTickx();
    //gPad->SetTicky();
    //frame->GetXaxis()->CenterTitle(1);
    //frame->GetYaxis()->CenterTitle(1);
    //frame->GetXaxis()->SetTitleSize(0.05);
    //frame->GetYaxis()->SetTitleSize(0.05);
    //frame->SetTitleOffset(1.1,"Y");
    //frame->SetTitleOffset(1.2,"X");
    //frame->GetXaxis()->SetTitle("KE_{T} (GeV)");
    //frame->GetYaxis()->SetTitle("v_{2}^{sig}");

    // Pull TGraph for Kshort and lambda
    TFile* file_pidv2 = TFile::Open("rootFiles/v2valuesRapidityPbPb.root");
    //TFile* file_pidv2 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/v2valuesRapidity_pPb_185_250_10_31_17.root");
    //TFile* file_pidv2_Omega = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/v2valuesRapidityOmegaCheck.root");


    //TGraphErrors* om8_v2_Fixed = (TGraphErrors*)file_pidv2_Omega->Get("v2omega");
    TGraphErrors* om8_v2 = (TGraphErrors*)file_pidv2->Get("v2omega");
    TGraphErrors* xi8_v2 = (TGraphErrors*)file_pidv2->Get("v2xi");
    TGraphErrors* ks8_v2 = (TGraphErrors*)file_pidv2->Get("v2kshort");
    TGraphErrors* la8_v2 = (TGraphErrors*)file_pidv2->Get("v2lambda");

    TGraphErrors* om8_v2kn = (TGraphErrors*)file_pidv2->Get("v2omega_ket_nq");
    TGraphErrors* xi8_v2kn = (TGraphErrors*)file_pidv2->Get("v2xi_ket_nq");
    TGraphErrors* ks8_v2kn = (TGraphErrors*)file_pidv2->Get("v2kshort_ket_nq");
    TGraphErrors* la8_v2kn = (TGraphErrors*)file_pidv2->Get("v2lambda_ket_nq");

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

    om8_v2->SetMarkerColor(kMagenta);
    om8_v2->SetMarkerStyle(29);
    om8_v2->SetMarkerSize(1.5);
    om8_v2->SetLineColor(kMagenta);

    //om8_v2_Fixed->SetMarkerColor(kBlack);
    //om8_v2_Fixed->SetMarkerStyle(29);
    //om8_v2_Fixed->SetMarkerSize(1.5);
    //om8_v2_Fixed->SetLineColor(kBlack);

    ks8_v2kn->SetMarkerColor(kRed);
    ks8_v2kn->SetMarkerStyle(20);
    ks8_v2kn->SetMarkerSize(1.5);
    ks8_v2kn->SetLineColor(kRed);

    xi8_v2kn->SetMarkerColor(kGreen+2);
    xi8_v2kn->SetMarkerStyle(21);
    xi8_v2kn->SetMarkerSize(1.5);
    xi8_v2kn->SetLineColor(kGreen+2);

    la8_v2kn->SetMarkerColor(kBlue-4);
    la8_v2kn->SetMarkerStyle(22);
    la8_v2kn->SetMarkerSize(1.5);
    la8_v2kn->SetLineColor(kBlue-4);

    om8_v2kn->SetMarkerColor(kMagenta);
    om8_v2kn->SetMarkerStyle(29);
    om8_v2kn->SetMarkerSize(1.5);
    om8_v2kn->SetLineColor(kMagenta);

    c1->cd();

    TLegend* leg = new TLegend(0.15,0.55,0.27,0.75);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    //leg->AddEntry(ha_v2, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}", "P");
    leg->AddEntry(ks8_v2, "K_{S}^{0}", "P");
    leg->AddEntry(la8_v2, "#Lambda / #bar{#Lambda}", "P");
    leg->AddEntry(xi8_v2, "#Xi^{#pm}", "P");
    leg->AddEntry(om8_v2, "#Omega^{#pm}", "P");
    leg->Draw();


    //ha_v2->Draw("PESAME");
    ks8_v2->Draw("P");
    la8_v2->Draw("P");
    xi8_v2->Draw("P");
    om8_v2->Draw("P");
    //om8_v2_Fixed->Draw("P");

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(62);
    tex->SetTextSize(0.05);
    //tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->DrawLatex(0.15,0.8,"CMS PbPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
    tex->SetTextSize(0.045);
    //tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
    // tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
    tex->SetTextFont(42);
    //tex->DrawLatex(0.40,0.24,"185 #leq N_{trk}^{offline} < 250");
    tex->DrawLatex(0.36,0.72,"Cent: 30% - 50%");
    /*tex->DrawLatex(0.15,0.74,"|y| < 1");*/
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");

    c2->cd();

    frame2 = c2->DrawFrame(0,-0.01,4.5,0.3);
    gPad->SetTickx();
    gPad->SetTicky();
    frame2->GetXaxis()->CenterTitle(1);
    frame2->GetYaxis()->CenterTitle(1);
    frame2->GetXaxis()->SetTitleSize(0.05);
    frame2->GetYaxis()->SetTitleSize(0.05);
    frame2->SetTitleOffset(1.1,"Y");
    frame2->SetTitleOffset(1.2,"X");
    frame2->GetXaxis()->SetTitle("KE_{T}/n_{q} (GeV)");
    frame2->GetYaxis()->SetTitle("v_{2}^{sig}/n_{q}");

    TLegend* legkn = new TLegend(0.15,0.55,0.27,0.75);
    legkn->SetFillColor(10);
    legkn->SetFillStyle(0);
    legkn->SetBorderSize(0);
    legkn->SetTextFont(42);
    legkn->SetTextSize(0.05);
    legkn->AddEntry(ks8_v2kn, "K_{S}^{0}", "P");
    legkn->AddEntry(la8_v2kn, "#Lambda / #bar{#Lambda}", "P");
    legkn->AddEntry(xi8_v2kn, "#Xi^{#pm}", "P");
    legkn->AddEntry(om8_v2kn, "#Omega^{#pm}", "P");
    legkn->Draw();

    ks8_v2kn->Draw("P");
    la8_v2kn->Draw("P");
    xi8_v2kn->Draw("P");
    om8_v2kn->Draw("P");

    // Draw Legend and write points into rootfile

    tex->SetTextFont(62);
    //tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->DrawLatex(0.15,0.8,"CMS PbPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
    tex->SetTextSize(0.045);
    //tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
    // tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
    tex->SetTextFont(42);
    //tex->DrawLatex(0.40,0.24,"185 #leq N_{trk}^{offline} < 250");
    tex->DrawLatex(0.33,0.72,"Cent: 30% - 50%");
    /*tex->DrawLatex(0.15,0.74,"|y| < 1");*/
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");

    c1->Print("v2SigRapidityPbPb.pdf");
    c1->Print("v2SigRapidityPbPb.png");
    //c1->Print("v2SigRapidityKET.pdf");
    c2->Print("v2SigRapidityDividednqPbPb.pdf");
    c2->Print("v2SigRapidityDividednqPbPb.png");
}

void Rap_v2obsbkg_pPb()
{
    TFile* f = new TFile("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/v2valuesRapidity_pPb_185_250_10_31_17.root");
    MITStyle();
    TCanvas* c1 = MakeCanvas("c1", "Individual");
    c1->SetLeftMargin(0.12);

    TCanvas* c2 = MakeCanvas("c2", "Combined");
    c2->cd();
    c2->SetLeftMargin(0.12);

    // draw the frame using a histogram frame
    TH1F* frame_co = c2->DrawFrame(0,-0.05,9,0.45);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_co->GetXaxis()->CenterTitle(1);
    frame_co->GetYaxis()->CenterTitle(1);
    frame_co->GetXaxis()->SetTitleSize(0.05);
    frame_co->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_co->GetYaxis()->SetTitle("v_{2}^{Obs,Bkg}");
    frame_co->GetYaxis()->SetTitleSize(0.05);
    frame_co->SetTitleOffset(1.1,"Y");
    frame_co->SetTitleOffset(1.2,"X");

    TLegend* leg_co_spec = new TLegend(0.15,0.55,0.27,0.75);
    leg_co_spec->SetFillColor(10);
    leg_co_spec->SetFillStyle(0);
    leg_co_spec->SetBorderSize(0);
    leg_co_spec->SetTextFont(42);
    leg_co_spec->SetTextSize(0.05);

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(62);
    tex->SetTextSize(0.05);

    std::string legendLabel = "";
    TLegend* leg_co_Label = new TLegend(0.75,0.75,0.85,0.85);
    leg_co_Label->SetFillColor(10);
    leg_co_Label->SetFillStyle(0);
    leg_co_Label->SetBorderSize(0);
    leg_co_Label->SetTextFont(42);
    leg_co_Label->SetTextSize(0.05);

    c1->cd();

    std::string kshortv2 = "";
    std::string lambdav2 = "";
    std::string cascadev2 = "";
    std::string yaxis = "";

    for(int i=0; i<2; i++)
    {
        TGraphErrors* xi8_v2;
        TGraphErrors* om8_v2;
        TGraphErrors* ks8_v2;
        TGraphErrors* la8_v2;

        if(i==0)
        {
            ks8_v2 = (TGraphErrors*)f->Get("v2obskshort");
            la8_v2 = (TGraphErrors*)f->Get("v2obslambda");
            xi8_v2 = (TGraphErrors*)f->Get("v2obsxi");
            om8_v2 = (TGraphErrors*)f->Get("v2obsomega");

            yaxis = "v_{2}^{obs}";
        }
        else
        {
            ks8_v2 = (TGraphErrors*)f->Get("v2bkgkshort");
            la8_v2 = (TGraphErrors*)f->Get("v2bkglambda");
            xi8_v2 = (TGraphErrors*)f->Get("v2bkgxi");
            om8_v2 = (TGraphErrors*)f->Get("v2bkgomega");

            yaxis = "v_{2}^{bkg}";

            leg_co_spec->Draw();
            leg_co_Label->Draw();
            tex->SetTextFont(62);
            tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
            tex->SetTextSize(0.045);
            tex->SetTextFont(42);
            tex->DrawLatex(0.40,0.24,"185 #leq N_{trk}^{offline} < 250");
            c1->cd();
        }

        ks8_v2->SetMarkerColor(kRed);
        ks8_v2->SetMarkerStyle(20);
        ks8_v2->SetMarkerSize(1.5);
        ks8_v2->SetLineColor(kRed);

        xi8_v2->SetMarkerColor(kGreen+2);
        xi8_v2->SetMarkerStyle(21);
        xi8_v2->SetMarkerSize(1.5);
        xi8_v2->SetLineColor(kGreen+2);

        om8_v2->SetMarkerColor(kMagenta);
        om8_v2->SetMarkerStyle(29);
        om8_v2->SetMarkerSize(1.5);
        om8_v2->SetLineColor(kMagenta);

        la8_v2->SetMarkerColor(kBlue-4);
        la8_v2->SetMarkerStyle(22);
        la8_v2->SetMarkerSize(1.5);
        la8_v2->SetLineColor(kBlue-4);

        // draw the frame using a histogram frame
        TH1F* frame = c1->DrawFrame(0,-0.05,9,0.45);
        //TH1F* frame = c1->DrawFrame(0,0.01,20,1);
        gPad->SetTickx();
        gPad->SetTicky();
        frame->GetXaxis()->CenterTitle(1);
        frame->GetYaxis()->CenterTitle(1);
        frame->GetXaxis()->SetTitleSize(0.05);
        frame->GetXaxis()->SetTitle("p_{T} (GeV)");
        frame->GetYaxis()->SetTitle(yaxis.c_str());
        frame->GetYaxis()->SetTitleSize(0.05);
        frame->SetTitleOffset(1.1,"Y");
        frame->SetTitleOffset(1.2,"X");

        //ha_v2->Draw("PESAME");
        ks8_v2->Draw("P");
        la8_v2->Draw("P");
        xi8_v2->Draw("P");
        om8_v2->Draw("P");

        TLegend* leg = new TLegend(0.15,0.55,0.27,0.75);
        leg->SetFillColor(10);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.05);
        leg->AddEntry(ks8_v2, "K_{S}^{0}", "P");
        leg->AddEntry(la8_v2, "#Lambda / #bar{#Lambda}", "P");
        leg->AddEntry(xi8_v2, "#Xi^{#pm}", "P");
        leg->AddEntry(om8_v2, "#Omega^{#pm}", "P");
        leg->Draw();

        tex->SetTextFont(62);
        tex->SetTextSize(0.05);
        tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
        tex->SetTextSize(0.045);
        //tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
        // tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
        tex->SetTextFont(42);
        tex->DrawLatex(0.40,0.24,"185 #leq N_{trk}^{offline} < 250");
        //tex->DrawLatex(0.15,0.74,"|y| < 1");
        //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");

        if(i==0) c1->Print("v2ObsRapidity.pdf");
        if(i==1) c1->Print("v2BkgRapidity.pdf");

        c2->cd();
        tex->SetTextFont(62);
        tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
        tex->SetTextSize(0.045);
        tex->SetTextFont(42);
        tex->DrawLatex(0.40,0.24,"185 #leq N_{trk}^{offline} < 250");
        leg->Draw();
        if(i==1)
        {
            ks8_v2->SetMarkerStyle(24);
            xi8_v2->SetMarkerStyle(25);
            om8_v2->SetMarkerStyle(30);
            la8_v2->SetMarkerStyle(26);
        }
        ks8_v2->Draw("P");
        la8_v2->Draw("P");
        xi8_v2->Draw("P");
        om8_v2->Draw("P");

        if(i==0)
        {
            leg_co_spec->AddEntry(ks8_v2, "K_{S}^{0}", "P");
            leg_co_spec->AddEntry(la8_v2, "#Lambda/#bar{#Lambda}", "P");
            leg_co_spec->AddEntry(xi8_v2, "#Xi^{#pm}", "P");
            leg_co_spec->AddEntry(om8_v2, "#Omega^{#pm}", "P");
            TGraphErrors* ks_clone = (TGraphErrors*)ks8_v2->Clone("ks8_v2_obs");
            ks_clone->SetMarkerStyle(20);
            ks_clone->SetMarkerColor(kBlack);
            legendLabel = "Obs";
            leg_co_Label->AddEntry(ks_clone,legendLabel.c_str(), "P");
            leg_co_Label->Draw();
        }
        else
        {
            legendLabel = "Bkg";
            TGraphErrors* ks_clone = (TGraphErrors*)ks8_v2->Clone("ks8_v2_bkg");
            ks_clone->SetMarkerStyle(24);
            ks_clone->SetMarkerColor(kBlack);
            leg_co_Label->AddEntry(ks_clone,legendLabel.c_str(), "P");
            leg_co_Label->Draw();
        }

        c1->cd();

        // Write points into rootfile
        //TFile out(v2RootFileName.c_str(),"UPDATE");
        //ks8_v2->Write(kshortv2.c_str(),TObject::kOverwrite);
        //la8_v2->Write(lambdav2.c_str(),TObject::kOverwrite);
        //xi8_v2->Write(cascadev2.c_str(),TObject::kOverwrite);
        //out.Close();
    }

    c2->Print("v2ObsBkgRapidity.pdf");
}

void Rap_v2obsbkg_PbPb()
{
    TFile* f = new TFile("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/v2valuesRapidityPbPb.root");
    MITStyle();
    TCanvas* c1 = MakeCanvas("c1", "Individual");
    c1->SetLeftMargin(0.12);

    TCanvas* c2 = MakeCanvas("c2", "Combined");
    c2->cd();
    c2->SetLeftMargin(0.12);

    // draw the frame using a histogram frame
    TH1F* frame_co = c2->DrawFrame(0,-0.05,9,0.75);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_co->GetXaxis()->CenterTitle(1);
    frame_co->GetYaxis()->CenterTitle(1);
    frame_co->GetXaxis()->SetTitleSize(0.05);
    frame_co->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_co->GetYaxis()->SetTitle("v_{2}^{Obs,Bkg}");
    frame_co->GetYaxis()->SetTitleSize(0.05);
    frame_co->SetTitleOffset(1.1,"Y");
    frame_co->SetTitleOffset(1.2,"X");

    TLegend* leg_co_spec = new TLegend(0.15,0.55,0.27,0.75);
    leg_co_spec->SetFillColor(10);
    leg_co_spec->SetFillStyle(0);
    leg_co_spec->SetBorderSize(0);
    leg_co_spec->SetTextFont(42);
    leg_co_spec->SetTextSize(0.05);

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(62);
    tex->SetTextSize(0.05);

    std::string legendLabel = "";
    TLegend* leg_co_Label = new TLegend(0.75,0.75,0.85,0.85);
    leg_co_Label->SetFillColor(10);
    leg_co_Label->SetFillStyle(0);
    leg_co_Label->SetBorderSize(0);
    leg_co_Label->SetTextFont(42);
    leg_co_Label->SetTextSize(0.05);

    c1->cd();

    std::string kshortv2 = "";
    std::string lambdav2 = "";
    std::string cascadev2 = "";
    std::string yaxis = "";

    for(int i=0; i<2; i++)
    {
        TGraphErrors* xi8_v2;
        TGraphErrors* om8_v2;
        TGraphErrors* ks8_v2;
        TGraphErrors* la8_v2;

        if(i==0)
        {
            ks8_v2 = (TGraphErrors*)f->Get("v2obskshort");
            la8_v2 = (TGraphErrors*)f->Get("v2obslambda");
            xi8_v2 = (TGraphErrors*)f->Get("v2obsxi");
            om8_v2 = (TGraphErrors*)f->Get("v2obsomega");

            yaxis = "v_{2}^{obs}";
        }
        else
        {
            ks8_v2 = (TGraphErrors*)f->Get("v2bkgkshort");
            la8_v2 = (TGraphErrors*)f->Get("v2bkglambda");
            xi8_v2 = (TGraphErrors*)f->Get("v2bkgxi");
            om8_v2 = (TGraphErrors*)f->Get("v2bkgomega");

            yaxis = "v_{2}^{bkg}";

            leg_co_spec->Draw();
            leg_co_Label->Draw();
            tex->SetTextFont(62);
            tex->DrawLatex(0.15,0.8,"CMS PbPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
            tex->SetTextSize(0.045);
            tex->SetTextFont(42);
            tex->DrawLatex(0.40,0.24,"Cent: 30 - 50%");
            c1->cd();
        }

        ks8_v2->SetMarkerColor(kRed);
        ks8_v2->SetMarkerStyle(20);
        ks8_v2->SetMarkerSize(1.5);
        ks8_v2->SetLineColor(kRed);

        xi8_v2->SetMarkerColor(kGreen+2);
        xi8_v2->SetMarkerStyle(21);
        xi8_v2->SetMarkerSize(1.5);
        xi8_v2->SetLineColor(kGreen+2);

        om8_v2->SetMarkerColor(kMagenta);
        om8_v2->SetMarkerStyle(29);
        om8_v2->SetMarkerSize(1.5);
        om8_v2->SetLineColor(kMagenta);

        la8_v2->SetMarkerColor(kBlue-4);
        la8_v2->SetMarkerStyle(22);
        la8_v2->SetMarkerSize(1.5);
        la8_v2->SetLineColor(kBlue-4);

        // draw the frame using a histogram frame
        TH1F* frame = c1->DrawFrame(0,-0.05,9,0.75);
        //TH1F* frame = c1->DrawFrame(0,0.01,20,1);
        gPad->SetTickx();
        gPad->SetTicky();
        frame->GetXaxis()->CenterTitle(1);
        frame->GetYaxis()->CenterTitle(1);
        frame->GetXaxis()->SetTitleSize(0.05);
        frame->GetXaxis()->SetTitle("p_{T} (GeV)");
        frame->GetYaxis()->SetTitle(yaxis.c_str());
        frame->GetYaxis()->SetTitleSize(0.05);
        frame->SetTitleOffset(1.1,"Y");
        frame->SetTitleOffset(1.2,"X");

        //ha_v2->Draw("PESAME");
        ks8_v2->Draw("P");
        la8_v2->Draw("P");
        xi8_v2->Draw("P");
        om8_v2->Draw("P");

        TLegend* leg = new TLegend(0.15,0.55,0.27,0.75);
        leg->SetFillColor(10);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.05);
        leg->AddEntry(ks8_v2, "K_{S}^{0}", "P");
        leg->AddEntry(la8_v2, "#Lambda / #bar{#Lambda}", "P");
        leg->AddEntry(xi8_v2, "#Xi^{#pm}", "P");
        leg->AddEntry(om8_v2, "#Omega^{#pm}", "P");
        leg->Draw();

        tex->SetTextFont(62);
        tex->SetTextSize(0.05);
        tex->DrawLatex(0.15,0.8,"CMS PbPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
        tex->SetTextSize(0.045);
        //tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
        // tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
        tex->SetTextFont(42);
        tex->DrawLatex(0.37,0.72,"Cent: 30 - 50%");
        //tex->DrawLatex(0.15,0.74,"|y| < 1");
        //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");

        if(i==0) c1->Print("v2ObsRapidityPbPb.pdf");
        if(i==1) c1->Print("v2BkgRapidityPbPb.pdf");

        c2->cd();
        tex->SetTextFont(62);
        tex->DrawLatex(0.15,0.8,"CMS PbPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
        tex->SetTextSize(0.045);
        tex->SetTextFont(42);
        tex->DrawLatex(0.34,0.72,"Cent: 30 - 50%");
        leg->Draw();
        if(i==1)
        {
            ks8_v2->SetMarkerStyle(24);
            xi8_v2->SetMarkerStyle(25);
            om8_v2->SetMarkerStyle(30);
            la8_v2->SetMarkerStyle(26);
        }
        ks8_v2->Draw("P");
        la8_v2->Draw("P");
        xi8_v2->Draw("P");
        om8_v2->Draw("P");

        if(i==0)
        {
            leg_co_spec->AddEntry(ks8_v2, "K_{S}^{0}", "P");
            leg_co_spec->AddEntry(la8_v2, "#Lambda/#bar{#Lambda}", "P");
            leg_co_spec->AddEntry(xi8_v2, "#Xi^{#pm}", "P");
            leg_co_spec->AddEntry(om8_v2, "#Omega^{#pm}", "P");
            TGraphErrors* ks_clone = (TGraphErrors*)ks8_v2->Clone("ks8_v2_obs");
            ks_clone->SetMarkerStyle(20);
            ks_clone->SetMarkerColor(kBlack);
            legendLabel = "Obs";
            leg_co_Label->AddEntry(ks_clone,legendLabel.c_str(), "P");
            leg_co_Label->Draw();
        }
        else
        {
            legendLabel = "Bkg";
            TGraphErrors* ks_clone = (TGraphErrors*)ks8_v2->Clone("ks8_v2_bkg");
            ks_clone->SetMarkerStyle(24);
            ks_clone->SetMarkerColor(kBlack);
            leg_co_Label->AddEntry(ks_clone,legendLabel.c_str(), "P");
            leg_co_Label->Draw();
        }

        c1->cd();

        // Write points into rootfile
        //TFile out(v2RootFileName.c_str(),"UPDATE");
        //ks8_v2->Write(kshortv2.c_str(),TObject::kOverwrite);
        //la8_v2->Write(lambdav2.c_str(),TObject::kOverwrite);
        //xi8_v2->Write(cascadev2.c_str(),TObject::kOverwrite);
        //out.Close();
    }

    c2->Print("v2ObsBkgRapidityPbPb.pdf");
}

void Rap_fsig_pPb()
{
    MITStyle();

    TFile* f = new TFile("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/v2valuesRapidity_pPb_185_250_10_31_17.root");
    TGraphErrors* v2xi     = (TGraphErrors*)f->Get("v2xi");
    TGraphErrors* v2kshort = (TGraphErrors*)f->Get("v2kshort");
    TGraphErrors* v2lambda = (TGraphErrors*)f->Get("v2lambda");
    TGraphErrors* v2omega  = (TGraphErrors*)f->Get("v2omega");


    std::vector<double> fsig_Xi8  = {0.959427 ,0.976239 ,0.979161 ,0.980678 ,0.980661 ,0.981534 ,0.981502 ,0.979289 ,0.979192};
    double* pTXi8  = v2xi->GetX();
    const int xi_npoints = fsig_Xi8.size();

    std::vector<double> fsig_Om8 = {0.721825, 0.807658, 0.861617, 0.908762, 0.937818, 0.970699, 0.966515, 0.96482};
    double* pTOm8 = v2omega->GetX();
    const int om_npoints = fsig_Om8.size();

    std::vector<double> fsig_Ks8  = {0.999476 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,0.999999 ,0.999999 ,0.999988 ,0.999997 ,0.992327};
    double* pTKs8  = v2kshort->GetX();
	const int ks_npoints = fsig_Ks8.size();

    std::vector<double> fsig_La8 = {0.99882 ,0.999987 ,1 ,0.999524 ,0.999632 ,0.999855 ,0.999698 ,0.998783 ,0.999771 ,0.997088};
    double* pTLa8  = v2lambda->GetX();
	const int la_npoints = fsig_La8.size();

    // Pull TGraph for Kshort and lambda

    TGraphErrors* xi8_fsig_ = new TGraphErrors(xi_npoints,pTXi8,&fsig_Xi8[0],0,0);
    TGraphErrors* om8_fsig_ = new TGraphErrors(om_npoints,pTOm8,&fsig_Om8[0],0,0);
    TGraphErrors* ks8_fsig_ = new TGraphErrors(ks_npoints,pTKs8,&fsig_Ks8[0],0,0);
	TGraphErrors* la8_fsig_ = new TGraphErrors(la_npoints,pTLa8,&fsig_La8[0],0,0);

    ks8_fsig_->SetMarkerColor(kRed);
    ks8_fsig_->SetMarkerStyle(20);
    ks8_fsig_->SetMarkerSize(1.5);
    ks8_fsig_->SetLineColor(kRed);

    xi8_fsig_->SetMarkerColor(kGreen+2);
    xi8_fsig_->SetMarkerStyle(21);
    xi8_fsig_->SetMarkerSize(1.5);
    xi8_fsig_->SetLineColor(kGreen+2);

    om8_fsig_->SetMarkerColor(kMagenta);
    om8_fsig_->SetMarkerStyle(29);
    om8_fsig_->SetMarkerSize(1.5);
    om8_fsig_->SetLineColor(kMagenta);

    la8_fsig_->SetMarkerColor(kBlue-4);
    la8_fsig_->SetMarkerStyle(22);
    la8_fsig_->SetMarkerSize(1.5);
    la8_fsig_->SetLineColor(kBlue-4);

    TCanvas* c1 = MakeCanvas("c1", "Plot");
    c1->cd();
    /*c1->SetLogy();*/
    c1->SetLeftMargin(0.12);

    // draw the frame using a histogram frame

    TH1F* frame = c1->DrawFrame(0,0.65,9,1.25);
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
    om8_fsig_->Draw("P");

    // Write points into rootfile
    TFile out("8TeVfsig_GraphpPb185-250.root","RECREATE");
    ks8_fsig_->Write("kshortfsig");
    la8_fsig_->Write("lambdafsig");
    xi8_fsig_->Write("cascadefsig");
    om8_fsig_->Write("omegafsig");

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
    leg->AddEntry(xi8_fsig_, "#Xi^{#pm}", "P");
    leg->AddEntry(om8_fsig_, "#Omega^{#pm}", "P");
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
    tex->DrawLatex(0.15,0.75,"185 #leq N_{trk}^{offline} < 250");
    /*tex->DrawLatex(0.15,0.74,"|y| < 1");*/
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");
    TLine* line = new TLine(0,1,9,1);
    line->SetLineStyle(2);
    line->Draw("same");

    c1->Print("fsig_Rapidity.pdf");

}
void Rap_fsig_PbPb()
{
    MITStyle();

    TFile* f = new TFile("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/v2valuesRapidityPbPb.root");
    TGraphErrors* v2xi     = (TGraphErrors*)f->Get("v2xi");
    TGraphErrors* v2kshort = (TGraphErrors*)f->Get("v2kshort");
    TGraphErrors* v2lambda = (TGraphErrors*)f->Get("v2lambda");
    TGraphErrors* v2omega  = (TGraphErrors*)f->Get("v2omega");


    std::vector<double> fsig_Xi8  = {0.869648 ,0.90817 ,0.939477 ,0.944172 ,0.940623 ,0.930704 ,0.936042 ,0.938471};
    double* pTXi8  = v2xi->GetX();
    const int xi_npoints = fsig_Xi8.size();

    std::vector<double> fsig_Om8 = {0.647487 ,0.774924 ,0.872547 ,0.911986 ,0.937495 ,0.95909};
    double* pTOm8 = v2omega->GetX();
    const int om_npoints = fsig_Om8.size();

    std::vector<double> fsig_Ks8  = {0.88371 ,0.921642 ,0.933393 ,0.944669 ,0.957631 ,0.958982 ,0.949047 ,0.947948 ,0.91845 ,0.915206 ,0.913798 ,0.90336 ,0.889312};
    double* pTKs8  = v2kshort->GetX();
	const int ks_npoints = fsig_Ks8.size();

    std::vector<double> fsig_La8 = {0.685977,0.851926 ,0.928584 ,0.944353 ,0.952072 ,0.951948 ,0.946417 ,0.931891 ,0.88348 ,0.872633};
    double* pTLa8  = v2lambda->GetX();
	const int la_npoints = fsig_La8.size();

    // Pull TGraph for Kshort and lambda

    TGraphErrors* xi8_fsig_ = new TGraphErrors(xi_npoints,pTXi8,&fsig_Xi8[0],0,0);
    TGraphErrors* om8_fsig_ = new TGraphErrors(om_npoints,pTOm8,&fsig_Om8[0],0,0);
    TGraphErrors* ks8_fsig_ = new TGraphErrors(ks_npoints,pTKs8,&fsig_Ks8[0],0,0);
	TGraphErrors* la8_fsig_ = new TGraphErrors(la_npoints,pTLa8,&fsig_La8[0],0,0);

    ks8_fsig_->SetMarkerColor(kRed);
    ks8_fsig_->SetMarkerStyle(20);
    ks8_fsig_->SetMarkerSize(1.5);
    ks8_fsig_->SetLineColor(kRed);

    xi8_fsig_->SetMarkerColor(kGreen+2);
    xi8_fsig_->SetMarkerStyle(21);
    xi8_fsig_->SetMarkerSize(1.5);
    xi8_fsig_->SetLineColor(kGreen+2);

    om8_fsig_->SetMarkerColor(kMagenta);
    om8_fsig_->SetMarkerStyle(29);
    om8_fsig_->SetMarkerSize(1.5);
    om8_fsig_->SetLineColor(kMagenta);

    la8_fsig_->SetMarkerColor(kBlue-4);
    la8_fsig_->SetMarkerStyle(22);
    la8_fsig_->SetMarkerSize(1.5);
    la8_fsig_->SetLineColor(kBlue-4);

    TCanvas* c1 = MakeCanvas("c1", "Plot");
    c1->cd();
    /*c1->SetLogy();*/
    c1->SetLeftMargin(0.12);

    // draw the frame using a histogram frame

    TH1F* frame = c1->DrawFrame(0,0.55,9,1.25);
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
    om8_fsig_->Draw("P");

    // Write points into rootfile
    TFile out("8TeVfsig_GraphPbPb_cent_30_50.root","RECREATE");
    ks8_fsig_->Write("kshortfsig");
    la8_fsig_->Write("lambdafsig");
    xi8_fsig_->Write("cascadefsig");
    om8_fsig_->Write("omegafsig");

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
    leg->AddEntry(xi8_fsig_, "#Xi^{#pm}", "P");
    leg->AddEntry(om8_fsig_, "#Omega^{#pm}", "P");
    leg->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(62);
    tex->SetTextSize(0.045);
    tex->DrawLatex(0.15,0.8,"CMS PbPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
    tex->SetTextSize(0.04);
    //tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
    // tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
    tex->SetTextFont(42);
    tex->DrawLatex(0.15,0.75,"Cent: 30 - 50%");
    /*tex->DrawLatex(0.15,0.74,"|y| < 1");*/
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");
    TLine* line = new TLine(0,1,9,1);
    line->SetLineStyle(2);
    line->Draw("same");

    c1->Print("fsig_Rapidity_PbPb.pdf");

}

void Rap_perisub()
{
    MITStyle();
    TCanvas* c1 = MakeCanvas("c1", "Plot");
    /*c1->SetLogy();*/
    c1->SetLeftMargin(0.12);

    TCanvas* c2 = MakeCanvas("c2", "Plot");
    c2->SetLeftMargin(0.12);

    TCanvas* c3 = MakeCanvas("c3", "Plot");
    c3->SetLeftMargin(0.12);

    TCanvas* c4 = new TCanvas("c4","Combined", 1200,900);
    c4->SetLeftMargin(0.12);
    c4->Divide(2,2);

    c1->cd();

    // draw the frame using a histogram frame
    TH1F* frame1;
    TH1F* frame2;
    TH1F* frame3;
    TH1F* frame4_1;
    TH1F* frame4_2;
    TH1F* frame4_3;
    TH1F* frame4_4;
    TH1F* frame5;

    frame1 = c1->DrawFrame(0,-0.01,9,0.35);
    gPad->SetTickx();
    gPad->SetTicky();
    frame1->GetXaxis()->CenterTitle(1);
    frame1->GetYaxis()->CenterTitle(1);
    frame1->GetXaxis()->SetTitleSize(0.05);
    frame1->GetYaxis()->SetTitleSize(0.05);
    frame1->SetTitleOffset(1.1,"Y");
    frame1->SetTitleOffset(1.2,"X");
    frame1->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame1->GetYaxis()->SetTitle("v_{2}^{sig}");
    const int ks_npoints = 13;
    const int la_npoints = 10;

    // Pull TGraph for Kshort and lambda
    //TFile* file_hadv2 = TFile::Open("lrgraphv2_v3_pPb_hadron_185-above.root");
    //TFile* file_pidv2 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/v2valuesRapidityPeripheralSubFix.root");
    //TFile* file_pidv2 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/V0v2perisubFix.root");
    //TFile* file_pidv2 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/V0v2perisub.root");
    //TFile* file_pidv2 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/V0v2perisubAltLongRange.root");
    //TFile* file_pidv2 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/V0v2perisubFixedZYAM.root");
    //TFile* file_pidv2 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/V0v2perisubEG1.root");
    TFile* file_pidv2 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/V0v2perisubXiFake.root");
    //TFile* file_pidv2 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/V0v2perisubXiHighNoNorm.root");

    //TGraphErrors* ha_v2 = (TGraphErrors*)file_hadv2->Get("hadronv2");
    //TGraphErrors* ha_v2 = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033pPb_v2.txt"),1,28,1.2);

    TGraphErrors* ks8_v2 = (TGraphErrors*)file_pidv2->Get("kshortv2true");
    TGraphErrors* ks8_v2sub = (TGraphErrors*)file_pidv2->Get("kshortv2truesub");
    TGraphErrors* ks8_v2Dsub = (TGraphErrors*)file_pidv2->Get("kshortv2trueDirectSub");
    TGraphErrors* la8_v2 = (TGraphErrors*)file_pidv2->Get("lambdav2true");
    TGraphErrors* la8_v2sub = (TGraphErrors*)file_pidv2->Get("lambdav2truesub");
    TGraphErrors* la8_v2Dsub = (TGraphErrors*)file_pidv2->Get("lambdav2trueDirectSub");
    TGraphErrors* xi8_v2 = (TGraphErrors*)file_pidv2->Get("xiv2true");
    TGraphErrors* xi8_v2sub = (TGraphErrors*)file_pidv2->Get("xiv2truesub");
    TGraphErrors* xi8_v2Dsub = (TGraphErrors*)file_pidv2->Get("xiv2trueDirectSub");

    TGraphErrors* DirectSubNass_ks     = (TGraphErrors*)file_pidv2->Get("DirectSubNass_ks");
    TGraphErrors* Nassoc_bkg_ks        = (TGraphErrors*)file_pidv2->Get("Nassoc_bkg_ks");
    TGraphErrors* Nassoc_obs_ks        = (TGraphErrors*)file_pidv2->Get("Nassoc_obs_ks");
    TGraphErrors* DirectSubNass_low_ks = (TGraphErrors*)file_pidv2->Get("DirectSubNass_low_ks");
    TGraphErrors* Nassoc_bkg_low_ks    = (TGraphErrors*)file_pidv2->Get("Nassoc_low_bkg_ks");
    TGraphErrors* Nassoc_obs_low_ks    = (TGraphErrors*)file_pidv2->Get("Nassoc_low_obs_ks");

    TGraphErrors* DirectSubNass_la     = (TGraphErrors*)file_pidv2->Get("DirectSubNass_la");
    TGraphErrors* Nassoc_bkg_la        = (TGraphErrors*)file_pidv2->Get("Nassoc_bkg_la");
    TGraphErrors* Nassoc_obs_la        = (TGraphErrors*)file_pidv2->Get("Nassoc_obs_la");
    TGraphErrors* DirectSubNass_low_la = (TGraphErrors*)file_pidv2->Get("DirectSubNass_low_la");
    TGraphErrors* Nassoc_bkg_low_la    = (TGraphErrors*)file_pidv2->Get("Nassoc_low_bkg_la");
    TGraphErrors* Nassoc_obs_low_la    = (TGraphErrors*)file_pidv2->Get("Nassoc_low_obs_la");

    TGraphErrors* DirectSubNass_xi     = (TGraphErrors*)file_pidv2->Get("DirectSubNass_xi");
    TGraphErrors* Nassoc_bkg_xi        = (TGraphErrors*)file_pidv2->Get("Nassoc_bkg_xi");
    TGraphErrors* Nassoc_obs_xi        = (TGraphErrors*)file_pidv2->Get("Nassoc_obs_xi");
    TGraphErrors* DirectSubNass_low_xi = (TGraphErrors*)file_pidv2->Get("DirectSubNass_low_xi");
    TGraphErrors* Nassoc_bkg_low_xi    = (TGraphErrors*)file_pidv2->Get("Nassoc_low_bkg_xi");
    TGraphErrors* Nassoc_obs_low_xi    = (TGraphErrors*)file_pidv2->Get("Nassoc_low_obs_xi");

    TGraphErrors* ks8_v2kn = (TGraphErrors*)file_pidv2->Get("kshortv2true_KET");
    TGraphErrors* ks8_v2kn_sub = (TGraphErrors*)file_pidv2->Get("kshortv2truesub_KET");
    TGraphErrors* ks8_v2kn_Dsub = (TGraphErrors*)file_pidv2->Get("kshortv2trueDirectSub_KET");
    TGraphErrors* la8_v2kn = (TGraphErrors*)file_pidv2->Get("lambdav2true_KET");
    TGraphErrors* la8_v2kn_sub = (TGraphErrors*)file_pidv2->Get("lambdav2truesub_KET");
    TGraphErrors* la8_v2kn_Dsub = (TGraphErrors*)file_pidv2->Get("lambdav2trueDirectSub_KET");
    TGraphErrors* xi8_v2kn = (TGraphErrors*)file_pidv2->Get("xiv2true_KET");
    TGraphErrors* xi8_v2kn_sub = (TGraphErrors*)file_pidv2->Get("xiv2truesub_KET");
    TGraphErrors* xi8_v2kn_Dsub = (TGraphErrors*)file_pidv2->Get("xiv2trueDirectSub_KET");

    double* ks8_v2Y = ks8_v2->GetY();
    double* ks8_v2EY = ks8_v2->GetEY();
    double* la8_v2Y = la8_v2->GetY();
    double* la8_v2EY = la8_v2->GetEY();
    double* xi8_v2Y = xi8_v2->GetY();
    double* xi8_v2EY = xi8_v2->GetEY();

    double* ks8_v2subY = ks8_v2sub->GetY();
    double* ks8_v2subEY = ks8_v2sub->GetEY();
    double* la8_v2subY = la8_v2sub->GetY();
    double* la8_v2subEY = la8_v2sub->GetEY();
    double* xi8_v2subY = xi8_v2sub->GetY();
    double* xi8_v2subEY = xi8_v2sub->GetEY();

    double* ks8_v2DsubY = ks8_v2Dsub->GetY();
    double* ks8_v2DsubEY = ks8_v2Dsub->GetEY();
    double* la8_v2DsubY = la8_v2Dsub->GetY();
    double* la8_v2DsubEY = la8_v2Dsub->GetEY();
    double* xi8_v2DsubY = xi8_v2Dsub->GetY();
    double* xi8_v2DsubEY = xi8_v2Dsub->GetEY();

    double* pt_ks = ks8_v2->GetX();
    double* pt_la = la8_v2->GetX();
    double* pt_xi = xi8_v2->GetX();


    ks8_v2->SetMarkerColor(kRed);
    ks8_v2->SetMarkerStyle(20);
    ks8_v2->SetMarkerSize(1.5);
    ks8_v2->SetLineColor(kRed);

    ks8_v2sub->SetMarkerColor(kRed);
    ks8_v2sub->SetMarkerStyle(24);
    ks8_v2sub->SetMarkerSize(1.5);
    ks8_v2sub->SetLineColor(kRed);

    ks8_v2Dsub->SetMarkerColor(kRed);
    ks8_v2Dsub->SetMarkerStyle(29);
    ks8_v2Dsub->SetMarkerSize(1.5);
    ks8_v2Dsub->SetLineColor(kRed);

    ks8_v2kn->SetMarkerColor(kRed);
    ks8_v2kn->SetMarkerStyle(20);
    ks8_v2kn->SetMarkerSize(1.5);
    ks8_v2kn->SetLineColor(kRed);

    ks8_v2kn_sub->SetMarkerColor(kRed);
    ks8_v2kn_sub->SetMarkerStyle(24);
    ks8_v2kn_sub->SetMarkerSize(1.5);
    ks8_v2kn_sub->SetLineColor(kRed);

    ks8_v2kn_Dsub->SetMarkerColor(kRed);
    ks8_v2kn_Dsub->SetMarkerStyle(29);
    ks8_v2kn_Dsub->SetMarkerSize(1.5);
    ks8_v2kn_Dsub->SetLineColor(kRed);

    la8_v2->SetMarkerColor(kBlue-4);
    la8_v2->SetMarkerStyle(20);
    la8_v2->SetMarkerSize(1.5);
    la8_v2->SetLineColor(kBlue-4);

    la8_v2sub->SetMarkerColor(kBlue-4);
    la8_v2sub->SetMarkerStyle(24);
    la8_v2sub->SetMarkerSize(1.5);
    la8_v2sub->SetLineColor(kBlue-4);

    la8_v2Dsub->SetMarkerColor(kBlue-4);
    la8_v2Dsub->SetMarkerStyle(29);
    la8_v2Dsub->SetMarkerSize(1.5);
    la8_v2Dsub->SetLineColor(kBlue-4);

    la8_v2kn->SetMarkerColor(kBlue-4);
    la8_v2kn->SetMarkerStyle(20);
    la8_v2kn->SetMarkerSize(1.5);
    la8_v2kn->SetLineColor(kBlue-4);

    la8_v2kn_sub->SetMarkerColor(kBlue-4);
    la8_v2kn_sub->SetMarkerStyle(24);
    la8_v2kn_sub->SetMarkerSize(1.5);
    la8_v2kn_sub->SetLineColor(kBlue-4);

    la8_v2kn_Dsub->SetMarkerColor(kBlue-4);
    la8_v2kn_Dsub->SetMarkerStyle(29);
    la8_v2kn_Dsub->SetMarkerSize(1.5);
    la8_v2kn_Dsub->SetLineColor(kBlue-4);

    xi8_v2->SetMarkerColor(kGreen-2);
    xi8_v2->SetMarkerStyle(20);
    xi8_v2->SetMarkerSize(1.5);
    xi8_v2->SetLineColor(kGreen-2);

    xi8_v2sub->SetMarkerColor(kGreen-2);
    xi8_v2sub->SetMarkerStyle(24);
    xi8_v2sub->SetMarkerSize(1.5);
    xi8_v2sub->SetLineColor(kGreen-2);

    xi8_v2Dsub->SetMarkerColor(kGreen-2);
    xi8_v2Dsub->SetMarkerStyle(29);
    xi8_v2Dsub->SetMarkerSize(1.5);
    xi8_v2Dsub->SetLineColor(kGreen-2);

    xi8_v2kn->SetMarkerColor(kGreen-2);
    xi8_v2kn->SetMarkerStyle(20);
    xi8_v2kn->SetMarkerSize(1.5);
    xi8_v2kn->SetLineColor(kGreen-2);

    xi8_v2kn_sub->SetMarkerColor(kGreen-2);
    xi8_v2kn_sub->SetMarkerStyle(24);
    xi8_v2kn_sub->SetMarkerSize(1.5);
    xi8_v2kn_sub->SetLineColor(kGreen-2);

    xi8_v2kn_Dsub->SetMarkerColor(kGreen-2);
    xi8_v2kn_Dsub->SetMarkerStyle(29);
    xi8_v2kn_Dsub->SetMarkerSize(1.5);
    xi8_v2kn_Dsub->SetLineColor(kGreen-2);



    c1->cd();

    TLegend* leg = new TLegend(0.15,0.60,0.27,0.85);
    leg->SetFillColor(10);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    //leg->AddEntry(ha_v2, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}", "P");
    leg->AddEntry(ks8_v2, "K_{S}^{0}", "P");
    leg->AddEntry(ks8_v2sub, "K_{S}^{0} Sub", "P");
    leg->AddEntry(ks8_v2Dsub, "K_{S}^{0} Direct Sub", "P");
    leg->AddEntry(la8_v2, "#Lambda / #bar{#Lambda}", "P");
    leg->AddEntry(la8_v2sub, "#Lambda / #bar{#Lambda} Sub", "P");
    leg->AddEntry(la8_v2Dsub, "#Lambda / #bar{#Lambda} Direct Sub", "P");
    leg->Draw();


    //ha_v2->Draw("PESAME");
    ks8_v2->Draw("P");
    la8_v2->Draw("P");
    ks8_v2sub->Draw("P");
    la8_v2sub->Draw("P");
    ks8_v2Dsub->Draw("P");
    la8_v2Dsub->Draw("P");

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(62);
    tex->SetTextSize(0.05);
    tex->DrawLatex(0.12,0.9,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->SetTextFont(42);
    tex->SetTextSize(0.04);
    tex->DrawLatex(0.40,0.82,"(185 #leq N_{trk}^{offline} < 250) - (0 < N_{trk}^{offline} < 20)");
    //tex->DrawLatex(0.36,0.72,"Cent: 30% - 50%");
    /*tex->DrawLatex(0.15,0.74,"|y| < 1");*/
    //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");

    c2->cd();

    frame2 = c2->DrawFrame(0,-0.01,7.5,0.3);
    gPad->SetTickx();
    gPad->SetTicky();
    frame2->GetXaxis()->CenterTitle(1);
    frame2->GetYaxis()->CenterTitle(1);
    frame2->GetXaxis()->SetTitleSize(0.05);
    frame2->GetYaxis()->SetTitleSize(0.05);
    frame2->SetTitleOffset(1.1,"Y");
    frame2->SetTitleOffset(1.2,"X");
    frame2->GetXaxis()->SetTitle("KE_{T} (GeV)");
    frame2->GetYaxis()->SetTitle("v_{2}^{sig}");

    TLegend* legkn = new TLegend(0.145,0.56,0.26,0.76);
    legkn->SetFillColor(10);
    legkn->SetFillStyle(0);
    legkn->SetBorderSize(0);
    legkn->SetTextFont(42);
    legkn->SetTextSize(0.05);
    legkn->AddEntry(ks8_v2kn, "K_{S}^{0}", "P");
    legkn->AddEntry(ks8_v2kn_sub, "K_{S}^{0} Sub", "P");
    legkn->AddEntry(la8_v2kn, "#Lambda / #bar{#Lambda}", "P");
    legkn->AddEntry(la8_v2kn_sub, "#Lambda / #bar{#Lambda} Sub", "P");
    legkn->Draw();

    ks8_v2kn->Draw("P");
    la8_v2kn->Draw("P");
    ks8_v2kn_sub->Draw("P");
    la8_v2kn_sub->Draw("P");

    // Draw Legend and write points into rootfile

    tex->SetTextFont(62);
    tex->SetTextSize(0.05);
    tex->DrawLatex(0.12,0.9,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->SetTextSize(0.045);
    tex->SetTextFont(42);
    tex->DrawLatex(0.46,0.73,"185 #leq N_{trk}^{offline} < 250");
    tex->DrawLatex(0.145,0.80,"Peripheral Subtraction 0 #leq N_{trk}^{offline} < 20");

    c1->Print("v2SigRapidityPeriSub.pdf");
    c1->Print("v2SigRapidityPeriSub.png");
    c2->Print("v2SigRapidityPeriSubKET.pdf");
    c2->Print("v2SigRapidityPeriSubKET.png");

    std::vector<double> ks8_v2ratio;
    std::vector<double> ks8_v2Eratio;
    std::vector<double> ks8_v2Dratio;
    std::vector<double> ks8_v2EDratio;

    std::vector<double> la8_v2ratio;
    std::vector<double> la8_v2Eratio;
    std::vector<double> la8_v2Dratio;
    std::vector<double> la8_v2EDratio;

    std::vector<double> xi8_v2ratio;
    std::vector<double> xi8_v2Eratio;
    std::vector<double> xi8_v2Dratio;
    std::vector<double> xi8_v2EDratio;

    //Ratio
    for(int i=0; i<ks8_v2->GetN(); i++)
    {
        ks8_v2ratio.push_back(ks8_v2subY[i]/ks8_v2Y[i]);
        ks8_v2Eratio.push_back(ks8_v2subEY[i]/ks8_v2EY[i]);
        ks8_v2Dratio.push_back(ks8_v2DsubY[i]/ks8_v2Y[i]);
        ks8_v2EDratio.push_back(ks8_v2DsubEY[i]/ks8_v2EY[i]);
    }
    for(int i=0; i<la8_v2->GetN(); i++)
    {
        la8_v2ratio.push_back(la8_v2subY[i]/la8_v2Y[i]);
        la8_v2Eratio.push_back(la8_v2subEY[i]/la8_v2Y[i]);
        la8_v2Dratio.push_back(la8_v2DsubY[i]/la8_v2Y[i]);
        la8_v2EDratio.push_back(la8_v2DsubEY[i]/la8_v2EY[i]);
    }
    for(int i=0; i<xi8_v2->GetN(); i++)
    {
        xi8_v2ratio.push_back(xi8_v2subY[i]/xi8_v2Y[i]);
        xi8_v2Eratio.push_back(xi8_v2subEY[i]/xi8_v2Y[i]);
        xi8_v2Dratio.push_back(xi8_v2DsubY[i]/xi8_v2Y[i]);
        xi8_v2EDratio.push_back(xi8_v2DsubEY[i]/xi8_v2EY[i]);
    }

    TGraphErrors* ks8_v2r = new TGraphErrors(ks8_v2->GetN(),pt_ks,&ks8_v2ratio[0],0,0);
    TGraphErrors* ks8_v2Dr = new TGraphErrors(ks8_v2->GetN(),pt_ks,&ks8_v2Dratio[0],0,0);
    TGraphErrors* la8_v2r = new TGraphErrors(la8_v2->GetN(),pt_la,&la8_v2ratio[0],0,0);
    TGraphErrors* la8_v2Dr = new TGraphErrors(la8_v2->GetN(),pt_la,&la8_v2Dratio[0],0,0);
    TGraphErrors* xi8_v2r = new TGraphErrors(xi8_v2->GetN(),pt_xi,&xi8_v2ratio[0],0,0);
    TGraphErrors* xi8_v2Dr = new TGraphErrors(xi8_v2->GetN(),pt_xi,&xi8_v2Dratio[0],0,0);

    c3->cd();

    //Ratios
    frame3 = c3->DrawFrame(0,0.05,9,1.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame3->GetXaxis()->CenterTitle(1);
    frame3->GetYaxis()->CenterTitle(1);
    frame3->GetXaxis()->SetTitleSize(0.05);
    frame3->GetYaxis()->SetTitleSize(0.05);
    frame3->SetTitleOffset(1.1,"Y");
    frame3->SetTitleOffset(1.2,"X");
    frame3->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame3->GetYaxis()->SetTitle("v_{2}^{#left{sub#right}}/v_{2}^{sig}");

    ks8_v2r->SetMarkerColor(kRed);
    ks8_v2r->SetMarkerStyle(20);
    ks8_v2r->SetMarkerSize(1.5);
    ks8_v2r->SetLineColor(kRed);

    ks8_v2Dr->SetMarkerColor(kRed);
    ks8_v2Dr->SetMarkerStyle(24);
    ks8_v2Dr->SetMarkerSize(1.5);
    ks8_v2Dr->SetLineColor(kRed);

    la8_v2r->SetMarkerColor(kBlue-4);
    la8_v2r->SetMarkerStyle(20);
    la8_v2r->SetMarkerSize(1.5);
    la8_v2r->SetLineColor(kBlue-4);

    la8_v2Dr->SetMarkerColor(kBlue-4);
    la8_v2Dr->SetMarkerStyle(24);
    la8_v2Dr->SetMarkerSize(1.5);
    la8_v2Dr->SetLineColor(kBlue-4);

    xi8_v2r->SetMarkerColor(kGreen-2);
    xi8_v2r->SetMarkerStyle(20);
    xi8_v2r->SetMarkerSize(1.5);
    xi8_v2r->SetLineColor(kGreen-2);

    xi8_v2Dr->SetMarkerColor(kGreen-2);
    xi8_v2Dr->SetMarkerStyle(24);
    xi8_v2Dr->SetMarkerSize(1.5);
    xi8_v2Dr->SetLineColor(kGreen-2);


    TLegend* legr = new TLegend(0.60,0.68,0.75,0.88);
    legr->SetFillColor(10);
    legr->SetFillStyle(0);
    legr->SetBorderSize(0);
    legr->SetTextFont(42);
    legr->SetTextSize(0.04);
    legr->AddEntry(ks8_v2r, "K_{S}^{0} Sub", "P");
    legr->AddEntry(ks8_v2Dr, "K_{S}^{0} Direct Sub", "P");
    legr->AddEntry(la8_v2r, "#Lambda / #bar{#Lambda} Sub", "P");
    legr->AddEntry(la8_v2Dr, "#Lambda / #bar{#Lambda} Direct Sub", "P");
    legr->AddEntry(xi8_v2r, "#Xi^{-} Sub", "P");
    legr->AddEntry(xi8_v2Dr, "#Xi^{-} Direct Sub", "P");
    legr->Draw();

    TLine* line = new TLine(0,1,9,1);
    line->Draw();

    ks8_v2r->Draw("P");
    ks8_v2Dr->Draw("P");
    la8_v2r->Draw("P");
    la8_v2Dr->Draw("P");
    xi8_v2r->Draw("P");
    xi8_v2Dr->Draw("P");

    c3->Print("v2SigRapidityPeriSubRatio.pdf");
    c3->Print("v2SigRapidityPeriSubRatio.png");

    frame4_1 = c4->cd(1)->DrawFrame(0,-0.01,9,0.3);
    gPad->SetTickx();
    gPad->SetTicky();
    frame4_1->GetXaxis()->CenterTitle(1);
    frame4_1->GetYaxis()->CenterTitle(1);
    frame4_1->GetXaxis()->SetTitleSize(0.05);
    frame4_1->GetYaxis()->SetTitleSize(0.05);
    frame4_1->SetTitleOffset(1.1,"Y");
    frame4_1->SetTitleOffset(1.2,"X");
    frame4_1->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame4_1->GetYaxis()->SetTitle("v_{2}^{sig}");

    TLegend* leg4_1 = new TLegend(0.20,0.65,0.34,0.90);
    leg4_1->SetFillColor(10);
    leg4_1->SetFillStyle(0);
    leg4_1->SetBorderSize(0);
    leg4_1->SetTextFont(42);
    leg4_1->SetTextSize(0.04);
    leg4_1->AddEntry(ks8_v2, "K_{S}^{0}", "P");
    leg4_1->AddEntry(ks8_v2sub, "K_{S}^{0} Sub", "P");
    leg4_1->AddEntry(ks8_v2Dsub, "K_{S}^{0} Direct Sub", "P");
    leg4_1->Draw();
    ks8_v2->Draw("P");
    ks8_v2sub->Draw("P");
    ks8_v2Dsub->Draw("P");

    frame4_2 = c4->cd(2)->DrawFrame(0,-0.01,9,0.3);
    gPad->SetTickx();
    gPad->SetTicky();
    frame4_2->GetXaxis()->CenterTitle(1);
    frame4_2->GetYaxis()->CenterTitle(1);
    frame4_2->GetXaxis()->SetTitleSize(0.05);
    frame4_2->GetYaxis()->SetTitleSize(0.05);
    frame4_2->SetTitleOffset(1.1,"Y");
    frame4_2->SetTitleOffset(1.2,"X");
    frame4_2->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame4_2->GetYaxis()->SetTitle("v_{2}^{sig}");

    TLegend* leg4_2 = new TLegend(0.20,0.65,0.34,0.90);
    leg4_2->SetFillColor(10);
    leg4_2->SetFillStyle(0);
    leg4_2->SetBorderSize(0);
    leg4_2->SetTextFont(42);
    leg4_2->SetTextSize(0.04);
    leg4_2->AddEntry(la8_v2, "#Lambda / #bar{#Lambda}", "P");
    leg4_2->AddEntry(la8_v2sub, "#Lambda / #bar{#Lambda} Sub", "P");
    leg4_2->AddEntry(la8_v2Dsub, "#Lambda / #bar{#Lambda} Direct Sub", "P");
    leg4_2->Draw();
    la8_v2->Draw("P");
    la8_v2sub->Draw("P");
    la8_v2Dsub->Draw("P");

    frame4_3 = c4->cd(3)->DrawFrame(0,-0.01,9,0.3);
    gPad->SetTickx();
    gPad->SetTicky();
    frame4_3->GetXaxis()->CenterTitle(1);
    frame4_3->GetYaxis()->CenterTitle(1);
    frame4_3->GetXaxis()->SetTitleSize(0.05);
    frame4_3->GetYaxis()->SetTitleSize(0.05);
    frame4_3->SetTitleOffset(1.1,"Y");
    frame4_3->SetTitleOffset(1.2,"X");
    frame4_3->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame4_3->GetYaxis()->SetTitle("v_{2}^{sig}");

    TLegend* leg4_3 = new TLegend(0.20,0.65,0.34,0.90);
    leg4_3->SetFillColor(10);
    leg4_3->SetFillStyle(0);
    leg4_3->SetBorderSize(0);
    leg4_3->SetTextFont(42);
    leg4_3->SetTextSize(0.04);
    leg4_3->AddEntry(xi8_v2, "#Xi^{-}", "P");
    leg4_3->AddEntry(xi8_v2sub, "#Xi^{-} Sub", "P");
    leg4_3->AddEntry(xi8_v2Dsub, "Xi^{-} Direct Sub", "P");
    leg4_3->Draw();
    xi8_v2->Draw("P");
    xi8_v2sub->Draw("P");
    xi8_v2Dsub->Draw("P");

    frame4_4 = c4->cd(4)->DrawFrame(0,0.05,9,1.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame4_3->GetXaxis()->CenterTitle(1);
    frame4_3->GetYaxis()->CenterTitle(1);
    frame4_3->GetXaxis()->SetTitleSize(0.05);
    frame4_3->GetYaxis()->SetTitleSize(0.05);
    frame4_3->SetTitleOffset(1.1,"Y");
    frame4_3->SetTitleOffset(1.2,"X");
    frame4_3->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame4_3->GetYaxis()->SetTitle("v_{2}^{sig}");

    legr->Draw();

    line->Draw();

    ks8_v2r->Draw("P");
    ks8_v2Dr->Draw("P");
    la8_v2r->Draw("P");
    la8_v2Dr->Draw("P");
    xi8_v2r->Draw("P");
    xi8_v2Dr->Draw("P");

    c4->Print("v2SigRapidityPeriSubCombined.pdf");
    c4->Print("v2SigRapidityPeriSubCombined.png");

    c3->cd();

    frame5 = c3->DrawFrame(0,-0.01,9,0.3);
    gPad->SetTickx();
    gPad->SetTicky();
    frame5->GetXaxis()->CenterTitle(1);
    frame5->GetYaxis()->CenterTitle(1);
    frame5->GetXaxis()->SetTitleSize(0.05);
    frame5->GetYaxis()->SetTitleSize(0.05);
    frame5->SetTitleOffset(1.1,"Y");
    frame5->SetTitleOffset(1.2,"X");
    frame5->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame5->GetYaxis()->SetTitle("v_{2}^{#left{sub#right}}");

    TGraphErrors* ks8_v2Dsub_co = (TGraphErrors*)ks8_v2Dsub->Clone();
    TGraphErrors* la8_v2Dsub_co = (TGraphErrors*)la8_v2Dsub->Clone();
    TGraphErrors* xi8_v2Dsub_co = (TGraphErrors*)xi8_v2Dsub->Clone();

    ks8_v2Dsub_co->SetMarkerStyle(20);
    la8_v2Dsub_co->SetMarkerStyle(22);
    xi8_v2Dsub_co->SetMarkerStyle(21);

    TLegend* leg5 = new TLegend(0.16,0.6,0.3,0.86);
    leg5->SetFillColor(10);
    leg5->SetFillStyle(0);
    leg5->SetBorderSize(0);
    leg5->SetTextFont(42);
    leg5->SetTextSize(0.04);
    leg5->AddEntry(ks8_v2Dsub_co, "K_{S}^{0} Direct Sub", "P");
    leg5->AddEntry(la8_v2Dsub_co, "#Lambda / #bar{#Lambda} Direct Sub", "P");
    leg5->AddEntry(xi8_v2Dsub_co, "Xi^{-} Direct Sub", "P");
    leg5->Draw();
    ks8_v2Dsub_co->Draw("P");
    la8_v2Dsub_co->Draw("P");
    xi8_v2Dsub_co->Draw("P");

    c3->Print("v2SigRapiditySubtracted.pdf");
    c3->Print("v2SigRapiditySubtracted.png");

}
void RapSys_perisub()
{
    TFile* file_pidv2 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/V0v2perisubFixedZYAM.root");

    //Nass
    TGraphErrors* DirectSubNass_ks     = (TGraphErrors*)file_pidv2->Get("DirectSubNass_ks");
    TGraphErrors* Nassoc_bkg_ks        = (TGraphErrors*)file_pidv2->Get("Nassoc_bkg_ks");
    TGraphErrors* Nassoc_obs_ks        = (TGraphErrors*)file_pidv2->Get("Nassoc_obs_ks");
    TGraphErrors* DirectSubNass_low_ks = (TGraphErrors*)file_pidv2->Get("DirectSubNass_low_ks");
    TGraphErrors* Nassoc_bkg_low_ks    = (TGraphErrors*)file_pidv2->Get("Nassoc_low_bkg_ks");
    TGraphErrors* Nassoc_obs_low_ks    = (TGraphErrors*)file_pidv2->Get("Nassoc_low_obs_ks");

    TGraphErrors* DirectSubNass_la     = (TGraphErrors*)file_pidv2->Get("DirectSubNass_la");
    TGraphErrors* Nassoc_bkg_la        = (TGraphErrors*)file_pidv2->Get("Nassoc_bkg_la");
    TGraphErrors* Nassoc_obs_la        = (TGraphErrors*)file_pidv2->Get("Nassoc_obs_la");
    TGraphErrors* DirectSubNass_low_la = (TGraphErrors*)file_pidv2->Get("DirectSubNass_low_la");
    TGraphErrors* Nassoc_bkg_low_la    = (TGraphErrors*)file_pidv2->Get("Nassoc_low_bkg_la");
    TGraphErrors* Nassoc_obs_low_la    = (TGraphErrors*)file_pidv2->Get("Nassoc_low_obs_la");

    TGraphErrors* DirectSubNass_xi     = (TGraphErrors*)file_pidv2->Get("DirectSubNass_xi");
    TGraphErrors* Nassoc_bkg_xi        = (TGraphErrors*)file_pidv2->Get("Nassoc_bkg_xi");
    TGraphErrors* Nassoc_obs_xi        = (TGraphErrors*)file_pidv2->Get("Nassoc_obs_xi");
    TGraphErrors* DirectSubNass_low_xi = (TGraphErrors*)file_pidv2->Get("DirectSubNass_low_xi");
    TGraphErrors* Nassoc_bkg_low_xi    = (TGraphErrors*)file_pidv2->Get("Nassoc_low_bkg_xi");
    TGraphErrors* Nassoc_obs_low_xi    = (TGraphErrors*)file_pidv2->Get("Nassoc_low_obs_xi");

    double* Nassoc_bkg_ks_Y = Nassoc_bkg_ks->GetY();
    double* Nassoc_obs_ks_Y = Nassoc_obs_ks->GetY();
    double* Nassoc_bkg_low_ks_Y = Nassoc_bkg_low_ks->GetY();
    double* Nassoc_obs_low_ks_Y = Nassoc_obs_low_ks->GetY();

    double* Nassoc_bkg_la_Y = Nassoc_bkg_la->GetY();
    double* Nassoc_obs_la_Y = Nassoc_obs_la->GetY();
    double* Nassoc_bkg_low_la_Y = Nassoc_bkg_low_la->GetY();
    double* Nassoc_obs_low_la_Y = Nassoc_obs_low_la->GetY();

    double* pt_ks = Nassoc_obs_ks->GetX();
    double* pt_la = Nassoc_obs_la->GetX();

    std::vector<double> Nassoc_high_ratio_ks;
    std::vector<double> Nassoc_high_err_ratio_ks;
    std::vector<double> Nassoc_low_ratio_ks;
    std::vector<double> Nassoc_low_err_ratio_ks;

    std::vector<double> Nassoc_high_ratio_la;
    std::vector<double> Nassoc_high_err_ratio_la;
    std::vector<double> Nassoc_low_ratio_la;
    std::vector<double> Nassoc_low_err_ratio_la;

    int numKs = Nassoc_obs_ks->GetN();
    int numLa = Nassoc_obs_la->GetN();

    //Do obs/bkg
    for(int i=0; i<Nassoc_bkg_ks->GetN(); i++)
    {
        Nassoc_high_ratio_ks.push_back(Nassoc_obs_ks_Y[i]/Nassoc_bkg_ks_Y[i]);

        Nassoc_low_ratio_ks.push_back(Nassoc_obs_low_ks_Y[i]/Nassoc_bkg_low_ks_Y[i]);
    }

    for(int i=0; i<Nassoc_bkg_la->GetN(); i++)
    {
        Nassoc_high_ratio_la.push_back(Nassoc_obs_la_Y[i]/Nassoc_bkg_la_Y[i]);

        Nassoc_low_ratio_la.push_back(Nassoc_obs_low_la_Y[i]/Nassoc_bkg_low_la_Y[i]);
    }

    //TGraphErrors* TGNassoc_high_ratio_ks = new TGraphErrors(numKs,pt_ks,&Nassoc_high_ratio_ks[0],0,0);
    //TGraphErrors* TGNassoc_low_ratio_ks = new TGraphErrors(numKs,pt_ks,&Nassoc_low_ratio_ks[0],0,0);
    //TGraphErrors* TGNassoc_high_ratio_la = new TGraphErrors(numLa,pt_la,&Nassoc_high_ratio_la[0],0,0);
    //TGraphErrors* TGNassoc_low_ratio_la = new TGraphErrors(numLa,pt_la,&Nassoc_low_ratio_la[0],0,0);

    TGraphErrors* TGNassoc_high_ratio_ks = TGDivideSameX(Nassoc_obs_ks,Nassoc_bkg_ks);
    TGraphErrors* TGNassoc_low_ratio_ks = TGDivideSameX(Nassoc_obs_low_ks,Nassoc_bkg_low_ks);
    TGraphErrors* TGNassoc_high_ratio_la = TGDivideSameX(Nassoc_obs_la,Nassoc_bkg_la);
    TGraphErrors* TGNassoc_low_ratio_la = TGDivideSameX(Nassoc_obs_low_la,Nassoc_bkg_low_la);

    TGNassoc_high_ratio_ks->SetMarkerColor(kRed);
    TGNassoc_high_ratio_ks->SetMarkerStyle(20);
    TGNassoc_high_ratio_ks->SetMarkerSize(1.2);
    TGNassoc_high_ratio_ks->SetLineColor(kRed);

    TGNassoc_low_ratio_ks->SetMarkerColor(kBlue);
    TGNassoc_low_ratio_ks->SetMarkerStyle(20);
    TGNassoc_low_ratio_ks->SetMarkerSize(1.2);
    TGNassoc_low_ratio_ks->SetLineColor(kBlue);

    TGNassoc_high_ratio_la->SetMarkerColor(kRed);
    TGNassoc_high_ratio_la->SetMarkerStyle(22);
    TGNassoc_high_ratio_la->SetMarkerSize(1.2);
    TGNassoc_high_ratio_la->SetLineColor(kRed);

    TGNassoc_low_ratio_la->SetMarkerColor(kBlue);
    TGNassoc_low_ratio_la->SetMarkerStyle(22);
    TGNassoc_low_ratio_la->SetMarkerSize(1.2);
    TGNassoc_low_ratio_la->SetLineColor(kBlue);

    DirectSubNass_ks     ->SetMarkerColor(kGreen-2);
    DirectSubNass_ks     ->SetMarkerStyle(24);
    DirectSubNass_ks     ->SetMarkerSize(1.5);
    DirectSubNass_ks     ->SetLineColor(kGreen-2);
    Nassoc_bkg_ks        ->SetMarkerColor(kBlue);
    Nassoc_bkg_ks        ->SetMarkerStyle(24);
    Nassoc_bkg_ks        ->SetMarkerSize(1.5);
    Nassoc_bkg_ks        ->SetLineColor(kBlue);
    Nassoc_obs_ks        ->SetMarkerColor(kRed);
    Nassoc_obs_ks        ->SetMarkerStyle(24);
    Nassoc_obs_ks        ->SetMarkerSize(1.5);
    Nassoc_obs_ks        ->SetLineColor(kRed);

    DirectSubNass_low_ks ->SetMarkerColor(kGreen-2);
    DirectSubNass_low_ks ->SetMarkerStyle(24);
    DirectSubNass_low_ks ->SetMarkerSize(1.5);
    DirectSubNass_low_ks ->SetLineColor(kGreen-2);
    Nassoc_bkg_low_ks    ->SetMarkerColor(kBlue);
    Nassoc_bkg_low_ks    ->SetMarkerStyle(24);
    Nassoc_bkg_low_ks    ->SetMarkerSize(1.5);
    Nassoc_bkg_low_ks    ->SetLineColor(kBlue);
    Nassoc_obs_low_ks    ->SetMarkerColor(kRed);
    Nassoc_obs_low_ks    ->SetMarkerStyle(24);
    Nassoc_obs_low_ks    ->SetMarkerSize(1.5);
    Nassoc_obs_low_ks    ->SetLineColor(kRed);

    DirectSubNass_la     ->SetMarkerColor(kGreen-2);
    DirectSubNass_la     ->SetMarkerStyle(26);
    DirectSubNass_la     ->SetMarkerSize(1.5);
    DirectSubNass_la     ->SetLineColor(kGreen-2);
    Nassoc_bkg_la        ->SetMarkerColor(kBlue);
    Nassoc_bkg_la        ->SetMarkerStyle(26);
    Nassoc_bkg_la        ->SetMarkerSize(1.5);
    Nassoc_bkg_la        ->SetLineColor(kBlue);
    Nassoc_obs_la        ->SetMarkerColor(kRed);
    Nassoc_obs_la        ->SetMarkerStyle(26);
    Nassoc_obs_la        ->SetMarkerSize(1.5);
    Nassoc_obs_la        ->SetLineColor(kRed);

    DirectSubNass_low_la ->SetMarkerColor(kGreen-2);
    DirectSubNass_low_la ->SetMarkerStyle(26);
    DirectSubNass_low_la ->SetMarkerSize(1.5);
    DirectSubNass_low_la ->SetLineColor(kGreen-2);
    Nassoc_bkg_low_la    ->SetMarkerColor(kBlue);
    Nassoc_bkg_low_la    ->SetMarkerStyle(26);
    Nassoc_bkg_low_la    ->SetMarkerSize(1.5);
    Nassoc_bkg_low_la    ->SetLineColor(kBlue);
    Nassoc_obs_low_la    ->SetMarkerColor(kRed);
    Nassoc_obs_low_la    ->SetMarkerStyle(26);
    Nassoc_obs_low_la    ->SetMarkerSize(1.5);
    Nassoc_obs_low_la    ->SetLineColor(kRed);

    //DirectSubNass_xi     
    //Nassoc_bkg_xi        
    //Nassoc_obs_xi        
    //DirectSubNass_low_xi 
    //Nassoc_bkg_low_xi    
    //Nassoc_obs_low_xi    
    //
    /*
    TH1F *h1 = new TH1F("h1","",100,-5,5);
    TH1F* h2 = new TH1F("h2","",100,-5,5);
    h1->FillRandom("gaus");
    h2->FillRandom("gaus");
    TCanvas* c12 = new TCanvas("c12","c12",800,800);
    TPad *padc121 = new TPad("padc121","padc121",0,0.3,1,1.0);
    padc121->SetBottomMargin(0);
    padc121->SetGridx();
    padc121->Draw();
    padc121->cd();
    //DirectSubNass_ks->Draw("P");
    //Nassoc_obs_ks->Draw("P");
    //Nassoc_bkg_ks->Draw("P");
    h1->Draw();
    h2->Draw("same");

    h1->GetYaxis()->SetLabelSize(0);
    TGaxis *axis = new TGaxis(-5,20,-5,220,20,220,510,"");
    axis->SetLabelFont(43);
    axis->SetLabelSize(15);
    axis->Draw();

    c12->cd();

    TPad* padc122 = new TPad("padc122","padc122",0,0.05,1,0.3);
    padc122->SetTopMargin(0);
    padc122->SetBottomMargin(0.2);
    padc122->SetGridx();
    padc122->Draw();
    padc122->cd();

    TH1F *h3 = (TH1F*)h1->Clone("h3");
    h3->SetLineColor(kBlack);
    h3->SetMinimum(0.8);
    h3->SetMaximum(1.35);
    h3->Sumw2();
    h3->SetStats(0);
    h3->Divide(h2);
    h3->SetMarkerStyle(21);
    h3->Draw("ep");
    */

    double bottomMargin = 0.10;

    TCanvas* Nass_ks = new TCanvas("Nass_ks","Nass_ks",1300,1000);
    Nass_ks->Divide(2,2);

    TCanvas* Nass_la = new TCanvas("Nass_la","Nass_la",1300,1000);
    Nass_la->Divide(2,2);
    TH1F* frame_Nass_ks;
    TH1F* frame_Nass_ratio_ks;
    TH1F* frame_Nass_low_ks;
    TH1F* frame_Nass_ratio_low_ks;
    TH1F* frame_Nass_la;
    TH1F* frame_Nass_ratio_la;
    TH1F* frame_Nass_low_la;
    TH1F* frame_Nass_ratio_low_la;

    TLine* line_ratio = new TLine(0,1,9,1);

    frame_Nass_ks = Nass_ks->cd(1)->DrawFrame(0,6.8,9,8.2);
    Nass_ks->cd(1)->SetLeftMargin(bottomMargin);
    Nass_ks->cd(1)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Nass_ks->GetXaxis()->CenterTitle(1);
    frame_Nass_ks->GetYaxis()->CenterTitle(1);
    frame_Nass_ks->GetXaxis()->SetTitleSize(0.05);
    frame_Nass_ks->GetYaxis()->SetTitleSize(0.05);
    frame_Nass_ks->SetTitleOffset(1.1,"Y");
    frame_Nass_ks->SetTitleOffset(1.2,"X");
    frame_Nass_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Nass_ks->GetYaxis()->SetTitle("HM K_{S}^{0} N_{assoc}");
    TLegend* leg_Nass_ks = new TLegend(0.69,0.60,0.81,0.85);
    leg_Nass_ks->SetFillColor(10);
    leg_Nass_ks->SetFillStyle(0);
    leg_Nass_ks->SetBorderSize(0);
    leg_Nass_ks->SetTextFont(42);
    leg_Nass_ks->SetTextSize(0.04);
    leg_Nass_ks->AddEntry(DirectSubNass_ks, "Direct Sub", "P");
    leg_Nass_ks->AddEntry(Nassoc_obs_ks, "Sub Obs", "P");
    leg_Nass_ks->AddEntry(Nassoc_bkg_ks, "Sub Bkg", "P");
    leg_Nass_ks->Draw();
    DirectSubNass_ks->Draw("P");
    Nassoc_obs_ks->Draw("P");
    Nassoc_bkg_ks->Draw("P");

    frame_Nass_ratio_ks = Nass_ks->cd(2)->DrawFrame(0,0.95,9,1.05);
    Nass_ks->cd(2)->SetLeftMargin(bottomMargin);
    Nass_ks->cd(2)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Nass_ratio_ks->GetXaxis()->CenterTitle(1);
    frame_Nass_ratio_ks->GetYaxis()->CenterTitle(1);
    frame_Nass_ratio_ks->GetXaxis()->SetTitleSize(0.05);
    frame_Nass_ratio_ks->GetYaxis()->SetTitleSize(0.05);
    frame_Nass_ratio_ks->SetTitleOffset(1.1,"Y");
    frame_Nass_ratio_ks->SetTitleOffset(1.2,"X");
    frame_Nass_ratio_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Nass_ratio_ks->GetYaxis()->SetTitle("HM K_{S}^{0} N_{assoc}^{obs}/N_{assoc}^{bkg} ");
    TGNassoc_high_ratio_ks->Draw("P");
    line_ratio->Draw("same");

    frame_Nass_low_ks = Nass_ks->cd(3)->DrawFrame(0,0.52,9,0.65);
    Nass_ks->cd(3)->SetLeftMargin(bottomMargin);
    Nass_ks->cd(3)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Nass_low_ks->GetXaxis()->CenterTitle(1);
    frame_Nass_low_ks->GetYaxis()->CenterTitle(1);
    frame_Nass_low_ks->GetXaxis()->SetTitleSize(0.05);
    frame_Nass_low_ks->GetYaxis()->SetTitleSize(0.05);
    frame_Nass_low_ks->SetTitleOffset(1.1,"Y");
    frame_Nass_low_ks->SetTitleOffset(1.2,"X");
    frame_Nass_low_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Nass_low_ks->GetYaxis()->SetTitle("MB K_{S}^{0} N_{assoc}");
    TLegend* leg_Nass_low_ks = new TLegend(0.69,0.60,0.81,0.85);
    leg_Nass_low_ks->SetFillColor(10);
    leg_Nass_low_ks->SetFillStyle(0);
    leg_Nass_low_ks->SetBorderSize(0);
    leg_Nass_low_ks->SetTextFont(42);
    leg_Nass_low_ks->SetTextSize(0.04);
    leg_Nass_low_ks->AddEntry(DirectSubNass_low_ks, "Direct Sub", "P");
    leg_Nass_low_ks->AddEntry(Nassoc_obs_low_ks, "Sub Obs", "P");
    leg_Nass_low_ks->AddEntry(Nassoc_bkg_low_ks, "Sub Bkg", "P");
    leg_Nass_low_ks->Draw();
    DirectSubNass_low_ks->Draw("P");
    Nassoc_obs_low_ks->Draw("P");
    Nassoc_bkg_low_ks->Draw("P");

    frame_Nass_ratio_low_ks = Nass_ks->cd(4)->DrawFrame(0,0.90,9,1.15);
    Nass_ks->cd(4)->SetLeftMargin(bottomMargin);
    Nass_ks->cd(4)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Nass_ratio_low_ks->GetXaxis()->CenterTitle(1);
    frame_Nass_ratio_low_ks->GetYaxis()->CenterTitle(1);
    frame_Nass_ratio_low_ks->GetXaxis()->SetTitleSize(0.05);
    frame_Nass_ratio_low_ks->GetYaxis()->SetTitleSize(0.05);
    frame_Nass_ratio_low_ks->SetTitleOffset(1.1,"Y");
    frame_Nass_ratio_low_ks->SetTitleOffset(1.2,"X");
    frame_Nass_ratio_low_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Nass_ratio_low_ks->GetYaxis()->SetTitle("MB K_{S}^{0} N_{assoc}^{obs}/N_{assoc}^{bkg} ");
    TGNassoc_low_ratio_ks->Draw("P");
    line_ratio->Draw("same");


    frame_Nass_la = Nass_la->cd(1)->DrawFrame(0,6.8,9,8.0);
    Nass_la->cd(1)->SetLeftMargin(bottomMargin);
    Nass_la->cd(1)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Nass_la->GetXaxis()->CenterTitle(1);
    frame_Nass_la->GetYaxis()->CenterTitle(1);
    frame_Nass_la->GetXaxis()->SetTitleSize(0.05);
    frame_Nass_la->GetYaxis()->SetTitleSize(0.05);
    frame_Nass_la->SetTitleOffset(1.1,"Y");
    frame_Nass_la->SetTitleOffset(1.2,"X");
    frame_Nass_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Nass_la->GetYaxis()->SetTitle("HM #Lambda/#bar{#Lambda} N_{assoc}");
    TLegend* leg_Nass_la = new TLegend(0.69,0.60,0.81,0.85);
    leg_Nass_la->SetFillColor(10);
    leg_Nass_la->SetFillStyle(0);
    leg_Nass_la->SetBorderSize(0);
    leg_Nass_la->SetTextFont(42);
    leg_Nass_la->SetTextSize(0.04);
    leg_Nass_la->AddEntry(DirectSubNass_la, "Direct Sub", "P");
    leg_Nass_la->AddEntry(Nassoc_obs_la, "Sub Obs", "P");
    leg_Nass_la->AddEntry(Nassoc_bkg_la, "Sub Bkg", "P");
    leg_Nass_la->Draw();
    DirectSubNass_la->Draw("P");
    Nassoc_obs_la->Draw("P");
    Nassoc_bkg_la->Draw("P");

    frame_Nass_ratio_la = Nass_la->cd(2)->DrawFrame(0,0.95,9,1.05);
    Nass_la->cd(2)->SetLeftMargin(bottomMargin);
    Nass_la->cd(2)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Nass_ratio_la->GetXaxis()->CenterTitle(1);
    frame_Nass_ratio_la->GetYaxis()->CenterTitle(1);
    frame_Nass_ratio_la->GetXaxis()->SetTitleSize(0.05);
    frame_Nass_ratio_la->GetYaxis()->SetTitleSize(0.05);
    frame_Nass_ratio_la->SetTitleOffset(1.1,"Y");
    frame_Nass_ratio_la->SetTitleOffset(1.2,"X");
    frame_Nass_ratio_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Nass_ratio_la->GetYaxis()->SetTitle("HM #Lambda/#bar{#Lambda} N_{assoc}^{obs}/N_{assoc}^{bkg} ");
    TGNassoc_high_ratio_la->Draw("P");
    line_ratio->Draw("same");

    frame_Nass_low_la = Nass_la->cd(3)->DrawFrame(0,0.52,9,0.65);
    Nass_la->cd(3)->SetLeftMargin(bottomMargin);
    Nass_la->cd(3)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Nass_low_la->GetXaxis()->CenterTitle(1);
    frame_Nass_low_la->GetYaxis()->CenterTitle(1);
    frame_Nass_low_la->GetXaxis()->SetTitleSize(0.05);
    frame_Nass_low_la->GetYaxis()->SetTitleSize(0.05);
    frame_Nass_low_la->SetTitleOffset(1.1,"Y");
    frame_Nass_low_la->SetTitleOffset(1.2,"X");
    frame_Nass_low_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Nass_low_la->GetYaxis()->SetTitle("MB #Lambda/#bar{#Lambda} N_{assoc}");
    TLegend* leg_Nass_low_la = new TLegend(0.69,0.60,0.81,0.85);
    leg_Nass_low_la->SetFillColor(10);
    leg_Nass_low_la->SetFillStyle(0);
    leg_Nass_low_la->SetBorderSize(0);
    leg_Nass_low_la->SetTextFont(42);
    leg_Nass_low_la->SetTextSize(0.04);
    leg_Nass_low_la->AddEntry(DirectSubNass_low_la, "Direct Sub", "P");
    leg_Nass_low_la->AddEntry(Nassoc_obs_low_la, "Sub Obs", "P");
    leg_Nass_low_la->AddEntry(Nassoc_bkg_low_la, "Sub Bkg", "P");
    leg_Nass_low_la->Draw();
    DirectSubNass_low_la->Draw("P");
    Nassoc_obs_low_la->Draw("P");
    Nassoc_bkg_low_la->Draw("P");

    frame_Nass_ratio_low_la = Nass_la->cd(4)->DrawFrame(0,0.9,9,1.1);
    Nass_la->cd(4)->SetLeftMargin(bottomMargin);
    Nass_la->cd(4)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Nass_ratio_low_la->GetXaxis()->CenterTitle(1);
    frame_Nass_ratio_low_la->GetYaxis()->CenterTitle(1);
    frame_Nass_ratio_low_la->GetXaxis()->SetTitleSize(0.05);
    frame_Nass_ratio_low_la->GetYaxis()->SetTitleSize(0.05);
    frame_Nass_ratio_low_la->SetTitleOffset(1.1,"Y");
    frame_Nass_ratio_low_la->SetTitleOffset(1.2,"X");
    frame_Nass_ratio_low_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Nass_ratio_low_la->GetYaxis()->SetTitle("MB #Lambda/#bar{#Lambda} N_{assoc}^{obs}/N_{assoc}^{bkg} ");
    TGNassoc_low_ratio_la->Draw("P");
    line_ratio->Draw("same");

    //Yield
    TGraphErrors* DirectSubYield_ks     = (TGraphErrors*)file_pidv2->Get("DirectSubYield_ks");
    TGraphErrors* Yield_bkg_ks        = (TGraphErrors*)file_pidv2->Get("YieldPlot_bkg_ks");
    TGraphErrors* Yield_obs_ks        = (TGraphErrors*)file_pidv2->Get("YieldPlot_obs_ks");
    TGraphErrors* DirectSubYield_low_ks = (TGraphErrors*)file_pidv2->Get("DirectSubYield_low_ks");
    TGraphErrors* Yield_bkg_low_ks    = (TGraphErrors*)file_pidv2->Get("YieldPlot_low_bkg_ks");
    TGraphErrors* Yield_obs_low_ks    = (TGraphErrors*)file_pidv2->Get("YieldPlot_low_obs_ks");

    TGraphErrors* DirectSubYield_la     = (TGraphErrors*)file_pidv2->Get("DirectSubYield_la");
    TGraphErrors* Yield_bkg_la        = (TGraphErrors*)file_pidv2->Get("YieldPlot_bkg_la");
    TGraphErrors* Yield_obs_la        = (TGraphErrors*)file_pidv2->Get("YieldPlot_obs_la");
    TGraphErrors* DirectSubYield_low_la = (TGraphErrors*)file_pidv2->Get("DirectSubYield_low_la");
    TGraphErrors* Yield_bkg_low_la    = (TGraphErrors*)file_pidv2->Get("YieldPlot_low_bkg_la");
    TGraphErrors* Yield_obs_low_la    = (TGraphErrors*)file_pidv2->Get("YieldPlot_low_obs_la");

    TGraphErrors* DirectSubYield_xi     = (TGraphErrors*)file_pidv2->Get("DirectSubYield_xi");
    TGraphErrors* Yield_bkg_xi        = (TGraphErrors*)file_pidv2->Get("YieldPlot_bkg_xi");
    TGraphErrors* Yield_obs_xi        = (TGraphErrors*)file_pidv2->Get("YieldPlot_obs_xi");
    TGraphErrors* DirectSubYield_low_xi = (TGraphErrors*)file_pidv2->Get("DirectSubYield_low_xi");
    TGraphErrors* Yield_bkg_low_xi    = (TGraphErrors*)file_pidv2->Get("YieldPlot_low_bkg_xi");
    TGraphErrors* Yield_obs_low_xi    = (TGraphErrors*)file_pidv2->Get("YieldPlot_low_obs_xi");

    double* Yield_bkg_ks_Y = Yield_bkg_ks->GetY();
    double* Yield_obs_ks_Y = Yield_obs_ks->GetY();
    double* Yield_bkg_low_ks_Y = Yield_bkg_low_ks->GetY();
    double* Yield_obs_low_ks_Y = Yield_obs_low_ks->GetY();

    double* Yield_bkg_la_Y = Yield_bkg_la->GetY();
    double* Yield_obs_la_Y = Yield_obs_la->GetY();
    double* Yield_bkg_low_la_Y = Yield_bkg_low_la->GetY();
    double* Yield_obs_low_la_Y = Yield_obs_low_la->GetY();

    std::vector<double> Yield_high_ratio_ks;
    std::vector<double> Yield_high_err_ratio_ks;
    std::vector<double> Yield_low_ratio_ks;
    std::vector<double> Yield_low_err_ratio_ks;

    std::vector<double> Yield_high_ratio_la;
    std::vector<double> Yield_high_err_ratio_la;
    std::vector<double> Yield_low_ratio_la;
    std::vector<double> Yield_low_err_ratio_la;

    //Do obs/bkg
    for(int i=0; i<Yield_bkg_ks->GetN(); i++)
    {
        Yield_high_ratio_ks.push_back(Yield_obs_ks_Y[i]/Yield_bkg_ks_Y[i]);
        Yield_low_ratio_ks.push_back(Yield_obs_low_ks_Y[i]/Yield_bkg_low_ks_Y[i]);
    }

    for(int i=0; i<Yield_bkg_la->GetN(); i++)
    {
        Yield_high_ratio_la.push_back(Yield_obs_la_Y[i]/Yield_bkg_la_Y[i]);
        Yield_low_ratio_la.push_back(Yield_obs_low_la_Y[i]/Yield_bkg_low_la_Y[i]);
    }

    TGraphErrors* TGYield_high_ratio_ks = new TGraphErrors(numKs,pt_ks,&Yield_high_ratio_ks[0],0,0);
    TGraphErrors* TGYield_low_ratio_ks = new TGraphErrors(numKs,pt_ks,&Yield_low_ratio_ks[0],0,0);
    TGraphErrors* TGYield_high_ratio_la = new TGraphErrors(numLa,pt_la,&Yield_high_ratio_la[0],0,0);
    TGraphErrors* TGYield_low_ratio_la = new TGraphErrors(numLa,pt_la,&Yield_low_ratio_la[0],0,0);

    TGYield_high_ratio_ks->SetMarkerColor(kRed);
    TGYield_high_ratio_ks->SetMarkerStyle(20);
    TGYield_high_ratio_ks->SetMarkerSize(1.2);
    TGYield_high_ratio_ks->SetLineColor(kRed);

    TGYield_low_ratio_ks->SetMarkerColor(kBlue);
    TGYield_low_ratio_ks->SetMarkerStyle(20);
    TGYield_low_ratio_ks->SetMarkerSize(1.2);
    TGYield_low_ratio_ks->SetLineColor(kBlue);

    TGYield_high_ratio_la->SetMarkerColor(kRed);
    TGYield_high_ratio_la->SetMarkerStyle(22);
    TGYield_high_ratio_la->SetMarkerSize(1.2);
    TGYield_high_ratio_la->SetLineColor(kRed);

    TGYield_low_ratio_la->SetMarkerColor(kBlue);
    TGYield_low_ratio_la->SetMarkerStyle(22);
    TGYield_low_ratio_la->SetMarkerSize(1.2);
    TGYield_low_ratio_la->SetLineColor(kBlue);

    DirectSubYield_ks     ->SetMarkerColor(kGreen-2);
    DirectSubYield_ks     ->SetMarkerStyle(24);
    DirectSubYield_ks     ->SetMarkerSize(1.5);
    DirectSubYield_ks     ->SetLineColor(kGreen-2);
    Yield_bkg_ks        ->SetMarkerColor(kBlue);
    Yield_bkg_ks        ->SetMarkerStyle(24);
    Yield_bkg_ks        ->SetMarkerSize(1.5);
    Yield_bkg_ks        ->SetLineColor(kBlue);
    Yield_obs_ks        ->SetMarkerColor(kRed);
    Yield_obs_ks        ->SetMarkerStyle(24);
    Yield_obs_ks        ->SetMarkerSize(1.5);
    Yield_obs_ks        ->SetLineColor(kRed);

    DirectSubYield_low_ks ->SetMarkerColor(kGreen-2);
    DirectSubYield_low_ks ->SetMarkerStyle(24);
    DirectSubYield_low_ks ->SetMarkerSize(1.5);
    DirectSubYield_low_ks ->SetLineColor(kGreen-2);
    Yield_bkg_low_ks    ->SetMarkerColor(kBlue);
    Yield_bkg_low_ks    ->SetMarkerStyle(24);
    Yield_bkg_low_ks    ->SetMarkerSize(1.5);
    Yield_bkg_low_ks    ->SetLineColor(kBlue);
    Yield_obs_low_ks    ->SetMarkerColor(kRed);
    Yield_obs_low_ks    ->SetMarkerStyle(24);
    Yield_obs_low_ks    ->SetMarkerSize(1.5);
    Yield_obs_low_ks    ->SetLineColor(kRed);

    DirectSubYield_la     ->SetMarkerColor(kGreen-2);
    DirectSubYield_la     ->SetMarkerStyle(26);
    DirectSubYield_la     ->SetMarkerSize(1.5);
    DirectSubYield_la     ->SetLineColor(kGreen-2);
    Yield_bkg_la        ->SetMarkerColor(kBlue);
    Yield_bkg_la        ->SetMarkerStyle(26);
    Yield_bkg_la        ->SetMarkerSize(1.5);
    Yield_bkg_la        ->SetLineColor(kBlue);
    Yield_obs_la        ->SetMarkerColor(kRed);
    Yield_obs_la        ->SetMarkerStyle(26);
    Yield_obs_la        ->SetMarkerSize(1.5);
    Yield_obs_la        ->SetLineColor(kRed);

    DirectSubYield_low_la ->SetMarkerColor(kGreen-2);
    DirectSubYield_low_la ->SetMarkerStyle(26);
    DirectSubYield_low_la ->SetMarkerSize(1.5);
    DirectSubYield_low_la ->SetLineColor(kGreen-2);
    Yield_bkg_low_la    ->SetMarkerColor(kBlue);
    Yield_bkg_low_la    ->SetMarkerStyle(26);
    Yield_bkg_low_la    ->SetMarkerSize(1.5);
    Yield_bkg_low_la    ->SetLineColor(kBlue);
    Yield_obs_low_la    ->SetMarkerColor(kRed);
    Yield_obs_low_la    ->SetMarkerStyle(26);
    Yield_obs_low_la    ->SetMarkerSize(1.5);
    Yield_obs_low_la    ->SetLineColor(kRed);

    TCanvas* Yield_ks = new TCanvas("Yield_ks","Yield_ks",1200,900);
    Yield_ks->Divide(2,2);

    TCanvas* Yield_la = new TCanvas("Yield_la","Yield_la",1200,900);
    Yield_la->Divide(2,2);
    TH1F* frame_Yield_ks;
    TH1F* frame_Yield_ratio_ks;
    TH1F* frame_Yield_low_ks;
    TH1F* frame_Yield_ratio_low_ks;
    TH1F* frame_Yield_la;
    TH1F* frame_Yield_ratio_la;
    TH1F* frame_Yield_low_la;
    TH1F* frame_Yield_ratio_low_la;

    frame_Yield_ks = Yield_ks->cd(1)->DrawFrame(0,0,9,1.3);
    Yield_ks->cd(1)->SetLeftMargin(bottomMargin);
    Yield_ks->cd(1)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Yield_ks->GetXaxis()->CenterTitle(1);
    frame_Yield_ks->GetYaxis()->CenterTitle(1);
    frame_Yield_ks->GetXaxis()->SetTitleSize(0.05);
    frame_Yield_ks->GetYaxis()->SetTitleSize(0.05);
    frame_Yield_ks->SetTitleOffset(1.1,"Y");
    frame_Yield_ks->SetTitleOffset(1.2,"X");
    frame_Yield_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Yield_ks->GetYaxis()->SetTitle("HM K_{S}^{0} Y");
    TLegend* leg_Yield_ks = new TLegend(0.17,0.63,0.29,0.88);
    leg_Yield_ks->SetFillColor(10);
    leg_Yield_ks->SetFillStyle(0);
    leg_Yield_ks->SetBorderSize(0);
    leg_Yield_ks->SetTextFont(42);
    leg_Yield_ks->SetTextSize(0.04);
    leg_Yield_ks->AddEntry(DirectSubYield_ks, "Direct Sub", "P");
    leg_Yield_ks->AddEntry(Yield_obs_ks, "Sub Obs", "P");
    leg_Yield_ks->AddEntry(Yield_bkg_ks, "Sub Bkg", "P");
    leg_Yield_ks->Draw();
    DirectSubYield_ks->Draw("P");
    Yield_obs_ks->Draw("P");
    Yield_bkg_ks->Draw("P");

    frame_Yield_ratio_ks = Yield_ks->cd(2)->DrawFrame(0,0.0,9,1.5);
    Yield_ks->cd(2)->SetLeftMargin(bottomMargin);
    Yield_ks->cd(2)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Yield_ratio_ks->GetXaxis()->CenterTitle(1);
    frame_Yield_ratio_ks->GetYaxis()->CenterTitle(1);
    frame_Yield_ratio_ks->GetXaxis()->SetTitleSize(0.05);
    frame_Yield_ratio_ks->GetYaxis()->SetTitleSize(0.05);
    frame_Yield_ratio_ks->SetTitleOffset(1.1,"Y");
    frame_Yield_ratio_ks->SetTitleOffset(1.2,"X");
    frame_Yield_ratio_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Yield_ratio_ks->GetYaxis()->SetTitle("HM K_{S}^{0} Y^{obs}/Y^{bkg} ");
    TGYield_high_ratio_ks->Draw("P");
    line_ratio->Draw("same");

    frame_Yield_low_ks = Yield_ks->cd(3)->DrawFrame(0,0.0,9,0.8);
    Yield_ks->cd(3)->SetLeftMargin(bottomMargin);
    Yield_ks->cd(3)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Yield_low_ks->GetXaxis()->CenterTitle(1);
    frame_Yield_low_ks->GetYaxis()->CenterTitle(1);
    frame_Yield_low_ks->GetXaxis()->SetTitleSize(0.05);
    frame_Yield_low_ks->GetYaxis()->SetTitleSize(0.05);
    frame_Yield_low_ks->SetTitleOffset(1.1,"Y");
    frame_Yield_low_ks->SetTitleOffset(1.2,"X");
    frame_Yield_low_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Yield_low_ks->GetYaxis()->SetTitle("MB K_{S}^{0} Y");
    TLegend* leg_Yield_low_ks = new TLegend(0.17,0.63,0.29,0.88);
    leg_Yield_low_ks->SetFillColor(10);
    leg_Yield_low_ks->SetFillStyle(0);
    leg_Yield_low_ks->SetBorderSize(0);
    leg_Yield_low_ks->SetTextFont(42);
    leg_Yield_low_ks->SetTextSize(0.04);
    leg_Yield_low_ks->AddEntry(DirectSubYield_low_ks, "Direct Sub", "P");
    leg_Yield_low_ks->AddEntry(Yield_obs_low_ks, "Sub Obs", "P");
    leg_Yield_low_ks->AddEntry(Yield_bkg_low_ks, "Sub Bkg", "P");
    leg_Yield_low_ks->Draw();
    DirectSubYield_low_ks->Draw("P");
    Yield_obs_low_ks->Draw("P");
    Yield_bkg_low_ks->Draw("P");

    frame_Yield_ratio_low_ks = Yield_ks->cd(4)->DrawFrame(0,0.00,9,1.5);
    Yield_ks->cd(4)->SetLeftMargin(bottomMargin);
    Yield_ks->cd(4)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Yield_ratio_low_ks->GetXaxis()->CenterTitle(1);
    frame_Yield_ratio_low_ks->GetYaxis()->CenterTitle(1);
    frame_Yield_ratio_low_ks->GetXaxis()->SetTitleSize(0.05);
    frame_Yield_ratio_low_ks->GetYaxis()->SetTitleSize(0.05);
    frame_Yield_ratio_low_ks->SetTitleOffset(1.1,"Y");
    frame_Yield_ratio_low_ks->SetTitleOffset(1.2,"X");
    frame_Yield_ratio_low_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Yield_ratio_low_ks->GetYaxis()->SetTitle("MB K_{S}^{0} Y^{obs}/Y^{bkg} ");
    TGYield_low_ratio_ks->Draw("P");
    line_ratio->Draw("same");


    frame_Yield_la = Yield_la->cd(1)->DrawFrame(0,0.0,9,1.3);
    Yield_la->cd(1)->SetLeftMargin(bottomMargin);
    Yield_la->cd(1)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Yield_la->GetXaxis()->CenterTitle(1);
    frame_Yield_la->GetYaxis()->CenterTitle(1);
    frame_Yield_la->GetXaxis()->SetTitleSize(0.05);
    frame_Yield_la->GetYaxis()->SetTitleSize(0.05);
    frame_Yield_la->SetTitleOffset(1.1,"Y");
    frame_Yield_la->SetTitleOffset(1.2,"X");
    frame_Yield_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Yield_la->GetYaxis()->SetTitle("HM #Lambda/#bar{#Lambda} Y");
    TLegend* leg_Yield_la = new TLegend(0.17,0.63,0.29,0.88);
    leg_Yield_la->SetFillColor(10);
    leg_Yield_la->SetFillStyle(0);
    leg_Yield_la->SetBorderSize(0);
    leg_Yield_la->SetTextFont(42);
    leg_Yield_la->SetTextSize(0.04);
    leg_Yield_la->AddEntry(DirectSubYield_la, "Direct Sub", "P");
    leg_Yield_la->AddEntry(Yield_obs_la, "Sub Obs", "P");
    leg_Yield_la->AddEntry(Yield_bkg_la, "Sub Bkg", "P");
    leg_Yield_la->Draw();
    DirectSubYield_la->Draw("P");
    Yield_obs_la->Draw("P");
    Yield_bkg_la->Draw("P");

    frame_Yield_ratio_la = Yield_la->cd(2)->DrawFrame(0,0.00,9,1.5);
    Yield_la->cd(2)->SetLeftMargin(bottomMargin);
    Yield_la->cd(2)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Yield_ratio_la->GetXaxis()->CenterTitle(1);
    frame_Yield_ratio_la->GetYaxis()->CenterTitle(1);
    frame_Yield_ratio_la->GetXaxis()->SetTitleSize(0.05);
    frame_Yield_ratio_la->GetYaxis()->SetTitleSize(0.05);
    frame_Yield_ratio_la->SetTitleOffset(1.1,"Y");
    frame_Yield_ratio_la->SetTitleOffset(1.2,"X");
    frame_Yield_ratio_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Yield_ratio_la->GetYaxis()->SetTitle("HM #Lambda/#bar{#Lambda} Y^{obs}/Y^{bkg} ");
    TGYield_high_ratio_la->Draw("P");
    line_ratio->Draw("same");

    frame_Yield_low_la = Yield_la->cd(3)->DrawFrame(0,0.00,9,1.5);
    Yield_la->cd(3)->SetLeftMargin(bottomMargin);
    Yield_la->cd(3)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Yield_low_la->GetXaxis()->CenterTitle(1);
    frame_Yield_low_la->GetYaxis()->CenterTitle(1);
    frame_Yield_low_la->GetXaxis()->SetTitleSize(0.05);
    frame_Yield_low_la->GetYaxis()->SetTitleSize(0.05);
    frame_Yield_low_la->SetTitleOffset(1.1,"Y");
    frame_Yield_low_la->SetTitleOffset(1.2,"X");
    frame_Yield_low_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Yield_low_la->GetYaxis()->SetTitle("MB #Lambda/#bar{#Lambda} Y");
    TLegend* leg_Yield_low_la = new TLegend(0.17,0.63,0.29,0.88);
    leg_Yield_low_la->SetFillColor(10);
    leg_Yield_low_la->SetFillStyle(0);
    leg_Yield_low_la->SetBorderSize(0);
    leg_Yield_low_la->SetTextFont(42);
    leg_Yield_low_la->SetTextSize(0.04);
    leg_Yield_low_la->AddEntry(DirectSubYield_low_la, "Direct Sub", "P");
    leg_Yield_low_la->AddEntry(Yield_obs_low_la, "Sub Obs", "P");
    leg_Yield_low_la->AddEntry(Yield_bkg_low_la, "Sub Bkg", "P");
    leg_Yield_low_la->Draw();
    DirectSubYield_low_la->Draw("P");
    Yield_obs_low_la->Draw("P");
    Yield_bkg_low_la->Draw("P");

    frame_Yield_ratio_low_la = Yield_la->cd(4)->DrawFrame(0,0.0,9,1.5);
    Yield_la->cd(4)->SetLeftMargin(bottomMargin);
    Yield_la->cd(4)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Yield_ratio_low_la->GetXaxis()->CenterTitle(1);
    frame_Yield_ratio_low_la->GetYaxis()->CenterTitle(1);
    frame_Yield_ratio_low_la->GetXaxis()->SetTitleSize(0.05);
    frame_Yield_ratio_low_la->GetYaxis()->SetTitleSize(0.05);
    frame_Yield_ratio_low_la->SetTitleOffset(1.1,"Y");
    frame_Yield_ratio_low_la->SetTitleOffset(1.2,"X");
    frame_Yield_ratio_low_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Yield_ratio_low_la->GetYaxis()->SetTitle("MB #Lambda/#bar{#Lambda} Y^{obs}/Y^{bkg} ");
    TGYield_low_ratio_la->Draw("P");
    line_ratio->Draw("same");

    //V2 obs and bkg
    TGraphErrors* V2_bkg_ks        = (TGraphErrors*)file_pidv2->Get("V2plot_bkg_ks");
    TGraphErrors* V2_obs_ks        = (TGraphErrors*)file_pidv2->Get("V2plot_obs_ks");
    TGraphErrors* V2_bkg_low_ks    = (TGraphErrors*)file_pidv2->Get("V2plot_bkg_low_ks");
    TGraphErrors* V2_obs_low_ks    = (TGraphErrors*)file_pidv2->Get("V2plot_obs_low_ks");

    TGraphErrors* V2_bkg_la        = (TGraphErrors*)file_pidv2->Get("V2plot_bkg_la");
    TGraphErrors* V2_obs_la        = (TGraphErrors*)file_pidv2->Get("V2plot_obs_la");
    TGraphErrors* V2_bkg_low_la    = (TGraphErrors*)file_pidv2->Get("V2plot_bkg_low_la");
    TGraphErrors* V2_obs_low_la    = (TGraphErrors*)file_pidv2->Get("V2plot_obs_low_la");

    TGraphErrors* V2_bkg_xi        = (TGraphErrors*)file_pidv2->Get("V2plot_bkg_xi");
    TGraphErrors* V2_obs_xi        = (TGraphErrors*)file_pidv2->Get("V2plot_obs_xi");
    TGraphErrors* V2_bkg_low_xi    = (TGraphErrors*)file_pidv2->Get("V2plot_bkg_low_xi");
    TGraphErrors* V2_obs_low_xi    = (TGraphErrors*)file_pidv2->Get("V2plot_obs_low_xi");

    double* V2_bkg_ks_Y = V2_bkg_ks->GetY();
    double* V2_bkg_ks_EY = V2_bkg_ks->GetEY();
    double* V2_obs_ks_Y = V2_obs_ks->GetY();
    double* V2_obs_ks_EY = V2_obs_ks->GetEY();
    double* V2_bkg_low_ks_Y = V2_bkg_low_ks->GetY();
    double* V2_bkg_low_ks_EY = V2_bkg_low_ks->GetEY();
    double* V2_obs_low_ks_Y = V2_obs_low_ks->GetY();
    double* V2_obs_low_ks_EY = V2_obs_low_ks->GetEY();

    double* V2_bkg_la_Y = V2_bkg_la->GetY();
    double* V2_bkg_la_EY = V2_bkg_la->GetEY();
    double* V2_obs_la_Y = V2_obs_la->GetY();
    double* V2_obs_la_EY = V2_obs_la->GetEY();
    double* V2_bkg_low_la_Y = V2_bkg_low_la->GetY();
    double* V2_bkg_low_la_EY = V2_bkg_low_la->GetEY();
    double* V2_obs_low_la_Y = V2_obs_low_la->GetY();
    double* V2_obs_low_la_EY = V2_obs_low_la->GetEY();

    std::vector<double> V2_high_ratio_ks;
    std::vector<double> V2_high_err_ratio_ks;
    std::vector<double> V2_low_ratio_ks;
    std::vector<double> V2_low_err_ratio_ks;

    std::vector<double> V2_high_ratio_la;
    std::vector<double> V2_high_err_ratio_la;
    std::vector<double> V2_low_ratio_la;
    std::vector<double> V2_low_err_ratio_la;

    //Do obs/bkg
    for(int i=0; i<V2_bkg_ks->GetN(); i++)
    {
        V2_high_ratio_ks.push_back(V2_obs_ks_Y[i]/V2_bkg_ks_Y[i]);
        V2_high_err_ratio_ks.push_back(V2_obs_ks_EY[i]/V2_bkg_ks_Y[i]);
        V2_low_ratio_ks.push_back(V2_obs_low_ks_Y[i]/V2_bkg_low_ks_Y[i]);
        V2_low_err_ratio_ks.push_back(V2_obs_low_ks_EY[i]/V2_bkg_low_ks_Y[i]);
    }

    for(int i=0; i<V2_bkg_la->GetN(); i++)
    {
        V2_high_ratio_la.push_back(V2_obs_la_Y[i]/V2_bkg_la_Y[i]);
        V2_high_err_ratio_la.push_back(V2_obs_la_EY[i]/V2_bkg_la_Y[i]);
        V2_low_ratio_la.push_back(V2_obs_low_la_Y[i]/V2_bkg_low_la_Y[i]);
        V2_low_err_ratio_la.push_back(V2_obs_low_la_EY[i]/V2_bkg_low_la_Y[i]);
    }

    TGraphErrors* TGV2_high_ratio_ks = new TGraphErrors(numKs,pt_ks,&V2_high_ratio_ks[0],0,&V2_high_err_ratio_ks[0]);
    TGraphErrors* TGV2_low_ratio_ks = new TGraphErrors(numKs,pt_ks,&V2_low_ratio_ks[0],0,&V2_low_err_ratio_ks[0]);
    TGraphErrors* TGV2_high_ratio_la = new TGraphErrors(numLa,pt_la,&V2_high_ratio_la[0],0,&V2_high_err_ratio_la[0]);
    TGraphErrors* TGV2_low_ratio_la = new TGraphErrors(numLa,pt_la,&V2_low_ratio_la[0],0,&V2_low_err_ratio_la[0]);

    TGV2_high_ratio_ks->SetMarkerColor(kRed);
    TGV2_high_ratio_ks->SetMarkerStyle(20);
    TGV2_high_ratio_ks->SetMarkerSize(1.2);
    TGV2_high_ratio_ks->SetLineColor(kRed);

    TGV2_low_ratio_ks->SetMarkerColor(kBlue);
    TGV2_low_ratio_ks->SetMarkerStyle(20);
    TGV2_low_ratio_ks->SetMarkerSize(1.2);
    TGV2_low_ratio_ks->SetLineColor(kBlue);

    TGV2_high_ratio_la->SetMarkerColor(kRed);
    TGV2_high_ratio_la->SetMarkerStyle(22);
    TGV2_high_ratio_la->SetMarkerSize(1.2);
    TGV2_high_ratio_la->SetLineColor(kRed);

    TGV2_low_ratio_la->SetMarkerColor(kBlue);
    TGV2_low_ratio_la->SetMarkerStyle(22);
    TGV2_low_ratio_la->SetMarkerSize(1.2);
    TGV2_low_ratio_la->SetLineColor(kBlue);

    V2_bkg_ks        ->SetMarkerColor(kBlue);
    V2_bkg_ks        ->SetMarkerStyle(24);
    V2_bkg_ks        ->SetMarkerSize(1.5);
    V2_bkg_ks        ->SetLineColor(kBlue);
    V2_obs_ks        ->SetMarkerColor(kRed);
    V2_obs_ks        ->SetMarkerStyle(24);
    V2_obs_ks        ->SetMarkerSize(1.5);
    V2_obs_ks        ->SetLineColor(kRed);

    V2_bkg_low_ks    ->SetMarkerColor(kBlue);
    V2_bkg_low_ks    ->SetMarkerStyle(24);
    V2_bkg_low_ks    ->SetMarkerSize(1.5);
    V2_bkg_low_ks    ->SetLineColor(kBlue);
    V2_obs_low_ks    ->SetMarkerColor(kRed);
    V2_obs_low_ks    ->SetMarkerStyle(24);
    V2_obs_low_ks    ->SetMarkerSize(1.5);
    V2_obs_low_ks    ->SetLineColor(kRed);

    V2_bkg_la        ->SetMarkerColor(kBlue);
    V2_bkg_la        ->SetMarkerStyle(26);
    V2_bkg_la        ->SetMarkerSize(1.5);
    V2_bkg_la        ->SetLineColor(kBlue);
    V2_obs_la        ->SetMarkerColor(kRed);
    V2_obs_la        ->SetMarkerStyle(26);
    V2_obs_la        ->SetMarkerSize(1.5);
    V2_obs_la        ->SetLineColor(kRed);

    V2_bkg_low_la    ->SetMarkerColor(kBlue);
    V2_bkg_low_la    ->SetMarkerStyle(26);
    V2_bkg_low_la    ->SetMarkerSize(1.5);
    V2_bkg_low_la    ->SetLineColor(kBlue);
    V2_obs_low_la    ->SetMarkerColor(kRed);
    V2_obs_low_la    ->SetMarkerStyle(26);
    V2_obs_low_la    ->SetMarkerSize(1.5);
    V2_obs_low_la    ->SetLineColor(kRed);

    TCanvas* V2_ks = new TCanvas("V2_ks","V2_ks",1200,900);
    V2_ks->Divide(2,2);

    TCanvas* V2_la = new TCanvas("V2_la","V2_la",1200,900);
    V2_la->Divide(2,2);
    TH1F* frame_V2_ks;
    TH1F* frame_V2_ratio_ks;
    TH1F* frame_V2_low_ks;
    TH1F* frame_V2_ratio_low_ks;
    TH1F* frame_V2_la;
    TH1F* frame_V2_ratio_la;
    TH1F* frame_V2_low_la;
    TH1F* frame_V2_ratio_low_la;

    frame_V2_ks = V2_ks->cd(1)->DrawFrame(0,-0.001,9,0.04);
    V2_ks->cd(1)->SetLeftMargin(bottomMargin);
    V2_ks->cd(1)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_V2_ks->GetXaxis()->CenterTitle(1);
    frame_V2_ks->GetYaxis()->CenterTitle(1);
    frame_V2_ks->GetXaxis()->SetTitleSize(0.05);
    frame_V2_ks->GetYaxis()->SetTitleSize(0.05);
    frame_V2_ks->SetTitleOffset(1.1,"Y");
    frame_V2_ks->SetTitleOffset(1.2,"X");
    frame_V2_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_V2_ks->GetYaxis()->SetTitle("HM K_{S}^{0} V_{2}");
    TLegend* leg_V2_ks = new TLegend(0.17,0.63,0.29,0.88);
    leg_V2_ks->SetFillColor(10);
    leg_V2_ks->SetFillStyle(0);
    leg_V2_ks->SetBorderSize(0);
    leg_V2_ks->SetTextFont(42);
    leg_V2_ks->SetTextSize(0.04);
    leg_V2_ks->AddEntry(V2_obs_ks, "Obs", "P");
    leg_V2_ks->AddEntry(V2_bkg_ks, "Bkg", "P");
    leg_V2_ks->Draw();
    V2_obs_ks->Draw("P");
    V2_bkg_ks->Draw("P");

    frame_V2_ratio_ks = V2_ks->cd(2)->DrawFrame(0,0.0,9,1.5);
    V2_ks->cd(2)->SetLeftMargin(bottomMargin);
    V2_ks->cd(2)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_V2_ratio_ks->GetXaxis()->CenterTitle(1);
    frame_V2_ratio_ks->GetYaxis()->CenterTitle(1);
    frame_V2_ratio_ks->GetXaxis()->SetTitleSize(0.05);
    frame_V2_ratio_ks->GetYaxis()->SetTitleSize(0.05);
    frame_V2_ratio_ks->SetTitleOffset(1.1,"Y");
    frame_V2_ratio_ks->SetTitleOffset(1.2,"X");
    frame_V2_ratio_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_V2_ratio_ks->GetYaxis()->SetTitle("HM K_{S}^{0} V_{2}^{obs}/V_{2}^{bkg} ");
    TGV2_high_ratio_ks->Draw("P");
    line_ratio->Draw("same");

    frame_V2_low_ks = V2_ks->cd(3)->DrawFrame(0,-0.01,9,0.1);
    V2_ks->cd(3)->SetLeftMargin(bottomMargin);
    V2_ks->cd(3)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_V2_low_ks->GetXaxis()->CenterTitle(1);
    frame_V2_low_ks->GetYaxis()->CenterTitle(1);
    frame_V2_low_ks->GetXaxis()->SetTitleSize(0.05);
    frame_V2_low_ks->GetYaxis()->SetTitleSize(0.05);
    frame_V2_low_ks->SetTitleOffset(1.1,"Y");
    frame_V2_low_ks->SetTitleOffset(1.2,"X");
    frame_V2_low_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_V2_low_ks->GetYaxis()->SetTitle("MB K_{S}^{0} V_{2}");
    TLegend* leg_V2_low_ks = new TLegend(0.17,0.63,0.29,0.88);
    leg_V2_low_ks->SetFillColor(10);
    leg_V2_low_ks->SetFillStyle(0);
    leg_V2_low_ks->SetBorderSize(0);
    leg_V2_low_ks->SetTextFont(42);
    leg_V2_low_ks->SetTextSize(0.04);
    leg_V2_low_ks->AddEntry(V2_obs_low_ks, "Obs", "P");
    leg_V2_low_ks->AddEntry(V2_bkg_low_ks, "Bkg", "P");
    leg_V2_low_ks->Draw();
    V2_obs_low_ks->Draw("P");
    V2_bkg_low_ks->Draw("P");

    frame_V2_ratio_low_ks = V2_ks->cd(4)->DrawFrame(0,0.00,9,1.5);
    V2_ks->cd(4)->SetLeftMargin(bottomMargin);
    V2_ks->cd(4)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_V2_ratio_low_ks->GetXaxis()->CenterTitle(1);
    frame_V2_ratio_low_ks->GetYaxis()->CenterTitle(1);
    frame_V2_ratio_low_ks->GetXaxis()->SetTitleSize(0.05);
    frame_V2_ratio_low_ks->GetYaxis()->SetTitleSize(0.05);
    frame_V2_ratio_low_ks->SetTitleOffset(1.1,"Y");
    frame_V2_ratio_low_ks->SetTitleOffset(1.2,"X");
    frame_V2_ratio_low_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_V2_ratio_low_ks->GetYaxis()->SetTitle("MB K_{S}^{0} V_{2}^{obs}/V_{2}^{bkg} ");
    TGV2_low_ratio_ks->Draw("P");
    line_ratio->Draw("same");


    frame_V2_la = V2_la->cd(1)->DrawFrame(0,-0.001,9,0.04);
    V2_la->cd(1)->SetLeftMargin(bottomMargin);
    V2_la->cd(1)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_V2_la->GetXaxis()->CenterTitle(1);
    frame_V2_la->GetYaxis()->CenterTitle(1);
    frame_V2_la->GetXaxis()->SetTitleSize(0.05);
    frame_V2_la->GetYaxis()->SetTitleSize(0.05);
    frame_V2_la->SetTitleOffset(1.1,"Y");
    frame_V2_la->SetTitleOffset(1.2,"X");
    frame_V2_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_V2_la->GetYaxis()->SetTitle("HM #Lambda/#bar{#Lambda} V_{2}");
    TLegend* leg_V2_la = new TLegend(0.17,0.63,0.29,0.88);
    leg_V2_la->SetFillColor(10);
    leg_V2_la->SetFillStyle(0);
    leg_V2_la->SetBorderSize(0);
    leg_V2_la->SetTextFont(42);
    leg_V2_la->SetTextSize(0.04);
    leg_V2_la->AddEntry(V2_obs_la, "Obs", "P");
    leg_V2_la->AddEntry(V2_bkg_la, "Bkg", "P");
    leg_V2_la->Draw();
    V2_obs_la->Draw("P");
    V2_bkg_la->Draw("P");

    frame_V2_ratio_la = V2_la->cd(2)->DrawFrame(0,0.00,9,1.5);
    V2_la->cd(2)->SetLeftMargin(bottomMargin);
    V2_la->cd(2)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_V2_ratio_la->GetXaxis()->CenterTitle(1);
    frame_V2_ratio_la->GetYaxis()->CenterTitle(1);
    frame_V2_ratio_la->GetXaxis()->SetTitleSize(0.05);
    frame_V2_ratio_la->GetYaxis()->SetTitleSize(0.05);
    frame_V2_ratio_la->SetTitleOffset(1.1,"Y");
    frame_V2_ratio_la->SetTitleOffset(1.2,"X");
    frame_V2_ratio_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_V2_ratio_la->GetYaxis()->SetTitle("HM #Lambda/#bar{#Lambda} V_{2}^{obs}/V_{2}^{bkg} ");
    TGV2_high_ratio_la->Draw("P");
    line_ratio->Draw("same");

    frame_V2_low_la = V2_la->cd(3)->DrawFrame(0,-0.01,9,0.1);
    V2_la->cd(3)->SetLeftMargin(bottomMargin);
    V2_la->cd(3)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_V2_low_la->GetXaxis()->CenterTitle(1);
    frame_V2_low_la->GetYaxis()->CenterTitle(1);
    frame_V2_low_la->GetXaxis()->SetTitleSize(0.05);
    frame_V2_low_la->GetYaxis()->SetTitleSize(0.05);
    frame_V2_low_la->SetTitleOffset(1.1,"Y");
    frame_V2_low_la->SetTitleOffset(1.2,"X");
    frame_V2_low_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_V2_low_la->GetYaxis()->SetTitle("MB #Lambda/#bar{#Lambda} V_{2}");
    TLegend* leg_V2_low_la = new TLegend(0.17,0.63,0.29,0.88);
    leg_V2_low_la->SetFillColor(10);
    leg_V2_low_la->SetFillStyle(0);
    leg_V2_low_la->SetBorderSize(0);
    leg_V2_low_la->SetTextFont(42);
    leg_V2_low_la->SetTextSize(0.04);
    leg_V2_low_la->AddEntry(V2_obs_low_la, "Obs", "P");
    leg_V2_low_la->AddEntry(V2_bkg_low_la, "Bkg", "P");
    leg_V2_low_la->Draw();
    V2_obs_low_la->Draw("P");
    V2_bkg_low_la->Draw("P");

    frame_V2_ratio_low_la = V2_la->cd(4)->DrawFrame(0,0.0,9,1.5);
    V2_la->cd(4)->SetLeftMargin(bottomMargin);
    V2_la->cd(4)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_V2_ratio_low_la->GetXaxis()->CenterTitle(1);
    frame_V2_ratio_low_la->GetYaxis()->CenterTitle(1);
    frame_V2_ratio_low_la->GetXaxis()->SetTitleSize(0.05);
    frame_V2_ratio_low_la->GetYaxis()->SetTitleSize(0.05);
    frame_V2_ratio_low_la->SetTitleOffset(1.1,"Y");
    frame_V2_ratio_low_la->SetTitleOffset(1.2,"X");
    frame_V2_ratio_low_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_V2_ratio_low_la->GetYaxis()->SetTitle("MB #Lambda/#bar{#Lambda} V_{2}^{obs}/V_{2}^{bkg} ");
    TGV2_low_ratio_la->Draw("P");
    line_ratio->Draw("same");

   Nass_ks->Print("Image/RapSys_perisub/Nass_ks.pdf");
   Nass_la->Print("Image/RapSys_perisub/Nass_la.pdf");
   
   Yield_ks->Print("Image/RapSys_perisub/Yield_ks.pdf");
   Yield_la->Print("Image/RapSys_perisub/Yield_la.pdf");

   V2_ks->Print("Image/RapSys_perisub/V2_ks.pdf");
   V2_la->Print("Image/RapSys_perisub/V2_la.pdf");

    /*
    //v2 obs and bkg
    TGraphErrors* v2_bkg_ks        = (TGraphErrors*)file_pidv2->Get("v2plot_bkg_ks");
    TGraphErrors* v2_obs_ks        = (TGraphErrors*)file_pidv2->Get("v2plot_obs_ks");
    TGraphErrors* v2_bkg_low_ks    = (TGraphErrors*)file_pidv2->Get("v2plot_bkg_low_ks");
    TGraphErrors* v2_obs_low_ks    = (TGraphErrors*)file_pidv2->Get("v2plot_obs_low_ks");

    TGraphErrors* v2_bkg_la        = (TGraphErrors*)file_pidv2->Get("v2plot_bkg_la");
    TGraphErrors* v2_obs_la        = (TGraphErrors*)file_pidv2->Get("v2plot_obs_la");
    TGraphErrors* v2_bkg_low_la    = (TGraphErrors*)file_pidv2->Get("v2plot_bkg_low_la");
    TGraphErrors* v2_obs_low_la    = (TGraphErrors*)file_pidv2->Get("v2plot_obs_low_la");

    TGraphErrors* v2_bkg_xi        = (TGraphErrors*)file_pidv2->Get("v2plot_bkg_xi");
    TGraphErrors* v2_obs_xi        = (TGraphErrors*)file_pidv2->Get("v2plot_obs_xi");
    TGraphErrors* v2_bkg_low_xi    = (TGraphErrors*)file_pidv2->Get("v2plot_bkg_low_xi");
    TGraphErrors* v2_obs_low_xi    = (TGraphErrors*)file_pidv2->Get("v2plot_obs_low_xi");

    double* v2_bkg_ks_Y = v2_bkg_ks->GetY();
    double* v2_bkg_ks_EY = v2_bkg_ks->GetEY();
    double* v2_obs_ks_Y = v2_obs_ks->GetY();
    double* v2_obs_ks_EY = v2_obs_ks->GetEY();
    double* v2_bkg_low_ks_Y = v2_bkg_low_ks->GetY();
    double* v2_bkg_low_ks_EY = v2_bkg_low_ks->GetEY();
    double* v2_obs_low_ks_Y = v2_obs_low_ks->GetY();
    double* v2_obs_low_ks_EY = v2_obs_low_ks->GetEY();

    double* v2_bkg_la_Y = v2_bkg_la->GetY();
    double* v2_bkg_la_EY = v2_bkg_la->GetEY();
    double* v2_obs_la_Y = v2_obs_la->GetY();
    double* v2_obs_la_EY = v2_obs_la->GetEY();
    double* v2_bkg_low_la_Y = v2_bkg_low_la->GetY();
    double* v2_bkg_low_la_EY = v2_bkg_low_la->GetEY();
    double* v2_obs_low_la_Y = v2_obs_low_la->GetY();
    double* v2_obs_low_la_EY = v2_obs_low_la->GetEY();

    std::vector<double> v2_high_ratio_ks;
    std::vector<double> v2_high_err_ratio_ks;
    std::vector<double> v2_low_ratio_ks;
    std::vector<double> v2_low_err_ratio_ks;

    std::vector<double> v2_high_ratio_la;
    std::vector<double> v2_high_err_ratio_la;
    std::vector<double> v2_low_ratio_la;
    std::vector<double> v2_low_err_ratio_la;

    //Do obs/bkg
    for(int i=0; i<v2_bkg_ks->GetN(); i++)
    {
        v2_high_ratio_ks.push_back(v2_obs_ks_Y[i]/v2_bkg_ks_Y[i]);
        v2_high_err_ratio_ks.push_back(v2_obs_ks_EY[i]/v2_bkg_ks_Y[i]);
        v2_low_ratio_ks.push_back(v2_obs_low_ks_Y[i]/v2_bkg_low_ks_Y[i]);
        v2_low_err_ratio_ks.push_back(v2_obs_low_ks_EY[i]/v2_bkg_low_ks_Y[i]);
    }

    for(int i=0; i<v2_bkg_la->GetN(); i++)
    {
        v2_high_ratio_la.push_back(v2_obs_la_Y[i]/v2_bkg_la_Y[i]);
        v2_high_err_ratio_la.push_back(v2_obs_la_EY[i]/v2_bkg_la_Y[i]);
        v2_low_ratio_la.push_back(v2_obs_low_la_Y[i]/v2_bkg_low_la_Y[i]);
        v2_low_err_ratio_la.push_back(v2_obs_low_la_EY[i]/v2_bkg_low_la_Y[i]);
    }

    TGraphErrors* TGv2_high_ratio_ks = new TGraphErrors(numKs,pt_ks,&v2_high_ratio_ks[0],0,&v2_high_err_ratio_ks[0]);
    TGraphErrors* TGv2_low_ratio_ks = new TGraphErrors(numKs,pt_ks,&v2_low_ratio_ks[0],0,&v2_low_err_ratio_ks[0]);
    TGraphErrors* TGv2_high_ratio_la = new TGraphErrors(numLa,pt_la,&v2_high_ratio_la[0],0,&v2_high_err_ratio_la[0]);
    TGraphErrors* TGv2_low_ratio_la = new TGraphErrors(numLa,pt_la,&v2_low_ratio_la[0],0,&v2_low_err_ratio_la[0]);

    TGv2_high_ratio_ks->SetMarkerColor(kRed);
    TGv2_high_ratio_ks->SetMarkerStyle(20);
    TGv2_high_ratio_ks->SetMarkerSize(1.2);
    TGv2_high_ratio_ks->SetLineColor(kRed);

    TGv2_low_ratio_ks->SetMarkerColor(kBlue);
    TGv2_low_ratio_ks->SetMarkerStyle(20);
    TGv2_low_ratio_ks->SetMarkerSize(1.2);
    TGv2_low_ratio_ks->SetLineColor(kBlue);

    TGv2_high_ratio_la->SetMarkerColor(kRed);
    TGv2_high_ratio_la->SetMarkerStyle(22);
    TGv2_high_ratio_la->SetMarkerSize(1.2);
    TGv2_high_ratio_la->SetLineColor(kRed);

    TGv2_low_ratio_la->SetMarkerColor(kBlue);
    TGv2_low_ratio_la->SetMarkerStyle(22);
    TGv2_low_ratio_la->SetMarkerSize(1.2);
    TGv2_low_ratio_la->SetLineColor(kBlue);

    v2_bkg_ks        ->SetMarkerColor(kBlue);
    v2_bkg_ks        ->SetMarkerStyle(24);
    v2_bkg_ks        ->SetMarkerSize(1.5);
    v2_bkg_ks        ->SetLineColor(kBlue);
    v2_obs_ks        ->SetMarkerColor(kRed);
    v2_obs_ks        ->SetMarkerStyle(24);
    v2_obs_ks        ->SetMarkerSize(1.5);
    v2_obs_ks        ->SetLineColor(kRed);

    v2_bkg_low_ks    ->SetMarkerColor(kBlue);
    v2_bkg_low_ks    ->SetMarkerStyle(24);
    v2_bkg_low_ks    ->SetMarkerSize(1.5);
    v2_bkg_low_ks    ->SetLineColor(kBlue);
    v2_obs_low_ks    ->SetMarkerColor(kRed);
    v2_obs_low_ks    ->SetMarkerStyle(24);
    v2_obs_low_ks    ->SetMarkerSize(1.5);
    v2_obs_low_ks    ->SetLineColor(kRed);

    v2_bkg_la        ->SetMarkerColor(kBlue);
    v2_bkg_la        ->SetMarkerStyle(26);
    v2_bkg_la        ->SetMarkerSize(1.5);
    v2_bkg_la        ->SetLineColor(kBlue);
    v2_obs_la        ->SetMarkerColor(kRed);
    v2_obs_la        ->SetMarkerStyle(26);
    v2_obs_la        ->SetMarkerSize(1.5);
    v2_obs_la        ->SetLineColor(kRed);

    v2_bkg_low_la    ->SetMarkerColor(kBlue);
    v2_bkg_low_la    ->SetMarkerStyle(26);
    v2_bkg_low_la    ->SetMarkerSize(1.5);
    v2_bkg_low_la    ->SetLineColor(kBlue);
    v2_obs_low_la    ->SetMarkerColor(kRed);
    v2_obs_low_la    ->SetMarkerStyle(26);
    v2_obs_low_la    ->SetMarkerSize(1.5);
    v2_obs_low_la    ->SetLineColor(kRed);

    TCanvas* v2_ks = new TCanvas("v2_ks","v2_ks",1200,900);
    v2_ks->Divide(2,2);

    TCanvas* v2_la = new TCanvas("v2_la","v2_la",1200,900);
    v2_la->Divide(2,2);
    TH1F* frame_v2_ks;
    TH1F* frame_v2_ratio_ks;
    TH1F* frame_v2_low_ks;
    TH1F* frame_v2_ratio_low_ks;
    TH1F* frame_v2_la;
    TH1F* frame_v2_ratio_la;
    TH1F* frame_v2_low_la;
    TH1F* frame_v2_ratio_low_la;

    frame_v2_ks = v2_ks->cd(1)->DrawFrame(0,-0.001,9,0.04);
    v2_ks->cd(1)->SetLeftMargin(bottomMargin);
    v2_ks->cd(1)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_v2_ks->GetXaxis()->CenterTitle(1);
    frame_v2_ks->GetYaxis()->CenterTitle(1);
    frame_v2_ks->GetXaxis()->SetTitleSize(0.05);
    frame_v2_ks->GetYaxis()->SetTitleSize(0.05);
    frame_v2_ks->SetTitleOffset(1.1,"Y");
    frame_v2_ks->SetTitleOffset(1.2,"X");
    frame_v2_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_v2_ks->GetYaxis()->SetTitle("HM K_{S}^{0} V_{2}");
    TLegend* leg_v2_ks = new TLegend(0.17,0.63,0.29,0.88);
    leg_v2_ks->SetFillColor(10);
    leg_v2_ks->SetFillStyle(0);
    leg_v2_ks->SetBorderSize(0);
    leg_v2_ks->SetTextFont(42);
    leg_v2_ks->SetTextSize(0.04);
    leg_v2_ks->AddEntry(v2_obs_ks, "Obs", "P");
    leg_v2_ks->AddEntry(v2_bkg_ks, "Bkg", "P");
    leg_v2_ks->Draw();
    v2_obs_ks->Draw("P");
    v2_bkg_ks->Draw("P");

    frame_v2_ratio_ks = v2_ks->cd(2)->DrawFrame(0,0.0,9,1.5);
    v2_ks->cd(2)->SetLeftMargin(bottomMargin);
    v2_ks->cd(2)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_v2_ratio_ks->GetXaxis()->CenterTitle(1);
    frame_v2_ratio_ks->GetYaxis()->CenterTitle(1);
    frame_v2_ratio_ks->GetXaxis()->SetTitleSize(0.05);
    frame_v2_ratio_ks->GetYaxis()->SetTitleSize(0.05);
    frame_v2_ratio_ks->SetTitleOffset(1.1,"Y");
    frame_v2_ratio_ks->SetTitleOffset(1.2,"X");
    frame_v2_ratio_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_v2_ratio_ks->GetYaxis()->SetTitle("HM K_{S}^{0} v_{2}^{obs}/v_{2}^{bkg} ");
    TGv2_high_ratio_ks->Draw("P");
    line_ratio->Draw("same");

    frame_v2_low_ks = v2_ks->cd(3)->DrawFrame(0,-0.01,9,0.1);
    v2_ks->cd(3)->SetLeftMargin(bottomMargin);
    v2_ks->cd(3)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_v2_low_ks->GetXaxis()->CenterTitle(1);
    frame_v2_low_ks->GetYaxis()->CenterTitle(1);
    frame_v2_low_ks->GetXaxis()->SetTitleSize(0.05);
    frame_v2_low_ks->GetYaxis()->SetTitleSize(0.05);
    frame_v2_low_ks->SetTitleOffset(1.1,"Y");
    frame_v2_low_ks->SetTitleOffset(1.2,"X");
    frame_v2_low_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_v2_low_ks->GetYaxis()->SetTitle("MB K_{S}^{0} v_{2}");
    TLegend* leg_v2_low_ks = new TLegend(0.17,0.63,0.29,0.88);
    leg_v2_low_ks->SetFillColor(10);
    leg_v2_low_ks->SetFillStyle(0);
    leg_v2_low_ks->SetBorderSize(0);
    leg_v2_low_ks->SetTextFont(42);
    leg_v2_low_ks->SetTextSize(0.04);
    leg_v2_low_ks->AddEntry(v2_obs_low_ks, "Obs", "P");
    leg_v2_low_ks->AddEntry(v2_bkg_low_ks, "Bkg", "P");
    leg_v2_low_ks->Draw();
    v2_obs_low_ks->Draw("P");
    v2_bkg_low_ks->Draw("P");

    frame_v2_ratio_low_ks = v2_ks->cd(4)->DrawFrame(0,0.00,9,1.5);
    v2_ks->cd(4)->SetLeftMargin(bottomMargin);
    v2_ks->cd(4)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_v2_ratio_low_ks->GetXaxis()->CenterTitle(1);
    frame_v2_ratio_low_ks->GetYaxis()->CenterTitle(1);
    frame_v2_ratio_low_ks->GetXaxis()->SetTitleSize(0.05);
    frame_v2_ratio_low_ks->GetYaxis()->SetTitleSize(0.05);
    frame_v2_ratio_low_ks->SetTitleOffset(1.1,"Y");
    frame_v2_ratio_low_ks->SetTitleOffset(1.2,"X");
    frame_v2_ratio_low_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_v2_ratio_low_ks->GetYaxis()->SetTitle("MB K_{S}^{0} v_{2}^{obs}/v_{2}^{bkg} ");
    TGv2_low_ratio_ks->Draw("P");
    line_ratio->Draw("same");


    frame_v2_la = v2_la->cd(1)->DrawFrame(0,-0.001,9,0.04);
    v2_la->cd(1)->SetLeftMargin(bottomMargin);
    v2_la->cd(1)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_v2_la->GetXaxis()->CenterTitle(1);
    frame_v2_la->GetYaxis()->CenterTitle(1);
    frame_v2_la->GetXaxis()->SetTitleSize(0.05);
    frame_v2_la->GetYaxis()->SetTitleSize(0.05);
    frame_v2_la->SetTitleOffset(1.1,"Y");
    frame_v2_la->SetTitleOffset(1.2,"X");
    frame_v2_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_v2_la->GetYaxis()->SetTitle("HM #Lambda/#bar{#Lambda} v_{2}");
    TLegend* leg_v2_la = new TLegend(0.17,0.63,0.29,0.88);
    leg_v2_la->SetFillColor(10);
    leg_v2_la->SetFillStyle(0);
    leg_v2_la->SetBorderSize(0);
    leg_v2_la->SetTextFont(42);
    leg_v2_la->SetTextSize(0.04);
    leg_v2_la->AddEntry(v2_obs_la, "Obs", "P");
    leg_v2_la->AddEntry(v2_bkg_la, "Bkg", "P");
    leg_v2_la->Draw();
    v2_obs_la->Draw("P");
    v2_bkg_la->Draw("P");

    frame_v2_ratio_la = v2_la->cd(2)->DrawFrame(0,0.00,9,1.5);
    v2_la->cd(2)->SetLeftMargin(bottomMargin);
    v2_la->cd(2)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_v2_ratio_la->GetXaxis()->CenterTitle(1);
    frame_v2_ratio_la->GetYaxis()->CenterTitle(1);
    frame_v2_ratio_la->GetXaxis()->SetTitleSize(0.05);
    frame_v2_ratio_la->GetYaxis()->SetTitleSize(0.05);
    frame_v2_ratio_la->SetTitleOffset(1.1,"Y");
    frame_v2_ratio_la->SetTitleOffset(1.2,"X");
    frame_v2_ratio_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_v2_ratio_la->GetYaxis()->SetTitle("HM #Lambda/#bar{#Lambda} v_{2}^{obs}/v_{2}^{bkg} ");
    TGv2_high_ratio_la->Draw("P");
    line_ratio->Draw("same");

    frame_v2_low_la = v2_la->cd(3)->DrawFrame(0,-0.01,9,0.1);
    v2_la->cd(3)->SetLeftMargin(bottomMargin);
    v2_la->cd(3)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_v2_low_la->GetXaxis()->CenterTitle(1);
    frame_v2_low_la->GetYaxis()->CenterTitle(1);
    frame_v2_low_la->GetXaxis()->SetTitleSize(0.05);
    frame_v2_low_la->GetYaxis()->SetTitleSize(0.05);
    frame_v2_low_la->SetTitleOffset(1.1,"Y");
    frame_v2_low_la->SetTitleOffset(1.2,"X");
    frame_v2_low_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_v2_low_la->GetYaxis()->SetTitle("MB #Lambda/#bar{#Lambda} v_{2}");
    TLegend* leg_v2_low_la = new TLegend(0.17,0.63,0.29,0.88);
    leg_v2_low_la->SetFillColor(10);
    leg_v2_low_la->SetFillStyle(0);
    leg_v2_low_la->SetBorderSize(0);
    leg_v2_low_la->SetTextFont(42);
    leg_v2_low_la->SetTextSize(0.04);
    leg_v2_low_la->AddEntry(v2_obs_low_la, "Obs", "P");
    leg_v2_low_la->AddEntry(v2_bkg_low_la, "Bkg", "P");
    leg_v2_low_la->Draw();
    v2_obs_low_la->Draw("P");
    v2_bkg_low_la->Draw("P");

    frame_v2_ratio_low_la = v2_la->cd(4)->DrawFrame(0,0.0,9,1.5);
    v2_la->cd(4)->SetLeftMargin(bottomMargin);
    v2_la->cd(4)->SetBottomMargin(bottomMargin);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_v2_ratio_low_la->GetXaxis()->CenterTitle(1);
    frame_v2_ratio_low_la->GetYaxis()->CenterTitle(1);
    frame_v2_ratio_low_la->GetXaxis()->SetTitleSize(0.05);
    frame_v2_ratio_low_la->GetYaxis()->SetTitleSize(0.05);
    frame_v2_ratio_low_la->SetTitleOffset(1.1,"Y");
    frame_v2_ratio_low_la->SetTitleOffset(1.2,"X");
    frame_v2_ratio_low_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_v2_ratio_low_la->GetYaxis()->SetTitle("MB #Lambda/#bar{#Lambda} v_{2}^{obs}/v_{2}^{bkg} ");
    TGv2_low_ratio_la->Draw("P");
    line_ratio->Draw("same");
    */
}
void RapSys_RecoCutsV0()
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

    c1_ks->Print("v2RecoSystematicsKshort.pdf");

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

    c1_la->Print("v2RecoSystematicsLambda.pdf");

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

    c1_co->Print("v2CombinedRecoSystematicsV0.pdf");

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

    TLine* LineRatio_min_ks = new TLine(0,0.995,8,0.995);
    LineRatio_min_ks->SetLineStyle(2);
    LineRatio_min_ks->Draw();

    TLine* LineRatio_max_ks = new TLine(0,1.005,8,1.005);
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

    TLine* LineRatio_min_la = new TLine(0,0.995,8,0.995);
    LineRatio_min_la->SetLineStyle(2);
    LineRatio_min_la->Draw();

    TLine* LineRatio_max_la = new TLine(0,1.005,8,1.005);
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

void RapSys_RecoCutsXi()
{
    MITStyle();
    gStyle->SetTitleAlign(33);
    const int xi_npoints = 8;
    double v2loose_xi[xi_npoints]  = {0.0373282 ,0.073848 ,0.0950396 ,0.120215 ,0.168178 ,0.185988 ,0.209509 ,0.189538};
    double pTloose_xi[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.474};
    double v2loose_xiE[xi_npoints] = {0.0107406 ,0.00562262 ,0.00481564 ,0.00380466 ,0.003604 ,0.00422626 ,0.00607231 ,0.0132913};

    double v2standard_xi[xi_npoints]  = {0.0371521 ,0.0669904 ,0.0877339 ,0.123959 ,0.163467 ,0.183918 ,0.209974 ,0.17582};
    double pTstandard_xi[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.474};
    double v2standard_xiE[xi_npoints] = {0.00932224 ,0.00484125 ,0.0041527 ,0.00326934 ,0.00309064 ,0.00361334 ,0.00517473};

    double v2tight_xi[xi_npoints]  = {0.0386033 ,0.0732162 ,0.08916 ,0.120257 ,0.166137 ,0.180869 ,0.212351 ,0.199091};
    double pTtight_xi[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.474};
    double v2tight_xiE[xi_npoints] = {0.01219 ,0.00620487 ,0.00526521 ,0.00414194 ,0.00391181 ,0.00455734 ,0.0065219, 0.0142387};


    // Pull TGraph for Xihort and lambda

    TGraphErrors* loose_v2_xi = new TGraphErrors(xi_npoints,pTloose_xi,v2loose_xi,0,v2loose_xiE);
    TGraphErrors* standard_v2_xi = new TGraphErrors(xi_npoints,pTstandard_xi,v2standard_xi,0,v2standard_xiE);
	TGraphErrors* tight_v2_xi = new TGraphErrors(xi_npoints,pTtight_xi,v2tight_xi,0,v2tight_xiE);

    standard_v2_xi->SetMarkerColor(kRed);
    standard_v2_xi->SetMarkerStyle(20);
    standard_v2_xi->SetMarkerSize(1);
    standard_v2_xi->SetLineColor(kRed);

    loose_v2_xi->SetMarkerColor(kRed);
    loose_v2_xi->SetMarkerStyle(25);
    loose_v2_xi->SetMarkerSize(1);
    loose_v2_xi->SetLineColor(kRed);

    tight_v2_xi->SetMarkerColor(kRed);
    tight_v2_xi->SetMarkerStyle(26);
    tight_v2_xi->SetMarkerSize(1);
    tight_v2_xi->SetLineColor(kRed);

    TCanvas* c1_xi = MakeCanvas("c1_xi", "Plot_xi");
    c1_xi->cd();
    /*c1_xi->SetLogy();*/
    c1_xi->SetLeftMargin(0.12);

    // draw the frame_xi using a histogram frame_xi

    TH1F* frame_xi = c1_xi->DrawFrame(0,-0.01,8,0.45);
    /*TH1F* frame_xi = c1_xi->DrawFrame(0,0.01,20,1);*/
    gPad->SetTickx();
    gPad->SetTicky();
    frame_xi->SetTitle("K_{S}^{0} Reconstruction Cuts");
    frame_xi->SetTitleSize(0.055,"t");
    frame_xi->GetXaxis()->CenterTitle(1);
    frame_xi->GetYaxis()->CenterTitle(1);
    frame_xi->GetXaxis()->SetTitleSize(0.05);
    frame_xi->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_xi->GetYaxis()->SetTitle("v_{2}^{sig}");
    frame_xi->GetYaxis()->SetTitleSize(0.05);
    frame_xi->SetTitleOffset(1.1,"Y");
    frame_xi->SetTitleOffset(1.2,"X");

    //ha_v2->Draw("PESAME");
    standard_v2_xi->Draw("P");
    tight_v2_xi->Draw("P");
    loose_v2_xi->Draw("P");

    TLegend* leg_xi = new TLegend(0.15,0.55,0.27,0.75);
    leg_xi->SetFillColor(10);
    leg_xi->SetFillStyle(0);
    leg_xi->SetBorderSize(0);
    leg_xi->SetTextFont(42);
    leg_xi->SetTextSize(0.03);
    leg_xi->AddEntry(standard_v2_xi, "Standard reconstruction", "P");
    leg_xi->AddEntry(tight_v2_xi, "Tight reconstruction", "P");
    leg_xi->AddEntry(loose_v2_xi, "Loose reconstruction", "P");
    leg_xi->Draw();

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(62);
    tex->SetTextSize(0.045);
    tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
    tex->SetTextSize(0.04);
    tex->SetTextFont(42);
    tex->DrawLatex(0.40,0.2,"185 #leq N_{trk}^{offline} < 250");
    TLine* line = new TLine(0,0,6,0);
    line->SetLineStyle(2);
    line->Draw("same");

    c1_xi->Print("v2RecoSystematicsXi.pdf");

    //Calculate Ratios
    std::vector<double> v2RatioLooseXi;
    std::vector<double> v2RatioTightXi;

    std::vector<double> v2ErrorRatioLooseXi;
    std::vector<double> v2ErrorRatioTightXi;

    //Xi
    for(int i=0; i<xi_npoints; i++)
    {
        v2RatioLooseXi.push_back(v2loose_xi[i]/v2standard_xi[i]);
        v2ErrorRatioLooseXi.push_back(v2loose_xiE[i]/v2standard_xi[i]);

        v2RatioTightXi.push_back(v2tight_xi[i]/v2standard_xi[i]);
        v2ErrorRatioTightXi.push_back(v2tight_xiE[i]/v2standard_xi[i]);
    }

    double* av2RatioLooseXi = &v2RatioLooseXi[0];
    double* av2RatioTightXi = &v2RatioTightXi[0];

    double* av2ErrorRatioLooseXi = &v2ErrorRatioLooseXi[0];
    double* av2ErrorRatioTightXi = &v2ErrorRatioTightXi[0];

    TGraphErrors* Ratioloose_v2_xi = new TGraphErrors(xi_npoints,pTloose_xi,av2RatioLooseXi,0,av2ErrorRatioLooseXi);
	TGraphErrors* Ratiotight_v2_xi = new TGraphErrors(xi_npoints,pTtight_xi,av2RatioTightXi,0,av2ErrorRatioTightXi);

    Ratioloose_v2_xi->SetMarkerColor(kRed);
    Ratioloose_v2_xi->SetMarkerStyle(20);
    Ratioloose_v2_xi->SetMarkerSize(1.5);
    Ratioloose_v2_xi->SetLineColor(kRed);

    Ratiotight_v2_xi->SetMarkerColor(kBlue-4);
    Ratiotight_v2_xi->SetMarkerStyle(21);
    Ratiotight_v2_xi->SetMarkerSize(1.5);
    Ratiotight_v2_xi->SetLineColor(kBlue-4);

    TCanvas* c1_ratio_xi = MakeCanvas("c1_ratio_xi", "Plot_ratio_xi");
    TH1F* frame_Ratio_xi = c1_ratio_xi->DrawFrame(0,0.85,8,1.2);
    /*TH1F* frame_co = c1_co->DrawFrame(0,0.01,20,1);*/
    gPad->SetTickx();
    gPad->SetTicky();
    /*frame_Ratio_xi->SetTitle("K_{S}^{-1} Reconstruction Cut Ratio");*/
    frame_Ratio_xi->GetXaxis()->CenterTitle(1);
    frame_Ratio_xi->GetYaxis()->CenterTitle(1);
    frame_Ratio_xi->GetXaxis()->SetTitleSize(0.05);
    frame_Ratio_xi->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Ratio_xi->GetYaxis()->SetTitle("v_{2}^{sig}");
    frame_Ratio_xi->GetYaxis()->SetTitleSize(0.05);
    frame_Ratio_xi->SetTitleOffset(1.2,"Y");
    frame_Ratio_xi->SetTitleOffset(1.2,"X");

    TLine* LineRatio_xi = new TLine(0,1,8,1);
    /*LineRatio_xi->SetLineStyle(2);*/
    LineRatio_xi->Draw();

    TLine* LineRatio_min_xi = new TLine(0,0.995,8,0.995);
    LineRatio_min_xi->SetLineStyle(2);
    LineRatio_min_xi->Draw();

    TLine* LineRatio_max_xi = new TLine(0,1.005,8,1.005);
    LineRatio_max_xi->SetLineStyle(2);
    LineRatio_max_xi->Draw();

    Ratioloose_v2_xi->Draw("P");
    Ratiotight_v2_xi->Draw("P");

    TLegend* leg_ratio_xi = new TLegend(0.45,0.2,0.57,0.3);
    leg_ratio_xi->SetFillColor(10);
    leg_ratio_xi->SetFillStyle(0);
    leg_ratio_xi->SetBorderSize(0);
    leg_ratio_xi->SetTextFont(42);
    leg_ratio_xi->SetTextSize(0.03);
    leg_ratio_xi->AddEntry(Ratioloose_v2_xi, "Loose #Xi^{#pm} Reconstruction", "P");
    leg_ratio_xi->AddEntry(Ratiotight_v2_xi, "Tight #Xi^{#pm} Reconstruction", "P");
    leg_ratio_xi->Draw();

    c1_ratio_xi->Print("v2RecoRatioSystematicsXi.pdf");

    // Write points into rootfile
    TFile out("V0ReconstructionSys.root","RECREATE");
    standard_v2_xi->Write("StandardXi");
    tight_v2_xi->Write("TightXi");
    loose_v2_xi->Write("LooseXi");
    Ratioloose_v2_xi->Write("RatioLooseXi");
    Ratiotight_v2_xi->Write("RatioTightXi");
}

void RapSys_Closure_Study()
{
    MITStyle();
    gStyle->SetTitleAlign(33);
    TVirtualFitter::SetMaxIterations( 300000 );
    const int ks_npoints = 13;
    const int la_npoints = 10;

    //TFile* file_pidv2_match_recoRef = TFile::Open("rootFiles/v2valuesRapidityClosure_MatchV0ClosureBpPb_09_20_17.root"); //Match V0 w/ reco ref
    TFile* file_pidv2_match_recoRef = TFile::Open("rootFiles/v2valuesRapidityClosure_V0ClosureGen_RecoRef_10_23_17.root"); //Gen V0 w/ reco ref
    TFile* file_pidv2_reco_genRef = TFile::Open("rootFiles/v2valuesRapidityClosure_V0ClosureReco_GenRef_10_24_17.root"); //Reco V0 w/ gen ref
    TFile* file_pidv2_reco_recoRef = TFile::Open("rootFiles/v2valuesRapidityClosure_V0CorrelationClosureTotal_08_28_17.root"); //Reco V0 w/ reco ref
    //TFile* file_pidv2_gen = TFile::Open("rootFiles/v2valuesRapidityClosure_V0CorrelationClosureGenTotal_08_28_17.root"); // Gen w/ gen ref
    TFile* file_pidv2_gen = TFile::Open("rootFiles/v2valuesRapidityClosure_V0ClosureGenNoStrange_10_25_17.root"); // Gen w/ gen ref




    // Pull TGraph for Kshort and lambda
    TGraphErrors* Reco_recoRef_v2_ks = (TGraphErrors*)file_pidv2_reco_recoRef->Get("v2kshort");
    TGraphErrors* RecoMatch_v2_ks = (TGraphErrors*)file_pidv2_match_recoRef->Get("v2kshort");
    TGraphErrors* Reco_genRef_v2_ks = (TGraphErrors*)file_pidv2_reco_genRef->Get("v2kshort");
    TGraphErrors* Gen_v2_ks =(TGraphErrors*)file_pidv2_gen->Get("v2kshort");

    Gen_v2_ks->SetMarkerColor(kRed);
    Gen_v2_ks->SetMarkerStyle(20);
    Gen_v2_ks->SetMarkerSize(1.5);
    Gen_v2_ks->SetLineColor(kRed);

    Reco_genRef_v2_ks->SetMarkerColor(kBlue-4);
    Reco_genRef_v2_ks->SetMarkerStyle(25);
    Reco_genRef_v2_ks->SetMarkerSize(1.5);
    Reco_genRef_v2_ks->SetLineColor(kBlue-4);

    RecoMatch_v2_ks->SetMarkerColor(kGreen);
    RecoMatch_v2_ks->SetMarkerStyle(26);
    RecoMatch_v2_ks->SetMarkerSize(1.5);
    RecoMatch_v2_ks->SetLineColor(kGreen);

    Reco_recoRef_v2_ks->SetMarkerColor(kMagenta);
    Reco_recoRef_v2_ks->SetMarkerStyle(27);
    Reco_recoRef_v2_ks->SetMarkerSize(1.5);
    Reco_recoRef_v2_ks->SetLineColor(kMagenta);

    TCanvas* c1_ks = MakeCanvas("c1_ks", "Plot_ks");
    c1_ks->cd();
    /*c1_ks->SetLogy();*/
    c1_ks->SetLeftMargin(0.12);

    // draw the frame_ks using a histogram frame_ks

    TH1F* frame_ks = c1_ks->DrawFrame(0,-0.08,8.5,1.00);
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

    TF1* ksFit = new TF1("ksFit","([0]/(1 + exp(-(x-[1])/[2])) - [3])*pol1(4) + [5]*pol2(6)",0,8);
    ksFit->SetParameter(0,1);
    ksFit->SetParameter(1,1);
    ksFit->SetParameter(2,1);
    ksFit->SetParameter(3,1);
    ksFit->SetParameter(4,1);
    ksFit->SetParameter(5,1);
    ksFit->SetNpx(250);
    ksFit->SetLineColor(kRed);
    ksFit->SetLineStyle(2);

    Gen_v2_ks->Fit("ksFit");


    /*double parhold[6];*/

    //ha_v2->Draw("PESAME");
    Gen_v2_ks->Draw("P");
    Reco_recoRef_v2_ks->Draw("P");
    RecoMatch_v2_ks->Draw("P");
    Reco_genRef_v2_ks->Draw("P");

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045);
    tex->SetTextFont(42);
    TLine* line = new TLine(0,0,6,0);
    line->SetLineStyle(2);
    line->Draw("same");

    TLegend* leg_ks = new TLegend(0.15,0.65,0.3,0.75);
    leg_ks->SetFillColor(10);
    leg_ks->SetFillStyle(0);
    leg_ks->SetBorderSize(0);
    leg_ks->SetTextFont(42);
    leg_ks->SetTextSize(0.03);
    leg_ks->AddEntry(Gen_v2_ks, "Gen K_{S}^{0} w/ Gen Ref", "P");
    leg_ks->AddEntry(Reco_recoRef_v2_ks, "Reco K_{S}^{0} w/ Reco Ref", "P");
    leg_ks->AddEntry(RecoMatch_v2_ks, "Gen K_{S}^{0} w/ Reco Ref", "P");
    leg_ks->AddEntry(Reco_genRef_v2_ks, "Reco K_{S}^{0} w/ Gen Ref", "P");
    leg_ks->Draw();

    tex->SetTextFont(62);
    tex->SetTextSize(0.045);
    tex->DrawLatex(0.15,0.8,"pPb EPOS");

    c1_ks->Print("v2ClosureSystematicsKshort.pdf");

    //Do lambda


    // Pull TGraph for Kshort and lambda

    TGraphErrors* Reco_recoRef_v2_la = (TGraphErrors*)file_pidv2_reco_recoRef->Get("v2lambda");
    TGraphErrors* RecoMatch_v2_la = (TGraphErrors*)file_pidv2_match_recoRef->Get("v2lambda");
    TGraphErrors* Reco_genRef_v2_la = (TGraphErrors*)file_pidv2_reco_genRef->Get("v2lambda");
    TGraphErrors* Gen_v2_la =(TGraphErrors*)file_pidv2_gen->Get("v2lambda");

    Gen_v2_la->SetMarkerColor(kGreen-2);
    Gen_v2_la->SetMarkerStyle(29);
    Gen_v2_la->SetMarkerSize(1.5);
    Gen_v2_la->SetLineColor(kGreen-2);

    Reco_genRef_v2_la->SetMarkerColor(kBlue-4);
    Reco_genRef_v2_la->SetMarkerStyle(25);
    Reco_genRef_v2_la->SetMarkerSize(1.5);
    Reco_genRef_v2_la->SetLineColor(kBlue-4);

    RecoMatch_v2_la->SetMarkerColor(kGreen);
    RecoMatch_v2_la->SetMarkerStyle(26);
    RecoMatch_v2_la->SetMarkerSize(1.5);
    RecoMatch_v2_la->SetLineColor(kGreen);

    Reco_recoRef_v2_la->SetMarkerColor(kMagenta);
    Reco_recoRef_v2_la->SetMarkerStyle(27);
    Reco_recoRef_v2_la->SetMarkerSize(1.5);
    Reco_recoRef_v2_la->SetLineColor(kMagenta);


    TCanvas* c1_la = MakeCanvas("c1_la", "Plot_la");
    c1_la->cd();
    /*c1_la->SetLogy();*/
    c1_la->SetLeftMargin(0.12);

    // draw the frame_la using a histogram frame_la

    TH1F* frame_la = c1_la->DrawFrame(0,-0.15,8.5,1.0);
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

    TF1* laFit = new TF1("laFit","([0]/(1 + exp(-(x-[1])/[2])) - [3])*pol1(4) + [5]*pol2(6)",0,8);
    laFit->SetParameter(0,1);
    laFit->SetParameter(1,1);
    laFit->SetParameter(2,1);
    laFit->SetParameter(3,1);
    laFit->SetParameter(4,1);
    laFit->SetParameter(5,1);
    laFit->SetNpx(250);
    laFit->SetLineColor(kBlue-4);
    laFit->SetLineStyle(2);

    Gen_v2_la->Fit("laFit");

    //ha_v2->Draw("PESAME");
    Gen_v2_la->Draw("P");
    Reco_recoRef_v2_la->Draw("P");
    RecoMatch_v2_la->Draw("P");
    Reco_genRef_v2_la->Draw("P");

    TLegend* leg_la = new TLegend(0.15,0.65,0.3,0.75);
    leg_la->SetFillColor(10);
    leg_la->SetFillStyle(0);
    leg_la->SetBorderSize(0);
    leg_la->SetTextFont(42);
    leg_la->SetTextSize(0.03);
    leg_la->AddEntry(Gen_v2_la, "Gen #Lambda/#bar{#Lambda} w/ Gen Ref", "P");
    leg_la->AddEntry(Reco_recoRef_v2_la, "Reco #Lambda/#bar{#Lambda} w/ Reco Ref", "P");
    leg_la->AddEntry(RecoMatch_v2_la, "Gen #Lambda/#bar{#Lambda} w/ Reco Ref", "P");
    leg_la->AddEntry(Reco_genRef_v2_la, "Reco #Lambda/#bar{#Lambda} w/ Gen Ref", "P" );
    leg_la->Draw();

    tex->SetTextFont(62);
    tex->SetTextSize(0.04);
    tex->DrawLatex(0.15,0.8,"pPb EPOS");
    tex->SetTextSize(0.04);
    tex->SetTextFont(42);
    line->Draw("same");

    c1_la->Print("v2ClosureSystematicsLambda.pdf");

    //Combined
    TCanvas* c1_co = MakeCanvas("c1_co", "Plot_co");
    c1_co->cd();
    /*c1_co->SetLogy();*/
    c1_co->SetLeftMargin(0.12);

    // draw the frame_co using a histogram frame_co

    TH1F* frame_co = c1_co->DrawFrame(0,-0.15,8.5,1.00);
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

    Gen_v2_ks->Draw("P");
    Reco_recoRef_v2_ks->Draw("P");
    RecoMatch_v2_ks->Draw("P");
    Reco_genRef_v2_ks->Draw("P");

    Gen_v2_la->Draw("P");
    Reco_recoRef_v2_la->Draw("P");
    RecoMatch_v2_la->Draw("P");
    Reco_genRef_v2_la->Draw("P");

    TLegend* leg_co1 = new TLegend(0.15,0.65,0.27,0.75);
    leg_co1->SetFillColor(10);
    leg_co1->SetFillStyle(0);
    leg_co1->SetBorderSize(0);
    leg_co1->SetTextFont(42);
    leg_co1->SetTextSize(0.04);
    leg_co1->AddEntry(Gen_v2_ks, "K_{S}^{0}", "P");
    leg_co1->AddEntry(Gen_v2_la, "#Lambda/#bar{#Lambda}", "P");
    leg_co1->Draw();

    TLegend* leg_co2 = new TLegend(0.6,0.75,0.85,0.9);
    leg_co2->SetFillColor(10);
    leg_co2->SetFillStyle(0);
    leg_co2->SetBorderSize(0);
    leg_co2->SetTextFont(42);
    leg_co2->SetTextSize(0.03);
    leg_co2->AddEntry(Gen_v2_la, "Gen w/ Gen Ref", "P");
    leg_co2->AddEntry(Reco_recoRef_v2_la, "Reco w/ Reco Ref", "P");
    leg_co2->AddEntry(RecoMatch_v2_la, "Gen w/ Reco Ref", "P");
    leg_co2->AddEntry(Reco_genRef_v2_la, "Reco w/ Gen Ref", "P");
    leg_co2->Draw();

    tex->SetTextFont(62);
    tex->SetTextSize(0.045);
    tex->DrawLatex(0.15,0.8,"pPb EPOS");
    tex->SetTextSize(0.04);
    tex->SetTextFont(42);
    line->Draw("same");

    c1_co->Print("v2CombinedClosureSystematics.pdf");

    //Calculate Ratios Reco/Gen
    std::vector<double> v2RatioKsReco;
    std::vector<double> v2RatioLaReco;
    std::vector<double> v2RatioKsRecoMatch;
    std::vector<double> v2RatioLaRecoMatch;

    std::vector<double> v2RatioKsGen;
    std::vector<double> v2RatioLaGen;

    std::vector<double> v2RatioKsGenRecoref;
    std::vector<double> v2RatioLaGenRecoref;

    std::vector<double> v2ErrorRatioKsReco;
    std::vector<double> v2ErrorRatioLaReco;
    std::vector<double> v2ErrorRatioKsRecoMatch;
    std::vector<double> v2ErrorRatioLaRecoMatch;

    std::vector<double> v2ErrorRatioKsGen;
    std::vector<double> v2ErrorRatioLaGen;

    std::vector<double> v2ErrorRatioKsGenRecoref;
    std::vector<double> v2ErrorRatioLaGenRecoref;

    double *aReco_recoRef_v2_ksX = Reco_recoRef_v2_ks->GetX();
    double *aReco_recoRef_v2_ksY = Reco_recoRef_v2_ks->GetY();
    double *aReco_recoRef_v2_ksEY = Reco_recoRef_v2_ks->GetEY();

    double *aRecoMatch_v2_ksX = RecoMatch_v2_ks->GetX();
    double *aRecoMatch_v2_ksY = RecoMatch_v2_ks->GetY();
    double *aRecoMatch_v2_ksEY = RecoMatch_v2_ks->GetEY();

    double *aGen_v2_ksX = Gen_v2_ks->GetX();
    double *aGen_v2_ksY = Gen_v2_ks->GetY();
    double *aGen_v2_ksEY = Gen_v2_ks->GetEY();

    double *aReco_genRef_v2_ksX = Reco_genRef_v2_ks->GetX();
    double *aReco_genRef_v2_ksY = Reco_genRef_v2_ks->GetY();
    double *aReco_genRef_v2_ksEY = Reco_genRef_v2_ks->GetEY();

    double *aReco_recoRef_v2_laX = Reco_recoRef_v2_la->GetX();
    double *aReco_recoRef_v2_laY = Reco_recoRef_v2_la->GetY();
    double *aReco_recoRef_v2_laEY = Reco_recoRef_v2_la->GetEY();

    double *aRecoMatch_v2_laX = RecoMatch_v2_la->GetX();
    double *aRecoMatch_v2_laY = RecoMatch_v2_la->GetY();
    double *aRecoMatch_v2_laEY = RecoMatch_v2_la->GetEY();

    double *aGen_v2_laX = Gen_v2_la->GetX();
    double *aGen_v2_laY = Gen_v2_la->GetY();
    double *aGen_v2_laEY = Gen_v2_la->GetEY();

    double *aReco_genRef_v2_laX = Reco_genRef_v2_la->GetX();
    double *aReco_genRef_v2_laY = Reco_genRef_v2_la->GetY();
    double *aReco_genRef_v2_laEY = Reco_genRef_v2_la->GetEY();


    //Kshort
    for(int i=0; i<Gen_v2_ks->GetN(); i++)
    {
        v2RatioKsReco.push_back(aReco_recoRef_v2_ksY[i]/ksFit->Eval(aReco_recoRef_v2_ksX[i]));
        v2ErrorRatioKsReco.push_back(aReco_recoRef_v2_ksEY[i]/ksFit->Eval(aReco_recoRef_v2_ksX[i]));

        v2RatioKsRecoMatch.push_back(aRecoMatch_v2_ksY[i]/ksFit->Eval(aRecoMatch_v2_ksX[i]));
        v2ErrorRatioKsRecoMatch.push_back(aRecoMatch_v2_ksEY[i]/ksFit->Eval(aRecoMatch_v2_ksX[i]));

        v2RatioKsGen.push_back(aGen_v2_ksY[i]/ksFit->Eval(aGen_v2_ksX[i]));
        v2ErrorRatioKsGen.push_back(aGen_v2_ksEY[i]/ksFit->Eval(aGen_v2_ksX[i]));

        v2RatioKsGenRecoref.push_back(aReco_genRef_v2_ksY[i]/ksFit->Eval(aReco_genRef_v2_ksX[i]));
        v2ErrorRatioKsGenRecoref.push_back(aReco_genRef_v2_ksEY[i]/ksFit->Eval(aReco_genRef_v2_ksX[i]));

        cout << "Ratio " << i << ": " << v2RatioKsReco[i] << endl;
    }

    //Lambda
    for(int i=0; i<Gen_v2_la->GetN(); i++)
    {
        v2RatioLaReco.push_back(aReco_recoRef_v2_laY[i]/laFit->Eval(aReco_recoRef_v2_laX[i]));
        v2ErrorRatioLaReco.push_back(aReco_recoRef_v2_laEY[i]/laFit->Eval(aReco_recoRef_v2_laX[i]));

        v2RatioLaRecoMatch.push_back(aRecoMatch_v2_laY[i]/laFit->Eval(aRecoMatch_v2_laX[i]));
        v2ErrorRatioLaRecoMatch.push_back(aRecoMatch_v2_laEY[i]/laFit->Eval(aRecoMatch_v2_laX[i]));

        v2RatioLaGen.push_back(aGen_v2_laY[i]/laFit->Eval(aGen_v2_laX[i]));
        v2ErrorRatioLaGen.push_back(aGen_v2_laEY[i]/laFit->Eval(aGen_v2_laX[i]));

        v2RatioLaGenRecoref.push_back(aReco_genRef_v2_laY[i]/laFit->Eval(aReco_genRef_v2_laX[i]));
        v2ErrorRatioLaGenRecoref.push_back(aReco_genRef_v2_laEY[i]/laFit->Eval(aReco_genRef_v2_laX[i]));
    }

    TGraphErrors* RatioReco_v2_la = new TGraphErrors(la_npoints,aReco_recoRef_v2_laX,&v2RatioLaReco[0],0,&v2ErrorRatioLaReco[0]);
    TGraphErrors* RatioRecoMatch_recoref_v2_la = new TGraphErrors(la_npoints,aReco_genRef_v2_laX,&v2RatioLaRecoMatch[0],0,&v2ErrorRatioLaRecoMatch[0]);
    TGraphErrors* RatioReco_genref_v2_la = new TGraphErrors(la_npoints,aReco_genRef_v2_laX,&v2RatioLaGenRecoref[0],0,&v2ErrorRatioLaGenRecoref[0]);
    TGraphErrors* RatioGen_v2_la = new TGraphErrors(la_npoints,aGen_v2_laX,&v2RatioLaGen[0],0,&v2ErrorRatioLaGen[0]);

    TGraphErrors* RatioReco_v2_ks = new TGraphErrors(ks_npoints,aReco_recoRef_v2_ksX,&v2RatioKsReco[0],0,&v2ErrorRatioKsReco[0]);
    TGraphErrors* RatioRecoMatch_recoref_v2_ks = new TGraphErrors(ks_npoints,aReco_genRef_v2_ksX,&v2RatioKsRecoMatch[0],0,&v2ErrorRatioKsRecoMatch[0]);
    TGraphErrors* RatioReco_genref_v2_ks = new TGraphErrors(ks_npoints,aReco_genRef_v2_ksX,&v2RatioKsGenRecoref[0],0,&v2ErrorRatioKsGenRecoref[0]);
    TGraphErrors* RatioGen_v2_ks = new TGraphErrors(ks_npoints,aGen_v2_ksX,&v2RatioKsGen[0],0,&v2ErrorRatioKsGen[0]);

    RatioReco_v2_ks->SetMarkerColor(kRed);
    RatioReco_v2_ks->SetMarkerStyle(25);
    RatioReco_v2_ks->SetMarkerSize(1.3);
    RatioReco_v2_ks->SetLineColor(kRed);

    RatioRecoMatch_recoref_v2_ks->SetMarkerColor(kRed);
    RatioRecoMatch_recoref_v2_ks->SetMarkerStyle(26);
    RatioRecoMatch_recoref_v2_ks->SetMarkerSize(1.3);
    RatioRecoMatch_recoref_v2_ks->SetLineColor(kRed);

    RatioReco_genref_v2_ks->SetMarkerColor(kRed);
    RatioReco_genref_v2_ks->SetMarkerStyle(27);
    RatioReco_genref_v2_ks->SetMarkerSize(1.3);
    RatioReco_genref_v2_ks->SetLineColor(kRed);

    RatioGen_v2_ks->SetMarkerColor(kRed);
    RatioGen_v2_ks->SetMarkerStyle(20);
    RatioGen_v2_ks->SetMarkerSize(1.3);
    RatioGen_v2_ks->SetLineColor(kRed);

    RatioReco_v2_la->SetMarkerColor(kBlue-4);
    RatioReco_v2_la->SetMarkerStyle(25);
    RatioReco_v2_la->SetMarkerSize(1.3);
    RatioReco_v2_la->SetLineColor(kBlue-4);

    RatioRecoMatch_recoref_v2_la->SetMarkerColor(kBlue-4);
    RatioRecoMatch_recoref_v2_la->SetMarkerStyle(26);
    RatioRecoMatch_recoref_v2_la->SetMarkerSize(1.3);
    RatioRecoMatch_recoref_v2_la->SetLineColor(kBlue-4);

    RatioReco_genref_v2_la->SetMarkerColor(kBlue-4);
    RatioReco_genref_v2_la->SetMarkerStyle(27);
    RatioReco_genref_v2_la->SetMarkerSize(1.3);
    RatioReco_genref_v2_la->SetLineColor(kBlue-4);

    RatioGen_v2_la->SetMarkerColor(kBlue-4);
    RatioGen_v2_la->SetMarkerStyle(20);
    RatioGen_v2_la->SetMarkerSize(1.3);
    RatioGen_v2_la->SetLineColor(kBlue-4);

    TCanvas* c1_ratio = MakeCanvas("c1_ratio", "Plot_ratio");
    TH1F* frame_Ratio = c1_ratio->DrawFrame(0,0.8,8.5,1.2);
    //TH1F* frame_Ratio = c1_ratio->DrawFrame(0,0,8.5,2);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Ratio->GetXaxis()->CenterTitle(1);
    frame_Ratio->GetYaxis()->CenterTitle(1);
    frame_Ratio->GetXaxis()->SetTitleSize(0.05);
    frame_Ratio->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Ratio->GetYaxis()->SetTitle("#frac{Combination v^{sig}_{2}}{Gen/Gen v^{sig}_{2}}");
    frame_Ratio->GetYaxis()->SetTitleSize(0.05);
    frame_Ratio->SetTitleOffset(1.2,"Y");
    frame_Ratio->SetTitleOffset(1.2,"X");

    TLine* LineRatio_ks = new TLine(0,1,8.5,1);
    LineRatio_ks->Draw();

    TLine* LineRatio_min_ks = new TLine(0,0.95,8.5,0.95);
    LineRatio_min_ks->SetLineStyle(2);
    LineRatio_min_ks->Draw();

    TLine* LineRatio_max_ks = new TLine(0,1.05,8.5,1.05);
    LineRatio_max_ks->SetLineStyle(2);
    LineRatio_max_ks->Draw();

    RatioReco_v2_ks->Draw("P");
    RatioRecoMatch_recoref_v2_ks->Draw("P");
    RatioReco_genref_v2_ks->Draw("P");
    RatioGen_v2_ks->Draw("P");
    RatioGen_v2_la->Draw("P");
    RatioReco_v2_la->Draw("P");
    RatioRecoMatch_recoref_v2_la->Draw("P");
    RatioReco_genref_v2_la->Draw("P");

    TLegend* leg_ratio_species = new TLegend(0.80,0.75,0.90,0.85);
    leg_ratio_species->SetFillColor(10);
    leg_ratio_species->SetFillStyle(0);
    leg_ratio_species->SetBorderSize(0);
    leg_ratio_species->SetTextFont(42);
    leg_ratio_species->SetTextSize(0.03);
    leg_ratio_species->AddEntry(RatioGen_v2_ks, "K_{S}^{0}", "P");
    leg_ratio_species->AddEntry(RatioGen_v2_la, "#Lambda/#bar{#Lambda}", "P");
    leg_ratio_species->Draw();

    TGraphErrors* RatioGen_v2_ks_clone = (TGraphErrors*)RatioGen_v2_ks->Clone("RatioGen_v2_ks_clone");
    RatioGen_v2_ks_clone->SetMarkerColor(kBlack);
    TGraphErrors* RatioReco_v2_ks_clone = (TGraphErrors*)RatioReco_v2_ks->Clone("RatioReco_v2_ks_clone");
    RatioReco_v2_ks_clone->SetMarkerColor(kBlack);
    TGraphErrors* RatioRecoMatch_recoref_v2_ks_clone = (TGraphErrors*)RatioRecoMatch_recoref_v2_ks->Clone("RatioRecoMatch_recoref_v2_ks_clone");
    RatioRecoMatch_recoref_v2_ks_clone->SetMarkerColor(kBlack);
    TGraphErrors* RatioReco_genref_v2_ks_clone = (TGraphErrors*)RatioReco_genref_v2_ks->Clone("RatioReco_genref_v2_ks_clone");
    RatioReco_genref_v2_ks_clone->SetMarkerColor(kBlack);
    TGraphErrors* refpoint = (TGraphErrors*)RatioGen_v2_ks->Clone("ref");
    refpoint->SetMarkerStyle(1);

    TLegend* leg_ratio_label = new TLegend(0.64,0.70,0.74,0.85);
    leg_ratio_label->SetFillColor(10);
    leg_ratio_label->SetFillStyle(0);
    leg_ratio_label->SetBorderSize(0);
    leg_ratio_label->SetTextFont(42);
    leg_ratio_label->SetTextSize(0.03);
    leg_ratio_label->AddEntry(refpoint, "V0 / Ref", "P");
    leg_ratio_label->AddEntry(RatioGen_v2_ks_clone, "Gen/Gen", "P");
    leg_ratio_label->AddEntry(RatioReco_v2_ks_clone, "Reco/Reco", "P");
    //leg_ratio_label->AddEntry(RatioRecoMatch_recoref_v2_ks_clone,"Match/Reco","P");
    leg_ratio_label->AddEntry(RatioRecoMatch_recoref_v2_ks_clone,"Gen/Reco","P");
    leg_ratio_label->AddEntry(RatioReco_genref_v2_ks_clone, "Reco/Gen", "P");
    leg_ratio_label->Draw();

    tex->SetTextSize(0.03);
    tex->DrawLatex(0.4,0.85,"pPb EPOS");
    tex->DrawLatex(0.4,0.80,"|#eta| > 1");
    tex->DrawLatex(0.4,0.75,"0.3 < p_{T}^{assoc} < 3.0 GeV");

    //c1_ratio->Print("v2RatioClosureSystematics.pdf");
    c1_ratio->Print("v2RatioClosureSystematics_All_Fix.pdf");
    //c1_ratio->Print("v2RatioClosureSystematics_CheckGen_RecoRef_with_Match_RecoRef.pdf");

     //Write points into rootfile
    //TFile out("V0ClosureSys.root","RECREATE");
    //Gen_v2_la->Write("GenLa");
    //Reco_v2_la->Write("RecoLa");
    //Gen_v2_ks->Write("GenKs");
    //Reco_v2_ks->Write("RecoKs");
    //RatioReco_v2_ks->Write("RatioRecoKs");
    //RatioReco_v2_la->Write("RatioRecoLa");

    TCanvas* c1_ratio_ks = MakeCanvas("ks_ratio", "Plot_ks_ratio");
    TH1F* frame_Ratio_ks = c1_ratio_ks->DrawFrame(0,0.8,8.5,1.2);
    //TH1F* frame_Ratio = c1_ratio->DrawFrame(0,0,8.5,2);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Ratio_ks->GetXaxis()->CenterTitle(1);
    frame_Ratio_ks->GetYaxis()->CenterTitle(1);
    frame_Ratio_ks->GetXaxis()->SetTitleSize(0.05);
    frame_Ratio_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Ratio_ks->GetYaxis()->SetTitle("#frac{Combination v^{sig}_{2}}{Gen/Gen v^{sig}_{2}}");
    frame_Ratio_ks->GetYaxis()->SetTitleSize(0.05);
    frame_Ratio_ks->SetTitleOffset(1.2,"Y");
    frame_Ratio_ks->SetTitleOffset(1.2,"X");

    RatioReco_v2_ks->Draw("P");
    RatioRecoMatch_recoref_v2_ks->Draw("P");
    RatioReco_genref_v2_ks->Draw("P");
    RatioGen_v2_ks->Draw("P");

    TLegend* leg_ratio_label_ks = new TLegend(0.64,0.70,0.74,0.85);
    leg_ratio_label_ks->SetFillColor(10);
    leg_ratio_label_ks->SetFillStyle(0);
    leg_ratio_label_ks->SetBorderSize(0);
    leg_ratio_label_ks->SetTextFont(42);
    leg_ratio_label_ks->SetTextSize(0.03);
    leg_ratio_label_ks->AddEntry(refpoint, "V0 / Ref", "P");
    leg_ratio_label_ks->AddEntry(RatioGen_v2_ks_clone, "Gen/Gen", "P");
    leg_ratio_label_ks->AddEntry(RatioReco_v2_ks_clone, "Reco/Reco", "P");
    leg_ratio_label_ks->AddEntry(RatioRecoMatch_recoref_v2_ks_clone,"Gen/Reco","P");
    leg_ratio_label_ks->AddEntry(RatioReco_genref_v2_ks_clone, "Reco/Gen", "P");
    leg_ratio_label_ks->Draw();

    TLine* line_ks = new TLine(0,1,8.5,1);
    line_ks->Draw();

    tex->SetTextSize(0.035);
    tex->DrawLatex(0.35,0.85,"K_{S}^{0} pPb EPOS");
    tex->DrawLatex(0.35,0.80,"|#eta| > 1");
    tex->DrawLatex(0.35,0.75,"0.3 < p_{T}^{assoc} < 3.0 GeV");

    c1_ratio_ks->Print("v2RatioClosureSystematics_Kshort_Fix.pdf");
    c1_ratio_ks->Print("v2RatioClosureSystematics_Kshort_Fix.png");

    TCanvas* c1_ratio_la = MakeCanvas("la_ratio", "Plot_la_ratio");
    TH1F* frame_Ratio_la = c1_ratio_la->DrawFrame(0,0.8,8.5,1.2);
    //TH1F* frame_Ratio = c1_ratio->DrawFrame(0,0,8.5,2);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Ratio_la->GetXaxis()->CenterTitle(1);
    frame_Ratio_la->GetYaxis()->CenterTitle(1);
    frame_Ratio_la->GetXaxis()->SetTitleSize(0.05);
    frame_Ratio_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Ratio_la->GetYaxis()->SetTitle("#frac{Combination v^{sig}_{2}}{Gen/Gen v^{sig}_{2}}");
    frame_Ratio_la->GetYaxis()->SetTitleSize(0.05);
    frame_Ratio_la->SetTitleOffset(1.2,"Y");
    frame_Ratio_la->SetTitleOffset(1.2,"X");

    RatioReco_v2_la->Draw("P");
    RatioRecoMatch_recoref_v2_la->Draw("P");
    RatioReco_genref_v2_la->Draw("P");
    RatioGen_v2_la->Draw("P");

    TLegend* leg_ratio_label_la = new TLegend(0.64,0.70,0.74,0.85);
    leg_ratio_label_la->SetFillColor(10);
    leg_ratio_label_la->SetFillStyle(0);
    leg_ratio_label_la->SetBorderSize(0);
    leg_ratio_label_la->SetTextFont(42);
    leg_ratio_label_la->SetTextSize(0.03);
    leg_ratio_label_la->AddEntry(refpoint, "V0 / Ref", "P");
    leg_ratio_label_la->AddEntry(RatioGen_v2_ks_clone, "Gen/Gen", "P");
    leg_ratio_label_la->AddEntry(RatioReco_v2_ks_clone, "Reco/Reco", "P");
    leg_ratio_label_la->AddEntry(RatioRecoMatch_recoref_v2_ks_clone,"Gen/Reco","P");
    leg_ratio_label_la->AddEntry(RatioReco_genref_v2_ks_clone, "Reco/Gen", "P");
    leg_ratio_label_la->Draw();

    TLine* line_la = new TLine(0,1,8.5,1);
    line_la->Draw();

    tex->SetTextSize(0.035);
    tex->DrawLatex(0.35,0.85,"#Lambda / #bar{#Lambda} pPb EPOS");
    tex->DrawLatex(0.35,0.80,"|#eta| > 1");
    tex->DrawLatex(0.35,0.75,"0.3 < p_{T}^{assoc} < 3.0 GeV");

    c1_ratio_la->Print("v2RatioClosureSystematics_Lambda_Fix.pdf");
    c1_ratio_la->Print("v2RatioClosureSystematics_Lambda_Fix.png");

}

void RapSys_Closure_plot()
{
    MITStyle();
    gStyle->SetTitleAlign(33);
    TVirtualFitter::SetMaxIterations( 300000 );
    const int ks_npoints = 13;
    const int la_npoints = 10;

    TFile* file_pidv2_reco_recoRef = TFile::Open("rootFiles/v2valuesRapidityClosure_V0CorrelationClosureTotal_08_28_17.root"); //Reco V0 w/ reco ref
    //TFile* file_pidv2_gen = TFile::Open("rootFiles/v2valuesRapidityClosure_V0CorrelationClosureGenTotal_08_28_17.root"); // Gen w/ gen ref
    TFile* file_pidv2_gen = TFile::Open("rootFiles/v2valuesRapidityClosure_V0ClosureGenNoStrange_10_25_17.root"); // Gen w/ gen ref




    // Pull TGraph for Kshort and lambda
    TGraphErrors* Reco_recoRef_v2_ks = (TGraphErrors*)file_pidv2_reco_recoRef->Get("v2kshort");
    TGraphErrors* Gen_v2_ks =(TGraphErrors*)file_pidv2_gen->Get("v2kshort");

    Gen_v2_ks->SetMarkerColor(kRed);
    Gen_v2_ks->SetMarkerStyle(20);
    Gen_v2_ks->SetMarkerSize(1.5);
    Gen_v2_ks->SetLineColor(kRed);

    Reco_recoRef_v2_ks->SetMarkerColor(kRed);
    Reco_recoRef_v2_ks->SetMarkerStyle(24);
    Reco_recoRef_v2_ks->SetMarkerSize(1.5);
    Reco_recoRef_v2_ks->SetLineColor(kRed);

    TCanvas* c1_ks = MakeCanvas("c1_ks", "Plot_ks");
    c1_ks->cd();
    /*c1_ks->SetLogy();*/
    c1_ks->SetLeftMargin(0.12);

    // draw the frame_ks using a histogram frame_ks

    TH1F* frame_ks = c1_ks->DrawFrame(0,-0.08,8.5,1.00);
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

    TF1* ksFit = new TF1("ksFit","([0]/(1 + exp(-(x-[1])/[2])) - [3])*pol1(4) + [5]*pol2(6)",0,8);
    ksFit->SetParameter(0,1);
    ksFit->SetParameter(1,1);
    ksFit->SetParameter(2,1);
    ksFit->SetParameter(3,1);
    ksFit->SetParameter(4,1);
    ksFit->SetParameter(5,1);
    ksFit->SetNpx(250);
    ksFit->SetLineColor(kRed);
    ksFit->SetLineStyle(2);

    Gen_v2_ks->Fit("ksFit");

    /*double parhold[6];*/

    //ha_v2->Draw("PESAME");
    Gen_v2_ks->Draw("P");
    Reco_recoRef_v2_ks->Draw("P");

    TLatex *tex = new TLatex();
    tex->SetNDC();
    tex->SetTextSize(0.045);
    tex->SetTextFont(42);
    TLine* line = new TLine(0,0,6,0);
    line->SetLineStyle(2);
    line->Draw("same");

    TLegend* leg_ks = new TLegend(0.15,0.65,0.3,0.75);
    leg_ks->SetFillColor(10);
    leg_ks->SetFillStyle(0);
    leg_ks->SetBorderSize(0);
    leg_ks->SetTextFont(42);
    leg_ks->SetTextSize(0.03);
    leg_ks->AddEntry(Gen_v2_ks, "Gen K_{S}^{0}", "P");
    leg_ks->AddEntry(Reco_recoRef_v2_ks, "Reco K_{S}^{0}", "P");
    leg_ks->Draw();

    tex->SetTextFont(62);
    tex->SetTextSize(0.045);
    tex->DrawLatex(0.15,0.8,"pPb EPOS");

    c1_ks->Print("v2ClosureSystematicsKshort.pdf");

    //Do lambda


    // Pull TGraph for Kshort and lambda

    TGraphErrors* Reco_recoRef_v2_la = (TGraphErrors*)file_pidv2_reco_recoRef->Get("v2lambda");
    TGraphErrors* Gen_v2_la =(TGraphErrors*)file_pidv2_gen->Get("v2lambda");

    Gen_v2_la->SetMarkerColor(kBlue);
    Gen_v2_la->SetMarkerStyle(22);
    Gen_v2_la->SetMarkerSize(1.5);
    Gen_v2_la->SetLineColor(kBlue);



    Reco_recoRef_v2_la->SetMarkerColor(kBlue);
    Reco_recoRef_v2_la->SetMarkerStyle(26);
    Reco_recoRef_v2_la->SetMarkerSize(1.5);
    Reco_recoRef_v2_la->SetLineColor(kBlue);


    TCanvas* c1_la = MakeCanvas("c1_la", "Plot_la");
    c1_la->cd();
    /*c1_la->SetLogy();*/
    c1_la->SetLeftMargin(0.12);

    // draw the frame_la using a histogram frame_la

    TH1F* frame_la = c1_la->DrawFrame(0,-0.15,8.5,1.0);
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

    TF1* laFit = new TF1("laFit","([0]/(1 + exp(-(x-[1])/[2])) - [3])*pol1(4) + [5]*pol2(6)",-1,9);
    laFit->SetParameter(0,1);
    laFit->SetParameter(1,1);
    laFit->SetParameter(2,1);
    laFit->SetParameter(3,1);
    laFit->SetParameter(4,1);
    laFit->SetParameter(5,1);
    laFit->SetNpx(250);
    laFit->SetLineColor(kBlue-4);
    laFit->SetLineStyle(2);

    Gen_v2_la->Fit("laFit","","",0.5,8.25);

    //ha_v2->Draw("PESAME");
    Gen_v2_la->Draw("P");
    Reco_recoRef_v2_la->Draw("P");

    TLegend* leg_la = new TLegend(0.15,0.65,0.3,0.75);
    leg_la->SetFillColor(10);
    leg_la->SetFillStyle(0);
    leg_la->SetBorderSize(0);
    leg_la->SetTextFont(42);
    leg_la->SetTextSize(0.03);
    leg_la->AddEntry(Gen_v2_la, "Gen #Lambda/#bar{#Lambda}", "P");
    leg_la->AddEntry(Reco_recoRef_v2_la, "Reco #Lambda/#bar{#Lambda}", "P");
    leg_la->Draw();

    tex->SetTextFont(62);
    tex->SetTextSize(0.04);
    tex->DrawLatex(0.15,0.8,"pPb EPOS");
    tex->SetTextSize(0.04);
    tex->SetTextFont(42);
    line->Draw("same");

    c1_la->Print("v2ClosureSystematicsLambda.pdf");

    //Combined
    TCanvas* c1_co = MakeCanvas("c1_co", "Plot_co");
    c1_co->cd();
    /*c1_co->SetLogy();*/
    c1_co->SetLeftMargin(0.12);

    // draw the frame_co using a histogram frame_co

    TH1F* frame_co = c1_co->DrawFrame(0,-0.15,8.5,1.00);
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

    Gen_v2_ks->Draw("P");
    Reco_recoRef_v2_ks->Draw("P");

    Gen_v2_la->Draw("P");
    Reco_recoRef_v2_la->Draw("P");

    TLegend* leg_co1 = new TLegend(0.15,0.65,0.27,0.75);
    leg_co1->SetFillColor(10);
    leg_co1->SetFillStyle(0);
    leg_co1->SetBorderSize(0);
    leg_co1->SetTextFont(42);
    leg_co1->SetTextSize(0.04);
    leg_co1->AddEntry(Gen_v2_ks, "K_{S}^{0}", "P");
    leg_co1->AddEntry(Gen_v2_la, "#Lambda/#bar{#Lambda}", "P");
    leg_co1->Draw();

    TLegend* leg_co2 = new TLegend(0.6,0.75,0.85,0.9);
    leg_co2->SetFillColor(10);
    leg_co2->SetFillStyle(0);
    leg_co2->SetBorderSize(0);
    leg_co2->SetTextFont(42);
    leg_co2->SetTextSize(0.03);
    leg_co2->AddEntry(Gen_v2_la, "Gen", "P");
    leg_co2->AddEntry(Reco_recoRef_v2_la, "Reco", "P");
    leg_co2->Draw();

    tex->SetTextFont(62);
    tex->SetTextSize(0.045);
    tex->DrawLatex(0.15,0.8,"pPb EPOS");
    tex->SetTextSize(0.04);
    tex->SetTextFont(42);
    line->Draw("same");

    c1_co->Print("v2CombinedClosureSystematics.pdf");

    //Calculate Ratios Reco/Gen
    std::vector<double> v2RatioKsReco;
    std::vector<double> v2RatioLaReco;
    std::vector<double> v2RatioKsRecoMatch;
    std::vector<double> v2RatioLaRecoMatch;

    std::vector<double> v2RatioKsGen;
    std::vector<double> v2RatioLaGen;

    std::vector<double> v2RatioKsGenRecoref;
    std::vector<double> v2RatioLaGenRecoref;

    std::vector<double> v2ErrorRatioKsReco;
    std::vector<double> v2ErrorRatioLaReco;
    std::vector<double> v2ErrorRatioKsRecoMatch;
    std::vector<double> v2ErrorRatioLaRecoMatch;

    std::vector<double> v2ErrorRatioKsGen;
    std::vector<double> v2ErrorRatioLaGen;

    std::vector<double> v2ErrorRatioKsGenRecoref;
    std::vector<double> v2ErrorRatioLaGenRecoref;

    double *aReco_recoRef_v2_ksX = Reco_recoRef_v2_ks->GetX();
    double *aReco_recoRef_v2_ksY = Reco_recoRef_v2_ks->GetY();
    double *aReco_recoRef_v2_ksEY = Reco_recoRef_v2_ks->GetEY();


    double *aGen_v2_ksX = Gen_v2_ks->GetX();
    double *aGen_v2_ksY = Gen_v2_ks->GetY();
    double *aGen_v2_ksEY = Gen_v2_ks->GetEY();


    double *aReco_recoRef_v2_laX = Reco_recoRef_v2_la->GetX();
    double *aReco_recoRef_v2_laY = Reco_recoRef_v2_la->GetY();
    double *aReco_recoRef_v2_laEY = Reco_recoRef_v2_la->GetEY();


    double *aGen_v2_laX = Gen_v2_la->GetX();
    double *aGen_v2_laY = Gen_v2_la->GetY();
    double *aGen_v2_laEY = Gen_v2_la->GetEY();



    //Kshort
    for(int i=0; i<Gen_v2_ks->GetN(); i++)
    {
        v2RatioKsReco.push_back(aReco_recoRef_v2_ksY[i]/ksFit->Eval(aReco_recoRef_v2_ksX[i]));
        v2ErrorRatioKsReco.push_back(aReco_recoRef_v2_ksEY[i]/ksFit->Eval(aReco_recoRef_v2_ksX[i]));


        v2RatioKsGen.push_back(aGen_v2_ksY[i]/ksFit->Eval(aGen_v2_ksX[i]));
        v2ErrorRatioKsGen.push_back(aGen_v2_ksEY[i]/ksFit->Eval(aGen_v2_ksX[i]));


        cout << "Ratio " << i << ": " << v2RatioKsReco[i] << endl;
    }

    //Lambda
    for(int i=0; i<Gen_v2_la->GetN(); i++)
    {
        v2RatioLaReco.push_back(aReco_recoRef_v2_laY[i]/laFit->Eval(aReco_recoRef_v2_laX[i]));
        v2ErrorRatioLaReco.push_back(aReco_recoRef_v2_laEY[i]/laFit->Eval(aReco_recoRef_v2_laX[i]));


        v2RatioLaGen.push_back(aGen_v2_laY[i]/laFit->Eval(aGen_v2_laX[i]));
        v2ErrorRatioLaGen.push_back(aGen_v2_laEY[i]/laFit->Eval(aGen_v2_laX[i]));

    }

    TGraphErrors* RatioReco_v2_la = new TGraphErrors(la_npoints,aReco_recoRef_v2_laX,&v2RatioLaReco[0],0,&v2ErrorRatioLaReco[0]);
    TGraphErrors* RatioGen_v2_la = new TGraphErrors(la_npoints,aGen_v2_laX,&v2RatioLaGen[0],0,&v2ErrorRatioLaGen[0]);

    TGraphErrors* RatioReco_v2_ks = new TGraphErrors(ks_npoints,aReco_recoRef_v2_ksX,&v2RatioKsReco[0],0,&v2ErrorRatioKsReco[0]);
    TGraphErrors* RatioGen_v2_ks = new TGraphErrors(ks_npoints,aGen_v2_ksX,&v2RatioKsGen[0],0,&v2ErrorRatioKsGen[0]);

    RatioReco_v2_ks->SetMarkerColor(kRed);
    RatioReco_v2_ks->SetMarkerStyle(25);
    RatioReco_v2_ks->SetMarkerSize(1.3);
    RatioReco_v2_ks->SetLineColor(kRed);


    RatioGen_v2_ks->SetMarkerColor(kRed);
    RatioGen_v2_ks->SetMarkerStyle(20);
    RatioGen_v2_ks->SetMarkerSize(1.3);
    RatioGen_v2_ks->SetLineColor(kRed);

    RatioReco_v2_la->SetMarkerColor(kBlue-4);
    RatioReco_v2_la->SetMarkerStyle(25);
    RatioReco_v2_la->SetMarkerSize(1.3);
    RatioReco_v2_la->SetLineColor(kBlue-4);

    RatioGen_v2_la->SetMarkerColor(kBlue-4);
    RatioGen_v2_la->SetMarkerStyle(20);
    RatioGen_v2_la->SetMarkerSize(1.3);
    RatioGen_v2_la->SetLineColor(kBlue-4);

    TCanvas* c1_ratio = MakeCanvas("c1_ratio", "Plot_ratio");
    TH1F* frame_Ratio = c1_ratio->DrawFrame(0,0.8,8.5,1.2);
    //TH1F* frame_Ratio = c1_ratio->DrawFrame(0,0,8.5,2);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Ratio->GetXaxis()->CenterTitle(1);
    frame_Ratio->GetYaxis()->CenterTitle(1);
    frame_Ratio->GetXaxis()->SetTitleSize(0.05);
    frame_Ratio->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Ratio->GetYaxis()->SetTitle("#frac{Reco v^{sig}_{2}}{Gen v^{sig}_{2}}");
    frame_Ratio->GetYaxis()->SetTitleSize(0.05);
    frame_Ratio->SetTitleOffset(1.2,"Y");
    frame_Ratio->SetTitleOffset(1.2,"X");

    TLine* LineRatio_ks = new TLine(0,1,8.5,1);
    LineRatio_ks->Draw();

    TLine* LineRatio_min_ks = new TLine(0,0.95,8.5,0.95);
    LineRatio_min_ks->SetLineStyle(2);
    LineRatio_min_ks->Draw();

    TLine* LineRatio_max_ks = new TLine(0,1.05,8.5,1.05);
    LineRatio_max_ks->SetLineStyle(2);
    LineRatio_max_ks->Draw();

    RatioReco_v2_ks->Draw("P");
    RatioGen_v2_ks->Draw("P");
    RatioGen_v2_la->Draw("P");
    RatioReco_v2_la->Draw("P");

    TLegend* leg_ratio_species = new TLegend(0.80,0.75,0.90,0.85);
    leg_ratio_species->SetFillColor(10);
    leg_ratio_species->SetFillStyle(0);
    leg_ratio_species->SetBorderSize(0);
    leg_ratio_species->SetTextFont(42);
    leg_ratio_species->SetTextSize(0.03);
    leg_ratio_species->AddEntry(RatioGen_v2_ks, "K_{S}^{0}", "P");
    leg_ratio_species->AddEntry(RatioGen_v2_la, "#Lambda/#bar{#Lambda}", "P");
    leg_ratio_species->Draw();

    TGraphErrors* RatioGen_v2_ks_clone = (TGraphErrors*)RatioGen_v2_ks->Clone("RatioGen_v2_ks_clone");
    RatioGen_v2_ks_clone->SetMarkerColor(kBlack);
    TGraphErrors* RatioReco_v2_ks_clone = (TGraphErrors*)RatioReco_v2_ks->Clone("RatioReco_v2_ks_clone");
    RatioReco_v2_ks_clone->SetMarkerColor(kBlack);
    TGraphErrors* refpoint = (TGraphErrors*)RatioGen_v2_ks->Clone("ref");
    refpoint->SetMarkerStyle(1);

    TLegend* leg_ratio_label = new TLegend(0.64,0.75,0.74,0.85);
    leg_ratio_label->SetFillColor(10);
    leg_ratio_label->SetFillStyle(0);
    leg_ratio_label->SetBorderSize(0);
    leg_ratio_label->SetTextFont(42);
    leg_ratio_label->SetTextSize(0.03);
    leg_ratio_label->AddEntry(RatioGen_v2_ks_clone, "Gen", "P");
    leg_ratio_label->AddEntry(RatioReco_v2_ks_clone, "Reco", "P");
    leg_ratio_label->Draw();

    tex->SetTextSize(0.03);
    tex->DrawLatex(0.4,0.85,"pPb EPOS");
    tex->DrawLatex(0.4,0.80,"|#eta| > 1");
    tex->DrawLatex(0.4,0.75,"0.3 < p_{T}^{assoc} < 3.0 GeV");

    //c1_ratio->Print("v2RatioClosureSystematics.pdf");
    c1_ratio->Print("v2RatioClosureSystematics_Combined.pdf");
    //c1_ratio->Print("v2RatioClosureSystematics_CheckGen_RecoRef_with_Match_RecoRef.pdf");

     //Write points into rootfile
    //TFile out("V0ClosureSys.root","RECREATE");
    //Gen_v2_la->Write("GenLa");
    //Reco_v2_la->Write("RecoLa");
    //Gen_v2_ks->Write("GenKs");
    //Reco_v2_ks->Write("RecoKs");
    //RatioReco_v2_ks->Write("RatioRecoKs");
    //RatioReco_v2_la->Write("RatioRecoLa");

    TCanvas* c1_ratio_ks = MakeCanvas("ks_ratio", "Plot_ks_ratio");
    TH1F* frame_Ratio_ks = c1_ratio_ks->DrawFrame(0,0.8,8.5,1.2);
    //TH1F* frame_Ratio = c1_ratio->DrawFrame(0,0,8.5,2);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Ratio_ks->GetXaxis()->CenterTitle(1);
    frame_Ratio_ks->GetYaxis()->CenterTitle(1);
    frame_Ratio_ks->GetXaxis()->SetTitleSize(0.05);
    frame_Ratio_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Ratio_ks->GetYaxis()->SetTitle("#frac{Reco v^{sig}_{2}}{Gen v^{sig}_{2}}");
    frame_Ratio_ks->GetYaxis()->SetTitleSize(0.05);
    frame_Ratio_ks->SetTitleOffset(1.2,"Y");
    frame_Ratio_ks->SetTitleOffset(1.2,"X");

    RatioReco_v2_ks->Draw("P");
    RatioGen_v2_ks->Draw("P");

    TLegend* leg_ratio_label_ks = new TLegend(0.64,0.70,0.74,0.85);
    leg_ratio_label_ks->SetFillColor(10);
    leg_ratio_label_ks->SetFillStyle(0);
    leg_ratio_label_ks->SetBorderSize(0);
    leg_ratio_label_ks->SetTextFont(42);
    leg_ratio_label_ks->SetTextSize(0.03);
    leg_ratio_label_ks->AddEntry(RatioGen_v2_ks_clone, "Gen K_{S}^{0}", "P");
    leg_ratio_label_ks->AddEntry(RatioReco_v2_ks_clone, "Reco K_{S}^{0}", "P");
    leg_ratio_label_ks->Draw();

    TLine* line_ks = new TLine(0,1,8.5,1);
    line_ks->Draw();

    tex->SetTextSize(0.035);
    tex->DrawLatex(0.35,0.85,"K_{S}^{0} pPb EPOS");
    tex->DrawLatex(0.35,0.80,"|#eta| > 1");
    tex->DrawLatex(0.35,0.75,"0.3 < p_{T}^{assoc} < 3.0 GeV");

    c1_ratio_ks->Print("v2RatioClosureSystematics_Kshort.pdf");
    c1_ratio_ks->Print("v2RatioClosureSystematics_Kshort.png");

    TCanvas* c1_ratio_la = MakeCanvas("la_ratio", "Plot_la_ratio");
    TH1F* frame_Ratio_la = c1_ratio_la->DrawFrame(0,0.8,8.5,1.2);
    //TH1F* frame_Ratio = c1_ratio->DrawFrame(0,0,8.5,2);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_Ratio_la->GetXaxis()->CenterTitle(1);
    frame_Ratio_la->GetYaxis()->CenterTitle(1);
    frame_Ratio_la->GetXaxis()->SetTitleSize(0.05);
    frame_Ratio_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_Ratio_la->GetYaxis()->SetTitle("#frac{Reco v^{sig}_{2}}{Gen v^{sig}_{2}}");
    frame_Ratio_la->GetYaxis()->SetTitleSize(0.05);
    frame_Ratio_la->SetTitleOffset(1.2,"Y");
    frame_Ratio_la->SetTitleOffset(1.2,"X");

    RatioReco_v2_la->Draw("P");
    RatioGen_v2_la->Draw("P");

    TLegend* leg_ratio_label_la = new TLegend(0.64,0.70,0.74,0.85);
    leg_ratio_label_la->SetFillColor(10);
    leg_ratio_label_la->SetFillStyle(0);
    leg_ratio_label_la->SetBorderSize(0);
    leg_ratio_label_la->SetTextFont(42);
    leg_ratio_label_la->SetTextSize(0.03);
    leg_ratio_label_la->AddEntry(RatioGen_v2_ks_clone, "Gen", "P");
    leg_ratio_label_la->AddEntry(RatioReco_v2_ks_clone, "Reco", "P");
    leg_ratio_label_la->Draw();

    TLine* line_la = new TLine(0,1,8.5,1);
    line_la->Draw();

    tex->SetTextSize(0.035);
    tex->DrawLatex(0.35,0.85,"#Lambda / #bar{#Lambda} pPb EPOS");
    tex->DrawLatex(0.35,0.80,"|#eta| > 1");
    tex->DrawLatex(0.35,0.75,"0.3 < p_{T}^{assoc} < 3.0 GeV");

    c1_ratio_la->Print("v2RatioClosureSystematics_Lambda.pdf");
    c1_ratio_la->Print("v2RatioClosureSystematics_Lambda.png");

}

void RapSys_etaGap()
{
    TFile* f_RapSys_etaGap_2 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Fit/v2Fit/data/v2valuesRapidity_etaGap_2.root");
    TFile* f_RapSys_etaGap_1p65 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Fit/v2Fit/data/v2valuesRapidity_etaGap_1p65.root");
    TFile* f_RapSys_etaGap_1p35 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Fit/v2Fit/data/v2valuesRapidity_etaGap_1p35.root");
    TFile* f_RapSys_etaGap_1p05 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Fit/v2Fit/data/v2valuesRapidity_etaGap_1p05.root");
    TFile* f_RapSys_etaGap_0p75 = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Fit/v2Fit/data/v2valuesRapidity_etaGap_0p75.root");
    MITStyle();
    TCanvas* c1 = MakeCanvas("c1", "Plot_ks");
    TCanvas* c2 = MakeCanvas("c2", "Plot_la");
    TCanvas* c3 = MakeCanvas("c3", "Ratio_ks");
    TCanvas* c4 = MakeCanvas("c4", "Ratio_la");
    TCanvas* c5 = MakeCanvas("c5", "Plot_xi");
    TCanvas* c6 = MakeCanvas("c6", "Plot_om");
    TCanvas* c7 = MakeCanvas("c7", "Ratio_xi");
    TCanvas* c8 = MakeCanvas("c8", "Ratio_om");
    c1->cd();
    /*c1->SetLogy();*/
    c1->SetLeftMargin(0.12);

    TGraphErrors* ks8_v2_2 = 0;
    TGraphErrors* la8_v2_2 = 0;
    TGraphErrors* xi8_v2_2 = 0;
    TGraphErrors* om8_v2_2 = 0;
    TGraphErrors* ks8_v2_1p65 = 0;
    TGraphErrors* la8_v2_1p65 = 0;
    TGraphErrors* xi8_v2_1p65 = 0;
    TGraphErrors* om8_v2_1p65 = 0;
    TGraphErrors* ks8_v2_1p35 = 0;
    TGraphErrors* la8_v2_1p35 = 0;
    TGraphErrors* xi8_v2_1p35 = 0;
    TGraphErrors* om8_v2_1p35 = 0;
    TGraphErrors* ks8_v2_1p05 = 0;
    TGraphErrors* la8_v2_1p05 = 0;
    TGraphErrors* xi8_v2_1p05 = 0;
    TGraphErrors* om8_v2_1p05 = 0;
    TGraphErrors* ks8_v2_0p75 = 0;
    TGraphErrors* la8_v2_0p75 = 0;
    TGraphErrors* xi8_v2_0p75 = 0;
    TGraphErrors* om8_v2_0p75 = 0;

    // draw the frame using a histogram frame
    TH1F* frame_ks;
    TH1F* frame_la;
    TH1F* frame_xi;
    TH1F* frame_om;
    TH1F* frame_ratio_ks;
    TH1F* frame_ratio_la;
    TH1F* frame_ratio_xi;
    TH1F* frame_ratio_om;

    ks8_v2_2 = (TGraphErrors*)f_RapSys_etaGap_2->Get("v2kshort");
    la8_v2_2 = (TGraphErrors*)f_RapSys_etaGap_2->Get("v2lambda");
    xi8_v2_2 = (TGraphErrors*)f_RapSys_etaGap_2->Get("v2xi");
    om8_v2_2 = (TGraphErrors*)f_RapSys_etaGap_2->Get("v2omega");
    ks8_v2_1p65 = (TGraphErrors*)f_RapSys_etaGap_1p65->Get("v2kshort");
    la8_v2_1p65 = (TGraphErrors*)f_RapSys_etaGap_1p65->Get("v2lambda");
    xi8_v2_1p65 = (TGraphErrors*)f_RapSys_etaGap_1p65->Get("v2xi");
    om8_v2_1p65 = (TGraphErrors*)f_RapSys_etaGap_1p65->Get("v2omega");
    ks8_v2_1p35 = (TGraphErrors*)f_RapSys_etaGap_1p35->Get("v2kshort");
    la8_v2_1p35 = (TGraphErrors*)f_RapSys_etaGap_1p35->Get("v2lambda");
    xi8_v2_1p35 = (TGraphErrors*)f_RapSys_etaGap_1p35->Get("v2xi");
    om8_v2_1p35 = (TGraphErrors*)f_RapSys_etaGap_1p35->Get("v2omega");
    ks8_v2_1p05 = (TGraphErrors*)f_RapSys_etaGap_1p05->Get("v2kshort");
    la8_v2_1p05 = (TGraphErrors*)f_RapSys_etaGap_1p05->Get("v2lambda");
    xi8_v2_1p05 = (TGraphErrors*)f_RapSys_etaGap_1p05->Get("v2xi");
    om8_v2_1p05 = (TGraphErrors*)f_RapSys_etaGap_1p05->Get("v2omega");
    ks8_v2_0p75 = (TGraphErrors*)f_RapSys_etaGap_0p75->Get("v2kshort");
    la8_v2_0p75 = (TGraphErrors*)f_RapSys_etaGap_0p75->Get("v2lambda");
    xi8_v2_0p75 = (TGraphErrors*)f_RapSys_etaGap_0p75->Get("v2xi");
    om8_v2_0p75 = (TGraphErrors*)f_RapSys_etaGap_0p75->Get("v2omega");

    double* ks_x = ks8_v2_2->GetX();
    double* ks_2_y = ks8_v2_2->GetY();
    double* ks_2_Ey = ks8_v2_2->GetEY();
    double* ks_1p65_y = ks8_v2_1p65->GetY();
    double* ks_1p65_Ey = ks8_v2_1p65->GetEY();
    double* ks_1p35_y = ks8_v2_1p35->GetY();
    double* ks_1p35_Ey = ks8_v2_1p35->GetEY();
    double* ks_1p05_y = ks8_v2_1p05->GetY();
    double* ks_1p05_Ey = ks8_v2_1p05->GetEY();
    double* ks_0p75_y = ks8_v2_0p75->GetY();
    double* ks_0p75_Ey = ks8_v2_0p75->GetEY();

    double* la_x = la8_v2_2->GetX();
    double* la_2_y = la8_v2_2->GetY();
    double* la_2_Ey = la8_v2_2->GetEY();
    double* la_1p65_y = la8_v2_1p65->GetY();
    double* la_1p65_Ey = la8_v2_1p65->GetEY();
    double* la_1p35_y = la8_v2_1p35->GetY();
    double* la_1p35_Ey = la8_v2_1p35->GetEY();
    double* la_1p05_y = la8_v2_1p05->GetY();
    double* la_1p05_Ey = la8_v2_1p05->GetEY();
    double* la_0p75_y = la8_v2_0p75->GetY();
    double* la_0p75_Ey = la8_v2_0p75->GetEY();

    double* xi_x = xi8_v2_2->GetX();
    double* xi_2_y = xi8_v2_2->GetY();
    double* xi_2_Ey = xi8_v2_2->GetEY();
    double* xi_1p65_y = xi8_v2_1p65->GetY();
    double* xi_1p65_Ey = xi8_v2_1p65->GetEY();
    double* xi_1p35_y = xi8_v2_1p35->GetY();
    double* xi_1p35_Ey = xi8_v2_1p35->GetEY();
    double* xi_1p05_y = xi8_v2_1p05->GetY();
    double* xi_1p05_Ey = xi8_v2_1p05->GetEY();
    double* xi_0p75_y = xi8_v2_0p75->GetY();
    double* xi_0p75_Ey = xi8_v2_0p75->GetEY();

    double* om_x = om8_v2_2->GetX();
    double* om_2_y = om8_v2_2->GetY();
    double* om_2_Ey = om8_v2_2->GetEY();
    double* om_1p65_y = om8_v2_1p65->GetY();
    double* om_1p65_Ey = om8_v2_1p65->GetEY();
    double* om_1p35_y = om8_v2_1p35->GetY();
    double* om_1p35_Ey = om8_v2_1p35->GetEY();
    double* om_1p05_y = om8_v2_1p05->GetY();
    double* om_1p05_Ey = om8_v2_1p05->GetEY();
    double* om_0p75_y = om8_v2_0p75->GetY();
    double* om_0p75_Ey = om8_v2_0p75->GetEY();

    int Npoints_ks = 0;
    int Npoints_la = 0;
    int Npoints_xi = 0;
    int Npoints_om = 0;
    Npoints_ks = ks8_v2_2->GetN();
    Npoints_la = la8_v2_2->GetN();
    Npoints_xi = xi8_v2_2->GetN();
    Npoints_om = om8_v2_2->GetN();


    std::vector<double> om_x_redo;
    std::vector<double> om_2_y_redo;
    std::vector<double> om_1p65_y_redo;
    std::vector<double> om_1p35_y_redo;
    std::vector<double> om_1p05_y_redo;
    std::vector<double> om_0p75_y_redo;

    std::vector<double> om_2_Ey_redo;
    std::vector<double> om_1p65_Ey_redo;
    std::vector<double> om_1p35_Ey_redo;
    std::vector<double> om_1p05_Ey_redo;
    std::vector<double> om_0p75_Ey_redo;

    for(int i=1; i<Npoints_om; i++)
    {
        om_x_redo.push_back(om_x[i]);
        om_2_y_redo.push_back(om_2_y[i]);
        om_1p65_y_redo.push_back(om_1p65_y[i]);
        om_1p35_y_redo.push_back(om_1p35_y[i]);
        om_1p05_y_redo.push_back(om_1p05_y[i]);
        om_0p75_y_redo.push_back(om_0p75_y[i]);

        om_2_Ey_redo.push_back(om_2_Ey[i]);
        om_1p65_Ey_redo.push_back(om_1p65_Ey[i]);
        om_1p35_Ey_redo.push_back(om_1p35_Ey[i]);
        om_1p05_Ey_redo.push_back(om_1p05_Ey[i]);
        om_0p75_Ey_redo.push_back(om_0p75_Ey[i]);
    }

    int Npoints_om_redo = Npoints_om-1;

    TGraphErrors* om8_v2_2_redo = new TGraphErrors(Npoints_om_redo,&om_x_redo[0],&om_2_y_redo[0],0,&om_2_Ey_redo[0]);
    TGraphErrors* om8_v2_1p65_redo = new TGraphErrors(Npoints_om_redo,&om_x_redo[0],&om_1p65_y_redo[0],0,&om_1p65_Ey_redo[0]);
    TGraphErrors* om8_v2_1p35_redo = new TGraphErrors(Npoints_om_redo,&om_x_redo[0],&om_1p35_y_redo[0],0,&om_1p35_Ey_redo[0]);
    TGraphErrors* om8_v2_1p05_redo = new TGraphErrors(Npoints_om_redo,&om_x_redo[0],&om_1p05_y_redo[0],0,&om_1p05_Ey_redo[0]);
    TGraphErrors* om8_v2_0p75_redo = new TGraphErrors(Npoints_om_redo,&om_x_redo[0],&om_0p75_y_redo[0],0,&om_0p75_Ey_redo[0]);

    c1->cd();
    frame_ks = c1->DrawFrame(0,-0.01,9,0.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_ks->SetTitle("K_{S}^{0}");
    frame_ks->GetXaxis()->CenterTitle(1);
    frame_ks->GetYaxis()->CenterTitle(1);
    frame_ks->GetXaxis()->SetTitleSize(0.05);
    frame_ks->GetYaxis()->SetTitleSize(0.05);
    frame_ks->SetTitleOffset(1.1,"Y");
    frame_ks->SetTitleOffset(1.2,"X");
    frame_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_ks->GetYaxis()->SetTitle("v_{2}^{sig}");

    c2->cd();
    frame_la = c2->DrawFrame(0,-0.01,9,0.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_la->SetTitle("#Lambda/#bar{#Lambda}");
    frame_la->GetXaxis()->CenterTitle(1);
    frame_la->GetYaxis()->CenterTitle(1);
    frame_la->GetXaxis()->SetTitleSize(0.05);
    frame_la->GetYaxis()->SetTitleSize(0.05);
    frame_la->SetTitleOffset(1.1,"Y");
    frame_la->SetTitleOffset(1.2,"X");
    frame_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_la->GetYaxis()->SetTitle("v_{2}^{sig}");

    c3->cd();
    frame_ratio_ks = c3->DrawFrame(0,0.90,9,1.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_ratio_ks->SetTitle("K_{S}^{0}");
    frame_ratio_ks->GetXaxis()->CenterTitle(1);
    frame_ratio_ks->GetYaxis()->CenterTitle(1);
    frame_ratio_ks->GetXaxis()->SetTitleSize(0.05);
    frame_ratio_ks->GetYaxis()->SetTitleSize(0.05);
    frame_ratio_ks->SetTitleOffset(1.1,"Y");
    frame_ratio_ks->SetTitleOffset(1.2,"X");
    frame_ratio_ks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_ratio_ks->GetYaxis()->SetTitle("v_{2}^{|#Delta#eta|}/v_{2}^{|#Delta#eta|>2}");

    c4->cd();
    frame_ratio_la = c4->DrawFrame(0,0.90,9,1.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_ratio_la->SetTitle("#Lambda/#bar{#Lambda}");
    frame_ratio_la->GetXaxis()->CenterTitle(1);
    frame_ratio_la->GetYaxis()->CenterTitle(1);
    frame_ratio_la->GetXaxis()->SetTitleSize(0.05);
    frame_ratio_la->GetYaxis()->SetTitleSize(0.05);
    frame_ratio_la->SetTitleOffset(1.1,"Y");
    frame_ratio_la->SetTitleOffset(1.2,"X");
    frame_ratio_la->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_ratio_la->GetYaxis()->SetTitle("v_{2}^{|#Delta#eta|}/v_{2}^{|#Delta#eta|>2}");

    c5->cd();
    frame_xi = c5->DrawFrame(0,-0.01,9,0.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_xi->SetTitle("#Xi^{#pm}");
    frame_xi->GetXaxis()->CenterTitle(1);
    frame_xi->GetYaxis()->CenterTitle(1);
    frame_xi->GetXaxis()->SetTitleSize(0.05);
    frame_xi->GetYaxis()->SetTitleSize(0.05);
    frame_xi->SetTitleOffset(1.1,"Y");
    frame_xi->SetTitleOffset(1.2,"X");
    frame_xi->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_xi->GetYaxis()->SetTitle("v_{2}^{sig}");

    c6->cd();
    frame_om = c6->DrawFrame(0,-0.01,9,0.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_om->SetTitle("#Omega^{#pm}");
    frame_om->GetXaxis()->CenterTitle(1);
    frame_om->GetYaxis()->CenterTitle(1);
    frame_om->GetXaxis()->SetTitleSize(0.05);
    frame_om->GetYaxis()->SetTitleSize(0.05);
    frame_om->SetTitleOffset(1.1,"Y");
    frame_om->SetTitleOffset(1.2,"X");
    frame_om->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_om->GetYaxis()->SetTitle("v_{2}^{sig}");

    c7->cd();
    frame_ratio_xi = c7->DrawFrame(0,0.90,9,1.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_ratio_xi->SetTitle("#Xi^{#pm}");
    frame_ratio_xi->GetXaxis()->CenterTitle(1);
    frame_ratio_xi->GetYaxis()->CenterTitle(1);
    frame_ratio_xi->GetXaxis()->SetTitleSize(0.05);
    frame_ratio_xi->GetYaxis()->SetTitleSize(0.05);
    frame_ratio_xi->SetTitleOffset(1.1,"Y");
    frame_ratio_xi->SetTitleOffset(1.2,"X");
    frame_ratio_xi->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_ratio_xi->GetYaxis()->SetTitle("v_{2}^{|#Delta#eta|}/v_{2}^{|#Delta#eta|>2}");

    c8->cd();
    frame_ratio_om = c8->DrawFrame(0,0.90,9,1.8);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_ratio_om->SetTitle("#Omega^{#pm}");
    frame_ratio_om->GetXaxis()->CenterTitle(1);
    frame_ratio_om->GetYaxis()->CenterTitle(1);
    frame_ratio_om->GetXaxis()->SetTitleSize(0.05);
    frame_ratio_om->GetYaxis()->SetTitleSize(0.05);
    frame_ratio_om->SetTitleOffset(1.1,"Y");
    frame_ratio_om->SetTitleOffset(1.2,"X");
    frame_ratio_om->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_ratio_om->GetYaxis()->SetTitle("v_{2}^{|#Delta#eta|}/v_{2}^{|#Delta#eta|>2}");

    std::vector<double> ratio_ks_1p65;
    std::vector<double> ratio_ks_1p35;
    std::vector<double> ratio_ks_1p05;
    std::vector<double> ratio_ks_0p75;

    std::vector<double> ratio_la_1p65;
    std::vector<double> ratio_la_1p35;
    std::vector<double> ratio_la_1p05;
    std::vector<double> ratio_la_0p75;

    std::vector<double> ratio_xi_1p65;
    std::vector<double> ratio_xi_1p35;
    std::vector<double> ratio_xi_1p05;
    std::vector<double> ratio_xi_0p75;

    std::vector<double> ratio_om_1p65;
    std::vector<double> ratio_om_1p35;
    std::vector<double> ratio_om_1p05;
    std::vector<double> ratio_om_0p75;

    std::vector<double> ratioE_ks_1p65;
    std::vector<double> ratioE_ks_1p35;
    std::vector<double> ratioE_ks_1p05;
    std::vector<double> ratioE_ks_0p75;

    std::vector<double> ratioE_la_1p65;
    std::vector<double> ratioE_la_1p35;
    std::vector<double> ratioE_la_1p05;
    std::vector<double> ratioE_la_0p75;

    std::vector<double> ratioE_xi_1p65;
    std::vector<double> ratioE_xi_1p35;
    std::vector<double> ratioE_xi_1p05;
    std::vector<double> ratioE_xi_0p75;

    std::vector<double> ratioE_om_1p65;
    std::vector<double> ratioE_om_1p35;
    std::vector<double> ratioE_om_1p05;
    std::vector<double> ratioE_om_0p75;

    std::vector<double> xcor_ks;
    std::vector<double> xcor_la;
    std::vector<double> xcor_xi;
    std::vector<double> xcor_om;

    for(int i=0; i<Npoints_ks; i++)
    {
        ratio_ks_1p65.push_back(ks_1p65_y[i]/ks_2_y[i] );
        ratio_ks_1p35.push_back(ks_1p35_y[i]/ks_2_y[i] );
        ratio_ks_1p05.push_back(ks_1p05_y[i]/ks_2_y[i] );
        ratio_ks_0p75.push_back(ks_0p75_y[i]/ks_2_y[i] );

        ratioE_ks_1p65.push_back(ks_1p65_Ey[i]/ks_2_y[i] );
        ratioE_ks_1p35.push_back(ks_1p35_Ey[i]/ks_2_y[i] );
        ratioE_ks_1p05.push_back(ks_1p05_Ey[i]/ks_2_y[i] );
        ratioE_ks_0p75.push_back(ks_0p75_Ey[i]/ks_2_y[i] );
    }
    for(int i=0; i<Npoints_la; i++)
    {
        ratio_la_1p65.push_back(la_1p65_y[i]/la_2_y[i] );
        ratio_la_1p35.push_back(la_1p35_y[i]/la_2_y[i] );
        ratio_la_1p05.push_back(la_1p05_y[i]/la_2_y[i] );
        ratio_la_0p75.push_back(la_0p75_y[i]/la_2_y[i] );

        ratioE_la_1p65.push_back(la_1p65_Ey[i]/la_2_y[i] );
        ratioE_la_1p35.push_back(la_1p35_Ey[i]/la_2_y[i] );
        ratioE_la_1p05.push_back(la_1p05_Ey[i]/la_2_y[i] );
        ratioE_la_0p75.push_back(la_0p75_Ey[i]/la_2_y[i] );
    }
    for(int i=0; i<Npoints_xi; i++)
    {
        ratio_xi_1p65.push_back(xi_1p65_y[i]/xi_2_y[i] );
        ratio_xi_1p35.push_back(xi_1p35_y[i]/xi_2_y[i] );
        ratio_xi_1p05.push_back(xi_1p05_y[i]/xi_2_y[i] );
        ratio_xi_0p75.push_back(xi_0p75_y[i]/xi_2_y[i] );

        ratioE_xi_1p65.push_back(xi_1p65_Ey[i]/xi_2_y[i] );
        ratioE_xi_1p35.push_back(xi_1p35_Ey[i]/xi_2_y[i] );
        ratioE_xi_1p05.push_back(xi_1p05_Ey[i]/xi_2_y[i] );
        ratioE_xi_0p75.push_back(xi_0p75_Ey[i]/xi_2_y[i] );
    }
    for(int i=0; i<Npoints_om; i++)
    {
        ratio_om_1p65.push_back(om_1p65_y[i]/om_2_y[i] );
        ratio_om_1p35.push_back(om_1p35_y[i]/om_2_y[i] );
        ratio_om_1p05.push_back(om_1p05_y[i]/om_2_y[i] );
        ratio_om_0p75.push_back(om_0p75_y[i]/om_2_y[i] );

        ratioE_om_1p65.push_back(om_1p65_Ey[i]/om_2_y[i] );
        ratioE_om_1p35.push_back(om_1p35_Ey[i]/om_2_y[i] );
        ratioE_om_1p05.push_back(om_1p05_Ey[i]/om_2_y[i] );
        ratioE_om_0p75.push_back(om_0p75_Ey[i]/om_2_y[i] );
    }


    TGraphErrors* ks8_v2_1p65_ratio = new TGraphErrors(Npoints_ks,ks_x,&ratio_ks_1p65[0],0,&ratioE_ks_1p65[0]);
    TGraphErrors* ks8_v2_1p35_ratio = new TGraphErrors(Npoints_ks,ks_x,&ratio_ks_1p35[0],0,&ratioE_ks_1p35[0]);
    TGraphErrors* ks8_v2_1p05_ratio = new TGraphErrors(Npoints_ks,ks_x,&ratio_ks_1p05[0],0,&ratioE_ks_1p05[0]);
    TGraphErrors* ks8_v2_0p75_ratio = new TGraphErrors(Npoints_ks,ks_x,&ratio_ks_0p75[0],0,&ratioE_ks_0p75[0]);

    TGraphErrors* la8_v2_1p65_ratio = new TGraphErrors(Npoints_la,la_x,&ratio_la_1p65[0],0,&ratioE_la_1p65[0]);
    TGraphErrors* la8_v2_1p35_ratio = new TGraphErrors(Npoints_la,la_x,&ratio_la_1p35[0],0,&ratioE_la_1p35[0]);
    TGraphErrors* la8_v2_1p05_ratio = new TGraphErrors(Npoints_la,la_x,&ratio_la_1p05[0],0,&ratioE_la_1p05[0]);
    TGraphErrors* la8_v2_0p75_ratio = new TGraphErrors(Npoints_la,la_x,&ratio_la_0p75[0],0,&ratioE_la_0p75[0]);

    TGraphErrors* xi8_v2_1p65_ratio = new TGraphErrors(Npoints_xi,xi_x,&ratio_xi_1p65[0],0,&ratioE_xi_1p65[0]);
    TGraphErrors* xi8_v2_1p35_ratio = new TGraphErrors(Npoints_xi,xi_x,&ratio_xi_1p35[0],0,&ratioE_xi_1p35[0]);
    TGraphErrors* xi8_v2_1p05_ratio = new TGraphErrors(Npoints_xi,xi_x,&ratio_xi_1p05[0],0,&ratioE_xi_1p05[0]);
    TGraphErrors* xi8_v2_0p75_ratio = new TGraphErrors(Npoints_xi,xi_x,&ratio_xi_0p75[0],0,&ratioE_xi_0p75[0]);

    TGraphErrors* om8_v2_1p65_ratio = new TGraphErrors(Npoints_om,om_x,&ratio_om_1p65[0],0,&ratioE_om_1p65[0]);
    TGraphErrors* om8_v2_1p35_ratio = new TGraphErrors(Npoints_om,om_x,&ratio_om_1p35[0],0,&ratioE_om_1p35[0]);
    TGraphErrors* om8_v2_1p05_ratio = new TGraphErrors(Npoints_om,om_x,&ratio_om_1p05[0],0,&ratioE_om_1p05[0]);
    TGraphErrors* om8_v2_0p75_ratio = new TGraphErrors(Npoints_om,om_x,&ratio_om_0p75[0],0,&ratioE_om_0p75[0]);

    ks8_v2_2->SetMarkerColor(kBlack);
    ks8_v2_2->SetMarkerStyle(20);
    ks8_v2_2->SetMarkerSize(1.5);
    ks8_v2_2->SetLineColor(kBlack);

    ks8_v2_1p65->SetMarkerColor(kRed+3);
    ks8_v2_1p65->SetMarkerStyle(24);
    ks8_v2_1p65->SetMarkerSize(1.5);
    ks8_v2_1p65->SetLineColor(kRed+3);

    ks8_v2_1p35->SetMarkerColor(kBlue);
    ks8_v2_1p35->SetMarkerStyle(24);
    ks8_v2_1p35->SetMarkerSize(1.5);
    ks8_v2_1p35->SetLineColor(kBlue);

    ks8_v2_1p05->SetMarkerColor(kGreen);
    ks8_v2_1p05->SetMarkerStyle(24);
    ks8_v2_1p05->SetMarkerSize(1.5);
    ks8_v2_1p05->SetLineColor(kGreen);

    ks8_v2_0p75->SetMarkerColor(kMagenta);
    ks8_v2_0p75->SetMarkerStyle(24);
    ks8_v2_0p75->SetMarkerSize(1.5);
    ks8_v2_0p75->SetLineColor(kMagenta);

    ks8_v2_1p65_ratio->SetMarkerColor(kRed);
    ks8_v2_1p65_ratio->SetMarkerStyle(24);
    ks8_v2_1p65_ratio->SetMarkerSize(1.5);
    ks8_v2_1p65_ratio->SetLineColor(kRed);

    ks8_v2_1p35_ratio->SetMarkerColor(kBlue);
    ks8_v2_1p35_ratio->SetMarkerStyle(24);
    ks8_v2_1p35_ratio->SetMarkerSize(1.5);
    ks8_v2_1p35_ratio->SetLineColor(kBlue);

    ks8_v2_1p05_ratio->SetMarkerColor(kGreen-2);
    ks8_v2_1p05_ratio->SetMarkerStyle(24);
    ks8_v2_1p05_ratio->SetMarkerSize(1.5);
    ks8_v2_1p05_ratio->SetLineColor(kGreen-2);

    ks8_v2_0p75_ratio->SetMarkerColor(kMagenta);
    ks8_v2_0p75_ratio->SetMarkerStyle(24);
    ks8_v2_0p75_ratio->SetMarkerSize(1.5);
    ks8_v2_0p75_ratio->SetLineColor(kMagenta);

    la8_v2_2->SetMarkerColor(kBlack);
    la8_v2_2->SetMarkerStyle(22);
    la8_v2_2->SetMarkerSize(1.5);
    la8_v2_2->SetLineColor(kBlack);

    la8_v2_1p65->SetMarkerColor(kRed+3);
    la8_v2_1p65->SetMarkerStyle(26);
    la8_v2_1p65->SetMarkerSize(1.5);
    la8_v2_1p65->SetLineColor(kRed+3);

    la8_v2_1p35->SetMarkerColor(kBlue);
    la8_v2_1p35->SetMarkerStyle(26);
    la8_v2_1p35->SetMarkerSize(1.5);
    la8_v2_1p35->SetLineColor(kBlue);

    la8_v2_1p05->SetMarkerColor(kGreen-2);
    la8_v2_1p05->SetMarkerStyle(26);
    la8_v2_1p05->SetMarkerSize(1.5);
    la8_v2_1p05->SetLineColor(kGreen-2);

    la8_v2_0p75->SetMarkerColor(kMagenta);
    la8_v2_0p75->SetMarkerStyle(26);
    la8_v2_0p75->SetMarkerSize(1.5);
    la8_v2_0p75->SetLineColor(kMagenta);

    la8_v2_1p65_ratio->SetMarkerColor(kRed);
    la8_v2_1p65_ratio->SetMarkerStyle(26);
    la8_v2_1p65_ratio->SetMarkerSize(1.5);
    la8_v2_1p65_ratio->SetLineColor(kRed);

    la8_v2_1p35_ratio->SetMarkerColor(kBlue);
    la8_v2_1p35_ratio->SetMarkerStyle(26);
    la8_v2_1p35_ratio->SetMarkerSize(1.5);
    la8_v2_1p35_ratio->SetLineColor(kBlue);

    la8_v2_1p05_ratio->SetMarkerColor(kGreen-2);
    la8_v2_1p05_ratio->SetMarkerStyle(26);
    la8_v2_1p05_ratio->SetMarkerSize(1.5);
    la8_v2_1p05_ratio->SetLineColor(kGreen-2);

    la8_v2_0p75_ratio->SetMarkerColor(kMagenta);
    la8_v2_0p75_ratio->SetMarkerStyle(26);
    la8_v2_0p75_ratio->SetMarkerSize(1.5);
    la8_v2_0p75_ratio->SetLineColor(kMagenta);

    xi8_v2_2->SetMarkerColor(kBlack);
    xi8_v2_2->SetMarkerStyle(20);
    xi8_v2_2->SetMarkerSize(1.5);
    xi8_v2_2->SetLineColor(kBlack);

    xi8_v2_1p65->SetMarkerColor(kRed+3);
    xi8_v2_1p65->SetMarkerStyle(24);
    xi8_v2_1p65->SetMarkerSize(1.5);
    xi8_v2_1p65->SetLineColor(kRed+3);

    xi8_v2_1p35->SetMarkerColor(kBlue);
    xi8_v2_1p35->SetMarkerStyle(24);
    xi8_v2_1p35->SetMarkerSize(1.5);
    xi8_v2_1p35->SetLineColor(kBlue);

    xi8_v2_1p05->SetMarkerColor(kGreen-2);
    xi8_v2_1p05->SetMarkerStyle(24);
    xi8_v2_1p05->SetMarkerSize(1.5);
    xi8_v2_1p05->SetLineColor(kGreen-2);

    xi8_v2_0p75->SetMarkerColor(kMagenta);
    xi8_v2_0p75->SetMarkerStyle(24);
    xi8_v2_0p75->SetMarkerSize(1.5);
    xi8_v2_0p75->SetLineColor(kMagenta);

    xi8_v2_1p65_ratio->SetMarkerColor(kRed);
    xi8_v2_1p65_ratio->SetMarkerStyle(24);
    xi8_v2_1p65_ratio->SetMarkerSize(1.5);
    xi8_v2_1p65_ratio->SetLineColor(kRed);

    xi8_v2_1p35_ratio->SetMarkerColor(kBlue);
    xi8_v2_1p35_ratio->SetMarkerStyle(24);
    xi8_v2_1p35_ratio->SetMarkerSize(1.5);
    xi8_v2_1p35_ratio->SetLineColor(kBlue);

    xi8_v2_1p05_ratio->SetMarkerColor(kGreen-2);
    xi8_v2_1p05_ratio->SetMarkerStyle(24);
    xi8_v2_1p05_ratio->SetMarkerSize(1.5);
    xi8_v2_1p05_ratio->SetLineColor(kGreen-2);

    xi8_v2_0p75_ratio->SetMarkerColor(kMagenta);
    xi8_v2_0p75_ratio->SetMarkerStyle(24);
    xi8_v2_0p75_ratio->SetMarkerSize(1.5);
    xi8_v2_0p75_ratio->SetLineColor(kMagenta);

    om8_v2_2_redo->SetMarkerColor(kBlack);
    om8_v2_2_redo->SetMarkerStyle(20);
    om8_v2_2_redo->SetMarkerSize(1.5);
    om8_v2_2_redo->SetLineColor(kBlack);

    om8_v2_1p65_redo->SetMarkerColor(kRed+3);
    om8_v2_1p65_redo->SetMarkerStyle(24);
    om8_v2_1p65_redo->SetMarkerSize(1.5);
    om8_v2_1p65_redo->SetLineColor(kRed+3);

    om8_v2_1p35_redo->SetMarkerColor(kBlue);
    om8_v2_1p35_redo->SetMarkerStyle(24);
    om8_v2_1p35_redo->SetMarkerSize(1.5);
    om8_v2_1p35_redo->SetLineColor(kBlue);

    om8_v2_1p05_redo->SetMarkerColor(kGreen-2);
    om8_v2_1p05_redo->SetMarkerStyle(24);
    om8_v2_1p05_redo->SetMarkerSize(1.5);
    om8_v2_1p05_redo->SetLineColor(kGreen-2);

    om8_v2_0p75_redo->SetMarkerColor(kMagenta);
    om8_v2_0p75_redo->SetMarkerStyle(24);
    om8_v2_0p75_redo->SetMarkerSize(1.5);
    om8_v2_0p75_redo->SetLineColor(kMagenta);

    om8_v2_1p65_ratio->SetMarkerColor(kRed);
    om8_v2_1p65_ratio->SetMarkerStyle(24);
    om8_v2_1p65_ratio->SetMarkerSize(1.5);
    om8_v2_1p65_ratio->SetLineColor(kRed);

    om8_v2_1p35_ratio->SetMarkerColor(kBlue);
    om8_v2_1p35_ratio->SetMarkerStyle(24);
    om8_v2_1p35_ratio->SetMarkerSize(1.5);
    om8_v2_1p35_ratio->SetLineColor(kBlue);

    om8_v2_1p05_ratio->SetMarkerColor(kGreen-2);
    om8_v2_1p05_ratio->SetMarkerStyle(24);
    om8_v2_1p05_ratio->SetMarkerSize(1.5);
    om8_v2_1p05_ratio->SetLineColor(kGreen-2);

    om8_v2_0p75_ratio->SetMarkerColor(kMagenta);
    om8_v2_0p75_ratio->SetMarkerStyle(24);
    om8_v2_0p75_ratio->SetMarkerSize(1.5);
    om8_v2_0p75_ratio->SetLineColor(kMagenta);

    TLegend* leg_ks = new TLegend(0.25,0.55,0.4,0.85);
    leg_ks->SetFillColor(10);
    leg_ks->SetFillStyle(0);
    leg_ks->SetBorderSize(0);
    leg_ks->SetTextFont(42);
    leg_ks->SetTextSize(0.05);
    leg_ks->AddEntry(ks8_v2_2, "|#Delta#eta| > 2.0", "P");
    leg_ks->AddEntry(ks8_v2_1p65, "|#Delta#eta| > 1.65", "P");
    leg_ks->AddEntry(ks8_v2_1p35, "|#Delta#eta| > 1.35", "P");
    leg_ks->AddEntry(ks8_v2_1p05, "|#Delta#eta| > 1.05", "P");
    leg_ks->AddEntry(ks8_v2_0p75, "|#Delta#eta| > 0.75", "P");
    //leg_ks->AddEntry(xi8_v2, "#Xi^{+}/ #Xi^{-}", "P");

    TLegend* leg_la = new TLegend(0.25,0.55,0.4,0.85);
    leg_la->SetFillColor(10);
    leg_la->SetFillStyle(0);
    leg_la->SetBorderSize(0);
    leg_la->SetTextFont(42);
    leg_la->SetTextSize(0.05);
    leg_la->AddEntry(la8_v2_2, "|#Delta#eta| > 2.0", "P");
    leg_la->AddEntry(la8_v2_1p65, "|#Delta#eta| > 1.65", "P");
    leg_la->AddEntry(la8_v2_1p35, "|#Delta#eta| > 1.35", "P");
    leg_la->AddEntry(la8_v2_1p05, "|#Delta#eta| > 1.05", "P");
    leg_la->AddEntry(la8_v2_0p75, "|#Delta#eta| > 0.75", "P");
    //leg_la->AddEntry(xi8_v2, "#Xi^{+}/ #Xi^{-}", "P");

    TLegend* leg_xi = new TLegend(0.25,0.55,0.40,0.85);
    leg_xi->SetFillColor(10);
    leg_xi->SetFillStyle(0);
    leg_xi->SetBorderSize(0);
    leg_xi->SetTextFont(42);
    leg_xi->SetTextSize(0.05);
    leg_xi->AddEntry(xi8_v2_2, "|#Delta#eta| > 2.0", "P");
    leg_xi->AddEntry(xi8_v2_1p65, "|#Delta#eta| > 1.65", "P");
    leg_xi->AddEntry(xi8_v2_1p35, "|#Delta#eta| > 1.35", "P");
    leg_xi->AddEntry(xi8_v2_1p05, "|#Delta#eta| > 1.05", "P");
    leg_xi->AddEntry(xi8_v2_0p75, "|#Delta#eta| > 0.75", "P");
    //leg_xi->AddEntry(xi8_v2, "#Xi^{+}/ #Xi^{-}", "P");

    TLegend* leg_om = new TLegend(0.25,0.55,0.40,0.85);
    leg_om->SetFillColor(10);
    leg_om->SetFillStyle(0);
    leg_om->SetBorderSize(0);
    leg_om->SetTextFont(42);
    leg_om->SetTextSize(0.05);
    leg_om->AddEntry(om8_v2_2_redo, "|#Delta#eta| > 2.0", "P");
    leg_om->AddEntry(om8_v2_1p65_redo, "|#Delta#eta| > 1.65", "P");
    leg_om->AddEntry(om8_v2_1p35_redo, "|#Delta#eta| > 1.35", "P");
    leg_om->AddEntry(om8_v2_1p05_redo, "|#Delta#eta| > 1.05", "P");
    leg_om->AddEntry(om8_v2_0p75_redo, "|#Delta#eta| > 0.75", "P");
    //leg_om->AddEntry(xi8_v2, "#Xi^{+}/ #Xi^{-}", "P");

    TLegend* leg_ks_ratio = new TLegend(0.72,0.59,0.87,0.89);
    leg_ks_ratio->SetFillColor(10);
    leg_ks_ratio->SetFillStyle(0);
    leg_ks_ratio->SetBorderSize(0);
    leg_ks_ratio->SetTextFont(42);
    leg_ks_ratio->SetTextSize(0.05);
    leg_ks_ratio->AddEntry(ks8_v2_1p65_ratio, "|#Delta#eta| > 1.65", "P");
    leg_ks_ratio->AddEntry(ks8_v2_1p35_ratio, "|#Delta#eta| > 1.35", "P");
    leg_ks_ratio->AddEntry(ks8_v2_1p05_ratio, "|#Delta#eta| > 1.05", "P");
    leg_ks_ratio->AddEntry(ks8_v2_0p75_ratio, "|#Delta#eta| > 0.75", "P");
    //leg_ks->AddEntry(xi8_v2, "#Xi^{+}/ #Xi^{-}", "P");

    TLegend* leg_la_ratio = new TLegend(0.72,0.59,0.87,0.89);
    leg_la_ratio->SetFillColor(10);
    leg_la_ratio->SetFillStyle(0);
    leg_la_ratio->SetBorderSize(0);
    leg_la_ratio->SetTextFont(42);
    leg_la_ratio->SetTextSize(0.05);
    leg_la_ratio->AddEntry(la8_v2_1p65_ratio, "|#Delta#eta| > 1.65", "P");
    leg_la_ratio->AddEntry(la8_v2_1p35_ratio, "|#Delta#eta| > 1.35", "P");
    leg_la_ratio->AddEntry(la8_v2_1p05_ratio, "|#Delta#eta| > 1.05", "P");
    leg_la_ratio->AddEntry(la8_v2_0p75_ratio, "|#Delta#eta| > 0.75", "P");

    TLegend* leg_xi_ratio = new TLegend(0.72,0.59,0.87,0.89);
    leg_xi_ratio->SetFillColor(10);
    leg_xi_ratio->SetFillStyle(0);
    leg_xi_ratio->SetBorderSize(0);
    leg_xi_ratio->SetTextFont(42);
    leg_xi_ratio->SetTextSize(0.05);
    leg_xi_ratio->AddEntry(xi8_v2_1p65_ratio, "|#Delta#eta| > 1.65", "P");
    leg_xi_ratio->AddEntry(xi8_v2_1p35_ratio, "|#Delta#eta| > 1.35", "P");
    leg_xi_ratio->AddEntry(xi8_v2_1p05_ratio, "|#Delta#eta| > 1.05", "P");
    leg_xi_ratio->AddEntry(xi8_v2_0p75_ratio, "|#Delta#eta| > 0.75", "P");
    //leg_xi->AddEntry(xi8_v2, "#Xi^{+}/ #Xi^{-}", "P");

    TLegend* leg_om_ratio = new TLegend(0.72,0.59,0.87,0.89);
    leg_om_ratio->SetFillColor(10);
    leg_om_ratio->SetFillStyle(0);
    leg_om_ratio->SetBorderSize(0);
    leg_om_ratio->SetTextFont(42);
    leg_om_ratio->SetTextSize(0.05);
    leg_om_ratio->AddEntry(om8_v2_1p65_ratio, "|#Delta#eta| > 1.65", "P");
    leg_om_ratio->AddEntry(om8_v2_1p35_ratio, "|#Delta#eta| > 1.35", "P");
    leg_om_ratio->AddEntry(om8_v2_1p05_ratio, "|#Delta#eta| > 1.05", "P");
    leg_om_ratio->AddEntry(om8_v2_0p75_ratio, "|#Delta#eta| > 0.75", "P");
    //leg_om->AddEntry(xi8_v2, "#Xi^{+}/ #Xi^{-}", "P");

    TLine* line = new TLine(0,1,9,1);
    line->SetLineStyle(2);

    c1->cd();
    ks8_v2_2->Draw("P");
    ks8_v2_1p65->Draw("P");
    ks8_v2_1p35->Draw("P");
    ks8_v2_1p05->Draw("P");
    ks8_v2_0p75->Draw("P");
    leg_ks->Draw("same");

    c2->cd();
    la8_v2_2->Draw("P");
    la8_v2_1p65->Draw("P");
    la8_v2_1p35->Draw("P");
    la8_v2_1p05->Draw("P");
    la8_v2_0p75->Draw("P");
    leg_la->Draw("same");
    //xi8_v2->Draw("P");

    c3->cd();
    ks8_v2_1p65_ratio->Draw("P");
    ks8_v2_1p35_ratio->Draw("P");
    ks8_v2_1p05_ratio->Draw("P");
    ks8_v2_0p75_ratio->Draw("P");
    leg_ks_ratio->Draw("same");
    line->Draw("same");

    c4->cd();
    la8_v2_1p65_ratio->Draw("P");
    la8_v2_1p35_ratio->Draw("P");
    la8_v2_1p05_ratio->Draw("P");
    la8_v2_0p75_ratio->Draw("P");
    leg_la_ratio->Draw("same");
    line->Draw("same");

    c5->cd();
    xi8_v2_2->Draw("P");
    xi8_v2_1p65->Draw("P");
    xi8_v2_1p35->Draw("P");
    xi8_v2_1p05->Draw("P");
    xi8_v2_0p75->Draw("P");
    leg_xi->Draw("same");

    c6->cd();
    om8_v2_2_redo->Draw("P");
    om8_v2_1p65_redo->Draw("P");
    om8_v2_1p35_redo->Draw("P");
    om8_v2_1p05_redo->Draw("P");
    om8_v2_0p75_redo->Draw("P");
    leg_om->Draw("same");

    c7->cd();
    xi8_v2_1p65_ratio->Draw("P");
    xi8_v2_1p35_ratio->Draw("P");
    xi8_v2_1p05_ratio->Draw("P");
    xi8_v2_0p75_ratio->Draw("P");
    leg_xi_ratio->Draw("same");
    line->Draw("same");

    c8->cd();
    om8_v2_1p65_ratio->Draw("P");
    om8_v2_1p35_ratio->Draw("P");
    om8_v2_1p05_ratio->Draw("P");
    om8_v2_0p75_ratio->Draw("P");
    leg_om_ratio->Draw("same");
    line->Draw("same");

    c1->Print("v2SigEtaGap_ks.pdf");
    c2->Print("v2SigEtaGap_la.pdf");
    c3->Print("v2SigEtaGap_ratio_ks.pdf");
    c4->Print("v2SigEtaGap_ratio_la.pdf");
    c5->Print("v2SigEtaGap_xi.pdf");
    c6->Print("v2SigEtaGap_om.pdf");
    c7->Print("v2SigEtaGap_ratio_xi.pdf");
    c8->Print("v2SigEtaGap_ratio_om.pdf");

    TCanvas* c9 = new TCanvas("Combined","Combined",1200,1000);
    c9->Divide(2,2);

    TH1F* frame_CoRatioks;
    TH1F* frame_CoRatiola;
    TH1F* frame_CoRatioxi;
    TH1F* frame_CoRatioom;


    frame_CoRatioks = c9->cd(1)->DrawFrame(0,0.90,9,1.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_CoRatioks->SetTitle("K_{S}^{0}");
    frame_CoRatioks->GetXaxis()->CenterTitle(1);
    frame_CoRatioks->GetYaxis()->CenterTitle(1);
    frame_CoRatioks->GetXaxis()->SetTitleSize(0.05);
    frame_CoRatioks->GetYaxis()->SetTitleSize(0.05);
    frame_CoRatioks->SetTitleOffset(1.1,"Y");
    frame_CoRatioks->SetTitleOffset(1.2,"X");
    frame_CoRatioks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_CoRatioks->GetYaxis()->SetTitle("v_{2}^{|#Delta#eta|}/v_{2}^{|#Delta#eta|>2}");

    ks8_v2_1p65_ratio->Draw("P");
    ks8_v2_1p35_ratio->Draw("P");
    ks8_v2_1p05_ratio->Draw("P");
    ks8_v2_0p75_ratio->Draw("P");
    leg_ks_ratio->Draw("same");
    line->Draw("same");

    frame_CoRatiola = c9->cd(2)->DrawFrame(0,0.90,9,1.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_CoRatiola->SetTitle("#Lambda/#bar{#Lambda}");
    frame_CoRatiola->GetXaxis()->CenterTitle(1);
    frame_CoRatiola->GetYaxis()->CenterTitle(1);
    frame_CoRatiola->GetXaxis()->SetTitleSize(0.05);
    frame_CoRatiola->GetYaxis()->SetTitleSize(0.05);
    frame_CoRatiola->SetTitleOffset(1.1,"Y");
    frame_CoRatiola->SetTitleOffset(1.2,"X");
    frame_CoRatiola->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_CoRatiola->GetYaxis()->SetTitle("v_{2}^{|#Delta#eta|}/v_{2}^{|#Delta#eta|>2}");
    la8_v2_1p65_ratio->Draw("P");
    la8_v2_1p35_ratio->Draw("P");
    la8_v2_1p05_ratio->Draw("P");
    la8_v2_0p75_ratio->Draw("P");
    leg_la_ratio->Draw("same");
    line->Draw("same");

    frame_CoRatioxi = c9->cd(3)->DrawFrame(0,0.90,9,1.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_CoRatioxi->SetTitle("#Xi^{#pm}");
    frame_CoRatioxi->GetXaxis()->CenterTitle(1);
    frame_CoRatioxi->GetYaxis()->CenterTitle(1);
    frame_CoRatioxi->GetXaxis()->SetTitleSize(0.05);
    frame_CoRatioxi->GetYaxis()->SetTitleSize(0.05);
    frame_CoRatioxi->SetTitleOffset(1.1,"Y");
    frame_CoRatioxi->SetTitleOffset(1.2,"X");
    frame_CoRatioxi->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_CoRatioxi->GetYaxis()->SetTitle("v_{2}^{|#Delta#eta|}/v_{2}^{|#Delta#eta|>2}");
    xi8_v2_1p65_ratio->Draw("P");
    xi8_v2_1p35_ratio->Draw("P");
    xi8_v2_1p05_ratio->Draw("P");
    xi8_v2_0p75_ratio->Draw("P");
    leg_xi_ratio->Draw("same");
    line->Draw("same");

    frame_CoRatioom = c9->cd(4)->DrawFrame(0,0.90,9,2.0);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_CoRatioom->SetTitle("#Omega^{#pm}");
    frame_CoRatioom->GetXaxis()->CenterTitle(1);
    frame_CoRatioom->GetYaxis()->CenterTitle(1);
    frame_CoRatioom->GetXaxis()->SetTitleSize(0.05);
    frame_CoRatioom->GetYaxis()->SetTitleSize(0.05);
    frame_CoRatioom->SetTitleOffset(1.1,"Y");
    frame_CoRatioom->SetTitleOffset(1.2,"X");
    frame_CoRatioom->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_CoRatioom->GetYaxis()->SetTitle("v_{2}^{|#Delta#eta|}/v_{2}^{|#Delta#eta|>2}");
    om8_v2_1p65_ratio->Draw("P");
    om8_v2_1p35_ratio->Draw("P");
    om8_v2_1p05_ratio->Draw("P");
    om8_v2_0p75_ratio->Draw("P");
    leg_om_ratio->Draw("same");
    line->Draw("same");

    c9->Print("v2SigEtaGap_ratio_Combined.pdf");
    c9->Print("v2SigEtaGap_ratio_Combined.png");

    TH1F* frame_coks;
    TH1F* frame_cola;
    TH1F* frame_coxi;
    TH1F* frame_coom;

    TCanvas* c10 = new TCanvas("C10","C10",1200,1000);
    c10->Divide(2,2);

    frame_coks = c10->cd(1)->DrawFrame(0,-0.01,9,0.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_coks->SetTitle("K_{S}^{0}");
    frame_coks->GetXaxis()->CenterTitle(1);
    frame_coks->GetYaxis()->CenterTitle(1);
    frame_coks->GetXaxis()->SetTitleSize(0.05);
    frame_coks->GetYaxis()->SetTitleSize(0.05);
    frame_coks->SetTitleOffset(1.1,"Y");
    frame_coks->SetTitleOffset(1.2,"X");
    frame_coks->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_coks->GetYaxis()->SetTitle("v_{2}^{sig}");
    ks8_v2_2->Draw("P");
    ks8_v2_1p65->Draw("P");
    ks8_v2_1p35->Draw("P");
    ks8_v2_1p05->Draw("P");
    ks8_v2_0p75->Draw("P");
    leg_ks->Draw("same");

    frame_cola = c10->cd(2)->DrawFrame(0,-0.01,9,0.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_cola->SetTitle("#Lambda / #bar{#Lambda}");
    frame_cola->GetXaxis()->CenterTitle(1);
    frame_cola->GetYaxis()->CenterTitle(1);
    frame_cola->GetXaxis()->SetTitleSize(0.05);
    frame_cola->GetYaxis()->SetTitleSize(0.05);
    frame_cola->SetTitleOffset(1.1,"Y");
    frame_cola->SetTitleOffset(1.2,"X");
    frame_cola->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_cola->GetYaxis()->SetTitle("v_{2}^{sig}");
    la8_v2_2->Draw("P");
    la8_v2_1p65->Draw("P");
    la8_v2_1p35->Draw("P");
    la8_v2_1p05->Draw("P");
    la8_v2_0p75->Draw("P");
    leg_la->Draw("same");

    frame_coxi = c10->cd(3)->DrawFrame(0,-0.01,9,0.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_coxi->SetTitle("#Xi^{#pm}");
    frame_coxi->GetXaxis()->CenterTitle(1);
    frame_coxi->GetYaxis()->CenterTitle(1);
    frame_coxi->GetXaxis()->SetTitleSize(0.05);
    frame_coxi->GetYaxis()->SetTitleSize(0.05);
    frame_coxi->SetTitleOffset(1.1,"Y");
    frame_coxi->SetTitleOffset(1.2,"X");
    frame_coxi->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_coxi->GetYaxis()->SetTitle("v_{2}^{sig}");
    xi8_v2_2->Draw("P");
    xi8_v2_1p65->Draw("P");
    xi8_v2_1p35->Draw("P");
    xi8_v2_1p05->Draw("P");
    xi8_v2_0p75->Draw("P");
    leg_xi->Draw("same");

    frame_coom = c10->cd(4)->DrawFrame(0,-0.01,9,0.5);
    gPad->SetTickx();
    gPad->SetTicky();
    frame_coom->SetTitle("#Omega^{#pm}");
    frame_coom->GetXaxis()->CenterTitle(1);
    frame_coom->GetYaxis()->CenterTitle(1);
    frame_coom->GetXaxis()->SetTitleSize(0.05);
    frame_coom->GetYaxis()->SetTitleSize(0.05);
    frame_coom->SetTitleOffset(1.1,"Y");
    frame_coom->SetTitleOffset(1.2,"X");
    frame_coom->GetXaxis()->SetTitle("p_{T} (GeV)");
    frame_coom->GetYaxis()->SetTitle("v_{2}^{sig}");
    om8_v2_2_redo->Draw("P");
    om8_v2_1p65_redo->Draw("P");
    om8_v2_1p35_redo->Draw("P");
    om8_v2_1p05_redo->Draw("P");
    om8_v2_0p75_redo->Draw("P");
    leg_om->Draw("same");

    c10->Print("v2SigRapidityCombined.pdf");
    c10->Print("v2SigRapidityCombined.png");

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

void CascadeRebin()
{
    TFile* f = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/PeripheralSubtractionMB.root");

    TH2D* MassPtXi = (TH2D*)f->Get("v0CasCorrelationRapidityPeriSub/MassPtXi");
    TH2D* MassPtOm = (TH2D*)f->Get("v0CasCorrelationRapidityPeriSub/MassPtOm");

    TCanvas* c_xi1 = new TCanvas("c_xi1","c_xi1",1200,900);
    TCanvas* c_xi2 = new TCanvas("c_xi2","c_xi2",1200,900);
    TCanvas* c_om1 = new TCanvas("c_om","c_om",1200,900);
    TCanvas* c_om2 = new TCanvas("c_om2","c_om2",1200,900);

    c_xi1->Divide(3,2);
    c_xi2->Divide(2,2);
    c_om1->Divide(3,2);
    c_om2->Divide(2,2);

    TH1D* hOm_1 [8];
    TH1D* hOm_2 [8];

    TH1D* hXi_1 [8];
    TH1D* hXi_2 [8];
    TLatex* tex = new TLatex();
    tex->SetTextFont(42);
    tex->SetTextSize(0.05);
    tex->SetNDC();

    std::vector<double> PtBinOm = {10,11,12,13,14,15,17};
    std::vector<double> PtBinOm2 = {10,18,23,30,72};

    hOm_1[0] = (TH1D*)MassPtOm->ProjectionX("om_1",PtBinOm[0]+1,PtBinOm[6]);
    hOm_1[1] = (TH1D*)MassPtOm->ProjectionX("om_2",PtBinOm[1]+1,PtBinOm[6]);
    hOm_1[2] = (TH1D*)MassPtOm->ProjectionX("om_3",PtBinOm[2]+1,PtBinOm[6]);
    hOm_1[3] = (TH1D*)MassPtOm->ProjectionX("om_4",PtBinOm[3]+1,PtBinOm[6]);
    hOm_1[4] = (TH1D*)MassPtOm->ProjectionX("om_5",PtBinOm[4]+1,PtBinOm[6]);
    hOm_1[5] = (TH1D*)MassPtOm->ProjectionX("om_6",PtBinOm[5]+1,PtBinOm[6]);

    hOm_2[0] = (TH1D*)MassPtOm->ProjectionX("om_7",PtBinOm2[0]+1,PtBinOm2[1]);
    hOm_2[1] = (TH1D*)MassPtOm->ProjectionX("om_8",PtBinOm2[1]+1,PtBinOm2[2]);
    hOm_2[2] = (TH1D*)MassPtOm->ProjectionX("om_9",PtBinOm2[2]+1,PtBinOm2[3]);
    hOm_2[3] = (TH1D*)MassPtOm->ProjectionX("om_10",PtBinOm2[3]+1,PtBinOm2[4]);

    std::vector<double> PtBinXi2 = {10,18,25,36,72};

    hXi_2[0] = (TH1D*)MassPtXi->ProjectionX("xi_7",PtBinXi2[0]+1,PtBinXi2[1]);
    hXi_2[1] = (TH1D*)MassPtXi->ProjectionX("xi_8",PtBinXi2[1]+1,PtBinXi2[2]);
    hXi_2[2] = (TH1D*)MassPtXi->ProjectionX("xi_9",PtBinXi2[2]+1,PtBinXi2[3]);
    hXi_2[3] = (TH1D*)MassPtXi->ProjectionX("xi_10",PtBinXi2[3]+1,PtBinXi2[4]);

    for(unsigned i=0; i<PtBinOm.size()-1; i++)
    {
        c_om1->cd(i+1);
        hOm_1[i]->SetMarkerStyle(20);
        hOm_1[i]->SetMarkerSize(0.5);
        hOm_1[i]->GetXaxis()->SetRangeUser(1.6,1.75);
        hOm_1[i]->Draw("P");
        //tex->DrawLatex(0.2,0.8,Form("%2.1f < Pt < %2.1f",(PtBinOm[i])/10,PtBinOm[i+1]/10));
        tex->DrawLatex(0.2,0.8,Form("%2.1f < Pt < %2.1f",(PtBinOm[i])/10,PtBinOm[6]/10));
    }

    //for(unsigned i=0; i<PtBinOm2.size()-1; i++)
    //{
        //c_om2->cd(i+1);
        //hOm_2[i]->SetMarkerStyle(20);
        //hOm_2[i]->SetMarkerSize(0.5);
        //hOm_2[i]->GetXaxis()->SetRangeUser(1.6,1.75);
        //hOm_2[i]->Draw("P");
        //tex->DrawLatex(0.2,0.8,Form("%2.1f < Pt < %2.1f",(PtBinOm[i])/10,PtBinOm[i+1]/10));
        //tex->DrawLatex(0.2,0.8,Form("%2.1f < Pt < %2.1f",(PtBinOm2[i])/10,PtBinOm2[i+1]/10));
    //}

    for(unsigned i=0; i<PtBinXi2.size()-1; i++)
    {
        c_xi2->cd(i+1);
        hXi_2[i]->SetMarkerStyle(20);
        hXi_2[i]->SetMarkerSize(0.5);
        hXi_2[i]->GetXaxis()->SetRangeUser(1.25,1.40);
        hXi_2[i]->Draw("P");
        //tex->DrawLatex(0.2,0.8,Form("%2.1f < Pt < %2.1f",(PtBinXi[i])/10,PtBinXi[i+1]/10));
        tex->DrawLatex(0.2,0.8,Form("%2.1f < Pt < %2.1f",(PtBinXi2[i])/10,PtBinXi2[i+1]/10));
    }

    ostringstream os;
    ostringstream osYield;
    for(unsigned i=0; i<PtBinOm2.size()-1; i++)
    {
        double s1 =0.003;
        double s2 = 0.003;
        RooRealVar x("x","mass",1.6,1.75);
        RooPlot* xframe_ = x.frame(150);
        RooRealVar mean("mean","mean",1.67,1.6,1.75);//Omega
        xframe_->GetYaxis()->SetTitle("Candidates / 0.001 GeV");
        RooDataHist data("data","dataset",x,hOm_2[i]);
        data.plotOn(xframe_,Name("data"));
        RooRealVar sigma1("sigma1","sigma1",s1,0.001,0.04);
        RooRealVar sigma2("sigma2","sigma2",s2,0.001,0.04);
        //RooRealVar sig1("sig1","signal1",100,-100,10000000);
        //RooRealVar sig2("sig2","signal2",100,-100,10000000);
        //RooRealVar qsig("qsig","qsig",50,0,1000000);
        RooRealVar sig1("sig1","signal1",35,-100,10000000);
        RooRealVar sig2("sig2","signal2",35,-100,10000000);
        RooRealVar qsig("qsig","qsig",100,0,1000000);
        RooRealVar alpha("alpha","alpha",0.5,0,2);
        RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
        RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
        RooRealVar ap("ap","ap",-0.1,-1,1);
        RooRealVar bp("bp","bp",-0.1,-1,1);
        RooRealVar cp("cp","cp",-0.1,-1,1);
        RooRealVar dp("dp","dp",-0.1,-1,1);
        RooChebychev background("background","background",x,RooArgList(ap,bp,cp,dp));
        RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,background),RooArgList(sig1,sig2,qsig));

        x.setRange("cut",1.625,1.75);

        RooFitResult* r_xi = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));

        double covQual = r_xi->covQual();
        double mean_xi = mean.getVal();

        double gaus1F_xi = sig1.getVal();
        double gaus2F_xi = sig2.getVal();
        double qsig_xi   = qsig.getVal();

        //set ranges for individual gaussian yield determination
        x.setRange("g1", mean.getVal() - 2*sigma1.getVal(), mean.getVal() + 2*sigma1.getVal());
        x.setRange("g2", mean.getVal() - 2*sigma2.getVal(), mean.getVal() + 2*sigma2.getVal());

        RooAbsReal* Intgaus1_yield_xi = gaus1.createIntegral(x,x,"g1");
        RooAbsReal* Intgaus2_yield_xi = gaus2.createIntegral(x,x,"g2");

        double gaus1_yield_xi = gaus1F_xi*Intgaus1_yield_xi->getVal();
        double gaus2_yield_xi = gaus2F_xi*Intgaus2_yield_xi->getVal();
        double gausTot_yield_xi = gaus1_yield_xi + gaus2_yield_xi;

        //cout << "Yield1: " << gaus1_yield_xi << endl;
        //cout << "Yield2: " << gaus2_yield_xi << endl;

        double rms_gaus1_sig_xi = gaus1_yield_xi/gausTot_yield_xi;
        double rms_gaus2_sig_xi = gaus2_yield_xi/gausTot_yield_xi;
        double rms_true_xi = TMath::Sqrt(rms_gaus1_sig_xi*sigma1.getVal()*sigma1.getVal() + rms_gaus2_sig_xi*sigma2.getVal()*sigma2.getVal());

        x.setRange("peak", mean.getVal() - 2*rms_true_xi, mean.getVal() + 2*rms_true_xi);
        RooAbsReal* Intgaus1_xi      = gaus1.createIntegral(x, x,  "peak");
        RooAbsReal* Intgaus2_xi      = gaus2.createIntegral(x, x, "peak");
        RooAbsReal* Intbackground_xi = background.createIntegral(x, x, "peak");

        double Intgaus1E_xi      = gaus1F_xi*Intgaus1_xi->getVal();
        double Intgaus2E_xi      = gaus2F_xi*Intgaus2_xi->getVal();
        double IntbackgroundE_xi = qsig_xi*Intbackground_xi->getVal();
        double totsig_xi         = Intgaus1E_xi + Intgaus2E_xi + IntbackgroundE_xi;
        double Yield_xi          = Intgaus1E_xi + Intgaus2E_xi;


        double Fsig_xi = Yield_xi/totsig_xi;

        double significance = Yield_xi/sqrt(totsig_xi);

        //cout << "Yield (xi): " << Yield_xi << endl;
        //cout << "Fsig (xi): " << Fsig_xi << endl;
        //cout << "std (xi): "  << rms_true_xi  << endl;
        //cout << "mass (xi): " << mean_xi << endl;

        //cout << "covQual (xi)" << covQual << endl;
        //cout << "Signal Sig (xi)" << significance << endl;


        sum.plotOn(xframe_,Name("sum"),NormRange("cut"),LineWidth(1),LineColor(kBlue));
        sum.plotOn(xframe_,Components(background),NormRange("cut"),LineStyle(kDashed),LineWidth(1),LineColor(kBlue));
        c_om2->cd(i+1);
        xframe_->Draw();
        TLine* t1 = new TLine(mean.getVal() - 2*rms_true_xi, 0, mean.getVal() - 2*rms_true_xi, gPad->GetUymax());
        TLine* t2 = new TLine(mean.getVal() + 2*rms_true_xi, 0, mean.getVal() + 2*rms_true_xi, gPad->GetUymax());
        t1->SetLineStyle(2);
        t1->SetLineColor(kGreen);
        t2->SetLineStyle(2);
        t2->SetLineColor(kGreen);
        t1->Draw("same");
        t2->Draw("same");

        double xpos = 0.64;
        double ypos = 0.85;
        double increment = 0.07;
        os << "Fsig: " << Fsig_xi;
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "Mean: " << std::setprecision(5) << mean_xi << " GeV" << std::setprecision(6);
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "#sigma :" << std::setprecision(2) << rms_true_xi << " GeV" << std::setprecision(6);
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "CovQual: " << covQual;
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        osYield << "S: " << std::setprecision(2) << Yield_xi;
        tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
        osYield.str(std::string());
        osYield << "B: " << std::setprecision(2) << IntbackgroundE_xi;
        tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
        osYield.str(std::string());
        osYield << "S-B: " << std::setprecision(2) << Yield_xi - IntbackgroundE_xi;
        tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
        osYield.str(std::string());
        os.str(std::string());
        os << "Sig: " << significance;
        tex->DrawLatex(xpos,0.85,os.str().c_str());
        os.str(std::string());
    }

    c_om2->Print("OmegaPtRebin_PeriSub.pdf");
    c_om2->Print("OmegaPtRebin_PeriSub.png");

}
