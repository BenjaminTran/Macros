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

std::string v2RootFileName = "v2V0sGraphs185_250D0Ana.root";

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

void Rap_v2sig()
{
    MITStyle();
    TCanvas* c1 = MakeCanvas("c1", "Plot");
    c1->cd();
    /*c1->SetLogy();*/
    c1->SetLeftMargin(0.12);

    // draw the frame using a histogram frame
    TH1F* frame;

    const int xi_npoints = 8;
    const int ks_npoints = 13;
    const int la_npoints = 10;

    std::vector<double> v2Xi8;
    std::vector<double> pTXi8;
    std::vector<double> v2Xi8E;

    std::vector<double> v2Ks8;
    std::vector<double> pTKs8;
    std::vector<double> v2Ks8E;

    std::vector<double> v2La8;
    std::vector<double> pTLa8;
    std::vector<double> v2La8E;

    for(int i=0; i<3; i++)
    {
        if(i==0) // v2 vs pt
        {
            v2Xi8.insert(v2Xi8.end(),{0.0371521 ,0.0669904 ,0.0877339 ,0.123959 ,0.163467 ,0.183918 ,0.209974 ,0.17582});
            pTXi8.insert(pTXi8.end(),{1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.474});
            v2Xi8E.insert(v2Xi8E.end(),{0.00932224 ,0.00484125 ,0.0041527 ,0.00326934 ,0.00309064 ,0.00361334 ,0.00517473 ,0.0113186});

            v2Ks8.insert(v2Ks8.end(),{0.0135543 ,0.0293041 ,0.0433881 ,0.0594498 ,0.082389 ,0.106184 ,0.124363 ,0.137749 ,0.146475 ,0.145182 ,0.135266 ,0.125505 ,0.124708});
            pTKs8.insert(pTKs8.end(),{0.3666, 0.5309, 0.711, 0.9046, 1.202, 1.591, 1.986, 2.465, 3.136, 4.008, 5.142, 6.431, 7.619, 9.142, 11.64, 16.86});
            v2Ks8E.insert(v2Ks8E.end(),{0.00217719 ,0.000522063 ,0.000305164 ,0.000248705 ,0.000161989 ,0.000173887 ,0.000208535 ,0.000225897 ,0.000298544 ,0.000448324 ,0.000705022 ,0.00146462 ,0.00192261});

            v2La8.insert(v2La8.end(),{0.0345669 ,0.0510091 ,0.0771498 ,0.106513 ,0.138504 ,0.171975 ,0.19385 ,0.202724 ,0.19716 ,0.178541});
            pTLa8.insert(pTLa8.end(),{0.9252, 1.224, 1.603, 1.995, 2.485, 3.156, 4.007, 5.116, 6.414, 7.588});
            v2La8E.insert(v2La8E.end(),{0.00106898 ,0.000437383 ,0.000387893 ,0.000396151 ,0.000360204 ,0.000398977 ,0.000561998 ,0.000940498 ,0.00223511 ,0.0032202});

            frame = c1->DrawFrame(0,-0.01,8,0.5);
            gPad->SetTickx();
            gPad->SetTicky();
            frame->GetXaxis()->CenterTitle(1);
            frame->GetYaxis()->CenterTitle(1);
            frame->GetXaxis()->SetTitleSize(0.05);
            frame->GetYaxis()->SetTitleSize(0.05);
            frame->SetTitleOffset(1.1,"Y");
            frame->SetTitleOffset(1.2,"X");
            frame->GetXaxis()->SetTitle("p_{T} (GeV)");
            frame->GetYaxis()->SetTitle("v_{2}^{sig}");
        }
        else if(i==1) // v2 vs KET
        {
            v2Xi8.clear();
            pTXi8.clear();
            v2Xi8E.clear();

            v2Ks8.clear();
            pTKs8.clear();
            v2Ks8E.clear();

            v2La8.clear();
            pTLa8.clear();
            v2La8E.clear();

            v2Xi8.insert(v2Xi8.end(),{0.0371521 ,0.0669904 ,0.0877339 ,0.123959 ,0.163467 ,0.183918 ,0.209974 ,0.17582});
            pTXi8.insert(pTXi8.end(),{0.515319 ,0.768532 ,1.07944 ,1.50413 ,2.12034 ,2.92705 ,3.98122 ,5.28277});
            v2Xi8E.insert(v2Xi8E.end(),{0.00932224 ,0.00484125 ,0.0041527 ,0.00326934 ,0.00309064 ,0.00361334 ,0.00517473 ,0.0113186});

            v2Ks8.insert(v2Ks8.end(),{0.0135543 ,0.0293041 ,0.0433881 ,0.0594498 ,0.082389 ,0.106184 ,0.124363 ,0.137749 ,0.146475 ,0.145182 ,0.135266 ,0.125505 ,0.124708});
            pTKs8.insert(pTKs8.end(),{0.120493 ,0.230933 ,0.370907 ,0.535241 ,0.804172 ,1.17001 ,1.55031 ,2.01776 ,2.67887 ,3.54188 ,4.67169 ,5.95357 ,7.13835});
            v2Ks8E.insert(v2Ks8E.end(),{0.00217719 ,0.000522063 ,0.000305164 ,0.000248705 ,0.000161989 ,0.000173887 ,0.000208535 ,0.000225897 ,0.000298544 ,0.000448324 ,0.000705022 ,0.00146462 ,0.00192261});

            v2La8.insert(v2La8.end(),{0.0345669 ,0.0510091 ,0.0771498 ,0.106513 ,0.138504 ,0.171975 ,0.19385 ,0.202724 ,0.19716 ,0.178541});
            pTLa8.insert(pTLa8.end(),{0.326395 ,0.537952 ,0.836813 ,1.16991 ,1.60772 ,2.23092 ,3.0437 ,4.1234 ,5.39811 ,6.55991});
            v2La8E.insert(v2La8E.end(),{0.00106898 ,0.000437383 ,0.000387893 ,0.000396151 ,0.000360204 ,0.000398977 ,0.000561998 ,0.000940498 ,0.00223511 ,0.0032202});

            frame = c1->DrawFrame(0,-0.01,8,0.5);
            gPad->SetTickx();
            gPad->SetTicky();
            frame->GetXaxis()->CenterTitle(1);
            frame->GetYaxis()->CenterTitle(1);
            frame->GetXaxis()->SetTitleSize(0.05);
            frame->GetYaxis()->SetTitleSize(0.05);
            frame->SetTitleOffset(1.1,"Y");
            frame->SetTitleOffset(1.2,"X");
            frame->GetXaxis()->SetTitle("KE_{T} (GeV)");
            frame->GetYaxis()->SetTitle("v_{2}^{sig}");

        }
        else // v2/nq vs KET/nq
        {
            v2Xi8.clear();
            pTXi8.clear();
            v2Xi8E.clear();

            v2Ks8.clear();
            pTKs8.clear();
            v2Ks8E.clear();

            v2La8.clear();
            pTLa8.clear();
            v2La8E.clear();

            v2Xi8.insert(v2Xi8.end(),{0.012384 ,0.0223301 ,0.0292446 ,0.0413197 ,0.0544889 ,0.0613061 ,0.0699914 ,0.0586068});
            pTXi8.insert(pTXi8.end(),{0.171773 ,0.256177 ,0.359813 ,0.501378 ,0.70678 ,0.975684 ,1.32707 ,1.76092});
            v2Xi8E.insert(v2Xi8E.end(),{0.00310741 ,0.00161375 ,0.00138423 ,0.00108978 ,0.00103021 ,0.00120445 ,0.00172491 ,0.00377288});

            v2Ks8.insert(v2Ks8.end(),{0.00677717 ,0.0146521 ,0.021694 ,0.0297249 ,0.0411945 ,0.0530921 ,0.0621816 ,0.0688744 ,0.0732373 ,0.0725908 ,0.0676328 ,0.0627523 ,0.062354});
            pTKs8.insert(pTKs8.end(),{0.0602465 ,0.115467 ,0.185454 ,0.26762 ,0.402086 ,0.585004 ,0.775155 ,1.00888 ,1.33943 ,1.77094 ,2.33584 ,2.97678 ,3.56917});
            v2Ks8E.insert(v2Ks8E.end(),{0.0010886 ,0.000261031 ,0.000152582 ,0.000124353 ,8.09944e-05 ,8.69433e-05 ,0.000104268 ,0.000112948 ,0.000149272 ,0.000224162 ,0.000352511 ,0.00073231 ,0.000961305});

            v2La8.insert(v2La8.end(),{0.0115223 ,0.017003 ,0.0257166 ,0.0355043 ,0.0461679 ,0.0573251 ,0.0646166 ,0.0675745 ,0.0657199 ,0.0595136});
            pTLa8.insert(pTLa8.end(),{0.108798 ,0.179317 ,0.278938 ,0.389971 ,0.535908 ,0.74364 ,1.01457 ,1.37447 ,1.79937 ,2.18664});
            v2La8E.insert(v2La8E.end(),{0.000356326 ,0.000145794 ,0.000129298 ,0.00013205 ,0.000120068 ,0.000132992 ,0.000187333 ,0.000313499 ,0.000745036 ,0.0010734});

            frame = c1->DrawFrame(0,-0.01,5,0.2);
            gPad->SetTickx();
            gPad->SetTicky();
            frame->GetXaxis()->CenterTitle(1);
            frame->GetYaxis()->CenterTitle(1);
            frame->GetXaxis()->SetTitleSize(0.05);
            frame->GetYaxis()->SetTitleSize(0.05);
            frame->SetTitleOffset(1.1,"Y");
            frame->SetTitleOffset(1.2,"X");
            frame->GetXaxis()->SetTitle("KE_{T}/n_{q} (GeV)");
            frame->GetYaxis()->SetTitle("v_{2}^{sig}/n_{q}");
        }

        // Pull TGraph for Kshort and lambda
        TFile* file_pPbv2 = TFile::Open("lrgraphv2_v3_pPb_185-220.root");
        TFile* file_hadv2 = TFile::Open("lrgraphv2_v3_pPb_hadron_185-above.root");

        TGraphErrors* ks_v2 = (TGraphErrors*)file_pPbv2->Get("kshortv2true");
        TGraphErrors* la_v2 = (TGraphErrors*)file_pPbv2->Get("lambdav2true");
        //TGraphErrors* ha_v2 = (TGraphErrors*)file_hadv2->Get("hadronv2");
        TGraphErrors* ha_v2 = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033pPb_v2.txt"),1,28,1.2);

        TGraphErrors* xi8_v2 = new TGraphErrors(xi_npoints,&pTXi8[0],&v2Xi8[0],0,&v2Xi8E[0]);
        TGraphErrors* ks8_v2 = new TGraphErrors(ks_npoints,&pTKs8[0],&v2Ks8[0],0,&v2Ks8E[0]);
        TGraphErrors* la8_v2 = new TGraphErrors(la_npoints,&pTLa8[0],&v2La8[0],0,&v2La8E[0]);

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


        //ha_v2->Draw("PESAME");
        ks8_v2->Draw("P");
        la8_v2->Draw("P");
        xi8_v2->Draw("P");

        std::string kshortv2 = "";
        std::string lambdav2 = "";
        std::string cascadev2 = "";

        // Draw Legend and write points into rootfile
        TFile* out = NULL;
        if(i==0)
        {
            kshortv2 = "kshortv2";
            lambdav2 = "lambdav2";
            cascadev2 = "cascadev2";
            out = new TFile(v2RootFileName.c_str(),"RECREATE");
        }
        else if(i==1)
        {
            kshortv2 = "kshortv2KET";
            lambdav2 = "lambdav2KET";
            cascadev2 = "cascadev2KET";
            out = new TFile(v2RootFileName.c_str(),"UPDATE");
        }
        else
        {
            kshortv2 = "kshortv2nq";
            lambdav2 = "lambdav2nq";
            cascadev2 = "cascadev2nq";
            out = new TFile(v2RootFileName.c_str(),"UPDATE");
        }
        ks8_v2->Write(kshortv2.c_str());
        la8_v2->Write(lambdav2.c_str());
        xi8_v2->Write(cascadev2.c_str());
        out->Close();


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

        if(i==0)
        {
            c1->Print("v2SigRapidity.pdf");
        }
        else if(i==1)
        {
            c1->Print("v2SigRapidityKET.pdf");
        }
        else
        {
            c1->Print("v2SigRapidityDividednq.pdf");
        }
    }
}

void Rap_v2obs()
{
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

    const int xi_npoints = 8;
    const int ks_npoints = 13;
    const int la_npoints = 10;

    std::string kshortv2 = "";
    std::string lambdav2 = "";
    std::string cascadev2 = "";
    std::string yaxis = "";

    std::vector<double> v2Xi8;
    std::vector<double> pTXi8;
    std::vector<double> v2Xi8E;

    std::vector<double> v2Ks8;
    std::vector<double> pTKs8;
    std::vector<double> v2Ks8E;

    std::vector<double> v2La8;
    std::vector<double> pTLa8;
    std::vector<double> v2La8E;

    for(int i=0; i<2; i++)
    {
        if(i==0)
        {
            kshortv2 = "kshortv2Obs";
            lambdav2 = "lambdav2Obs";
            cascadev2 = "cascadev2Obs";
            yaxis = "v_{2}^{obs}";
            v2Xi8.insert(v2Xi8.end(),{0.0397521 ,0.0672835 ,0.0879851 ,0.123906 ,0.163679 ,0.184638 ,0.210454 ,0.177214});
            pTXi8.insert(pTXi8.end(),{1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.938});
            v2Xi8E.insert(v2Xi8E.end(),{0.00885641 ,0.00470611 ,0.0040499 ,0.00319483 ,0.0030202 ,0.0035308 ,0.00505887 ,0.0110131});

            v2Ks8.insert(v2Ks8.end(),{0.0135496 ,0.0293041 ,0.0433881 ,0.0594498 ,0.082389 ,0.106184 ,0.124363 ,0.137749 ,0.146475 ,0.145182 ,0.135266 ,0.125505 ,0.12497});
            pTKs8.insert(pTKs8.end(),{0.3666, 0.5309, 0.711, 0.9046, 1.202, 1.591, 1.986, 2.465, 3.136, 4.008, 5.142, 6.431, 7.619});
            v2Ks8E.insert(v2Ks8E.end(),{0.00217605 ,0.000522193 ,0.000305651 ,0.000249826 ,0.000165266 ,0.000178936 ,0.000214314 ,0.000232438 ,0.000304168 ,0.000452022 ,0.000707061 ,0.00146546 ,0.00190789});

            v2La8.insert(v2La8.end(),{0.0346226 ,0.0510097 ,0.0771498 ,0.106527 ,0.13851 ,0.171976 ,0.19385 ,0.202711 ,0.197155 ,0.178511});
            pTLa8.insert(pTLa8.end(),{0.9252, 1.224, 1.603, 1.995, 2.485, 3.156, 4.007, 5.116, 6.414, 7.588});
            v2La8E.insert(v2La8E.end(),{0.00106781 ,0.000437847 ,0.000389104 ,0.00039822 ,0.000364256 ,0.000404734 ,0.000567087 ,0.000942798 ,0.00223597 ,0.00321154});
        }
        else
        {
            v2Xi8.clear();
            pTXi8.clear();
            v2Xi8E.clear();

            v2Ks8.clear();
            pTKs8.clear();
            v2Ks8E.clear();

            v2La8.clear();
            pTLa8.clear();
            v2La8E.clear();

            kshortv2 = "kshortv2Bkg";
            lambdav2 = "lambdav2Bkg";
            cascadev2 = "cascadev2Bkg";
            yaxis = "v_{2}^{bkg}";

            v2Xi8.insert(v2Xi8.end(),{0.0936965 ,0.0782133 ,0.0985179 ,0.121518 ,0.173158 ,0.216723 ,0.232391 ,0.230779});
            pTXi8.insert(pTXi8.end(),{1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.938});
            v2Xi8E.insert(v2Xi8E.end(),{0.00749852 ,0.00950025 ,0.00868329 ,0.00770257 ,0.00664219 ,0.00701412 ,0.010755 ,0.0195874});

            v2Ks8.insert(v2Ks8.end(),{0.00449492 ,0.0473127 ,0.0585964 ,0.0687393 ,0.0931208 ,0.114931 ,0.132041 ,0.149003 ,0.166686 ,0.172154 ,0.171293 ,0.158914 ,0.158798});
            pTKs8.insert(pTKs8.end(),{0.3666, 0.5309, 0.711, 0.9046, 1.202, 1.591, 1.986, 2.465, 3.136, 4.008, 5.142, 6.431, 7.619});
            v2Ks8E.insert(v2Ks8E.end(),{0.0390969 ,0.00255086 ,0.00181359 ,0.00163701 ,0.000986848 ,0.000960547 ,0.0010765 ,0.00101774 ,0.00108115 ,0.00148885 ,0.00198329 ,0.00404482 ,0.00497703});

            v2La8.insert(v2La8.end(),{0.0817952 ,0.0973297 ,0.117934 ,0.136615 ,0.156659 ,0.178532 ,0.193261 ,0.191957 ,0.177981 ,0.168276});
            pTLa8.insert(pTLa8.end(),{0.9252, 1.224, 1.603, 1.995, 2.485, 3.156, 4.007, 5.116, 6.414, 7.588});
            v2La8E.insert(v2La8E.end(),{0.000541466 ,0.000463843 ,0.000595759 ,0.00077596 ,0.000850846 ,0.00108663 ,0.00159224 ,0.00261702 ,0.00559253 ,0.00709033});

            c2->cd();

            leg_co_spec->Draw();
            leg_co_Label->Draw();
            tex->SetTextFont(62);
            tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
            tex->SetTextSize(0.045);
            tex->SetTextFont(42);
            tex->DrawLatex(0.40,0.24,"185 #leq N_{trk}^{offline} < 250");
            c1->cd();
        }

        TGraphErrors* xi8_v2 = new TGraphErrors(xi_npoints,&pTXi8[0],&v2Xi8[0],0,&v2Xi8E[0]);
        TGraphErrors* ks8_v2 = new TGraphErrors(ks_npoints,&pTKs8[0],&v2Ks8[0],0,&v2Ks8E[0]);
        TGraphErrors* la8_v2 = new TGraphErrors(la_npoints,&pTLa8[0],&v2La8[0],0,&v2La8E[0]);

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

        // draw the frame using a histogram frame
        TH1F* frame = c1->DrawFrame(0,-0.05,9,0.45);
        /*TH1F* frame = c1->DrawFrame(0,0.01,20,1);*/
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
        leg->AddEntry(xi8_v2, "#Xi^{+}/#Xi^{-}", "P");
        leg->Draw();

        tex->SetTextFont(62);
        tex->DrawLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
        tex->SetTextSize(0.045);
        //tex->DrawLatex(0.15,0.72, "L_{#lower[-0.25]{int}} = #color[38]{35} nb^{#font[122]{\55}1}, #color[46]{62} nb^{#font[122]{\55}1}");
        // tex->DrawLatex(0.23,0.72, "L_{#lower[-0.25]{int}} = #color[kOrange+8]{35} nb^{#font[122]{\55}1}, 62 nb^{#font[122]{\55}1}");
        tex->SetTextFont(42);
        tex->DrawLatex(0.40,0.24,"185 #leq N_{trk}^{offline} < 250");
        /*tex->DrawLatex(0.15,0.74,"|y| < 1");*/
        //tex->DrawLatex(0.4,0.7, "L_{#lower[-0.25]{int}} = 35 nb^{#font[122]{\55}1}, 185 nb^{#font[122]{\55}1}");

        if(i==0) c1->Print("v2ObsRapidity.pdf");
        if(i==1) c1->Print("v2BkgRapidity.pdf");

        c2->cd();
        if(i==1)
        {
            ks8_v2->SetMarkerStyle(24);
            xi8_v2->SetMarkerStyle(25);
            la8_v2->SetMarkerStyle(26);
        }
        ks8_v2->Draw("P");
        la8_v2->Draw("P");
        xi8_v2->Draw("P");

        if(i==0)
        {
            leg_co_spec->AddEntry(ks8_v2, "K_{S}^{0}", "P");
            leg_co_spec->AddEntry(la8_v2, "#Lambda/#bar{#Lambda}", "P");
            leg_co_spec->AddEntry(xi8_v2, "#Xi^{+}/#Xi^{-}", "P");
            legendLabel = "Obs";
            leg_co_Label->AddEntry(ks8_v2,legendLabel.c_str(), "P");
        }
        else
        {
            legendLabel = "Bkg";
            TGraphErrors* ks_clone = (TGraphErrors*)ks8_v2->Clone("ks8_v2_bkg");
            ks_clone->SetMarkerStyle(24);
            leg_co_Label->AddEntry(ks_clone,legendLabel.c_str(), "P");
        }

        c1->cd();

        // Write points into rootfile
        TFile out(v2RootFileName.c_str(),"UPDATE");
        ks8_v2->Write(kshortv2.c_str(),TObject::kOverwrite);
        la8_v2->Write(lambdav2.c_str(),TObject::kOverwrite);
        xi8_v2->Write(cascadev2.c_str(),TObject::kOverwrite);
        out.Close();
    }

    c2->Print("v2ObsBkgRapidity.pdf");
}

void Rap_v2bkg()
{
    MITStyle();

    const int xi_npoints = 8;
    double v2Xi8[xi_npoints]  = {0.0936965 ,0.0782133 ,0.0985179 ,0.121518 ,0.173158 ,0.216723 ,0.232391 ,0.230779};
    double pTXi8[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.938};
    double v2Xi8E[xi_npoints] = {0.00749852 ,0.00950025 ,0.00868329 ,0.00770257 ,0.00664219 ,0.00701412 ,0.010755 ,0.0195874};

	const int ks_npoints = 13;
    double v2Ks8[ks_npoints]  = {0.00449492 ,0.0473127 ,0.0585964 ,0.0687393 ,0.0931208 ,0.114931 ,0.132041 ,0.149003 ,0.166686 ,0.172154 ,0.171293 ,0.158914 ,0.158798};
    double pTKs8[ks_npoints]  = {0.3666, 0.5309, 0.711, 0.9046, 1.202, 1.591, 1.986, 2.465, 3.136, 4.008, 5.142, 6.431, 7.619};
    double v2Ks8E[ks_npoints] = {0.0390969 ,0.00255086 ,0.00181359 ,0.00163701 ,0.000986848 ,0.000960547 ,0.0010765 ,0.00101774 ,0.00108115 ,0.00148885 ,0.00198329 ,0.00404482 ,0.00497703};

	const int la_npoints = 10;
    double v2La8[la_npoints]  = {0.0817952 ,0.0973297 ,0.117934 ,0.136615 ,0.156659 ,0.178532 ,0.193261 ,0.191957 ,0.177981 ,0.168276};
    double pTLa8[la_npoints]  = {0.9252, 1.224, 1.603, 1.995, 2.485, 3.156, 4.007, 5.116, 6.414, 7.588};
    double v2La8E[la_npoints] = {0.000541466 ,0.000463843 ,0.000595759 ,0.00077596 ,0.000850846 ,0.00108663 ,0.00159224 ,0.00261702 ,0.00559253 ,0.00709033};

    TGraphErrors* xi8_v2 = new TGraphErrors(xi_npoints,pTXi8,v2Xi8,0,v2Xi8E);
    TGraphErrors* ks8_v2 = new TGraphErrors(ks_npoints,pTKs8,v2Ks8,0,v2Ks8E);
	TGraphErrors* la8_v2 = new TGraphErrors(la_npoints,pTLa8,v2La8,0,v2La8E);

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

    TH1F* frame = c1->DrawFrame(0,-0.05,9,0.45);
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

    ks8_v2->Draw("P");
    la8_v2->Draw("P");
    xi8_v2->Draw("P");

    // Write points into rootfile
    TFile out(v2RootFileName.c_str(),"UPDATE");
    ks8_v2->Write("kshortv2bkg",TObject::kOverwrite);
    la8_v2->Write("lambdav2bkg",TObject::kOverwrite);
    xi8_v2->Write("cascadev2bkg",TObject::kOverwrite);
    out.Close();

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

    const int xi_npoints = 8;
    double fsig_Xi8[xi_npoints]  = {0.954019 ,0.973881 ,0.976705 ,0.97829 ,0.978074 ,0.978057 ,0.978603 ,0.974644};
    double pTXi8[xi_npoints]  = {1.267, 1.62, 2.008, 2.501, 3.173, 4.029, 5.055, 6.938};

	const int ks_npoints = 13;
    double fsig_Ks8[ks_npoints]  = {0.999476 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,0.999999 ,0.999999 ,0.999988 ,0.999997 ,0.992327};
    double pTKs8[ks_npoints]  = {0.3666, 0.5309, 0.711, 0.9046, 1.202, 1.591, 1.986, 2.465, 3.136, 4.008, 5.142, 6.431, 7.619};

	const int la_npoints = 10;
    double fsig_La8[la_npoints]  = {0.99882 ,0.999987 ,1 ,0.999524 ,0.999632 ,0.999855 ,0.999698 ,0.998783 ,0.999771 ,0.997088};
    double pTLa8[la_npoints]  = {0.9252, 1.224, 1.603, 1.995, 2.485, 3.156, 4.007, 5.116, 6.414, 7.588};

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

    TH1F* frame = c1->DrawFrame(0,0.85,9,1.15);
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
    TLine* line = new TLine(0,1,9,1);
    line->SetLineStyle(2);
    line->Draw("same");

    c1->Print("fsig_Rapidity.pdf");

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

void RapSys_Closure()
{
    MITStyle();
    gStyle->SetTitleAlign(33);
    TVirtualFitter::SetMaxIterations( 300000 );
    const int ks_npoints = 13;
    const int la_npoints = 10;
    std::vector<double> PtMeanReco_ks;
    std::vector<double> PtMeanReco_la;
    std::vector<double> PtMeanRecoMatch_ks;
    std::vector<double> PtMeanRecoMatch_la;

    TFile* f_Reco = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0CorrelationClosureTotal_08_28_2017.root");
    TFile* f_Match = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/MatchV0ClosureBpPb_09_20_17.root"); //1 | Forgot to cut on abs(dpt/pt_gen)

    //Determine Mean pt points for Ks
    for(int i=0; i<ks_npoints; i++)
    {
        TH1D* pt_ks = (TH1D*)f_Reco->Get(Form("v0CorrelationRapidityMC/Ptkshort_pt%d",i));
        TH1D* pt_bkg_ks = (TH1D*)f_Reco->Get(Form("v0CorrelationRapidityMC/Ptkshort_bkg_pt%d",i));
        double PtMean = 0;
        int n = 0;
        for(int j=pt_ks->FindFirstBinAbove(0,1); j<pt_ks->FindLastBinAbove(0,1); j++)
        {
            PtMean += pt_ks->GetBinCenter(j)*(pt_ks->GetBinContent(j));
            n += pt_ks->GetBinContent(j);
        }
        for(int j=pt_bkg_ks->FindFirstBinAbove(0,1); j<pt_bkg_ks->FindLastBinAbove(0,1); j++)
        {
            PtMean += pt_bkg_ks->GetBinCenter(j)*(pt_bkg_ks->GetBinContent(j));
            n += pt_bkg_ks->GetBinContent(j);
        }
        PtMeanReco_ks.push_back(PtMean/n);
        cout << i << ": " << PtMean/n << endl;
    }
    //Determine Mean pt points for La
    for(int i=0; i<la_npoints; i++)
    {
        TH1D* pt_la = (TH1D*)f_Reco->Get(Form("v0CorrelationRapidityMC/Ptlambda_pt%d",i));
        TH1D* pt_bkg_la = (TH1D*)f_Reco->Get(Form("v0CorrelationRapidityMC/Ptlambda_bkg_pt%d",i));
        double PtMean = 0;
        int n = 0;
        for(int j=pt_la->FindFirstBinAbove(0,1); j<pt_la->FindLastBinAbove(0,1); j++)
        {
            PtMean += pt_la->GetBinCenter(j)*(pt_la->GetBinContent(j));
            n += pt_la->GetBinContent(j);
        }
        for(int j=pt_bkg_la->FindFirstBinAbove(0,1); j<pt_bkg_la->FindLastBinAbove(0,1); j++)
        {
            PtMean += pt_bkg_la->GetBinCenter(j)*(pt_bkg_la->GetBinContent(j));
            n += pt_bkg_la->GetBinContent(j);
        }
        PtMeanReco_la.push_back(PtMean/n);
    }

    //MATCHED
    //Determine Mean pt points for Ks
    for(int i=0; i<ks_npoints; i++)
    {
        TH1D* pt_ks = (TH1D*)f_Match->Get(Form("v0CorrelationRapidityMatchMC/Ptkshort_pt%d",i));
        TH1D* pt_bkg_ks = (TH1D*)f_Match->Get(Form("v0CorrelationRapidityMatchMC/Ptkshort_bkg_pt%d",i));
        double PtMean = 0;
        int n = 0;
        for(int j=pt_ks->FindFirstBinAbove(0,1); j<pt_ks->FindLastBinAbove(0,1); j++)
        {
            PtMean += pt_ks->GetBinCenter(j)*(pt_ks->GetBinContent(j));
            n += pt_ks->GetBinContent(j);
        }
        for(int j=pt_bkg_ks->FindFirstBinAbove(0,1); j<pt_bkg_ks->FindLastBinAbove(0,1); j++)
        {
            PtMean += pt_bkg_ks->GetBinCenter(j)*(pt_bkg_ks->GetBinContent(j));
            n += pt_bkg_ks->GetBinContent(j);
        }
        PtMeanRecoMatch_ks.push_back(PtMean/n);
        cout << i << ": " << PtMean/n << endl;
    }
    //Determine Mean pt points for La
    for(int i=0; i<la_npoints; i++)
    {
        TH1D* pt_la = (TH1D*)f_Match->Get(Form("v0CorrelationRapidityMatchMC/Ptlambda_pt%d",i));
        TH1D* pt_bkg_la = (TH1D*)f_Match->Get(Form("v0CorrelationRapidityMatchMC/Ptlambda_bkg_pt%d",i));
        double PtMean = 0;
        int n = 0;
        for(int j=pt_la->FindFirstBinAbove(0,1); j<pt_la->FindLastBinAbove(0,1); j++)
        {
            PtMean += pt_la->GetBinCenter(j)*(pt_la->GetBinContent(j));
            n += pt_la->GetBinContent(j);
        }
        for(int j=pt_bkg_la->FindFirstBinAbove(0,1); j<pt_bkg_la->FindLastBinAbove(0,1); j++)
        {
            PtMean += pt_bkg_la->GetBinCenter(j)*(pt_bkg_la->GetBinContent(j));
            n += pt_bkg_la->GetBinContent(j);
        }
        PtMeanRecoMatch_la.push_back(PtMean/n);
    }

    double v2Reco_ks[ks_npoints]  = {0.0722071 ,0.147723 ,0.222303 ,0.294967 ,0.369263 ,0.445429 ,0.491016 ,0.521116 ,0.533933 ,0.519577 ,0.464798 ,0.381484 ,0.290978};
    double* pTReco_ks = &PtMeanReco_ks[0];
    double v2Reco_ksE[ks_npoints] = {0.00437059 ,0.00115752 ,0.000754802 ,0.000676901 ,0.000484584 ,0.000557721 ,0.0006959 ,0.000776143 ,0.00105449 ,0.00166647 ,0.00308999 ,0.0067008 ,0.00956059};

    double v2RecoMatch_ks[ks_npoints]  = {0.0553694 ,0.145913 ,0.220285 ,0.294232 ,0.369009 ,0.446266 ,0.490768 ,0.521643 ,0.532904 ,0.520321 ,0.473552 ,0.392277 ,0.287408};
    double* pTRecoMatch_ks = &PtMeanRecoMatch_ks[0];
    double v2RecoMatch_ksE[ks_npoints] = {0.00365939 ,0.00126627 ,0.000890618 ,0.000837354 ,0.000628735 ,0.000747037 ,0.000940407 ,0.00105765 ,0.00144895 ,0.00229439 ,0.00427065 ,0.0090483 ,0.0130823};

    double v2Gen_ks[ks_npoints]  = {0.052706 ,0.14045 ,0.229782 ,0.307633 ,0.384587 ,0.467894 ,0.519706 ,0.551816 ,0.565919 ,0.554509 ,0.49801 ,0.406751 ,0.301587};
    double pTGen_ks[ks_npoints]  = {0.303, 0.4971, 0.6951, 0.8948, 1.178, 1.578, 1.979, 2.456, 3.129, 4, 5.131, 6.427, 7.61};
    double v2Gen_ksE[ks_npoints] = {0.000121163 ,0.00012872 ,0.000144935 ,0.000167038 ,0.00015137 ,0.000204732 ,0.000277958 ,0.000327538 ,0.000456264 ,0.000720603 ,0.00119137 ,0.00266441 ,0.00361488};

    double v2Reco_la[la_npoints]  = {0.191419 ,0.257195 ,0.359238 ,0.449797 ,0.531741 ,0.592744 ,0.63408 ,0.638523 ,0.649031 ,0.630982};
    double* pTReco_la = &PtMeanReco_la[0];
    double v2Reco_laE[la_npoints] = {0.00263104 ,0.00117218 ,0.00113757 ,0.00124031 ,0.0011602 ,0.00128685 ,0.00180476 ,0.00313557 ,0.0111194 ,0.0282399};

    double v2RecoMatch_la[la_npoints]  = {0.186242 ,0.252737 ,0.359515 ,0.449549 ,0.52772 ,0.591169 ,0.629537 ,0.637112 ,0.616482 ,0.544056};
    double* pTRecoMatch_la = &PtMeanRecoMatch_la[0];
    double v2RecoMatch_laE[la_npoints] = {0.00275591 ,0.00145987 ,0.00143559 ,0.00156058 ,0.00149181 ,0.00170974 ,0.00243067 ,0.00416321 ,0.0155256 ,0.0435596};

    double v2Gen_la[la_npoints]  = {0.176432 ,0.261996 ,0.374751 ,0.473161 ,0.557898 ,0.629242 ,0.66937 ,0.679275 ,0.663517 ,0.628933};
    double pTGen_la[la_npoints]  = {0.8981, 1.189, 1.586, 1.985, 2.466, 3.137, 4.002, 5.118, 6.412, 7.568};
    double v2Gen_laE[la_npoints] = {0.000204122 ,0.000165959 ,0.000198333 ,0.000243174 ,0.000260547 ,0.000333439 ,0.000501804 ,0.000835003 ,0.0019952 ,0.00308793};

    // Pull TGraph for Kshort and lambda
    TGraphErrors* Reco_v2_ks = new TGraphErrors(ks_npoints,pTReco_ks,v2Reco_ks,0,v2Reco_ksE);
    TGraphErrors* RecoMatch_v2_ks = new TGraphErrors(ks_npoints,pTRecoMatch_ks,v2RecoMatch_ks,0,v2RecoMatch_ksE);
    TGraphErrors* Gen_v2_ks = new TGraphErrors(ks_npoints,pTGen_ks,v2Gen_ks,0,v2Gen_ksE);

    Gen_v2_ks->SetMarkerColor(kRed);
    Gen_v2_ks->SetMarkerStyle(20);
    Gen_v2_ks->SetMarkerSize(1.5);
    Gen_v2_ks->SetLineColor(kRed);

    Reco_v2_ks->SetMarkerColor(kRed);
    Reco_v2_ks->SetMarkerStyle(25);
    Reco_v2_ks->SetMarkerSize(1.5);
    Reco_v2_ks->SetLineColor(kRed);

    RecoMatch_v2_ks->SetMarkerColor(kGreen);
    RecoMatch_v2_ks->SetMarkerStyle(25);
    RecoMatch_v2_ks->SetMarkerSize(1.5);
    RecoMatch_v2_ks->SetLineColor(kGreen);

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
    Reco_v2_ks->Draw("P");
    RecoMatch_v2_ks->Draw("P");

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
    leg_ks->AddEntry(Reco_v2_ks, "Reco K_{S}^{0}", "P");
    leg_ks->AddEntry(RecoMatch_v2_ks, "Reco Matched K_{S}^{0}", "P");
    leg_ks->Draw();

    tex->SetTextFont(62);
    tex->SetTextSize(0.045);
    tex->DrawLatex(0.15,0.8,"pPb EPOS");

    c1_ks->Print("v2ClosureSystematicsKshort.pdf");

    //Do lambda


    // Pull TGraph for Kshort and lambda

    TGraphErrors* Reco_v2_la = new TGraphErrors(la_npoints,pTReco_la,v2Reco_la,0,v2Reco_laE);
    TGraphErrors* RecoMatch_v2_la = new TGraphErrors(la_npoints,pTRecoMatch_la,v2RecoMatch_la,0,v2RecoMatch_laE);
    TGraphErrors* Gen_v2_la = new TGraphErrors(la_npoints,pTGen_la,v2Gen_la,0,v2Gen_laE);

    Gen_v2_la->SetMarkerColor(kBlue-4);
    Gen_v2_la->SetMarkerStyle(20);
    Gen_v2_la->SetMarkerSize(1.5);
    Gen_v2_la->SetLineColor(kBlue-4);

    Reco_v2_la->SetMarkerColor(kBlue-4);
    Reco_v2_la->SetMarkerStyle(25);
    Reco_v2_la->SetMarkerSize(1.5);
    Reco_v2_la->SetLineColor(kBlue-4);

    RecoMatch_v2_la->SetMarkerColor(kGreen);
    RecoMatch_v2_la->SetMarkerStyle(25);
    RecoMatch_v2_la->SetMarkerSize(1.5);
    RecoMatch_v2_la->SetLineColor(kGreen);


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
    Reco_v2_la->Draw("P");
    RecoMatch_v2_la->Draw("P");

    TLegend* leg_la = new TLegend(0.15,0.65,0.3,0.75);
    leg_la->SetFillColor(10);
    leg_la->SetFillStyle(0);
    leg_la->SetBorderSize(0);
    leg_la->SetTextFont(42);
    leg_la->SetTextSize(0.03);
    leg_la->AddEntry(Gen_v2_la, "Gen #Lambda/#bar{#Lambda}", "P");
    leg_la->AddEntry(Reco_v2_la, "Reco #Lambda/#bar{#Lambda}", "P");
    leg_la->AddEntry(RecoMatch_v2_la, "Reco Match #Lambda/#bar{#Lambda}", "P");
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
    Reco_v2_ks->Draw("P");
    RecoMatch_v2_ks->Draw("P");

    Gen_v2_la->Draw("P");
    Reco_v2_la->Draw("P");
    RecoMatch_v2_la->Draw("P");

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
    leg_co2->AddEntry(Reco_v2_la, "Reco", "P");
    leg_co2->AddEntry(RecoMatch_v2_la, "Reco Match", "P");
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

    std::vector<double> v2ErrorRatioKsReco;
    std::vector<double> v2ErrorRatioLaReco;
    std::vector<double> v2ErrorRatioKsRecoMatch;
    std::vector<double> v2ErrorRatioLaRecoMatch;

    std::vector<double> v2ErrorRatioKsGen;
    std::vector<double> v2ErrorRatioLaGen;

    //Kshort
    for(int i=0; i<ks_npoints; i++)
    {
        v2RatioKsReco.push_back(v2Reco_ks[i]/ksFit->Eval(pTReco_ks[i]));
        v2ErrorRatioKsReco.push_back(v2Reco_ksE[i]/ksFit->Eval(pTReco_ks[i]));

        v2RatioKsRecoMatch.push_back(v2RecoMatch_ks[i]/ksFit->Eval(pTRecoMatch_ks[i]));
        v2ErrorRatioKsRecoMatch.push_back(v2RecoMatch_ksE[i]/ksFit->Eval(pTRecoMatch_ks[i]));

        v2RatioKsGen.push_back(v2Gen_ks[i]/ksFit->Eval(pTReco_ks[i]));
        v2ErrorRatioKsGen.push_back(v2Gen_ksE[i]/ksFit->Eval(pTReco_ks[i]));

        cout << "Ratio " << i << ": " << v2RatioKsReco[i] << endl;
    }

    //Lambda
    for(int i=0; i<la_npoints; i++)
    {
        v2RatioLaReco.push_back(v2Reco_la[i]/laFit->Eval(pTReco_la[i]));
        v2ErrorRatioLaReco.push_back(v2Reco_laE[i]/laFit->Eval(pTReco_la[i]));
        cout << i << endl;

        v2RatioLaRecoMatch.push_back(v2RecoMatch_la[i]/laFit->Eval(pTRecoMatch_la[i]));
        v2ErrorRatioLaRecoMatch.push_back(v2RecoMatch_laE[i]/laFit->Eval(pTRecoMatch_la[i]));
        cout << i << endl;

        v2RatioLaGen.push_back(v2Gen_la[i]/laFit->Eval(pTReco_la[i]));
        v2ErrorRatioLaGen.push_back(v2Gen_laE[i]/laFit->Eval(pTReco_la[i]));
    }

    double* av2RatioKsReco = &v2RatioKsReco[0];
    double* av2RatioLaReco = &v2RatioLaReco[0];
    double* av2RatioKsRecoMatch = &v2RatioKsRecoMatch[0];
    double* av2RatioLaRecoMatch = &v2RatioLaRecoMatch[0];
    double* av2RatioKsGen = &v2RatioKsGen[0];
    double* av2RatioLaGen = &v2RatioLaGen[0];

    double* av2ErrorRatioKsReco = &v2ErrorRatioKsReco[0];
    double* av2ErrorRatioLaReco = &v2ErrorRatioLaReco[0];
    double* av2ErrorRatioKsRecoMatch = &v2ErrorRatioKsRecoMatch[0];
    double* av2ErrorRatioLaRecoMatch = &v2ErrorRatioLaRecoMatch[0];
    double* av2ErrorRatioKsGen = &v2ErrorRatioKsGen[0];
    double* av2ErrorRatioLaGen = &v2ErrorRatioLaGen[0];

    TGraphErrors* RatioReco_v2_la = new TGraphErrors(la_npoints,pTReco_la,av2RatioLaReco,0,av2ErrorRatioLaReco);
    TGraphErrors* RatioRecoMatch_v2_la = new TGraphErrors(la_npoints,pTRecoMatch_la,av2RatioLaRecoMatch,0,av2ErrorRatioLaRecoMatch);
    TGraphErrors* RatioGen_v2_la = new TGraphErrors(la_npoints,pTGen_la,av2RatioLaGen,0,av2ErrorRatioLaGen);


    TGraphErrors* RatioReco_v2_ks = new TGraphErrors(ks_npoints,pTReco_ks,av2RatioKsReco,0,av2ErrorRatioKsReco);
    TGraphErrors* RatioRecoMatch_v2_ks = new TGraphErrors(ks_npoints,pTRecoMatch_ks,av2RatioKsRecoMatch,0,av2ErrorRatioKsRecoMatch);
    TGraphErrors* RatioGen_v2_ks = new TGraphErrors(ks_npoints,pTGen_ks,av2RatioKsGen,0,av2ErrorRatioKsGen);

    RatioReco_v2_ks->SetMarkerColor(kRed);
    RatioReco_v2_ks->SetMarkerStyle(25);
    RatioReco_v2_ks->SetMarkerSize(1.3);
    RatioReco_v2_ks->SetLineColor(kRed);

    RatioRecoMatch_v2_ks->SetMarkerColor(kRed);
    RatioRecoMatch_v2_ks->SetMarkerStyle(26);
    RatioRecoMatch_v2_ks->SetMarkerSize(1.3);
    RatioRecoMatch_v2_ks->SetLineColor(kRed);

    RatioGen_v2_ks->SetMarkerColor(kRed);
    RatioGen_v2_ks->SetMarkerStyle(20);
    RatioGen_v2_ks->SetMarkerSize(1.3);
    RatioGen_v2_ks->SetLineColor(kRed);

    RatioReco_v2_la->SetMarkerColor(kBlue-4);
    RatioReco_v2_la->SetMarkerStyle(25);
    RatioReco_v2_la->SetMarkerSize(1.3);
    RatioReco_v2_la->SetLineColor(kBlue-4);

    RatioRecoMatch_v2_la->SetMarkerColor(kBlue-4);
    RatioRecoMatch_v2_la->SetMarkerStyle(26);
    RatioRecoMatch_v2_la->SetMarkerSize(1.3);
    RatioRecoMatch_v2_la->SetLineColor(kBlue-4);

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
    frame_Ratio->GetYaxis()->SetTitle("Reco/Gen");
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
    RatioReco_v2_la->Draw("P");
    RatioRecoMatch_v2_ks->Draw("P");
    RatioRecoMatch_v2_la->Draw("P");
    RatioGen_v2_ks->Draw("P");
    RatioGen_v2_la->Draw("P");

    TLegend* leg_ratio_species = new TLegend(0.70,0.75,0.80,0.85);
    leg_ratio_species->SetFillColor(10);
    leg_ratio_species->SetFillStyle(0);
    leg_ratio_species->SetBorderSize(0);
    leg_ratio_species->SetTextFont(42);
    leg_ratio_species->SetTextSize(0.03);
    leg_ratio_species->AddEntry(RatioGen_v2_ks, "K_{S}^{0}", "P");
    leg_ratio_species->AddEntry(RatioGen_v2_la, "#Lambda/#bar{#Lambda}", "P");
    leg_ratio_species->Draw();

    TLegend* leg_ratio_label = new TLegend(0.80,0.75,0.9,0.85);
    leg_ratio_label->SetFillColor(10);
    leg_ratio_label->SetFillStyle(0);
    leg_ratio_label->SetBorderSize(0);
    leg_ratio_label->SetTextFont(42);
    leg_ratio_label->SetTextSize(0.03);
    leg_ratio_label->AddEntry(RatioGen_v2_ks, "Gen", "P");
    leg_ratio_label->AddEntry(RatioReco_v2_ks, "Reco", "P");
    leg_ratio_label->Draw();

    tex->SetTextSize(0.03);
    tex->DrawLatex(0.4,0.85,"pPb EPOS");
    tex->DrawLatex(0.4,0.80,"|#eta| > 2");
    tex->DrawLatex(0.4,0.75,"0.3 < p_{T}^{assoc} < 3.0 GeV");

    c1_ratio->Print("v2RatioClosureSystematics.pdf");

    // Write points into rootfile
    TFile out("V0ClosureSys.root","RECREATE");
    Gen_v2_la->Write("GenLa");
    Reco_v2_la->Write("RecoLa");
    Gen_v2_ks->Write("GenKs");
    Reco_v2_ks->Write("RecoKs");
    RatioReco_v2_ks->Write("RatioRecoKs");
    RatioReco_v2_la->Write("RatioRecoLa");
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

    std::vector<double> xcor_ks;
    std::vector<double> xcor_la;

    std::vector<double> ratio_xi_1p65;
    std::vector<double> ratio_xi_1p35;
    std::vector<double> ratio_xi_1p05;
    std::vector<double> ratio_xi_0p75;

    std::vector<double> ratio_om_1p65;
    std::vector<double> ratio_om_1p35;
    std::vector<double> ratio_om_1p05;
    std::vector<double> ratio_om_0p75;

    std::vector<double> xcor_xi;
    std::vector<double> xcor_om;

    double x_2=-99;
    double y_2=-99;
    double x_1p65=-99;
    double y_1p65=-99;
    double x_1p35=-99;
    double y_1p35=-99;
    double x_1p05=-99;
    double y_1p05=-99;
    double x_0p75=-99;
    double y_0p75=-99;

    int Npoints_ks = 0;
    int Npoints_la = 0;
    int Npoints_xi = 0;
    int Npoints_om = 0;
    Npoints_ks = ks8_v2_2->GetN();
    Npoints_la = la8_v2_2->GetN();
    Npoints_xi = xi8_v2_2->GetN();
    Npoints_om = om8_v2_2->GetN();

    for(int i=0; i<Npoints_ks; i++)
    {
        ks8_v2_2->GetPoint(i,x_2,y_2);
        ks8_v2_1p65->GetPoint(i,x_1p65,y_1p65);
        ks8_v2_1p35->GetPoint(i,x_1p35,y_1p35);
        ks8_v2_1p05->GetPoint(i,x_1p05,y_1p05);
        ks8_v2_0p75->GetPoint(i,x_0p75,y_0p75);

        ratio_ks_1p65.push_back(y_1p65/y_2 );
        ratio_ks_1p35.push_back(y_1p35/y_2 );
        ratio_ks_1p05.push_back(y_1p05/y_2 );
        ratio_ks_0p75.push_back(y_0p75/y_2 );

        xcor_ks.push_back(x_2);

    }
    for(int i=0; i<Npoints_la; i++)
    {
        la8_v2_2->GetPoint(i,x_2,y_2);
        la8_v2_1p65->GetPoint(i,x_1p65,y_1p65);
        la8_v2_1p35->GetPoint(i,x_1p35,y_1p35);
        la8_v2_1p05->GetPoint(i,x_1p05,y_1p05);
        la8_v2_0p75->GetPoint(i,x_0p75,y_0p75);

        ratio_la_1p65.push_back(y_1p65/y_2 );
        ratio_la_1p35.push_back(y_1p35/y_2 );
        ratio_la_1p05.push_back(y_1p05/y_2 );
        ratio_la_0p75.push_back(y_0p75/y_2 );

        xcor_la.push_back(x_2);
    }
    for(int i=0; i<Npoints_xi; i++)
    {
        xi8_v2_2->GetPoint(i,x_2,y_2);
        xi8_v2_1p65->GetPoint(i,x_1p65,y_1p65);
        xi8_v2_1p35->GetPoint(i,x_1p35,y_1p35);
        xi8_v2_1p05->GetPoint(i,x_1p05,y_1p05);
        xi8_v2_0p75->GetPoint(i,x_0p75,y_0p75);

        ratio_xi_1p65.push_back(y_1p65/y_2 );
        ratio_xi_1p35.push_back(y_1p35/y_2 );
        ratio_xi_1p05.push_back(y_1p05/y_2 );
        ratio_xi_0p75.push_back(y_0p75/y_2 );

        xcor_xi.push_back(x_2);

    }
    for(int i=0; i<Npoints_om; i++)
    {
        om8_v2_2->GetPoint(i,x_2,y_2);
        om8_v2_1p65->GetPoint(i,x_1p65,y_1p65);
        om8_v2_1p35->GetPoint(i,x_1p35,y_1p35);
        om8_v2_1p05->GetPoint(i,x_1p05,y_1p05);
        om8_v2_0p75->GetPoint(i,x_0p75,y_0p75);

        ratio_om_1p65.push_back(y_1p65/y_2 );
        ratio_om_1p35.push_back(y_1p35/y_2 );
        ratio_om_1p05.push_back(y_1p05/y_2 );
        ratio_om_0p75.push_back(y_0p75/y_2 );

        xcor_om.push_back(x_2);

    }


    TGraphErrors* ks8_v2_1p65_ratio = new TGraphErrors(Npoints_ks,&xcor_ks[0],&ratio_ks_1p65[0],0,0);
    TGraphErrors* ks8_v2_1p35_ratio = new TGraphErrors(Npoints_ks,&xcor_ks[0],&ratio_ks_1p35[0],0,0);
    TGraphErrors* ks8_v2_1p05_ratio = new TGraphErrors(Npoints_ks,&xcor_ks[0],&ratio_ks_1p05[0],0,0);
    TGraphErrors* ks8_v2_0p75_ratio = new TGraphErrors(Npoints_ks,&xcor_ks[0],&ratio_ks_0p75[0],0,0);

    TGraphErrors* la8_v2_1p65_ratio = new TGraphErrors(Npoints_la,&xcor_la[0],&ratio_la_1p65[0],0,0);
    TGraphErrors* la8_v2_1p35_ratio = new TGraphErrors(Npoints_la,&xcor_la[0],&ratio_la_1p35[0],0,0);
    TGraphErrors* la8_v2_1p05_ratio = new TGraphErrors(Npoints_la,&xcor_la[0],&ratio_la_1p05[0],0,0);
    TGraphErrors* la8_v2_0p75_ratio = new TGraphErrors(Npoints_la,&xcor_la[0],&ratio_la_0p75[0],0,0);

    TGraphErrors* xi8_v2_1p65_ratio = new TGraphErrors(Npoints_xi,&xcor_xi[0],&ratio_xi_1p65[0],0,0);
    TGraphErrors* xi8_v2_1p35_ratio = new TGraphErrors(Npoints_xi,&xcor_xi[0],&ratio_xi_1p35[0],0,0);
    TGraphErrors* xi8_v2_1p05_ratio = new TGraphErrors(Npoints_xi,&xcor_xi[0],&ratio_xi_1p05[0],0,0);
    TGraphErrors* xi8_v2_0p75_ratio = new TGraphErrors(Npoints_xi,&xcor_xi[0],&ratio_xi_0p75[0],0,0);

    TGraphErrors* om8_v2_1p65_ratio = new TGraphErrors(Npoints_om,&xcor_om[0],&ratio_om_1p65[0],0,0);
    TGraphErrors* om8_v2_1p35_ratio = new TGraphErrors(Npoints_om,&xcor_om[0],&ratio_om_1p35[0],0,0);
    TGraphErrors* om8_v2_1p05_ratio = new TGraphErrors(Npoints_om,&xcor_om[0],&ratio_om_1p05[0],0,0);
    TGraphErrors* om8_v2_0p75_ratio = new TGraphErrors(Npoints_om,&xcor_om[0],&ratio_om_0p75[0],0,0);

    ks8_v2_2->SetMarkerColor(kRed+4);
    ks8_v2_2->SetMarkerStyle(20);
    ks8_v2_2->SetMarkerSize(1.5);
    ks8_v2_2->SetLineColor(kRed+4);

    ks8_v2_1p65->SetMarkerColor(kRed+3);
    ks8_v2_1p65->SetMarkerStyle(24);
    ks8_v2_1p65->SetMarkerSize(1.5);
    ks8_v2_1p65->SetLineColor(kRed+3);

    ks8_v2_1p35->SetMarkerColor(kRed+2);
    ks8_v2_1p35->SetMarkerStyle(24);
    ks8_v2_1p35->SetMarkerSize(1.5);
    ks8_v2_1p35->SetLineColor(kRed+2);

    ks8_v2_1p05->SetMarkerColor(kRed+1);
    ks8_v2_1p05->SetMarkerStyle(24);
    ks8_v2_1p05->SetMarkerSize(1.5);
    ks8_v2_1p05->SetLineColor(kRed+1);

    ks8_v2_0p75->SetMarkerColor(kRed);
    ks8_v2_0p75->SetMarkerStyle(24);
    ks8_v2_0p75->SetMarkerSize(1.5);
    ks8_v2_0p75->SetLineColor(kRed);

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
    //xi8_v2->SetMarkerColor(kGreen+2);
    //xi8_v2->SetMarkerStyle(21);
    //xi8_v2->SetMarkerSize(1.5);
    //xi8_v2->SetLineColor(kGreen+2);

    la8_v2_2->SetMarkerColor(kBlue+4);
    la8_v2_2->SetMarkerStyle(22);
    la8_v2_2->SetMarkerSize(1.5);
    la8_v2_2->SetLineColor(kBlue+4);

    la8_v2_1p65->SetMarkerColor(kBlue+3);
    la8_v2_1p65->SetMarkerStyle(26);
    la8_v2_1p65->SetMarkerSize(1.5);
    la8_v2_1p65->SetLineColor(kBlue+3);

    la8_v2_1p35->SetMarkerColor(kBlue+2);
    la8_v2_1p35->SetMarkerStyle(26);
    la8_v2_1p35->SetMarkerSize(1.5);
    la8_v2_1p35->SetLineColor(kBlue+2);

    la8_v2_1p05->SetMarkerColor(kBlue+1);
    la8_v2_1p05->SetMarkerStyle(26);
    la8_v2_1p05->SetMarkerSize(1.5);
    la8_v2_1p05->SetLineColor(kBlue+1);

    la8_v2_0p75->SetMarkerColor(kBlue);
    la8_v2_0p75->SetMarkerStyle(26);
    la8_v2_0p75->SetMarkerSize(1.5);
    la8_v2_0p75->SetLineColor(kBlue);

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

    xi8_v2_2->SetMarkerColor(kRed+4);
    xi8_v2_2->SetMarkerStyle(20);
    xi8_v2_2->SetMarkerSize(1.5);
    xi8_v2_2->SetLineColor(kRed+4);

    xi8_v2_1p65->SetMarkerColor(kRed+3);
    xi8_v2_1p65->SetMarkerStyle(24);
    xi8_v2_1p65->SetMarkerSize(1.5);
    xi8_v2_1p65->SetLineColor(kRed+3);

    xi8_v2_1p35->SetMarkerColor(kRed+2);
    xi8_v2_1p35->SetMarkerStyle(24);
    xi8_v2_1p35->SetMarkerSize(1.5);
    xi8_v2_1p35->SetLineColor(kRed+2);

    xi8_v2_1p05->SetMarkerColor(kRed+1);
    xi8_v2_1p05->SetMarkerStyle(24);
    xi8_v2_1p05->SetMarkerSize(1.5);
    xi8_v2_1p05->SetLineColor(kRed+1);

    xi8_v2_0p75->SetMarkerColor(kRed);
    xi8_v2_0p75->SetMarkerStyle(24);
    xi8_v2_0p75->SetMarkerSize(1.5);
    xi8_v2_0p75->SetLineColor(kRed);

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

    om8_v2_2->SetMarkerColor(kRed+4);
    om8_v2_2->SetMarkerStyle(20);
    om8_v2_2->SetMarkerSize(1.5);
    om8_v2_2->SetLineColor(kRed+4);

    om8_v2_1p65->SetMarkerColor(kRed+3);
    om8_v2_1p65->SetMarkerStyle(24);
    om8_v2_1p65->SetMarkerSize(1.5);
    om8_v2_1p65->SetLineColor(kRed+3);

    om8_v2_1p35->SetMarkerColor(kRed+2);
    om8_v2_1p35->SetMarkerStyle(24);
    om8_v2_1p35->SetMarkerSize(1.5);
    om8_v2_1p35->SetLineColor(kRed+2);

    om8_v2_1p05->SetMarkerColor(kRed+1);
    om8_v2_1p05->SetMarkerStyle(24);
    om8_v2_1p05->SetMarkerSize(1.5);
    om8_v2_1p05->SetLineColor(kRed+1);

    om8_v2_0p75->SetMarkerColor(kRed);
    om8_v2_0p75->SetMarkerStyle(24);
    om8_v2_0p75->SetMarkerSize(1.5);
    om8_v2_0p75->SetLineColor(kRed);

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

    TLegend* leg_om = new TLegend(0.25,0.55,0.4,0.85);
    leg_om->SetFillColor(10);
    leg_om->SetFillStyle(0);
    leg_om->SetBorderSize(0);
    leg_om->SetTextFont(42);
    leg_om->SetTextSize(0.05);
    leg_om->AddEntry(om8_v2_2, "|#Delta#eta| > 2.0", "P");
    leg_om->AddEntry(om8_v2_1p65, "|#Delta#eta| > 1.65", "P");
    leg_om->AddEntry(om8_v2_1p35, "|#Delta#eta| > 1.35", "P");
    leg_om->AddEntry(om8_v2_1p05, "|#Delta#eta| > 1.05", "P");
    leg_om->AddEntry(om8_v2_0p75, "|#Delta#eta| > 0.75", "P");
    //leg_om->AddEntry(xi8_v2, "#Xi^{+}/ #Xi^{-}", "P");

    TLegend* leg_ks_ratio = new TLegend(0.25,0.55,0.4,0.85);
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

    TLegend* leg_la_ratio = new TLegend(0.25,0.55,0.4,0.85);
    leg_la_ratio->SetFillColor(10);
    leg_la_ratio->SetFillStyle(0);
    leg_la_ratio->SetBorderSize(0);
    leg_la_ratio->SetTextFont(42);
    leg_la_ratio->SetTextSize(0.05);
    leg_la_ratio->AddEntry(la8_v2_1p65_ratio, "|#Delta#eta| > 1.65", "P");
    leg_la_ratio->AddEntry(la8_v2_1p35_ratio, "|#Delta#eta| > 1.35", "P");
    leg_la_ratio->AddEntry(la8_v2_1p05_ratio, "|#Delta#eta| > 1.05", "P");
    leg_la_ratio->AddEntry(la8_v2_0p75_ratio, "|#Delta#eta| > 0.75", "P");

    TLegend* leg_xi_ratio = new TLegend(0.25,0.55,0.4,0.85);
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

    TLegend* leg_om_ratio = new TLegend(0.65,0.55,0.77,0.85);
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
    om8_v2_2->Draw("P");
    om8_v2_1p65->Draw("P");
    om8_v2_1p35->Draw("P");
    om8_v2_1p05->Draw("P");
    om8_v2_0p75->Draw("P");
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
