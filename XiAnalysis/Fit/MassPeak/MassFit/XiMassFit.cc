//Includes
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <TROOT.h>
#include <TStyle.h>

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooChi2Var.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "TString.h"
#include "TGaxis.h"

#include <vector>

void XiMassFit()
{
    //Initializers
    using namespace RooFit;
    gStyle->SetMarkerSize(0.5);
    TGaxis::SetMaxDigits(3);
    RooMsgService::instance().setStreamStatus(0,kFALSE);
    RooMsgService::instance().setStreamStatus(1,kFALSE);
    std::ostringstream os;
    std::ostringstream osYield;
    std::ofstream myfile;

    TH1D* massxi;
    TH2D* MassXi;
    std::vector<RooPlot*> xframe;
    std::vector<double> mass_xi;
    std::vector<double> std_xi;
    std::vector<double> fsig_xi;
    std::vector<double> covQual_xi;
    bool doPbPb = true;
    bool doMB = false;
    bool pPb = false;

    //std::vector<double> pxi = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,72, 73,85, 86,100, 101,200, 201,300};
    //std::vector<double> pxi = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,100, 101,200};//, 201,300};
    std::vector<double> pxi = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,72, 73,100}; //pPb
    //std::vector<double> pxi = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,72}; //PbPb

    TCanvas* cc1 = new TCanvas("cc1","cc1",1200,1200);
    cc1->Divide(3,3);

    TCanvas* c_xi = new TCanvas("c_xi","c_xi",800,600);

    //File Creation
    myfile.open("XiPeakParam.txt");
    TFile* file;
    //file = new TFile("/Volumes/MacHD/Users/blt1/research/CascadeV2pPb/RootFiles/Flow/CasCutLoose/CasCutLooseJL40.root");
    //file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/Thesis/XiAnalysisCorrelationPtCut8TeVPD1_4_ForFinal.root");
    //file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/Composites/V0CasMassPtPD11_16.root");
    //file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/PeripheralSubtractionMB.root"); //MB
    //if(!doPbPb) file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationRapidityTotal_08_20_2017.root"); //pPb old
    if(doPbPb) file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/V0CasCorrelationPbPbTotal_10_30_17.root"); //PbPb
    else if(doMB) file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MBCorr/XiMB_0_35_1_02_18.root");
    else if(pPb) file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationHM_11_07_17.root"); //pPb

    //file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/Composites/V0CasMassPtPD5JL12.root"); //only one PD

    //MassXi = (TH2D*)file->Get("XiMassPt/MassPt");
    if(doMB) MassXi = (TH2D*)file->Get("v0CasCorrelationRapidityPeriSub/MassPtXi"); //MB
    else if(doPbPb) MassXi = (TH2D*)file->Get("v0CasCorrelationRapidityPbPb/MassPtXi"); //PbPb
    else if(pPb)   MassXi = (TH2D*)file->Get("v0CasCorrelationRapidity/MassPtXi"); //pPb
    //MassXi = (TH2D*)file->Get("t/MassPt");

    //Fit
    int pxicounter = 0; //for correct bin counting
    int hbincounter =1; //histogram bin counting
    //for(unsigned i=0; i<1; i++)
    for(unsigned i=0; i<pxi.size(); i++)
    {
        TCanvas* cc2 = new TCanvas("cc2","",600,450);
        int index = (i+2)/2;
        massxi = (TH1D*)MassXi->ProjectionX("massxi", pxi[i],pxi[i+1]);
        massxi->GetXaxis()->SetRangeUser(1.25,1.40);

        gStyle->SetOptTitle(kFALSE);

        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextFont(62);
        tex->SetTextSize(0.04);
        //tex->SetTextAlign(10);

        RooRealVar x("x","mass",1.25,1.40);
        RooPlot* xframe_ = x.frame(150);
        xframe_->GetXaxis()->SetTitle("Invariant mass (GeV)");
        xframe_->GetYaxis()->SetTitle("Candidates / 1 MeV");
        xframe_->GetXaxis()->CenterTitle(1);
        xframe_->GetYaxis()->CenterTitle(1);
        xframe_->GetXaxis()->SetTickSize(0.02);
        xframe_->GetYaxis()->SetTickSize(0.02);
        xframe_->GetXaxis()->SetNdivisions(407);
        xframe_->GetYaxis()->SetNdivisions(410);
        xframe_->GetXaxis()->SetTitleSize(0.06);
        xframe_->GetYaxis()->SetTitleSize(0.06);
        xframe_->GetYaxis()->SetTitleOffset(1.05);
        xframe_->GetXaxis()->SetTitleOffset(0.5);
        //xframe_->GetXaxis()->SetLabelSize(xframe_->GetXaxis()->GetLabelSize()*2.0);
        xframe_->GetYaxis()->SetLabelSize(0.1);
        xframe_->GetXaxis()->SetLabelSize(0.1);
        //xframe_->GetYaxis()->SetLabelSize(xframe_->GetYaxis()->GetLabelSize()*2.0);
        RooDataHist data("data","dataset",x,massxi);
        data.plotOn(xframe_,Name("data"));
        RooRealVar mean("mean","mean",1.32,1.29,1.33);
        RooRealVar sigma1("sigma1","sigma1",0.004,0.001,0.04);
        RooRealVar sigma2("sigma2","sigma2",0.006,0.001,0.04);
        RooRealVar sig1("sig1","signal1",2200,0,1000000000);
        RooRealVar sig2("sig2","signal2",1500,0,1000000000);
        RooRealVar qsig("qsig","qsig",7000,0,1000000000);
        RooRealVar alpha("alpha","alpha",1,0,10);
        RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
        RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
        RooRealVar ap("ap","ap",0.1,-1,1);
        RooRealVar bp("bp","bp",-0.1,-1,1);
        RooRealVar cp("cp","cp",-0.1,-1,1);
        RooRealVar dp("dp","dp",-0.1,-1,1);
        RooChebychev background("background","background",x,RooArgList(ap,bp,cp,dp));
        //RooGenericPdf background("background", "x - (1.115683 + 0.13957018)^alpha", RooArgList(x,alpha));
        RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,background),RooArgList(sig1,sig2,qsig));

        //x.setRange("cut",1.285,1.375);
        x.setRange("cut",1.26,1.39);

        RooFitResult* r_xi = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));
        //RooChi2Var chi2_xiVar("chi2_xi","chi2",sum,data);

        double covQual = r_xi->covQual();
        double mean_xi = mean.getVal();

        double gaus1F_xi = sig1.getVal();
        double gaus2F_xi = sig2.getVal();
        double qsig_xi   = qsig.getVal();
        double Norm = (gaus2F_xi + gaus1F_xi + qsig_xi)*0.001;

        //set ranges for individual gaussian yield determination
        x.setRange("g1", mean.getVal() - 2*sigma1.getVal(), mean.getVal() + 2*sigma1.getVal());
        x.setRange("g2", mean.getVal() - 2*sigma2.getVal(), mean.getVal() + 2*sigma2.getVal());

        RooAbsReal* Intgaus1_yield_xi = gaus1.createIntegral(x,x,"g1");
        RooAbsReal* Intgaus2_yield_xi = gaus2.createIntegral(x,x,"g2");

        double gaus1_yield_xi = gaus1F_xi*Intgaus1_yield_xi->getVal();
        double gaus2_yield_xi = gaus2F_xi*Intgaus2_yield_xi->getVal();
        double gausTot_yield_xi = gaus1_yield_xi + gaus2_yield_xi;

        cout << "Yield1: " << gaus1_yield_xi << endl;
        cout << "Yield2: " << gaus2_yield_xi << endl;

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

        mass_xi.push_back(mean_xi);
        std_xi.push_back(rms_true_xi);
        fsig_xi.push_back(Fsig_xi);
        covQual_xi.push_back(covQual);

        cout << "Yield (xi): " << Yield_xi << endl;
        cout << "Fsig (xi): " << Fsig_xi << endl;
        cout << "std (xi): "  << rms_true_xi  << endl;
        cout << "mass (xi): " << mean_xi << endl;

        cout << "covQual (xi)" << covQual << endl;

        sum.plotOn(xframe_,Name("sum"),NormRange("cut"),LineWidth(2),LineColor(kBlue));
        sum.plotOn(xframe_,Components(background),NormRange("cut"),LineStyle(kDashed),LineWidth(2),LineColor(kBlue));

        RooArgSet s(ap,bp,cp,dp,mean,qsig,sig1,sig2,sigma1);
        s.add(sigma2);
        TF1* func = (TF1*)sum.asTF(RooArgList(x),RooArgList(s),x);
        std::vector<double> pull_xi;
        std::vector<double> pull_xi_x;
        std::vector<double> ratio_xi;
        for(int j=1; j<131; j++)
        {
            double data_point = massxi->GetBinContent(10+j);
            double data_error = massxi->GetBinError(10+j);
            double fit_point = Norm*func->Eval(massxi->GetBinCenter(10+j));
            if(data_point == 0) data_error = 1;
            pull_xi.push_back((data_point - fit_point)/data_error);
            pull_xi_x.push_back(massxi->GetBinCenter(10+j));
            ratio_xi.push_back(data_point/fit_point);
        }

        TGraphErrors* TGpull = new TGraphErrors(130,&pull_xi_x[0],&pull_xi[0],0,0);
        TGraphErrors* TGratio= new TGraphErrors(130,&pull_xi_x[0],&ratio_xi[0],0,0);

        double chisquare = 0;

        for(int j=0; j<pull_xi.size(); j++)
        {
            chisquare+=TMath::Power(pull_xi[j],2);
        }

        int ndf = 130 - (10 - 1);

        cc1->cd(index);
        gPad->SetBottomMargin(0.15); //gives more space for titles
        gPad->SetLeftMargin(0.15);
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        xframe.push_back(xframe_);
        xframe_->Draw();
        cc1->Update();
        if(i==8)
        {
            RooPlot* xframe_sxi = x.frame(150);
            xframe_sxi->GetXaxis()->SetTitle("m_{#Lambda#pi} (GeV/c^{2})");
            xframe_sxi->GetYaxis()->SetTitle("Candidates / 1 MeV");
            xframe_sxi->GetXaxis()->CenterTitle(1);
            xframe_sxi->GetYaxis()->CenterTitle(1);
            xframe_sxi->GetXaxis()->SetTickSize(0.02);
            xframe_sxi->GetYaxis()->SetTickSize(0.02);
            xframe_sxi->GetXaxis()->SetNdivisions(407);
            xframe_sxi->GetYaxis()->SetNdivisions(410);
            xframe_sxi->GetXaxis()->SetTitleSize(0.05);
            xframe_sxi->GetYaxis()->SetTitleSize(0.05);
            xframe_sxi->GetYaxis()->SetTitleOffset(0.9);
            xframe_sxi->GetXaxis()->SetTitleOffset(0.9);
            //xframe_sxi->GetXaxis()->SetLabelSize(xframe_sxi->GetXaxis()->GetLabelSize()*2.0);
            xframe_sxi->GetYaxis()->SetLabelSize(0.04);
            xframe_sxi->GetXaxis()->SetLabelSize(0.04);
            data.plotOn(xframe_sxi,Name("data"));
            sum.plotOn(xframe_sxi,Name("sum"),NormRange("cut"),LineWidth(2),LineColor(kBlue));
            sum.plotOn(xframe_sxi,Name("poly"),Components(background),NormRange("cut"),LineStyle(kDashed),LineWidth(2),LineColor(kBlue));
            c_xi->cd();
            gPad->SetTickx();
            gPad->SetTicky();
            xframe_sxi->Draw();
            c_xi->Update();
            TLine* t1_xi = new TLine(mean.getVal() - 2*rms_true_xi, 0, mean.getVal() - 2*rms_true_xi, gPad->GetUymax());
            TLine* t2_xi = new TLine(mean.getVal() + 2*rms_true_xi, 0, mean.getVal() + 2*rms_true_xi, gPad->GetUymax());
            t1_xi->SetLineStyle(2);
            t1_xi->SetLineColor(kGreen+1);
            t2_xi->SetLineStyle(2);
            t2_xi->SetLineColor(kGreen+1);
            t1_xi->Draw();
            t2_xi->Draw();

            double xstart_xi = 0.13;
            double ystart_xi = 0.82;
            double xpos = xstart_xi;
            double ypos = ystart_xi;
            double increment = 0.07;
            tex->SetTextSize(0.05);
            os << "CMS pPb #sqrt{S_{NN}} = 8.16 TeV";
            tex->DrawLatex(0.513,0.918,os.str().c_str());
            os.str(std::string());
            tex->SetTextSize(0.045);
            os << "185 #leq N_{trk}^{offline} < 250 ";
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            os << (pxi[i]-1)/10 << " < p_{T} < " << pxi[i+1]/10  << " GeV/c";
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            os << "|y| < 1";
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());

            TLegend* leg = new TLegend(0.58,0.61,0.75,0.8);
            leg->SetFillColor(10);
            leg->SetFillStyle(0);
            leg->SetBorderSize(0);
            leg->SetTextFont(62);
            leg->SetTextSize(0.045);
            leg->AddEntry(xframe_sxi->findObject("data"),"Data","P");
            leg->AddEntry(xframe_sxi->findObject("sum"),"Fit","l");
            leg->AddEntry("poly","Background","l");
            leg->AddEntry(t1_xi,"#pm 2#sigma","l");
            leg->Draw();

            c_xi->Print("Image/XiMassFit/XiPlotForZhenyu.pdf");
        }

        cc1->cd(index);
        double chi2_xi = xframe_->chiSquare("sum","data",4);

        TLine* t1 = new TLine(mean.getVal() - 2*rms_true_xi, 0, mean.getVal() - 2*rms_true_xi, gPad->GetUymax());
        TLine* t2 = new TLine(mean.getVal() + 2*rms_true_xi, 0, mean.getVal() + 2*rms_true_xi, gPad->GetUymax());
        t1->SetLineStyle(2);
        t1->SetLineColor(kGreen);
        t2->SetLineStyle(2);
        t2->SetLineColor(kGreen);
        t1->Draw("same");
        t2->Draw("same");
        double xpos = 0.55;
        double ypos = 0.85;
        double increment = 0.07;
        if(i==0)
        {
            if(!doPbPb)
            {
                os << "CMS pPb";
                tex->SetTextSize(0.06);
                tex->DrawLatex(0.17,ypos-increment,os.str().c_str());
                tex->SetTextSize(0.04);
                os.str(std::string());
                os << "185 #leq N_{trk}^{offline} < 250";
                tex->DrawLatex(0.17,ypos-2*increment,os.str().c_str());
                os.str(std::string());
            }
            else
            {
                os << "CMS PbPb";
                tex->SetTextSize(0.06);
                tex->DrawLatex(0.17,ypos-increment,os.str().c_str());
                tex->SetTextSize(0.04);
                os.str(std::string());
                os << "Centrality 30 - 50%";
                tex->DrawLatex(0.17,ypos-2*increment,os.str().c_str());
                os.str(std::string());
            }
        }
        os  << (pxi[i]-1)/10 << " < P_{t} < "  << pxi[i+1]/10 << " GeV";
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
        os << "Fsig: " << Fsig_xi;
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        //os << "#chi^{2}/ndf: " << std::setprecision(3) << chi2_xi << std::setprecision(6);
        //tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        //os.str(std::string());
        osYield << "Yield: " << std::setprecision(2) << Yield_xi;
        tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
        osYield.str(std::string());
        os.str(std::string());


        //tex->SetTextSize(tex->GetTextSize()*0.95);

        cc2->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        xframe_->GetXaxis()->SetTitleOffset(1);
        xframe_->GetXaxis()->SetTitleSize(xframe_->GetXaxis()->GetTitleSize()*0.8);
        //xframe_->GetYaxis()->SetTitleSize(xframe_->GetYaxis()->GetTitleSize()*1.3);
        xframe_->GetXaxis()->SetLabelSize(xframe_->GetXaxis()->GetLabelSize()*0.5);
        xframe_->GetYaxis()->SetLabelSize(xframe_->GetYaxis()->GetLabelSize()*0.5);
        TPad* pad1 = new TPad("pad1","top pad",0.0,0.4,1.0,1.0);
        TPad* pad2 = new TPad("pad3","middle pad",0.0,0.0,1.0,0.25);
        TPad* pad3 = new TPad("pad2","bottom pad",0.0,0.25,1.0,0.4);
        pad1->SetTopMargin(0.1);
        pad1->SetBottomMargin(0.0);
        pad1->SetRightMargin(0.038);
        pad1->SetLeftMargin(0.17);
        pad1->Draw();

        pad2->SetFrameFillStyle(4000);
        pad2->SetTopMargin(0.0);
        pad2->SetBottomMargin(0.37);
        pad2->SetRightMargin(0.038);
        pad2->SetLeftMargin(0.17);
        pad2->Draw();

        pad3->SetFrameFillStyle(4000);
        pad3->SetTopMargin(0.0);
        pad3->SetBottomMargin(0.0);
        pad3->SetRightMargin(0.038);
        pad3->SetLeftMargin(0.17);
        pad3->Draw();
        pad1->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        xframe_->Draw();
        cc2->Update();

        t1->Draw("same");
        t2->Draw("same");
        xpos = 0.72;
        ypos = 0.85;
        if(i==0)
        {
            if(!doPbPb)
            {
                os << "CMS pPb";
                tex->SetTextSize(0.06);
                tex->DrawLatex(0.25,0.78,os.str().c_str());
                tex->SetTextSize(0.04);
                os.str(std::string());
                os << "185 #leq N_{trk}^{offline} < 250";
                tex->DrawLatex(0.25,0.71,os.str().c_str());
                os.str(std::string());
            }
            else
            {
                os << "CMS PbPb";
                tex->SetTextSize(0.06);
                tex->DrawLatex(0.25,0.78,os.str().c_str());
                tex->SetTextSize(0.04);
                os.str(std::string());
                os << "Centrality 30 - 50%";
                tex->DrawLatex(0.25,0.71,os.str().c_str());
                os.str(std::string());
            }
        }
        os << "P_{t} Bin: " << (pxi[i]-1)/10 << " - " << pxi[i+1]/10;
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "Mean: " << mean_xi;
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "#sigma :" << rms_true_xi;
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "CovQual: " << covQual;
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        osYield << "Yield: " << std::setprecision(2) << Yield_xi;
        tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
        osYield.str(std::string());
        os << "#chi^{2}/ndf: " << chisquare << "/" << ndf;
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());

        pad2->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        TH1F* frame = pad2->cd()->DrawFrame(1.25,-6,1.40,6);
        frame->GetXaxis()->SetTitle("#Lambda #pi Invariant mass (GeV)");
        frame->GetYaxis()->SetTitle("Pull");
        frame->GetXaxis()->CenterTitle(1);
        frame->GetYaxis()->CenterTitle(1);
        frame->GetXaxis()->SetTickSize(0.02);
        frame->GetYaxis()->SetTickSize(0.02);
        frame->GetXaxis()->SetNdivisions(410);
        frame->GetYaxis()->SetNdivisions(407);
        frame->GetXaxis()->SetTitleSize(0.12);
        frame->GetYaxis()->SetTitleSize(0.12);
        frame->GetYaxis()->SetTitleOffset(0.5);
        frame->GetXaxis()->SetTitleOffset(1.2);
        frame->GetYaxis()->SetLabelSize(0.1);
        frame->GetXaxis()->SetLabelSize(0.1);
        TGpull->SetMarkerStyle(20);
        TGpull->Draw("P");
        TLine* line = new TLine(1.25, 0, 1.40, 0);
        line->SetLineStyle(2);
        line->Draw("same");

        pad3->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        double low = 0.1;
        double high = 2.75;
        //if(doPbPb)
        //{
            //low = 0.5;
            //high= 1.5;
        //}
        frame = pad3->cd()->DrawFrame(1.25,low,1.40,high);
        frame->GetYaxis()->SetTitle("Ratio data/fit");
        frame->GetXaxis()->CenterTitle(1);
        frame->GetYaxis()->CenterTitle(1);
        frame->GetXaxis()->SetTickSize(0.02);
        frame->GetYaxis()->SetTickSize(0.02);
        frame->GetXaxis()->SetNdivisions(410);
        frame->GetYaxis()->SetNdivisions(407);
        frame->GetXaxis()->SetTitleSize(0.12);
        frame->GetYaxis()->SetTitleSize(0.12);
        frame->GetYaxis()->SetTitleOffset(0.5);
        frame->GetXaxis()->SetTitleOffset(1.2);
        frame->GetYaxis()->SetLabelSize(0.1);
        frame->GetXaxis()->SetLabelSize(0.1);
        TGratio->SetMarkerStyle(20);
        TGratio->Draw("P");
        TLine* line2 = new TLine(1.25, 1, 1.40, 1);
        line2->SetLineStyle(2);
        line2->Draw("same");

        hbincounter++;
        if(!doPbPb) cc2->Print(Form("Image/XiMassFit/XiMassFit_Pull_pPb%d.pdf",i/2));
        else if(doMB) cc2->Print(Form("Image/XiMassFit/XiMassFit_Pull_MB%d.pdf",i/2));
        else cc2->Print(Form("Image/XiMassFit/XiMassFit_Pull_PbPb%d.pdf",i/2));
        //if(i==0) cc2->Print("XiMassFitInd.pdf(","pdf");
        //else if(i < pxi.size() - 2) cc2->Print("XiMassFitInd.pdf","pdf");
        //else cc2->Print("XiMassFitInd.pdf)","pdf");
        i++; //to access correct bins
    }
    cc1->Print("Image/XiMassFit/XiMassFitComposite.pdf");
    cc1->Print("Image/XiMassFit/XiMassFitComposite.png");

    //Output
    pxicounter = 0;
    for(unsigned i=0; i<mass_xi.size(); i++)
    {
        cout <<  "====================" << endl;
        cout << "Pt Bin: " << (pxi[pxicounter]-1)/10 << " - " << pxi[pxicounter+1]/10 << endl;
        cout <<  "====================" << endl;
        cout << "Mass_xi: "   << mass_xi[i]    << endl;
        cout << "Fsig_xi: "   << fsig_xi[i]    << endl;
        cout << "std_xi: "    << std_xi[i]     << endl;
        cout << "covQual_xi " << covQual_xi[i] << endl;

        myfile <<  "====================" << "\n";
        myfile << "Pt Bin: " << (pxi[pxicounter]-1)/10 << " - " << pxi[pxicounter+1]/10 << "\n";
        myfile <<  "====================" << "\n";
        myfile << "Mass_xi: "   << mass_xi[i]    << "\n";
        myfile << "Fsig_xi: "   << fsig_xi[i]    << "\n";
        myfile << "std_xi: "    << std_xi[i]     << "\n";
        myfile << "covQual_xi " << covQual_xi[i] << "\n";

        pxicounter+=2;
    }
}
