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
#include "TH3.h"
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
#include "RooChi2Var.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "TString.h"
#include "TGaxis.h"

#include <vector>

void V0MassFitFromHist()
{
    //Initializers
    using namespace RooFit;
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetMarkerSize(0.5);
    TGaxis::SetMaxDigits(3);
    RooMsgService::instance().setStreamStatus(0,kFALSE);
    RooMsgService::instance().setStreamStatus(1,kFALSE);
    std::ostringstream os;
    std::ostringstream osYield;
    std::ofstream myfile;

    TH1D* massks[18];
    TH1D* massla[18];
    TH2D* MassKs;
    TH2D* MassLa;
    TH3D* MassKs3D;
    TH3D* MassLa3D;
    bool doPbPb = false;
    bool pPb = false;
    bool mc = false;
    bool mb = true;
    bool lambda;

    std::vector<RooPlot*> Xframe_Ks;
    std::vector<double> mass_ks;
    std::vector<double> std_ks;
    std::vector<double> fsig_ks;
    std::vector<double> covQual_ks;

    std::vector<RooPlot*> Xframe_La;
    std::vector<double> mass_la;
    std::vector<double> std_la;
    std::vector<double> fsig_la;
    std::vector<double> covQual_la;

    //int pTksLength = 26; // the number of bins to be fitted is half of this number
    //std::vector<double> pks = {3,4, 5,6, 7,8, 9,10, 11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,70, 71,85, 86,100};//, 201,250, 251,300};
    //std::vector<double> pla = {0,0, 0,0, 0,0, 9,10, 11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,70, 71,85, 86,100};//, 201,250, 251,300};
    std::vector<double> pks = {0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0};//pPb
    std::vector<double> pla = {0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0}; //pPb
    //std::vector<double> pks = {0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5};//PbPb
    //std::vector<double> pla = {0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5}; //PbPb

    TCanvas* Composite_Ks = new TCanvas("Composite_Ks","",1000,1200);
    Composite_Ks->Divide(3,5);

    TCanvas* Composite_La = new TCanvas("Composite_La","",1000,1200);
    Composite_La->Divide(3,4);


    TCanvas* c_ks = new TCanvas("c_ks","c_ks",800,600);
    TCanvas* c_la = new TCanvas("c_la","c_la",800,600);


    //File Creation
    myfile.open("V0PeakParam.txt");
    TFile* file = 0;
    //Get root files with mass distributions
    if(doPbPb) file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/V0CasCorrelationPbPbTotal_10_30_17.root");
    else if(mc) file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/All/MCMassPtTotal_08_23_2017.root");
    else if(mb) file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MBCorr/V0MB_0_35_1_02_18.root");
    else if(pPb) file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityCorrectMultB_09_19_17.root");

    if(mc)
    {
        MassKs3D = (TH3D*)file->Get("MassPtRapidityMC/KsMassPtRap");
        MassLa3D = (TH3D*)file->Get("MassPtRapidityMC/LaMassPtRap");
        MassKs = (TH2D*)MassKs3D->Project3D("yx");
        MassLa = (TH2D*)MassLa3D->Project3D("yx");
    }



    //Do Fit and calculate values
    int pkscounter  = 0; //for correct bin counting
    int placounter  = 0;
    int hbincounter = 1;
    int Rebin = 1;
    for(unsigned i=0; i<14; i++)
    //for(unsigned i=0; i<1; i++)
    {
        TCanvas* cc1 = new TCanvas("cc1","cc1",600,600);
        TCanvas* cc2 = new TCanvas("cc2","cc2",600,600);
        //i = pks.size()-2;
        if(i >= pla.size()-1) lambda = false;
        else lambda= true;
        if(doPbPb)
        {
            massks[i] = (TH1D*)file->Get(Form("v0CasCorrelationRapidityPbPb/masskshort_pt%d",i));
            massla[i] = (TH1D*)file->Get(Form("v0CasCorrelationRapidityPbPb/masslambda_pt%d",i));
        }
        else if(mc)
        {
            massks[i] = (TH1D*)MassKs->ProjectionX(Form("massks%d",i),pks[i]*10,pks[i+1]*10);
            massla[i] = (TH1D*)MassLa->ProjectionX(Form("massla%d",i),pla[i]*10,pla[i+1]*10);
        }
        else if(mb)
        {
            massks[i]= (TH1D*)file->Get(Form("v0CasCorrelationRapidityPeriSub/masskshort_pt%d",i)); //mb
            massla[i]= (TH1D*)file->Get(Form("v0CasCorrelationRapidityPeriSub/masslambda_pt%d",i)); //mb
        }
        else if(pPb)
        {
            massks[i] = (TH1D*)file->Get(Form("v0CorrelationRapidity/masskshort_pt%d",i));
            massla[i] = (TH1D*)file->Get(Form("v0CorrelationRapidity/masslambda_pt%d",i));
        }
        massks[i]->Rebin(1);
        massla[i]->Rebin(1);
        massks[i]->GetXaxis()->SetRangeUser(0.43,0.565);
        massla[i]->GetXaxis()->SetRangeUser(1.08,1.16);


        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextFont(42);
        tex->SetTextSize(0.05);

        //tex->SetTextSize(tex->GetTextSize()*0.95);

        //kshort
        RooRealVar x("x","mass",0.43,0.565);
        RooDataHist data("data","dataset",x,massks[i]);
        RooRealVar mean("mean","mean",0.5,0.49,0.51);
        RooRealVar sigma1("sigma1","sigma1",0.005,0.001,0.01);
        RooRealVar sigma2("sigma2","sigma2",0.005,0.0001,0.01);
        RooRealVar sig1("sig1","signal1",4e4,0,100000000);
        RooRealVar sig2("sig2","signal2",5e4,0,100000000);
        //RooRealVar sig1("sig1","signal1",1e4,0,100000000);
        //RooRealVar sig2("sig2","signal2",1e4,0,100000000);
        RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
        RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
        RooRealVar ap("ap","ap",0,-5,5);
        RooRealVar bp("bp","bp",0.1,-5,5);
        RooRealVar cp("cp","cp",-0.1,-5,5);
        RooRealVar dp("dp","dp",0.1,-5,5);
        //RooRealVar ap("ap","ap",1,-100,100);
        //RooRealVar bp("bp","bp",1,-100,100);
        //RooRealVar cp("cp","cp",1,-100,100);
        //RooRealVar dp("dp","dp",1,-100,100);
        RooChebychev poly("poly","poly",x,RooArgList(ap,bp,cp,dp));
        //RooPolynomial poly("poly","poly",x,RooArgList(ap,bp,cp,dp));
        RooRealVar polysig("polysig","polysig",4e6,0,10000000);
        RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));

        x.setRange("cut",0.43,0.565);
        //x.setRange("cut",0.45,0.54);

        RooFitResult* r_ks = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));

        double covQuality_ks = r_ks->covQual();
        //double covQuality_ks = 0;
        double mean_ks = mean.getVal();

        double gaus1F_ks = sig1.getVal();
        double gaus2F_ks = sig2.getVal();
        double polyF_ks  = polysig.getVal();
        double Norm = (gaus2F_ks + gaus1F_ks + polyF_ks)*0.0005;

        //set ranges for individual gaussian yield determination
        x.setRange("g1", mean.getVal() - 2*sigma1.getVal(), mean.getVal() + 2*sigma1.getVal());
        x.setRange("g2", mean.getVal() - 2*sigma2.getVal(), mean.getVal() + 2*sigma2.getVal());

        RooAbsReal* Intgaus1_yield_ks = gaus1.createIntegral(x, x, "g1");
        RooAbsReal* Intgaus2_yield_ks = gaus2.createIntegral(x, x, "g2");

        double gaus1_yield_ks = gaus1F_ks*Intgaus1_yield_ks->getVal();
        double gaus2_yield_ks = gaus2F_ks*Intgaus2_yield_ks->getVal();
        double gausTot_yield_ks = gaus1_yield_ks + gaus2_yield_ks;

        double rms_gaus1_sig_ks = gaus1_yield_ks/gausTot_yield_ks;
        double rms_gaus2_sig_ks = gaus2_yield_ks/gausTot_yield_ks;

        double rms_true_ks = TMath::Sqrt(rms_gaus1_sig_ks*sigma1.getVal()*sigma1.getVal() + rms_gaus2_sig_ks*sigma2.getVal()*sigma2.getVal());

        x.setRange("peak", mean.getVal() - 2*rms_true_ks, mean.getVal() + 2*rms_true_ks);


        RooAbsReal* Intgaus1_ks = gaus1.createIntegral(x, x, "peak");
        RooAbsReal* Intgaus2_ks = gaus2.createIntegral(x, x, "peak");
        RooAbsReal* Intpoly_ks  = poly.createIntegral(x, x, "peak");

        double Intgaus1E_ks = gaus1F_ks*Intgaus1_ks->getVal();
        double Intgaus2E_ks = gaus2F_ks*Intgaus2_ks->getVal();
        double IntpolyE_ks  = polyF_ks*Intpoly_ks->getVal();
        double totsig_ks    = Intgaus1E_ks + Intgaus2E_ks + IntpolyE_ks;
        double yield_ks       = Intgaus1E_ks + Intgaus2E_ks;

        double Fsig_ks = yield_ks/totsig_ks;

        mass_ks.push_back(mean_ks);
        std_ks.push_back(rms_true_ks);
        fsig_ks.push_back(Fsig_ks);
        covQual_ks.push_back(covQuality_ks);

        cout << "Yield (ks): "<< yield_ks << endl;
        cout << "Fsig (ks): " << Fsig_ks << endl;
        cout << "std (ks): " << rms_true_ks << endl;
        cout << "mass (ks): " << mean_ks << endl;
        cout << "covQual (ks)" << covQuality_ks << endl;

        RooPlot* xframe_ks = x.frame(270/Rebin);
        xframe_ks->GetXaxis()->SetTitle("#pi^{+} #pi^{-} Invariant mass (GeV)");
        xframe_ks->GetYaxis()->SetTitle("Candidates / 0.5 MeV");
        xframe_ks->GetXaxis()->CenterTitle(1);
        xframe_ks->GetYaxis()->CenterTitle(1);
        xframe_ks->GetXaxis()->SetTickSize(0.02);
        xframe_ks->GetYaxis()->SetTickSize(0.02);
        xframe_ks->GetXaxis()->SetNdivisions(407);
        xframe_ks->GetYaxis()->SetNdivisions(410);
        xframe_ks->GetXaxis()->SetTitleSize(0.04);
        xframe_ks->GetYaxis()->SetTitleSize(0.03);
        xframe_ks->GetYaxis()->SetTitleOffset(0.1);
        xframe_ks->GetXaxis()->SetTitleOffset(0.1);
        //xframe_ks->GetXaxis()->SetLabelSize(xframe_ks->GetXaxis()->GetLabelSize()*2.0);
        xframe_ks->GetYaxis()->SetLabelSize(0.02);
        xframe_ks->GetXaxis()->SetLabelSize(0.02);
        //xframe_ks->GetYaxis()->SetLabelSize(xframe_ks->GetYaxis()->GetLabelSize()*2.0);
        data.plotOn(xframe_ks,Name("data"));
        sum.plotOn(xframe_ks,Name("sum"),NormRange("cut"),LineWidth(2),LineColor(kBlue));
        sum.plotOn(xframe_ks,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(2),LineColor(kBlue));

        RooArgSet s(ap,bp,cp,dp,mean,polysig,sig1,sig2,sigma1);
        s.add(sigma2);
        TF1* func_ks = (TF1*)sum.asTF(RooArgList(x),RooArgList(s),x);
        std::vector<double> pull_ks;
        std::vector<double> pull_ks_x;
        std::vector<double> ratio_ks;
        int start = 860;
        if(mc) start = 0;
        for(int j=1; j<(270/Rebin + 1); j++)
        //for(int j=1; j<136; j++)
        {
            double data_point = massks[i]->GetBinContent(start/Rebin+j);
            double data_error = massks[i]->GetBinError(start/Rebin+j);
            double fit_point = Norm*func_ks->Eval(massks[i]->GetBinCenter(start/Rebin+j));
            if(data_point == 0) data_error = 1;
            pull_ks.push_back((data_point - fit_point)/data_error);
            pull_ks_x.push_back(massks[i]->GetBinCenter(start/Rebin+j));
            ratio_ks.push_back(data_point/fit_point);
        }

        TGraphErrors* TGpull = new TGraphErrors(270/Rebin,&pull_ks_x[0],&pull_ks[0],0,0);
        TGraphErrors* TGratio= new TGraphErrors(270/Rebin,&pull_ks_x[0],&ratio_ks[0],0,0);

        double chisquare = 0;

        for(int j=0; j<pull_ks.size(); j++)
        {
            chisquare+=TMath::Power(pull_ks[j],2);
        }

        int ndf = 270/Rebin - (10 - 1);

        /*
        if(i==8)
        {
            RooPlot* xframe_sks = x.frame(270);
            xframe_sks->GetXaxis()->SetTitle("#pi^{+}#pi^{-} Invariant Mass (GeV/c^{2})");
            xframe_sks->GetYaxis()->SetTitle("Candidates / 0.5 MeV");
            xframe_sks->GetXaxis()->CenterTitle(1);
            xframe_sks->GetYaxis()->CenterTitle(1);
            xframe_sks->GetXaxis()->SetTickSize(0.02);
            xframe_sks->GetYaxis()->SetTickSize(0.02);
            xframe_sks->GetXaxis()->SetNdivisions(407);
            xframe_sks->GetYaxis()->SetNdivisions(410);
            xframe_sks->GetXaxis()->SetTitleSize(0.05);
            xframe_sks->GetYaxis()->SetTitleSize(0.05);
            xframe_sks->GetYaxis()->SetTitleOffset(0.9);
            xframe_sks->GetXaxis()->SetTitleOffset(0.9);
            //xframe_sks->GetXaxis()->SetLabelSize(xframe_sks->GetXaxis()->GetLabelSize()*2.0);
            xframe_sks->GetYaxis()->SetLabelSize(0.04);
            xframe_sks->GetXaxis()->SetLabelSize(0.04);
            data.plotOn(xframe_sks,Name("data"));
            sum.plotOn(xframe_sks,Name("sum"),NormRange("cut"),LineWidth(2),LineColor(kBlue));
            sum.plotOn(xframe_sks,Name("poly"),Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(2),LineColor(kBlue-4));
            c_ks->cd();
            gPad->SetTickx();
            gPad->SetTicky();
            xframe_sks->Draw();
            c_ks->Update();
            TLine* t1_ks = new TLine(mean.getVal() - 2*rms_true_ks, 0, mean.getVal() - 2*rms_true_ks, gPad->GetUymax());
            TLine* t2_ks = new TLine(mean.getVal() + 2*rms_true_ks, 0, mean.getVal() + 2*rms_true_ks, gPad->GetUymax());
            t1_ks->SetLineStyle(2);
            t1_ks->SetLineColor(kGreen+1);
            t2_ks->SetLineStyle(2);
            t2_ks->SetLineColor(kGreen+1);
            t1_ks->Draw();
            t2_ks->Draw();

            double xstart_ks = 0.15;
            double ystart_ks = 0.89;
            double xpos = xstart_ks;
            double ypos = ystart_ks;
            double increment = 0.07;
            os << "CMS pPb #sqrt{S_{NN}} = 8.16 TeV";
            tex->DrawLatex(0.514,0.918,os.str().c_str());
            os.str(std::string());
            tex->SetTextSize(0.06);
            tex->DrawLatex(0.20,ypos-=increment,os.str().c_str());
            os.str(std::string());
            tex->SetTextSize(0.045);
            os << "185 #leq N_{trk}^{offline} < 250 ";
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            os << pks[i] << " < p_{T} < " << pks[i+1]  << " GeV/c";
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            os << "|y| < 1";
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());

            TLegend* leg = new TLegend(0.65,0.61,0.8,0.8);
            leg->SetFillColor(10);
            leg->SetFillStyle(0);
            leg->SetBorderSize(0);
            leg->SetTextFont(62);
            leg->SetTextSize(0.045);
            leg->AddEntry(xframe_ks->findObject("data"),"Data","P");
            leg->AddEntry(xframe_ks->findObject("sum"),"Fit","l");
            leg->AddEntry("poly","Background","l");
            leg->AddEntry(t1_ks,"#pm 2#sigma","l");
            leg->Draw();

            c_ks->Print("KsPlotForZhenyu.pdf");
        }
        */
        cc1->cd();

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
        xframe_ks->Draw();
        cc1->Update();

        TLine* t1_ks = new TLine(mean.getVal() - 2*rms_true_ks, 0, mean.getVal() - 2*rms_true_ks, gPad->GetUymax());
        TLine* t2_ks = new TLine(mean.getVal() + 2*rms_true_ks, 0, mean.getVal() + 2*rms_true_ks, gPad->GetUymax());
        t1_ks->SetLineStyle(2);
        t1_ks->SetLineColor(kGreen);
        t2_ks->SetLineStyle(2);
        t2_ks->SetLineColor(kGreen);
        t1_ks->Draw("same");
        t2_ks->Draw("same");
        if(i==0)
        {
            if(!doPbPb && !mc)
            {
                os << "CMS pPb";
                tex->SetTextSize(0.06);
                tex->DrawLatex(0.22,0.78,os.str().c_str());
                tex->SetTextSize(0.05);
                os.str(std::string());
                os << "185 #leq N_{trk}^{offline} < 250";
                tex->DrawLatex(0.22,0.71,os.str().c_str());
                os.str(std::string());
                os << "|y| < 1";
                tex->DrawLatex(0.22,0.64,os.str().c_str());
                os.str(std::string());
            }
            else if(mc)
            {
                os << "EPOS pPb";
                tex->SetTextSize(0.06);
                tex->DrawLatex(0.22,0.78,os.str().c_str());
                os.str(std::string());
                os << "|y| < 1";
                tex->DrawLatex(0.22,0.71,os.str().c_str());
                os.str(std::string());
            }
            else
            {
                os << "CMS PbPb";
                tex->SetTextSize(0.06);
                tex->DrawLatex(0.22,0.78,os.str().c_str());
                tex->SetTextSize(0.04);
                os.str(std::string());
                os << "Centrality 30 - 50%";
                tex->DrawLatex(0.22,0.71,os.str().c_str());
                os.str(std::string());
                os << "|y| < 1";
                tex->DrawLatex(0.22,0.64,os.str().c_str());
                os.str(std::string());
            }
        }
        double xstart_ks = 0.65;
        double ystart_ks = 0.85;
        double xpos = xstart_ks;
        double ypos = ystart_ks;
        double increment = 0.07;
        os << pks[i] << " < P_{t} < " << pks[i+1]  << " GeV";
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "Mean: " << std::setprecision(4) << mean_ks << " GeV" << std::setprecision(6);
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "#sigma :" << std::setprecision(2) << rms_true_ks << " GeV" << std::setprecision(6);
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "CovQual: " << covQuality_ks;
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        osYield << "Yield: " << std::setprecision(2) << yield_ks;
        tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
        osYield.str(std::string());
        os << "#chi^{2}/ndf: " << chisquare << "/" << ndf;
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        //osYield << "fsig: " << std::setprecision(6) << Fsig_ks;
        //tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
        //osYield.str(std::string());

        pad2->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        TH1F* frame_ks = pad2->cd()->DrawFrame(0.43,-9,0.565,9);
        frame_ks->GetXaxis()->SetTitle("#pi^{+} #pi^{-} Invariant mass (GeV)");
        frame_ks->GetYaxis()->SetTitle("Pull");
        frame_ks->GetXaxis()->CenterTitle(1);
        frame_ks->GetYaxis()->CenterTitle(1);
        frame_ks->GetXaxis()->SetTickSize(0.02);
        frame_ks->GetYaxis()->SetTickSize(0.02);
        frame_ks->GetXaxis()->SetNdivisions(410);
        frame_ks->GetYaxis()->SetNdivisions(407);
        frame_ks->GetXaxis()->SetTitleSize(0.12);
        frame_ks->GetYaxis()->SetTitleSize(0.12);
        frame_ks->GetYaxis()->SetTitleOffset(0.5);
        frame_ks->GetXaxis()->SetTitleOffset(1.2);
        frame_ks->GetYaxis()->SetLabelSize(0.1);
        frame_ks->GetXaxis()->SetLabelSize(0.1);
        TGpull->SetMarkerStyle(20);
        TGpull->Draw("P");
        TLine* line = new TLine(0.43, 0, 0.565, 0);
        line->SetLineStyle(2);
        line->Draw("same");

        pad3->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        double low = 0.1;
        double high = 2.75;
        if(doPbPb)
        {
            low = 0.5;
            high= 1.5;
        }
        frame_ks = pad3->cd()->DrawFrame(0.43,low,0.565,high);
        frame_ks->GetYaxis()->SetTitle("Ratio data/fit");
        frame_ks->GetXaxis()->CenterTitle(1);
        frame_ks->GetYaxis()->CenterTitle(1);
        frame_ks->GetXaxis()->SetTickSize(0.02);
        frame_ks->GetYaxis()->SetTickSize(0.02);
        frame_ks->GetXaxis()->SetNdivisions(410);
        frame_ks->GetYaxis()->SetNdivisions(407);
        frame_ks->GetXaxis()->SetTitleSize(0.12);
        frame_ks->GetYaxis()->SetTitleSize(0.12);
        frame_ks->GetYaxis()->SetTitleOffset(0.5);
        frame_ks->GetXaxis()->SetTitleOffset(1.2);
        frame_ks->GetYaxis()->SetLabelSize(0.1);
        frame_ks->GetXaxis()->SetLabelSize(0.1);
        TGratio->SetMarkerStyle(20);
        TGratio->Draw("P");
        TLine* line2 = new TLine(0.43, 1, 0.565, 1);
        line2->SetLineStyle(2);
        line2->Draw("same");

        Composite_Ks->cd(i+1);
        gPad->SetBottomMargin(0.15);
        gPad->SetLeftMargin(0.15);
        gPad->SetTickx();
        gPad->SetTicky();
        xframe_ks->GetXaxis()->SetNdivisions(507);
        xframe_ks->GetXaxis()->SetTitleSize(0.07);
        xframe_ks->GetYaxis()->SetTitleSize(0.07);
        xframe_ks->GetYaxis()->SetTitleOffset(0.95);
        xframe_ks->GetXaxis()->SetTitleOffset(1.1);
        xframe_ks->GetYaxis()->SetLabelSize(0.07);
        xframe_ks->GetXaxis()->SetLabelSize(0.07);
        xframe_ks->Draw();
        Composite_Ks->Update();

        t1_ks->Draw("same");
        t2_ks->Draw("same");
        xpos = xstart_ks;
        ypos = ystart_ks;
        increment = 0.07;
        if(i==0)
        {
            if(!doPbPb && !mc)
            {
                os << "CMS pPb";
                tex->SetTextSize(0.06);
                tex->DrawLatex(0.17,ypos-increment,os.str().c_str());
                tex->SetTextSize(0.05);
                os.str(std::string());
                os << "185 #leq N_{trk}^{offline} < 250";
                tex->DrawLatex(0.17,ypos-2*increment,os.str().c_str());
                os.str(std::string());
            }
            else if(mc)
            {
                os << "EPOS pPb";
                tex->SetTextSize(0.06);
                tex->DrawLatex(0.17,ypos-increment,os.str().c_str());
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
        os << pks[i] << " < P_{t} < " << pks[i+1] << " GeV";
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "Mean: " << std::setprecision(4) << mean_ks << " GeV" << std::setprecision(6);
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "#sigma :" << std::setprecision(2) << rms_true_ks << " GeV" << std::setprecision(6);
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "CovQual: " << covQuality_ks;
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        //os << "#chi^{2}/ndf: " << std::setprecision(3) << chi2_ks << std::setprecision(6);
        //tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        //os.str(std::string());
        osYield << "Yield: " << std::setprecision(2) << yield_ks;
        tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
        osYield.str(std::string());
        osYield << "fsig: " << std::setprecision(6) << Fsig_ks;
        tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
        osYield.str(std::string());
        if(pPb) cc1->Print(Form("KsMassFit_Pull_pPb%d.pdf",i));
        else if(mc) cc1->Print(Form("KsMassFit_Pull_MC%d.pdf",i));
        else if(doPbPb) cc1->Print(Form("KsMassFit_Pull_PbPb%d.pdf",i));
        else if(mb) cc1->Print(Form("KsMassFit_Pull_MB%d.pdf",i));

        if(lambda)
        {
            //lambda
            RooRealVar x("x","mass",1.08,1.16);
            RooPlot* xframe_la = x.frame(160);
            RooDataHist data("data","dataset",x,massla[i]);
            RooRealVar mean("mean","mean",1.115,1.11,1.12);
            RooRealVar sigma1("sigma1","sigma1",0.005,0.001,0.01);
            RooRealVar sigma2("sigma2","sigma2",0.005,0.001,0.01);
            //RooRealVar sig1("sig1","signal1",2e6,0,10000000); //pPb
            //RooRealVar sig2("sig2","signal2",8e6,0,10000000);
            //RooRealVar sig1("sig1","signal1",1e6,0,10000000); //pbpb
            //RooRealVar sig2("sig2","signal2",5e6,0,10000000);
            RooRealVar sig1("sig1","signal1",1e4,0,10000000); //mc
            RooRealVar sig2("sig2","signal2",5e4,0,10000000);
            RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
            RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
            RooRealVar ap("ap" , "ap" , 0    , -5 , 2);
            RooRealVar bp("bp" , "bp" , 0.1  , -5 , 2);
            RooRealVar cp("cp" , "cp" , -0.1 , -5 , 2);
            RooRealVar dp("dp" , "dp" , 0.1  , -5 , 2);
            //RooRealVar ap("ap","ap",10,0,1000000);
            //RooRealVar bp("bp","bp",10,0,1000000);
            //RooRealVar cp("cp","cp",10,0,1000000);
            //RooRealVar dp("dp","dp",10,0,1000000);
            //RooPolynomial poly("poly","poly",x,RooArgList(ap,bp,cp,dp));
            RooChebychev poly("poly","poly",x,RooArgList(ap,bp,cp,dp));
            //RooRealVar polysig("polysig","polysig",8.0e6,0,1e7);
            RooRealVar polysig("polysig","polysig",1.0e5,0,1e7); //mc
            //RooRealVar polysig("polysig","polysig",2e3,1e3,1e6);
            RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));

            if(i>=7 && i<10 && pPb)
            {
                sig1.setVal(1e6);
                sig2.setVal(1e6);
                polysig.setVal(2.3e5);
            }
            if(i==9 && mb)
            {
                sig1.setVal(1e4);
                sig2.setVal(1e4);
                polysig.setVal(2.3e4);
            }
            if(i==10 && mb)
            {
                sig1.setVal(0.5e4);
                sig2.setVal(0.5e4);
                polysig.setVal(1.3e4);
            }
            if((i==8) && mc)
            {
                sig1.setVal(1e4);
                sig2.setVal(1e4);
                polysig.setVal(2e4);
            }
            if((i==10) && mc)
            {
                sig1.setVal(10);
                sig2.setVal(10);
                polysig.setVal(8e2);
            }
            if(i==10 && doPbPb)
            {
                sig1.setVal(5e5);
                sig2.setVal(5e5);
                polysig.setVal(5e5);
            }
            if(i==0 && doPbPb)
            {
                sig1.setVal(1e6);
                sig2.setVal(4e5);
                polysig.setVal(4e6);
            }
            if(i==8 && doPbPb)
            {
                sig1.setVal(3e5);
                sig2.setVal(3e5);
                polysig.setVal(8e6);
            }
            if(i==9 && doPbPb)
            {
                sig1.setVal(2e3);
                sig2.setVal(2e3);
                polysig.setVal(4.5e6);
            }
            if(!doPbPb || mc) x.setRange("cut",1.080,1.159);
            //else x.setRange("cut",1.0920,1.1565);
            else x.setRange("cut",1.091,1.159);
            //if(i == pks.size()-2) polysig.setVal(4e7);

            RooFitResult* r_la = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));
            cout << "UID:" << i << endl;
            RooChi2Var chi2_laVar("chi2_laVar","chi2",sum,data);

            double covQuality_la = r_la->covQual();
            double mean_la = mean.getVal();

            double gaus1F_la = sig1.getVal();
            double gaus2F_la = sig2.getVal();
            double polyF_la  = polysig.getVal();
            Norm = (gaus1F_la + gaus2F_la + polyF_la)*0.0005;

            //set ranges for individual gaussian yield determination
            x.setRange("g1", mean.getVal() - 2*sigma1.getVal(), mean.getVal() + 2*sigma1.getVal());
            x.setRange("g2", mean.getVal() - 2*sigma2.getVal(), mean.getVal() + 2*sigma2.getVal());

            RooAbsReal* Intgaus1_yield_la = gaus1.createIntegral(x, x, "g1");
            RooAbsReal* Intgaus2_yield_la = gaus2.createIntegral(x, x, "g2");

            double gaus1_yield_la = gaus1F_la*Intgaus1_yield_la->getVal();
            double gaus2_yield_la = gaus2F_la*Intgaus2_yield_la->getVal();
            double gausTot_yield_la = gaus1_yield_la + gaus2_yield_la;

            double rms_gaus1_sig_la = gaus1_yield_la/gausTot_yield_la;
            double rms_gaus2_sig_la = gaus2_yield_la/gausTot_yield_la;

            double rms_true_la = TMath::Sqrt(rms_gaus1_sig_la*sigma1.getVal()*sigma1.getVal() + rms_gaus2_sig_la*sigma2.getVal()*sigma2.getVal());

            x.setRange("peak", mean.getVal() - 2*rms_true_la, mean.getVal() + 2*rms_true_la);
            RooAbsReal* Intpoly_la  = poly.createIntegral(x, x,"peak");
            RooAbsReal* Intgaus1_la = gaus1.createIntegral(x,x,"peak");
            RooAbsReal* Intgaus2_la = gaus2.createIntegral(x,x,"peak");


            double Intgaus1E_la = gaus1F_la*Intgaus1_la->getVal();
            double Intgaus2E_la = gaus2F_la*Intgaus2_la->getVal();
            double IntpolyE_la  = polyF_la*Intpoly_la->getVal();
            double totsig_la    = Intgaus1E_la + Intgaus2E_la + IntpolyE_la;
            double yield_la       = Intgaus1E_la + Intgaus2E_la;

            double Fsig_la = yield_la/totsig_la;

            mass_la.push_back(mean_la);
            std_la.push_back(rms_true_la);
            fsig_la.push_back(Fsig_la);
            covQual_la.push_back(r_la->covQual());

            cout << "Yield (la):" << yield_la << endl;
            cout << "Fsig (la): " << Fsig_la << endl;
            cout << "std (la): " << rms_true_la << endl;
            cout << "covQual (la): " << covQuality_la << endl;

            xframe_la->GetXaxis()->SetTitle("Invariant mass (GeV)");
            xframe_la->GetYaxis()->SetTitle("Candidates / 0.5 MeV");
            xframe_la->GetXaxis()->CenterTitle(1);
            xframe_la->GetYaxis()->CenterTitle(1);
            xframe_la->GetXaxis()->SetTickSize(0.02);
            xframe_la->GetYaxis()->SetTickSize(0.02);
            xframe_la->GetXaxis()->SetNdivisions(407);
            xframe_la->GetYaxis()->SetNdivisions(410);
            xframe_la->GetXaxis()->SetTitleSize(0.07);
            xframe_la->GetYaxis()->SetTitleSize(0.06);
            xframe_la->GetYaxis()->SetTitleOffset(0.95);
            xframe_la->GetXaxis()->SetTitleOffset(0.5);
            //xframe_la->GetXaxis()->SetLabelSize(xframe_la_->GetXaxis()->GetLabelSize()*2.0);
            xframe_la->GetYaxis()->SetLabelSize(0.08);
            xframe_la->GetXaxis()->SetLabelSize(0.08);
            //xframe_la->GetYaxis()->SetLabelSize(xframe_la_->GetYaxis()->GetLabelSize()*2.0);
            data.plotOn(xframe_la,Name("data"));
            sum.plotOn(xframe_la,Name("sum"),NormRange("cut"),LineWidth(2),LineColor(kRed));
            sum.plotOn(xframe_la,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(2),LineColor(kRed));

            RooArgSet s(ap,bp,cp,dp,mean,polysig,sig1,sig2,sigma1);
            s.add(sigma2);
            TF1* func = (TF1*)sum.asTF(RooArgList(x),RooArgList(s),x);
            std::vector<double> pull_la;
            std::vector<double> pull_la_x;
            std::vector<double> ratio_la;
            int nbins=0;
            int start = 0;
            if(!doPbPb && !mc)
            {
                nbins = 158;
                start = 1160;
            }
            else
            {
                nbins = 136;
                start = 1183;
            }
            if(mc)
            {
                nbins = 158;
                start = 0;
            }
            for(int j=0; j<nbins; j++)
            {
                double data_point = massla[i]->GetBinContent(start+j);
                double data_error = massla[i]->GetBinError(start+j);
                double fit_point = Norm*func->Eval(massla[i]->GetBinCenter(start+j));
                if(data_point == 0) data_error = 1;
                pull_la.push_back((data_point - fit_point)/data_error);
                pull_la_x.push_back(massla[i]->GetBinCenter(start+j));
                ratio_la.push_back(data_point/fit_point);
            }
            ndf = nbins - (10 - 1);

            TGraphErrors* TGpull = new TGraphErrors(nbins,&pull_la_x[0],&pull_la[0],0,0);
            TGraphErrors* TGratio= new TGraphErrors(nbins,&pull_la_x[0],&ratio_la[0],0,0);

            chisquare = 0;

            for(int j=0; j<pull_la.size(); j++)
            {
                chisquare+=TMath::Power(pull_la[j],2);
            }


            /*
            if(i==5)
            {
                RooPlot* xframe_sla = x.frame(160);
                xframe_sla->GetXaxis()->SetTitle("m_{p#pi^{-}} (GeV/c^{2})");
                xframe_sla->GetYaxis()->SetTitle("Candidates / 0.5 MeV");
                xframe_sla->GetXaxis()->CenterTitle(1);
                xframe_sla->GetYaxis()->CenterTitle(1);
                xframe_sla->GetXaxis()->SetTickSize(0.02);
                xframe_sla->GetYaxis()->SetTickSize(0.02);
                xframe_sla->GetXaxis()->SetNdivisions(407);
                xframe_sla->GetYaxis()->SetNdivisions(410);
                xframe_sla->GetXaxis()->SetTitleSize(0.05);
                xframe_sla->GetYaxis()->SetTitleSize(0.05);
                xframe_sla->GetYaxis()->SetTitleOffset(0.9);
                xframe_sla->GetXaxis()->SetTitleOffset(0.9);
                //xframe_sla->GetXaxis()->SetLabelSize(xframe_sla->GetXaxis()->GetLabelSize()*2.0);
                xframe_sla->GetYaxis()->SetLabelSize(0.04);
                xframe_sla->GetXaxis()->SetLabelSize(0.04);
                data.plotOn(xframe_sla,Name("data"));
                sum.plotOn(xframe_sla,Name("sum"),NormRange("cut"),LineWidth(2),LineColor(kRed));
                sum.plotOn(xframe_sla,Name("poly"),Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(2),LineColor(kRed));
                c_la->cd();
                gPad->SetTickx();
                gPad->SetTicky();
                xframe_sla->Draw();
                c_la->Update();
                TLine* t1_la = new TLine(mean.getVal() - 2*rms_true_la, 0, mean.getVal() - 2*rms_true_la, gPad->GetUymax());
                TLine* t2_la = new TLine(mean.getVal() + 2*rms_true_la, 0, mean.getVal() + 2*rms_true_la, gPad->GetUymax());
                t1_la->SetLineStyle(2);
                t1_la->SetLineColor(kGreen+1);
                t2_la->SetLineStyle(2);
                t2_la->SetLineColor(kGreen+1);
                t1_la->Draw();
                t2_la->Draw();

                double xstart_la = 0.13;
                double ystart_la = 0.89;
                double xpos = xstart_la;
                double ypos = ystart_la;
                double increment = 0.07;
                os << "CMS pPb #sqrt{S_{NN}} = 8.16 TeV";
                tex->DrawLatex(0.514,0.918,os.str().c_str());
                os.str(std::string());
                tex->SetTextSize(0.06);
                tex->DrawLatex(0.20,ypos-=increment,os.str().c_str());
                os.str(std::string());
                tex->SetTextSize(0.045);
                os << "185 #leq N_{trk}^{offline} < 250 ";
                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                os.str(std::string());
                os << pla[i] << " < p_{T} < " << pla[i+1]  << " GeV/c";
                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                os.str(std::string());
                os << "|y| < 1";
                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                os.str(std::string());

                TLegend* leg = new TLegend(0.63,0.61,0.78,0.8);
                leg->SetFillColor(10);
                leg->SetFillStyle(0);
                leg->SetBorderSize(0);
                leg->SetTextFont(62);
                leg->SetTextSize(0.045);
                leg->AddEntry(xframe_la->findObject("data"),"Data","P");
                leg->AddEntry(xframe_la->findObject("sum"),"Fit","l");
                leg->AddEntry("poly","Background","l");
                leg->AddEntry(t1_la,"#pm 2#sigma","l");
                leg->Draw();

                c_la->Print("LaPlotForZhenyu.pdf");
            }
        */
            cc2->cd();
            gPad->SetTickx();
            gPad->SetTicky();
            xframe_la->GetXaxis()->SetTitleOffset(1);
            xframe_la->GetXaxis()->SetTitleSize(xframe_la->GetXaxis()->GetTitleSize()*0.8);
            xframe_la->GetXaxis()->SetLabelSize(xframe_la->GetXaxis()->GetLabelSize()*0.5);
            xframe_la->GetYaxis()->SetLabelSize(xframe_la->GetYaxis()->GetLabelSize()*0.5);
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
            xframe_la->Draw();
            cc2->Update();

            TLine* t1_la = new TLine(mean.getVal() - 2*rms_true_la, 0, mean.getVal() - 2*rms_true_la, gPad->GetUymax());
            TLine* t2_la = new TLine(mean.getVal() + 2*rms_true_la, 0, mean.getVal() + 2*rms_true_la, gPad->GetUymax());
            t1_la->SetLineStyle(2);
            t1_la->SetLineColor(kGreen);
            t2_la->SetLineStyle(2);
            t2_la->SetLineColor(kGreen);
            t1_la->Draw("same");
            t2_la->Draw("same");
            if(i==0)
            {
                if(!doPbPb && !mc)
                {
                    os << "CMS pPb";
                    tex->SetTextSize(0.06);
                    tex->DrawLatex(0.22,0.78,os.str().c_str());
                    tex->SetTextSize(0.04);
                    os.str(std::string());
                    os << "185 #leq N_{trk}^{offline} < 250";
                    tex->DrawLatex(0.22,0.71,os.str().c_str());
                    os.str(std::string());
                    os << "|y| < 1";
                    tex->DrawLatex(0.22,0.64,os.str().c_str());
                    os.str(std::string());
                }
                else if(mc)
                {
                    os << "EPOS pPb";
                    tex->SetTextSize(0.06);
                    tex->DrawLatex(0.22,0.78,os.str().c_str());
                    os.str(std::string());
                    os << "|y| < 1";
                    tex->DrawLatex(0.22,0.71,os.str().c_str());
                    os.str(std::string());
                }
                else
                {
                    os << "CMS PbPb";
                    tex->SetTextSize(0.06);
                    tex->DrawLatex(0.22,0.78,os.str().c_str());
                    tex->SetTextSize(0.04);
                    os.str(std::string());
                    os << "Centrality 30 - 50%";
                    tex->DrawLatex(0.22,0.71,os.str().c_str());
                    os.str(std::string());
                    os << "|y| < 1";
                    tex->DrawLatex(0.22,0.64,os.str().c_str());
                    os.str(std::string());
                }
            }
            double xstart_la = 0.60;
            double ystart_la = 0.85;
            xpos = xstart_la;
            ypos = ystart_la;
            increment = 0.07;
            os << pla[i] << " < P_{t} < " << pla[i+1] << " GeV";
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            os << "Mean: " << std::setprecision(4) << mean_la << " GeV" << std::setprecision(6);
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            os << "#sigma :" << std::setprecision(2) << rms_true_la << " GeV" << std::setprecision(6);
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            os << "CovQual: " << covQuality_la;
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            osYield << "Yield: " << std::setprecision(2) << yield_la;
            tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
            osYield.str(std::string());
            os << "#chi^{2}/ndf: " << chisquare << "/" << ndf;
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            //osYield << "fsig: " << std::setprecision(6) << Fsig_la;
            //tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
            //osYield.str(std::string());

            pad2->cd();
            gPad->SetTickx();
            gPad->SetTicky();
            TH1F* frame_la = pad2->cd()->DrawFrame(1.08,-9,1.16,9);
            frame_la->GetXaxis()->SetTitle("#pi p Invariant mass (GeV)");
            frame_la->GetYaxis()->SetTitle("Pull");
            frame_la->GetXaxis()->CenterTitle(1);
            frame_la->GetYaxis()->CenterTitle(1);
            frame_la->GetXaxis()->SetTickSize(0.02);
            frame_la->GetYaxis()->SetTickSize(0.02);
            frame_la->GetXaxis()->SetNdivisions(410);
            frame_la->GetYaxis()->SetNdivisions(407);
            frame_la->GetXaxis()->SetTitleSize(0.12);
            frame_la->GetYaxis()->SetTitleSize(0.12);
            frame_la->GetYaxis()->SetTitleOffset(0.5);
            frame_la->GetXaxis()->SetTitleOffset(1.2);
            frame_la->GetYaxis()->SetLabelSize(0.1);
            frame_la->GetXaxis()->SetLabelSize(0.1);
            TGpull->SetMarkerStyle(20);
            TGpull->Draw("P");
            TLine* line = new TLine(1.08, 0, 1.16, 0);
            line->SetLineStyle(2);
            line->Draw("same");

            pad3->cd();
            gPad->SetTickx();
            gPad->SetTicky();
            low = 0.1;
            high = 2.75;
            if(doPbPb)
            {
                low = 0.5;
                high= 1.5;
            }
            frame_la = pad3->cd()->DrawFrame(1.08,low,1.160,high);
            frame_la->GetYaxis()->SetTitle("Ratio data/fit");
            frame_la->GetXaxis()->CenterTitle(1);
            frame_la->GetYaxis()->CenterTitle(1);
            frame_la->GetXaxis()->SetTickSize(0.02);
            frame_la->GetYaxis()->SetTickSize(0.02);
            frame_la->GetXaxis()->SetNdivisions(410);
            frame_la->GetYaxis()->SetNdivisions(407);
            frame_la->GetXaxis()->SetTitleSize(0.12);
            frame_la->GetYaxis()->SetTitleSize(0.12);
            frame_la->GetYaxis()->SetTitleOffset(0.5);
            frame_la->GetXaxis()->SetTitleOffset(1.2);
            frame_la->GetYaxis()->SetLabelSize(0.1);
            frame_la->GetXaxis()->SetLabelSize(0.1);
            TGratio->SetMarkerStyle(20);
            TGratio->Draw("P");
            TLine* line2 = new TLine(1.08, 1, 1.16, 1);
            line2->SetLineStyle(2);
            line2->Draw("same");

            Composite_La->cd(i+1);
            gPad->SetBottomMargin(0.15);
            gPad->SetLeftMargin(0.15);
            gPad->SetTickx();
            gPad->SetTicky();
            xframe_la->GetXaxis()->SetNdivisions(507);
            xframe_la->GetXaxis()->SetTitleSize(0.06);
            xframe_la->GetYaxis()->SetTitleSize(0.06);
            xframe_la->GetYaxis()->SetTitleOffset(1.1);
            xframe_la->GetXaxis()->SetTitleOffset(1.1);
            xframe_la->GetYaxis()->SetLabelSize(0.06);
            xframe_la->GetXaxis()->SetLabelSize(0.06);
            xframe_la->Draw();
            Composite_La->Update();

            t1_la->Draw("same");
            t2_la->Draw("same");
            xpos = xstart_la;
            ypos = ystart_la;
            increment = 0.07;
            if(i==10)
            {
                xpos = 0.65;
            }
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
            os << pla[i] << " < P_{t} < " << pla[i+1] << " GeV";
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            os << "Mean: " << std::setprecision(4) << mean_la << " GeV" << std::setprecision(6);
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            os << "#sigma :" << std::setprecision(2) << rms_true_la << " GeV" << std::setprecision(6);
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            os << "CovQual: " << covQuality_la;
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            //os << "#chi^{2}/ndf: " << std::setprecision(3) << chi2_la << std::setprecision(6);
            //tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            //os.str(std::string());
            osYield << "Yield: " << std::setprecision(2) << yield_la;
            tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
            osYield.str(std::string());
            osYield << "fsig: " << std::setprecision(6) << Fsig_la;
            tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
            osYield.str(std::string());
        }

        if(pPb) cc2->Print(Form("LaMassFit_Pull_pPb%d.pdf",i));
        else if(mc)cc2->Print(Form("LaMassFit_Pull_MC%d.pdf",i));
        else if(doPbPb) cc2->Print(Form("LaMassFit_Pull_PbPb%d.pdf",i));
        else if(mb) cc2->Print(Form("LaMassFit_Pull_MB%d.pdf",i));
        //if(i==0) cc1->Print("V0MassFitInd.pdf(","pdf");
        //else if(i < pks.size() - 2) cc1->Print("V0MassFitInd.pdf","pdf");
        //else cc1->Print("V0MassFitInd.pdf)","pdf");
        hbincounter++;
    }
    if(doPbPb)
    {
        Composite_Ks->Print("KsMassFitCompositePbPb.pdf");
        Composite_La->Print("LaMassFitCompositePbPb.pdf");
    }
    else
    {
        Composite_Ks->Print("KsMassFitComposite.pdf");
        Composite_La->Print("LaMassFitComposite.pdf");
    }
    //Output
    pkscounter = 0;
    myfile << "KSHORT KSHORT KSHORT\n";
    for(unsigned i=0; i<mass_ks.size(); i++)
    {
        cout <<  "====================" << endl;
        cout << "Pt Bin: " << (pks[pkscounter]-1)/10 << " - " << pks[pkscounter+1]/10 << endl;
        cout <<  "====================" << endl;
        cout << "Mass_ks: "   << mass_ks[i]    << endl;
        cout << "Fsig_ks: "   << fsig_ks[i]    << endl;
        cout << "std_ks: "    << std_ks[i]     << endl;
        cout << "covQual_ks " << covQual_ks[i] << endl;

        myfile <<  "====================" << "\n";
        myfile << "Pt Bin: " << (pks[pkscounter]-1)/10 << " - " << pks[pkscounter+1]/10 << "\n";
        myfile <<  "====================" << "\n";
        myfile << "Mass_ks: "   << mass_ks[i]    << "\n";
        myfile << "Fsig_ks: "   << fsig_ks[i]    << "\n";
        myfile << "std_ks: "    << std_ks[i]     << "\n";
        myfile << "covQual_ks " << covQual_ks[i] << "\n";
        pkscounter+=2;
    }

    placounter=6;
    myfile << "LAMBDA LAMBDA LAMBDA\n";
    for(unsigned i=0; i<mass_la.size(); i++)
    {
        cout <<  "====================" << endl;
        cout << "Pt Bin: " << (pla[placounter]-1)/10 << " - " << pla[placounter+1]/10 << endl;
        cout <<  "====================" << endl;
        cout << "Mass_la: "   << mass_la[i]    << endl;
        cout << "Fsig_la: "   << fsig_la[i]    << endl;
        cout << "std_la: "    << std_la[i]     << endl;
        cout << "covQual_la " << covQual_la[i] << endl;

        myfile <<  "====================" << "\n";
        myfile << "Pt Bin: " << (pla[placounter]-1)/10 << " - " << pla[placounter+1]/10 << "\n";
        myfile <<  "====================" << "\n";
        myfile << "Mass_la: "   << mass_la[i]    << "\n";
        myfile << "Fsig_la: "   << fsig_la[i]    << "\n";
        myfile << "std_la: "    << std_la[i]     << "\n";
        myfile << "covQual_la " << covQual_la[i] << "\n";
        placounter+=2;
    }
}
