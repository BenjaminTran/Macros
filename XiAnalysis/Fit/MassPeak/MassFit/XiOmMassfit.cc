#include <iostream>

#include "TF1.h"
#include "TH1.h"
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
//#include "RooFit.h"
//#include "RooRealVar.h"
//#include "RooPlot.h"
//#include "RooDataHist.h"
//#include "RooGaussian.h"
//#include "RooPolynomial.h"
//#include "RooAddPdf.h"
#include "TString.h"

#include <vector>

void XiOmMassfit()
{
    TH1::SetDefaultSumw2();
    RooMsgService::instance().setStreamStatus(0,kFALSE);
    RooMsgService::instance().setStreamStatus(1,kFALSE);
    TGaxis::SetMaxDigits(3);
    gStyle->SetMarkerSize(0.8);
    TFile* f1_Xi;
    TFile* f1_Om;
    using namespace RooFit;
    gStyle->SetOptTitle(kFALSE);

    bool do_pPb = true;
    bool do_PbPb = !do_pPb;
    bool doxi = true;
    bool doom = true;

    if(do_pPb)
    {
        f1_Xi = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationHM_11_07_17.root");
        f1_Om = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/OmCorr/OmegaHMRapidity_Total_12_04_17.root");
    }
    else if(do_PbPb)
    {
        f1_Xi = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/V0CasCorrelationPbPbTotal_10_30_17.root");
        f1_Om = f1_Xi;
    }

    TCanvas* cc1 = new TCanvas("cc1","cc1",700,550);
    TCanvas* cc2 = new TCanvas("cc2","cc2",700,550);
    //cc1->Divide(2,1);
    
    TH2D* MassXi;
    TH2D* MassOm;
    TH1D* massxi;
    TH1D* massom;
    
    TLatex* tex = new TLatex();
    TLatex* tex_pid = new TLatex();
    tex->SetNDC();
    tex_pid->SetNDC();
    
    char label_energy[2][200]={"CMS #scale[1.1]{#font[12]{Preliminary}}","pPb #sqrt{s_{NN}} = 8.16 TeV"};
    char label_n[200]={"185 #leq N^{offline}_{trk} < 250"};
    char label_pid[2][200]={"#Xi^{-}","#Omega^{-}"};
    char label_mean[2][200]={"Mean: 1.3228 GeV","Mean: 1.6727 GeV"};
    char label_pt[2][200] = {"1 < p_{T} < 3 GeV","1.5 < p_{T} < 3 GeV"};
    char label_sigma[2][200]={"Average #sigma: 0.0043 GeV","Average #sigma: 0.0044 GeV"};
    char label_cms[2][200]={"Preliminary","L_{int} = 162 nb^{-1}"};

    char label_energy_PbPb[2][200]={"CMS #scale[1.1]{#font[12]{Preliminary}}","PbPb #sqrt{s_{NN}} = 5.02 TeV"};
    char label_n_PbPb[200]={"30-50% Centrality"};
    char label_mean_PbPb[2][200]={"Mean: 1.3227 GeV","Mean: 1.6727 GeV"};
    char label_sigma_PbPb[2][200]={"Average #sigma: 0.0042 GeV","Average #sigma: 0.0059 GeV"};
    
    tex->SetTextSize(tex->GetTextSize()*0.75);
    tex_pid->SetTextSize(tex_pid->GetTextSize());

    if(do_pPb)
    {
        MassXi = (TH2D*)f1_Xi->Get("v0CasCorrelationRapidity/MassPtXi");
        MassOm = (TH2D*)f1_Om->Get("v0CasCorrelationRapidity/MassPtOm");
    }
    if(do_PbPb)
    {
        MassXi = (TH2D*)f1_Xi->Get("v0CasCorrelationRapidityPbPb/MassPtXi");
        MassOm = (TH2D*)f1_Om->Get("v0CasCorrelationRapidityPbPb/MassPtOm");
    }

    massxi = (TH1D*)MassXi->ProjectionX("massxi", 11,30);
    massom = (TH1D*)MassOm->ProjectionX("massom", 16,30);
    
    if(doxi)
    {
        //kshort
        //massxi = (TH1D*)f1->Get(Form("demo/massxihort_pt%d",0));
        RooRealVar x("x","mass",1.25,1.40);
        RooPlot* xframe = x.frame(150);
        RooDataHist data("data","dataset",x,massxi);
        data.plotOn(xframe,Name("data"));
        xframe->GetXaxis()->SetTitle("#Lambda#pi^{-} and charge conjugate invariant mass (GeV/c^{2})");
        xframe->GetYaxis()->SetTitle("Candidates / 1 MeV");
        xframe->GetXaxis()->CenterTitle(1);
        xframe->GetYaxis()->CenterTitle(1);
        xframe->GetXaxis()->SetTitleSize(xframe->GetXaxis()->GetTitleSize()*1.3);
        xframe->GetYaxis()->SetTitleSize(xframe->GetYaxis()->GetTitleSize()*1.2);
        xframe->GetYaxis()->SetTitleOffset(1);
        //xframe->GetYaxis()->SetRangeUser(0,20e6);
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

        x.setRange("cut",1.26,1.39);

        sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));


        sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(2),LineColor(kBlue));
        sum.plotOn(xframe,Components(background),NormRange("cut"),LineStyle(kDashed),LineWidth(2),LineColor(kBlue));
        cc1->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        xframe->Draw();
        if(do_pPb)
        {
            tex->DrawLatex(0.55,0.80,label_mean[0]);
            tex->DrawLatex(0.55,0.74,label_sigma[0]);
        }
        else if(do_PbPb)
        {
            tex->DrawLatex(0.55,0.80,label_mean[0]);
            tex->DrawLatex(0.55,0.74,label_sigma_PbPb[0]);
        }
    }
    if(doom)
    {
        //lambda
        //massom = (TH1D*)f1->Get(Form("demo/massommbda_pt%d",0));
        RooRealVar x("x","mass",1.60,1.75);
        RooDataHist data("data","dataset",x,massom);
        RooPlot* xframe = x.frame(150);
        xframe->GetXaxis()->SetTitle("#Lambda K + charge conjugate invariant mass (GeV/c^{2})");
        xframe->GetYaxis()->SetTitle("Candidates / 1 MeV");
        xframe->GetXaxis()->CenterTitle(1);
        xframe->GetYaxis()->CenterTitle(1);
        xframe->GetXaxis()->SetTitleSize(xframe->GetXaxis()->GetTitleSize()*1.3);
        xframe->GetYaxis()->SetTitleSize(xframe->GetYaxis()->GetTitleSize()*1.2);
        xframe->GetYaxis()->SetTitleOffset(1);
        //xframe->GetYaxis()->SetRangeUser(0,12e6);
        data.plotOn(xframe,Name("data"));
        RooRealVar mean("mean","mean",1.67,1.65,1.69);
        RooRealVar sigma1("sigma1","sigma1",0.004,0.001,0.04);
        RooRealVar sigma2("sigma2","sigma2",0.005,0.001,0.04);
        RooRealVar sig1("sig1","signal1",45000,0,1000000);
        RooRealVar sig2("sig2","signal2",45000,0,1000000);
        RooRealVar qsig("qsig","qsig",52000,0,1000000);
        //RooRealVar sig1("sig1","signal1",1000,-100,1000000);
        //RooRealVar sig2("sig2","signal2",1000,-100,1000000);
        //RooRealVar qsig("qsig","qsig",1000,0,1000000);
        RooRealVar alpha("alpha","alpha",0.001,-1,10);
        RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
        RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
        RooRealVar ap("ap","ap",-0.5,-1,1);
        RooRealVar bp("bp","bp",-0.1,-1,1);
        RooRealVar cp("cp","cp",-0.1,-1,1);
        RooRealVar dp("dp","dp",-0.1,-1,1);
        RooChebychev background("background","background",x,RooArgList(ap,bp,cp,dp));
        //RooExtendPdf egaus1("egaus1","egaus1",gaus1,sig1);
        //RooExtendPdf egaus2("egaus2","egaus2",gaus2,sig2);
        //RooExtendPdf ebackground("ebackground","ebackground",background,qsig);
        //RooGenericPdf background("background", "x - (1.115683 + 0.493677)^alpha", RooArgList(x,alpha));
        RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,background),RooArgList(sig1,sig2,qsig));
        //RooAddPdf sum("sum","sum",RooArgList(egaus1,egaus2,ebackground));

        x.setRange("cut",1.62,1.75);

        sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));

        sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(2),LineColor(kRed));
        sum.plotOn(xframe,Components(background),NormRange("cut"),LineStyle(kDashed),LineWidth(2),LineColor(kRed));
        cc2->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        xframe->Draw();
        if(do_pPb)
        {
            tex->DrawLatex(0.55,0.80,label_mean[1]);
            tex->DrawLatex(0.55,0.74,label_sigma[1]);
        }
        else if(do_PbPb)
        {
            tex->DrawLatex(0.55,0.80,label_mean[1]);
            tex->DrawLatex(0.55,0.74,label_sigma_PbPb[1]);
        }

        cc1->cd();
        double x_start = 0.139;
        double y_start = 0.83+0.06;
        double increment = 0.06;
        if(do_pPb)
        {
            tex->DrawLatex(x_start,y_start-=increment,label_energy[0]);
            tex->DrawLatex(x_start,y_start-=increment,label_energy[1]);
            //tex->DrawLatex(x_start,y_start-=increment,label_cms[1]);
            tex->DrawLatex(x_start,y_start-=increment,label_n);
            tex->DrawLatex(x_start,y_start-=increment,label_pt[0]);
            tex_pid->DrawLatex(x_start,y_start-=increment,label_pid[0]);
        }
        else if(do_PbPb)
        {
            tex->DrawLatex(x_start,y_start-=increment,label_energy_PbPb[0]);
            tex->DrawLatex(x_start,y_start-=increment,label_energy_PbPb[1]);
            //tex->DrawLatex(x_start,y_start-=increment,label_cms[1]);
            tex->DrawLatex(x_start,y_start-=increment,label_n_PbPb);
            tex->DrawLatex(x_start,y_start-=increment,label_pt[0]);
            tex_pid->DrawLatex(x_start,y_start-=increment,label_pid[0]);
        }

        cc2->cd();
        x_start = 0.139;
        y_start = 0.83+0.06;
        if(do_pPb)
        {
            tex->DrawLatex(x_start,y_start-=increment,label_energy[0]);
            tex->DrawLatex(x_start,y_start-=increment,label_energy[1]);
            //tex->DrawLatex(x_start,y_start-=increment,label_cms[1]);
            tex->DrawLatex(x_start,y_start-=increment,label_n);
            tex->DrawLatex(x_start,y_start-=increment,label_pt[1]);
            tex_pid->DrawLatex(x_start,y_start-=increment,label_pid[1]);
            //tex->DrawLatex(0.22,0.70,label_cms[0]);
        }
        else if(do_PbPb)
        {
            tex->DrawLatex(x_start,y_start-=increment,label_energy_PbPb[0]);
            tex->DrawLatex(x_start,y_start-=increment,label_energy_PbPb[1]);
            //tex->DrawLatex(x_start,y_start-=increment,label_cms[1]);
            tex->DrawLatex(x_start,y_start-=increment,label_n_PbPb);
            tex->DrawLatex(x_start,y_start-=increment,label_pt[1]);
            tex_pid->DrawLatex(x_start,y_start-=increment,label_pid[1]);
        }
    }

    if(do_pPb)
    {
        cc1->Print("Image/XiOmMassfit/massfit_pPb_Xi.pdf");
        cc2->Print("Image/XiOmMassfit/massfit_pPb_Om.pdf");
    }
    else if(do_PbPb)
    {
        cc1->Print("Image/XiOmMassfit/massfit_PbPb_Xi.pdf");
        cc2->Print("Image/XiOmMassfit/massfit_PbPb_Om.pdf");
    }
}
