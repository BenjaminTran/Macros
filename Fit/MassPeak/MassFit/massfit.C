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

void massfit()
{
    gStyle->SetMarkerSize(0.8);
    
    TFile* file = new TFile("ppCorr_all_V0_pt13.root");

    TH1D* massks;
    TH1D* massla;

    massks = (TH1D*)file->Get("pp_HM105_GplusPP/masskshort_pt0");
    massla = (TH1D*)file->Get("pp_HM105_GplusPP/masslambda_pt0");

    using namespace RooFit;
    gStyle->SetOptTitle(kFALSE);

    TCanvas* cc1 = new TCanvas("cc1","cc1",1200,450);
    cc1->Divide(2,1);
    
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    
    char label_energy[200]={"CMS pp #sqrt{s} = 13 TeV"};
    char label_n[200]={"105 #leq N^{offline}_{trk} < 150"};
    char label_pid[2][200]={"K^{0}_{S}","#Lambda/#bar{#Lambda}"};
    char label_mean[2][200]={"Mean: 0.4976 GeV","Mean: 1.1159 GeV"};
    char label_sigma[2][200]={"Average #sigma: 0.0067 GeV","Average #sigma: 0.0031 GeV"};
    char label_cms[2][200]={"Preliminary","L = 0.7 pb^{-1}"};
    
    tex->SetTextSize(tex->GetTextSize()*0.95);
    
    //kshort
        RooRealVar x("x","mass",0.43,0.565);
        RooDataHist data("data","dataset",x,massks);
        RooPlot* xframe = x.frame(270);
        xframe->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");
        xframe->GetYaxis()->SetTitle("Candidates / 0.0005 GeV");
        xframe->GetXaxis()->CenterTitle(1);
        xframe->GetYaxis()->CenterTitle(1);
        xframe->GetXaxis()->SetTitleSize(xframe->GetXaxis()->GetTitleSize()*1.4);
    xframe->GetXaxis()->SetTitleOffset(1);
        xframe->GetYaxis()->SetTitleSize(xframe->GetYaxis()->GetTitleSize()*1.3);
    //xframe->GetYaxis()->SetTitleOffset(1);
        data->plotOn(xframe,Name("data"));
        RooRealVar mean("mean","mean",0.50,0.49,0.51);
        RooRealVar sigma1("sigma1","sigma1",0.01,0.001,0.04);
        RooRealVar sigma2("sigma2","sigma2",0.01,0.001,0.04);
        RooRealVar sig1("sig1","signal1",10,0,10000000);
        RooRealVar sig2("sig2","signal2",10,0,10000000);
        RooRealVar a("a","a",0,-100000,100000);
        RooRealVar b("b","b",0,-100000,100000);
        RooRealVar cp("cp","cp",0,-100000,100000);
        RooRealVar d("d","d",0,-100000,100000);
        RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
        RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
        RooPolynomial poly("poly","poly",x,RooArgList(a,b,cp,d));
        RooRealVar polysig("polysig","polysig",10,0,10000000);
        RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));
        
        x->setRange("cut",0.44,0.56);

        sum->fitTo(data,Range("cut"));
        sum->fitTo(data,Range("cut"));
        sum->fitTo(data,Range("cut"));
        sum->fitTo(data,Range("cut"));
        sum->fitTo(data,Range("cut"));
        sum->fitTo(data,Range("cut"));
        
        
        sum->plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(0.5),LineColor(kBlue));
        sum->plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(0.5),LineColor(kBlue));
        cc1->cd(1);
        xframe->Draw();
    tex->DrawLatex(0.59,0.87,label_mean[0]);
    tex->DrawLatex(0.59,0.81,label_sigma[0]);
    
    
    //lambda
    RooRealVar x("x","mass",1.08,1.16);
    RooDataHist data("data","dataset",x,massla);
    RooPlot* xframe = x.frame(160);
    xframe->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");
    xframe->GetYaxis()->SetTitle("Candidates / 0.0005GeV");
    xframe->GetXaxis()->CenterTitle(1);
    xframe->GetYaxis()->CenterTitle(1);
    xframe->GetXaxis()->SetTitleSize(xframe->GetXaxis()->GetTitleSize()*1.4);
    xframe->GetYaxis()->SetTitleSize(xframe->GetYaxis()->GetTitleSize()*1.3);
    xframe->GetXaxis()->SetTitleOffset(1);
    data->plotOn(xframe,Name("data"));
    RooRealVar mean("mean","mean",1.115,1.11,1.12);
    RooRealVar sigma1("sigma1","sigma1",0.005,0.001,0.01);
    RooRealVar sigma2("sigma2","sigma2",0.005,0.001,0.01);
    RooRealVar sig1("sig1","signal1",10,0,10000000);
    RooRealVar sig2("sig2","signal2",10,0,10000000);
    RooRealVar a("a","a",0,-100000,100000);
    RooRealVar b("b","b",0,-100000,100000);
    RooRealVar cp("cp","cp",0,-100000,100000);
    RooRealVar d("d","d",0,-100000,100000);
    RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
    RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
    RooPolynomial poly("poly","poly",x,RooArgList(a,b,cp,d));
    RooRealVar polysig("polysig","polysig",10,0,10000000);
    RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));
    
    x->setRange("cut",1.09,1.14);
    
    sum->fitTo(data,Range("cut"));
    sum->fitTo(data,Range("cut"));
    sum->fitTo(data,Range("cut"));
    sum->fitTo(data,Range("cut"));
    sum->fitTo(data,Range("cut"));
    sum->fitTo(data,Range("cut"));
    
    sum->plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(0.5),LineColor(kRed));
    sum->plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(0.5),LineColor(kRed));
    cc1->cd(2);
    xframe->Draw();
    tex->DrawLatex(0.59,0.87,label_mean[1]);
    tex->DrawLatex(0.59,0.81,label_sigma[1]);
    
    cc1->cd(1);
    tex->DrawLatex(0.22,0.87,label_energy);
    tex->DrawLatex(0.22,0.81,label_n);
    tex->DrawLatex(0.22,0.75,label_pid[0]);
    tex->DrawLatex(0.22,0.69,"1 < p_{T} < 3 GeV/c");
    tex->DrawLatex(0.22,0.63,"Preliminary");
    //tex->DrawLatex(0.15,0.58,label_cms[1]);

    cc1->cd(2);
    tex->DrawLatex(0.22,0.87,label_energy);
    tex->DrawLatex(0.22,0.81,label_n);
    tex->DrawLatex(0.22,0.75,label_pid[1]);
    tex->DrawLatex(0.22,0.69,"1 < p_{T} < 3 GeV/c");
    tex->DrawLatex(0.22,0.63,"Preliminary");
    //tex->DrawLatex(0.15,0.58,label_cms[1]);

    cc1->Print("massfit.pdf");
    cc1->Print("massfit.gif");
}
