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

void MassFitV0()
{
    TH1::SetDefaultSumw2();
    RooMsgService::instance().setStreamStatus(0,kFALSE);
    RooMsgService::instance().setStreamStatus(1,kFALSE);
    TGaxis::SetMaxDigits(3);
    gStyle->SetMarkerSize(0.8);

    bool do_pPb = true;
    bool do_PbPb = !do_pPb;
    bool doks = true;
    bool dola = true;


    TFile* f1;

    if(do_pPb)
    {
        f1 = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0Correlation_TrkEff_1_11_18.root");
    }
    else if (do_PbPb)
    {
        f1 = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/V0CasCorrelationPbPbTotal_10_30_17.root");
    }
    using namespace RooFit;
    gStyle->SetOptTitle(kFALSE);


    TCanvas* cc1 = new TCanvas("cc1","cc1",700,550);
    TCanvas* cc2 = new TCanvas("cc2","cc2",700,550);
    //cc1->Divide(2,1);
    
    TH1D* massks;
    TH1D* massla;
    
    TLatex* tex = new TLatex();
    TLatex* tex_pid = new TLatex();
    tex->SetNDC();
    tex_pid->SetNDC();
    
    char label_energy[2][200]={"CMS #scale[1.1]{#font[12]{Preliminary}}","pPb #sqrt{s_{NN}} = 8.16 TeV"};
    char label_n[200]={"185 #leq N^{offline}_{trk} < 250"};
    char label_pid[2][200]={"K^{0}_{S}","#Lambda/#bar{#Lambda}"};
    char label_mean[2][200]={"Mean: 0.4976 GeV","Mean: 1.1159 GeV"};
    char label_pt[200] = {"1 < p_{T} < 3 GeV"};
    char label_sigma[2][200]={"Average #sigma: 0.0047 GeV","Average #sigma: 0.0018 GeV"};
    char label_cms[2][200]={"Preliminary","L_{int} = 162 nb^{-1}"};
    
    tex->SetTextSize(tex->GetTextSize()*0.75);

    std::vector<TH1D*> massks_h;
    std::vector<TH1D*> massla_h;

    std::vector<double> pks = {1.0,1.4,1.8,2.2,2.8};
    std::vector<double> pla = {1.0,1.4,1.8,2.2,2.8};

    for(int i=4; i<(4+pks.size()); i++)
    {
        if(do_pPb) massks_h.push_back((TH1D*)f1->Get(Form("v0CasCorrelationRapidity/masskshort_pt%d",i)));
        else if(do_PbPb) massks_h.push_back((TH1D*)f1->Get(Form("v0CasCorrelationRapidityPbPb/masskshort_pt%d",i)));
    }
    for(int i=1; i<(1+pla.size()); i++)
    {
        if(do_pPb) massla_h.push_back((TH1D*)f1->Get(Form("v0CasCorrelationRapidity/masslambda_pt%d",i)));
        else if(do_PbPb) massla_h.push_back((TH1D*)f1->Get(Form("v0CasCorrelationRapidityPbPb/masslambda_pt%d",i)));
    }

    for(int i=1; i<massks_h.size(); i++)
    {
        massks_h[0]->Add(massks_h[i]);
        massla_h[0]->Add(massla_h[i]);
    }

    massks = massks_h[0];
    massla = massla_h[0];
    cout << "Setup" << endl;
    
    if(doks)
    {
        //kshort
        //massks = (TH1D*)f1->Get(Form("demo/masskshort_pt%d",0));
        RooRealVar x("x","mass",0.43,0.565);
        RooPlot* xframe = x.frame(270);
        RooDataHist data("data","dataset",x,massks);
        data.plotOn(xframe,Name("data"));
        xframe->GetXaxis()->SetTitle("#pi^{+}#pi^{-} invariant mass (GeV/c^{2})");
        xframe->GetYaxis()->SetTitle("Candidates / 0.5 MeV");
        xframe->GetXaxis()->CenterTitle(1);
        xframe->GetYaxis()->CenterTitle(1);
        xframe->GetXaxis()->SetTitleSize(xframe->GetXaxis()->GetTitleSize()*1.3);
        xframe->GetYaxis()->SetTitleSize(xframe->GetYaxis()->GetTitleSize()*1.2);
        xframe->GetYaxis()->SetTitleOffset(1);
        xframe->GetYaxis()->SetRangeUser(0,20e6);
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

        x.setRange("cut",0.44,0.56);

        sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));


        sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(2),LineColor(kBlue));
        sum.plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(2),LineColor(kBlue));
        cc1->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        xframe->Draw();
        tex->DrawLatex(0.55,0.77,label_mean[0]);
        tex->DrawLatex(0.55,0.71,label_sigma[0]);
    }
    cout << "Did Ks" << endl;
    if(dola)
    {
        //lambda
        //massla = (TH1D*)f1->Get(Form("demo/masslambda_pt%d",0));
        RooRealVar x("x","mass",1.08,1.16);
        RooDataHist data("data","dataset",x,massla);
        RooPlot* xframe = x.frame(160);
        xframe->GetXaxis()->SetTitle("p#pi^{-} + charge conjugate invariant mass (GeV/c^{2})");
        xframe->GetYaxis()->SetTitle("Candidates / 0.5 MeV");
        xframe->GetXaxis()->CenterTitle(1);
        xframe->GetYaxis()->CenterTitle(1);
        xframe->GetXaxis()->SetTitleSize(xframe->GetXaxis()->GetTitleSize()*1.3);
        xframe->GetYaxis()->SetTitleSize(xframe->GetYaxis()->GetTitleSize()*1.2);
        xframe->GetYaxis()->SetTitleOffset(1);
        xframe->GetYaxis()->SetRangeUser(0,12e6);
        data.plotOn(xframe,Name("data"));
        RooRealVar mean("mean","mean",1.115,1.11,1.12);
        RooRealVar sigma1("sigma1","sigma1",0.002,0.001,0.01);
        RooRealVar sigma2("sigma2","sigma2",0.002,0.001,0.01);
        RooRealVar sig1("sig1","signal1",8e7,0,100000000); //pPb
        RooRealVar sig2("sig2","signal2",8e7,0,100000000);
        //RooRealVar sig1("sig1","signal1",1e6,0,10000000); //pbpb
        //RooRealVar sig2("sig2","signal2",5e6,0,10000000);
        RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
        RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
        RooRealVar ap("ap" , "ap" , 0    , -5 , 2);
        RooRealVar bp("bp" , "bp" , 0.1  , -5 , 2);
        RooRealVar cp("cp" , "cp" , -0.1 , -5 , 2);
        RooRealVar dp("dp" , "dp" , 0.1  , -5 , 2);
        RooChebychev poly("poly","poly",x,RooArgList(ap,bp,cp,dp));
        RooRealVar polysig("polysig","polysig",1.0e7,0,1e9);
        RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));

        x.setRange("cut",1.085,1.155);

        sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));

        sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(2),LineColor(kRed));
        sum.plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(2),LineColor(kRed));
        cc2->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        xframe->Draw();
        tex->DrawLatex(0.55,0.77,label_mean[1]);
        tex->DrawLatex(0.55,0.71,label_sigma[1]);

        cc1->cd();
        double x_start = 0.139;
        double y_start = 0.83+0.06;
        double increment = 0.06;
        tex->DrawLatex(x_start,y_start-=increment,label_energy[0]);
        tex->DrawLatex(x_start,y_start-=increment,label_energy[1]);
        tex->DrawLatex(x_start,y_start-=increment,label_cms[1]);
        tex->DrawLatex(x_start,y_start-=increment,label_n);
        tex->DrawLatex(x_start,y_start-=increment,label_pt);
        //tex->DrawLatex(x_start,y_start-=increment,label_cms[0]);
        tex->DrawLatex(x_start,y_start-=increment,label_pid[0]);

        cc2->cd();
        x_start = 0.139;
        y_start = 0.83+0.06;
        tex->DrawLatex(x_start,y_start-=increment,label_energy[0]);
        tex->DrawLatex(x_start,y_start-=increment,label_energy[1]);
        tex->DrawLatex(x_start,y_start-=increment,label_cms[1]);
        tex->DrawLatex(x_start,y_start-=increment,label_n);
        tex->DrawLatex(x_start,y_start-=increment,label_pt);
        tex_pid->DrawLatex(x_start,y_start-=increment,label_pid[1]);
        //tex->DrawLatex(0.22,0.70,label_cms[0]);
    }

    cc1->Print("Image/MassFitV0/massfit_ks.pdf");
    cc2->Print("Image/MassFitV0/massfit_la.pdf");
}
