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

void Efficiency()
{
    //Initializers
    std::string FileName = "";
    FileName = "V0Efficiency.txt";
    using namespace RooFit;
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetMarkerSize(0.5);
    TGaxis::SetMaxDigits(3);
    RooMsgService::instance().setStreamStatus(0,kFALSE);
    RooMsgService::instance().setStreamStatus(1,kFALSE);
    std::ostringstream os;
    std::ostringstream osYield;
    std::ofstream myfile;

    TH1D* massks;
    TH1D* massla;
    //TH1D* massks_gen;
    //TH1D* massla_gen;
    TH2D* MassKs;
    TH2D* MassLa;
    TH2D* MassKs_gen;
    TH2D* MassLa_gen;

    TCanvas* Composite_Ks[10];
    TCanvas* Composite_La[10];
    //TCanvas* Composite_Ks_Gen[7];
    //TCanvas* Composite_La_Gen[7];
    TCanvas* cc1 = new TCanvas("cc1","cc1",1000,450);
    cc1->Divide(2,1);
    bool lambda;

    std::vector<RooPlot*> Xframe_Ks;
    // rapidity index, values by pt bin
    std::map<int,std::vector<double> > mass_ks_;
    std::map<int,std::vector<double> > std_ks_;
    std::map<int,std::vector<double> > fsig_ks_;
    std::map<int,std::vector<double> > covQual_ks_;
    std::map<int,std::vector<double> > yield_ks_;
    std::map<int,std::vector<double> > yield_ks_gen;

    std::vector<RooPlot*> Xframe_La;
    std::map<int,std::vector<double> > mass_la_;
    std::map<int,std::vector<double> > std_la_;
    std::map<int,std::vector<double> > fsig_la_;
    std::map<int,std::vector<double> > covQual_la_;
    std::map<int,std::vector<double> > yield_la_;
    std::map<int,std::vector<double> > yield_la_gen;

    std::map<int,std::vector<double> > efficiency_ks;
    std::map<int,std::vector<double> > efficiency_la;

    //double pks[] = {3,4, 5,6, 7,8, 9,10, 11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,90, 91,120};
    //double pla[] = {0,0, 0,0, 0,0, 9,10, 11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,90, 91,120};
    //std::vector<double> pks = {1,2, 3,4, 5,6, 7,8, 9,10, 11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,70};//, 71,85, 86,100, 101,150};//, 151,200};//, 201,250, 251,300};
    //std::vector<double> pla = {1,2, 3,4, 5,6, 7,8, 9,10, 11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,70};//, 71,85, 86,100, 101,150};//, 151,200};//, 201,250, 251,300};
    //std::vector<double> pks = {1,1, 2,2, 3,3, 4,4, 5,5, 6,6, 7,7, 8,8, 9,9, 10,11, 12,13, 14,15, 16,17, 18,19, 20,22, 23,25, 26,28, 29,31, 32,34, 35,37, 38,40, 41,44, 45,50, 51,60, 61,70, 71,100};
    //std::vector<double> pla = {1,5, 6,8, 9,11, 12,14, 15,17, 18,20, 21,23, 24,26, 27,29, 30,32, 33,36, 37,40, 41,44, 45,48, 49,52, 53,56, 57,64, 65,80, 81,100};
    //std::vector<double> rap_bin = {2,4, 5,7, 8,10, 11,12, 13,15, 16,18, 19,21};
    std::vector<double> pks = {0,5,6,7,8,9,11,13,15,17,19,22,25,28,31,34,37,40,44,50,60,70,100};
    std::vector<double> pla = {0,8,11,14,17,20,23,26,29,32,36,40,44,48,52,56,64,80,100};
    std::vector<double> rap_bin = {2,6, 7,10, 11,14, 15,21};//, 19,21, 16,18, 19,21};
    int numKsBins = pks.size()-1;
    int numLaBins = pla.size()-1;
    int numRapBins = rap_bin.size()/2;
    //int numRapBins = 1;
    int pages=8;

    for(int i=0; i<pages; i++)
    {
        Composite_Ks[i] = new TCanvas(Form("Composite_Ks_%d",i),Form("Composite_Ks_RapBin_%d",i),1000,1200);
        Composite_Ks[i]->Divide(3,5);

        Composite_La[i] = new TCanvas(Form("Composite_La_%d",i),Form("Composite_La_RapBin_%d",i),1000,1200);
        Composite_La[i]->Divide(3,5);

        //Composite_Ks_Gen[i] = new TCanvas(Form("Composite_Ks_Gen_%d",i),Form("Composite_Ks_Gen_RapBin_%d"),1000,1200);
        //Composite_Ks_Gen[i]->Divide(4,4);

        //Composite_La_Gen[i] = new TCanvas(Form("Composite_La_Gen_%d",i),Form("Composite_La_Gen_RapBin_%d"),1000,1200);
        //Composite_La_Gen[i]->Divide(4,4);
    }

    //File Creation
    myfile.open(FileName.c_str());
    TFile* f1 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/MCMassPtTotal_08_23_2017.root");

    TH3D* LaMassPtRap = (TH3D*)f1->Get("MassPtRapidityMC/LaMassPtRap");
    TH3D* KsMassPtRap = (TH3D*)f1->Get("MassPtRapidityMC/KsMassPtRap");

    MassKs = (TH2D*)KsMassPtRap->Project3D("yx");
    MassLa = (TH2D*)LaMassPtRap->Project3D("yx");

    TH3D* LaMassPtRap_Gen = (TH3D*)f1->Get("MassPtRapidityMC/LaMassPtRap_Gen");
    TH3D* KsMassPtRap_Gen = (TH3D*)f1->Get("MassPtRapidityMC/KsMassPtRap_Gen");

    MassKs_gen = (TH2D*)KsMassPtRap_Gen->Project3D("yx");
    MassLa_gen = (TH2D*)LaMassPtRap_Gen->Project3D("yx");

    //Fit
    int hcounter_ks = -1;
    int hcounter_la = -1;
    int index = 1;
    double xstart_ks = 0.65;
    double ystart_ks = 0.85;
    double xstart_la = 0.60;
    double ystart_la = 0.85;
    double xpos = xstart_ks;
    double ypos = ystart_ks;
    double increment = 0.07;
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(62);
    tex->SetTextSize(0.05);
    //for(unsigned j=0; j<rap_bin.size(); j+=2)
    for(int j=0; j<numRapBins*2; j+=2)
    {
        for(int i=0; i<numKsBins; i++)
        //for(unsigned i=0; i<60; i+=2)
        {
            cout << i << endl;
            if((double)i/15 == 1 || i==0)
            {
                hcounter_ks++;
                index = 1;
            }
            lambda=true;

            massks = (TH1D*)KsMassPtRap->ProjectionX(Form("massks_%d",(int)(j*pks.size()+i)/2), pks[i]+1,pks[i+1],rap_bin[j],rap_bin[j+1]);

            //tex->SetTextSize(tex->GetTextSize()*0.95);

            //kshort
            RooRealVar x("x","mass",0.43,0.565);
            RooDataHist data("data","dataset",x,massks);
            RooRealVar mean("mean","mean",0.50,0.49,0.51);
            RooRealVar sigma1("sigma1","sigma1",0.01,0.001,0.04);
            RooRealVar sigma2("sigma2","sigma2",0.01,0.001,0.04);
            RooRealVar sig1("sig1","signal1",10,0,1000000000);
            RooRealVar sig2("sig2","signal2",10,0,1000000000);
            RooRealVar a("a","a",0,-100000,100000);
            RooRealVar b("b","b",0,-100000,100000);
            RooRealVar cp("cp","cp",0,-100000,100000);
            RooRealVar d("d","d",0,-100000,100000);
            RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
            RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
            RooPolynomial poly("poly","poly",x,RooArgList(a,b,cp,d));
            RooRealVar polysig("polysig","polysig",10,0,1000000000);
            RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));

            x.setRange("cut",0.44,0.56);

            RooFitResult* r_ks = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));

            double covQuality_ks = r_ks->covQual();
            //double covQuality_ks = 1;
            double mean_ks = mean.getVal();

            double gaus1F_ks = sig1.getVal();
            double gaus2F_ks = sig2.getVal();
            double polyF_ks  = poly.getVal();

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

            //mass_ks_[j/2].push_back(mean_ks);
            //std_ks_[j/2].push_back(rms_true_ks);
            //fsig_ks_[j/2].push_back(Fsig_ks);
            //covQual_ks_[j/2].push_back(covQuality_ks);
            yield_ks_[j/2].push_back(yield_ks);

            cout << "Yield (ks): "<< yield_ks << endl;
            cout << "Fsig (ks): " << Fsig_ks << endl;
            cout << "std (ks): " << rms_true_ks << endl;
            cout << "mass (ks): " << mean_ks << endl;
            cout << "covQual (ks)" << covQuality_ks << endl;

            RooPlot* xframe_ks = x.frame(270);
            xframe_ks->GetXaxis()->SetTitle("Invariant mass (GeV)");
            xframe_ks->GetYaxis()->SetTitle("Candidates / 0.0005 GeV");
            xframe_ks->GetXaxis()->CenterTitle(1);
            xframe_ks->GetYaxis()->CenterTitle(1);
            xframe_ks->GetXaxis()->SetTickSize(0.02);
            xframe_ks->GetYaxis()->SetTickSize(0.02);
            xframe_ks->GetXaxis()->SetNdivisions(407);
            xframe_ks->GetYaxis()->SetNdivisions(410);
            xframe_ks->GetXaxis()->SetTitleSize(0.06);
            xframe_ks->GetYaxis()->SetTitleSize(0.06);
            xframe_ks->GetYaxis()->SetTitleOffset(0.85);
            xframe_ks->GetXaxis()->SetTitleOffset(0.5);
            //xframe_ks->GetXaxis()->SetLabelSize(xframe_ks->GetXaxis()->GetLabelSize()*2.0);
            xframe_ks->GetYaxis()->SetLabelSize(0.1);
            xframe_ks->GetXaxis()->SetLabelSize(0.1);
            //xframe_ks->GetYaxis()->SetLabelSize(xframe_ks->GetYaxis()->GetLabelSize()*2.0);
            data.plotOn(xframe_ks,Name("data"));
            sum.plotOn(xframe_ks,Name("sum"),NormRange("cut"),LineWidth(1),LineColor(kBlue));
            sum.plotOn(xframe_ks,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(1),LineColor(kBlue));
            cc1->cd(1);
            gPad->SetTickx();
            gPad->SetTicky();
            Xframe_Ks.push_back(xframe_ks);
            xframe_ks->Draw();
            cc1->Update();
            //double chi2_ks = xframe_ks->chiSquare("sum","data",7);

            TLine* t1_ks = new TLine(mean.getVal() - 2*rms_true_ks, 0, mean.getVal() - 2*rms_true_ks, gPad->GetUymax());
            TLine* t2_ks = new TLine(mean.getVal() + 2*rms_true_ks, 0, mean.getVal() + 2*rms_true_ks, gPad->GetUymax());
            t1_ks->SetLineStyle(2);
            t1_ks->SetLineColor(kGreen);
            t2_ks->SetLineStyle(2);
            t2_ks->SetLineColor(kGreen);
            t1_ks->Draw("same");
            t2_ks->Draw("same");
            xpos = xstart_ks;
            ypos = ystart_ks;
            increment = 0.07;
            os << (pks[i])/10 << " < P_{t} < " << pks[i+1]/10  << " GeV";
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

            Composite_Ks[hcounter_ks]->cd(index);
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
            Composite_Ks[hcounter_ks]->Update();

            t1_ks->Draw("same");
            t2_ks->Draw("same");
            xpos = xstart_ks;
            ypos = ystart_ks;
            increment = 0.07;
            if(i==0)
            {
                os << "CMS pPb";
                tex->SetTextSize(0.06);
                tex->DrawLatex(0.18,ypos-increment,os.str().c_str());
                tex->SetTextSize(0.05);
                os.str(std::string());
                os << "185 #leq N_{trk}^{offline} < 250";
                tex->DrawLatex(0.18,ypos-2*increment,os.str().c_str());
                os.str(std::string());
                os << (-1.1 + (rap_bin[j]-1)/10) << " < y < " << (-1.1 + (rap_bin[j+1])/10);
                tex->DrawLatex(0.18,ypos-3*increment,os.str().c_str());
                os.str(std::string());
            }
            os << (pks[i])/10 << " < P_{t} < " << pks[i+1]/10 << " GeV";
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
            index++;
        }

        for(int i=0; i<numLaBins; i++)
        {
            if((double)i/15 == 1 || i==0)
            {
                hcounter_la++;
                index = 1;
            }

            massla = (TH1D*)LaMassPtRap->ProjectionX(Form("massla_%d",(int)(j*pla.size()+i)/2), pla[i]+1,pla[i+1],rap_bin[j],rap_bin[j+1]);

            //lambda
            RooRealVar x("x","mass",1.08,1.155);
            RooPlot* xframe_la = x.frame(160);
            RooDataHist data("data","dataset",x,massla);
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

            x.setRange("cut",1.1,1.14);

            RooFitResult* r_la = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));

            double covQuality_la = r_la->covQual();
            //double covQuality_la = 1;
            double mean_la = mean.getVal();

            double gaus1F_la = sig1.getVal();
            double gaus2F_la = sig2.getVal();
            double polyF_la  = poly.getVal();

            //set ranges for individual gaussian yield determination
            x.setRange("g1", mean.getVal() - 2*sigma1.getVal(), mean.getVal() + 2*sigma1.getVal());
            x.setRange("g2", mean.getVal() - 2*sigma2.getVal(), mean.getVal() + 2*sigma2.getVal());

            RooAbsReal* Intgaus1_yield_la_ = gaus1.createIntegral(x, x, "g1");
            RooAbsReal* Intgaus2_yield_la_ = gaus2.createIntegral(x, x, "g2");

            double gaus1_yield_la_ = gaus1F_la*Intgaus1_yield_la_->getVal();
            double gaus2_yield_la_ = gaus2F_la*Intgaus2_yield_la_->getVal();
            double gausTot_yield_la_ = gaus1_yield_la_ + gaus2_yield_la_;

            double rms_gaus1_sig_la = gaus1_yield_la_/gausTot_yield_la_;
            double rms_gaus2_sig_la = gaus2_yield_la_/gausTot_yield_la_;

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

            //mass_la_[j/2].push_back(mean_la);
            //std_la_[j/2].push_back(rms_true_la);
            //fsig_la_[j/2].push_back(Fsig_la);
            //covQual_la_[j/2].push_back(covQuality_la);
            //covQual_la_[j/2].push_back(1);
            yield_la_[j/2].push_back(yield_la);

            cout << "Yield (la):" << yield_la << endl;
            cout << "Fsig (la): " << Fsig_la << endl;
            cout << "std (la): " << rms_true_la << endl;
            cout << "covQual (la): " << covQuality_la << endl;
            cout << (j*pks.size() + i)/2 << endl;

            xframe_la->GetXaxis()->SetTitle("Invariant mass (GeV)");
            xframe_la->GetYaxis()->SetTitle("Candidates / 0.0005 GeV");
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
            sum.plotOn(xframe_la,Name("sum"),NormRange("cut"),LineWidth(1),LineColor(kRed));
            sum.plotOn(xframe_la,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(1),LineColor(kRed));
            cc1->cd(2);
            gPad->SetTickx();
            gPad->SetTicky();
            Xframe_La.push_back(xframe_la);
            xframe_la->Draw();
            cc1->Update();
            //double chi2_la = xframe_la->chiSquare("sum","data",7);

            TLine* t1_la = new TLine(mean.getVal() - 2*rms_true_la, 0, mean.getVal() - 2*rms_true_la, gPad->GetUymax());
            TLine* t2_la = new TLine(mean.getVal() + 2*rms_true_la, 0, mean.getVal() + 2*rms_true_la, gPad->GetUymax());
            t1_la->SetLineStyle(2);
            t1_la->SetLineColor(kGreen);
            t2_la->SetLineStyle(2);
            t2_la->SetLineColor(kGreen);
            t1_la->Draw("same");
            t2_la->Draw("same");
            xpos = xstart_la;
            ypos = ystart_la;
            increment = 0.07;
            os << (pla[i])/10 << " < P_{t} < " << pla[i+1]/10 << " GeV";
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

            Composite_La[hcounter_la]->cd(index);
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
            Composite_La[hcounter_la]->Update();

            t1_la->Draw("same");
            t2_la->Draw("same");
            xpos = xstart_la;
            ypos = ystart_la;
            increment = 0.07;
            if(i==numLaBins)
            {
                xpos = 0.65;
            }
            if(i==2)
            {
                os << "CMS pPb";
                tex->SetTextSize(0.06);
                tex->DrawLatex(0.18,ypos-increment,os.str().c_str());
                tex->SetTextSize(0.05);
                os.str(std::string());
                os << "185 #leq N_{trk}^{offline} < 250";
                tex->DrawLatex(0.18,ypos-2*increment,os.str().c_str());
                os.str(std::string());
                os << (-1.1 + (rap_bin[j]-1)/10) << " < y < " << (-1.1 + (rap_bin[j+1])/10);
                tex->DrawLatex(0.18,ypos-3*increment,os.str().c_str());
                os.str(std::string());
            }
            os << (pla[i])/10 << " < P_{t} < " << pla[i+1]/10 << " GeV";
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
            index++;
        }

        //if(i==0) cc1->Print("RECOV0MassFitInd.pdf(","pdf");
        //else if(i < pks.size() - 2) cc1->Print("RECOV0MassFitInd.pdf","pdf");
        //else cc1->Print("RECOV0MassFitInd.pdf)","pdf");

    }

    for(int j=0; j<numRapBins; j++)
    {
        Composite_Ks[j]->Print(Form("RECOKsMassFitComposite_%d.pdf",j));
        Composite_La[j]->Print(Form("RECOLaMassFitComposite_%d.pdf",j));
    }

    for(unsigned j=0; j<rap_bin.size(); j+=2)
    {
        for(unsigned i=0; i<numKsBins; i++)
        {
            TH1D* massks_gen = (TH1D*)KsMassPtRap_Gen->ProjectionX(Form("massks_gen_%d",(int)(j*pks.size()+i)/2), pks[i]+1,pks[i+1],rap_bin[j],rap_bin[j+1]);
            yield_ks_gen[j/2].push_back(massks_gen->GetBinContent(massks_gen->GetMaximumBin()));
        }
        for(unsigned i=0; i<numLaBins; i++)
        {
            TH1D* massla_gen = (TH1D*)LaMassPtRap_Gen->ProjectionX(Form("massla_gen_%d",(int)(j*pla.size()+i)/2), pla[i]+1,pla[i+1],rap_bin[j],rap_bin[j+1]);
            yield_la_gen[j/2].push_back(massla_gen->GetBinContent(massla_gen->GetMaximumBin()));
        }
    }
    //Composite_Ks_Gen->Print("Composite_Ks_Gen.pdf");
    //Composite_La_Gen->Print("Composite_La_Gen.pdf");

    //Calculate efficiency
    for(int j=0; j<numRapBins; j++)
    {
        for(unsigned i=0; i<yield_ks_[j].size(); i++)
        {
            efficiency_ks[j].push_back(yield_ks_[j][i]/yield_ks_gen[j][i]);
        }
        for(unsigned i=0; i<yield_la_[j].size(); i++)
        {
            efficiency_la[j].push_back(yield_la_[j][i]/yield_la_gen[j][i]);
        }
    }

    //Output
/*
    int pkscounter = 0;
    myfile << "KSHORT KSHORT KSHORT\n";
    for(unsigned i=0; i<mass_ks_[j].size(); i++)
    {
        cout <<  "====================" << endl;
        cout << "Pt Bin: " << (pks[pkscounter]-1)/10 << " - " << pks[pkscounter+1]/10 << endl;
        cout <<  "====================" << endl;
        cout << "Mass_ks: "   << mass_ks_[i]    << endl;
        cout << "Fsig_ks: "   << fsig_ks_[i]    << endl;
        cout << "std_ks_: "    << std_ks_[i]     << endl;
        cout << "covQual_ks_: " << covQual_ks_[i] << endl;

        myfile <<  "====================" << "\n";
        myfile << "Pt Bin: " << (pks[pkscounter]-1)/10 << " - " << pks[pkscounter+1]/10 << "\n";
        myfile <<  "====================" << "\n";
        myfile << "Mass_ks: "   << mass_ks_[i]    << "\n";
        myfile << "Fsig_ks: "   << fsig_ks_[i]    << "\n";
        myfile << "std_ks_: "    << std_ks_[i]     << "\n";
        myfile << "covQual_ks_: " << covQual_ks_[i] << "\n";
        myfile << "Yield:" << yield_ks_[i] << "\n";
        pkscounter+=2;
    }

    int placounter=0;
    myfile << "LAMBDA LAMBDA LAMBDA\n";
    for(unsigned i=0; i<mass_la_[j].size(); i++)
    {
        cout <<  "====================" << endl;
        cout << "Pt Bin: " << (pla[placounter]-1)/10 << " - " << pla[placounter+1]/10 << endl;
        cout <<  "====================" << endl;
        cout << "Mass_la: "   << mass_la_[i]    << endl;
        cout << "Fsig_la: "   << fsig_la_[i]    << endl;
        cout << "std_la_: "    << std_la_[i]     << endl;
        cout << "covQual_la_ " << covQual_la_[i] << endl;

        myfile <<  "====================" << "\n";
        myfile << "Pt Bin: " << (pla[placounter]-1)/10 << " - " << pla[placounter+1]/10 << "\n";
        myfile <<  "====================" << "\n";
        myfile << "Mass_la: "   << mass_la_[i]    << "\n";
        myfile << "Fsig_la: "   << fsig_la_[i]    << "\n";
        myfile << "std_la_: "    << std_la_[i]     << "\n";
        myfile << "covQual_la_ " << covQual_la_[i] << "\n";
        myfile << "Yield:" << yield_la_[i] << "\n";
        placounter+=2;
    }
    */

    //myfile << "Yield Gen Kshort\n";
    //for(unsigned i=0; i<yield_ks_gen.size(); i++)
        //myfile << yield_ks_gen[i] << "\n";

    //myfile << "Yield Gen Lambda\n";
    //for(unsigned i=0; i<yield_la_gen.size(); i++)
        //myfile << yield_la_gen[i] << "\n";


    myfile << "Efficiency Kshort\n";
    for(int j=0; j<numRapBins; j++)
    {
        myfile << "Rap_bin " << j << "\n";
        for(unsigned i=0; i<efficiency_ks[j].size(); i++)
            myfile << efficiency_ks[j][i] << "\n";
    }

    myfile << "Efficiency Lambda\n";
    for(int j=0; j<numRapBins; j++)
    {
        myfile << "Rap_bin " << j << "\n";
        for(unsigned i=0; i<efficiency_la[j].size(); i++)
            myfile << efficiency_la[j][i] << "\n";
    }

    const int rbinrap = 4;

    double Rebin_ks[pks.size()];
    double Rebin_la[pla.size()];
    double Rebin_rap[] = {-1.0,-0.5,-0.1,0.3,1.0};

    for(unsigned i=0; i<pks.size(); i++)
    {
        Rebin_ks[i] = pks[i]/10;
        cout << Rebin_ks[i] << ", ";
    }
    cout << endl;
    for(unsigned i=0; i<pla.size(); i++)
    {
        Rebin_la[i] = pla[i]/10;
        cout << Rebin_la[i] << ", ";
    }
    cout << endl;

    TH2D* Effhisto_ks = new TH2D("EffHistoKs","EffHistoKs",rbinrap,Rebin_rap,numKsBins,Rebin_ks);
    TH2D* Effhisto_la = new TH2D("EffHistoLa","EffHistoLa",rbinrap,Rebin_rap,numLaBins,Rebin_la);

    for(int j=0; j<numRapBins; j++)
    {
        for(unsigned i=0; i<efficiency_ks[j].size(); i++)
            Effhisto_ks->Fill(Rebin_rap[j],Rebin_ks[i]+0.01,efficiency_ks[j][i]);
    }

    for(int j=0; j<numRapBins; j++)
    {
        for(unsigned i=0; i<efficiency_la[j].size(); i++)
            Effhisto_la->Fill(Rebin_rap[j],Rebin_la[i]+0.01,efficiency_la[j][i]);
    }

    TCanvas* Eff = new TCanvas("Eff","",1200,800);
    Eff->Divide(2,1);
    Eff->cd(1);
    Effhisto_ks->Draw("Lego2");
    Eff->cd(2);
    Effhisto_la->Draw("Lego2");

    TFile histos("Effhisto.root","RECREATE");
    Effhisto_ks->Write();
    Effhisto_la->Write();

    std::vector<double> eta_trg = {-0.9,-0.8,-0.4,-0.2,-0.1,0.2,0.5,0.9};
    std::vector<double> pt_trg = {0.1,1.2, 3.4, 4.0, 6.7, 8.6, 1.5, 2.0};
    for(unsigned i=0; i<eta_trg.size(); i++)
    {
        cout << Effhisto_ks->GetBinContent(Effhisto_ks->FindBin(eta_trg[i],pt_trg[i])) << endl;
        cout << Effhisto_la->GetBinContent(Effhisto_la->FindBin(eta_trg[i],pt_trg[i])) << endl;
    }
}

void Plotter()
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 25;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);

    TFile* f = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/Efficiency/Effhisto.root");

    TH2D* effhistks = (TH2D*)f->Get("EffHistoKs");
    TH2D* effhistla = (TH2D*)f->Get("EffHistoLa");

    effhistks->SetTitle("Efficiency of K_{S}^{0}");
    effhistks->GetZaxis()->SetRangeUser(0,0.3);
    effhistks->GetZaxis()->SetNdivisions(404);
    effhistks->GetYaxis()->SetNdivisions(505);
    effhistks->GetXaxis()->SetNdivisions(505);
    effhistks->GetYaxis()->SetTitle("P_{T} (GeV/c)");
    effhistks->GetXaxis()->SetTitle("y");
    effhistks->GetYaxis()->CenterTitle(1);
    effhistks->GetXaxis()->CenterTitle(1);
    effhistks->GetXaxis()->SetTitleOffset(1.5);
    effhistks->SetStats(kFALSE);

    effhistla->SetTitle("Efficiency of #Lambda");
    effhistla->GetZaxis()->SetRangeUser(0,0.15);
    effhistla->GetZaxis()->SetNdivisions(404);
    effhistla->GetYaxis()->SetNdivisions(505);
    effhistla->GetXaxis()->SetNdivisions(505);
    effhistla->GetYaxis()->SetTitle("P_{T} (GeV/c)");
    effhistla->GetXaxis()->SetTitle("y");
    effhistla->GetYaxis()->CenterTitle(1);
    effhistla->GetXaxis()->CenterTitle(1);
    effhistla->GetXaxis()->SetTitleOffset(1.5);
    effhistla->SetStats(kFALSE);

    TCanvas* cc1 = new TCanvas("cc1","Efficiency", 1600,800);
    cc1->Divide(2,1);
    cc1->cd(1);
    effhistks->Draw("Lego2 FB");
    cc1->cd(2);
    effhistla->Draw("Lego2 FB");

    cc1->Print("Effhistos.pdf");
}
