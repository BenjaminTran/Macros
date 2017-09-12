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

void V0MassFit()
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

    TH1D* massks;
    TH1D* massla;
    TH2D* MassKs;
    TH2D* MassLa;
    TH3D* MassKs3D;
    TH3D* MassLa3D;
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
    std::vector<double> pks = {3,4, 5,6, 7,8, 9,10, 11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,70, 71,85};//, 86,100};//, 201,250, 251,300};
    std::vector<double> pla = {0,0, 0,0, 0,0, 9,10, 11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,70, 71,85};//, 86,100};//, 201,250, 251,300};

    TCanvas* Composite_Ks = new TCanvas("Composite_Ks","",1000,1200);
    Composite_Ks->Divide(3,5);

    TCanvas* Composite_La = new TCanvas("Composite_La","",1000,1200);
    Composite_La->Divide(3,4);

    //File Creation
    myfile.open("V0PeakParam.txt");
    //TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/Ksla/kslaMassPtJL1.root");
    //TFile* file_ks = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/Composites/V0CasMassPtPD11_16.root");
    TFile* file_ks = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/All/MCMassPtTotal_08_23_2017.root"); //Gen
    //TFile* file_la = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/La/LaMassPtPD1_6JL17_22.root");
    TFile* file_la = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/All/MCMassPtTotal_08_23_2017.root"); //Gen


    //MassKs = (TH2D*)file_ks->Get("MassPt/KsMassPt");
    //MassLa = (TH2D*)file_la->Get("LaMassPt/LaMassPt");
    MassKs3D = (TH3D*)file_ks->Get("MassPtRapidityMC/KsMassPtRap");
    MassLa3D = (TH3D*)file_la->Get("MassPtRapidityMC/LaMassPtRap");
    MassKs = (TH2D*)MassKs3D->Project3D("yx");
    MassLa = (TH2D*)MassLa3D->Project3D("yx");

    //Fit
    int pkscounter  = 0; //for correct bin counting
    int placounter  = 0;
    int hbincounter = 1;
    for(unsigned i=0; i<pks.size(); i++)
    //for(unsigned i=0; i<10; i++)
    {
        //i = pks.size()-2;
        int index = (i+2)/2;
        if(pla[i] == 0) lambda = false;
        else lambda = true;
        massks = (TH1D*)MassKs->ProjectionX("massks", pks[i],pks[i+1]);
        massla = (TH1D*)MassLa->ProjectionX("massla", pla[i],pla[i+1]);

        TCanvas* cc1 = new TCanvas("cc1","cc1",1000,450);
        cc1->Divide(2,1);

        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextFont(62);
        tex->SetTextSize(0.05);

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

        mass_ks.push_back(mean_ks);
        std_ks.push_back(rms_true_ks);
        fsig_ks.push_back(Fsig_ks);
        covQual_ks.push_back(covQuality_ks);

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
        double chi2_ks = xframe_ks->chiSquare("sum","data",7);

        TLine* t1_ks = new TLine(mean.getVal() - 2*rms_true_ks, 0, mean.getVal() - 2*rms_true_ks, gPad->GetUymax());
        TLine* t2_ks = new TLine(mean.getVal() + 2*rms_true_ks, 0, mean.getVal() + 2*rms_true_ks, gPad->GetUymax());
        t1_ks->SetLineStyle(2);
        t1_ks->SetLineColor(kGreen);
        t2_ks->SetLineStyle(2);
        t2_ks->SetLineColor(kGreen);
        t1_ks->Draw("same");
        t2_ks->Draw("same");
        double xstart_ks = 0.65;
        double ystart_ks = 0.85;
        double xpos = xstart_ks;
        double ypos = ystart_ks;
        double increment = 0.07;
        os << (pks[i]-1)/10 << " < P_{t} < " << pks[i+1]/10  << " GeV";
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

        Composite_Ks->cd(index);
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
            os << "CMS pPb";
            tex->SetTextSize(0.06);
            tex->DrawLatex(0.18,ypos-increment,os.str().c_str());
            tex->SetTextSize(0.05);
            os.str(std::string());
            os << "185 #leq N_{trk}^{offline} < 250";
            tex->DrawLatex(0.18,ypos-2*increment,os.str().c_str());
            os.str(std::string());
        }
        os << (pks[i]-1)/10 << " < P_{t} < " << pks[i+1]/10 << " GeV";
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

        if(lambda)
        {
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
            RooChi2Var chi2_laVar("chi2_laVar","chi2",sum,data);

            double covQuality_la = r_la->covQual();
            double mean_la = mean.getVal();

            double gaus1F_la = sig1.getVal();
            double gaus2F_la = sig2.getVal();
            double polyF_la  = poly.getVal();

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
            double chi2_la = xframe_la->chiSquare("sum","data",7);

            TLine* t1_la = new TLine(mean.getVal() - 2*rms_true_la, 0, mean.getVal() - 2*rms_true_la, gPad->GetUymax());
            TLine* t2_la = new TLine(mean.getVal() + 2*rms_true_la, 0, mean.getVal() + 2*rms_true_la, gPad->GetUymax());
            t1_la->SetLineStyle(2);
            t1_la->SetLineColor(kGreen);
            t2_la->SetLineStyle(2);
            t2_la->SetLineColor(kGreen);
            t1_la->Draw("same");
            t2_la->Draw("same");
            double xstart_la = 0.60;
            double ystart_la = 0.85;
            xpos = xstart_la;
            ypos = ystart_la;
            increment = 0.07;
            os << (pks[i]-1)/10 << " < P_{t} < " << pks[i+1]/10 << " GeV";
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

            Composite_La->cd(index-3);
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
            if(i==pks.size()-2)
            {
                xpos = 0.65;
            }
            if(i==6)
            {
                os << "CMS pPb";
                tex->SetTextSize(0.06);
                tex->DrawLatex(0.18,ypos-increment,os.str().c_str());
                tex->SetTextSize(0.05);
                os.str(std::string());
                os << "185 #leq N_{trk}^{offline} < 250";
                tex->DrawLatex(0.18,ypos-2*increment,os.str().c_str());
                os.str(std::string());
            }
            os << (pks[i]-1)/10 << " < P_{t} < " << pks[i+1]/10 << " GeV";
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
        }

        if(i==0) cc1->Print("V0MassFitInd.pdf(","pdf");
        else if(i < pks.size() - 2) cc1->Print("V0MassFitInd.pdf","pdf");
        else cc1->Print("V0MassFitInd.pdf)","pdf");
        i++; //to access correct bins
        hbincounter++;
    }
    Composite_Ks->Print("KsMassFitComposite.pdf");
    Composite_La->Print("LaMassFitComposite.pdf");

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
