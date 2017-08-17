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
    TGaxis::SetMaxDigits(2);
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

    //std::vector<double> pxi = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,72, 73,85, 86,100, 101,200, 201,300};
    //std::vector<double> pxi = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,100, 101,200, 201,300};
    std::vector<double> pxi = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,72, 73,100, 101,200, 201,300};

    TCanvas* cc1 = new TCanvas("cc1","cc1",1000,1200);
    cc1->Divide(3,4);

    //File Creation
    myfile.open("XiPeakParam.txt");
    //TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/CascadeV2pPb/RootFiles/Flow/CasCutLoose/CasCutLooseJL40.root");
    //TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/Thesis/XiAnalysisCorrelationPtCut8TeVPD1_4_ForFinal.root");
    TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/Composites/V0CasMassPtPD11_16.root");

    //MassXi = (TH2D*)file->Get("XiMassPt/MassPt");
    //MassXi = (TH2D*)file->Get("xiCorrelation/MassPt");
    MassXi = (TH2D*)file->Get("MassPt/MassPt");

    //Fit
    int pxicounter = 0; //for correct bin counting
    int hbincounter =1; //histogram bin counting
    for(unsigned i=0; i<pxi.size(); i++)
    {
        TCanvas* cc2 = new TCanvas("cc2","",600,450);
        int index = (i+2)/2;
        massxi = (TH1D*)MassXi->ProjectionX("massxi", pxi[i],pxi[i+1]);

        gStyle->SetOptTitle(kFALSE);

        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextFont(42);

        RooRealVar x("x","mass",1.26,1.4);
        RooDataHist data("data","dataset",x,massxi);
        RooRealVar mean("mean","mean",1.32,1.29,1.33);
        RooRealVar sigma1("sigma1","sigma1",0.004,0.001,0.04);
        RooRealVar sigma2("sigma2","sigma2",0.006,0.001,0.04);
        RooRealVar sig1("sig1","signal1",1200,0,1000000000);
        RooRealVar sig2("sig2","signal2",1000,0,1000000000);
        RooRealVar qsig("qsig","qsig",600,0,1000000000);
        RooRealVar alpha("alpha","alpha",1,0,10);
        RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
        RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
        RooGenericPdf background("background", "x - (1.115683 + 0.13957018)^alpha", RooArgList(x,alpha));
        RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,background),RooArgList(sig1,sig2,qsig));

        x.setRange("cut",1.285,1.375);

        RooFitResult* r_xi = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));
        RooChi2Var chi2_xiVar("chi2_xi","chi2",sum,data);

        double chi2_xi = chi2_xiVar.getVal();
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

        RooPlot* xframe_ = x.frame(150);
        xframe_->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");
        xframe_->GetYaxis()->SetTitle("Candidates / 0.001 GeV");
        xframe_->GetXaxis()->CenterTitle(1);
        xframe_->GetYaxis()->CenterTitle(1);
        xframe_->GetXaxis()->SetNdivisions(507);
        xframe_->GetXaxis()->SetTitleSize(0.06);
        xframe_->GetYaxis()->SetTitleSize(0.055);
        xframe_->GetYaxis()->SetTitleOffset(1);
        xframe_->GetXaxis()->SetTitleOffset(1.15);
        //xframe_->GetXaxis()->SetLabelSize(xframe_->GetXaxis()->GetLabelSize()*2.0);
        xframe_->GetYaxis()->SetLabelSize(0.12);
        xframe_->GetXaxis()->SetLabelSize(0.13);
        //xframe_->GetYaxis()->SetLabelSize(xframe_->GetYaxis()->GetLabelSize()*2.0);
        data.plotOn(xframe_,Name("data"));
        sum.plotOn(xframe_,Name("sum"),NormRange("cut"),LineWidth(1),LineColor(kRed));
        sum.plotOn(xframe_,Components(background),NormRange("cut"),LineStyle(kDashed),LineWidth(1),LineColor(kBlue));
        cc1->cd(index);
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        xframe.push_back(xframe_);
        xframe_->Draw();
        cc1->Update();

        TLine* t1 = new TLine(mean.getVal() - 2*rms_true_xi, 0, mean.getVal() - 2*rms_true_xi, gPad->GetUymax());
        TLine* t2 = new TLine(mean.getVal() + 2*rms_true_xi, 0, mean.getVal() + 2*rms_true_xi, gPad->GetUymax());
        t1->SetLineStyle(2);
        t1->SetLineColor(kGreen);
        t2->SetLineStyle(2);
        t2->SetLineColor(kGreen);
        t1->Draw("same");
        t2->Draw("same");
        os << "P_{t} Bin: " << (pxi[i]-1)/10 << " - " << pxi[i+1]/10;
        tex->DrawLatex(0.15,0.8,os.str().c_str());
        os.str(std::string());
        os << "Mean: " << mean_xi;
        tex->DrawLatex(0.15,0.75,os.str().c_str());
        os.str(std::string());
        os << "#sigma :" << rms_true_xi;
        tex->DrawLatex(0.15,0.70,os.str().c_str());
        os.str(std::string());
        os << "CovQual: " << covQual;
        tex->DrawLatex(0.15,0.65,os.str().c_str());
        os.str(std::string());
        osYield << "Yield: " << std::setprecision(2) << Yield_xi;
        tex->DrawLatex(0.15,0.60,osYield.str().c_str());
        osYield.str(std::string());
        os.str(std::string());
        //os << "#chi^{2}/ndf: " << chi2_xi;
        //tex->DrawLatex(0.15,0.60,os.str().c_str());
        //os.str(std::string());


        tex->SetTextSize(tex->GetTextSize()*0.95);

        cc2->cd();
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        xframe_->GetXaxis()->SetTitleOffset(1);
        xframe_->GetXaxis()->SetTitleSize(xframe_->GetXaxis()->GetTitleSize()*0.8);
        //xframe_->GetYaxis()->SetTitleSize(xframe_->GetYaxis()->GetTitleSize()*1.3);
        xframe_->GetXaxis()->SetLabelSize(xframe_->GetXaxis()->GetLabelSize()*0.5);
        xframe_->GetYaxis()->SetLabelSize(xframe_->GetYaxis()->GetLabelSize()*0.5);
        xframe_->Draw();
        cc2->Update();

        t1->Draw("same");
        t2->Draw("same");
        os << "P_{t} Bin: " << (pxi[i]-1)/10 << " - " << pxi[i+1]/10;
        tex->DrawLatex(0.15,0.8,os.str().c_str());
        os.str(std::string());
        os << "Mean: " << mean_xi;
        tex->DrawLatex(0.15,0.75,os.str().c_str());
        os.str(std::string());
        os << "#sigma :" << rms_true_xi;
        tex->DrawLatex(0.15,0.70,os.str().c_str());
        os.str(std::string());
        os << "CovQual: " << covQual;
        tex->DrawLatex(0.15,0.65,os.str().c_str());
        os.str(std::string());
        osYield << "Yield: " << std::setprecision(2) << Yield_xi;
        tex->DrawLatex(0.15,0.60,osYield.str().c_str());
        osYield.str(std::string());
        //os << "#chi^{2}/ndf: " << chi2_xi;
        //tex->DrawLatex(0.15,0.60,os.str().c_str());
        //os.str(std::string());

        hbincounter++;
        if(i==0) cc2->Print("XiMassFitInd.pdf(","pdf");
        else if(i < pxi.size() - 2) cc2->Print("XiMassFitInd.pdf","pdf");
        else cc2->Print("XiMassFitInd.pdf)","pdf");
        i++; //to access correct bins
    }
    cc1->Print("XiMassFitComposite.pdf");

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
