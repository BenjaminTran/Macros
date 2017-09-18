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

void OmMassFit()
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

    //std::vector<double> pxi = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,72, 73,85, 86,100, 101,200, 201,300};
    //std::vector<double> pxi = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,100, 101,200};//, 201,300};
    //std::vector<double> pxi = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,72, 73,100, 101,200};//, 201,300};
    std::vector<double> pxi = {11,15, 16,19, 20,23, 24,27, 28,33, 34,41, 42,50 ,51,60, 61,80};//, 81,100, 101,200};//, 201,300};

    TCanvas* cc1 = new TCanvas("cc1","cc1",1200,1200);
    cc1->Divide(3,3);

    //File Creation
    myfile.open("OmPeakParam.txt");
    //TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/CascadeV2pPb/RootFiles/Flow/CasCutLoose/CasCutLooseJL40.root");
    //TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/Thesis/XiAnalysisCorrelationPtCut8TeVPD1_4_ForFinal.root");
    //TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/Composites/V0CasMassPtPD5JL12.root"); //only one PD
    TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/Omega/OmMassPt.root"); //only one PD

    //MassXi = (TH2D*)file->Get("XiMassPt/MassPt");
    //MassXi = (TH2D*)file->Get("xiCorrelation/MassPt");
    MassXi = (TH2D*)file->Get("hom_Default");

    //Fit
    int pxicounter = 0; //for correct bin counting
    int hbincounter =1; //histogram bin counting
    //for(unsigned i=0; i<6; i++)
    for(unsigned i=0; i<pxi.size(); i++)
    {
        TCanvas* cc2 = new TCanvas("cc2","",600,450);
        int index = (i+2)/2;
        massxi = (TH1D*)MassXi->ProjectionX("massxi", pxi[i],pxi[i+1])->Rebin(2);

        gStyle->SetOptTitle(kFALSE);

        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextFont(62);
        tex->SetTextSize(0.04);
        //tex->SetTextAlign(10);

        RooRealVar x("x","mass",1.62,1.73);
        RooPlot* xframe_ = x.frame(50);
        xframe_->GetXaxis()->SetTitle("#Lambda K Invariant mass (GeV)");
        xframe_->GetYaxis()->SetTitle("Candidates / 0.002 GeV");
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
        RooRealVar mean("mean","mean",1.67,1.29,1.99);
        RooRealVar sigma1("sigma1","sigma1",0.004,0.001,0.04);
        RooRealVar sigma2("sigma2","sigma2",0.006,0.001,0.04);
        RooRealVar sig1("sig1","signal1",12,0,1000000000);
        RooRealVar sig2("sig2","signal2",10,0,1000000000);
        RooRealVar qsig("qsig","qsig",500,0,1000000000);
        RooRealVar alpha("alpha","alpha",0.001,-1,10);
        RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
        RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
        //RooRealVar a("a","a",0,-100000,100000);
        //RooRealVar b("b","b",0,-100000,100000);
        //RooRealVar cp("cp","cp",0,-100000,100000);
        //RooRealVar d("d","d",0,-100000,100000);
        //RooPolynomial background("poly","poly",x,RooArgList(a,b,cp,d));
        //RooPolynomial background("poly","poly",x,RooArgList(a,b,cp));
        //RooRealVar qsig("polysig","polysig",10,0,1000000000);
        RooGenericPdf background("background", "x - (1.115683 + 0.493677)^alpha", RooArgList(x,alpha));
        RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,background),RooArgList(sig1,sig2,qsig));

        x.setRange("cut",1.645,1.7);

        RooFitResult* r_xi = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));
        //RooChi2Var chi2_xiVar("chi2_xi","chi2",sum,data);

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

        sum.plotOn(xframe_,Name("sum"),NormRange("cut"),LineWidth(1),LineColor(kBlue));
        sum.plotOn(xframe_,Components(background),NormRange("cut"),LineStyle(kDashed),LineWidth(1),LineColor(kBlue));
        cc1->cd(index);
        gPad->SetBottomMargin(0.15); //gives more space for titles
        gPad->SetLeftMargin(0.15);
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        xframe.push_back(xframe_);
        xframe_->Draw();
        cc1->Update();
        double chi2_xi = xframe_->chiSquare("sum","data",4);

        TLine* t1 = new TLine(mean.getVal() - 2*rms_true_xi, 0, mean.getVal() - 2*rms_true_xi, gPad->GetUymax());
        TLine* t2 = new TLine(mean.getVal() + 2*rms_true_xi, 0, mean.getVal() + 2*rms_true_xi, gPad->GetUymax());
        t1->SetLineStyle(2);
        t1->SetLineColor(kGreen);
        t2->SetLineStyle(2);
        t2->SetLineColor(kGreen);
        t1->Draw("same");
        t2->Draw("same");
        double xpos = 0.60;
        double ypos = 0.85;
        double increment = 0.07;
        if(i==0)
        {
            os << "CMS pPb";
            tex->SetTextSize(0.06);
            tex->DrawLatex(0.17,ypos-increment,os.str().c_str());
            tex->SetTextSize(0.04);
            os.str(std::string());
            os << "185 #leq N_{trk}^{offline} < 220";
            tex->DrawLatex(0.17,ypos-2*increment,os.str().c_str());
            os.str(std::string());
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
        //os << "#chi^{2}/ndf: " << std::setprecision(3) << chi2_xi << std::setprecision(6);
        //tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        //os.str(std::string());
        osYield << "Yield: " << std::setprecision(2) << Yield_xi;
        tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
        osYield.str(std::string());
        os.str(std::string());


        //tex->SetTextSize(tex->GetTextSize()*0.95);

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
        tex->DrawLatex(xpos,0.8,os.str().c_str());
        os.str(std::string());
        os << "Mean: " << mean_xi;
        tex->DrawLatex(xpos,0.75,os.str().c_str());
        os.str(std::string());
        os << "#sigma :" << rms_true_xi;
        tex->DrawLatex(xpos,0.70,os.str().c_str());
        os.str(std::string());
        os << "CovQual: " << covQual;
        tex->DrawLatex(xpos,0.65,os.str().c_str());
        os.str(std::string());
        osYield << "Yield: " << std::setprecision(2) << Yield_xi;
        tex->DrawLatex(xpos,0.60,osYield.str().c_str());
        osYield.str(std::string());
        //os << "#chi^{2}/ndf: " << chi2_xi;
        //tex->DrawLatex(xpos,0.60,os.str().c_str());
        //os.str(std::string());

        hbincounter++;
        if(i==0) cc2->Print("OmMassFitInd.pdf(","pdf");
        else if(i < pxi.size() - 2) cc2->Print("OmMassFitInd.pdf","pdf");
        else cc2->Print("OmMassFitInd.pdf)","pdf");
        i++; //to access correct bins
    }
    cc1->Print("OmMassFitComposite.pdf");
    cc1->Print("OmMassFitComposite.png");

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
