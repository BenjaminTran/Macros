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
    //double pks[] = {3,4, 5,6, 7,8, 9,10, 11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,90, 91,120};
    //double pla[] = {0,0, 0,0, 0,0, 9,10, 11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,90, 91,120};
    std::vector<double> pks = {3,4, 5,6, 7,8, 9,10, 11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,70, 71,85, 86,100, 101,150, 151,200, 201,250, 251,300};
    std::vector<double> pla = {0,0, 0,0, 0,0, 9,10, 11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,70, 71,85, 86,100, 101,150, 151,200, 201,250, 251,300};

    TCanvas* Composite_Ks = new TCanvas("Composite_Ks","",1000,1600);
    Composite_Ks->Divide(3,6);

    TCanvas* Composite_La = new TCanvas("Composite_La","",1000,1200);
    Composite_La->Divide(3,5);

    //File Creation
    myfile.open("V0PeakParam.txt");
    //TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/Ksla/kslaMassPtJL1.root");
    TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/Composites/V0CasMassPtPD11_16.root");

    MassKs = (TH2D*)file->Get("MassPt/KsMassPt");
    MassLa = (TH2D*)file->Get("MassPt/LaMassPt");

    //Fit
    int pkscounter  = 0; //for correct bin counting
    int placounter  = 0;
    int hbincounter = 1;
    for(unsigned i=0; i<pks.size(); i++)
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
        tex->SetTextFont(42);

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
        RooChi2Var chi2_ksVar("chi2_ks","chi2",sum,data);

        double chi2_ks = chi2_ksVar.getVal();
        double covQuality_ks = r_ks->covQual();
        double mean_ks = mean.getVal();
        double rms_ks  = TMath::Sqrt(0.5*sigma1.getVal()*sigma1.getVal() + 0.5*sigma2.getVal()*sigma2.getVal());

        x.setRange("peak", mean.getVal() - 2*rms_ks, mean.getVal() + 2*rms_ks);

        double gaus1F_ks = sig1.getVal();
        double gaus2F_ks = sig2.getVal();
        double polyF_ks  = poly.getVal();

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
        std_ks.push_back(rms_ks);
        fsig_ks.push_back(Fsig_ks);
        covQual_ks.push_back(covQuality_ks);

        cout << "Yield (ks): "<< yield_ks << endl;
        cout << "Fsig (ks): " << Fsig_ks << endl;
        cout << "std (ks): " << rms_ks << endl;
        cout << "mass (ks): " << mean_ks << endl;
        cout << "covQual (ks)" << covQuality_ks << endl;

        RooPlot* xframe = x.frame(270);
        xframe->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");
        xframe->GetYaxis()->SetTitle("Candidates / 0.0005 GeV");
        xframe->GetXaxis()->CenterTitle(1);
        xframe->GetYaxis()->CenterTitle(1);
        xframe->GetXaxis()->SetNdivisions(507);
        xframe->GetXaxis()->SetTitleSize(0.05);
        xframe->GetYaxis()->SetTitleSize(0.05);
        xframe->GetYaxis()->SetTitleOffset(1.1);
        xframe->GetXaxis()->SetTitleOffset(1);
        xframe->GetYaxis()->SetLabelSize(0.05);
        xframe->GetXaxis()->SetLabelSize(0.05);
        data.plotOn(xframe,Name("data"));
        sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(1),LineColor(kBlue));
        sum.plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(1),LineColor(kBlue));
        cc1->cd(1);
        gPad->SetTickx();
        gPad->SetTicky();
        Xframe_Ks.push_back(xframe);
        xframe->Draw();

        os << "Pt Bin: " << (pks[i]-1)/10 << " - " << pks[i+1]/10;
        tex->DrawLatex(0.15,0.8,os.str().c_str());
        os.str(std::string());
        os << "Mean: " << mean_ks;
        tex->DrawLatex(0.15,0.75,os.str().c_str());
        os.str(std::string());
        os << "#sigma :" << rms_ks;
        tex->DrawLatex(0.15,0.70,os.str().c_str());
        os.str(std::string());
        os << "CovQual: " << covQuality_ks;
        tex->DrawLatex(0.15,0.65,os.str().c_str());
        os.str(std::string());
        osYield << "Yield: " << std::setprecision(2) << yield_ks;
        tex->DrawLatex(0.15,0.60,osYield.str().c_str());
        osYield.str(std::string());
        //os << "#chi^{2}/ndf: " << chi2_ks;
        //tex->DrawLatex(0.15,0.60,os.str().c_str());
        //os.str(std::string());

        Composite_Ks->cd(index);
        gPad->SetTickx();
        gPad->SetTicky();
        xframe->GetXaxis()->SetNdivisions(507);
        xframe->GetXaxis()->SetTitleSize(0.05);
        xframe->GetYaxis()->SetTitleSize(0.05);
        xframe->GetYaxis()->SetTitleOffset(1.1);
        xframe->GetXaxis()->SetTitleOffset(1);
        xframe->GetYaxis()->SetLabelSize(0.045);
        xframe->GetXaxis()->SetLabelSize(0.045);
        xframe->Draw();

        os << "Pt Bin: " << (pks[i]-1)/10 << " - " << pks[i+1]/10;
        tex->DrawLatex(0.15,0.8,os.str().c_str());
        os.str(std::string());
        os << "Mean: " << mean_ks;
        tex->DrawLatex(0.15,0.75,os.str().c_str());
        os.str(std::string());
        os << "#sigma :" << rms_ks;
        tex->DrawLatex(0.15,0.70,os.str().c_str());
        os.str(std::string());
        os << "CovQual: " << covQuality_ks;
        tex->DrawLatex(0.15,0.65,os.str().c_str());
        os.str(std::string());
        osYield << "Yield: " << std::setprecision(2) << yield_ks;
        tex->DrawLatex(0.15,0.60,osYield.str().c_str());
        osYield.str(std::string());
        //os << "#chi^{2}/ndf: " << chi2_ks;
        //tex->DrawLatex(0.15,0.60,os.str().c_str());
        //os.str(std::string());

        if(lambda)
        {
            //lambda
            RooRealVar x("x","mass",1.08,1.155);
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

            x.setRange("cut",1.09,1.14);

            RooFitResult* r_la = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));
            RooChi2Var chi2_laVar("chi2_laVar","chi2",sum,data);

            double chi2_la = chi2_laVar.getVal();
            double covQuality_la = r_la->covQual();
            double mean_la = mean.getVal();
            double rms_la  = TMath::Sqrt(0.5*sigma1.getVal()*sigma1.getVal() + 0.5*sigma2.getVal()*sigma2.getVal());

            x.setRange("peak", mean.getVal() - 2*rms_la, mean.getVal() + 2*rms_la);

            double gaus1F_la = sig1.getVal();
            double gaus2F_la = sig2.getVal();
            double polyF_la  = poly.getVal();

            RooAbsReal* Intgaus1_la = gaus1.createIntegral(x, x,  "peak");
            RooAbsReal* Intgaus2_la = gaus2.createIntegral(x, x, "peak");
            RooAbsReal* Intpoly_la  = poly.createIntegral(x, x, "peak");

            double Intgaus1E_la = gaus1F_la*Intgaus1_la->getVal();
            double Intgaus2E_la = gaus2F_la*Intgaus2_la->getVal();
            double IntpolyE_la  = polyF_la*Intpoly_la->getVal();
            double totsig_la    = Intgaus1E_la + Intgaus2E_la + IntpolyE_la;
            double yield_la       = Intgaus1E_la + Intgaus2E_la;

            double Fsig_la = yield_la/totsig_la;

            mass_la.push_back(mean_la);
            std_la.push_back(rms_la);
            fsig_la.push_back(Fsig_la);
            covQual_la.push_back(r_la->covQual());

            cout << "Yield (la):" << yield_la << endl;
            cout << "Fsig (la): " << Fsig_la << endl;
            cout << "std (la): " << rms_la << endl;
            cout << "covQual (la)" << covQuality_la << endl;

            RooPlot* xframe = x.frame(160);
            xframe->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");
            xframe->GetYaxis()->SetTitle("Candidates / 0.0005GeV");
            xframe->GetXaxis()->CenterTitle(1);
            xframe->GetYaxis()->CenterTitle(1);
            xframe->GetXaxis()->SetNdivisions(507);
            xframe->GetXaxis()->SetTitleSize(0.05);
            xframe->GetYaxis()->SetTitleSize(0.05);
            xframe->GetYaxis()->SetTitleOffset(1);
            xframe->GetXaxis()->SetTitleOffset(1);
            xframe->GetYaxis()->SetLabelSize(0.05);
            xframe->GetXaxis()->SetLabelSize(0.05);
            data.plotOn(xframe,Name("data"));
            sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(1),LineColor(kRed));
            sum.plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(1),LineColor(kRed));
            cc1->cd(2);
            gPad->SetTickx();
            gPad->SetTicky();
            Xframe_La.push_back(xframe);
            xframe->Draw();

            os << "Pt Bin: " << (pla[i]-1)/10 << " - " << pla[i+1]/10;
            tex->DrawLatex(0.15,0.8,os.str().c_str());
            os.str(std::string());
            os << "Mean: " << mean_la;
            tex->DrawLatex(0.15,0.75,os.str().c_str());
            os.str(std::string());
            os << "#sigma :" << rms_la;
            tex->DrawLatex(0.15,0.70,os.str().c_str());
            os.str(std::string());
            os << "CovQual: " << covQuality_la;
            tex->DrawLatex(0.15,0.65,os.str().c_str());
            os.str(std::string());
            osYield << "Yield: " << std::setprecision(2) << yield_la;
            tex->DrawLatex(0.15,0.60,osYield.str().c_str());
            osYield.str(std::string());
            //os << "#chi^{2}/ndf: " << chi2_la;
            //tex->DrawLatex(0.15,0.60,os.str().c_str());
            //os.str(std::string());

            Composite_La->cd(index-3);
            gPad->SetTickx();
            gPad->SetTicky();
            xframe->GetXaxis()->SetNdivisions(507);
            xframe->GetXaxis()->SetTitleSize(0.05);
            xframe->GetYaxis()->SetTitleSize(0.05);
            xframe->GetYaxis()->SetTitleOffset(1.1);
            xframe->GetXaxis()->SetTitleOffset(1);
            xframe->GetYaxis()->SetLabelSize(0.05);
            xframe->GetXaxis()->SetLabelSize(0.05);
            xframe->Draw();

            os << "Pt Bin: " << (pks[i]-1)/10 << " - " << pks[i+1]/10;
            tex->DrawLatex(0.15,0.8,os.str().c_str());
            os.str(std::string());
            os << "Mean: " << mean_la;
            tex->DrawLatex(0.15,0.75,os.str().c_str());
            os.str(std::string());
            os << "#sigma :" << rms_la;
            tex->DrawLatex(0.15,0.70,os.str().c_str());
            os.str(std::string());
            os << "CovQual: " << covQuality_la;
            tex->DrawLatex(0.15,0.65,os.str().c_str());
            os.str(std::string());
            osYield << "Yield: " << std::setprecision(2) << yield_la;
            tex->DrawLatex(0.15,0.60,osYield.str().c_str());
            osYield.str(std::string());
            //os << "#chi^{2}/ndf: " << chi2_la;
            //tex->DrawLatex(0.15,0.60,os.str().c_str());
            //os.str(std::string());
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
