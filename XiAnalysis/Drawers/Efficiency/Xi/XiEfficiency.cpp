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
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "TString.h"
#include "TGaxis.h"

#include <vector>

void XiEfficiency()
{
    //Initializers
    std::string FileName = "";
    FileName = "XiEfficiency.txt";
    using namespace RooFit;
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetMarkerSize(0.5);
    TGaxis::SetMaxDigits(3);
    RooMsgService::instance().setStreamStatus(0,kFALSE);
    RooMsgService::instance().setStreamStatus(1,kFALSE);
    std::ostringstream os;
    std::ostringstream osYield;
    std::ofstream myfile;

    TH1D* massxi;
    //TH1D* massxi_gen;
    TH2D* MassXi;
    TH2D* MassXi_gen;

    TCanvas* Composite_Xi[8];
    //TCanvas* Composite_Xi_Gen[7];
    TCanvas* cc1 = new TCanvas("cc1","cc1",1000,450);
    cc1->Divide(2,1);

    std::vector<RooPlot*> Xframe_Xi;
    // rapidity index, values by pt bin
    std::map<int,std::vector<double> > yield_xi_;
    std::map<int,std::vector<double> > yield_xi_gen;


    std::map<int,std::vector<double> > efficiency_xi;

    std::vector<double> pxi = {11,13, 14,15, 16,17, 18,19, 20,21, 22,23, 24,25, 26,27, 28,29, 30,32, 33,34, 35,37, 38,40, 41,44, 45,60, 61,100};
    std::vector<double> rap_bin = {2,11,12,21};
    int numRapBins = rap_bin.size()/2;
    int pages=4;

    for(int i=0; i<pages; i++)
    {
        Composite_Xi[i] = new TCanvas(Form("Composite_Xi_%d",i),Form("Composite_Xi_RapBin_%d",i),1000,1200);
        Composite_Xi[i]->Divide(3,5);

        //Composite_Xi_Gen[i] = new TCanvas(Form("Composite_Xi_Gen_%d",i),Form("Composite_Xi_Gen_RapBin_%d"),1000,1200);
        //Composite_Xi_Gen[i]->Divide(4,4);
    }

    //File Creation
    myfile.open(FileName.c_str());
    TFile* f1 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/MCMassPtTotal_08_23_2017.root");

    TH3D* XiMassPtRap = (TH3D*)f1->Get("MassPtRapidityMC/XiMassPtRap");

    MassXi = (TH2D*)XiMassPtRap->Project3D("yx");

    TH3D* XiMassPtRap_Gen = (TH3D*)f1->Get("MassPtRapidityMC/XiMassPtRap_Gen");

    MassXi_gen = (TH2D*)XiMassPtRap_Gen->Project3D("yx");

    //Fit
    int pxicounter  = 0; //for correct bin counting
    int hcounter_xi = -1;
    int index = 1;
    double xstart_xi = 0.65;
    double ystart_xi = 0.85;
    double xpos = xstart_xi;
    double ypos = ystart_xi;
    double increment = 0.07;
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(62);
    tex->SetTextSize(0.05);
    for(unsigned j=0; j<rap_bin.size(); j+=2)
    //for(unsigned j=0; j<2; j+=2)
    {
        for(unsigned i=0; i<pxi.size(); i+=2)
        //for(unsigned i=0; i<60; i+=2)
        {
            cout << i << endl;
            if((double)i/30 == 1 || i==0) 
            {
                hcounter_xi++;
                index = 1;
            }

            massxi = (TH1D*)XiMassPtRap->ProjectionX(Form("massxi_%d",(j*pxi.size()+i)/2), pxi[i],pxi[i+1],rap_bin[j],rap_bin[j+1]);

            //tex->SetTextSize(tex->GetTextSize()*0.95);

            RooRealVar x("x","mass",1.26,1.4);
            RooPlot* xframe_xi = x.frame(150);
            xframe_xi->GetXaxis()->SetTitle("Invariant mass (GeV)");
            xframe_xi->GetYaxis()->SetTitle("Candidates / 0.001 GeV");
            xframe_xi->GetXaxis()->CenterTitle(1);
            xframe_xi->GetYaxis()->CenterTitle(1);
            xframe_xi->GetXaxis()->SetTickSize(0.02);
            xframe_xi->GetYaxis()->SetTickSize(0.02);
            xframe_xi->GetXaxis()->SetNdivisions(407);
            xframe_xi->GetYaxis()->SetNdivisions(410);
            xframe_xi->GetXaxis()->SetTitleSize(0.06);
            xframe_xi->GetYaxis()->SetTitleSize(0.06);
            xframe_xi->GetYaxis()->SetTitleOffset(1.05);
            xframe_xi->GetXaxis()->SetTitleOffset(0.5);
            //xframe_xi->GetXaxis()->SetLabelSize(xframe_xi->GetXaxis()->GetLabelSize()*2.0);
            xframe_xi->GetYaxis()->SetLabelSize(0.1);
            xframe_xi->GetXaxis()->SetLabelSize(0.1);
            //xframe_xi->GetYaxis()->SetLabelSize(xframe_xi->GetYaxis()->GetLabelSize()*2.0);
            RooDataHist data("data","dataset",x,massxi);
            data.plotOn(xframe_xi,Name("data"));
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
            //RooChi2Var chi2_xiVar("chi2_xi","chi2",sum,data);

            //double covQual = r_xi->covQual();
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
            double yield_xi          = Intgaus1E_xi + Intgaus2E_xi;


            double Fsig_xi = yield_xi/totsig_xi;

            yield_xi_[j/2].push_back(yield_xi);

            sum.plotOn(xframe_xi,Name("sum"),NormRange("cut"),LineWidth(1),LineColor(kBlue));
            sum.plotOn(xframe_xi,Components(background),NormRange("cut"),LineStyle(kDashed),LineWidth(1),LineColor(kBlue));
            Composite_Xi[hcounter_xi]->cd(index);
            gPad->SetBottomMargin(0.15); //gives more space for titles
            gPad->SetLeftMargin(0.15);
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            Xframe_Xi.push_back(xframe_xi);
            xframe_xi->Draw();
            Composite_Xi[hcounter_xi]->Update();

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
                os << "CMS pPb";
                tex->SetTextSize(0.06);
                tex->DrawLatex(0.17,ypos-increment,os.str().c_str());
                tex->SetTextSize(0.04);
                os.str(std::string());
                os << "185 #leq N_{trk}^{offline} < 250";
                tex->DrawLatex(0.17,ypos-2*increment,os.str().c_str());
                os.str(std::string());
                os << (-1.1 + (rap_bin[j]-1)/10) << " < y < " << (-1.1 + (rap_bin[j+1])/10);
                tex->DrawLatex(0.18,ypos-3*increment,os.str().c_str());
                os.str(std::string());
            }
            os  << (pxi[i]-1)/10 << " < P_{t} < "  << pxi[i+1]/10 << " GeV";
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            osYield << "Yield: " << std::setprecision(2) << yield_xi;
            tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
            osYield.str(std::string());
            os.str(std::string());
            index++;
        }
    }

    for(int j=0; j<pages; j++)
    {
        Composite_Xi[j]->Print(Form("RECOXiMassFitComposite_%d.pdf",j));
    }

    for(unsigned j=0; j<rap_bin.size(); j+=2)
    {
        for(unsigned i=0; i<pxi.size(); i+=2)
        {
            TH1D* massxi_gen = (TH1D*)XiMassPtRap_Gen->ProjectionX(Form("massxi_gen_%d",(int)(j*pxi.size()+i)/2), pxi[i],pxi[i+1],rap_bin[j],rap_bin[j+1]);
            yield_xi_gen[j/2].push_back(massxi_gen->GetBinContent(massxi_gen->GetMaximumBin()));
        }
    }
    //Composite_Xi_Gen->Print("Composite_Xi_Gen.pdf");

    //Calculate efficiency
    for(unsigned j=0; j<numRapBins; j++)
    {
        for(unsigned i=0; i<yield_xi_[j].size(); i++)
        {
            efficiency_xi[j].push_back(yield_xi_[j][i]/yield_xi_gen[j][i]);
        }
    }

    myfile << "Efficiency Xi\n";
    for(int j=0; j<numRapBins; j++)
    {
        myfile << "Rap_bin " << j << "\n";
        for(unsigned i=0; i<efficiency_xi[j].size(); i++)
            myfile << efficiency_xi[j][i] << "\n";
    }

    //const int rbinks = 26; // This is the number of BINS
    //const int rbinla = 19;
    //const int rbinrap = 4;

    //std::vector<double> vRebin_ks = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.3,1.5,1.7,1.9,2.2,2.5,2.8,3.1,3.4,3.7,4.0,4.4,5.0,6.0,7.0,10.0};
    //double Rebin_ks[] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.3,1.5,1.7,1.9,2.2,2.5,2.8,3.1,3.4,3.7,4.0,4.4,5.0,6.0,7.0,10.0};
    //std::vector<double> vRebin_la = {0.0,0.5,0.8,1.1,1.4,1.7,2.0,2.3,2.6,2.9,3.2,3.6,4.0,4.4,4.8,5.2,5.6,6.4,8.0,10.0};
    //double Rebin_la[] = {0.0,0.5,0.8,1.1,1.4,1.7,2.0,2.3,2.6,2.9,3.2,3.6,4.0,4.4,4.8,5.2,5.6,6.4,8.0,10.0};
    //std::vector<double> vRebin_rap = {-1.0,-0.5,-0.1,0.3,1.0};
    //double Rebin_rap[] = {-1.0,-0.5,-0.1,0.3,1.0};

    //std::vector<double> ptValues_ks;
    //std::vector<double> ptValues_la;
    //std::vector<double> RapValues;
    //for(unsigned i=0; i<vRebin_ks.size(); i++)
    //{
        //ptValues_ks.push_back(Rebin_ks[i]+0.01);
    //}

    //for(unsigned i=0; i<vRebin_la.size(); i++)
    //{
        //ptValues_la.push_back(Rebin_la[i]+0.01);
    //}

    //for(int i=0; i<vRebin_rap.size(); i++)
    //{
        //RapValues.push_back(Rebin_rap[i]);
    //}
    //TH2D* Effhisto_ks = new TH2D("EffHistoXi","EffHistoXi",rbinrap,Rebin_rap,rbinks,Rebin_ks);
    //TH2D* Effhisto_la = new TH2D("EffHistoLa","EffHistoLa",rbinrap,Rebin_rap,rbinla,Rebin_la);

    //for(int j=0; j<4; j++)
    //{
        //for(unsigned i=0; i<efficiency_xi[j].size(); i++)
            //Effhisto_ks->Fill(RapValues[j],ptValues_ks[i],efficiency_xi[j][i]);
    //}

    //for(int j=0; j<4; j++)
    //{
        //for(unsigned i=0; i<efficiency_la[j].size(); i++)
            //Effhisto_la->Fill(RapValues[j],ptValues_la[i],efficiency_la[j][i]);
    //}
    //TCanvas* Eff = new TCanvas("Eff","",1200,800);
    //Eff->Divide(2,1);
    //Eff->cd(1);
    //Effhisto_ks->Draw("Lego2");
    //Eff->cd(2);
    //Effhisto_la->Draw("Lego2");

    //Eff->Print("EffhistoXi.pdf");

    //TFile histos("EffhistoXi.root","RECREATE");
    //Effhisto_ks->Write();
    //Effhisto_la->Write();

    //std::vector<double> eta_trg = {-0.9,-0.8,-0.4,-0.2,-0.1,0.2,0.5,0.9};
    //std::vector<double> pt_trg = {0.5,1.2, 3.4, 4.0, 6.7, 8.6, 1.5, 2.0};
    //for(unsigned i=0; i<eta_trg.size(); i++)
    //{
        //cout << Effhisto_ks->GetBinContent(Effhisto_ks->FindBin(eta_trg[i],pt_trg[i])) << endl;
    //}


}
