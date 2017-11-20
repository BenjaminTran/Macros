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

    TH1D* massom;
    TH2D* MassXi;
    TH3D* MassPtRapXi;
    std::vector<RooPlot*> xframe;
    std::vector<double> mass_om;
    std::vector<double> std_om;
    std::vector<double> fsig_om;
    std::vector<double> covQual_om;
    bool doRap = true;
    bool doPbPb = true;

    //std::vector<double> pom = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,72, 73,85, 86,100, 101,200, 201,300};
    //std::vector<double> pom = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,100, 101,200};//, 201,300};
    //std::vector<double> pom = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,72, 73,100, 101,200};//, 201,300};
    //std::vector<double> pom = {11,15, 16,19, 20,23, 24,27, 28,33, 34,41, 42,50 ,51,60, 61,72, 73,100};//, 81,100, 101,200};//, 201,300};
    std::vector<double> pom = {16,18, 19,22, 23,28, 29,36, 37,46, 47,60 ,61,72, 73,100}; //pPb
    //std::vector<double> pom = {16,18, 19,22, 23,28, 29,36, 37,50, 51,80}; //pPb MB 0-35
    //std::vector<double> pom = {11,18, 19,23, 24,30, 31,100}; //pPb MB 0-20
    //std::vector<double> pom = {16,18, 19,22, 23,28, 29,36, 37,46, 47,60};; //PbPb

    TCanvas* cc1 = new TCanvas("cc1","cc1",1600,900);
    if(!doPbPb)cc1->Divide(3,3);
    else cc1->Divide(3,2);

    TCanvas* c_om_NF = new TCanvas("c_om_NF","c_om_NoFit",1600,900);
    c_om_NF->Divide(3,3);
    TCanvas* c_om = new TCanvas("c_om","c_om",800,600);

    //File Creation
    myfile.open("OmPeakParam_10_23_17.txt");
    TFile* file = nullptr;
    //TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/CascadeV2pPb/RootFiles/Flow/CasCutLoose/CasCutLooseJL40.root");
    //TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/Thesis/XiAnalysisCorrelationPtCut8TeVPD1_4_ForFinal.root");
    //TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/Composites/V0CasMassPtPD5JL12.root"); //only one PD
    //file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MBCorr/XiOmegaMB_0_N_20_Partial_11_8_17.root"); //MB 0-20 Full Stats
    //file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/MinBias/Comparison/XiOmegaCompPD1_0_20_11_10_17.root"); //MB one PD for comp to 0-35
    //file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/MinBias/Comparison/XiOmegaMassPtMBPD1_0_35_11_10_17.root"); //MB one PD for comp to 0-35
    //file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/MinBias/Comparison/OmegaHMPD1Comparison_11_14_17.root"); //HM one PD for comp to 0-35
    if(doRap && !doPbPb) file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/OmCorr/OmCorrelationRapidityTotal_09_24_17.root"); //pPb Full Stats
    else if(doRap && doPbPb) file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/V0CasCorrelationPbPbTotal_10_30_17.root"); //PbPb Full Stats
    //else  file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MassPt/Omega/OmMassPt.root"); //only one PD

    //MassXi = (TH2D*)file->Get("XiMassPt/MassPt");
    //MassXi = (TH2D*)file->Get("omCorrelation/MassPt");
    //MassXi = (TH2D*)file->Get("v0CasCorrelationRapidityPeriSub/MassPtOm");//pPb MB
    //MassPtRapXi = (TH3D*)file->Get("MassPtRapidityMB/OmMassPt");
    //MassXi = (TH2D*)MassPtRapXi->Project3D("yx");
    if(doRap && !doPbPb)MassXi = (TH2D*)file->Get("omCorrelationRapidity/MassPt");//pPb
    else if(doRap && doPbPb)MassXi = (TH2D*)file->Get("v0CasCorrelationRapidityPbPb/MassPtOm"); //PbPb

    //Fit
    int pomcounter = 0; //for correct bin counting
    int hbincounter =1; //histogram bin counting
    //for(unsigned i=0; i<6; i++)
    for(unsigned i=0; i<pom.size(); i++)
    {
        double xpos = 0.60;
        double ypos = 0.85;
        double increment = 0.07;
        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextFont(62);
        tex->SetTextSize(0.04);
        //tex->SetTextAlign(10);
        TCanvas* cc2 = new TCanvas("cc2","",600,450);
        int index = (i+2)/2;
        massom = (TH1D*)MassXi->ProjectionX("massom", pom[i],pom[i+1]);
        massom->GetXaxis()->SetRangeUser(1.60, 1.75);

        //TH1D* massom_clone = (TH1D*)massom->Clone("massom_clone");

        //c_om_NF->cd(index);
        //massom_clone->Draw();

        os  << (pom[i]-1)/10 << " < P_{t} < "  << pom[i+1]/10 << " GeV";
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());

        gStyle->SetOptTitle(kFALSE);


        RooRealVar x("x","mass",1.60,1.75);
        RooPlot* xframe_ = x.frame(150);
        xframe_->GetXaxis()->SetTitle("#Lambda K Invariant mass (GeV)");
        xframe_->GetYaxis()->SetTitle("Candidates / 1 MeV");
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
        RooDataHist data("data","dataset",x,massom);
        data.plotOn(xframe_,Name("data"));
        RooRealVar mean("mean","mean",1.67,1.65,1.69);
        RooRealVar sigma1("sigma1","sigma1",0.004,0.001,0.04);
        RooRealVar sigma2("sigma2","sigma2",0.005,0.001,0.04);
        RooRealVar sig1("sig1","signal1",4500,0,1000000);
        RooRealVar sig2("sig2","signal2",4500,0,1000000);
        RooRealVar qsig("qsig","qsig",22000,0,1000000);
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
        //RooGenericPdf background("background", "x - (1.115683 + 0.493677)^alpha", RooArgList(x,alpha));
        RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,background),RooArgList(sig1,sig2,qsig));

        if(!doRap)
        {
            //if(i==2 || i==4)
                //x.setRange("cut",1.65,1.694);
            //else if(i==pom.size()-2)
                //x.setRange("cut",1.645,1.71);
            //else
                //x.setRange("cut",1.645,1.7);
            x.setRange("cut",1.62,1.75);
        }
        else
        {
            x.setRange("cut",1.61,1.75);
        }

        RooFitResult* r_om = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));
        //RooChi2Var chi2_omVar("chi2_om","chi2",sum,data);

        double covQual = r_om->covQual();
        double mean_om = mean.getVal();

        double gaus1F_om = sig1.getVal();
        double gaus2F_om = sig2.getVal();
        double qsig_om   = qsig.getVal();

        //set ranges for individual gaussian yield determination
        x.setRange("g1", mean.getVal() - 2*sigma1.getVal(), mean.getVal() + 2*sigma1.getVal());
        x.setRange("g2", mean.getVal() - 2*sigma2.getVal(), mean.getVal() + 2*sigma2.getVal());

        RooAbsReal* Intgaus1_yield_om = gaus1.createIntegral(x,x,"g1");
        RooAbsReal* Intgaus2_yield_om = gaus2.createIntegral(x,x,"g2");

        double gaus1_yield_om = gaus1F_om*Intgaus1_yield_om->getVal();
        double gaus2_yield_om = gaus2F_om*Intgaus2_yield_om->getVal();
        double gausTot_yield_om = gaus1_yield_om + gaus2_yield_om;

        cout << "Yield1: " << gaus1_yield_om << endl;
        cout << "Yield2: " << gaus2_yield_om << endl;

        double rms_gaus1_sig_om = gaus1_yield_om/gausTot_yield_om;
        double rms_gaus2_sig_om = gaus2_yield_om/gausTot_yield_om;
        double rms_true_om = TMath::Sqrt(rms_gaus1_sig_om*sigma1.getVal()*sigma1.getVal() + rms_gaus2_sig_om*sigma2.getVal()*sigma2.getVal());

        x.setRange("peak", mean.getVal() - 2*rms_true_om, mean.getVal() + 2*rms_true_om);
        RooAbsReal* Intgaus1_om      = gaus1.createIntegral(x, x,  "peak");
        RooAbsReal* Intgaus2_om      = gaus2.createIntegral(x, x, "peak");
        RooAbsReal* Intbackground_om = background.createIntegral(x, x, "peak");

        double Intgaus1E_om      = gaus1F_om*Intgaus1_om->getVal();
        double Intgaus2E_om      = gaus2F_om*Intgaus2_om->getVal();
        double IntbackgroundE_om = qsig_om*Intbackground_om->getVal();
        double totsig_om         = Intgaus1E_om + Intgaus2E_om + IntbackgroundE_om;
        double Yield_om          = Intgaus1E_om + Intgaus2E_om;


        double Fsig_om = Yield_om/totsig_om;

        mass_om.push_back(mean_om);
        std_om.push_back(rms_true_om);
        fsig_om.push_back(Fsig_om);
        covQual_om.push_back(covQual);

        cout << "Yield (om): " << Yield_om << endl;
        cout << "Fsig (om): " << Fsig_om << endl;
        cout << "std (om): "  << rms_true_om  << endl;
        cout << "mass (om): " << mean_om << endl;

        cout << "covQual (om)" << covQual << endl;

        sum.plotOn(xframe_,Name("sum"),NormRange("cut"),LineWidth(2),LineColor(kBlue));
        sum.plotOn(xframe_,Components(background),NormRange("cut"),LineStyle(kDashed),LineWidth(2),LineColor(kBlue));
        cc1->cd(index);
        gPad->SetBottomMargin(0.15); //gives more space for titles
        gPad->SetLeftMargin(0.15);
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        xframe.push_back(xframe_);
        xframe_->Draw();
        cc1->Update();
        /*
        if(i==6)
        {
            RooPlot* xframe_som = x.frame(150);
            xframe_som->GetXaxis()->SetTitle("m_{#Lambda K} (GeV/c^{2})");
            xframe_som->GetYaxis()->SetTitle("Candidates / 1 MeV");
            xframe_som->GetXaxis()->CenterTitle(1);
            xframe_som->GetYaxis()->CenterTitle(1);
            xframe_som->GetXaxis()->SetTickSize(0.02);
            xframe_som->GetYaxis()->SetTickSize(0.02);
            xframe_som->GetXaxis()->SetNdivisions(407);
            xframe_som->GetYaxis()->SetNdivisions(410);
            xframe_som->GetXaxis()->SetTitleSize(0.05);
            xframe_som->GetYaxis()->SetTitleSize(0.05);
            xframe_som->GetYaxis()->SetTitleOffset(0.9);
            xframe_som->GetXaxis()->SetTitleOffset(0.9);
            //xframe_som->GetXaxis()->SetLabelSize(xframe_som->GetXaxis()->GetLabelSize()*2.0);
            xframe_som->GetYaxis()->SetLabelSize(0.04);
            xframe_som->GetXaxis()->SetLabelSize(0.04);
            data.plotOn(xframe_som,Name("data"));
            sum.plotOn(xframe_som,Name("sum"),NormRange("cut"),LineWidth(2),LineColor(kBlue));
            sum.plotOn(xframe_som,Name("poly"),Components(background),NormRange("cut"),LineStyle(kDashed),LineWidth(2),LineColor(kBlue));
            c_om->cd();
            gPad->SetTickx();
            gPad->SetTicky();
            xframe_som->Draw();
            c_om->Update();
            TLine* t1_om = new TLine(mean.getVal() - 2*rms_true_om, 0, mean.getVal() - 2*rms_true_om, gPad->GetUymax());
            TLine* t2_om = new TLine(mean.getVal() + 2*rms_true_om, 0, mean.getVal() + 2*rms_true_om, gPad->GetUymax());
            t1_om->SetLineStyle(2);
            t1_om->SetLineColor(kGreen+1);
            t2_om->SetLineStyle(2);
            t2_om->SetLineColor(kGreen+1);
            t1_om->Draw();
            t2_om->Draw();

            double xstart_om = 0.13;
            double ystart_om = 0.82;
            double xpos = xstart_om;
            double ypos = ystart_om;
            double increment = 0.07;
            tex->SetTextSize(0.05);
            os << "CMS pPb #sqrt{S_{NN}} = 8.16 TeV";
            tex->DrawLatex(0.513,0.918,os.str().c_str());
            os.str(std::string());
            tex->SetTextSize(0.045);
            os << "185 #leq N_{trk}^{offline} < 250 ";
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            os << (pom[i]-1)/10 << " < p_{T} < " << pom[i+1]/10  << " GeV/c";
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());
            os << "|y| < 1";
            tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
            os.str(std::string());

            TLegend* leg = new TLegend(0.58,0.61,0.75,0.8);
            leg->SetFillColor(10);
            leg->SetFillStyle(0);
            leg->SetBorderSize(0);
            leg->SetTextFont(62);
            leg->SetTextSize(0.045);
            leg->AddEntry(xframe_som->findObject("data"),"Data","P");
            leg->AddEntry(xframe_som->findObject("sum"),"Fit","l");
            leg->AddEntry("poly","Background","l");
            leg->AddEntry(t1_om,"#pm 2#sigma","l");
            leg->Draw();

            c_om->Print("OmPlotForZhenyu.pdf");
        }
        cc1->cd();
        */
        double chi2_om = xframe_->chiSquare("sum","data",4);

        TLine* t1 = new TLine(mean.getVal() - 2*rms_true_om, 0, mean.getVal() - 2*rms_true_om, gPad->GetUymax());
        TLine* t2 = new TLine(mean.getVal() + 2*rms_true_om, 0, mean.getVal() + 2*rms_true_om, gPad->GetUymax());
        t1->SetLineStyle(2);
        t1->SetLineColor(kGreen);
        t2->SetLineStyle(2);
        t2->SetLineColor(kGreen);
        t1->Draw("same");
        t2->Draw("same");
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
        os  << (pom[i]-1)/10 << " < P_{t} < "  << pom[i+1]/10 << " GeV";
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "Mean: " << std::setprecision(5) << mean_om << " GeV" << std::setprecision(6);
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "#sigma :" << std::setprecision(2) << rms_true_om << " GeV" << std::setprecision(6);
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        os << "CovQual: " << covQual;
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        os.str(std::string());
        //os << "#chi^{2}/ndf: " << std::setprecision(3) << chi2_om << std::setprecision(6);
        //tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
        //os.str(std::string());
        osYield << "Yield: " << std::setprecision(2) << Yield_om;
        tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
        osYield.str(std::string());
        //os.str(std::string());
        os << "Fsig: " << Fsig_om;
        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
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
        os << "P_{t} Bin: " << (pom[i]-1)/10 << " - " << pom[i+1]/10;
        tex->DrawLatex(xpos,0.8,os.str().c_str());
        os.str(std::string());
        os << "Mean: " << mean_om;
        tex->DrawLatex(xpos,0.75,os.str().c_str());
        os.str(std::string());
        os << "#sigma :" << rms_true_om;
        tex->DrawLatex(xpos,0.70,os.str().c_str());
        os.str(std::string());
        os << "CovQual: " << covQual;
        tex->DrawLatex(xpos,0.65,os.str().c_str());
        os.str(std::string());
        osYield << "Yield: " << std::setprecision(2) << Yield_om;
        tex->DrawLatex(xpos,0.60,osYield.str().c_str());
        osYield.str(std::string());
        //os << "#chi^{2}/ndf: " << chi2_om;
        //tex->DrawLatex(xpos,0.60,os.str().c_str());
        //os.str(std::string());


        hbincounter++;
        if(i==0) cc2->Print("OmMassFitInd.pdf(","pdf");
        else if(i < pom.size() - 2) cc2->Print("OmMassFitInd.pdf","pdf");
        else cc2->Print("OmMassFitInd.pdf)","pdf");
        i++; //to access correct bins
    }
    if(!doRap)
    {
        cc1->Print("OmMassFitComposite.pdf");
        cc1->Print("OmMassFitComposite.png");
    }
    else
    {
        cc1->Print("OmMassFitCompositeD0Ana.pdf");
        cc1->Print("OmMassFitCompositeD0Ana.png");
    }
    //c_om_NF->Print("OmMassFitHM185_250.pdf");

    //Output
    pomcounter = 0;
    for(unsigned i=0; i<mass_om.size(); i++)
    {
        cout <<  "====================" << endl;
        cout << "Pt Bin: " << (pom[pomcounter]-1)/10 << " - " << pom[pomcounter+1]/10 << endl;
        cout <<  "====================" << endl;
        cout << "Mass_om: "   << mass_om[i]    << endl;
        cout << "Fsig_om: "   << fsig_om[i]    << endl;
        cout << "std_om: "    << std_om[i]     << endl;
        cout << "covQual_om " << covQual_om[i] << endl;

        myfile <<  "====================" << "\n";
        myfile << "Pt Bin: " << (pom[pomcounter]-1)/10 << " - " << pom[pomcounter+1]/10 << "\n";
        myfile <<  "====================" << "\n";
        myfile << "Mass_om: "   << mass_om[i]    << "\n";
        myfile << "Fsig_om: "   << fsig_om[i]    << "\n";
        myfile << "std_om: "    << std_om[i]     << "\n";
        myfile << "covQual_om " << covQual_om[i] << "\n";

        pomcounter+=2;
    }
}
