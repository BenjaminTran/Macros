//#include "makeMultiPanelCanvas.C"
#include"MITStyleOrgiinal.C"
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
#include "TBox.h"

#include <vector>

void plot_v2_ncq_sub()
{
    gStyle->SetOptTitle(0);
    gStyle->SetErrorX(0.);
    
    TFile* file[5];
    
    TFile* file[0] = TFile::Open("V0v2perisub_Default_FullXiHM.root");
    TFile* file[1] = TFile::Open("2PC_v2vspt_fromfit_d0_HM185_250_combine_deta2_sub.root");
    //TFile* file[2] = TFile::Open("SP_v2vspt_fromfit_d0_HM185_250_combine_Pbside.root");
    
    TFile ofile("NCQratio_8TeVpPb_sub.root","RECREATE");

    TGraphErrors* v2ksplot = file[0]->Get("kshortv2truesub");
    TGraphErrors* v2laplot = file[0]->Get("lambdav2truesub");
    TGraphErrors* v2xiplot = file[0]->Get("xiv2truesub");
    TGraphErrors* v2ksplot_KET = file[0]->Get("kshortv2truesub_KET");
    TGraphErrors* v2laplot_KET = file[0]->Get("lambdav2truesub_KET");
    TGraphErrors* v2xiplot_KET = file[0]->Get("xiv2truesub_KET");
    //TGraphErrors* v2omplot = file[0]->Get("v2omega");
    TGraphErrors* v2d2PCplot = file[1]->Get("v2vspt");
    //TGraphErrors* v2dSPplot = file[2]->Get("v2vspt");

    //make NCQ TGraph
    double ptks[100],ptla[100],ptxi[100],v2ks[100],v2la[100],v2xi[100],v2kse[100],v2lae[100],v2xie[100],KET_ks[100],KET_la[100],KET_xi[100];
    double v2d2PC[100], v2dSP[100], v2d2PCe[100], v2dSPe[100], KET_d2PC[100], KET_dSP[100];
    
    int Nks = v2ksplot->GetN();
    int Nla = v2laplot->GetN();
    int Nxi = v2xiplot->GetN();
    //int Nd = v2d2PCplot->GetN();
    
    double ksMass = 0.497;
    double laMass = 1.1156;
    double xiMass = 1.322;
    
    for(int i=0;i<Nks;i++)
    {
        v2ksplot_KET->GetPoint(i,ptks[i],v2ks[i]);
        v2ks[i] = v2ks[i]/2.0;
        v2kse[i] = v2ksplot->GetErrorY(i)/2.0;
        KET_ks[i] = ptks[i]/2.0;
    }
    for(int i=0;i<Nla;i++)
    {
        v2laplot_KET->GetPoint(i,ptla[i],v2la[i]);
        v2la[i] = v2la[i]/3.0;
        v2lae[i] = v2laplot->GetErrorY(i)/3.0;
        KET_la[i] = ptla[i]/3.0;
    }
    for(int i=0;i<Nxi;i++)
    {
        v2xiplot_KET->GetPoint(i,ptxi[i],v2xi[i]);
        v2xi[i] = v2xi[i]/3.0;
        v2xie[i] = v2xiplot->GetErrorY(i)/3.0;
        KET_xi[i] = ptxi[i]/3.0;
    }

    TGraphErrors* D0v22PCncq = file[1].Get("v2vsKET_ncq");
    //TGraphErrors* D0v2SPncq = file[2].Get("v2vsKET_ncq");
    TGraphErrors* ksNCQ = new TGraphErrors(Nks,KET_ks,v2ks,0,v2kse);
    TGraphErrors* laNCQ = new TGraphErrors(Nla,KET_la,v2la,0,v2lae);
    TGraphErrors* xiNCQ = new TGraphErrors(Nxi,KET_xi,v2xi,0,v2xie);
    
    //TGraphErrors* ksNCQ = file[0]->Get("v2kshort_ket_nq");
    //TGraphErrors* laNCQ = file[0]->Get("v2lambda_ket_nq");
    //TGraphErrors* xiNCQ = file[0]->Get("v2xi_ket_nq");
    //TGraphErrors* omNCQ = file[0]->Get("v2omega_ket_nq");
    
    //make NCQ fit&ratio
    fitfunc_v2QNS = new TF1("fitfunc_v2QNS","([0]/(1+exp(-(x-[1])/[2]))-[3])*pol1(4)",0,3.5);
    fitfunc_v2QNS->SetParameters(0.5,-0.1,0.2,0.2,0,0);
    ksNCQ->Fit(fitfunc_v2QNS,"NOB",0,3.5);
    fitfunc_v2QNS->SetLineStyle(2);
    fitfunc_v2QNS->SetLineWidth(1);

    TGraphErrors* ksNCQratio = new TGraphErrors(ksNCQ->GetN());
    for(int nn=0;nn<ksNCQ->GetN();nn++)
    {
        double xx, yy, xx_err, yy_err;
        ksNCQ->GetPoint(nn,xx,yy);
        xx_err = ksNCQ->GetErrorX(nn);
        yy_err = ksNCQ->GetErrorY(nn);
        ksNCQratio->SetPoint(nn,xx,yy/fitfunc_v2QNS->Eval(xx));
        ksNCQratio->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_v2QNS->Eval(xx));
    }
    TGraphErrors* laNCQratio = new TGraphErrors(laNCQ->GetN());
    for(int nn=0;nn<laNCQ->GetN();nn++)
    {
        double xx, yy, xx_err, yy_err;
        laNCQ->GetPoint(nn,xx,yy);
        xx_err = laNCQ->GetErrorX(nn);
        yy_err = laNCQ->GetErrorY(nn);
        laNCQratio->SetPoint(nn,xx,yy/fitfunc_v2QNS->Eval(xx));
        laNCQratio->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_v2QNS->Eval(xx));
    }
    TGraphErrors* xiNCQratio = new TGraphErrors(xiNCQ->GetN());
    for(int nn=0;nn<xiNCQ->GetN();nn++)
    {
        double xx, yy, xx_err, yy_err;
        xiNCQ->GetPoint(nn,xx,yy);
        xx_err = xiNCQ->GetErrorX(nn);
        yy_err = xiNCQ->GetErrorY(nn);
        xiNCQratio->SetPoint(nn,xx,yy/fitfunc_v2QNS->Eval(xx));
        xiNCQratio->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_v2QNS->Eval(xx));
    }
    /*TGraphErrors* omNCQratio = new TGraphErrors(omNCQ->GetN());
    for(int nn=0;nn<omNCQ->GetN();nn++)
    {
        double xx, yy, xx_err, yy_err;
        omNCQ->GetPoint(nn,xx,yy);
        xx_err = omNCQ->GetErrorX(nn);
        yy_err = omNCQ->GetErrorY(nn);
        omNCQratio->SetPoint(nn,xx,yy/fitfunc_v2QNS->Eval(xx));
        omNCQratio->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_v2QNS->Eval(xx));
    }*/
    TGraphErrors* d2PCNCQratio = new TGraphErrors(D0v22PCncq->GetN());
    for(int nn=0;nn<D0v22PCncq->GetN();nn++)
    {
        double xx, yy, xx_err, yy_err;
        D0v22PCncq->GetPoint(nn,xx,yy);
        xx_err = D0v22PCncq->GetErrorX(nn);
        yy_err = D0v22PCncq->GetErrorY(nn);
        d2PCNCQratio->SetPoint(nn,xx,yy/fitfunc_v2QNS->Eval(xx));
        d2PCNCQratio->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_v2QNS->Eval(xx));
    }
    /*TGraphErrors* dSPNCQratio = new TGraphErrors(D0v2SPncq->GetN());
    for(int nn=0;nn<D0v2SPncq->GetN();nn++)
    {
        double xx, yy, xx_err, yy_err;
        D0v2SPncq->GetPoint(nn,xx,yy);
        xx_err = D0v2SPncq->GetErrorX(nn);
        yy_err = D0v2SPncq->GetErrorY(nn);
        dSPNCQratio->SetPoint(nn,xx,yy/fitfunc_v2QNS->Eval(xx));
        dSPNCQratio->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_v2QNS->Eval(xx));
    }*/
   
    //set line width
    v2ksplot->SetLineWidth(1);
    v2laplot->SetLineWidth(1);
    v2xiplot->SetLineWidth(1);
    //v2omplot->SetLineWidth(1);
    v2d2PCplot->SetLineWidth(1);
    //v2dSPplot->SetLineWidth(1);
    D0v22PCncq->SetLineWidth(1);
    //D0v2SPncq->SetLineWidth(1);
    ksNCQ->SetLineWidth(1);
    laNCQ->SetLineWidth(1);
    xiNCQ->SetLineWidth(1);
    //omNCQ->SetLineWidth(1);
    ksNCQratio->SetLineWidth(1);
    laNCQratio->SetLineWidth(1);
    xiNCQratio->SetLineWidth(1);
    d2PCNCQratio->SetLineWidth(1);
    //dSPNCQratio->SetLineWidth(1);
    
    
    //set color
    v2ksplot->SetMarkerColor(4);
    v2ksplot->SetLineColor(4);
    v2ksplot->SetMarkerStyle(25);
    ksNCQ->SetMarkerColor(4);
    ksNCQ->SetLineColor(4);
    ksNCQ->SetMarkerStyle(25);
    ksNCQratio->SetMarkerColor(4);
    ksNCQratio->SetLineColor(4);
    ksNCQratio->SetMarkerStyle(25);
    
    v2laplot->SetMarkerColor(1);
    v2laplot->SetLineColor(1);
    v2laplot->SetMarkerStyle(24);
    v2laplot->SetMarkerSize(1.4);
    laNCQ->SetMarkerColor(1);
    laNCQ->SetLineColor(1);
    laNCQ->SetMarkerStyle(24);
    laNCQ->SetMarkerSize(1.4);
    laNCQratio->SetMarkerColor(1);
    laNCQratio->SetLineColor(1);
    laNCQratio->SetMarkerStyle(24);
    laNCQratio->SetMarkerSize(1.4);

    v2xiplot->SetMarkerColor(kGreen+4);
    v2xiplot->SetLineColor(kGreen+4);
    v2xiplot->SetMarkerStyle(28);
    v2xiplot->SetMarkerSize(1.8);
    xiNCQ->SetMarkerColor(kGreen+4);
    xiNCQ->SetLineColor(kGreen+4);
    xiNCQ->SetMarkerStyle(28);
    xiNCQ->SetMarkerSize(1.8);
    xiNCQratio->SetMarkerColor(kGreen+4);
    xiNCQratio->SetLineColor(kGreen+4);
    xiNCQratio->SetMarkerStyle(28);
    xiNCQratio->SetMarkerSize(1.8);
    
    v2d2PCplot->SetMarkerColor(2);
    v2d2PCplot->SetLineColor(2);
    v2d2PCplot->SetMarkerStyle(20);
    v2d2PCplot->SetMarkerSize(1.5);
    D0v22PCncq->SetMarkerColor(2);
    D0v22PCncq->SetLineColor(2);
    D0v22PCncq->SetMarkerStyle(20);
    D0v22PCncq->SetMarkerSize(1.5);
    d2PCNCQratio->SetMarkerColor(2);
    d2PCNCQratio->SetLineColor(2);
    d2PCNCQratio->SetMarkerStyle(20);
    d2PCNCQratio->SetMarkerSize(1.5);
    
    /*v2dSPplot->SetMarkerColor(2);
    v2dSPplot->SetLineColor(2);
    v2dSPplot->SetMarkerStyle(21);
    D0v2SPncq->SetMarkerColor(2);
    D0v2SPncq->SetLineColor(2);
    D0v2SPncq->SetMarkerStyle(21);
    dSPNCQratio->SetMarkerColor(2);
    dSPNCQratio->SetLineColor(2);
    dSPNCQratio->SetMarkerStyle(21);*/
    
    /*v2omplot->SetMarkerColor(kGreen+2);
    v2omplot->SetLineColor(kGreen+2);
    v2omplot->SetMarkerStyle(27);
    omNCQ->SetMarkerColor(kGreen+2);
    omNCQ->SetLineColor(kGreen+2);
    omNCQ->SetMarkerStyle(27);
    omNCQratio->SetMarkerColor(kGreen+2);
    omNCQratio->SetLineColor(kGreen+2);
    omNCQratio->SetMarkerStyle(27);
    omNCQratio->SetMarkerSize(1.4);
    omNCQ->SetMarkerSize(1.4);
    v2omplot->SetMarkerSize(1.4);*/
    
    //drop points
    laNCQ->RemovePoint(10);
    //laNCQ->RemovePoint(12);
    laNCQratio->RemovePoint(10);
    //laNCQratio->RemovePoint(12);
    ksNCQ->RemovePoint(11);
    ksNCQratio->RemovePoint(11);
    
    //xiNCQ->RemovePoint(8);
    //xiNCQ->RemovePoint(9);
    //xiNCQratio->RemovePoint(8);
    //xiNCQratio->RemovePoint(9);
    
    //omNCQ->RemovePoint(0);
    //omNCQratio->RemovePoint(0);
    //v2omplot->RemovePoint(0);
    
    //make hist
    TH1D* hist = new TH1D("hist","",80,0,9.9);
    hist->SetLineWidth(0);

    fixedFontHist1D(hist,2.4,2.2);
    hist->SetLineStyle(9);
    hist->GetXaxis()->SetTitleSize(hist->GetXaxis()->GetTitleSize()*1.);
    hist->GetYaxis()->SetTitleSize(hist->GetYaxis()->GetTitleSize()*1.6);
    hist->GetYaxis()->SetTitleOffset(1.8);
    hist->GetXaxis()->SetLabelSize(hist->GetXaxis()->GetLabelSize()*1.1);
    hist->GetYaxis()->SetLabelSize(hist->GetYaxis()->GetLabelSize()*1.2);
    hist->SetXTitle("p_{T} (GeV/c)");
    hist->SetYTitle("v_{2}^{sub}");
    hist->GetXaxis()->CenterTitle(1);
    hist->GetYaxis()->CenterTitle(1);
    hist->SetMinimum(0.001);
    hist->SetMaximum(0.24);
    hist->GetXaxis()->SetRangeUser(0,8.5);
    
    TH1D* hist1 = new TH1D("hist1","",80,-0.1,6);
    hist1->SetLineWidth(0);

    fixedFontHist1D(hist1,2.4,2.2);
    hist1->SetLineStyle(9);
    hist1->GetXaxis()->CenterTitle(1);
    hist1->GetYaxis()->CenterTitle(1);
    hist1->GetXaxis()->SetTitleSize(hist1->GetXaxis()->GetTitleSize()*1.);
    hist1->GetYaxis()->SetTitleSize(hist1->GetYaxis()->GetTitleSize()*1.5);
    hist1->GetYaxis()->SetTitleOffset(2.0);
    hist1->GetXaxis()->SetLabelSize(hist1->GetXaxis()->GetLabelSize()*1.2);
    hist1->GetYaxis()->SetLabelSize(hist1->GetYaxis()->GetLabelSize()*1.2);
    hist1->SetXTitle("KE_{T}/n_{q} (GeV)");
    hist1->SetYTitle("v_{2}^{sub}/n_{q}");
    hist1->SetMinimum(0.001);
    hist1->SetMaximum(0.084);
    hist1->GetXaxis()->SetRangeUser(0,3.2);

    TH1D* hist2 = new TH1D("hist2","",80,-0.1,6);
    hist2->SetLineWidth(0);

    fixedFontHist1D(hist2,2.4,2.2);
    hist2->GetXaxis()->CenterTitle(1);
    hist2->GetYaxis()->CenterTitle(1);
    hist2->GetXaxis()->SetTitleSize(hist2->GetXaxis()->GetTitleSize()*1.);
    hist2->GetYaxis()->SetTitleSize(hist2->GetYaxis()->GetTitleSize()*1.2);
    hist2->GetXaxis()->SetTitleOffset(hist2->GetXaxis()->GetTitleOffset()*2.5);
    hist2->GetYaxis()->SetTitleOffset(2.4);
    hist2->GetXaxis()->SetLabelSize(hist2->GetXaxis()->GetLabelSize()*1.2);
    hist2->GetYaxis()->SetLabelSize(hist2->GetYaxis()->GetLabelSize()*1.2);
    hist2->SetXTitle("KE_{T}/n_{q} (GeV)");
    hist2->SetYTitle("Data/Fit");
    hist2->SetMinimum(0.1);
    hist2->SetMaximum(1.9);
    hist2->GetXaxis()->SetRangeUser(0,3.2);
    
    TCanvas *c2 = new TCanvas("c2","c2",1,1,960,750);
    c2->Range(0,0,1,1);
    TPad* pad[3];
    pad[0] = new TPad("pad0", "pad0",0.0,0.54125,0.5,1);
    pad[1] = new TPad("pad2", "pad2",0.0,0.17625,0.5,0.54125);
    pad[2] = new TPad("pad4", "pad4",0.0,0.0,0.5,0.17625);
    
    for(int i=0; i<3; i++)
    {
        pad[i]->SetLeftMargin(0);
        pad[i]->SetRightMargin(0);
        pad[i]->SetTopMargin(0);
        pad[i]->SetBottomMargin(0);
        pad[i]->Draw();
    }

    pad[0]->SetLeftMargin(0.16);
    pad[1]->SetLeftMargin(0.16);
    pad[2]->SetLeftMargin(0.16);
    pad[0]->SetRightMargin(0.05);
    pad[1]->SetRightMargin(0.05);
    pad[2]->SetRightMargin(0.05);
    
    pad[0]->SetTopMargin(0.1);
    pad[1]->SetTopMargin(0.1);
    pad[0]->SetBottomMargin(0.15);
    pad[2]->SetBottomMargin(0.42);
    
    TPad* pad_PbPb[3];
    pad_PbPb[0] = new TPad("pad0", "pad0",0.5,0.54125,1,1);
    pad_PbPb[1] = new TPad("pad2", "pad2",0.5,0.17625,1,0.54125);
    pad_PbPb[2] = new TPad("pad4", "pad4",0.5,0.0,1,0.17625);
    
    for(int i=0; i<3; i++)
    {
        pad_PbPb[i]->SetLeftMargin(0);
        pad_PbPb[i]->SetRightMargin(0);
        pad_PbPb[i]->SetTopMargin(0);
        pad_PbPb[i]->SetBottomMargin(0);
        pad_PbPb[i]->Draw();
    }
    
    pad_PbPb[0]->SetLeftMargin(0.16);
    pad_PbPb[1]->SetLeftMargin(0.16);
    pad_PbPb[2]->SetLeftMargin(0.16);
    pad_PbPb[0]->SetRightMargin(0.05);
    pad_PbPb[1]->SetRightMargin(0.05);
    pad_PbPb[2]->SetRightMargin(0.05);
    
    pad_PbPb[0]->SetTopMargin(0.1);
    pad_PbPb[1]->SetTopMargin(0.1);
    pad_PbPb[0]->SetBottomMargin(0.15);
    pad_PbPb[2]->SetBottomMargin(0.42);
    
    //185-250
    pad[0]->cd();
    hist->Draw();
    
    double percent_ks = 0.052;
    double percent_la = 0.052;
    double percent_xi = 0.052;
    double percent_om = 0.052;
    double percent_d2PC_lowpt = 0.109;
    double percent_d2PC_highpt = 0.109;
    double percent_dSP = 0.0;
    
    for(int i=0;i<13;i++)
    {
        double xh,yh,xhsub,yhsub;
        v2ksplot->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;

        //yherr = fabs(yh)*percent_ks;
        yhsuberr = fabs(yhsub)*percent_ks;

    }
    
    for(int i=0;i<10;i++)
    {
        double xh,yh,xhsub,yhsub;
        v2laplot->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;

        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_la;

    }
    
    for(int i=0;i<9;i++)
    {
        double xh,yh,xhsub,yhsub;
        v2xiplot->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_xi;
        
    }

    /*for(int i=0;i<8;i++)
    {
        double xh,yh,xhsub,yhsub;
        v2omplot->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_om;
        
    }*/
    
    for(int i=0;i<2;i++)
    {
        double xh,yh,xhsub,yhsub;
        v2d2PCplot->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_d2PC_lowpt;
        
    }
    for(int i=2;i<6;i++)
    {
        double xh,yh,xhsub,yhsub;
        v2d2PCplot->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_d2PC_highpt;
        
    }


    /*for(int i=0;i<6;i++)
    {
        double xh,yh,xhsub,yhsub;
        v2dSPplot->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_dSP;
        
    }*/
    
    v2ksplot->Draw("PESAME");
    v2laplot->Draw("PESAME");
    v2xiplot->Draw("PESAME");
    //v2omplot->Draw("PESAME");
    v2d2PCplot->Draw("PESAME");
    //v2dSPplot->Draw("PESAME");
    
    //110-150 sub NCQ
    pad[1]->cd();
    hist1->Draw();
    
    for(int i=0;i<12;i++)
    {
        double xh,yh,xhsub,yhsub;
        ksNCQ->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_ks;
        yhsuberr = fabs(yhsub)*percent_ks;
        
    }
    
    for(int i=0;i<10;i++)
    {
        double xh,yh,xhsub,yhsub;
        laNCQ->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_la;
        
    }
    
    for(int i=0;i<9;i++)
    {
        double xh,yh,xhsub,yhsub;
        xiNCQ->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_xi;
        
    }

    /*for(int i=0;i<8;i++)
    {
        double xh,yh,xhsub,yhsub;
        omNCQ->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_om;
        
    }*/

    for(int i=0;i<2;i++)
    {
        double xh,yh,xhsub,yhsub;
        D0v22PCncq->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_d2PC_lowpt;
        
    }
    for(int i=2;i<6;i++)
    {
        double xh,yh,xhsub,yhsub;
        D0v22PCncq->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_d2PC_highpt;
        
    }
    
    /*for(int i=0;i<6;i++)
    {
        double xh,yh,xhsub,yhsub;
        D0v2SPncq->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_dSP;
        
    }*/

    ksNCQ->Draw("PESAME");
    laNCQ->Draw("PESAME");
    xiNCQ->Draw("PESAME");
    //omNCQ->Draw("PESAME");
    D0v22PCncq->Draw("PESAME");
    //D0v2SPncq->Draw("PESAME");
    fitfunc_v2QNS->Draw("LSAME");
    
    //110-150 sub NCQ ratio
    pad[2]->cd();
    hist2->Draw();
    
    for(int i=0;i<12;i++)
    {
        double xh,yh,xhsub,yhsub;
        ksNCQratio->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_ks;
        yhsuberr = fabs(yhsub)*percent_ks;
        
    }
    
    for(int i=0;i<10;i++)
    {
        double xh,yh,xhsub,yhsub;
        laNCQratio->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_la;
        
    }
    
    for(int i=0;i<9;i++)
    {
        double xh,yh,xhsub,yhsub;
        xiNCQratio->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_xi;
        
    }
    
    /*for(int i=0;i<8;i++)
    {
        double xh,yh,xhsub,yhsub;
        omNCQratio->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_om;
        
    }*/
    
    for(int i=0;i<2;i++)
    {
        double xh,yh,xhsub,yhsub;
        d2PCNCQratio->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_d2PC_lowpt;
        
    }
    for(int i=2;i<6;i++)
    {
        double xh,yh,xhsub,yhsub;
        d2PCNCQratio->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_d2PC_highpt;
        
    }
    
    /*for(int i=0;i<6;i++)
    {
        double xh,yh,xhsub,yhsub;
        dSPNCQratio->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_dSP;
        
    }*/

    
    TLine* l = new TLine(0,1,3.2,1);
    l->SetLineStyle(7);
    l->Draw("LSAME");
    ksNCQratio->Draw("PESAME");
    laNCQratio->Draw("PESAME");
    //xiNCQratio->Draw("PESAME");
    //omNCQratio->Draw("PESAME");
    d2PCNCQratio->Draw("PESAME");
    //dSPNCQratio->Draw("PESAME");
    
    ksNCQratio->SetName("kshort");
    laNCQratio->SetName("lambda");
    
    ksNCQratio->Write();
    laNCQratio->Write();

    //return;
    
    TLatex *tex= new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.055);
    pad[0]->cd();
    tex->DrawLatex(0.16,0.93,"#font[61]{CMS} #it{Preliminary}");
    tex->DrawLatex(0.74,0.93,"pPb 8.16TeV");
    pad[1]->cd();
    tex->SetTextSize(0.069);
    tex->DrawLatex(0.16,0.93,"#font[61]{CMS} #it{Preliminary}");
    tex->DrawLatex(0.74,0.93,"pPb 8.16TeV");
    tex->SetTextSize(0.055);
    //pad[1]->cd();
    //tex->SetTextSize(0.052);
    pad[0]->cd();
    tex->DrawLatex(0.18,0.82,"185 #leq N_{trk}^{offline} < 250, |y| < 1");
    pad[1]->cd();
    tex->SetTextSize(0.069);
    tex->DrawLatex(0.18,0.8,"185 #leq N_{trk}^{offline} < 250, |y| < 1");
    tex->SetTextSize(0.055);
    //tex->SetTextSize(0.062);
    //tex->DrawLatex(0.66,0.75,"|#Delta#eta| > 2");
    
    TLegend* leg1 = new TLegend(0.18,0.63,0.5,0.78);
    leg1->SetMargin(0.2);
    leg1->SetFillColor(10);
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0.035);
    leg1->SetTextFont(42);
    leg1->SetTextSize(0.06);
    leg1->AddEntry(v2d2PCplot,"D^{0}","p");
    //leg1->AddEntry(v2d2PCplot,"D^{0}, 2PC","p");
    //leg1->AddEntry(v2dSPplot,"D^{0}, SP","p");
    leg1->AddEntry(v2ksplot,"K^{0}_{S}","p");
    //leg1->AddEntry(v2omplot,"#Omega^{#pm}","p");
    
    TLegend* leg2 = new TLegend(0.29,0.63,0.64,0.78);
    leg2->SetMargin(0.2);
    leg2->SetFillColor(10);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0.035);
    leg2->SetTextFont(42);
    leg2->SetTextSize(0.065);
    leg2->AddEntry(v2laplot,"#Lambda/#bar{#Lambda}","p");
    leg2->AddEntry(v2xiplot,"#Xi^{#pm}","p");

    TLegend* leg11 = new TLegend(0.68,0.65,1,0.85);
    leg11->SetMargin(0.2);
    leg11->SetFillColor(10);
    leg11->SetFillStyle(0);
    leg11->SetBorderSize(0.035);
    leg11->SetTextFont(42);
    leg11->SetTextSize(0.065);
    leg11->AddEntry(v2d2PCplot,"D^{0}","p");
    //leg1->AddEntry(v2d2PCplot,"D^{0}, 2PC","p");
    //leg1->AddEntry(v2dSPplot,"D^{0}, SP","p");
    leg11->AddEntry(v2ksplot,"K^{0}_{S}","p");
    //leg11->AddEntry(v2omplot,"#Omega^{#pm}","p");

    TLegend* leg22 = new TLegend(0.79,0.65,1.11,0.85);
    leg22->SetMargin(0.2);
    leg22->SetFillColor(10);
    leg22->SetFillStyle(0);
    leg22->SetBorderSize(0.035);
    leg22->SetTextFont(42);
    leg22->SetTextSize(0.069);
    leg22->AddEntry(v2laplot,"#Lambda/#bar{#Lambda}","p");
    leg22->AddEntry(v2xiplot,"#Xi^{#pm}","p");


    TLegend *leg1221 = new TLegend(0.54,0.03,0.84,0.13);
    leg1221->SetFillColor(10);
    leg1221->SetFillStyle(0);
    leg1221->SetBorderSize(0.035);
    leg1221->SetTextFont(42);
    leg1221->SetTextSize(0.064);
    leg1221->AddEntry(fitfunc_v2QNS,"Polynomial fits to K_{S}^{0}","L");
    //leg1221->Draw();

    pad[0]->cd();
    leg1->Draw();
    leg2->Draw();
    
    pad[1]->cd();
    //leg1->Draw();
    leg11->Draw();
    leg22->Draw();
    leg1221->Draw();
    
    //============================START PbPb===============================
    TFile* file[2] = TFile::Open("2PC_5TeVv2sigGraphv2_PbPbCent30-50FullStats_deta1p05.root");
    TFile* file[3] = TFile::Open("vn_finalcombinedfit_vnvsmass_MBtrig_SP_cent30to50_poly3bkg_floatwidth_Bfeeddownsys_Alice0_data1_effcorrected0.root");
    //TFile* file[2] = TFile::Open("SP_v2vspt_fromfit_d0_HM185_250_combine_Pbside.root");
    
    TGraphErrors* v2ksplot_PbPb = file[2]->Get("v2kshort");
    TGraphErrors* v2laplot_PbPb = file[2]->Get("v2lambda");
    TGraphErrors* v2xiplot_PbPb = file[2]->Get("v2xi");
    TGraphErrors* v2omplot_PbPb = file[2]->Get("v2omega");
    TGraphErrors* v2d2PCplot_PbPb = file[3]->Get("gr_v2_pt_sys");
    //TGraphErrors* v2dSPplot = file[2]->Get("v2vspt");
    v2d2PCplot_PbPb->RemovePoint(0);
    v2d2PCplot_PbPb->RemovePoint(0);
    
    for(int i=0;i<7;i++)
    {
        v2d2PCplot_PbPb->SetPointError(i,0,v2d2PCplot_PbPb->GetErrorY(i+1));
    }
    
    //make NCQ TGraph
    double ptks[100],ptla[100],ptxi[100],v2ks[100],v2la[100],v2xi[100],v2kse[100],v2lae[100],v2xie[100],KET_ks[100],KET_la[100],KET_xi[100];
    double v2d2PC[100], v2dSP[100], v2d2PCe[100], v2dSPe[100], KET_d2PC[100], KET_dSP[100];
    
    double dMass = 1.865;
    double dv2ncq[10];
    double dv2ncqe[10];
    double dKET[10];
    double dpt[10];
    
    for(int i=0;i<8;i++)
    {
        v2d2PCplot_PbPb->GetPoint(i,dpt[i],dv2ncq[i]);
        dv2ncq[i] = dv2ncq[i]/2.0;
        dv2ncqe[i] = v2d2PCplot_PbPb->GetErrorY(i)/2.0;
        dKET[i] = (sqrt(dMass**2 + dpt[i]**2) - dMass)/2.0;
    }
    
    int Nks = v2ksplot->GetN();
    int Nla = v2laplot->GetN();
    int Nxi = v2xiplot->GetN();
    int Nd = v2d2PCplot->GetN();
    
    double ksMass = 0.497;
    double laMass = 1.1156;
    double xiMass = 1.322;
    
    for(int i=0;i<Nks;i++)
    {
        v2ksplot->GetPoint(i,ptks[i],v2ks[i]);
        v2ks[i] = v2ks[i]/2.0;
        v2kse[i] = v2ksplot->GetErrorY(i)/2.0;
        KET_ks[i] = (sqrt(ksMass**2 + ptks[i]**2) - ksMass)/2.0;
    }
    for(int i=0;i<Nla;i++)
    {
        v2laplot->GetPoint(i,ptla[i],v2la[i]);
        v2la[i] = v2la[i]/3.0;
        v2lae[i] = v2laplot->GetErrorY(i)/3.0;
        KET_la[i] = (sqrt(laMass**2 + ptla[i]**2) - laMass)/3.0;
    }
    for(int i=0;i<Nxi;i++)
    {
        v2xiplot->GetPoint(i,ptxi[i],v2xi[i]);
        v2xi[i] = v2xi[i]/3.0;
        v2xie[i] = v2xiplot->GetErrorY(i)/3.0;
        KET_xi[i] = (sqrt(xiMass**2 + ptxi[i]**2) - xiMass)/3.0;
    }
    
    TGraphErrors* D0v22PCncq_PbPb = new TGraphErrors(8,dKET,dv2ncq,0,dv2ncqe);
    
    //TGraphErrors* D0v22PCncq = file[1].Get("v2vsKET_ncq");
    //TGraphErrors* D0v2SPncq = file[2].Get("v2vsKET_ncq");
    TGraphErrors* ksNCQ_PbPb = file[2]->Get("v2kshort_ket_nq");
    TGraphErrors* laNCQ_PbPb = file[2]->Get("v2lambda_ket_nq");
    TGraphErrors* xiNCQ_PbPb = file[2]->Get("v2xi_ket_nq");
    TGraphErrors* omNCQ_PbPb = file[2]->Get("v2omega_ket_nq");
    
    //make NCQ fit&ratio
    //make NCQ fit&ratio
    fitfunc_v2QNS_PbPb = new TF1("fitfunc_v2QNS_PbPb","([0]/(1+exp(-(x-[1])/[2]))-[3])*pol1(4)",0,3.5);
    fitfunc_v2QNS_PbPb->SetParameters(0.5,-0.1,0.2,0.2,0,0);
    ksNCQ_PbPb->Fit(fitfunc_v2QNS_PbPb,"NOB",0,3.5);
    fitfunc_v2QNS_PbPb->SetLineStyle(2);
    fitfunc_v2QNS_PbPb->SetLineWidth(1);
    
    TGraphErrors* ksNCQratio_PbPb = new TGraphErrors(ksNCQ_PbPb->GetN());
    for(int nn=0;nn<ksNCQ_PbPb->GetN();nn++)
    {
        double xx, yy, xx_err, yy_err;
        ksNCQ_PbPb->GetPoint(nn,xx,yy);
        xx_err = ksNCQ_PbPb->GetErrorX(nn);
        yy_err = ksNCQ_PbPb->GetErrorY(nn);
        ksNCQratio_PbPb->SetPoint(nn,xx,yy/fitfunc_v2QNS_PbPb->Eval(xx));
        ksNCQratio_PbPb->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_v2QNS_PbPb->Eval(xx));
    }
    TGraphErrors* laNCQratio_PbPb = new TGraphErrors(laNCQ_PbPb->GetN());
    for(int nn=0;nn<laNCQ_PbPb->GetN();nn++)
    {
        double xx, yy, xx_err, yy_err;
        laNCQ_PbPb->GetPoint(nn,xx,yy);
        xx_err = laNCQ_PbPb->GetErrorX(nn);
        yy_err = laNCQ_PbPb->GetErrorY(nn);
        laNCQratio_PbPb->SetPoint(nn,xx,yy/fitfunc_v2QNS_PbPb->Eval(xx));
        laNCQratio_PbPb->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_v2QNS_PbPb->Eval(xx));
    }
    TGraphErrors* xiNCQratio_PbPb = new TGraphErrors(xiNCQ_PbPb->GetN());
    for(int nn=0;nn<xiNCQ_PbPb->GetN();nn++)
    {
        double xx, yy, xx_err, yy_err;
        xiNCQ_PbPb->GetPoint(nn,xx,yy);
        xx_err = xiNCQ_PbPb->GetErrorX(nn);
        yy_err = xiNCQ_PbPb->GetErrorY(nn);
        xiNCQratio_PbPb->SetPoint(nn,xx,yy/fitfunc_v2QNS_PbPb->Eval(xx));
        xiNCQratio_PbPb->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_v2QNS_PbPb->Eval(xx));
    }
    TGraphErrors* omNCQratio_PbPb = new TGraphErrors(omNCQ_PbPb->GetN());
    for(int nn=0;nn<omNCQ_PbPb->GetN();nn++)
    {
        double xx, yy, xx_err, yy_err;
        omNCQ_PbPb->GetPoint(nn,xx,yy);
        xx_err = omNCQ_PbPb->GetErrorX(nn);
        yy_err = omNCQ_PbPb->GetErrorY(nn);
        omNCQratio_PbPb->SetPoint(nn,xx,yy/fitfunc_v2QNS_PbPb->Eval(xx));
        omNCQratio_PbPb->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_v2QNS_PbPb->Eval(xx));
    }
    TGraphErrors* d2PCNCQratio_PbPb = new TGraphErrors(D0v22PCncq_PbPb->GetN());
    for(int nn=0;nn<D0v22PCncq_PbPb->GetN();nn++)
    {
        double xx, yy, xx_err, yy_err;
        D0v22PCncq_PbPb->GetPoint(nn,xx,yy);
        xx_err = D0v22PCncq_PbPb->GetErrorX(nn);
        yy_err = D0v22PCncq_PbPb->GetErrorY(nn);
        d2PCNCQratio_PbPb->SetPoint(nn,xx,yy/fitfunc_v2QNS_PbPb->Eval(xx));
        d2PCNCQratio_PbPb->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_v2QNS_PbPb->Eval(xx));
    }
    
    //set line width
    v2ksplot_PbPb->SetLineWidth(1);
    v2laplot_PbPb->SetLineWidth(1);
    v2xiplot_PbPb->SetLineWidth(1);
    v2omplot_PbPb->SetLineWidth(1);
    v2d2PCplot_PbPb->SetLineWidth(1);
    D0v22PCncq_PbPb->SetLineWidth(1);
    ksNCQ_PbPb->SetLineWidth(1);
    laNCQ_PbPb->SetLineWidth(1);
    xiNCQ_PbPb->SetLineWidth(1);
    omNCQ_PbPb->SetLineWidth(1);
    ksNCQratio_PbPb->SetLineWidth(1);
    laNCQratio_PbPb->SetLineWidth(1);
    xiNCQratio_PbPb->SetLineWidth(1);
    d2PCNCQratio_PbPb->SetLineWidth(1);
    
    //set color
    v2ksplot_PbPb->SetMarkerColor(4);
    v2ksplot_PbPb->SetLineColor(4);
    v2ksplot_PbPb->SetMarkerStyle(25);
    ksNCQ_PbPb->SetMarkerColor(4);
    ksNCQ_PbPb->SetLineColor(4);
    ksNCQ_PbPb->SetMarkerStyle(25);
    ksNCQratio_PbPb->SetMarkerColor(4);
    ksNCQratio_PbPb->SetLineColor(4);
    ksNCQratio_PbPb->SetMarkerStyle(25);
    
    v2laplot_PbPb->SetMarkerColor(1);
    v2laplot_PbPb->SetLineColor(1);
    v2laplot_PbPb->SetMarkerStyle(24);
    v2laplot_PbPb->SetMarkerSize(1.4);
    laNCQ_PbPb->SetMarkerColor(1);
    laNCQ_PbPb->SetLineColor(1);
    laNCQ_PbPb->SetMarkerStyle(24);
    laNCQ_PbPb->SetMarkerSize(1.4);
    laNCQratio_PbPb->SetMarkerColor(1);
    laNCQratio_PbPb->SetLineColor(1);
    laNCQratio_PbPb->SetMarkerStyle(24);
    laNCQratio_PbPb->SetMarkerSize(1.4);
    
    v2xiplot_PbPb->SetMarkerColor(1);
    v2xiplot_PbPb->SetLineColor(1);
    v2xiplot_PbPb->SetMarkerStyle(28);
    v2xiplot_PbPb->SetMarkerSize(1.8);
    xiNCQ_PbPb->SetMarkerColor(kGreen+4);
    xiNCQ_PbPb->SetLineColor(kGreen+4);
    xiNCQ_PbPb->SetMarkerStyle(28);
    xiNCQ_PbPb->SetMarkerSize(1.8);
    xiNCQratio_PbPb->SetMarkerColor(kGreen+4);
    xiNCQratio_PbPb->SetLineColor(kGreen+4);
    xiNCQratio_PbPb->SetMarkerStyle(28);
    xiNCQratio_PbPb->SetMarkerSize(1.8);
    
    v2d2PCplot_PbPb->SetMarkerColor(2);
    v2d2PCplot_PbPb->SetLineColor(2);
    v2d2PCplot_PbPb->SetMarkerStyle(20);
    v2d2PCplot_PbPb->SetMarkerSize(1.5);
    D0v22PCncq_PbPb->SetMarkerColor(2);
    D0v22PCncq_PbPb->SetLineColor(2);
    D0v22PCncq_PbPb->SetMarkerStyle(20);
    D0v22PCncq_PbPb->SetMarkerSize(1.5);
    d2PCNCQratio_PbPb->SetMarkerColor(2);
    d2PCNCQratio_PbPb->SetLineColor(2);
    d2PCNCQratio_PbPb->SetMarkerStyle(20);
    d2PCNCQratio_PbPb->SetMarkerSize(1.5);
    
    v2omplot_PbPb->SetMarkerColor(kGreen+2);
    v2omplot_PbPb->SetLineColor(kGreen+2);
    v2omplot_PbPb->SetMarkerStyle(27);
    v2omplot_PbPb->SetMarkerSize(1.4);
    omNCQ_PbPb->SetMarkerColor(kGreen+2);
    omNCQ_PbPb->SetLineColor(kGreen+2);
    omNCQ_PbPb->SetMarkerStyle(27);
    omNCQ_PbPb->SetMarkerStyle(27);
    omNCQ_PbPb->SetMarkerSize(1.4);
    omNCQratio_PbPb->SetMarkerColor(kGreen+2);
    omNCQratio_PbPb->SetLineColor(kGreen+2);
    omNCQratio_PbPb->SetMarkerStyle(27);
    omNCQratio_PbPb->SetMarkerSize(1.4);
    
    //drop points
    //laNCQ->RemovePoint(10);
    //laNCQ->RemovePoint(12);
    //laNCQratio->RemovePoint(10);
    //laNCQratio->RemovePoint(12);
    
    //xiNCQ->RemovePoint(8);
    //xiNCQ->RemovePoint(9);
    //xiNCQratio->RemovePoint(8);
    //xiNCQratio->RemovePoint(9);
    
    //make hist
    TH1D* hist_PbPb = new TH1D("hist_PbPb","",80,0,9.9);
    hist_PbPb->SetLineWidth(0);
    
    fixedFontHist1D(hist_PbPb,2.4,2.2);
    hist_PbPb->SetLineStyle(9);
    hist_PbPb->GetXaxis()->SetTitleSize(hist_PbPb->GetXaxis()->GetTitleSize()*1.);
    hist_PbPb->GetYaxis()->SetTitleSize(hist_PbPb->GetYaxis()->GetTitleSize()*1.6);
    hist_PbPb->GetYaxis()->SetTitleOffset(1.8);
    hist_PbPb->GetXaxis()->SetLabelSize(hist_PbPb->GetXaxis()->GetLabelSize()*1.1);
    hist_PbPb->GetYaxis()->SetLabelSize(hist_PbPb->GetYaxis()->GetLabelSize()*1.2);
    hist_PbPb->SetXTitle("p_{T} (GeV/c)");
    hist_PbPb->SetYTitle("v_{2}");
    hist_PbPb->GetXaxis()->CenterTitle(1);
    hist_PbPb->GetYaxis()->CenterTitle(1);
    hist_PbPb->SetMinimum(0.001);
    hist_PbPb->SetMaximum(0.49);
    hist_PbPb->GetXaxis()->SetRangeUser(0,8.5);
    
    TH1D* hist1_PbPb = new TH1D("hist1_PbPb","",80,-0.1,6);
    hist1_PbPb->SetLineWidth(0);
    
    fixedFontHist1D(hist1_PbPb,2.4,2.2);
    hist1_PbPb->SetLineStyle(9);
    hist1_PbPb->GetXaxis()->CenterTitle(1);
    hist1_PbPb->GetYaxis()->CenterTitle(1);
    hist1_PbPb->GetXaxis()->SetTitleSize(hist1_PbPb->GetXaxis()->GetTitleSize()*1.);
    hist1_PbPb->GetYaxis()->SetTitleSize(hist1_PbPb->GetYaxis()->GetTitleSize()*1.5);
    hist1_PbPb->GetYaxis()->SetTitleOffset(2.0);
    hist1_PbPb->GetXaxis()->SetLabelSize(hist1_PbPb->GetXaxis()->GetLabelSize()*1.2);
    hist1_PbPb->GetYaxis()->SetLabelSize(hist1_PbPb->GetYaxis()->GetLabelSize()*1.2);
    hist1_PbPb->SetXTitle("KE_{T}/n_{q} (GeV)");
    hist1_PbPb->SetYTitle("v_{2}/n_{q}");
    hist1_PbPb->SetMinimum(0.001);
    hist1_PbPb->SetMaximum(0.155);
    hist1_PbPb->GetXaxis()->SetRangeUser(0,3.2);
    
    TH1D* hist2_PbPb = new TH1D("hist2_PbPb","",80,-0.1,6);
    hist2_PbPb->SetLineWidth(0);
    
    fixedFontHist1D(hist2_PbPb,2.4,2.2);
    hist2_PbPb->GetXaxis()->CenterTitle(1);
    hist2_PbPb->GetYaxis()->CenterTitle(1);
    hist2_PbPb->GetXaxis()->SetTitleSize(hist2_PbPb->GetXaxis()->GetTitleSize()*1.);
    hist2_PbPb->GetYaxis()->SetTitleSize(hist2_PbPb->GetYaxis()->GetTitleSize()*1.2);
    hist2_PbPb->GetXaxis()->SetTitleOffset(hist2_PbPb->GetXaxis()->GetTitleOffset()*2.5);
    hist2_PbPb->GetYaxis()->SetTitleOffset(2.4);
    hist2_PbPb->GetXaxis()->SetLabelSize(hist2_PbPb->GetXaxis()->GetLabelSize()*1.2);
    hist2_PbPb->GetYaxis()->SetLabelSize(hist2_PbPb->GetYaxis()->GetLabelSize()*1.2);
    hist2_PbPb->SetXTitle("KE_{T}/n_{q} (GeV)");
    hist2_PbPb->SetYTitle("Data/Fit");
    hist2_PbPb->SetMinimum(0.1);
    hist2_PbPb->SetMaximum(1.9);
    hist2_PbPb->GetXaxis()->SetRangeUser(0,3.2);
    
    //185-250
    pad_PbPb[0]->cd();
    hist_PbPb->Draw();
    
    double percent_ks = 0.052;
    double percent_la = 0.052;
    double percent_xi = 0.052;
    double percent_om = 0.052;
    double percent_d2PC_lowpt = 0.109;
    double percent_d2PC_highpt = 0.109;
    double percent_dSP = 0.0;
    
    double Dsyst[6] = {0.0679,0.0157,0.0157,0.0157,0.0157,0.0157};
    double Dsyst_rel[6];
    
    for(int i=0;i<13;i++)
    {
        double xh,yh,xhsub,yhsub;
        v2ksplot_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_ks;
        yhsuberr = fabs(yhsub)*percent_ks;
        
    }
    
    for(int i=0;i<10;i++)
    {
        double xh,yh,xhsub,yhsub;
        v2laplot_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_la;
        
    }
    
    for(int i=0;i<8;i++)
    {
        double xh,yh,xhsub,yhsub;
        v2xiplot_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_xi;
        
    }
    for(int i=0;i<6;i++)
    {
        double xh,yh,xhsub,yhsub;
        v2omplot_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_om;
        
    }
    
    
    for(int i=0;i<6;i++)
    {
        double xh,yh,xhsub,yhsub;
        v2d2PCplot_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = Dsyst[i];
        Dsyst_rel[i] = Dsyst[i]/yhsub;
        
    }
    /*for(int i=2;i<8;i++)
    {
        double xh,yh,xhsub,yhsub;
        v2d2PCplot_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_d2PC_highpt;
        
    }*/
    
    
    v2ksplot_PbPb->Draw("PESAME");
    v2laplot_PbPb->Draw("PESAME");
    v2xiplot_PbPb->Draw("PESAME");
    v2omplot_PbPb->Draw("PESAME");
    v2d2PCplot_PbPb->Draw("PESAME");
    
    //110-150 sub NCQ
    pad_PbPb[1]->cd();
    hist1_PbPb->Draw();
    
    for(int i=0;i<13;i++)
    {
        double xh,yh,xhsub,yhsub;
        ksNCQ_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_ks;
        yhsuberr = fabs(yhsub)*percent_ks;
        
    }
    
    for(int i=0;i<10;i++)
    {
        double xh,yh,xhsub,yhsub;
        laNCQ_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_la;
        
    }
    
    for(int i=0;i<8;i++)
    {
        double xh,yh,xhsub,yhsub;
        xiNCQ_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_xi;
        
    }
    
    for(int i=0;i<6;i++)
    {
        double xh,yh,xhsub,yhsub;
        omNCQ_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_om;
        
    }
    
    for(int i=0;i<6;i++)
    {
        double xh,yh,xhsub,yhsub;
        D0v22PCncq_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = Dsyst[i]/2.0;
        
    }
    /*for(int i=2;i<8;i++)
    {
        double xh,yh,xhsub,yhsub;
        D0v22PCncq_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_d2PC_highpt;
        
    }*/
    
    ksNCQ_PbPb->Draw("PESAME");
    laNCQ_PbPb->Draw("PESAME");
    xiNCQ_PbPb->Draw("PESAME");
    omNCQ_PbPb->Draw("PESAME");
    D0v22PCncq_PbPb->Draw("PESAME");
    //D0v2SPncq->Draw("PESAME");
    fitfunc_v2QNS_PbPb->Draw("LSAME");
    
    //110-150 sub NCQ ratio
    pad_PbPb[2]->cd();
    hist2_PbPb->Draw();
    
    for(int i=1;i<13;i++)
    {
        double xh,yh,xhsub,yhsub;
        ksNCQratio_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_ks;
        yhsuberr = fabs(yhsub)*percent_ks;
        
    }
    
    for(int i=0;i<10;i++)
    {
        double xh,yh,xhsub,yhsub;
        laNCQratio_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_la;
        
    }
    
    for(int i=0;i<8;i++)
    {
        double xh,yh,xhsub,yhsub;
        xiNCQratio_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_xi;
        
    }
    for(int i=0;i<6;i++)
    {
        double xh,yh,xhsub,yhsub;
        omNCQratio_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_om;
        
    }
    
    for(int i=0;i<6;i++)
    {
        double xh,yh,xhsub,yhsub;
        d2PCNCQratio_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*Dsyst_rel[i];
        
    }
    /*for(int i=2;i<8;i++)
    {
        double xh,yh,xhsub,yhsub;
        d2PCNCQratio_PbPb->GetPoint(i,xh,yhsub);
        double yherr,yhsuberr;
        
        //yherr = fabs(yh)*percent_la;
        yhsuberr = fabs(yhsub)*percent_d2PC_highpt;
        
    }*/
    
    TLine* l = new TLine(0,1,3.2,1);
    l->SetLineStyle(7);
    l->Draw("LSAME");
    ksNCQratio_PbPb->Draw("PESAME");
    laNCQratio_PbPb->Draw("PESAME");
    xiNCQratio_PbPb->Draw("PESAME");
    omNCQratio_PbPb->Draw("PESAME");
    d2PCNCQratio_PbPb->Draw("PESAME");
    //dSPNCQratio->Draw("PESAME");
    
    //return;
    
    TLatex *tex= new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.055);
    pad_PbPb[0]->cd();
    tex->DrawLatex(0.16,0.93,"#font[61]{CMS} #it{Preliminary}");
    tex->DrawLatex(0.72,0.93,"PbPb 5.02TeV");
    pad_PbPb[1]->cd();
    tex->SetTextSize(0.069);
    tex->DrawLatex(0.16,0.93,"#font[61]{CMS} #it{Preliminary}");
    tex->DrawLatex(0.72,0.93,"PbPb 5.02TeV");
    tex->SetTextSize(0.055);
    //pad[1]->cd();
    //tex->SetTextSize(0.052);
    pad_PbPb[0]->cd();
    tex->DrawLatex(0.18,0.82,"Centrality 30-50%, |y| < 1");
    pad_PbPb[1]->cd();
    tex->SetTextSize(0.069);
    tex->DrawLatex(0.18,0.8,"Centrality 30-50%, |y| < 1");
    tex->SetTextSize(0.055);
    //tex->SetTextSize(0.062);
    //tex->DrawLatex(0.66,0.75,"|#Delta#eta| > 2");
    
    TLegend* leg1_PbPb = new TLegend(0.18,0.555,0.5,0.78);
    leg1_PbPb->SetMargin(0.2);
    leg1_PbPb->SetFillColor(10);
    leg1_PbPb->SetFillStyle(0);
    leg1_PbPb->SetBorderSize(0.035);
    leg1_PbPb->SetTextFont(42);
    leg1_PbPb->SetTextSize(0.06);
    leg1_PbPb->AddEntry(v2d2PCplot_PbPb,"D^{0}","p");
    //leg1->AddEntry(v2d2PCplot,"D^{0}, 2PC","p");
    //leg1->AddEntry(v2dSPplot,"D^{0}, SP","p");
    leg1_PbPb->AddEntry(v2ksplot_PbPb,"K^{0}_{S}","p");
    leg1_PbPb->AddEntry(v2omplot_PbPb,"#Omega^{#pm}","p");
    
    TLegend* leg2_PbPb = new TLegend(0.29,0.63,0.64,0.78);
    leg2_PbPb->SetMargin(0.2);
    leg2_PbPb->SetFillColor(10);
    leg2_PbPb->SetFillStyle(0);
    leg2_PbPb->SetBorderSize(0.035);
    leg2_PbPb->SetTextFont(42);
    leg2_PbPb->SetTextSize(0.065);
    leg2_PbPb->AddEntry(v2laplot_PbPb,"#Lambda/#bar{#Lambda}","p");
    leg2_PbPb->AddEntry(v2xiplot_PbPb,"#Xi^{#pm}","p");
    
    TLegend* leg11_PbPb = new TLegend(0.68,0.55,1,0.85);
    leg11_PbPb->SetMargin(0.2);
    leg11_PbPb->SetFillColor(10);
    leg11_PbPb->SetFillStyle(0);
    leg11_PbPb->SetBorderSize(0.035);
    leg11_PbPb->SetTextFont(42);
    leg11_PbPb->SetTextSize(0.065);
    leg11_PbPb->AddEntry(v2d2PCplot_PbPb,"D^{0}","p");
    //leg1->AddEntry(v2d2PCplot,"D^{0}, 2PC","p");
    //leg1->AddEntry(v2dSPplot,"D^{0}, SP","p");
    leg11_PbPb->AddEntry(v2ksplot_PbPb,"K^{0}_{S}","p");
    leg11_PbPb->AddEntry(v2omplot_PbPb,"#Omega^{#pm}","p");
    
    TLegend* leg22_PbPb = new TLegend(0.79,0.65,1.11,0.85);
    leg22_PbPb->SetMargin(0.2);
    leg22_PbPb->SetFillColor(10);
    leg22_PbPb->SetFillStyle(0);
    leg22_PbPb->SetBorderSize(0.035);
    leg22_PbPb->SetTextFont(42);
    leg22_PbPb->SetTextSize(0.069);
    leg22_PbPb->AddEntry(v2laplot_PbPb,"#Lambda/#bar{#Lambda}","p");
    leg22_PbPb->AddEntry(v2xiplot_PbPb,"#Xi^{#pm}","p");
    
    
    TLegend *leg1221_PbPb = new TLegend(0.54,0.03,0.84,0.13);
    leg1221_PbPb->SetFillColor(10);
    leg1221_PbPb->SetFillStyle(0);
    leg1221_PbPb->SetBorderSize(0.035);
    leg1221_PbPb->SetTextFont(42);
    leg1221_PbPb->SetTextSize(0.064);
    leg1221_PbPb->AddEntry(fitfunc_v2QNS_PbPb,"Polynomial fits to K_{S}^{0}","L");
    //leg1221->Draw();
    
    pad_PbPb[0]->cd();
    leg1_PbPb->Draw();
    leg2_PbPb->Draw();
    
    pad_PbPb[1]->cd();
    //leg1->Draw();
    leg11_PbPb->Draw();
    leg22_PbPb->Draw();
    leg1221_PbPb->Draw();
    
    c2->Print("v2vspt_NCQ_deta2_sub.pdf");
    //c2->Print("v2v3_pt_NCQsub_totalsyst.gif");
}
