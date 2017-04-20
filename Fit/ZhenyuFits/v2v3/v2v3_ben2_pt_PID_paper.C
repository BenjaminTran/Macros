#include "GetGraphFromFile.C"
#include "makeMultiPanelCanvas.C"
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

#include <vector>

void v2v3_ben2_pt_PID_paper()
{

// Results of charged hadrons from HIN-13-002
   TGraphErrors* pPb_v2[7];
   TGraphErrors* PbPb_v2[7];
   TGraphErrors* pPb_v3[7];
   TGraphErrors* PbPb_v3[7];

    TFile* file_pPb_hadron[4];
    file_pPb_hadron[0] = TFile::Open("../MultiStudyV2/rootfile/pPb/v15/lrgraphv2_v3_pPb_hadron_0-35.root");
    file_pPb_hadron[1] = TFile::Open("../MultiStudyV2/rootfile/pPb/v15/lrgraphv2_v3_pPb_hadron_35-60.root");
    file_pPb_hadron[2] = TFile::Open("../MultiStudyV2/rootfile/pPb/v15/lrgraphv2_v3_pPb_hadron_60-120.root");
    file_pPb_hadron[3] = TFile::Open("../MultiStudyV2/rootfile/pPb/v15/lrgraphv2_v3_v4_pPb_hadron_185-above.root");

    TGraphErrors* pPb_v2[0] = (TGraphErrors*)file_pPb_hadron[0]->Get("hadronv2");
    TGraphErrors* pPb_v2[1] = (TGraphErrors*)file_pPb_hadron[1]->Get("hadronv2");
    TGraphErrors* pPb_v2[2] = (TGraphErrors*)file_pPb_hadron[2]->Get("hadronv2");
    
    pPb_v2[0]->SetMarkerStyle(28);
    pPb_v2[1]->SetMarkerStyle(28);
    pPb_v2[2]->SetMarkerStyle(28);
    
    TFile* file_PbPb_hadron[4];
    file_PbPb_hadron[0] = TFile::Open("../MultiStudyV2/rootfile/PbPb/v15/lrgraphv2_v3_PbPb_hadron_0-35.root");
    file_PbPb_hadron[1] = TFile::Open("../MultiStudyV2/rootfile/PbPb/v15/lrgraphv2_v3_PbPb_hadron_35-60.root");
    file_PbPb_hadron[2] = TFile::Open("../MultiStudyV2/rootfile/PbPb/v15/lrgraphv2_v3_PbPb_hadron_60-120.root");
    file_PbPb_hadron[3] = TFile::Open("../MultiStudyV2/rootfile/PbPb/v15/lrgraphv2_v3_v4_PbPb_hadron_185-above.root");
    
    TGraphErrors* PbPb_v2[0] = (TGraphErrors*)file_PbPb_hadron[0]->Get("hadronv2");
    TGraphErrors* PbPb_v2[1] = (TGraphErrors*)file_PbPb_hadron[1]->Get("hadronv2");
    TGraphErrors* PbPb_v2[2] = (TGraphErrors*)file_PbPb_hadron[2]->Get("hadronv2");
    
    PbPb_v2[0]->SetMarkerStyle(28);
    PbPb_v2[1]->SetMarkerStyle(28);
    PbPb_v2[2]->SetMarkerStyle(28);

    
   pPb_v2[3] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n120_150_ptass033pPb_v2.txt"),1,28,1.2);
   pPb_v2[4] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n150_185_ptass033pPb_v2.txt"),1,28,1.2);
   pPb_v2[5] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033pPb_v2.txt"),1,28,1.2);
   pPb_v2[6] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n220_260_ptass033pPb_v2.txt"),1,28,1.2);

   PbPb_v2[3] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n120_150_ptass033_v2.txt"),1,28,1.2);
   PbPb_v2[4] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n150_185_ptass033_v2.txt"),1,28,1.2);
   PbPb_v2[5] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033_v2.txt"),1,28,1.2);
   PbPb_v2[6] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n220_260_ptass033_v2.txt"),1,28,1.2);

   pPb_v3[3] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n120_150_ptass033pPb_v3.txt"),1,28,1.2);
   pPb_v3[4] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n150_185_ptass033pPb_v3.txt"),1,28,1.2);
   //pPb_v3[5] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033pPb_v3.txt"),1,28,1.2);
   pPb_v3[6] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n220_260_ptass033pPb_v3.txt"),1,28,1.2);
    
    TGraphErrors* pPb_v3[5] = (TGraphErrors*)file_pPb_hadron[3]->Get("hadronv3");
    pPb_v3[5]->SetMarkerStyle(28);
    
   PbPb_v3[3] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n120_150_ptass033_v3.txt"),1,28,1.2);
   PbPb_v3[4] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n150_185_ptass033_v3.txt"),1,28,1.2);
   PbPb_v3[5] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n185_220_ptass033_v3.txt"),1,28,1.2);
   PbPb_v3[6] = (TGraphErrors*)GetGraphWithSymmYErrorsFromFile(Form("data/%s","n220_260_ptass033_v3.txt"),1,28,1.2);

   TFile* file_pPb[8];
   TFile* file_PbPb[8];

   file_pPb[0] = TFile::Open("../MultiStudyV2/rootfile/pPb/v15/lrgraphv2_v3_pPb_0-35.root");
   file_pPb[1] = TFile::Open("../MultiStudyV2/rootfile/pPb/v15/lrgraphv2_v3_pPb_35-60.root");
   file_pPb[2] = TFile::Open("../MultiStudyV2/rootfile/pPb/v15/lrgraphv2_v3_pPb_60-120.root");
   file_pPb[3] = TFile::Open("../MultiStudyV2/rootfile/pPb/v15/lrgraphv2_v3_pPb_120-150.root");
   file_pPb[4] = TFile::Open("../MultiStudyV2/rootfile/pPb/v15/lrgraphv2_v3_pPb_150-185.root");
   file_pPb[5] = TFile::Open("../MultiStudyV2/rootfile/pPb/v15/lrgraphv2_v3_pPb_185-220.root");
   file_pPb[6] = TFile::Open("../MultiStudyV2/rootfile/pPb/v15/lrgraphv2_v3_pPb_220-260.root");
   file_pPb[7] = TFile::Open("../MultiStudyV2/rootfile/pPb/v15/lrgraphv2_v3_pPb_185-above.root");

   file_PbPb[0] = TFile::Open("../MultiStudyV2/rootfile/PbPb/v15/lrgraphv2_v3_PbPb_0-35.root");
   file_PbPb[1] = TFile::Open("../MultiStudyV2/rootfile/PbPb/v15/lrgraphv2_v3_PbPb_35-60.root");
   file_PbPb[2] = TFile::Open("../MultiStudyV2/rootfile/PbPb/v15/lrgraphv2_v3_PbPb_60-120.root");
   file_PbPb[3] = TFile::Open("../MultiStudyV2/rootfile/PbPb/v15/lrgraphv2_v3_PbPb_120-150.root");
   file_PbPb[4] = TFile::Open("../MultiStudyV2/rootfile/PbPb/v15/lrgraphv2_v3_PbPb_150-185.root");
   file_PbPb[5] = TFile::Open("../MultiStudyV2/rootfile/PbPb/v15/lrgraphv2_v3_PbPb_185-220.root");
   file_PbPb[6] = TFile::Open("../MultiStudyV2/rootfile/PbPb/v15/lrgraphv2_v3_PbPb_220-260.root");
   file_PbPb[7] = TFile::Open("../MultiStudyV2/rootfile/PbPb/v15/lrgraphv2_v3_PbPb_220-260.root");

   TGraphErrors* ks_pPb_v2[8];
   TGraphErrors* lambda_pPb_v2[8];
   TGraphErrors* ks_PbPb_v2[8];
   TGraphErrors* lambda_PbPb_v2[8];  

   TGraphErrors* ks_pPb_v3[8];
   TGraphErrors* lambda_pPb_v3[8];
   TGraphErrors* ks_PbPb_v3[8];
   TGraphErrors* lambda_PbPb_v3[8];

   TGraphErrors* ks_pPb_v2QNS[8];
   TGraphErrors* lambda_pPb_v2QNS[8];
   TGraphErrors* ks_PbPb_v2QNS[8];
   TGraphErrors* lambda_PbPb_v2QNS[8];

   TGraphErrors* ks_pPb_v3QNS[8];
   TGraphErrors* lambda_pPb_v3QNS[8];
   TGraphErrors* ks_PbPb_v3QNS[8];
   TGraphErrors* lambda_PbPb_v3QNS[8];

   TMultiGraph*  pPb_v2QNS[8];
   TMultiGraph*  PbPb_v2QNS[8];
   TMultiGraph*  pPb_v3QNS[8];
   TMultiGraph*  PbPb_v3QNS[8];
   TF1*          fitfunc_pPb_v2QNS[8];
   TF1*          fitfunc_PbPb_v2QNS[8];
   TF1*          fitfunc_pPb_v3QNS[8];
   TF1*          fitfunc_PbPb_v3QNS[8];

   TGraphErrors* ks_pPb_v2QNS_ratio[8];
   TGraphErrors* lambda_pPb_v2QNS_ratio[8];
   TGraphErrors* ks_PbPb_v2QNS_ratio[8];
   TGraphErrors* lambda_PbPb_v2QNS_ratio[8];

   TGraphErrors* ks_pPb_v3QNS_ratio[8];
   TGraphErrors* lambda_pPb_v3QNS_ratio[8];
   TGraphErrors* ks_PbPb_v3QNS_ratio[8];
   TGraphErrors* lambda_PbPb_v3QNS_ratio[8];

   for(int i=0;i<8;i++)
   {
     ks_pPb_v2[i] = (TGraphErrors*)file_pPb[i]->Get("kshortv2true");
     ks_PbPb_v2[i] = (TGraphErrors*)file_PbPb[i]->Get("kshortv2true");
     lambda_pPb_v2[i] = (TGraphErrors*)file_pPb[i]->Get("lambdav2true");
     lambda_PbPb_v2[i] = (TGraphErrors*)file_PbPb[i]->Get("lambdav2true");

     ks_pPb_v3[i] = (TGraphErrors*)file_pPb[i]->Get("kshortv3true");
     ks_PbPb_v3[i] = (TGraphErrors*)file_PbPb[i]->Get("kshortv3true");
     lambda_pPb_v3[i] = (TGraphErrors*)file_pPb[i]->Get("lambdav3true");
     lambda_PbPb_v3[i] = (TGraphErrors*)file_PbPb[i]->Get("lambdav3true");

     ks_pPb_v2QNS[i] = (TGraphErrors*)file_pPb[i]->Get("kshortv2trueQNS");
     ks_PbPb_v2QNS[i] = (TGraphErrors*)file_PbPb[i]->Get("kshortv2trueQNS");
     lambda_pPb_v2QNS[i] = (TGraphErrors*)file_pPb[i]->Get("lambdav2trueQNS");
     lambda_PbPb_v2QNS[i] = (TGraphErrors*)file_PbPb[i]->Get("lambdav2trueQNS");

     ks_pPb_v3QNS[i] = (TGraphErrors*)file_pPb[i]->Get("kshortv3trueQNS");
     ks_PbPb_v3QNS[i] = (TGraphErrors*)file_PbPb[i]->Get("kshortv3trueQNS");
     lambda_pPb_v3QNS[i] = (TGraphErrors*)file_pPb[i]->Get("lambdav3trueQNS");
     lambda_PbPb_v3QNS[i] = (TGraphErrors*)file_PbPb[i]->Get("lambdav3trueQNS");

     pPb_v2QNS[i] = new TMultiGraph();
     pPb_v2QNS[i]->Add(ks_pPb_v2QNS[i]);
     pPb_v2QNS[i]->Add(lambda_pPb_v2QNS[i]);
     PbPb_v2QNS[i] = new TMultiGraph();
     PbPb_v2QNS[i]->Add(ks_PbPb_v2QNS[i]);
     PbPb_v2QNS[i]->Add(lambda_PbPb_v2QNS[i]);
     pPb_v3QNS[i] = new TMultiGraph();
     pPb_v3QNS[i]->Add(ks_pPb_v3QNS[i]);
     pPb_v3QNS[i]->Add(lambda_pPb_v3QNS[i]);
     PbPb_v3QNS[i] = new TMultiGraph();
     PbPb_v3QNS[i]->Add(ks_PbPb_v3QNS[i]);
     PbPb_v3QNS[i]->Add(lambda_PbPb_v3QNS[i]);

     // QNS scaling ratio fit

//     fitfunc_pPb_v2QNS[i] = new TF1(Form("fitfunc_pPb_v2QNS_%d",i),"pol3",0,2.0);
//     fitfunc_pPb_v2QNS[i]->FixParameter(0,0);
     fitfunc_pPb_v2QNS[i] = new TF1(Form("fitfunc_pPb_v2QNS_%d",i),"([0]/(1+exp(-(x-[1])/[2]))-[3])*pol1(4)",0,2.0);
     fitfunc_pPb_v2QNS[i]->SetParameters(0.5,-0.1,0.2,0.2,0,0);
//     fitfunc_pPb_v2QNS[i]->SetParLimits(1,-10,0);
//     fitfunc_pPb_v2QNS[i]->FixParameter(4,0);
     ks_pPb_v2QNS[i]->Fit(fitfunc_pPb_v2QNS[i],"NOB",0,2.0);
cout<<"pPb: n="<<i<<endl;
     fitfunc_pPb_v2QNS[i]->SetLineStyle(2);
     fitfunc_pPb_v2QNS[i]->SetLineWidth(1);
     ks_pPb_v2QNS_ratio[i] = new TGraphErrors(ks_pPb_v2QNS[i]->GetN());
     for(int nn=0;nn<ks_pPb_v2QNS_ratio[i]->GetN();nn++)
     {
       double xx, yy, xx_err, yy_err;
       ks_pPb_v2QNS[i]->GetPoint(nn,xx,yy);
       xx_err = ks_pPb_v2QNS[i]->GetErrorX(nn);
       yy_err = ks_pPb_v2QNS[i]->GetErrorY(nn);
       ks_pPb_v2QNS_ratio[i]->SetPoint(nn,xx,yy/fitfunc_pPb_v2QNS[i]->Eval(xx));
       ks_pPb_v2QNS_ratio[i]->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_pPb_v2QNS[i]->Eval(xx));          
     }
     lambda_pPb_v2QNS_ratio[i] = new TGraphErrors(lambda_pPb_v2QNS[i]->GetN());
     for(int nn=0;nn<lambda_pPb_v2QNS_ratio[i]->GetN();nn++)
     {
       double xx, yy, xx_err, yy_err;
       lambda_pPb_v2QNS[i]->GetPoint(nn,xx,yy);
       xx_err = lambda_pPb_v2QNS[i]->GetErrorX(nn);
       yy_err = lambda_pPb_v2QNS[i]->GetErrorY(nn);
       lambda_pPb_v2QNS_ratio[i]->SetPoint(nn,xx,yy/fitfunc_pPb_v2QNS[i]->Eval(xx));
       lambda_pPb_v2QNS_ratio[i]->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_pPb_v2QNS[i]->Eval(xx));
     }

//     fitfunc_PbPb_v2QNS[i] = new TF1(Form("fitfunc_PbPb_v2QNS_%d",i),"pol3",0,2.0);
     fitfunc_PbPb_v2QNS[i] = new TF1(Form("fitfunc_PbPb_v2QNS_%d",i),"([0]/(1+exp(-(x-[1])**2/[2]))-[3])*pol1(4)",0,2.0);
     fitfunc_PbPb_v2QNS[i]->SetParameters(0.2,-0.1,0.2,0.2);
//     fitfunc_PbPb_v2QNS[i]->SetParLimits(1,-10,0);
//     fitfunc_PbPb_v2QNS[i]->FixParameter(4,0);
     if(i==5) ks_PbPb_v2QNS[i]->Fit(fitfunc_PbPb_v2QNS[i],"NO",0,1.8);
     else if(i==4) ks_PbPb_v2QNS[i]->Fit(fitfunc_PbPb_v2QNS[i],"NO",0,1.5); 
     else ks_PbPb_v2QNS[i]->Fit(fitfunc_PbPb_v2QNS[i],"NO",0,1.8);

     fitfunc_PbPb_v2QNS[i]->SetLineStyle(2);
     fitfunc_PbPb_v2QNS[i]->SetLineWidth(1);
     ks_PbPb_v2QNS_ratio[i] = new TGraphErrors(ks_PbPb_v2QNS[i]->GetN());
     for(int nn=0;nn<ks_PbPb_v2QNS_ratio[i]->GetN();nn++)
     {
       double xx, yy, xx_err, yy_err;
       ks_PbPb_v2QNS[i]->GetPoint(nn,xx,yy);
       xx_err = ks_PbPb_v2QNS[i]->GetErrorX(nn);
       yy_err = ks_PbPb_v2QNS[i]->GetErrorY(nn);
       ks_PbPb_v2QNS_ratio[i]->SetPoint(nn,xx,yy/fitfunc_PbPb_v2QNS[i]->Eval(xx));
       ks_PbPb_v2QNS_ratio[i]->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_PbPb_v2QNS[i]->Eval(xx));
     }
     lambda_PbPb_v2QNS_ratio[i] = new TGraphErrors(lambda_PbPb_v2QNS[i]->GetN());
     for(int nn=0;nn<lambda_PbPb_v2QNS_ratio[i]->GetN();nn++)
     {
       double xx, yy, xx_err, yy_err;
       lambda_PbPb_v2QNS[i]->GetPoint(nn,xx,yy);
       xx_err = lambda_PbPb_v2QNS[i]->GetErrorX(nn);
       yy_err = lambda_PbPb_v2QNS[i]->GetErrorY(nn);
       lambda_PbPb_v2QNS_ratio[i]->SetPoint(nn,xx,yy/fitfunc_PbPb_v2QNS[i]->Eval(xx));
       lambda_PbPb_v2QNS_ratio[i]->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_PbPb_v2QNS[i]->Eval(xx));
     }

//     fitfunc_pPb_v3QNS[i] = new TF1(Form("fitfunc_pPb_v3QNS_%d",i),"pol3",0,2.0);
//     fitfunc_pPb_v3QNS[i]->FixParameter(0,0);
     fitfunc_pPb_v3QNS[i] = new TF1(Form("fitfunc_pPb_v3QNS_%d",i),"([0]/(1+exp(-(x-[1])**3.4/[2]))-[3])*pol1(4)",0.0,2.0);
     fitfunc_pPb_v3QNS[i]->SetParameters(0.2,-0.1,0.2,0.2);
//     fitfunc_pPb_v3QNS[i]->SetParLimits(1,-10,0);
//     fitfunc_pPb_v3QNS[i]->FixParameter(4,0.00000);
     ks_pPb_v3QNS[i]->Fit(fitfunc_pPb_v3QNS[i],"NO","",0.0,1.9);
     fitfunc_pPb_v3QNS[i]->SetLineStyle(2);
     fitfunc_pPb_v3QNS[i]->SetLineWidth(1);
     ks_pPb_v3QNS_ratio[i] = new TGraphErrors(ks_pPb_v3QNS[i]->GetN());
     for(int nn=0;nn<ks_pPb_v3QNS_ratio[i]->GetN();nn++)
     {
       double xx, yy, xx_err, yy_err;
       ks_pPb_v3QNS[i]->GetPoint(nn,xx,yy);
       xx_err = ks_pPb_v3QNS[i]->GetErrorX(nn);
       yy_err = ks_pPb_v3QNS[i]->GetErrorY(nn);
       ks_pPb_v3QNS_ratio[i]->SetPoint(nn,xx,yy/fitfunc_pPb_v3QNS[i]->Eval(xx));
       ks_pPb_v3QNS_ratio[i]->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_pPb_v3QNS[i]->Eval(xx));
     }
     lambda_pPb_v3QNS_ratio[i] = new TGraphErrors(lambda_pPb_v3QNS[i]->GetN());
     for(int nn=0;nn<lambda_pPb_v3QNS_ratio[i]->GetN();nn++)
     {
       double xx, yy, xx_err, yy_err;
       lambda_pPb_v3QNS[i]->GetPoint(nn,xx,yy);
       xx_err = lambda_pPb_v3QNS[i]->GetErrorX(nn);
       yy_err = lambda_pPb_v3QNS[i]->GetErrorY(nn);
       lambda_pPb_v3QNS_ratio[i]->SetPoint(nn,xx,yy/fitfunc_pPb_v3QNS[i]->Eval(xx));
       lambda_pPb_v3QNS_ratio[i]->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_pPb_v3QNS[i]->Eval(xx));
     }

     fitfunc_PbPb_v3QNS[i] = new TF1(Form("fitfunc_PbPb_v3QNS_%d",i),"pol3",0,2.0);
     fitfunc_PbPb_v3QNS[i]->FixParameter(0,0);
     ks_PbPb_v3QNS[i]->Fit(fitfunc_PbPb_v3QNS[i],"RNOB");
     fitfunc_PbPb_v3QNS[i]->SetLineStyle(2);
     fitfunc_PbPb_v3QNS[i]->SetLineWidth(1);
     ks_PbPb_v3QNS_ratio[i] = new TGraphErrors(ks_PbPb_v3QNS[i]->GetN());
     for(int nn=0;nn<ks_PbPb_v3QNS_ratio[i]->GetN();nn++)
     {
       double xx, yy, xx_err, yy_err;
       ks_PbPb_v3QNS[i]->GetPoint(nn,xx,yy);
       xx_err = ks_PbPb_v3QNS[i]->GetErrorX(nn);
       yy_err = ks_PbPb_v3QNS[i]->GetErrorY(nn);
       ks_PbPb_v3QNS_ratio[i]->SetPoint(nn,xx,yy/fitfunc_PbPb_v3QNS[i]->Eval(xx));
       ks_PbPb_v3QNS_ratio[i]->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_PbPb_v3QNS[i]->Eval(xx));
     }
     lambda_PbPb_v3QNS_ratio[i] = new TGraphErrors(lambda_PbPb_v3QNS[i]->GetN());
     for(int nn=0;nn<lambda_PbPb_v3QNS_ratio[i]->GetN();nn++)
     {
       double xx, yy, xx_err, yy_err;
       lambda_PbPb_v3QNS[i]->GetPoint(nn,xx,yy);
       xx_err = lambda_PbPb_v3QNS[i]->GetErrorX(nn);
       yy_err = lambda_PbPb_v3QNS[i]->GetErrorY(nn);
       lambda_PbPb_v3QNS_ratio[i]->SetPoint(nn,xx,yy/fitfunc_PbPb_v3QNS[i]->Eval(xx));
       lambda_PbPb_v3QNS_ratio[i]->SetPointError(nn,xx_err,yy_err/yy*yy/fitfunc_PbPb_v3QNS[i]->Eval(xx));
     }

     ks_pPb_v2[i]->SetMarkerColor(4);
     ks_PbPb_v2[i]->SetMarkerColor(4);
     lambda_pPb_v2[i]->SetMarkerColor(2);
     lambda_PbPb_v2[i]->SetMarkerColor(2);
       
       ks_pPb_v2[i]->SetLineColor(4);
       ks_PbPb_v2[i]->SetLineColor(4);
       lambda_pPb_v2[i]->SetLineColor(2);
       lambda_PbPb_v2[i]->SetLineColor(2);
       
     ks_pPb_v2[i]->SetMarkerStyle(21);
     ks_PbPb_v2[i]->SetMarkerStyle(21);
     lambda_pPb_v2[i]->SetMarkerStyle(20);
     lambda_PbPb_v2[i]->SetMarkerStyle(20); 
     ks_pPb_v2[i]->SetLineWidth(1);
     ks_PbPb_v2[i]->SetLineWidth(1);
     lambda_pPb_v2[i]->SetLineWidth(1);
     lambda_PbPb_v2[i]->SetLineWidth(1);

     ks_pPb_v3[i]->SetMarkerColor(4);
     ks_PbPb_v3[i]->SetMarkerColor(4);
     lambda_pPb_v3[i]->SetMarkerColor(2);
     lambda_PbPb_v3[i]->SetMarkerColor(2);
       
       ks_pPb_v3[i]->SetLineColor(4);
       ks_PbPb_v3[i]->SetLineColor(4);
       lambda_pPb_v3[i]->SetLineColor(2);
       lambda_PbPb_v3[i]->SetLineColor(2);

     ks_pPb_v3[i]->SetMarkerStyle(21);
     ks_PbPb_v3[i]->SetMarkerStyle(21);
     lambda_pPb_v3[i]->SetMarkerStyle(20);
     lambda_PbPb_v3[i]->SetMarkerStyle(20);
     ks_pPb_v3[i]->SetLineWidth(1);
     ks_PbPb_v3[i]->SetLineWidth(1);
     lambda_pPb_v3[i]->SetLineWidth(1);
     lambda_PbPb_v3[i]->SetLineWidth(1);

     ks_pPb_v2QNS[i]->SetMarkerColor(4);
     ks_PbPb_v2QNS[i]->SetMarkerColor(4);
     lambda_pPb_v2QNS[i]->SetMarkerColor(2);
     lambda_PbPb_v2QNS[i]->SetMarkerColor(2);
       
       ks_pPb_v2QNS[i]->SetLineColor(4);
       ks_PbPb_v2QNS[i]->SetLineColor(4);
       lambda_pPb_v2QNS[i]->SetLineColor(2);
       lambda_PbPb_v2QNS[i]->SetLineColor(2);

     ks_pPb_v2QNS[i]->SetMarkerStyle(21);
     ks_PbPb_v2QNS[i]->SetMarkerStyle(21);
     lambda_pPb_v2QNS[i]->SetMarkerStyle(20);
     lambda_PbPb_v2QNS[i]->SetMarkerStyle(20);
     ks_pPb_v2QNS[i]->SetLineWidth(1);
     ks_PbPb_v2QNS[i]->SetLineWidth(1);
     lambda_pPb_v2QNS[i]->SetLineWidth(1);
     lambda_PbPb_v2QNS[i]->SetLineWidth(1);

     ks_pPb_v3QNS[i]->SetMarkerColor(4);
     ks_PbPb_v3QNS[i]->SetMarkerColor(4);
     lambda_pPb_v3QNS[i]->SetMarkerColor(2);
     lambda_PbPb_v3QNS[i]->SetMarkerColor(2);
       
       ks_pPb_v3QNS[i]->SetLineColor(4);
       ks_PbPb_v3QNS[i]->SetLineColor(4);
       lambda_pPb_v3QNS[i]->SetLineColor(2);
       lambda_PbPb_v3QNS[i]->SetLineColor(2);
       
     ks_pPb_v3QNS[i]->SetMarkerStyle(21);
     ks_PbPb_v3QNS[i]->SetMarkerStyle(21);
     lambda_pPb_v3QNS[i]->SetMarkerStyle(20);
     lambda_PbPb_v3QNS[i]->SetMarkerStyle(20);
     ks_pPb_v3QNS[i]->SetLineWidth(1);
     ks_PbPb_v3QNS[i]->SetLineWidth(1);
     lambda_pPb_v3QNS[i]->SetLineWidth(1);
     lambda_PbPb_v3QNS[i]->SetLineWidth(1);

     ks_pPb_v2QNS_ratio[i]->SetMarkerColor(4);
     ks_PbPb_v2QNS_ratio[i]->SetMarkerColor(4);
     lambda_pPb_v2QNS_ratio[i]->SetMarkerColor(2);
     lambda_PbPb_v2QNS_ratio[i]->SetMarkerColor(2);
       
       ks_pPb_v2QNS_ratio[i]->SetLineColor(4);
       ks_PbPb_v2QNS_ratio[i]->SetLineColor(4);
       lambda_pPb_v2QNS_ratio[i]->SetLineColor(2);
       lambda_PbPb_v2QNS_ratio[i]->SetLineColor(2);

     ks_pPb_v2QNS_ratio[i]->SetMarkerStyle(21);
     ks_PbPb_v2QNS_ratio[i]->SetMarkerStyle(21);
     lambda_pPb_v2QNS_ratio[i]->SetMarkerStyle(20);
     lambda_PbPb_v2QNS_ratio[i]->SetMarkerStyle(20);
     ks_pPb_v2QNS_ratio[i]->SetLineWidth(1);
     ks_PbPb_v2QNS_ratio[i]->SetLineWidth(1);
     lambda_pPb_v2QNS_ratio[i]->SetLineWidth(1);
     lambda_PbPb_v2QNS_ratio[i]->SetLineWidth(1);

     ks_pPb_v3QNS_ratio[i]->SetMarkerColor(4);
     ks_PbPb_v3QNS_ratio[i]->SetMarkerColor(4);
     lambda_pPb_v3QNS_ratio[i]->SetMarkerColor(2);
     lambda_PbPb_v3QNS_ratio[i]->SetMarkerColor(2);
       
       ks_pPb_v3QNS_ratio[i]->SetLineColor(4);
       ks_PbPb_v3QNS_ratio[i]->SetLineColor(4);
       lambda_pPb_v3QNS_ratio[i]->SetLineColor(2);
       lambda_PbPb_v3QNS_ratio[i]->SetLineColor(2);

     ks_pPb_v3QNS_ratio[i]->SetMarkerStyle(21);
     ks_PbPb_v3QNS_ratio[i]->SetMarkerStyle(21);
     lambda_pPb_v3QNS_ratio[i]->SetMarkerStyle(20);
     lambda_PbPb_v3QNS_ratio[i]->SetMarkerStyle(20);
     ks_pPb_v3QNS_ratio[i]->SetLineWidth(1);
     ks_PbPb_v3QNS_ratio[i]->SetLineWidth(1);
     lambda_pPb_v3QNS_ratio[i]->SetLineWidth(1);
     lambda_PbPb_v3QNS_ratio[i]->SetLineWidth(1);
   }
  
   TCanvas *c1 = new TCanvas("c1","c1",1,1,1200,720);
   c1->Range(0,0,1,1);
   TPad* pad1[12];
   pad1[0] = new TPad("pad10", "pad10",0.0,0.53125,0.3,1);
   pad1[1] = new TPad("pad11", "pad11",0.3,0.53125,0.533,1);
   pad1[2] = new TPad("pad12", "pad12",0.533,0.53125,0.766,1);
   pad1[3] = new TPad("pad13", "pad13",0.766,0.53125,1.0,1);
   pad1[4] = new TPad("pad14", "pad14",0.0,0.17625,0.3,0.53125);
   pad1[5] = new TPad("pad15", "pad15",0.3,0.17625,0.533,0.53125);
   pad1[6] = new TPad("pad16", "pad16",0.533,0.17625,0.766,0.53125);
   pad1[7] = new TPad("pad17", "pad17",0.766,0.17625,1.0,0.53125);
   pad1[8] = new TPad("pad18", "pad18",0.0,0.0,0.3,0.17625);
   pad1[9] = new TPad("pad19", "pad19",0.3,0.0,0.533,0.17625);
   pad1[10] = new TPad("pad110", "pad110",0.533,0.0,0.766,0.17625);
   pad1[11] = new TPad("pad111", "pad111",0.766,0.0,1.0,0.17625);

   for(int i=3; i>=0; i--)
   {
     pad1[i]->SetLeftMargin(0.0);
     pad1[i]->SetRightMargin(0);
     pad1[i]->SetTopMargin(0.0);
     pad1[i]->SetBottomMargin(0);
     pad1[i]->Draw();
   }
   for(int i=7; i>=4; i--)
   {
     pad1[i]->SetLeftMargin(0.0);
     pad1[i]->SetRightMargin(0);
     pad1[i]->SetTopMargin(0.0);
     pad1[i]->SetBottomMargin(0);
     pad1[i]->Draw();
   }
   for(int i=11; i>=8; i--)
   {
     pad1[i]->SetLeftMargin(0.0);
     pad1[i]->SetRightMargin(0);
     pad1[i]->SetTopMargin(0.0);
     pad1[i]->SetBottomMargin(0);
     pad1[i]->Draw();
   }

   pad1[0]->SetLeftMargin(0.2);
   pad1[4]->SetLeftMargin(0.2);
   pad1[8]->SetLeftMargin(0.2);
   pad1[3]->SetRightMargin(0.005);
   pad1[7]->SetRightMargin(0.005);
   pad1[11]->SetRightMargin(0.005);
   pad1[0]->SetTopMargin(0.02);
   pad1[1]->SetTopMargin(0.02);
   pad1[2]->SetTopMargin(0.02);
   pad1[3]->SetTopMargin(0.02);
   pad1[0]->SetBottomMargin(0.18);
   pad1[1]->SetBottomMargin(0.18);
   pad1[2]->SetBottomMargin(0.18);
   pad1[3]->SetBottomMargin(0.18);
   pad1[8]->SetBottomMargin(0.42);
   pad1[9]->SetBottomMargin(0.42);
   pad1[10]->SetBottomMargin(0.42);
   pad1[11]->SetBottomMargin(0.42);

   TCanvas *c2 = new TCanvas("c2","c2",1,1,1200,720);
   c2->Range(0,0,1,1);
   TPad* pad[12];

   pad[0] = new TPad("pad0", "pad0",0.0,0.53125,0.3,1);
   pad[1] = new TPad("pad1", "pad1",0.3,0.53125,0.533,1);
   pad[2] = new TPad("pad2", "pad2",0.533,0.53125,0.766,1);
   pad[3] = new TPad("pad3", "pad3",0.766,0.53125,1.0,1);
   pad[4] = new TPad("pad4", "pad4",0.0,0.17625,0.3,0.53125);
   pad[5] = new TPad("pad5", "pad5",0.3,0.17625,0.533,0.53125);
   pad[6] = new TPad("pad6", "pad6",0.533,0.17625,0.766,0.53125);
   pad[7] = new TPad("pad7", "pad7",0.766,0.17625,1.0,0.53125);
   pad[8] = new TPad("pad8", "pad8",0.0,0.0,0.3,0.17625);
   pad[9] = new TPad("pad9", "pad9",0.3,0.0,0.533,0.17625);
   pad[10] = new TPad("pad10", "pad10",0.533,0.0,0.766,0.17625);
   pad[11] = new TPad("pad11", "pad11",0.766,0.0,1.0,0.17625);
 
   for(int i=3; i>=0; i--)
   {
     pad[i]->SetLeftMargin(0.0);
     pad[i]->SetRightMargin(0);
     pad[i]->SetTopMargin(0.0);
     pad[i]->SetBottomMargin(0);
     pad[i]->Draw();
   }
   for(int i=7; i>=4; i--)
   {
     pad[i]->SetLeftMargin(0.0);
     pad[i]->SetRightMargin(0);
     pad[i]->SetTopMargin(0.0);
     pad[i]->SetBottomMargin(0);
     pad[i]->Draw();
   }
   for(int i=11; i>=8; i--)
   {
     pad[i]->SetLeftMargin(0.0);
     pad[i]->SetRightMargin(0);
     pad[i]->SetTopMargin(0.0);
     pad[i]->SetBottomMargin(0.0);
     pad[i]->Draw();
   }

   pad[0]->SetLeftMargin(0.2);
   pad[4]->SetLeftMargin(0.2);
   pad[8]->SetLeftMargin(0.2);  
   pad[3]->SetRightMargin(0.005);
   pad[7]->SetRightMargin(0.005);
   pad[11]->SetRightMargin(0.005);
   pad[0]->SetTopMargin(0.02);
   pad[1]->SetTopMargin(0.02);
   pad[2]->SetTopMargin(0.02);
   pad[3]->SetTopMargin(0.02);
   pad[0]->SetBottomMargin(0.18);
   pad[1]->SetBottomMargin(0.18);
   pad[2]->SetBottomMargin(0.18);
   pad[3]->SetBottomMargin(0.18);
   pad[8]->SetBottomMargin(0.42);
   pad[9]->SetBottomMargin(0.42);
   pad[10]->SetBottomMargin(0.42);
   pad[11]->SetBottomMargin(0.42);

   TCanvas *c3v2 = new TCanvas("c3v2","c3v2",1,1,370*1.2,720*1.2);
   c3v2->Range(0,0,1,1);
   TPad* pad3v2[3];

   pad3v2[0] = new TPad("pad3v20", "pad3v20",0.0,0.5,1.0,1.0);
   pad3v2[1] = new TPad("pad3v21", "pad3v21",0.0,0.15,1.0,0.5);
   pad3v2[2] = new TPad("pad3v22", "pad3v22",0.0,0.0,1.0,0.15);

   for(int i=0; i<3; i++)
   {
     pad3v2[i]->SetLeftMargin(0.0);
     pad3v2[i]->SetRightMargin(0);
     pad3v2[i]->SetTopMargin(0.0);
     pad3v2[i]->SetBottomMargin(0);
     pad3v2[i]->Draw();
   }
   pad3v2[0]->SetLeftMargin(0.16);
   pad3v2[1]->SetLeftMargin(0.16);
   pad3v2[2]->SetLeftMargin(0.16);
   pad3v2[0]->SetRightMargin(0.015);
   pad3v2[1]->SetRightMargin(0.015);
   pad3v2[2]->SetRightMargin(0.015);
   pad3v2[0]->SetTopMargin(0.015);
   pad3v2[1]->SetTopMargin(0.015);
   pad3v2[0]->SetBottomMargin(0.12);
   pad3v2[2]->SetBottomMargin(0.4);

   TCanvas *c3v3 = new TCanvas("c3v3","c3v3",1,1,370*1.2,720*1.2);
   c3v3->Range(0,0,1,1);
   TPad* pad3v3[3];

   pad3v3[0] = new TPad("pad3v30", "pad3v30",0.0,0.5,1.0,1.0);
   pad3v3[1] = new TPad("pad3v31", "pad3v31",0.0,0.15,1.0,0.5);
   pad3v3[2] = new TPad("pad3v32", "pad3v32",0.0,0.0,1.0,0.15);

   for(int i=0; i<3; i++)
   {
     pad3v3[i]->SetLeftMargin(0.0);
     pad3v3[i]->SetRightMargin(0);
     pad3v3[i]->SetTopMargin(0.0);
     pad3v3[i]->SetBottomMargin(0);
     pad3v3[i]->Draw();
   }
   pad3v3[0]->SetLeftMargin(0.16);
   pad3v3[1]->SetLeftMargin(0.16);
   pad3v3[2]->SetLeftMargin(0.16);
   pad3v3[0]->SetRightMargin(0.015);
   pad3v3[1]->SetRightMargin(0.015);
   pad3v3[2]->SetRightMargin(0.015);
   pad3v3[0]->SetTopMargin(0.015);
   pad3v3[1]->SetTopMargin(0.015);
   pad3v3[0]->SetBottomMargin(0.12);
   pad3v3[2]->SetBottomMargin(0.4);

   TCanvas *c4v3 = new TCanvas("c4v3","c4v3",1,1,740*1.2,340*1.2);
   c4v3->Range(0,0,1,1);
   TPad* pad4v3[3];

   pad4v3[0] = new TPad("pad4v30", "pad4v30",0.0,0.0,0.5,1);
   pad4v3[1] = new TPad("pad4v31", "pad4v31",0.5,0.3,1.0,1.0);
   pad4v3[2] = new TPad("pad4v32", "pad4v32",0.5,0.0,1.0,0.3);

   for(int i=0; i<3; i++)
   {
     pad4v3[i]->SetLeftMargin(0.0);
     pad4v3[i]->SetRightMargin(0);
     pad4v3[i]->SetTopMargin(0.0);
     pad4v3[i]->SetBottomMargin(0);
     pad4v3[i]->Draw();
   }
   pad4v3[0]->SetLeftMargin(0.16);
   pad4v3[1]->SetLeftMargin(0.2);
   pad4v3[2]->SetLeftMargin(0.2);
   pad4v3[0]->SetRightMargin(0.015);
   pad4v3[1]->SetRightMargin(0.005);
   pad4v3[2]->SetRightMargin(0.005);
   pad4v3[0]->SetTopMargin(0.015);
   pad4v3[1]->SetTopMargin(0.015);
   pad4v3[0]->SetBottomMargin(0.12);
   pad4v3[2]->SetBottomMargin(0.4);
/*
   TCanvas *c3 = new TCanvas("c3","c3",1,1,1200,700);
   c3->Range(0,0,1,1);
   TPad* pad3[12];

   pad3[0] = new TPad("pad30", "pad30",0.0,0.53125,0.3,1);
   pad3[1] = new TPad("pad31", "pad31",0.3,0.53125,0.533,1);
   pad3[2] = new TPad("pad32", "pad32",0.533,0.53125,0.766,1);
   pad3[3] = new TPad("pad33", "pad33",0.766,0.53125,1.0,1);
   pad3[4] = new TPad("pad34", "pad34",0.0,0.17625,0.3,0.53125);
   pad3[5] = new TPad("pad35", "pad35",0.3,0.17625,0.533,0.53125);
   pad3[6] = new TPad("pad36", "pad36",0.533,0.17625,0.766,0.53125);
   pad3[7] = new TPad("pad37", "pad37",0.766,0.17625,1.0,0.53125);
   pad3[8] = new TPad("pad38", "pad38",0.0,0.0,0.3,0.17625);
   pad3[9] = new TPad("pad39", "pad39",0.3,0.0,0.533,0.17625);
   pad3[10] = new TPad("pad310", "pad310",0.533,0.0,0.766,0.17625);
   pad3[11] = new TPad("pad311", "pad311",0.766,0.0,1.0,0.17625);

   for(int i=0; i<12; i++)
   {
     pad3[i]->SetLeftMargin(0.0);
     pad3[i]->SetRightMargin(0);
     pad3[i]->SetTopMargin(0.0);
     pad3[i]->SetBottomMargin(0);
     pad3[i]->Draw();
   }
   pad3[0]->SetLeftMargin(0.2);
   pad3[4]->SetLeftMargin(0.2);
   pad3[8]->SetLeftMargin(0.2);
   pad3[3]->SetRightMargin(0.005);
   pad3[7]->SetRightMargin(0.005);
   pad3[11]->SetRightMargin(0.005);
   pad3[0]->SetTopMargin(0.02);
   pad3[1]->SetTopMargin(0.02);
   pad3[2]->SetTopMargin(0.02);
   pad3[3]->SetTopMargin(0.02);
   pad3[0]->SetBottomMargin(0.18);
   pad3[1]->SetBottomMargin(0.18);
   pad3[2]->SetBottomMargin(0.18);
   pad3[3]->SetBottomMargin(0.18);
   pad3[8]->SetBottomMargin(0.4);
   pad3[9]->SetBottomMargin(0.4);
   pad3[10]->SetBottomMargin(0.4);
   pad3[11]->SetBottomMargin(0.4);

  TCanvas *c4 = new TCanvas("c4","c4",1,1,1160,580);
  makeMultiPanelCanvas(c4,4,2,0.01,0.01,0.15,0.15,0.01);
*/
  TCanvas *c4v2 = new TCanvas("c4v2","c4v2",1,1,910,580);
  makeMultiPanelCanvas(c4v2,3,2,0.01,0.01,0.16,0.16,0.01);
/*
  TCanvas *c1QNS = new TCanvas("c1QNS","c1QNS",1,1,880,580);
  makeMultiPanelCanvas(c1QNS,3,2,0.01,0.01,0.18,0.15,0.01);

  TCanvas *c2QNS = new TCanvas("c2QNS","c2QNS",1,1,1160,580);
  makeMultiPanelCanvas(c2QNS,4,2,0.01,0.01,0.18,0.15,0.01);

  TCanvas *c3QNS = new TCanvas("c3QNS","c3QNS",1,1,880,580);
  makeMultiPanelCanvas(c3QNS,3,2,0.01,0.01,0.18,0.15,0.01);

  TCanvas *c4QNS = new TCanvas("c4QNS","c4QNS",1,1,1160,580);
  makeMultiPanelCanvas(c4QNS,4,2,0.01,0.01,0.18,0.15,0.01);
*/
  TH1D* hist12 = new TH1D("hist12","",100,-0.3,5.9);
  hist12->SetLineWidth(0);
  hist12->SetXTitle("p_{T} (GeV)");
  hist12->SetYTitle("v_{2}{2, |#Delta#eta|>2}");
  hist12->SetMinimum(-0.015);
  hist12->SetMaximum(0.315);
  fixedFontHist1D(hist12,1.5,2.4);
  hist12->GetXaxis()->CenterTitle(1);
  hist12->GetYaxis()->CenterTitle(1);
    hist12->GetXaxis()->SetTitleSize(hist12->GetXaxis()->GetTitleSize()*1.3);
  hist12->GetYaxis()->SetTitleSize(hist12->GetYaxis()->GetTitleSize()*1.3);
    hist12->GetXaxis()->SetLabelSize(hist12->GetXaxis()->GetLabelSize()*1.2);
    hist12->GetYaxis()->SetLabelSize(hist12->GetYaxis()->GetLabelSize()*1.2);
    hist12->GetXaxis()->SetTitleOffset(hist12->GetXaxis()->GetTitleOffset()*1.1);
  hist12->GetYaxis()->SetTitleOffset(hist12->GetYaxis()->GetTitleOffset()*1.1);

  TH1D* hist = new TH1D("hist","",100,-0.12,5.9);
  hist->SetLineWidth(0);
    hist->SetLineStyle(7);
  hist->SetXTitle("p_{T} (GeV)");
  hist->SetYTitle("v_{2}");
  hist->SetMinimum(-0.015);
  hist->SetMaximum(0.315);
  fixedFontHist1D(hist,1.65,2.1);
  hist->GetXaxis()->CenterTitle(1);
  hist->GetYaxis()->CenterTitle(1);
    hist->GetXaxis()->SetTitleSize(hist->GetXaxis()->GetTitleSize()*1.3);
  hist->GetYaxis()->SetTitleSize(hist->GetYaxis()->GetTitleSize()*1.3);
    hist->GetXaxis()->SetLabelSize(hist->GetXaxis()->GetLabelSize()*1.2);
    hist->GetYaxis()->SetLabelSize(hist->GetYaxis()->GetLabelSize()*1.2);

  TH1D* hist_PbPb = new TH1D("hist_PbPb","",100,-0.12,5.9);
  hist_PbPb->SetLineWidth(0);
    hist_PbPb->SetLineStyle(7);
  hist_PbPb->SetXTitle("p_{T} (GeV)");
  hist_PbPb->SetYTitle("v_{2}");
  hist_PbPb->SetMinimum(-0.015);
  hist_PbPb->SetMaximum(0.315);
  fixedFontHist1D(hist_PbPb,1.65,2.1);
  hist_PbPb->GetXaxis()->CenterTitle(1);
  hist_PbPb->GetYaxis()->CenterTitle(1);
    hist_PbPb->GetXaxis()->SetTitleSize(hist_PbPb->GetXaxis()->GetTitleSize()*1.3);
  hist_PbPb->GetYaxis()->SetTitleSize(hist_PbPb->GetYaxis()->GetTitleSize()*1.3);
    hist_PbPb->GetXaxis()->SetLabelSize(hist_PbPb->GetXaxis()->GetLabelSize()*1.2);
    hist_PbPb->GetYaxis()->SetLabelSize(hist_PbPb->GetYaxis()->GetLabelSize()*1.2);

  TH1D* histlow = new TH1D("histlow","",100,-0.12,5.9);
  histlow->SetLineWidth(0);
    histlow->SetLineStyle(7);
  histlow->SetXTitle("p_{T} (GeV)");
  histlow->SetYTitle("v_{2}");
  histlow->SetMinimum(-0.04);
  histlow->SetMaximum(0.599);
  fixedFontHist1D(histlow,1.5,2.1);
  histlow->GetXaxis()->CenterTitle(1);
  histlow->GetYaxis()->CenterTitle(1);
    histlow->GetXaxis()->SetTitleSize(histlow->GetXaxis()->GetTitleSize()*1.3);
  histlow->GetYaxis()->SetTitleSize(histlow->GetYaxis()->GetTitleSize()*1.3);
  histlow->GetYaxis()->SetTitleOffset(histlow->GetYaxis()->GetTitleOffset()*0.8);
    histlow->GetXaxis()->SetLabelSize(histlow->GetXaxis()->GetLabelSize()*1.2);
    histlow->GetYaxis()->SetLabelSize(histlow->GetYaxis()->GetLabelSize()*1.2);


  TH1D* hist1 = new TH1D("hist1","",100,-0.12,5.9);
  hist1->SetLineWidth(0);
  hist1->SetXTitle("p_{T} (GeV)");
  hist1->SetYTitle("v_{2}");
  hist1->SetMinimum(-0.045);
  hist1->SetMaximum(0.48);
  fixedFontHist1D(hist1,1.8,1.5);
  hist1->GetXaxis()->CenterTitle(1);
  hist1->GetYaxis()->CenterTitle(1);
    hist1->GetXaxis()->SetTitleSize(hist1->GetXaxis()->GetTitleSize()*1.3);
  hist1->GetYaxis()->SetTitleSize(hist1->GetYaxis()->GetTitleSize()*1.3);
    hist1->GetXaxis()->SetLabelSize(hist1->GetXaxis()->GetLabelSize()*1.2);
    hist1->GetYaxis()->SetLabelSize(hist1->GetYaxis()->GetLabelSize()*1.2);

  TH1D* hist2 = new TH1D("hist2","",100,-0.12,5.9);
  hist2->SetLineWidth(0);
    hist2->SetLineStyle(9);
  hist2->SetXTitle("p_{T} (GeV)");
  hist2->SetYTitle("v_{3}");
  hist2->SetMinimum(-0.005);
  hist2->SetMaximum(0.149);
  fixedFontHist1D(hist2,1.5,2.4);
  hist2->GetXaxis()->CenterTitle(1);
  hist2->GetYaxis()->CenterTitle(1);
    hist2->GetXaxis()->SetTitleSize(hist2->GetXaxis()->GetTitleSize()*1.3);
  hist2->GetYaxis()->SetTitleSize(hist2->GetYaxis()->GetTitleSize()*1.3);
    hist2->GetXaxis()->SetLabelSize(hist2->GetXaxis()->GetLabelSize()*1.2);
    hist2->GetYaxis()->SetLabelSize(hist2->GetYaxis()->GetLabelSize()*1.2);


  TH1D* hist3 = new TH1D("hist3","",100,-0.12,5.9);
  hist3->SetLineWidth(0);
  hist3->SetXTitle("p_{T} (GeV)");
  hist3->SetYTitle("v_{3}");
  hist3->SetMinimum(-0.005);
  hist3->SetMaximum(0.099);
  fixedFontHist1D(hist3,0.8,1.);
  hist3->GetXaxis()->CenterTitle(1);
  hist3->GetYaxis()->CenterTitle(1);
    hist3->GetXaxis()->SetTitleSize(hist3->GetXaxis()->GetTitleSize()*1.3);
  hist3->GetYaxis()->SetTitleSize(hist3->GetYaxis()->GetTitleSize()*1.3);
    hist3->GetXaxis()->SetLabelSize(hist3->GetXaxis()->GetLabelSize()*1.2);
    hist3->GetYaxis()->SetLabelSize(hist3->GetYaxis()->GetLabelSize()*1.2);

  TH1D* histQNS12 = new TH1D("histQNS12","",100,-0.12,2.12);
  histQNS12->SetLineWidth(0);
  histQNS12->SetXTitle("KE_{T}/n_{q} (GeV)");
  histQNS12->SetYTitle("v_{2}/n_{q}");
  histQNS12->SetMinimum(-0.008);
  histQNS12->SetMaximum(0.129);
  fixedFontHist1D(histQNS12,1.8,2.4);
  histQNS12->GetXaxis()->CenterTitle(1);
  histQNS12->GetYaxis()->CenterTitle(1);
    histQNS12->GetXaxis()->SetTitleSize(histQNS12->GetXaxis()->GetTitleSize()*1.3);
  histQNS12->GetYaxis()->SetTitleSize(histQNS12->GetYaxis()->GetTitleSize()*1.3);
    histQNS12->GetXaxis()->SetLabelSize(histQNS12->GetXaxis()->GetLabelSize()*1.2);
    histQNS12->GetYaxis()->SetLabelSize(histQNS12->GetYaxis()->GetLabelSize()*1.2);
    histQNS12->GetXaxis()->SetTitleOffset(histQNS12->GetXaxis()->GetTitleOffset()*1.1);
  histQNS12->GetYaxis()->SetTitleOffset(histQNS12->GetYaxis()->GetTitleOffset()*1.1);

  TH1D* histQNS = new TH1D("histQNS","",100,-0.12,2.12);
  histQNS->SetLineWidth(0);
    histQNS->SetLineStyle(7);
  histQNS->SetXTitle("KE_{T}/n_{q} (GeV)");
  histQNS->SetYTitle("v_{2}/n_{q}");
  histQNS->SetMinimum(-0.008);
  histQNS->SetMaximum(0.129);
  fixedFontHist1D(histQNS,1.9,2.3);
  histQNS->GetXaxis()->CenterTitle(1);
  histQNS->GetYaxis()->CenterTitle(1);
    histQNS->GetXaxis()->SetTitleSize(histQNS->GetXaxis()->GetTitleSize()*1.3);
  histQNS->GetYaxis()->SetTitleSize(histQNS->GetYaxis()->GetTitleSize()*1.3);
    histQNS->GetXaxis()->SetLabelSize(histQNS->GetXaxis()->GetLabelSize()*1.2);
    histQNS->GetYaxis()->SetLabelSize(histQNS->GetYaxis()->GetLabelSize()*1.2);

  TH1D* histQNS_PbPb = new TH1D("histQNS_PbPb","",100,-0.12,2.12);
  histQNS_PbPb->SetLineWidth(0);
    histQNS_PbPb->SetLineStyle(7);
  histQNS_PbPb->SetXTitle("KE_{T}/n_{q} (GeV)");
  histQNS_PbPb->SetYTitle("v_{2}/n_{q}");
  histQNS_PbPb->SetMinimum(-0.008);
  histQNS_PbPb->SetMaximum(0.129);
  fixedFontHist1D(histQNS_PbPb,1.9,2.3);
  histQNS_PbPb->GetXaxis()->CenterTitle(1);
  histQNS_PbPb->GetYaxis()->CenterTitle(1);
    histQNS_PbPb->GetXaxis()->SetTitleSize(histQNS_PbPb->GetXaxis()->GetTitleSize()*1.3);
  histQNS_PbPb->GetYaxis()->SetTitleSize(histQNS_PbPb->GetYaxis()->GetTitleSize()*1.3);
    histQNS_PbPb->GetXaxis()->SetLabelSize(histQNS_PbPb->GetXaxis()->GetLabelSize()*1.2);
    histQNS_PbPb->GetYaxis()->SetLabelSize(histQNS_PbPb->GetYaxis()->GetLabelSize()*1.2);

  TH1D* hist1QNS = new TH1D("hist1QNS","",100,-0.12,2.12);
  hist1QNS->SetLineWidth(0);
  hist1QNS->SetXTitle("KE_{T}/n_{q} (GeV)");
  hist1QNS->SetYTitle("v_{2}/n_{q}");
  hist1QNS->SetMinimum(-0.008);
  hist1QNS->SetMaximum(0.22);
  fixedFontHist1D(hist1QNS,1.8,2.2);
  hist1QNS->GetXaxis()->CenterTitle(1);
  hist1QNS->GetYaxis()->CenterTitle(1);
    hist1QNS->GetXaxis()->SetTitleSize(hist1QNS->GetXaxis()->GetTitleSize()*1.3);
  hist1QNS->GetYaxis()->SetTitleSize(hist1QNS->GetYaxis()->GetTitleSize()*1.3);
    hist1QNS->GetXaxis()->SetLabelSize(hist1QNS->GetXaxis()->GetLabelSize()*1.2);
    hist1QNS->GetYaxis()->SetLabelSize(hist1QNS->GetYaxis()->GetLabelSize()*1.2);

  TH1D* hist2QNS = new TH1D("hist2QNS","",100,-0.12,2.12);
  hist2QNS->SetLineWidth(0);
    hist2QNS->SetLineStyle(9);
  hist2QNS->SetXTitle("KE_{T}/n_{q} (GeV)");
  hist2QNS->SetYTitle("v_{3}/n_{q}");
  hist2QNS->SetMinimum(-0.003);
  hist2QNS->SetMaximum(0.049);
  fixedFontHist1D(hist2QNS,1.8,2.4);
  hist2QNS->GetXaxis()->CenterTitle(1);
  hist2QNS->GetYaxis()->CenterTitle(1);
    hist2QNS->GetXaxis()->SetTitleSize(hist2QNS->GetXaxis()->GetTitleSize()*1.3);
  hist2QNS->GetYaxis()->SetTitleSize(hist2QNS->GetYaxis()->GetTitleSize()*1.3);
    hist2QNS->GetXaxis()->SetLabelSize(hist2QNS->GetXaxis()->GetLabelSize()*1.2);
    hist2QNS->GetYaxis()->SetLabelSize(hist2QNS->GetYaxis()->GetLabelSize()*1.2);

  TH1D* hist3QNS = new TH1D("hist3QNS","",100,-0.12,2.12);
  hist3QNS->SetLineWidth(0);
  hist3QNS->SetXTitle("KE_{T}/n_{q} (GeV)");
  hist3QNS->SetYTitle("v_{3}/n_{q}");
  hist3QNS->SetMinimum(-0.003);
  hist3QNS->SetMaximum(0.22);
  fixedFontHist1D(hist3QNS,0.9,1.0);
  hist3QNS->GetXaxis()->CenterTitle(1);
  hist3QNS->GetYaxis()->CenterTitle(1);
    hist3QNS->GetXaxis()->SetTitleSize(hist3QNS->GetXaxis()->GetTitleSize()*1.3);
  hist3QNS->GetYaxis()->SetTitleSize(hist3QNS->GetYaxis()->GetTitleSize()*1.3);
    hist3QNS->GetXaxis()->SetLabelSize(hist3QNS->GetXaxis()->GetLabelSize()*1.2);
    hist3QNS->GetYaxis()->SetLabelSize(hist3QNS->GetYaxis()->GetLabelSize()*1.2);

  TH1D* histratio = new TH1D("histratio","",100,-0.12,2.12);
  histratio->SetLineWidth(0);
  histratio->SetXTitle("KE_{T}/n_{q} (GeV)");
  histratio->SetYTitle("Data/Fit");
  histratio->SetMinimum(0.51);
  histratio->SetMaximum(1.49);
  fixedFontHist1D(histratio,4.7,2.65);
  histratio->GetXaxis()->CenterTitle(1);
  histratio->GetYaxis()->CenterTitle(1);
    histratio->GetXaxis()->SetTitleSize(histratio->GetXaxis()->GetTitleSize()*1.3);
  histratio->GetYaxis()->SetTitleSize(histratio->GetYaxis()->GetTitleSize()*1.1);
    histratio->GetXaxis()->SetLabelSize(histratio->GetXaxis()->GetLabelSize()*1.2);
    histratio->GetYaxis()->SetLabelSize(histratio->GetYaxis()->GetLabelSize()*1.2);
  TLine* lline = new TLine(-0.1,1.0,2.1,1.0);
  lline->SetLineStyle(2);

  TH1D* histratio_v3 = new TH1D("histratio_v3","",100,-0.12,2.12);
  histratio_v3->SetLineWidth(0);
  histratio_v3->SetXTitle("KE_{T}/n_{q} (GeV)");
  histratio_v3->SetYTitle("Data/Fit");
  histratio_v3->SetMinimum(0.46);
  histratio_v3->SetMaximum(1.59);
  fixedFontHist1D(histratio_v3,5.1,2.8);
  histratio_v3->GetXaxis()->CenterTitle(1);
  histratio_v3->GetYaxis()->CenterTitle(1);
    histratio_v3->GetXaxis()->SetTitleSize(histratio_v3->GetXaxis()->GetTitleSize()*1.3);
  histratio_v3->GetYaxis()->SetTitleSize(histratio_v3->GetYaxis()->GetTitleSize()*1.1);
    histratio_v3->GetXaxis()->SetLabelSize(histratio_v3->GetXaxis()->GetLabelSize()*1.2);
    histratio_v3->GetYaxis()->SetLabelSize(histratio_v3->GetYaxis()->GetLabelSize()*1.2);

  TH1D* histratio_v2single = new TH1D("histratio_v2single","",100,-0.12,2.12);
  histratio_v2single->SetLineWidth(0);
  histratio_v2single->SetXTitle("KE_{T}/n_{q} (GeV)");
  histratio_v2single->SetYTitle("Data/Fit");
  histratio_v2single->SetMinimum(0.46);
  histratio_v2single->SetMaximum(1.59);
  fixedFontHist1D(histratio_v2single,5.1,2.8);
  histratio_v2single->GetXaxis()->CenterTitle(1);
  histratio_v2single->GetYaxis()->CenterTitle(1);
    histratio_v2single->GetXaxis()->SetTitleSize(histratio_v2single->GetXaxis()->GetTitleSize()*1.3);
  histratio_v2single->GetYaxis()->SetTitleSize(histratio_v2single->GetYaxis()->GetTitleSize()*1.1);
    histratio_v2single->GetXaxis()->SetLabelSize(histratio_v2single->GetXaxis()->GetLabelSize()*1.2);
    histratio_v2single->GetYaxis()->SetLabelSize(histratio_v2single->GetYaxis()->GetLabelSize()*1.2);

  TLine* lline = new TLine(-0.1,1.0,2.1,1.0);
  lline->SetLineStyle(2);
/*
  for(int i=0;i<3;i++)
  {
    c1QNS->cd(i+1);
    hist1QNS->Draw();
    ks_PbPb_v2QNS[i]->Draw("PESAME");
    lambda_PbPb_v2QNS[i]->Draw("PESAME");
    c1QNS->cd(i+1+3);
    hist1QNS->Draw();
    ks_pPb_v2QNS[i]->Draw("PESAME");
    lambda_pPb_v2QNS[i]->Draw("PESAME");

    c3QNS->cd(i+1);
    hist3QNS->Draw();
    ks_PbPb_v3QNS[i]->Draw("PESAME");
    lambda_PbPb_v3QNS[i]->Draw("PESAME");
    c3QNS->cd(i+1+3);
    hist3QNS->Draw();
    ks_pPb_v3QNS[i]->Draw("PESAME");
    lambda_pPb_v3QNS[i]->Draw("PESAME");
  }
*/

  double syst_pPb = 0.069;
  double syst_PbPb = 0.066;
  double syst_pPbH = 0.039;
  double syst_PbPbH = 0.03;

//begin add
  TBox* box;
//end add
  for(int i=0;i<4;i++)
  {
    pad[i]->cd();
    hist->Draw();
//begin add
    int n = ks_pPb_v2[i+3]->GetN();
    for(int j=0; j<n; j++)
    {
      double xk,yk,xl,yl,xh,yh;
      if(lambda_pPb_v2[i+3]->GetN()-1 >= j){
        lambda_pPb_v2[i+3]->GetPoint(j,xl,yl);
        double percent1L = syst_pPb;
      }
      if(pPb_v2[i+3]->GetN()-1 >= j)
      {
        pPb_v2[i+3]->GetPoint(j,xh,yh);
        double percentH = syst_pPbH;
      }
      ks_pPb_v2[i+3]->GetPoint(j,xk,yk);
      double percent1K = syst_pPb;
      double yLerr, yKerr, yHerr;

      yKerr=fabs(yk)*percent1K;
      yLerr=fabs(yl)*percent1L;
      yHerr=fabs(yh)*percentH;

      if(lambda_pPb_v2[i+3]->GetN()-1 >= j)
      {
        box = new TBox(xl-.15,yl-yLerr,xl+.15,yl+yLerr);
        box->SetFillColor(17);
        box->SetFillStyle(1001);
        box->SetLineWidth(0);
        box->Draw("SAME");
      }
      if(pPb_v2[i+3]->GetN()-1 >= j){
        box = new TBox(xh-.15,yh-yHerr,xh+.15,yh+yHerr);
        box->SetFillColor(17);
        box->SetFillStyle(1001);
        box->SetLineWidth(0);
        box->Draw("SAME");
      }
      box = new TBox(xk-.15,yk-yKerr,xk+.15,yk+yKerr);
      box->SetFillColor(17);
      box->SetFillStyle(1001);
      box->SetLineWidth(0);
      box->Draw("SAME");
    }
//end add

    ks_pPb_v2[i+3]->Draw("PESAME");
    lambda_pPb_v2[i+3]->Draw("PESAME");
    pPb_v2[i+3]->Draw("PESAME");
    
    
    pad[i+4]->cd();
    histQNS->Draw();
//begin add
    int n = ks_pPb_v2QNS[i+3]->GetN()-1;

    for(int j=0; j<n; j++)
    {
      double xk,yk,xl,yl;
      if(lambda_pPb_v2QNS[i+3]->GetN()-1 >= j){
        lambda_pPb_v2QNS[i+3]->GetPoint(j,xl,yl);
        double percent1L = syst_pPb;
      }
      ks_pPb_v2QNS[i+3]->GetPoint(j,xk,yk);
      double percent1K = syst_pPb;
      double yLerr, yKerr;

      yKerr=fabs(yk)*percent1K;
      yLerr=fabs(yl)*percent1L;

      if(lambda_pPb_v2QNS[i+3]->GetN()-1 >= j){
        box = new TBox(xl-.049,yl-yLerr,xl+.049,yl+yLerr);
        box->SetFillColor(17);
        box->SetFillStyle(1001);
        box->SetLineWidth(0);
        box->Draw("SAME");
      }
		
      box = new TBox(xk-.049,yk-yKerr,xk+.049,yk+yKerr);
      box->SetFillColor(17);
      box->SetFillStyle(1001);
      box->SetLineWidth(0);
      box->Draw("SAME");
    }
//end add

    ks_pPb_v2QNS[i+3]->Draw("PESAME");
    lambda_pPb_v2QNS[i+3]->Draw("PESAME");
    fitfunc_pPb_v2QNS[i+3]->Draw("LSAME");


    pad[i+8]->cd();
    histratio->Draw();
    lline->Draw("same");
//begin add
    int n = ks_pPb_v2QNS_ratio[i+3]->GetN()-1;

    for(int j=0; j<n; j++)
    {
      double xk,yk,xl,yl;
      if(lambda_pPb_v2QNS_ratio[i+3]->GetN()-1 >= j){
        lambda_pPb_v2QNS_ratio[i+3]->GetPoint(j,xl,yl);
        double percent1L = syst_pPb;
      }
      ks_pPb_v2QNS_ratio[i+3]->GetPoint(j,xk,yk);
      double percent1K = syst_pPb;
      double yLerr, yKerr;

      yKerr=fabs(yk)*percent1K;
      yLerr=fabs(yl)*percent1L;
		
      if(lambda_pPb_v2QNS_ratio[i+3]->GetN()-1 >= j){
        box = new TBox(xl-.049,yl-yLerr,xl+.049,yl+yLerr);
        box->SetFillColor(17);
        box->SetFillStyle(1001);
        box->SetLineWidth(0);
        box->Draw("SAME");
      }
		
      box = new TBox(xk-.049,yk-yKerr,xk+.049,yk+yKerr);
      box->SetFillColor(17);
      box->SetFillStyle(1001);
      box->SetLineWidth(0);
      box->Draw("SAME");
    }
//end add    
    ks_pPb_v2QNS_ratio[i+3]->Draw("PESAME");
    lambda_pPb_v2QNS_ratio[i+3]->Draw("PESAME");


    pad1[i]->cd();
    hist_PbPb->Draw();
//begin add
    int n = ks_PbPb_v2[i+3]->GetN();

    for(int j=0; j<n; j++)
    {
      double xk,yk,xl,yl,xh,yh;
      if(lambda_PbPb_v2[i+3]->GetN()-1 >= j){
        lambda_PbPb_v2[i+3]->GetPoint(j,xl,yl);
        double percent1L = syst_PbPb;
      }
      if(PbPb_v2[i+3]->GetN()-1 >= j)
      {
        PbPb_v2[i+3]->GetPoint(j,xh,yh);
        double percentH = syst_PbPbH;
      } 
      ks_PbPb_v2[i+3]->GetPoint(j,xk,yk);
      double percent1K = syst_PbPb;
      double yLerr, yKerr, yHerr;

      yKerr=fabs(yk)*percent1K;
      yLerr=fabs(yl)*percent1L;

      if(PbPb_v2[i+3]->GetN()-1 >= j){
        yHerr=fabs(yh)*percentH;
      }

      if(lambda_PbPb_v2[i+3]->GetN()-1 >= j){
        box = new TBox(xl-.15,yl-yLerr,xl+.15,yl+yLerr);
        box->SetFillColor(17);
        box->SetFillStyle(1001);
        box->SetLineWidth(0);
        box->Draw("Fsame");
      }
      if(PbPb_v2[i+3]->GetN()-1 >= j){
        box = new TBox(xh-.15,yh-yHerr,xh+.15,yh+yHerr);
        box->SetFillColor(17);
        box->SetFillStyle(1001);
        box->SetLineWidth(0);
        box->Draw("SAME");
      }

      box = new TBox(xk-.15,yk-yKerr,xk+.15,yk+yKerr);
      box->SetFillColor(17);
      box->SetFillStyle(1001);
      box->SetLineWidth(0);
      box->Draw("Fsame");

    }
//end add
    ks_PbPb_v2[i+3]->Draw("PESAME");
    lambda_PbPb_v2[i+3]->Draw("PESAME");
    PbPb_v2[i+3]->Draw("PESAME");

    
    pad1[i+4]->cd();
    histQNS_PbPb->Draw();
//begin add
    int n = ks_PbPb_v2QNS[i+3]->GetN()-1;
    
    for(int j=0; j<n; j++)
    {
      double xk,yk,xl,yl;
      if(lambda_PbPb_v2QNS[i+3]->GetN()-1 >= j){
        lambda_PbPb_v2QNS[i+3]->GetPoint(j,xl,yl);
        double percent1L = syst_PbPb;
      }
      ks_PbPb_v2QNS[i+3]->GetPoint(j,xk,yk);
      double percent1K = syst_PbPb;
      double yLerr, yKerr;

      yKerr=yk*percent1K;
      yLerr=yl*percent1L;

      if(lambda_PbPb_v2QNS[i+3]->GetN()-1 >= j){
        box = new TBox(xl-.049,yl-yLerr,xl+.049,yl+yLerr);
        box->SetFillColor(17);
        box->SetFillStyle(1001);
        box->SetLineWidth(0);
        box->Draw("SAME");
      }

      box = new TBox(xk-.049,yk-yKerr,xk+.049,yk+yKerr);
      box->SetFillColor(17);
      box->SetFillStyle(1001);
      box->SetLineWidth(0);
      box->Draw("SAME");
    }
//end add

    ks_PbPb_v2QNS[i+3]->Draw("PESAME");
    lambda_PbPb_v2QNS[i+3]->Draw("PESAME");
    fitfunc_PbPb_v2QNS[i+3]->Draw("LSAME");
    pad1[i+8]->cd();
    histratio->Draw();
    lline->Draw("same");
//begin add
    int n = ks_PbPb_v2QNS_ratio[i+3]->GetN()-1;

    for(int j=0; j<n; j++)
    {
      double xk,yk,xl,yl;
      if(lambda_PbPb_v2QNS_ratio[i+3]->GetN()-1 >= j){
        lambda_PbPb_v2QNS_ratio[i+3]->GetPoint(j,xl,yl);
        double percent1L = syst_PbPb;
      }
      ks_PbPb_v2QNS_ratio[i+3]->GetPoint(j,xk,yk);
      double percent1K = syst_PbPb;
      double yLerr, yKerr;

      yKerr=fabs(yk)*percent1K;
      yLerr=fabs(yl)*percent1L;

      if(lambda_PbPb_v2QNS_ratio[i+3]->GetN()-1 >= j){
        box = new TBox(xl-.049,yl-yLerr,xl+.049,yl+yLerr);
        box->SetFillColor(17);
        box->SetFillStyle(1001);
        box->SetLineWidth(0);
        box->Draw("SAME");
      }
		
      box = new TBox(xk-.049,yk-yKerr,xk+.049,yk+yKerr);
      box->SetFillColor(17);
      box->SetFillStyle(1001);
      box->SetLineWidth(0);
      box->Draw("SAME");
    }
//end add
    ks_PbPb_v2QNS_ratio[i+3]->Draw("PESAME");
    lambda_PbPb_v2QNS_ratio[i+3]->Draw("PESAME");

/*
    pad3[i]->cd();
    hist2->Draw();
    ks_pPb_v3[i+3]->Draw("PESAME");
    lambda_pPb_v3[i+3]->Draw("PESAME");
    pPb_v3[i+3]->Draw("PESAME");
    pad3[i+4]->cd();
    hist2QNS->Draw();
    ks_pPb_v3QNS[i+3]->Draw("PESAME");
    lambda_pPb_v3QNS[i+3]->Draw("PESAME");
    fitfunc_pPb_v3QNS[i+3]->Draw("LSAME");
    pad3[i+8]->cd();
    histratio->Draw();
    lline->Draw("same");
    ks_pPb_v3QNS_ratio[i+3]->Draw("PESAME");
    lambda_pPb_v3QNS_ratio[i+3]->Draw("PESAME");

    c4->cd(i+1);
    hist2->Draw();
    ks_PbPb_v3[i+3]->Draw("PESAME");
    lambda_PbPb_v3[i+3]->Draw("PESAME");
    c4->cd(i+1+4);
    hist2->Draw();
    ks_pPb_v3[i+3]->Draw("PESAME");
    lambda_pPb_v3[i+3]->Draw("PESAME");

    c2QNS->cd(i+1);
    histQNS->Draw();
    ks_PbPb_v2QNS[i+3]->Draw("PESAME");
    lambda_PbPb_v2QNS[i+3]->Draw("PESAME");
    c2QNS->cd(i+1+4);
    histQNS->Draw();
    ks_pPb_v2QNS[i+3]->Draw("PESAME");
    lambda_pPb_v2QNS[i+3]->Draw("PESAME");

    c4QNS->cd(i+1);
    hist2QNS->Draw();
    ks_PbPb_v3QNS[i+3]->Draw("PESAME");
    lambda_PbPb_v3QNS[i+3]->Draw("PESAME");
    c4QNS->cd(i+1+4);
    hist2QNS->Draw();
    ks_pPb_v3QNS[i+3]->Draw("PESAME");
    lambda_pPb_v3QNS[i+3]->Draw("PESAME");
*/
  }

  for(int i=2;i>=0;i--)
  {
    c4v2->cd(i+1);
    c4v2->GetPad(i+1)->cd();
    histlow->Draw();
//begin add
    int n = ks_PbPb_v2[i]->GetN();

	  for(int j=0; j<n; j++)
    {
		  double xk,yk,xl,yl,xh,yh;
		  if(lambda_PbPb_v2[i]->GetN()-1 >= j){
			  lambda_PbPb_v2[i]->GetPoint(j,xl,yl);
			  double percent1L = syst_PbPb;
		  }
		  ks_PbPb_v2[i]->GetPoint(j,xk,yk);
		  double percent1K = syst_PbPb;
		  double yLerr, yKerr;
        
        PbPb_v2[i]->GetPoint(j,xh,yh);
        double percent1H = syst_PbPbH;
        double yHerr;
		
		  yKerr=fabs(yk)*percent1K;
		  yLerr=fabs(yl)*percent1L;
        yHerr=fabs(yh)*percent1H;

		  if(lambda_PbPb_v2[i]->GetN()-1 >= j){
			  box = new TBox(xl-.15,yl-yLerr,xl+.15,yl+yLerr);
			  box->SetFillColor(17);
			  box->SetFillStyle(1001);
			  box->SetLineWidth(0);
			  box->Draw("SAME");
		  }
		
		  box = new TBox(xk-.15,yk-yKerr,xk+.15,yk+yKerr);
		  box->SetFillColor(17);
		  box->SetFillStyle(1001);
		  box->SetLineWidth(0);
		  box->Draw("SAME");

        box = new TBox(xh-.15,yh-yHerr,xh+.15,yh+yHerr);
        box->SetFillColor(17);
        box->SetFillStyle(1001);
        box->SetLineWidth(0);
        box->Draw("SAME");
    }
//end add
    ks_PbPb_v2[i]->Draw("PESAME");
    lambda_PbPb_v2[i]->Draw("PESAME");
      PbPb_v2[i]->Draw("PESAME");
    c4v2->cd(i+1+3);
    c4v2->GetPad(i+1+3)->Draw();
    histlow->Draw();
//begin add
    int n = ks_pPb_v2[i]->GetN();

	  for(int j=0; j<n; j++)
    {
		  double xk,yk,xl,yl,xh,yh;
		  if(lambda_pPb_v2[i]->GetN()-1 >= j){
			  lambda_pPb_v2[i]->GetPoint(j,xl,yl);
			  double percent1L = syst_pPb;
		  }
		  ks_pPb_v2[i]->GetPoint(j,xk,yk);
		  double percent1K = syst_pPb;
		  double yLerr, yKerr;
        
        pPb_v2[i]->GetPoint(j,xh,yh);
        double percent1H = syst_pPbH;
        double yHerr;
		
			  yKerr=fabs(yk)*percent1K;
			  yLerr=fabs(yl)*percent1L;
        yHerr=fabs(yh)*percent1H;

		  if(lambda_pPb_v2[i]->GetN()-1 >= j){
			  box = new TBox(xl-.15,yl-yLerr,xl+.15,yl+yLerr);
			  box->SetFillColor(17);
			  box->SetFillStyle(1001);
			  box->SetLineWidth(0);
			  box->Draw("SAME");
		  }
		
		  box = new TBox(xk-.15,yk-yKerr,xk+.15,yk+yKerr);
		  box->SetFillColor(17);
		  box->SetFillStyle(1001);
		  box->SetLineWidth(0);
		  box->Draw("SAME");
        
        box = new TBox(xh-.15,yh-yHerr,xh+.15,yh+yHerr);
        box->SetFillColor(17);
        box->SetFillStyle(1001);
        box->SetLineWidth(0);
        box->Draw("SAME");
	  }
//end add
    ks_pPb_v2[i]->Draw("PESAME");
    lambda_pPb_v2[i]->Draw("PESAME");
      pPb_v2[i]->Draw("PESAME");
  }

  pad3v3[0]->cd();
  hist2->Draw();
//begin add
  int n = ks_pPb_v3[7]->GetN();
    int nh = pPb_v3[5]->GetN();
    
    for(int j=0; j<nh; j++)
    {
        double xh,yh;
        pPb_v3[5]->GetPoint(j,xh,yh);
        double percentH = syst_pPbH;
        
        double yHerr;
        yHerr=fabs(yh)*percentH;
        
        box = new TBox(xh-.15,yh-yHerr,xh+.15,yh+yHerr);
        box->SetFillColor(17);
        box->SetFillStyle(1001);
        box->SetLineWidth(0);
        box->Draw("SAME");
    }

  for(int j=0; j<n; j++)
  {
    double xk,yk,xl,yl;
    if(lambda_pPb_v3[7]->GetN()-1 >= j){
      lambda_pPb_v3[7]->GetPoint(j,xl,yl);
      double percent1L = syst_pPb;
    }
    /*if(pPb_v3[5]->GetN()-1 >= j)
    {
      pPb_v3[5]->GetPoint(j,xh,yh);
      double percentH = syst_pPbH;
    }*/
    ks_pPb_v3[7]->GetPoint(j,xk,yk);
    double percent1K = syst_pPb;
    double yLerr, yKerr;

      yKerr=fabs(yk)*percent1K;
      yLerr=fabs(yl)*percent1L;
      //yHerr=fabs(yh)*percentH;

    if(lambda_pPb_v3[7]->GetN()-1 >= j){
      box = new TBox(xl-.15,yl-yLerr,xl+.15,yl+yLerr);
      box->SetFillColor(17);
      box->SetFillStyle(1001);
      box->SetLineWidth(0);
      box->Draw("SAME");
    }
	  /*if(pPb_v3[5]->GetN()-1 >= j){
      box = new TBox(xh-.15,yh-yHerr,xh+.15,yh+yHerr);
      box->SetFillColor(17);
      box->SetFillStyle(1001);
      box->SetLineWidth(0);
      box->Draw("SAME");
      }*/
    box = new TBox(xk-.15,yk-yKerr,xk+.15,yk+yKerr);
    box->SetFillColor(17);
    box->SetFillStyle(1001);
    box->SetLineWidth(0);
    box->Draw("SAME");
  }
//end add

  pPb_v3[5]->Draw("PESAME");
  ks_pPb_v3[7]->Draw("PESAME");
  lambda_pPb_v3[7]->Draw("PESAME");
  pad3v3[1]->cd();
  hist2QNS->Draw();
//begin add
  int n = ks_pPb_v3QNS[7]->GetN()-1;

  for(int j=0; j<n; j++)
  {
    double xk,yk,xl,yl;
    if(lambda_pPb_v3QNS[7]->GetN()-1 >= j){
      lambda_pPb_v3QNS[7]->GetPoint(j,xl,yl);
      double percent1L = syst_pPb;
    }
    ks_pPb_v3QNS[7]->GetPoint(j,xk,yk);
    double percent1K = syst_pPb;
    double yLerr, yKerr;

      yKerr=fabs(yk)*percent1K;
      yLerr=fabs(yl)*percent1L;

    if(lambda_pPb_v3QNS[7]->GetN()-1 >= j){
      box = new TBox(xl-.049,yl-yLerr,xl+.049,yl+yLerr);
      box->SetFillColor(17);
      box->SetFillStyle(1001);
      box->SetLineWidth(0);
      box->Draw("SAME");
    }
		
    box = new TBox(xk-.049,yk-yKerr,xk+.049,yk+yKerr);
    box->SetFillColor(17);
    box->SetFillStyle(1001);
    box->SetLineWidth(0);
    box->Draw("SAME");
  }

//end add

  ks_pPb_v3QNS[7]->Draw("PESAME");
  lambda_pPb_v3QNS[7]->Draw("PESAME");  
  fitfunc_pPb_v3QNS[7]->Draw("LSAME");


  pad3v3[2]->cd();
  histratio_v3->Draw();
  lline->Draw("same");
//begin add
  int n = ks_pPb_v3QNS_ratio[7]->GetN()-1;

  for(int j=0; j<n; j++)
  {
    double xk,yk,xl,yl;
    if(lambda_pPb_v3QNS_ratio[7]->GetN()-1 >= j){
      lambda_pPb_v3QNS_ratio[7]->GetPoint(j,xl,yl);
      double percent1L = syst_pPb;
    }
    ks_pPb_v3QNS_ratio[7]->GetPoint(j,xk,yk);
    double percent1K = syst_pPb;
    double yLerr, yKerr;
  
      yKerr=fabs(yk)*percent1K;
      yLerr=fabs(yl)*percent1L;

    if(lambda_pPb_v3QNS_ratio[7]->GetN()-1 >= j){
      box = new TBox(xl-.049,yl-yLerr,xl+.049,yl+yLerr);
      box->SetFillColor(17);
      box->SetFillStyle(1001);
      box->SetLineWidth(0);
      box->Draw("SAME");
    }
		
    box = new TBox(xk-.049,yk-yKerr,xk+.049,yk+yKerr);
    box->SetFillColor(17);
    box->SetFillStyle(1001);
    box->SetLineWidth(0);
    box->Draw("SAME");
  }
//end add

  ks_pPb_v3QNS_ratio[7]->Draw("PESAME");
  lambda_pPb_v3QNS_ratio[7]->Draw("PESAME");

  pad3v2[0]->cd();
  hist12->Draw();
//begin add
  int n = ks_PbPb_v2[3]->GetN();

	for(int j=0; j<n; j++)
  {
		double xk,yk,xl,yl,xh,yh;
		if(lambda_PbPb_v2[3]->GetN()-1 >= j){
			lambda_PbPb_v2[3]->GetPoint(j,xl,yl);
			double percent1L = syst_PbPb;
		}
    if(PbPb_v2[3]->GetN()-1 >= j)
    {
      PbPb_v2[3]->GetPoint(j,xh,yh);
      double percentH = syst_PbPbH;
    }
		ks_PbPb_v2[3]->GetPoint(j,xk,yk);
		double percent1K = syst_PbPb;
		double yLerr, yKerr, yHerr;
		
			yKerr=fabs(yk)*percent1K;
			yLerr=fabs(yl)*percent1L;
                        yHerr=fabs(yh)*percentH;

		if(lambda_PbPb_v2[3]->GetN()-1 >= j){
			box = new TBox(xl-.15,yl-yLerr,xl+.15,yl+yLerr);
			box->SetFillColor(17);
			box->SetFillStyle(1001);
			box->SetLineWidth(0);
			box->Draw("SAME");
		}
		if(PbPb_v2[3]->GetN()-1 >= j){
      box = new TBox(xh-.15,yh-yHerr,xh+.15,yh+yHerr);
      box->SetFillColor(17);
      box->SetFillStyle(1001);
      box->SetLineWidth(0);
      box->Draw("SAME");
    }
		box = new TBox(xk-.15,yk-yKerr,xk+.15,yk+yKerr);
		box->SetFillColor(17);
		box->SetFillStyle(1001);
		box->SetLineWidth(0);
		box->Draw("SAME");
	}
//end add
  PbPb_v2[3]->Draw("PESAME");
  ks_PbPb_v2[3]->Draw("PESAME");
  lambda_PbPb_v2[3]->Draw("PESAME");

  pad3v2[1]->cd();
  histQNS12->Draw();
//begin add
  int n = ks_pPb_v2QNS[5]->GetN()-1;

	for(int j=0; j<n; j++)
  {
		double xk,yk,xl,yl;
		if(lambda_pPb_v2QNS[3]->GetN()-1 >= j){
			lambda_pPb_v2QNS[3]->GetPoint(j,xl,yl);
			double percent1L = syst_pPb;
		}
		ks_pPb_v2QNS[3]->GetPoint(j,xk,yk);
		double percent1K = syst_pPb;
		double yLerr, yKerr;
		
			yKerr=fabs(yk)*percent1K;
			yLerr=fabs(yl)*percent1L;

		if(lambda_pPb_v2QNS[3]->GetN()-1 >= j){
			box = new TBox(xl-.049,yl-yLerr,xl+.049,yl+yLerr);
			box->SetFillColor(17);
			box->SetFillStyle(1001);
			box->SetLineWidth(0);
			box->Draw("PESAME");
		}
		
		box = new TBox(xk-.049,yk-yKerr,xk+.049,yk+yKerr);
		box->SetFillColor(17);
		box->SetFillStyle(1001);
		box->SetLineWidth(0);
		box->Draw("PESAME");
	}
//end add
  ks_pPb_v2QNS[3]->Draw("PESAME");
  lambda_pPb_v2QNS[3]->Draw("PESAME");
  fitfunc_pPb_v2QNS[3]->Draw("LSAME");


  pad3v2[2]->cd();
  histratio_v2single->Draw();
  lline->Draw("same");
//begin add
  int n = ks_pPb_v2QNS_ratio[3]->GetN()-1;

	for(int j=0; j<n; j++)
  {
		double xk,yk,xl,yl;
		if(lambda_pPb_v2QNS_ratio[3]->GetN()-1 >= j){
			lambda_pPb_v2QNS_ratio[3]->GetPoint(j,xl,yl);
			double percent1L = syst_pPb;
		}
		ks_pPb_v2QNS_ratio[3]->GetPoint(j,xk,yk);
		double percent1K = syst_pPb;
		double yLerr, yKerr;
		
			yKerr=fabs(yk)*percent1K;
			yLerr=fabs(yl)*percent1L;

		if(lambda_pPb_v2QNS_ratio[3]->GetN()-1 >= j){
			box = new TBox(xl-.049,yl-yLerr,xl+.049,yl+yLerr);
			box->SetFillColor(17);
			box->SetFillStyle(1001);
			box->SetLineWidth(0);
			box->Draw("SAME");
		}
		
		box = new TBox(xk-.049,yk-yKerr,xk+.049,yk+yKerr);
		box->SetFillColor(17);
		box->SetFillStyle(1001);
		box->SetLineWidth(0);
		box->Draw("SAME");
	}
//end add
  ks_pPb_v2QNS_ratio[3]->Draw("PESAME");
  lambda_pPb_v2QNS_ratio[3]->Draw("PESAME");

  char label_energy[2][200]={"CMS PbPb #sqrt{s_{NN}} = 2.76 TeV","CMS pPb #sqrt{s_{NN}} = 5.02 TeV"};
  char label_energy1[2][200]={"#sqrt{s_{NN}} = 2.76 TeV, L_{int} = 2.3 #mub^{#font[122]{\55}1}","#sqrt{s_{NN}} = 5.02 TeV, L_{int} = 35 nb^{#font[122]{\55}1}"};
  char label_n[4][200]={"120 #leq N^{offline}_{trk} < 150","150 #leq N^{offline}_{trk} < 185","185 #leq N^{offline}_{trk} < 220","220 #leq N^{offline}_{trk} < 260"};
  char label_n1[3][200]={"N^{offline}_{trk} < 35","35 #leq N^{offline}_{trk} < 60","60 #leq N^{offline}_{trk} < 120"};
  char label_frac[4][200]={"(0.5-2.5%)","(0.06-0.5%)","(0.006-0.06%)","(0.0003-0.006%)"};
  char label_frac1[3][200]={"(48-100%)","(24.5-48%)","(2.5-24.5%)"};
  char label_frac_PbPb[4][200]={"(67#pm3%)","(64#pm3%)","(62#pm2%)","(59#pm2%)"};
  char label_frac1_PbPb[3][200]={"(90#pm5%)","(80#pm4%)","(73#pm4%)"};

  TLatex *tex2= new TLatex();
  tex2->SetNDC();
    tex2->SetTextFont(42);
  pad[0]->cd();
  tex2->SetTextSize(tex2->GetTextSize()*1.2);
    tex2->DrawLatex(0.235,0.89,"CMS pPb #sqrt{s_{NN}} = 5.02 TeV");
    //tex2->SetTextSize(tex2->GetTextSize()/1.4*1.2);
  tex2->DrawLatex(0.235,0.80,"L_{int} = 35 nb^{#font[122]{\55}1}");
  tex2->DrawLatex(0.55,0.32,label_n[0]);
  tex2->DrawLatex(0.55,0.24,label_frac[0]);
  pad[1]->cd();
    tex2->SetTextSize(tex2->GetTextSize()*1.2);
  tex2->DrawLatex(0.4,0.32,label_n[1]);
  tex2->DrawLatex(0.4,0.24,label_frac[1]);
  pad[2]->cd();
  tex2->DrawLatex(0.4,0.32,label_n[2]);
  tex2->DrawLatex(0.4,0.24,label_frac[2]);
  pad[3]->cd();
  tex2->DrawLatex(0.4,0.32,label_n[3]);
  tex2->DrawLatex(0.4,0.24,label_frac[3]);

  pad[0]->cd();
  TLegend *leg1 = new TLegend(0.22,0.52,0.46,0.75);
  leg1->SetFillColor(10);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0.035);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.07);
  leg1->AddEntry(ks_pPb_v2[0],"K_{S}^{0}","p");
  leg1->AddEntry(lambda_pPb_v2[0],"#Lambda/#bar{#Lambda}","P");
  leg1->AddEntry(pPb_v2[3],"h^{#pm}","P");
  leg1->Draw();
  pad[4]->cd();
  TLegend *leg2 = new TLegend(0.43,0.1,0.88,0.285);
  leg2->SetFillColor(10);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0.035);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.073);
    leg2->AddEntry(fitfunc_pPb_v2QNS[3]," "," ");
  leg2->AddEntry(fitfunc_pPb_v2QNS[3],"Polynomial fits to K_{S}^{0}","L");
  leg2->Draw();

  TLatex *tex1= new TLatex();
  tex1->SetNDC();
    tex1->SetTextFont(42);
  pad1[0]->cd();
  tex1->SetTextSize(tex1->GetTextSize()*1.2);
    tex1->DrawLatex(0.235,0.89,"CMS PbPb #sqrt{s_{NN}} = 2.76 TeV");
    //tex1->SetTextSize(tex1->GetTextSize()/1.4*1.3);
  tex1->DrawLatex(0.235,0.80,"L_{int} = 2.3 #mub^{#font[122]{\55}1}");
  tex1->DrawLatex(0.55,0.32,label_n[0]);
  tex1->DrawLatex(0.55,0.24,label_frac_PbPb[0]);
  tex1->SetTextSize(tex1->GetTextSize()/1.25*1.35);
  pad1[1]->cd();
  tex1->DrawLatex(0.4,0.32,label_n[1]);
  tex1->DrawLatex(0.4,0.24,label_frac_PbPb[1]);
  pad1[2]->cd();
  tex1->DrawLatex(0.4,0.32,label_n[2]);
  tex1->DrawLatex(0.4,0.24,label_frac_PbPb[2]);
  pad1[3]->cd();
  tex1->DrawLatex(0.4,0.32,label_n[3]);
  tex1->DrawLatex(0.4,0.24,label_frac_PbPb[3]);

  pad1[0]->cd();
  TLegend *leg11 = new TLegend(0.22,0.52,0.46,0.75);
  leg11->SetFillColor(10);
  leg11->SetFillStyle(0);
  leg11->SetBorderSize(0.035);
  leg11->SetTextFont(42);
  leg11->SetTextSize(0.06);
  leg11->AddEntry(ks_PbPb_v2[0],"K_{S}^{0}","p");
  leg11->AddEntry(lambda_PbPb_v2[0],"#Lambda/#bar{#Lambda}","P");
  leg11->AddEntry(PbPb_v2[3],"h^{#pm}","P");
  leg11->Draw();
  pad1[4]->cd();
  TLegend *leg12 = new TLegend(0.43,0.1,0.88,0.285);
  leg12->SetFillColor(10);
  leg12->SetFillStyle(0);
  leg12->SetBorderSize(0.035);
  leg12->SetTextFont(42);
  leg12->SetTextSize(0.073);
    leg12->AddEntry(fitfunc_pPb_v2QNS[3]," "," ");
  leg12->AddEntry(fitfunc_PbPb_v2QNS[3],"Polynomial fits to K_{S}^{0}","L");
  leg12->Draw();

  TLatex *tex4= new TLatex();
  tex4->SetNDC();
    tex4->SetTextFont(42);
    TLatex *tex41= new TLatex();
    tex41->SetNDC();
    tex41->SetTextFont(42);
  
  c4v2->cd(1);
  tex4->SetTextSize(tex4->GetTextSize()*1.5);
  tex4->DrawLatex(0.19,0.885,"CMS PbPb #sqrt{s_{NN}} = 2.76 TeV");
    //tex4->SetTextSize(tex4->GetTextSize()/1.5*1.4);
  tex4->DrawLatex(0.19,0.775,"L_{int} = 2.3 #mub^{#font[122]{\55}1}");
  tex4->DrawLatex(0.61,0.185,label_n1[0]);
  tex4->DrawLatex(0.61,0.095,label_frac1_PbPb[0]);
  c4v2->cd(2);
  tex4->DrawLatex(0.51,0.785,label_n1[1]);
  tex4->DrawLatex(0.51,0.695,label_frac1_PbPb[1]);
  c4v2->cd(3);
  tex4->DrawLatex(0.46,0.785,label_n1[2]);
  tex4->DrawLatex(0.46,0.695,label_frac1_PbPb[2]);
  c4v2->cd(4);
  tex41->SetTextSize(tex41->GetTextSize()*1.4);
    tex41->DrawLatex(0.19,0.905,"CMS pPb #sqrt{s_{NN}} = 5.02 TeV");
    //tex41->SetTextSize(tex41->GetTextSize()/1.5*1.45);
  tex41->DrawLatex(0.19,0.795,"L_{int} = 35 nb^{#font[122]{\55}1}");
    tex41->SetTextSize(tex41->GetTextSize()*1);
  tex41->DrawLatex(0.61,0.305,label_n1[0]);
  tex41->DrawLatex(0.61,0.23,label_frac1[0]);
  c4v2->cd(5);
  tex41->DrawLatex(0.51,0.80,label_n1[1]);
  tex41->DrawLatex(0.51,0.725,label_frac1[1]);
  c4v2->cd(6);
  tex41->DrawLatex(0.46,0.80,label_n1[2]);
  tex41->DrawLatex(0.46,0.725,label_frac1[2]);

  c4v2->cd(1);
  TLegend *leg4 = new TLegend(0.17,0.41,0.50,0.72);
  leg4->SetFillColor(10);
  leg4->SetFillStyle(0);
  leg4->SetBorderSize(0.035);
  leg4->SetTextFont(42);
  leg4->SetTextSize(0.09);
  leg4->AddEntry(ks_pPb_v2[0],"K_{S}^{0}","p");
  leg4->AddEntry(lambda_pPb_v2[0],"#Lambda/#bar{#Lambda}","P");
    leg4->AddEntry(pPb_v2[0],"h^{#pm}","P");
  leg4->Draw();

  pad3v3[0]->cd();
  TLegend *leg111 = new TLegend(0.17,0.56,0.51,0.79);
  leg111->SetFillColor(10);
  leg111->SetFillStyle(0);
  leg111->SetBorderSize(0.035);
  leg111->SetTextFont(42);
  leg111->SetTextSize(0.055);
  leg111->AddEntry(ks_pPb_v3[5],"K_{S}^{0}","p");
  leg111->AddEntry(lambda_pPb_v3[5],"#Lambda/#bar{#Lambda}","P");
  leg111->AddEntry(pPb_v3[5],"h^{#pm}","P");
  leg111->Draw();
  TLatex *tex5 = new TLatex();
  tex5->SetNDC();
    tex5->SetTextFont(42);
  //tex5->SetTextSize(tex5->GetTextSize()*1.2);
    tex5->DrawLatex(0.2,0.915,"CMS pPb #sqrt{s_{NN}} = 5.02 TeV");
  tex5->DrawLatex(0.2,0.84,"L_{int} = 35 nb^{#font[122]{\55}1}");
  //tex5->SetTextSize(tex5->GetTextSize()/1.2);
  tex5->DrawLatex(0.55,0.775,"185 #leq N^{offline}_{trk} < 350");
  tex5->DrawLatex(0.55,0.71,"(0-0.06%)");
  pad3v3[1]->cd();
  TLegend *leg122 = new TLegend(0.46,0.08,0.8,0.265);
  leg122->SetFillColor(10);
  leg122->SetFillStyle(0);
  leg122->SetBorderSize(0.035);
  leg122->SetTextFont(42);
  leg122->SetTextSize(0.072);
    leg122->AddEntry(fitfunc_pPb_v3QNS[5]," "," ");
  leg122->AddEntry(fitfunc_pPb_v3QNS[5],"Polynomial fit to K_{S}^{0}","L");
  leg122->Draw();

  pad3v2[0]->cd();
  TLegend *leg1111 = new TLegend(0.18,0.56,0.51,0.79);
  leg1111->SetFillColor(10);
  leg1111->SetFillStyle(0);
  leg1111->SetBorderSize(0.035);
  leg1111->SetTextFont(42);
  leg1111->SetTextSize(0.065);
  leg1111->AddEntry(ks_pPb_v2[5],"K_{S}^{0}","p");
  leg1111->AddEntry(lambda_pPb_v2[5],"#Lambda/#bar{#Lambda}","P");
  leg1111->AddEntry(pPb_v2[5],"h^{#pm}","P");
  leg1111->Draw();
  TLatex *tex6 = new TLatex();
  tex6->SetNDC();
    tex6->SetTextFont(42);
  tex6->SetTextSize(tex6->GetTextSize()*0.95);
    tex5->DrawLatex(0.2,0.92,"CMS");
  tex6->DrawLatex(0.2,0.85,Form("PbPb %s",label_energy1[0]));
  tex6->SetTextSize(tex6->GetTextSize()*1.0);
  tex6->DrawLatex(0.5,0.26,"120 #leq N^{offline}_{trk} < 150");
  tex6->DrawLatex(0.5,0.195,"(67#pm3%)");
  pad3v2[1]->cd();
  TLegend *leg1221 = new TLegend(0.46,0.08,0.9,0.265);
  leg1221->SetFillColor(10);
  leg1221->SetFillStyle(0);
  leg1221->SetBorderSize(0.035);
  leg1221->SetTextFont(42);
  leg1221->SetTextSize(0.065);
  leg1221->AddEntry(fitfunc_pPb_v2QNS[5],"Polynomial fits to K_{S}^{0}","L");
  leg1221->Draw();

    //return;
  /*c1->Print("./paperplots/v2_PID_highN_PbPb.gif");
  c1->Print("./paperplots/v2_PID_highN_PbPb.pdf");
  c2->Print("./paperplots/v2_PID_highN_pPb.gif");
  c2->Print("./paperplots/v2_PID_highN_pPb.pdf");

  c4v2->Print("./paperplots/v2_PID_lowN.gif");
  c4v2->Print("./paperplots/v2_PID_lowN.pdf");
  c3v3->Print("./paperplots/v3_PID_highN.gif");
  c3v3->Print("./paperplots/v3_PID_highN.pdf");
  c3v2->Print("v2_PID_highN_single.gif");
  c3v2->Print("v2_PID_highN_single.pdf");
  c3v2->Print("v2_PID_highN_single.C");*/
    
    

}
