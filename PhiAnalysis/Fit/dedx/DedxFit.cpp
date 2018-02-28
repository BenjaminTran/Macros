//Includes
#include <TLatex.h>
#include <TStyle.h>
#include "TColor.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TMathText.h"
#include "TPad.h"
#include <TString.h>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <fstream>

void DedxFit()
{
    TFile *f = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/PhiDeDx/v3/Phi_DeDx_pPb_v3_PD1_020918.root");

    TH2D* h_Dedx;
    TCanvas* c = new TCanvas("c","c",800,600);

    f->GetObject("PhiSelector/Dedx_harm",h_Dedx);
    h_Dedx->SetMaximum(1.9e6);
    h_Dedx->SetContour(99);
    h_Dedx->GetXaxis()->SetTitle("p");

    //TF1* f_Dedx_top = new TF1("f_Dedx","1.1*TMath::Power(1.019/x,2) + 3.1",0.4,5);
    //TF1* f_Dedx_top1 = new TF1("f_Dedx","0.55*(TMath::Power(2.15/(x+0.02),2) - 2*TMath::Power(2.2/(x-0.02),1)) + 4.3",0.3,5); // makes the bottom come up instead of lay flat. Must be due to sigma_2 > sigma_1
    //TF1* f_Dedx_top1 = new TF1("f_Dedx","0.55*(TMath::Power(2.15/(x),2) - 2*TMath::Power(1.17/(x+0.17),1)) + 3.7",0.3,5);
    //TF1* f_Dedx_top2 = new TF1("f_Dedx","0.85*(TMath::Power(0.95/(x),4) - 2.4*TMath::Power(0.6/(x),2)) + 3.6",0.3,5);
    //TF1* f_Dedx_top1 = new TF1("f_Dedx","0.85*(TMath::Power(0.95/(x),2)) + 3.6",0.3,5);
    TF1* f_Dedx_top2 = new TF1("f_Dedx","0.55*(TMath::Power(1.62/(x),2) - 2.95*TMath::Power(0.6/(x),1)) + 3.6",0.1,5); //magenta
    TF1* f_Dedx_top3 = new TF1("f_Dedx","0.55*(TMath::Power(1.62/x,2) - 2*TMath::Power(0.6/x,1)) + 3.6",0.1,5);

    //TF1* f_Dedx_bot = new TF1("f_Dedx","1.1*TMath::Power(1.019/x,2.2) + 1.1",0.4,5);
    //TF1* f_Dedx_bot1 = new TF1("f_Dedx","0.5*(TMath::Power(0.8/x,4) - 2*TMath::Power(0.6/x,2)) + 3",0.4,5);
    //TF1* f_Dedx_bot2 = new TF1("f_Dedx","0.5*(TMath::Power(0.5/x,4) - 1*TMath::Power(0.50/x,2)) + 2.7",0.1,5); //Loose
    TF1* f_Dedx_bot2 = new TF1("f_Dedx","0.55*(TMath::Power(1.15/x,2) - 1.7*TMath::Power(0.6/x,1)) + 2.5",0.1,5);
    //TF1* f_Dedx_bot3 = new TF1("f_Dedx","0.55*(TMath::Power(1.15/x,2) - 2*TMath::Power(0.6/x,1)) + 3",0.1,5);
    TF1* f_Dedx_bot3 = new TF1("f_Dedx","0.55*(TMath::Power(1.15/x,2) - 1.7*TMath::Power(0.6/x,1)) + 2.7",0.1,5); //tight
    //TF1* f_Dedx_bot4 = new TF1("f_Dedx","0.5*(TMath::Power(0.76/x,4) - 2*TMath::Power(0.5/x,2)) + 3",0.4,5);

    //f_Dedx_top1->SetLineColor(6);
    f_Dedx_top2->SetLineColor(6);
    f_Dedx_top3->SetLineColor(2);
    //f_Dedx_bot1->SetLineColor(2);
    f_Dedx_bot2->SetLineColor(4);
    f_Dedx_bot3->SetLineColor(3);
    //f_Dedx_bot4->SetLineColor(3);

    c->cd();
    h_Dedx->Draw("cont0 SCAT");
    //f_Dedx_top->Draw("SAME");
    //f_Dedx_top1->Draw("SAME");
    f_Dedx_top2->Draw("SAME");
    f_Dedx_top3->Draw("SAME");
    //f_Dedx_bot4->Draw("SAME");
    f_Dedx_bot3->Draw("SAME");
    f_Dedx_bot2->Draw("SAME");
    //f_Dedx_bot1->Draw("SAME");

    c->Print("Image/DeDxPrelimFit_v1.pdf");
    c->Print("Image/DeDxPrelimFit_v1.png");

}
