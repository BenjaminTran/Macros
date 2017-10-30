#include <iostream>
#include <vector>
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"

using namespace std;

struct ParticleData{
    std::string fn = "V0v2perisub.root";
    std::string fn_V0 = "v0CorrelationRapidity";
    std::string fn_Xi = "xiCorrelationRapidity";
    std::string fn_v0cas = "v0CasCorrelationRapidityPeriSub";
    std::vector<double> fsig_ks = {0.989638, 0.992109, 0.992861, 0.992984, 0.992159, 0.989499, 0.986445, 0.982275, 0.976946, 0.970546, 0.964845, 0.958563, 0.969557};
    std::vector<double> fsig_la = {0.908266, 0.96779, 0.979912, 0.981899, 0.982888, 0.982854, 0.980766, 0.97569, 0.97569, 0.964048};
    std::vector<double> fsig_xi = {0.959427 ,0.976239 ,0.979161 ,0.980678 ,0.980661 ,0.981534 ,0.981502 ,0.979289 ,0.979192};
    std::vector<double> fsig_ks_MB = {0.999 ,0.999 ,0.995 ,0.996 ,0.995 ,0.992 ,0.990 ,0.987 ,0.983 ,0.980 ,0.976 ,0.965 ,0.965};
    std::vector<double> fsig_la_MB = {0.990,0.995,0.995,0.999,0.992,0.989,0.985 ,0.998,0.998 ,0.995};
    std::vector<double> PtBin_ks = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5};//, 10.0, 15.0, 20.0}; 
    std::vector<double> PtBin_la = {0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5};//, 10.0, 15.0, 20.0}; 
    std::vector<double> PtBin_xi = {1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.2,10.0};
    TFile *f_perisub = TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/PeripheralSubtractionMB.root");
    TFile *f_V0 = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityCorrectMultB_09_19_17.root");
    TFile* f_Xi = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationRapidityTotal_08_20_2017.root");
    TFile *f_low_ref = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/HadCorr/Combine_MB0_corr_ref.root");
    TFile *f_high_ref = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/HadCorr/XiCorrelationRapidityTotal_08_20_2017.root");
} V0;

void PeriSubObs(ParticleData PD)
{
    TH1::SetDefaultSumw2();

    double BW2D = (9.9/33)*((2-1.0/16)*3.1416/31);

    TLatex* tex = new TLatex();
    tex->SetNDC();


    TCanvas* csPeaksr_low_ref = new TCanvas("csPeaksr_low_ref","csPeaksr_low_ref",600,600);
    TCanvas* csPeaklr_low_ref = new TCanvas("csPeaklr_low_ref","csPeaklr_low_ref",600,600);
    TCanvas* csPeaksr_high_ref = new TCanvas("csPeaksr_high_ref","csPeaksr_high_ref",600,600);
    TCanvas* csPeaklr_high_ref = new TCanvas("csPeaklr_high_ref","csPeaklr_high_ref",600,600);
    TCanvas* csPeaksr_low_ks = new TCanvas("csPeaksr_low_ks","csPeaksr_low_ks",1200,900);
    TCanvas* csPeaklr_low_ks = new TCanvas("csPeaklr_low_ks","csPeaklr_low_ks",1200,900);
    TCanvas* csPeaksr_ks = new TCanvas("csPeaksr_ks","csPeaksr_ks",1200,900);
    TCanvas* csPeaklr_ks = new TCanvas("csPeaklr_ks","csPeaklr_ks",1200,900);
    TCanvas* csPeaksr_low_la = new TCanvas("csPeaksr_low_la","csPeaksr_low_la",1200,900);
    TCanvas* csPeaklr_low_la = new TCanvas("csPeaklr_low_la","csPeaklr_low_la",1200,900);
    TCanvas* csPeaksr_la = new TCanvas("csPeaksr_la","csPeaksr_la",1200,900);
    TCanvas* csPeaklr_la = new TCanvas("csPeaklr_la","csPeaklr_la",1200,900);
    TCanvas* csPeaksr_low_xi = new TCanvas("csPeaksr_low_xi","csPeaksr_low_xi",1200,900);
    TCanvas* csPeaklr_low_xi = new TCanvas("csPeaklr_low_xi","csPeaklr_low_xi",1200,900);
    TCanvas* csPeaksr_xi = new TCanvas("csPeaksr_xi","csPeaksr_xi",1200,900);
    TCanvas* csPeaklr_xi = new TCanvas("csPeaklr_xi","csPeaklr_xi",1200,900);
    TCanvas* c = new TCanvas("c","c",400,400);

    csPeaksr_low_ks->Divide(4,4);
    csPeaklr_low_ks->Divide(4,4);
    csPeaksr_ks->Divide(4,4);
    csPeaklr_ks->Divide(4,4);
    csPeaksr_low_la->Divide(4,3);
    csPeaklr_low_la->Divide(4,3);
    csPeaksr_la->Divide(4,3);
    csPeaklr_la->Divide(4,3);
    csPeaksr_low_xi->Divide(4,3);
    csPeaklr_low_xi->Divide(4,3);
    csPeaksr_xi->Divide(4,3);
    csPeaklr_xi->Divide(4,3);

    //Low N ref
    //Containers
    std::vector<double> Nassoc_low_ref;
    std::vector<double> Jyieldsr_low_ref;
    std::vector<double> Jyieldlr_low_ref;
    double Jyieldsr_err_low_ref;
    double Jyieldlr_err_low_ref;
    std::vector<double> V2Values_low_ref;
    std::vector<double> V2Values_err_low_ref;
    std::vector<double> V3Values_low_ref;
    std::vector<double> V3Values_err_low_ref;
    std::vector<double> Bz_low_ref;
    std::vector<double> nEvent_low_ref;

    std::vector<double> JyieldSub_low_ref;
    std::vector<double> JyieldSub_err_low_ref;


    TH1D* hsPeaksr_low_ref;
    TH1D* hbPeaksr_low_ref;
    TH1D* hsPeaklr_low_ref;
    TH1D* hbPeaklr_low_ref;
    TH1D* V2lrs_low_ref;
    TH1D* V2lrb_low_ref;
    TH2D* hbackgroundPeak_low_ref;
    TH2D* hsignalPeak_low_ref;

    //Calculate Nassoc, Jet yield, Low N

    PD.f_low_ref->GetObject("pPbCorr/background",hbackgroundPeak_low_ref);
    PD.f_low_ref->GetObject("pPbCorr/signal",hsignalPeak_low_ref);
    TH1D* mult_low_ref = (TH1D*) PD.f_low_ref->Get("pPbCorr/mult");

    hbPeaksr_low_ref = hbackgroundPeak_low_ref->ProjectionY("hbPeaksr_low_ref", 14, 20);
    hsPeaksr_low_ref = hsignalPeak_low_ref->ProjectionY("hsPeaksr_low_ref", 14, 20);

    TF1* quadFit1 = new TF1("quadFit1","[0]*x^2+[1]*x+[2]",0.6,2.2);
    quadFit1->SetParameters(1,1,1);

    //nEvent_low_ref.push_back(mult_low_ref->Integral(0,100000));
    nEvent_low_ref.push_back(mult_low_ref->Integral(3,100000));
    Bz_low_ref.push_back(hbackgroundPeak_low_ref->GetBinContent(hbackgroundPeak_low_ref->FindBin(0,0)));

    hsPeaksr_low_ref->Divide(hbPeaksr_low_ref);
    hsPeaksr_low_ref->Scale(Bz_low_ref[0]/nEvent_low_ref[0]/BW2D);

    hsPeaksr_low_ref->Fit("quadFit1","R");
    hsPeaksr_low_ref->Fit("quadFit1","R");
    hsPeaksr_low_ref->Fit("quadFit1","R");


    double minVal_sr = quadFit1->GetMinimum(0.6,2.2);
    double minVal_srX = quadFit1->GetMinimumX(0.6,2.2);
    TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    minConst_sr->SetParameter(0,-minVal_sr);
    TH1D* hsPeaksr_zeroed_low_ref = (TH1D*)hsPeaksr_low_ref->Clone();
    hsPeaksr_zeroed_low_ref->Add(minConst_sr);
    csPeaksr_low_ref->cd();
    hsPeaksr_low_ref->Draw();
    double xcoor = 0.52;
    double ycoor = 0.90;
    double increment = 0.07;
    tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
    tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",0.3,3.0));
    tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
    c->cd();
    Jyieldsr_low_ref.push_back(hsPeaksr_zeroed_low_ref->IntegralAndError(hsPeaksr_zeroed_low_ref->FindBin(0.0),hsPeaksr_zeroed_low_ref->FindBin(minVal_srX),Jyieldsr_err_low_ref,"width"));
    double bin0yield = hsPeaksr_zeroed_low_ref->GetBinContent(hsPeaksr_zeroed_low_ref->FindBin(0.0))*0.19635;
    Jyieldsr_low_ref[0] = Jyieldsr_low_ref[0]*2 - bin0yield;

    hbPeaklr_low_ref = hbackgroundPeak_low_ref->ProjectionY("hbPeaklr_low_ref",1,14);
    TH1D* ahbPeaklr_low_ref = hbackgroundPeak_low_ref->ProjectionY("ahbPeaklr_low_ref",20,33);
    hsPeaklr_low_ref = hsignalPeak_low_ref->ProjectionY("hsPeaklr_low_ref",1,14);
    TH1D* ahsPeaklr_low_ref = hsignalPeak_low_ref->ProjectionY("ahsPeaklr_low_ref",20,33);

    hbPeaklr_low_ref->Add(ahbPeaklr_low_ref);
    hsPeaklr_low_ref->Add(ahsPeaklr_low_ref);
    hsPeaklr_low_ref->Divide(hbPeaklr_low_ref);
    hsPeaklr_low_ref->Scale(Bz_low_ref[0]/nEvent_low_ref[0]/BW2D);

    TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
    quadFit2->SetParameters(1,1,1);

    hsPeaklr_low_ref->Fit("quadFit2","R");
    hsPeaklr_low_ref->Fit("quadFit2","R");
    hsPeaklr_low_ref->Fit("quadFit2","R");

    double minVal_lr = quadFit2->GetMinimum(0.6,2.2);
    double minVal_lrX = quadFit2->GetMinimumX(0.6,2.2);
    TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    minConst_lr->SetParameter(0,-minVal_lr);
    TH1D* hsPeaklr_zeroed_low_ref = (TH1D*)hsPeaklr_low_ref->Clone();
    hsPeaklr_zeroed_low_ref->Add(minConst_lr);
    csPeaklr_low_ref->cd();
    hsPeaklr_low_ref->Draw();
    xcoor = 0.52;
    ycoor = 0.90;
    increment = 0.07;
    tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
    tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",0.3,3.0));
    tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
    c->cd();
    Jyieldlr_low_ref.push_back(hsPeaklr_zeroed_low_ref->IntegralAndError(hsPeaklr_zeroed_low_ref->FindBin(0.0),hsPeaklr_zeroed_low_ref->FindBin(minVal_lrX),Jyieldlr_err_low_ref,"width"));
    bin0yield = hsPeaklr_zeroed_low_ref->GetBinContent(hsPeaklr_zeroed_low_ref->FindBin(0.0))*0.19635;
    Jyieldlr_low_ref[0] = Jyieldlr_low_ref[0]*2 - bin0yield;

    JyieldSub_low_ref.push_back(Jyieldsr_low_ref[0] - Jyieldlr_low_ref[0]);
    JyieldSub_err_low_ref.push_back(sqrt(Jyieldsr_err_low_ref*Jyieldsr_err_low_ref + Jyieldlr_err_low_ref*Jyieldlr_err_low_ref));

    TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit1->SetParNames("N","V1","V2","V3");
    fit1->SetParameters(10,1,1,1);
    fit1->SetLineColor(2);

    V2lrs_low_ref = hsignalPeak_low_ref->ProjectionY("V2lrs_low_ref",1,14);
    TH1D* aV2lrs_low_ref = hsignalPeak_low_ref->ProjectionY("aV2lrs_low_ref",20,33);
    V2lrb_low_ref = hbackgroundPeak_low_ref->ProjectionY("V2lrb_low_ref",1,14);
    TH1D* aV2lrb_low_ref = hbackgroundPeak_low_ref->ProjectionY("aV2lrb_low_ref",20,33);
    V2lrs_low_ref->Add(aV2lrs_low_ref);
    V2lrb_low_ref->Add(aV2lrb_low_ref);
    V2lrs_low_ref->Divide(V2lrb_low_ref);
    V2lrs_low_ref->Scale(Bz_low_ref[0]/nEvent_low_ref[0]/BW2D);

    V2lrs_low_ref->Fit("fit1","R");
    V2lrs_low_ref->Fit("fit1","R");
    V2lrs_low_ref->Fit("fit1","R");

    V2Values_low_ref.push_back(fit1->GetParameter(2));
    V2Values_err_low_ref.push_back(fit1->GetParError(2));
    V3Values_low_ref.push_back(fit1->GetParameter(3));
    V3Values_err_low_ref.push_back(fit1->GetParError(3));

    double v2_low_ref = sqrt(V2Values_low_ref[0]);
    double v2e_low_ref = sqrt(V2Values_low_ref[0])*(V2Values_err_low_ref[0]/V2Values_low_ref[0])/2;
    double v3_low_ref = sqrt(V3Values_low_ref[0]);
    double v3e_low_ref = sqrt(V3Values_low_ref[0])*(V3Values_err_low_ref[0]/V3Values_low_ref[0])/2;

    Nassoc_low_ref.push_back(fit1->GetParameter(0));
    c->cd();

    //High N ref
    //Containers
    std::vector<double> Nassoc_high_ref;
    std::vector<double> Jyieldsr_high_ref;
    std::vector<double> Jyieldlr_high_ref;
    double Jyieldsr_err_high_ref;
    double Jyieldlr_err_high_ref;
    std::vector<double> V2Values_high_ref;
    std::vector<double> V2Values_err_high_ref;
    std::vector<double> V3Values_high_ref;
    std::vector<double> V3Values_err_high_ref;
    std::vector<double> Bz_high_ref;
    std::vector<double> nEvent_high_ref;

    std::vector<double> JyieldSub_high_ref;
    std::vector<double> JyieldSub_err_high_ref;


    TH1D* hsPeaksr_high_ref;;
    TH1D* hbPeaksr_high_ref;
    TH1D* hsPeaklr_high_ref;
    TH1D* hbPeaklr_high_ref;
    TH1D* V2lrs_high_ref;
    TH1D* V2lrb_high_ref;
    TH2D* hbackgroundPeak_high_ref;
    TH2D* hsignalPeak_high_ref;

    //Calculate Nassoc, Jet yield, High N

    PD.f_high_ref->GetObject("xiCorrelationRapidity/BackgroundHad",hbackgroundPeak_high_ref);
    PD.f_high_ref->GetObject("xiCorrelationRapidity/SignalHad",hsignalPeak_high_ref);
    TH1D* mult_high_ref = (TH1D*) PD.f_high_ref->Get("xiCorrelationRapidity/nTrk");

    hbPeaksr_high_ref = hbackgroundPeak_high_ref->ProjectionY("hbPeaksr_high_ref", 14, 20);
    hsPeaksr_high_ref = hsignalPeak_high_ref->ProjectionY("hsPeaksr_high_ref", 14, 20);

    TF1* quadFit12 = new TF1("quadFit12","[0]*x^2+[1]*x+[2]",0.6,2.2);
    quadFit12->SetParameters(1,1,1);

    //nEvent_high_ref.push_back(mult_high_ref->Integral(0,100000));
    nEvent_high_ref.push_back(mult_high_ref->Integral(3,100000));
    Bz_high_ref.push_back(hbackgroundPeak_high_ref->GetBinContent(hbackgroundPeak_high_ref->FindBin(0,0)));

    hsPeaksr_high_ref->Divide(hbPeaksr_high_ref);
    hsPeaksr_high_ref->Scale(Bz_high_ref[0]/nEvent_high_ref[0]/BW2D);

    hsPeaksr_high_ref->Fit("quadFit12","R");
    hsPeaksr_high_ref->Fit("quadFit12","R");
    hsPeaksr_high_ref->Fit("quadFit12","R");

    minVal_sr = quadFit12->GetMinimum(0.6,2.2);
    minVal_srX = quadFit12->GetMinimumX(0.6,2.2);
    TF1* minConst_sr2 = new TF1("minConst_sr2","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    minConst_sr2->SetParameter(0,-minVal_sr);
    TH1D* hsPeaksr_zeroed_high_ref = (TH1D*)hsPeaksr_high_ref->Clone();
    hsPeaksr_zeroed_high_ref->Add(minConst_sr2);
    csPeaksr_high_ref->cd();
    hsPeaksr_high_ref->Draw();
    xcoor = 0.52;
    ycoor = 0.90;
    increment = 0.07;
    tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
    tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",0.3,3.0));
    tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
    c->cd();
    Jyieldsr_high_ref.push_back(hsPeaksr_zeroed_high_ref->IntegralAndError(hsPeaksr_zeroed_high_ref->FindBin(0.0),hsPeaksr_zeroed_high_ref->FindBin(minVal_srX),Jyieldsr_err_high_ref,"width"));
    bin0yield = hsPeaksr_zeroed_high_ref->GetBinContent(hsPeaksr_zeroed_high_ref->FindBin(0.0))*0.19635;
    Jyieldsr_high_ref[0] = Jyieldsr_high_ref[0]*2 - bin0yield;

    hbPeaklr_high_ref = hbackgroundPeak_high_ref->ProjectionY("hbPeaklr_high_ref",1,14);
    TH1D* ahbPeaklr_high_ref = hbackgroundPeak_high_ref->ProjectionY("ahbPeaklr_high_ref",20,33);
    hsPeaklr_high_ref = hsignalPeak_high_ref->ProjectionY("hsPeaklr_high_ref",1,14);
    TH1D* ahsPeaklr_high_ref = hsignalPeak_high_ref->ProjectionY("ahsPeaklr_high_ref",20,33);

    hbPeaklr_high_ref->Add(ahbPeaklr_high_ref);
    hsPeaklr_high_ref->Add(ahsPeaklr_high_ref);
    hsPeaklr_high_ref->Divide(hbPeaklr_high_ref);
    hsPeaklr_high_ref->Scale(Bz_high_ref[0]/nEvent_high_ref[0]/BW2D);

    TF1* quadFit21 = new TF1("quadFit21","[0]*x^2+[1]*x+[2]",0.6,2.2);
    quadFit21->SetParameters(1,1,1);

    hsPeaklr_high_ref->Fit("quadFit21","R");
    hsPeaklr_high_ref->Fit("quadFit21","R");
    hsPeaklr_high_ref->Fit("quadFit21","R");

    minVal_lr = quadFit21->GetMinimum(0.6,2.2);
    minVal_lrX = quadFit21->GetMinimumX(0.6,2.2);
    TF1* minConst_lr2 = new TF1("minConst_lr2","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    minConst_lr2->SetParameter(0,-minVal_lr);
    TH1D* hsPeaklr_zeroed_high_ref = (TH1D*)hsPeaklr_high_ref->Clone();
    hsPeaklr_zeroed_high_ref->Add(minConst_lr2);
    csPeaklr_high_ref->cd();
    hsPeaklr_high_ref->Draw();
    xcoor = 0.52;
    ycoor = 0.90;
    increment = 0.07;
    tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
    tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",0.3,3.0));
    tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
    c->cd();
    Jyieldlr_high_ref.push_back(hsPeaklr_zeroed_high_ref->IntegralAndError(hsPeaklr_zeroed_high_ref->FindBin(0.0),hsPeaklr_zeroed_high_ref->FindBin(minVal_lrX),Jyieldlr_err_high_ref,"width"));
    bin0yield = hsPeaklr_zeroed_high_ref->GetBinContent(hsPeaklr_zeroed_high_ref->FindBin(0.0))*0.19635;
    Jyieldlr_high_ref[0] = Jyieldlr_high_ref[0]*2 - bin0yield;

    JyieldSub_high_ref.push_back(Jyieldsr_high_ref[0] - Jyieldlr_high_ref[0]);
    JyieldSub_err_high_ref.push_back(sqrt(Jyieldsr_err_high_ref*Jyieldsr_err_high_ref + Jyieldlr_err_high_ref*Jyieldlr_err_high_ref));

    TF1* fit2 = new TF1("fit2","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit2->SetParNames("N","V1","V2","V3");
    fit2->SetParameters(10,1,1,1);
    fit2->SetLineColor(2);

    V2lrs_high_ref = hsignalPeak_high_ref->ProjectionY("V2lrs_high_ref",1,14);
    TH1D* aV2lrs_high_ref = hsignalPeak_high_ref->ProjectionY("aV2lrs_high_ref",20,33);
    V2lrb_high_ref = hbackgroundPeak_high_ref->ProjectionY("V2lrb_high_ref",1,14);
    TH1D* aV2lrb_high_ref = hbackgroundPeak_high_ref->ProjectionY("aV2lrb_high_ref",20,33);
    V2lrs_high_ref->Add(aV2lrs_high_ref);
    V2lrb_high_ref->Add(aV2lrb_high_ref);
    V2lrs_high_ref->Divide(V2lrb_high_ref);
    V2lrs_high_ref->Scale(Bz_high_ref[0]/nEvent_high_ref[0]/BW2D);

    V2lrs_high_ref->Fit("fit2","R");
    V2lrs_high_ref->Fit("fit2","R");
    V2lrs_high_ref->Fit("fit2","R");

    V2Values_high_ref.push_back(fit2->GetParameter(2));
    V2Values_err_high_ref.push_back(fit2->GetParError(2));
    V3Values_high_ref.push_back(fit2->GetParameter(3));
    V3Values_err_high_ref.push_back(fit2->GetParError(3));

    double v2_high_ref = sqrt(V2Values_high_ref[0]);
    double v2e_high_ref = sqrt(V2Values_high_ref[0])*(V2Values_err_high_ref[0]/V2Values_high_ref[0])/2;
    double v3_high_ref = sqrt(V3Values_high_ref[0]);
    double v3e_high_ref = sqrt(V3Values_high_ref[0])*(V3Values_err_high_ref[0]/V3Values_high_ref[0])/2;

    Nassoc_high_ref.push_back(fit2->GetParameter(0));
    c->cd();

    std::string PathBackground_low_ks;
    std::string PathSignal_low_ks;
    std::string PathMult_low_ks;
    std::string PathBackground_high_ks;
    std::string PathSignal_high_ks;
    std::string PathMult_high_ks;

    std::string PathBackground_low_la;
    std::string PathSignal_low_la;
    std::string PathMult_low_la;
    std::string PathBackground_high_la;
    std::string PathSignal_high_la;
    std::string PathMult_high_la;

    std::string PathBackground_low_xi;
    std::string PathSignal_low_xi;
    std::string PathMult_low_xi;
    std::string PathBackground_high_xi;
    std::string PathSignal_high_xi;

    int numPtBins_ks = PD.PtBin_ks.size()-1;
    int numPtBins_la = PD.PtBin_la.size()-1;
    int numPtBins_xi = PD.PtBin_xi.size()-1;


    for(int j=0; j<2; j++)
    {
        if(j == 0)
        {
            PathBackground_low_ks  = PD.fn_v0cas + "/backgroundkshort_pt%d";
            PathSignal_low_ks = PD.fn_v0cas + "/signalkshort_pt%d";
            PathMult_low_ks = PD.fn_v0cas + "/mult_ks_pt%d";
            PathBackground_high_ks  = PD.fn_V0 + "/backgroundkshort_pt%d";
            PathSignal_high_ks = PD.fn_V0 + "/signalkshort_pt%d";
            PathMult_high_ks = PD.fn_V0 + "/mult_ks_pt%d";

            PathBackground_low_la  = PD.fn_v0cas + "/backgroundlambda_pt%d";
            PathSignal_low_la = PD.fn_v0cas + "/signallambda_pt%d";
            PathMult_low_la = PD.fn_v0cas + "/mult_la_pt%d";
            PathBackground_high_la  = PD.fn_V0 + "/backgroundlambda_pt%d";
            PathSignal_high_la = PD.fn_V0 + "/signallambda_pt%d";
            PathMult_high_la = PD.fn_V0 + "/mult_la_pt%d";

            PathBackground_low_xi  = PD.fn_v0cas + "/BackgroundXiPeak_pt%d";
            PathSignal_low_xi = PD.fn_v0cas + "/SignalXiPeak_pt%d";
            PathMult_low_xi = PD.fn_v0cas + "/mult_xi_pt%d";
            PathBackground_high_xi  = PD.fn_Xi + "/BackgroundPeak_pt%d";
            PathSignal_high_xi = PD.fn_Xi + "/SignalPeak_pt%d";
        }
        if(j == 1)
        {
            PathBackground_low_ks  = PD.fn_v0cas + "/backgroundkshort_bkg_pt%d";
            PathSignal_low_ks = PD.fn_v0cas + "/signalkshort_bkg_pt%d";
            PathMult_low_ks = PD.fn_v0cas + "/mult_ks_bkg_pt%d";
            PathBackground_high_ks  = PD.fn_V0 + "/backgroundkshort_bkg_pt%d";
            PathSignal_high_ks = PD.fn_V0 + "/signalkshort_bkg_pt%d";
            PathMult_high_ks = PD.fn_V0 + "/mult_ks_bkg_pt%d";

            PathBackground_low_la  = PD.fn_v0cas + "/backgroundlambda_bkg_pt%d";
            PathSignal_low_la = PD.fn_v0cas + "/signallambda_bkg_pt%d";
            PathMult_low_la = PD.fn_v0cas + "/mult_la_bkg_pt%d";
            PathBackground_high_la  = PD.fn_V0 + "/backgroundlambda_bkg_pt%d";
            PathSignal_high_la = PD.fn_V0 + "/signallambda_bkg_pt%d";
            PathMult_high_la = PD.fn_V0 + "/mult_la_bkg_pt%d";

            PathBackground_low_xi  = PD.fn_v0cas + "/BackgroundXiSide_pt%d";
            PathSignal_low_xi = PD.fn_v0cas + "/SignalXiSide_pt%d";
            PathMult_low_xi = PD.fn_v0cas + "/mult_xi_bkg_pt%d";
            PathBackground_high_xi  = PD.fn_Xi + "/BackgroundSide_pt%d";
            PathSignal_high_xi = PD.fn_Xi + "/SignalSide_pt%d";
        }
        //Low N
        //KSHORT
        //Containers
        std::vector<double> Nassoc_low_ks;
        std::vector<double> Jyieldsr_low_ks;
        std::vector<double> Jyieldlr_low_ks;
        std::vector<double> Jyieldsr_err_low_ks(18);
        std::vector<double> Jyieldlr_err_low_ks(18);
        std::vector<double> V2Values_low_ks;
        std::vector<double> V2Values_err_low_ks;
        std::vector<double> V3Values_low_ks;
        std::vector<double> V3Values_err_low_ks;
        std::vector<double> Bz_low_ks;
        std::vector<double> nEvent_low_ks;

        std::vector<double> JyieldSub_low_ks;
        std::vector<double> JyieldSub_err_low_ks;

        int arraySize_ks = PD.PtBin_ks.size();

        TH1D* hsPeaksr_low_ks[arraySize_ks];
        TH1D* hbPeaksr_low_ks[arraySize_ks];
        TH1D* hsPeaklr_low_ks[arraySize_ks];
        TH1D* hbPeaklr_low_ks[arraySize_ks];
        TH1D* V2lrs_low_ks[arraySize_ks];
        TH1D* V2lrb_low_ks[arraySize_ks];
        TH2D* hbackgroundPeak_low_ks[arraySize_ks];
        TH2D* hsignalPeak_low_ks[arraySize_ks];


        //Calculate Nassoc, Jet yield, Peak region Low N

        for(int i=0; i<numPtBins_ks; i++)
        {
            PD.f_perisub->GetObject(Form(PathBackground_low_ks.c_str(),i),hbackgroundPeak_low_ks[i]);
            PD.f_perisub->GetObject(Form(PathSignal_low_ks.c_str(),i),hsignalPeak_low_ks[i]);
            TH1D* mult_low_ks = (TH1D*) PD.f_perisub->Get(Form(PathMult_low_ks.c_str(),i));

            hbPeaksr_low_ks[i] = hbackgroundPeak_low_ks[i]->ProjectionY("hbPeaksr_low_ks", 14, 20);
            hsPeaksr_low_ks[i] = hsignalPeak_low_ks[i]->ProjectionY("hsPeaksr_low_ks", 14, 20);

            TF1* quadFit13 = new TF1("quadFit13","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit13->SetParameters(1,1,1);


            //nEvent_low_ks.push_back(mult_low_ks->Integral(0,100000));
            nEvent_low_ks.push_back(mult_low_ks->Integral(2,100000));
            Bz_low_ks.push_back(hbackgroundPeak_low_ks[i]->GetBinContent(hbackgroundPeak_low_ks[i]->FindBin(0,0)));

            hsPeaksr_low_ks[i]->Divide(hbPeaksr_low_ks[i]);
            hsPeaksr_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);

            hsPeaksr_low_ks[i]->Fit("quadFit13","R");
            hsPeaksr_low_ks[i]->Fit("quadFit13","R");
            hsPeaksr_low_ks[i]->Fit("quadFit13","R");

            double minVal_sr = quadFit13->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit13->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            TH1D* hsPeaksr_zeroed_low_ks = (TH1D*)hsPeaksr_low_ks[i]->Clone();
            hsPeaksr_zeroed_low_ks->Add(minConst_sr);
            csPeaksr_low_ks->cd(i+1);
            hsPeaksr_low_ks[i]->Draw();
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_ks[i],PD.PtBin_ks[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_low_ks.push_back(hsPeaksr_zeroed_low_ks->IntegralAndError(hsPeaksr_zeroed_low_ks->FindBin(0),hsPeaksr_zeroed_low_ks->FindBin(minVal_srX),Jyieldsr_err_low_ks[i],"width"));
            double bin0yield = hsPeaksr_zeroed_low_ks->GetBinContent(hsPeaksr_zeroed_low_ks->FindBin(0.0))*0.19635;
            Jyieldsr_low_ks[i] = Jyieldsr_low_ks[i]*2 - bin0yield;

            hbPeaklr_low_ks[i] = hbackgroundPeak_low_ks[i]->ProjectionY("hbPeaklr_low_ks",1,14);
            TH1D* ahbPeaklr_low_ks = hbackgroundPeak_low_ks[i]->ProjectionY("ahbPeaklr_low_ks",20,33);
            hsPeaklr_low_ks[i] = hsignalPeak_low_ks[i]->ProjectionY("hsPeaklr_low_ks",1,14);
            TH1D* ahsPeaklr_low_ks = hsignalPeak_low_ks[i]->ProjectionY("ahsPeaklr_low_ks",20,33);

            hbPeaklr_low_ks[i]->Add(ahbPeaklr_low_ks);
            hsPeaklr_low_ks[i]->Add(ahsPeaklr_low_ks);
            hsPeaklr_low_ks[i]->Divide(hbPeaklr_low_ks[i]);
            hsPeaklr_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);

            TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit2->SetParameters(1,1,1);
            csPeaklr_low_ks->cd(i+1);

            hsPeaklr_low_ks[i]->Fit("quadFit2","R");
            //hsPeaklr_low_ks[i]->Fit("quadFit2","R");
            //hsPeaklr_low_ks[i]->Fit("quadFit2","R");

            double minVal_lr = quadFit2->GetMinimum(0.6,2.2);
            double minVal_lrX = quadFit2->GetMinimumX(0.6,2.2);
            TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_lr->SetParameter(0,-minVal_lr);
            TH1D* hsPeaklr_zeroed_low_ks = (TH1D*)hsPeaklr_low_ks[i]->Clone();
            hsPeaklr_zeroed_low_ks->Add(minConst_lr);
            hsPeaklr_low_ks[i]->Draw();
            xcoor = 0.52;
            ycoor = 0.90;
            increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_ks[i],PD.PtBin_ks[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldlr_low_ks.push_back(hsPeaklr_zeroed_low_ks->IntegralAndError(hsPeaklr_zeroed_low_ks->FindBin(0.0),hsPeaklr_zeroed_low_ks->FindBin(minVal_lrX),Jyieldlr_err_low_ks[i],"width"));
            bin0yield = hsPeaklr_zeroed_low_ks->GetBinContent(hsPeaklr_zeroed_low_ks->FindBin(0.0))*0.19635;
            Jyieldlr_low_ks[i] = Jyieldlr_low_ks[i]*2 - bin0yield;

            JyieldSub_low_ks.push_back(Jyieldsr_low_ks[i] - Jyieldlr_low_ks[i]);
            JyieldSub_err_low_ks.push_back(sqrt(Jyieldsr_err_low_ks[i]*Jyieldsr_err_low_ks[i] + Jyieldlr_err_low_ks[i]*Jyieldlr_err_low_ks[i]));

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_low_ks[i] = hsignalPeak_low_ks[i]->ProjectionY(Form("V2lrs_low_ks%d",i),1,14);
            TH1D* aV2lrs_low_ks = hsignalPeak_low_ks[i]->ProjectionY("aV2lrs_low_ks",20,33);
            V2lrb_low_ks[i] = hbackgroundPeak_low_ks[i]->ProjectionY(Form("V2lrb_low_ks%d",i),1,14);
            TH1D* aV2lrb_low_ks = hbackgroundPeak_low_ks[i]->ProjectionY("aV2lrb_low_ks",20,33);
            V2lrs_low_ks[i]->Add(aV2lrs_low_ks);
            V2lrb_low_ks[i]->Add(aV2lrb_low_ks);
            V2lrs_low_ks[i]->Divide(V2lrb_low_ks[i]);
            V2lrs_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);

            V2lrs_low_ks[i]->Fit("fit1","R");
            //V2lrs_low_ks[i]->Fit("fit1","R");
            //V2lrs_low_ks[i]->Fit("fit1","R");

            V2Values_low_ks.push_back(fit1->GetParameter(2));
            V2Values_err_low_ks.push_back(fit1->GetParError(2));
            V3Values_low_ks.push_back(fit1->GetParameter(3));
            V3Values_err_low_ks.push_back(fit1->GetParError(3));

            Nassoc_low_ks.push_back(fit1->GetParameter(0));
            c->cd();
        }

        //Lambda N low
        //Containers
        std::vector<double> Nassoc_low_la;
        std::vector<double> Jyieldsr_low_la;
        std::vector<double> Jyieldlr_low_la;
        std::vector<double> Jyieldsr_err_low_la(18);
        std::vector<double> Jyieldlr_err_low_la(18);
        std::vector<double> V2Values_low_la;
        std::vector<double> V2Values_err_low_la;
        std::vector<double> V3Values_low_la;
        std::vector<double> V3Values_err_low_la;
        std::vector<double> Bz_low_la;
        std::vector<double> nEvent_low_la;

        std::vector<double> JyieldSub_low_la;
        std::vector<double> JyieldSub_err_low_la;

        int arraySize_la = PD.PtBin_la.size();

        TH1D* hsPeaksr_low_la[arraySize_la];
        TH1D* hbPeaksr_low_la[arraySize_la];
        TH1D* hsPeaklr_low_la[arraySize_la];
        TH1D* hbPeaklr_low_la[arraySize_la];
        TH1D* V2lrs_low_la[arraySize_la];
        TH1D* V2lrb_low_la[arraySize_la];
        TH2D* hbackgroundPeak_low_la[arraySize_la];
        TH2D* hsignalPeak_low_la[arraySize_la];

        //Calculate Nassoc, Jet yield, Peak region Low N

        for(int i=0; i<numPtBins_la; i++)
        {
            PD.f_perisub->GetObject(Form(PathBackground_low_la.c_str(),i),hbackgroundPeak_low_la[i]);
            PD.f_perisub->GetObject(Form(PathSignal_low_la.c_str(),i),hsignalPeak_low_la[i]);
            TH1D* mult_low_la = (TH1D*) PD.f_perisub->Get(Form(PathMult_low_la.c_str(),i));

            hbPeaksr_low_la[i] = hbackgroundPeak_low_la[i]->ProjectionY("hbPeaksr_low_la", 14, 20);
            hsPeaksr_low_la[i] = hsignalPeak_low_la[i]->ProjectionY("hsPeaksr_low_la", 14, 20);

            TF1* quadFit14 = new TF1("quadFit14","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit14->SetParameters(1,1,1);

            //nEvent_low_la.push_back(mult_low_la->Integral(0,100000));
            nEvent_low_la.push_back(mult_low_la->Integral(2,100000));
            Bz_low_la.push_back(hbackgroundPeak_low_la[i]->GetBinContent(hbackgroundPeak_low_la[i]->FindBin(0,0)));

            hsPeaksr_low_la[i]->Divide(hbPeaksr_low_la[i]);
            hsPeaksr_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);

            csPeaksr_low_la->cd(i+1);

            hsPeaksr_low_la[i]->Fit("quadFit14","R");
            //hsPeaksr_low_la[i]->Fit("quadFit14","R");
            //hsPeaksr_low_la[i]->Fit("quadFit14","R");

            double minVal_sr = quadFit14->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit14->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            TH1D* hsPeaksr_zeroed_low_la = (TH1D*)hsPeaksr_low_la[i]->Clone();
            hsPeaksr_zeroed_low_la->Add(minConst_sr);
            hsPeaksr_low_la[i]->Draw();
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_la[i],PD.PtBin_la[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_low_la.push_back(hsPeaksr_zeroed_low_la->IntegralAndError(hsPeaksr_zeroed_low_la->FindBin(0.0),hsPeaksr_zeroed_low_la->FindBin(minVal_srX),Jyieldsr_err_low_la[i],"width"));
            double bin0yield = hsPeaksr_zeroed_low_la->GetBinContent(hsPeaksr_zeroed_low_la->FindBin(0.0))*0.19635;
            Jyieldsr_low_la[i] = Jyieldsr_low_la[i]*2 - bin0yield;

            hbPeaklr_low_la[i] = hbackgroundPeak_low_la[i]->ProjectionY("hbPeaklr_low_la",1,14);
            TH1D* ahbPeaklr_low_la = hbackgroundPeak_low_la[i]->ProjectionY("ahbPeaklr_low_la",20,33);
            hsPeaklr_low_la[i] = hsignalPeak_low_la[i]->ProjectionY("hsPeaklr_low_la",1,14);
            TH1D* ahsPeaklr_low_la = hsignalPeak_low_la[i]->ProjectionY("ahsPeaklr_low_la",20,33);

            hbPeaklr_low_la[i]->Add(ahbPeaklr_low_la);
            hsPeaklr_low_la[i]->Add(ahsPeaklr_low_la);
            hsPeaklr_low_la[i]->Divide(hbPeaklr_low_la[i]);
            hsPeaklr_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);

            TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit2->SetParameters(1,1,1);

            hsPeaklr_low_la[i]->Fit("quadFit2","R");
            //hsPeaklr_low_la[i]->Fit("quadFit2","R");
            //hsPeaklr_low_la[i]->Fit("quadFit2","R");

            double minVal_lr = quadFit2->GetMinimum(0.6,2.2);
            double minVal_lrX = quadFit2->GetMinimumX(0.6,2.2);
            TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_lr->SetParameter(0,-minVal_lr);
            TH1D* hsPeaklr_zeroed_low_la = (TH1D*)hsPeaklr_low_la[i]->Clone();
            hsPeaklr_zeroed_low_la->Add(minConst_lr);
            csPeaklr_low_la->cd(i+1);
            hsPeaklr_low_la[i]->Draw();
            xcoor = 0.52;
            ycoor = 0.90;
            increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_la[i],PD.PtBin_la[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldlr_low_la.push_back(hsPeaklr_zeroed_low_la->IntegralAndError(hsPeaklr_zeroed_low_la->FindBin(0.0),hsPeaklr_zeroed_low_la->FindBin(minVal_lrX),Jyieldlr_err_low_la[i],"width"));
            bin0yield = hsPeaklr_zeroed_low_la->GetBinContent(hsPeaklr_zeroed_low_la->FindBin(0.0))*0.19635;
            Jyieldlr_low_la[i] = Jyieldlr_low_la[i]*2 - bin0yield;

            JyieldSub_low_la.push_back(Jyieldsr_low_la[i] - Jyieldlr_low_la[i]);
            JyieldSub_err_low_la.push_back(sqrt(Jyieldsr_err_low_la[i]*Jyieldsr_err_low_la[i] + Jyieldlr_err_low_la[i]*Jyieldlr_err_low_la[i]));

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_low_la[i] = hsignalPeak_low_la[i]->ProjectionY(Form("V2lrs_low_la%d",i),1,14);
            TH1D* aV2lrs_low_la = hsignalPeak_low_la[i]->ProjectionY("aV2lrs_low_la",20,33);
            V2lrb_low_la[i] = hbackgroundPeak_low_la[i]->ProjectionY(Form("V2lrb_low_la%d",i),1,14);
            TH1D* aV2lrb_low_la = hbackgroundPeak_low_la[i]->ProjectionY("aV2lrb_low_la",20,33);
            V2lrs_low_la[i]->Add(aV2lrs_low_la);
            V2lrb_low_la[i]->Add(aV2lrb_low_la);
            V2lrs_low_la[i]->Divide(V2lrb_low_la[i]);
            V2lrs_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);

            V2lrs_low_la[i]->Fit("fit1","R");
            //V2lrs_low_la[i]->Fit("fit1","R");
            //V2lrs_low_la[i]->Fit("fit1","R");

            V2Values_low_la.push_back(fit1->GetParameter(2));
            V2Values_err_low_la.push_back(fit1->GetParError(2));
            V3Values_low_la.push_back(fit1->GetParameter(3));
            V3Values_err_low_la.push_back(fit1->GetParError(3));

            Nassoc_low_la.push_back(fit1->GetParameter(0));
            c->cd();
        }

        //Xi Low N
        std::vector<double> Nassoc_low_xi;
        std::vector<double> Jyieldsr_low_xi;
        std::vector<double> Jyieldlr_low_xi;
        std::vector<double> Jyieldsr_err_low_xi(9);
        std::vector<double> Jyieldlr_err_low_xi(9);
        std::vector<double> V2Values_low_xi;
        std::vector<double> V2Values_err_low_xi;
        std::vector<double> V3Values_low_xi;
        std::vector<double> V3Values_err_low_xi;
        std::vector<double> Bz_low_xi;
        std::vector<double> nEvent_low_xi;

        std::vector<double> JyieldSub_low_xi;
        std::vector<double> JyieldSub_err_low_xi;

        int arraySize_xi = PD.PtBin_xi.size();

        TH1D* hsPeaksr_low_xi[arraySize_xi];
        TH1D* hbPeaksr_low_xi[arraySize_xi];
        TH1D* hsPeaklr_low_xi[arraySize_xi];
        TH1D* hbPeaklr_low_xi[arraySize_xi];
        TH1D* V2lrs_low_xi[arraySize_xi];
        TH1D* V2lrb_low_xi[arraySize_xi];
        TH2D* hbackgroundPeak_low_xi[arraySize_xi];
        TH2D* hsignalPeak_low_xi[arraySize_xi];


        //Calculate Nassoc, Jet yield, Peak region Low N

        for(int i=0; i<numPtBins_xi; i++)
        {
            PD.f_perisub->GetObject(Form(PathBackground_low_xi.c_str(),i),hbackgroundPeak_low_xi[i]);
            PD.f_perisub->GetObject(Form(PathSignal_low_xi.c_str(),i),hsignalPeak_low_xi[i]);
            TH1D* mult_low_xi = (TH1D*) PD.f_perisub->Get(Form(PathMult_low_xi.c_str(),i));

            hbPeaksr_low_xi[i] = hbackgroundPeak_low_xi[i]->ProjectionY("hbPeaksr_low_xi", 14, 20);
            hsPeaksr_low_xi[i] = hsignalPeak_low_xi[i]->ProjectionY("hsPeaksr_low_xi", 14, 20);

            TF1* quadFit13 = new TF1("quadFit13","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit13->SetParameters(1,1,1);


            //nEvent_low_xi.push_back(mult_low_xi->Integral(0,100000));
            //nEvent_low_xi.push_back(mult_low_xi->Integral(2,100000));
            Bz_low_xi.push_back(hbackgroundPeak_low_xi[i]->GetBinContent(hbackgroundPeak_low_xi[i]->FindBin(0,0)));

            hsPeaksr_low_xi[i]->Divide(hbPeaksr_low_xi[i]);
            //hsPeaksr_low_xi[i]->Scale(Bz_low_xi[i]/nEvent_low_xi[i]/BW2D);
            hsPeaksr_low_xi[i]->Scale(Bz_low_xi[i]/BW2D);

            csPeaksr_low_xi->cd(i+1);

            hsPeaksr_low_xi[i]->Fit("quadFit13","R");
            //hsPeaksr_low_xi[i]->Fit("quadFit13","R");
            //hsPeaksr_low_xi[i]->Fit("quadFit13","R");

            double minVal_sr = quadFit13->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit13->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            TH1D* hsPeaksr_zeroed_low_xi = (TH1D*)hsPeaksr_low_xi[i]->Clone();
            hsPeaksr_zeroed_low_xi->Add(minConst_sr);
            hsPeaksr_low_xi[i]->Draw();
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_xi[i],PD.PtBin_xi[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_low_xi.push_back(hsPeaksr_zeroed_low_xi->IntegralAndError(hsPeaksr_zeroed_low_xi->FindBin(0),hsPeaksr_zeroed_low_xi->FindBin(minVal_srX),Jyieldsr_err_low_xi[i],"width"));
            double bin0yield = hsPeaksr_zeroed_low_xi->GetBinContent(hsPeaksr_zeroed_low_xi->FindBin(0.0))*0.19635;
            Jyieldsr_low_xi[i] = Jyieldsr_low_xi[i]*2 - bin0yield;

            hbPeaklr_low_xi[i] = hbackgroundPeak_low_xi[i]->ProjectionY("hbPeaklr_low_xi",1,14);
            TH1D* ahbPeaklr_low_xi = hbackgroundPeak_low_xi[i]->ProjectionY("ahbPeaklr_low_xi",20,33);
            hsPeaklr_low_xi[i] = hsignalPeak_low_xi[i]->ProjectionY("hsPeaklr_low_xi",1,14);
            TH1D* ahsPeaklr_low_xi = hsignalPeak_low_xi[i]->ProjectionY("ahsPeaklr_low_xi",20,33);

            hbPeaklr_low_xi[i]->Add(ahbPeaklr_low_xi);
            hsPeaklr_low_xi[i]->Add(ahsPeaklr_low_xi);
            hsPeaklr_low_xi[i]->Divide(hbPeaklr_low_xi[i]);
            //hsPeaklr_low_xi[i]->Scale(Bz_low_xi[i]/BW2D);

            csPeaklr_low_xi->cd(i+1);

            TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit2->SetParameters(1,1,1);

            hsPeaklr_low_xi[i]->Fit("quadFit2","R");
            //hsPeaklr_low_xi[i]->Fit("quadFit2","R");
            //hsPeaklr_low_xi[i]->Fit("quadFit2","R");

            double minVal_lr = quadFit2->GetMinimum(0.6,2.2);
            double minVal_lrX = quadFit2->GetMinimumX(0.6,2.2);
            TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_lr->SetParameter(0,-minVal_lr);
            TH1D* hsPeaklr_zeroed_low_xi = (TH1D*)hsPeaklr_low_xi[i]->Clone();
            hsPeaklr_zeroed_low_xi->Add(minConst_lr);
            hsPeaklr_low_xi[i]->Draw();
            xcoor = 0.52;
            ycoor = 0.90;
            increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_xi[i],PD.PtBin_xi[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldlr_low_xi.push_back(hsPeaklr_zeroed_low_xi->IntegralAndError(hsPeaklr_zeroed_low_xi->FindBin(0.0),hsPeaklr_zeroed_low_xi->FindBin(minVal_lrX),Jyieldlr_err_low_xi[i],"width"));
            bin0yield = hsPeaklr_zeroed_low_xi->GetBinContent(hsPeaklr_zeroed_low_xi->FindBin(0.0))*0.19635;
            Jyieldlr_low_xi[i] = Jyieldlr_low_xi[i]*2 - bin0yield;

            JyieldSub_low_xi.push_back(Jyieldsr_low_xi[i] - Jyieldlr_low_xi[i]);
            JyieldSub_err_low_xi.push_back(sqrt(Jyieldsr_err_low_xi[i]*Jyieldsr_err_low_xi[i] + Jyieldlr_err_low_xi[i]*Jyieldlr_err_low_xi[i]));

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_low_xi[i] = hsignalPeak_low_xi[i]->ProjectionY(Form("V2lrs_low_xi%d",i),1,14);
            TH1D* aV2lrs_low_xi = hsignalPeak_low_xi[i]->ProjectionY("aV2lrs_low_xi",20,33);
            V2lrb_low_xi[i] = hbackgroundPeak_low_xi[i]->ProjectionY(Form("V2lrb_low_xi%d",i),1,14);
            TH1D* aV2lrb_low_xi = hbackgroundPeak_low_xi[i]->ProjectionY("aV2lrb_low_xi",20,33);
            V2lrs_low_xi[i]->Add(aV2lrs_low_xi);
            V2lrb_low_xi[i]->Add(aV2lrb_low_xi);
            V2lrs_low_xi[i]->Divide(V2lrb_low_xi[i]);
            //V2lrs_low_xi[i]->Scale(Bz_low_xi[i]/nEvent_low_xi[i]/BW2D);
            V2lrs_low_xi[i]->Scale(Bz_low_xi[i]/BW2D);

            V2lrs_low_xi[i]->Fit("fit1","R");
            //V2lrs_low_xi[i]->Fit("fit1","R");
            //V2lrs_low_xi[i]->Fit("fit1","R");

            V2Values_low_xi.push_back(fit1->GetParameter(2));
            V2Values_err_low_xi.push_back(fit1->GetParError(2));
            V3Values_low_xi.push_back(fit1->GetParameter(3));
            V3Values_err_low_xi.push_back(fit1->GetParError(3));

            Nassoc_low_xi.push_back(fit1->GetParameter(0));
            c->cd();
        }

        //High N
        //KSHORT
        //Containers
        std::vector<double> Nassoc_ks;
        std::vector<double> Jyieldsr_ks;
        std::vector<double> Jyieldlr_ks;
        std::vector<double> Jyieldsr_err_ks(18);
        std::vector<double> Jyieldlr_err_ks(18);
        std::vector<double> V2Values_ks;
        std::vector<double> V2Values_err_ks;
        std::vector<double> V3Values_ks;
        std::vector<double> V3Values_err_ks;
        std::vector<double> Bz_ks;
        std::vector<double> nEvent_ks;

        std::vector<double> JyieldSub_ks;
        std::vector<double> JyieldSub_err_ks;

        TH1D* hsPeaksr_ks[arraySize_ks];
        TH1D* hbPeaksr_ks[arraySize_ks];
        TH1D* hsPeaklr_ks[arraySize_ks];
        TH1D* hbPeaklr_ks[arraySize_ks];
        TH1D* V2lrs_ks[arraySize_ks];
        TH1D* V2lrb_ks[arraySize_ks];
        TH2D* hbackgroundPeak_ks[arraySize_ks];
        TH2D* hsignalPeak_ks[arraySize_ks];

        //Calculate Nassoc, Jet yield, Peak region
        //Jet Yield

        for(int i=0; i<numPtBins_ks; i++)
        {
            PD.f_V0->GetObject(Form(PathBackground_high_ks.c_str(),i),hbackgroundPeak_ks[i]);
            PD.f_V0->GetObject(Form(PathSignal_high_ks.c_str(),i),hsignalPeak_ks[i]);
            TH1D* mult_ks = (TH1D*) PD.f_V0->Get(Form(PathMult_high_ks.c_str(),i));

            hbPeaksr_ks[i] = hbackgroundPeak_ks[i]->ProjectionY("hbPeaksr_ks", 14, 20);
            hsPeaksr_ks[i] = hsignalPeak_ks[i]->ProjectionY("hsPeaksr_ks", 14, 20);

            TF1* quadFit15 = new TF1("quadFit15","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit15->SetParameters(1,1,1);

            //nEvent_ks.push_back(mult_ks->Integral(0,100000));
            nEvent_ks.push_back(mult_ks->Integral(2,100000));
            Bz_ks.push_back(hbackgroundPeak_ks[i]->GetBinContent(hbackgroundPeak_ks[i]->FindBin(0,0)));

            hsPeaksr_ks[i]->Divide(hbPeaksr_ks[i]);
            hsPeaksr_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);

            csPeaksr_ks->cd(i+1);

            hsPeaksr_ks[i]->Fit("quadFit15","R");
            //hsPeaksr_ks[i]->Fit("quadFit15","R");
            //hsPeaksr_ks[i]->Fit("quadFit15","R");

            double minVal_sr = quadFit15->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit15->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            TH1D* hsPeaksr_zeroed_ks = (TH1D*)hsPeaksr_ks[i]->Clone();
            hsPeaksr_zeroed_ks->Add(minConst_sr);
            hsPeaksr_ks[i]->Draw();
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_ks[i],PD.PtBin_ks[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_ks.push_back(hsPeaksr_zeroed_ks->IntegralAndError(hsPeaksr_zeroed_ks->FindBin(0.0),hsPeaksr_zeroed_ks->FindBin(minVal_srX),Jyieldsr_err_ks[i],"width"));
            double bin0yield = hsPeaksr_zeroed_ks->GetBinContent(hsPeaksr_zeroed_ks->FindBin(0.0))*0.19635;
            Jyieldsr_ks[i] = Jyieldsr_ks[i]*2 - bin0yield;

            hbPeaklr_ks[i] = hbackgroundPeak_ks[i]->ProjectionY("hbPeaklr_ks",1,14);
            TH1D* ahbPeaklr_ks = hbackgroundPeak_ks[i]->ProjectionY("ahbPeaklr_ks",20,33);
            hsPeaklr_ks[i] = hsignalPeak_ks[i]->ProjectionY("hsPeaklr_ks",1,14);
            TH1D* ahsPeaklr_ks = hsignalPeak_ks[i]->ProjectionY("ahsPeaklr_ks",20,33);

            hbPeaklr_ks[i]->Add(ahbPeaklr_ks);
            hsPeaklr_ks[i]->Add(ahsPeaklr_ks);
            hsPeaklr_ks[i]->Divide(hbPeaklr_ks[i]);
            hsPeaklr_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);

            csPeaklr_ks->cd(i+1);

            TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit2->SetParameters(1,1,1);

            hsPeaklr_ks[i]->Fit("quadFit2","R");
            //hsPeaklr_ks[i]->Fit("quadFit2","R");
            //hsPeaklr_ks[i]->Fit("quadFit2","R");

            double minVal_lr = quadFit2->GetMinimum(0.6,2.2);
            double minVal_lrX = quadFit2->GetMinimumX(0.6,2.2);
            TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_lr->SetParameter(0,-minVal_lr);
            TH1D* hsPeaklr_zeroed_ks = (TH1D*)hsPeaklr_ks[i]->Clone();
            hsPeaklr_zeroed_ks->Add(minConst_lr);
            hsPeaklr_ks[i]->Draw();
            xcoor = 0.52;
            ycoor = 0.90;
            increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_ks[i],PD.PtBin_ks[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldlr_ks.push_back(hsPeaklr_zeroed_ks->IntegralAndError(hsPeaklr_zeroed_ks->FindBin(0.0),hsPeaklr_zeroed_ks->FindBin(minVal_lrX),Jyieldlr_err_ks[i],"width"));
            bin0yield = hsPeaklr_zeroed_ks->GetBinContent(hsPeaklr_zeroed_ks->FindBin(0.0))*0.19635;
            Jyieldlr_ks[i] = Jyieldlr_ks[i]*2 - bin0yield;

            JyieldSub_ks.push_back(Jyieldsr_ks[i] - Jyieldlr_ks[i]);
            JyieldSub_err_ks.push_back(sqrt(Jyieldsr_err_ks[i]*Jyieldsr_err_ks[i] + Jyieldlr_err_ks[i]*Jyieldlr_err_ks[i]));

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_ks[i] = hsignalPeak_ks[i]->ProjectionY(Form("V2lrs_ks%d",i),1,14);
            TH1D* aV2lrs_ks = hsignalPeak_ks[i]->ProjectionY("aV2lrs_ks",20,33);
            V2lrb_ks[i] = hbackgroundPeak_ks[i]->ProjectionY(Form("V2lrb_ks%d",i),1,14);
            TH1D* aV2lrb_ks = hbackgroundPeak_ks[i]->ProjectionY("aV2lrb_ks",20,33);
            V2lrs_ks[i]->Add(aV2lrs_ks);
            V2lrb_ks[i]->Add(aV2lrb_ks);
            V2lrs_ks[i]->Divide(V2lrb_ks[i]);
            V2lrs_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);

            V2lrs_ks[i]->Fit("fit1","R");
            //V2lrs_ks[i]->Fit("fit1","R");
            //V2lrs_ks[i]->Fit("fit1","R");

            V2Values_ks.push_back(fit1->GetParameter(2));
            V2Values_err_ks.push_back(fit1->GetParError(2));
            V3Values_ks.push_back(fit1->GetParameter(3));
            V3Values_err_ks.push_back(fit1->GetParError(3));

            Nassoc_ks.push_back(fit1->GetParameter(0));
            c->cd();
        }

        //LAMBDA
        std::vector<double> Nassoc_la;
        std::vector<double> Jyieldsr_la;
        std::vector<double> Jyieldlr_la;
        std::vector<double> Jyieldsr_err_la(18);
        std::vector<double> Jyieldlr_err_la(18);
        std::vector<double> V2Values_la;
        std::vector<double> V2Values_err_la;
        std::vector<double> V3Values_la;
        std::vector<double> V3Values_err_la;
        std::vector<double> Bz_la;
        std::vector<double> nEvent_la;

        std::vector<double> JyieldSub_la;
        std::vector<double> JyieldSub_err_la;

        TH1D* hsPeaksr_la[arraySize_la];
        TH1D* hbPeaksr_la[arraySize_la];
        TH1D* hsPeaklr_la[arraySize_la];
        TH1D* hbPeaklr_la[arraySize_la];
        TH1D* V2lrs_la[arraySize_la];
        TH1D* V2lrb_la[arraySize_la];
        TH2D* hbackgroundPeak_la[arraySize_la];
        TH2D* hsignalPeak_la[arraySize_la];

        //Calculate Nassoc, Jet yield, Peak region
        //Jet Yield

        for(int i=0; i<numPtBins_la; i++)
        {
            PD.f_V0->GetObject(Form(PathBackground_high_la.c_str(),i),hbackgroundPeak_la[i]);
            PD.f_V0->GetObject(Form(PathSignal_high_la.c_str(),i),hsignalPeak_la[i]);
            //hbackgroundPeak_la[i] = (TH2D*)PD.f_V0->Get(Form((fn_V0 + "/backgroundlambda_pt%d").c_str(),i));
            //hsignalPeak_la[i] = (TH2D*)f_V0->Get(Form((fn_V0 + "/signallambda_pt%d").c_str(),i));
            TH1D* mult_la = (TH1D*) PD.f_V0->Get(Form(PathMult_high_la.c_str(),i));

            hbPeaksr_la[i] = hbackgroundPeak_la[i]->ProjectionY("hbPeaksr_la", 14, 20);
            hsPeaksr_la[i] = hsignalPeak_la[i]->ProjectionY("hsPeaksr_la", 14, 20);

            TF1* quadFit16 = new TF1("quadFit16","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit16->SetParameters(1,1,1);

            //nEvent_la.push_back(mult_la->Integral(0,10000));
            nEvent_la.push_back(mult_la->Integral(2,10000));
            Bz_la.push_back(hbackgroundPeak_la[i]->GetBinContent(hbackgroundPeak_la[i]->FindBin(0,0)));

            hsPeaksr_la[i]->Divide(hbPeaksr_la[i]);
            hsPeaksr_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);

            csPeaksr_la->cd(i+1);

            TH1D* hsPeaksr_zeroed_la = (TH1D*)hsPeaksr_la[i]->Clone();

            hsPeaksr_la[i]->Fit("quadFit16","R");
            //hsPeaksr_la[i]->Fit("quadFit16","R");
            //hsPeaksr_la[i]->Fit("quadFit16","R");

            double minVal_sr = quadFit16->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit16->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            hsPeaksr_zeroed_la->Add(minConst_sr);
            hsPeaksr_la[i]->Draw();
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%4.1f<p_{T}^{trg}<%4.1f",PD.PtBin_la[i],PD.PtBin_la[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_la.push_back(hsPeaksr_zeroed_la->IntegralAndError(hsPeaksr_zeroed_la->FindBin(0.0),hsPeaksr_zeroed_la->FindBin(minVal_srX),Jyieldsr_err_la[i],"width"));
            double bin0yield = hsPeaksr_zeroed_la->GetBinContent(hsPeaksr_zeroed_la->FindBin(0.0))*0.19635;
            Jyieldsr_la[i] = Jyieldsr_la[i]*2 - bin0yield;

            hbPeaklr_la[i] = hbackgroundPeak_la[i]->ProjectionY("hbPeaklr_la",1,14);
            TH1D* ahbPeaklr_la = hbackgroundPeak_la[i]->ProjectionY("ahbPeaklr_la",20,33);
            hsPeaklr_la[i] = hsignalPeak_la[i]->ProjectionY("hsPeaklr_la",1,14);
            TH1D* ahsPeaklr_la = hsignalPeak_la[i]->ProjectionY("ahsPeaklr_la",20,33);

            hbPeaklr_la[i]->Add(ahbPeaklr_la);
            hsPeaklr_la[i]->Add(ahsPeaklr_la);
            hsPeaklr_la[i]->Divide(hbPeaklr_la[i]);
            hsPeaklr_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);

            csPeaklr_la->cd(i+1);

            TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit2->SetParameters(1,1,1);

            hsPeaklr_la[i]->Fit("quadFit2","R");
            //hsPeaklr_la[i]->Fit("quadFit2","R");
            //hsPeaklr_la[i]->Fit("quadFit2","R");

            double minVal_lr = quadFit2->GetMinimum(0.6,2.2);
            double minVal_lrX = quadFit2->GetMinimumX(0.6,2.2);
            TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_lr->SetParameter(0,-minVal_lr);
            TH1D* hsPeaklr_zeroed_la = (TH1D*)hsPeaklr_la[i]->Clone();
            hsPeaklr_zeroed_la->Add(minConst_lr);
            hsPeaklr_la[i]->Draw();
            xcoor = 0.52;
            ycoor = 0.90;
            increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_la[i],PD.PtBin_la[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldlr_la.push_back(hsPeaklr_zeroed_la->IntegralAndError(hsPeaklr_zeroed_la->FindBin(0.0),hsPeaklr_zeroed_la->FindBin(minVal_lrX),Jyieldlr_err_la[i],"width"));
            bin0yield = hsPeaklr_zeroed_la->GetBinContent(hsPeaklr_zeroed_la->FindBin(0.0))*0.19635;
            Jyieldlr_la[i] = Jyieldlr_la[i]*2 - bin0yield;

            JyieldSub_la.push_back(Jyieldsr_la[i] - Jyieldlr_la[i]);
            JyieldSub_err_la.push_back(sqrt(Jyieldsr_err_la[i]*Jyieldsr_err_la[i] + Jyieldlr_err_la[i]*Jyieldlr_err_la[i]));

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_la[i] = hsignalPeak_la[i]->ProjectionY(Form("V2lrs_la%d",i),1,14);
            TH1D* aV2lrs_la = hsignalPeak_la[i]->ProjectionY("aV2lrs_la",20,33);
            V2lrb_la[i] = hbackgroundPeak_la[i]->ProjectionY(Form("V2lrb_la%d",i),1,14);
            TH1D* aV2lrb_la = hbackgroundPeak_la[i]->ProjectionY("aV2lrb_la",20,33);
            V2lrs_la[i]->Add(aV2lrs_la);
            V2lrb_la[i]->Add(aV2lrb_la);
            V2lrs_la[i]->Divide(V2lrb_la[i]);
            V2lrs_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);

            V2lrs_la[i]->Fit("fit1","R");
            //V2lrs_la[i]->Fit("fit1","R");
            //V2lrs_la[i]->Fit("fit1","R");

            V2Values_la.push_back(fit1->GetParameter(2));
            V2Values_err_la.push_back(fit1->GetParError(2));
            V3Values_la.push_back(fit1->GetParameter(3));
            V3Values_err_la.push_back(fit1->GetParError(3));

            Nassoc_la.push_back(fit1->GetParameter(0));
            c->cd();
        }

        //Xi High N
        std::vector<double> Nassoc_xi;
        std::vector<double> Jyieldsr_xi;
        std::vector<double> Jyieldlr_xi;
        std::vector<double> Jyieldsr_err_xi(18);
        std::vector<double> Jyieldlr_err_xi(18);
        std::vector<double> V2Values_xi;
        std::vector<double> V2Values_err_xi;
        std::vector<double> V3Values_xi;
        std::vector<double> V3Values_err_xi;
        std::vector<double> Bz_xi;
        std::vector<double> nEvent_xi;

        std::vector<double> JyieldSub_xi;
        std::vector<double> JyieldSub_err_xi;

        TH1D* hsPeaksr_xi[arraySize_xi];
        TH1D* hbPeaksr_xi[arraySize_xi];
        TH1D* hsPeaklr_xi[arraySize_xi];
        TH1D* hbPeaklr_xi[arraySize_xi];
        TH1D* V2lrs_xi[arraySize_xi];
        TH1D* V2lrb_xi[arraySize_xi];
        TH2D* hbackgroundPeak_xi[arraySize_xi];
        TH2D* hsignalPeak_xi[arraySize_xi];

        //Calculate Nassoc, Jet yield, Peak region
        //Jet Yield

        for(int i=0; i<numPtBins_xi; i++)
        {
            PD.f_Xi->GetObject(Form(PathBackground_high_xi.c_str(),i),hbackgroundPeak_xi[i]);
            PD.f_Xi->GetObject(Form(PathSignal_high_xi.c_str(),i),hsignalPeak_xi[i]);
            //hbackgroundPeak_xi[i] = (TH2D*)PD.f_V0->Get(Form((fn_V0 + "/backgroundlambda_pt%d").c_str(),i));
            //hsignalPeak_xi[i] = (TH2D*)f_V0->Get(Form((fn_V0 + "/signallambda_pt%d").c_str(),i));
            //TH1D* mult_xi = (TH1D*) PD.f_Xi->Get(Form((PD.fn_Xi + "/mult_xi_pt%d").c_str(),i));

            hbPeaksr_xi[i] = hbackgroundPeak_xi[i]->ProjectionY("hbPeaksr_xi", 14, 20);
            hsPeaksr_xi[i] = hsignalPeak_xi[i]->ProjectionY("hsPeaksr_xi", 14, 20);

            TF1* quadFit16 = new TF1("quadFit16","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit16->SetParameters(1,1,1);

            //nEvent_xi.push_back(mult_xi->Integral(0,10000));
            //nEvent_xi.push_back(mult_xi->Integral(2,10000));
            Bz_xi.push_back(hbackgroundPeak_xi[i]->GetBinContent(hbackgroundPeak_xi[i]->FindBin(0,0)));

            hsPeaksr_xi[i]->Divide(hbPeaksr_xi[i]);
            //hsPeaksr_xi[i]->Scale(Bz_xi[i]/nEvent_xi[i]/BW2D);
            hsPeaksr_xi[i]->Scale(Bz_xi[i]/BW2D);

            csPeaksr_xi->cd(i+1);

            TH1D* hsPeaksr_zeroed_xi = (TH1D*)hsPeaksr_xi[i]->Clone();

            hsPeaksr_xi[i]->Fit("quadFit16","R");
            //hsPeaksr_xi[i]->Fit("quadFit16","R");
            //hsPeaksr_xi[i]->Fit("quadFit16","R");

            double minVal_sr = quadFit16->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit16->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            hsPeaksr_zeroed_xi->Add(minConst_sr);
            hsPeaksr_xi[i]->Draw();
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%4.1f<p_{T}^{trg}<%4.1f",PD.PtBin_xi[i],PD.PtBin_xi[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_xi.push_back(hsPeaksr_zeroed_xi->IntegralAndError(hsPeaksr_zeroed_xi->FindBin(0.0),hsPeaksr_zeroed_xi->FindBin(minVal_srX),Jyieldsr_err_xi[i],"width"));
            double bin0yield = hsPeaksr_zeroed_xi->GetBinContent(hsPeaksr_zeroed_xi->FindBin(0.0))*0.19635;
            Jyieldsr_xi[i] = Jyieldsr_xi[i]*2 - bin0yield;

            hbPeaklr_xi[i] = hbackgroundPeak_xi[i]->ProjectionY("hbPeaklr_xi",1,14);
            TH1D* ahbPeaklr_xi = hbackgroundPeak_xi[i]->ProjectionY("ahbPeaklr_xi",20,33);
            hsPeaklr_xi[i] = hsignalPeak_xi[i]->ProjectionY("hsPeaklr_xi",1,14);
            TH1D* ahsPeaklr_xi = hsignalPeak_xi[i]->ProjectionY("ahsPeaklr_xi",20,33);

            hbPeaklr_xi[i]->Add(ahbPeaklr_xi);
            hsPeaklr_xi[i]->Add(ahsPeaklr_xi);
            hsPeaklr_xi[i]->Divide(hbPeaklr_xi[i]);
            //hsPeaklr_xi[i]->Scale(Bz_xi[i]/nEvent_xi[i]/BW2D);
            hsPeaklr_xi[i]->Scale(Bz_xi[i]/BW2D);

            csPeaklr_xi->cd(i+1);

            TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit2->SetParameters(1,1,1);

            hsPeaklr_xi[i]->Fit("quadFit2","R");
            //hsPeaklr_xi[i]->Fit("quadFit2","R");
            //hsPeaklr_xi[i]->Fit("quadFit2","R");

            double minVal_lr = quadFit2->GetMinimum(0.6,2.2);
            double minVal_lrX = quadFit2->GetMinimumX(0.6,2.2);
            TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_lr->SetParameter(0,-minVal_lr);
            TH1D* hsPeaklr_zeroed_xi = (TH1D*)hsPeaklr_xi[i]->Clone();
            hsPeaklr_zeroed_xi->Add(minConst_lr);
            hsPeaklr_xi[i]->Draw();
            xcoor = 0.52;
            ycoor = 0.90;
            increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_xi[i],PD.PtBin_xi[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldlr_xi.push_back(hsPeaklr_zeroed_xi->IntegralAndError(hsPeaklr_zeroed_xi->FindBin(0.0),hsPeaklr_zeroed_xi->FindBin(minVal_lrX),Jyieldlr_err_xi[i],"width"));
            bin0yield = hsPeaklr_zeroed_xi->GetBinContent(hsPeaklr_zeroed_xi->FindBin(0.0))*0.19635;
            Jyieldlr_xi[i] = Jyieldlr_xi[i]*2 - bin0yield;

            JyieldSub_xi.push_back(Jyieldsr_xi[i] - Jyieldlr_xi[i]);
            JyieldSub_err_xi.push_back(sqrt(Jyieldsr_err_xi[i]*Jyieldsr_err_xi[i] + Jyieldlr_err_xi[i]*Jyieldlr_err_xi[i]));

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_xi[i] = hsignalPeak_xi[i]->ProjectionY(Form("V2lrs_xi%d",i),1,14);
            TH1D* aV2lrs_xi = hsignalPeak_xi[i]->ProjectionY("aV2lrs_xi",20,33);
            V2lrb_xi[i] = hbackgroundPeak_xi[i]->ProjectionY(Form("V2lrb_xi%d",i),1,14);
            TH1D* aV2lrb_xi = hbackgroundPeak_xi[i]->ProjectionY("aV2lrb_xi",20,33);
            V2lrs_xi[i]->Add(aV2lrs_xi);
            V2lrb_xi[i]->Add(aV2lrb_xi);
            V2lrs_xi[i]->Divide(V2lrb_xi[i]);
            //V2lrs_xi[i]->Scale(Bz_xi[i]/nEvent_xi[i]/BW2D);
            V2lrs_xi[i]->Scale(Bz_xi[i]/BW2D);

            V2lrs_xi[i]->Fit("fit1","R");
            //V2lrs_xi[i]->Fit("fit1","R");
            //V2lrs_xi[i]->Fit("fit1","R");

            V2Values_xi.push_back(fit1->GetParameter(2));
            V2Values_err_xi.push_back(fit1->GetParError(2));
            V3Values_xi.push_back(fit1->GetParameter(3));
            V3Values_err_xi.push_back(fit1->GetParError(3));

            Nassoc_xi.push_back(fit1->GetParameter(0));
            c->cd();
        }

        // Reference Vn after corrections
        double V2sub_ref = V2Values_high_ref[0] - V2Values_low_ref[0]*Nassoc_low_ref[0]/Nassoc_high_ref[0]*Jyieldsr_high_ref[0]/Jyieldsr_low_ref[0];
        double V2sube_ref = sqrt(TMath::Power(V2Values_err_high_ref[0],2) + TMath::Power(V2Values_err_low_ref[0]*Nassoc_low_ref[0]/Nassoc_high_ref[0]*Jyieldsr_high_ref[0]/Jyieldsr_low_ref[0],2));

        double V3sub_ref = V3Values_high_ref[0] - V3Values_low_ref[0]*Nassoc_low_ref[0]/Nassoc_high_ref[0]*Jyieldsr_high_ref[0]/Jyieldsr_low_ref[0];
        double V3sube_ref = sqrt(TMath::Power(V3Values_err_high_ref[0],2) + TMath::Power(V3Values_err_low_ref[0]*Nassoc_low_ref[0]/Nassoc_high_ref[0]*Jyieldsr_high_ref[0]/Jyieldsr_low_ref[0],2));


        double v2sub_ref = sqrt(V2sub_ref);
        double v2sube_ref = sqrt(V2sub_ref)*(V2sube_ref/V2sub_ref)/2;

        double v3sub_ref = sqrt(V3sub_ref);
        double v3sube_ref = sqrt(V3sub_ref)*(V3sube_ref/V3sub_ref)/2;



        //Ks Vn correction
        std::vector<double> V2sub_ks;
        std::vector<double> V2sube_ks;

        std::vector<double> V3sub_ks;
        std::vector<double> V3sube_ks;

        std::vector<double> jetYfactor_ks;

        for(unsigned i=0; i<numPtBins_ks; i++){
            V2sub_ks.push_back(V2Values_ks[i] - V2Values_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*JyieldSub_ks[i]/JyieldSub_low_ks[i]);
            V2sube_ks.push_back(sqrt(TMath::Power(V2Values_err_ks[i],2) + sqrt(TMath::Power(V2Values_err_low_ks[i]/V2Values_low_ks[i],2) + TMath::Power(JyieldSub_err_ks[i]/JyieldSub_ks[i],2) + TMath::Power(JyieldSub_err_low_ks[i]/JyieldSub_low_ks[i],2))*V2Values_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*JyieldSub_ks[i]/JyieldSub_low_ks[i]*sqrt(TMath::Power(V2Values_err_low_ks[i]/V2Values_low_ks[i],2) + TMath::Power(JyieldSub_err_ks[i]/JyieldSub_ks[i],2) + TMath::Power(JyieldSub_err_low_ks[i]/JyieldSub_low_ks[i],2))*V2Values_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*JyieldSub_ks[i]/JyieldSub_low_ks[i]));

            V3sub_ks.push_back(V3Values_ks[i] - V3Values_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*JyieldSub_ks[i]/JyieldSub_low_ks[i]);
            V3sube_ks.push_back(sqrt(TMath::Power(V3Values_err_ks[i],2) + sqrt(TMath::Power(V3Values_err_low_ks[i]/V3Values_low_ks[i],2) + TMath::Power(JyieldSub_err_ks[i]/JyieldSub_ks[i],2) + TMath::Power(JyieldSub_err_low_ks[i]/JyieldSub_low_ks[i],2))*V3Values_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*JyieldSub_ks[i]/JyieldSub_low_ks[i]*sqrt(TMath::Power(V3Values_err_low_ks[i]/V3Values_low_ks[i],2) + TMath::Power(JyieldSub_err_ks[i]/JyieldSub_ks[i],2) + TMath::Power(JyieldSub_err_low_ks[i]/JyieldSub_low_ks[i],2))*V3Values_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*JyieldSub_ks[i]/JyieldSub_low_ks[i]));

            jetYfactor_ks.push_back(JyieldSub_ks[i]/JyieldSub_low_ks[i]);
        }

        //La Vn correction
        std::vector<double> V2sub_la;
        std::vector<double> V2sube_la;

        std::vector<double> V3sub_la;
        std::vector<double> V3sube_la;

        std::vector<double> jetYfactor_la;

        for(unsigned i=0; i<numPtBins_la; i++){
            V2sub_la.push_back(V2Values_la[i] - V2Values_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*JyieldSub_la[i]/JyieldSub_low_la[i]);
            V2sube_la.push_back(sqrt(TMath::Power(V2Values_err_la[i],2) + sqrt(TMath::Power(V2Values_err_low_la[i]/V2Values_low_la[i],2) + TMath::Power(JyieldSub_err_la[i]/JyieldSub_la[i],2) + TMath::Power(JyieldSub_err_low_la[i]/JyieldSub_low_la[i],2))*V2Values_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*JyieldSub_la[i]/JyieldSub_low_la[i]*sqrt(TMath::Power(V2Values_err_low_la[i]/V2Values_low_la[i],2) + TMath::Power(JyieldSub_err_la[i]/JyieldSub_la[i],2) + TMath::Power(JyieldSub_err_low_la[i]/JyieldSub_low_la[i],2))*V2Values_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*JyieldSub_la[i]/JyieldSub_low_la[i]));

            V3sub_la.push_back(V3Values_la[i] - V3Values_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*JyieldSub_la[i]/JyieldSub_low_la[i]);
            V3sube_la.push_back(sqrt(TMath::Power(V3Values_err_la[i],2) + sqrt(TMath::Power(V3Values_err_low_la[i]/V3Values_low_la[i],2) + TMath::Power(JyieldSub_err_la[i]/JyieldSub_la[i],2) + TMath::Power(JyieldSub_err_low_la[i]/JyieldSub_low_la[i],2))*V3Values_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*JyieldSub_la[i]/JyieldSub_low_la[i]*sqrt(TMath::Power(V3Values_err_low_la[i]/V3Values_low_la[i],2) + TMath::Power(JyieldSub_err_la[i]/JyieldSub_la[i],2) + TMath::Power(JyieldSub_err_low_la[i]/JyieldSub_low_la[i],2))*V3Values_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*JyieldSub_la[i]/JyieldSub_low_la[i]));

            jetYfactor_la.push_back(JyieldSub_la[i]/JyieldSub_low_la[i]);
        }

        //Xi Vn correction
        std::vector<double> V2sub_xi;
        std::vector<double> V2sube_xi;

        std::vector<double> V3sub_xi;
        std::vector<double> V3sube_xi;

        std::vector<double> jetYfactor_xi;

        for(unsigned i=0; i<numPtBins_xi; i++){
            V2sub_xi.push_back(V2Values_xi[i] - V2Values_low_xi[i]*Nassoc_low_xi[i]/Nassoc_xi[i]*JyieldSub_xi[i]/JyieldSub_low_xi[i]);
            V2sube_xi.push_back(sqrt(TMath::Power(V2Values_err_xi[i],2) + sqrt(TMath::Power(V2Values_err_low_xi[i]/V2Values_low_xi[i],2) + TMath::Power(JyieldSub_err_xi[i]/JyieldSub_xi[i],2) + TMath::Power(JyieldSub_err_low_xi[i]/JyieldSub_low_xi[i],2))*V2Values_low_xi[i]*Nassoc_low_xi[i]/Nassoc_xi[i]*JyieldSub_xi[i]/JyieldSub_low_xi[i]*sqrt(TMath::Power(V2Values_err_low_xi[i]/V2Values_low_xi[i],2) + TMath::Power(JyieldSub_err_xi[i]/JyieldSub_xi[i],2) + TMath::Power(JyieldSub_err_low_xi[i]/JyieldSub_low_xi[i],2))*V2Values_low_xi[i]*Nassoc_low_xi[i]/Nassoc_xi[i]*JyieldSub_xi[i]/JyieldSub_low_xi[i]));

            V3sub_xi.push_back(V3Values_xi[i] - V3Values_low_xi[i]*Nassoc_low_xi[i]/Nassoc_xi[i]*JyieldSub_xi[i]/JyieldSub_low_xi[i]);
            V3sube_xi.push_back(sqrt(TMath::Power(V3Values_err_xi[i],2) + sqrt(TMath::Power(V3Values_err_low_xi[i]/V3Values_low_xi[i],2) + TMath::Power(JyieldSub_err_xi[i]/JyieldSub_xi[i],2) + TMath::Power(JyieldSub_err_low_xi[i]/JyieldSub_low_xi[i],2))*V3Values_low_xi[i]*Nassoc_low_xi[i]/Nassoc_xi[i]*JyieldSub_xi[i]/JyieldSub_low_xi[i]*sqrt(TMath::Power(V3Values_err_low_xi[i]/V3Values_low_xi[i],2) + TMath::Power(JyieldSub_err_xi[i]/JyieldSub_xi[i],2) + TMath::Power(JyieldSub_err_low_xi[i]/JyieldSub_low_xi[i],2))*V3Values_low_xi[i]*Nassoc_low_xi[i]/Nassoc_xi[i]*JyieldSub_xi[i]/JyieldSub_low_xi[i]));

            jetYfactor_xi.push_back(JyieldSub_xi[i]/JyieldSub_low_xi[i]);
        }


        //Ks vn calculation
        std::vector<double> v2sub_ks;
        std::vector<double> v2sube_ks;
        std::vector<double> v2_ks;
        std::vector<double> v2e_ks;
        std::vector<double> v2_low_ks;
        std::vector<double> v2e_low_ks;

        std::vector<double> v3sub_ks;
        std::vector<double> v3sube_ks;
        std::vector<double> v3_ks;
        std::vector<double> v3e_ks;
        std::vector<double> v3_low_ks;
        std::vector<double> v3e_low_ks;

        for(unsigned i=0; i<numPtBins_ks; i++){
            v2sub_ks.push_back(V2sub_ks[i]/v2sub_ref);
            v2sube_ks.push_back(fabs(sqrt(V2sube_ks[i]/V2sub_ks[i]*V2sube_ks[i]/V2sub_ks[i] + v2sube_ref/v2sub_ref*v2sube_ref/v2sub_ref)*v2sub_ks[i]));

            v2_ks.push_back(V2Values_ks[i]/v2_high_ref);
            v2e_ks.push_back(sqrt(TMath::Power(V2Values_err_ks[i]/V2Values_ks[i],2) + TMath::Power(v2e_high_ref/v2_high_ref,2))*v2_ks[i]);

            v2_low_ks.push_back(V2Values_low_ks[i]/v2_low_ref);
            v2e_low_ks.push_back(sqrt(TMath::Power(V2Values_err_low_ks[i]/V2Values_low_ks[i],2) + TMath::Power(v2e_low_ref/v2_low_ref,2))*v2_low_ks[i]);

            v3sub_ks.push_back(V3sub_ks[i]/v3sub_ref);
            v3sube_ks.push_back(fabs(sqrt(V3sube_ks[i]/V3sub_ks[i]*V3sube_ks[i]/V3sub_ks[i] + v3sube_ref/v3sub_ref*v3sube_ref/v3sub_ref))*v3sub_ks[i]);

            v3_ks.push_back(V3Values_ks[i]/v3_high_ref);
            v3e_ks.push_back(sqrt(TMath::Power(V3Values_err_ks[i]/V3Values_ks[i],2) + TMath::Power(v3e_high_ref/v3_high_ref,2))*v3_ks[i]);

            v3_low_ks.push_back(V3Values_low_ks[i]/v3_low_ref);
            v3e_low_ks.push_back(sqrt(TMath::Power(V3Values_err_low_ks[i]/V3Values_low_ks[i],2) + TMath::Power(v3e_low_ref/v3_low_ref,2))*v3_low_ks[i]);
        }

        //La vn calculation
        std::vector<double> v2sub_la;
        std::vector<double> v2sube_la;
        std::vector<double> v2_la;
        std::vector<double> v2e_la;
        std::vector<double> v2_low_la;
        std::vector<double> v2e_low_la;

        std::vector<double> v3sub_la;
        std::vector<double> v3sube_la;
        std::vector<double> v3_la;
        std::vector<double> v3e_la;
        std::vector<double> v3_low_la;
        std::vector<double> v3e_low_la;

        for(unsigned i=0; i<numPtBins_la; i++){
            v2sub_la.push_back(V2sub_la[i]/v2sub_ref);
            v2sube_la.push_back(fabs(sqrt(TMath::Power(V2sube_la[i]/V2sub_la[i],2) + TMath::Power(v2sube_ref/v2sub_ref,2))*v2sub_la[i]));

            v2_la.push_back(V2Values_la[i]/v2_high_ref);
            v2e_la.push_back(sqrt(TMath::Power(V2Values_err_la[i]/V2Values_la[i],2) + TMath::Power(v2e_high_ref/v2_high_ref,2))*v2_la[i]);

            v2_low_la.push_back(V2Values_low_la[i]/v2_low_ref);
            v2e_low_la.push_back(sqrt(TMath::Power(V2Values_err_low_la[i]/V2Values_low_la[i],2) + TMath::Power(v2e_low_ref/v2_low_ref,2))*v2_low_la[i]);

            v3sub_la.push_back(V3sub_la[i]/v3sub_ref);
            v3sube_la.push_back(fabs(sqrt(TMath::Power(V3sube_la[i]/V3sub_la[i],2) + TMath::Power(v3sube_ref/v3sub_ref,2)))*v3sub_la[i]);

            v3_la.push_back(V3Values_la[i]/v3_high_ref);
            v3e_la.push_back(sqrt(TMath::Power(V3Values_err_la[i]/V3Values_la[i],2) + TMath::Power(v3e_high_ref/v3_high_ref,2))*v3_la[i]);

            v3_low_la.push_back(V3Values_low_la[i]/v3_low_ref);
            v3e_low_la.push_back(sqrt(TMath::Power(V3Values_err_low_la[i]/V3Values_low_la[i],2) + TMath::Power(v3e_low_ref/v3_low_ref,2))*v3_low_la[i]);
        }

        //Xi vn calculation
        std::vector<double> v2sub_xi;
        std::vector<double> v2sube_xi;
        std::vector<double> v2_xi;
        std::vector<double> v2e_xi;
        std::vector<double> v2_low_xi;
        std::vector<double> v2e_low_xi;

        std::vector<double> v3sub_xi;
        std::vector<double> v3sube_xi;
        std::vector<double> v3_xi;
        std::vector<double> v3e_xi;
        std::vector<double> v3_low_xi;
        std::vector<double> v3e_low_xi;

        for(unsigned i=0; i<numPtBins_xi; i++){
            v2sub_xi.push_back(V2sub_xi[i]/v2sub_ref);
            v2sube_xi.push_back(fabs(sqrt(TMath::Power(V2sube_xi[i]/V2sub_xi[i],2) + TMath::Power(v2sube_ref/v2sub_ref,2)))*v2sub_xi[i]);

            v2_xi.push_back(V2Values_xi[i]/v2_high_ref);
            v2e_xi.push_back(sqrt(TMath::Power(V2Values_err_xi[i]/V2Values_xi[i],2) + TMath::Power(v2e_high_ref/v2_high_ref,2))*v2_xi[i]);

            v2_low_xi.push_back(V2Values_low_xi[i]/v2_low_ref);
            v2e_low_xi.push_back(sqrt(TMath::Power(V2Values_err_low_xi[i]/V2Values_low_xi[i],2) + TMath::Power(v2e_low_ref/v2_low_ref,2))*v2_low_xi[i]);

            v3sub_xi.push_back(V3sub_xi[i]/v3sub_ref);
            v3sube_xi.push_back(fabs(sqrt(TMath::Power(V3sube_xi[i]/V3sub_xi[i],2) + TMath::Power(v3sube_ref/v3sub_ref,2)))*v3sub_xi[i]);

            v3_xi.push_back(V3Values_xi[i]/v3_high_ref);
            v3e_xi.push_back(sqrt(TMath::Power(V3Values_err_xi[i]/V3Values_xi[i],2) + TMath::Power(v3e_high_ref/v3_high_ref,2))*v3_xi[i]);

            v3_low_xi.push_back(V3Values_low_xi[i]/v3_low_ref);
            v3e_low_xi.push_back(sqrt(TMath::Power(V3Values_err_low_xi[i]/V3Values_low_xi[i],2) + TMath::Power(v3e_low_ref/v3_low_ref,2))*v3_low_xi[i]);
        }

        //Produce RootFiles
        TFile output(PD.fn.c_str(),"UPDATE");

        //Get Mean Pt and KET of High N
        std::vector<double> pt_ks;
        std::vector<double> Ket_ks;
        for(unsigned i=0; i<numPtBins_ks; i++){
            TH1D* hpt = (TH1D*)PD.f_V0->Get(Form((PD.fn_V0 + "/Ptkshort_pt%d").c_str(),i));
            TH1D* hket = (TH1D*)PD.f_V0->Get(Form((PD.fn_V0 + "/KETkshort_pt%d").c_str(),i));

            pt_ks.push_back(hpt->GetMean(1));
            Ket_ks.push_back(hket->GetMean(1));
        }

        std::vector<double> pt_la;
        std::vector<double> Ket_la;
        for(unsigned i=0; i<numPtBins_la; i++){
            TH1D* hpt = (TH1D*)PD.f_V0->Get(Form((PD.fn_V0 + "/Ptlambda_pt%d").c_str(),i));
            TH1D* hket = (TH1D*)PD.f_V0->Get(Form((PD.fn_V0 + "/KETlambda_pt%d").c_str(),i));

            pt_la.push_back(hpt->GetMean(1));
            Ket_la.push_back(hket->GetMean(1));
        }

        std::vector<double> pt_xi;
        std::vector<double> Ket_xi;
        for(unsigned i=0; i<numPtBins_xi; i++){
            TH1D* hpt = (TH1D*)PD.f_Xi->Get(Form((PD.fn_Xi + "/Pt_xi_pt%d").c_str(),i));
            TH1D* hket = (TH1D*)PD.f_Xi->Get(Form((PD.fn_Xi + "/KET_xi_pt%d").c_str(),i));

            pt_xi.push_back(hpt->GetMean(1));
            Ket_xi.push_back(hket->GetMean(1));
        }

        std::vector<double> perisubfactor_ks;
        std::vector<double> perisubfactore_ks;
        for(int i=0; i<numPtBins_ks; i++){
            perisubfactor_ks.push_back(V2Values_low_ks[i]*Nassoc_low_ks[i]/JyieldSub_low_ks[i]);
            perisubfactore_ks.push_back(sqrt(TMath::Power(V2Values_err_low_ks[i]/V2Values_low_ks[i],2) + TMath::Power(JyieldSub_err_low_ks[i]/JyieldSub_low_ks[i],2))*perisubfactor_ks[i]);
        }
        std::vector<double> perisubfactor_la;
        std::vector<double> perisubfactore_la;
        for(int i=0; i<numPtBins_la; i++){
            perisubfactor_la.push_back(V2Values_low_la[i]*Nassoc_low_la[i]/JyieldSub_low_la[i]);
            perisubfactore_la.push_back(sqrt(TMath::Power(V2Values_err_low_la[i]/V2Values_low_la[i],2) + TMath::Power(JyieldSub_err_low_la[i]/JyieldSub_low_la[i],2))*perisubfactor_la[i]);
        }

        std::vector<double> perisubfactor_xi;
        std::vector<double> perisubfactore_xi;
        for(int i=0; i<numPtBins_xi; i++){
            perisubfactor_xi.push_back(V2Values_low_xi[i]*Nassoc_low_xi[i]/JyieldSub_low_xi[i]);
            perisubfactore_xi.push_back(sqrt(TMath::Power(V2Values_err_low_xi[i]/V2Values_low_xi[i],2) + TMath::Power(JyieldSub_err_low_xi[i]/JyieldSub_low_xi[i],2))*perisubfactor_xi[i]);
        }

        //Obs v2 values
        TGraphErrors* v2plot_ks        = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&v2_ks[0]      ,0,&v2e_ks[0]);
        TGraphErrors* v2plot_KET_ks    = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&v2_ks[0]      ,0,&v2e_ks[0]);
        TGraphErrors* v2subplot_ks     = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&v2sub_ks[0]   ,0,&v2sube_ks[0]);
        TGraphErrors* v2subplot_KET_ks = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&v2sub_ks[0]   ,0,&v2sube_ks[0]);
        TGraphErrors* v3plot_ks        = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&v3_ks[0]      ,0,&v3e_ks[0]);
        TGraphErrors* v3plot_KET_ks    = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&v3_ks[0]      ,0,&v3e_ks[0]);
        TGraphErrors* v3subplot_ks     = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&v3sub_ks[0]   ,0,&v3sube_ks[0]);
        TGraphErrors* v3subplot_KET_ks = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&v3sub_ks[0]   ,0,&v3sube_ks[0]);

        TGraphErrors* v2plot_la        = new TGraphErrors(numPtBins_la,&pt_la[0] ,&v2_la[0]      ,0,&v2e_la[0]);
        TGraphErrors* v2plot_KET_la    = new TGraphErrors(numPtBins_la,&Ket_la[0],&v2_la[0]      ,0,&v2e_la[0]);
        TGraphErrors* v2subplot_la     = new TGraphErrors(numPtBins_la,&pt_la[0] ,&v2sub_la[0]   ,0,&v2sube_la[0]);
        TGraphErrors* v2subplot_KET_la = new TGraphErrors(numPtBins_la,&Ket_la[0],&v2sub_la[0]   ,0,&v2sube_la[0]);
        TGraphErrors* v3plot_la        = new TGraphErrors(numPtBins_la,&pt_la[0] ,&v3_la[0]      ,0,&v3e_la[0]);
        TGraphErrors* v3plot_KET_la    = new TGraphErrors(numPtBins_la,&Ket_la[0],&v3_la[0]      ,0,&v3e_la[0]);
        TGraphErrors* v3subplot_la     = new TGraphErrors(numPtBins_la,&pt_la[0] ,&v3sub_la[0]   ,0,&v3sube_la[0]);
        TGraphErrors* v3subplot_KET_la = new TGraphErrors(numPtBins_la,&Ket_la[0],&v3sub_la[0]   ,0,&v3sube_la[0]);

        TGraphErrors* v2plot_xi        = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&v2_xi[0]      ,0,&v2e_xi[0]);
        TGraphErrors* v2plot_KET_xi    = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&v2_xi[0]      ,0,&v2e_xi[0]);
        TGraphErrors* v2subplot_xi     = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&v2sub_xi[0]   ,0,&v2sube_xi[0]);
        TGraphErrors* v2subplot_KET_xi = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&v2sub_xi[0]   ,0,&v2sube_xi[0]);
        TGraphErrors* v3plot_xi        = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&v3_xi[0]      ,0,&v3e_xi[0]);
        TGraphErrors* v3plot_KET_xi    = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&v3_xi[0]      ,0,&v3e_xi[0]);
        TGraphErrors* v3subplot_xi     = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&v3sub_xi[0]   ,0,&v3sube_xi[0]);
        TGraphErrors* v3subplot_KET_xi = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&v3sub_xi[0]   ,0,&v3sube_xi[0]);

        //Obs V2 values
        TGraphErrors* V2plot_ks        = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&V2Values_ks[0],0,&V2Values_err_ks[0]);
        TGraphErrors* V2plot_KET_ks    = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&V2Values_ks[0],0,&V2Values_err_ks[0]);
        TGraphErrors* V2subplot_ks     = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&V2sub_ks[0]   ,0,&V2sube_ks[0]);
        TGraphErrors* V2subplot_KET_ks = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&V2sub_ks[0]   ,0,&V2sube_ks[0]);
        TGraphErrors* V3plot_ks        = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&V3Values_ks[0],0,&V3Values_err_ks[0]);
        TGraphErrors* V3plot_KET_ks    = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&V3Values_ks[0],0,&V3Values_err_ks[0]);
        TGraphErrors* V3subplot_ks     = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&V3sub_ks[0]   ,0,&V3sube_ks[0]);
        TGraphErrors* V3subplot_KET_ks = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&V3sub_ks[0]   ,0,&V3sube_ks[0]);

        TGraphErrors* V2plot_la        = new TGraphErrors(numPtBins_la,&pt_la[0] ,&V2Values_la[0],0,&V2Values_err_la[0]);
        TGraphErrors* V2plot_KET_la    = new TGraphErrors(numPtBins_la,&Ket_la[0],&V2Values_la[0],0,&V2Values_err_la[0]);
        TGraphErrors* V2subplot_la     = new TGraphErrors(numPtBins_la,&pt_la[0] ,&V2sub_la[0]   ,0,&V2sube_la[0]);
        TGraphErrors* V2subplot_KET_la = new TGraphErrors(numPtBins_la,&Ket_la[0],&V2sub_la[0]   ,0,&V2sube_la[0]);
        TGraphErrors* V3plot_la        = new TGraphErrors(numPtBins_la,&pt_la[0] ,&V3Values_la[0],0,&V3Values_err_la[0]);
        TGraphErrors* V3plot_KET_la    = new TGraphErrors(numPtBins_la,&Ket_la[0],&V3Values_la[0],0,&V3Values_err_la[0]);
        TGraphErrors* V3subplot_la     = new TGraphErrors(numPtBins_la,&pt_la[0] ,&V3sub_la[0]   ,0,&V3sube_la[0]);
        TGraphErrors* V3subplot_KET_la = new TGraphErrors(numPtBins_la,&Ket_la[0],&V3sub_la[0]   ,0,&V3sube_la[0]);

        TGraphErrors* V2plot_xi        = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&V2Values_xi[0],0,&V2Values_err_xi[0]);
        TGraphErrors* V2plot_KET_xi    = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&V2Values_xi[0],0,&V2Values_err_xi[0]);
        TGraphErrors* V2subplot_xi     = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&V2sub_xi[0]   ,0,&V2sube_xi[0]);
        TGraphErrors* V2subplot_KET_xi = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&V2sub_xi[0]   ,0,&V2sube_xi[0]);
        TGraphErrors* V3plot_xi        = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&V3Values_xi[0],0,&V3Values_err_xi[0]);
        TGraphErrors* V3plot_KET_xi    = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&V3Values_xi[0],0,&V3Values_err_xi[0]);
        TGraphErrors* V3subplot_xi     = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&V3sub_xi[0]   ,0,&V3sube_xi[0]);
        TGraphErrors* V3subplot_KET_xi = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&V3sub_xi[0]   ,0,&V3sube_xi[0]);

        //Yields
        TGraphErrors* srYieldPlot_ks       = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Jyieldsr_ks[0]     ,0,&Jyieldsr_err_ks[0]);
        TGraphErrors* lrYieldPlot_ks       = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Jyieldlr_ks[0]     ,0,&Jyieldlr_err_ks[0]);
        TGraphErrors* subYieldPlot_ks      = new TGraphErrors(numPtBins_ks,&pt_ks[0],&JyieldSub_ks[0]    ,0,&JyieldSub_err_ks[0]);
        TGraphErrors* perisubfactorplot_ks = new TGraphErrors(numPtBins_ks,&pt_ks[0],&perisubfactor_ks[0],0,&perisubfactore_ks[0]);

        TGraphErrors* srYieldPlot_la       = new TGraphErrors(numPtBins_la,&pt_la[0],&Jyieldsr_la[0]     ,0,&Jyieldsr_err_la[0]);
        TGraphErrors* lrYieldPlot_la       = new TGraphErrors(numPtBins_la,&pt_la[0],&Jyieldlr_la[0]     ,0,&Jyieldlr_err_la[0]);
        TGraphErrors* subYieldPlot_la      = new TGraphErrors(numPtBins_la,&pt_la[0],&JyieldSub_la[0]    ,0,&JyieldSub_err_la[0]);
        TGraphErrors* perisubfactorplot_la = new TGraphErrors(numPtBins_la,&pt_la[0],&perisubfactor_la[0],0,&perisubfactore_la[0]);

        TGraphErrors* srYieldPlot_xi       = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Jyieldsr_xi[0]     ,0,&Jyieldsr_err_xi[0]);
        TGraphErrors* lrYieldPlot_xi       = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Jyieldlr_xi[0]     ,0,&Jyieldlr_err_xi[0]);
        TGraphErrors* subYieldPlot_xi      = new TGraphErrors(numPtBins_xi,&pt_xi[0],&JyieldSub_xi[0]    ,0,&JyieldSub_err_xi[0]);
        TGraphErrors* perisubfactorplot_xi = new TGraphErrors(numPtBins_xi,&pt_xi[0],&perisubfactor_xi[0],0,&perisubfactore_xi[0]);

        //For Direct subtraction method
        TGraphErrors* Nasslow_ks = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Nassoc_low_ks[0],0,0);
        TGraphErrors* Nass_ks = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Nassoc_ks[0],0,0);
        TGraphErrors* subYieldPlotLow_ks = new TGraphErrors(numPtBins_ks,&pt_ks[0],&JyieldSub_low_ks[0],0,&JyieldSub_err_low_ks[0]);

        TGraphErrors* Nasslow_la = new TGraphErrors(numPtBins_la,&pt_la[0],&Nassoc_low_la[0],0,0);
        TGraphErrors* Nass_la = new TGraphErrors(numPtBins_la,&pt_la[0],&Nassoc_la[0],0,0);
        TGraphErrors* subYieldPlotLow_la = new TGraphErrors(numPtBins_la,&pt_la[0],&JyieldSub_low_la[0],0,&JyieldSub_err_low_la[0]);

        TGraphErrors* Nasslow_xi = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Nassoc_low_xi[0],0,0);
        TGraphErrors* Nass_xi = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Nassoc_xi[0],0,0);
        TGraphErrors* subYieldPlotLow_xi = new TGraphErrors(numPtBins_xi,&pt_xi[0],&JyieldSub_low_xi[0],0,&JyieldSub_err_low_xi[0]);

        TGraphErrors* bz_ks = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Bz_ks[0],0,0);
        TGraphErrors* bz_la = new TGraphErrors(numPtBins_la,&pt_la[0],&Bz_la[0],0,0);
        TGraphErrors* bz_xi = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Bz_xi[0],0,0);

        std::string region;
        if(j == 0) region = "obs";
        if(j == 1) region = "bkg";

        v2plot_ks       ->Write(("v2plot_" + region + "_ks"       ).c_str(),TObject::kOverwrite);
        v2plot_KET_ks   ->Write(("v2plot_KET_" + region + "_ks"   ).c_str(),TObject::kOverwrite);
        v2subplot_ks    ->Write(("v2subplot_" + region + "_ks"    ).c_str(),TObject::kOverwrite);
        v2subplot_KET_ks->Write(("v2subplot_KET_" + region + "_ks").c_str(),TObject::kOverwrite);
        v3plot_ks       ->Write(("v3plot_" + region + "_ks"       ).c_str(),TObject::kOverwrite);
        v3plot_KET_ks   ->Write(("v3plot_KET_" + region + "_ks"   ).c_str(),TObject::kOverwrite);
        v3subplot_ks    ->Write(("v3subplot_" + region + "_ks"    ).c_str(),TObject::kOverwrite);
        v3subplot_KET_ks->Write(("v3subplot_KET_" + region + "_ks").c_str(),TObject::kOverwrite);

        v2plot_la       ->Write(("v2plot_" + region + "_la"       ).c_str(),TObject::kOverwrite);
        v2plot_KET_la   ->Write(("v2plot_KET_" + region + "_la"   ).c_str(),TObject::kOverwrite);
        v2subplot_la    ->Write(("v2subplot_" + region + "_la"    ).c_str(),TObject::kOverwrite);
        v2subplot_KET_la->Write(("v2subplot_KET_" + region + "_la").c_str(),TObject::kOverwrite);
        v3plot_la       ->Write(("v3plot_" + region + "_la"       ).c_str(),TObject::kOverwrite);
        v3plot_KET_la   ->Write(("v3plot_KET_" + region + "_la"   ).c_str(),TObject::kOverwrite);
        v3subplot_la    ->Write(("v3subplot_" + region + "_la"    ).c_str(),TObject::kOverwrite);
        v3subplot_KET_la->Write(("v3subplot_KET_" + region + "_la").c_str(),TObject::kOverwrite);

        v2plot_xi       ->Write(("v2plot_" + region + "_xi"       ).c_str(),TObject::kOverwrite);
        v2plot_KET_xi   ->Write(("v2plot_KET_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        v2subplot_xi    ->Write(("v2subplot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);
        v2subplot_KET_xi->Write(("v2subplot_KET_" + region + "_xi").c_str(),TObject::kOverwrite);
        v3plot_xi       ->Write(("v3plot_" + region + "_xi"       ).c_str(),TObject::kOverwrite);
        v3plot_KET_xi   ->Write(("v3plot_KET_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        v3subplot_xi    ->Write(("v3subplot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);
        v3subplot_KET_xi->Write(("v3subplot_KET_" + region + "_xi").c_str(),TObject::kOverwrite);

        V2plot_ks       ->Write(("V2plot_" + region + "_ks"       ).c_str(),TObject::kOverwrite);
        V2plot_KET_ks   ->Write(("V2plot_KET_" + region + "_ks"   ).c_str(),TObject::kOverwrite);
        V2subplot_ks    ->Write(("V2subplot_" + region + "_ks"    ).c_str(),TObject::kOverwrite);
        V2subplot_KET_ks->Write(("V2subplot_KET_" + region + "_ks").c_str(),TObject::kOverwrite);
        V3plot_ks       ->Write(("V3plot_" + region + "_ks"       ).c_str(),TObject::kOverwrite);
        V3plot_KET_ks   ->Write(("V3plot_KET_" + region + "_ks"   ).c_str(),TObject::kOverwrite);
        V3subplot_ks    ->Write(("V3subplot_" + region + "_ks"    ).c_str(),TObject::kOverwrite);
        V3subplot_KET_ks->Write(("V3subplot_KET_" + region + "_ks").c_str(),TObject::kOverwrite);

        V2plot_la       ->Write(("V2plot_" + region + "_la"       ).c_str(),TObject::kOverwrite);
        V2plot_KET_la   ->Write(("V2plot_KET_" + region + "_la"   ).c_str(),TObject::kOverwrite);
        V2subplot_la    ->Write(("V2subplot_" + region + "_la"    ).c_str(),TObject::kOverwrite);
        V2subplot_KET_la->Write(("V2subplot_KET_" + region + "_la").c_str(),TObject::kOverwrite);
        V3plot_la       ->Write(("V3plot_" + region + "_la"       ).c_str(),TObject::kOverwrite);
        V3plot_KET_la   ->Write(("V3plot_KET_" + region + "_la"   ).c_str(),TObject::kOverwrite);
        V3subplot_la    ->Write(("V3subplot_" + region + "_la"    ).c_str(),TObject::kOverwrite);
        V3subplot_KET_la->Write(("V3subplot_KET_" + region + "_la").c_str(),TObject::kOverwrite);

        V2plot_xi       ->Write(("V2plot_" + region + "_xi"       ).c_str(),TObject::kOverwrite);
        V2plot_KET_xi   ->Write(("V2plot_KET_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        V2subplot_xi    ->Write(("V2subplot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);
        V2subplot_KET_xi->Write(("V2subplot_KET_" + region + "_xi").c_str(),TObject::kOverwrite);
        V3plot_xi       ->Write(("V3plot_" + region + "_xi"       ).c_str(),TObject::kOverwrite);
        V3plot_KET_xi   ->Write(("V3plot_KET_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        V3subplot_xi    ->Write(("V3subplot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);
        V3subplot_KET_xi->Write(("V3subplot_KET_" + region + "_xi").c_str(),TObject::kOverwrite);

        perisubfactorplot_ks->Write(("perisubfactor_" + region + "_ks").c_str(),TObject::kOverwrite);
        perisubfactorplot_la->Write(("perisubfactor_" + region + "_la").c_str(),TObject::kOverwrite);
        perisubfactorplot_xi->Write(("perisubfactor_" + region + "_xi").c_str(),TObject::kOverwrite);

        Nass_ks           ->Write(("Nassoc_" + region + "_ks"       ).c_str(),TObject::kOverwrite);
        Nasslow_ks        ->Write(("Nassoc_low_" + region + "_ks"   ).c_str(),TObject::kOverwrite);
        subYieldPlotLow_ks->Write(("YieldPlot_low_" + region + "_ks").c_str(),TObject::kOverwrite);
        subYieldPlot_ks   ->Write(("YieldPlot_" + region + "_ks"    ).c_str(),TObject::kOverwrite);

        Nass_la           ->Write(("Nassoc_" + region + "_la"       ).c_str(),TObject::kOverwrite);
        Nasslow_la        ->Write(("Nassoc_low_" + region + "_la"   ).c_str(),TObject::kOverwrite);
        subYieldPlotLow_la->Write(("YieldPlot_low_" + region + "_la").c_str(),TObject::kOverwrite);
        subYieldPlot_la   ->Write(("YieldPlot_" + region + "_la"    ).c_str(),TObject::kOverwrite);

        Nass_xi           ->Write(("Nassoc_" + region + "_xi"       ).c_str(),TObject::kOverwrite);
        Nasslow_xi        ->Write(("Nassoc_low_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        subYieldPlotLow_xi->Write(("YieldPlot_low_" + region + "_xi").c_str(),TObject::kOverwrite);
        subYieldPlot_xi   ->Write(("YieldPlot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);

        bz_ks->Write(("bz_" + region + "_ks").c_str(),TObject::kOverwrite);
        bz_la->Write(("bz_" + region + "_la").c_str(),TObject::kOverwrite);
        bz_xi->Write(("bz_" + region + "_xi").c_str(),TObject::kOverwrite);
    }
}

/*
void PeriSubBkg(ParticleData PD)
{
    TH1::SetDefaultSumw2();

    //TFile *f_perisub = TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/PeripheralSubtractionMB.root");
    //TFile *f_V0 = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityCorrectMultB_09_19_17.root");
    //TFile *f_low_ref = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/HadCorr/Combine_MB0_corr_ref.root");
    //TFile *f_high_ref = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/HadCorr/XiCorrelationRapidityTotal_08_20_2017.root");

    double BW2D = (9.9/33)*((2-1.0/16)*3.1416/31);

    TLatex* tex = new TLatex();
    tex->SetNDC();

    TCanvas* csSidesr_low_ref  = new TCanvas("csSidesr_low_ref","csSidesr_low_ref",600,600);
    TCanvas* csSidelr_low_ref  = new TCanvas("csSidelr_low_ref","csSidelr_low_ref",600,600);
    TCanvas* csSidesr_high_ref = new TCanvas("csSidesr_high_ref","csSidesr_high_ref",600,600);
    TCanvas* csSidelr_high_ref = new TCanvas("csSidelr_high_ref","csSidelr_high_ref",600,600);
    TCanvas* csSidesr_low_ks   = new TCanvas("csSidesr_low_ks","csSidesr_low_ks",1200,900);
    TCanvas* csSidelr_low_ks   = new TCanvas("csSidelr_low_ks","csSidelr_low_ks",1200,900);
    TCanvas* csSidesr_ks       = new TCanvas("csSidesr_ks","csSidesr_ks",1200,900);
    TCanvas* csSidelr_ks       = new TCanvas("csSidelr_ks","csSidelr_ks",1200,900);
    TCanvas* csSidesr_low_la   = new TCanvas("csSidesr_low_la","csSidesr_low_la",1200,900);
    TCanvas* csSidelr_low_la   = new TCanvas("csSidelr_low_la","csSidelr_low_la",1200,900);
    TCanvas* csSidesr_la       = new TCanvas("csSidesr_la","csSidesr_la",1200,900);
    TCanvas* csSidelr_la       = new TCanvas("csSidelr_la","csSidelr_la",1200,900);
    TCanvas* csSidesr_low_xi   = new TCanvas("csSidesr_low_xi","csSidesr_low_xi",1200,900);
    TCanvas* csSidelr_low_xi   = new TCanvas("csSidelr_low_xi","csSidelr_low_xi",1200,900);
    TCanvas* csSidesr_xi       = new TCanvas("csSidesr_xi","csSidesr_xi",1200,900);
    TCanvas* csSidelr_xi       = new TCanvas("csSidelr_xi","csSidelr_xi",1200,900);
    TCanvas* c                 = new TCanvas("c","c",400,400);

    csSidesr_low_ks->Divide(4,4);
    csSidelr_low_ks->Divide(4,4);
    csSidesr_low_la->Divide(4,3);
    csSidelr_low_la->Divide(4,3);
    csSidesr_ks->Divide(4,4);
    csSidelr_ks->Divide(4,4);
    csSidesr_la->Divide(4,3);
    csSidelr_la->Divide(4,3);

    //Low N ref
    //Containers
    std::vector<double> Nassoc_low_ref;
    std::vector<double> Jyieldsr_low_ref;
    std::vector<double> Jyieldlr_low_ref;
    double Jyieldsr_err_low_ref;
    double Jyieldlr_err_low_ref;
    std::vector<double> V2Values_low_ref;
    std::vector<double> V2Values_err_low_ref;
    std::vector<double> V3Values_low_ref;
    std::vector<double> V3Values_err_low_ref;
    std::vector<double> Bz_low_ref;
    std::vector<double> nEvent_low_ref;

    std::vector<double> JyieldSub_low_ref;
    std::vector<double> JyieldSub_err_low_ref;


    TH1D* hsSidesr_low_ref;
    TH1D* hbSidesr_low_ref;
    TH1D* hsSidelr_low_ref;
    TH1D* hbSidelr_low_ref;
    TH1D* V2lrs_low_ref;
    TH1D* V2lrb_low_ref;
    TH2D* hbackgroundSide_low_ref;
    TH2D* hsignalSide_low_ref;

    //Calculate Nassoc, Jet yield, Low N

    PD.f_low_ref->GetObject("pPbCorr/background",hbackgroundSide_low_ref);
    PD.f_low_ref->GetObject("pPbCorr/signal",hsignalSide_low_ref);
    TH1D* mult_low_ref = (TH1D*) PD.f_low_ref->Get("pPbCorr/mult");

    hbSidesr_low_ref = hbackgroundSide_low_ref->ProjectionY("hbSidesr_low_ref", 14, 20);
    hsSidesr_low_ref = hsignalSide_low_ref->ProjectionY("hsSidesr_low_ref", 14, 20);

    TF1* quadFit1 = new TF1("quadFit1","[0]*x^2+[1]*x+[2]",0.6,2.2);
    quadFit1->SetParameters(1,1,1);

    nEvent_low_ref.push_back(mult_low_ref->Integral(0,100000));
    nEvent_low_ref.push_back(mult_low_ref->Integral(3,100000));
    Bz_low_ref.push_back(hbackgroundSide_low_ref->GetBinContent(hbackgroundSide_low_ref->FindBin(0,0)));

    hsSidesr_low_ref->Divide(hbSidesr_low_ref);
    hsSidesr_low_ref->Scale(Bz_low_ref[0]/nEvent_low_ref[0]/BW2D);

    hsSidesr_low_ref->Fit("quadFit1","R");
    hsSidesr_low_ref->Fit("quadFit1","R");
    hsSidesr_low_ref->Fit("quadFit1","R");


    double minVal_sr = quadFit1->GetMinimum(0.6,2.2);
    double minVal_srX = quadFit1->GetMinimumX(0.6,2.2);
    TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    minConst_sr->SetParameter(0,-minVal_sr);
    TH1D* hsSidesr_zeroed_low_ref = (TH1D*)hsSidesr_low_ref->Clone();
    hsSidesr_zeroed_low_ref->Add(minConst_sr);
    csSidesr_low_ref->cd();
    hsSidesr_low_ref->Draw();
    double xcoor = 0.52;
    double ycoor = 0.90;
    double increment = 0.07;
    tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
    tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",0.3,3.0));
    tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
    c->cd();
    Jyieldsr_low_ref.push_back(hsSidesr_zeroed_low_ref->IntegralAndError(hsSidesr_zeroed_low_ref->FindBin(0.0),hsSidesr_zeroed_low_ref->FindBin(minVal_srX),Jyieldsr_err_low_ref,"width"));
    double bin0yield = hsSidesr_zeroed_low_ref->GetBinContent(hsSidesr_zeroed_low_ref->FindBin(0.0))*0.19635;
    Jyieldsr_low_ref[0] = Jyieldsr_low_ref[0]*2 - bin0yield;

    hbSidelr_low_ref = hbackgroundSide_low_ref->ProjectionY("hbSidelr_low_ref",1,14);
    TH1D* ahbSidelr_low_ref = hbackgroundSide_low_ref->ProjectionY("ahbSidelr_low_ref",20,33);
    hsSidelr_low_ref = hsignalSide_low_ref->ProjectionY("hsSidelr_low_ref",1,14);
    TH1D* ahsSidelr_low_ref = hsignalSide_low_ref->ProjectionY("ahsSidelr_low_ref",20,33);

    hbSidelr_low_ref->Add(ahbSidelr_low_ref);
    hsSidelr_low_ref->Add(ahsSidelr_low_ref);
    hsSidelr_low_ref->Divide(hbSidelr_low_ref);
    hsSidelr_low_ref->Scale(Bz_low_ref[0]/nEvent_low_ref[0]/BW2D);

    TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
    quadFit2->SetParameters(1,1,1);

    hsSidelr_low_ref->Fit("quadFit2","R");
    hsSidelr_low_ref->Fit("quadFit2","R");
    hsSidelr_low_ref->Fit("quadFit2","R");

    double minVal_lr = quadFit2->GetMinimum(0.6,2.2);
    double minVal_lrX = quadFit2->GetMinimumX(0.6,2.2);
    TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    minConst_lr->SetParameter(0,-minVal_lr);
    TH1D* hsSidelr_zeroed_low_ref = (TH1D*)hsSidelr_low_ref->Clone();
    hsSidelr_zeroed_low_ref->Add(minConst_lr);
    csSidelr_low_ref->cd();
    hsSidelr_low_ref->Draw();
    xcoor = 0.52;
    ycoor = 0.90;
    increment = 0.07;
    tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
    tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",0.3,3.0));
    tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
    c->cd();
    Jyieldlr_low_ref.push_back(hsSidelr_zeroed_low_ref->IntegralAndError(hsSidelr_zeroed_low_ref->FindBin(0.0),hsSidelr_zeroed_low_ref->FindBin(minVal_srX),Jyieldlr_err_low_ref,"width"));
    bin0yield = hsSidelr_zeroed_low_ref->GetBinContent(hsSidelr_zeroed_low_ref->FindBin(0.0))*0.19635;
    Jyieldlr_low_ref[0] = Jyieldlr_low_ref[0]*2 - bin0yield;

    JyieldSub_low_ref.push_back(Jyieldsr_low_ref[0] - Jyieldlr_low_ref[0]);
    JyieldSub_err_low_ref.push_back(sqrt(Jyieldsr_err_low_ref*Jyieldsr_err_low_ref + Jyieldlr_err_low_ref*Jyieldlr_err_low_ref));

    TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit1->SetParNames("N","V1","V2","V3");
    fit1->SetParameters(10,1,1,1);
    fit1->SetLineColor(2);

    V2lrs_low_ref = hsignalSide_low_ref->ProjectionY("V2lrs_low_ref",1,14);
    TH1D* aV2lrs_low_ref = hsignalSide_low_ref->ProjectionY("aV2lrs_low_ref",20,33);
    V2lrb_low_ref = hbackgroundSide_low_ref->ProjectionY("V2lrb_low_ref",1,14);
    TH1D* aV2lrb_low_ref = hbackgroundSide_low_ref->ProjectionY("aV2lrb_low_ref",20,33);
    V2lrs_low_ref->Add(aV2lrs_low_ref);
    V2lrb_low_ref->Add(aV2lrb_low_ref);
    V2lrs_low_ref->Divide(V2lrb_low_ref);
    V2lrs_low_ref->Scale(Bz_low_ref[0]/nEvent_low_ref[0]/BW2D);

    V2lrs_low_ref->Fit("fit1","R");
    V2lrs_low_ref->Fit("fit1","R");
    V2lrs_low_ref->Fit("fit1","R");

    V2Values_low_ref.push_back(fit1->GetParameter(2));
    V2Values_err_low_ref.push_back(fit1->GetParError(2));
    V3Values_low_ref.push_back(fit1->GetParameter(3));
    V3Values_err_low_ref.push_back(fit1->GetParError(3));

    double v2_low_ref = sqrt(V2Values_low_ref[0]);
    double v2e_low_ref = sqrt(V2Values_low_ref[0])*(V2Values_err_low_ref[0]/V2Values_low_ref[0])/2;
    double v3_low_ref = sqrt(V3Values_low_ref[0]);
    double v3e_low_ref = sqrt(V3Values_low_ref[0])*(V3Values_err_low_ref[0]/V3Values_low_ref[0])/2;

    Nassoc_low_ref.push_back(fit1->GetParameter(0));
    c->cd();

    //High N ref
    //Containers
    std::vector<double> Nassoc_high_ref;
    std::vector<double> Jyieldsr_high_ref;
    std::vector<double> Jyieldlr_high_ref;
    double Jyieldsr_err_high_ref;
    double Jyieldlr_err_high_ref;
    std::vector<double> V2Values_high_ref;
    std::vector<double> V2Values_err_high_ref;
    std::vector<double> V3Values_high_ref;
    std::vector<double> V3Values_err_high_ref;
    std::vector<double> Bz_high_ref;
    std::vector<double> nEvent_high_ref;

    std::vector<double> JyieldSub_high_ref;
    std::vector<double> JyieldSub_err_high_ref;

    TH1D* hsSidesr_high_ref;;
    TH1D* hbSidesr_high_ref;
    TH1D* hsSidelr_high_ref;
    TH1D* hbSidelr_high_ref;
    TH1D* V2lrs_high_ref;
    TH1D* V2lrb_high_ref;
    TH2D* hbackgroundSide_high_ref;
    TH2D* hsignalSide_high_ref;

    //Calculate Nassoc, Jet yield, High N

    PD.f_high_ref->GetObject("xiCorrelationRapidity/BackgroundHad",hbackgroundSide_high_ref);
    PD.f_high_ref->GetObject("xiCorrelationRapidity/SignalHad",hsignalSide_high_ref);
    TH1D* mult_high_ref = (TH1D*) PD.f_high_ref->Get("xiCorrelationRapidity/nTrk");

    hbSidesr_high_ref = hbackgroundSide_high_ref->ProjectionY("hbSidesr_high_ref", 14, 20);
    hsSidesr_high_ref = hsignalSide_high_ref->ProjectionY("hsSidesr_high_ref", 14, 20);

    TF1* quadFit12 = new TF1("quadFit12","[0]*x^2+[1]*x+[2]",0.6,2.2);
    quadFit12->SetParameters(1,1,1);

    nEvent_high_ref.push_back(mult_high_ref->Integral(0,100000));
    nEvent_high_ref.push_back(mult_high_ref->Integral(3,100000));
    Bz_high_ref.push_back(hbackgroundSide_high_ref->GetBinContent(hbackgroundSide_high_ref->FindBin(0,0)));

    hsSidesr_high_ref->Divide(hbSidesr_high_ref);
    hsSidesr_high_ref->Scale(Bz_high_ref[0]/nEvent_high_ref[0]/BW2D);

    hsSidesr_high_ref->Fit("quadFit12","R");
    hsSidesr_high_ref->Fit("quadFit12","R");
    hsSidesr_high_ref->Fit("quadFit12","R");

    minVal_sr = quadFit12->GetMinimum(0.6,2.2);
    minVal_srX = quadFit12->GetMinimumX(0.6,2.2);
    TF1* minConst_sr2 = new TF1("minConst_sr2","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    minConst_sr2->SetParameter(0,-minVal_sr);
    TH1D* hsSidesr_zeroed_high_ref = (TH1D*)hsSidesr_high_ref->Clone();
    hsSidesr_zeroed_high_ref->Add(minConst_sr2);
    csSidesr_high_ref->cd();
    hsSidesr_high_ref->Draw();
    xcoor = 0.52;
    ycoor = 0.90;
    increment = 0.07;
    tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
    tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",0.3,3.0));
    tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
    c->cd();
    Jyieldsr_high_ref.push_back(hsSidesr_zeroed_high_ref->IntegralAndError(hsSidesr_zeroed_high_ref->FindBin(0.0),hsSidesr_zeroed_high_ref->FindBin(minVal_srX),Jyieldsr_err_high_ref,"width"));
    bin0yield = hsSidesr_zeroed_high_ref->GetBinContent(hsSidesr_zeroed_high_ref->FindBin(0.0))*0.19635;
    Jyieldsr_high_ref[0] = Jyieldsr_high_ref[0]*2 - bin0yield;

    hbSidelr_high_ref = hbackgroundSide_high_ref->ProjectionY("hbSidelr_high_ref",1,14);
    TH1D* ahbSidelr_high_ref = hbackgroundSide_high_ref->ProjectionY("ahbSidelr_high_ref",20,33);
    hsSidelr_high_ref = hsignalSide_high_ref->ProjectionY("hsSidelr_high_ref",1,14);
    TH1D* ahsSidelr_high_ref = hsignalSide_high_ref->ProjectionY("ahsSidelr_high_ref",20,33);

    hbSidelr_high_ref->Add(ahbSidelr_high_ref);
    hsSidelr_high_ref->Add(ahsSidelr_high_ref);
    hsSidelr_high_ref->Divide(hbSidelr_high_ref);
    hsSidelr_high_ref->Scale(Bz_high_ref[0]/nEvent_high_ref[0]/BW2D);

    TF1* quadFit21 = new TF1("quadFit21","[0]*x^2+[1]*x+[2]",0.6,2.2);
    quadFit21->SetParameters(1,1,1);

    hsSidelr_high_ref->Fit("quadFit21","R");
    hsSidelr_high_ref->Fit("quadFit21","R");
    hsSidelr_high_ref->Fit("quadFit21","R");

    minVal_lr = quadFit21->GetMinimum(0.6,2.2);
    minVal_lrX = quadFit21->GetMinimumX(0.6,2.2);
    TF1* minConst_lr2 = new TF1("minConst_lr2","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    minConst_lr2->SetParameter(0,-minVal_lr);
    TH1D* hsSidelr_zeroed_high_ref = (TH1D*)hsSidelr_high_ref->Clone();
    hsSidelr_zeroed_high_ref->Add(minConst_lr2);
    csSidelr_high_ref->cd();
    hsSidelr_high_ref->Draw();
    xcoor = 0.52;
    ycoor = 0.90;
    increment = 0.07;
    tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
    tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",0.3,3.0));
    tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
    c->cd();
    Jyieldlr_high_ref.push_back(hsSidelr_zeroed_high_ref->IntegralAndError(hsSidelr_zeroed_high_ref->FindBin(0.0),hsSidelr_zeroed_high_ref->FindBin(minVal_srX),Jyieldlr_err_high_ref,"width"));
    bin0yield = hsSidelr_zeroed_high_ref->GetBinContent(hsSidelr_zeroed_high_ref->FindBin(0.0))*0.19635;
    Jyieldlr_high_ref[0] = Jyieldlr_high_ref[0]*2 - bin0yield;

    JyieldSub_high_ref.push_back(Jyieldsr_high_ref[0] - Jyieldlr_high_ref[0]);
    JyieldSub_err_high_ref.push_back(sqrt(Jyieldsr_err_high_ref*Jyieldsr_err_high_ref + Jyieldlr_err_high_ref*Jyieldlr_err_high_ref));

    TF1* fit2 = new TF1("fit2","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit2->SetParNames("N","V1","V2","V3");
    fit2->SetParameters(10,1,1,1);
    fit2->SetLineColor(2);

    V2lrs_high_ref = hsignalSide_high_ref->ProjectionY("V2lrs_high_ref",1,14);
    TH1D* aV2lrs_high_ref = hsignalSide_high_ref->ProjectionY("aV2lrs_high_ref",20,33);
    V2lrb_high_ref = hbackgroundSide_high_ref->ProjectionY("V2lrb_high_ref",1,14);
    TH1D* aV2lrb_high_ref = hbackgroundSide_high_ref->ProjectionY("aV2lrb_high_ref",20,33);
    V2lrs_high_ref->Add(aV2lrs_high_ref);
    V2lrb_high_ref->Add(aV2lrb_high_ref);
    V2lrs_high_ref->Divide(V2lrb_high_ref);
    V2lrs_high_ref->Scale(Bz_high_ref[0]/nEvent_high_ref[0]/BW2D);

    V2lrs_high_ref->Fit("fit2","R");
    V2lrs_high_ref->Fit("fit2","R");
    V2lrs_high_ref->Fit("fit2","R");

    V2Values_high_ref.push_back(fit2->GetParameter(2));
    V2Values_err_high_ref.push_back(fit2->GetParError(2));
    V3Values_high_ref.push_back(fit2->GetParameter(3));
    V3Values_err_high_ref.push_back(fit2->GetParError(3));

    double v2_high_ref = sqrt(V2Values_high_ref[0]);
    double v2e_high_ref = sqrt(V2Values_high_ref[0])*(V2Values_err_high_ref[0]/V2Values_high_ref[0])/2;
    double v3_high_ref = sqrt(V3Values_high_ref[0]);
    double v3e_high_ref = sqrt(V3Values_high_ref[0])*(V3Values_err_high_ref[0]/V3Values_high_ref[0])/2;

    Nassoc_high_ref.push_back(fit2->GetParameter(0));
    c->cd();

    //Low N
//KSHORT
    //Containers
    std::vector<double> Nassoc_low_ks;
    std::vector<double> Jyieldsr_low_ks;
    std::vector<double> Jyieldlr_low_ks;
    std::vector<double> Jyieldsr_err_low_ks(18);
    std::vector<double> Jyieldlr_err_low_ks(18);
    std::vector<double> V2Values_low_ks;
    std::vector<double> V2Values_err_low_ks;
    std::vector<double> V3Values_low_ks;
    std::vector<double> V3Values_err_low_ks;
    std::vector<double> Bz_low_ks;
    std::vector<double> nEvent_low_ks;

    std::vector<double> JyieldSub_low_ks;
    std::vector<double> JyieldSub_err_low_ks;

    int arraySize_ks = PD.PtBin_ks.size();

    TH1D* hsSidesr_low_ks[arraySize_ks];
    TH1D* hbSidesr_low_ks[arraySize_ks];
    TH1D* hsSidelr_low_ks[arraySize_ks];
    TH1D* hbSidelr_low_ks[arraySize_ks];
    TH1D* V2lrs_low_ks[arraySize_ks];
    TH1D* V2lrb_low_ks[arraySize_ks];
    TH2D* hbackgroundSide_low_ks[arraySize_ks];
    TH2D* hsignalSide_low_ks[arraySize_ks];


    //Calculate Nassoc, Jet yield, Side region Low N

    for(int i=0; i<PD.PtBin_ks.size(); i++)
    {
        PD.f_perisub->GetObject(Form((PD.fn_v0cas + "/backgroundkshort_bkg_pt%d").c_str(),i),hbackgroundSide_low_ks[i]);
        PD.f_perisub->GetObject(Form((PD.fn_v0cas + "/signalkshort_bkg_pt%d").c_str(),i),hsignalSide_low_ks[i]);
        TH1D* mult_low_ks = (TH1D*) PD.f_perisub->Get(Form((PD.fn_v0cas + "/mult_ks_bkg_pt%d").c_str(),i));

        hbSidesr_low_ks[i] = hbackgroundSide_low_ks[i]->ProjectionY("hbSidesr_low_ks", 14, 20);
        hsSidesr_low_ks[i] = hsignalSide_low_ks[i]->ProjectionY("hsSidesr_low_ks", 14, 20);

        TF1* quadFit13 = new TF1("quadFit13","[0]*x^2+[1]*x+[2]",0.6,2.2);
        quadFit13->SetParameters(1,1,1);


        nEvent_low_ks.push_back(mult_low_ks->Integral(0,100000));
        //nEvent_low_ks.push_back(mult_low_ks->Integral(2,100000));
        Bz_low_ks.push_back(hbackgroundSide_low_ks[i]->GetBinContent(hbackgroundSide_low_ks[i]->FindBin(0,0)));

        hsSidesr_low_ks[i]->Divide(hbSidesr_low_ks[i]);
        hsSidesr_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);

        hsSidesr_low_ks[i]->Fit("quadFit13","R");
        hsSidesr_low_ks[i]->Fit("quadFit13","R");
        hsSidesr_low_ks[i]->Fit("quadFit13","R");

        double minVal_sr = quadFit13->GetMinimum(0.6,2.2);
        double minVal_srX = quadFit13->GetMinimumX(0.6,2.2);
        TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        minConst_sr->SetParameter(0,-minVal_sr);
        TH1D* hsSidesr_zeroed_low_ks = (TH1D*)hsSidesr_low_ks[i]->Clone();
        hsSidesr_zeroed_low_ks->Add(minConst_sr);
        csSidesr_low_ks->cd(i+1);
        hsSidesr_low_ks[i]->Draw();
        double xcoor = 0.52;
        double ycoor = 0.90;
        double increment = 0.07;
        tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
        tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_ks[i],PD.PtBin_ks[i+1]));
        tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
        c->cd();
        Jyieldsr_low_ks.push_back(hsSidesr_zeroed_low_ks->IntegralAndError(hsSidesr_zeroed_low_ks->FindBin(0),hsSidesr_zeroed_low_ks->FindBin(minVal_srX),Jyieldsr_err_low_ks[i],"width"));
        double bin0yield = hsSidesr_zeroed_low_ks->GetBinContent(hsSidesr_zeroed_low_ks->FindBin(0.0))*0.19635;
        Jyieldsr_low_ks[i] = Jyieldsr_low_ks[i]*2 - bin0yield;

        hbSidelr_low_ks[i] = hbackgroundSide_low_ks[i]->ProjectionY("hbSidelr_low_ks",1,14);
        TH1D* ahbSidelr_low_ks = hbackgroundSide_low_ks[i]->ProjectionY("ahbSidelr_low_ks",20,33);
        hsSidelr_low_ks[i] = hsignalSide_low_ks[i]->ProjectionY("hsSidelr_low_ks",1,14);
        TH1D* ahsSidelr_low_ks = hsignalSide_low_ks[i]->ProjectionY("ahsSidelr_low_ks",20,33);

        hbSidelr_low_ks[i]->Add(ahbSidelr_low_ks);
        hsSidelr_low_ks[i]->Add(ahsSidelr_low_ks);
        hsSidelr_low_ks[i]->Divide(hbSidelr_low_ks[i]);
        hsSidelr_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);

        TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
        quadFit2->SetParameters(1,1,1);


        hsSidelr_low_ks[i]->Fit("quadFit2","R");
        hsSidelr_low_ks[i]->Fit("quadFit2","R");
        hsSidelr_low_ks[i]->Fit("quadFit2","R");

        double minVal_lr = quadFit2->GetMinimum(0.6,2.2);
        double minVal_lrX = quadFit2->GetMinimumX(0.6,2.2);
        TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        minConst_lr->SetParameter(0,-minVal_lr);
        TH1D* hsSidelr_zeroed_low_ks = (TH1D*)hsSidelr_low_ks[i]->Clone();
        hsSidelr_zeroed_low_ks->Add(minConst_lr);
        csSidelr_low_ks->cd(i+1);
        hsSidelr_low_ks[i]->Draw();
        xcoor = 0.52;
        ycoor = 0.90;
        increment = 0.07;
        tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
        tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_ks[i],PD.PtBin_ks[i+1]));
        tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
        c->cd();
        Jyieldlr_low_ks.push_back(hsSidelr_zeroed_low_ks->IntegralAndError(hsSidelr_zeroed_low_ks->FindBin(0.0),hsSidelr_zeroed_low_ks->FindBin(minVal_srX),Jyieldlr_err_low_ks[i],"width"));
        bin0yield = hsSidelr_zeroed_low_ks->GetBinContent(hsSidelr_zeroed_low_ks->FindBin(0.0))*0.19635;
        Jyieldlr_low_ks[i] = Jyieldlr_low_ks[i]*2 - bin0yield;

        JyieldSub_low_ks.push_back(Jyieldsr_low_ks[i] - Jyieldlr_low_ks[i]);
        JyieldSub_err_low_ks.push_back(sqrt(Jyieldsr_err_low_ks[i]*Jyieldsr_err_low_ks[i] + Jyieldlr_err_low_ks[i]*Jyieldlr_err_low_ks[i]));

        TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        fit1->SetParNames("N","V1","V2","V3");
        fit1->SetParameters(10,1,1,1);
        fit1->SetLineColor(2);

        V2lrs_low_ks[i] = hsignalSide_low_ks[i]->ProjectionY(Form("V2lrs_low_ks%d",i),1,14);
        TH1D* aV2lrs_low_ks = hsignalSide_low_ks[i]->ProjectionY("aV2lrs_low_ks",20,33);
        V2lrb_low_ks[i] = hbackgroundSide_low_ks[i]->ProjectionY(Form("V2lrb_low_ks%d",i),1,14);
        TH1D* aV2lrb_low_ks = hbackgroundSide_low_ks[i]->ProjectionY("aV2lrb_low_ks",20,33);
        V2lrs_low_ks[i]->Add(aV2lrs_low_ks);
        V2lrb_low_ks[i]->Add(aV2lrb_low_ks);
        V2lrs_low_ks[i]->Divide(V2lrb_low_ks[i]);
        V2lrs_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);

        V2lrs_low_ks[i]->Fit("fit1","R");
        V2lrs_low_ks[i]->Fit("fit1","R");
        V2lrs_low_ks[i]->Fit("fit1","R");

        V2Values_low_ks.push_back(fit1->GetParameter(2));
        V2Values_err_low_ks.push_back(fit1->GetParError(2));
        V3Values_low_ks.push_back(fit1->GetParameter(3));
        V3Values_err_low_ks.push_back(fit1->GetParError(3));

        Nassoc_low_ks.push_back(fit1->GetParameter(0));
        c->cd();
    }

//Lambda N low
    //Containers
    std::vector<double> Nassoc_low_la;
    std::vector<double> Jyieldsr_low_la;
    std::vector<double> Jyieldlr_low_la;
    std::vector<double> Jyieldsr_err_low_la(18);
    std::vector<double> Jyieldlr_err_low_la(18);
    std::vector<double> V2Values_low_la;
    std::vector<double> V2Values_err_low_la;
    std::vector<double> V3Values_low_la;
    std::vector<double> V3Values_err_low_la;
    std::vector<double> Bz_low_la;
    std::vector<double> nEvent_low_la;

    std::vector<double> JyieldSub_low_la;
    std::vector<double> JyieldSub_err_low_la;

    int arraySize_la = PD.PtBin_la.size();

    TH1D* hsSidesr_low_la[arraySize_la];
    TH1D* hbSidesr_low_la[arraySize_la];
    TH1D* hsSidelr_low_la[arraySize_la];
    TH1D* hbSidelr_low_la[arraySize_la];
    TH1D* V2lrs_low_la[arraySize_la];
    TH1D* V2lrb_low_la[arraySize_la];
    TH2D* hbackgroundSide_low_la[arraySize_la];
    TH2D* hsignalSide_low_la[arraySize_la];

    //Calculate Nassoc, Jet yield, Side region Low N

    for(int i=0; i<PD.PtBin_la.size(); i++)
    {
        PD.f_perisub->GetObject(Form((PD.fn_v0cas + "/backgroundlambda_bkg_pt%d").c_str(),i),hbackgroundSide_low_la[i]);
        PD.f_perisub->GetObject(Form((PD.fn_v0cas + "/signallambda_bkg_pt%d").c_str(),i),hsignalSide_low_la[i]);
        TH1D* mult_low_la = (TH1D*) PD.f_perisub->Get(Form((PD.fn_v0cas + "/mult_la_bkg_pt%d").c_str(),i));

        hbSidesr_low_la[i] = hbackgroundSide_low_la[i]->ProjectionY("hbSidesr_low_la", 14, 20);
        hsSidesr_low_la[i] = hsignalSide_low_la[i]->ProjectionY("hsSidesr_low_la", 14, 20);

        TF1* quadFit14 = new TF1("quadFit14","[0]*x^2+[1]*x+[2]",0.6,2.2);
        quadFit14->SetParameters(1,1,1);

        nEvent_low_la.push_back(mult_low_la->Integral(0,100000));
        //nEvent_low_la.push_back(mult_low_la->Integral(2,100000));
        Bz_low_la.push_back(hbackgroundSide_low_la[i]->GetBinContent(hbackgroundSide_low_la[i]->FindBin(0,0)));

        hsSidesr_low_la[i]->Divide(hbSidesr_low_la[i]);
        hsSidesr_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);

        hsSidesr_low_la[i]->Fit("quadFit14","R");
        hsSidesr_low_la[i]->Fit("quadFit14","R");
        hsSidesr_low_la[i]->Fit("quadFit14","R");

        double minVal_sr = quadFit14->GetMinimum(0.6,2.2);
        double minVal_srX = quadFit14->GetMinimumX(0.6,2.2);
        TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        minConst_sr->SetParameter(0,-minVal_sr);
        TH1D* hsSidesr_zeroed_low_la = (TH1D*)hsSidesr_low_la[i]->Clone();
        hsSidesr_zeroed_low_la->Add(minConst_sr);
        csSidesr_low_la->cd(i+1);
        hsSidesr_low_la[i]->Draw();
        double xcoor = 0.52;
        double ycoor = 0.90;
        double increment = 0.07;
        tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
        tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_la[i],PD.PtBin_la[i+1]));
        tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
        c->cd();
        Jyieldsr_low_la.push_back(hsSidesr_zeroed_low_la->IntegralAndError(hsSidesr_zeroed_low_la->FindBin(0.0),hsSidesr_zeroed_low_la->FindBin(minVal_srX),Jyieldsr_err_low_la[i],"width"));
        double bin0yield = hsSidesr_zeroed_low_la->GetBinContent(hsSidesr_zeroed_low_la->FindBin(0.0))*0.19635;
        Jyieldsr_low_la[i] = Jyieldsr_low_la[i]*2 - bin0yield;

        hbSidelr_low_la[i] = hbackgroundSide_low_la[i]->ProjectionY("hbSidelr_low_la",1,14);
        TH1D* ahbSidelr_low_la = hbackgroundSide_low_la[i]->ProjectionY("ahbSidelr_low_la",20,33);
        hsSidelr_low_la[i] = hsignalSide_low_la[i]->ProjectionY("hsSidelr_low_la",1,14);
        TH1D* ahsSidelr_low_la = hsignalSide_low_la[i]->ProjectionY("ahsSidelr_low_la",20,33);

        hbSidelr_low_la[i]->Add(ahbSidelr_low_la);
        hsSidelr_low_la[i]->Add(ahsSidelr_low_la);
        hsSidelr_low_la[i]->Divide(hbSidelr_low_la[i]);
        hsSidelr_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);

        TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
        quadFit2->SetParameters(1,1,1);

        hsSidelr_low_la[i]->Fit("quadFit2","R");
        hsSidelr_low_la[i]->Fit("quadFit2","R");
        hsSidelr_low_la[i]->Fit("quadFit2","R");

        double minVal_lr = quadFit2->GetMinimum(0.6,2.2);
        double minVal_lrX = quadFit2->GetMinimumX(0.6,2.2);
        TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        minConst_lr->SetParameter(0,-minVal_lr);
        TH1D* hsSidelr_zeroed_low_la = (TH1D*)hsSidelr_low_la[i]->Clone();
        hsSidelr_zeroed_low_la->Add(minConst_lr);
        csSidelr_low_la->cd(i+1);
        hsSidelr_low_la[i]->Draw();
        xcoor = 0.52;
        ycoor = 0.90;
        increment = 0.07;
        tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
        tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_la[i],PD.PtBin_la[i+1]));
        tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
        c->cd();
        Jyieldlr_low_la.push_back(hsSidelr_zeroed_low_la->IntegralAndError(hsSidelr_zeroed_low_la->FindBin(0.0),hsSidelr_zeroed_low_la->FindBin(minVal_srX),Jyieldlr_err_low_la[i],"width"));
        bin0yield = hsSidelr_zeroed_low_la->GetBinContent(hsSidelr_zeroed_low_la->FindBin(0.0))*0.19635;
        Jyieldlr_low_la[i] = Jyieldlr_low_la[i]*2 - bin0yield;

        JyieldSub_low_la.push_back(Jyieldsr_low_la[i] - Jyieldlr_low_la[i]);
        JyieldSub_err_low_la.push_back(sqrt(Jyieldsr_err_low_la[i]*Jyieldsr_err_low_la[i] + Jyieldlr_err_low_la[i]*Jyieldlr_err_low_la[i]));

        TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        fit1->SetParNames("N","V1","V2","V3");
        fit1->SetParameters(10,1,1,1);
        fit1->SetLineColor(2);

        V2lrs_low_la[i] = hsignalSide_low_la[i]->ProjectionY(Form("V2lrs_low_la%d",i),1,14);
        TH1D* aV2lrs_low_la = hsignalSide_low_la[i]->ProjectionY("aV2lrs_low_la",20,33);
        V2lrb_low_la[i] = hbackgroundSide_low_la[i]->ProjectionY(Form("V2lrb_low_la%d",i),1,14);
        TH1D* aV2lrb_low_la = hbackgroundSide_low_la[i]->ProjectionY("aV2lrb_low_la",20,33);
        V2lrs_low_la[i]->Add(aV2lrs_low_la);
        V2lrb_low_la[i]->Add(aV2lrb_low_la);
        V2lrs_low_la[i]->Divide(V2lrb_low_la[i]);
        V2lrs_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);

        V2lrs_low_la[i]->Fit("fit1","R");
        V2lrs_low_la[i]->Fit("fit1","R");
        V2lrs_low_la[i]->Fit("fit1","R");

        V2Values_low_la.push_back(fit1->GetParameter(2));
        V2Values_err_low_la.push_back(fit1->GetParError(2));
        V3Values_low_la.push_back(fit1->GetParameter(3));
        V3Values_err_low_la.push_back(fit1->GetParError(3));

        Nassoc_low_la.push_back(fit1->GetParameter(0));
        c->cd();
    }

    //High N
//KSHORT
    //Containers
    std::vector<double> Nassoc_ks;
    std::vector<double> Jyieldsr_ks;
    std::vector<double> Jyieldlr_ks;
    std::vector<double> Jyieldsr_err_ks(18);
    std::vector<double> Jyieldlr_err_ks(18);
    std::vector<double> V2Values_ks;
    std::vector<double> V2Values_err_ks;
    std::vector<double> V3Values_ks;
    std::vector<double> V3Values_err_ks;
    std::vector<double> Bz_ks;
    std::vector<double> nEvent_ks;

    std::vector<double> JyieldSub_ks;
    std::vector<double> JyieldSub_err_ks;

    //std::vector<double> PtBin_ks = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5};//, 10.0, 15.0, 20.0}; 

    TH1D* hsSidesr_ks[arraySize_ks];
    TH1D* hbSidesr_ks[arraySize_ks];
    TH1D* hsSidelr_ks[arraySize_ks];
    TH1D* hbSidelr_ks[arraySize_ks];
    TH1D* V2lrs_ks[arraySize_ks];
    TH1D* V2lrb_ks[arraySize_ks];
    TH2D* hbackgroundSide_ks[arraySize_ks];
    TH2D* hsignalSide_ks[arraySize_ks];

    //Calculate Nassoc, Jet yield, Side region
    //Jet Yield

    for(int i=0; i<PD.PtBin_ks.size(); i++)
    {
        PD.f_V0->GetObject(Form((PD.fn_V0 + "/backgroundkshort_bkg_pt%d").c_str(),i),hbackgroundSide_ks[i]);
        PD.f_V0->GetObject(Form((PD.fn_V0 + "/signalkshort_bkg_pt%d").c_str(),i),hsignalSide_ks[i]);
        TH1D* mult_ks = (TH1D*) PD.f_V0->Get(Form((PD.fn_V0 + "/mult_ks_bkg_pt%d").c_str(),i));

        hbSidesr_ks[i] = hbackgroundSide_ks[i]->ProjectionY("hbSidesr_ks", 14, 20);
        hsSidesr_ks[i] = hsignalSide_ks[i]->ProjectionY("hsSidesr_ks", 14, 20);

        TF1* quadFit15 = new TF1("quadFit15","[0]*x^2+[1]*x+[2]",0.6,2.2);
        quadFit15->SetParameters(1,1,1);

        nEvent_ks.push_back(mult_ks->Integral(0,100000));
        //nEvent_ks.push_back(mult_ks->Integral(2,100000));
        Bz_ks.push_back(hbackgroundSide_ks[i]->GetBinContent(hbackgroundSide_ks[i]->FindBin(0,0)));

        hsSidesr_ks[i]->Divide(hbSidesr_ks[i]);
        hsSidesr_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);

        hsSidesr_ks[i]->Fit("quadFit15","R");
        hsSidesr_ks[i]->Fit("quadFit15","R");
        hsSidesr_ks[i]->Fit("quadFit15","R");

        double minVal_sr = quadFit15->GetMinimum(0.6,2.2);
        double minVal_srX = quadFit15->GetMinimumX(0.6,2.2);
        TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        minConst_sr->SetParameter(0,-minVal_sr);
        TH1D* hsSidesr_zeroed_ks = (TH1D*)hsSidesr_ks[i]->Clone();
        hsSidesr_zeroed_ks->Add(minConst_sr);
        csSidesr_ks->cd(i+1);
        hsSidesr_ks[i]->Draw();
        double xcoor = 0.52;
        double ycoor = 0.90;
        double increment = 0.07;
        tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
        tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_ks[i],PD.PtBin_ks[i+1]));
        tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
        c->cd();
        Jyieldsr_ks.push_back(hsSidesr_zeroed_ks->IntegralAndError(hsSidesr_zeroed_ks->FindBin(0.0),hsSidesr_zeroed_ks->FindBin(minVal_srX),Jyieldsr_err_ks[i],"width"));
        double bin0yield = hsSidesr_zeroed_ks->GetBinContent(hsSidesr_zeroed_ks->FindBin(0.0))*0.19635;
        Jyieldsr_ks[i] = Jyieldsr_ks[i]*2 - bin0yield;

        hbSidelr_ks[i] = hbackgroundSide_ks[i]->ProjectionY("hbSidelr_ks",1,14);
        TH1D* ahbSidelr_ks = hbackgroundSide_ks[i]->ProjectionY("ahbSidelr_ks",20,33);
        hsSidelr_ks[i] = hsignalSide_ks[i]->ProjectionY("hsSidelr_ks",1,14);
        TH1D* ahsSidelr_ks = hsignalSide_ks[i]->ProjectionY("ahsSidelr_ks",20,33);

        hbSidelr_ks[i]->Add(ahbSidelr_ks);
        hsSidelr_ks[i]->Add(ahsSidelr_ks);
        hsSidelr_ks[i]->Divide(hbSidelr_ks[i]);
        hsSidelr_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);

        TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
        quadFit2->SetParameters(1,1,1);

        hsSidelr_ks[i]->Fit("quadFit2","R");
        hsSidelr_ks[i]->Fit("quadFit2","R");
        hsSidelr_ks[i]->Fit("quadFit2","R");

        double minVal_lr = quadFit2->GetMinimum(0.6,2.2);
        double minVal_lrX = quadFit2->GetMinimumX(0.6,2.2);
        TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        minConst_lr->SetParameter(0,-minVal_lr);
        TH1D* hsSidelr_zeroed_ks = (TH1D*)hsSidelr_ks[i]->Clone();
        hsSidelr_zeroed_ks->Add(minConst_lr);
        csSidelr_ks->cd(i+1);
        hsSidelr_ks[i]->Draw();
        xcoor = 0.52;
        ycoor = 0.90;
        increment = 0.07;
        tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
        tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_ks[i],PD.PtBin_ks[i+1]));
        tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
        c->cd();
        Jyieldlr_ks.push_back(hsSidelr_zeroed_ks->IntegralAndError(hsSidelr_zeroed_ks->FindBin(0.0),hsSidelr_zeroed_ks->FindBin(minVal_srX),Jyieldlr_err_ks[i],"width"));
        bin0yield = hsSidelr_zeroed_ks->GetBinContent(hsSidelr_zeroed_ks->FindBin(0.0))*0.19635;
        Jyieldlr_ks[i] = Jyieldlr_ks[i]*2 - bin0yield;

        JyieldSub_ks.push_back(Jyieldsr_ks[i] - Jyieldlr_ks[i]);
        JyieldSub_err_ks.push_back(sqrt(Jyieldsr_err_ks[i]*Jyieldsr_err_ks[i] + Jyieldlr_err_ks[i]*Jyieldlr_err_ks[i]));

        TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        fit1->SetParNames("N","V1","V2","V3");
        fit1->SetParameters(10,1,1,1);
        fit1->SetLineColor(2);

        V2lrs_ks[i] = hsignalSide_ks[i]->ProjectionY(Form("V2lrs_ks%d",i),1,14);
        TH1D* aV2lrs_ks = hsignalSide_ks[i]->ProjectionY("aV2lrs_ks",20,33);
        V2lrb_ks[i] = hbackgroundSide_ks[i]->ProjectionY(Form("V2lrb_ks%d",i),1,14);
        TH1D* aV2lrb_ks = hbackgroundSide_ks[i]->ProjectionY("aV2lrb_ks",20,33);
        V2lrs_ks[i]->Add(aV2lrs_ks);
        V2lrb_ks[i]->Add(aV2lrb_ks);
        V2lrs_ks[i]->Divide(V2lrb_ks[i]);
        V2lrs_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);

        V2lrs_ks[i]->Fit("fit1","R");
        V2lrs_ks[i]->Fit("fit1","R");
        V2lrs_ks[i]->Fit("fit1","R");

        V2Values_ks.push_back(fit1->GetParameter(2));
        V2Values_err_ks.push_back(fit1->GetParError(2));
        V3Values_ks.push_back(fit1->GetParameter(3));
        V3Values_err_ks.push_back(fit1->GetParError(3));

        Nassoc_ks.push_back(fit1->GetParameter(0));
        c->cd();
    }

//LAMBDA
    std::vector<double> Nassoc_la;
    std::vector<double> Jyieldsr_la;
    std::vector<double> Jyieldlr_la;
    std::vector<double> Jyieldsr_err_la(18);
    std::vector<double> Jyieldlr_err_la(18);
    std::vector<double> V2Values_la;
    std::vector<double> V2Values_err_la;
    std::vector<double> V3Values_la;
    std::vector<double> V3Values_err_la;
    std::vector<double> Bz_la;
    std::vector<double> nEvent_la;

    std::vector<double> JyieldSub_la;
    std::vector<double> JyieldSub_err_la;

    //std::vector<double> PD.PtBin_la = {0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5};//, 10.0, 15.0, 20.0}; 

    TH1D* hsSidesr_la[arraySize_la];
    TH1D* hbSidesr_la[arraySize_la];
    TH1D* hsSidelr_la[arraySize_la];
    TH1D* hbSidelr_la[arraySize_la];
    TH1D* V2lrs_la[arraySize_la];
    TH1D* V2lrb_la[arraySize_la];
    TH2D* hbackgroundSide_la[arraySize_la];
    TH2D* hsignalSide_la[arraySize_la];

    //Calculate Nassoc, Jet yield, Side region
    //Jet Yield

    for(int i=0; i<PD.PtBin_la.size(); i++)
    {
        PD.f_V0->GetObject(Form((PD.fn_V0 + "/backgroundlambda_bkg_pt%d").c_str(),i),hbackgroundSide_la[i]);
        PD.f_V0->GetObject(Form((PD.fn_V0 + "/signallambda_bkg_pt%d").c_str(),i),hsignalSide_la[i]);
        TH1D* mult_la = (TH1D*) PD.f_V0->Get(Form((PD.fn_V0 + "/mult_la_bkg_pt%d").c_str(),i));

        hbSidesr_la[i] = hbackgroundSide_la[i]->ProjectionY("hbSidesr_la", 14, 20);
        hsSidesr_la[i] = hsignalSide_la[i]->ProjectionY("hsSidesr_la", 14, 20);

        TF1* quadFit16 = new TF1("quadFit16","[0]*x^2+[1]*x+[2]",0.6,2.2);
        quadFit16->SetParameters(1,1,1);

        nEvent_la.push_back(mult_la->Integral(0,10000));
        //nEvent_la.push_back(mult_la->Integral(2,10000));
        Bz_la.push_back(hbackgroundSide_la[i]->GetBinContent(hbackgroundSide_la[i]->FindBin(0,0)));

        hsSidesr_la[i]->Divide(hbSidesr_la[i]);
        hsSidesr_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);

        hsSidesr_la[i]->Fit("quadFit16","R");
        hsSidesr_la[i]->Fit("quadFit16","R");
        hsSidesr_la[i]->Fit("quadFit16","R");

        double minVal_sr = quadFit16->GetMinimum(0.6,2.2);
        double minVal_srX = quadFit16->GetMinimumX(0.6,2.2);
        TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        minConst_sr->SetParameter(0,-minVal_sr);
        TH1D* hsSidesr_zeroed_la = (TH1D*)hsSidesr_la[i]->Clone();
        hsSidesr_zeroed_la->Add(minConst_sr);
        csSidesr_la->cd(i+1);
        hsSidesr_la[i]->Draw();
        double xcoor = 0.52;
        double ycoor = 0.90;
        double increment = 0.07;
        tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
        tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_la[i],PD.PtBin_la[i+1]));
        tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
        c->cd();
        Jyieldsr_la.push_back(hsSidesr_zeroed_la->IntegralAndError(hsSidesr_zeroed_la->FindBin(0.0),hsSidesr_zeroed_la->FindBin(minVal_srX),Jyieldsr_err_la[i],"width"));
        double bin0yield = hsSidesr_zeroed_la->GetBinContent(hsSidesr_zeroed_la->FindBin(0.0))*0.19635;
        Jyieldsr_la[i] = Jyieldsr_la[i]*2 - bin0yield;

        hbSidelr_la[i] = hbackgroundSide_la[i]->ProjectionY("hbSidelr_la",1,14);
        TH1D* ahbSidelr_la = hbackgroundSide_la[i]->ProjectionY("ahbSidelr_la",20,33);
        hsSidelr_la[i] = hsignalSide_la[i]->ProjectionY("hsSidelr_la",1,14);
        TH1D* ahsSidelr_la = hsignalSide_la[i]->ProjectionY("ahsSidelr_la",20,33);

        hbSidelr_la[i]->Add(ahbSidelr_la);
        hsSidelr_la[i]->Add(ahsSidelr_la);
        hsSidelr_la[i]->Divide(hbSidelr_la[i]);
        hsSidelr_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);

        TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
        quadFit2->SetParameters(1,1,1);

        hsSidelr_la[i]->Fit("quadFit2","R");
        hsSidelr_la[i]->Fit("quadFit2","R");
        hsSidelr_la[i]->Fit("quadFit2","R");

        double minVal_lr = quadFit2->GetMinimum(0.6,2.2);
        double minVal_lrX = quadFit2->GetMinimumX(0.6,2.2);
        TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        minConst_lr->SetParameter(0,-minVal_lr);
        TH1D* hsSidelr_zeroed_la = (TH1D*)hsSidelr_la[i]->Clone();
        hsSidelr_zeroed_la->Add(minConst_lr);
        csSidelr_la->cd(i+1);
        hsSidelr_la[i]->Draw();
        xcoor = 0.52;
        ycoor = 0.90;
        increment = 0.07;
        tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}<20");
        tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_la[i],PD.PtBin_la[i+1]));
        tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
        c->cd();
        Jyieldlr_la.push_back(hsSidelr_zeroed_la->IntegralAndError(hsSidelr_zeroed_la->FindBin(0.0),hsSidelr_zeroed_la->FindBin(minVal_srX),Jyieldlr_err_la[i],"width"));
        bin0yield = hsSidelr_zeroed_la->GetBinContent(hsSidelr_zeroed_la->FindBin(0.0))*0.19635;
        Jyieldlr_la[i] = Jyieldlr_la[i]*2 - bin0yield;

        JyieldSub_la.push_back(Jyieldsr_la[i] - Jyieldlr_la[i]);
        JyieldSub_err_la.push_back(sqrt(Jyieldsr_err_la[i]*Jyieldsr_err_la[i] + Jyieldlr_err_la[i]*Jyieldlr_err_la[i]));

        TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        fit1->SetParNames("N","V1","V2","V3");
        fit1->SetParameters(10,1,1,1);
        fit1->SetLineColor(2);

        V2lrs_la[i] = hsignalSide_la[i]->ProjectionY(Form("V2lrs_la%d",i),1,14);
        TH1D* aV2lrs_la = hsignalSide_la[i]->ProjectionY("aV2lrs_la",20,33);
        V2lrb_la[i] = hbackgroundSide_la[i]->ProjectionY(Form("V2lrb_la%d",i),1,14);
        TH1D* aV2lrb_la = hbackgroundSide_la[i]->ProjectionY("aV2lrb_la",20,33);
        V2lrs_la[i]->Add(aV2lrs_la);
        V2lrb_la[i]->Add(aV2lrb_la);
        V2lrs_la[i]->Divide(V2lrb_la[i]);
        V2lrs_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);

        V2lrs_la[i]->Fit("fit1","R");
        V2lrs_la[i]->Fit("fit1","R");
        V2lrs_la[i]->Fit("fit1","R");

        V2Values_la.push_back(fit1->GetParameter(2));
        V2Values_err_la.push_back(fit1->GetParError(2));
        V3Values_la.push_back(fit1->GetParameter(3));
        V3Values_err_la.push_back(fit1->GetParError(3));

        Nassoc_la.push_back(fit1->GetParameter(0));
        c->cd();
    }

    // Reference Vn after corrections
    double V2sub_ref = V2Values_high_ref[0] - V2Values_low_ref[0]*Nassoc_low_ref[0]/Nassoc_high_ref[0]*Jyieldsr_high_ref[0]/Jyieldsr_low_ref[0];
    double V2sube_ref = sqrt(TMath::Power(V2Values_err_high_ref[0],2) + TMath::Power(V2Values_err_low_ref[0]*Nassoc_low_ref[0]/Nassoc_high_ref[0]*Jyieldsr_high_ref[0]/Jyieldsr_low_ref[0],2));

    double V3sub_ref = V3Values_high_ref[0] - V3Values_low_ref[0]*Nassoc_low_ref[0]/Nassoc_high_ref[0]*Jyieldsr_high_ref[0]/Jyieldsr_low_ref[0];
    double V3sube_ref = sqrt(TMath::Power(V3Values_err_high_ref[0],2) + TMath::Power(V3Values_err_low_ref[0]*Nassoc_low_ref[0]/Nassoc_high_ref[0]*Jyieldsr_high_ref[0]/Jyieldsr_low_ref[0],2));


    double v2sub_ref = sqrt(V2sub_ref);
    double v2sube_ref = sqrt(V2sub_ref)*(V2sube_ref/V2sub_ref)/2;

    double v3sub_ref = sqrt(V3sub_ref);
    double v3sube_ref = sqrt(V3sub_ref)*(V3sube_ref/V3sub_ref)/2;



    //Ks Vn correction
    std::vector<double> V2sub_ks;
    std::vector<double> V2sube_ks;

    std::vector<double> V3sub_ks;
    std::vector<double> V3sube_ks;

    std::vector<double> jetYfactor_ks;

    for(unsigned i=0; i<PD.PtBin_ks.size(); i++){
        V2sub_ks.push_back(V2Values_ks[i] - V2Values_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*JyieldSub_ks[i]/JyieldSub_low_ks[i]);
        V2sube_ks.push_back(sqrt(TMath::Power(V2Values_err_ks[i],2) + sqrt(TMath::Power(V2Values_err_low_ks[i]/V2Values_low_ks[i],2) + TMath::Power(JyieldSub_err_ks[i]/JyieldSub_ks[i],2) + TMath::Power(JyieldSub_err_low_ks[i]/JyieldSub_low_ks[i],2))*V2Values_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*JyieldSub_ks[i]/JyieldSub_low_ks[i]*sqrt(TMath::Power(V2Values_err_low_ks[i]/V2Values_low_ks[i],2) + TMath::Power(JyieldSub_err_ks[i]/JyieldSub_ks[i],2) + TMath::Power(JyieldSub_err_low_ks[i]/JyieldSub_low_ks[i],2))*V2Values_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*JyieldSub_ks[i]/JyieldSub_low_ks[i]));

        V3sub_ks.push_back(V3Values_ks[i] - V3Values_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*JyieldSub_ks[i]/JyieldSub_low_ks[i]);
        V3sube_ks.push_back(sqrt(TMath::Power(V3Values_err_ks[i],2) + sqrt(TMath::Power(V3Values_err_low_ks[i]/V3Values_low_ks[i],2) + TMath::Power(JyieldSub_err_ks[i]/JyieldSub_ks[i],2) + TMath::Power(JyieldSub_err_low_ks[i]/JyieldSub_low_ks[i],2))*V3Values_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*JyieldSub_ks[i]/JyieldSub_low_ks[i]*sqrt(TMath::Power(V3Values_err_low_ks[i]/V3Values_low_ks[i],2) + TMath::Power(JyieldSub_err_ks[i]/JyieldSub_ks[i],2) + TMath::Power(JyieldSub_err_low_ks[i]/JyieldSub_low_ks[i],2))*V3Values_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*JyieldSub_ks[i]/JyieldSub_low_ks[i]));

        jetYfactor_ks.push_back(JyieldSub_ks[i]/JyieldSub_low_ks[i]);
    }

    //La Vn correction
    std::vector<double> V2sub_la;
    std::vector<double> V2sube_la;

    std::vector<double> V3sub_la;
    std::vector<double> V3sube_la;

    std::vector<double> jetYfactor_la;

    for(unsigned i=0; i<PD.PtBin_la.size(); i++){
        V2sub_la.push_back(V2Values_la[i] - V2Values_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*JyieldSub_la[i]/JyieldSub_low_la[i]);
        V2sube_la.push_back(sqrt(TMath::Power(V2Values_err_la[i],2) + sqrt(TMath::Power(V2Values_err_low_la[i]/V2Values_low_la[i],2) + TMath::Power(JyieldSub_err_la[i]/JyieldSub_la[i],2) + TMath::Power(JyieldSub_err_low_la[i]/JyieldSub_low_la[i],2))*V2Values_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*JyieldSub_la[i]/JyieldSub_low_la[i]*sqrt(TMath::Power(V2Values_err_low_la[i]/V2Values_low_la[i],2) + TMath::Power(JyieldSub_err_la[i]/JyieldSub_la[i],2) + TMath::Power(JyieldSub_err_low_la[i]/JyieldSub_low_la[i],2))*V2Values_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*JyieldSub_la[i]/JyieldSub_low_la[i]));

        V3sub_la.push_back(V3Values_la[i] - V3Values_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*JyieldSub_la[i]/JyieldSub_low_la[i]);
        V3sube_la.push_back(sqrt(TMath::Power(V3Values_err_la[i],2) + sqrt(TMath::Power(V3Values_err_low_la[i]/V3Values_low_la[i],2) + TMath::Power(JyieldSub_err_la[i]/JyieldSub_la[i],2) + TMath::Power(JyieldSub_err_low_la[i]/JyieldSub_low_la[i],2))*V3Values_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*JyieldSub_la[i]/JyieldSub_low_la[i]*sqrt(TMath::Power(V3Values_err_low_la[i]/V3Values_low_la[i],2) + TMath::Power(JyieldSub_err_la[i]/JyieldSub_la[i],2) + TMath::Power(JyieldSub_err_low_la[i]/JyieldSub_low_la[i],2))*V3Values_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*JyieldSub_la[i]/JyieldSub_low_la[i]));

        jetYfactor_la.push_back(JyieldSub_la[i]/JyieldSub_low_la[i]);
    }


    //Ks vn calculation
    std::vector<double> v2sub_ks;
    std::vector<double> v2sube_ks;
    std::vector<double> v2_ks;
    std::vector<double> v2e_ks;
    std::vector<double> v2_low_ks;
    std::vector<double> v2e_low_ks;

    std::vector<double> v3sub_ks;
    std::vector<double> v3sube_ks;
    std::vector<double> v3_ks;
    std::vector<double> v3e_ks;
    std::vector<double> v3_low_ks;
    std::vector<double> v3e_low_ks;

    for(unsigned i=0; i<PD.PtBin_ks.size(); i++){
        v2sub_ks.push_back(V2sub_ks[i]/v2sub_ref);
        v2sube_ks.push_back(fabs(sqrt(V2sube_ks[i]/V2sub_ks[i]*V2sube_ks[i]/V2sub_ks[i] + v2sube_ref/v2sub_ref*v2sube_ref/v2sub_ref)*v2sub_ks[i]));

        v2_ks.push_back(V2Values_ks[i]/v2_high_ref);
        v2e_ks.push_back(sqrt(TMath::Power(V2Values_err_ks[i]/v2_high_ref,2) + TMath::Power(v2e_high_ref/v2_high_ref,2)*v2_ks[i]));

        v2_low_ks.push_back(V2Values_low_ks[i]/v2_low_ref);
        v2e_low_ks.push_back(sqrt(TMath::Power(V2Values_err_low_ks[i]/v2_low_ref,2) + TMath::Power(v2e_low_ref/v2_low_ref,2))*v2_low_ks[i]);

        v3sub_ks.push_back(V3sub_ks[i]/v3sub_ref);
        v3sube_ks.push_back(fabs(sqrt(V3sube_ks[i]/V3sub_ks[i]*V3sube_ks[i]/V3sub_ks[i] + v3sube_ref/v3sub_ref*v3sube_ref/v3sub_ref)*v3sub_ks[i]));

        v3_ks.push_back(V3Values_ks[i]/v3_high_ref);
        v3e_ks.push_back(sqrt(TMath::Power(V3Values_err_ks[i]/v3_high_ref,2) + TMath::Power(v3e_high_ref/v3_high_ref,2)*v3_ks[i]));

        v3_low_ks.push_back(V3Values_low_ks[i]/v3_low_ref);
        v3e_low_ks.push_back(sqrt(TMath::Power(V3Values_err_low_ks[i]/v3_low_ref,2) + TMath::Power(v3e_low_ref/v3_low_ref,2))*v3_low_ks[i]);
        }

    //La vn calculation
    std::vector<double> v2sub_la;
    std::vector<double> v2sube_la;
    std::vector<double> v2_la;
    std::vector<double> v2e_la;
    std::vector<double> v2_low_la;
    std::vector<double> v2e_low_la;

    std::vector<double> v3sub_la;
    std::vector<double> v3sube_la;
    std::vector<double> v3_la;
    std::vector<double> v3e_la;
    std::vector<double> v3_low_la;
    std::vector<double> v3e_low_la;

    for(unsigned i=0; i<PD.PtBin_la.size(); i++){
        v2sub_la.push_back(V2sub_la[i]/v2sub_ref);
        v2sube_la.push_back(fabs(sqrt(V2sube_la[i]/V2sub_la[i]*V2sube_la[i]/V2sub_la[i] + v2sube_ref/v2sub_ref*v2sube_ref/v2sub_ref)*v2sub_la[i]));

        v2_la.push_back(V2Values_la[i]/v2_high_ref);
        v2e_la.push_back(sqrt(TMath::Power(V2Values_err_la[i]/v2_high_ref,2) + TMath::Power(v2e_high_ref/v2_high_ref,2)*v2_la[i]));

        v2_low_la.push_back(V2Values_low_la[i]/v2_low_ref);
        v2e_low_la.push_back(sqrt(TMath::Power(V2Values_err_low_la[i]/v2_low_ref,2) + TMath::Power(v2e_low_ref/v2_low_ref,2))*v2_low_la[i]);

        v3sub_la.push_back(V3sub_la[i]/v3sub_ref);
        v3sube_la.push_back(fabs(sqrt(V3sube_la[i]/V3sub_la[i]*V3sube_la[i]/V3sub_la[i] + v3sube_ref/v3sub_ref*v3sube_ref/v3sub_ref)*v3sub_la[i]));

        v3_la.push_back(V3Values_la[i]/v3_high_ref);
        v3e_la.push_back(sqrt(TMath::Power(V3Values_err_la[i]/v3_high_ref,2) + TMath::Power(v3e_high_ref/v3_high_ref,2)*v3_la[i]));

       v3_low_la.push_back(V3Values_low_la[i]/v3_low_ref);
       v3e_low_la.push_back(sqrt(TMath::Power(V3Values_err_low_la[i]/v3_low_ref,2) + TMath::Power(v3e_low_ref/v3_low_ref,2))*v3_low_la[i]);
    }

    //Produce RootFiles
    TFile output(PD.fn.c_str(),"UPDATE");

    //Get Mean Pt and KET of High N
    std::vector<double> pt_ks;
    std::vector<double> Ket_ks;
    for(unsigned i=0; i<PD.PtBin_ks.size(); i++){
        TH1D* hpt = (TH1D*)PD.f_V0->Get(Form((PD.fn_V0 + "/Ptkshort_pt%d").c_str(),i));
        TH1D* hket = (TH1D*)PD.f_V0->Get(Form((PD.fn_V0 + "/KETkshort_pt%d").c_str(),i));

        pt_ks.push_back(hpt->GetMean(1));
        Ket_ks.push_back(hket->GetMean(1));
    }

    std::vector<double> pt_la;
    std::vector<double> Ket_la;
    for(unsigned i=0; i<PD.PtBin_la.size(); i++){
        TH1D* hpt = (TH1D*)PD.f_V0->Get(Form((PD.fn_V0 + "/Ptlambda_pt%d").c_str(),i));
        TH1D* hket = (TH1D*)PD.f_V0->Get(Form((PD.fn_V0 + "/KETlambda_pt%d").c_str(),i));

        pt_la.push_back(hpt->GetMean(1));
        Ket_la.push_back(hket->GetMean(1));
    }

    std::vector<double> perisubfactor_ks;
    std::vector<double> perisubfactore_ks;
    for(int i=0; i<arraySize_ks; i++){
        perisubfactor_ks.push_back(V2Values_low_ks[i]*Nassoc_low_ks[i]/JyieldSub_low_ks[i]);
        perisubfactore_ks.push_back(sqrt(TMath::Power(V2Values_err_low_ks[i]/V2Values_low_ks[i],2) + TMath::Power(JyieldSub_err_low_ks[i]/JyieldSub_low_ks[i],2))*perisubfactor_ks[i]);
    }
    std::vector<double> perisubfactor_la;
    std::vector<double> perisubfactore_la;
    for(int i=0; i<arraySize_la; i++){
        perisubfactor_la.push_back(V2Values_low_la[i]*Nassoc_low_la[i]/JyieldSub_low_la[i]);
        perisubfactore_la.push_back(sqrt(TMath::Power(V2Values_err_low_la[i]/V2Values_low_la[i],2) + TMath::Power(JyieldSub_err_low_la[i]/JyieldSub_low_la[i],2))*perisubfactor_la[i]);
    }

    

    //Obs v2 values
    TGraphErrors* v2plot_ks        = new TGraphErrors(arraySize_ks,&pt_ks[0] ,&v2_ks[0]      ,0,&v2e_ks[0]);
    TGraphErrors* v2plot_KET_ks    = new TGraphErrors(arraySize_ks,&Ket_ks[0],&v2_ks[0]      ,0,&v2e_ks[0]);
    TGraphErrors* v2subplot_ks     = new TGraphErrors(arraySize_ks,&pt_ks[0] ,&v2sub_ks[0]   ,0,&v2sube_ks[0]);
    TGraphErrors* v2subplot_KET_ks = new TGraphErrors(arraySize_ks,&Ket_ks[0],&v2sub_ks[0]   ,0,&v2sube_ks[0]);
    TGraphErrors* v3plot_ks        = new TGraphErrors(arraySize_ks,&pt_ks[0] ,&v3_ks[0]      ,0,&v3e_ks[0]);
    TGraphErrors* v3plot_KET_ks    = new TGraphErrors(arraySize_ks,&Ket_ks[0],&v3_ks[0]      ,0,&v3e_ks[0]);
    TGraphErrors* v3subplot_ks     = new TGraphErrors(arraySize_ks,&pt_ks[0] ,&v3sub_ks[0]   ,0,&v3sube_ks[0]);
    TGraphErrors* v3subplot_KET_ks = new TGraphErrors(arraySize_ks,&Ket_ks[0],&v3sub_ks[0]   ,0,&v3sube_ks[0]);

    TGraphErrors* v2plot_la        = new TGraphErrors(arraySize_la,&pt_la[0] ,&v2_la[0]      ,0,&v2e_la[0]);
    TGraphErrors* v2plot_KET_la    = new TGraphErrors(arraySize_la,&Ket_la[0],&v2_la[0]      ,0,&v2e_la[0]);
    TGraphErrors* v2subplot_la     = new TGraphErrors(arraySize_la,&pt_la[0] ,&v2sub_la[0]   ,0,&v2sube_la[0]);
    TGraphErrors* v2subplot_KET_la = new TGraphErrors(arraySize_la,&Ket_la[0],&v2sub_la[0]   ,0,&v2sube_la[0]);
    TGraphErrors* v3plot_la        = new TGraphErrors(arraySize_la,&pt_la[0] ,&v3_la[0]      ,0,&v3e_la[0]);
    TGraphErrors* v3plot_KET_la    = new TGraphErrors(arraySize_la,&Ket_la[0],&v3_la[0]      ,0,&v3e_la[0]);
    TGraphErrors* v3subplot_la     = new TGraphErrors(arraySize_la,&pt_la[0] ,&v3sub_la[0]   ,0,&v3sube_la[0]);
    TGraphErrors* v3subplot_KET_la = new TGraphErrors(arraySize_la,&Ket_la[0],&v3sub_la[0]   ,0,&v3sube_la[0]);

    //Obs V2 values
    TGraphErrors* V2plot_ks        = new TGraphErrors(arraySize_ks,&pt_ks[0] ,&V2Values_ks[0],0,&V2Values_err_ks[0]);
    TGraphErrors* V2plot_KET_ks    = new TGraphErrors(arraySize_ks,&Ket_ks[0],&V2Values_ks[0],0,&V2Values_err_ks[0]);
    TGraphErrors* V2subplot_ks     = new TGraphErrors(arraySize_ks,&pt_ks[0] ,&V2sub_ks[0]   ,0,&V2sube_ks[0]);
    TGraphErrors* V2subplot_KET_ks = new TGraphErrors(arraySize_ks,&Ket_ks[0],&V2sub_ks[0]   ,0,&V2sube_ks[0]);
    TGraphErrors* V3plot_ks        = new TGraphErrors(arraySize_ks,&pt_ks[0] ,&V3Values_ks[0],0,&V3Values_err_ks[0]);
    TGraphErrors* V3plot_KET_ks    = new TGraphErrors(arraySize_ks,&Ket_ks[0],&V3Values_ks[0],0,&V3Values_err_ks[0]);
    TGraphErrors* V3subplot_ks     = new TGraphErrors(arraySize_ks,&pt_ks[0] ,&V3sub_ks[0]   ,0,&V3sube_ks[0]);
    TGraphErrors* V3subplot_KET_ks = new TGraphErrors(arraySize_ks,&Ket_ks[0],&V3sub_ks[0]   ,0,&V3sube_ks[0]);

    TGraphErrors* V2plot_la        = new TGraphErrors(arraySize_la,&pt_la[0] ,&V2Values_la[0],0,&V2Values_err_la[0]);
    TGraphErrors* V2plot_KET_la    = new TGraphErrors(arraySize_la,&Ket_la[0],&V2Values_la[0],0,&V2Values_err_la[0]);
    TGraphErrors* V2subplot_la     = new TGraphErrors(arraySize_la,&pt_la[0] ,&V2sub_la[0]   ,0,&V2sube_la[0]);
    TGraphErrors* V2subplot_KET_la = new TGraphErrors(arraySize_la,&Ket_la[0],&V2sub_la[0]   ,0,&V2sube_la[0]);
    TGraphErrors* V3plot_la        = new TGraphErrors(arraySize_la,&pt_la[0] ,&V3Values_la[0],0,&V3Values_err_la[0]);
    TGraphErrors* V3plot_KET_la    = new TGraphErrors(arraySize_la,&Ket_la[0],&V3Values_la[0],0,&V3Values_err_la[0]);
    TGraphErrors* V3subplot_la     = new TGraphErrors(arraySize_la,&pt_la[0] ,&V3sub_la[0]   ,0,&V3sube_la[0]);
    TGraphErrors* V3subplot_KET_la = new TGraphErrors(arraySize_la,&Ket_la[0],&V3sub_la[0]   ,0,&V3sube_la[0]);

    //Yields
    TGraphErrors* srYieldPlot_ks       = new TGraphErrors(arraySize_ks,&pt_ks[0],&Jyieldsr_ks[0]     ,0,&Jyieldsr_err_ks[0]);
    TGraphErrors* lrYieldPlot_ks       = new TGraphErrors(arraySize_ks,&pt_ks[0],&Jyieldlr_ks[0]     ,0,&Jyieldlr_err_ks[0]);
    TGraphErrors* subYieldPlot_ks      = new TGraphErrors(arraySize_ks,&pt_ks[0],&JyieldSub_ks[0]    ,0,&JyieldSub_err_ks[0]);
    TGraphErrors* perisubfactorplot_ks = new TGraphErrors(arraySize_ks,&pt_ks[0],&perisubfactor_ks[0],0,&perisubfactore_ks[0]);

    TGraphErrors* srYieldPlot_la       = new TGraphErrors(arraySize_la,&pt_la[0],&Jyieldsr_la[0]     ,0,&Jyieldsr_err_la[0]);
    TGraphErrors* lrYieldPlot_la       = new TGraphErrors(arraySize_la,&pt_la[0],&Jyieldlr_la[0]     ,0,&Jyieldlr_err_la[0]);
    TGraphErrors* subYieldPlot_la      = new TGraphErrors(arraySize_la,&pt_la[0],&JyieldSub_la[0]    ,0,&JyieldSub_err_la[0]);
    TGraphErrors* perisubfactorplot_la = new TGraphErrors(arraySize_la,&pt_la[0],&perisubfactor_la[0],0,&perisubfactore_la[0]);

    //For Direct subtraction method
    TGraphErrors* Nasslow_ks = new TGraphErrors(arraySize_ks,&pt_ks[0],&Nassoc_low_ks[0],0,0);
    TGraphErrors* Nasslow_la = new TGraphErrors(arraySize_la,&pt_la[0],&Nassoc_low_la[0],0,0);
    TGraphErrors* Nass_ks = new TGraphErrors(arraySize_ks,&pt_ks[0],&Nassoc_ks[0],0,0);
    TGraphErrors* Nass_la = new TGraphErrors(arraySize_la,&pt_la[0],&Nassoc_la[0],0,0);
    TGraphErrors* subYieldPlot_low_ks = new TGraphErrors(arraySize_ks,&pt_ks[0],&JyieldSub_low_ks[0],0,&JyieldSub_err_low_ks[0]);
    TGraphErrors* subYieldPlot_low_la = new TGraphErrors(arraySize_la,&pt_la[0],&JyieldSub_low_la[0],0,&JyieldSub_err_low_la[0]);

    TGraphErrors* bz_ks = new TGraphErrors(arraySize_ks,&pt_ks[0],&Bz_ks[0],0,0);
    TGraphErrors* bz_la = new TGraphErrors(arraySize_la,&pt_la[0],&Bz_la[0],0,0);

    v2plot_ks       ->Write("v2plot_bkg_ks"       ,TObject::kOverwrite);
    v2plot_KET_ks   ->Write("v2plot_KET_bkg_ks"   ,TObject::kOverwrite);
    v2subplot_ks    ->Write("v2subplot_bkg_ks"    ,TObject::kOverwrite);
    v2subplot_KET_ks->Write("v2subplot_KET_bkg_ks",TObject::kOverwrite);
    v3plot_ks       ->Write("v3plot_bkg_ks"       ,TObject::kOverwrite);
    v3plot_KET_ks   ->Write("v3plot_KET_bkg_ks"   ,TObject::kOverwrite);
    v3subplot_ks    ->Write("v3subplot_bkg_ks"    ,TObject::kOverwrite);
    v3subplot_KET_ks->Write("v3subplot_KET_bkg_ks",TObject::kOverwrite);

    v2plot_la       ->Write("v2plot_bkg_la"       ,TObject::kOverwrite);
    v2plot_KET_la   ->Write("v2plot_KET_bkg_la"   ,TObject::kOverwrite);
    v2subplot_la    ->Write("v2subplot_bkg_la"    ,TObject::kOverwrite);
    v2subplot_KET_la->Write("v2subplot_KET_bkg_la",TObject::kOverwrite);
    v3plot_la       ->Write("v3plot_bkg_la"       ,TObject::kOverwrite);
    v3plot_KET_la   ->Write("v3plot_KET_bkg_la"   ,TObject::kOverwrite);
    v3subplot_la    ->Write("v3subplot_bkg_la"    ,TObject::kOverwrite);
    v3subplot_KET_la->Write("v3subplot_KET_bkg_la",TObject::kOverwrite);

    V2plot_ks       ->Write("V2plot_bkg_ks"       ,TObject::kOverwrite);
    V2plot_KET_ks   ->Write("V2plot_KET_bkg_ks"   ,TObject::kOverwrite);
    V2subplot_ks    ->Write("V2subplot_bkg_ks"    ,TObject::kOverwrite);
    V2subplot_KET_ks->Write("V2subplot_KET_bkg_ks",TObject::kOverwrite);
    V3plot_ks       ->Write("V3plot_bkg_ks"       ,TObject::kOverwrite);
    V3plot_KET_ks   ->Write("V3plot_KET_bkg_ks"   ,TObject::kOverwrite);
    V3subplot_ks    ->Write("V3subplot_bkg_ks"    ,TObject::kOverwrite);
    V3subplot_KET_ks->Write("V3subplot_KET_bkg_ks",TObject::kOverwrite);

    V2plot_la       ->Write("V2plot_bkg_la"       ,TObject::kOverwrite);
    V2plot_KET_la   ->Write("V2plot_KET_bkg_la"   ,TObject::kOverwrite);
    V2subplot_la    ->Write("V2subplot_bkg_la"    ,TObject::kOverwrite);
    V2subplot_KET_la->Write("V2subplot_KET_bkg_la",TObject::kOverwrite);
    V3plot_la       ->Write("V3plot_bkg_la"       ,TObject::kOverwrite);
    V3plot_KET_la   ->Write("V3plot_KET_bkg_la"   ,TObject::kOverwrite);
    V3subplot_la    ->Write("V3subplot_bkg_la"    ,TObject::kOverwrite);
    V3subplot_KET_la->Write("V3subplot_KET_bkg_la",TObject::kOverwrite);

    perisubfactorplot_ks->Write("perisubfactor_bkg_ks",TObject::kOverwrite);
    perisubfactorplot_la->Write("perisubfactor_bkg_la",TObject::kOverwrite);

    Nass_ks        ->Write("Nassoc_bkg_ks"      ,TObject::kOverwrite);
    Nasslow_ks     ->Write("Nassoc_low_bkg_ks"  ,TObject::kOverwrite);
    Nass_la        ->Write("Nassoc_bkg_la"      ,TObject::kOverwrite);
    Nasslow_la     ->Write("Nassoc_low_bkg_la"  ,TObject::kOverwrite);
    subYieldPlot_low_ks->Write("YieldPlot_low_bkg_ks",TObject::kOverwrite);
    subYieldPlot_low_la->Write("YieldPlot_low_bkg_la",TObject::kOverwrite);
    subYieldPlot_ks->Write("YieldPlot_bkg_ks"   ,TObject::kOverwrite);
    subYieldPlot_la->Write("YieldPlot_bkg_la"   ,TObject::kOverwrite);

    bz_ks ->Write("bz_bkg_ks",TObject::kOverwrite);
    bz_la ->Write("bz_bkg_la",TObject::kOverwrite);
}
*/

void PeriSubSigIndirect(ParticleData PD)
{
    TH1::SetDefaultSumw2();

    TFile* f = TFile::Open(PD.fn.c_str());

    TGraphErrors* v2obs_ks;
    TGraphErrors* v2obs_sub_ks;
    TGraphErrors* v2bkg_ks;
    TGraphErrors* v2bkg_sub_ks;

    TGraphErrors* v2obs_KET_ks;
    TGraphErrors* v2obs_KET_sub_ks;
    TGraphErrors* v2bkg_KET_ks;
    TGraphErrors* v2bkg_KET_sub_ks;

    TGraphErrors* v2obs_la;
    TGraphErrors* v2obs_sub_la;
    TGraphErrors* v2bkg_la;
    TGraphErrors* v2bkg_sub_la;

    TGraphErrors* v2obs_KET_la;
    TGraphErrors* v2obs_KET_sub_la;
    TGraphErrors* v2bkg_KET_la;
    TGraphErrors* v2bkg_KET_sub_la;

    TGraphErrors* v2obs_xi;
    TGraphErrors* v2obs_sub_xi;
    TGraphErrors* v2bkg_xi;
    TGraphErrors* v2bkg_sub_xi;

    TGraphErrors* v2obs_KET_xi;
    TGraphErrors* v2obs_KET_sub_xi;
    TGraphErrors* v2bkg_KET_xi;
    TGraphErrors* v2bkg_KET_sub_xi;

    f->GetObject("v2plot_obs_ks"   ,v2obs_ks);
    f->GetObject("v2subplot_obs_ks",v2obs_sub_ks);
    f->GetObject("v2plot_bkg_ks"   ,v2bkg_ks);
    f->GetObject("v2subplot_bkg_ks",v2bkg_sub_ks);

    f->GetObject("v2plot_KET_obs_ks"   ,v2obs_KET_ks);
    f->GetObject("v2subplot_KET_obs_ks",v2obs_KET_sub_ks);
    f->GetObject("v2plot_KET_bkg_ks"   ,v2bkg_KET_ks);
    f->GetObject("v2subplot_KET_bkg_ks",v2bkg_KET_sub_ks);

    f->GetObject("v2plot_obs_la"   ,v2obs_la);
    f->GetObject("v2subplot_obs_la",v2obs_sub_la);
    f->GetObject("v2plot_bkg_la"   ,v2bkg_la);
    f->GetObject("v2subplot_bkg_la",v2bkg_sub_la);

    f->GetObject("v2plot_KET_obs_la"   ,v2obs_KET_la);
    f->GetObject("v2subplot_KET_obs_la",v2obs_KET_sub_la);
    f->GetObject("v2plot_KET_bkg_la"   ,v2bkg_KET_la);
    f->GetObject("v2subplot_KET_bkg_la",v2bkg_KET_sub_la);

    f->GetObject("v2plot_obs_xi"   ,v2obs_xi);
    f->GetObject("v2subplot_obs_xi",v2obs_sub_xi);
    f->GetObject("v2plot_bkg_xi"   ,v2bkg_xi);
    f->GetObject("v2subplot_bkg_xi",v2bkg_sub_xi);

    f->GetObject("v2plot_KET_obs_xi"   ,v2obs_KET_xi);
    f->GetObject("v2subplot_KET_obs_xi",v2obs_KET_sub_xi);
    f->GetObject("v2plot_KET_bkg_xi"   ,v2bkg_KET_xi);
    f->GetObject("v2subplot_KET_bkg_xi",v2bkg_KET_sub_xi);

    double* v2_obs_ks     = v2obs_ks->GetY();
    double* v2_obs_err_ks = v2obs_ks->GetEY();
    double* v2_bkg_ks     = v2bkg_ks->GetY();
    double* v2_bkg_err_ks = v2bkg_ks->GetEY();

    double* v2_obs_sub_ks     = v2obs_sub_ks->GetY();
    double* v2_obs_err_sub_ks = v2obs_sub_ks->GetEY();
    double* v2_bkg_sub_ks     = v2bkg_sub_ks->GetY();
    double* v2_bkg_err_sub_ks = v2bkg_sub_ks->GetEY();

    double* v2_obs_la     = v2obs_la->GetY();
    double* v2_obs_err_la = v2obs_la->GetEY();
    double* v2_bkg_la     = v2bkg_la->GetY();
    double* v2_bkg_err_la = v2bkg_la->GetEY();

    double* v2_obs_sub_la     = v2obs_sub_la->GetY();
    double* v2_obs_err_sub_la = v2obs_sub_la->GetEY();
    double* v2_bkg_sub_la     = v2bkg_sub_la->GetY();
    double* v2_bkg_err_sub_la = v2bkg_sub_la->GetEY();

    double* v2_obs_xi     = v2obs_xi->GetY();
    double* v2_obs_err_xi = v2obs_xi->GetEY();
    double* v2_bkg_xi     = v2bkg_xi->GetY();
    double* v2_bkg_err_xi = v2bkg_xi->GetEY();

    double* v2_obs_sub_xi     = v2obs_sub_xi->GetY();
    double* v2_obs_err_sub_xi = v2obs_sub_xi->GetEY();
    double* v2_bkg_sub_xi     = v2bkg_sub_xi->GetY();
    double* v2_bkg_err_sub_xi = v2bkg_sub_xi->GetEY();

    double* pt_ks = v2obs_ks->GetX();
    double* pt_la = v2obs_la->GetX();
    double* pt_xi = v2obs_xi->GetX();

    double* KET_ks = v2obs_KET_ks->GetX();
    double* KET_la = v2obs_KET_la->GetX();
    double* KET_xi = v2obs_KET_xi->GetX();

    std::vector<double> bkgfrac_ks;
    std::vector<double> v2true_ks;
    std::vector<double> v2true_err_ks;
    std::vector<double> v2true_sub_ks;
    std::vector<double> v2true_sub_err_ks;

    std::vector<double> bkgfrac_la;
    std::vector<double> v2true_la;
    std::vector<double> v2true_err_la;
    std::vector<double> v2true_sub_la;
    std::vector<double> v2true_sub_err_la;

    std::vector<double> bkgfrac_xi;
    std::vector<double> v2true_xi;
    std::vector<double> v2true_err_xi;
    std::vector<double> v2true_sub_xi;
    std::vector<double> v2true_sub_err_xi;

    for(unsigned i=0; i<PD.PtBin_ks.size(); i++)
    {
        bkgfrac_ks.push_back(1 - PD.fsig_ks[i]);
        v2true_ks.push_back((v2_obs_ks[i] - v2_bkg_ks[i]*bkgfrac_ks[i])/PD.fsig_ks[i]);
        v2true_err_ks.push_back(sqrt(TMath::Power((v2_bkg_err_ks[i]*bkgfrac_ks[i]),2) + TMath::Power(v2_obs_err_ks[i],2))/PD.fsig_ks[i]);
        v2true_sub_ks.push_back((v2_obs_sub_ks[i] - v2_bkg_sub_ks[i]*bkgfrac_ks[i])/PD.fsig_ks[i]);
        v2true_sub_err_ks.push_back(sqrt(TMath::Power(v2_bkg_err_sub_ks[i]*bkgfrac_ks[i],2) + TMath::Power(v2_obs_err_sub_ks[i],2))/PD.fsig_ks[i]);
    }

    for(unsigned i=0; i<PD.PtBin_la.size(); i++)
    {
        bkgfrac_la.push_back(1 - PD.fsig_la[i]);
        v2true_la.push_back((v2_obs_la[i] - v2_bkg_la[i]*bkgfrac_la[i])/PD.fsig_la[i]);
        v2true_err_la.push_back(sqrt(TMath::Power((v2_bkg_err_la[i]*bkgfrac_la[i]),2) + TMath::Power(v2_obs_err_la[i],2))/PD.fsig_la[i]);
        v2true_sub_la.push_back((v2_obs_sub_la[i] - v2_bkg_sub_la[i]*bkgfrac_la[i])/PD.fsig_la[i]);
        v2true_sub_err_la.push_back(sqrt(TMath::Power(v2_bkg_err_sub_la[i]*bkgfrac_la[i],2) + TMath::Power(v2_obs_err_sub_la[i],2))/PD.fsig_la[i]);
    }

    for(unsigned i=0; i<PD.PtBin_xi.size(); i++)
    {
        bkgfrac_xi.push_back(1 - PD.fsig_xi[i]);
        v2true_xi.push_back((v2_obs_xi[i] - v2_bkg_xi[i]*bkgfrac_xi[i])/PD.fsig_xi[i]);
        v2true_err_xi.push_back(sqrt(TMath::Power((v2_bkg_err_xi[i]*bkgfrac_xi[i]),2) + TMath::Power(v2_obs_err_xi[i],2))/PD.fsig_xi[i]);
        v2true_sub_xi.push_back((v2_obs_sub_xi[i] - v2_bkg_sub_xi[i]*bkgfrac_xi[i])/PD.fsig_xi[i]);
        v2true_sub_err_xi.push_back(sqrt(TMath::Power(v2_bkg_err_sub_xi[i]*bkgfrac_xi[i],2) + TMath::Power(v2_obs_err_sub_xi[i],2))/PD.fsig_xi[i]);
    }

    int numPtBins_ks = PD.PtBin_ks.size()-1;
    int numPtBins_la = PD.PtBin_la.size()-1;
    int numPtBins_xi = PD.PtBin_xi.size()-1;

    TGraphErrors* ksv2true    = new TGraphErrors(numPtBins_ks,pt_ks,&v2true_ks[0]    ,0,&v2true_err_ks[0]);
    TGraphErrors* ksv2truesub = new TGraphErrors(numPtBins_ks,pt_ks,&v2true_sub_ks[0],0,&v2true_sub_err_ks[0]);

    TGraphErrors* ksv2true_KET    = new TGraphErrors(numPtBins_ks,KET_ks,&v2true_ks[0]    ,0,&v2true_err_ks[0]);
    TGraphErrors* ksv2truesub_KET = new TGraphErrors(numPtBins_ks,KET_ks,&v2true_sub_ks[0],0,&v2true_sub_err_ks[0]);

    TGraphErrors* lav2true    = new TGraphErrors(numPtBins_la,pt_la,&v2true_la[0]    ,0,&v2true_err_la[0]);
    TGraphErrors* lav2truesub = new TGraphErrors(numPtBins_la,pt_la,&v2true_sub_la[0],0,&v2true_sub_err_la[0]);

    TGraphErrors* lav2true_KET    = new TGraphErrors(numPtBins_la,KET_la,&v2true_la[0]    ,0,&v2true_err_la[0]);
    TGraphErrors* lav2truesub_KET = new TGraphErrors(numPtBins_la,KET_la,&v2true_sub_la[0],0,&v2true_sub_err_la[0]);

    TGraphErrors* xiv2true    = new TGraphErrors(numPtBins_xi,pt_xi,&v2true_xi[0]    ,0,&v2true_err_xi[0]);
    TGraphErrors* xiv2truesub = new TGraphErrors(numPtBins_xi,pt_xi,&v2true_sub_xi[0],0,&v2true_sub_err_xi[0]);

    TGraphErrors* xiv2true_KET    = new TGraphErrors(numPtBins_xi,KET_xi,&v2true_xi[0]    ,0,&v2true_err_xi[0]);
    TGraphErrors* xiv2truesub_KET = new TGraphErrors(numPtBins_xi,KET_xi,&v2true_sub_xi[0],0,&v2true_sub_err_xi[0]);

    TFile output(PD.fn.c_str(),"UPDATE");

    ksv2true       ->Write("kshortv2true"       ,TObject::kOverwrite);
    ksv2true_KET   ->Write("kshortv2true_KET"   ,TObject::kOverwrite);
    ksv2truesub    ->Write("kshortv2truesub"    ,TObject::kOverwrite);
    ksv2truesub_KET->Write("kshortv2truesub_KET",TObject::kOverwrite);

    lav2true       ->Write("lambdav2true"       ,TObject::kOverwrite);
    lav2truesub    ->Write("lambdav2truesub"    ,TObject::kOverwrite);
    lav2true_KET   ->Write("lambdav2true_KET"   ,TObject::kOverwrite);
    lav2truesub_KET->Write("lambdav2truesub_KET",TObject::kOverwrite);

    xiv2true       ->Write("xiv2true"       ,TObject::kOverwrite);
    xiv2truesub    ->Write("xiv2truesub"    ,TObject::kOverwrite);
    xiv2true_KET   ->Write("xiv2true_KET"   ,TObject::kOverwrite);
    xiv2truesub_KET->Write("xiv2truesub_KET",TObject::kOverwrite);
}

void PeriSubSigDirect(ParticleData PD)
{
    TH1::SetDefaultSumw2();

    TFile* f = TFile::Open(PD.fn.c_str());
    TFile* f_MB = TFile::Open("/Volumes/MacHD/Users/blt1/research/Macros/Drawers/TGraph/rootFiles/V0Casv2LowMult.root");

    TGraphErrors* ksv2true;
    TGraphErrors* ksv2true_KET;
    TGraphErrors* ksv2_MB;

    TGraphErrors* lav2true;
    TGraphErrors* lav2true_KET;
    TGraphErrors* lav2_MB;

    TGraphErrors* xiv2true;
    TGraphErrors* xiv2true_KET;
    TGraphErrors* xiv2_MB;

    TGraphErrors* Nass_obs_ks;
    TGraphErrors* Nass_low_obs_ks;
    TGraphErrors* subYieldPlot_obs_ks;
    TGraphErrors* subYieldPlot_low_obs_ks;

    TGraphErrors* Nass_bkg_ks;
    TGraphErrors* Nass_low_bkg_ks;
    TGraphErrors* subYieldPlot_bkg_ks;
    TGraphErrors* subYieldPlot_low_bkg_ks;

    TGraphErrors* Nass_obs_la;
    TGraphErrors* Nass_low_obs_la;
    TGraphErrors* subYieldPlot_obs_la;
    TGraphErrors* subYieldPlot_low_obs_la;

    TGraphErrors* Nass_low_bkg_la;
    TGraphErrors* Nass_bkg_la;
    TGraphErrors* subYieldPlot_bkg_la;
    TGraphErrors* subYieldPlot_low_bkg_la;

    TGraphErrors* Nass_obs_xi;
    TGraphErrors* Nass_low_obs_xi;
    TGraphErrors* subYieldPlot_obs_xi;
    TGraphErrors* subYieldPlot_low_obs_xi;

    TGraphErrors* Nass_low_bkg_xi;
    TGraphErrors* Nass_bkg_xi;
    TGraphErrors* subYieldPlot_bkg_xi;
    TGraphErrors* subYieldPlot_low_bkg_xi;

    f->GetObject("kshortv2true", ksv2true);
    f->GetObject("kshortv2true_KET", ksv2true_KET);
    f_MB->GetObject("v2kshort", ksv2_MB);

    f->GetObject("lambdav2true", lav2true);
    f->GetObject("lambdav2true_KET", lav2true_KET);
    f_MB->GetObject("v2lambda", lav2_MB);

    f->GetObject("xiv2true", xiv2true);
    f->GetObject("xiv2true_KET", xiv2true_KET);
    f_MB->GetObject("v2xi", xiv2_MB);

    f->GetObject("Nassoc_obs_ks", Nass_obs_ks);
    f->GetObject("Nassoc_low_obs_ks", Nass_low_obs_ks);
    f->GetObject("YieldPlot_low_obs_ks", subYieldPlot_low_obs_ks);
    f->GetObject("YieldPlot_obs_ks", subYieldPlot_obs_ks);

    f->GetObject("Nassoc_bkg_ks", Nass_bkg_ks);
    f->GetObject("Nassoc_low_bkg_ks", Nass_low_bkg_ks);
    f->GetObject("YieldPlot_low_bkg_ks", subYieldPlot_low_bkg_ks);
    f->GetObject("YieldPlot_bkg_ks", subYieldPlot_bkg_ks);

    f->GetObject("Nassoc_obs_la", Nass_obs_la);
    f->GetObject("Nassoc_low_obs_la", Nass_low_obs_la);
    f->GetObject("YieldPlot_low_obs_la", subYieldPlot_low_obs_la);
    f->GetObject("YieldPlot_obs_la", subYieldPlot_obs_la);

    f->GetObject("Nassoc_bkg_la", Nass_bkg_la);
    f->GetObject("Nassoc_low_bkg_la", Nass_low_bkg_la);
    f->GetObject("YieldPlot_low_bkg_la", subYieldPlot_low_bkg_la);
    f->GetObject("YieldPlot_bkg_la", subYieldPlot_bkg_la);

    f->GetObject("Nassoc_obs_xi", Nass_obs_xi);
    f->GetObject("Nassoc_low_obs_xi", Nass_low_obs_xi);
    f->GetObject("YieldPlot_low_obs_xi", subYieldPlot_low_obs_xi);
    f->GetObject("YieldPlot_obs_xi", subYieldPlot_obs_xi);

    f->GetObject("Nassoc_bkg_xi", Nass_bkg_xi);
    f->GetObject("Nassoc_low_bkg_xi", Nass_low_bkg_xi);
    f->GetObject("YieldPlot_low_bkg_xi", subYieldPlot_low_bkg_xi);
    f->GetObject("YieldPlot_bkg_xi", subYieldPlot_bkg_xi);

    double* Nassoc_obs_ks            = Nass_obs_ks->GetY();
    double* Nassoc_low_obs_ks        = Nass_low_obs_ks->GetY();
    double* YieldPlot_obs_ks         = subYieldPlot_obs_ks->GetY();
    double* YieldPlot_obs_err_ks     = subYieldPlot_obs_ks->GetEY();
    double* YieldPlot_low_obs_ks     = subYieldPlot_low_obs_ks->GetY();
    double* YieldPlot_low_obs_err_ks = subYieldPlot_low_obs_ks->GetEY();

    double* Nassoc_bkg_ks            = Nass_bkg_ks->GetY();
    double* Nassoc_low_bkg_ks        = Nass_low_bkg_ks->GetY();
    double* YieldPlot_bkg_ks         = subYieldPlot_bkg_ks->GetY();
    double* YieldPlot_bkg_err_ks     = subYieldPlot_bkg_ks->GetEY();
    double* YieldPlot_low_bkg_ks     = subYieldPlot_low_bkg_ks->GetY();
    double* YieldPlot_low_bkg_err_ks = subYieldPlot_low_bkg_ks->GetEY();

    double* Nassoc_obs_la            = Nass_obs_la->GetY();
    double* Nassoc_low_obs_la        = Nass_low_obs_la->GetY();
    double* YieldPlot_obs_la         = subYieldPlot_obs_la->GetY();
    double* YieldPlot_obs_err_la     = subYieldPlot_obs_la->GetEY();
    double* YieldPlot_low_obs_la     = subYieldPlot_low_obs_la->GetY();
    double* YieldPlot_low_obs_err_la = subYieldPlot_low_obs_la->GetEY();

    double* Nassoc_bkg_la            = Nass_bkg_la->GetY();
    double* Nassoc_low_bkg_la        = Nass_low_bkg_la->GetY();
    double* YieldPlot_low_bkg_la     = subYieldPlot_low_bkg_la->GetY();
    double* YieldPlot_low_bkg_err_la = subYieldPlot_low_bkg_la->GetEY();
    double* YieldPlot_bkg_la         = subYieldPlot_bkg_la->GetY();
    double* YieldPlot_bkg_err_la     = subYieldPlot_bkg_la->GetEY();

    double* Nassoc_obs_xi            = Nass_obs_xi->GetY();
    double* Nassoc_low_obs_xi        = Nass_low_obs_xi->GetY();
    double* YieldPlot_obs_xi         = subYieldPlot_obs_xi->GetY();
    double* YieldPlot_obs_err_xi     = subYieldPlot_obs_xi->GetEY();
    double* YieldPlot_low_obs_xi     = subYieldPlot_low_obs_xi->GetY();
    double* YieldPlot_low_obs_err_xi = subYieldPlot_low_obs_xi->GetEY();

    double* Nassoc_bkg_xi            = Nass_bkg_xi->GetY();
    double* Nassoc_low_bkg_xi        = Nass_low_bkg_xi->GetY();
    double* YieldPlot_low_bkg_xi     = subYieldPlot_low_bkg_xi->GetY();
    double* YieldPlot_low_bkg_err_xi = subYieldPlot_low_bkg_xi->GetEY();
    double* YieldPlot_bkg_xi         = subYieldPlot_bkg_xi->GetY();
    double* YieldPlot_bkg_err_xi     = subYieldPlot_bkg_xi->GetEY();

    double* ksv2true_Y = ksv2true->GetY();
    double* ksv2true_err = ksv2true->GetEY();
    double* lav2true_Y = lav2true->GetY();
    double* lav2true_err = lav2true->GetEY();
    double* xiv2true_Y = xiv2true->GetY();
    double* xiv2true_err = xiv2true->GetEY();

    double* ksv2_MB_Y = ksv2_MB->GetY();
    double* ksv2_MB_err = ksv2_MB->GetEY();
    double* lav2_MB_Y = lav2_MB->GetY();
    double* lav2_MB_err = lav2_MB->GetEY();
    double* xiv2_MB_Y = xiv2_MB->GetY();
    double* xiv2_MB_err = xiv2_MB->GetEY();

    double* pt_ks = ksv2true->GetX();
    double* pt_la = lav2true->GetX();
    double* pt_xi = xiv2true->GetX();

    double* KET_ks = ksv2true_KET->GetX();
    double* KET_la = lav2true_KET->GetX();
    double* KET_xi = xiv2true_KET->GetX();

    std::vector<double> bkgfrac_ks;
    std::vector<double> bkgfrac_MB_ks;
    std::vector<double> DirectSubYield_ks;
    std::vector<double> DirectSubYield_err_ks;
    std::vector<double> DirectSubYield_low_ks;
    std::vector<double> DirectSubYield_low_err_ks;
    std::vector<double> DirectSubNass_ks;
    std::vector<double> DirectSubNass_err_ks;
    std::vector<double> DirectSubNass_low_ks;
    std::vector<double> DirectSubNass_low_err_ks;

    std::vector<double> bkgfrac_la;
    std::vector<double> bkgfrac_MB_la;
    std::vector<double> DirectSubYield_la;
    std::vector<double> DirectSubYield_err_la;
    std::vector<double> DirectSubYield_low_la;
    std::vector<double> DirectSubYield_low_err_la;
    std::vector<double> DirectSubNass_la;
    std::vector<double> DirectSubNass_err_la;
    std::vector<double> DirectSubNass_low_la;
    std::vector<double> DirectSubNass_low_err_la;

    std::vector<double> bkgfrac_xi;
    std::vector<double> bkgfrac_MB_xi;
    std::vector<double> DirectSubYield_xi;
    std::vector<double> DirectSubYield_err_xi;
    std::vector<double> DirectSubYield_low_xi;
    std::vector<double> DirectSubYield_low_err_xi;
    std::vector<double> DirectSubNass_xi;
    std::vector<double> DirectSubNass_err_xi;
    std::vector<double> DirectSubNass_low_xi;
    std::vector<double> DirectSubNass_low_err_xi;

    std::vector<double> v2true_DirectSub_ks;
    std::vector<double> v2true_DirectSub_err_ks;

    std::vector<double> v2true_DirectSub_la;
    std::vector<double> v2true_DirectSub_err_la;

    std::vector<double> v2true_DirectSub_xi;
    std::vector<double> v2true_DirectSub_err_xi;

    for(unsigned i=0; i<PD.PtBin_ks.size(); i++)
    {
        bkgfrac_ks.push_back(1 - PD.fsig_ks[i]);
        bkgfrac_MB_ks.push_back(1 - PD.fsig_ks_MB[i]);
        DirectSubYield_ks.push_back((YieldPlot_obs_ks[i] - bkgfrac_MB_ks[i]*YieldPlot_bkg_ks[i])/PD.fsig_ks_MB[i]);
        DirectSubYield_err_ks.push_back(sqrt(TMath::Power((YieldPlot_bkg_err_ks[i]*bkgfrac_MB_ks[i]),2) + TMath::Power(YieldPlot_obs_err_ks[i],2))/PD.fsig_ks_MB[i]);
        DirectSubYield_low_ks.push_back((YieldPlot_low_obs_ks[i] - bkgfrac_MB_ks[i]*YieldPlot_low_bkg_ks[i])/PD.fsig_ks_MB[i]);
        DirectSubYield_low_err_ks.push_back(sqrt(TMath::Power((YieldPlot_low_bkg_err_ks[i]*bkgfrac_MB_ks[i]),2) + TMath::Power(YieldPlot_low_obs_err_ks[i],2))/PD.fsig_ks_MB[i]);
        DirectSubNass_ks.push_back((Nassoc_obs_ks[i] - bkgfrac_MB_ks[i]*Nassoc_bkg_ks[i])/PD.fsig_ks_MB[i]);
        DirectSubNass_low_ks.push_back((Nassoc_low_obs_ks[i] - bkgfrac_MB_ks[i]*Nassoc_low_bkg_ks[i])/PD.fsig_ks_MB[i]);

        v2true_DirectSub_ks.push_back(ksv2true_Y[i] - ksv2_MB_Y[i]*DirectSubNass_low_ks[i]/DirectSubNass_ks[i]*DirectSubYield_ks[i]/DirectSubYield_low_ks[i]);
        v2true_DirectSub_err_ks.push_back(sqrt(TMath::Power(ksv2true_err[i],2) + (TMath::Power(ksv2_MB_err[i]/ksv2_MB_Y[i],2) + TMath::Power(DirectSubYield_err_ks[i]/DirectSubYield_ks[i],2) + TMath::Power(DirectSubYield_low_err_ks[i]/DirectSubYield_low_ks[i],2))*TMath::Power(ksv2_MB_Y[i]*DirectSubNass_low_ks[i]/DirectSubNass_ks[i]*DirectSubYield_ks[i]/DirectSubYield_low_ks[i],2)));
    }

    for(unsigned i=0; i<PD.PtBin_la.size(); i++)
    {
        bkgfrac_la.push_back(1 - PD.fsig_la[i]);
        bkgfrac_MB_la.push_back(1 - PD.fsig_la_MB[i]);
        DirectSubYield_la.push_back((YieldPlot_obs_la[i] - bkgfrac_MB_la[i]*YieldPlot_bkg_la[i])/PD.fsig_la_MB[i]);
        DirectSubYield_err_la.push_back(sqrt(TMath::Power((YieldPlot_bkg_err_la[i]*bkgfrac_MB_la[i]),2) + TMath::Power(YieldPlot_obs_err_la[i],2))/PD.fsig_la_MB[i]);
        DirectSubYield_low_la.push_back((YieldPlot_low_obs_la[i] - bkgfrac_MB_la[i]*YieldPlot_low_bkg_la[i])/PD.fsig_la_MB[i]);
        DirectSubYield_low_err_la.push_back(sqrt(TMath::Power((YieldPlot_low_bkg_err_la[i]*bkgfrac_MB_la[i]),2) + TMath::Power(YieldPlot_low_obs_err_la[i],2))/PD.fsig_la_MB[i]);
        DirectSubNass_la.push_back((Nassoc_obs_la[i] - bkgfrac_MB_la[i]*Nassoc_bkg_la[i])/PD.fsig_la_MB[i]);
        DirectSubNass_low_la.push_back((Nassoc_low_obs_la[i] - bkgfrac_MB_la[i]*Nassoc_low_bkg_la[i])/PD.fsig_la_MB[i]);

        v2true_DirectSub_la.push_back(lav2true_Y[i] - lav2_MB_Y[i]*DirectSubNass_low_la[i]/DirectSubNass_la[i]*DirectSubYield_la[i]/DirectSubYield_low_la[i]);
        v2true_DirectSub_err_la.push_back(sqrt(TMath::Power(lav2true_err[i],2) + (TMath::Power(lav2_MB_err[i]/lav2_MB_Y[i],2) + TMath::Power(DirectSubYield_err_la[i]/DirectSubYield_la[i],2) + TMath::Power(DirectSubYield_low_err_la[i]/DirectSubYield_low_la[i],2))*TMath::Power(lav2_MB_Y[i]*DirectSubNass_low_la[i]/DirectSubNass_la[i]*DirectSubYield_la[i]/DirectSubYield_low_la[i],2)));
    }

    for(unsigned i=0; i<PD.PtBin_xi.size(); i++)
    {
        bkgfrac_xi.push_back(1 - PD.fsig_xi[i]);
        DirectSubYield_xi.push_back((YieldPlot_obs_xi[i] - bkgfrac_xi[i]*YieldPlot_bkg_xi[i])/PD.fsig_xi[i]);
        DirectSubYield_err_xi.push_back(sqrt(TMath::Power((YieldPlot_bkg_err_xi[i]*bkgfrac_xi[i]),2) + TMath::Power(YieldPlot_obs_err_xi[i],2))/PD.fsig_xi[i]);
        DirectSubYield_low_xi.push_back((YieldPlot_low_obs_xi[i] - bkgfrac_xi[i]*YieldPlot_low_bkg_xi[i])/PD.fsig_xi[i]);
        DirectSubYield_low_err_xi.push_back(sqrt(TMath::Power((YieldPlot_low_bkg_err_xi[i]*bkgfrac_xi[i]),2) + TMath::Power(YieldPlot_low_obs_err_xi[i],2))/PD.fsig_xi[i]);
        DirectSubNass_xi.push_back((Nassoc_obs_xi[i] - bkgfrac_xi[i]*Nassoc_bkg_xi[i])/PD.fsig_xi[i]);
        DirectSubNass_low_xi.push_back((Nassoc_low_obs_xi[i] - bkgfrac_xi[i]*Nassoc_low_bkg_xi[i])/PD.fsig_xi[i]);

        v2true_DirectSub_xi.push_back(xiv2true_Y[i] - xiv2_MB_Y[i]*DirectSubNass_low_xi[i]/DirectSubNass_xi[i]*DirectSubYield_xi[i]/DirectSubYield_low_xi[i]);
        v2true_DirectSub_err_xi.push_back(sqrt(TMath::Power(xiv2true_err[i],2) + (TMath::Power(xiv2_MB_err[i]/xiv2_MB_Y[i],2) + TMath::Power(DirectSubYield_err_xi[i]/DirectSubYield_xi[i],2) + TMath::Power(DirectSubYield_low_err_xi[i]/DirectSubYield_low_xi[i],2))*TMath::Power(xiv2_MB_Y[i]*DirectSubNass_low_xi[i]/DirectSubNass_xi[i]*DirectSubYield_xi[i]/DirectSubYield_low_xi[i],2)));
    }

    int numPtBins_ks = PD.PtBin_ks.size()-1;
    int numPtBins_la = PD.PtBin_la.size()-1;
    int numPtBins_xi = PD.PtBin_xi.size()-1;

    TGraphErrors* ksv2truesub = new TGraphErrors(numPtBins_ks,pt_ks,&v2true_DirectSub_ks[0],0,&v2true_DirectSub_err_ks[0]);
    TGraphErrors* lav2truesub = new TGraphErrors(numPtBins_la,pt_la,&v2true_DirectSub_la[0],0,&v2true_DirectSub_err_la[0]);
    TGraphErrors* xiv2truesub = new TGraphErrors(numPtBins_xi,pt_xi,&v2true_DirectSub_xi[0],0,&v2true_DirectSub_err_xi[0]);

    TGraphErrors* ksv2truesub_KET = new TGraphErrors(numPtBins_ks,KET_ks,&v2true_DirectSub_ks[0],0,&v2true_DirectSub_err_ks[0]);
    TGraphErrors* lav2truesub_KET = new TGraphErrors(numPtBins_la,KET_la,&v2true_DirectSub_la[0],0,&v2true_DirectSub_err_la[0]);
    TGraphErrors* xiv2truesub_KET = new TGraphErrors(numPtBins_xi,KET_xi,&v2true_DirectSub_xi[0],0,&v2true_DirectSub_err_xi[0]);

    TGraphErrors* TGDirectSubYield_ks = new TGraphErrors(numPtBins_ks,pt_ks,&DirectSubYield_ks[0],0,&DirectSubYield_err_ks[0]);
    TGraphErrors* TGDirectSubYield_low_ks = new TGraphErrors(numPtBins_ks,pt_ks,&DirectSubYield_low_ks[0],0,&DirectSubYield_low_err_ks[0]);
    TGraphErrors* TGDirectSubYield_la = new TGraphErrors(numPtBins_la,pt_la,&DirectSubYield_la[0],0,&DirectSubYield_err_la[0]);
    TGraphErrors* TGDirectSubYield_low_la = new TGraphErrors(numPtBins_la,pt_la,&DirectSubYield_low_la[0],0,&DirectSubYield_low_err_la[0]);
    TGraphErrors* TGDirectSubYield_xi = new TGraphErrors(numPtBins_xi,pt_xi,&DirectSubYield_xi[0],0,&DirectSubYield_err_xi[0]);
    TGraphErrors* TGDirectSubYield_low_xi = new TGraphErrors(numPtBins_xi,pt_xi,&DirectSubYield_low_xi[0],0,&DirectSubYield_low_err_xi[0]);

    TFile output(PD.fn.c_str(),"UPDATE");

    ksv2truesub->Write("kshortv2trueDirectSub",TObject::kOverwrite);
    lav2truesub->Write("lambdav2trueDirectSub",TObject::kOverwrite);
    xiv2truesub->Write("xiv2trueDirectSub",TObject::kOverwrite);

    ksv2truesub_KET->Write("kshortv2trueDirectSub_KET",TObject::kOverwrite);
    lav2truesub_KET->Write("lambdav2trueDirectSub_KET",TObject::kOverwrite);
    xiv2truesub_KET->Write("xiv2trueDirectSub_KET",TObject::kOverwrite);

    TGDirectSubYield_ks->Write("DirectSubYield_ks",TObject::kOverwrite);
    TGDirectSubYield_low_ks->Write("DirectSubYield_low_ks",TObject::kOverwrite);
    TGDirectSubYield_la->Write("DirectSubYield_la",TObject::kOverwrite);
    TGDirectSubYield_low_la->Write("DirectSubYield_low_la",TObject::kOverwrite);
    TGDirectSubYield_xi->Write("DirectSubYield_la",TObject::kOverwrite);
    TGDirectSubYield_low_xi->Write("DirectSubYield_low_la",TObject::kOverwrite);
}

int main()
{
    PeriSubObs(V0);
    //PeriSubBkg(V0);
    PeriSubSigIndirect(V0);
    PeriSubSigDirect(V0);
    return 1;
}
