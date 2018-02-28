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

void perisubrefalltestV0bkg()
{
    TH1::SetDefaultSumw2();

    TFile *_file0 = TFile::Open("../MB10_20/ppCorr_pt033.root");
    TFile *_file1 = TFile::Open("ppCorr_pt033.root");

    TFile *_file2 = TFile::Open("../MB10_20/ppCorr_all_V0_pt033_7TeVEff.root");
    TFile *_file3 = TFile::Open("ppCorr_all_V0_pt033_7TeVEff.root");
    
    TLatex* tex = new TLatex();
    tex->SetNDC();
    
    char Nplot[13][200] = {"0.2","0.4","0.6","0.8","1.0","1.4","1.8","2.2","2.8","3.6","4.6","6.0","9.0"};
 
    //=============create canvas for fitting==============
    
    TCanvas* c1 = new TCanvas("lowNsr","lowNsr",1200,900);
    TCanvas* c2 = new TCanvas("lowNlr","lowNlr",1200,900);
    TCanvas* c3 = new TCanvas("highNsr","highNsr",1200,900);
    TCanvas* c4 = new TCanvas("highNlr","highNlr",1200,900);
    TCanvas* c5 = new TCanvas("lowFit","lowFit",1200,900);
    TCanvas* c6 = new TCanvas("highFit","highFit",1200,900);

    c1->Divide(4,3);
    c2->Divide(4,3);
    c3->Divide(4,3);
    c4->Divide(4,3);
    c5->Divide(4,3);
    c6->Divide(4,3);

    TCanvas* c = new TCanvas("c","c",400,400);
    
    c->cd();
   
    double BW2D = (9.9/33)*((2-1.0/16)*3.1416/31);
    double errfactor = sqrt(1.0);
    
    //==============High N ref Vn and jet yield=========================
    TH2D* signal_ref = _file1.Get("pp_HM105_GplusPP/signal");
    TH2D* background_ref = _file1.Get("pp_HM105_GplusPP/background");
    TH1D* mult_ref = _file1.Get("pp_HM105_GplusPP/mult");
    
    TF1* fit = new TF1("fit","[0]*x^2+[1]*x+[2]",0.6,2.2);
    fit->SetParameters(1,1,1);
    
    double nEvent_ref = mult_ref->Integral(3,10000);
    double Bz_ref = background_ref->GetBinContent(background_ref->FindBin(0,0));
    
    TH1D* srsks = signal_ref->ProjectionY("srsks",14,20);
    TH1D* srbks = background_ref->ProjectionY("srbks",14,20);
    srsks->Divide(srbks);
    srsks->Scale(Bz_ref/nEvent_ref/BW2D);
    
    srsks->Fit("fit","R");
    srsks->Fit("fit","R");
    srsks->Fit("fit","R");
    srsks->Fit("fit","R");
    srsks->Fit("fit","R");
    
    double srmks = fit->GetMinimum(0.6,2.2);
    double srmksx = fit->GetMinimumX(0.6,2.2);
    TF1* mfsrks = new TF1("mfsrks","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    mfsrks->SetParameter(0,-srmks);
    srsks->Add(mfsrks);
    double srye_ref;
    double sry_ref = srsks->IntegralAndError(srsks->FindBin(0.0),srsks->FindBin(srmksx),srye_ref,"width");
    double bin0yield = srsks->GetBinContent(srsks->FindBin(0.0))*0.19635;
    sry_ref = sry_ref*2 - bin0yield;
    
TF1* fit10 = new TF1("fit10","[0]*x^2+[1]*x+[2]",0,2.0);
    fit10->SetParameters(1,1,1);

    TH1D* alrs = signal_ref->ProjectionY("alrsks22",1,10);
    TH1D* alrs1 = signal_ref->ProjectionY("alrs1",24,33);
    TH1D* alrb = background_ref->ProjectionY("alrb",1,10);
    TH1D* alrb1 = background_ref->ProjectionY("alrb1",24,33);
    alrs->Add(alrs1);
    alrb->Add(alrb1);
    alrs->Divide(alrb);
    alrs->Scale(Bz_ref/nEvent_ref/BW2D);
    
    alrs->Fit("fit10","R");
    alrs->Fit("fit10","R");
    alrs->Fit("fit10","R");
    alrs->Fit("fit10","R");
    alrs->Fit("fit10","R");
    
    double lrm = fit10->GetMinimum(0,2.0);
    double lrmksx = fit10->GetMinimumX(0,2.0);
    TF1* mflr = new TF1("mflr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    mflr->SetParameter(0,-lrm);
    alrs->Add(mflr);
    double lrye_ref;
    double lry_ref = alrs->IntegralAndError(alrs->FindBin(0.0),alrs->FindBin(lrmksx),lrye_ref,"width");
    double bin0yield = alrs->GetBinContent(alrs->FindBin(0.0))*0.19635;
    lry_ref = lry_ref*2 - bin0yield;
    
    double suby_ref = sry_ref - lry_ref;
    double subye_ref = sqrt(srye_ref*srye_ref+lrye_ref*lrye_ref);
    
    /*TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit1->SetParNames("N","V1","V2","V3","V4");
    fit1->SetParameters(10,1,1,1,1);*/
    
    TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit1->SetParNames("N","V1","V2","V3");
    fit1->SetParameters(10,1,1,1);
    
    TH1D* lrs = signal_ref->ProjectionY("lrs",1,10);
    TH1D* lrs1 = signal_ref->ProjectionY("lrs1",24,33);
    TH1D* lrb = background_ref->ProjectionY("lrb",1,10);
    TH1D* lrb1 = background_ref->ProjectionY("lrb1",24,33);
    lrs->Add(lrs1);
    lrb->Add(lrb1);
    lrs->Divide(lrb);
lrs->Scale(Bz_ref/nEvent_ref/BW2D);
    lrs->Fit("fit1","R");
    lrs->Fit("fit1","R");
    lrs->Fit("fit1","R");
    lrs->Fit("fit1","R");
    
    double V2_ref = fit1->GetParameter(2);
    double V2e_ref = fit1->GetParError(2)*errfactor;

    double V3_ref = fit1->GetParameter(3);
    double V3e_ref = fit1->GetParError(3)*errfactor;

    double V1_ref = fit1->GetParameter(1);
    double V1e_ref = fit1->GetParError(1)*errfactor;
    
    double v2_ref = sqrt(V2_ref);
    double v2e_ref = sqrt(V2_ref)*(V2e_ref/V2_ref)/2;
    double v3_ref = sqrt(V3_ref);
    double v3e_ref = sqrt(V3_ref)*(V3e_ref/V3_ref)/2;
    double v1_ref = sqrt(V1_ref);
    double v1e_ref = sqrt(V1_ref)*(V1e_ref/V1_ref)/2;
    
    double Nassoc_ref_fit = fit1->GetParameter(0);
    
    //=======================Low N ref Vn and jet yield===================
    TH2D* signal_ref_low = _file0.Get("pp_MB10_GplusPP/signal");
    TH2D* background_ref_low = _file0.Get("pp_MB10_GplusPP/background");
    TH1D* mult_ref_low = _file0.Get("pp_MB10_GplusPP/mult");
    
    TF1* fit = new TF1("fit","[0]*x^2+[1]*x+[2]",0.6,2.2);
    fit->SetParameters(1,1,1);
    
    double nEvent_ref_low = mult_ref_low->Integral(3,10000);
    double Bz_ref_low = background_ref_low->GetBinContent(background_ref_low->FindBin(0,0));
    
    TH1D* srsks_low = signal_ref_low->ProjectionY("srsks_low",14,20);
    TH1D* srbks_low = background_ref_low->ProjectionY("srbks_low",14,20);
    srsks_low->Divide(srbks_low);
    srsks_low->Scale(Bz_ref_low/nEvent_ref_low/BW2D);
    
    srsks_low->Fit("fit","R");
    srsks_low->Fit("fit","R");
    srsks_low->Fit("fit","R");
    srsks_low->Fit("fit","R");
    srsks_low->Fit("fit","R");
    
    TCanvas* c22 = new TCanvas;
    c22->cd();
    srsks_low->Draw();
    c->cd();
    
    double srmks_low = fit->GetMinimum(0.6,2.2);
    double srmksx_low = fit->GetMinimumX(0.6,2.2);
    TF1* mfsrks_low = new TF1("mfsrks_low","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    mfsrks_low->SetParameter(0,-srmks_low);
    srsks_low->Add(mfsrks_low);
    double srye_ref_low;
    double sry_ref_low = srsks_low->IntegralAndError(srsks_low->FindBin(0.0),srsks_low->FindBin(srmksx_low),srye_ref_low,"width");
    double bin0yield = srsks_low->GetBinContent(srsks_low->FindBin(0.0))*0.19635;
    sry_ref_low = sry_ref_low*2 - bin0yield;
    
    TF1* fit10 = new TF1("fit10","[0]*x^2+[1]*x+[2]",0,2.0);
    fit10->SetParameters(1,1,1);
    
    TH1D* alrsks = signal_ref_low->ProjectionY("alrsks",1,10);
    TH1D* alrs1 = signal_ref_low->ProjectionY("alrs1",24,33);
    TH1D* alrb = background_ref_low->ProjectionY("alrb",1,10);
    TH1D* alrb1 = background_ref_low->ProjectionY("alrb1",24,33);
    alrsks->Add(alrs1);
    alrb->Add(alrb1);
    alrsks->Divide(alrb);
    alrsks->Scale(Bz_ref_low/nEvent_ref_low/BW2D);
    
    alrsks->Fit("fit10","R");
    alrsks->Fit("fit10","R");
    alrsks->Fit("fit10","R");
    alrsks->Fit("fit10","R");
    alrsks->Fit("fit10","R");
    
    double lrm = fit10->GetMinimum(0,2.0);
	double lrmksx_low = fit10->GetMinimumX(0,2.0);
    TF1* mflr = new TF1("mflr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    mflr->SetParameter(0,-lrm);
    alrsks->Add(mflr);
    double lrye_ref_low;
    double lry_ref_low = alrsks->IntegralAndError(alrsks->FindBin(0.0),alrsks->FindBin(lrmksx_low),lrye_ref_low,"width");
    double bin0yield = alrsks->GetBinContent(alrsks->FindBin(0.0))*0.19635;
    lry_ref_low = lry_ref_low*2 -bin0yield;
    
    double suby_ref_low = sry_ref_low - lry_ref_low;
    //double suby_ref_low = sry_ref_low;
    double subye_ref_low = sqrt(srye_ref_low*srye_ref_low+lrye_ref_low*lrye_ref_low);
    
    /*TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit1->SetParNames("N","V1","V2","V3","V4");
    fit1->SetParameters(10,1,1,1,1);*/
    
    TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit1->SetParNames("N","V1","V2","V3");
    fit1->SetParameters(10,1,1,1);
    
    TH1D* lrs = signal_ref_low->ProjectionY("lrs",1,10);
    TH1D* lrs1 = signal_ref_low->ProjectionY("lrs1",24,33);
    TH1D* lrb = background_ref_low->ProjectionY("lrb",1,10);
    TH1D* lrb1 = background_ref_low->ProjectionY("lrb1",24,33);
    lrs->Add(lrs1);
    lrb->Add(lrb1);
    lrs->Divide(lrb);
    lrs->Scale(Bz_ref_low/nEvent_ref_low/BW2D);
    lrs->Fit("fit1","R");
    lrs->Fit("fit1","R");
    lrs->Fit("fit1","R");
    lrs->Fit("fit1","R");
    
    double V1_ref_low = fit1->GetParameter(1);
    double V1e_ref_low = fit1->GetParError(1)*errfactor;
    
    double v1_ref_low = sqrt(V1_ref_low);
    double v1e_ref_low = sqrt(V1_ref_low)*(V1e_ref_low/V1_ref_low)/2;
    
    double V2_ref_low = fit1->GetParameter(2);
    double V2e_ref_low = fit1->GetParError(2)*errfactor;
    
    double v2_ref_low = sqrt(V2_ref_low);
    double v2e_ref_low = sqrt(V2_ref_low)*(V2e_ref_low/V2_ref_low)/2;

    double V3_ref_low = fit1->GetParameter(3);
    double V3e_ref_low = fit1->GetParError(3)*errfactor;
    
    double v3_ref_low = sqrt(V3_ref_low);
    double v3e_ref_low = sqrt(V3_ref_low)*(V3e_ref_low/V3_ref_low)/2;

    double Nassoc_ref_fit_low = fit1->GetParameter(0);
    double Nassoc_ref_low;
    TH1D* mult_assoc_low_ref;
    _file0->GetObject("pp_MB10_GplusPP/mult_assoc",mult_assoc_low_ref);
    //mult_assoc_low_ref->GetXaxis()->SetRangeUser(2,300);
    Nassoc_ref_low = mult_assoc_low_ref->GetMean(1);
    
    //====================Low N pT differential Vn and jet yield for obs================
    //===============KS==================
    int nEvent_low_ks[13];
    double Bz_low_ks[13];
    double Nassoc_low_ks[13];
    double Nassoc_fit_low_ks[13];
    
    double srye_low_ks[13];
    double sry_low_ks[13];
    double lrye_low_ks[13];
    double lry_low_ks[13];
    double subye_low_ks[13];
    double suby_low_ks[13];
    
    double V2_low_ks[13];
    double V2e_low_ks[13];
    
    double V3_low_ks[13];
    double V3e_low_ks[13];

    double V1_low_ks[13];
    double V1e_low_ks[13];
    
    TH2D* signal_low_ks[13];
    TH2D* background_low_ks[13];
    
    TH1D* mult_low_ks[13];
    TH1D* mult_assoc_low_ks[13];
    
    TH1D* srs_low_ks[13];
    TH1D* alrs_low_ks[13];

    TH1D* lrb_low_ks[13];
    TH1D* lrs_low_ks[13];
    
    for(int i=0;i<9;i++){
        //sr Y
        _file2->GetObject(Form("pp_MB10_GplusPP/signalkshort_bkg_pt%d",i),signal_low_ks[i]);
        _file2->GetObject(Form("pp_MB10_GplusPP/backgroundkshort_bkg_pt%d",i),background_low_ks[i]);
        _file2->GetObject(Form("pp_MB10_GplusPP/mult_ks_bkg_pt%d",i),mult_low_ks[i]);
        _file2->GetObject(Form("pp_MB10_GplusPP/mult_ass",i),mult_assoc_low_ks[i]);
        Nassoc_low_ks[i] = mult_assoc_low_ks[i]->GetMean(1);

        TF1* fit2 = new TF1("fit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
        fit2->SetParameters(1,1,1);
        fit2->SetLineColor(2);
        TF1* fit = new TF1("fit","[0]*x^2+[1]*x+[2]",0.0,2.0);
        fit->SetParameters(1,1,1);
        fit->SetLineColor(2);

        nEvent_low_ks[i] = mult_low_ks[i]->Integral(2,10000);
        Bz_low_ks[i] = background_low_ks[i]->GetBinContent(background_low_ks[i]->FindBin(0,0));
        c->cd();
        (TH1D*)srs_low_ks[i] = signal_low_ks[i]->ProjectionY(Form("srslow_ks%d",i),14,20);
        TH1D* srb = background_low_ks[i]->ProjectionY("srb",14,20);
        srs_low_ks[i]->Divide(srb);
        srs_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);
        
        srs_low_ks[i]->Fit("fit2","R");
        srs_low_ks[i]->Fit("fit2","R");
        srs_low_ks[i]->Fit("fit2","R");
        srs_low_ks[i]->Fit("fit2","R");
        srs_low_ks[i]->Fit("fit2","R");
        
        double srm = fit2->GetMinimum(0.6,2.2);
	double srmx = fit2->GetMinimumX(0.6,2.2);
        TF1* mfsr = new TF1("mfsr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        mfsr->SetParameter(0,-srm);
        TH1D* srslow = srs_low_ks[i]->Clone();
        srslow->Add(mfsr);
        c1->cd(i+1);
        srs_low_ks[i]->Draw();
        tex->DrawLatex(0.52,0.88,"CMS pp 7TeV");
        tex->DrawLatex(0.52,0.82,"10<N_{trk}^{offline}<20");
        tex->DrawLatex(0.52,0.74,Form("%s<p_{T}^{trg}<%s",Nplot[i],Nplot[i+1]));
        tex->DrawLatex(0.52,0.67,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(0.52,0.59,"|#Delta#eta|<1");
        c->cd();
        sry_low_ks[i] = srslow->IntegralAndError(srslow->FindBin(0.0),srslow->FindBin(srmx),srye_low_ks[i],"width");
	double bin0yield = srslow->GetBinContent(srslow->FindBin(0.0))*0.19635;
	sry_low_ks[i] = sry_low_ks[i]*2 - bin0yield;
        //cout<<"111"<<endl;
        //lr Y
        (TH1D*)alrs_low_ks[i] = signal_low_ks[i]->ProjectionY(Form("alrslow_ks%d",i),1,10);
        alrs1 = signal_low_ks[i]->ProjectionY("alrs1",24,33);
        alrb = background_low_ks[i]->ProjectionY("alrb",1,10);
        alrb1 = background_low_ks[i]->ProjectionY("alrb1",24,33);
        alrs_low_ks[i]->Add(alrs1);
        alrb->Add(alrb1);
        alrs_low_ks[i]->Divide(alrb);
        alrs_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);
        
        alrs_low_ks[i]->Fit("fit","R");
        alrs_low_ks[i]->Fit("fit","R");
        alrs_low_ks[i]->Fit("fit","R");
        alrs_low_ks[i]->Fit("fit","R");
        alrs_low_ks[i]->Fit("fit","R");
        
        double lrm = fit->GetMinimum(0.0,2.0);
        double lrmx = fit->GetMinimumX(0.0,2.0);
        TF1* mflr = new TF1("mflr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        mflr->SetParameter(0,-lrm);
        TH1D* alrslow = alrs_low_ks[i]->Clone();
        alrslow->Add(mflr);
        if(i==11) proj_low = alrslow->Clone();
        c2->cd(i+1);
        alrs_low_ks[i]->Draw();
        tex->DrawLatex(0.22,0.88,"CMS pp 7TeV");
        tex->DrawLatex(0.22,0.82,"10<N_{trk}^{offline}<20");
        tex->DrawLatex(0.22,0.74,Form("%s<p_{T}^{trg}<%s",Nplot[i],Nplot[i+1]));
        tex->DrawLatex(0.22,0.67,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(0.22,0.59,"|#Delta#eta|>2");
        c->cd();
        lry_low_ks[i] = alrslow->IntegralAndError(alrslow->FindBin(0.0),alrslow->FindBin(lrmx),lrye_low_ks[i],"width");
	double bin0yield = alrslow->GetBinContent(alrslow->FindBin(0.0))*0.19635;
	lry_low_ks[i] = lry_low_ks[i]*2 - bin0yield;
        
        //sub Y
        suby_low_ks[i] = sry_low_ks[i] - lry_low_ks[i];
        //suby_low[i] = sry_low[i];
        subye_low_ks[i] = sqrt(srye_low_ks[i]*srye_low_ks[i]+lrye_low_ks[i]*lrye_low_ks[i]);
        c5->cd(i+1);

        /*TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        fit1->SetParNames("N","V1","V2","V3","V4");
        fit1->SetParameters(10,1,1,1,1);
        fit1->SetLineColor(2);*/
        
        TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        fit1->SetParNames("N","V1","V2","V3");
        fit1->SetParameters(10,1,1,1);
        fit1->SetLineColor(2);
        
        //lr V2 obs
        lrs_low_ks[i] = signal_low_ks[i]->ProjectionY(Form("lrs_ks%d",i),1,10);
        TH1D* lrsa = signal_low_ks[i]->ProjectionY("lrsa",24,33);
        lrb_low_ks[i] = background_low_ks[i]->ProjectionY(Form("lrb_ks%d",i),1,10);
        TH1D* lrba = background_low_ks[i]->ProjectionY("lrba",24,33);
        lrs_low_ks[i]->Add(lrsa);
        lrb_low_ks[i]->Add(lrba);
        lrs_low_ks[i]->Divide(lrb_low_ks[i]);
        lrs_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);
        //lrs_low[i]->Scale(1.0/4.2);
        lrs_low_ks[i]->Fit("fit1","R");
        lrs_low_ks[i]->Fit("fit1","R");
        lrs_low_ks[i]->Fit("fit1","R");
        lrs_low_ks[i]->Fit("fit1","R");
        
        V2_low_ks[i] = fit1->GetParameter(2);
        V2e_low_ks[i] = fit1->GetParError(2)*errfactor;
        V3_low_ks[i] = fit1->GetParameter(3);
        V3e_low_ks[i] = fit1->GetParError(3)*errfactor;
        V1_low_ks[i] = fit1->GetParameter(1);
        V1e_low_ks[i] = fit1->GetParError(1)*errfactor;

        Nassoc_fit_low_ks[i] = fit1->GetParameter(0);
        c->cd();
    }
    
    //===============LA==================
    int nEvent_low_la[13];
    double Bz_low_la[13];
    double Nassoc_low_la[13];
    double Nassoc_fit_low_la[13];
    
    double srye_low_la[13];
    double sry_low_la[13];
    double lrye_low_la[13];
    double lry_low_la[13];
    double subye_low_la[13];
    double suby_low_la[13];
    
    double V2_low_la[13];
    double V2e_low_la[13];
    
    double V3_low_la[13];
    double V3e_low_la[13];
    
    double V1_low_la[13];
    double V1e_low_la[13];
    
    TH2D* signal_low_la[13];
    TH2D* background_low_la[13];
    
    TH1D* mult_low_la[13];
    TH1D* mult_assoc_low_la[13];
    
    TH1D* srs_low_la[13];
    TH1D* alrs_low_la[13];
    
    TH1D* lrb_low_la[13];
    TH1D* lrs_low_la[13];
    
    for(int i=0;i<9;i++){
        //sr Y
        _file2->GetObject(Form("pp_MB10_GplusPP/signallambda_bkg_pt%d",i),signal_low_la[i]);
        _file2->GetObject(Form("pp_MB10_GplusPP/backgroundlambda_bkg_pt%d",i),background_low_la[i]);
        _file2->GetObject(Form("pp_MB10_GplusPP/mult_la_bkg_pt%d",i),mult_low_la[i]);
        _file2->GetObject(Form("pp_MB10_GplusPP/mult_ass",i),mult_assoc_low_la[i]);
        Nassoc_low_la[i] = mult_assoc_low_la[i]->GetMean(1);
        
        TF1* fit2 = new TF1("fit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
        fit2->SetParameters(1,1,1);
        fit2->SetLineColor(2);
        TF1* fit = new TF1("fit","[0]*x^2+[1]*x+[2]",0.0,2.0);
        fit->SetParameters(1,1,1);
        fit->SetLineColor(2);
        
        nEvent_low_la[i] = mult_low_la[i]->Integral(2,10000);
        Bz_low_la[i] = background_low_la[i]->GetBinContent(background_low_la[i]->FindBin(0,0));
        c->cd();
        (TH1D*)srs_low_la[i] = signal_low_la[i]->ProjectionY(Form("srslow_la%d",i),14,20);
        TH1D* srb = background_low_la[i]->ProjectionY("srb",14,20);
        srs_low_la[i]->Divide(srb);
        srs_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);
        
        srs_low_la[i]->Fit("fit2","R");
        srs_low_la[i]->Fit("fit2","R");
        srs_low_la[i]->Fit("fit2","R");
        srs_low_la[i]->Fit("fit2","R");
        srs_low_la[i]->Fit("fit2","R");
        
        double srm = fit2->GetMinimum(0.6,2.2);
        double srmx = fit2->GetMinimumX(0.6,2.2);
        TF1* mfsr = new TF1("mfsr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        mfsr->SetParameter(0,-srm);
        TH1D* srslow = (TH1D*)srs_low_la[i]->Clone();
        srslow->Add(mfsr);
        c1->cd(i+1);
        srs_low_la[i]->Draw();
        tex->DrawLatex(0.52,0.88,"CMS pp 7TeV");
        tex->DrawLatex(0.52,0.82,"10<N_{trk}^{offline}<20");
        tex->DrawLatex(0.52,0.74,Form("%s<p_{T}^{trg}<%s",Nplot[i],Nplot[i+1]));
        tex->DrawLatex(0.52,0.67,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(0.52,0.59,"|#Delta#eta|<1");
        c->cd();
        sry_low_la[i] = srslow->IntegralAndError(srslow->FindBin(0.0),srslow->FindBin(srmx),srye_low_la[i],"width");
        double bin0yield = srslow->GetBinContent(srslow->FindBin(0.0))*0.19635;
        sry_low_la[i] = sry_low_la[i]*2 - bin0yield;
        //cout<<"111"<<endl;
        //lr Y
        (TH1D*)alrs_low_la[i] = signal_low_la[i]->ProjectionY(Form("alrslow_la%d",i),1,10);
        alrs1 = signal_low_la[i]->ProjectionY("alrs1",24,33);
        alrb = background_low_la[i]->ProjectionY("alrb",1,10);
        alrb1 = background_low_la[i]->ProjectionY("alrb1",24,33);
        alrs_low_la[i]->Add(alrs1);
        alrb->Add(alrb1);
        alrs_low_la[i]->Divide(alrb);
        alrs_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);
        
        alrs_low_la[i]->Fit("fit","R");
        alrs_low_la[i]->Fit("fit","R");
        alrs_low_la[i]->Fit("fit","R");
        alrs_low_la[i]->Fit("fit","R");
        alrs_low_la[i]->Fit("fit","R");
        
        double lrm = fit->GetMinimum(0.0,2.0);
        double lrmx = fit->GetMinimumX(0.0,2.0);
        TF1* mflr = new TF1("mflr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        mflr->SetParameter(0,-lrm);
        TH1D* alrslow = (TH1D*)alrs_low_la[i]->Clone();
        alrslow->Add(mflr);
        c2->cd(i+1);
        alrs_low_la[i]->Draw();
        tex->DrawLatex(0.22,0.88,"CMS pp 7TeV");
        tex->DrawLatex(0.22,0.82,"10<N_{trk}^{offline}<20");
        tex->DrawLatex(0.22,0.74,Form("%s<p_{T}^{trg}<%s",Nplot[i],Nplot[i+1]));
        tex->DrawLatex(0.22,0.67,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(0.22,0.59,"|#Delta#eta|>2");
        c->cd();
        lry_low_la[i] = alrslow->IntegralAndError(alrslow->FindBin(0.0),alrslow->FindBin(lrmx),lrye_low_la[i],"width");
        double bin0yield = alrslow->GetBinContent(alrslow->FindBin(0.0))*0.19635;
        lry_low_la[i] = lry_low_la[i]*2 - bin0yield;
        
        //sub Y
        suby_low_la[i] = sry_low_la[i] - lry_low_la[i];
        //suby_low[i] = sry_low[i];
        subye_low_la[i] = sqrt(srye_low_la[i]*srye_low_la[i]+lrye_low_la[i]*lrye_low_la[i]);
        c5->cd(i+1);
        
        /*TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        fit1->SetParNames("N","V1","V2","V3","V4");
        fit1->SetParameters(10,1,1,1,1);
        fit1->SetLineColor(2);*/
        
        TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        fit1->SetParNames("N","V1","V2","V3");
        fit1->SetParameters(10,1,1,1);
        fit1->SetLineColor(2);
        
        //lr V2 obs
        lrs_low_la[i] = signal_low_la[i]->ProjectionY(Form("lrs_la%d",i),1,10);
        TH1D* lrsa = signal_low_la[i]->ProjectionY("lrsa",24,33);
        lrb_low_la[i] = background_low_la[i]->ProjectionY(Form("lrb_la%d",i),1,10);
        TH1D* lrba = background_low_la[i]->ProjectionY("lrba",24,33);
        lrs_low_la[i]->Add(lrsa);
        lrb_low_la[i]->Add(lrba);
        lrs_low_la[i]->Divide(lrb_low_la[i]);
        lrs_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);
        //lrs_low[i]->Scale(1.0/4.2);
        lrs_low_la[i]->Fit("fit1","R");
        lrs_low_la[i]->Fit("fit1","R");
        lrs_low_la[i]->Fit("fit1","R");
        lrs_low_la[i]->Fit("fit1","R");
        
        V2_low_la[i] = fit1->GetParameter(2);
        V2e_low_la[i] = fit1->GetParError(2)*errfactor;
        V3_low_la[i] = fit1->GetParameter(3);
        V3e_low_la[i] = fit1->GetParError(3)*errfactor;
        V1_low_la[i] = fit1->GetParameter(1);
        V1e_low_la[i] = fit1->GetParError(1)*errfactor;
        
        Nassoc_fit_low_la[i] = fit1->GetParameter(0);
        c->cd();
    }
    
    //==========================High N pT differential Vn and jet yield for obs
    double Nassoc_ref;
    TH1D* mult_assoc_ref;
    _file1->GetObject("pp_HM105_GplusPP/mult_assoc",mult_assoc_ref);
    //mult_assoc_ref->GetXaxis()->SetRangeUser(2,300);
    Nassoc_ref = mult_assoc_ref->GetMean(1);
    
    //===============KS==================
    int nEvent_ks[13];
    double Bz_ks[13];
    double Nassoc_ks[13];
    double Nassoc_fit_ks[13];
    
    double srye_ks[13];
    double sry_ks[13];
    double lrye_ks[13];
    double lry_ks[13];
    double subye_ks[13];
    double suby_ks[13];
    
    double V2_ks[13];
    double V2e_ks[13];
    
    double V3_ks[13];
    double V3e_ks[13];
    
    double V1_ks[13];
    double V1e_ks[13];
    
    TH2D* signal_ks[13];
    TH2D* background_ks[13];
    
    TH1D* mult_ks[13];
    TH1D* mult_assoc_ks[13];
    
    TH1D* srs_ks[13];
    TH1D* alrs_ks[13];
    
    TH1D* lrb_ks[13];
    TH1D* lrs_ks[13];
    
    for(int i=0;i<9;i++){
        //sr Y
        _file3->GetObject(Form("pp_HM105_GplusPP/signalkshort_bkg_pt%d",i),signal_ks[i]);
        _file3->GetObject(Form("pp_HM105_GplusPP/backgroundkshort_bkg_pt%d",i),background_ks[i]);
        _file3->GetObject(Form("pp_HM105_GplusPP/mult_ks_bkg_pt%d",i),mult_ks[i]);
        _file3->GetObject(Form("pp_HM105_GplusPP/mult_ass",i),mult_assoc_ks[i]);
        Nassoc_ks[i] = mult_assoc_ks[i]->GetMean(1);
        
        TF1* fit2 = new TF1("fit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
        fit2->SetParameters(1,1,1);
        fit2->SetLineColor(2);
        TF1* fit = new TF1("fit","[0]*x^2+[1]*x+[2]",0.0,2.0);
        fit->SetParameters(1,1,1);
        fit->SetLineColor(2);
        
        nEvent_ks[i] = mult_ks[i]->Integral(2,10000);
        Bz_ks[i] = background_ks[i]->GetBinContent(background_ks[i]->FindBin(0,0));
        c->cd();
        (TH1D*)srs_ks[i] = signal_ks[i]->ProjectionY(Form("srsks%d",i),14,20);
        TH1D* srb = background_ks[i]->ProjectionY("srb",14,20);
        srs_ks[i]->Divide(srb);
        srs_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);
        
        srs_ks[i]->Fit("fit2","R");
        srs_ks[i]->Fit("fit2","R");
        srs_ks[i]->Fit("fit2","R");
        srs_ks[i]->Fit("fit2","R");
        srs_ks[i]->Fit("fit2","R");
        
        double srm = fit2->GetMinimum(0.6,2.2);
        double srmx = fit2->GetMinimumX(0.6,2.2);
        TF1* mfsr = new TF1("mfsr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        mfsr->SetParameter(0,-srm);
        TH1D* srslow = (TH1D*)srs_ks[i]->Clone();
        srslow->Add(mfsr);
        c1->cd(i+1);
        srs_ks[i]->Draw();
        tex->DrawLatex(0.52,0.88,"CMS pp 7TeV");
        tex->DrawLatex(0.52,0.82,"10<N_{trk}^{offline}<20");
        tex->DrawLatex(0.52,0.74,Form("%s<p_{T}^{trg}<%s",Nplot[i],Nplot[i+1]));
        tex->DrawLatex(0.52,0.67,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(0.52,0.59,"|#Delta#eta|<1");
        c->cd();
        sry_ks[i] = srslow->IntegralAndError(srslow->FindBin(0.0),srslow->FindBin(srmx),srye_ks[i],"width");
        double bin0yield = srslow->GetBinContent(srslow->FindBin(0.0))*0.19635;
        sry_ks[i] = sry_ks[i]*2 - bin0yield;
        //cout<<"111"<<endl;
        //lr Y
        (TH1D*)alrs_ks[i] = signal_ks[i]->ProjectionY(Form("alrsks%d",i),1,10);
        alrs1 = signal_ks[i]->ProjectionY("alrs1",24,33);
        alrb = background_ks[i]->ProjectionY("alrb",1,10);
        alrb1 = background_ks[i]->ProjectionY("alrb1",24,33);
        alrs_ks[i]->Add(alrs1);
        alrb->Add(alrb1);
        alrs_ks[i]->Divide(alrb);
        alrs_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);
        
        alrs_ks[i]->Fit("fit","R");
        alrs_ks[i]->Fit("fit","R");
        alrs_ks[i]->Fit("fit","R");
        alrs_ks[i]->Fit("fit","R");
        alrs_ks[i]->Fit("fit","R");
        
        double lrm = fit->GetMinimum(0.0,2.0);
        double lrmx = fit->GetMinimumX(0.0,2.0);
        TF1* mflr = new TF1("mflr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        mflr->SetParameter(0,-lrm);
        TH1D* alrslow = (TH1D*)alrs_ks[i]->Clone();
        alrslow->Add(mflr);
        c2->cd(i+1);
        alrs_ks[i]->Draw();
        tex->DrawLatex(0.22,0.88,"CMS pp 7TeV");
        tex->DrawLatex(0.22,0.82,"10<N_{trk}^{offline}<20");
        tex->DrawLatex(0.22,0.74,Form("%s<p_{T}^{trg}<%s",Nplot[i],Nplot[i+1]));
        tex->DrawLatex(0.22,0.67,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(0.22,0.59,"|#Delta#eta|>2");
        c->cd();
        lry_ks[i] = alrslow->IntegralAndError(alrslow->FindBin(0.0),alrslow->FindBin(lrmx),lrye_ks[i],"width");
        double bin0yield = alrslow->GetBinContent(alrslow->FindBin(0.0))*0.19635;
        lry_ks[i] = lry_ks[i]*2 - bin0yield;
        
        //sub Y
        suby_ks[i] = sry_ks[i] - lry_ks[i];
        //suby_low[i] = sry_low[i];
        subye_ks[i] = sqrt(srye_ks[i]*srye_ks[i]+lrye_ks[i]*lrye_ks[i]);
        c5->cd(i+1);
        
        /*TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        fit1->SetParNames("N","V1","V2","V3","V4");
        fit1->SetParameters(10,1,1,1,1);
        fit1->SetLineColor(2);*/
        
        TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        fit1->SetParNames("N","V1","V2","V3");
        fit1->SetParameters(10,1,1,1);
        fit1->SetLineColor(2);
        
        //lr V2 obs
        lrs_ks[i] = signal_ks[i]->ProjectionY(Form("lrs_ks%d",i),1,10);
        TH1D* lrsa = signal_ks[i]->ProjectionY("lrsa",24,33);
        lrb_ks[i] = background_ks[i]->ProjectionY(Form("lrb_ks%d",i),1,10);
        TH1D* lrba = background_ks[i]->ProjectionY("lrba",24,33);
        lrs_ks[i]->Add(lrsa);
        lrb_ks[i]->Add(lrba);
        lrs_ks[i]->Divide(lrb_ks[i]);
        lrs_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);
        //lrs_low[i]->Scale(1.0/4.2);
        lrs_ks[i]->Fit("fit1","R");
        lrs_ks[i]->Fit("fit1","R");
        lrs_ks[i]->Fit("fit1","R");
        lrs_ks[i]->Fit("fit1","R");
        
        V2_ks[i] = fit1->GetParameter(2);
        V2e_ks[i] = fit1->GetParError(2)*errfactor;
        V3_ks[i] = fit1->GetParameter(3);
        V3e_ks[i] = fit1->GetParError(3)*errfactor;
        V1_ks[i] = fit1->GetParameter(1);
        V1e_ks[i] = fit1->GetParError(1)*errfactor;
        
        Nassoc_fit_ks[i] = fit1->GetParameter(0);
        c->cd();
    }
    
    //===============LA==================
    int nEvent_la[13];
    double Bz_la[13];
    double Nassoc_la[13];
    double Nassoc_fit_la[13];
    
    double srye_la[13];
    double sry_la[13];
    double lrye_la[13];
    double lry_la[13];
    double subye_la[13];
    double suby_la[13];
    
    double V2_la[13];
    double V2e_la[13];
    
    double V3_la[13];
    double V3e_la[13];
    
    double V1_la[13];
    double V1e_la[13];
    
    TH2D* signal_la[13];
    TH2D* background_la[13];
    
    TH1D* mult_la[13];
    TH1D* mult_assoc_la[13];
    
    TH1D* srs_la[13];
    TH1D* alrs_la[13];
    
    TH1D* lrb_la[13];
    TH1D* lrs_la[13];
    
    for(int i=0;i<9;i++){
        //sr Y
        _file3->GetObject(Form("pp_HM105_GplusPP/signallambda_bkg_pt%d",i),signal_la[i]);
        _file3->GetObject(Form("pp_HM105_GplusPP/backgroundlambda_bkg_pt%d",i),background_la[i]);
        _file3->GetObject(Form("pp_HM105_GplusPP/mult_la_bkg_pt%d",i),mult_la[i]);
        _file3->GetObject(Form("pp_HM105_GplusPP/mult_ass",i),mult_assoc_la[i]);
        Nassoc_la[i] = mult_assoc_la[i]->GetMean(1);
        
        TF1* fit2 = new TF1("fit2","[0]*x^2+[1]*x+[2]",0.6,2.2);
        fit2->SetParameters(1,1,1);
        fit2->SetLineColor(2);
        TF1* fit = new TF1("fit","[0]*x^2+[1]*x+[2]",0.0,2.0);
        fit->SetParameters(1,1,1);
        fit->SetLineColor(2);
        
        nEvent_la[i] = mult_la[i]->Integral(2,10000);
        Bz_la[i] = background_la[i]->GetBinContent(background_la[i]->FindBin(0,0));
        c->cd();
        (TH1D*)srs_la[i] = signal_la[i]->ProjectionY(Form("srsla%d",i),14,20);
        TH1D* srb = background_la[i]->ProjectionY("srb",14,20);
        srs_la[i]->Divide(srb);
        srs_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);
        
        srs_la[i]->Fit("fit2","R");
        srs_la[i]->Fit("fit2","R");
        srs_la[i]->Fit("fit2","R");
        srs_la[i]->Fit("fit2","R");
        srs_la[i]->Fit("fit2","R");
        
        double srm = fit2->GetMinimum(0.6,2.2);
        double srmx = fit2->GetMinimumX(0.6,2.2);
        TF1* mfsr = new TF1("mfsr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        mfsr->SetParameter(0,-srm);
        TH1D* srslow = (TH1D*)srs_la[i]->Clone();
        srslow->Add(mfsr);
        c1->cd(i+1);
        srs_la[i]->Draw();
        tex->DrawLatex(0.52,0.88,"CMS pp 7TeV");
        tex->DrawLatex(0.52,0.82,"10<N_{trk}^{offline}<20");
        tex->DrawLatex(0.52,0.74,Form("%s<p_{T}^{trg}<%s",Nplot[i],Nplot[i+1]));
        tex->DrawLatex(0.52,0.67,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(0.52,0.59,"|#Delta#eta|<1");
        c->cd();
        sry_la[i] = srslow->IntegralAndError(srslow->FindBin(0.0),srslow->FindBin(srmx),srye_la[i],"width");
        double bin0yield = srslow->GetBinContent(srslow->FindBin(0.0))*0.19635;
        sry_la[i] = sry_la[i]*2 - bin0yield;
        //cout<<"111"<<endl;
        //lr Y
        (TH1D*)alrs_la[i] = signal_la[i]->ProjectionY(Form("alrsla%d",i),1,10);
        alrs1 = signal_la[i]->ProjectionY("alrs1",24,33);
        alrb = background_la[i]->ProjectionY("alrb",1,10);
        alrb1 = background_la[i]->ProjectionY("alrb1",24,33);
        alrs_la[i]->Add(alrs1);
        alrb->Add(alrb1);
        alrs_la[i]->Divide(alrb);
        alrs_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);
        
        alrs_la[i]->Fit("fit","R");
        alrs_la[i]->Fit("fit","R");
        alrs_la[i]->Fit("fit","R");
        alrs_la[i]->Fit("fit","R");
        alrs_la[i]->Fit("fit","R");
        
        double lrm = fit->GetMinimum(0.0,2.0);
        double lrmx = fit->GetMinimumX(0.0,2.0);
        TF1* mflr = new TF1("mflr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
        mflr->SetParameter(0,-lrm);
        TH1D* alrslow = (TH1D*)alrs_la[i]->Clone();
        alrslow->Add(mflr);
        c2->cd(i+1);
        alrs_la[i]->Draw();
        tex->DrawLatex(0.22,0.88,"CMS pp 7TeV");
        tex->DrawLatex(0.22,0.82,"10<N_{trk}^{offline}<20");
        tex->DrawLatex(0.22,0.74,Form("%s<p_{T}^{trg}<%s",Nplot[i],Nplot[i+1]));
        tex->DrawLatex(0.22,0.67,"0.3<p_{T}^{assoc}<3GeV");
        tex->DrawLatex(0.22,0.59,"|#Delta#eta|>2");
        c->cd();
        lry_la[i] = alrslow->IntegralAndError(alrslow->FindBin(0.0),alrslow->FindBin(lrmx),lrye_la[i],"width");
        double bin0yield = alrslow->GetBinContent(alrslow->FindBin(0.0))*0.19635;
        lry_la[i] = lry_la[i]*2 - bin0yield;
        
        //sub Y
        suby_la[i] = sry_la[i] - lry_la[i];
        //suby_low[i] = sry_low[i];
        subye_la[i] = sqrt(srye_la[i]*srye_la[i]+lrye_la[i]*lrye_la[i]);
        c5->cd(i+1);
        
        /*TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        fit1->SetParNames("N","V1","V2","V3","V4");
        fit1->SetParameters(10,1,1,1,1);
        fit1->SetLineColor(2);*/
        
        TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
        fit1->SetParNames("N","V1","V2","V3");
        fit1->SetParameters(10,1,1,1);
        fit1->SetLineColor(2);
        
        //lr V2 obs
        lrs_la[i] = signal_la[i]->ProjectionY(Form("lrs_la%d",i),1,10);
        TH1D* lrsa = signal_la[i]->ProjectionY("lrsa",24,33);
        lrb_la[i] = background_la[i]->ProjectionY(Form("lrb_la%d",i),1,10);
        TH1D* lrba = background_la[i]->ProjectionY("lrba",24,33);
        lrs_la[i]->Add(lrsa);
        lrb_la[i]->Add(lrba);
        lrs_la[i]->Divide(lrb_la[i]);
        lrs_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);
        //lrs_low[i]->Scale(1.0/4.2);
        lrs_la[i]->Fit("fit1","R");
        lrs_la[i]->Fit("fit1","R");
        lrs_la[i]->Fit("fit1","R");
        lrs_la[i]->Fit("fit1","R");
        
        V2_la[i] = fit1->GetParameter(2);
        V2e_la[i] = fit1->GetParError(2)*errfactor;
        V3_la[i] = fit1->GetParameter(3);
        V3e_la[i] = fit1->GetParError(3)*errfactor;
        V1_la[i] = fit1->GetParameter(1);
        V1e_la[i] = fit1->GetParError(1)*errfactor;
        
        Nassoc_fit_la[i] = fit1->GetParameter(0);
        c->cd();
    }
    
    //==================ref Vn after jet correction=====================
    double V2sub_ref = V2_ref - V2_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low;
    double V2sube_ref = sqrt(V2e_ref*V2e_ref + V2e_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low*V2e_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low);
    
    double V3sub_ref = V3_ref - V3_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low;
    double V3sube_ref = sqrt(V3e_ref*V3e_ref + V3e_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low*V3e_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low);
    
    double V1sub_ref = V1_ref - V1_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low;
    double V1sube_ref = sqrt(V1e_ref*V1e_ref + V1e_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low*V1e_ref_low*Nassoc_ref_fit_low/Nassoc_ref_fit*sry_ref/sry_ref_low);
    
    double v2sub_ref = sqrt(V2sub_ref);
    double v2sube_ref = sqrt(V2sub_ref)*(V2sube_ref/V2sub_ref)/2;

    double v3sub_ref = sqrt(V3sub_ref);
    double v3sube_ref = sqrt(V3sub_ref)*(V3sube_ref/V3sub_ref)/2;
    
    double v1sub_ref = sqrt(V1sub_ref);
    double v1sube_ref = sqrt(V1sub_ref)*(V1sube_ref/V1sub_ref)/2;
    
    //=================pT differential Vn for obs after jet correction=====================
    //==========KS============
    double V2sub_ks[13];
    double V2sube_ks[13];
    
    double V3sub_ks[13];
    double V3sube_ks[13];

    double V1sub_ks[13];
    double V1sube_ks[13];

    double jetYfactor_ks[13];
    
    for(int i=0;i<9;i++)
    {
	Nassoc_low_ks[i] = Nassoc_fit_low_ks[i];
	Nassoc_ks[i] = Nassoc_fit_ks[i];
        V2sub_ks[i] = V2_ks[i] - V2_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*suby_ks[i]/suby_low_ks[i];
        V2sube_ks[i] = sqrt(V2e_ks[i]*V2e_ks[i] + sqrt(V2e_low_ks[i]/V2_low_ks[i]*V2e_low_ks[i]/V2_low_ks[i] + subye_ks[i]/suby_ks[i]*subye_ks[i]/suby_ks[i] + subye_low_ks[i]/suby_low_ks[i]*subye_low_ks[i]/suby_low_ks[i])*V2_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*suby_ks[i]/suby_low_ks[i]*sqrt(V2e_low_ks[i]/V2_low_ks[i]*V2e_low_ks[i]/V2_low_ks[i] + subye_ks[i]/suby_ks[i]*subye_ks[i]/suby_ks[i] + subye_low_ks[i]/suby_low_ks[i]*subye_low_ks[i]/suby_low_ks[i])*V2_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*suby_ks[i]/suby_low_ks[i]);
        
        V3sub_ks[i] = V3_ks[i] - V3_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*suby_ks[i]/suby_low_ks[i];
        V3sube_ks[i] = sqrt(V3e_ks[i]*V3e_ks[i] + sqrt(V3e_low_ks[i]/V3_low_ks[i]*V3e_low_ks[i]/V3_low_ks[i] + subye_ks[i]/suby_ks[i]*subye_ks[i]/suby_ks[i] + subye_low_ks[i]/suby_low_ks[i]*subye_low_ks[i]/suby_low_ks[i])*V3_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*suby_ks[i]/suby_low_ks[i]*sqrt(V3e_low_ks[i]/V3_low_ks[i]*V3e_low_ks[i]/V3_low_ks[i] + subye_ks[i]/suby_ks[i]*subye_ks[i]/suby_ks[i] + subye_low_ks[i]/suby_low_ks[i]*subye_low_ks[i]/suby_low_ks[i])*V3_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*suby_ks[i]/suby_low_ks[i]);

        V1sub_ks[i] = V1_ks[i] - V1_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*suby_ks[i]/suby_low_ks[i];
        V1sube_ks[i] = sqrt(V1e_ks[i]*V1e_ks[i] + sqrt(V1e_low_ks[i]/V1_low_ks[i]*V1e_low_ks[i]/V1_low_ks[i] + subye_ks[i]/suby_ks[i]*subye_ks[i]/suby_ks[i] + subye_low_ks[i]/suby_low_ks[i]*subye_low_ks[i]/suby_low_ks[i])*V1_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*suby_ks[i]/suby_low_ks[i]*sqrt(V1e_low_ks[i]/V1_low_ks[i]*V1e_low_ks[i]/V1_low_ks[i] + subye_ks[i]/suby_ks[i]*subye_ks[i]/suby_ks[i] + subye_low_ks[i]/suby_low_ks[i]*subye_low_ks[i]/suby_low_ks[i])*V1_low_ks[i]*Nassoc_low_ks[i]/Nassoc_ks[i]*suby_ks[i]/suby_low_ks[i]);
        
        jetYfactor_ks[i] = suby_ks[i]/suby_low_ks[i];
        
    }
    
    //==========LA============
    double V2sub_la[13];
    double V2sube_la[13];
    
    double V3sub_la[13];
    double V3sube_la[13];
    
    double V1sub_la[13];
    double V1sube_la[13];
    
    double jetYfactor_la[13];
    
    for(int i=0;i<9;i++)
    {
        Nassoc_low_la[i] = Nassoc_fit_low_la[i];
        Nassoc_la[i] = Nassoc_fit_la[i];
        V2sub_la[i] = V2_la[i] - V2_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*suby_la[i]/suby_low_la[i];
        V2sube_la[i] = sqrt(V2e_la[i]*V2e_la[i] + sqrt(V2e_low_la[i]/V2_low_la[i]*V2e_low_la[i]/V2_low_la[i] + subye_la[i]/suby_la[i]*subye_la[i]/suby_la[i] + subye_low_la[i]/suby_low_la[i]*subye_low_la[i]/suby_low_la[i])*V2_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*suby_la[i]/suby_low_la[i]*sqrt(V2e_low_la[i]/V2_low_la[i]*V2e_low_la[i]/V2_low_la[i] + subye_la[i]/suby_la[i]*subye_la[i]/suby_la[i] + subye_low_la[i]/suby_low_la[i]*subye_low_la[i]/suby_low_la[i])*V2_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*suby_la[i]/suby_low_la[i]);
        
        V3sub_la[i] = V3_la[i] - V3_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*suby_la[i]/suby_low_la[i];
        V3sube_la[i] = sqrt(V3e_la[i]*V3e_la[i] + sqrt(V3e_low_la[i]/V3_low_la[i]*V3e_low_la[i]/V3_low_la[i] + subye_la[i]/suby_la[i]*subye_la[i]/suby_la[i] + subye_low_la[i]/suby_low_la[i]*subye_low_la[i]/suby_low_la[i])*V3_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*suby_la[i]/suby_low_la[i]*sqrt(V3e_low_la[i]/V3_low_la[i]*V3e_low_la[i]/V3_low_la[i] + subye_la[i]/suby_la[i]*subye_la[i]/suby_la[i] + subye_low_la[i]/suby_low_la[i]*subye_low_la[i]/suby_low_la[i])*V3_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*suby_la[i]/suby_low_la[i]);
        
        V1sub_la[i] = V1_la[i] - V1_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*suby_la[i]/suby_low_la[i];
        V1sube_la[i] = sqrt(V1e_la[i]*V1e_la[i] + sqrt(V1e_low_la[i]/V1_low_la[i]*V1e_low_la[i]/V1_low_la[i] + subye_la[i]/suby_la[i]*subye_la[i]/suby_la[i] + subye_low_la[i]/suby_low_la[i]*subye_low_la[i]/suby_low_la[i])*V1_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*suby_la[i]/suby_low_la[i]*sqrt(V1e_low_la[i]/V1_low_la[i]*V1e_low_la[i]/V1_low_la[i] + subye_la[i]/suby_la[i]*subye_la[i]/suby_la[i] + subye_low_la[i]/suby_low_la[i]*subye_low_la[i]/suby_low_la[i])*V1_low_la[i]*Nassoc_low_la[i]/Nassoc_la[i]*suby_la[i]/suby_low_la[i]);
        
        jetYfactor_la[i] = suby_la[i]/suby_low_la[i];
        
    }

    //=========================pT differential vn for obs=================
    //=========KS===========
    double v2sub_ks[13];
    double v2sube_ks[13];
    double v2_ks[13];
    double v2e_ks[13];
    double v2_low_ks[13];
    double v2e_low_ks[13];
    
    double v3sub_ks[13];
    double v3sube_ks[13];
    double v3_ks[13];
    double v3e_ks[13];
    double v3_low_ks[13];
    double v3e_low_ks[13];
    
    double v1sub_ks[13];
    double v1sube_ks[13];
    double v1_ks[13];
    double v1e_ks[13];
    double v1_low_ks[13];
    double v1e_low_ks[13];
    
    for(int i=0;i<9;i++)
    {
        v2sub_ks[i] = V2sub_ks[i]/v2sub_ref;
        v2sube_ks[i] = fabs(sqrt(V2sube_ks[i]/V2sub_ks[i]*V2sube_ks[i]/V2sub_ks[i] + v2sube_ref/v2sub_ref*v2sube_ref/v2sub_ref)*v2sub_ks[i]);
        
        v2_ks[i] = V2_ks[i]/v2_ref;
        v2e_ks[i] = sqrt(V2e_ks[i]/V2_ks[i]*V2e_ks[i]/V2_ks[i] + v2e_ref/v2_ref*v2e_ref/v2_ref)*v2_ks[i];
        v2_low_ks[i] = V2_low_ks[i]/v2_ref_low;
        v2e_low_ks[i] = sqrt(V2e_low_ks[i]/V2_low_ks[i]*V2e_low_ks[i]/V2_low_ks[i] + v2e_ref_low/v2_ref_low*v2e_ref_low/v2_ref_low)*v2_low_ks[i];
        
        v3sub_ks[i] = V3sub_ks[i]/v3sub_ref;
        v3sube_ks[i] = fabs(sqrt(V3sube_ks[i]/V3sub_ks[i]*V3sube_ks[i]/V3sub_ks[i] + v3sube_ref/v3sub_ref*v3sube_ref/v3sub_ref)*v3sub_ks[i]);
        
        v3_ks[i] = V3_ks[i]/v3_ref;
        v3e_ks[i] = sqrt(V3e_ks[i]/V3_ks[i]*V3e_ks[i]/V3_ks[i] + v3e_ref/v3_ref*v3e_ref/v3_ref)*v3_ks[i];
        v3_low_ks[i] = V3_low_ks[i]/v3_ref_low;
        v3e_low_ks[i] = sqrt(V3e_low_ks[i]/V3_low_ks[i]*V3e_low_ks[i]/V3_low_ks[i] + v3e_ref_low/v3_ref_low*v3e_ref_low/v3_ref_low)*v3_low_ks[i];

        v1sub_ks[i] = V1sub_ks[i]/v1sub_ref;
        v1sube_ks[i] = fabs(sqrt(V1sube_ks[i]/V1sub_ks[i]*V1sube_ks[i]/V1sub_ks[i] + v1sube_ref/v1sub_ref*v1sube_ref/v1sub_ref)*v1sub_ks[i]);
        
        v1_ks[i] = V1_ks[i]/v1_ref;
        v1e_ks[i] = sqrt(V1e_ks[i]/V1_ks[i]*V1e_ks[i]/V1_ks[i] + v1e_ref/v1_ref*v1e_ref/v1_ref)*v1_ks[i];
        v1_low_ks[i] = V1_low_ks[i]/v1_ref_low;
        v1e_low_ks[i] = sqrt(V1e_low_ks[i]/V1_low_ks[i]*V1e_low_ks[i]/V1_low_ks[i] + v1e_ref_low/v1_ref_low*v1e_ref_low/v1_ref_low)*v1_low_ks[i];
    }
    
    //=========LA===========
    double v2sub_la[13];
    double v2sube_la[13];
    double v2_la[13];
    double v2e_la[13];
    double v2_low_la[13];
    double v2e_low_la[13];
    
    double v3sub_la[13];
    double v3sube_la[13];
    double v3_la[13];
    double v3e_la[13];
    double v3_low_la[13];
    double v3e_low_la[13];
    
    double v1sub_la[13];
    double v1sube_la[13];
    double v1_la[13];
    double v1e_la[13];
    double v1_low_la[13];
    double v1e_low_la[13];
    
    for(int i=0;i<9;i++)
    {
        v2sub_la[i] = V2sub_la[i]/v2sub_ref;
        v2sube_la[i] = fabs(sqrt(V2sube_la[i]/V2sub_la[i]*V2sube_la[i]/V2sub_la[i] + v2sube_ref/v2sub_ref*v2sube_ref/v2sub_ref)*v2sub_la[i]);
        
        v2_la[i] = V2_la[i]/v2_ref;
        v2e_la[i] = sqrt(V2e_la[i]/V2_la[i]*V2e_la[i]/V2_la[i] + v2e_ref/v2_ref*v2e_ref/v2_ref)*v2_la[i];
        v2_low_la[i] = V2_low_la[i]/v2_ref_low;
        v2e_low_la[i] = sqrt(V2e_low_la[i]/V2_low_la[i]*V2e_low_la[i]/V2_low_la[i] + v2e_ref_low/v2_ref_low*v2e_ref_low/v2_ref_low)*v2_low_la[i];
        
        v3sub_la[i] = V3sub_la[i]/v3sub_ref;
        v3sube_la[i] = fabs(sqrt(V3sube_la[i]/V3sub_la[i]*V3sube_la[i]/V3sub_la[i] + v3sube_ref/v3sub_ref*v3sube_ref/v3sub_ref)*v3sub_la[i]);
        
        v3_la[i] = V3_la[i]/v3_ref;
        v3e_la[i] = sqrt(V3e_la[i]/V3_la[i]*V3e_la[i]/V3_la[i] + v3e_ref/v3_ref*v3e_ref/v3_ref)*v3_la[i];
        v3_low_la[i] = V3_low_la[i]/v3_ref_low;
        v3e_low_la[i] = sqrt(V3e_low_la[i]/V3_low_la[i]*V3e_low_la[i]/V3_low_la[i] + v3e_ref_low/v3_ref_low*v3e_ref_low/v3_ref_low)*v3_low_la[i];
        
        v1sub_la[i] = V1sub_la[i]/v1sub_ref;
        v1sube_la[i] = fabs(sqrt(V1sube_la[i]/V1sub_la[i]*V1sube_la[i]/V1sub_la[i] + v1sube_ref/v1sub_ref*v1sube_ref/v1sub_ref)*v1sub_la[i]);
        
        v1_la[i] = V1_la[i]/v1_ref;
        v1e_la[i] = sqrt(V1e_la[i]/V1_la[i]*V1e_la[i]/V1_la[i] + v1e_ref/v1_ref*v1e_ref/v1_ref)*v1_la[i];
        v1_low_la[i] = V1_low_la[i]/v1_ref_low;
        v1e_low_la[i] = sqrt(V1e_low_la[i]/V1_low_la[i]*V1e_low_la[i]/V1_low_la[i] + v1e_ref_low/v1_ref_low*v1e_ref_low/v1_ref_low)*v1_low_la[i];
        
        v2sub_la[0] = -999;
        v2la[0] = -999;
    }
    
    //======================Get pT of each bin===================
    //===========KS============
    TH1D* hPt_ks[13];
    
    double pt_ks[13];
    
    for(int i=0;i<9;i++)
    {
        _file3->GetObject(Form("pp_HM105_GplusPP/Ptkshort_bkg_pt%d",i),hPt_ks[i]);

        pt_ks[i] = hPt_ks[i]->GetMean(1);
        
    }
    //===========LA============
    TH1D* hPt_la[13];
    
    double pt_la[13];
    
    for(int i=0;i<9;i++)
    {
        _file3->GetObject(Form("pp_HM105_GplusPP/Ptlambda_bkg_pt%d",i),hPt_la[i]);
        
        pt_la[i] = hPt_la[i]->GetMean(1);
        pt_la[0] = 0;
    }

    //======================Get KET of each bin===================
    //===========KS============
    TH1D* hKET_ks[13];
    
    double KET_ks[13];
    
    for(int i=0;i<9;i++)
    {
        _file3->GetObject(Form("pp_HM105_GplusPP/KETkshort_bkg_pt%d",i),hKET_ks[i]);
        
        KET_ks[i] = hKET_ks[i]->GetMean(1);
        
    }
    //===========LA============
    TH1D* hKET_la[13];
    
    double KET_la[13];
    
    for(int i=0;i<9;i++)
    {
        _file3->GetObject(Form("pp_HM105_GplusPP/KETlambda_bkg_pt%d",i),hKET_la[i]);
        
        KET_la[i] = hKET_la[i]->GetMean(1);
        KET_la[0] = 0;
        
    }

    //======================Get perisub factor V2*N/Yjet===============
    //=======KS===========
	double perisubfactor_ks[13];
	double perisubfactore_ks[13];
	for(int i=0;i<9;i++)
	{
		perisubfactor_ks[i] = V2_low_ks[i]*Nassoc_low_ks[i]/suby_low_ks[i];
		perisubfactore_ks[i] = sqrt((V2e_low_ks[i]/V2_low_ks[i])*(V2e_low_ks[i]/V2_low_ks[i]) + subye_low_ks[i]/suby_low_ks[i]*subye_low_ks[i]/suby_low_ks[i])*V2_low_ks[i]*Nassoc_low_ks[i]/suby_low_ks[i];
	}
    //=======LA===========
    double perisubfactor_la[13];
    double perisubfactore_la[13];
    for(int i=0;i<9;i++)
    {
        perisubfactor_la[i] = V2_low_la[i]*Nassoc_low_la[i]/suby_low_la[i];
        perisubfactore_la[i] = sqrt((V2e_low_la[i]/V2_low_la[i])*(V2e_low_la[i]/V2_low_la[i]) + subye_low_la[i]/suby_low_la[i]*subye_low_la[i]/suby_low_la[i])*V2_low_la[i]*Nassoc_low_la[i]/suby_low_la[i];
    }

    //========================output file==============================
    TFile ofile("V0v2_vspt_sub1020_NassFit_highpt_pt033_bkg_3terms_7TeVEff.root","RECREATE");
    
    //========================create TGraph for obs====================
    //========KS=========
    TGraphErrors* V2plot_ks = new TGraphErrors(9,pt_ks,V2_ks,0,V2e_ks);
    TGraphErrors* V2_subplot_ks = new TGraphErrors(9,pt_ks,V2sub_ks,0,V2sube_ks);
    TGraphErrors* V2plot_ks_KET = new TGraphErrors(11,KET_ks,V2_ks,0,V2e_ks);
    TGraphErrors* V2_subplot_ks_KET = new TGraphErrors(11,KET_ks,V2sub_ks,0,V2sube_ks);
    TGraphErrors* V3plot_ks = new TGraphErrors(9,pt_ks,V3_ks,0,V3e_ks);
    TGraphErrors* V3_subplot_ks = new TGraphErrors(9,pt_ks,V3sub_ks,0,V3sube_ks);
    TGraphErrors* V1plot_ks = new TGraphErrors(9,pt_ks,V1_ks,0,V2e_ks);
    TGraphErrors* V1_subplot_ks = new TGraphErrors(9,pt_ks,V1sub_ks,0,V1sube_ks);
    
    V2plot_ks->SetName("kshortV2_bkg_GplusPP");
    V2_subplot_ks->SetName("kshortV2sub_bkg_GplusPP");
    
    V2plot_ks->Write();
    V2_subplot_ks->Write();
    
    V2plot_ks_KET->SetName("kshortV2_bkg_GplusPP_KET");
    V2_subplot_ks_KET->SetName("kshortV2sub_bkg_GplusPP_KET");
    
    V2plot_ks_KET->Write();
    V2_subplot_ks_KET->Write();

    V3plot_ks->SetName("kshortV3_bkg_GplusPP");
    V3_subplot_ks->SetName("kshortV3sub_bkg_GplusPP");
    
    V3plot_ks->Write();
    V3_subplot_ks->Write();
    
    V1plot_ks->SetName("kshortV1_bkg_GplusPP");
    V1_subplot_ks->SetName("kshortV1sub_bkg_GplusPP");
    
    V1plot_ks->Write();
    V1_subplot_ks->Write();

    TGraphErrors* v2plot_ks = new TGraphErrors(9,pt_ks,v2_ks,0,v2e_ks);
    TGraphErrors* v2_subplot_ks = new TGraphErrors(9,pt_ks,v2sub_ks,0,v2sube_ks);
    TGraphErrors* v2plot_ks_KET = new TGraphErrors(11,KET_ks,v2_ks,0,v2e_ks);
    TGraphErrors* v2_subplot_ks_KET = new TGraphErrors(11,KET_ks,v2sub_ks,0,v2sube_ks);
    TGraphErrors* v3plot_ks = new TGraphErrors(9,pt_ks,v3_ks,0,v3e_ks);
    TGraphErrors* v3_subplot_ks = new TGraphErrors(9,pt_ks,v3sub_ks,0,v3sube_ks);
    TGraphErrors* v1plot_ks = new TGraphErrors(9,pt_ks,v1_ks,0,v1e_ks);
    TGraphErrors* v1_subplot_ks = new TGraphErrors(9,pt_ks,v1sub_ks,0,v1sube_ks);

    v2plot_ks->SetName("kshortv2_bkg_GplusPP");
    v2_subplot_ks->SetName("kshortv2sub_bkg_GplusPP");

    v2plot_ks->Write();
    v2_subplot_ks->Write();
    
    v2plot_ks_KET->SetName("kshortv2_bkg_GplusPP_KET");
    v2_subplot_ks_KET->SetName("kshortv2sub_bkg_GplusPP_KET");
    
    v2plot_ks_KET->Write();
    v2_subplot_ks_KET->Write();
    
    v3plot_ks->SetName("kshortv3_bkg_GplusPP");
    v3_subplot_ks->SetName("kshortv3sub_bkg_GplusPP");
    
    v3plot_ks->Write();
    v3_subplot_ks->Write();

    v1plot_ks->SetName("kshortv1_bkg_GplusPP");
    v1_subplot_ks->SetName("kshortv1sub_bkg_GplusPP");
    
    v1plot_ks->Write();
    v1_subplot_ks->Write();

    TGraphErrors* sryplot_ks = new TGraphErrors(9,pt_ks,sry_ks,0,srye_ks);
    TGraphErrors* lryplot_ks = new TGraphErrors(9,pt_ks,lry_ks,0,lrye_ks);
    TGraphErrors* subyplot_ks = new TGraphErrors(9,pt_ks,suby_ks,0,subye_ks);
	TGraphErrors* perisubfactorplot_ks = new TGraphErrors(9,pt_ks,perisubfactor_ks,0,perisubfactore_ks);

	sryplot_ks->SetName("sry_ks_bkg_GplusPP");
	lryplot_ks->SetName("lry_ks_bkg_GplusPP");
	subyplot_ks->SetName("suby_ks_bkg_GplusPP");
	perisubfactorplot_ks->SetName("perisubfactor_ks_bkg_GplusPP");

	sryplot_ks->Write();
	lryplot_ks->Write();
	subyplot_ks->Write();
	perisubfactorplot_ks->Write();

    //========LA=========
    TGraphErrors* V2plot_la = new TGraphErrors(9,pt_la,V2_la,0,V2e_la);
    TGraphErrors* V2_subplot_la = new TGraphErrors(9,pt_la,V2sub_la,0,V2sube_la);
    TGraphErrors* V2plot_la_KET = new TGraphErrors(11,KET_la,V2_la,0,V2e_la);
    TGraphErrors* V2_subplot_la_KET = new TGraphErrors(11,KET_la,V2sub_la,0,V2sube_la);
    TGraphErrors* V3plot_la = new TGraphErrors(9,pt_la,V3_la,0,V3e_la);
    TGraphErrors* V3_subplot_la = new TGraphErrors(9,pt_la,V3sub_la,0,V3sube_la);
    TGraphErrors* V1plot_la = new TGraphErrors(9,pt_la,V1_la,0,V2e_la);
    TGraphErrors* V1_subplot_la = new TGraphErrors(9,pt_la,V1sub_la,0,V1sube_la);
    
    V2plot_la->SetName("lambdaV2_bkg_GplusPP");
    V2_subplot_la->SetName("lambdaV2sub_bkg_GplusPP");
    
    V2plot_la->Write();
    V2_subplot_la->Write();
    
    V2plot_la_KET->SetName("lambdaV2_bkg_GplusPP_KET");
    V2_subplot_la_KET->SetName("lambdaV2sub_bkg_GplusPP_KET");
    
    V2plot_la_KET->Write();
    V2_subplot_la_KET->Write();
    
    V3plot_la->SetName("lambdaV3_bkg_GplusPP");
    V3_subplot_la->SetName("lambdaV3sub_bkg_GplusPP");
    
    V3plot_la->Write();
    V3_subplot_la->Write();
    
    V1plot_la->SetName("lambdaV1_bkg_GplusPP");
    V1_subplot_la->SetName("lambdaV1sub_bkg_GplusPP");
    
    V1plot_la->Write();
    V1_subplot_la->Write();
    
    TGraphErrors* v2plot_la = new TGraphErrors(9,pt_la,v2_la,0,v2e_la);
    TGraphErrors* v2_subplot_la = new TGraphErrors(9,pt_la,v2sub_la,0,v2sube_la);
    TGraphErrors* v2plot_la_KET = new TGraphErrors(11,KET_la,v2_la,0,v2e_la);
    TGraphErrors* v2_subplot_la_KET = new TGraphErrors(11,KET_la,v2sub_la,0,v2sube_la);
    TGraphErrors* v3plot_la = new TGraphErrors(9,pt_la,v3_la,0,v3e_la);
    TGraphErrors* v3_subplot_la = new TGraphErrors(9,pt_la,v3sub_la,0,v3sube_la);
    TGraphErrors* v1plot_la = new TGraphErrors(9,pt_la,v1_la,0,v1e_la);
    TGraphErrors* v1_subplot_la = new TGraphErrors(9,pt_la,v1sub_la,0,v1sube_la);
    
    v2plot_la->SetName("lambdav2_bkg_GplusPP");
    v2_subplot_la->SetName("lambdav2sub_bkg_GplusPP");
    
    v2plot_la->Write();
    v2_subplot_la->Write();
    
    v2plot_la_KET->SetName("lambdav2_bkg_GplusPP_KET");
    v2_subplot_la_KET->SetName("lambdav2sub_bkg_GplusPP_KET");
    
    v2plot_la_KET->Write();
    v2_subplot_la_KET->Write();
    
    v3plot_la->SetName("lambdav3_bkg_GplusPP");
    v3_subplot_la->SetName("lambdav3sub_bkg_GplusPP");
    
    v3plot_la->Write();
    v3_subplot_la->Write();
    
    v1plot_la->SetName("lambdav1_bkg_GplusPP");
    v1_subplot_la->SetName("lambdav1sub_bkg_GplusPP");
    
    v1plot_la->Write();
    v1_subplot_la->Write();
    
    TGraphErrors* sryplot_la = new TGraphErrors(9,pt_la,sry_la,0,srye_la);
    TGraphErrors* lryplot_la = new TGraphErrors(9,pt_la,lry_la,0,lrye_la);
    TGraphErrors* subyplot_la = new TGraphErrors(9,pt_la,suby_la,0,subye_la);
    TGraphErrors* perisubfactorplot_la = new TGraphErrors(9,pt_la,perisubfactor_la,0,perisubfactore_la);
    
    sryplot_la->SetName("sry_la_bkg_GplusPP");
    lryplot_la->SetName("lry_la_bkg_GplusPP");
    subyplot_la->SetName("suby_la_bkg_GplusPP");
    perisubfactorplot_la->SetName("perisubfactor_la_bkg_GplusPP");
    
    sryplot_la->Write();
    lryplot_la->Write();
    subyplot_la->Write();
    perisubfactorplot_la->Write();
    
    return;
    
    for(int i=0;i<9;i++)
    {
        cout<<"jet factor: "<<jetYfactor[i]<<endl;
    }
    
    cout<<"Nassoc: "<<Nassoc<<"  Nassoc_low: "<<Nassoc_low<<endl;
    
cout<<"jet yield: high:"<<suby_ref<<" low:"<<suby_ref_low<<endl;
    cout<<"sr yield: high:"<<sry_ref<<" low:"<<sry_ref_low<<endl;
    cout<<"lr yield: high:"<<lry_ref<<" low:"<<lry_ref_low<<endl;

    
    cout<<"Ref V2: before sub:"<<V2_ref<<" after sub:"<<V2sub_ref<<"  V2low:"<<V2_ref_low<<endl;
cout<<"Ref v2: before sub:"<<v2_ref<<" after sub:"<<v2sub_ref<<"  v2low:"<<v2_ref_low<<endl;
    cout<<"V2 vs pT before sub:"<<endl;
    for(int i=0;i<9;i++)
    {
        cout<<v2[i]<<"   "<<v2e[i]<<endl;
	//cout<<"sry:"<<sry[i]<<endl;
	//cout<<"sry_low:"<<sry_low[i]<<endl;
	//cout<<"lry:"<<lry[i]<<endl;
	//cout<<"lry_low:"<<lry_low[i]<<endl;
	cout<<"Nassfit:"<<Nassoc[i]<<endl;
	cout<<"Nassfit_low:"<<Nassoc_low[i]<<endl;
    }
    cout<<"V2 vs pT after sub:"<<endl;
    for(int i=0;i<9;i++)
    {
        cout<<v2sub[i]<<"   "<<v2sube[i]<<endl;
    }
    cout<<"V2_low:"<<endl;
    for(int i=0;i<9;i++)
    {
        cout<<v2_low[i]<<"   "<<v2e_low[i]<<endl;
    }
    
}
