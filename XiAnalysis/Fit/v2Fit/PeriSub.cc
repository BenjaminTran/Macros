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
    //std::string fn = "V0v2perisubFixedWindow.root";
    //std::string fn = "FitRootFiles/PeriSub/V0v2perisub_Default_AllStrange_EG1_0_35_CorrectRef_CorrectGap_FullStats_5percentReject_1_11_18.root"; USing Old V0 data
    //std::string fn = "FitRootFiles/PeriSub/V0v2perisub_Default_AllStrange_EG1_0_35_CorrectRef_CorrectGap_FullStats_5percentReject_1_10_18.root";
    //std::string fn = "FitRootFiles/PeriSub/V0v2perisub_Default_AllStrange_EG1_0_35_CorrectRef_CorrectGap_FullStats_8p5-10GeV_1_10_18.root";
    std::string fn = "FitRootFiles/PeriSub/V0v2perisub_Default_EG1_0_35_V0wTrkEff_01_11_18.root";
    //std::string fn = "FitRootFiles/PeriSub/V0v2perisub_Default_EG1_0_35_CorrectRef_InCorrectGap_12_05_17.root";
    //std::string fn = "FitRootFiles/PeriSub/V0v2perisub_Default_EG1_0_35_JetPeak1p64_12_17_17.root";
    //std::string fn = "FitRootFiles/PeriSub/V0v2perisub_FixedWindow1p7_EG1_0_35_12_17_17.root";
    //std::string fn = "FitRootFiles/PeriSub/V0v2perisub_FixedWindow1p2_EG1_0_35_12_17_17.root";
    //std::string fn = "FitRootFiles/PeriSub/V0v2perisub_FixedWindow0p5_EG1_0_35_12_17_17.root";
    //std::string fn = "JustForPrintingNassoc.root";
    //std::string fn_V0 = "v0CorrelationRapidity"; //For original HM V0
    std::string fn_V0 = "v0CasCorrelationRapidity";
    //std::string fn_Xi = "xiCorrelationRapidity";
    std::string fn_Xi = "v0CasCorrelationRapidity";
    std::string fn_Om = "v0CasCorrelationRapidity";
    std::string fn_v0cas = "v0CasCorrelationRapidityPeriSub";
    // For 8.5-10 Pt Bin
    //std::vector<double> fsig_ks    = {0.989638 , 0.992109 , 0.992861 , 0.992984 , 0.992159 , 0.989499 , 0.986445 , 0.982275 , 0.976946 , 0.970546 , 0.964845 , 0.958563 , 0.969557 , 0.949};
    //std::vector<double> fsig_ks_MB = {0.997    , 0.997    , 0.993    , 0.997    , 0.995    , 0.993    , 0.991    , 0.989    , 0.986    , 0.982    , 0.978    , 0.971    , 0.968, 0.965};
    //std::vector<double> fsig_la    = {0.908266 , 0.96779 , 0.979912 , 0.981899 , 0.982888 , 0.982854 , 0.980766 , 0.97569 , 0.97569 , 0.964048, 0.945};
    //std::vector<double> fsig_la_MB = {0.996    , 0.994   , 0.995    , 0.994    , 0.992    , 0.990    , 0.987    , 0.981   , 0.973   , 0.963 , 0.952};
    //std::vector<double> PtBin_ks   = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5, 10.0};
    //std::vector<double> PtBin_la   = {0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5, 10.0};

    // For 5% rejection
    //std::vector<double> fsig_ks = {0.9918, 0.9928, 0.9934, 0.9936, 0.9929, 0.9906, 0.9879, 0.9841, 0.9793, 0.9736, 0.9685, 0.9633, 0.9577};
    //std::vector<double> fsig_la = {0.8979, 0.9977, 0.9811, 0.9833, 0.9835, 0.9847, 0.9831, 0.9785, 0.9702, 0.9609};

    // Default
    std::vector<double> fsig_ks    = {0.989638 , 0.992109 , 0.992861 , 0.992984 , 0.992159 , 0.989499 , 0.986445 , 0.982275 , 0.976946 , 0.970546 , 0.964845 , 0.958563 , 0.969557};
    std::vector<double> fsig_ks_MB = {0.997    , 0.997    , 0.993    , 0.997    , 0.995    , 0.993    , 0.991    , 0.989    , 0.986    , 0.982    , 0.978    , 0.971    , 0.968};
    std::vector<double> fsig_la    = {0.908266 , 0.96779 , 0.979912 , 0.981899 , 0.982888 , 0.982854 , 0.980766 , 0.97569 , 0.97569 , 0.964048};
    std::vector<double> fsig_la_MB = {0.996    , 0.994   , 0.995    , 0.994    , 0.992    , 0.990    , 0.987    , 0.981   , 0.973   , 0.963};
    std::vector<double> fsig_xi    = {0.959427 , 0.976239 , 0.979161 , 0.980678 , 0.980661 , 0.981534 , 0.981502 , 0.979289 , 0.979192};
    std::vector<double> fsig_xi_MB = {0.967  , 0.977 , 0.980 , 0.983   , 0.982 , 0.980 , 0.983 , 0.981 , 0.971};
    std::vector<double> fsig_om    = {0.778226 ,0.874337, 0.922863, 0.953094, 0.966833};// ,0.976012}; //For Peripheral Sub
    std::vector<double> fsig_om_MB = {0.896448 ,0.941455 ,0.956581 ,0.966235 ,0.978242};// ,0.976012}; //For Peripheral Sub
    std::vector<double> PtBin_ks   = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5};//, 10.0, 15.0, 20.0};
    std::vector<double> PtBin_la   = {0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5};//, 10.0, 15.0, 20.0};
    std::vector<double> PtBin_xi   = {1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.2,10.0};
    std::vector<double> PtBin_om   = {1.5, 2.2, 2.8, 3.6, 5.0, 8.0};//, 20.0}; // PbPb
    TFile *f_perisub      = TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/MBCorr/V0MB_0_35_1_02_18.root");
    //TFile *f_perisub      = TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/MBCorr/V0MB_w8p5-10pt_1_08_18.root");
    //TFile *f_perisub = TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/MBCorr/PeripheralSubtractionMB_0_n_20_V0Only.root");
    //TFile *f_perisub_xi   = TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/MBCorr/XiOmegaMB_0_N_20_Partial_11_8_17.root");
    TFile *f_perisub_xi   = TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/MBCorr/XiMB_0_35_1_02_18.root");
    //TFile *f_perisub_xi   = TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/MBCorr/XiMB_partial_0_35_12_05_17.root");
    TFile *f_perisub_om   = TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/MBCorr/OmegaMB_Total_0_35_MergedBin1-2_11_28_17.root");
    //TFile *f_V0           = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityCorrectMultB_09_19_17.root"); Old No trk eff DO NOT USE
    TFile *f_V0           = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0Correlation_TrkEff_1_11_18.root");
    //TFile *f_V0           = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0_5percentDump_ARC3_1_10_18.root"); //Throw away 5% of data as requested by ARC3
    //TFile *f_V0           = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0_OLD_5percentDump_ARC3_1_11_18.root"); //Throw away 5% of data as requested by ARC3 (FROM OLD FILE)
    //TFile* f_Xi         = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationRapidityTotal_08_20_2017.root");
    TFile* f_Xi           = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationHM_11_07_17.root");
    //TFile *f_low_ref      = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/HadCorr/Combine_MB0_corr_ref.root"); //0-20
    TFile *f_low_ref      = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/HadCorr/Combine_MB0_ref.root"); //0-35
    TFile* f_Om = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/OmCorr/OmegaHMRapidity_Total_12_04_17.root");
    //TFile *f_high_ref   = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/HadCorr/XiCorrelationRapidityTotal_08_20_2017.root");
    TFile *f_high_ref     = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/HadCorr/Combine_HM185_corr_ref_PUrej.root");
    int binlow_V2 = 13;
    int binhigh_V2 = 21;
    // Only need to worry about fixd window if using fixed window function
    double fw_low = -1.2; //Fixed window low
    double fw_high= 1.2;
    double sr_low = 14; //14 = -1, 12 = -1.64
    double sr_high = 20; // 20 = 1, 22 = 1.64
} V0;

void PeriSubFixedWindow(ParticleData PD)
{
    TH1::SetDefaultSumw2();

    double BW2D = (9.9/33)*((2-1.0/16)*3.1416/31);

    TLatex* tex = new TLatex();
    tex->SetNDC();

    TCanvas* c_sr_low_ref;
    TCanvas* c_lr_low_ref;
    TCanvas* c_sr_high_ref;
    TCanvas* c_lr_high_ref;
    TCanvas* c_sr_low_ks[2];
    TCanvas* c_lr_low_ks[2];
    TCanvas* c_sr_ks[2];
    TCanvas* c_lr_ks[2];
    TCanvas* c_sr_low_la[2];
    TCanvas* c_lr_low_la[2];
    TCanvas* c_sr_la[2];
    TCanvas* c_lr_la[2];
    TCanvas* c_sr_low_xi[2];
    TCanvas* c_lr_low_xi[2];
    TCanvas* c_sr_xi[2];
    TCanvas* c_lr_xi[2];
    TCanvas* c_sr_low_om[2];
    TCanvas* c_lr_low_om[2];
    TCanvas* c_sr_om[2];
    TCanvas* c_lr_om[2];

    TCanvas* c_low_2Dks_1[2];
    TCanvas* c_low_2Dks_2[2];

    TCanvas* c_high_2Dks_1[2];
    TCanvas* c_high_2Dks_2[2];

    TCanvas* c_low_2Dla_1[2];
    TCanvas* c_low_2Dla_2[2];

    TCanvas* c_high_2Dla_1[2];
    TCanvas* c_high_2Dla_2[2];

    TCanvas* c_low_2Dxi_1[2];
    TCanvas* c_low_2Dxi_2[2];

    TCanvas* c_high_2Dxi_1[2];
    TCanvas* c_high_2Dxi_2[2];

    TCanvas* c_low_2Dom_1[2];
    TCanvas* c_low_2Dom_2[2];

    TCanvas* c_high_2Dom_1[2];
    TCanvas* c_high_2Dom_2[2];

    c_sr_low_ref  = new TCanvas("csPeaksr_low_ref" , "cssr_low_ref"  , 600  , 600);
    c_lr_low_ref  = new TCanvas("csPeaklr_low_ref" , "cslr_low_ref"  , 600  , 600);
    c_sr_high_ref = new TCanvas("csPeaksr_high_ref", "cssr_high_ref" , 600  , 600);
    c_lr_high_ref = new TCanvas("csPeaklr_high_ref", "cslr_high_ref" , 600  , 600);
    for(int j=0; j<2; j++)
    {
        std::string region;
        if( j==0 ) region = "Peak";
        if( j==1 ) region = "Side";
        c_sr_low_ks[j]   = new TCanvas(("cs" + region + "sr_low_ks"   ).c_str(), ("cs" + region + "sr_low_ks"   ).c_str(), 1200 , 900);
        c_lr_low_ks[j]   = new TCanvas(("cs" + region + "lr_low_ks"   ).c_str(), ("cs" + region + "lr_low_ks"   ).c_str(), 1200 , 900);
        c_sr_ks[j]       = new TCanvas(("cs" + region + "sr_ks"       ).c_str(), ("cs" + region + "sr_ks"       ).c_str(), 1200 , 900);
        c_lr_ks[j]       = new TCanvas(("cs" + region + "lr_ks"       ).c_str(), ("cs" + region + "lr_ks"       ).c_str(), 1200 , 900);
        c_sr_low_la[j]   = new TCanvas(("cs" + region + "sr_low_la"   ).c_str(), ("cs" + region + "sr_low_la"   ).c_str(), 1200 , 900);
        c_lr_low_la[j]   = new TCanvas(("cs" + region + "lr_low_la"   ).c_str(), ("cs" + region + "lr_low_la"   ).c_str(), 1200 , 900);
        c_sr_la[j]       = new TCanvas(("cs" + region + "sr_la"       ).c_str(), ("cs" + region + "sr_la"       ).c_str(), 1200 , 900);
        c_lr_la[j]       = new TCanvas(("cs" + region + "lr_la"       ).c_str(), ("cs" + region + "lr_la"       ).c_str(), 1200 , 900);
        c_sr_low_xi[j]   = new TCanvas(("cs" + region + "sr_low_xi"   ).c_str(), ("cs" + region + "sr_low_xi"   ).c_str(), 1200 , 900);
        c_lr_low_xi[j]   = new TCanvas(("cs" + region + "lr_low_xi"   ).c_str(), ("cs" + region + "lr_low_xi"   ).c_str(), 1200 , 900);
        c_sr_xi[j]       = new TCanvas(("cs" + region + "sr_xi"       ).c_str(), ("cs" + region + "sr_xi"       ).c_str(), 1200 , 900);
        c_lr_xi[j]       = new TCanvas(("cs" + region + "lr_xi"       ).c_str(), ("cs" + region + "lr_xi"       ).c_str(), 1200 , 900);
        c_sr_low_om[j]   = new TCanvas(("cs" + region + "sr_low_om"   ).c_str(), ("cs" + region + "sr_low_om"   ).c_str(), 1200 , 900);
        c_lr_low_om[j]   = new TCanvas(("cs" + region + "lr_low_om"   ).c_str(), ("cs" + region + "lr_low_om"   ).c_str(), 1200 , 900);
        c_sr_om[j]       = new TCanvas(("cs" + region + "sr_om"       ).c_str(), ("cs" + region + "sr_om"       ).c_str(), 1200 , 900);
        c_lr_om[j]       = new TCanvas(("cs" + region + "lr_om"       ).c_str(), ("cs" + region + "lr_om"       ).c_str(), 1200 , 900);

        c_low_2Dks_1[j] = new TCanvas(("cs" + region + "_low_2Dks_1").c_str(),("cs" + region + "_low_2Dks_1").c_str(),1200,900);
        c_low_2Dks_2[j] = new TCanvas(("cs" + region + "_low_2Dks_2").c_str(),("cs" + region + "_low_2Dks_2").c_str(),1200,900);

        c_high_2Dks_1[j] = new TCanvas(("cs" + region + "_high_2Dks_1").c_str(),("cs" + region + "_high_2Dks_1").c_str(),1200,900);
        c_high_2Dks_2[j] = new TCanvas(("cs" + region + "_high_2Dks_2").c_str(),("cs" + region + "_high_2Dks_2").c_str(),1200,900);

        c_low_2Dla_1[j] = new TCanvas(("cs" + region + "_low_2Dla_1").c_str(),("cs" + region + "_low_2Dla_1").c_str(),1200,900);
        c_low_2Dla_2[j] = new TCanvas(("cs" + region + "_low_2Dla_2").c_str(),("cs" + region + "_low_2Dla_2").c_str(),1200,900);

        c_high_2Dla_1[j] = new TCanvas(("cs" + region + "_high_2Dla_1").c_str(),("cs" + region + "_high_2Dla_1").c_str(),1200,900);
        c_high_2Dla_2[j] = new TCanvas(("cs" + region + "_high_2Dla_2").c_str(),("cs" + region + "_high_2Dla_2").c_str(),1200,900);

        c_low_2Dxi_1[j] = new TCanvas(("cs" + region + "_low_2Dxi_1").c_str(),("cs" + region + "_low_2Dxi_1").c_str(),1200,900);
        c_low_2Dxi_2[j] = new TCanvas(("cs" + region + "_low_2Dxi_2").c_str(),("cs" + region + "_low_2Dxi_2").c_str(),1200,900);

        c_high_2Dxi_1[j] = new TCanvas(("cs" + region + "_high_2Dxi_1").c_str(),("cs" + region + "_high_2Dxi_1").c_str(),1200,900);
        c_high_2Dxi_2[j] = new TCanvas(("cs" + region + "_high_2Dxi_2").c_str(),("cs" + region + "_high_2Dxi_2").c_str(),1200,900);

        c_low_2Dom_1[j] = new TCanvas(("cs" + region + "_low_2Dom_1").c_str(),("cs" + region + "_low_2Dom_1").c_str(),1200,900);
        c_low_2Dom_2[j] = new TCanvas(("cs" + region + "_low_2Dom_2").c_str(),("cs" + region + "_low_2Dom_2").c_str(),1200,900);

        c_high_2Dom_1[j] = new TCanvas(("cs" + region + "_high_2Dom_1").c_str(),("cs" + region + "_high_2Dom_1").c_str(),1200,900);
        c_high_2Dom_2[j] = new TCanvas(("cs" + region + "_high_2Dom_2").c_str(),("cs" + region + "_high_2Dom_2").c_str(),1200,900);

        c_sr_low_ks[j]->Divide(4,4);
        c_lr_low_ks[j]->Divide(4,4);
        c_sr_ks[j]->Divide(4,4);
        c_lr_ks[j]->Divide(4,4);
        c_sr_low_la[j]->Divide(4,3);
        c_lr_low_la[j]->Divide(4,3);
        c_sr_la[j]->Divide(4,3);
        c_lr_la[j]->Divide(4,3);
        c_sr_low_xi[j]->Divide(4,3);
        c_lr_low_xi[j]->Divide(4,3);
        c_sr_xi[j]->Divide(4,3);
        c_lr_xi[j]->Divide(4,3);
        c_sr_low_om[j]->Divide(4,3);
        c_lr_low_om[j]->Divide(4,3);
        c_sr_om[j]->Divide(4,3);
        c_lr_om[j]->Divide(4,3);

        c_low_2Dks_1[j]  -> Divide(4,2);
        c_low_2Dks_2[j]  -> Divide(4,2);
        c_high_2Dks_1[j] -> Divide(4,2);
        c_high_2Dks_2[j] -> Divide(4,2);

        c_low_2Dla_1[j]  -> Divide(4,2);
        c_low_2Dla_2[j]  -> Divide(4,2);
        c_high_2Dla_1[j] -> Divide(4,2);
        c_high_2Dla_2[j] -> Divide(4,2);

        c_low_2Dxi_1[j]  -> Divide(4,2);
        c_low_2Dxi_2[j]  -> Divide(4,2);
        c_high_2Dxi_1[j] -> Divide(4,2);
        c_high_2Dxi_2[j] -> Divide(4,2);

        c_low_2Dom_1[j]  -> Divide(4,2);
        c_low_2Dom_2[j]  -> Divide(4,2);
        c_high_2Dom_1[j] -> Divide(4,2);
        c_high_2Dom_2[j] -> Divide(4,2);
    }
    TCanvas* c = new TCanvas("c","c",400,400);



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

    hbPeaksr_low_ref = hbackgroundPeak_low_ref->ProjectionY("hbPeaksr_low_ref", PD.sr_low, PD.sr_high);
    hsPeaksr_low_ref = hsignalPeak_low_ref->ProjectionY("hsPeaksr_low_ref", PD.sr_low, PD.sr_high);

    //nEvent_low_ref.push_back(mult_low_ref->Integral(0,100000));
    nEvent_low_ref.push_back(mult_low_ref->Integral(3,100000));
    Bz_low_ref.push_back(hbackgroundPeak_low_ref->GetBinContent(hbackgroundPeak_low_ref->FindBin(0,0)));

    hsPeaksr_low_ref->Divide(hbPeaksr_low_ref);
    hsPeaksr_low_ref->Scale(Bz_low_ref[0]/nEvent_low_ref[0]/BW2D);

    hbPeaklr_low_ref = hbackgroundPeak_low_ref->ProjectionY("hbPeaklr_low_ref",1,10);
    TH1D* ahbPeaklr_low_ref = hbackgroundPeak_low_ref->ProjectionY("ahbPeaklr_low_ref",24,33);
    hsPeaklr_low_ref = hsignalPeak_low_ref->ProjectionY("hsPeaklr_low_ref",1,10);
    TH1D* ahsPeaklr_low_ref = hsignalPeak_low_ref->ProjectionY("ahsPeaklr_low_ref",24,33);

    hbPeaklr_low_ref->Add(ahbPeaklr_low_ref);
    hsPeaklr_low_ref->Add(ahsPeaklr_low_ref);
    hsPeaklr_low_ref->Divide(hbPeaklr_low_ref);
    hsPeaklr_low_ref->Scale(Bz_low_ref[0]/nEvent_low_ref[0]/BW2D);

    hsPeaksr_low_ref->Add(hsPeaklr_low_ref,-1);

    c_lr_low_ref->cd();
    TH2D* hbackgroundPeak_low_ref_clone = (TH2D*)hbackgroundPeak_low_ref->Clone("hbackgroundPeak_low_ref_clone");
    TH2D* hsignalPeak_low_ref_clone = (TH2D*)hsignalPeak_low_ref->Clone("hsignalPeak_low_ref_clone");
    hsignalPeak_low_ref_clone->Divide(hbackgroundPeak_low_ref_clone);
    hsignalPeak_low_ref_clone->Draw("Surf1");

    c->cd();

    TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.0,2.0);
    quadFit2->SetParameters(1,1,1);

    hsPeaksr_low_ref->Fit("quadFit2","R");
    hsPeaksr_low_ref->Fit("quadFit2","R");
    hsPeaksr_low_ref->Fit("quadFit2","R");

    double minVal_sr = quadFit2->GetMinimum(0.0,2.0);
    double minVal_srX = quadFit2->GetMinimumX(0.0,2.0);
    TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    minConst_sr->SetParameter(0,-minVal_sr);
    TH1D* hsPeaksr_zeroed_low_ref = (TH1D*)hsPeaksr_low_ref->Clone();
    hsPeaksr_zeroed_low_ref->Add(minConst_sr);
    c_sr_low_ref->cd();
    hsPeaksr_low_ref->Draw();
    double xcoor = 0.52;
    double ycoor = 0.90;
    double increment = 0.07;
    tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
    tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",0.3,3.0));
    tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|>2");
    c->cd();
    Jyieldsr_low_ref.push_back(hsPeaksr_zeroed_low_ref->IntegralAndError(hsPeaksr_zeroed_low_ref->FindBin(PD.fw_low),hsPeaksr_zeroed_low_ref->FindBin(PD.fw_high),Jyieldsr_err_low_ref,"width"));
    //bin0yield = hsPeaklr_zeroed_low_ref->GetBinContent(hsPeaklr_zeroed_low_ref->FindBin(0.0))*0.19635;
    //Jyieldlr_low_ref[0] = Jyieldlr_low_ref[0]*2 - bin0yield;

    JyieldSub_low_ref.push_back(Jyieldsr_low_ref[0]);
    JyieldSub_err_low_ref.push_back(Jyieldsr_err_low_ref);

    TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit1->SetParNames("N","V1","V2","V3");
    fit1->SetParameters(10,1,1,1);
    fit1->SetLineColor(2);

    V2lrs_low_ref = hsignalPeak_low_ref            -> ProjectionY("V2lrs_low_ref",1,PD.binlow_V2);
    TH1D* aV2lrs_low_ref = hsignalPeak_low_ref     -> ProjectionY("aV2lrs_low_ref",PD.binhigh_V2,33);
    V2lrb_low_ref = hbackgroundPeak_low_ref        -> ProjectionY("V2lrb_low_ref",1,PD.binlow_V2);
    TH1D* aV2lrb_low_ref = hbackgroundPeak_low_ref -> ProjectionY("aV2lrb_low_ref",PD.binhigh_V2,33);
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

    //PD.f_high_ref->GetObject("xiCorrelationRapidity/BackgroundHad",hbackgroundPeak_high_ref);
    //PD.f_high_ref->GetObject("xiCorrelationRapidity/SignalHad",hsignalPeak_high_ref);
    //TH1D* mult_high_ref = (TH1D*) PD.f_high_ref->Get("xiCorrelationRapidity/nTrk");
    PD.f_high_ref->GetObject("pPbCorr/background",hbackgroundPeak_high_ref);
    PD.f_high_ref->GetObject("pPbCorr/signal",hsignalPeak_high_ref);
    TH1D* mult_high_ref = (TH1D*) PD.f_high_ref->Get("pPbCorr/mult");

    hbPeaksr_high_ref = hbackgroundPeak_high_ref->ProjectionY("hbPeaksr_high_ref", PD.sr_low, PD.sr_high);
    hsPeaksr_high_ref = hsignalPeak_high_ref->ProjectionY("hsPeaksr_high_ref", PD.sr_low, PD.sr_high);

    //nEvent_high_ref.push_back(mult_high_ref->Integral(0,100000));
    nEvent_high_ref.push_back(mult_high_ref->Integral(3,100000));
    Bz_high_ref.push_back(hbackgroundPeak_high_ref->GetBinContent(hbackgroundPeak_high_ref->FindBin(0,0)));

    hsPeaksr_high_ref->Divide(hbPeaksr_high_ref);
    hsPeaksr_high_ref->Scale(Bz_high_ref[0]/nEvent_high_ref[0]/BW2D);

    hbPeaklr_high_ref = hbackgroundPeak_high_ref->ProjectionY("hbPeaklr_high_ref",1,10);
    TH1D* ahbPeaklr_high_ref = hbackgroundPeak_high_ref->ProjectionY("ahbPeaklr_high_ref",24,33);
    hsPeaklr_high_ref = hsignalPeak_high_ref->ProjectionY("hsPeaklr_high_ref",1,10);
    TH1D* ahsPeaklr_high_ref = hsignalPeak_high_ref->ProjectionY("ahsPeaklr_high_ref",24,33);

    hbPeaklr_high_ref->Add(ahbPeaklr_high_ref);
    hsPeaklr_high_ref->Add(ahsPeaklr_high_ref);
    hsPeaklr_high_ref->Divide(hbPeaklr_high_ref);
    hsPeaklr_high_ref->Scale(Bz_high_ref[0]/nEvent_high_ref[0]/BW2D);

    hsPeaksr_high_ref->Add(hsPeaklr_high_ref,-1);

    TF1* quadFit21 = new TF1("quadFit21","[0]*x^2+[1]*x+[2]",0.0,2.0);
    quadFit21->SetParameters(1,1,1);

    hsPeaksr_high_ref->Fit("quadFit21","R");
    hsPeaksr_high_ref->Fit("quadFit21","R");
    hsPeaksr_high_ref->Fit("quadFit21","R");

    minVal_sr = quadFit21->GetMinimum(0.6,2.0);
    minVal_srX = quadFit21->GetMinimumX(0.6,2.0);
    TF1* minConst_sr2 = new TF1("minConst_sr2","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    minConst_sr2->SetParameter(0,-minVal_sr);
    TH1D* hsPeaksr_zeroed_high_ref = (TH1D*)hsPeaksr_high_ref->Clone();
    hsPeaksr_zeroed_high_ref->Add(minConst_sr2);
    c_sr_high_ref->cd();
    hsPeaksr_high_ref->Draw();
    xcoor = 0.52;
    ycoor = 0.90;
    increment = 0.07;
    tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
    tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",0.3,3.0));
    tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|>2");
    c->cd();
    Jyieldsr_high_ref.push_back(hsPeaksr_zeroed_high_ref->IntegralAndError(hsPeaksr_zeroed_high_ref->FindBin(PD.fw_low),hsPeaksr_zeroed_high_ref->FindBin(PD.fw_high),Jyieldsr_err_high_ref,"width"));
    //bin0yield = hsPeaksr_zeroed_high_ref->GetBinContent(hsPeaksr_zeroed_high_ref->FindBin(0.0))*0.19635;
    //Jyieldsr_high_ref[0] = Jyieldsr_high_ref[0]*2 - bin0yield;

    JyieldSub_high_ref.push_back(Jyieldsr_high_ref[0]);
    JyieldSub_err_high_ref.push_back(Jyieldsr_err_high_ref);

    TF1* fit2 = new TF1("fit2","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit2->SetParNames("N","V1","V2","V3");
    fit2->SetParameters(10,1,1,1);
    fit2->SetLineColor(2);

    V2lrs_high_ref = hsignalPeak_high_ref->ProjectionY("V2lrs_high_ref",1,PD.binlow_V2);
    TH1D* aV2lrs_high_ref = hsignalPeak_high_ref->ProjectionY("aV2lrs_high_ref",PD.binhigh_V2,33);
    V2lrb_high_ref = hbackgroundPeak_high_ref->ProjectionY("V2lrb_high_ref",1,PD.binlow_V2);
    TH1D* aV2lrb_high_ref = hbackgroundPeak_high_ref->ProjectionY("aV2lrb_high_ref",PD.binhigh_V2,33);
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
    std::string PathMult_high_xi;

    std::string PathBackground_low_om;
    std::string PathSignal_low_om;
    std::string PathMult_low_om;
    std::string PathBackground_high_om;
    std::string PathSignal_high_om;
    std::string PathMult_high_om;

    int numPtBins_ks = PD.PtBin_ks.size()-1;
    int numPtBins_la = PD.PtBin_la.size()-1;
    int numPtBins_xi = PD.PtBin_xi.size()-1;
    int numPtBins_om = PD.PtBin_om.size()-1;

    std::string region_label;


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
            PathBackground_high_xi  = PD.fn_Xi + "/BackgroundXiPeak_pt%d";
            PathSignal_high_xi = PD.fn_Xi + "/SignalXiPeak_pt%d";
            PathMult_high_xi = PD.fn_Xi + "/mult_xi_pt%d";

            PathBackground_low_om  = PD.fn_v0cas + "/BackgroundOmPeak_pt%d";
            PathSignal_low_om = PD.fn_v0cas + "/SignalOmPeak_pt%d";
            PathMult_low_om = PD.fn_v0cas + "/mult_om_pt%d";
            PathBackground_high_om  = PD.fn_Om + "/BackgroundOmPeak_pt%d";
            PathSignal_high_om = PD.fn_Om + "/SignalOmPeak_pt%d";
            PathMult_high_om = PD.fn_Om + "/mult_om_pt%d";

            region_label = "Peak";
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
            PathBackground_high_xi  = PD.fn_Xi + "/BackgroundXiSide_pt%d";
            PathSignal_high_xi = PD.fn_Xi + "/SignalXiSide_pt%d";
            PathMult_high_xi = PD.fn_Xi + "/mult_xi_bkg_pt%d";

            PathBackground_low_om  = PD.fn_v0cas + "/BackgroundOmSide_pt%d";
            PathSignal_low_om = PD.fn_v0cas + "/SignalOmSide_pt%d";
            PathMult_low_om = PD.fn_v0cas + "/mult_om_bkg_pt%d";
            PathBackground_high_om  = PD.fn_Om + "/BackgroundOmSide_pt%d";
            PathSignal_high_om = PD.fn_Om + "/SignalOmSide_pt%d";
            PathMult_high_om = PD.fn_Om + "/mult_om_bkg_pt%d";

            region_label = "Side";
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

        TH1D* hsPeaksr_low_ks[3*arraySize_ks];
        TH1D* hbPeaksr_low_ks[3*arraySize_ks];
        TH1D* hsPeaklr_low_ks[3*arraySize_ks];
        TH1D* hbPeaklr_low_ks[3*arraySize_ks];
        TH1D* V2lrs_low_ks[3*arraySize_ks];
        TH1D* V2lrb_low_ks[3*arraySize_ks];
        TH2D* hbackgroundPeak_low_ks[3*arraySize_ks];
        TH2D* hsignalPeak_low_ks[3*arraySize_ks];


        //Calculate Nassoc, Jet yield, Peak region Low N

        for(int i=0; i<numPtBins_ks; i++)
        {
            PD.f_perisub->GetObject(Form(PathBackground_low_ks.c_str(),i),hbackgroundPeak_low_ks[i]);
            PD.f_perisub->GetObject(Form(PathSignal_low_ks.c_str(),i),hsignalPeak_low_ks[i]);
            TH1D* mult_low_ks = (TH1D*) PD.f_perisub->Get(Form(PathMult_low_ks.c_str(),i));

            hbPeaksr_low_ks[i] = hbackgroundPeak_low_ks[i]->ProjectionY(Form("hbPeaksr_low_ks%d",i), PD.sr_low, PD.sr_high);
            hsPeaksr_low_ks[i] = hsignalPeak_low_ks[i]->ProjectionY(Form(("hs"+region_label+"sr_low_ks%d").c_str(),i), PD.sr_low, PD.sr_high);

            TF1* quadFit13 = new TF1("quadFit13","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit13->SetParameters(1,1,1);

            if(i<8 && j==0)
            {
                c_low_2Dks_1[j]->cd(i+1);
                //TH2D* hPeak_low_2Dks = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_ks.c_str(),i));
                TH2D* hPeak_low_2Dks = (TH2D*)hsignalPeak_low_ks[i]->Clone();
                hPeak_low_2Dks->SetName(Form("hPeak_low_2Dks%d",i));
                hPeak_low_2Dks->Divide(hbackgroundPeak_low_ks[i]);

                hPeak_low_2Dks->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dks->Draw("SURF1");
                //hbackgroundPeak_low_ks[i]->Draw();
            }
            if(i>=8 && j==0)
            {
                c_low_2Dks_2[j]->cd(i-8+1);
                //TH2D* hPeak_low_2Dks = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_ks.c_str(),i));
                TH2D* hPeak_low_2Dks = (TH2D*)hsignalPeak_low_ks[i]->Clone();
                hPeak_low_2Dks->SetName(Form("hPeak_low_2Dks%d",i));
                hPeak_low_2Dks->Divide(hbackgroundPeak_low_ks[i]);
                hPeak_low_2Dks->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dks->Draw("SURF1");
            }

            if(i<8 && j==1)
            {
                c_low_2Dks_1[j]->cd(i+1);
                //TH2D* hPeak_low_2Dks = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_ks.c_str(),i));
                TH2D* hPeak_low_2Dks = (TH2D*)hsignalPeak_low_ks[i]->Clone();
                hPeak_low_2Dks->SetName(Form("hPeak_low_2Dks%d",i));
                hPeak_low_2Dks->Divide(hbackgroundPeak_low_ks[i]);

                hPeak_low_2Dks->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dks->Draw("SURF1");
                //hbackgroundPeak_low_ks[i]->Draw();
            }
            if(i>=8 && j==1)
            {
                c_low_2Dks_2[j]->cd(i-8+1);
                //TH2D* hPeak_low_2Dks = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_ks.c_str(),i));
                TH2D* hPeak_low_2Dks = (TH2D*)hsignalPeak_low_ks[i]->Clone();
                hPeak_low_2Dks->SetName(Form("hPeak_low_2Dks%d",i));
                hPeak_low_2Dks->Divide(hbackgroundPeak_low_ks[i]);
                hPeak_low_2Dks->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dks->Draw("SURF1");
            }
            c->cd();


            //nEvent_low_ks.push_back(mult_low_ks->Integral(0,100000));
            nEvent_low_ks.push_back(mult_low_ks->Integral(2,100000));
            Bz_low_ks.push_back(hbackgroundPeak_low_ks[i]->GetBinContent(hbackgroundPeak_low_ks[i]->FindBin(0,0)));

            hsPeaksr_low_ks[i]->Divide(hbPeaksr_low_ks[i]);
            hsPeaksr_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);
            //hsPeaksr_low_ks[i]->Scale(Bz_low_ks[i]/BW2D);

            hbPeaklr_low_ks[i] = hbackgroundPeak_low_ks[i]->ProjectionY("hbPeaklr_low_ks",1,10);
            TH1D* ahbPeaklr_low_ks = hbackgroundPeak_low_ks[i]->ProjectionY("ahbPeaklr_low_ks",24,33);
            c_lr_low_ks[j]->cd(i+1);
            hsPeaklr_low_ks[i] = hsignalPeak_low_ks[i]->ProjectionY("hsPeaklr_low_ks",1,10);
            TH1D* ahsPeaklr_low_ks = hsignalPeak_low_ks[i]->ProjectionY("ahsPeaklr_low_ks",24,33);

            hbPeaklr_low_ks[i]->Add(ahbPeaklr_low_ks);
            hsPeaklr_low_ks[i]->Add(ahsPeaklr_low_ks);
            hsPeaklr_low_ks[i]->Divide(hbPeaklr_low_ks[i]);
            hsPeaklr_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);

            hsPeaksr_low_ks[i]->Add(hsPeaklr_low_ks[i],-1);

            hsPeaksr_low_ks[i]->Fit("quadFit13","R");
            hsPeaksr_low_ks[i]->Fit("quadFit13","R");
            hsPeaksr_low_ks[i]->Fit("quadFit13","R");

            c_sr_low_ks[j]->cd(i+1);
            double minVal_sr = quadFit13->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit13->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            TH1D* hsPeaksr_zeroed_low_ks = (TH1D*)hsPeaksr_low_ks[i]->Clone();
            hsPeaksr_zeroed_low_ks->Add(minConst_sr);
            hsPeaksr_low_ks[i]->Draw();
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_ks[i],PD.PtBin_ks[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_low_ks.push_back(hsPeaksr_zeroed_low_ks->IntegralAndError(hsPeaksr_zeroed_low_ks->FindBin(PD.fw_low),hsPeaksr_zeroed_low_ks->FindBin(PD.fw_high),Jyieldsr_err_low_ks[i],"width"));
            //double bin0yield = hsPeaksr_zeroed_low_ks->GetBinContent(hsPeaksr_zeroed_low_ks->FindBin(0.0))*0.19635;
            //Jyieldsr_low_ks[i] = Jyieldsr_low_ks[i]*2 - bin0yield;

            //hsPeaklr_low_ks[i]->Scale(Bz_low_ks[i]/BW2D);

            JyieldSub_low_ks.push_back(Jyieldsr_low_ks[i]);
            JyieldSub_err_low_ks.push_back(Jyieldsr_err_low_ks[i]);
            c->cd();

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_low_ks[i] = hsignalPeak_low_ks[i]->ProjectionY(Form("V2lrs_low_ks%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_low_ks = hsignalPeak_low_ks[i]->ProjectionY("aV2lrs_low_ks",PD.binhigh_V2,33);
            V2lrb_low_ks[i] = hbackgroundPeak_low_ks[i]->ProjectionY(Form("V2lrb_low_ks%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_low_ks = hbackgroundPeak_low_ks[i]->ProjectionY("aV2lrb_low_ks",PD.binhigh_V2,33);
            V2lrs_low_ks[i]->Add(aV2lrs_low_ks);
            V2lrb_low_ks[i]->Add(aV2lrb_low_ks);
            V2lrs_low_ks[i]->Divide(V2lrb_low_ks[i]);
            V2lrs_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);
            //V2lrs_low_ks[i]->Scale(Bz_low_ks[i]/BW2D);

            V2lrs_low_ks[i]->Fit("fit1","R");
            //V2lrs_low_ks[i]->Fit("fit1","R");
            //V2lrs_low_ks[i]->Fit("fit1","R");

            V2Values_low_ks.push_back(fit1->GetParameter(2));
            V2Values_err_low_ks.push_back(fit1->GetParError(2));
            V3Values_low_ks.push_back(fit1->GetParameter(3));
            V3Values_err_low_ks.push_back(fit1->GetParError(3));

            Nassoc_low_ks.push_back(fit1->GetParameter(0));
            c->cd();
            //if(i==4) return;
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

            hbPeaksr_low_la[i] = hbackgroundPeak_low_la[i]->ProjectionY("hbPeaksr_low_la", PD.sr_low, PD.sr_high);
            hsPeaksr_low_la[i] = hsignalPeak_low_la[i]->ProjectionY(("hs"+region_label+"sr_low_la").c_str(), PD.sr_low, PD.sr_high);

            if(i<8 && j==0)
            {
                c_low_2Dla_1[j]->cd(i+1);
                //TH2D* hPeak_low_2Dla = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_la.c_str(),i));
                TH2D* hPeak_low_2Dla = (TH2D*)hsignalPeak_low_la[i]->Clone();
                hPeak_low_2Dla->SetName(Form("hPeak_low_2Dla%d",i));
                hPeak_low_2Dla->Divide(hbackgroundPeak_low_la[i]);

                hPeak_low_2Dla->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dla->Draw("SURF1");
                //hbackgroundPeak_low_la[i]->Draw();
            }
            if(i>=8 && j==0)
            {
                c_low_2Dla_2[j]->cd(i-8+1);
                //TH2D* hPeak_low_2Dla = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_la.c_str(),i));
                TH2D* hPeak_low_2Dla = (TH2D*)hsignalPeak_low_la[i]->Clone();
                hPeak_low_2Dla->SetName(Form("hPeak_low_2Dla%d",i));
                hPeak_low_2Dla->Divide(hbackgroundPeak_low_la[i]);
                hPeak_low_2Dla->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dla->Draw("SURF1");
            }

            if(i<8 && j==1)
            {
                c_low_2Dla_1[j]->cd(i+1);
                //TH2D* hPeak_low_2Dla = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_la.c_str(),i));
                TH2D* hPeak_low_2Dla = (TH2D*)hsignalPeak_low_la[i]->Clone();
                hPeak_low_2Dla->SetName(Form("hSide_low_2Dla%d",i));
                hPeak_low_2Dla->Divide(hbackgroundPeak_low_la[i]);

                hPeak_low_2Dla->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dla->Draw("SURF1");
                //hbackgroundPeak_low_la[i]->Draw();
            }
            if(i>=8 && j==1)
            {
                c_low_2Dla_2[j]->cd(i-8+1);
                //TH2D* hPeak_low_2Dla = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_la.c_str(),i));
                TH2D* hPeak_low_2Dla = (TH2D*)hsignalPeak_low_la[i]->Clone();
                hPeak_low_2Dla->SetName(Form("hSide_low_2Dla%d",i));
                hPeak_low_2Dla->Divide(hbackgroundPeak_low_la[i]);
                hPeak_low_2Dla->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dla->Draw("SURF1");
            }

            TF1* quadFit14 = new TF1("quadFit14","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit14->SetParameters(1,1,1);

            //nEvent_low_la.push_back(mult_low_la->Integral(0,100000));
            nEvent_low_la.push_back(mult_low_la->Integral(2,100000));
            Bz_low_la.push_back(hbackgroundPeak_low_la[i]->GetBinContent(hbackgroundPeak_low_la[i]->FindBin(0,0)));

            hsPeaksr_low_la[i]->Divide(hbPeaksr_low_la[i]);
            hsPeaksr_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);
            //hsPeaksr_low_la[i]->Scale(Bz_low_la[i]/BW2D);

            hbPeaklr_low_la[i] = hbackgroundPeak_low_la[i]->ProjectionY("hbPeaklr_low_la",1,10);
            TH1D* ahbPeaklr_low_la = hbackgroundPeak_low_la[i]->ProjectionY("ahbPeaklr_low_la",24,33);
            hsPeaklr_low_la[i] = hsignalPeak_low_la[i]->ProjectionY("hsPeaklr_low_la",1,10);
            TH1D* ahsPeaklr_low_la = hsignalPeak_low_la[i]->ProjectionY("ahsPeaklr_low_la",24,33);

            hbPeaklr_low_la[i]->Add(ahbPeaklr_low_la);
            hsPeaklr_low_la[i]->Add(ahsPeaklr_low_la);
            hsPeaklr_low_la[i]->Divide(hbPeaklr_low_la[i]);
            hsPeaklr_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);
            //hsPeaklr_low_la[i]->Scale(Bz_low_la[i]/BW2D);

            hsPeaksr_low_la[i]->Add(hsPeaklr_low_la[i],-1);

            c_sr_low_la[j]->cd(i+1);

            hsPeaksr_low_la[i]->Fit("quadFit14","R");
            hsPeaksr_low_la[i]->Fit("quadFit14","R");
            hsPeaksr_low_la[i]->Fit("quadFit14","R");

            double minVal_sr = quadFit14->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit14->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            TH1D* hsPeaksr_zeroed_low_la = (TH1D*)hsPeaksr_low_la[i]->Clone();
            TH1D* hsPeaksr_drawed_low_la = (TH1D*)hsPeaksr_low_la[i]->Clone();
            hsPeaksr_drawed_low_la->Draw();
            hsPeaksr_zeroed_low_la->Add(minConst_sr);
            //hsPeaksr_low_la[i]->Draw();
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_la[i],PD.PtBin_la[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_low_la.push_back(hsPeaksr_zeroed_low_la->IntegralAndError(hsPeaksr_zeroed_low_la->FindBin(PD.fw_low),hsPeaksr_zeroed_low_la->FindBin(PD.fw_high),Jyieldsr_err_low_la[i],"width"));
            //double bin0yield = hsPeaksr_zeroed_low_la->GetBinContent(hsPeaksr_zeroed_low_la->FindBin(0.0))*0.19635;
            //Jyieldsr_low_la[i] = Jyieldsr_low_la[i]*2 - bin0yield;

            JyieldSub_low_la.push_back(Jyieldsr_low_la[i]);
            JyieldSub_err_low_la.push_back(Jyieldsr_err_low_la[i]);

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_low_la[i] = hsignalPeak_low_la[i]->ProjectionY(Form("V2lrs_low_la%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_low_la = hsignalPeak_low_la[i]->ProjectionY("aV2lrs_low_la",PD.binhigh_V2,33);
            V2lrb_low_la[i] = hbackgroundPeak_low_la[i]->ProjectionY(Form("V2lrb_low_la%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_low_la = hbackgroundPeak_low_la[i]->ProjectionY("aV2lrb_low_la",PD.binhigh_V2,33);
            V2lrs_low_la[i]->Add(aV2lrs_low_la);
            V2lrb_low_la[i]->Add(aV2lrb_low_la);
            V2lrs_low_la[i]->Divide(V2lrb_low_la[i]);
            V2lrs_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);
            //V2lrs_low_la[i]->Scale(Bz_low_la[i]/BW2D);

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
            PD.f_perisub_xi->GetObject(Form(PathBackground_low_xi.c_str(),i),hbackgroundPeak_low_xi[i]);
            PD.f_perisub_xi->GetObject(Form(PathSignal_low_xi.c_str(),i),hsignalPeak_low_xi[i]);
            TH1D* mult_low_xi = (TH1D*) PD.f_perisub_xi->Get(Form(PathMult_low_xi.c_str(),i));

            hbPeaksr_low_xi[i] = hbackgroundPeak_low_xi[i]->ProjectionY("hbPeaksr_low_xi", PD.sr_low, PD.sr_high);
            hsPeaksr_low_xi[i] = hsignalPeak_low_xi[i]->ProjectionY(Form(("hs"+region_label+"sr_low_xi%d").c_str(),i), PD.sr_low, PD.sr_high);

            TF1* quadFit13 = new TF1("quadFit13","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit13->SetParameters(1,1,1);


            //nEvent_low_xi.push_back(mult_low_xi->Integral(0,100000));
            nEvent_low_xi.push_back(mult_low_xi->Integral(2,100000));
            Bz_low_xi.push_back(hbackgroundPeak_low_xi[i]->GetBinContent(hbackgroundPeak_low_xi[i]->FindBin(0,0)));

            hsPeaksr_low_xi[i]->Divide(hbPeaksr_low_xi[i]);
            hsPeaksr_low_xi[i]->Scale(Bz_low_xi[i]/nEvent_low_xi[i]/BW2D);
            //hsPeaksr_low_xi[i]->Scale(Bz_low_xi[i]/BW2D);

            hbPeaklr_low_xi[i] = hbackgroundPeak_low_xi[i]->ProjectionY("hbPeaklr_low_xi",1,10);
            TH1D* ahbPeaklr_low_xi = hbackgroundPeak_low_xi[i]->ProjectionY("ahbPeaklr_low_xi",24,33);
            hsPeaklr_low_xi[i] = hsignalPeak_low_xi[i]->ProjectionY("hsPeaklr_low_xi",1,10);
            TH1D* ahsPeaklr_low_xi = hsignalPeak_low_xi[i]->ProjectionY("ahsPeaklr_low_xi",24,33);

            hbPeaklr_low_xi[i]->Add(ahbPeaklr_low_xi);
            hsPeaklr_low_xi[i]->Add(ahsPeaklr_low_xi);
            hsPeaklr_low_xi[i]->Divide(hbPeaklr_low_xi[i]);
            hsPeaklr_low_xi[i]->Scale(Bz_low_xi[i]/nEvent_low_xi[i]/BW2D);
            //hsPeaklr_low_xi[i]->Scale(Bz_low_xi[i]/BW2D);

            hsPeaksr_low_xi[i]->Add(hsPeaklr_low_xi[i],-1);

            c_sr_low_xi[j]->cd(i+1);

            hsPeaksr_low_xi[i]->Fit("quadFit13","R");
            hsPeaksr_low_xi[i]->Fit("quadFit13","R");
            hsPeaksr_low_xi[i]->Fit("quadFit13","R");

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
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_xi[i],PD.PtBin_xi[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_low_xi.push_back(hsPeaksr_zeroed_low_xi->IntegralAndError(hsPeaksr_zeroed_low_xi->FindBin(PD.fw_low),hsPeaksr_zeroed_low_xi->FindBin(PD.fw_high),Jyieldsr_err_low_xi[i],"width"));
            //double bin0yield = hsPeaksr_zeroed_low_xi->GetBinContent(hsPeaksr_zeroed_low_xi->FindBin(0.0))*0.19635;
            //Jyieldsr_low_xi[i] = Jyieldsr_low_xi[i]*2 - bin0yield;

            JyieldSub_low_xi.push_back(Jyieldsr_low_xi[i]);
            JyieldSub_err_low_xi.push_back(Jyieldsr_err_low_xi[i]);

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_low_xi[i] = hsignalPeak_low_xi[i]->ProjectionY(Form("V2lrs_low_xi%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_low_xi = hsignalPeak_low_xi[i]->ProjectionY("aV2lrs_low_xi",PD.binhigh_V2,33);
            V2lrb_low_xi[i] = hbackgroundPeak_low_xi[i]->ProjectionY(Form("V2lrb_low_xi%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_low_xi = hbackgroundPeak_low_xi[i]->ProjectionY("aV2lrb_low_xi",PD.binhigh_V2,33);
            V2lrs_low_xi[i]->Add(aV2lrs_low_xi);
            V2lrb_low_xi[i]->Add(aV2lrb_low_xi);
            V2lrs_low_xi[i]->Divide(V2lrb_low_xi[i]);
            V2lrs_low_xi[i]->Scale(Bz_low_xi[i]/nEvent_low_xi[i]/BW2D);
            //V2lrs_low_xi[i]->Scale(Bz_low_xi[i]/BW2D);

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

        //Om Low N
        std::vector<double> Nassoc_low_om;
        std::vector<double> Jyieldsr_low_om;
        std::vector<double> Jyieldlr_low_om;
        std::vector<double> Jyieldsr_err_low_om(9);
        std::vector<double> Jyieldlr_err_low_om(9);
        std::vector<double> V2Values_low_om;
        std::vector<double> V2Values_err_low_om;
        std::vector<double> V3Values_low_om;
        std::vector<double> V3Values_err_low_om;
        std::vector<double> Bz_low_om;
        std::vector<double> nEvent_low_om;

        std::vector<double> JyieldSub_low_om;
        std::vector<double> JyieldSub_err_low_om;

        int arraySize_om = PD.PtBin_om.size();

        TH1D* hsPeaksr_low_om[arraySize_om];
        TH1D* hbPeaksr_low_om[arraySize_om];
        TH1D* hsPeaklr_low_om[arraySize_om];
        TH1D* hbPeaklr_low_om[arraySize_om];
        TH1D* V2lrs_low_om[arraySize_om];
        TH1D* V2lrb_low_om[arraySize_om];
        TH2D* hbackgroundPeak_low_om[arraySize_om];
        TH2D* hsignalPeak_low_om[arraySize_om];


        //Calculate Nassoc, Jet yield, Peak region Low N

        for(int i=0; i<numPtBins_om; i++)
        {
            PD.f_perisub_om->GetObject(Form(PathBackground_low_om.c_str(),i),hbackgroundPeak_low_om[i]);
            PD.f_perisub_om->GetObject(Form(PathSignal_low_om.c_str(),i),hsignalPeak_low_om[i]);
            TH1D* mult_low_om = (TH1D*) PD.f_perisub_om->Get(Form(PathMult_low_om.c_str(),i));

            hbPeaksr_low_om[i] = hbackgroundPeak_low_om[i]->ProjectionY(Form(("hb"+region_label+"sr_low_om%d").c_str(),i), PD.sr_low, PD.sr_high);
            hsPeaksr_low_om[i] = hsignalPeak_low_om[i]->ProjectionY(Form(("hs"+region_label+"sr_low_om%d").c_str(),i), PD.sr_low, PD.sr_high);

            TF1* quadFit13 = new TF1("quadFit13","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit13->SetParameters(1,1,1);


            //nEvent_low_om.push_back(mult_low_om->Integral(0,100000));
            nEvent_low_om.push_back(mult_low_om->Integral(2,100000));
            Bz_low_om.push_back(hbackgroundPeak_low_om[i]->GetBinContent(hbackgroundPeak_low_om[i]->FindBin(0,0)));

            hsPeaksr_low_om[i]->Divide(hbPeaksr_low_om[i]);
            hsPeaksr_low_om[i]->Scale(Bz_low_om[i]/nEvent_low_om[i]/BW2D);

            hbPeaklr_low_om[i] = hbackgroundPeak_low_om[i]->ProjectionY("hbPeaklr_low_om",1,10);
            TH1D* ahbPeaklr_low_om = hbackgroundPeak_low_om[i]->ProjectionY("ahbPeaklr_low_om",24,33);
            hsPeaklr_low_om[i] = hsignalPeak_low_om[i]->ProjectionY("hsPeaklr_low_om",1,10);
            TH1D* ahsPeaklr_low_om = hsignalPeak_low_om[i]->ProjectionY("ahsPeaklr_low_om",24,33);

            hbPeaklr_low_om[i]->Add(ahbPeaklr_low_om);
            hsPeaklr_low_om[i]->Add(ahsPeaklr_low_om);
            hsPeaklr_low_om[i]->Divide(hbPeaklr_low_om[i]);
            hsPeaklr_low_om[i]->Scale(Bz_low_om[i]/nEvent_low_om[i]/BW2D);

            hsPeaksr_low_om[i]->Add(hsPeaklr_low_om[i],-1);

            c_sr_low_om[j]->cd(i+1);

            hsPeaksr_low_om[i]->Fit("quadFit13","R");
            hsPeaksr_low_om[i]->Fit("quadFit13","R");
            hsPeaksr_low_om[i]->Fit("quadFit13","R");

            double minVal_sr = quadFit13->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit13->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            TH1D* hsPeaksr_zeroed_low_om = (TH1D*)hsPeaksr_low_om[i]->Clone();
            hsPeaksr_zeroed_low_om->Add(minConst_sr);
            hsPeaksr_low_om[i]->Draw();
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_om[i],PD.PtBin_om[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_low_om.push_back(hsPeaksr_zeroed_low_om->IntegralAndError(hsPeaksr_zeroed_low_om->FindBin(PD.fw_low),hsPeaksr_zeroed_low_om->FindBin(PD.fw_high),Jyieldsr_err_low_om[i],"width"));
            //double bin0yield = hsPeaksr_zeroed_low_om->GetBinContent(hsPeaksr_zeroed_low_om->FindBin(0.0))*0.19635;
            //Jyieldsr_low_om[i] = Jyieldsr_low_om[i]*2 - bin0yield;

            JyieldSub_low_om.push_back(Jyieldsr_low_om[i]);
            JyieldSub_err_low_om.push_back(Jyieldsr_err_low_om[i]);

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_low_om[i] = hsignalPeak_low_om[i]->ProjectionY(Form("V2lrs_low_om%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_low_om = hsignalPeak_low_om[i]->ProjectionY("aV2lrs_low_om",PD.binhigh_V2,33);
            V2lrb_low_om[i] = hbackgroundPeak_low_om[i]->ProjectionY(Form("V2lrb_low_om%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_low_om = hbackgroundPeak_low_om[i]->ProjectionY("aV2lrb_low_om",PD.binhigh_V2,33);
            V2lrs_low_om[i]->Add(aV2lrs_low_om);
            V2lrb_low_om[i]->Add(aV2lrb_low_om);
            V2lrs_low_om[i]->Divide(V2lrb_low_om[i]);
            V2lrs_low_om[i]->Scale(Bz_low_om[i]/nEvent_low_om[i]/BW2D);
            //V2lrs_low_om[i]->Scale(Bz_low_om[i]/BW2D);

            V2lrs_low_om[i]->Fit("fit1","R");
            //V2lrs_low_om[i]->Fit("fit1","R");
            //V2lrs_low_om[i]->Fit("fit1","R");

            V2Values_low_om.push_back(fit1->GetParameter(2));
            V2Values_err_low_om.push_back(fit1->GetParError(2));
            V3Values_low_om.push_back(fit1->GetParameter(3));
            V3Values_err_low_om.push_back(fit1->GetParError(3));

            Nassoc_low_om.push_back(fit1->GetParameter(0));
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

            hbPeaksr_ks[i] = hbackgroundPeak_ks[i]->ProjectionY("hbPeaksr_ks", PD.sr_low, PD.sr_high);
            hsPeaksr_ks[i] = hsignalPeak_ks[i]->ProjectionY(Form(("hs"+region_label+"Peaksr_ks%d").c_str(),i), PD.sr_low, PD.sr_high);

            if(i<8 && j==0)
            {
                c_high_2Dks_1[j]->cd(i+1);
                //TH2D* hPeak_high_2Dks = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_ks.c_str(),i));
                TH2D* hPeak_high_2Dks = (TH2D*)hsignalPeak_ks[i]->Clone();
                hPeak_high_2Dks->SetName(Form("hPeak_high_2Dks%d",i));
                hPeak_high_2Dks->Divide(hbackgroundPeak_ks[i]);

                hPeak_high_2Dks->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dks->Draw("SURF1");
                //hbackgroundPeak_high_ks[i]->Draw();
            }
            if(i>=8 && j==0)
            {
                c_high_2Dks_2[j]->cd(i-8+1);
                //TH2D* hPeak_high_2Dks = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_ks.c_str(),i));
                TH2D* hPeak_high_2Dks = (TH2D*)hsignalPeak_ks[i]->Clone();
                hPeak_high_2Dks->SetName(Form("hPeak_high_2Dks%d",i));
                hPeak_high_2Dks->Divide(hbackgroundPeak_ks[i]);
                hPeak_high_2Dks->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dks->Draw("SURF1");
            }

            if(i<8 && j==1)
            {
                c_high_2Dks_1[j]->cd(i+1);
                //TH2D* hPeak_high_2Dks = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_ks.c_str(),i));
                TH2D* hPeak_high_2Dks = (TH2D*)hsignalPeak_ks[i]->Clone();
                hPeak_high_2Dks->SetName(Form("hSide_high_2Dks%d",i));
                hPeak_high_2Dks->Divide(hbackgroundPeak_ks[i]);

                hPeak_high_2Dks->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dks->Draw("SURF1");
                //hbackgroundPeak_high_ks[i]->Draw();
            }
            if(i>=8 && j==1)
            {
                c_high_2Dks_2[j]->cd(i-8+1);
                //TH2D* hPeak_high_2Dks = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_ks.c_str(),i));
                TH2D* hPeak_high_2Dks = (TH2D*)hsignalPeak_ks[i]->Clone();
                hPeak_high_2Dks->SetName(Form("hPeak_high_2Dks%d",i));
                hPeak_high_2Dks->Divide(hbackgroundPeak_ks[i]);
                hPeak_high_2Dks->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dks->Draw("SURF1");
            }

            TF1* quadFit15 = new TF1("quadFit15","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit15->SetParameters(1,1,1);

            //nEvent_ks.push_back(mult_ks->Integral(0,100000));
            nEvent_ks.push_back(mult_ks->Integral(2,100000));
            Bz_ks.push_back(hbackgroundPeak_ks[i]->GetBinContent(hbackgroundPeak_ks[i]->FindBin(0,0)));

            hsPeaksr_ks[i]->Divide(hbPeaksr_ks[i]);
            hsPeaksr_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);
            //hsPeaksr_ks[i]->Scale(Bz_ks[i]/BW2D);

            hbPeaklr_ks[i] = hbackgroundPeak_ks[i]->ProjectionY("hbPeaklr_ks",1,10);
            TH1D* ahbPeaklr_ks = hbackgroundPeak_ks[i]->ProjectionY("ahbPeaklr_ks",24,33);
            hsPeaklr_ks[i] = hsignalPeak_ks[i]->ProjectionY("hsPeaklr_ks",1,10);
            TH1D* ahsPeaklr_ks = hsignalPeak_ks[i]->ProjectionY("ahsPeaklr_ks",24,33);

            hbPeaklr_ks[i]->Add(ahbPeaklr_ks);
            hsPeaklr_ks[i]->Add(ahsPeaklr_ks);
            hsPeaklr_ks[i]->Divide(hbPeaklr_ks[i]);
            hsPeaklr_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);
            //hsPeaklr_ks[i]->Scale(Bz_ks[i]/BW2D);

            hsPeaksr_ks[i]->Add(hsPeaklr_ks[i],-1);

            c_sr_ks[j]->cd(i+1);

            hsPeaksr_ks[i]->Fit("quadFit15","R");
            hsPeaksr_ks[i]->Fit("quadFit15","R");
            hsPeaksr_ks[i]->Fit("quadFit15","R");

            double minVal_sr = quadFit15->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit15->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            TH1D* hsPeaksr_zeroed_ks = (TH1D*)hsPeaksr_ks[i]->Clone();
            TH1D* hsPeaksr_drawed_ks = (TH1D*)hsPeaksr_ks[i]->Clone();
            hsPeaksr_drawed_ks->Draw();
            hsPeaksr_zeroed_ks->Add(minConst_sr);
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_ks[i],PD.PtBin_ks[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_ks.push_back(hsPeaksr_zeroed_ks->IntegralAndError(hsPeaksr_zeroed_ks->FindBin(PD.fw_low),hsPeaksr_zeroed_ks->FindBin(PD.fw_high),Jyieldsr_err_ks[i],"width"));
            //double bin0yield = hsPeaksr_zeroed_ks->GetBinContent(hsPeaksr_zeroed_ks->FindBin(0.0))*0.19635;
            //Jyieldsr_ks[i] = Jyieldsr_ks[i]*2 - bin0yield;


            JyieldSub_ks.push_back(Jyieldsr_ks[i]);
            JyieldSub_err_ks.push_back(Jyieldsr_err_ks[i]);

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_ks[i] = hsignalPeak_ks[i]->ProjectionY(Form("V2lrs_ks%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_ks = hsignalPeak_ks[i]->ProjectionY("aV2lrs_ks",PD.binhigh_V2,33);
            V2lrb_ks[i] = hbackgroundPeak_ks[i]->ProjectionY(Form("V2lrb_ks%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_ks = hbackgroundPeak_ks[i]->ProjectionY("aV2lrb_ks",PD.binhigh_V2,33);
            V2lrs_ks[i]->Add(aV2lrs_ks);
            V2lrb_ks[i]->Add(aV2lrb_ks);
            V2lrs_ks[i]->Divide(V2lrb_ks[i]);
            V2lrs_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);
            //V2lrs_ks[i]->Scale(Bz_ks[i]/BW2D);

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
            if(i<8 && j==0)
            {
                c_high_2Dla_1[j]->cd(i+1);
                //TH2D* hPeak_high_2Dla = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_la.c_str(),i));
                TH2D* hPeak_high_2Dla = (TH2D*)hsignalPeak_la[i]->Clone();
                hPeak_high_2Dla->SetName(Form("hPeak_high_2Dla%d",i));
                hPeak_high_2Dla->Divide(hbackgroundPeak_la[i]);

                hPeak_high_2Dla->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dla->Draw("SURF1");
                //hbackgroundPeak_high_la[i]->Draw();
            }
            if(i>=8 && j==0)
            {
                c_high_2Dla_2[j]->cd(i-8+1);
                //TH2D* hPeak_high_2Dla = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_la.c_str(),i));
                TH2D* hPeak_high_2Dla = (TH2D*)hsignalPeak_la[i]->Clone();
                hPeak_high_2Dla->SetName(Form("hPeak_high_2Dla%d",i));
                hPeak_high_2Dla->Divide(hbackgroundPeak_la[i]);
                hPeak_high_2Dla->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dla->Draw("SURF1");
            }

            if(i<8 && j==1)
            {
                c_high_2Dla_1[j]->cd(i+1);
                //TH2D* hPeak_high_2Dla = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_la.c_str(),i));
                TH2D* hPeak_high_2Dla = (TH2D*)hsignalPeak_la[i]->Clone();
                hPeak_high_2Dla->SetName(Form("hSide_high_2Dla%d",i));
                hPeak_high_2Dla->Divide(hbackgroundPeak_la[i]);

                hPeak_high_2Dla->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dla->Draw("SURF1");
                //hbackgroundPeak_high_la[i]->Draw();
            }
            if(i>=8 && j==1)
            {
                c_high_2Dla_2[j]->cd(i-8+1);
                //TH2D* hPeak_high_2Dla = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_la.c_str(),i));
                TH2D* hPeak_high_2Dla = (TH2D*)hsignalPeak_la[i]->Clone();
                hPeak_high_2Dla->SetName(Form("hSide_high_2Dla%d",i));
                hPeak_high_2Dla->Divide(hbackgroundPeak_la[i]);
                hPeak_high_2Dla->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dla->Draw("SURF1");
            }

            hbPeaksr_la[i] = hbackgroundPeak_la[i]->ProjectionY("hbPeaksr_la", PD.sr_low, PD.sr_high);
            hsPeaksr_la[i] = hsignalPeak_la[i]->ProjectionY(Form(("hs"+region_label+"Peaksr_la%d").c_str(),i), PD.sr_low, PD.sr_high);

            TF1* quadFit16 = new TF1("quadFit16","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit16->SetParameters(1,1,1);

            //nEvent_la.push_back(mult_la->Integral(0,10000));
            nEvent_la.push_back(mult_la->Integral(2,10000));
            Bz_la.push_back(hbackgroundPeak_la[i]->GetBinContent(hbackgroundPeak_la[i]->FindBin(0,0)));

            hsPeaksr_la[i]->Divide(hbPeaksr_la[i]);
            hsPeaksr_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);
            //hsPeaksr_la[i]->Scale(Bz_la[i]/BW2D);

            hbPeaklr_la[i] = hbackgroundPeak_la[i]->ProjectionY("hbPeaklr_la",1,10);
            TH1D* ahbPeaklr_la = hbackgroundPeak_la[i]->ProjectionY("ahbPeaklr_la",24,33);
            hsPeaklr_la[i] = hsignalPeak_la[i]->ProjectionY("hsPeaklr_la",1,10);
            TH1D* ahsPeaklr_la = hsignalPeak_la[i]->ProjectionY("ahsPeaklr_la",24,33);

            hbPeaklr_la[i]->Add(ahbPeaklr_la);
            hsPeaklr_la[i]->Add(ahsPeaklr_la);
            hsPeaklr_la[i]->Divide(hbPeaklr_la[i]);
            hsPeaklr_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);
            //hsPeaklr_la[i]->Scale(Bz_la[i]/BW2D);

            hsPeaksr_la[i]->Add(hsPeaklr_la[i],-1);

            c_sr_la[j]->cd(i+1);

            TH1D* hsPeaksr_zeroed_la = (TH1D*)hsPeaksr_la[i]->Clone();

            hsPeaksr_la[i]->Fit("quadFit16","R");
            hsPeaksr_la[i]->Fit("quadFit16","R");
            hsPeaksr_la[i]->Fit("quadFit16","R");
            TH1D* hsPeaksr_drawed_la = (TH1D*)hsPeaksr_la[i]->Clone();

            double minVal_sr = quadFit16->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit16->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            hsPeaksr_drawed_la->Draw();
            hsPeaksr_zeroed_la->Add(minConst_sr);
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%4.1f<p_{T}^{trg}<%4.1f",PD.PtBin_la[i],PD.PtBin_la[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_la.push_back(hsPeaksr_zeroed_la->IntegralAndError(hsPeaksr_zeroed_la->FindBin(PD.fw_low),hsPeaksr_zeroed_la->FindBin(PD.fw_high),Jyieldsr_err_la[i],"width"));
            //double bin0yield = hsPeaksr_zeroed_la->GetBinContent(hsPeaksr_zeroed_la->FindBin(0.0))*0.19635;
            //Jyieldsr_la[i] = Jyieldsr_la[i]*2 - bin0yield;

            JyieldSub_la.push_back(Jyieldsr_la[i]);
            JyieldSub_err_la.push_back(Jyieldsr_err_la[i]);

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_la[i] = hsignalPeak_la[i]->ProjectionY(Form("V2lrs_la%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_la = hsignalPeak_la[i]->ProjectionY("aV2lrs_la",PD.binhigh_V2,33);
            V2lrb_la[i] = hbackgroundPeak_la[i]->ProjectionY(Form("V2lrb_la%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_la = hbackgroundPeak_la[i]->ProjectionY("aV2lrb_la",PD.binhigh_V2,33);
            V2lrs_la[i]->Add(aV2lrs_la);
            V2lrb_la[i]->Add(aV2lrb_la);
            V2lrs_la[i]->Divide(V2lrb_la[i]);
            V2lrs_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);
            //V2lrs_la[i]->Scale(Bz_la[i]/BW2D);

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
            //TH1D* mult_xi = (TH1D*) PD.f_Xi->Get((PD.fn_Xi + "/nEvtCut").c_str());
            TH1D* mult_xi = (TH1D*) PD.f_Xi->Get(Form(PathMult_high_xi.c_str(),i));

            hbPeaksr_xi[i] = hbackgroundPeak_xi[i]->ProjectionY("hbPeaksr_xi", PD.sr_low, PD.sr_high);
            hsPeaksr_xi[i] = hsignalPeak_xi[i]->ProjectionY(Form(("hs"+region_label+"Peaksr_xi%d").c_str(),i), PD.sr_low, PD.sr_high);

            TF1* quadFit16 = new TF1("quadFit16","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit16->SetParameters(1,1,1);

            //nEvent_xi.push_back(mult_xi->Integral(1,10000));
            nEvent_xi.push_back(mult_xi->Integral(2,10000));
            Bz_xi.push_back(hbackgroundPeak_xi[i]->GetBinContent(hbackgroundPeak_xi[i]->FindBin(0,0)));

            hsPeaksr_xi[i]->Divide(hbPeaksr_xi[i]);
            hsPeaksr_xi[i]->Scale(Bz_xi[i]/nEvent_xi[i]/BW2D);
            //hsPeaksr_xi[i]->Scale(Bz_xi[i]/BW2D);

            hbPeaklr_xi[i] = hbackgroundPeak_xi[i]->ProjectionY("hbPeaklr_xi",1,10);
            TH1D* ahbPeaklr_xi = hbackgroundPeak_xi[i]->ProjectionY("ahbPeaklr_xi",24,33);
            hsPeaklr_xi[i] = hsignalPeak_xi[i]->ProjectionY("hsPeaklr_xi",1,10);
            TH1D* ahsPeaklr_xi = hsignalPeak_xi[i]->ProjectionY("ahsPeaklr_xi",24,33);

            hbPeaklr_xi[i]->Add(ahbPeaklr_xi);
            hsPeaklr_xi[i]->Add(ahsPeaklr_xi);
            hsPeaklr_xi[i]->Divide(hbPeaklr_xi[i]);
            hsPeaklr_xi[i]->Scale(Bz_xi[i]/nEvent_xi[i]/BW2D);
            //hsPeaklr_xi[i]->Scale(Bz_xi[i]/BW2D);


            hsPeaksr_xi[i]->Add(hsPeaklr_xi[i],-1);

            c_sr_xi[j]->cd(i+1);

            TH1D* hsPeaksr_zeroed_xi = (TH1D*)hsPeaksr_xi[i]->Clone();

            hsPeaksr_xi[i]->Fit("quadFit16","R");
            hsPeaksr_xi[i]->Fit("quadFit16","R");
            hsPeaksr_xi[i]->Fit("quadFit16","R");

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
            tex->DrawLatex(xcoor,ycoor-=increment,"185<N_{trk}^{offline}<250");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%4.1f<p_{T}^{trg}<%4.1f",PD.PtBin_xi[i],PD.PtBin_xi[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_xi.push_back(hsPeaksr_zeroed_xi->IntegralAndError(hsPeaksr_zeroed_xi->FindBin(PD.fw_low),hsPeaksr_zeroed_xi->FindBin(PD.fw_high),Jyieldsr_err_xi[i],"width"));
            //double bin0yield = hsPeaksr_zeroed_xi->GetBinContent(hsPeaksr_zeroed_xi->FindBin(0.0))*0.19635;
            //Jyieldsr_xi[i] = Jyieldsr_xi[i]*2 - bin0yield;

            JyieldSub_xi.push_back(Jyieldsr_xi[i]);
            JyieldSub_err_xi.push_back(Jyieldsr_err_xi[i]);

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_xi[i] = hsignalPeak_xi[i]->ProjectionY(Form("V2lrs_xi%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_xi = hsignalPeak_xi[i]->ProjectionY("aV2lrs_xi",PD.binhigh_V2,33);
            V2lrb_xi[i] = hbackgroundPeak_xi[i]->ProjectionY(Form("V2lrb_xi%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_xi = hbackgroundPeak_xi[i]->ProjectionY("aV2lrb_xi",PD.binhigh_V2,33);
            V2lrs_xi[i]->Add(aV2lrs_xi);
            V2lrb_xi[i]->Add(aV2lrb_xi);
            V2lrs_xi[i]->Divide(V2lrb_xi[i]);
            V2lrs_xi[i]->Scale(Bz_xi[i]/nEvent_xi[i]/BW2D);
            //V2lrs_xi[i]->Scale(Bz_xi[i]/BW2D);

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

        //Om High N
        std::vector<double> Nassoc_om;
        std::vector<double> Jyieldsr_om;
        std::vector<double> Jyieldlr_om;
        std::vector<double> Jyieldsr_err_om(18);
        std::vector<double> Jyieldlr_err_om(18);
        std::vector<double> V2Values_om;
        std::vector<double> V2Values_err_om;
        std::vector<double> V3Values_om;
        std::vector<double> V3Values_err_om;
        std::vector<double> Bz_om;
        std::vector<double> nEvent_om;

        std::vector<double> JyieldSub_om;
        std::vector<double> JyieldSub_err_om;

        TH1D* hsPeaksr_om[arraySize_om];
        TH1D* hbPeaksr_om[arraySize_om];
        TH1D* hsPeaklr_om[arraySize_om];
        TH1D* hbPeaklr_om[arraySize_om];
        TH1D* V2lrs_om[arraySize_om];
        TH1D* V2lrb_om[arraySize_om];
        TH2D* hbackgroundPeak_om[arraySize_om];
        TH2D* hsignalPeak_om[arraySize_om];

        //Calculate Nassoc, Jet yield, Peak region
        //Jet Yield

        for(int i=0; i<numPtBins_om; i++)
        {
            PD.f_Om->GetObject(Form(PathBackground_high_om.c_str(),i),hbackgroundPeak_om[i]);
            PD.f_Om->GetObject(Form(PathSignal_high_om.c_str(),i),hsignalPeak_om[i]);
            TH1D* mult_om = (TH1D*) PD.f_Om->Get(Form(PathMult_high_om.c_str(),i));

            hbPeaksr_om[i] = hbackgroundPeak_om[i]->ProjectionY("hbPeaksr_om", PD.sr_low, PD.sr_high);
            hsPeaksr_om[i] = hsignalPeak_om[i]->ProjectionY(Form(("hs"+region_label+"Peaksr_om%d").c_str(),i), PD.sr_low, PD.sr_high);

            TF1* quadFit16 = new TF1("quadFit16","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit16->SetParameters(1,1,1);

            //nEvent_om.push_back(mult_om->Integral(1,10000));
            nEvent_om.push_back(mult_om->Integral(2,10000));
            Bz_om.push_back(hbackgroundPeak_om[i]->GetBinContent(hbackgroundPeak_om[i]->FindBin(0,0)));

            hsPeaksr_om[i]->Divide(hbPeaksr_om[i]);
            hsPeaksr_om[i]->Scale(Bz_om[i]/nEvent_om[i]/BW2D);
            //hsPeaksr_om[i]->Scale(Bz_om[i]/BW2D);

            hbPeaklr_om[i] = hbackgroundPeak_om[i]->ProjectionY("hbPeaklr_om",1,10);
            TH1D* ahbPeaklr_om = hbackgroundPeak_om[i]->ProjectionY("ahbPeaklr_om",24,33);
            hsPeaklr_om[i] = hsignalPeak_om[i]->ProjectionY("hsPeaklr_om",1,10);
            TH1D* ahsPeaklr_om = hsignalPeak_om[i]->ProjectionY("ahsPeaklr_om",24,33);

            hbPeaklr_om[i]->Add(ahbPeaklr_om);
            hsPeaklr_om[i]->Add(ahsPeaklr_om);
            hsPeaklr_om[i]->Divide(hbPeaklr_om[i]);
            hsPeaklr_om[i]->Scale(Bz_om[i]/nEvent_om[i]/BW2D);
            //hsPeaklr_om[i]->Scale(Bz_om[i]/BW2D);


            hsPeaksr_om[i]->Add(hsPeaklr_om[i],-1);

            c_sr_om[j]->cd(i+1);

            TH1D* hsPeaksr_zeroed_om = (TH1D*)hsPeaksr_om[i]->Clone();

            hsPeaksr_om[i]->Fit("quadFit16","R");
            hsPeaksr_om[i]->Fit("quadFit16","R");
            hsPeaksr_om[i]->Fit("quadFit16","R");

            double minVal_sr = quadFit16->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit16->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            hsPeaksr_zeroed_om->Add(minConst_sr);
            hsPeaksr_om[i]->Draw();
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"185<N_{trk}^{offline}<250");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%4.1f<p_{T}^{trg}<%4.1f",PD.PtBin_om[i],PD.PtBin_om[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_om.push_back(hsPeaksr_zeroed_om->IntegralAndError(hsPeaksr_zeroed_om->FindBin(PD.fw_low),hsPeaksr_zeroed_om->FindBin(PD.fw_high),Jyieldsr_err_om[i],"width"));
            //double bin0yield = hsPeaksr_zeroed_om->GetBinContent(hsPeaksr_zeroed_om->FindBin(0.0))*0.19635;
            //Jyieldsr_om[i] = Jyieldsr_om[i]*2 - bin0yield;

            JyieldSub_om.push_back(Jyieldsr_om[i]);
            JyieldSub_err_om.push_back(Jyieldsr_err_om[i]);

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_om[i] = hsignalPeak_om[i]->ProjectionY(Form("V2lrs_om%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_om = hsignalPeak_om[i]->ProjectionY("aV2lrs_om",PD.binhigh_V2,33);
            V2lrb_om[i] = hbackgroundPeak_om[i]->ProjectionY(Form("V2lrb_om%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_om = hbackgroundPeak_om[i]->ProjectionY("aV2lrb_om",PD.binhigh_V2,33);
            V2lrs_om[i]->Add(aV2lrs_om);
            V2lrb_om[i]->Add(aV2lrb_om);
            V2lrs_om[i]->Divide(V2lrb_om[i]);
            V2lrs_om[i]->Scale(Bz_om[i]/nEvent_om[i]/BW2D);
            //V2lrs_om[i]->Scale(Bz_om[i]/BW2D);

            V2lrs_om[i]->Fit("fit1","R");
            //V2lrs_om[i]->Fit("fit1","R");
            //V2lrs_om[i]->Fit("fit1","R");

            V2Values_om.push_back(fit1->GetParameter(2));
            V2Values_err_om.push_back(fit1->GetParError(2));
            V3Values_om.push_back(fit1->GetParameter(3));
            V3Values_err_om.push_back(fit1->GetParError(3));

            Nassoc_om.push_back(fit1->GetParameter(0));
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

        //Om Vn correction
        std::vector<double> V2sub_om;
        std::vector<double> V2sube_om;

        std::vector<double> V3sub_om;
        std::vector<double> V3sube_om;

        std::vector<double> jetYfactor_om;

        for(unsigned i=0; i<numPtBins_om; i++){
            V2sub_om.push_back(V2Values_om[i] - V2Values_low_om[i]*Nassoc_low_om[i]/Nassoc_om[i]*JyieldSub_om[i]/JyieldSub_low_om[i]);
            V2sube_om.push_back(sqrt(TMath::Power(V2Values_err_om[i],2) + sqrt(TMath::Power(V2Values_err_low_om[i]/V2Values_low_om[i],2) + TMath::Power(JyieldSub_err_om[i]/JyieldSub_om[i],2) + TMath::Power(JyieldSub_err_low_om[i]/JyieldSub_low_om[i],2))*V2Values_low_om[i]*Nassoc_low_om[i]/Nassoc_om[i]*JyieldSub_om[i]/JyieldSub_low_om[i]*sqrt(TMath::Power(V2Values_err_low_om[i]/V2Values_low_om[i],2) + TMath::Power(JyieldSub_err_om[i]/JyieldSub_om[i],2) + TMath::Power(JyieldSub_err_low_om[i]/JyieldSub_low_om[i],2))*V2Values_low_om[i]*Nassoc_low_om[i]/Nassoc_om[i]*JyieldSub_om[i]/JyieldSub_low_om[i]));

            V3sub_om.push_back(V3Values_om[i] - V3Values_low_om[i]*Nassoc_low_om[i]/Nassoc_om[i]*JyieldSub_om[i]/JyieldSub_low_om[i]);
            V3sube_om.push_back(sqrt(TMath::Power(V3Values_err_om[i],2) + sqrt(TMath::Power(V3Values_err_low_om[i]/V3Values_low_om[i],2) + TMath::Power(JyieldSub_err_om[i]/JyieldSub_om[i],2) + TMath::Power(JyieldSub_err_low_om[i]/JyieldSub_low_om[i],2))*V3Values_low_om[i]*Nassoc_low_om[i]/Nassoc_om[i]*JyieldSub_om[i]/JyieldSub_low_om[i]*sqrt(TMath::Power(V3Values_err_low_om[i]/V3Values_low_om[i],2) + TMath::Power(JyieldSub_err_om[i]/JyieldSub_om[i],2) + TMath::Power(JyieldSub_err_low_om[i]/JyieldSub_low_om[i],2))*V3Values_low_om[i]*Nassoc_low_om[i]/Nassoc_om[i]*JyieldSub_om[i]/JyieldSub_low_om[i]));

            jetYfactor_om.push_back(JyieldSub_om[i]/JyieldSub_low_om[i]);
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
            v3sube_ks.push_back(fabs(sqrt(V3sube_ks[i]/V3sub_ks[i]*V3sube_ks[i]/V3sub_ks[i] + v3sube_ref/v3sub_ref*v3sube_ref/v3sub_ref)*v3sub_ks[i]));

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
            v3sube_la.push_back(fabs(sqrt(TMath::Power(V3sube_la[i]/V3sub_la[i],2) + TMath::Power(v3sube_ref/v3sub_ref,2))*v3sub_la[i]));

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
            v3sube_xi.push_back(fabs(sqrt(TMath::Power(V3sube_xi[i]/V3sub_xi[i],2) + TMath::Power(v3sube_ref/v3sub_ref,2))*v3sub_xi[i]));

            v3_xi.push_back(V3Values_xi[i]/v3_high_ref);
            v3e_xi.push_back(sqrt(TMath::Power(V3Values_err_xi[i]/V3Values_xi[i],2) + TMath::Power(v3e_high_ref/v3_high_ref,2))*v3_xi[i]);

            v3_low_xi.push_back(V3Values_low_xi[i]/v3_low_ref);
            v3e_low_xi.push_back(sqrt(TMath::Power(V3Values_err_low_xi[i]/V3Values_low_xi[i],2) + TMath::Power(v3e_low_ref/v3_low_ref,2))*v3_low_xi[i]);
        }

        //Om vn calculation
        std::vector<double> v2sub_om;
        std::vector<double> v2sube_om;
        std::vector<double> v2_om;
        std::vector<double> v2e_om;
        std::vector<double> v2_low_om;
        std::vector<double> v2e_low_om;

        std::vector<double> v3sub_om;
        std::vector<double> v3sube_om;
        std::vector<double> v3_om;
        std::vector<double> v3e_om;
        std::vector<double> v3_low_om;
        std::vector<double> v3e_low_om;

        for(unsigned i=0; i<numPtBins_om; i++){
            v2sub_om.push_back(V2sub_om[i]/v2sub_ref);
            v2sube_om.push_back(fabs(sqrt(TMath::Power(V2sube_om[i]/V2sub_om[i],2) + TMath::Power(v2sube_ref/v2sub_ref,2)))*v2sub_om[i]);

            v2_om.push_back(V2Values_om[i]/v2_high_ref);
            v2e_om.push_back(sqrt(TMath::Power(V2Values_err_om[i]/V2Values_om[i],2) + TMath::Power(v2e_high_ref/v2_high_ref,2))*v2_om[i]);

            v2_low_om.push_back(V2Values_low_om[i]/v2_low_ref);
            v2e_low_om.push_back(sqrt(TMath::Power(V2Values_err_low_om[i]/V2Values_low_om[i],2) + TMath::Power(v2e_low_ref/v2_low_ref,2))*v2_low_om[i]);

            v3sub_om.push_back(V3sub_om[i]/v3sub_ref);
            v3sube_om.push_back(fabs(sqrt(TMath::Power(V3sube_om[i]/V3sub_om[i],2) + TMath::Power(v3sube_ref/v3sub_ref,2))*v3sub_om[i]));

            v3_om.push_back(V3Values_om[i]/v3_high_ref);
            v3e_om.push_back(sqrt(TMath::Power(V3Values_err_om[i]/V3Values_om[i],2) + TMath::Power(v3e_high_ref/v3_high_ref,2))*v3_om[i]);

            v3_low_om.push_back(V3Values_low_om[i]/v3_low_ref);
            v3e_low_om.push_back(sqrt(TMath::Power(V3Values_err_low_om[i]/V3Values_low_om[i],2) + TMath::Power(v3e_low_ref/v3_low_ref,2))*v3_low_om[i]);
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

        std::vector<double> pt_om;
        std::vector<double> Ket_om;
        for(unsigned i=0; i<numPtBins_om; i++){
            TH1D* hpt = (TH1D*)PD.f_Om->Get(Form((PD.fn_Om + "/Pt_om_pt%d").c_str(),i));
            TH1D* hket = (TH1D*)PD.f_Om->Get(Form((PD.fn_Om + "/KET_om_pt%d").c_str(),i));

            pt_om.push_back(hpt->GetMean(1));
            Ket_om.push_back(hket->GetMean(1));
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

        std::vector<double> perisubfactor_om;
        std::vector<double> perisubfactore_om;
        for(int i=0; i<numPtBins_om; i++){
            perisubfactor_om.push_back(V2Values_low_om[i]*Nassoc_low_om[i]/JyieldSub_low_om[i]);
            perisubfactore_om.push_back(sqrt(TMath::Power(V2Values_err_low_om[i]/V2Values_low_om[i],2) + TMath::Power(JyieldSub_err_low_om[i]/JyieldSub_low_om[i],2))*perisubfactor_om[i]);
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
        TGraphErrors* v2plot_low_xi    = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&v2_low_xi[0]      ,0,&v2e_low_xi[0]);
        TGraphErrors* v2plot_KET_xi    = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&v2_xi[0]      ,0,&v2e_xi[0]);
        TGraphErrors* v2subplot_xi     = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&v2sub_xi[0]   ,0,&v2sube_xi[0]);
        TGraphErrors* v2subplot_KET_xi = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&v2sub_xi[0]   ,0,&v2sube_xi[0]);
        TGraphErrors* v3plot_xi        = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&v3_xi[0]      ,0,&v3e_xi[0]);
        TGraphErrors* v3plot_KET_xi    = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&v3_xi[0]      ,0,&v3e_xi[0]);
        TGraphErrors* v3subplot_xi     = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&v3sub_xi[0]   ,0,&v3sube_xi[0]);
        TGraphErrors* v3subplot_KET_xi = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&v3sub_xi[0]   ,0,&v3sube_xi[0]);

        TGraphErrors* v2plot_om        = new TGraphErrors(numPtBins_om,&pt_om[0] ,&v2_om[0]      ,0,&v2e_om[0]);
        TGraphErrors* v2plot_low_om    = new TGraphErrors(numPtBins_om,&pt_om[0] ,&v2_low_om[0]      ,0,&v2e_low_om[0]);
        TGraphErrors* v2plot_KET_om    = new TGraphErrors(numPtBins_om,&Ket_om[0],&v2_om[0]      ,0,&v2e_om[0]);
        TGraphErrors* v2subplot_om     = new TGraphErrors(numPtBins_om,&pt_om[0] ,&v2sub_om[0]   ,0,&v2sube_om[0]);
        TGraphErrors* v2subplot_KET_om = new TGraphErrors(numPtBins_om,&Ket_om[0],&v2sub_om[0]   ,0,&v2sube_om[0]);
        TGraphErrors* v3plot_om        = new TGraphErrors(numPtBins_om,&pt_om[0] ,&v3_om[0]      ,0,&v3e_om[0]);
        TGraphErrors* v3plot_KET_om    = new TGraphErrors(numPtBins_om,&Ket_om[0],&v3_om[0]      ,0,&v3e_om[0]);
        TGraphErrors* v3subplot_om     = new TGraphErrors(numPtBins_om,&pt_om[0] ,&v3sub_om[0]   ,0,&v3sube_om[0]);
        TGraphErrors* v3subplot_KET_om = new TGraphErrors(numPtBins_om,&Ket_om[0],&v3sub_om[0]   ,0,&v3sube_om[0]);

        //Obs V2 values
        TGraphErrors* V2plot_ks        = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&V2Values_ks[0]    ,0,&V2Values_err_ks[0]);
        TGraphErrors* V2plot_low_ks    = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&V2Values_low_ks[0],0,&V2Values_err_low_ks[0]);
        TGraphErrors* V2plot_KET_ks    = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&V2Values_ks[0]    ,0,&V2Values_err_ks[0]);
        TGraphErrors* V2subplot_ks     = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&V2sub_ks[0]       ,0,&V2sube_ks[0]);
        TGraphErrors* V2subplot_KET_ks = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&V2sub_ks[0]       ,0,&V2sube_ks[0]);
        TGraphErrors* V3plot_ks        = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&V3Values_ks[0]    ,0,&V3Values_err_ks[0]);
        TGraphErrors* V3plot_KET_ks    = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&V3Values_ks[0]    ,0,&V3Values_err_ks[0]);
        TGraphErrors* V3subplot_ks     = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&V3sub_ks[0]       ,0,&V3sube_ks[0]);
        TGraphErrors* V3subplot_KET_ks = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&V3sub_ks[0]       ,0,&V3sube_ks[0]);

        TGraphErrors* V2plot_la        = new TGraphErrors(numPtBins_la,&pt_la[0] ,&V2Values_la[0]    ,0,&V2Values_err_la[0]);
        TGraphErrors* V2plot_low_la    = new TGraphErrors(numPtBins_la,&pt_la[0] ,&V2Values_low_la[0],0,&V2Values_err_low_la[0]);
        TGraphErrors* V2plot_KET_la    = new TGraphErrors(numPtBins_la,&Ket_la[0],&V2Values_la[0]    ,0,&V2Values_err_la[0]);
        TGraphErrors* V2subplot_la     = new TGraphErrors(numPtBins_la,&pt_la[0] ,&V2sub_la[0]       ,0,&V2sube_la[0]);
        TGraphErrors* V2subplot_KET_la = new TGraphErrors(numPtBins_la,&Ket_la[0],&V2sub_la[0]       ,0,&V2sube_la[0]);
        TGraphErrors* V3plot_la        = new TGraphErrors(numPtBins_la,&pt_la[0] ,&V3Values_la[0]    ,0,&V3Values_err_la[0]);
        TGraphErrors* V3plot_KET_la    = new TGraphErrors(numPtBins_la,&Ket_la[0],&V3Values_la[0]    ,0,&V3Values_err_la[0]);
        TGraphErrors* V3subplot_la     = new TGraphErrors(numPtBins_la,&pt_la[0] ,&V3sub_la[0]       ,0,&V3sube_la[0]);
        TGraphErrors* V3subplot_KET_la = new TGraphErrors(numPtBins_la,&Ket_la[0],&V3sub_la[0]       ,0,&V3sube_la[0]);

        TGraphErrors* V2plot_xi        = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&V2Values_xi[0]    ,0,&V2Values_err_xi[0]);
        TGraphErrors* V2plot_low_xi    = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&V2Values_low_xi[0],0,&V2Values_err_low_xi[0]);
        TGraphErrors* V2plot_KET_xi    = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&V2Values_xi[0]    ,0,&V2Values_err_xi[0]);
        TGraphErrors* V2subplot_xi     = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&V2sub_xi[0]       ,0,&V2sube_xi[0]);
        TGraphErrors* V2subplot_KET_xi = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&V2sub_xi[0]       ,0,&V2sube_xi[0]);
        TGraphErrors* V3plot_xi        = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&V3Values_xi[0]    ,0,&V3Values_err_xi[0]);
        TGraphErrors* V3plot_KET_xi    = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&V3Values_xi[0]    ,0,&V3Values_err_xi[0]);
        TGraphErrors* V3subplot_xi     = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&V3sub_xi[0]       ,0,&V3sube_xi[0]);
        TGraphErrors* V3subplot_KET_xi = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&V3sub_xi[0]       ,0,&V3sube_xi[0]);

        TGraphErrors* V2plot_om        = new TGraphErrors(numPtBins_om,&pt_om[0] ,&V2Values_om[0]    ,0,&V2Values_err_om[0]);
        TGraphErrors* V2plot_low_om    = new TGraphErrors(numPtBins_om,&pt_om[0] ,&V2Values_low_om[0],0,&V2Values_err_low_om[0]);
        TGraphErrors* V2plot_KET_om    = new TGraphErrors(numPtBins_om,&Ket_om[0],&V2Values_om[0]    ,0,&V2Values_err_om[0]);
        TGraphErrors* V2subplot_om     = new TGraphErrors(numPtBins_om,&pt_om[0] ,&V2sub_om[0]       ,0,&V2sube_om[0]);
        TGraphErrors* V2subplot_KET_om = new TGraphErrors(numPtBins_om,&Ket_om[0],&V2sub_om[0]       ,0,&V2sube_om[0]);
        TGraphErrors* V3plot_om        = new TGraphErrors(numPtBins_om,&pt_om[0] ,&V3Values_om[0]    ,0,&V3Values_err_om[0]);
        TGraphErrors* V3plot_KET_om    = new TGraphErrors(numPtBins_om,&Ket_om[0],&V3Values_om[0]    ,0,&V3Values_err_om[0]);
        TGraphErrors* V3subplot_om     = new TGraphErrors(numPtBins_om,&pt_om[0] ,&V3sub_om[0]       ,0,&V3sube_om[0]);
        TGraphErrors* V3subplot_KET_om = new TGraphErrors(numPtBins_om,&Ket_om[0],&V3sub_om[0]       ,0,&V3sube_om[0]);

        //Yields
        TGraphErrors* srYieldPlot_ks       = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Jyieldsr_ks[0]     ,0,&Jyieldsr_err_ks[0]);
        //TGraphErrors* lrYieldPlot_ks       = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Jyieldlr_ks[0]     ,0,&Jyieldlr_err_ks[0]);
        TGraphErrors* srYieldPlot_low_ks       = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Jyieldsr_low_ks[0]     ,0,&Jyieldsr_err_low_ks[0]);
        //TGraphErrors* lrYieldPlot_low_ks       = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Jyieldlr_low_ks[0]     ,0,&Jyieldlr_err_low_ks[0]);
        TGraphErrors* subYieldPlot_ks      = new TGraphErrors(numPtBins_ks,&pt_ks[0],&JyieldSub_ks[0]    ,0,&JyieldSub_err_ks[0]);
        TGraphErrors* perisubfactorplot_ks = new TGraphErrors(numPtBins_ks,&pt_ks[0],&perisubfactor_ks[0],0,&perisubfactore_ks[0]);

        TGraphErrors* srYieldPlot_la       = new TGraphErrors(numPtBins_la,&pt_la[0],&Jyieldsr_la[0]     ,0,&Jyieldsr_err_la[0]);
        //TGraphErrors* lrYieldPlot_la       = new TGraphErrors(numPtBins_la,&pt_la[0],&Jyieldlr_la[0]     ,0,&Jyieldlr_err_la[0]);
        TGraphErrors* srYieldPlot_low_la       = new TGraphErrors(numPtBins_la,&pt_la[0],&Jyieldsr_low_la[0]     ,0,&Jyieldsr_err_low_la[0]);
        //TGraphErrors* lrYieldPlot_low_la       = new TGraphErrors(numPtBins_la,&pt_la[0],&Jyieldlr_low_la[0]     ,0,&Jyieldlr_err_low_la[0]);
        TGraphErrors* subYieldPlot_la      = new TGraphErrors(numPtBins_la,&pt_la[0],&JyieldSub_la[0]    ,0,&JyieldSub_err_la[0]);
        TGraphErrors* perisubfactorplot_la = new TGraphErrors(numPtBins_la,&pt_la[0],&perisubfactor_la[0],0,&perisubfactore_la[0]);

        TGraphErrors* srYieldPlot_xi       = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Jyieldsr_xi[0]     ,0,&Jyieldsr_err_xi[0]);
        //TGraphErrors* lrYieldPlot_xi       = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Jyieldlr_xi[0]     ,0,&Jyieldlr_err_xi[0]);
        TGraphErrors* srYieldPlot_low_xi       = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Jyieldsr_low_xi[0]     ,0,&Jyieldsr_err_low_xi[0]);
        //TGraphErrors* lrYieldPlot_low_xi       = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Jyieldlr_low_xi[0]     ,0,&Jyieldlr_err_low_xi[0]);
        TGraphErrors* subYieldPlot_xi      = new TGraphErrors(numPtBins_xi,&pt_xi[0],&JyieldSub_xi[0]    ,0,&JyieldSub_err_xi[0]);
        TGraphErrors* perisubfactorplot_xi = new TGraphErrors(numPtBins_xi,&pt_xi[0],&perisubfactor_xi[0],0,&perisubfactore_xi[0]);

        TGraphErrors* srYieldPlot_om       = new TGraphErrors(numPtBins_om,&pt_om[0],&Jyieldsr_om[0]     ,0,&Jyieldsr_err_om[0]);
        //TGraphErrors* lrYieldPlot_om       = new TGraphErrors(numPtBins_om,&pt_om[0],&Jyieldlr_om[0]     ,0,&Jyieldlr_err_om[0]);
        TGraphErrors* srYieldPlot_low_om       = new TGraphErrors(numPtBins_om,&pt_om[0],&Jyieldsr_low_om[0]     ,0,&Jyieldsr_err_low_om[0]);
        //TGraphErrors* lrYieldPlot_low_om       = new TGraphErrors(numPtBins_om,&pt_om[0],&Jyieldlr_low_om[0]     ,0,&Jyieldlr_err_low_om[0]);
        TGraphErrors* subYieldPlot_om      = new TGraphErrors(numPtBins_om,&pt_om[0],&JyieldSub_om[0]    ,0,&JyieldSub_err_om[0]);
        TGraphErrors* perisubfactorplot_om = new TGraphErrors(numPtBins_om,&pt_om[0],&perisubfactor_om[0],0,&perisubfactore_om[0]);

        //For Direct subtraction method
        TGraphErrors* Nasslow_ks         = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Nassoc_low_ks[0]   ,0,0);
        TGraphErrors* Nass_ks            = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Nassoc_ks[0]       ,0,0);
        TGraphErrors* subYieldPlotLow_ks = new TGraphErrors(numPtBins_ks,&pt_ks[0],&JyieldSub_low_ks[0],0,&JyieldSub_err_low_ks[0]);

        TGraphErrors* Nasslow_la         = new TGraphErrors(numPtBins_la,&pt_la[0],&Nassoc_low_la[0]   ,0,0);
        TGraphErrors* Nass_la            = new TGraphErrors(numPtBins_la,&pt_la[0],&Nassoc_la[0]       ,0,0);
        TGraphErrors* subYieldPlotLow_la = new TGraphErrors(numPtBins_la,&pt_la[0],&JyieldSub_low_la[0],0,&JyieldSub_err_low_la[0]);

        TGraphErrors* Nasslow_xi         = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Nassoc_low_xi[0]   ,0,0);
        TGraphErrors* Nass_xi            = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Nassoc_xi[0]       ,0,0);
        TGraphErrors* subYieldPlotLow_xi = new TGraphErrors(numPtBins_xi,&pt_xi[0],&JyieldSub_low_xi[0],0,&JyieldSub_err_low_xi[0]);

        TGraphErrors* Nasslow_om         = new TGraphErrors(numPtBins_om,&pt_om[0],&Nassoc_low_om[0]   ,0,0);
        TGraphErrors* Nass_om            = new TGraphErrors(numPtBins_om,&pt_om[0],&Nassoc_om[0]       ,0,0);
        TGraphErrors* subYieldPlotLow_om = new TGraphErrors(numPtBins_om,&pt_om[0],&JyieldSub_low_om[0],0,&JyieldSub_err_low_om[0]);

        TGraphErrors* bz_ks = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Bz_ks[0],0,0);
        TGraphErrors* bz_la = new TGraphErrors(numPtBins_la,&pt_la[0],&Bz_la[0],0,0);
        TGraphErrors* bz_xi = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Bz_xi[0],0,0);
        TGraphErrors* bz_om = new TGraphErrors(numPtBins_om,&pt_om[0],&Bz_om[0],0,0);

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
        v2plot_low_xi       ->Write(("v2plot_" + region + "_low_xi"       ).c_str(),TObject::kOverwrite);
        v2plot_KET_xi   ->Write(("v2plot_KET_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        v2subplot_xi    ->Write(("v2subplot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);
        v2subplot_KET_xi->Write(("v2subplot_KET_" + region + "_xi").c_str(),TObject::kOverwrite);
        v3plot_xi       ->Write(("v3plot_" + region + "_xi"       ).c_str(),TObject::kOverwrite);
        v3plot_KET_xi   ->Write(("v3plot_KET_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        v3subplot_xi    ->Write(("v3subplot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);
        v3subplot_KET_xi->Write(("v3subplot_KET_" + region + "_xi").c_str(),TObject::kOverwrite);

        v2plot_om       ->Write(("v2plot_" + region + "_om"       ).c_str(),TObject::kOverwrite);
        v2plot_low_om       ->Write(("v2plot_" + region + "_low_om"       ).c_str(),TObject::kOverwrite);
        v2plot_KET_om   ->Write(("v2plot_KET_" + region + "_om"   ).c_str(),TObject::kOverwrite);
        v2subplot_om    ->Write(("v2subplot_" + region + "_om"    ).c_str(),TObject::kOverwrite);
        v2subplot_KET_om->Write(("v2subplot_KET_" + region + "_om").c_str(),TObject::kOverwrite);
        v3plot_om       ->Write(("v3plot_" + region + "_om"       ).c_str(),TObject::kOverwrite);
        v3plot_KET_om   ->Write(("v3plot_KET_" + region + "_om"   ).c_str(),TObject::kOverwrite);
        v3subplot_om    ->Write(("v3subplot_" + region + "_om"    ).c_str(),TObject::kOverwrite);
        v3subplot_KET_om->Write(("v3subplot_KET_" + region + "_om").c_str(),TObject::kOverwrite);

        V2plot_ks        -> Write(("V2plot_" + region + "_ks"       ).c_str(),TObject::kOverwrite);
        V2plot_low_ks    -> Write(("V2plot_" + region + "_low_ks"   ).c_str(),TObject::kOverwrite);
        V2plot_KET_ks    -> Write(("V2plot_KET_" + region + "_ks"   ).c_str(),TObject::kOverwrite);
        V2subplot_ks     -> Write(("V2subplot_" + region + "_ks"    ).c_str(),TObject::kOverwrite);
        V2subplot_KET_ks -> Write(("V2subplot_KET_" + region + "_ks").c_str(),TObject::kOverwrite);
        V3plot_ks        -> Write(("V3plot_" + region + "_ks"       ).c_str(),TObject::kOverwrite);
        V3plot_KET_ks    -> Write(("V3plot_KET_" + region + "_ks"   ).c_str(),TObject::kOverwrite);
        V3subplot_ks     -> Write(("V3subplot_" + region + "_ks"    ).c_str(),TObject::kOverwrite);
        V3subplot_KET_ks -> Write(("V3subplot_KET_" + region + "_ks").c_str(),TObject::kOverwrite);

        V2plot_la        -> Write(("V2plot_" + region + "_la"       ).c_str(),TObject::kOverwrite);
        V2plot_low_la    -> Write(("V2plot_" + region + "_low_la"   ).c_str(),TObject::kOverwrite);
        V2plot_KET_la    -> Write(("V2plot_KET_" + region + "_la"   ).c_str(),TObject::kOverwrite);
        V2subplot_la     -> Write(("V2subplot_" + region + "_la"    ).c_str(),TObject::kOverwrite);
        V2subplot_KET_la -> Write(("V2subplot_KET_" + region + "_la").c_str(),TObject::kOverwrite);
        V3plot_la        -> Write(("V3plot_" + region + "_la"       ).c_str(),TObject::kOverwrite);
        V3plot_KET_la    -> Write(("V3plot_KET_" + region + "_la"   ).c_str(),TObject::kOverwrite);
        V3subplot_la     -> Write(("V3subplot_" + region + "_la"    ).c_str(),TObject::kOverwrite);
        V3subplot_KET_la -> Write(("V3subplot_KET_" + region + "_la").c_str(),TObject::kOverwrite);

        V2plot_xi        -> Write(("V2plot_" + region + "_xi"       ).c_str(),TObject::kOverwrite);
        V2plot_low_xi    -> Write(("V2plot_" + region + "_low_xi"   ).c_str(),TObject::kOverwrite);
        V2plot_KET_xi    -> Write(("V2plot_KET_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        V2subplot_xi     -> Write(("V2subplot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);
        V2subplot_KET_xi -> Write(("V2subplot_KET_" + region + "_xi").c_str(),TObject::kOverwrite);
        V3plot_xi        -> Write(("V3plot_" + region + "_xi"       ).c_str(),TObject::kOverwrite);
        V3plot_KET_xi    -> Write(("V3plot_KET_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        V3subplot_xi     -> Write(("V3subplot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);
        V3subplot_KET_xi -> Write(("V3subplot_KET_" + region + "_xi").c_str(),TObject::kOverwrite);

        V2plot_om        -> Write(("V2plot_" + region + "_om"       ).c_str(),TObject::kOverwrite);
        V2plot_low_om    -> Write(("V2plot_" + region + "_low_om"   ).c_str(),TObject::kOverwrite);
        V2plot_KET_om    -> Write(("V2plot_KET_" + region + "_om"   ).c_str(),TObject::kOverwrite);
        V2subplot_om     -> Write(("V2subplot_" + region + "_om"    ).c_str(),TObject::kOverwrite);
        V2subplot_KET_om -> Write(("V2subplot_KET_" + region + "_om").c_str(),TObject::kOverwrite);
        V3plot_om        -> Write(("V3plot_" + region + "_om"       ).c_str(),TObject::kOverwrite);
        V3plot_KET_om    -> Write(("V3plot_KET_" + region + "_om"   ).c_str(),TObject::kOverwrite);
        V3subplot_om     -> Write(("V3subplot_" + region + "_om"    ).c_str(),TObject::kOverwrite);
        V3subplot_KET_om -> Write(("V3subplot_KET_" + region + "_om").c_str(),TObject::kOverwrite);

        perisubfactorplot_ks->Write(("perisubfactor_" + region + "_ks").c_str(),TObject::kOverwrite);
        perisubfactorplot_la->Write(("perisubfactor_" + region + "_la").c_str(),TObject::kOverwrite);
        perisubfactorplot_xi->Write(("perisubfactor_" + region + "_xi").c_str(),TObject::kOverwrite);
        perisubfactorplot_om->Write(("perisubfactor_" + region + "_om").c_str(),TObject::kOverwrite);

        Nass_ks           ->Write(("Nassoc_" + region + "_ks"       ).c_str(),TObject::kOverwrite);
        Nasslow_ks        ->Write(("Nassoc_low_" + region + "_ks"   ).c_str(),TObject::kOverwrite);
        subYieldPlotLow_ks->Write(("YieldPlot_low_" + region + "_ks").c_str(),TObject::kOverwrite);
        subYieldPlot_ks   ->Write(("YieldPlot_" + region + "_ks"    ).c_str(),TObject::kOverwrite);
        srYieldPlot_ks->Write(("YieldsrPlot_" + region + "_ks").c_str(),TObject::kOverwrite);
        //lrYieldPlot_ks->Write(("YieldlrPlot_" + region + "_ks").c_str(),TObject::kOverwrite);
        srYieldPlot_low_ks->Write(("YieldsrPlot_" + region + "low_ks").c_str(),TObject::kOverwrite);
        //lrYieldPlot_low_ks->Write(("YieldlrPlot_" + region + "low_ks").c_str(),TObject::kOverwrite);

        Nass_la           ->Write(("Nassoc_" + region + "_la"       ).c_str(),TObject::kOverwrite);
        Nasslow_la        ->Write(("Nassoc_low_" + region + "_la"   ).c_str(),TObject::kOverwrite);
        subYieldPlotLow_la->Write(("YieldPlot_low_" + region + "_la").c_str(),TObject::kOverwrite);
        subYieldPlot_la   ->Write(("YieldPlot_" + region + "_la"    ).c_str(),TObject::kOverwrite);
        srYieldPlot_la->Write(("YieldsrPlot_" + region + "_la").c_str(),TObject::kOverwrite);
        //lrYieldPlot_la->Write(("YieldlrPlot_" + region + "_la").c_str(),TObject::kOverwrite);
        srYieldPlot_low_la->Write(("YieldsrPlot_" + region + "low_la").c_str(),TObject::kOverwrite);
        //lrYieldPlot_low_la->Write(("YieldlrPlot_" + region + "low_la").c_str(),TObject::kOverwrite);

        Nass_xi           ->Write(("Nassoc_" + region + "_xi"       ).c_str(),TObject::kOverwrite);
        Nasslow_xi        ->Write(("Nassoc_low_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        subYieldPlotLow_xi->Write(("YieldPlot_low_" + region + "_xi").c_str(),TObject::kOverwrite);
        subYieldPlot_xi   ->Write(("YieldPlot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);
        srYieldPlot_xi->Write(("YieldsrPlot_" + region + "_xi").c_str(),TObject::kOverwrite);
        //lrYieldPlot_xi->Write(("YieldlrPlot_" + region + "_xi").c_str(),TObject::kOverwrite);
        srYieldPlot_low_xi->Write(("YieldsrPlot_" + region + "low_xi").c_str(),TObject::kOverwrite);
        //lrYieldPlot_low_xi->Write(("YieldlrPlot_" + region + "low_xi").c_str(),TObject::kOverwrite);

        Nass_om           ->Write(("Nassoc_" + region + "_om"       ).c_str(),TObject::kOverwrite);
        Nasslow_om        ->Write(("Nassoc_low_" + region + "_om"   ).c_str(),TObject::kOverwrite);
        subYieldPlotLow_om->Write(("YieldPlot_low_" + region + "_om").c_str(),TObject::kOverwrite);
        subYieldPlot_om   ->Write(("YieldPlot_" + region + "_om"    ).c_str(),TObject::kOverwrite);
        srYieldPlot_om->Write(("YieldsrPlot_" + region + "_om").c_str(),TObject::kOverwrite);
        //lrYieldPlot_om->Write(("YieldlrPlot_" + region + "_om").c_str(),TObject::kOverwrite);
        srYieldPlot_low_om->Write(("YieldsrPlot_" + region + "low_om").c_str(),TObject::kOverwrite);
        //lrYieldPlot_low_om->Write(("YieldlrPlot_" + region + "low_om").c_str(),TObject::kOverwrite);

        bz_ks->Write(("bz_" + region + "_ks").c_str(),TObject::kOverwrite);
        bz_la->Write(("bz_" + region + "_la").c_str(),TObject::kOverwrite);
        bz_xi->Write(("bz_" + region + "_xi").c_str(),TObject::kOverwrite);
        bz_om->Write(("bz_" + region + "_om").c_str(),TObject::kOverwrite);

        //c_lr_low_Fourier_xi[j]->Print(("Image/PeriSub/Fourier_" + region + "_low_xiEG2_2terms.pdf").c_str());
        c_sr_low_xi[j]->Print(("Image/PeriSub/Yield_sr_" + region + "_low_xi.pdf").c_str());
        //c_lr_low_xi[j]->Print(("Image/PeriSub/Yield_lr_" + region + "_low_xi.pdf").c_str());

        //c_lr_low_Fourier_om[j]->Print(("Image/PeriSub/Fourier_" + region + "_low_omEG2_2terms.pdf").c_str());
        c_sr_low_om[j]->Print(("Image/PeriSub/Yield_sr_" + region + "_low_om.pdf").c_str());
        //c_lr_low_om[j]->Print(("Image/PeriSub/Yield_lr_" + region + "_low_om.pdf").c_str());

    }
}

void PeriSubObs(ParticleData PD)
{
    TH1::SetDefaultSumw2();

    double BW2D = (9.9/33)*((2-1.0/16)*3.1416/31);

    TLatex* tex = new TLatex();
    tex->SetNDC();


    //TCanvas* csPeaksr_low_ref = new TCanvas("csPeaksr_low_ref","csPeaksr_low_ref",600,600);
    //TCanvas* csPeaklr_low_ref = new TCanvas("csPeaklr_low_ref","csPeaklr_low_ref",600,600);
    //TCanvas* csPeaksr_high_ref = new TCanvas("csPeaksr_high_ref","csPeaksr_high_ref",600,600);
    //TCanvas* csPeaklr_high_ref = new TCanvas("csPeaklr_high_ref","csPeaklr_high_ref",600,600);
    //TCanvas* csPeaksr_low_ks = new TCanvas("csPeaksr_low_ks","csPeaksr_low_ks",1200,900);
    //TCanvas* csPeaklr_low_ks = new TCanvas("csPeaklr_low_ks","csPeaklr_low_ks",1200,900);
    //TCanvas* csPeaksr_ks = new TCanvas("csPeaksr_ks","csPeaksr_ks",1200,900);
    //TCanvas* csPeaklr_ks = new TCanvas("csPeaklr_ks","csPeaklr_ks",1200,900);
    //TCanvas* csPeaksr_low_la = new TCanvas("csPeaksr_low_la","csPeaksr_low_la",1200,900);
    //TCanvas* csPeaklr_low_la = new TCanvas("csPeaklr_low_la","csPeaklr_low_la",1200,900);
    //TCanvas* csPeaksr_la = new TCanvas("csPeaksr_la","csPeaksr_la",1200,900);
    //TCanvas* csPeaklr_la = new TCanvas("csPeaklr_la","csPeaklr_la",1200,900);
    //TCanvas* csPeaksr_low_xi = new TCanvas("csPeaksr_low_xi","csPeaksr_low_xi",1200,900);
    //TCanvas* csPeaklr_low_xi = new TCanvas("csPeaklr_low_xi","csPeaklr_low_xi",1200,900);
    //TCanvas* csPeaksr_xi = new TCanvas("csPeaksr_xi","csPeaksr_xi",1200,900);
    //TCanvas* csPeaklr_xi = new TCanvas("csPeaklr_xi","csPeaklr_xi",1200,900);

    TCanvas* c_sr_low_ref;
    TCanvas* c_lr_low_ref;
    TCanvas* c_sr_high_ref;
    TCanvas* c_lr_high_ref;
    TCanvas* c_sr_low_ks[2];
    TCanvas* c_lr_low_ks[2];
    TCanvas* c_sr_ks[2];
    TCanvas* c_lr_ks[2];
    TCanvas* c_sr_low_la[2];
    TCanvas* c_lr_low_la[2];
    TCanvas* c_sr_la[2];
    TCanvas* c_lr_la[2];
    TCanvas* c_sr_low_xi[2];
    TCanvas* c_lr_low_xi[2];
    TCanvas* c_sr_xi[2];
    TCanvas* c_lr_xi[2];
    TCanvas* c_sr_low_om[2];
    TCanvas* c_lr_low_om[2];
    TCanvas* c_sr_om[2];
    TCanvas* c_lr_om[2];

    TCanvas* c_lr_low_Fourier_xi[2];
    TCanvas* c_lr_high_Fourier_xi[2];

    TCanvas* c_low_2Dks_1[2];
    TCanvas* c_low_2Dks_2[2];

    TCanvas* c_high_2Dks_1[2];
    TCanvas* c_high_2Dks_2[2];

    TCanvas* c_low_2Dla_1[2];
    TCanvas* c_low_2Dla_2[2];

    TCanvas* c_high_2Dla_1[2];
    TCanvas* c_high_2Dla_2[2];

    TCanvas* c_low_2Dxi_1[2];
    TCanvas* c_low_2Dxi_2[2];

    TCanvas* c_high_2Dxi_1[2];
    TCanvas* c_high_2Dxi_2[2];

    TCanvas* c_low_2Dom_1[2];
    TCanvas* c_low_2Dom_2[2];

    TCanvas* c_high_2Dom_1[2];
    TCanvas* c_high_2Dom_2[2];

    c_sr_low_ref  = new TCanvas("csPeaksr_low_ref" , "cssr_low_ref"  , 600  , 600);
    c_lr_low_ref  = new TCanvas("csPeaklr_low_ref" , "cslr_low_ref"  , 600  , 600);
    c_sr_high_ref = new TCanvas("csPeaksr_high_ref", "cssr_high_ref" , 600  , 600);
    c_lr_high_ref = new TCanvas("csPeaklr_high_ref", "cslr_high_ref" , 600  , 600);
    for(int j=0; j<2; j++)
    {
        std::string region;
        if( j==0 ) region = "Peak";
        if( j==1 ) region = "Side";
        c_sr_low_ks[j]   = new TCanvas(("cs" + region + "sr_low_ks"   ).c_str(), ("cs" + region + "sr_low_ks"   ).c_str(), 1200 , 900);
        c_lr_low_ks[j]   = new TCanvas(("cs" + region + "lr_low_ks"   ).c_str(), ("cs" + region + "lr_low_ks"   ).c_str(), 1200 , 900);
        c_sr_ks[j]       = new TCanvas(("cs" + region + "sr_ks"       ).c_str(), ("cs" + region + "sr_ks"       ).c_str(), 1200 , 900);
        c_lr_ks[j]       = new TCanvas(("cs" + region + "lr_ks"       ).c_str(), ("cs" + region + "lr_ks"       ).c_str(), 1200 , 900);
        c_sr_low_la[j]   = new TCanvas(("cs" + region + "sr_low_la"   ).c_str(), ("cs" + region + "sr_low_la"   ).c_str(), 1200 , 900);
        c_lr_low_la[j]   = new TCanvas(("cs" + region + "lr_low_la"   ).c_str(), ("cs" + region + "lr_low_la"   ).c_str(), 1200 , 900);
        c_sr_la[j]       = new TCanvas(("cs" + region + "sr_la"       ).c_str(), ("cs" + region + "sr_la"       ).c_str(), 1200 , 900);
        c_lr_la[j]       = new TCanvas(("cs" + region + "lr_la"       ).c_str(), ("cs" + region + "lr_la"       ).c_str(), 1200 , 900);
        c_sr_low_xi[j]   = new TCanvas(("cs" + region + "sr_low_xi"   ).c_str(), ("cs" + region + "sr_low_xi"   ).c_str(), 1200 , 900);
        c_lr_low_xi[j]   = new TCanvas(("cs" + region + "lr_low_xi"   ).c_str(), ("cs" + region + "lr_low_xi"   ).c_str(), 1200 , 900);
        c_sr_xi[j]       = new TCanvas(("cs" + region + "sr_xi"       ).c_str(), ("cs" + region + "sr_xi"       ).c_str(), 1200 , 900);
        c_lr_xi[j]       = new TCanvas(("cs" + region + "lr_xi"       ).c_str(), ("cs" + region + "lr_xi"       ).c_str(), 1200 , 900);
        c_sr_low_om[j]   = new TCanvas(("cs" + region + "sr_low_om"   ).c_str(), ("cs" + region + "sr_low_om"   ).c_str(), 1200 , 900);
        c_lr_low_om[j]   = new TCanvas(("cs" + region + "lr_low_om"   ).c_str(), ("cs" + region + "lr_low_om"   ).c_str(), 1200 , 900);
        c_sr_om[j]       = new TCanvas(("cs" + region + "sr_om"       ).c_str(), ("cs" + region + "sr_om"       ).c_str(), 1200 , 900);
        c_lr_om[j]       = new TCanvas(("cs" + region + "lr_om"       ).c_str(), ("cs" + region + "lr_om"       ).c_str(), 1200 , 900);

        c_lr_low_Fourier_xi[j] = new TCanvas(("cs" + region + "lr_low_Fourier_xi").c_str(),("cs" + region + "lr_low_Fourier_xi").c_str(), 1200,900);
        c_lr_high_Fourier_xi[j] = new TCanvas(("cs" + region + "lr_high_Fourier_xi").c_str(),("cs" + region + "lr_high_Fourier_xi").c_str(), 1200,900);

        c_low_2Dks_1[j] = new TCanvas(("cs" + region + "_low_2Dks_1").c_str(),("cs" + region + "_low_2Dks_1").c_str(),1200,900);
        c_low_2Dks_2[j] = new TCanvas(("cs" + region + "_low_2Dks_2").c_str(),("cs" + region + "_low_2Dks_2").c_str(),1200,900);

        c_high_2Dks_1[j] = new TCanvas(("cs" + region + "_high_2Dks_1").c_str(),("cs" + region + "_high_2Dks_1").c_str(),1200,900);
        c_high_2Dks_2[j] = new TCanvas(("cs" + region + "_high_2Dks_2").c_str(),("cs" + region + "_high_2Dks_2").c_str(),1200,900);

        c_low_2Dla_1[j] = new TCanvas(("cs" + region + "_low_2Dla_1").c_str(),("cs" + region + "_low_2Dla_1").c_str(),1200,900);
        c_low_2Dla_2[j] = new TCanvas(("cs" + region + "_low_2Dla_2").c_str(),("cs" + region + "_low_2Dla_2").c_str(),1200,900);

        c_high_2Dla_1[j] = new TCanvas(("cs" + region + "_high_2Dla_1").c_str(),("cs" + region + "_high_2Dla_1").c_str(),1200,900);
        c_high_2Dla_2[j] = new TCanvas(("cs" + region + "_high_2Dla_2").c_str(),("cs" + region + "_high_2Dla_2").c_str(),1200,900);

        c_low_2Dxi_1[j] = new TCanvas(("cs" + region + "_low_2Dxi_1").c_str(),("cs" + region + "_low_2Dxi_1").c_str(),1200,900);
        c_low_2Dxi_2[j] = new TCanvas(("cs" + region + "_low_2Dxi_2").c_str(),("cs" + region + "_low_2Dxi_2").c_str(),1200,900);

        c_high_2Dxi_1[j] = new TCanvas(("cs" + region + "_high_2Dxi_1").c_str(),("cs" + region + "_high_2Dxi_1").c_str(),1200,900);
        c_high_2Dxi_2[j] = new TCanvas(("cs" + region + "_high_2Dxi_2").c_str(),("cs" + region + "_high_2Dxi_2").c_str(),1200,900);

        c_low_2Dom_1[j] = new TCanvas(("cs" + region + "_low_2Dom_1").c_str(),("cs" + region + "_low_2Dom_1").c_str(),1200,900);
        c_low_2Dom_2[j] = new TCanvas(("cs" + region + "_low_2Dom_2").c_str(),("cs" + region + "_low_2Dom_2").c_str(),1200,900);

        c_high_2Dom_1[j] = new TCanvas(("cs" + region + "_high_2Dom_1").c_str(),("cs" + region + "_high_2Dom_1").c_str(),1200,900);
        c_high_2Dom_2[j] = new TCanvas(("cs" + region + "_high_2Dom_2").c_str(),("cs" + region + "_high_2Dom_2").c_str(),1200,900);

        c_sr_low_ks[j]->Divide(4,4);
        c_lr_low_ks[j]->Divide(4,4);
        c_sr_ks[j]->Divide(4,4);
        c_lr_ks[j]->Divide(4,4);
        c_sr_low_la[j]->Divide(4,3);
        c_lr_low_la[j]->Divide(4,3);
        c_sr_la[j]->Divide(4,3);
        c_lr_la[j]->Divide(4,3);
        c_sr_low_xi[j]->Divide(4,3);
        c_lr_low_xi[j]->Divide(4,3);
        c_sr_xi[j]->Divide(4,3);
        c_lr_xi[j]->Divide(4,3);
        c_sr_low_om[j]->Divide(4,3);
        c_lr_low_om[j]->Divide(4,3);
        c_sr_om[j]->Divide(4,3);
        c_lr_om[j]->Divide(4,3);

        c_lr_low_Fourier_xi[j]->Divide(4,3);
        c_lr_high_Fourier_xi[j]->Divide(4,3);

        c_low_2Dks_1[j]  -> Divide(4,2);
        c_low_2Dks_2[j]  -> Divide(4,2);
        c_high_2Dks_1[j] -> Divide(4,2);
        c_high_2Dks_2[j] -> Divide(4,2);

        c_low_2Dla_1[j]  -> Divide(4,2);
        c_low_2Dla_2[j]  -> Divide(4,2);
        c_high_2Dla_1[j] -> Divide(4,2);
        c_high_2Dla_2[j] -> Divide(4,2);

        c_low_2Dxi_1[j]  -> Divide(4,2);
        c_low_2Dxi_2[j]  -> Divide(4,2);
        c_high_2Dxi_1[j] -> Divide(4,2);
        c_high_2Dxi_2[j] -> Divide(4,2);

        c_low_2Dom_1[j]  -> Divide(4,2);
        c_low_2Dom_2[j]  -> Divide(4,2);
        c_high_2Dom_1[j] -> Divide(4,2);
        c_high_2Dom_2[j] -> Divide(4,2);
    }

    //TCanvas* csPeak_low_2Dks_1 = new TCanvas("csPeak_low_2Dks_1","csPeak_low_2Dks_1",1200,900);
    //TCanvas* csPeak_low_2Dks_2 = new TCanvas("csPeak_low_2Dks_2","csPeak_low_2Dks_2",1200,900);
    //TCanvas* csSide_low_2Dks_1 = new TCanvas("csSide_low_2Dks_1","csSide_low_2Dks_1",1200,900);
    //TCanvas* csSide_low_2Dks_2 = new TCanvas("csSide_low_2Dks_2","csSide_low_2Dks_2",1200,900);

    //TCanvas* csPeak_high_2Dks_1 = new TCanvas("csPeak_high_2Dks_1","csPeak_high_2Dks_1",1200,900);
    //TCanvas* csPeak_high_2Dks_2 = new TCanvas("csPeak_high_2Dks_2","csPeak_high_2Dks_2",1200,900);
    //TCanvas* csSide_high_2Dks_1 = new TCanvas("csSide_high_2Dks_1","csSide_high_2Dks_1",1200,900);
    //TCanvas* csSide_high_2Dks_2 = new TCanvas("csSide_high_2Dks_2","csSide_high_2Dks_2",1200,900);

    //TCanvas* csPeak_low_2Dla_1 = new TCanvas("csPeak_low_2Dla_1","csPeak_low_2Dla_1",1200,900);
    //TCanvas* csPeak_low_2Dla_2 = new TCanvas("csPeak_low_2Dla_2","csPeak_low_2Dla_2",1200,900);
    //TCanvas* csSide_low_2Dla_1 = new TCanvas("csSide_low_2Dla_1","csSide_low_2Dla_1",1200,900);
    //TCanvas* csSide_low_2Dla_2 = new TCanvas("csSide_low_2Dla_2","csSide_low_2Dla_2",1200,900);

    //TCanvas* csPeak_high_2Dla_1 = new TCanvas("csPeak_high_2Dla_1","csPeak_high_2Dla_1",1200,900);
    //TCanvas* csPeak_high_2Dla_2 = new TCanvas("csPeak_high_2Dla_2","csPeak_high_2Dla_2",1200,900);
    //TCanvas* csSide_high_2Dla_1 = new TCanvas("csSide_high_2Dla_1","csSide_high_2Dla_1",1200,900);
    //TCanvas* csSide_high_2Dla_2 = new TCanvas("csSide_high_2Dla_2","csSide_high_2Dla_2",1200,900);
    TCanvas* c = new TCanvas("c","c",400,400);



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

    hbPeaksr_low_ref = hbackgroundPeak_low_ref->ProjectionY("hbPeaksr_low_ref", PD.sr_low, PD.sr_high);
    hsPeaksr_low_ref = hsignalPeak_low_ref->ProjectionY("hsPeaksr_low_ref", PD.sr_low, PD.sr_high);

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
    c_sr_low_ref->cd();
    hsPeaksr_low_ref->Draw();
    double xcoor = 0.52;
    double ycoor = 0.90;
    double increment = 0.07;
    tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
    tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",0.3,3.0));
    tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
    c->cd();
    Jyieldsr_low_ref.push_back(hsPeaksr_zeroed_low_ref->IntegralAndError(hsPeaksr_zeroed_low_ref->FindBin(0.0),hsPeaksr_zeroed_low_ref->FindBin(minVal_srX),Jyieldsr_err_low_ref,"width"));
    double bin0yield = hsPeaksr_zeroed_low_ref->GetBinContent(hsPeaksr_zeroed_low_ref->FindBin(0.0))*0.19635;
    Jyieldsr_low_ref[0] = Jyieldsr_low_ref[0]*2 - bin0yield;

    hbPeaklr_low_ref = hbackgroundPeak_low_ref->ProjectionY("hbPeaklr_low_ref",1,10);
    TH1D* ahbPeaklr_low_ref = hbackgroundPeak_low_ref->ProjectionY("ahbPeaklr_low_ref",24,33);
    hsPeaklr_low_ref = hsignalPeak_low_ref->ProjectionY("hsPeaklr_low_ref",1,10);
    TH1D* ahsPeaklr_low_ref = hsignalPeak_low_ref->ProjectionY("ahsPeaklr_low_ref",24,33);

    hbPeaklr_low_ref->Add(ahbPeaklr_low_ref);
    hsPeaklr_low_ref->Add(ahsPeaklr_low_ref);
    hsPeaklr_low_ref->Divide(hbPeaklr_low_ref);
    hsPeaklr_low_ref->Scale(Bz_low_ref[0]/nEvent_low_ref[0]/BW2D);

    TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.0,2.0);
    quadFit2->SetParameters(1,1,1);

    hsPeaklr_low_ref->Fit("quadFit2","R");
    hsPeaklr_low_ref->Fit("quadFit2","R");
    hsPeaklr_low_ref->Fit("quadFit2","R");

    double minVal_lr = quadFit2->GetMinimum(0.0,2.0);
    double minVal_lrX = quadFit2->GetMinimumX(0.0,2.0);
    TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    minConst_lr->SetParameter(0,-minVal_lr);
    TH1D* hsPeaklr_zeroed_low_ref = (TH1D*)hsPeaklr_low_ref->Clone();
    hsPeaklr_zeroed_low_ref->Add(minConst_lr);
    c_lr_low_ref->cd();
    hsPeaklr_low_ref->Draw();
    xcoor = 0.52;
    ycoor = 0.90;
    increment = 0.07;
    tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
    tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",0.3,3.0));
    tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|>2");
    c->cd();
    Jyieldlr_low_ref.push_back(hsPeaklr_zeroed_low_ref->IntegralAndError(hsPeaklr_zeroed_low_ref->FindBin(0.0),hsPeaklr_zeroed_low_ref->FindBin(minVal_lrX),Jyieldlr_err_low_ref,"width"));
    bin0yield = hsPeaklr_zeroed_low_ref->GetBinContent(hsPeaklr_zeroed_low_ref->FindBin(0.0))*0.19635;
    Jyieldlr_low_ref[0] = Jyieldlr_low_ref[0]*2 - bin0yield;

    JyieldSub_low_ref.push_back(Jyieldsr_low_ref[0] - Jyieldlr_low_ref[0]);
    JyieldSub_err_low_ref.push_back(sqrt(Jyieldsr_err_low_ref*Jyieldsr_err_low_ref + Jyieldlr_err_low_ref*Jyieldlr_err_low_ref));

    //TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
    fit1->SetParNames("N","V1","V2","V3");
    fit1->SetParameters(10,1,1,1);
    fit1->SetLineColor(2);

    V2lrs_low_ref = hsignalPeak_low_ref->ProjectionY("V2lrs_low_ref",1,PD.binlow_V2);
    TH1D* aV2lrs_low_ref = hsignalPeak_low_ref->ProjectionY("aV2lrs_low_ref",PD.binhigh_V2,33);
    V2lrb_low_ref = hbackgroundPeak_low_ref->ProjectionY("V2lrb_low_ref",1,PD.binlow_V2);
    TH1D* aV2lrb_low_ref = hbackgroundPeak_low_ref->ProjectionY("aV2lrb_low_ref",PD.binhigh_V2,33);
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

    //PD.f_high_ref->GetObject("xiCorrelationRapidity/BackgroundHad",hbackgroundPeak_high_ref);
    //PD.f_high_ref->GetObject("xiCorrelationRapidity/SignalHad",hsignalPeak_high_ref);
    //TH1D* mult_high_ref = (TH1D*) PD.f_high_ref->Get("xiCorrelationRapidity/nTrk");
    PD.f_high_ref->GetObject("pPbCorr/background",hbackgroundPeak_high_ref);
    PD.f_high_ref->GetObject("pPbCorr/signal",hsignalPeak_high_ref);
    TH1D* mult_high_ref = (TH1D*) PD.f_high_ref->Get("pPbCorr/mult");

    hbPeaksr_high_ref = hbackgroundPeak_high_ref->ProjectionY("hbPeaksr_high_ref", PD.sr_low, PD.sr_high);
    hsPeaksr_high_ref = hsignalPeak_high_ref->ProjectionY("hsPeaksr_high_ref", PD.sr_low, PD.sr_high);

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
    c_sr_high_ref->cd();
    hsPeaksr_high_ref->Draw();
    xcoor = 0.52;
    ycoor = 0.90;
    increment = 0.07;
    tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
    tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",0.3,3.0));
    tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
    c->cd();
    Jyieldsr_high_ref.push_back(hsPeaksr_zeroed_high_ref->IntegralAndError(hsPeaksr_zeroed_high_ref->FindBin(0.0),hsPeaksr_zeroed_high_ref->FindBin(minVal_srX),Jyieldsr_err_high_ref,"width"));
    bin0yield = hsPeaksr_zeroed_high_ref->GetBinContent(hsPeaksr_zeroed_high_ref->FindBin(0.0))*0.19635;
    Jyieldsr_high_ref[0] = Jyieldsr_high_ref[0]*2 - bin0yield;

    hbPeaklr_high_ref = hbackgroundPeak_high_ref->ProjectionY("hbPeaklr_high_ref",1,10);
    TH1D* ahbPeaklr_high_ref = hbackgroundPeak_high_ref->ProjectionY("ahbPeaklr_high_ref",24,33);
    hsPeaklr_high_ref = hsignalPeak_high_ref->ProjectionY("hsPeaklr_high_ref",1,10);
    TH1D* ahsPeaklr_high_ref = hsignalPeak_high_ref->ProjectionY("ahsPeaklr_high_ref",24,33);

    hbPeaklr_high_ref->Add(ahbPeaklr_high_ref);
    hsPeaklr_high_ref->Add(ahsPeaklr_high_ref);
    hsPeaklr_high_ref->Divide(hbPeaklr_high_ref);
    hsPeaklr_high_ref->Scale(Bz_high_ref[0]/nEvent_high_ref[0]/BW2D);

    TF1* quadFit21 = new TF1("quadFit21","[0]*x^2+[1]*x+[2]",0.0,2.0);
    quadFit21->SetParameters(1,1,1);

    hsPeaklr_high_ref->Fit("quadFit21","R");
    hsPeaklr_high_ref->Fit("quadFit21","R");
    hsPeaklr_high_ref->Fit("quadFit21","R");

    minVal_lr = quadFit21->GetMinimum(0.6,2.0);
    minVal_lrX = quadFit21->GetMinimumX(0.6,2.0);
    TF1* minConst_lr2 = new TF1("minConst_lr2","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
    minConst_lr2->SetParameter(0,-minVal_lr);
    TH1D* hsPeaklr_zeroed_high_ref = (TH1D*)hsPeaklr_high_ref->Clone();
    hsPeaklr_zeroed_high_ref->Add(minConst_lr2);
    c_lr_high_ref->cd();
    hsPeaklr_high_ref->Draw();
    xcoor = 0.52;
    ycoor = 0.90;
    increment = 0.07;
    tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
    tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",0.3,3.0));
    tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
    tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|>2");
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

    V2lrs_high_ref = hsignalPeak_high_ref->ProjectionY("V2lrs_high_ref",1,PD.binlow_V2);
    TH1D* aV2lrs_high_ref = hsignalPeak_high_ref->ProjectionY("aV2lrs_high_ref",PD.binhigh_V2,33);
    V2lrb_high_ref = hbackgroundPeak_high_ref->ProjectionY("V2lrb_high_ref",1,PD.binlow_V2);
    TH1D* aV2lrb_high_ref = hbackgroundPeak_high_ref->ProjectionY("aV2lrb_high_ref",PD.binhigh_V2,33);
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
    std::string PathMult_high_xi;

    std::string PathBackground_low_om;
    std::string PathSignal_low_om;
    std::string PathMult_low_om;
    std::string PathBackground_high_om;
    std::string PathSignal_high_om;
    std::string PathMult_high_om;

    std::string region_label;

    int numPtBins_ks = PD.PtBin_ks.size()-1;
    int numPtBins_la = PD.PtBin_la.size()-1;
    int numPtBins_xi = PD.PtBin_xi.size()-1;
    int numPtBins_om = PD.PtBin_om.size()-1;


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
            PathBackground_high_xi  = PD.fn_Xi + "/BackgroundXiPeak_pt%d";
            PathSignal_high_xi = PD.fn_Xi + "/SignalXiPeak_pt%d";
            PathMult_high_xi = PD.fn_Xi + "/mult_xi_pt%d";

            PathBackground_low_om  = PD.fn_v0cas + "/BackgroundOmPeak_pt%d";
            PathSignal_low_om = PD.fn_v0cas + "/SignalOmPeak_pt%d";
            PathMult_low_om = PD.fn_v0cas + "/mult_om_pt%d";
            PathBackground_high_om  = PD.fn_Om + "/BackgroundOmPeak_pt%d";
            PathSignal_high_om = PD.fn_Om + "/SignalOmPeak_pt%d";
            PathMult_high_om = PD.fn_Om + "/mult_om_pt%d";

            region_label = "peak";
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
            PathBackground_high_xi  = PD.fn_Xi + "/BackgroundXiSide_pt%d";
            PathSignal_high_xi = PD.fn_Xi + "/SignalXiSide_pt%d";
            PathMult_high_xi = PD.fn_Xi + "/mult_xi_bkg_pt%d";

            PathBackground_low_om  = PD.fn_v0cas + "/BackgroundOmSide_pt%d";
            PathSignal_low_om = PD.fn_v0cas + "/SignalOmSide_pt%d";
            PathMult_low_om = PD.fn_v0cas + "/mult_om_bkg_pt%d";
            PathBackground_high_om  = PD.fn_Om + "/BackgroundOmSide_pt%d";
            PathSignal_high_om = PD.fn_Om + "/SignalOmSide_pt%d";
            PathMult_high_om = PD.fn_Om + "/mult_om_bkg_pt%d";

            region_label = "side";

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

        TH1D* hsPeaksr_low_ks[3*arraySize_ks];
        TH1D* hbPeaksr_low_ks[3*arraySize_ks];
        TH1D* hsPeaklr_low_ks[3*arraySize_ks];
        TH1D* hbPeaklr_low_ks[3*arraySize_ks];
        TH1D* V2lrs_low_ks[3*arraySize_ks];
        TH1D* V2lrb_low_ks[3*arraySize_ks];
        TH2D* hbackgroundPeak_low_ks[3*arraySize_ks];
        TH2D* hsignalPeak_low_ks[3*arraySize_ks];


        //Calculate Nassoc, Jet yield, Peak region Low N

        for(int i=0; i<numPtBins_ks; i++)
        {
            PD.f_perisub->GetObject(Form(PathBackground_low_ks.c_str(),i),hbackgroundPeak_low_ks[i]);
            PD.f_perisub->GetObject(Form(PathSignal_low_ks.c_str(),i),hsignalPeak_low_ks[i]);
            TH1D* mult_low_ks = (TH1D*) PD.f_perisub->Get(Form(PathMult_low_ks.c_str(),i));

            hbPeaksr_low_ks[i] = hbackgroundPeak_low_ks[i]->ProjectionY(Form(("hb"+region_label+"sr_low_ks%d").c_str(),i), PD.sr_low, PD.sr_high);
            hsPeaksr_low_ks[i] = hsignalPeak_low_ks[i]->ProjectionY(Form(("hs"+region_label+"sr_low_ks%d").c_str(),i), PD.sr_low, PD.sr_high);

            TF1* quadFit13 = new TF1("quadFit13","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit13->SetParameters(1,1,1);

            if(i<8)
            {
                c_low_2Dks_1[j]->cd(i+1);
                //TH2D* hPeak_low_2Dks = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_ks.c_str(),i));
                TH2D* hPeak_low_2Dks = (TH2D*)hsignalPeak_low_ks[i]->Clone();
                hPeak_low_2Dks->SetName(Form(("h"+region_label+"_low_2Dks%d").c_str(),i));
                hPeak_low_2Dks->Divide(hbackgroundPeak_low_ks[i]);

                hPeak_low_2Dks->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dks->Draw("SURF1");
                //hbackgroundPeak_low_ks[i]->Draw();
            }
            if(i>=8)
            {
                c_low_2Dks_2[j]->cd(i-8+1);
                //TH2D* hPeak_low_2Dks = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_ks.c_str(),i));
                TH2D* hPeak_low_2Dks = (TH2D*)hsignalPeak_low_ks[i]->Clone();
                hPeak_low_2Dks->SetName(Form(("h"+region_label+"_low_2Dks%d").c_str(),i));
                hPeak_low_2Dks->Divide(hbackgroundPeak_low_ks[i]);
                hPeak_low_2Dks->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dks->Draw("SURF1");
            }
            c->cd();


            //nEvent_low_ks.push_back(mult_low_ks->Integral(0,100000));
            nEvent_low_ks.push_back(mult_low_ks->Integral(2,100000));
            Bz_low_ks.push_back(hbackgroundPeak_low_ks[i]->GetBinContent(hbackgroundPeak_low_ks[i]->FindBin(0,0)));

            hsPeaksr_low_ks[i]->Divide(hbPeaksr_low_ks[i]);
            hsPeaksr_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);
            //hsPeaksr_low_ks[i]->Scale(Bz_low_ks[i]/BW2D);

            hsPeaksr_low_ks[i]->Fit("quadFit13","R");
            hsPeaksr_low_ks[i]->Fit("quadFit13","R");
            hsPeaksr_low_ks[i]->Fit("quadFit13","R");

            c_sr_low_ks[j]->cd(i+1);
            double minVal_sr = quadFit13->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit13->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            TH1D* hsPeaksr_zeroed_low_ks = (TH1D*)hsPeaksr_low_ks[i]->Clone();
            TH1D* hsPeaksr_drawed_low_ks = (TH1D*)hsPeaksr_low_ks[i]->Clone();
            hsPeaksr_drawed_low_ks->Draw();
            hsPeaksr_zeroed_low_ks->Add(minConst_sr);
            //hsPeaksr_low_ks[i]->Draw();
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_ks[i],PD.PtBin_ks[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_low_ks.push_back(hsPeaksr_zeroed_low_ks->IntegralAndError(hsPeaksr_zeroed_low_ks->FindBin(0),hsPeaksr_zeroed_low_ks->FindBin(minVal_srX),Jyieldsr_err_low_ks[i],"width"));
            double bin0yield = hsPeaksr_zeroed_low_ks->GetBinContent(hsPeaksr_zeroed_low_ks->FindBin(0.0))*0.19635;
            Jyieldsr_low_ks[i] = Jyieldsr_low_ks[i]*2 - bin0yield;

            hbPeaklr_low_ks[i] = hbackgroundPeak_low_ks[i]->ProjectionY("hbPeaklr_low_ks",1,10);
            TH1D* ahbPeaklr_low_ks = hbackgroundPeak_low_ks[i]->ProjectionY("ahbPeaklr_low_ks",24,33);
            c_lr_low_ks[j]->cd(i+1);
            hsPeaklr_low_ks[i] = hsignalPeak_low_ks[i]->ProjectionY("hsPeaklr_low_ks",1,10);
            TH1D* ahsPeaklr_low_ks = hsignalPeak_low_ks[i]->ProjectionY("ahsPeaklr_low_ks",24,33);

            hbPeaklr_low_ks[i]->Add(ahbPeaklr_low_ks);
            hsPeaklr_low_ks[i]->Add(ahsPeaklr_low_ks);
            hsPeaklr_low_ks[i]->Divide(hbPeaklr_low_ks[i]);
            hsPeaklr_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);
            //hsPeaklr_low_ks[i]->Scale(Bz_low_ks[i]/BW2D);

            TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.0,2.0);
            quadFit2->SetParameters(1,1,1);

            hsPeaklr_low_ks[i]->Fit("quadFit2","R");
            //hsPeaklr_low_ks[i]->Fit("quadFit2","R");
            //hsPeaklr_low_ks[i]->Fit("quadFit2","R");

            double minVal_lr = quadFit2->GetMinimum(0.0,2.0);
            double minVal_lrX = quadFit2->GetMinimumX(0.0,2.0);
            TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_lr->SetParameter(0,-minVal_lr);
            //hsPeaklr_low_ks[i]->Draw();
            TH1D* hsPeaklr_zeroed_low_ks = (TH1D*)hsPeaklr_low_ks[i]->Clone();
            TH1D* hsPeaklr_drawed_low_ks = (TH1D*)hsPeaklr_low_ks[i]->Clone();
            hsPeaklr_drawed_low_ks->Draw();
            hsPeaklr_zeroed_low_ks->Add(minConst_lr);
            xcoor = 0.52;
            ycoor = 0.90;
            increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_ks[i],PD.PtBin_ks[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|>2");
            Jyieldlr_low_ks.push_back(hsPeaklr_zeroed_low_ks->IntegralAndError(hsPeaklr_zeroed_low_ks->FindBin(0.0),hsPeaklr_zeroed_low_ks->FindBin(minVal_lrX),Jyieldlr_err_low_ks[i],"width"));
            bin0yield = hsPeaklr_zeroed_low_ks->GetBinContent(hsPeaklr_zeroed_low_ks->FindBin(0.0))*0.19635;
            Jyieldlr_low_ks[i] = Jyieldlr_low_ks[i]*2 - bin0yield;

            JyieldSub_low_ks.push_back(Jyieldsr_low_ks[i] - Jyieldlr_low_ks[i]);
            JyieldSub_err_low_ks.push_back(sqrt(Jyieldsr_err_low_ks[i]*Jyieldsr_err_low_ks[i] + Jyieldlr_err_low_ks[i]*Jyieldlr_err_low_ks[i]));
            c->cd();

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_low_ks[i] = hsignalPeak_low_ks[i]->ProjectionY(Form("V2lrs_low_ks%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_low_ks = hsignalPeak_low_ks[i]->ProjectionY("aV2lrs_low_ks",PD.binhigh_V2,33);
            V2lrb_low_ks[i] = hbackgroundPeak_low_ks[i]->ProjectionY(Form("V2lrb_low_ks%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_low_ks = hbackgroundPeak_low_ks[i]->ProjectionY("aV2lrb_low_ks",PD.binhigh_V2,33);
            V2lrs_low_ks[i]->Add(aV2lrs_low_ks);
            V2lrb_low_ks[i]->Add(aV2lrb_low_ks);
            V2lrs_low_ks[i]->Divide(V2lrb_low_ks[i]);
            V2lrs_low_ks[i]->Scale(Bz_low_ks[i]/nEvent_low_ks[i]/BW2D);
            //V2lrs_low_ks[i]->Scale(Bz_low_ks[i]/BW2D);

            V2lrs_low_ks[i]->Fit("fit1","R");
            V2lrs_low_ks[i]->Fit("fit1","R");
            V2lrs_low_ks[i]->Fit("fit1","R");

            V2Values_low_ks.push_back(fit1->GetParameter(2));
            V2Values_err_low_ks.push_back(fit1->GetParError(2));
            V3Values_low_ks.push_back(fit1->GetParameter(3));
            V3Values_err_low_ks.push_back(fit1->GetParError(3));

            Nassoc_low_ks.push_back(fit1->GetParameter(0));
            c->cd();
            //if(i==4) return;
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

            hbPeaksr_low_la[i] = hbackgroundPeak_low_la[i]->ProjectionY(("hb"+region_label+"sr_low_la").c_str(), PD.sr_low, PD.sr_high);
            hsPeaksr_low_la[i] = hsignalPeak_low_la[i]->ProjectionY(("hs"+region_label+"sr_low_la").c_str(), PD.sr_low, PD.sr_high);

            if(i<8)
            {
                c_low_2Dla_1[j]->cd(i+1);
                //TH2D* hPeak_low_2Dla = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_la.c_str(),i));
                TH2D* hPeak_low_2Dla = (TH2D*)hsignalPeak_low_la[i]->Clone();
                hPeak_low_2Dla->SetName(Form(("h"+region_label+"_low_2Dla%d").c_str(),i));
                hPeak_low_2Dla->Divide(hbackgroundPeak_low_la[i]);

                hPeak_low_2Dla->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dla->Draw("SURF1");
                //hbackgroundPeak_low_la[i]->Draw();
            }
            if(i>=8)
            {
                c_low_2Dla_2[j]->cd(i-8+1);
                //TH2D* hPeak_low_2Dla = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_la.c_str(),i));
                TH2D* hPeak_low_2Dla = (TH2D*)hsignalPeak_low_la[i]->Clone();
                hPeak_low_2Dla->SetName(Form(("h"+region_label+"_low_2Dla%d").c_str(),i));
                hPeak_low_2Dla->Divide(hbackgroundPeak_low_la[i]);
                hPeak_low_2Dla->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dla->Draw("SURF1");
            }

            TF1* quadFit14 = new TF1("quadFit14","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit14->SetParameters(1,1,1);

            //nEvent_low_la.push_back(mult_low_la->Integral(0,100000));
            nEvent_low_la.push_back(mult_low_la->Integral(2,100000));
            Bz_low_la.push_back(hbackgroundPeak_low_la[i]->GetBinContent(hbackgroundPeak_low_la[i]->FindBin(0,0)));

            hsPeaksr_low_la[i]->Divide(hbPeaksr_low_la[i]);
            hsPeaksr_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);
            //hsPeaksr_low_la[i]->Scale(Bz_low_la[i]/BW2D);

            c_sr_low_la[j]->cd(i+1);

            hsPeaksr_low_la[i]->Fit("quadFit14","R");
            //hsPeaksr_low_la[i]->Fit("quadFit14","R");
            //hsPeaksr_low_la[i]->Fit("quadFit14","R");

            double minVal_sr = quadFit14->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit14->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            TH1D* hsPeaksr_zeroed_low_la = (TH1D*)hsPeaksr_low_la[i]->Clone();
            TH1D* hsPeaksr_drawed_low_la = (TH1D*)hsPeaksr_low_la[i]->Clone();
            hsPeaksr_drawed_low_la->Draw();
            hsPeaksr_zeroed_low_la->Add(minConst_sr);
            //hsPeaksr_low_la[i]->Draw();
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_la[i],PD.PtBin_la[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_low_la.push_back(hsPeaksr_zeroed_low_la->IntegralAndError(hsPeaksr_zeroed_low_la->FindBin(0.0),hsPeaksr_zeroed_low_la->FindBin(minVal_srX),Jyieldsr_err_low_la[i],"width"));
            double bin0yield = hsPeaksr_zeroed_low_la->GetBinContent(hsPeaksr_zeroed_low_la->FindBin(0.0))*0.19635;
            Jyieldsr_low_la[i] = Jyieldsr_low_la[i]*2 - bin0yield;

            hbPeaklr_low_la[i] = hbackgroundPeak_low_la[i]->ProjectionY("hbPeaklr_low_la",1,10);
            TH1D* ahbPeaklr_low_la = hbackgroundPeak_low_la[i]->ProjectionY("ahbPeaklr_low_la",24,33);
            hsPeaklr_low_la[i] = hsignalPeak_low_la[i]->ProjectionY("hsPeaklr_low_la",1,10);
            TH1D* ahsPeaklr_low_la = hsignalPeak_low_la[i]->ProjectionY("ahsPeaklr_low_la",24,33);

            hbPeaklr_low_la[i]->Add(ahbPeaklr_low_la);
            hsPeaklr_low_la[i]->Add(ahsPeaklr_low_la);
            hsPeaklr_low_la[i]->Divide(hbPeaklr_low_la[i]);
            hsPeaklr_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);
            //hsPeaklr_low_la[i]->Scale(Bz_low_la[i]/BW2D);

            TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.0,2.0);
            quadFit2->SetParameters(1,1,1);

            hsPeaklr_low_la[i]->Fit("quadFit2","R");
            //hsPeaklr_low_la[i]->Fit("quadFit2","R");
            //hsPeaklr_low_la[i]->Fit("quadFit2","R");

            c_lr_low_la[j]->cd(i+1);
            double minVal_lr = quadFit2->GetMinimum(0.0,2.0);
            double minVal_lrX = quadFit2->GetMinimumX(0.0,2.0);
            TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_lr->SetParameter(0,-minVal_lr);
            TH1D* hsPeaklr_zeroed_low_la = (TH1D*)hsPeaklr_low_la[i]->Clone();
            TH1D* hsPeaklr_drawed_low_la = (TH1D*)hsPeaklr_low_la[i]->Clone();
            hsPeaklr_drawed_low_la->Draw();
            hsPeaklr_zeroed_low_la->Add(minConst_lr);
            xcoor = 0.52;
            ycoor = 0.90;
            increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_la[i],PD.PtBin_la[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|>2");
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

            V2lrs_low_la[i] = hsignalPeak_low_la[i]->ProjectionY(Form("V2lrs_low_la%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_low_la = hsignalPeak_low_la[i]->ProjectionY("aV2lrs_low_la",PD.binhigh_V2,33);
            V2lrb_low_la[i] = hbackgroundPeak_low_la[i]->ProjectionY(Form("V2lrb_low_la%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_low_la = hbackgroundPeak_low_la[i]->ProjectionY("aV2lrb_low_la",PD.binhigh_V2,33);
            V2lrs_low_la[i]->Add(aV2lrs_low_la);
            V2lrb_low_la[i]->Add(aV2lrb_low_la);
            V2lrs_low_la[i]->Divide(V2lrb_low_la[i]);
            V2lrs_low_la[i]->Scale(Bz_low_la[i]/nEvent_low_la[i]/BW2D);
            //V2lrs_low_la[i]->Scale(Bz_low_la[i]/BW2D);

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
            PD.f_perisub_xi->GetObject(Form(PathBackground_low_xi.c_str(),i),hbackgroundPeak_low_xi[i]);
            PD.f_perisub_xi->GetObject(Form(PathSignal_low_xi.c_str(),i),hsignalPeak_low_xi[i]);
            TH1D* mult_low_xi = (TH1D*) PD.f_perisub_xi->Get(Form(PathMult_low_xi.c_str(),i));

            hbPeaksr_low_xi[i] = hbackgroundPeak_low_xi[i]->ProjectionY(Form(("hb"+region_label+"sr_low_xi%d").c_str(),i), PD.sr_low, PD.sr_high);
            hsPeaksr_low_xi[i] = hsignalPeak_low_xi[i]->ProjectionY(Form(("hs"+region_label+"sr_low_xi%d").c_str(),i), PD.sr_low, PD.sr_high);

            TF1* quadFit13 = new TF1("quadFit13","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit13->SetParameters(1,1,1);

            if(i<8)
            {
                c_low_2Dxi_1[j]->cd(i+1);
                //TH2D* hPeak_low_2Dxi = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_xi.c_str(),i));
                TH2D* hPeak_low_2Dxi = (TH2D*)hsignalPeak_low_xi[i]->Clone();
                hPeak_low_2Dxi->SetName(Form(("h"+region_label+"_low_2Dxi%d").c_str(),i));
                hPeak_low_2Dxi->Divide(hbackgroundPeak_low_xi[i]);

                hPeak_low_2Dxi->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dxi->Draw("SURF1");
                //hbackgroundPeak_low_xi[i]->Draw();
            }
            if(i>=8)
            {
                c_low_2Dxi_2[j]->cd(i-8+1);
                //TH2D* hPeak_low_2Dxi = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_xi.c_str(),i));
                TH2D* hPeak_low_2Dxi = (TH2D*)hsignalPeak_low_xi[i]->Clone();
                hPeak_low_2Dxi->SetName(Form(("h"+region_label+"_low_2Dxi%d").c_str(),i));
                hPeak_low_2Dxi->Divide(hbackgroundPeak_low_xi[i]);
                hPeak_low_2Dxi->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dxi->Draw("SURF1");
            }


            //nEvent_low_xi.push_back(mult_low_xi->Integral(0,100000));
            nEvent_low_xi.push_back(mult_low_xi->Integral(2,100000));
            Bz_low_xi.push_back(hbackgroundPeak_low_xi[i]->GetBinContent(hbackgroundPeak_low_xi[i]->FindBin(0,0)));

            hsPeaksr_low_xi[i]->Divide(hbPeaksr_low_xi[i]);
            hsPeaksr_low_xi[i]->Scale(Bz_low_xi[i]/nEvent_low_xi[i]/BW2D);
            //hsPeaksr_low_xi[i]->Scale(Bz_low_xi[i]/BW2D);

            c_sr_low_xi[j]->cd(i+1);

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
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_xi[i],PD.PtBin_xi[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_low_xi.push_back(hsPeaksr_zeroed_low_xi->IntegralAndError(hsPeaksr_zeroed_low_xi->FindBin(0),hsPeaksr_zeroed_low_xi->FindBin(minVal_srX),Jyieldsr_err_low_xi[i],"width"));
            double bin0yield = hsPeaksr_zeroed_low_xi->GetBinContent(hsPeaksr_zeroed_low_xi->FindBin(0.0))*0.19635;
            Jyieldsr_low_xi[i] = Jyieldsr_low_xi[i]*2 - bin0yield;

            hbPeaklr_low_xi[i] = hbackgroundPeak_low_xi[i]->ProjectionY("hbPeaklr_low_xi",1,10);
            TH1D* ahbPeaklr_low_xi = hbackgroundPeak_low_xi[i]->ProjectionY("ahbPeaklr_low_xi",24,33);
            hsPeaklr_low_xi[i] = hsignalPeak_low_xi[i]->ProjectionY(Form("hsPeaklr_low_xi%d",i),1,10);
            TH1D* ahsPeaklr_low_xi = hsignalPeak_low_xi[i]->ProjectionY("ahsPeaklr_low_xi",24,33);

            hbPeaklr_low_xi[i]->Add(ahbPeaklr_low_xi);
            hsPeaklr_low_xi[i]->Add(ahsPeaklr_low_xi);
            hsPeaklr_low_xi[i]->Divide(hbPeaklr_low_xi[i]);
            hsPeaklr_low_xi[i]->Scale(Bz_low_xi[i]/nEvent_low_xi[i]/BW2D);
            //hsPeaklr_low_xi[i]->Scale(Bz_low_xi[i]/BW2D);

            c_lr_low_xi[j]->cd(i+1);

            TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.0,2.0);
            quadFit2->SetParameters(1,1,1);

            hsPeaklr_low_xi[i]->Fit("quadFit2","R");
            //hsPeaklr_low_xi[i]->Fit("quadFit2","R");
            //hsPeaklr_low_xi[i]->Fit("quadFit2","R");

            double minVal_lr = quadFit2->GetMinimum(0.0,2.0);
            double minVal_lrX = quadFit2->GetMinimumX(0.0,2.0);
            TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_lr->SetParameter(0,-minVal_lr);
            TH1D* hsPeaklr_zeroed_low_xi = (TH1D*)hsPeaklr_low_xi[i]->Clone();
            hsPeaklr_zeroed_low_xi->Add(minConst_lr);
            hsPeaklr_low_xi[i]->Draw();
            xcoor = 0.52;
            ycoor = 0.90;
            increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_xi[i],PD.PtBin_xi[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|>2");
            Jyieldlr_low_xi.push_back(hsPeaklr_zeroed_low_xi->IntegralAndError(hsPeaklr_zeroed_low_xi->FindBin(0.0),hsPeaklr_zeroed_low_xi->FindBin(minVal_lrX),Jyieldlr_err_low_xi[i],"width"));
            bin0yield = hsPeaklr_zeroed_low_xi->GetBinContent(hsPeaklr_zeroed_low_xi->FindBin(0.0))*0.19635;
            Jyieldlr_low_xi[i] = Jyieldlr_low_xi[i]*2 - bin0yield;

            JyieldSub_low_xi.push_back(Jyieldsr_low_xi[i] - Jyieldlr_low_xi[i]);
            JyieldSub_err_low_xi.push_back(sqrt(Jyieldsr_err_low_xi[i]*Jyieldsr_err_low_xi[i] + Jyieldlr_err_low_xi[i]*Jyieldlr_err_low_xi[i]));

            //TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            //TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_low_xi[i] = hsignalPeak_low_xi[i]->ProjectionY(Form("V2lrs_low_xi%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_low_xi = hsignalPeak_low_xi[i]->ProjectionY("aV2lrs_low_xi",PD.binhigh_V2,33);
            V2lrb_low_xi[i] = hbackgroundPeak_low_xi[i]->ProjectionY(Form("V2lrb_low_xi%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_low_xi = hbackgroundPeak_low_xi[i]->ProjectionY("aV2lrb_low_xi",PD.binhigh_V2,33);
            V2lrs_low_xi[i]->Add(aV2lrs_low_xi);
            V2lrb_low_xi[i]->Add(aV2lrb_low_xi);
            V2lrs_low_xi[i]->Divide(V2lrb_low_xi[i]);
            V2lrs_low_xi[i]->Scale(Bz_low_xi[i]/nEvent_low_xi[i]/BW2D);
            //V2lrs_low_xi[i]->Scale(Bz_low_xi[i]/BW2D);
            
            c_lr_low_Fourier_xi[j]->cd(i+1);

            V2lrs_low_xi[i]->Fit("fit1","R");
            V2lrs_low_xi[i]->Fit("fit1","R");
            V2lrs_low_xi[i]->Fit("fit1","R");

            V2lrs_low_xi[i]->Draw();

            V2Values_low_xi.push_back(fit1->GetParameter(2));
            V2Values_err_low_xi.push_back(fit1->GetParError(2));
            V3Values_low_xi.push_back(fit1->GetParameter(3));
            V3Values_err_low_xi.push_back(fit1->GetParError(3));

            Nassoc_low_xi.push_back(fit1->GetParameter(0));
            c->cd();
        }

        //Omega Low N
        std::vector<double> Nassoc_low_om;
        std::vector<double> Jyieldsr_low_om;
        std::vector<double> Jyieldlr_low_om;
        std::vector<double> Jyieldsr_err_low_om(9);
        std::vector<double> Jyieldlr_err_low_om(9);
        std::vector<double> V2Values_low_om;
        std::vector<double> V2Values_err_low_om;
        std::vector<double> V3Values_low_om;
        std::vector<double> V3Values_err_low_om;
        std::vector<double> Bz_low_om;
        std::vector<double> nEvent_low_om;

        std::vector<double> JyieldSub_low_om;
        std::vector<double> JyieldSub_err_low_om;

        int arraySize_om = PD.PtBin_om.size();

        TH1D* hsPeaksr_low_om[arraySize_om];
        TH1D* hbPeaksr_low_om[arraySize_om];
        TH1D* hsPeaklr_low_om[arraySize_om];
        TH1D* hbPeaklr_low_om[arraySize_om];
        TH1D* V2lrs_low_om[arraySize_om];
        TH1D* V2lrb_low_om[arraySize_om];
        TH2D* hbackgroundPeak_low_om[arraySize_om];
        TH2D* hsignalPeak_low_om[arraySize_om];


        //Calculate Nassoc, Jet yield, Peak region Low N

        for(int i=0; i<numPtBins_om; i++)
        {
            PD.f_perisub_om->GetObject(Form(PathBackground_low_om.c_str(),i),hbackgroundPeak_low_om[i]);
            PD.f_perisub_om->GetObject(Form(PathSignal_low_om.c_str(),i),hsignalPeak_low_om[i]);
            TH1D* mult_low_om = (TH1D*) PD.f_perisub_om->Get(Form(PathMult_low_om.c_str(),i));

            hbPeaksr_low_om[i] = hbackgroundPeak_low_om[i]->ProjectionY(Form(("hb"+region_label+"sr_low_om%d").c_str(),i), PD.sr_low, PD.sr_high);
            hsPeaksr_low_om[i] = hsignalPeak_low_om[i]->ProjectionY(Form(("hs"+region_label+"sr_low_om%d").c_str(),i), PD.sr_low, PD.sr_high);

            TF1* quadFit13 = new TF1("quadFit13","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit13->SetParameters(1,1,1);

            if(i<8)
            {
                c_low_2Dom_1[j]->cd(i+1);
                //TH2D* hPeak_low_2Dom = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_om.c_str(),i));
                TH2D* hPeak_low_2Dom = (TH2D*)hsignalPeak_low_om[i]->Clone();
                hPeak_low_2Dom->SetName(Form(("h"+region_label+"_low_2Dom%d").c_str(),i));
                hPeak_low_2Dom->Divide(hbackgroundPeak_low_om[i]);

                hPeak_low_2Dom->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dom->Draw("SURF1");
                //hbackgroundPeak_low_om[i]->Draw();
            }
            if(i>=8)
            {
                c_low_2Dom_2[j]->cd(i-8+1);
                //TH2D* hPeak_low_2Dom = (TH2D*)PD.f_perisub->Get(Form(PathSignal_low_om.c_str(),i));
                TH2D* hPeak_low_2Dom = (TH2D*)hsignalPeak_low_om[i]->Clone();
                hPeak_low_2Dom->SetName(Form(("h"+region_label+"_low_2Dom%d").c_str(),i));
                hPeak_low_2Dom->Divide(hbackgroundPeak_low_om[i]);
                hPeak_low_2Dom->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_low_2Dom->Draw("SURF1");
            }

            nEvent_low_om.push_back(mult_low_om->Integral(2,100000));
            //nEvent_low_om.push_back(mult_low_om->Integral(2,100000));
            Bz_low_om.push_back(hbackgroundPeak_low_om[i]->GetBinContent(hbackgroundPeak_low_om[i]->FindBin(0,0)));

            hsPeaksr_low_om[i]->Divide(hbPeaksr_low_om[i]);
            hsPeaksr_low_om[i]->Scale(Bz_low_om[i]/nEvent_low_om[i]/BW2D);
            //hsPeaksr_low_om[i]->Scale(Bz_low_om[i]/BW2D);

            c_sr_low_om[j]->cd(i+1);

            hsPeaksr_low_om[i]->Fit("quadFit13","R");
            //hsPeaksr_low_om[i]->Fit("quadFit13","R");
            //hsPeaksr_low_om[i]->Fit("quadFit13","R");

            double minVal_sr = quadFit13->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit13->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            TH1D* hsPeaksr_zeroed_low_om = (TH1D*)hsPeaksr_low_om[i]->Clone();
            hsPeaksr_zeroed_low_om->Add(minConst_sr);
            hsPeaksr_low_om[i]->Draw();
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_om[i],PD.PtBin_om[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_low_om.push_back(hsPeaksr_zeroed_low_om->IntegralAndError(hsPeaksr_zeroed_low_om->FindBin(0),hsPeaksr_zeroed_low_om->FindBin(minVal_srX),Jyieldsr_err_low_om[i],"width"));
            double bin0yield = hsPeaksr_zeroed_low_om->GetBinContent(hsPeaksr_zeroed_low_om->FindBin(0.0))*0.19635;
            Jyieldsr_low_om[i] = Jyieldsr_low_om[i]*2 - bin0yield;

            hbPeaklr_low_om[i] = hbackgroundPeak_low_om[i]->ProjectionY("hbPeaklr_low_om",1,10);
            TH1D* ahbPeaklr_low_om = hbackgroundPeak_low_om[i]->ProjectionY("ahbPeaklr_low_om",24,33);
            hsPeaklr_low_om[i] = hsignalPeak_low_om[i]->ProjectionY(Form("hsPeaklr_low_om%d",i),1,10);
            TH1D* ahsPeaklr_low_om = hsignalPeak_low_om[i]->ProjectionY("ahsPeaklr_low_om",24,33);

            hbPeaklr_low_om[i]->Add(ahbPeaklr_low_om);
            hsPeaklr_low_om[i]->Add(ahsPeaklr_low_om);
            hsPeaklr_low_om[i]->Divide(hbPeaklr_low_om[i]);
            hsPeaklr_low_om[i]->Scale(Bz_low_om[i]/nEvent_low_om[i]/BW2D);
            //hsPeaklr_low_om[i]->Scale(Bz_low_om[i]/BW2D);

            c_lr_low_om[j]->cd(i+1);

            TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.0,2.0);
            quadFit2->SetParameters(1,1,1);

            hsPeaklr_low_om[i]->Fit("quadFit2","R");
            //hsPeaklr_low_om[i]->Fit("quadFit2","R");
            //hsPeaklr_low_om[i]->Fit("quadFit2","R");

            double minVal_lr = quadFit2->GetMinimum(0.0,2.0);
            double minVal_lrX = quadFit2->GetMinimumX(0.0,2.0);
            TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_lr->SetParameter(0,-minVal_lr);
            TH1D* hsPeaklr_zeroed_low_om = (TH1D*)hsPeaklr_low_om[i]->Clone();
            hsPeaklr_zeroed_low_om->Add(minConst_lr);
            hsPeaklr_low_om[i]->Draw();
            xcoor = 0.52;
            ycoor = 0.90;
            increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_om[i],PD.PtBin_om[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|>2");
            Jyieldlr_low_om.push_back(hsPeaklr_zeroed_low_om->IntegralAndError(hsPeaklr_zeroed_low_om->FindBin(0.0),hsPeaklr_zeroed_low_om->FindBin(minVal_lrX),Jyieldlr_err_low_om[i],"width"));
            bin0yield = hsPeaklr_zeroed_low_om->GetBinContent(hsPeaklr_zeroed_low_om->FindBin(0.0))*0.19635;
            Jyieldlr_low_om[i] = Jyieldlr_low_om[i]*2 - bin0yield;

            JyieldSub_low_om.push_back(Jyieldsr_low_om[i] - Jyieldlr_low_om[i]);
            JyieldSub_err_low_om.push_back(sqrt(Jyieldsr_err_low_om[i]*Jyieldsr_err_low_om[i] + Jyieldlr_err_low_om[i]*Jyieldlr_err_low_om[i]));
            c->cd();

            //TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            //TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_low_om[i] = hsignalPeak_low_om[i]->ProjectionY(Form("V2lrs_low_om%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_low_om = hsignalPeak_low_om[i]->ProjectionY("aV2lrs_low_om",PD.binhigh_V2,33);
            V2lrb_low_om[i] = hbackgroundPeak_low_om[i]->ProjectionY(Form("V2lrb_low_om%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_low_om = hbackgroundPeak_low_om[i]->ProjectionY("aV2lrb_low_om",PD.binhigh_V2,33);
            V2lrs_low_om[i]->Add(aV2lrs_low_om);
            V2lrb_low_om[i]->Add(aV2lrb_low_om);
            V2lrs_low_om[i]->Divide(V2lrb_low_om[i]);
            V2lrs_low_om[i]->Scale(Bz_low_om[i]/nEvent_low_om[i]/BW2D);
            //V2lrs_low_om[i]->Scale(Bz_low_om[i]/BW2D);

            V2lrs_low_om[i]->Fit("fit1","R");
            V2lrs_low_om[i]->Fit("fit1","R");
            V2lrs_low_om[i]->Fit("fit1","R");

            V2lrs_low_om[i]->Draw();

            V2Values_low_om.push_back(fit1->GetParameter(2));
            V2Values_err_low_om.push_back(fit1->GetParError(2));
            V3Values_low_om.push_back(fit1->GetParameter(3));
            V3Values_err_low_om.push_back(fit1->GetParError(3));

            Nassoc_low_om.push_back(fit1->GetParameter(0));
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

            hbPeaksr_ks[i] = hbackgroundPeak_ks[i]->ProjectionY(("hb"+region_label+"sr_ks").c_str(), PD.sr_low, PD.sr_high);
            hsPeaksr_ks[i] = hsignalPeak_ks[i]->ProjectionY(("hs"+region_label+"sr_ks").c_str(), PD.sr_low, PD.sr_high);

            if(i<8)
            {
                c_high_2Dks_1[j]->cd(i+1);
                //TH2D* hPeak_high_2Dks = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_ks.c_str(),i));
                TH2D* hPeak_high_2Dks = (TH2D*)hsignalPeak_ks[i]->Clone();
                hPeak_high_2Dks->SetName(Form(("h"+region_label+"_high_2Dks%d").c_str(),i));
                hPeak_high_2Dks->Divide(hbackgroundPeak_ks[i]);

                hPeak_high_2Dks->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dks->Draw("SURF1");
                //hbackgroundPeak_high_ks[i]->Draw();
            }
            if(i>=8)
            {
                c_high_2Dks_2[j]->cd(i-8+1);
                //TH2D* hPeak_high_2Dks = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_ks.c_str(),i));
                TH2D* hPeak_high_2Dks = (TH2D*)hsignalPeak_ks[i]->Clone();
                hPeak_high_2Dks->SetName(Form(("h"+region_label+"_high_2Dks%d").c_str(),i));
                hPeak_high_2Dks->Divide(hbackgroundPeak_ks[i]);
                hPeak_high_2Dks->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dks->Draw("SURF1");
            }

            TF1* quadFit15 = new TF1("quadFit15","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit15->SetParameters(1,1,1);

            //nEvent_ks.push_back(mult_ks->Integral(0,100000));
            nEvent_ks.push_back(mult_ks->Integral(2,100000));
            Bz_ks.push_back(hbackgroundPeak_ks[i]->GetBinContent(hbackgroundPeak_ks[i]->FindBin(0,0)));

            hsPeaksr_ks[i]->Divide(hbPeaksr_ks[i]);
            hsPeaksr_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);
            //hsPeaksr_ks[i]->Scale(Bz_ks[i]/BW2D);

            c_sr_ks[j]->cd(i+1);

            hsPeaksr_ks[i]->Fit("quadFit15","R");
            //hsPeaksr_ks[i]->Fit("quadFit15","R");
            //hsPeaksr_ks[i]->Fit("quadFit15","R");

            double minVal_sr = quadFit15->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit15->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            TH1D* hsPeaksr_zeroed_ks = (TH1D*)hsPeaksr_ks[i]->Clone();
            TH1D* hsPeaksr_drawed_ks = (TH1D*)hsPeaksr_ks[i]->Clone();
            hsPeaksr_drawed_ks->Draw();
            hsPeaksr_zeroed_ks->Add(minConst_sr);
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_ks[i],PD.PtBin_ks[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_ks.push_back(hsPeaksr_zeroed_ks->IntegralAndError(hsPeaksr_zeroed_ks->FindBin(0.0),hsPeaksr_zeroed_ks->FindBin(minVal_srX),Jyieldsr_err_ks[i],"width"));
            double bin0yield = hsPeaksr_zeroed_ks->GetBinContent(hsPeaksr_zeroed_ks->FindBin(0.0))*0.19635;
            Jyieldsr_ks[i] = Jyieldsr_ks[i]*2 - bin0yield;

            hbPeaklr_ks[i] = hbackgroundPeak_ks[i]->ProjectionY("hbPeaklr_ks",1,10);
            TH1D* ahbPeaklr_ks = hbackgroundPeak_ks[i]->ProjectionY("ahbPeaklr_ks",24,33);
            hsPeaklr_ks[i] = hsignalPeak_ks[i]->ProjectionY("hsPeaklr_ks",1,10);
            TH1D* ahsPeaklr_ks = hsignalPeak_ks[i]->ProjectionY("ahsPeaklr_ks",24,33);

            hbPeaklr_ks[i]->Add(ahbPeaklr_ks);
            hsPeaklr_ks[i]->Add(ahsPeaklr_ks);
            hsPeaklr_ks[i]->Divide(hbPeaklr_ks[i]);
            hsPeaklr_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);
            //hsPeaklr_ks[i]->Scale(Bz_ks[i]/BW2D);

            c_lr_ks[j]->cd(i+1);

            TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.0,2.0);
            quadFit2->SetParameters(1,1,1);

            hsPeaklr_ks[i]->Fit("quadFit2","R");
            //hsPeaklr_ks[i]->Fit("quadFit2","R");
            //hsPeaklr_ks[i]->Fit("quadFit2","R");

            double minVal_lr = quadFit2->GetMinimum(0.0,2.0);
            double minVal_lrX = quadFit2->GetMinimumX(0.0,2.0);
            TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_lr->SetParameter(0,-minVal_lr);
            TH1D* hsPeaklr_zeroed_ks = (TH1D*)hsPeaklr_ks[i]->Clone();
            TH1D* hsPeaklr_drawed_ks = (TH1D*)hsPeaklr_ks[i]->Clone();
            hsPeaklr_drawed_ks->Draw();
            hsPeaklr_zeroed_ks->Add(minConst_lr);
            xcoor = 0.52;
            ycoor = 0.90;
            increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_ks[i],PD.PtBin_ks[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|>2");
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

            V2lrs_ks[i] = hsignalPeak_ks[i]->ProjectionY(Form("V2lrs_ks%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_ks = hsignalPeak_ks[i]->ProjectionY("aV2lrs_ks",PD.binhigh_V2,33);
            V2lrb_ks[i] = hbackgroundPeak_ks[i]->ProjectionY(Form("V2lrb_ks%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_ks = hbackgroundPeak_ks[i]->ProjectionY("aV2lrb_ks",PD.binhigh_V2,33);
            V2lrs_ks[i]->Add(aV2lrs_ks);
            V2lrb_ks[i]->Add(aV2lrb_ks);
            V2lrs_ks[i]->Divide(V2lrb_ks[i]);
            V2lrs_ks[i]->Scale(Bz_ks[i]/nEvent_ks[i]/BW2D);
            //V2lrs_ks[i]->Scale(Bz_ks[i]/BW2D);

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

            hbPeaksr_la[i] = hbackgroundPeak_la[i]->ProjectionY(("hb"+region_label+"sr_la").c_str(), PD.sr_low, PD.sr_high);
            hsPeaksr_la[i] = hsignalPeak_la[i]->ProjectionY(("hs"+region_label+"sr_la").c_str(), PD.sr_low, PD.sr_high);

            if(i<8)
            {
                c_high_2Dla_1[j]->cd(i+1);
                //TH2D* hPeak_high_2Dla = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_la.c_str(),i));
                TH2D* hPeak_high_2Dla = (TH2D*)hsignalPeak_la[i]->Clone();
                hPeak_high_2Dla->SetName(Form(("h"+region_label+"_high_2Dla%d").c_str(),i));
                hPeak_high_2Dla->Divide(hbackgroundPeak_la[i]);

                hPeak_high_2Dla->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dla->Draw("SURF1");
                //hbackgroundPeak_high_la[i]->Draw();
            }
            if(i>=8)
            {
                c_high_2Dla_2[j]->cd(i-8+1);
                //TH2D* hPeak_high_2Dla = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_la.c_str(),i));
                TH2D* hPeak_high_2Dla = (TH2D*)hsignalPeak_la[i]->Clone();
                hPeak_high_2Dla->SetName(Form(("h"+region_label+"_high_2Dla%d").c_str(),i));
                hPeak_high_2Dla->Divide(hbackgroundPeak_la[i]);
                hPeak_high_2Dla->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dla->Draw("SURF1");
            }


            TF1* quadFit16 = new TF1("quadFit16","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit16->SetParameters(1,1,1);

            //nEvent_la.push_back(mult_la->Integral(0,10000));
            nEvent_la.push_back(mult_la->Integral(2,10000));
            Bz_la.push_back(hbackgroundPeak_la[i]->GetBinContent(hbackgroundPeak_la[i]->FindBin(0,0)));

            hsPeaksr_la[i]->Divide(hbPeaksr_la[i]);
            hsPeaksr_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);
            //hsPeaksr_la[i]->Scale(Bz_la[i]/BW2D);

            c_sr_la[j]->cd(i+1);

            TH1D* hsPeaksr_zeroed_la = (TH1D*)hsPeaksr_la[i]->Clone();

            hsPeaksr_la[i]->Fit("quadFit16","R");
            //hsPeaksr_la[i]->Fit("quadFit16","R");
            //hsPeaksr_la[i]->Fit("quadFit16","R");
            TH1D* hsPeaksr_drawed_la = (TH1D*)hsPeaksr_la[i]->Clone();

            double minVal_sr = quadFit16->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit16->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            hsPeaksr_drawed_la->Draw();
            hsPeaksr_zeroed_la->Add(minConst_sr);
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%4.1f<p_{T}^{trg}<%4.1f",PD.PtBin_la[i],PD.PtBin_la[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_la.push_back(hsPeaksr_zeroed_la->IntegralAndError(hsPeaksr_zeroed_la->FindBin(0.0),hsPeaksr_zeroed_la->FindBin(minVal_srX),Jyieldsr_err_la[i],"width"));
            double bin0yield = hsPeaksr_zeroed_la->GetBinContent(hsPeaksr_zeroed_la->FindBin(0.0))*0.19635;
            Jyieldsr_la[i] = Jyieldsr_la[i]*2 - bin0yield;

            hbPeaklr_la[i] = hbackgroundPeak_la[i]->ProjectionY("hbPeaklr_la",1,10);
            TH1D* ahbPeaklr_la = hbackgroundPeak_la[i]->ProjectionY("ahbPeaklr_la",24,33);
            hsPeaklr_la[i] = hsignalPeak_la[i]->ProjectionY("hsPeaklr_la",1,10);
            TH1D* ahsPeaklr_la = hsignalPeak_la[i]->ProjectionY("ahsPeaklr_la",24,33);

            hbPeaklr_la[i]->Add(ahbPeaklr_la);
            hsPeaklr_la[i]->Add(ahsPeaklr_la);
            hsPeaklr_la[i]->Divide(hbPeaklr_la[i]);
            hsPeaklr_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);
            //hsPeaklr_la[i]->Scale(Bz_la[i]/BW2D);

            c_lr_la[j]->cd(i+1);

            TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.5,2.0);
            quadFit2->SetParameters(1,1,1);

            hsPeaklr_la[i]->Fit("quadFit2","R");
            //hsPeaklr_la[i]->Fit("quadFit2","R");
            //hsPeaklr_la[i]->Fit("quadFit2","R");

            double minVal_lr = quadFit2->GetMinimum(0.5,2.0);
            double minVal_lrX = quadFit2->GetMinimumX(0.5,2.0);
            TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_lr->SetParameter(0,-minVal_lr);
            TH1D* hsPeaklr_zeroed_la = (TH1D*)hsPeaklr_la[i]->Clone();
            TH1D* hsPeaklr_drawed_la = (TH1D*)hsPeaklr_la[i]->Clone();
            hsPeaklr_drawed_la->Draw();
            hsPeaklr_zeroed_la->Add(minConst_lr);
            xcoor = 0.52;
            ycoor = 0.90;
            increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"0<N_{trk}^{offline}35");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_la[i],PD.PtBin_la[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|>2");
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

            V2lrs_la[i] = hsignalPeak_la[i]->ProjectionY(Form("V2lrs_la%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_la = hsignalPeak_la[i]->ProjectionY("aV2lrs_la",PD.binhigh_V2,33);
            V2lrb_la[i] = hbackgroundPeak_la[i]->ProjectionY(Form("V2lrb_la%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_la = hbackgroundPeak_la[i]->ProjectionY("aV2lrb_la",PD.binhigh_V2,33);
            V2lrs_la[i]->Add(aV2lrs_la);
            V2lrb_la[i]->Add(aV2lrb_la);
            V2lrs_la[i]->Divide(V2lrb_la[i]);
            V2lrs_la[i]->Scale(Bz_la[i]/nEvent_la[i]/BW2D);
            //V2lrs_la[i]->Scale(Bz_la[i]/BW2D);

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
            //TH1D* mult_xi = (TH1D*) PD.f_Xi->Get((PD.fn_Xi + "/nEvtCut").c_str());
            TH1D* mult_xi = (TH1D*) PD.f_Xi->Get(Form(PathMult_high_xi.c_str(),i));

            hbPeaksr_xi[i] = hbackgroundPeak_xi[i]->ProjectionY(Form(("hb"+region_label+"sr_xi%d").c_str(),i), PD.sr_low, PD.sr_high);
            hsPeaksr_xi[i] = hsignalPeak_xi[i]->ProjectionY(Form(("hs"+region_label+"sr_xi%d").c_str(),i), PD.sr_low, PD.sr_high);

            if(i<8)
            {
                c_high_2Dxi_1[j]->cd(i+1);
                //TH2D* hPeak_high_2Dxi = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_xi.c_str(),i));
                TH2D* hPeak_high_2Dxi = (TH2D*)hsignalPeak_xi[i]->Clone();
                hPeak_high_2Dxi->SetName(Form(("h"+region_label+"_high_2Dxi%d").c_str(),i));
                hPeak_high_2Dxi->Divide(hbackgroundPeak_xi[i]);

                hPeak_high_2Dxi->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dxi->Draw("SURF1");
                //hbackgroundPeak_high_xi[i]->Draw();
            }
            if(i>=8)
            {
                c_high_2Dxi_2[j]->cd(i-8+1);
                //TH2D* hPeak_high_2Dxi = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_xi.c_str(),i));
                TH2D* hPeak_high_2Dxi = (TH2D*)hsignalPeak_xi[i]->Clone();
                hPeak_high_2Dxi->SetName(Form(("h"+region_label+"_high_2Dxi%d").c_str(),i));
                hPeak_high_2Dxi->Divide(hbackgroundPeak_xi[i]);
                hPeak_high_2Dxi->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dxi->Draw("SURF1");
            }

            TF1* quadFit16 = new TF1("quadFit16","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit16->SetParameters(1,1,1);

            //nEvent_xi.push_back(mult_xi->Integral(1,10000));
            nEvent_xi.push_back(mult_xi->Integral(2,10000));
            Bz_xi.push_back(hbackgroundPeak_xi[i]->GetBinContent(hbackgroundPeak_xi[i]->FindBin(0,0)));

            hsPeaksr_xi[i]->Divide(hbPeaksr_xi[i]);
            hsPeaksr_xi[i]->Scale(Bz_xi[i]/nEvent_xi[i]/BW2D);
            //hsPeaksr_xi[i]->Scale(Bz_xi[i]/BW2D);

            c_sr_xi[j]->cd(i+1);

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
            tex->DrawLatex(xcoor,ycoor-=increment,"185<N_{trk}^{offline}<250");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%4.1f<p_{T}^{trg}<%4.1f",PD.PtBin_xi[i],PD.PtBin_xi[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_xi.push_back(hsPeaksr_zeroed_xi->IntegralAndError(hsPeaksr_zeroed_xi->FindBin(0.0),hsPeaksr_zeroed_xi->FindBin(minVal_srX),Jyieldsr_err_xi[i],"width"));
            double bin0yield = hsPeaksr_zeroed_xi->GetBinContent(hsPeaksr_zeroed_xi->FindBin(0.0))*0.19635;
            Jyieldsr_xi[i] = Jyieldsr_xi[i]*2 - bin0yield;

            hbPeaklr_xi[i] = hbackgroundPeak_xi[i]->ProjectionY("hbPeaklr_xi",1,10);
            TH1D* ahbPeaklr_xi = hbackgroundPeak_xi[i]->ProjectionY("ahbPeaklr_xi",24,33);
            hsPeaklr_xi[i] = hsignalPeak_xi[i]->ProjectionY(Form("hsPeaklr_xi%d",i),1,10);
            TH1D* ahsPeaklr_xi = hsignalPeak_xi[i]->ProjectionY("ahsPeaklr_xi",24,33);

            hbPeaklr_xi[i]->Add(ahbPeaklr_xi);
            hsPeaklr_xi[i]->Add(ahsPeaklr_xi);
            hsPeaklr_xi[i]->Divide(hbPeaklr_xi[i]);
            hsPeaklr_xi[i]->Scale(Bz_xi[i]/nEvent_xi[i]/BW2D);
            //hsPeaklr_xi[i]->Scale(Bz_xi[i]/BW2D);

            c_lr_xi[j]->cd(i+1);

            TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.0,2.0);
            quadFit2->SetParameters(1,1,1);

            hsPeaklr_xi[i]->Fit("quadFit2","R");
            //hsPeaklr_xi[i]->Fit("quadFit2","R");
            //hsPeaklr_xi[i]->Fit("quadFit2","R");

            double minVal_lr = quadFit2->GetMinimum(0.0,2.0);
            double minVal_lrX = quadFit2->GetMinimumX(0.0,2.0);
            TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_lr->SetParameter(0,-minVal_lr);
            TH1D* hsPeaklr_zeroed_xi = (TH1D*)hsPeaklr_xi[i]->Clone(Form("hsPeaklr_zeroes_xi%d",i));
            hsPeaklr_zeroed_xi->Add(minConst_lr);
            hsPeaklr_xi[i]->Draw();
            xcoor = 0.52;
            ycoor = 0.90;
            increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"185 < N_{trk}^{offline} < 250");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_xi[i],PD.PtBin_xi[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|>2");
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

            c_lr_high_Fourier_xi[j]->cd(i+1);

            V2lrs_xi[i] = hsignalPeak_xi[i]->ProjectionY(Form("V2lrs_xi%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_xi = hsignalPeak_xi[i]->ProjectionY("aV2lrs_xi",PD.binhigh_V2,33);
            V2lrb_xi[i] = hbackgroundPeak_xi[i]->ProjectionY(Form("V2lrb_xi%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_xi = hbackgroundPeak_xi[i]->ProjectionY("aV2lrb_xi",PD.binhigh_V2,33);
            V2lrs_xi[i]->Add(aV2lrs_xi);
            V2lrb_xi[i]->Add(aV2lrb_xi);
            V2lrs_xi[i]->Divide(V2lrb_xi[i]);
            V2lrs_xi[i]->Scale(Bz_xi[i]/nEvent_xi[i]/BW2D);
            //V2lrs_xi[i]->Scale(Bz_xi[i]/BW2D);

            V2lrs_xi[i]->Fit("fit1","R");
            V2lrs_xi[i]->Fit("fit1","R");
            V2lrs_xi[i]->Fit("fit1","R");

            V2lrs_xi[i]->Draw();

            V2Values_xi.push_back(fit1->GetParameter(2));
            V2Values_err_xi.push_back(fit1->GetParError(2));
            V3Values_xi.push_back(fit1->GetParameter(3));
            V3Values_err_xi.push_back(fit1->GetParError(3));

            Nassoc_xi.push_back(fit1->GetParameter(0));
            c->cd();
        }

        //Om High N
        std::vector<double> Nassoc_om;
        std::vector<double> Jyieldsr_om;
        std::vector<double> Jyieldlr_om;
        std::vector<double> Jyieldsr_err_om(18);
        std::vector<double> Jyieldlr_err_om(18);
        std::vector<double> V2Values_om;
        std::vector<double> V2Values_err_om;
        std::vector<double> V3Values_om;
        std::vector<double> V3Values_err_om;
        std::vector<double> Bz_om;
        std::vector<double> nEvent_om;

        std::vector<double> JyieldSub_om;
        std::vector<double> JyieldSub_err_om;

        TH1D* hsPeaksr_om[arraySize_om];
        TH1D* hbPeaksr_om[arraySize_om];
        TH1D* hsPeaklr_om[arraySize_om];
        TH1D* hbPeaklr_om[arraySize_om];
        TH1D* V2lrs_om[arraySize_om];
        TH1D* V2lrb_om[arraySize_om];
        TH2D* hbackgroundPeak_om[arraySize_om];
        TH2D* hsignalPeak_om[arraySize_om];
        for(int i=0; i<numPtBins_om; i++)
        {
            PD.f_Om->GetObject(Form(PathBackground_high_om.c_str(),i),hbackgroundPeak_om[i]);
            PD.f_Om->GetObject(Form(PathSignal_high_om.c_str(),i),hsignalPeak_om[i]);
            TH1D* mult_om = (TH1D*) PD.f_Om->Get(Form(PathMult_high_om.c_str(),i));

            hbPeaksr_om[i] = hbackgroundPeak_om[i]->ProjectionY(Form(("hb"+region_label+"sr_om%d").c_str(),i), PD.sr_low, PD.sr_high);
            hsPeaksr_om[i] = hsignalPeak_om[i]->ProjectionY(Form(("hs"+region_label+"sr_om%d").c_str(),i), PD.sr_low, PD.sr_high);

            if(i<8)
            {
                c_high_2Dom_1[j]->cd(i+1);
                //TH2D* hPeak_high_2Dom = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_om.c_str(),i));
                TH2D* hPeak_high_2Dom = (TH2D*)hsignalPeak_om[i]->Clone();
                hPeak_high_2Dom->SetName(Form(("h"+region_label+"_high_2Dom%d").c_str(),i));
                hPeak_high_2Dom->Divide(hbackgroundPeak_om[i]);

                hPeak_high_2Dom->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dom->Draw("SURF1");
                //hbackgroundPeak_high_om[i]->Draw();
            }
            if(i>=8)
            {
                c_high_2Dom_2[j]->cd(i-8+1);
                //TH2D* hPeak_high_2Dom = (TH2D*)PD.f_V0->Get(Form(PathSignal_high_om.c_str(),i));
                TH2D* hPeak_high_2Dom = (TH2D*)hsignalPeak_om[i]->Clone();
                hPeak_high_2Dom->SetName(Form(("h"+region_label+"_high_2Dom%d").c_str(),i));
                hPeak_high_2Dom->Divide(hbackgroundPeak_om[i]);
                hPeak_high_2Dom->GetXaxis()->SetRangeUser(-3,3 );
                hPeak_high_2Dom->Draw("SURF1");
            }

            TF1* quadFit16 = new TF1("quadFit16","[0]*x^2+[1]*x+[2]",0.6,2.2);
            quadFit16->SetParameters(1,1,1);

            //nEvent_om.push_back(mult_om->Integral(1,10000));
            nEvent_om.push_back(mult_om->Integral(2,10000));
            Bz_om.push_back(hbackgroundPeak_om[i]->GetBinContent(hbackgroundPeak_om[i]->FindBin(0,0)));

            hsPeaksr_om[i]->Divide(hbPeaksr_om[i]);
            hsPeaksr_om[i]->Scale(Bz_om[i]/nEvent_om[i]/BW2D);

            c_sr_om[j]->cd(i+1);

            TH1D* hsPeaksr_zeroed_om = (TH1D*)hsPeaksr_om[i]->Clone();

            hsPeaksr_om[i]->Fit("quadFit16","R");
            //hsPeaksr_om[i]->Fit("quadFit16","R");
            //hsPeaksr_om[i]->Fit("quadFit16","R");

            double minVal_sr = quadFit16->GetMinimum(0.6,2.2);
            double minVal_srX = quadFit16->GetMinimumX(0.6,2.2);
            TF1* minConst_sr = new TF1("minConst_sr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_sr->SetParameter(0,-minVal_sr);
            hsPeaksr_zeroed_om->Add(minConst_sr);
            hsPeaksr_om[i]->Draw();
            double xcoor = 0.52;
            double ycoor = 0.90;
            double increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"185<N_{trk}^{offline}<250");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%4.1f<p_{T}^{trg}<%4.1f",PD.PtBin_om[i],PD.PtBin_om[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|<1");
            c->cd();
            Jyieldsr_om.push_back(hsPeaksr_zeroed_om->IntegralAndError(hsPeaksr_zeroed_om->FindBin(0.0),hsPeaksr_zeroed_om->FindBin(minVal_srX),Jyieldsr_err_om[i],"width"));
            double bin0yield = hsPeaksr_zeroed_om->GetBinContent(hsPeaksr_zeroed_om->FindBin(0.0))*0.19635;
            Jyieldsr_om[i] = Jyieldsr_om[i]*2 - bin0yield;

            hbPeaklr_om[i] = hbackgroundPeak_om[i]->ProjectionY("hbPeaklr_om",1,10);
            TH1D* ahbPeaklr_om = hbackgroundPeak_om[i]->ProjectionY("ahbPeaklr_om",24,33);
            hsPeaklr_om[i] = hsignalPeak_om[i]->ProjectionY(Form("hsPeaklr_om%d",i),1,10);
            TH1D* ahsPeaklr_om = hsignalPeak_om[i]->ProjectionY("ahsPeaklr_om",24,33);

            hbPeaklr_om[i]->Add(ahbPeaklr_om);
            hsPeaklr_om[i]->Add(ahsPeaklr_om);
            hsPeaklr_om[i]->Divide(hbPeaklr_om[i]);
            hsPeaklr_om[i]->Scale(Bz_om[i]/nEvent_om[i]/BW2D);

            c_lr_om[j]->cd(i+1);

            TF1* quadFit2 = new TF1("quadFit2","[0]*x^2+[1]*x+[2]",0.0,2.0);
            quadFit2->SetParameters(1,1,1);

            hsPeaklr_om[i]->Fit("quadFit2","R");
            //hsPeaklr_om[i]->Fit("quadFit2","R");
            //hsPeaklr_om[i]->Fit("quadFit2","R");

            double minVal_lr = quadFit2->GetMinimum(0.0,2.0);
            double minVal_lrX = quadFit2->GetMinimumX(0.0,2.0);
            TF1* minConst_lr = new TF1("minConst_lr","[0]",-(0.5-1.0/32)*3.1416,(1.5-1.0/32)*3.1416);
            minConst_lr->SetParameter(0,-minVal_lr);
            TH1D* hsPeaklr_zeroed_om = (TH1D*)hsPeaklr_om[i]->Clone(Form("hsPeaklr_zeroes_om%d",i));
            hsPeaklr_zeroed_om->Add(minConst_lr);
            hsPeaklr_om[i]->Draw();
            xcoor = 0.52;
            ycoor = 0.90;
            increment = 0.07;
            tex->DrawLatex(xcoor,ycoor-=increment,"CMS pPb 8.16 TeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"185 < N_{trk}^{offline} < 250");
            tex->DrawLatex(xcoor,ycoor-=increment,Form("%f<p_{T}^{trg}<%f",PD.PtBin_om[i],PD.PtBin_om[i+1]));
            tex->DrawLatex(xcoor,ycoor-=increment,"0.3<p_{T}^{assoc}<3GeV");
            tex->DrawLatex(xcoor,ycoor-=increment,"|#Delta#eta|>2");
            c->cd();
            Jyieldlr_om.push_back(hsPeaklr_zeroed_om->IntegralAndError(hsPeaklr_zeroed_om->FindBin(0.0),hsPeaklr_zeroed_om->FindBin(minVal_lrX),Jyieldlr_err_om[i],"width"));
            bin0yield = hsPeaklr_zeroed_om->GetBinContent(hsPeaklr_zeroed_om->FindBin(0.0))*0.19635;
            Jyieldlr_om[i] = Jyieldlr_om[i]*2 - bin0yield;

            JyieldSub_om.push_back(Jyieldsr_om[i] - Jyieldlr_om[i]);
            JyieldSub_err_om.push_back(sqrt(Jyieldsr_err_om[i]*Jyieldsr_err_om[i] + Jyieldlr_err_om[i]*Jyieldlr_err_om[i]));

            TF1* fit1 = new TF1("fit1","[0]*(1.0+2.0*[1]*cos(x)+2.0*[2]*cos(2.0*x)+2.0*[3]*cos(3.0*x)+2.0*[4]*cos(4.0*x))",-0.5*TMath::Pi(),1.5*TMath::Pi());
            fit1->SetParNames("N","V1","V2","V3");
            fit1->SetParameters(10,1,1,1);
            fit1->SetLineColor(2);

            V2lrs_om[i] = hsignalPeak_om[i]->ProjectionY(Form("V2lrs_om%d",i),1,PD.binlow_V2);
            TH1D* aV2lrs_om = hsignalPeak_om[i]->ProjectionY("aV2lrs_om",PD.binhigh_V2,33);
            V2lrb_om[i] = hbackgroundPeak_om[i]->ProjectionY(Form("V2lrb_om%d",i),1,PD.binlow_V2);
            TH1D* aV2lrb_om = hbackgroundPeak_om[i]->ProjectionY("aV2lrb_om",PD.binhigh_V2,33);
            V2lrs_om[i]->Add(aV2lrs_om);
            V2lrb_om[i]->Add(aV2lrb_om);
            V2lrs_om[i]->Divide(V2lrb_om[i]);
            V2lrs_om[i]->Scale(Bz_om[i]/nEvent_om[i]/BW2D);

            V2lrs_om[i]->Fit("fit1","R");
            V2lrs_om[i]->Fit("fit1","R");
            V2lrs_om[i]->Fit("fit1","R");

            V2lrs_om[i]->Draw();

            V2Values_om.push_back(fit1->GetParameter(2));
            V2Values_err_om.push_back(fit1->GetParError(2));
            V3Values_om.push_back(fit1->GetParameter(3));
            V3Values_err_om.push_back(fit1->GetParError(3));

            Nassoc_om.push_back(fit1->GetParameter(0));
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

        //Om Vn correction
        std::vector<double> V2sub_om;
        std::vector<double> V2sube_om;

        std::vector<double> V3sub_om;
        std::vector<double> V3sube_om;

        std::vector<double> jetYfactor_om;

        for(unsigned i=0; i<numPtBins_om; i++){
            V2sub_om.push_back(V2Values_om[i] - V2Values_low_om[i]*Nassoc_low_om[i]/Nassoc_om[i]*JyieldSub_om[i]/JyieldSub_low_om[i]);
            V2sube_om.push_back(sqrt(TMath::Power(V2Values_err_om[i],2) + sqrt(TMath::Power(V2Values_err_low_om[i]/V2Values_low_om[i],2) + TMath::Power(JyieldSub_err_om[i]/JyieldSub_om[i],2) + TMath::Power(JyieldSub_err_low_om[i]/JyieldSub_low_om[i],2))*V2Values_low_om[i]*Nassoc_low_om[i]/Nassoc_om[i]*JyieldSub_om[i]/JyieldSub_low_om[i]*sqrt(TMath::Power(V2Values_err_low_om[i]/V2Values_low_om[i],2) + TMath::Power(JyieldSub_err_om[i]/JyieldSub_om[i],2) + TMath::Power(JyieldSub_err_low_om[i]/JyieldSub_low_om[i],2))*V2Values_low_om[i]*Nassoc_low_om[i]/Nassoc_om[i]*JyieldSub_om[i]/JyieldSub_low_om[i]));

            V3sub_om.push_back(V3Values_om[i] - V3Values_low_om[i]*Nassoc_low_om[i]/Nassoc_om[i]*JyieldSub_om[i]/JyieldSub_low_om[i]);
            V3sube_om.push_back(sqrt(TMath::Power(V3Values_err_om[i],2) + sqrt(TMath::Power(V3Values_err_low_om[i]/V3Values_low_om[i],2) + TMath::Power(JyieldSub_err_om[i]/JyieldSub_om[i],2) + TMath::Power(JyieldSub_err_low_om[i]/JyieldSub_low_om[i],2))*V3Values_low_om[i]*Nassoc_low_om[i]/Nassoc_om[i]*JyieldSub_om[i]/JyieldSub_low_om[i]*sqrt(TMath::Power(V3Values_err_low_om[i]/V3Values_low_om[i],2) + TMath::Power(JyieldSub_err_om[i]/JyieldSub_om[i],2) + TMath::Power(JyieldSub_err_low_om[i]/JyieldSub_low_om[i],2))*V3Values_low_om[i]*Nassoc_low_om[i]/Nassoc_om[i]*JyieldSub_om[i]/JyieldSub_low_om[i]));

            jetYfactor_om.push_back(JyieldSub_om[i]/JyieldSub_low_om[i]);
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
            v3sube_la.push_back(fabs(sqrt(TMath::Power(V3sube_la[i]/V3sub_la[i],2) + TMath::Power(v3sube_ref/v3sub_ref,2))*v3sub_la[i]));

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
            v2sube_xi.push_back(fabs(sqrt(TMath::Power(V2sube_xi[i]/V2sub_xi[i],2) + TMath::Power(v2sube_ref/v2sub_ref,2))*v2sub_xi[i]));

            v2_xi.push_back(V2Values_xi[i]/v2_high_ref);
            v2e_xi.push_back(sqrt(TMath::Power(V2Values_err_xi[i]/V2Values_xi[i],2) + TMath::Power(v2e_high_ref/v2_high_ref,2))*v2_xi[i]);

            v2_low_xi.push_back(V2Values_low_xi[i]/v2_low_ref);
            v2e_low_xi.push_back(sqrt(TMath::Power(V2Values_err_low_xi[i]/V2Values_low_xi[i],2) + TMath::Power(v2e_low_ref/v2_low_ref,2))*v2_low_xi[i]);

            v3sub_xi.push_back(V3sub_xi[i]/v3sub_ref);
            v3sube_xi.push_back(fabs(sqrt(TMath::Power(V3sube_xi[i]/V3sub_xi[i],2) + TMath::Power(v3sube_ref/v3sub_ref,2))*v3sub_xi[i]));

            v3_xi.push_back(V3Values_xi[i]/v3_high_ref);
            v3e_xi.push_back(sqrt(TMath::Power(V3Values_err_xi[i]/V3Values_xi[i],2) + TMath::Power(v3e_high_ref/v3_high_ref,2))*v3_xi[i]);

            v3_low_xi.push_back(V3Values_low_xi[i]/v3_low_ref);
            v3e_low_xi.push_back(sqrt(TMath::Power(V3Values_err_low_xi[i]/V3Values_low_xi[i],2) + TMath::Power(v3e_low_ref/v3_low_ref,2))*v3_low_xi[i]);
        }

        //Om vn calculation
        std::vector<double> v2sub_om;
        std::vector<double> v2sube_om;
        std::vector<double> v2_om;
        std::vector<double> v2e_om;
        std::vector<double> v2_low_om;
        std::vector<double> v2e_low_om;

        std::vector<double> v3sub_om;
        std::vector<double> v3sube_om;
        std::vector<double> v3_om;
        std::vector<double> v3e_om;
        std::vector<double> v3_low_om;
        std::vector<double> v3e_low_om;

        for(unsigned i=0; i<numPtBins_om; i++){
            v2sub_om.push_back(V2sub_om[i]/v2sub_ref);
            v2sube_om.push_back(fabs(sqrt(TMath::Power(V2sube_om[i]/V2sub_om[i],2) + TMath::Power(v2sube_ref/v2sub_ref,2)))*v2sub_om[i]);

            v2_om.push_back(V2Values_om[i]/v2_high_ref);
            v2e_om.push_back(sqrt(TMath::Power(V2Values_err_om[i]/V2Values_om[i],2) + TMath::Power(v2e_high_ref/v2_high_ref,2))*v2_om[i]);

            v2_low_om.push_back(V2Values_low_om[i]/v2_low_ref);
            v2e_low_om.push_back(sqrt(TMath::Power(V2Values_err_low_om[i]/V2Values_low_om[i],2) + TMath::Power(v2e_low_ref/v2_low_ref,2))*v2_low_om[i]);

            v3sub_om.push_back(V3sub_om[i]/v3sub_ref);
            v3sube_om.push_back(fabs(sqrt(TMath::Power(V3sube_om[i]/V3sub_om[i],2) + TMath::Power(v3sube_ref/v3sub_ref,2))*v3sub_om[i]));

            v3_om.push_back(V3Values_om[i]/v3_high_ref);
            v3e_om.push_back(sqrt(TMath::Power(V3Values_err_om[i]/V3Values_om[i],2) + TMath::Power(v3e_high_ref/v3_high_ref,2))*v3_om[i]);

            v3_low_om.push_back(V3Values_low_om[i]/v3_low_ref);
            v3e_low_om.push_back(sqrt(TMath::Power(V3Values_err_low_om[i]/V3Values_low_om[i],2) + TMath::Power(v3e_low_ref/v3_low_ref,2))*v3_low_om[i]);
        }

        //Produce RootFiles
        std::string command = "";
        if(j==0) command = "RECREATE";
        else command = "UPDATE";
        TFile output(PD.fn.c_str(),command.c_str());

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

        std::vector<double> pt_om;
        std::vector<double> Ket_om;
        for(unsigned i=0; i<numPtBins_om; i++){
            TH1D* hpt = (TH1D*)PD.f_Om->Get(Form((PD.fn_Om + "/Pt_om_pt%d").c_str(),i));
            TH1D* hket = (TH1D*)PD.f_Om->Get(Form((PD.fn_Om + "/KET_om_pt%d").c_str(),i));

            pt_om.push_back(hpt->GetMean(1));
            Ket_om.push_back(hket->GetMean(1));
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

        std::vector<double> perisubfactor_om;
        std::vector<double> perisubfactore_om;
        for(int i=0; i<numPtBins_om; i++){
            perisubfactor_om.push_back(V2Values_low_om[i]*Nassoc_low_om[i]/JyieldSub_low_om[i]);
            perisubfactore_om.push_back(sqrt(TMath::Power(V2Values_err_low_om[i]/V2Values_low_om[i],2) + TMath::Power(JyieldSub_err_low_om[i]/JyieldSub_low_om[i],2))*perisubfactor_om[i]);
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
        TGraphErrors* v2plot_low_xi    = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&v2_low_xi[0]      ,0,&v2e_low_xi[0]);
        TGraphErrors* v2plot_KET_xi    = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&v2_xi[0]      ,0,&v2e_xi[0]);
        TGraphErrors* v2subplot_xi     = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&v2sub_xi[0]   ,0,&v2sube_xi[0]);
        TGraphErrors* v2subplot_KET_xi = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&v2sub_xi[0]   ,0,&v2sube_xi[0]);
        TGraphErrors* v3plot_xi        = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&v3_xi[0]      ,0,&v3e_xi[0]);
        TGraphErrors* v3plot_KET_xi    = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&v3_xi[0]      ,0,&v3e_xi[0]);
        TGraphErrors* v3subplot_xi     = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&v3sub_xi[0]   ,0,&v3sube_xi[0]);
        TGraphErrors* v3subplot_KET_xi = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&v3sub_xi[0]   ,0,&v3sube_xi[0]);

        TGraphErrors* v2plot_om        = new TGraphErrors(numPtBins_om,&pt_om[0] ,&v2_om[0]      ,0,&v2e_om[0]);
        TGraphErrors* v2plot_low_om    = new TGraphErrors(numPtBins_om,&pt_om[0] ,&v2_low_om[0]      ,0,&v2e_low_om[0]);
        TGraphErrors* v2plot_KET_om    = new TGraphErrors(numPtBins_om,&Ket_om[0],&v2_om[0]      ,0,&v2e_om[0]);
        TGraphErrors* v2subplot_om     = new TGraphErrors(numPtBins_om,&pt_om[0] ,&v2sub_om[0]   ,0,&v2sube_om[0]);
        TGraphErrors* v2subplot_KET_om = new TGraphErrors(numPtBins_om,&Ket_om[0],&v2sub_om[0]   ,0,&v2sube_om[0]);
        TGraphErrors* v3plot_om        = new TGraphErrors(numPtBins_om,&pt_om[0] ,&v3_om[0]      ,0,&v3e_om[0]);
        TGraphErrors* v3plot_KET_om    = new TGraphErrors(numPtBins_om,&Ket_om[0],&v3_om[0]      ,0,&v3e_om[0]);
        TGraphErrors* v3subplot_om     = new TGraphErrors(numPtBins_om,&pt_om[0] ,&v3sub_om[0]   ,0,&v3sube_om[0]);
        TGraphErrors* v3subplot_KET_om = new TGraphErrors(numPtBins_om,&Ket_om[0],&v3sub_om[0]   ,0,&v3sube_om[0]);

        //Obs V2 values
        TGraphErrors* V2plot_ks        = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&V2Values_ks[0]    ,0,&V2Values_err_ks[0]);
        TGraphErrors* V2plot_low_ks    = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&V2Values_low_ks[0],0,&V2Values_err_low_ks[0]);
        TGraphErrors* V2plot_KET_ks    = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&V2Values_ks[0]    ,0,&V2Values_err_ks[0]);
        TGraphErrors* V2subplot_ks     = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&V2sub_ks[0]       ,0,&V2sube_ks[0]);
        TGraphErrors* V2subplot_KET_ks = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&V2sub_ks[0]       ,0,&V2sube_ks[0]);
        TGraphErrors* V3plot_ks        = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&V3Values_ks[0]    ,0,&V3Values_err_ks[0]);
        TGraphErrors* V3plot_KET_ks    = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&V3Values_ks[0]    ,0,&V3Values_err_ks[0]);
        TGraphErrors* V3subplot_ks     = new TGraphErrors(numPtBins_ks,&pt_ks[0] ,&V3sub_ks[0]       ,0,&V3sube_ks[0]);
        TGraphErrors* V3subplot_KET_ks = new TGraphErrors(numPtBins_ks,&Ket_ks[0],&V3sub_ks[0]       ,0,&V3sube_ks[0]);

        TGraphErrors* V2plot_la        = new TGraphErrors(numPtBins_la,&pt_la[0] ,&V2Values_la[0]    ,0,&V2Values_err_la[0]);
        TGraphErrors* V2plot_low_la    = new TGraphErrors(numPtBins_la,&pt_la[0] ,&V2Values_low_la[0],0,&V2Values_err_low_la[0]);
        TGraphErrors* V2plot_KET_la    = new TGraphErrors(numPtBins_la,&Ket_la[0],&V2Values_la[0]    ,0,&V2Values_err_la[0]);
        TGraphErrors* V2subplot_la     = new TGraphErrors(numPtBins_la,&pt_la[0] ,&V2sub_la[0]       ,0,&V2sube_la[0]);
        TGraphErrors* V2subplot_KET_la = new TGraphErrors(numPtBins_la,&Ket_la[0],&V2sub_la[0]       ,0,&V2sube_la[0]);
        TGraphErrors* V3plot_la        = new TGraphErrors(numPtBins_la,&pt_la[0] ,&V3Values_la[0]    ,0,&V3Values_err_la[0]);
        TGraphErrors* V3plot_KET_la    = new TGraphErrors(numPtBins_la,&Ket_la[0],&V3Values_la[0]    ,0,&V3Values_err_la[0]);
        TGraphErrors* V3subplot_la     = new TGraphErrors(numPtBins_la,&pt_la[0] ,&V3sub_la[0]       ,0,&V3sube_la[0]);
        TGraphErrors* V3subplot_KET_la = new TGraphErrors(numPtBins_la,&Ket_la[0],&V3sub_la[0]       ,0,&V3sube_la[0]);

        TGraphErrors* V2plot_xi        = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&V2Values_xi[0]    ,0,&V2Values_err_xi[0]);
        TGraphErrors* V2plot_low_xi    = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&V2Values_low_xi[0],0,&V2Values_err_low_xi[0]);
        TGraphErrors* V2plot_KET_xi    = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&V2Values_xi[0]    ,0,&V2Values_err_xi[0]);
        TGraphErrors* V2subplot_xi     = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&V2sub_xi[0]       ,0,&V2sube_xi[0]);
        TGraphErrors* V2subplot_KET_xi = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&V2sub_xi[0]       ,0,&V2sube_xi[0]);
        TGraphErrors* V3plot_xi        = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&V3Values_xi[0]    ,0,&V3Values_err_xi[0]);
        TGraphErrors* V3plot_KET_xi    = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&V3Values_xi[0]    ,0,&V3Values_err_xi[0]);
        TGraphErrors* V3subplot_xi     = new TGraphErrors(numPtBins_xi,&pt_xi[0] ,&V3sub_xi[0]       ,0,&V3sube_xi[0]);
        TGraphErrors* V3subplot_KET_xi = new TGraphErrors(numPtBins_xi,&Ket_xi[0],&V3sub_xi[0]       ,0,&V3sube_xi[0]);

        TGraphErrors* V2plot_om        = new TGraphErrors(numPtBins_om,&pt_om[0] ,&V2Values_om[0]    ,0,&V2Values_err_om[0]);
        TGraphErrors* V2plot_low_om    = new TGraphErrors(numPtBins_om,&pt_om[0] ,&V2Values_low_om[0],0,&V2Values_err_low_om[0]);
        TGraphErrors* V2plot_KET_om    = new TGraphErrors(numPtBins_om,&Ket_om[0],&V2Values_om[0]    ,0,&V2Values_err_om[0]);
        TGraphErrors* V2subplot_om     = new TGraphErrors(numPtBins_om,&pt_om[0] ,&V2sub_om[0]       ,0,&V2sube_om[0]);
        TGraphErrors* V2subplot_KET_om = new TGraphErrors(numPtBins_om,&Ket_om[0],&V2sub_om[0]       ,0,&V2sube_om[0]);
        TGraphErrors* V3plot_om        = new TGraphErrors(numPtBins_om,&pt_om[0] ,&V3Values_om[0]    ,0,&V3Values_err_om[0]);
        TGraphErrors* V3plot_KET_om    = new TGraphErrors(numPtBins_om,&Ket_om[0],&V3Values_om[0]    ,0,&V3Values_err_om[0]);
        TGraphErrors* V3subplot_om     = new TGraphErrors(numPtBins_om,&pt_om[0] ,&V3sub_om[0]       ,0,&V3sube_om[0]);
        TGraphErrors* V3subplot_KET_om = new TGraphErrors(numPtBins_om,&Ket_om[0],&V3sub_om[0]       ,0,&V3sube_om[0]);

        //Yields
        TGraphErrors* srYieldPlot_ks       = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Jyieldsr_ks[0]     ,0,&Jyieldsr_err_ks[0]);
        TGraphErrors* lrYieldPlot_ks       = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Jyieldlr_ks[0]     ,0,&Jyieldlr_err_ks[0]);
        TGraphErrors* srYieldPlot_low_ks   = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Jyieldsr_low_ks[0]     ,0,&Jyieldsr_err_low_ks[0]);
        TGraphErrors* lrYieldPlot_low_ks   = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Jyieldlr_low_ks[0]     ,0,&Jyieldlr_err_low_ks[0]);
        TGraphErrors* subYieldPlot_ks      = new TGraphErrors(numPtBins_ks,&pt_ks[0],&JyieldSub_ks[0]    ,0,&JyieldSub_err_ks[0]);
        TGraphErrors* perisubfactorplot_ks = new TGraphErrors(numPtBins_ks,&pt_ks[0],&perisubfactor_ks[0],0,&perisubfactore_ks[0]);

        TGraphErrors* srYieldPlot_la       = new TGraphErrors(numPtBins_la,&pt_la[0],&Jyieldsr_la[0]     ,0,&Jyieldsr_err_la[0]);
        TGraphErrors* lrYieldPlot_la       = new TGraphErrors(numPtBins_la,&pt_la[0],&Jyieldlr_la[0]     ,0,&Jyieldlr_err_la[0]);
        TGraphErrors* srYieldPlot_low_la       = new TGraphErrors(numPtBins_la,&pt_la[0],&Jyieldsr_low_la[0]     ,0,&Jyieldsr_err_low_la[0]);
        TGraphErrors* lrYieldPlot_low_la       = new TGraphErrors(numPtBins_la,&pt_la[0],&Jyieldlr_low_la[0]     ,0,&Jyieldlr_err_low_la[0]);
        TGraphErrors* subYieldPlot_la      = new TGraphErrors(numPtBins_la,&pt_la[0],&JyieldSub_la[0]    ,0,&JyieldSub_err_la[0]);
        TGraphErrors* perisubfactorplot_la = new TGraphErrors(numPtBins_la,&pt_la[0],&perisubfactor_la[0],0,&perisubfactore_la[0]);

        TGraphErrors* srYieldPlot_xi       = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Jyieldsr_xi[0]     ,0,&Jyieldsr_err_xi[0]);
        TGraphErrors* lrYieldPlot_xi       = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Jyieldlr_xi[0]     ,0,&Jyieldlr_err_xi[0]);
        TGraphErrors* srYieldPlot_low_xi       = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Jyieldsr_low_xi[0]     ,0,&Jyieldsr_err_low_xi[0]);
        TGraphErrors* lrYieldPlot_low_xi       = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Jyieldlr_low_xi[0]     ,0,&Jyieldlr_err_low_xi[0]);
        TGraphErrors* subYieldPlot_xi      = new TGraphErrors(numPtBins_xi,&pt_xi[0],&JyieldSub_xi[0]    ,0,&JyieldSub_err_xi[0]);
        TGraphErrors* perisubfactorplot_xi = new TGraphErrors(numPtBins_xi,&pt_xi[0],&perisubfactor_xi[0],0,&perisubfactore_xi[0]);

        TGraphErrors* srYieldPlot_om       = new TGraphErrors(numPtBins_om,&pt_om[0],&Jyieldsr_om[0]     ,0,&Jyieldsr_err_om[0]);
        TGraphErrors* lrYieldPlot_om       = new TGraphErrors(numPtBins_om,&pt_om[0],&Jyieldlr_om[0]     ,0,&Jyieldlr_err_om[0]);
        TGraphErrors* srYieldPlot_low_om       = new TGraphErrors(numPtBins_om,&pt_om[0],&Jyieldsr_low_om[0]     ,0,&Jyieldsr_err_low_om[0]);
        TGraphErrors* lrYieldPlot_low_om       = new TGraphErrors(numPtBins_om,&pt_om[0],&Jyieldlr_low_om[0]     ,0,&Jyieldlr_err_low_om[0]);
        TGraphErrors* subYieldPlot_om      = new TGraphErrors(numPtBins_om,&pt_om[0],&JyieldSub_om[0]    ,0,&JyieldSub_err_om[0]);
        TGraphErrors* perisubfactorplot_om = new TGraphErrors(numPtBins_om,&pt_om[0],&perisubfactor_om[0],0,&perisubfactore_om[0]);

        //For Direct subtraction method
        TGraphErrors* Nasslow_ks         = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Nassoc_low_ks[0]   ,0,0);
        TGraphErrors* Nass_ks            = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Nassoc_ks[0]       ,0,0);
        TGraphErrors* subYieldPlotLow_ks = new TGraphErrors(numPtBins_ks,&pt_ks[0],&JyieldSub_low_ks[0],0,&JyieldSub_err_low_ks[0]);

        TGraphErrors* Nasslow_la         = new TGraphErrors(numPtBins_la,&pt_la[0],&Nassoc_low_la[0]   ,0,0);
        TGraphErrors* Nass_la            = new TGraphErrors(numPtBins_la,&pt_la[0],&Nassoc_la[0]       ,0,0);
        TGraphErrors* subYieldPlotLow_la = new TGraphErrors(numPtBins_la,&pt_la[0],&JyieldSub_low_la[0],0,&JyieldSub_err_low_la[0]);

        TGraphErrors* Nasslow_xi         = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Nassoc_low_xi[0]   ,0,0);
        TGraphErrors* Nass_xi            = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Nassoc_xi[0]       ,0,0);
        TGraphErrors* subYieldPlotLow_xi = new TGraphErrors(numPtBins_xi,&pt_xi[0],&JyieldSub_low_xi[0],0,&JyieldSub_err_low_xi[0]);

        TGraphErrors* Nasslow_om         = new TGraphErrors(numPtBins_om,&pt_om[0],&Nassoc_low_om[0]   ,0,0);
        TGraphErrors* Nass_om            = new TGraphErrors(numPtBins_om,&pt_om[0],&Nassoc_om[0]       ,0,0);
        TGraphErrors* subYieldPlotLow_om = new TGraphErrors(numPtBins_om,&pt_om[0],&JyieldSub_low_om[0],0,&JyieldSub_err_low_om[0]);

        TGraphErrors* bz_ks = new TGraphErrors(numPtBins_ks,&pt_ks[0],&Bz_ks[0],0,0);
        TGraphErrors* bz_la = new TGraphErrors(numPtBins_la,&pt_la[0],&Bz_la[0],0,0);
        TGraphErrors* bz_xi = new TGraphErrors(numPtBins_xi,&pt_xi[0],&Bz_xi[0],0,0);
        TGraphErrors* bz_om = new TGraphErrors(numPtBins_om,&pt_om[0],&Bz_om[0],0,0);

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
        v2plot_low_xi       ->Write(("v2plot_" + region + "_low_xi"       ).c_str(),TObject::kOverwrite);
        v2plot_KET_xi   ->Write(("v2plot_KET_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        v2subplot_xi    ->Write(("v2subplot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);
        v2subplot_KET_xi->Write(("v2subplot_KET_" + region + "_xi").c_str(),TObject::kOverwrite);
        v3plot_xi       ->Write(("v3plot_" + region + "_xi"       ).c_str(),TObject::kOverwrite);
        v3plot_KET_xi   ->Write(("v3plot_KET_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        v3subplot_xi    ->Write(("v3subplot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);
        v3subplot_KET_xi->Write(("v3subplot_KET_" + region + "_xi").c_str(),TObject::kOverwrite);

        v2plot_om       ->Write(("v2plot_" + region + "_om"       ).c_str(),TObject::kOverwrite);
        v2plot_low_om       ->Write(("v2plot_" + region + "_low_om"       ).c_str(),TObject::kOverwrite);
        v2plot_KET_om   ->Write(("v2plot_KET_" + region + "_om"   ).c_str(),TObject::kOverwrite);
        v2subplot_om    ->Write(("v2subplot_" + region + "_om"    ).c_str(),TObject::kOverwrite);
        v2subplot_KET_om->Write(("v2subplot_KET_" + region + "_om").c_str(),TObject::kOverwrite);
        v3plot_om       ->Write(("v3plot_" + region + "_om"       ).c_str(),TObject::kOverwrite);
        v3plot_KET_om   ->Write(("v3plot_KET_" + region + "_om"   ).c_str(),TObject::kOverwrite);
        v3subplot_om    ->Write(("v3subplot_" + region + "_om"    ).c_str(),TObject::kOverwrite);
        v3subplot_KET_om->Write(("v3subplot_KET_" + region + "_om").c_str(),TObject::kOverwrite);

        V2plot_ks        -> Write(("V2plot_" + region + "_ks"       ).c_str(),TObject::kOverwrite);
        V2plot_low_ks    -> Write(("V2plot_" + region + "_low_ks"   ).c_str(),TObject::kOverwrite);
        V2plot_KET_ks    -> Write(("V2plot_KET_" + region + "_ks"   ).c_str(),TObject::kOverwrite);
        V2subplot_ks     -> Write(("V2subplot_" + region + "_ks"    ).c_str(),TObject::kOverwrite);
        V2subplot_KET_ks -> Write(("V2subplot_KET_" + region + "_ks").c_str(),TObject::kOverwrite);
        V3plot_ks        -> Write(("V3plot_" + region + "_ks"       ).c_str(),TObject::kOverwrite);
        V3plot_KET_ks    -> Write(("V3plot_KET_" + region + "_ks"   ).c_str(),TObject::kOverwrite);
        V3subplot_ks     -> Write(("V3subplot_" + region + "_ks"    ).c_str(),TObject::kOverwrite);
        V3subplot_KET_ks -> Write(("V3subplot_KET_" + region + "_ks").c_str(),TObject::kOverwrite);

        V2plot_la        -> Write(("V2plot_" + region + "_la"       ).c_str(),TObject::kOverwrite);
        V2plot_low_la    -> Write(("V2plot_" + region + "_low_la"   ).c_str(),TObject::kOverwrite);
        V2plot_KET_la    -> Write(("V2plot_KET_" + region + "_la"   ).c_str(),TObject::kOverwrite);
        V2subplot_la     -> Write(("V2subplot_" + region + "_la"    ).c_str(),TObject::kOverwrite);
        V2subplot_KET_la -> Write(("V2subplot_KET_" + region + "_la").c_str(),TObject::kOverwrite);
        V3plot_la        -> Write(("V3plot_" + region + "_la"       ).c_str(),TObject::kOverwrite);
        V3plot_KET_la    -> Write(("V3plot_KET_" + region + "_la"   ).c_str(),TObject::kOverwrite);
        V3subplot_la     -> Write(("V3subplot_" + region + "_la"    ).c_str(),TObject::kOverwrite);
        V3subplot_KET_la -> Write(("V3subplot_KET_" + region + "_la").c_str(),TObject::kOverwrite);

        V2plot_xi        -> Write(("V2plot_" + region + "_xi"       ).c_str(),TObject::kOverwrite);
        V2plot_low_xi    -> Write(("V2plot_" + region + "_low_xi"   ).c_str(),TObject::kOverwrite);
        V2plot_KET_xi    -> Write(("V2plot_KET_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        V2subplot_xi     -> Write(("V2subplot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);
        V2subplot_KET_xi -> Write(("V2subplot_KET_" + region + "_xi").c_str(),TObject::kOverwrite);
        V3plot_xi        -> Write(("V3plot_" + region + "_xi"       ).c_str(),TObject::kOverwrite);
        V3plot_KET_xi    -> Write(("V3plot_KET_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        V3subplot_xi     -> Write(("V3subplot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);
        V3subplot_KET_xi -> Write(("V3subplot_KET_" + region + "_xi").c_str(),TObject::kOverwrite);

        V2plot_om        -> Write(("V2plot_" + region + "_om"       ).c_str(),TObject::kOverwrite);
        V2plot_low_om    -> Write(("V2plot_" + region + "_low_om"   ).c_str(),TObject::kOverwrite);
        V2plot_KET_om    -> Write(("V2plot_KET_" + region + "_om"   ).c_str(),TObject::kOverwrite);
        V2subplot_om     -> Write(("V2subplot_" + region + "_om"    ).c_str(),TObject::kOverwrite);
        V2subplot_KET_om -> Write(("V2subplot_KET_" + region + "_om").c_str(),TObject::kOverwrite);
        V3plot_om        -> Write(("V3plot_" + region + "_om"       ).c_str(),TObject::kOverwrite);
        V3plot_KET_om    -> Write(("V3plot_KET_" + region + "_om"   ).c_str(),TObject::kOverwrite);
        V3subplot_om     -> Write(("V3subplot_" + region + "_om"    ).c_str(),TObject::kOverwrite);
        V3subplot_KET_om -> Write(("V3subplot_KET_" + region + "_om").c_str(),TObject::kOverwrite);

        perisubfactorplot_ks->Write(("perisubfactor_" + region + "_ks").c_str(),TObject::kOverwrite);
        perisubfactorplot_la->Write(("perisubfactor_" + region + "_la").c_str(),TObject::kOverwrite);
        perisubfactorplot_xi->Write(("perisubfactor_" + region + "_xi").c_str(),TObject::kOverwrite);
        perisubfactorplot_om->Write(("perisubfactor_" + region + "_om").c_str(),TObject::kOverwrite);

        Nass_ks           ->Write(("Nassoc_" + region + "_ks"       ).c_str(),TObject::kOverwrite);
        Nasslow_ks        ->Write(("Nassoc_low_" + region + "_ks"   ).c_str(),TObject::kOverwrite);
        subYieldPlotLow_ks->Write(("YieldPlot_low_" + region + "_ks").c_str(),TObject::kOverwrite);
        subYieldPlot_ks   ->Write(("YieldPlot_" + region + "_ks"    ).c_str(),TObject::kOverwrite);
        srYieldPlot_ks->Write(("YieldsrPlot_" + region + "_ks").c_str(),TObject::kOverwrite);
        lrYieldPlot_ks->Write(("YieldlrPlot_" + region + "_ks").c_str(),TObject::kOverwrite);
        srYieldPlot_low_ks->Write(("YieldsrPlot_" + region + "low_ks").c_str(),TObject::kOverwrite);
        lrYieldPlot_low_ks->Write(("YieldlrPlot_" + region + "low_ks").c_str(),TObject::kOverwrite);

        Nass_la           ->Write(("Nassoc_" + region + "_la"       ).c_str(),TObject::kOverwrite);
        Nasslow_la        ->Write(("Nassoc_low_" + region + "_la"   ).c_str(),TObject::kOverwrite);
        subYieldPlotLow_la->Write(("YieldPlot_low_" + region + "_la").c_str(),TObject::kOverwrite);
        subYieldPlot_la   ->Write(("YieldPlot_" + region + "_la"    ).c_str(),TObject::kOverwrite);
        srYieldPlot_la->Write(("YieldsrPlot_" + region + "_la").c_str(),TObject::kOverwrite);
        lrYieldPlot_la->Write(("YieldlrPlot_" + region + "_la").c_str(),TObject::kOverwrite);
        srYieldPlot_low_la->Write(("YieldsrPlot_" + region + "low_la").c_str(),TObject::kOverwrite);
        lrYieldPlot_low_la->Write(("YieldlrPlot_" + region + "low_la").c_str(),TObject::kOverwrite);

        Nass_xi           ->Write(("Nassoc_" + region + "_xi"       ).c_str(),TObject::kOverwrite);
        Nasslow_xi        ->Write(("Nassoc_low_" + region + "_xi"   ).c_str(),TObject::kOverwrite);
        subYieldPlotLow_xi->Write(("YieldPlot_low_" + region + "_xi").c_str(),TObject::kOverwrite);
        subYieldPlot_xi   ->Write(("YieldPlot_" + region + "_xi"    ).c_str(),TObject::kOverwrite);
        srYieldPlot_xi->Write(("YieldsrPlot_" + region + "_xi").c_str(),TObject::kOverwrite);
        lrYieldPlot_xi->Write(("YieldlrPlot_" + region + "_xi").c_str(),TObject::kOverwrite);
        srYieldPlot_low_xi->Write(("YieldsrPlot_" + region + "low_xi").c_str(),TObject::kOverwrite);
        lrYieldPlot_low_xi->Write(("YieldlrPlot_" + region + "low_xi").c_str(),TObject::kOverwrite);

        Nass_om           ->Write(("Nassoc_" + region + "_om"       ).c_str(),TObject::kOverwrite);
        Nasslow_om        ->Write(("Nassoc_low_" + region + "_om"   ).c_str(),TObject::kOverwrite);
        subYieldPlotLow_om->Write(("YieldPlot_low_" + region + "_om").c_str(),TObject::kOverwrite);
        subYieldPlot_om   ->Write(("YieldPlot_" + region + "_om"    ).c_str(),TObject::kOverwrite);
        srYieldPlot_om->Write(("YieldsrPlot_" + region + "_om").c_str(),TObject::kOverwrite);
        lrYieldPlot_om->Write(("YieldlrPlot_" + region + "_om").c_str(),TObject::kOverwrite);
        srYieldPlot_low_om->Write(("YieldsrPlot_" + region + "low_om").c_str(),TObject::kOverwrite);
        lrYieldPlot_low_om->Write(("YieldlrPlot_" + region + "low_om").c_str(),TObject::kOverwrite);

        bz_ks->Write(("bz_" + region + "_ks").c_str(),TObject::kOverwrite);
        bz_la->Write(("bz_" + region + "_la").c_str(),TObject::kOverwrite);
        bz_xi->Write(("bz_" + region + "_xi").c_str(),TObject::kOverwrite);
        bz_om->Write(("bz_" + region + "_om").c_str(),TObject::kOverwrite);

        c_lr_low_Fourier_xi[j]->Print(("Image/PeriSub/Fourier_" + region + "_low_xiEG2_2terms.pdf").c_str());
        c_sr_low_xi[j]->Print(("Image/PeriSub/Yield_sr_" + region + "_low_xi.pdf").c_str());
        c_lr_low_xi[j]->Print(("Image/PeriSub/Yield_lr_" + region + "_low_xi.pdf").c_str());

        //c_lr_low_Fourier_om[j]->Print(("Image/PeriSub/Fourier_" + region + "_low_omEG2_2terms.pdf").c_str());
        //c_sr_low_om[j]->Print(("Image/PeriSub/Yield_sr_" + region + "_low_om.pdf").c_str());
        //c_lr_low_om[j]->Print(("Image/PeriSub/Yield_lr_" + region + "_low_om.pdf").c_str());

        c_sr_low_ks[j] ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_sr_low_ks[j]->GetTitle() + ".pdf").c_str() );
        c_lr_low_ks[j] ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_lr_low_ks[j]->GetTitle() + ".pdf").c_str() );
        c_sr_ks[j]     ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_sr_ks[j]->GetTitle()     + ".pdf").c_str() );
        c_lr_ks[j]     ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_lr_ks[j]->GetTitle()     + ".pdf").c_str() );
        c_sr_low_la[j] ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_sr_low_la[j]->GetTitle() + ".pdf").c_str() );
        c_lr_low_la[j] ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_lr_low_la[j]->GetTitle() + ".pdf").c_str() );
        c_sr_la[j]     ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_sr_la[j]->GetTitle()     + ".pdf").c_str() );
        c_lr_la[j]     ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_lr_la[j]->GetTitle()     + ".pdf").c_str() );
        c_sr_low_xi[j] ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_sr_low_xi[j]->GetTitle() + ".pdf").c_str() );
        c_lr_low_xi[j] ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_lr_low_xi[j]->GetTitle() + ".pdf").c_str() );
        c_sr_xi[j]     ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_sr_xi[j]->GetTitle()     + ".pdf").c_str() );
        c_lr_xi[j]     ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_lr_xi[j]->GetTitle()     + ".pdf").c_str() );
        c_sr_low_om[j] ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_sr_low_om[j]->GetTitle() + ".pdf").c_str() );
        c_lr_low_om[j] ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_lr_low_om[j]->GetTitle() + ".pdf").c_str() );
        c_sr_om[j]     ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_sr_om[j]->GetTitle()     + ".pdf").c_str() );
        c_lr_om[j]     ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_lr_om[j]->GetTitle()     + ".pdf").c_str() );
    }
    c_sr_low_ref ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_sr_low_ref ->GetTitle()+".pdf").c_str());
    c_lr_low_ref ->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_lr_low_ref ->GetTitle()+".pdf").c_str());
    c_sr_high_ref->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_sr_high_ref->GetTitle()+".pdf").c_str());
    c_lr_high_ref->Print(("Image/PeriSub/JetPeakFitting/"+(std::string)c_lr_high_ref->GetTitle()+".pdf").c_str());
}

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

    TGraphErrors* v2obs_om;
    TGraphErrors* v2obs_sub_om;
    TGraphErrors* v2bkg_om;
    TGraphErrors* v2bkg_sub_om;

    TGraphErrors* v2obs_KET_om;
    TGraphErrors* v2obs_KET_sub_om;
    TGraphErrors* v2bkg_KET_om;
    TGraphErrors* v2bkg_KET_sub_om;

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

    f->GetObject("v2plot_obs_om"   ,v2obs_om);
    f->GetObject("v2subplot_obs_om",v2obs_sub_om);
    f->GetObject("v2plot_bkg_om"   ,v2bkg_om);
    f->GetObject("v2subplot_bkg_om",v2bkg_sub_om);

    f->GetObject("v2plot_KET_obs_om"   ,v2obs_KET_om);
    f->GetObject("v2subplot_KET_obs_om",v2obs_KET_sub_om);
    f->GetObject("v2plot_KET_bkg_om"   ,v2bkg_KET_om);
    f->GetObject("v2subplot_KET_bkg_om",v2bkg_KET_sub_om);

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

    double* v2_obs_om     = v2obs_om->GetY();
    double* v2_obs_err_om = v2obs_om->GetEY();
    double* v2_bkg_om     = v2bkg_om->GetY();
    double* v2_bkg_err_om = v2bkg_om->GetEY();

    double* v2_obs_sub_om     = v2obs_sub_om->GetY();
    double* v2_obs_err_sub_om = v2obs_sub_om->GetEY();
    double* v2_bkg_sub_om     = v2bkg_sub_om->GetY();
    double* v2_bkg_err_sub_om = v2bkg_sub_om->GetEY();

    double* pt_ks = v2obs_ks->GetX();
    double* pt_la = v2obs_la->GetX();
    double* pt_xi = v2obs_xi->GetX();
    double* pt_om = v2obs_om->GetX();

    double* KET_ks = v2obs_KET_ks->GetX();
    double* KET_la = v2obs_KET_la->GetX();
    double* KET_xi = v2obs_KET_xi->GetX();
    double* KET_om = v2obs_KET_om->GetX();

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

    std::vector<double> bkgfrac_om;
    std::vector<double> v2true_om;
    std::vector<double> v2true_err_om;
    std::vector<double> v2true_sub_om;
    std::vector<double> v2true_sub_err_om;

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

    for(unsigned i=0; i<PD.PtBin_om.size(); i++)
    {
        bkgfrac_om.push_back(1 - PD.fsig_om[i]);
        v2true_om.push_back((v2_obs_om[i] - v2_bkg_om[i]*bkgfrac_om[i])/PD.fsig_om[i]);
        v2true_err_om.push_back(sqrt(TMath::Power((v2_bkg_err_om[i]*bkgfrac_om[i]),2) + TMath::Power(v2_obs_err_om[i],2))/PD.fsig_om[i]);
        v2true_sub_om.push_back((v2_obs_sub_om[i] - v2_bkg_sub_om[i]*bkgfrac_om[i])/PD.fsig_om[i]);
        v2true_sub_err_om.push_back(sqrt(TMath::Power(v2_bkg_err_sub_om[i]*bkgfrac_om[i],2) + TMath::Power(v2_obs_err_sub_om[i],2))/PD.fsig_om[i]);
    }

    int numPtBins_ks = PD.PtBin_ks.size()-1;
    int numPtBins_la = PD.PtBin_la.size()-1;
    int numPtBins_xi = PD.PtBin_xi.size()-1;
    int numPtBins_om = PD.PtBin_om.size()-1;

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

    TGraphErrors* omv2true    = new TGraphErrors(numPtBins_om,pt_om,&v2true_om[0]    ,0,&v2true_err_om[0]);
    TGraphErrors* omv2truesub = new TGraphErrors(numPtBins_om,pt_om,&v2true_sub_om[0],0,&v2true_sub_err_om[0]);

    TGraphErrors* omv2true_KET    = new TGraphErrors(numPtBins_om,KET_om,&v2true_om[0]    ,0,&v2true_err_om[0]);
    TGraphErrors* omv2truesub_KET = new TGraphErrors(numPtBins_om,KET_om,&v2true_sub_om[0],0,&v2true_sub_err_om[0]);

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

    omv2true       ->Write("omv2true"       ,TObject::kOverwrite);
    omv2truesub    ->Write("omv2truesub"    ,TObject::kOverwrite);
    omv2true_KET   ->Write("omv2true_KET"   ,TObject::kOverwrite);
    omv2truesub_KET->Write("omv2truesub_KET",TObject::kOverwrite);
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
        DirectSubYield_ks.push_back((YieldPlot_obs_ks[i] - bkgfrac_ks[i]*YieldPlot_bkg_ks[i])/PD.fsig_ks[i]);
        DirectSubYield_err_ks.push_back(sqrt(TMath::Power((YieldPlot_bkg_err_ks[i]*bkgfrac_ks[i]),2) + TMath::Power(YieldPlot_obs_err_ks[i],2))/PD.fsig_ks[i]);
        DirectSubYield_low_ks.push_back((YieldPlot_low_obs_ks[i] - bkgfrac_MB_ks[i]*YieldPlot_low_bkg_ks[i])/PD.fsig_ks_MB[i]);
        DirectSubYield_low_err_ks.push_back(sqrt(TMath::Power((YieldPlot_low_bkg_err_ks[i]*bkgfrac_MB_ks[i]),2) + TMath::Power(YieldPlot_low_obs_err_ks[i],2))/PD.fsig_ks_MB[i]);
        DirectSubNass_ks.push_back((Nassoc_obs_ks[i] - bkgfrac_ks[i]*Nassoc_bkg_ks[i])/PD.fsig_ks[i]);
        DirectSubNass_low_ks.push_back((Nassoc_low_obs_ks[i] - bkgfrac_MB_ks[i]*Nassoc_low_bkg_ks[i])/PD.fsig_ks_MB[i]);

        v2true_DirectSub_ks.push_back(ksv2true_Y[i] - ksv2_MB_Y[i]*DirectSubNass_low_ks[i]/DirectSubNass_ks[i]*DirectSubYield_ks[i]/DirectSubYield_low_ks[i]);
        v2true_DirectSub_err_ks.push_back(sqrt(TMath::Power(ksv2true_err[i],2) + (TMath::Power(ksv2_MB_err[i]/ksv2_MB_Y[i],2) + TMath::Power(DirectSubYield_err_ks[i]/DirectSubYield_ks[i],2) + TMath::Power(DirectSubYield_low_err_ks[i]/DirectSubYield_low_ks[i],2))*TMath::Power(ksv2_MB_Y[i]*DirectSubNass_low_ks[i]/DirectSubNass_ks[i]*DirectSubYield_ks[i]/DirectSubYield_low_ks[i],2)));
    }

    for(unsigned i=0; i<PD.PtBin_la.size(); i++)
    {
        cout<< "Y obs " << i << ": " << YieldPlot_obs_la[i] << endl;
        cout<< " Y bkg " << i << ": " << YieldPlot_bkg_la[i] << endl;
        bkgfrac_la.push_back(1 - PD.fsig_la[i]);
        bkgfrac_MB_la.push_back(1 - PD.fsig_la_MB[i]);
        cout << "bkgfrac " << i << ": " << bkgfrac_la[i] << endl;
        cout << "bkgfrac MB" << i << ": " <<  bkgfrac_MB_la[i] << endl;
        DirectSubYield_la.push_back((YieldPlot_obs_la[i] - bkgfrac_la[i]*YieldPlot_bkg_la[i])/PD.fsig_la[i]);
        DirectSubYield_err_la.push_back(sqrt(TMath::Power((YieldPlot_bkg_err_la[i]*bkgfrac_la[i]),2) + TMath::Power(YieldPlot_obs_err_la[i],2))/PD.fsig_la[i]);
        DirectSubYield_low_la.push_back((YieldPlot_low_obs_la[i] - bkgfrac_MB_la[i]*YieldPlot_low_bkg_la[i])/PD.fsig_la_MB[i]);
        DirectSubYield_low_err_la.push_back(sqrt(TMath::Power((YieldPlot_low_bkg_err_la[i]*bkgfrac_MB_la[i]),2) + TMath::Power(YieldPlot_low_obs_err_la[i],2))/PD.fsig_la_MB[i]);
        DirectSubNass_la.push_back((Nassoc_obs_la[i] - bkgfrac_la[i]*Nassoc_bkg_la[i])/PD.fsig_la[i]);
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
    TGraphErrors* TGDirectSubNass_ks = new TGraphErrors(numPtBins_ks,pt_ks,&DirectSubNass_ks[0],0,0);
    TGraphErrors* TGDirectSubNass_low_ks = new TGraphErrors(numPtBins_ks,pt_ks,&DirectSubNass_low_ks[0],0,0);
    TGraphErrors* TGDirectSubNass_la = new TGraphErrors(numPtBins_la,pt_la,&DirectSubNass_la[0],0,0);
    TGraphErrors* TGDirectSubNass_low_la = new TGraphErrors(numPtBins_la,pt_la,&DirectSubNass_low_la[0],0,0);
    TGraphErrors* TGDirectSubNass_xi = new TGraphErrors(numPtBins_xi,pt_xi,&DirectSubNass_xi[0],0,0);
    TGraphErrors* TGDirectSubNass_low_xi = new TGraphErrors(numPtBins_xi,pt_xi,&DirectSubNass_low_xi[0],0,0);

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
    TGDirectSubYield_xi->Write("DirectSubYield_xi",TObject::kOverwrite);
    TGDirectSubYield_low_xi->Write("DirectSubYield_low_xi",TObject::kOverwrite);
    TGDirectSubNass_ks->Write("DirectSubNass_ks",TObject::kOverwrite);
    TGDirectSubNass_low_ks->Write("DirectSubNass_low_ks",TObject::kOverwrite);
    TGDirectSubNass_la->Write("DirectSubNass_la",TObject::kOverwrite);
    TGDirectSubNass_low_la->Write("DirectSubNass_low_la",TObject::kOverwrite);
    TGDirectSubNass_xi->Write("DirectSubNass_xi",TObject::kOverwrite);
    TGDirectSubNass_low_xi->Write("DirectSubNass_low_xi",TObject::kOverwrite);
    for(int i =0; i<PD.PtBin_la.size()-1; i++)
    {
        cout << "DirectSubYield_la " << i << ": " << DirectSubYield_la[i] << endl;
    }
}

void PrintV2(ParticleData PD)
{
    TFile* f_corr = TFile::Open("FitRootFiles/PeriSub/V0v2perisub_Default_AllStrange_EG1_0_35_CorrectRef_CorrectGap_12_06_17.root");
    TFile* f_incorr = TFile::Open("FitRootFiles/PeriSub/V0v2perisub_Default_EG1_0_35_CorrectRef_InCorrectGap_12_05_17.root");
    TFile* f_corr_20 = TFile::Open("FitRootFiles/PeriSub/V0v2perisub_Default_EG1_0_20_CorrectRef_CorrectGap_12_06_17.root");
    TFile* f_incorr_20 = TFile::Open("FitRootFiles/PeriSub/V0v2perisub_Default_EG1_0_20_CorrectRef_InCorrectGap_12_05_17.root");

    TGraphErrors* Nasslow_obs_ks_corr = (TGraphErrors*)f_corr->Get("Nassoc_low_obs_ks");
    TGraphErrors* Nasslow_bkg_ks_corr = (TGraphErrors*)f_corr->Get("Nassoc_low_bkg_ks");
    TGraphErrors* Nasslow_obs_la_corr = (TGraphErrors*)f_corr->Get("Nassoc_low_obs_la");
    TGraphErrors* Nasslow_bkg_la_corr = (TGraphErrors*)f_corr->Get("Nassoc_low_bkg_la");
    TGraphErrors* Nasslow_obs_xi_corr = (TGraphErrors*)f_corr->Get("Nassoc_low_obs_xi");
    TGraphErrors* Nasslow_bkg_xi_corr = (TGraphErrors*)f_corr->Get("Nassoc_low_bkg_xi");
    TGraphErrors* Nasslow_obs_om_corr = (TGraphErrors*)f_corr->Get("Nassoc_low_obs_om");
    TGraphErrors* Nasslow_bkg_om_corr = (TGraphErrors*)f_corr->Get("Nassoc_low_bkg_om");

    TGraphErrors* V2low_obs_ks_corr = (TGraphErrors*)f_corr->Get("V2plot_obs_low_ks");
    TGraphErrors* V2low_bkg_ks_corr = (TGraphErrors*)f_corr->Get("V2plot_bkg_low_ks");
    TGraphErrors* V2low_obs_la_corr = (TGraphErrors*)f_corr->Get("V2plot_obs_low_la");
    TGraphErrors* V2low_bkg_la_corr = (TGraphErrors*)f_corr->Get("V2plot_bkg_low_la");
    TGraphErrors* V2low_obs_xi_corr = (TGraphErrors*)f_corr->Get("V2plot_obs_low_xi");
    TGraphErrors* V2low_bkg_xi_corr = (TGraphErrors*)f_corr->Get("V2plot_bkg_low_xi");
    TGraphErrors* V2low_obs_om_corr = (TGraphErrors*)f_corr->Get("V2plot_obs_low_om");
    TGraphErrors* V2low_bkg_om_corr = (TGraphErrors*)f_corr->Get("V2plot_bkg_low_om");

    TGraphErrors* Nasslow_obs_ks_incorr = (TGraphErrors*)f_incorr->Get("Nassoc_low_obs_ks");
    TGraphErrors* Nasslow_bkg_ks_incorr = (TGraphErrors*)f_incorr->Get("Nassoc_low_bkg_ks");
    TGraphErrors* Nasslow_obs_la_incorr = (TGraphErrors*)f_incorr->Get("Nassoc_low_obs_la");
    TGraphErrors* Nasslow_bkg_la_incorr = (TGraphErrors*)f_incorr->Get("Nassoc_low_bkg_la");
    TGraphErrors* Nasslow_obs_om_incorr = (TGraphErrors*)f_incorr->Get("Nassoc_low_obs_om");
    TGraphErrors* Nasslow_bkg_om_incorr = (TGraphErrors*)f_incorr->Get("Nassoc_low_bkg_om");

    TGraphErrors* V2low_obs_ks_incorr = (TGraphErrors*)f_incorr->Get("V2plot_obs_low_ks");
    TGraphErrors* V2low_bkg_ks_incorr = (TGraphErrors*)f_incorr->Get("V2plot_bkg_low_ks");
    TGraphErrors* V2low_obs_la_incorr = (TGraphErrors*)f_incorr->Get("V2plot_obs_low_la");
    TGraphErrors* V2low_bkg_la_incorr = (TGraphErrors*)f_incorr->Get("V2plot_bkg_low_la");
    TGraphErrors* V2low_obs_om_incorr = (TGraphErrors*)f_incorr->Get("V2plot_obs_low_om");
    TGraphErrors* V2low_bkg_om_incorr = (TGraphErrors*)f_incorr->Get("V2plot_bkg_low_om");

    TGraphErrors* Nasslow_obs_ks_corr_20 = (TGraphErrors*)f_corr_20->Get("Nassoc_low_obs_ks");
    TGraphErrors* Nasslow_bkg_ks_corr_20 = (TGraphErrors*)f_corr_20->Get("Nassoc_low_bkg_ks");
    TGraphErrors* Nasslow_obs_la_corr_20 = (TGraphErrors*)f_corr_20->Get("Nassoc_low_obs_la");
    TGraphErrors* Nasslow_bkg_la_corr_20 = (TGraphErrors*)f_corr_20->Get("Nassoc_low_bkg_la");
    TGraphErrors* Nasslow_obs_xi_corr_20 = (TGraphErrors*)f_corr_20->Get("Nassoc_low_obs_xi");
    TGraphErrors* Nasslow_bkg_xi_corr_20 = (TGraphErrors*)f_corr_20->Get("Nassoc_low_bkg_xi");

    TGraphErrors* V2low_obs_ks_corr_20 = (TGraphErrors*)f_corr_20->Get("V2plot_obs_low_ks");
    TGraphErrors* V2low_bkg_ks_corr_20 = (TGraphErrors*)f_corr_20->Get("V2plot_bkg_low_ks");
    TGraphErrors* V2low_obs_la_corr_20 = (TGraphErrors*)f_corr_20->Get("V2plot_obs_low_la");
    TGraphErrors* V2low_bkg_la_corr_20 = (TGraphErrors*)f_corr_20->Get("V2plot_bkg_low_la");
    TGraphErrors* V2low_obs_xi_corr_20 = (TGraphErrors*)f_corr_20->Get("V2plot_obs_low_xi");
    TGraphErrors* V2low_bkg_xi_corr_20 = (TGraphErrors*)f_corr_20->Get("V2plot_bkg_low_xi");

    TGraphErrors* Nasslow_obs_ks_incorr_20 = (TGraphErrors*)f_incorr_20->Get("Nassoc_low_obs_ks");
    TGraphErrors* Nasslow_bkg_ks_incorr_20 = (TGraphErrors*)f_incorr_20->Get("Nassoc_low_bkg_ks");
    TGraphErrors* Nasslow_obs_la_incorr_20 = (TGraphErrors*)f_incorr_20->Get("Nassoc_low_obs_la");
    TGraphErrors* Nasslow_bkg_la_incorr_20 = (TGraphErrors*)f_incorr_20->Get("Nassoc_low_bkg_la");

    TGraphErrors* V2low_obs_ks_incorr_20 = (TGraphErrors*)f_incorr_20->Get("V2plot_obs_low_ks");
    TGraphErrors* V2low_bkg_ks_incorr_20 = (TGraphErrors*)f_incorr_20->Get("V2plot_bkg_low_ks");
    TGraphErrors* V2low_obs_la_incorr_20 = (TGraphErrors*)f_incorr_20->Get("V2plot_obs_low_la");
    TGraphErrors* V2low_bkg_la_incorr_20 = (TGraphErrors*)f_incorr_20->Get("V2plot_bkg_low_la");

    Nasslow_obs_ks_corr->SetMarkerColor(kRed);
    Nasslow_obs_ks_corr->SetMarkerStyle(20);
    Nasslow_obs_la_corr->SetMarkerColor(kBlue);
    Nasslow_obs_la_corr->SetMarkerStyle(22);
    Nasslow_obs_xi_corr->SetMarkerColor(kGreen+2);
    Nasslow_obs_xi_corr->SetMarkerStyle(21);
    Nasslow_obs_om_corr->SetMarkerColor(kMagenta);
    Nasslow_obs_om_corr->SetMarkerStyle(29);
    Nasslow_obs_ks_corr_20->SetMarkerColor(kRed);
    Nasslow_obs_ks_corr_20->SetMarkerStyle(24);
    Nasslow_obs_la_corr_20->SetMarkerColor(kBlue);
    Nasslow_obs_la_corr_20->SetMarkerStyle(26);
    Nasslow_obs_xi_corr_20->SetMarkerColor(kGreen+2);
    Nasslow_obs_xi_corr_20->SetMarkerStyle(25);

    Nasslow_bkg_ks_corr->SetMarkerColor(kRed);
    Nasslow_bkg_ks_corr->SetMarkerStyle(20);
    Nasslow_bkg_la_corr->SetMarkerColor(kBlue);
    Nasslow_bkg_la_corr->SetMarkerStyle(22);
    Nasslow_bkg_xi_corr->SetMarkerColor(kGreen+2);
    Nasslow_bkg_xi_corr->SetMarkerStyle(21);
    Nasslow_bkg_om_corr->SetMarkerColor(kMagenta);
    Nasslow_bkg_om_corr->SetMarkerStyle(29);
    Nasslow_bkg_ks_corr_20->SetMarkerColor(kRed);
    Nasslow_bkg_ks_corr_20->SetMarkerStyle(24);
    Nasslow_bkg_la_corr_20->SetMarkerColor(kBlue);
    Nasslow_bkg_la_corr_20->SetMarkerStyle(26);
    Nasslow_bkg_xi_corr_20->SetMarkerColor(kGreen+2);
    Nasslow_bkg_xi_corr_20->SetMarkerStyle(25);


    TCanvas* CompNass = new TCanvas("CompNass","CompNass",1200,600);
    CompNass->Divide(1,2);
    TH1F* Frame1 = CompNass->cd(1)->DrawFrame(0,0,9,2);
    Frame1->GetYaxis()->SetTitle("Nassoc Obs");
    TLegend* leg1 = new TLegend(0.20,0.60,0.27,0.8);
    leg1->AddEntry(Nasslow_obs_ks_corr,"Ks 0-35", "P");
    leg1->AddEntry(Nasslow_obs_la_corr,"La 0-35", "P");
    leg1->AddEntry(Nasslow_obs_ks_corr_20,"Ks 0-20", "P");
    leg1->AddEntry(Nasslow_obs_la_corr_20,"La 0-20", "P");
    leg1->Draw();
    Nasslow_obs_ks_corr->Draw("P");
    Nasslow_obs_la_corr->Draw("PSAME");
    Nasslow_obs_ks_corr_20->Draw("PSAME");
    Nasslow_obs_la_corr_20->Draw("PSAME");

    TH1F* Frame2 = CompNass->cd(2)->DrawFrame(0,0,9,2);
    Frame2->GetYaxis()->SetTitle("Nassoc Bkg");
    TLegend* leg2 = new TLegend(0.20,0.60,0.27,0.8);
    leg2->AddEntry(Nasslow_bkg_ks_corr,"Ks 0-35", "P");
    leg2->AddEntry(Nasslow_bkg_la_corr,"La 0-35", "P");
    leg2->AddEntry(Nasslow_bkg_ks_corr_20,"Ks 0-20", "P");
    leg2->AddEntry(Nasslow_bkg_la_corr_20,"La 0-20", "P");
    leg2->Draw();
    Nasslow_bkg_ks_corr->Draw("P");
    Nasslow_bkg_la_corr->Draw("PSAME");
    Nasslow_bkg_ks_corr_20->Draw("PSAME");
    Nasslow_bkg_la_corr_20->Draw("PSAME");

    TFile* f_print = new TFile("Nassoc_V0CasOmComparison_20_v_35.root","RECREATE");
    Nasslow_obs_ks_corr->Write("Kshort_obs_35");
    Nasslow_obs_ks_corr_20->Write("Kshort_obs_20");
    Nasslow_obs_la_corr->Write("Lambda_obs_35");
    Nasslow_obs_la_corr_20->Write("Lambda_obs_20");
    Nasslow_obs_xi_corr->Write("Xi_obs_35");
    Nasslow_obs_xi_corr_20->Write("Xi_obs_20");
    Nasslow_obs_om_corr->Write("Omega_obs_35");
    Nasslow_bkg_ks_corr->Write("Kshort_bkg_35");
    Nasslow_bkg_ks_corr_20->Write("Kshort_bkg_20");
    Nasslow_bkg_la_corr->Write("Lambda_bkg_35");
    Nasslow_bkg_la_corr_20->Write("Lambda_bkg_20");
    Nasslow_bkg_xi_corr->Write("Xi_bkg_35");
    Nasslow_bkg_xi_corr_20->Write("Xi_bkg_20");
    Nasslow_bkg_om_corr->Write("Omega_bkg_35");


    double* Nasslow_obs_ks_Y_corr = Nasslow_obs_ks_corr->GetY();
    double* Nasslow_bkg_ks_Y_corr = Nasslow_bkg_ks_corr->GetY();
    double* Nasslow_obs_la_Y_corr = Nasslow_obs_la_corr->GetY();
    double* Nasslow_bkg_la_Y_corr = Nasslow_bkg_la_corr->GetY();
    double* Nasslow_obs_om_Y_corr = Nasslow_obs_om_corr->GetY();
    double* Nasslow_bkg_om_Y_corr = Nasslow_bkg_om_corr->GetY();

    double* V2low_obs_ks_Y_corr = V2low_obs_ks_corr->GetY();
    double* V2low_bkg_ks_Y_corr = V2low_bkg_ks_corr->GetY();
    double* V2low_obs_la_Y_corr = V2low_obs_la_corr->GetY();
    double* V2low_bkg_la_Y_corr = V2low_bkg_la_corr->GetY();
    double* V2low_obs_om_Y_corr = V2low_obs_om_corr->GetY();
    double* V2low_bkg_om_Y_corr = V2low_bkg_om_corr->GetY();

    double* Nasslow_obs_ks_Y_incorr = Nasslow_obs_ks_incorr->GetY();
    double* Nasslow_bkg_ks_Y_incorr = Nasslow_bkg_ks_incorr->GetY();
    double* Nasslow_obs_la_Y_incorr = Nasslow_obs_la_incorr->GetY();
    double* Nasslow_bkg_la_Y_incorr = Nasslow_bkg_la_incorr->GetY();
    double* Nasslow_obs_om_Y_incorr = Nasslow_obs_om_incorr->GetY();
    double* Nasslow_bkg_om_Y_incorr = Nasslow_bkg_om_incorr->GetY();

    double* V2low_obs_ks_Y_incorr = V2low_obs_ks_incorr->GetY();
    double* V2low_bkg_ks_Y_incorr = V2low_bkg_ks_incorr->GetY();
    double* V2low_obs_la_Y_incorr = V2low_obs_la_incorr->GetY();
    double* V2low_bkg_la_Y_incorr = V2low_bkg_la_incorr->GetY();
    double* V2low_obs_om_Y_incorr = V2low_obs_om_incorr->GetY();
    double* V2low_bkg_om_Y_incorr = V2low_bkg_om_incorr->GetY();

    double* Nasslow_obs_ks_Y_corr_20 = Nasslow_obs_ks_corr_20->GetY();
    double* Nasslow_bkg_ks_Y_corr_20 = Nasslow_bkg_ks_corr_20->GetY();
    double* Nasslow_obs_la_Y_corr_20 = Nasslow_obs_la_corr_20->GetY();
    double* Nasslow_bkg_la_Y_corr_20 = Nasslow_bkg_la_corr_20->GetY();

    double* V2low_obs_ks_Y_corr_20 = V2low_obs_ks_corr_20->GetY();
    double* V2low_bkg_ks_Y_corr_20 = V2low_bkg_ks_corr_20->GetY();
    double* V2low_obs_la_Y_corr_20 = V2low_obs_la_corr_20->GetY();
    double* V2low_bkg_la_Y_corr_20 = V2low_bkg_la_corr_20->GetY();

    double* Nasslow_obs_ks_Y_incorr_20 = Nasslow_obs_ks_incorr_20->GetY();
    double* Nasslow_bkg_ks_Y_incorr_20 = Nasslow_bkg_ks_incorr_20->GetY();
    double* Nasslow_obs_la_Y_incorr_20 = Nasslow_obs_la_incorr_20->GetY();
    double* Nasslow_bkg_la_Y_incorr_20 = Nasslow_bkg_la_incorr_20->GetY();

    double* V2low_obs_ks_Y_incorr_20 = V2low_obs_ks_incorr_20->GetY();
    double* V2low_bkg_ks_Y_incorr_20 = V2low_bkg_ks_incorr_20->GetY();
    double* V2low_obs_la_Y_incorr_20 = V2low_obs_la_incorr_20->GetY();
    double* V2low_bkg_la_Y_incorr_20 = V2low_bkg_la_incorr_20->GetY();

    cout << "Low N Kshort obs Nassoc (0-35 Correct gap | Incorrect gap | difference || 0-20 ...)" << endl;
    for(int i=0; i<Nasslow_obs_ks_corr->GetN(); i++)
    {
        cout << "Pt bin " << i << ": " << std::left << std::setw(10) << Nasslow_obs_ks_Y_corr[i] << std::left << std::setw(10) << Nasslow_obs_ks_Y_incorr[i] << std::left << std::setw(10) << Nasslow_obs_ks_Y_corr[i] - Nasslow_obs_ks_Y_incorr[i] <<  "\t || \t" << std::left << std::setw(10) << Nasslow_obs_ks_Y_corr_20[i] << std::left << std::setw(10) << Nasslow_obs_ks_Y_incorr_20[i] << std::left << std::setw(10) << Nasslow_obs_ks_Y_corr_20[i] - Nasslow_obs_ks_Y_incorr_20[i] << endl;
    }

    cout << "Low N Kshort bkg Nassoc" << endl;
    for(int i=0; i<Nasslow_bkg_ks_corr->GetN(); i++)
    {
        cout << "Pt bin " << i << ": " << std::left << std::setw(10) << Nasslow_bkg_ks_Y_corr[i] << std::left << std::setw(10) << Nasslow_bkg_ks_Y_incorr[i] << std::left << std::setw(10) << Nasslow_bkg_ks_Y_corr[i] - Nasslow_bkg_ks_Y_incorr[i] << "\t || \t" << std::left << std::setw(10) << Nasslow_bkg_ks_Y_corr_20[i] << std::left << std::setw(10) << Nasslow_bkg_ks_Y_incorr_20[i] << std::left << std::setw(10) << Nasslow_bkg_ks_Y_corr_20[i] - Nasslow_bkg_ks_Y_incorr_20[i] << endl;
    }

    cout << "Low N Lambda obs Nassoc" << endl;
    for(int i=0; i<Nasslow_obs_la_corr->GetN(); i++)
    {
        cout << "Pt bin " << i << ": " << std::left << std::setw(10) << Nasslow_obs_la_Y_corr[i] << std::left << std::setw(10) << Nasslow_obs_la_Y_incorr[i] << std::left << std::setw(10) << Nasslow_obs_la_Y_corr[i] - Nasslow_obs_la_Y_incorr[i] <<  "\t || \t" << std::left << std::setw(10) << Nasslow_obs_la_Y_corr_20[i] << std::left << std::setw(10) << Nasslow_obs_la_Y_incorr_20[i] << std::left << std::setw(10) << Nasslow_obs_la_Y_corr_20[i] - Nasslow_obs_la_Y_incorr_20[i] << endl;
    }

    cout << "Low N Lambda bkg Nassoc" << endl;
    for(int i=0; i<Nasslow_bkg_la_corr->GetN(); i++)
    {
        cout << "Pt bin " << i << ": " << std::left << std::setw(10) << Nasslow_bkg_la_Y_corr[i] << std::left << std::setw(10) << Nasslow_bkg_la_Y_incorr[i] << std::left << std::setw(10) << Nasslow_bkg_la_Y_corr[i] - Nasslow_bkg_la_Y_incorr[i] << "\t || \t" << std::left << std::setw(10) << Nasslow_bkg_la_Y_corr_20[i] << std::left << std::setw(10) << Nasslow_bkg_la_Y_incorr_20[i] << std::left << std::setw(10) << Nasslow_bkg_la_Y_corr_20[i] - Nasslow_bkg_la_Y_incorr_20[i] << endl;
    }

    cout << "Low N Omega obs Nassoc" << endl;
    for(int i=0; i<Nasslow_obs_om_corr->GetN(); i++)
    {
        cout << "Pt bin " << i << ": " << std::left << std::setw(15) << Nasslow_obs_om_Y_corr[i] << std::left << std::setw(15) << Nasslow_bkg_om_Y_incorr[i] << std::left << std::setw(15) << Nasslow_obs_om_Y_corr[i] - Nasslow_bkg_om_Y_incorr[i] <<  endl;
    }

    cout << "Low N Omega bkg Nassoc" << endl;
    for(int i=0; i<Nasslow_bkg_om_corr->GetN(); i++)
    {
        cout << "Pt bin " << i << ": " << std::left << std::setw(15) << Nasslow_bkg_om_Y_corr[i] << std::left << std::setw(15) << Nasslow_bkg_om_Y_incorr[i] << std::left << std::setw(15) << Nasslow_bkg_om_Y_corr[i] - Nasslow_bkg_om_Y_incorr[i] << endl;
    }

    //V2
    cout << "Low N Kshort obs V2" << endl;
    for(int i=0; i<V2low_obs_ks_corr->GetN(); i++)
    {
        cout << "Pt bin " << i << ": " << std::left << std::setw(15) << V2low_obs_ks_Y_corr[i] << std::left << std::setw(15) << V2low_obs_ks_Y_incorr[i] << std::left << std::setw(15) << V2low_obs_ks_Y_corr[i] - V2low_obs_ks_Y_incorr[i] << "\t || \t" << std::left << std::setw(15) << V2low_obs_ks_Y_corr_20[i] << std::left << std::setw(15) << V2low_obs_ks_Y_incorr_20[i] << std::left << std::setw(15) << V2low_obs_ks_Y_corr_20[i] - V2low_obs_ks_Y_incorr_20[i] << endl;
    }

    cout << "Low N Kshort bkg V2" << endl;
    for(int i=0; i<V2low_bkg_ks_corr->GetN(); i++)
    {
        cout << "Pt bin " << i << ": " << std::left << std::setw(15) << V2low_bkg_ks_Y_corr[i] << std::left << std::setw(15) << V2low_bkg_ks_Y_incorr[i] << std::left << std::setw(15) << V2low_bkg_ks_Y_corr[i] - V2low_bkg_ks_Y_incorr[i] << "\t || \t" << std::left << std::setw(15) << V2low_bkg_ks_Y_corr_20[i] << std::left << std::setw(15) << V2low_bkg_ks_Y_incorr_20[i] << std::left << std::setw(15) << V2low_bkg_ks_Y_corr_20[i] - V2low_bkg_ks_Y_incorr_20[i] << endl;
    }

    cout << "Low N Lambda obs V2" << endl;
    for(int i=0; i<V2low_obs_la_corr->GetN(); i++)
    {
        cout << "Pt bin " << i << ": " << std::left << std::setw(15) << V2low_obs_la_Y_corr[i] << std::left << std::setw(15) << V2low_obs_la_Y_incorr[i] << std::left << std::setw(15) << V2low_obs_la_Y_corr[i] - V2low_obs_la_Y_incorr[i] << "\t || \t" << std::left << std::setw(15) << V2low_obs_la_Y_corr_20[i] << std::left << std::setw(15) << V2low_obs_la_Y_incorr_20[i] << std::left << std::setw(15) << V2low_obs_la_Y_corr_20[i] - V2low_obs_la_Y_incorr_20[i] << endl;
    }

    cout << "Low N Lambda bkg V2" << endl;
    for(int i=0; i<V2low_bkg_la_corr->GetN(); i++)
    {
        cout << "Pt bin "<< i << ": " << std::left << std::setw(15) << V2low_bkg_la_Y_corr[i] << std::left << std::setw(15) << V2low_bkg_la_Y_incorr[i] << std::left << std::setw(15) << V2low_bkg_la_Y_corr[i] - V2low_bkg_la_Y_incorr[i] << "\t || \t" << std::left << std::setw(15) << V2low_bkg_la_Y_corr_20[i] << std::left << std::setw(15) << V2low_bkg_la_Y_incorr_20[i] << std::left << std::setw(15) << V2low_bkg_la_Y_corr_20[i] - V2low_bkg_la_Y_incorr_20[i] << endl;
    }

    cout << "Low N Omega obs V2" << endl;
    for(int i=0; i<V2low_obs_om_corr->GetN(); i++)
    {
        cout << "Pt bin " << i << ": " << V2low_obs_om_Y_corr[i] << "\t" << V2low_obs_om_Y_incorr[i] << "\t" << V2low_obs_om_Y_corr[i] - V2low_obs_om_Y_incorr[i] << endl;
    }

    cout << "Low N Omega bkg V2" << endl;
    for(int i=0; i<V2low_bkg_om_corr->GetN(); i++)
    {
        cout << "Pt bin " << i << ": " << V2low_bkg_om_Y_corr[i] << "\t" << V2low_bkg_om_Y_incorr[i] << "\t" << V2low_bkg_om_Y_corr[i] - V2low_bkg_om_Y_incorr[i] << endl;
    }

}

int main()
{
    //PeriSubFixedWindow(V0);
    PeriSubObs(V0);
    PeriSubSigIndirect(V0);
    PeriSubSigDirect(V0);
    PrintV2(V0);
    return 1;
}
