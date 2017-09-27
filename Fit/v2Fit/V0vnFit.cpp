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

#define PI 3.1416

Double_t FourierHad(Double_t *x, Double_t *par)
{
    //Double_t xx1 = par[0]/(2*PI);
    Double_t xx1 = par[0];
    Double_t xx2 = 1 + 2*(par[1]*TMath::Cos(x[0]) + par[2]*TMath::Cos(2*x[0])  + par[3]*TMath::Cos(3*x[0])) + par[4]*TMath::Cos(4*x[0]);// + par[5]*TMath::Cos(5*x[0]) + par[6]*TMath::Cos(6*x[0]) + par[7]*TMath::Cos(7*x[0]));
    return xx1*xx2;
}

void OutputVnValues(int degree, std::string region, std::string V0IDname, std::vector<double> vnvalues, std::vector<double> ptbin, std::string file)
{
    std::ostringstream output;
    std::ofstream myfile;
    myfile.open(file.c_str(),std::ios_base::app);

    output << V0IDname << " V" << degree << " values " << region << "\n";
    cout << "==========================================================" << endl;
    cout << output.str();
    cout << "==========================================================" << endl;

    myfile << output.str();

    int PtBinCounter=0;
    for(std::vector<double>::iterator it = vnvalues.begin(); it != vnvalues.end(); ++it)
    {
        cout << ptbin[PtBinCounter] << " < Pt =< " << ptbin[PtBinCounter + 1] << ": " << *it << endl;
        myfile << *it << "\n";
        PtBinCounter++;
    }
    //myfile.close();
}

void OutputVnErrors(int degree, std::string region, std::string V0IDname, std::vector<double> vnerrors, std::vector<double> ptbin, std::string file)
{
    std::ostringstream output;
    std::ofstream myfile;
    myfile.open(file.c_str(),std::ios_base::app);

    output << V0IDname << " V" << degree << " errors " << region << "\n";
    cout << "==========================================================" << endl;
    cout << output.str();
    cout << "==========================================================" << endl;

    myfile << output.str();

    int PtBinCounter=0;
    for(std::vector<double>::iterator it = vnerrors.begin(); it != vnerrors.end(); ++it)
    {
        cout << ptbin[PtBinCounter] << " < Pt =< " << ptbin[PtBinCounter + 1] << ": " << *it << endl;
        myfile << *it << "\n";
        PtBinCounter++;
    }
}

void OutputVnH(int degree, std::vector<double> vnvalues, std::vector<double> vnerrors, std::string file)
{
    std::ostringstream output;
    std::ofstream myfile;
    myfile.open(file.c_str(),std::ios_base::app);

    output << "Hadron V" << degree << " value THEN error\n";

    cout << output.str();
    myfile << output.str();

    cout << vnvalues[degree] << endl;
    cout << vnerrors[degree] << endl;
    myfile << vnvalues[degree] << "\n";
    myfile << vnerrors[degree] << "\n";
}

void vnCalculate(int degree, std::string V0IDname, std::vector<double> vnvalues_peak, std::vector<double> vnerrors_peak, std::vector<double> vnvalues_side, std::vector<double> vnerrors_side, std::vector<double> vnvalues_h, std::vector<double> vnerrors_h, std::vector<double> fsig, std::string file)
{
    std::ostringstream output;
    std::ofstream myfile;
    myfile.open(file.c_str(),std::ios_base::app);

    std::vector<double> sigvalues;
    std::vector<double> sigerrors;
    std::vector<double> vnObs;
    std::vector<double> vnBkg;
    std::vector<double> vnRefError;

    output << V0IDname << " signal v" << degree << " values\n";
    cout << output.str();

    for(unsigned i=0; i<fsig.size(); i++)
    {
        cout << "Pt Bin " << i+1 << endl;
        double vnObs_ = vnvalues_peak[i]/TMath::Sqrt(vnvalues_h[degree]);
        double vnBkg_ = vnvalues_side[i]/TMath::Sqrt(vnvalues_h[degree]);
        vnObs.push_back(vnObs_);
        vnBkg.push_back(vnBkg_);
        double sig = (vnObs_ - (1 - fsig[i])*vnBkg_)/fsig[i];

        double vnObsError = vnObs_*TMath::Sqrt(TMath::Power(vnerrors_peak[i]/vnvalues_peak[i],2) + TMath::Power(0.5*vnerrors_h[degree]/vnvalues_h[degree],2));
        double vnBkgError = vnBkg_*TMath::Sqrt(TMath::Power(vnerrors_side[i]/vnvalues_side[i],2) + TMath::Power(0.5*vnerrors_h[degree]/vnvalues_h[degree],2));
        double sigError = TMath::Sqrt(vnObsError*vnObsError + TMath::Power(vnBkgError*(1-fsig[i]),2))/fsig[i];

        cout << "Sig: " << sig << endl;
        cout << "Error:" << sigError << endl;

        vnRefError.push_back(0.5*vnerrors_h[degree]/TMath::Sqrt(vnvalues_h[degree]));

        sigvalues.push_back(sig);
        sigerrors.push_back(sigError);
    }

    int divFactor = 1;
    if(V0IDname == "Kshort") divFactor = 2;
    if(V0IDname == "Lambda") divFactor = 3;

    myfile << output.str();
    output.str(std::string());
    for(unsigned i=0; i<sigvalues.size(); i++) myfile << sigvalues[i] << "\n";

    myfile << " signal v2/nq\n";
    for(unsigned i=0; i<sigvalues.size(); i++) myfile << sigvalues[i]/divFactor << "\n";

    output << V0IDname << " signal v" << degree << " errors\n";
    myfile << output.str();
    output.str(std::string());
    for(unsigned i=0; i<sigerrors.size(); i++) myfile << sigerrors[i] << "\n";

    output << V0IDname << " signal v" << degree << "/nq errors\n";
    myfile << output.str();
    output.str(std::string());
    for(unsigned i=0; i<sigerrors.size(); i++) myfile << sigerrors[i]/divFactor << "\n";

    output << V0IDname << " observed v" << degree << " values\n";
    myfile << output.str();
    output.str(std::string());
    for(unsigned i=0; i<vnvalues_peak.size(); i++) myfile << vnvalues_peak[i]/TMath::Sqrt(vnvalues_h[degree]) << "\n";

    output << V0IDname << " observed v" << degree << " errors\n";
    myfile << output.str();
    output.str(std::string());
    for(unsigned i=0; i<vnvalues_peak.size(); i++) myfile << vnObs[i]*TMath::Sqrt(TMath::Power(vnerrors_peak[i]/vnvalues_peak[i],2) + TMath::Power(vnRefError[i]/vnvalues_h[degree],2)) << "\n";

    output << V0IDname << " background v" << degree << " values\n";
    myfile << output.str();
    output.str(std::string());
    for(unsigned i=0; i<vnvalues_side.size(); i++) myfile << vnvalues_side[i]/TMath::Sqrt(vnvalues_h[degree]) << "\n";

    output << V0IDname << " background v" << degree << " errors\n";
    myfile << output.str();
    output.str(std::string());
    for(unsigned i=0; i<vnvalues_side.size(); i++) myfile << vnObs[i]*TMath::Sqrt(TMath::Power(vnerrors_side[i]/vnvalues_side[i],2) + TMath::Power(vnRefError[i]/vnvalues_h[degree],2)) << "\n";
}

std::map<std::string, std::vector<double> > vnCalculateMap(int degree, std::string V0IDname, std::vector<double> vnvalues_peak, std::vector<double> vnerrors_peak, std::vector<double> vnvalues_side, std::vector<double> vnerrors_side, std::vector<double> vnvalues_h, std::vector<double> vnerrors_h, std::vector<double> fsig)
{
    std::map<std::string, std::vector<double> > returnContainer;
    std::vector<double> sigvalues;
    std::vector<double> sigvalues_nq;
    std::vector<double> sigerrors;
    std::vector<double> sigerrors_nq;

    std::vector<double> vnObs;
    std::vector<double> vnObs_errors;
    std::vector<double> vnObs_nq;
    std::vector<double> vnObs_errors_nq;

    std::vector<double> vnBkg;
    std::vector<double> vnBkg_errors;
    std::vector<double> vnBkg_nq;
    std::vector<double> vnBkg_errors_nq;
    std::vector<double> vnRefError;

    for(unsigned i=0; i<fsig.size(); i++)
    {
        double vnObs_ = vnvalues_peak[i]/TMath::Sqrt(vnvalues_h[degree]);
        double vnBkg_ = vnvalues_side[i]/TMath::Sqrt(vnvalues_h[degree]);
        vnObs.push_back(vnObs_);
        vnBkg.push_back(vnBkg_);
        double sig = (vnObs_ - (1 - fsig[i])*vnBkg_)/fsig[i];

        double vnObsError = vnObs_*TMath::Sqrt(TMath::Power(vnerrors_peak[i]/vnvalues_peak[i],2) + TMath::Power(0.5*vnerrors_h[degree]/vnvalues_h[degree],2));
        double vnBkgError = vnBkg_*TMath::Sqrt(TMath::Power(vnerrors_side[i]/vnvalues_side[i],2) + TMath::Power(0.5*vnerrors_h[degree]/vnvalues_h[degree],2));
        double sigError = TMath::Sqrt(vnObsError*vnObsError + TMath::Power(vnBkgError*(1-fsig[i]),2))/fsig[i];

        vnRefError.push_back(0.5*vnerrors_h[degree]/TMath::Sqrt(vnvalues_h[degree]));

        sigvalues.push_back(sig);
        sigerrors.push_back(sigError);
        vnObs_errors.push_back(vnObsError);
        vnBkg_errors.push_back(vnBkgError);
    }

    int divFactor = 1;
    if(V0IDname == "Kshort") divFactor = 2;
    if(V0IDname == "Lambda") divFactor = 3;

    for(unsigned i=0; i<sigvalues.size(); i++   ) sigvalues_nq.push_back(sigvalues[i]/divFactor);
    for(unsigned i=0; i<sigerrors.size(); i++   ) sigerrors_nq.push_back(sigerrors[i]/divFactor);
    for(unsigned i=0; i<vnObs.size(); i++       ) vnObs_nq.push_back(vnObs[i]/divFactor);
    for(unsigned i=0; i<vnBkg.size(); i++       ) vnBkg_nq.push_back(vnBkg[i]/divFactor);
    for(unsigned i=0; i<vnObs_errors.size(); i++) vnObs_errors_nq.push_back(vnObs_errors[i]/divFactor);
    for(unsigned i=0; i<vnBkg_errors.size(); i++) vnBkg_errors_nq.push_back(vnBkg_errors[i]/divFactor);

    returnContainer["sig"]           = sigvalues;
    returnContainer["sig_errors"]    = sigerrors;
    returnContainer["sig_nq"]        = sigvalues_nq;
    returnContainer["sig_errors_nq"] = sigerrors_nq;
    returnContainer["obs"]           = vnObs;
    returnContainer["obs_errors"]    = vnObs_errors;
    returnContainer["obs_nq"]        = vnObs_nq;
    returnContainer["obs_errors_nq"] = vnObs_errors_nq;
    returnContainer["bkg"]           = vnBkg;
    returnContainer["bkg_errors"]    = vnBkg_errors;
    returnContainer["bkg_nq"]        = vnBkg_nq;
    returnContainer["bkg_errors_nq"] = vnBkg_errors_nq;

    return returnContainer;
}

std::vector<double> AvgX(TFile* file, std::string branch, std::string branch_bkg, int numPtBins)
{
    std::vector<double> AvgXcoor;
    branch +="%d";
    branch_bkg += "%d";
    for(int i=0; i<numPtBins; i++)
    {
        TH1D* hX = (TH1D*)file->Get(Form(branch.c_str(),i));
        TH1D* hX_bkg = (TH1D*)file->Get(Form(branch_bkg.c_str(),i));

        int nEntries = 0;
        double XTotal = 0;
        for(int j=hX->FindFirstBinAbove(0,1); j<=hX->FindLastBinAbove(0,1); j++)
        {
            double nX = hX->GetBinContent(j);
            double X = nX*(hX->GetBinCenter(j));
            nEntries+=nX;
            XTotal += X;
        }
        for(int j=hX_bkg->FindFirstBinAbove(0,1); j<=hX_bkg->FindLastBinAbove(0,1); j++)
        {
            double nX_bkg = hX_bkg->GetBinContent(j);
            double X_bkg = nX_bkg*(hX_bkg->GetBinCenter(j));
            nEntries += nX_bkg;
            XTotal += X_bkg;
        }
        AvgXcoor.push_back(XTotal/nEntries);
    }

    return AvgXcoor;
}

void vnGraph(std::map<std::string,std::vector<double> > returnContainer, std::vector<double> AvgX_pt, std::vector<double> AvgX_ket, std::string V0ID, std::string fn)
{
    TFile* out = new TFile(fn.c_str(),"UPDATE");
    int divFactor = 1;
    std::vector<double> AvgX_ket_nq = AvgX_ket;
    if(V0ID == "Kshort")
    {
        divFactor = 2;
        for(unsigned i=0; i<AvgX_ket.size(); i++)
        {
            AvgX_ket_nq[i] = AvgX_ket[i]/divFactor;
        }
    }

    if(V0ID == "Lambda")
    {
        divFactor = 3;
        for(unsigned i=0; i<AvgX_ket.size(); i++)
        {
            AvgX_ket_nq[i] = AvgX_ket[i]/divFactor;
        }
    }

    TGraphErrors* v2           = new TGraphErrors(returnContainer["sig"].size(),&AvgX_pt[0],&(returnContainer["sig"])[0],0,&(returnContainer["sig_errors"])[0]);
    TGraphErrors* v2_nq        = new TGraphErrors(returnContainer["sig_nq"].size(),&AvgX_pt[0],&(returnContainer["sig_nq"])[0],0,&(returnContainer["sig_errors_nq"])[0]);
    TGraphErrors* v2_ket       = new TGraphErrors(returnContainer["sig"].size(),&AvgX_ket[0],&(returnContainer["sig"])[0],0,&(returnContainer["sig_errors"])[0]);
    TGraphErrors* v2_ket_nq    = new TGraphErrors(returnContainer["sig_nq"].size(),&AvgX_ket_nq[0],&(returnContainer["sig_nq"])[0],0,&(returnContainer["sig_errors_nq"])[0]);
    TGraphErrors* v2obs        = new TGraphErrors(returnContainer["obs"].size(),&AvgX_pt[0],&(returnContainer["obs"])[0],0,&(returnContainer["obs_errors"])[0]);
    TGraphErrors* v2obs_nq     = new TGraphErrors(returnContainer["obs_nq"].size(),&AvgX_pt[0],&(returnContainer["obs_nq"])[0],0,&(returnContainer["obs_errors_nq"])[0]);
    TGraphErrors* v2obs_ket    = new TGraphErrors(returnContainer["obs"].size(),&AvgX_ket[0],&(returnContainer["obs"])[0],0,&(returnContainer["obs_errors"])[0]);
    TGraphErrors* v2obs_ket_nq = new TGraphErrors(returnContainer["obs_nq"].size(),&AvgX_ket_nq[0],&(returnContainer["obs_nq"])[0],0,&(returnContainer["obs_errors_nq"])[0]);
    TGraphErrors* v2bkg        = new TGraphErrors(returnContainer["bkg"].size(),&AvgX_pt[0],&(returnContainer["bkg"])[0],0,&(returnContainer["bkg_errors"])[0]);
    TGraphErrors* v2bkg_nq     = new TGraphErrors(returnContainer["bkg_nq"].size(),&AvgX_pt[0],&(returnContainer["bkg_nq"])[0],0,&(returnContainer["bkg_errors_nq"])[0]);
    TGraphErrors* v2bkg_ket    = new TGraphErrors(returnContainer["bkg"].size(),&AvgX_ket[0],&(returnContainer["bkg"])[0],0,&(returnContainer["bkg_errors"])[0]);
    TGraphErrors* v2bkg_ket_nq = new TGraphErrors(returnContainer["bkg_nq"].size(),&AvgX_ket_nq[0],&(returnContainer["bkg_nq"])[0],0,&(returnContainer["bkg_errors_nq"])[0]);

    if(V0ID == "Kshort")
    {
        v2           -> Write("v2kshort",TObject::kOverwrite);
        v2_nq        -> Write("v2kshort_nq",TObject::kOverwrite);
        v2_ket       -> Write("v2kshort_ket",TObject::kOverwrite);
        v2_ket_nq    -> Write("v2kshort_ket_nq",TObject::kOverwrite);
        v2obs        -> Write("v2obskshort",TObject::kOverwrite);
        v2obs_nq     -> Write("v2obskshort_nq",TObject::kOverwrite);
        v2obs_ket    -> Write("v2obskshort_ket",TObject::kOverwrite);
        v2obs_ket_nq -> Write("v2obskshort_ket_nq",TObject::kOverwrite);
        v2bkg        -> Write("v2bkgkshort",TObject::kOverwrite);
        v2bkg_nq     -> Write("v2bkgkshort_nq",TObject::kOverwrite);
        v2bkg_ket    -> Write("v2bkgkshort_ket",TObject::kOverwrite);
        v2bkg_ket_nq -> Write("v2bkgkshort_ket_nq",TObject::kOverwrite);
    }
    if(V0ID == "Lambda")
    {
        v2           -> Write("v2lambda",TObject::kOverwrite);
        v2_nq        -> Write("v2lambda_nq",TObject::kOverwrite);
        v2_ket       -> Write("v2lambda_ket",TObject::kOverwrite);
        v2_ket_nq    -> Write("v2lambda_ket_nq",TObject::kOverwrite);
        v2obs        -> Write("v2obslambda",TObject::kOverwrite);
        v2obs_nq     -> Write("v2obslambda_nq",TObject::kOverwrite);
        v2obs_ket    -> Write("v2obslambda_ket",TObject::kOverwrite);
        v2obs_ket_nq -> Write("v2obslambda_ket_nq",TObject::kOverwrite);
        v2bkg        -> Write("v2bkglambda",TObject::kOverwrite);
        v2bkg_nq     -> Write("v2bkglambda_nq",TObject::kOverwrite);
        v2bkg_ket    -> Write("v2bkglambda_ket",TObject::kOverwrite);
        v2bkg_ket_nq -> Write("v2bkglambda_ket_nq",TObject::kOverwrite);
    }

    out->Close();
}


void V0vnFit()
{
    double maxFact = 2; //2
    double minFact = 0.005; //0.005
    double mini = -0.005;
    double maxi = 0.06;
    int numFourierParams = 5;
    bool Peak = true;
	//bool Peak = false;
	//Aesthetics

    //TLatex labels
    std::ostringstream os; // stringstream for making dynamic TLatex labels
    double SNN = 8.16;
    int Lint = 35;
    int Nmin = 185;
    int Nmax = 220;
    int pTassMin = 1;
    int pTassMax = 3;
    int longRange = 2;

    //For Enabling TLatex labels
    bool publish = true;
    //bool publish = kFALSE;

    gStyle->SetOptFit(1111);
    gStyle->SetErrorX(0); //removes horizontal error bars

    //HIST APPEARANCE
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    TGaxis::SetMaxDigits(3); // Forces exponents after n number of digits

    //Font and Size
    gStyle->SetTextSize(20);
    gStyle->SetTextFont(42); //2=times-bold-r-normal, 2=precision for TLatex to work

	//Files
	//TFiles
    //TFile *f = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationJL7_8.root");
    //TFile *f = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationTotal_08_20_2017.root");
    //TFile *f = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityTotal_08_21_2017.root");
    //TFile *f = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/LooseAndTight/V0CorrelationTightMCTotal_08_23_2017.root");

    TFile *f = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityCorrectMultB_09_19_17.root");
    TFile *fhad = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationRapidityTotal_08_20_2017.root" );
    //TFile *f = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0CorrelationRapidityClosureNoEff_09_11_17.root"); //No Eff closure for D0 study
    //TFile *f = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0CorrelationClosureTotal_08_28_2017.root"); //Normal Closure for D0 study
    //TFile *f = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/MatchV0ClosureBpPb_09_20_17.root"); //Match Closure for D0 study
    //TFile *fhad = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0CorrelationClosureHadronRecoFix_09_10_17.root"); //For vn of hadron in closure

	//Txt files
	ofstream vnPeak;
	ofstream vnSide;
    ofstream vnH;
    ofstream vnCalculator;
    std::string vnPeakName = "vnPeakMatch.txt";
    std::string vnSideName = "vnPeakMatch.txt";
    std::string vnHName = "vnHadron.txt";
    std::string vnCalculatorName = "vnSignalMatch.txt";
    std::string branchname_ks = "v0CorrelationRapidity/Ptkshort_pt";
    std::string branchname_ks_bkg = "v0CorrelationRapidity/Ptkshort_bkg_pt";
    std::string branchname_la = "v0CorrelationRapidity/Ptlambda_pt";
    std::string branchname_la_bkg = "v0CorrelationRapidity/Ptlambda_bkg_pt";
    std::string branchname_ket_ks = "v0CorrelationRapidity/KETkshort_pt";
    std::string branchname_ket_ks_bkg = "v0CorrelationRapidity/KETkshort_bkg_pt";
    std::string branchname_ket_la = "v0CorrelationRapidity/KETlambda_pt";
    std::string branchname_ket_la_bkg = "v0CorrelationRapidity/KETlambda_bkg_pt";
    std::string graphName = "v2valuesRapidity_etaGap_0p75.root";
	vnPeak.open(vnPeakName);
	vnSide.open(vnSideName);
    vnH.open(vnHName);
    vnCalculator.open(vnCalculatorName);

    TVirtualFitter::SetMaxIterations(300000);
    TH1::SetDefaultSumw2();
    std::vector<double> PtBin_ks = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5};//, 10.0, 15.0, 20.0}; //if number of bins changes make sure you change numPtBins
    std::vector<double> PtBin_la = {0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5};//, 10.0, 15.0, 20.0}; //if number of bins changes make sure you change numPtBins
    int numPtBins_ks = PtBin_ks.size() - 1;
    int numPtBins_la = PtBin_la.size() - 1;
    TH1D* dPhiFourierPeak_ks[numPtBins_ks];
    TH1D* dPhiFourierSide_ks[numPtBins_ks];
    TH1D* dPhiFourierPeak_la[numPtBins_la];
    TH1D* dPhiFourierSide_la[numPtBins_la];

    TH1D* dPhiPeak_ks[numPtBins_ks];
    TH1D* dPhiSide_ks[numPtBins_ks];
    TH1D* dPhiPeak_la[numPtBins_la];
    TH1D* dPhiSide_la[numPtBins_la];

    TF1* FourierFit_ks[numPtBins_ks];
    TF1* FourierFit_la[numPtBins_la];
    std::map<int, std::vector<double> > vnValues_ks_peak;
    std::map<int, std::vector<double> > vnErrors_ks_peak;

    std::map<int,std::vector<double> > vnValues_ks_side;
    std::map<int,std::vector<double> > vnErrors_ks_side;

    std::map<int, std::vector<double> > vnValues_la_peak;
    std::map<int, std::vector<double> > vnValues_la_side;

    std::map<int, std::vector<double> > vnErrors_la_peak;
    std::map<int, std::vector<double> > vnErrors_la_side;

    /**
     * Initialized first two elements because I want to keep the index consistently starting at 2
     **/
    std::vector<double> vnValues_h = {-999,-999};
    std::vector<double> vnErrors_h = {-999,-999};

    //Fsig for vn calculations
    //std::vector<double> fsig_ks = {0.999666 ,0.999977 ,0.999972 ,0.999988 ,0.999998 ,0.999999 ,0.999999 ,0.999992 ,0.999994 ,0.999955, 0.999975};
    //std::vector<double> fsig_la = {0.988877 ,0.9967 ,0.99754 ,0.998939 ,0.999954 ,0.999951 ,0.999992 ,0.999842};
    //std::vector<double> fsig_ks = {0.999476 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,0.999999 ,0.999999 ,0.999988 ,0.999997 ,0.992327};// ,0.99836};// ,0.890106 ,0.481433};
    //std::vector<double> fsig_la = {0.99882 ,0.999987 ,1 ,0.999524 ,0.999632 ,0.999855 ,0.999698 ,0.998783 ,0.999771 ,0.997088};// ,0.92718};// ,0.990913 ,0.421011};
    std::vector<double> fsig_ks = {0.991217 ,1 ,1 ,1 ,1 ,1 ,1 ,0.999988 ,0.99998 ,0.998263 ,0.942699 ,0.999827 ,0.99997};// ,0.99836};// ,0.890106 ,0.481433}; //MC
    std::vector<double> fsig_la = {0.971604 ,0.985982 ,0.991208 ,0.989859 ,0.997168 ,0.996337 ,0.982485 ,0.953336 ,0.781606 ,0.596769};// ,0.92718};// ,0.990913 ,0.421011}; //MC


    if((PtBin_ks.size()-1 != fsig_ks.size()) || (PtBin_la.size()-1 != fsig_la.size()))
    {
        cout << "something is wrong with number of pt bins or number of fsig" << endl;
        return;
    }

    TLatex* ltx2 = new TLatex();
    ltx2->SetTextSize(0.045);
    ltx2->SetNDC(kTRUE);

    //FITTING
	//KSHORT
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
	cout << "KSHORT KSHORT KSHORT"                                                             << endl;
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
    for(int i=0; i<numPtBins_ks; i++)
    {


        //================================================================================
        //Peak Calculations
        //================================================================================
        dPhiPeak_ks[i] = new TH1D(Form("dPhiPeak_ks%d",i), "K_{S}^{0} - h^{#pm} ", 31, -(0.5 -1.0/32)*PI, (1.5 - 1.0/32)*PI);
        TH1D *dPhiHad = new TH1D("dPhiHad", "h^{#pm}- h^{#pm} ", 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        //Pull 2D Histograms
        TH2D *hbackgroundPeak = (TH2D*) f->Get(Form("v0CorrelationRapidity/backgroundkshort_pt%d",i));
        TH2D *hsignalPeak     = (TH2D*) f->Get(Form("v0CorrelationRapidity/signalkshort_pt%d",i));
        TH2D *hBackgroundHad  = (TH2D*) fhad->Get("xiCorrelationRapidity/BackgroundHad");
        TH2D *hSignalHad      = (TH2D*) fhad->Get("xiCorrelationRapidity/SignalHad");
        //TH2D *hBackgroundHad  = (TH2D*) fhad->Get("HadronCorrelation/BackgroundHadReco");
        //TH2D *hSignalHad      = (TH2D*) fhad->Get("HadronCorrelation/SignalHadReco");

        //Project Phi

        // For projecting both shoulders
        //TH1D* hbPhiTotPeak = hbackgroundPeak->ProjectionY("PhiBkgTotPeak", 1, 10);
        //TH1D* hbPhiOthPeak = hbackgroundPeak->ProjectionY("PhiBkgOthPeak", 23, -1);
        //TH1D* hsPhiTotPeak = hsignalPeak->ProjectionY("PhiSigTotPeak", 1, 10);
        //TH1D* hsPhiOthPeak = hsignalPeak->ProjectionY("PhiSigOthPeak", 23, -1);
        //TH1D* hbHadPhiTot = hBackgroundHad->ProjectionY("PhiBkgHadTot", 1, 10);
        //TH1D* hbHadPhiOth = hBackgroundHad->ProjectionY("PhiBkgHadOth", 23, -1);
        //TH1D* hsHadPhiTot = hSignalHad->ProjectionY("PhiSigHadTot", 1, 10);
        //TH1D* hsHadPhiOth = hSignalHad->ProjectionY("PhiSigHadOth", 23, -1);

        int binlow = 14;
        int binhigh = 19;
        TH1D* hbPhiTotPeak = hbackgroundPeak->ProjectionY("PhiBkgTotPeak", 1, binlow);
        TH1D* hbPhiOthPeak = hbackgroundPeak->ProjectionY("PhiBkgOthPeak", binhigh, -1);
        TH1D* hsPhiTotPeak = hsignalPeak->ProjectionY("PhiSigTotPeak", 1, binlow);
        TH1D* hsPhiOthPeak = hsignalPeak->ProjectionY("PhiSigOthPeak", binhigh, -1);
        TH1D* hbHadPhiTot = hBackgroundHad->ProjectionY("PhiBkgHadTot", 1, binlow);
        TH1D* hbHadPhiOth = hBackgroundHad->ProjectionY("PhiBkgHadOth", binhigh, -1);
        TH1D* hsHadPhiTot = hSignalHad->ProjectionY("PhiSigHadTot", 1, binlow);
        TH1D* hsHadPhiOth = hSignalHad->ProjectionY("PhiSigHadOth", binhigh, -1);

        hbPhiTotPeak->Add(hbPhiOthPeak);
        hsPhiTotPeak->Add(hsPhiOthPeak);

        hbHadPhiTot->Add(hbHadPhiOth);
        hsHadPhiTot->Add(hsHadPhiOth);

        //Divide
        dPhiPeak_ks[i]->Divide(hsPhiTotPeak, hbPhiTotPeak);
        dPhiHad->Divide(hsHadPhiTot, hbHadPhiTot);

        //Clone histograms for display without fit functions
        dPhiFourierPeak_ks[i] = (TH1D*)dPhiPeak_ks[i]->Clone();
        TH1D* dPhiHadFourier = (TH1D*)dPhiHad->Clone();

        FourierFit_ks[i] = new TF1(Form("FourierFit_ks%d",i), FourierHad, -1.5, 5, numFourierParams);
        FourierFit_ks[i]->SetNpx(250);
        FourierFit_ks[i]->SetParNames("Scale", "V_{1}", "V_{2}", "V_{3}");

        TCanvas *c4_ks_peak = new TCanvas("c4_ks_peak", "Fourier peak trg_ks", 800,800);
        c4_ks_peak->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        dPhiFourierPeak_ks[i]->Fit(Form("FourierFit_ks%d",i));
        dPhiFourierPeak_ks[i]->SetStats(kFALSE);
        for(int j=2; j<numFourierParams; j++)
        {
            vnValues_ks_peak[j].push_back(FourierFit_ks[i]->GetParameter(j));
            vnErrors_ks_peak[j].push_back(FourierFit_ks[i]->GetParError(j));
        }

        TF1 *FourierFitHad = new TF1("FourierFitHad", FourierHad, -1.5, 5, numFourierParams);
        FourierFitHad->SetNpx(250);
        FourierFitHad->SetParNames("Scale", "V_{1}", "V_{2}", "V_{3}");

        TCanvas *c5_ks = new TCanvas("c5_ks", "c5_ks", 800,800);
        c5_ks->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        dPhiHadFourier->Fit("FourierFitHad");
        dPhiHadFourier->SetStats(kFALSE);
        if(i==0)
        {
            for(int j=2; j<numFourierParams; j++)
            {
                vnValues_h.push_back(FourierFitHad->GetParameter(j));
                vnErrors_h.push_back(FourierFitHad->GetParError(j));
            }
        }

        double maxBinContent = dPhiFourierPeak_ks[i]->GetBinContent(dPhiFourierPeak_ks[i]->GetMaximumBin());
        double minBinContent = dPhiFourierPeak_ks[i]->GetBinContent(dPhiFourierPeak_ks[i]->GetMinimumBin());
        double minRange = minBinContent - minFact*minBinContent;
        double maxRange = minRange + maxFact*(maxBinContent - minBinContent);

        //ZYAM FITS
        TCanvas *c2_ks_peak = new TCanvas("c2_ks_peak", "ZYAM Peak Trg_ks", 800,800);
        c2_ks_peak->cd();

        gPad->SetTickx();
        gPad->SetTicky();
        dPhiPeak_ks[i]->SetMarkerStyle(21);
        dPhiPeak_ks[i]->SetMarkerColor(4);
        dPhiPeak_ks[i]->SetTitleOffset(2, "Y");
        dPhiPeak_ks[i]->SetTitle("Peak");
        dPhiPeak_ks[i]->GetYaxis()->SetRangeUser(mini , maxi);
        dPhiPeak_ks[i]->GetYaxis()->SetTitleSize(0.03);
        dPhiPeak_ks[i]->GetYaxis()->CenterTitle(true);
        dPhiPeak_ks[i]->GetYaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} ");
        dPhiPeak_ks[i]->SetTitleOffset(1.5, "X");
        dPhiPeak_ks[i]->GetXaxis()->SetTitleSize(0.035);
        dPhiPeak_ks[i]->GetXaxis()->CenterTitle(true);
        dPhiPeak_ks[i]->GetXaxis()->SetTitle("#Delta#phi (radians)");
        //dPhiPeak_ks[i]->Draw("E1");
        dPhiPeak_ks[i]->Fit("pol2","","", 0.4,2.4);
        dPhiPeak_ks[i]->SetStats(!publish);

        TLatex *ltx3 = new TLatex();
        ltx3->SetTextSize(0.035);
        ltx3->SetNDC(kTRUE);
        ltx3->SetTextFont(42);

        if(publish)
        {
            os << "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = " << fixed << std::setprecision(2) << SNN << " TeV";
            ltx3->DrawLatex(0.2, 0.82, os.str().c_str());
            os.str(std::string());
            os << "L_{#lower[-0.25]{int}} = " << Lint << " nb^{-1}";
            ltx3->DrawLatex(0.2, 0.74, os.str().c_str());
            os.str(std::string());
            os << Nmin << "  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< " << Nmax;
            ltx3->DrawLatex(0.2, 0.67, os.str().c_str());
            os.str(std::string());
            os << pTassMin << " < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < " << pTassMax << " GeV";
            ltx3->DrawLatex(0.2, 0.60, os.str().c_str());
            os.str(std::string());
            os << "Long range (|#Delta#eta| > " <<  longRange << ")";
            ltx3->DrawLatex(0.2, 0.53, os.str().c_str());
            os.str(std::string());
            ltx2->DrawLatex(0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }

        maxBinContent = dPhiHadFourier->GetBinContent(dPhiHadFourier->GetMaximumBin());
        minBinContent = dPhiHadFourier->GetBinContent(dPhiHadFourier->GetMinimumBin());
        minRange = minBinContent - minFact*minBinContent;
        maxRange = minRange + maxFact*(maxBinContent - minBinContent);

        TCanvas *c3_ks = new TCanvas("c3_ks", "ZYAM Hadronic", 800,800);
        c3_ks->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        dPhiHad->SetMarkerStyle(34);
        dPhiHad->SetMarkerSize(1.5);
        dPhiHad->SetTitleOffset(2, "Y");
        dPhiHad->SetTitle("");
        dPhiHad->GetYaxis()->SetRangeUser(mini , maxi);
        dPhiHad->GetYaxis()->SetTitleSize(0.03);
        dPhiHad->GetYaxis()->CenterTitle(true);
        dPhiHad->GetYaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi}");
        dPhiHad->SetTitleOffset(1.5, "X");
        dPhiHad->GetXaxis()->SetTitleSize(0.035);
        dPhiHad->GetXaxis()->CenterTitle(true);
        dPhiHad->GetXaxis()->SetTitle("#Delta#phi (radians)");
        dPhiHad->Fit("pol2","","", 0.4,2);
        dPhiHad->SetStats(!publish);

        if(publish)
        {
            ltx3->DrawLatex(0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
            ltx3->DrawLatex(0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}");
            ltx3->DrawLatex(0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
            ltx3->DrawLatex(0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
            ltx3->DrawLatex(0.2, 0.53, "Long range (|#Delta#eta| > 2)");
            ltx2->DrawLatex(0.7, 0.74, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }

        TF1 *dPhiFitPeak = dPhiPeak_ks[i]->GetFunction("pol2");
        TF1 *dPhiHadFit = dPhiHad->GetFunction("pol2");

        double dPhiFitMinPeak = dPhiFitPeak->GetMinimum();
        double dPhiHadFitMin = dPhiHadFit->GetMinimum();

        for(int j = 1; j < 32; j++)
        {
            dPhiFourierPeak_ks[i]->AddBinContent(j, -dPhiFitMinPeak);
        }

        for(int j = 1; j < 32; j++)
        {
            dPhiHadFourier->AddBinContent(j, -dPhiHadFitMin);
        }

        c4_ks_peak->cd();
        dPhiFourierPeak_ks[i]->SetMarkerStyle(21);
        dPhiFourierPeak_ks[i]->SetMarkerColor(4);
        //dPhiFourierPeak_ks[i]->Draw("E1");
        dPhiFourierPeak_ks[i]->SetStats(kFALSE);
        dPhiFourierPeak_ks[i]->Fit(Form("FourierFit_ks%d",i));
        os << "Peak " << PtBin_ks[i] << "_Pt_" << PtBin_ks[i+1];
        dPhiFourierPeak_ks[i]->SetTitle(os.str().c_str());
        os.str(std::string());
        dPhiFourierPeak_ks[i]->SetTitleOffset(2, "Y");
        dPhiFourierPeak_ks[i]->GetYaxis()->CenterTitle(true);
        dPhiFourierPeak_ks[i]->GetYaxis()->SetTitleSize(0.03);
        dPhiFourierPeak_ks[i]->GetYaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}");
        //dPhiFourierPeak_ks[i]->GetYaxis()->SetRangeUser(-0.0004, 0.008);
        dPhiFourierPeak_ks[i]->GetYaxis()->SetRangeUser(mini, maxi);
        dPhiFourierPeak_ks[i]->SetTitleOffset(1.5, "X");
        dPhiFourierPeak_ks[i]->GetXaxis()->SetTitleSize(0.035);
        dPhiFourierPeak_ks[i]->GetXaxis()->CenterTitle(true);
        dPhiFourierPeak_ks[i]->GetXaxis()->SetTitle("#Delta#phi (radians)");


        if(publish)
        {
            ltx3->DrawLatex(0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
            ltx3->DrawLatex(0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}");
            ltx3->DrawLatex(0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
            ltx3->DrawLatex(0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
            ltx3->DrawLatex(0.2, 0.53, "Long range (|#Delta#eta| > 2)");
            ltx2->DrawLatex(0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }

        c5_ks->cd();
        dPhiHadFourier->SetMarkerStyle(34);
        dPhiHadFourier->SetMarkerSize(1.5);
        dPhiHadFourier->Draw("E1");
        dPhiHadFourier->SetStats(kFALSE);
        dPhiHadFourier->Fit("FourierFitHad");
        dPhiHadFourier->SetTitle("");
        dPhiHadFourier->SetTitleOffset(2, "Y");
        dPhiHadFourier->GetYaxis()->CenterTitle(true);
        dPhiHadFourier->GetYaxis()->SetTitleSize(0.03);
        dPhiHadFourier->GetYaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}");
        dPhiHadFourier->GetYaxis()->SetRangeUser(-0.0004 , 0.008);
        dPhiHadFourier->SetTitleOffset(1.5, "X");
        dPhiHadFourier->GetXaxis()->SetTitleSize(0.035);
        dPhiHadFourier->GetXaxis()->CenterTitle(true);
        dPhiHadFourier->GetXaxis()->SetTitle("#Delta#phi (radians)");

        if(publish)
        {
            ltx3->DrawLatex(0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
            ltx3->DrawLatex(0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}");
            ltx3->DrawLatex(0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
            ltx3->DrawLatex(0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
            ltx3->DrawLatex(0.2, 0.53, "Long range (|#Delta#eta| > 2)");
            ltx2->DrawLatex(0.7, 0.74, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }

        //================================================================================
        //Sideband Calculations
        //================================================================================
        dPhiSide_ks[i] = new TH1D(Form("dPhiSide_ks%d",i), "K_{S}^{0} - h^{#pm} ", 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        TH2D *hbackgroundSide = (TH2D*) f->Get(Form("v0CorrelationRapidity/backgroundkshort_bkg_pt%d",i));
        TH2D *hsignalSide     = (TH2D*) f->Get(Form("v0CorrelationRapidity/signalkshort_bkg_pt%d",i));

        //Project Phi
        TH1D* hbPhiTotSide = hbackgroundSide->ProjectionY("PhiBkgTot", 0, 10);
        TH1D* hbPhiOthSide = hbackgroundSide->ProjectionY("PhiBkgOthPeak", 23, -1);
        TH1D* hsPhiTotSide = hsignalSide->ProjectionY("PhiSigTot", 0, 10);
        TH1D* hsPhiOthSide = hsignalSide->ProjectionY("PhiSigOthPeak", 23, -1);

        hbPhiTotSide->Add(hbPhiOthSide);
        hsPhiTotSide->Add(hsPhiOthSide);

        //Divide
        dPhiSide_ks[i]->Divide(hsPhiTotSide, hbPhiTotSide);

        //Clone histograms for display without fit functions
        dPhiFourierSide_ks[i] = (TH1D*)dPhiSide_ks[i]->Clone();

        FourierFit_ks[i] = new TF1(Form("FourierFit_ks%d",i) , FourierHad, -1.5, 5, numFourierParams);
        FourierFit_ks[i]->SetNpx(250);
        FourierFit_ks[i]->SetParNames("Scale", "V_{1}", "V_{2}", "V_{3}");

        TCanvas *c4_ks_side = new TCanvas("c4_ks_side", "Fourier side trg_ks", 800,800);
        c4_ks_side->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        //dPhiFourierSide_ks->Fit("FourierFit_ks","","",0,PI);
        dPhiFourierSide_ks[i]->Fit(Form("FourierFit_ks%d",i));
        dPhiFourierSide_ks[i]->SetStats(kFALSE);
        for(int j=2; j<numFourierParams; j++)
        {
            vnValues_ks_side[j].push_back(FourierFit_ks[i]->GetParameter(j));
            vnErrors_ks_side[j].push_back(FourierFit_ks[i]->GetParError(j));
        }

        maxBinContent = dPhiFourierSide_ks[i]->GetBinContent(dPhiFourierSide_ks[i]->GetMaximumBin());
        minBinContent = dPhiFourierSide_ks[i]->GetBinContent(dPhiFourierSide_ks[i]->GetMinimumBin());
        minRange = minBinContent - minFact*minBinContent;
        maxRange = minRange + maxFact*(maxBinContent - minBinContent);

        //ZYAM FITS
        TCanvas *c2_ks_side = new TCanvas("c2_ks_side", "ZYAM Side trg_ks", 800,800);
        c2_ks_side->cd();

        gPad->SetTickx();
        gPad->SetTicky();
        dPhiSide_ks[i]->SetMarkerStyle(21);
        dPhiSide_ks[i]->SetMarkerColor(4);
        dPhiSide_ks[i]->SetTitleOffset(2, "Y");
        dPhiSide_ks[i]->SetTitle("Sideband");
        dPhiSide_ks[i]->GetYaxis()->SetRangeUser(mini , maxi);
        dPhiSide_ks[i]->GetYaxis()->SetTitleSize(0.03);
        dPhiSide_ks[i]->GetYaxis()->CenterTitle(true);
        dPhiSide_ks[i]->GetYaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} ");
        dPhiSide_ks[i]->SetTitleOffset(1.5, "X");
        dPhiSide_ks[i]->GetXaxis()->SetTitleSize(0.035);
        dPhiSide_ks[i]->GetXaxis()->CenterTitle(true);
        dPhiSide_ks[i]->GetXaxis()->SetTitle("#Delta#phi (radians)");
        dPhiSide_ks[i]->Fit("pol2","","", 0.4,2.4);
        dPhiSide_ks[i]->SetStats(!publish);

        if(publish)
        {
            os << "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = " << fixed << std::setprecision(2) << SNN << " TeV";
            ltx3->DrawLatex(0.2, 0.82, os.str().c_str());
            os.str(std::string());
            os << "L_{#lower[-0.25]{int}} = " << Lint << " nb^{-1}";
            ltx3->DrawLatex(0.2, 0.74, os.str().c_str());
            os.str(std::string());
            os << Nmin << "  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< " << Nmax;
            ltx3->DrawLatex(0.2, 0.67, os.str().c_str());
            os.str(std::string());
            os << pTassMin << " < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < " << pTassMax << " GeV";
            ltx3->DrawLatex(0.2, 0.60, os.str().c_str());
            os.str(std::string());
            os << "Long range (|#Delta#eta| > " <<  longRange << ")";
            ltx3->DrawLatex(0.2, 0.53, os.str().c_str());
            os.str(std::string());
            ltx2->DrawLatex(0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }

        TF1 *dPhiFitSide = dPhiSide_ks[i]->GetFunction("pol2");

        double dPhiFitMinSide = dPhiFitSide->GetMinimum();

        for(int j = 1; j < 32; j++)
        {
            dPhiFourierSide_ks[i]->AddBinContent(j, -dPhiFitMinSide);
        }

        c4_ks_side->cd();
        dPhiFourierSide_ks[i]->SetMarkerStyle(21);
        dPhiFourierSide_ks[i]->SetMarkerColor(4);
        //dPhiFourierSide_ks[i]->Draw("E1");
        dPhiFourierSide_ks[i]->SetStats(kFALSE);
        dPhiFourierSide_ks[i]->Fit(Form("FourierFit_ks%d",i));
        os << "SideBand " << PtBin_ks[i] << "_Pt_" << PtBin_ks[i+1];
        dPhiFourierSide_ks[i]->SetTitle(os.str().c_str());
        os.str(std::string());
        dPhiFourierSide_ks[i]->SetTitleOffset(2, "Y");
        dPhiFourierSide_ks[i]->GetYaxis()->CenterTitle(true);
        dPhiFourierSide_ks[i]->GetYaxis()->SetTitleSize(0.03);
        dPhiFourierSide_ks[i]->GetYaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}");
        //dPhiFourierSide_ks[i]->GetYaxis()->SetRangeUser(-0.0004, 0.008);
        dPhiFourierSide_ks[i]->GetYaxis()->SetRangeUser(mini, maxi);
        dPhiFourierSide_ks[i]->SetTitleOffset(1.5, "X");
        dPhiFourierSide_ks[i]->GetXaxis()->SetTitleSize(0.035);
        dPhiFourierSide_ks[i]->GetXaxis()->CenterTitle(true);
        dPhiFourierSide_ks[i]->GetXaxis()->SetTitle("#Delta#phi (radians)");

        if(publish)
        {
            ltx3->DrawLatex(0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
            ltx3->DrawLatex(0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}");
            ltx3->DrawLatex(0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
            ltx3->DrawLatex(0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
            ltx3->DrawLatex(0.2, 0.53, "Long range (|#Delta#eta| > 2)");
            ltx2->DrawLatex(0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }
    }

	//LAMBDA
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
	cout << "LAMBDA LAMBDA LAMBDA"                                                             << endl;
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
	cout << "================================================================================" << endl;
    for(int i=0; i<numPtBins_la; i++)
    {
        //================================================================================
        //Peak Calculations
        //================================================================================
        dPhiPeak_la[i] = new TH1D(Form("dPhiPeak_la%d",i), "K_{S}^{0} - h^{#pm} ", 31, -(0.5 -1.0/32)*PI, (1.5 - 1.0/32)*PI);
        TH1D *dPhiHad = new TH1D("dPhiHad", "h^{#pm}- h^{#pm} ", 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        //Pull 2D Histograms
        TH2D *hbackgroundPeak = (TH2D*) f->Get(Form("v0CorrelationRapidity/backgroundlambda_pt%d",i));
        TH2D *hsignalPeak     = (TH2D*) f->Get(Form("v0CorrelationRapidity/signallambda_pt%d",i));
        TH2D *hBackgroundHad  = (TH2D*) fhad->Get("xiCorrelationRapidity/BackgroundHad");
        TH2D *hSignalHad      = (TH2D*) fhad->Get("xiCorrelationRapidity/SignalHad");
        //TH2D *hBackgroundHad  = (TH2D*) fhad->Get("HadronCorrelation/BackgroundHadReco");
        //TH2D *hSignalHad      = (TH2D*) fhad->Get("HadronCorrelation/SignalHadReco");

        //Project Phi

        // For projecting both shoulders
        //TH1D* hbPhiTotPeak = hbackgroundPeak->ProjectionY("PhiBkgTotPeak", 1, 10);
        //TH1D* hbPhiOthPeak = hbackgroundPeak->ProjectionY("PhiBkgOthPeak", 23, -1);
        //TH1D* hsPhiTotPeak = hsignalPeak->ProjectionY("PhiSigTotPeak", 1, 10);
        //TH1D* hsPhiOthPeak = hsignalPeak->ProjectionY("PhiSigOthPeak", 23, -1);
        //TH1D* hbHadPhiTot = hBackgroundHad->ProjectionY("PhiBkgHadTot", 1, 10);
        //TH1D* hbHadPhiOth = hBackgroundHad->ProjectionY("PhiBkgHadOth", 23, -1);
        //TH1D* hsHadPhiTot = hSignalHad->ProjectionY("PhiSigHadTot", 1, 10);
        //TH1D* hsHadPhiOth = hSignalHad->ProjectionY("PhiSigHadOth", 23, -1);

        int binlow = 14;
        int binhigh = 19;
        TH1D* hbPhiTotPeak = hbackgroundPeak->ProjectionY("PhiBkgTotPeak", 1, binlow);
        TH1D* hbPhiOthPeak = hbackgroundPeak->ProjectionY("PhiBkgOthPeak", binhigh, -1);
        TH1D* hsPhiTotPeak = hsignalPeak->ProjectionY("PhiSigTotPeak", 1, binlow);
        TH1D* hsPhiOthPeak = hsignalPeak->ProjectionY("PhiSigOthPeak", binhigh, -1);
        TH1D* hbHadPhiTot = hBackgroundHad->ProjectionY("PhiBkgHadTot", 1, binlow);
        TH1D* hbHadPhiOth = hBackgroundHad->ProjectionY("PhiBkgHadOth", binhigh, -1);
        TH1D* hsHadPhiTot = hSignalHad->ProjectionY("PhiSigHadTot", 1, binlow);
        TH1D* hsHadPhiOth = hSignalHad->ProjectionY("PhiSigHadOth", binhigh, -1);

        hbPhiTotPeak->Add(hbPhiOthPeak);
        hsPhiTotPeak->Add(hsPhiOthPeak);

        hbHadPhiTot->Add(hbHadPhiOth);
        hsHadPhiTot->Add(hsHadPhiOth);

        //Divide
        dPhiPeak_la[i]->Divide(hsPhiTotPeak, hbPhiTotPeak);
        dPhiHad->Divide(hsHadPhiTot, hbHadPhiTot);

        //Clone histograms for display without fit functions
        dPhiFourierPeak_la[i] = (TH1D*)dPhiPeak_la[i]->Clone();
        TH1D* dPhiHadFourier = (TH1D*)dPhiHad->Clone();

        FourierFit_la[i] = new TF1(Form("FourierFit_la%d",i), FourierHad, -1.5, 5, numFourierParams);
        FourierFit_la[i]->SetNpx(250);
        FourierFit_la[i]->SetParNames("Scale", "V_{1}", "V_{2}", "V_{3}");

        TCanvas *c4_la_peak = new TCanvas("c4_la_peak", "Fourier peak trg_la", 800,800);
        c4_la_peak->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        dPhiFourierPeak_la[i]->Fit(Form("FourierFit_la%d",i));
        dPhiFourierPeak_la[i]->SetStats(kFALSE);
        for(int j=2; j<numFourierParams; j++)
        {
            vnValues_la_peak[j].push_back(FourierFit_la[i]->GetParameter(j));
            vnErrors_la_peak[j].push_back(FourierFit_la[i]->GetParError(j));
        }

        TF1 *FourierFitHad = new TF1("FourierFitHadLa", FourierHad, -1.5, 5, numFourierParams);
        FourierFitHad->SetNpx(250);
        FourierFitHad->SetParNames("Scale", "V_{1}", "V_{2}", "V_{3}");

        TCanvas *c5_la = new TCanvas("c5_la", "c5_la", 800,800);
        c5_la->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        dPhiHadFourier->Fit("FourierFitHadLa");
        dPhiHadFourier->SetStats(kFALSE);

        double maxBinContent = dPhiFourierPeak_la[i]->GetBinContent(dPhiFourierPeak_la[i]->GetMaximumBin());
        double minBinContent = dPhiFourierPeak_la[i]->GetBinContent(dPhiFourierPeak_la[i]->GetMinimumBin());
        double minRange = minBinContent - minFact*minBinContent;
        double maxRange = minRange + maxFact*(maxBinContent - minBinContent);

        //ZYAM FITS
        TCanvas *c2_la_peak = new TCanvas("c2_la_peak", "ZYAM Peak Trg_la", 800,800);
        c2_la_peak->cd();

        gPad->SetTickx();
        gPad->SetTicky();
        dPhiPeak_la[i]->SetMarkerStyle(21);
        dPhiPeak_la[i]->SetMarkerColor(4);
        dPhiPeak_la[i]->SetTitleOffset(2, "Y");
        dPhiPeak_la[i]->SetTitle("Peak");
        dPhiPeak_la[i]->GetYaxis()->SetRangeUser(mini , maxi);
        dPhiPeak_la[i]->GetYaxis()->SetTitleSize(0.03);
        dPhiPeak_la[i]->GetYaxis()->CenterTitle(true);
        dPhiPeak_la[i]->GetYaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} ");
        dPhiPeak_la[i]->SetTitleOffset(1.5, "X");
        dPhiPeak_la[i]->GetXaxis()->SetTitleSize(0.035);
        dPhiPeak_la[i]->GetXaxis()->CenterTitle(true);
        dPhiPeak_la[i]->GetXaxis()->SetTitle("#Delta#phi (radians)");
        //dPhiPeak_la[i]->Draw("E1");
        dPhiPeak_la[i]->Fit("pol2","","", 0.4,2.4);
        dPhiPeak_la[i]->SetStats(!publish);

        TLatex *ltx3 = new TLatex();
        ltx3->SetTextSize(0.035);
        ltx3->SetNDC(kTRUE);
        ltx3->SetTextFont(42);

        if(publish)
        {
            os << "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = " << fixed << std::setprecision(2) << SNN << " TeV";
            ltx3->DrawLatex(0.2, 0.82, os.str().c_str());
            os.str(std::string());
            os << "L_{#lower[-0.25]{int}} = " << Lint << " nb^{-1}";
            ltx3->DrawLatex(0.2, 0.74, os.str().c_str());
            os.str(std::string());
            os << Nmin << "  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< " << Nmax;
            ltx3->DrawLatex(0.2, 0.67, os.str().c_str());
            os.str(std::string());
            os << pTassMin << " < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < " << pTassMax << " GeV";
            ltx3->DrawLatex(0.2, 0.60, os.str().c_str());
            os.str(std::string());
            os << "Long range (|#Delta#eta| > " <<  longRange << ")";
            ltx3->DrawLatex(0.2, 0.53, os.str().c_str());
            os.str(std::string());
            ltx2->DrawLatex(0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }

        maxBinContent = dPhiHadFourier->GetBinContent(dPhiHadFourier->GetMaximumBin());
        minBinContent = dPhiHadFourier->GetBinContent(dPhiHadFourier->GetMinimumBin());
        minRange = minBinContent - minFact*minBinContent;
        maxRange = minRange + maxFact*(maxBinContent - minBinContent);

        TCanvas *c3_la = new TCanvas("c3_la", "ZYAM Hadronic", 800,800);
        c3_la->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        dPhiHad->SetMarkerStyle(34);
        dPhiHad->SetMarkerSize(1.5);
        dPhiHad->SetTitleOffset(2, "Y");
        dPhiHad->SetTitle("");
        dPhiHad->GetYaxis()->SetRangeUser(mini , maxi);
        dPhiHad->GetYaxis()->SetTitleSize(0.03);
        dPhiHad->GetYaxis()->CenterTitle(true);
        dPhiHad->GetYaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi}");
        dPhiHad->SetTitleOffset(1.5, "X");
        dPhiHad->GetXaxis()->SetTitleSize(0.035);
        dPhiHad->GetXaxis()->CenterTitle(true);
        dPhiHad->GetXaxis()->SetTitle("#Delta#phi (radians)");
        dPhiHad->Fit("pol2","","", 0.4,2);
        dPhiHad->SetStats(!publish);

        if(publish)
        {
            ltx3->DrawLatex(0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
            ltx3->DrawLatex(0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}");
            ltx3->DrawLatex(0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
            ltx3->DrawLatex(0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
            ltx3->DrawLatex(0.2, 0.53, "Long range (|#Delta#eta| > 2)");
            ltx2->DrawLatex(0.7, 0.74, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }

        TF1 *dPhiFitPeak = dPhiPeak_la[i]->GetFunction("pol2");
        TF1 *dPhiHadFit = dPhiHad->GetFunction("pol2");

        double dPhiFitMinPeak = dPhiFitPeak->GetMinimum();
        double dPhiHadFitMin = dPhiHadFit->GetMinimum();

        for(int j = 1; j < 32; j++)
        {
            dPhiFourierPeak_la[i]->AddBinContent(j, -dPhiFitMinPeak);
        }

        for(int j = 1; j < 32; j++)
        {
            dPhiHadFourier->AddBinContent(j, -dPhiHadFitMin);
        }

        c4_la_peak->cd();
        dPhiFourierPeak_la[i]->SetMarkerStyle(21);
        dPhiFourierPeak_la[i]->SetMarkerColor(4);
        //dPhiFourierPeak_la[i]->Draw("E1");
        dPhiFourierPeak_la[i]->SetStats(kFALSE);
        dPhiFourierPeak_la[i]->Fit(Form("FourierFit_la%d",i));
        os << "Peak " << PtBin_la[i] << "_Pt_" << PtBin_la[i+1];
        dPhiFourierPeak_la[i]->SetTitle(os.str().c_str());
        os.str(std::string());
        dPhiFourierPeak_la[i]->SetTitleOffset(2, "Y");
        dPhiFourierPeak_la[i]->GetYaxis()->CenterTitle(true);
        dPhiFourierPeak_la[i]->GetYaxis()->SetTitleSize(0.03);
        dPhiFourierPeak_la[i]->GetYaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}");
        //dPhiFourierPeak_la[i]->GetYaxis()->SetRangeUser(-0.0004, 0.008);
        dPhiFourierPeak_la[i]->GetYaxis()->SetRangeUser(mini, maxi);
        dPhiFourierPeak_la[i]->SetTitleOffset(1.5, "X");
        dPhiFourierPeak_la[i]->GetXaxis()->SetTitleSize(0.035);
        dPhiFourierPeak_la[i]->GetXaxis()->CenterTitle(true);
        dPhiFourierPeak_la[i]->GetXaxis()->SetTitle("#Delta#phi (radians)");


        if(publish)
        {
            ltx3->DrawLatex(0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
            ltx3->DrawLatex(0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}");
            ltx3->DrawLatex(0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
            ltx3->DrawLatex(0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
            ltx3->DrawLatex(0.2, 0.53, "Long range (|#Delta#eta| > 2)");
            ltx2->DrawLatex(0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }

        c5_la->cd();
        dPhiHadFourier->SetMarkerStyle(34);
        dPhiHadFourier->SetMarkerSize(1.5);
        dPhiHadFourier->Draw("E1");
        dPhiHadFourier->SetStats(kFALSE);
        dPhiHadFourier->Fit("FourierFitHad");
        dPhiHadFourier->SetTitle("");
        dPhiHadFourier->SetTitleOffset(2, "Y");
        dPhiHadFourier->GetYaxis()->CenterTitle(true);
        dPhiHadFourier->GetYaxis()->SetTitleSize(0.03);
        dPhiHadFourier->GetYaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}");
        dPhiHadFourier->GetYaxis()->SetRangeUser(-0.0004 , 0.008);
        dPhiHadFourier->SetTitleOffset(1.5, "X");
        dPhiHadFourier->GetXaxis()->SetTitleSize(0.035);
        dPhiHadFourier->GetXaxis()->CenterTitle(true);
        dPhiHadFourier->GetXaxis()->SetTitle("#Delta#phi (radians)");

        if(publish)
        {
            ltx3->DrawLatex(0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
            ltx3->DrawLatex(0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}");
            ltx3->DrawLatex(0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
            ltx3->DrawLatex(0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
            ltx3->DrawLatex(0.2, 0.53, "Long range (|#Delta#eta| > 2)");
            ltx2->DrawLatex(0.7, 0.74, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }

        //================================================================================
        //Sideband Calculations
        //================================================================================
        dPhiSide_la[i] = new TH1D(Form("dPhiSide_la%d",i), "K_{S}^{0} - h^{#pm} ", 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
        TH2D *hbackgroundSide = (TH2D*) f->Get(Form("v0CorrelationRapidity/backgroundlambda_bkg_pt%d",i));
        TH2D *hsignalSide     = (TH2D*) f->Get(Form("v0CorrelationRapidity/signallambda_bkg_pt%d",i));

        //Project Phi
        TH1D* hbPhiTotSide = hbackgroundSide->ProjectionY("PhiBkgTot", 0, 10);
        TH1D* hbPhiOthSide = hbackgroundSide->ProjectionY("PhiBkgOthPeak", 23, -1);
        TH1D* hsPhiTotSide = hsignalSide->ProjectionY("PhiSigTot", 0, 10);
        TH1D* hsPhiOthSide = hsignalSide->ProjectionY("PhiSigOthPeak", 23, -1);

        hbPhiTotSide->Add(hbPhiOthSide);
        hsPhiTotSide->Add(hsPhiOthSide);

        //Divide
        dPhiSide_la[i]->Divide(hsPhiTotSide, hbPhiTotSide);

        //Clone histograms for display without fit functions
        dPhiFourierSide_la[i] = (TH1D*)dPhiSide_la[i]->Clone();

        FourierFit_la[i] = new TF1(Form("FourierFit_la%d",i) , FourierHad, -1.5, 5, numFourierParams);
        FourierFit_la[i]->SetNpx(250);
        FourierFit_la[i]->SetParNames("Scale", "V_{1}", "V_{2}", "V_{3}");

        TCanvas *c4_la_side = new TCanvas("c4_la_side", "Fourier side trg_la", 800,800);
        c4_la_side->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        //dPhiFourierSide_la->Fit("FourierFit_la","","",0,PI);
        dPhiFourierSide_la[i]->Fit(Form("FourierFit_la%d",i));
        dPhiFourierSide_la[i]->SetStats(kFALSE);
        for(int j=2; j<numFourierParams; j++)
        {
            vnValues_la_side[j].push_back(FourierFit_la[i]->GetParameter(j));
            vnErrors_la_side[j].push_back(FourierFit_la[i]->GetParError(j));
        }

        maxBinContent = dPhiFourierSide_la[i]->GetBinContent(dPhiFourierSide_la[i]->GetMaximumBin());
        minBinContent = dPhiFourierSide_la[i]->GetBinContent(dPhiFourierSide_la[i]->GetMinimumBin());
        minRange = minBinContent - minFact*minBinContent;
        maxRange = minRange + maxFact*(maxBinContent - minBinContent);

        //ZYAM FITS
        TCanvas *c2_la_side = new TCanvas("c2_la_side", "ZYAM Side trg_la", 800,800);
        c2_la_side->cd();

        gPad->SetTickx();
        gPad->SetTicky();
        dPhiSide_la[i]->SetMarkerStyle(21);
        dPhiSide_la[i]->SetMarkerColor(4);
        dPhiSide_la[i]->SetTitleOffset(2, "Y");
        dPhiSide_la[i]->SetTitle("Sideband");
        dPhiSide_la[i]->GetYaxis()->SetRangeUser(mini , maxi);
        dPhiSide_la[i]->GetYaxis()->SetTitleSize(0.03);
        dPhiSide_la[i]->GetYaxis()->CenterTitle(true);
        dPhiSide_la[i]->GetYaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} ");
        dPhiSide_la[i]->SetTitleOffset(1.5, "X");
        dPhiSide_la[i]->GetXaxis()->SetTitleSize(0.035);
        dPhiSide_la[i]->GetXaxis()->CenterTitle(true);
        dPhiSide_la[i]->GetXaxis()->SetTitle("#Delta#phi (radians)");
        dPhiSide_la[i]->Fit("pol2","","", 0.4,2.4);
        dPhiSide_la[i]->SetStats(!publish);

        if(publish)
        {
            os << "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = " << fixed << std::setprecision(2) << SNN << " TeV";
            ltx3->DrawLatex(0.2, 0.82, os.str().c_str());
            os.str(std::string());
            os << "L_{#lower[-0.25]{int}} = " << Lint << " nb^{-1}";
            ltx3->DrawLatex(0.2, 0.74, os.str().c_str());
            os.str(std::string());
            os << Nmin << "  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< " << Nmax;
            ltx3->DrawLatex(0.2, 0.67, os.str().c_str());
            os.str(std::string());
            os << pTassMin << " < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < " << pTassMax << " GeV";
            ltx3->DrawLatex(0.2, 0.60, os.str().c_str());
            os.str(std::string());
            os << "Long range (|#Delta#eta| > " <<  longRange << ")";
            ltx3->DrawLatex(0.2, 0.53, os.str().c_str());
            os.str(std::string());
            ltx2->DrawLatex(0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }

        TF1 *dPhiFitSide = dPhiSide_la[i]->GetFunction("pol2");

        double dPhiFitMinSide = dPhiFitSide->GetMinimum();

        for(int j = 1; j < 32; j++)
        {
            dPhiFourierSide_la[i]->AddBinContent(j, -dPhiFitMinSide);
        }

        c4_la_side->cd();
        dPhiFourierSide_la[i]->SetMarkerStyle(21);
        dPhiFourierSide_la[i]->SetMarkerColor(4);
        //dPhiFourierSide_la[i]->Draw("E1");
        dPhiFourierSide_la[i]->SetStats(kFALSE);
        dPhiFourierSide_la[i]->Fit(Form("FourierFit_la%d",i));
        os << "SideBand " << PtBin_la[i] << "_Pt_" << PtBin_la[i+1];
        dPhiFourierSide_la[i]->SetTitle(os.str().c_str());
        os.str(std::string());
        dPhiFourierSide_la[i]->SetTitleOffset(2, "Y");
        dPhiFourierSide_la[i]->GetYaxis()->CenterTitle(true);
        dPhiFourierSide_la[i]->GetYaxis()->SetTitleSize(0.03);
        dPhiFourierSide_la[i]->GetYaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}");
        //dPhiFourierSide_la[i]->GetYaxis()->SetRangeUser(-0.0004, 0.008);
        dPhiFourierSide_la[i]->GetYaxis()->SetRangeUser(mini, maxi);
        dPhiFourierSide_la[i]->SetTitleOffset(1.5, "X");
        dPhiFourierSide_la[i]->GetXaxis()->SetTitleSize(0.035);
        dPhiFourierSide_la[i]->GetXaxis()->CenterTitle(true);
        dPhiFourierSide_la[i]->GetXaxis()->SetTitle("#Delta#phi (radians)");

        if(publish)
        {
            ltx3->DrawLatex(0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
            ltx3->DrawLatex(0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}");
            ltx3->DrawLatex(0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
            ltx3->DrawLatex(0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
            ltx3->DrawLatex(0.2, 0.53, "Long range (|#Delta#eta| > 2)");
            ltx2->DrawLatex(0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }
    }

    //Output peak values and errors Kshort
    for(int i=2; i<numFourierParams; i++)
    {
        OutputVnValues(i,"Peak","Kshort",vnValues_ks_peak[i],PtBin_ks,vnPeakName);
        OutputVnErrors(i,"Peak","Kshort",vnErrors_ks_peak[i],PtBin_ks,vnPeakName);
    }
    //Output side values and errors Kshort
    for(int i=2; i<numFourierParams; i++)
    {
        OutputVnValues(i,"Side","Kshort",vnValues_ks_side[i],PtBin_ks,vnSideName);
        OutputVnErrors(i,"Side","Kshort",vnErrors_ks_side[i],PtBin_ks,vnSideName);
    }

    //Output peak values and errors Lambda
    //for(int i=2; i<numFourierParams; i++)
    for(int i=2; i<numFourierParams; i++)
    {
        OutputVnValues(i,"Peak","Lambda",vnValues_la_peak[i],PtBin_la,vnPeakName);
        OutputVnErrors(i,"Peak","Lambda",vnErrors_la_peak[i],PtBin_la,vnPeakName);
    }
    //Output side values and errors Lambda
    for(int i=2; i<numFourierParams; i++)
    {
        OutputVnValues(i,"Side","Lambda",vnValues_la_side[i],PtBin_la,vnSideName);
        OutputVnErrors(i,"Side","Lambda",vnErrors_la_side[i],PtBin_la,vnSideName);
    }

    //Output hadron vn values and errors
    for(int i=2; i<numFourierParams; i++)
    {
        OutputVnH(i,vnValues_h,vnErrors_h,"vnHadron.txt");
    }

    TCanvas* Fourier_ks_peak = new TCanvas("Fourier_ks_peak", "Fourier_ks_peak", 1600,800);
    Fourier_ks_peak->Divide(5,3);
    for(int i=0; i<numPtBins_ks; i++){
        Fourier_ks_peak->cd(i+1);
        gPad->SetTickx();
        gPad->SetTicky();
        dPhiFourierPeak_ks[i]->Draw("E1");
    }
    TCanvas* Fourier_ks_side = new TCanvas("Fourier_ks_side","Fourier_ks_side",1600,800);
    Fourier_ks_side->Divide(5,3);
    for(int i=0; i<numPtBins_ks; i++){
        Fourier_ks_side->cd(i+1);
        gPad->SetTickx();
        gPad->SetTicky();
        dPhiFourierSide_ks[i]->Draw("E1");
    }

    TCanvas* Fourier_la_peak = new TCanvas("Fourier_la_peak", "Fourier_la_peak", 1600,800);
    Fourier_la_peak->Divide(5,3);
    for(int i=0; i<numPtBins_la; i++){
        Fourier_la_peak->cd(i+1);
        gPad->SetTickx();
        gPad->SetTicky();
        dPhiFourierPeak_la[i]->Draw("E1");
    }
    TCanvas* Fourier_la_side = new TCanvas("Fourier_la_side","Fourier_la_side",1600,800);
    Fourier_la_side->Divide(5,3);
    for(int i=0; i<numPtBins_la; i++){
        Fourier_la_side->cd(i+1);
        gPad->SetTickx();
        gPad->SetTicky();
        dPhiFourierSide_la[i]->Draw("E1");
    }

    //Calculate flow
    for(int i=2; i<numFourierParams; i++) vnCalculate(i,"Kshort",vnValues_ks_peak[i],vnErrors_ks_peak[i],vnValues_ks_side[i],vnErrors_ks_side[i],vnValues_h,vnErrors_h,fsig_ks,vnCalculatorName);
    for(int i=2; i<numFourierParams; i++) vnCalculate(i,"Lambda",vnValues_la_peak[i],vnErrors_la_peak[i],vnValues_la_side[i],vnErrors_la_side[i],vnValues_h,vnErrors_h,fsig_la,vnCalculatorName);

    std::map<std::string, std::vector<double> > results_ks = vnCalculateMap(2,"Kshort",vnValues_ks_peak[2],vnErrors_ks_peak[2],vnValues_ks_side[2],vnErrors_ks_side[2],vnValues_h,vnErrors_h,fsig_ks);
    std::map<std::string, std::vector<double> > results_la = vnCalculateMap(2,"Lambda",vnValues_la_peak[2],vnErrors_la_peak[2],vnValues_la_side[2],vnErrors_la_side[2],vnValues_h,vnErrors_h,fsig_la);

    std::vector<double> AvgX_ks = AvgX(f,branchname_ks,branchname_ks_bkg,numPtBins_ks);
    std::vector<double> AvgX_la = AvgX(f,branchname_la,branchname_la_bkg,numPtBins_la);

    std::vector<double> AvgX_ket_ks = AvgX(f,branchname_ket_ks,branchname_ket_ks_bkg,numPtBins_ks);
    std::vector<double> AvgX_ket_la = AvgX(f,branchname_ket_la,branchname_ket_la_bkg,numPtBins_la);

    vnGraph(results_ks,AvgX_ks,AvgX_ket_ks,"Kshort",graphName);
    vnGraph(results_la,AvgX_la,AvgX_ket_la,"Lambda",graphName);

    //Output Publication plots
	if(publish)
	{
    //1D correlation functions
	//Kshort
    if(Peak){
        TCanvas* PubFourier_ks = new TCanvas("PubFourier_ks", "Pub_ks", 800,800);
        PubFourier_ks->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        TH1D* dPhiFourierPeakCopy_ks = (TH1D*)dPhiFourierPeak_ks[4]->Clone();
        dPhiFourierPeakCopy_ks->SetTitle("Peak");
        dPhiFourierPeakCopy_ks->Draw("E1");

        TLatex *ltx3 = new TLatex();
        ltx3->SetTextSize(0.035);
        ltx3->SetNDC(kTRUE);
        ltx3->SetTextFont(42);
        if(publish)
        {
            ltx3->DrawLatex(0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}");
            ltx3->DrawLatex(0.2, 0.74,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
            ltx3->DrawLatex(0.2, 0.67, "0.3 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
            ltx3->DrawLatex(0.2, 0.60, "2.8 < p_{T}^{#Xi} < 3.6 GeV");
            ltx3->DrawLatex(0.2, 0.53, "Long range (|#Delta#eta| > 2)");
            ltx3->SetTextSize(0.045);
            ltx3->DrawLatex(0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }
    }
    else{
        TCanvas* PubFourier_ks = new TCanvas("PubFourier_ks", "Pub_ks", 800,800);
        PubFourier_ks->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        TH1D* dPhiFourierSideCopy_ks = (TH1D*)dPhiFourierSide_ks[4]->Clone();
        dPhiFourierSideCopy_ks->SetTitle("SideBand");
        dPhiFourierSideCopy_ks->Draw("E1");


        TLatex *ltx3 = new TLatex();
        ltx3->SetTextSize(0.035);
        ltx3->SetNDC(kTRUE);
        ltx3->SetTextFont(42);
        if(publish)
        {
            ltx3->DrawLatex(0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}");
            ltx3->DrawLatex(0.2, 0.74,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
            ltx3->DrawLatex(0.2, 0.67, "0.3 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
            ltx3->DrawLatex(0.2, 0.60, "2.8 < p_{T}^{#Xi} < 3.6 GeV");
            ltx3->DrawLatex(0.2, 0.53, "Long range (|#Delta#eta| > 2)");
            ltx3->SetTextSize(0.045);
            ltx3->DrawLatex(0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }
    }


	// Lambda
    if(Peak){
        TCanvas* PubFourier_la = new TCanvas("PubFourier_la", "Pub_la", 800,800);
        PubFourier_la->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        TH1D* dPhiFourierPeakCopy_la = (TH1D*)dPhiFourierPeak_la[4]->Clone();
        dPhiFourierPeakCopy_la->SetTitle("Peak");
        dPhiFourierPeakCopy_la->Draw("E1");

        TLatex *ltx3 = new TLatex();
        ltx3->SetTextSize(0.035);
        ltx3->SetNDC(kTRUE);
        ltx3->SetTextFont(42);
        if(publish)
        {
            ltx3->DrawLatex(0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}");
            ltx3->DrawLatex(0.2, 0.74,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
            ltx3->DrawLatex(0.2, 0.67, "0.3 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
            ltx3->DrawLatex(0.2, 0.60, "2.8 < p_{T}^{#Xi} < 3.6 GeV");
            ltx3->DrawLatex(0.2, 0.53, "Long range (|#Delta#eta| > 2)");
            ltx3->SetTextSize(0.045);
            ltx3->DrawLatex(0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }
    }
    else{
        TCanvas* PubFourier_la = new TCanvas("PubFourier_la", "Pub_la", 800,800);
        PubFourier_la->cd();
        gPad->SetTickx();
        gPad->SetTicky();
        TH1D* dPhiFourierSideCopy_la = (TH1D*)dPhiFourierSide_la[4]->Clone();
        dPhiFourierSideCopy_la->SetTitle("SideBand");
        dPhiFourierSideCopy_la->Draw("E1");


        TLatex *ltx3 = new TLatex();
        ltx3->SetTextSize(0.035);
        ltx3->SetNDC(kTRUE);
        ltx3->SetTextFont(42);
        if(publish)
        {
            ltx3->DrawLatex(0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}");
            ltx3->DrawLatex(0.2, 0.74,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
            ltx3->DrawLatex(0.2, 0.67, "0.3 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
            ltx3->DrawLatex(0.2, 0.60, "2.8 < p_{T}^{#Xi} < 3.6 GeV");
            ltx3->DrawLatex(0.2, 0.53, "Long range (|#Delta#eta| > 2)");
            ltx3->SetTextSize(0.045);
            ltx3->DrawLatex(0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
        }
    }
	}

/*
    //2D Correlation function 1-3 GeV associated
    TLatex *ltx0 = new TLatex();
    ltx0->SetTextSize(0.031);
    ltx0->SetNDC(kTRUE);
    ltx0->SetTextFont(42);

	//Kshort
    TCanvas* TwoDCorrelation_ks = new TCanvas("TwoDCorrelation_ks", "", 1000, 1000);
    TwoDCorrelation_ks->SetLeftMargin(0.2);

    TH2D* Signal_ks = (TH2D*)f->Get("v0CorrelationRapidity/signalkshort_pt2");
    TH2D* Background_ks = (TH2D*)f->Get("v0CorrelationRapidity/backgroundkshort_pt2");

    TGaxis::SetMaxDigits(1);

    TH2D* Correlation_ks = (TH2D*)Signal_ks->Clone();
    Correlation_ks->Divide(Background_ks);
    Correlation_ks->GetXaxis()->SetRangeUser(-4.0,4.0);
    Correlation_ks->GetYaxis()->SetRangeUser(-PI/2.0,4.5);
    Correlation_ks->GetXaxis()->SetTitle("#Delta#eta");
    Correlation_ks->GetXaxis()->SetTitleOffset(1.4);
    Correlation_ks->GetXaxis()->CenterTitle(true);
    Correlation_ks->GetYaxis()->SetTitle("#Delta#phi (radians)");
    Correlation_ks->GetYaxis()->SetTitleOffset(1.4);
    Correlation_ks->GetYaxis()->CenterTitle(true);
    Correlation_ks->GetZaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{d^{2}N^{pair}}{d#Delta#eta d#Delta#phi}");
    Correlation_ks->GetZaxis()->SetTitleOffset(2.3);
    Correlation_ks->GetZaxis()->CenterTitle(true);
    Correlation_ks->GetXaxis()->SetNdivisions(405);
    Correlation_ks->GetYaxis()->SetNdivisions(405);
    Correlation_ks->GetZaxis()->SetNdivisions(4);
    Correlation_ks->SetTitle("");
    Correlation_ks->SetStats(kFALSE);
    Correlation_ks->Scale(10);

    const Int_t NRGBs = 5;
    const Int_t NCont = 20;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    //gStyle->SetNumberContours(90);

    TwoDCorrelation_ks->cd();
    Correlation_ks->Draw("SURF1 FB ");


    ltx0->DrawLatex(0.05, 0.95, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}");
    ltx0->DrawLatex(0.05, 0.88, "185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
    ltx0->DrawLatex(0.05, 0.81, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
    ltx0->DrawLatex(0.05, 0.75, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{trig}}} < 3 GeV");
    ltx0->SetTextSize(0.04);
    ltx0->DrawLatex(0.85, 0.88, "K_{S}^{0}#kern[-0.3]{#lower[0.02]{{}^{#pm}}}- h^{#pm}");


	//Lambda
    TCanvas* TwoDCorrelation_la = new TCanvas("TwoDCorrelation_la", "", 1000, 1000);
    TwoDCorrelation_la->SetLeftMargin(0.2);

    TH2D* Signal_la = (TH2D*)f->Get("v0CorrelationRapidity/signallambda_pt2");
    TH2D* Background_la = (TH2D*)f->Get("v0CorrelationRapidity/backgroundlambda_pt2");

    TGaxis::SetMaxDigits(1);

    TH2D* Correlation_la = (TH2D*)Signal_la->Clone();
    Correlation_la->Divide(Background_la);
    Correlation_la->GetXaxis()->SetRangeUser(-4.0,4.0);
    Correlation_la->GetYaxis()->SetRangeUser(-PI/2.0,4.5);
    Correlation_la->GetXaxis()->SetTitle("#Delta#eta");
    Correlation_la->GetXaxis()->SetTitleOffset(1.4);
    Correlation_la->GetXaxis()->CenterTitle(true);
    Correlation_la->GetYaxis()->SetTitle("#Delta#phi (radians)");
    Correlation_la->GetYaxis()->SetTitleOffset(1.4);
    Correlation_la->GetYaxis()->CenterTitle(true);
    Correlation_la->GetZaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{d^{2}N^{pair}}{d#Delta#eta d#Delta#phi}");
    Correlation_la->GetZaxis()->SetTitleOffset(2.3);
    Correlation_la->GetZaxis()->CenterTitle(true);
    Correlation_la->GetXaxis()->SetNdivisions(405);
    Correlation_la->GetYaxis()->SetNdivisions(405);
    Correlation_la->GetZaxis()->SetNdivisions(4);
    Correlation_la->SetTitle("");
    Correlation_la->SetStats(kFALSE);
    Correlation_la->Scale(10);

    TwoDCorrelation_la->cd();
    Correlation_la->Draw("SURF1 FB ");


    ltx0->DrawLatex(0.05, 0.95, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}");
    ltx0->DrawLatex(0.05, 0.88, "185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
    ltx0->DrawLatex(0.05, 0.81, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
    ltx0->DrawLatex(0.05, 0.75, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{trig}}} < 3 GeV");
    ltx0->SetTextSize(0.04);
    ltx0->DrawLatex(0.85, 0.88, "#Lambda#kern[-0.3]{#lower[0.02]{{}^{#pm}}}- h^{#pm}");
    */
}
