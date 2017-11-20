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
    Double_t xx1 = par[0]/(2*PI);
    //Double_t xx1 = par[0];
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
    if(V0IDname == "Lambda" || V0IDname == "Xi" || V0IDname == "Omega") divFactor = 3;

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
    if(V0IDname == "Lambda" || V0IDname == "Xi" || V0IDname == "Omega") divFactor = 3;

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

std::map<std::string, std::vector<double> > vnCalculateMapGen(int degree, std::string V0IDname, std::vector<double> vnvalues_peak, std::vector<double> vnerrors_peak, std::vector<double> vnvalues_h, std::vector<double> vnerrors_h, std::vector<double> fsig)
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
        vnObs.push_back(vnObs_);
        double sig = vnObs_;

        double vnObsError = vnObs_*TMath::Sqrt(TMath::Power(vnerrors_peak[i]/vnvalues_peak[i],2) + TMath::Power(0.5*vnerrors_h[degree]/vnvalues_h[degree],2));
        double sigError = vnObsError;

        vnRefError.push_back(0.5*vnerrors_h[degree]/TMath::Sqrt(vnvalues_h[degree]));

        sigvalues.push_back(sig);
        sigerrors.push_back(sigError);
        vnObs_errors.push_back(vnObsError);
    }

    int divFactor = 1;
    if(V0IDname == "Kshort") divFactor = 2;
    if(V0IDname == "Lambda" || V0IDname == "Xi" || V0IDname == "Omega") divFactor = 3;

    for(unsigned i=0; i<sigvalues.size(); i++   ) sigvalues_nq.push_back(sigvalues[i]/divFactor);
    for(unsigned i=0; i<sigerrors.size(); i++   ) sigerrors_nq.push_back(sigerrors[i]/divFactor);
    for(unsigned i=0; i<vnObs.size(); i++       ) vnObs_nq.push_back(vnObs[i]/divFactor);
    for(unsigned i=0; i<vnObs_errors.size(); i++) vnObs_errors_nq.push_back(vnObs_errors[i]/divFactor);

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

    if(V0ID == "Lambda" || V0ID == "Xi" || V0ID == "Omega")
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
    if(V0ID == "Xi")
    {
        v2           -> Write("v2xi",TObject::kOverwrite);
        v2_nq        -> Write("v2xi_nq",TObject::kOverwrite);
        v2_ket       -> Write("v2xi_ket",TObject::kOverwrite);
        v2_ket_nq    -> Write("v2xi_ket_nq",TObject::kOverwrite);
        v2obs        -> Write("v2obsxi",TObject::kOverwrite);
        v2obs_nq     -> Write("v2obsxi_nq",TObject::kOverwrite);
        v2obs_ket    -> Write("v2obsxi_ket",TObject::kOverwrite);
        v2obs_ket_nq -> Write("v2obsxi_ket_nq",TObject::kOverwrite);
        v2bkg        -> Write("v2bkgxi",TObject::kOverwrite);
        v2bkg_nq     -> Write("v2bkgxi_nq",TObject::kOverwrite);
        v2bkg_ket    -> Write("v2bkgxi_ket",TObject::kOverwrite);
        v2bkg_ket_nq -> Write("v2bkgxi_ket_nq",TObject::kOverwrite);
    }
    if(V0ID == "Omega")
    {
        v2           -> Write("v2omega",TObject::kOverwrite);
        v2_nq        -> Write("v2omega_nq",TObject::kOverwrite);
        v2_ket       -> Write("v2omega_ket",TObject::kOverwrite);
        v2_ket_nq    -> Write("v2omega_ket_nq",TObject::kOverwrite);
        v2obs        -> Write("v2obsomega",TObject::kOverwrite);
        v2obs_nq     -> Write("v2obsomega_nq",TObject::kOverwrite);
        v2obs_ket    -> Write("v2obsomega_ket",TObject::kOverwrite);
        v2obs_ket_nq -> Write("v2obsomega_ket_nq",TObject::kOverwrite);
        v2bkg        -> Write("v2bkgomega",TObject::kOverwrite);
        v2bkg_nq     -> Write("v2bkgomega_nq",TObject::kOverwrite);
        v2bkg_ket    -> Write("v2bkgomega_ket",TObject::kOverwrite);
        v2bkg_ket_nq -> Write("v2bkgomega_ket_nq",TObject::kOverwrite);
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
    //TFile *f = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/LooseAndTight/V0CorrelationTightMCTotal_08_23_2017.root");

    //TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/PeripheralSubtractionMB.root");
    TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/V0CasCorrelationPbPbTotal_10_30_17.root"); //PbPb
    TFile *f_Xi = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/V0CasCorrelationPbPbTotal_10_30_17.root"); //PbPb
    TFile *f_Om = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/V0CasCorrelationPbPbTotal_10_30_17.root"); //PbPb
    TFile *fhad = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/HadCorr/PbPb_Cent3050_ref.root" ); //PbPb
    //TFile *f_Xi = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationRapidityTotal_08_20_2017.root" );
    //TFile *f = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityCorrectMultB_09_19_17.root"); //For V0s only pPb Directory not called v0CasCorrelationRapidity Dont use
    //
    //Newest Files
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityHM185_250.root"); //For V0s only pPb
    //TFile *f_Xi = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationHM_11_07_17.root" ); //pPb
    //TFile *f_Om = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/OmCorr/OmegaCorrelationRapidityHM185_250.root" ); //pPb
    //TFile *fhad = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/HadCorr/Combine_HM185_corr_ref_PUrej.root" ); //pPb
    //
    //TFile *f = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/OmCorr/OmegaCorrelationRapidityFixedSideband_84percentStats_10_30_17.root" );
    //TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/V0CasCorrelationRapidityMerged.root"); // Merged pPb file
    //TFile *fhad = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationRapidityTotal_08_20_2017.root" ); //pPb
    //
    //TFile *f_Xi = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/PeripheralSubtractionMB.root"); //MB PeriSub BadMult Xi
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/PeripheralSubtractionMB.root"); //MB PeriSub
    //TFile *f_Xi = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MBCorr/XiOmegaMB_0_N_20_Partial_11_8_17.root"); //MB PeriSub
    //TFile *f_Om = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MBCorr/XiOmegaMB_0_N_20_Partial_11_8_17.root"); //MB PeriSub
    //TFile *fhad = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/HadCorr/Combine_MB0_corr_ref.root" ); //MB peripheral subtraction
    //
    //Closure
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0CorrelationRapidityClosureNoEff_09_11_17.root"); //No Eff closure for D0 study
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0CorrelationClosureTotal_08_28_2017.root"); //Normal Closure Reco level for D0 study
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0CorrelationClosureGenTotal_08_28_2017.root"); //Normal Closure Gen level for D0 study
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/MatchV0ClosureBpPb_09_20_17.root"); //Match Closure for D0 study using reco ref
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0ClosureReco_GenRef_10_24_17.root"); //Closure normal reco with genRef
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0ClosureGenNoStrange_10_25_17.root"); //Closure Gen reco with genRef NoStrange
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0ClosureGen_RecoRef_10_23_17.root"); //Closure Gen w/ reco ref to compare with Match w/ reco ref for D0 study
    //TFile *fhad = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0CorrelationClosureHadronRecoFix_09_10_17.root"); //For vn of hadron in closure (USE FOR GEN ONLY / RECO HAS NO EFF CORRECTION)
    //
    //Newer Closure Files 11_8_17
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0ClosureReco_GenRef_NoV0Eff_11_8_17.root"); //Closure normal reco No V0 Eff with genRef
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0ClosureReco_GenRef_11_8_17.root"); //Closure normal reco with genRef
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0ClosureGen_RecoRef_11_8_17.root"); //Closure Gen w/ reco ref
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0ClosureReco_RecoRef_11_8_17.root"); //Normal Closure Reco level for D0 study
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0ClosureReco_RecoRef_NoV0Eff_Total_11_8_17.root"); //Normal Closure Reco level No V0 Eff for D0 study
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0ClosureGen_GenRef_AllParticles_Total_11_8_17.root"); //Normal Closure Gen level for D0 study
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0ClosureReco_RecoRef_top20percent_Total_11_14_17.root"); //Normal Closure Reco level top 20% for D0 study
    //TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0ClosureGen_GenRef_top20percent_Total_11_14_17.root"); //Normal Closure Gen level top 20% for D0 study
    //TFile *fhad = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0ClosureRefReco_wEff_11_8_17.root"); //For vn of hadron in closure (FOR RECO ONLY)
    //TFile *fhad = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/V0CorrelationClosureHadronRecoFix_09_10_17.root"); //For vn of hadron in closure (USE FOR GEN ONLY / RECO HAS NO EFF CORRECTION)

	//Txt files
	ofstream vnPeak;
	ofstream vnSide;
    ofstream vnH;
    ofstream vnCalculator;
    bool doKs = true;
    bool doLa = true;
    bool doXi = true;
    bool doOm = true;
    bool doGen = false;
    //std::string fn = "v0CasCorrelationRapidityPeriSub";
    //std::string fn = "v0CorrelationRapidity";
    //std::string fn = "xiCorrelationRapidity";
    //std::string fn = "v0CasCorrelationRapidity";
    std::string fn = "v0CasCorrelationRapidityPbPb";
    //std::string fn = "v0CorrelationRapidityMatchMC";
    //std::string fn = "v0CorrelationRapidityMC";//Normal Reco Closure
    //std::string fn = "v0CorrelationMCRapidity";//Normal Gen Closure
    //std::string fn = "v0CorrelationRapidityMCGenRef";//Normal Reco w/ genRef Closure
    //std::string fn = "v0CorrelationMCRapidityRecoRef";//Normal gen w/ RecoRef Closure
    //std::string fn = "v0CorrelationRapidityMatchMC";//Match Reco w/ RecoRef Closure
    //std::string vnPeakName = "vnPeakPeriSub.txt";
    //std::string vnSideName = "vnPeakPeriSub.txt";
    //std::string vnPeakName = "vnPeak.txt";
    //std::string vnSideName = "vnPeak.txt";
    //std::string vnPeakName = "vnPeak_V0ClosureGenRef_10_24_17.txt";
    std::string vnPeakName = "vnPeak_V0ClosureGenNoStrange_10_25_17.txt";
    //std::string vnPeakName = "vnPeak_V0ClosureGen_RecoRef_10_23_17.txt";
    std::string vnSideName = vnPeakName;
    //std::string vnHName = "vnHadronMB.txt";
    //std::string vnHName = "vnHadron.txt";
    std::string vnHName = "vnHadron.txt";
    //std::string vnCalculatorName = "vnSignalPeriSub.txt";
    //std::string vnCalculatorName = "vnSignal.txt";
    std::string vnCalculatorName = "vnSignal_V0ClosureGenNoStrange_10_25_17.txt"; //Normal Reco w/ GenRef
    //std::string vnCalculatorName = "vnSignal_V0ClosureGen_RecoRef_10_23_17.txt"; //Normal Reco w/ GenRef
    std::string branchname_ks         = fn + "/Ptkshort_pt";
    std::string branchname_ks_bkg     = fn + "/Ptkshort_bkg_pt";
    std::string branchname_la         = fn + "/Ptlambda_pt";
    std::string branchname_la_bkg     = fn + "/Ptlambda_bkg_pt";
    std::string branchname_xi         = fn + "/Pt_xi_pt";
    std::string branchname_xi_bkg     = fn + "/Pt_xi_bkg_pt";
    std::string branchname_ket_ks     = fn + "/KETkshort_pt";
    std::string branchname_ket_ks_bkg = fn + "/KETkshort_bkg_pt";
    std::string branchname_ket_la     = fn + "/KETlambda_pt";
    std::string branchname_ket_la_bkg = fn + "/KETlambda_bkg_pt";
    std::string branchname_ket_xi     = fn + "/KET_xi_pt";
    std::string branchname_ket_xi_bkg = fn + "/KET_xi_bkg_pt";
    std::string branchname_om         = fn + "/Pt_om_pt";
    std::string branchname_om_bkg     = fn + "/Pt_om_bkg_pt";
    std::string branchname_ket_om     = fn + "/KET_om_pt";
    std::string branchname_ket_om_bkg = fn + "/KET_om_bkg_pt";
    //std::string graphName = "FitRootFiles/v2valuesRapidityClosureReco_RecoRef_Total_11_8_17.root";
    std::string graphName = "FitRootFiles/PbPb/v2valuesRapidity_Cent_30_50perc_EtaGap2_11_9_17.root"; //PbPb
    //std::string graphName = "FitRootFiles/HM/v2valuesRapidityHM_185_250_EtaGap2_11_9_17.root"; //HM
    //std::string graphName = "FitRootFiles/MB/v2valuesRapidityMB_0_20_EtaGap2_11_9_17.root"; //MB
    //std::string graphName = "v2valuesRapidityPeripheralSub.root";
    //std::string graphName = "FitRootFiles/Closure/v2valuesRapidityClosure_V0CorrelationClosureGenTotal_Top20percent_11_14_17.root";//Normal Gen Closure
    //std::string graphName = "v2valuesRapidityClosure_V0ClosureReco_GenRef_10_24_17.root";//Normal Reco w/ GenRef
    //std::string graphName = "v2valuesRapidityClosure_V0ClosureGenNoStrange_10_25_17.root";//Normal Reco w/ GenRef
    //std::string graphName = "v2valuesRapidityClosure_MatchV0ClosureBpPb_09_20_17.root";//Match Reco w/ GenRef
    //std::string graphName = "v2valuesRapidityClosure_V0ClosureGen_RecoRef_10_23_17.root";//Gen w/ RecoRef
    int binlow = 10;
    int binhigh = 24;
	vnPeak.open(vnPeakName);
	vnSide.open(vnSideName);
    vnH.open(vnHName);
    vnCalculator.open(vnCalculatorName);

    TVirtualFitter::SetMaxIterations(300000);
    TH1::SetDefaultSumw2();
    std::vector<double> PtBin_ks = {0.2, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5};//, 10.0, 15.0, 20.0}; 
    std::vector<double> PtBin_la = {0.8, 1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.0, 8.5};//, 10.0, 15.0, 20.0};
    //std::vector<double> PtBin_xi = {1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10.0};//, 20.0};
    std::vector<double> PtBin_xi = {1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2};//, 20.0}; // PbPb
    //std::vector<double> PtBin_om = {1.0, 1.5, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10.0};//, 20.0}; //pPb
    //std::vector<double> PtBin_om = {1.5, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10.0};//, 20.0}; //pPb removed 1.0-1.5
    //std::vector<double> PtBin_om = {1.0, 1.8, 2.3, 3.0, 10.0};//pPb MB
    //std::vector<double> PtBin_om = {1.5, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10.0};//, 20.0}; // For Peripheral Subtraction removing 1-1.5
    std::vector<double> PtBin_om = {1.5, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0};//, 20.0}; // PbPb
    int numPtBins_ks = PtBin_ks.size() - 1;
    int numPtBins_la = PtBin_la.size() - 1;
    int numPtBins_xi = PtBin_xi.size() - 1;
    int numPtBins_om = PtBin_om.size() - 1;
    TH1D* dPhiFourierPeak_ks[numPtBins_ks];
    TH1D* dPhiFourierSide_ks[numPtBins_ks];
    TH1D* dPhiFourierPeak_la[numPtBins_la];
    TH1D* dPhiFourierSide_la[numPtBins_la];
    TH1D* dPhiFourierPeak_xi[numPtBins_xi];
    TH1D* dPhiFourierSide_xi[numPtBins_xi];
    TH1D* dPhiFourierPeak_om[numPtBins_om];
    TH1D* dPhiFourierSide_om[numPtBins_om];

    TH1D* dPhiPeak_ks[numPtBins_ks];
    TH1D* dPhiSide_ks[numPtBins_ks];
    TH1D* dPhiPeak_la[numPtBins_la];
    TH1D* dPhiSide_la[numPtBins_la];
    TH1D* dPhiPeak_xi[numPtBins_xi];
    TH1D* dPhiSide_xi[numPtBins_xi];
    TH1D* dPhiPeak_om[numPtBins_om];
    TH1D* dPhiSide_om[numPtBins_om];

    TF1* FourierFit_ks[numPtBins_ks];
    TF1* FourierFit_la[numPtBins_la];
    TF1* FourierFit_xi[numPtBins_xi];
    TF1* FourierFit_om[numPtBins_om];

    std::map<int,double> FourierScale_ks;
    std::map<int,std::vector<double> > vnValues_ks_peak;
    std::map<int,std::vector<double> > vnValues_ks_side;

    std::map<int,std::vector<double> > vnErrors_ks_peak;
    std::map<int,std::vector<double> > vnErrors_ks_side;

    std::map<int,double> FourierScale_la;
    std::map<int,std::vector<double> > vnValues_la_peak;
    std::map<int,std::vector<double> > vnValues_la_side;

    std::map<int,std::vector<double> > vnErrors_la_peak;
    std::map<int,std::vector<double> > vnErrors_la_side;

    std::vector<double> v2values_xi_peak;
    std::vector<double> v2values_xi_side;

    std::vector<double> v2errors_xi_peak;
    std::vector<double> v2errors_xi_side;

    std::vector<double> v2values_om_peak;
    std::vector<double> v2values_om_side;

    std::vector<double> v2errors_om_peak;
    std::vector<double> v2errors_om_side;

    /**
     * Initialized first two elements because I want to keep the index consistently starting at 2
     **/
    std::vector<double> vnValues_h = {-999,-999};
    std::vector<double> vnErrors_h = {-999,-999};

    //Fsig for vn calculations
    //std::vector<double> fsig_ks = {0.999666 ,0.999977 ,0.999972 ,0.999988 ,0.999998 ,0.999999 ,0.999999 ,0.999992 ,0.999994 ,0.999955, 0.999975};
    //std::vector<double> fsig_la = {0.988877 ,0.9967 ,0.99754 ,0.998939 ,0.999954 ,0.999951 ,0.999992 ,0.999842};
    //std::vector<double> fsig_ks = {0.999476 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,0.999999 ,0.999999 ,0.999988 ,0.999997 ,0.992327};// ,0.99836};// ,0.890106 ,0.481433}; //Rapidity pPb partial stats
    //std::vector<double> fsig_la = {0.99882 ,0.999987 ,1 ,0.999524 ,0.999632 ,0.999855 ,0.999698 ,0.998783 ,0.999771 ,0.997088};// ,0.92718};// ,0.990913 ,0.421011}; //Rapidity pPb partial stats
    //std::vector<double> fsig_xi = {0.959427 ,0.976239 ,0.979161 ,0.980678 ,0.980661 ,0.981534 ,0.981502 ,0.979289 ,0.979192};// ,0.976012}; //Rapidity pPb Partial Stats
    //std::vector<double> fsig_om = {0.519579 ,0.664389 ,0.757589 ,0.820397 ,0.864778 ,0.903763 ,0.933867 ,0.928403 ,0.949529};// ,0.976012}; Rapidity pPb w/ 1.0-1.5
    //std::vector<double> fsig_om = {0.721825, 0.807658, 0.861617, 0.908762, 0.937818, 0.970699, 0.966515, 0.96482};// ,0.976012}; Rapidity pPb removed 1.0-1.5 Partial Stats
    //std::vector<double> fsig_ks = {0.991743 ,0.992798 ,0.993457 ,0.993601 ,0.992867 ,0.990667 ,0.98788 ,0.984159 ,0.979121 ,0.973607 ,0.968452 ,0.962748 ,0.957971};// ,0.99836};// ,0.890106 ,0.481433}; //Rapidity pPb full stats
    //std::vector<double> fsig_la = {0.899098,0.969946 ,0.981187 ,0.983357 ,0.98448 ,0.984754 ,0.983104 ,0.97844 ,0.970099 ,0.960906};// ,0.92718};// ,0.990913 ,0.421011}; //Rapidity pPb full stats
    //std::vector<double> fsig_xi = {0.946441 ,0.969308 ,0.972483 ,0.973308 ,0.973584 ,0.972577 ,0.972438 ,0.968581 ,0.966521};// ,0.976012}; //Rapidity pPb Full Stats
    //std::vector<double> fsig_om = {0.722342, 0.806663, 0.874069, 0.922878, 0.950436, 0.965793, 0.970866, 0.96041};// ,0.976012}; Rapidity pPb removed 1.0-1.5 Full Stats
    //std::vector<double> fsig_om = {0.817472, 0.911707, 0.943929, 0.962094};// ,0.976012}; Rapidity pPb removed 1.0-1.5
    //std::vector<double> fsig_om = {0.664389 ,0.757589 ,0.820397 ,0.864778 ,0.903763 ,0.933867 ,0.928403 ,0.949529};// ,0.976012}; //For Peripheral Sub
    //std::vector<double> fsig_ks = {0.88371 ,0.921642 ,0.933393 ,0.944669 ,0.957631 ,0.958982 ,0.949047 ,0.947948 ,0.91845 ,0.915206 ,0.913798 ,0.90336 ,0.889312}; //PbPb partial stats
    //std::vector<double> fsig_la = {0.685977,0.851926 ,0.928584 ,0.944353 ,0.952072 ,0.951948 ,0.946417 ,0.931891 ,0.88348 ,0.872633}; //PbPb Partial Stats 
    //std::vector<double> fsig_xi = {0.869648 ,0.90817 ,0.939477 ,0.944172 ,0.940623 ,0.930704 ,0.936042 ,0.938471};//,0.976012}; //PbPb Partial Stats
    //std::vector<double> fsig_om = {0.647487 ,0.774924 ,0.872547 ,0.911986 ,0.937495 ,0.95909};// ,0.976012}; //For PbPb partial Stat
    std::vector<double> fsig_ks = {0.892205,0.92336 ,0.935079 ,0.946728 ,0.961131 ,0.960816 ,0.951312 ,0.935589 ,0.922555,0.91897 ,0.918127 ,0.908964,0.897616}; //PbPb full stats
    std::vector<double> fsig_la = {0.685943,0.855173 ,0.928891 ,0.944278 ,0.952467 ,0.951831 ,0.946294 ,0.931161 ,0.89215 ,0.839763}; //PbPb Full Stats
    std::vector<double> fsig_xi = {0.877281 ,0.911338,0.942308 ,0.947788 ,0.943301 ,0.936111 ,0.935527 ,0.933834};//,0.976012}; //PbPb Full Stats
    std::vector<double> fsig_om = {0.659704 ,0.797872 ,0.887173 ,0.925122 ,0.942805 ,0.95212};// ,0.976012}; //For PbPb Full Stat
    //std::vector<double> fsig_ks = {0.991217 ,1 ,1 ,1 ,1 ,1 ,1 ,0.999988 ,0.99998 ,0.998263 ,0.942699 ,0.999827 ,0.99997};// ,0.99836};// ,0.890106 ,0.481433}; //MC
    //std::vector<double> fsig_la = {0.971604 ,0.985982 ,0.991208 ,0.989859 ,0.997168 ,0.996337 ,0.982485 ,0.953336 ,0.781606 ,0.596769};// ,0.92718};// ,0.990913 ,0.421011}; //MC


    if((PtBin_ks.size()-1 != fsig_ks.size()) || (PtBin_la.size()-1 != fsig_la.size()) || (PtBin_xi.size()-1 != fsig_xi.size()) || PtBin_om.size()-1 != fsig_om.size())
    {
        cout << "something is wrong with number of pt bins or number of fsig" << endl;
        return;
    }

    TLatex* ltx2 = new TLatex();
    ltx2->SetTextSize(0.045);
    ltx2->SetNDC(kTRUE);

    //Vn hadron
    //TH2D *hBackgroundHad  = (TH2D*) fhad->Get("pPbCorr/background"); //pPb & MB
    //TH2D *hSignalHad      = (TH2D*) fhad->Get("pPbCorr/signal"); //pPb & MB
    TH2D *hBackgroundHad  = (TH2D*) fhad->Get("PbPbCorr/background"); //PbPb
    TH2D *hSignalHad      = (TH2D*) fhad->Get("PbPbCorr/signal"); //PbPb
    //TH2D *hBackgroundHad  = (TH2D*) fhad->Get("HadronCorrelation/BackgroundHadReco"); //Closure
    //TH2D *hSignalHad      = (TH2D*) fhad->Get("HadronCorrelation/SignalHadReco"); //Closure
    //TH2D *hBackgroundHad  = (TH2D*) fhad->Get("HadronCorrelation/BackgroundHad"); //Closure Gen Ref
    //TH2D *hSignalHad      = (TH2D*) fhad->Get("HadronCorrelation/SignalHad"); // Closure Gen Ref
    TH1D *dPhiHad = new TH1D("dPhiHad", "h^{#pm}- h^{#pm} ", 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
    TH1D* hbHadPhiTot = hBackgroundHad->ProjectionY("PhiBkgHadTot", 1, binlow);
    TH1D* hbHadPhiOth = hBackgroundHad->ProjectionY("PhiBkgHadOth", binhigh, -1);
    TH1D* hsHadPhiTot = hSignalHad->ProjectionY("PhiSigHadTot", 1, binlow);
    TH1D* hsHadPhiOth = hSignalHad->ProjectionY("PhiSigHadOth", binhigh, -1);

    hbHadPhiTot->Add(hbHadPhiOth);
    hsHadPhiTot->Add(hsHadPhiOth);

    dPhiHad->Divide(hsHadPhiTot, hbHadPhiTot);
    TH1D* dPhiHadFourier = (TH1D*)dPhiHad->Clone();

    TF1 *FourierFitHad = new TF1("FourierFitHad", FourierHad, -1.5, 5, numFourierParams);
    FourierFitHad->SetNpx(250);
    FourierFitHad->SetParNames("Scale", "V_{1}", "V_{2}", "V_{3}");
    dPhiHadFourier->Fit("FourierFitHad");
    dPhiHadFourier->SetStats(kFALSE);
    for(int j=2; j<numFourierParams; j++)
    {
        vnValues_h.push_back(FourierFitHad->GetParameter(j));
        vnErrors_h.push_back(FourierFitHad->GetParError(j));
    }
    //vnValues_h.push_back(6.8e-2); //For reasonable Closure results V2 would have to be around this 
    //vnErrors_h.push_back(2e-6);

    if(doKs)
    {
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
            //Pull 2D Histograms
            TH2D *hbackgroundPeak = (TH2D*) f_V0->Get(Form((fn + "/backgroundkshort_pt%d").c_str(),i));
            TH2D *hsignalPeak     = (TH2D*) f_V0->Get(Form((fn + "/signalkshort_pt%d").c_str(),i));

            //Project Phi

            // For projecting both shoulders
            TH1D* hbPhiTotPeak = hbackgroundPeak->ProjectionY("PhiBkgTotPeak", 1, binlow);
            TH1D* hbPhiOthPeak = hbackgroundPeak->ProjectionY("PhiBkgOthPeak", binhigh, -1);
            TH1D* hsPhiTotPeak = hsignalPeak->ProjectionY("PhiSigTotPeak", 1, binlow);
            TH1D* hsPhiOthPeak = hsignalPeak->ProjectionY("PhiSigOthPeak", binhigh, -1);

            hbPhiTotPeak->Add(hbPhiOthPeak);
            hsPhiTotPeak->Add(hsPhiOthPeak);

            //Divide
            dPhiPeak_ks[i]->Divide(hsPhiTotPeak, hbPhiTotPeak);

            //Clone histograms for display without fit functions
            dPhiFourierPeak_ks[i] = (TH1D*)dPhiPeak_ks[i]->Clone();

            FourierFit_ks[i] = new TF1(Form("FourierFit_ks%d",i), FourierHad, -1.5, 5, numFourierParams);
            FourierFit_ks[i]->SetNpx(250);
            FourierFit_ks[i]->SetParNames("Scale", "V_{1}", "V_{2}", "V_{3}");

            dPhiFourierPeak_ks[i]->Fit(Form("FourierFit_ks%d",i));
            dPhiFourierPeak_ks[i]->SetStats(kFALSE);
            FourierScale_ks[i] = FourierFit_ks[i]->GetParameter(0);
            for(int j=2; j<numFourierParams; j++)
            {
                vnValues_ks_peak[j].push_back(FourierFit_ks[i]->GetParameter(j));
                vnErrors_ks_peak[j].push_back(FourierFit_ks[i]->GetParError(j));
            }


            //TCanvas *c5_ks = new TCanvas("c5_ks", "c5_ks", 800,800);
            //c5_ks->cd();
            //gPad->SetTickx();
            //gPad->SetTicky();

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

            //maxBinContent = dPhiHadFourier->GetBinContent(dPhiHadFourier->GetMaximumBin());
            //minBinContent = dPhiHadFourier->GetBinContent(dPhiHadFourier->GetMinimumBin());
            //minRange = minBinContent - minFact*minBinContent;
            //maxRange = minRange + maxFact*(maxBinContent - minBinContent);

            //TCanvas *c3_ks = new TCanvas("c3_ks", "ZYAM Hadronic", 800,800);
            //c3_ks->cd();
            //gPad->SetTickx();
            //gPad->SetTicky();
            //dPhiHad->SetMarkerStyle(34);
            //dPhiHad->SetMarkerSize(1.5);
            //dPhiHad->SetTitleOffset(2, "Y");
            //dPhiHad->SetTitle("");
            //dPhiHad->GetYaxis()->SetRangeUser(mini , maxi);
            //dPhiHad->GetYaxis()->SetTitleSize(0.03);
            //dPhiHad->GetYaxis()->CenterTitle(true);
            //dPhiHad->GetYaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi}");
            //dPhiHad->SetTitleOffset(1.5, "X");
            //dPhiHad->GetXaxis()->SetTitleSize(0.035);
            //dPhiHad->GetXaxis()->CenterTitle(true);
            //dPhiHad->GetXaxis()->SetTitle("#Delta#phi (radians)");
            //dPhiHad->Fit("pol2","","", 0.4,2);
            //dPhiHad->SetStats(!publish);

            //if(publish)
            //{
                //ltx3->DrawLatex(0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
                //ltx3->DrawLatex(0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}");
                //ltx3->DrawLatex(0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
                //ltx3->DrawLatex(0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
                //ltx3->DrawLatex(0.2, 0.53, "Long range (|#Delta#eta| > 2)");
                //ltx2->DrawLatex(0.7, 0.74, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
            //}

            TF1 *dPhiFitPeak = dPhiPeak_ks[i]->GetFunction("pol2");
            //TF1 *dPhiHadFit = dPhiHad->GetFunction("pol2");

            double dPhiFitMinPeak = dPhiFitPeak->GetMinimum();
            //double dPhiHadFitMin = dPhiHadFit->GetMinimum();

            for(int j = 1; j < 32; j++)
            {
                dPhiFourierPeak_ks[i]->AddBinContent(j, -dPhiFitMinPeak);
            }

            //for(int j = 1; j < 32; j++)
            //{
                //dPhiHadFourier->AddBinContent(j, -dPhiHadFitMin);
            //}

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

            //dPhiHadFourier->SetMarkerStyle(34);
            //dPhiHadFourier->SetMarkerSize(1.5);
            //dPhiHadFourier->Draw("E1");
            //dPhiHadFourier->SetStats(kFALSE);
            //dPhiHadFourier->Fit("FourierFitHad");
            //dPhiHadFourier->SetTitle("");
            //dPhiHadFourier->SetTitleOffset(2, "Y");
            //dPhiHadFourier->GetYaxis()->CenterTitle(true);
            //dPhiHadFourier->GetYaxis()->SetTitleSize(0.03);
            //dPhiHadFourier->GetYaxis()->SetTitle("#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}");
            //dPhiHadFourier->GetYaxis()->SetRangeUser(-0.0004 , 0.008);
            //dPhiHadFourier->SetTitleOffset(1.5, "X");
            //dPhiHadFourier->GetXaxis()->SetTitleSize(0.035);
            //dPhiHadFourier->GetXaxis()->CenterTitle(true);
            //dPhiHadFourier->GetXaxis()->SetTitle("#Delta#phi (radians)");

            //if(publish)
            //{
                //ltx3->DrawLatex(0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV");
                //ltx3->DrawLatex(0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}");
                //ltx3->DrawLatex(0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220");
                //ltx3->DrawLatex(0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV");
                //ltx3->DrawLatex(0.2, 0.53, "Long range (|#Delta#eta| > 2)");
                //ltx2->DrawLatex(0.7, 0.74, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}");
            //}

            if(!doGen)
            {
                //================================================================================
                //Sideband Calculations
                //================================================================================
                dPhiSide_ks[i] = new TH1D(Form("dPhiSide_ks%d",i), "K_{S}^{0} - h^{#pm} ", 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
                TH2D *hbackgroundSide = (TH2D*) f_V0->Get(Form((fn + "/backgroundkshort_bkg_pt%d").c_str(),i));
                TH2D *hsignalSide     = (TH2D*) f_V0->Get(Form((fn + "/signalkshort_bkg_pt%d").c_str(),i));

                //Project Phi
                TH1D* hbPhiTotSide = hbackgroundSide->ProjectionY("PhiBkgTot", 0, binlow);
                TH1D* hbPhiOthSide = hbackgroundSide->ProjectionY("PhiBkgOthPeak", binhigh, -1);
                TH1D* hsPhiTotSide = hsignalSide->ProjectionY("PhiSigTot", 0, binlow);
                TH1D* hsPhiOthSide = hsignalSide->ProjectionY("PhiSigOthPeak", binhigh, -1);

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
        }
        cout << "1" << endl;
        TCanvas* Fourier_ks_peak = new TCanvas("Fourier_ks_peak", "Fourier_ks_peak", 1600,800);
        Fourier_ks_peak->Divide(5,3);
        for(int i=0; i<numPtBins_ks; i++){
            Fourier_ks_peak->cd(i+1);
            gPad->SetTickx();
            gPad->SetTicky();
            dPhiFourierPeak_ks[i]->Draw("E1");
            cout << i << endl;
        }
        if(!doGen)
        {
            TCanvas* Fourier_ks_side = new TCanvas("Fourier_ks_side","Fourier_ks_side",1600,800);
            Fourier_ks_side->Divide(5,3);
            for(int i=0; i<numPtBins_ks; i++){
                Fourier_ks_side->cd(i+1);
                gPad->SetTickx();
                gPad->SetTicky();
                dPhiFourierSide_ks[i]->Draw("E1");
            }
        }

        std::map<std::string, std::vector<double> > results_ks;
        std::map<std::string, std::vector<double> > gen_results_ks;
        if(doGen)
        {
            std::vector<double> AvgX_ks;
            std::vector<double> AvgX_ket_ks;
            std::string branch_pt = branchname_ks += "%d";
            std::string branch_ket = branchname_ket_ks += "%d";
            for(int i=0; i<numPtBins_ks; i++)
            {
                AvgX_ks.push_back(((TH1D*)f_V0->Get(Form(branch_pt.c_str(),i)))->GetMean(1));
                AvgX_ket_ks.push_back(((TH1D*)f_V0->Get(Form((branch_ket).c_str(),i)))->GetMean(1));
            }
            gen_results_ks = vnCalculateMapGen(2,"Kshort",vnValues_ks_peak[2],vnErrors_ks_peak[2],vnValues_h,vnErrors_h,fsig_ks);
            vnGraph(gen_results_ks,AvgX_ks,AvgX_ket_ks,"Kshort",graphName);
        }
        else
        {
            std::vector<double> AvgX_ks = AvgX(f_V0,branchname_ks,branchname_ks_bkg,numPtBins_ks);
            std::vector<double> AvgX_ket_ks = AvgX(f_V0,branchname_ket_ks,branchname_ket_ks_bkg,numPtBins_ks);
            for(int i=2; i<numFourierParams; i++) vnCalculate(i,"Kshort",vnValues_ks_peak[i],vnErrors_ks_peak[i],vnValues_ks_side[i],vnErrors_ks_side[i],vnValues_h,vnErrors_h,fsig_ks,vnCalculatorName);
            results_ks = vnCalculateMap(2,"Kshort",vnValues_ks_peak[2],vnErrors_ks_peak[2],vnValues_ks_side[2],vnErrors_ks_side[2],vnValues_h,vnErrors_h,fsig_ks);
            vnGraph(results_ks,AvgX_ks,AvgX_ket_ks,"Kshort",graphName);
        }

    }

    if(doLa)
    {
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
            //Pull 2D Histograms
            TH2D *hbackgroundPeak = (TH2D*) f_V0->Get(Form((fn + "/backgroundlambda_pt%d").c_str(),i));
            TH2D *hsignalPeak     = (TH2D*) f_V0->Get(Form((fn + "/signallambda_pt%d").c_str(),i));

            //Project Phi

            TH1D* hbPhiTotPeak = hbackgroundPeak->ProjectionY("PhiBkgTotPeak", 1, binlow);
            TH1D* hbPhiOthPeak = hbackgroundPeak->ProjectionY("PhiBkgOthPeak", binhigh, -1);
            TH1D* hsPhiTotPeak = hsignalPeak->ProjectionY("PhiSigTotPeak", 1, binlow);
            TH1D* hsPhiOthPeak = hsignalPeak->ProjectionY("PhiSigOthPeak", binhigh, -1);

            hbPhiTotPeak->Add(hbPhiOthPeak);
            hsPhiTotPeak->Add(hsPhiOthPeak);

            //Divide
            dPhiPeak_la[i]->Divide(hsPhiTotPeak, hbPhiTotPeak);

            //Clone histograms for display without fit functions
            dPhiFourierPeak_la[i] = (TH1D*)dPhiPeak_la[i]->Clone();

            FourierFit_la[i] = new TF1(Form("FourierFit_la%d",i), FourierHad, -1.5, 5, numFourierParams);
            FourierFit_la[i]->SetNpx(250);
            FourierFit_la[i]->SetParNames("Scale", "V_{1}", "V_{2}", "V_{3}");

            TCanvas *c4_la_peak = new TCanvas("c4_la_peak", "Fourier peak trg_la", 800,800);
            c4_la_peak->cd();
            gPad->SetTickx();
            gPad->SetTicky();
            dPhiFourierPeak_la[i]->Fit(Form("FourierFit_la%d",i));
            dPhiFourierPeak_la[i]->SetStats(kFALSE);
            FourierScale_la[i] = FourierFit_la[i]->GetParameter(0);
            for(int j=2; j<numFourierParams; j++)
            {
                vnValues_la_peak[j].push_back(FourierFit_la[i]->GetParameter(j));
                vnErrors_la_peak[j].push_back(FourierFit_la[i]->GetParError(j));
            }

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

            TF1 *dPhiFitPeak = dPhiPeak_la[i]->GetFunction("pol2");

            double dPhiFitMinPeak = dPhiFitPeak->GetMinimum();

            for(int j = 1; j < 32; j++)
            {
                dPhiFourierPeak_la[i]->AddBinContent(j, -dPhiFitMinPeak);
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


            if(!doGen)
            {
                //================================================================================
                //Sideband Calculations
                //================================================================================
                dPhiSide_la[i] = new TH1D(Form("dPhiSide_la%d",i), "K_{S}^{0} - h^{#pm} ", 31, -(0.5 - 1.0/32)*PI, (1.5 - 1.0/32)*PI);
                TH2D *hbackgroundSide = (TH2D*) f_V0->Get(Form((fn + "/backgroundlambda_bkg_pt%d").c_str(),i));
                TH2D *hsignalSide     = (TH2D*) f_V0->Get(Form((fn + "/signallambda_bkg_pt%d").c_str(),i));

                //Project Phi
                TH1D* hbPhiTotSide = hbackgroundSide->ProjectionY("PhiBkgTot", 0, binlow);
                TH1D* hbPhiOthSide = hbackgroundSide->ProjectionY("PhiBkgOthPeak", binhigh, -1);
                TH1D* hsPhiTotSide = hsignalSide->ProjectionY("PhiSigTot", 0, binlow);
                TH1D* hsPhiOthSide = hsignalSide->ProjectionY("PhiSigOthPeak", binhigh, -1);

                hbPhiTotSide->Add(hbPhiOthSide);
                hsPhiTotSide->Add(hsPhiOthSide);

                //Divide
                dPhiSide_la[i]->Divide(hsPhiTotSide, hbPhiTotSide);

                //Clone histograms for display without fit functions
                dPhiFourierSide_la[i] = (TH1D*)dPhiSide_la[i]->Clone();

                FourierFit_la[i] = new TF1(Form("FourierFit_la%d",i) , FourierHad, -1.5, 5, numFourierParams);
                FourierFit_la[i]->SetNpx(250);
                FourierFit_la[i]->SetParNames("Scale", "V_{1}", "V_{2}", "V_{3}");

                dPhiFourierSide_la[i]->Fit(Form("FourierFit_la%d",i));
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
                dPhiSide_la[i]->Fit("pol2","","", 0.4,2.4);
                dPhiSide_la[i]->SetStats(!publish);

                TF1 *dPhiFitSide = dPhiSide_la[i]->GetFunction("pol2");

                double dPhiFitMinSide = dPhiFitSide->GetMinimum();

                for(int j = 1; j < 32; j++)
                {
                    dPhiFourierSide_la[i]->AddBinContent(j, -dPhiFitMinSide);
                }

                dPhiFourierSide_la[i]->Fit(Form("FourierFit_la%d",i));
            }
        }

        TCanvas* Fourier_la_peak = new TCanvas("Fourier_la_peak", "Fourier_la_peak", 1600,800);
        Fourier_la_peak->Divide(5,3);
        for(int i=0; i<numPtBins_la; i++){
            Fourier_la_peak->cd(i+1);
            gPad->SetTickx();
            gPad->SetTicky();
            dPhiFourierPeak_la[i]->Draw("E1");
        }
        if(!doGen)
        {
            TCanvas* Fourier_la_side = new TCanvas("Fourier_la_side","Fourier_la_side",1600,800);
            Fourier_la_side->Divide(5,3);
            for(int i=0; i<numPtBins_la; i++){
                Fourier_la_side->cd(i+1);
                gPad->SetTickx();
                gPad->SetTicky();
                dPhiFourierSide_la[i]->Draw("E1");
            }
        }

        std::map<std::string, std::vector<double> > results_la;
        std::map<std::string, std::vector<double> > gen_results_la;
        if(doGen)
        {
            std::vector<double> AvgX_la;
            std::vector<double> AvgX_ket_la;
            std::string branch_pt = branchname_la += "%d";
            std::string branch_ket = branchname_ket_la += "%d";
            for(int i=0; i<numPtBins_la; i++)
            {
                AvgX_la.push_back(((TH1D*)f_V0->Get(Form(branch_pt.c_str(),i)))->GetMean(1));
                AvgX_ket_la.push_back(((TH1D*)f_V0->Get(Form((branch_ket).c_str(),i)))->GetMean(1));
            }
            gen_results_la = vnCalculateMapGen(2,"Lambda",vnValues_la_peak[2],vnErrors_la_peak[2],vnValues_h,vnErrors_h,fsig_la);
            vnGraph(gen_results_la,AvgX_la,AvgX_ket_la,"Lambda",graphName);
        }
        else
        {
            std::vector<double> AvgX_la = AvgX(f_V0,branchname_la,branchname_la_bkg,numPtBins_la);
            std::vector<double> AvgX_ket_la = AvgX(f_V0,branchname_ket_la,branchname_ket_la_bkg,numPtBins_la);
            for(int i=2; i<numFourierParams; i++) vnCalculate(i,"Lambda",vnValues_la_peak[i],vnErrors_la_peak[i],vnValues_la_side[i],vnErrors_la_side[i],vnValues_h,vnErrors_h,fsig_la,vnCalculatorName);
            results_la = vnCalculateMap(2,"Lambda",vnValues_la_peak[2],vnErrors_la_peak[2],vnValues_la_side[2],vnErrors_la_side[2],vnValues_h,vnErrors_h,fsig_la);
            vnGraph(results_la,AvgX_la,AvgX_ket_la,"Lambda",graphName);
        }
    }

    if(doXi)
    {
        //XI
        cout << "================================================================================" << endl;
        cout << "================================================================================" << endl;
        cout << "================================================================================" << endl;
        cout << "XI XI XI"                                                                         << endl;
        cout << "================================================================================" << endl;
        cout << "================================================================================" << endl;
        cout << "================================================================================" << endl;
        for( int i=0; i<numPtBins_xi; i++ )
        {
            dPhiPeak_xi[i] = new TH1D( Form( "dPhiPeak_xi%d",i ), "#Xi - h^{#pm} ", 31, -( 0.5 -
                        1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI  );
            //Pull 2D Histograms
            TH2D *hbackgroundPeak = (TH2D*) f_Xi->Get( Form((fn + "/BackgroundXiPeak_pt%d").c_str(),i ) );
            TH2D *hsignalPeak     = (TH2D*) f_Xi->Get( Form((fn + "/SignalXiPeak_pt%d").c_str(),i ) );
            //TH2D *hbackgroundPeak = (TH2D*) f_Xi->Get( Form((fn + "/BackgroundPeak_pt%d").c_str(),i ) );
            //TH2D *hsignalPeak     = (TH2D*) f_Xi->Get( Form((fn + "/SignalPeak_pt%d").c_str(),i ) );

            TH1::SetDefaultSumw2(  );
            //Project Phi

            // For projecting both shoulders
            TH1D* hbPhiTotPeak = hbackgroundPeak->ProjectionY( "PhiBkgTotPeak", 1, binlow );
            TH1D* hbPhiOthPeak = hbackgroundPeak->ProjectionY( "PhiBkgOthPeak", binhigh, -1 );
            TH1D* hsPhiTotPeak = hsignalPeak->ProjectionY( "PhiSigTotPeak", 1, binlow );
            TH1D* hsPhiOthPeak = hsignalPeak->ProjectionY( "PhiSigOthPeak", binhigh, -1 );

            hbPhiTotPeak->Add( hbPhiOthPeak );
            hsPhiTotPeak->Add( hsPhiOthPeak );

            //Divide
            dPhiPeak_xi[i]->Divide( hsPhiTotPeak, hbPhiTotPeak );

            //Clone histograms for display without fit functions
            dPhiFourierPeak_xi[i] = ( TH1D* )dPhiPeak_xi[i]->Clone(  );

            FourierFit_xi[i] = new TF1( Form( "FourierFit_xi%d",i ), FourierHad, -1.5, 5, numFourierParams );
            FourierFit_xi[i]->SetNpx( 250 );
            FourierFit_xi[i]->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *FourierPeak = new TCanvas( "FourierPeak", "Fourier Peak", 800,800 );
            FourierPeak->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            //dPhiFourierPeak_xi->Fit( "FourierFit_xi","","",0,PI );
            dPhiFourierPeak_xi[i]->Fit( Form( "FourierFit_xi%d",i ) );
            dPhiFourierPeak_xi[i]->SetStats( kFALSE );
            v2values_xi_peak.push_back( FourierFit_xi[i]->GetParameter( 2 ) );
            v2errors_xi_peak.push_back( FourierFit_xi[i]->GetParError( 2 ) );
            cout << "---------------------------------" << endl;
            cout << "Peak V2 for xi-h is " << FourierFit_xi[i]->GetParameter( 2 ) << endl;
            cout << "---------------------------------" << endl;

            double maxBinContent = dPhiFourierPeak_xi[i]->GetBinContent( dPhiFourierPeak_xi[i]->GetMaximumBin(  ) );
            double minBinContent = dPhiFourierPeak_xi[i]->GetBinContent( dPhiFourierPeak_xi[i]->GetMinimumBin(  ) );
            double minRange = minBinContent - 0.005*minBinContent;
            double maxRange = minRange + 2*( maxBinContent - minBinContent );

            //ZYAM FITS
            TCanvas *ZYAMFitPeak = new TCanvas( "ZYAMFitPeak", "ZYAM Fit Peak", 800,800 );
            ZYAMFitPeak->cd(  );

            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiPeak_xi[i]->SetMarkerStyle( 21 );
            dPhiPeak_xi[i]->SetMarkerColor( 4 );
            dPhiPeak_xi[i]->SetTitleOffset( 2, "Y" );
            dPhiPeak_xi[i]->SetTitle( "Peak" );
            dPhiPeak_xi[i]->GetYaxis(  )->SetRangeUser( minRange , maxRange );
            dPhiPeak_xi[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiPeak_xi[i]->GetYaxis(  )->CenterTitle( true );
            dPhiPeak_xi[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} " );
            dPhiPeak_xi[i]->SetTitleOffset( 1.5, "X" );
            dPhiPeak_xi[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiPeak_xi[i]->GetXaxis(  )->CenterTitle( true );
            dPhiPeak_xi[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
            //dPhiPeak_xi[i]->Draw( "E1" );
            dPhiPeak_xi[i]->Fit( "pol2","","", 0.4,2.4 );
            dPhiPeak_xi[i]->SetStats( !publish );

            TLatex *ltx3 = new TLatex(  );
            ltx3->SetTextSize( 0.035 );
            ltx3->SetNDC( kTRUE );
            ltx3->SetTextFont( 42 );

            if( publish )
            {
                os << "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = " << fixed << std::setprecision( 2 ) << SNN << " TeV";
                ltx3->DrawLatex( 0.2, 0.82, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                os << "L_{#lower[-0.25]{int}} = " << Lint << " nb^{-1}";
                ltx3->DrawLatex( 0.2, 0.74, os.str( ).c_str(  ) );
                os.str( std::string(  ) );
                os << Nmin << "  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< " << Nmax;
                ltx3->DrawLatex( 0.2, 0.67, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                os << pTassMin << " < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < " << pTassMax << " GeV";
                ltx3->DrawLatex( 0.2, 0.60, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                os << "Long range (|#Delta#eta| > " <<  longRange << ")";
                ltx3->DrawLatex( 0.2, 0.53, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );

                //ltx3->DrawLatex( 0.2, 0.82, "CMS pPb  #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                //ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                //ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                //ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                //ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                //ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }


            TF1 *dPhiFitPeak = dPhiPeak_xi[i]->GetFunction( "pol2" );

            double dPhiFitMinPeak = dPhiFitPeak->GetMinimum(  );

            for( int j = 1; j < 32; j++ )
            {
                dPhiFourierPeak_xi[i]->AddBinContent( j, -dPhiFitMinPeak );
            }

            FourierPeak->cd(  );
            dPhiFourierPeak_xi[i]->SetMarkerStyle( 21 );
            dPhiFourierPeak_xi[i]->SetMarkerColor( 4 );
            //dPhiFourierPeak_xi[i]->Draw( "E1" );
            dPhiFourierPeak_xi[i]->SetStats( kFALSE );
            dPhiFourierPeak_xi[i]->Fit( Form( "FourierFit_xi%d",i ) );
            os << "Peak " << PtBin_xi[i] << "_Pt_" << PtBin_xi[i+1];
            dPhiFourierPeak_xi[i]->SetTitle( os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            dPhiFourierPeak_xi[i]->SetTitleOffset( 2, "Y" );
            dPhiFourierPeak_xi[i]->GetYaxis(  )->CenterTitle( true );
            dPhiFourierPeak_xi[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiFourierPeak_xi[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
            dPhiFourierPeak_xi[i]->GetYaxis(  )->SetRangeUser( -0.0004, 0.008 );
            dPhiFourierPeak_xi[i]->SetTitleOffset( 1.5, "X" );
            dPhiFourierPeak_xi[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiFourierPeak_xi[i]->GetXaxis(  )->CenterTitle( true );
            dPhiFourierPeak_xi[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );


            if( publish )
            {
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }

            //side
            dPhiSide_xi[i] = new TH1D( Form( "dPhiSide%d",i ), "#Xi - h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI  );
            //TH2D *hbackgroundSide = (TH2D*) f_Xi->Get( Form( (fn + "/BackgroundSide_pt%d").c_str(),i ) );
            //TH2D *hsignalSide     = (TH2D*) f_Xi->Get( Form( (fn + "/SignalSide_pt%d").c_str(),i ) );
            TH2D *hbackgroundSide = (TH2D*) f_Xi->Get( Form( (fn + "/BackgroundXiSide_pt%d").c_str(),i ) );
            TH2D *hsignalSide     = (TH2D*) f_Xi->Get( Form( (fn + "/SignalXiSide_pt%d").c_str(),i ) );

            TH1::SetDefaultSumw2(  );

            //Project Phi
            TH1D* hbPhiTotSide = hbackgroundSide->ProjectionY( "PhiBkgTot", 0, binlow );
            TH1D* hbPhiOthSide = hbackgroundSide->ProjectionY( "PhiBkgOthPeak", binhigh, -1 );
            TH1D* hsPhiTotSide = hsignalSide->ProjectionY( "PhiSigTot", 0, binlow );
            TH1D* hsPhiOthSide = hsignalSide->ProjectionY( "PhiSigOthPeak", binhigh, -1 );

            hbPhiTotSide->Add( hbPhiOthSide );
            hsPhiTotSide->Add( hsPhiOthSide );

            //Divide
            dPhiSide_xi[i]->Divide( hsPhiTotSide, hbPhiTotSide );

            //Clone histograms for display without fit functions
            dPhiFourierSide_xi[i] = ( TH1D* )dPhiSide_xi[i]->Clone(  );

            FourierFit_xi[i] = new TF1( Form( "FourierFit_xi%d",i ) , FourierHad, -1.5, 5, numFourierParams );
            FourierFit_xi[i]->SetNpx( 250 );
            FourierFit_xi[i]->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *FourierSide = new TCanvas( "FourierSide", "Fourier Side", 800,800 );
            FourierSide->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiFourierSide_xi[i]->Fit( Form( "FourierFit_xi%d",i ) );
            dPhiFourierSide_xi[i]->SetStats( kFALSE );
            v2values_xi_side.push_back( FourierFit_xi[i]->GetParameter( 2 ) );
            v2errors_xi_side.push_back( FourierFit_xi[i]->GetParError( 2 ) );

            maxBinContent = dPhiFourierSide_xi[i]->GetBinContent( dPhiFourierSide_xi[i]->GetMaximumBin(  ) );
            minBinContent = dPhiFourierSide_xi[i]->GetBinContent( dPhiFourierSide_xi[i]->GetMinimumBin(  ) );
            minRange = minBinContent - 0.005*minBinContent;
            maxRange = minRange + 2*( maxBinContent - minBinContent );


            //ZYAM FITS
            TCanvas *ZYAMFitSide = new TCanvas( "ZYAMFitSide", "ZYAM Fit Side", 800,800 );
            ZYAMFitSide->cd(  );

            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiSide_xi[i]->SetMarkerStyle( 21 );
            dPhiSide_xi[i]->SetMarkerColor( 4 );
            dPhiSide_xi[i]->SetTitleOffset( 2, "Y" );
            dPhiSide_xi[i]->SetTitle( "Sideband" );
            dPhiSide_xi[i]->GetYaxis(  )->SetRangeUser( minRange , maxRange );
            dPhiSide_xi[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiSide_xi[i]->GetYaxis(  )->CenterTitle( true );
            dPhiSide_xi[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} " );
            dPhiSide_xi[i]->SetTitleOffset( 1.5, "X" );
            dPhiSide_xi[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiSide_xi[i]->GetXaxis(  )->CenterTitle( true );
            dPhiSide_xi[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
            dPhiSide_xi[i]->Fit( "pol2","","", 0.4,2.4 );
            dPhiSide_xi[i]->SetStats( !publish );

            if( publish )
            {
                os << "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = " << fixed << std::setprecision( 2 ) << SNN << " TeV";
                ltx3->DrawLatex( 0.2, 0.82, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                os << "L_{#lower[-0.25]{int}} = " << Lint << " nb^{-1}";
                ltx3->DrawLatex( 0.2, 0.74, os.str( ).c_str(  ) );
                os.str( std::string(  ) );
                os << Nmin << "  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< " << Nmax;
                ltx3->DrawLatex( 0.2, 0.67, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                os << pTassMin << " < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < " << pTassMax << " GeV";
                ltx3->DrawLatex( 0.2, 0.60, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                os << "Long range (|#Delta#eta| > " <<  longRange << ")";
                ltx3->DrawLatex( 0.2, 0.53, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );

                //ltx3->DrawLatex( 0.2, 0.82, "CMS pPb  #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                //ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                //ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                //ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                //ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                //ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }

            TF1 *dPhiFitSide = dPhiSide_xi[i]->GetFunction( "pol2" );

            double dPhiFitMinSide = dPhiFitSide->GetMinimum(  );

            for( int j = 1; j < 32; j++ )
            {
                dPhiFourierSide_xi[i]->AddBinContent( j, -dPhiFitMinSide );
            }


            FourierSide->cd(  );
            dPhiFourierSide_xi[i]->SetMarkerStyle( 21 );
            dPhiFourierSide_xi[i]->SetMarkerColor( 4 );
            //dPhiFourierSide_xi[i]->Draw( "E1" );
            dPhiFourierSide_xi[i]->SetStats( kFALSE );
            dPhiFourierSide_xi[i]->Fit( Form( "FourierFit_xi%d",i ) );
            os << "SideBand " << PtBin_xi[i] << "_Pt_" << PtBin_xi[i+1];
            dPhiFourierSide_xi[i]->SetTitle( os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            dPhiFourierSide_xi[i]->SetTitleOffset( 2, "Y" );
            dPhiFourierSide_xi[i]->GetYaxis(  )->CenterTitle( true );
            dPhiFourierSide_xi[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiFourierSide_xi[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
            dPhiFourierSide_xi[i]->GetYaxis(  )->SetRangeUser( -0.0004, 0.008 );
            dPhiFourierSide_xi[i]->SetTitleOffset( 1.5, "X" );
            dPhiFourierSide_xi[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiFourierSide_xi[i]->GetXaxis(  )->CenterTitle( true );
            dPhiFourierSide_xi[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );

            if( publish )
            {
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }
        }
        vnCalculate(2,"Xi",v2values_xi_peak,v2errors_xi_peak,v2values_xi_side,v2errors_xi_side,vnValues_h,vnErrors_h,fsig_xi,vnCalculatorName);
        std::map<std::string, std::vector<double> > results_xi = vnCalculateMap(2,"Xi",v2values_xi_peak,v2errors_xi_peak,v2values_xi_side,v2errors_xi_side,vnValues_h,vnErrors_h,fsig_xi);
        std::vector<double> AvgX_xi = AvgX(f_Xi,branchname_xi,branchname_xi_bkg,numPtBins_xi);
        std::vector<double> AvgX_ket_xi = AvgX(f_Xi,branchname_ket_xi,branchname_ket_xi_bkg,numPtBins_xi);
        vnGraph(results_xi,AvgX_xi,AvgX_ket_xi,"Xi",graphName);
        OutputVnValues(2,"Peak","Xi",v2values_xi_peak,PtBin_xi,vnPeakName);
        TCanvas* Fourier_xi_peak = new TCanvas("Fourier_xi_peak", "Fourier_xi_peak", 1600,800);
        Fourier_xi_peak->Divide(5,3);
        for(int i=0; i<numPtBins_xi; i++){
            Fourier_xi_peak->cd(i+1);
            gPad->SetTickx();
            gPad->SetTicky();
            dPhiFourierPeak_xi[i]->Draw("E1");
        }
        TCanvas* Fourier_xi_side = new TCanvas("Fourier_xi_side","Fourier_xi_side",1600,800);
        Fourier_xi_side->Divide(5,3);
        for(int i=0; i<numPtBins_xi; i++){
            Fourier_xi_side->cd(i+1);
            gPad->SetTickx();
            gPad->SetTicky();
            dPhiFourierSide_xi[i]->Draw("E1");
        }
    }

    if(doOm)
    {
        //OMEGA
        cout << "================================================================================" << endl;
        cout << "================================================================================" << endl;
        cout << "================================================================================" << endl;
        cout << "OMEGA OMEGA OMEGA"                                                                << endl;
        cout << "================================================================================" << endl;
        cout << "================================================================================" << endl;
        cout << "================================================================================" << endl;
        for( int i=0; i<numPtBins_om; i++ )
        {

            dPhiPeak_om[i] = new TH1D( Form( "dPhiPeak_om%d",i ), "#Xi - h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI  );
            TH1D *dPhiHad = new TH1D( "dPhiHad", "h^{#pm}- h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
            //Pull 2D Histograms
            TH2D *hbackgroundPeak = (TH2D*) f_Om->Get( Form((fn + "/BackgroundOmPeak_pt%d").c_str(),i ) );
            TH2D *hsignalPeak     = (TH2D*) f_Om->Get( Form((fn + "/SignalOmPeak_pt%d").c_str(),i ) );

            TH1::SetDefaultSumw2(  );
            //Project Phi

            // For projecting both shoulders
            TH1D* hbPhiTotPeak = hbackgroundPeak->ProjectionY( "PhiBkgTotPeak", 1, binlow );
            TH1D* hbPhiOthPeak = hbackgroundPeak->ProjectionY( "PhiBkgOthPeak", binhigh, -1 );
            TH1D* hsPhiTotPeak = hsignalPeak->ProjectionY( "PhiSigTotPeak", 1, binlow );
            TH1D* hsPhiOthPeak = hsignalPeak->ProjectionY( "PhiSigOthPeak", binhigh, -1 );

            hbPhiTotPeak->Add( hbPhiOthPeak );
            hsPhiTotPeak->Add( hsPhiOthPeak );

            //Divide
            dPhiPeak_om[i]->Divide( hsPhiTotPeak, hbPhiTotPeak );

            //Clone histograms for display without fit functions
            dPhiFourierPeak_om[i] = ( TH1D* )dPhiPeak_om[i]->Clone(  );

            FourierFit_om[i] = new TF1( Form( "FourierFit_om%d",i ), FourierHad, -1.5, 5, numFourierParams );
            FourierFit_om[i]->SetNpx( 250 );
            FourierFit_om[i]->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *FourierPeak = new TCanvas( "FourierPeak", "Fourier Peak", 800,800 );
            FourierPeak->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            //dPhiFourierPeak_om->Fit( "FourierFit_om","","",0,PI );
            dPhiFourierPeak_om[i]->Fit( Form( "FourierFit_om%d",i ) );
            dPhiFourierPeak_om[i]->SetStats( kFALSE );
            v2values_om_peak.push_back( FourierFit_om[i]->GetParameter( 2 ) );
            v2errors_om_peak.push_back( FourierFit_om[i]->GetParError( 2 ) );

            double maxBinContent = dPhiFourierPeak_om[i]->GetBinContent( dPhiFourierPeak_om[i]->GetMaximumBin(  ) );
            double minBinContent = dPhiFourierPeak_om[i]->GetBinContent( dPhiFourierPeak_om[i]->GetMinimumBin(  ) );
            double minRange = minBinContent - 0.005*minBinContent;
            double maxRange = minRange + 2*( maxBinContent - minBinContent );

            //ZYAM FITS
            TCanvas *ZYAMFitPeak = new TCanvas( "ZYAMFitPeak", "ZYAM Fit Peak", 800,800 );
            ZYAMFitPeak->cd(  );

            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiPeak_om[i]->SetMarkerStyle( 21 );
            dPhiPeak_om[i]->SetMarkerColor( 4 );
            dPhiPeak_om[i]->SetTitleOffset( 2, "Y" );
            dPhiPeak_om[i]->SetTitle( "Peak" );
            dPhiPeak_om[i]->GetYaxis(  )->SetRangeUser( minRange , maxRange );
            dPhiPeak_om[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiPeak_om[i]->GetYaxis(  )->CenterTitle( true );
            dPhiPeak_om[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} " );
            dPhiPeak_om[i]->SetTitleOffset( 1.5, "X" );
            dPhiPeak_om[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiPeak_om[i]->GetXaxis(  )->CenterTitle( true );
            dPhiPeak_om[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
            //dPhiPeak_om[i]->Draw( "E1" );
            dPhiPeak_om[i]->Fit( "pol2","","", 0.4,2.4 );
            dPhiPeak_om[i]->SetStats( !publish );

            TLatex *ltx3 = new TLatex(  );
            ltx3->SetTextSize( 0.035 );
            ltx3->SetNDC( kTRUE );
            ltx3->SetTextFont( 42 );

            if( publish )
            {
                os << "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = " << fixed << std::setprecision( 2 ) << SNN << " TeV";
                ltx3->DrawLatex( 0.2, 0.82, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                os << "L_{#lower[-0.25]{int}} = " << Lint << " nb^{-1}";
                ltx3->DrawLatex( 0.2, 0.74, os.str( ).c_str(  ) );
                os.str( std::string(  ) );
                os << Nmin << "  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< " << Nmax;
                ltx3->DrawLatex( 0.2, 0.67, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                os << pTassMin << " < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < " << pTassMax << " GeV";
                ltx3->DrawLatex( 0.2, 0.60, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                os << "Long range (|#Delta#eta| > " <<  longRange << ")";
                ltx3->DrawLatex( 0.2, 0.53, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );

                //ltx3->DrawLatex( 0.2, 0.82, "CMS pPb  #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                //ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                //ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                //ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                //ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                //ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }


            TF1 *dPhiFitPeak = dPhiPeak_om[i]->GetFunction( "pol2" );

            double dPhiFitMinPeak = dPhiFitPeak->GetMinimum(  );

            for( int j = 1; j < 32; j++ )
            {
                dPhiFourierPeak_om[i]->AddBinContent( j, -dPhiFitMinPeak );
            }

            FourierPeak->cd(  );
            dPhiFourierPeak_om[i]->SetMarkerStyle( 21 );
            dPhiFourierPeak_om[i]->SetMarkerColor( 4 );
            //dPhiFourierPeak_om[i]->Draw( "E1" );
            dPhiFourierPeak_om[i]->SetStats( kFALSE );
            dPhiFourierPeak_om[i]->Fit( Form( "FourierFit_om%d",i ) );
            os << "Peak " << PtBin_om[i] << "_Pt_" << PtBin_om[i+1];
            dPhiFourierPeak_om[i]->SetTitle( os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            dPhiFourierPeak_om[i]->SetTitleOffset( 2, "Y" );
            dPhiFourierPeak_om[i]->GetYaxis(  )->CenterTitle( true );
            dPhiFourierPeak_om[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiFourierPeak_om[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
            dPhiFourierPeak_om[i]->GetYaxis(  )->SetRangeUser( -0.0004, 0.008 );
            dPhiFourierPeak_om[i]->SetTitleOffset( 1.5, "X" );
            dPhiFourierPeak_om[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiFourierPeak_om[i]->GetXaxis(  )->CenterTitle( true );
            dPhiFourierPeak_om[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );


            if( publish )
            {
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }

            //side
            dPhiSide_om[i] = new TH1D( Form( "dPhiSide_om%d",i ), "#Xi - h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI  );
            TH2D *hbackgroundSide = (TH2D*) f_Om->Get( Form( (fn + "/BackgroundOmSide_pt%d").c_str(),i ) );
            TH2D *hsignalSide     = (TH2D*) f_Om->Get( Form( (fn + "/SignalOmSide_pt%d").c_str(),i ) );

            TH1::SetDefaultSumw2(  );

            //Project Phi
            TH1D* hbPhiTotSide = hbackgroundSide->ProjectionY( "PhiBkgTot", 0, binlow );
            TH1D* hbPhiOthSide = hbackgroundSide->ProjectionY( "PhiBkgOthPeak", binhigh, -1 );
            TH1D* hsPhiTotSide = hsignalSide->ProjectionY( "PhiSigTot", 0, binlow );
            TH1D* hsPhiOthSide = hsignalSide->ProjectionY( "PhiSigOthPeak", binhigh, -1 );

            hbPhiTotSide->Add( hbPhiOthSide );
            hsPhiTotSide->Add( hsPhiOthSide );

            hbHadPhiTot->Add( hbHadPhiOth );
            hsHadPhiTot->Add( hsHadPhiOth );

            //Divide
            dPhiSide_om[i]->Divide( hsPhiTotSide, hbPhiTotSide );

            //Clone histograms for display without fit functions
            dPhiFourierSide_om[i] = ( TH1D* )dPhiSide_om[i]->Clone(  );

            FourierFit_om[i] = new TF1( Form( "FourierFit_om%d",i ) , FourierHad, -1.5, 5, numFourierParams);
            FourierFit_om[i]->SetNpx( 250 );
            FourierFit_om[i]->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *FourierSide = new TCanvas( "FourierSide", "Fourier Side", 800,800 );
            FourierSide->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiFourierSide_om[i]->Fit( Form( "FourierFit_om%d",i ) );
            dPhiFourierSide_om[i]->SetStats( kFALSE );
            v2values_om_side.push_back( FourierFit_om[i]->GetParameter( 2 ) );
            v2errors_om_side.push_back( FourierFit_om[i]->GetParError( 2 ) );

            maxBinContent = dPhiFourierSide_om[i]->GetBinContent( dPhiFourierSide_om[i]->GetMaximumBin(  ) );
            minBinContent = dPhiFourierSide_om[i]->GetBinContent( dPhiFourierSide_om[i]->GetMinimumBin(  ) );
            minRange = minBinContent - 0.005*minBinContent;
            maxRange = minRange + 2*( maxBinContent - minBinContent );


            //ZYAM FITS
            TCanvas *ZYAMFitSide = new TCanvas( "ZYAMFitSide", "ZYAM Fit Side", 800,800 );
            ZYAMFitSide->cd(  );

            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiSide_om[i]->SetMarkerStyle( 21 );
            dPhiSide_om[i]->SetMarkerColor( 4 );
            dPhiSide_om[i]->SetTitleOffset( 2, "Y" );
            dPhiSide_om[i]->SetTitle( "Sideband" );
            dPhiSide_om[i]->GetYaxis(  )->SetRangeUser( minRange , maxRange );
            dPhiSide_om[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiSide_om[i]->GetYaxis(  )->CenterTitle( true );
            dPhiSide_om[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} " );
            dPhiSide_om[i]->SetTitleOffset( 1.5, "X" );
            dPhiSide_om[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiSide_om[i]->GetXaxis(  )->CenterTitle( true );
            dPhiSide_om[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
            //dPhiSide_om[i]->Draw( "E1" );
            dPhiSide_om[i]->Fit( "pol2","","", 0.4,2.4 );
            dPhiSide_om[i]->SetStats( !publish );

            if( publish )
            {
                os << "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = " << fixed << std::setprecision( 2 ) << SNN << " TeV";
                ltx3->DrawLatex( 0.2, 0.82, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                os << "L_{#lower[-0.25]{int}} = " << Lint << " nb^{-1}";
                ltx3->DrawLatex( 0.2, 0.74, os.str( ).c_str(  ) );
                os.str( std::string(  ) );
                os << Nmin << "  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< " << Nmax;
                ltx3->DrawLatex( 0.2, 0.67, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                os << pTassMin << " < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < " << pTassMax << " GeV";
                ltx3->DrawLatex( 0.2, 0.60, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                os << "Long range (|#Delta#eta| > " <<  longRange << ")";
                ltx3->DrawLatex( 0.2, 0.53, os.str(  ).c_str(  ) );
                os.str( std::string(  ) );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );

                //ltx3->DrawLatex( 0.2, 0.82, "CMS pPb  #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                //ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                //ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                //ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                //ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                //ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }

            TF1 *dPhiFitSide = dPhiSide_om[i]->GetFunction( "pol2" );

            double dPhiFitMinSide = dPhiFitSide->GetMinimum(  );

            for( int j = 1; j < 32; j++ )
            {
                dPhiFourierSide_om[i]->AddBinContent( j, -dPhiFitMinSide );
            }


            FourierSide->cd(  );
            dPhiFourierSide_om[i]->SetMarkerStyle( 21 );
            dPhiFourierSide_om[i]->SetMarkerColor( 4 );
            //dPhiFourierSide_om[i]->Draw( "E1" );
            dPhiFourierSide_om[i]->SetStats( kFALSE );
            dPhiFourierSide_om[i]->Fit( Form( "FourierFit_om%d",i ) );
            os << "SideBand " << PtBin_om[i] << "_Pt_" << PtBin_om[i+1];
            dPhiFourierSide_om[i]->SetTitle( os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            dPhiFourierSide_om[i]->SetTitleOffset( 2, "Y" );
            dPhiFourierSide_om[i]->GetYaxis(  )->CenterTitle( true );
            dPhiFourierSide_om[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiFourierSide_om[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
            dPhiFourierSide_om[i]->GetYaxis(  )->SetRangeUser( -0.0004, 0.008 );
            dPhiFourierSide_om[i]->SetTitleOffset( 1.5, "X" );
            dPhiFourierSide_om[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiFourierSide_om[i]->GetXaxis(  )->CenterTitle( true );
            dPhiFourierSide_om[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );

            if( publish )
            {
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }
            cout << "iteration: " << i << endl;
        }
        vnCalculate(2,"Omega",v2values_om_peak,v2errors_om_peak,v2values_om_side,v2errors_om_side,vnValues_h,vnErrors_h,fsig_om,vnCalculatorName);
        std::map<std::string, std::vector<double> > results_om = vnCalculateMap(2,"Omega",v2values_om_peak,v2errors_om_peak,v2values_om_side,v2errors_om_side,vnValues_h,vnErrors_h,fsig_om);
        std::vector<double> AvgX_om = AvgX(f_Om,branchname_om,branchname_om_bkg,numPtBins_om);
        std::vector<double> AvgX_ket_om = AvgX(f_Om,branchname_ket_om,branchname_ket_om_bkg,numPtBins_om);
        vnGraph(results_om,AvgX_om,AvgX_ket_om,"Omega",graphName);
        TCanvas* Fourier_om_peak = new TCanvas("Fourier_om_peak", "Fourier_om_peak", 1600,800);
        Fourier_om_peak->Divide(5,3);
        for(int i=0; i<numPtBins_om; i++){
            Fourier_om_peak->cd(i+1);
            gPad->SetTickx();
            gPad->SetTicky();
            dPhiFourierPeak_om[i]->Draw("E1");
        }
        TCanvas* Fourier_om_side = new TCanvas("Fourier_om_side","Fourier_om_side",1600,800);
        Fourier_om_side->Divide(5,3);
        for(int i=0; i<numPtBins_om; i++){
            Fourier_om_side->cd(i+1);
            gPad->SetTickx();
            gPad->SetTicky();
            dPhiFourierSide_om[i]->Draw("E1");
        }
    }
}
