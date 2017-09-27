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
#include "TMathText.h"
#include "TPad.h"
#include <TString.h>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>

#define PI 3.1416

Double_t FourierHad( Double_t *x, Double_t *par )
{
    //Double_t xx1 = par[0]/(2*PI);
    Double_t xx1 = par[0];
    Double_t xx2 = 1 + 2*(par[1]*TMath::Cos( x[0] ) + par[2]*TMath::Cos( 2*x[0] )  + par[3]*TMath::Cos( 3*x[0] ) );
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

    myfile << output.str();
    output.str(std::string());
    for(unsigned i=0; i<sigvalues.size(); i++) myfile << sigvalues[i] << "\n";
    
    myfile << "v2/nq values\n";
    for(unsigned i=0; i<sigvalues.size(); i++) myfile << sigvalues[i]/3 << "\n";

    output << V0IDname << " signal v" << degree << " errors\n";
    myfile << output.str();
    output.str(std::string());
    for(unsigned i=0; i<sigerrors.size(); i++) myfile << sigerrors[i] << "\n";

    output << V0IDname << " signal v" << degree << "/nq errors\n";
    myfile << output.str();
    output.str(std::string());
    for(unsigned i=0; i<sigerrors.size(); i++) myfile << sigerrors[i]/3 << "\n";

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
    if(V0IDname == "Xi") divFactor = 3;

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

    if(V0ID == "Xi")
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

    out->Close();
}

void Xiv2Fit(  )
{
    //TLatex
    std::ostringstream os; // stringstream for making dynamic TLatex labels
    double SNN = 8.16;
    int Lint = 35;
    int Nmin = 185;
    int Nmax = 220;
    int pTassMin = 1;
    int pTassMax = 3;
    int longRange = 2;

    //For Enabling TLatex labels
    //Bool_t publish = kTRUE;
    Bool_t publish = kFALSE;

    gStyle->SetOptFit( 1111 );
    gStyle->SetErrorX( 0 ); //removes horizontal error bars

    //HIST APPEARANCE
    gStyle->SetPadLeftMargin( 0.15 );
    gStyle->SetPadBottomMargin( 0.15 );
    TGaxis::SetMaxDigits( 3 ); // Forces exponents after n number of digits

    //Font and Size
    gStyle->SetTextSize( 20 );
    gStyle->SetTextFont( 42 ); //2=times-bold-r-normal, 2=precision for TLatex to work

    std::string Xiv2PeakName = "Xiv2PeakLoose.txt";
    std::string Xiv2SideName = "Xiv2SideLoose.txt";
    std::string Xiv2CalculatorName = "Xiv2SignalLoose.txt";
    std::string branchname_xi = "xiCorrelationRapidity/Pt_xi_pt";
    std::string branchname_xi_bkg = "xiCorrelationRapidity/Pt_xi_bkg_pt";
    std::string branchname_ket_xi = "xiCorrelationRapidity/KET_xi_pt";
    std::string branchname_ket_xi_bkg = "xiCorrelationRapidity/KET_xi_bkg_pt";
    std::string graphName = "v2valuesRapidity_etaGap_0p75.root";
    int binlow = 14;
    int binhigh = 19;

    std::ofstream Xiv2Peak;
    std::ofstream Xiv2Side;
    std::ofstream Xiv2Calculator;
    Xiv2Peak.open(Xiv2PeakName.c_str());
    Xiv2Side.open(Xiv2SideName.c_str());
    Xiv2Calculator.open(Xiv2CalculatorName.c_str());

    //TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationPD1-6reverseJL10-15_08_15_2017.root" );
    //TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationpPbPD1-6_08_15_2017.root" );
    //TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationRapidityLoose_08_30_2017.root" );
    //TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationRapidityTight_08_30_2017.root" );
    TFile *f = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationRapidityTotal_08_20_2017.root" );
    TFile *fhad = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationRapidityTotal_08_20_2017.root" );

    TVirtualFitter::SetMaxIterations( 300000 );
    TH1::SetDefaultSumw2(  );
    //std::vector<double> PtBin = {1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 8.5, 10.0, 20.0};
    //std::vector<double> fsig_xi = {0.954019 ,0.973881 ,0.976705 ,0.97829 ,0.978074 ,0.978057 ,0.978603 ,0.974644 ,0.975063 ,0.97418 ,0.976012};
    //std::vector<double> PtBin = {1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 10.0, 20.0};
    //std::vector<double> fsig_xi = {0.954019 ,0.973881 ,0.976705 ,0.97829 ,0.978074 ,0.978057 ,0.978603 ,0.974693 ,0.976012};
    std::vector<double> PtBin = {1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10.0};//, 20.0};
    std::vector<double> fsig_xi = {0.954019 ,0.973881 ,0.976705 ,0.97829 ,0.978074 ,0.978057 ,0.978603 ,0.974644 ,0.974909};// ,0.976012};

    int numPtBins = PtBin.size()-1;
    TH1D* dPhiFourierPeak[numPtBins];
    TH1D* dPhiFourierSide[numPtBins];

    TH1D* dPhiPeak[numPtBins];
    TH1D* dPhiSide[numPtBins];

    TF1* FourierFitXi[numPtBins];
    std::vector<double> v2values_peak;
    std::vector<double> v2values_side;
    std::vector<double> v2errors_peak;
    std::vector<double> v2errors_side;
    std::vector<double> v2value_h; //Need to use vector for vnCalculate function
    std::vector<double> v2error_h;
    std::vector<double> AvgKetXi;


    TLatex* ltx2 = new TLatex(  );
    ltx2->SetTextSize( 0.045 );
    ltx2->SetNDC( kTRUE );

    //FITTING FOR V2
    //
    //Define divided hist
    bool Peak = true;
    //bool Peak = false;
    for( int i=0; i<numPtBins; i++ )
    {
        //================================================================================
        //KET Calculations
        //================================================================================
        TH1D* hKetXi = (TH1D*)f->Get(Form("xiCorrelationRapidity/KET_xi_pt%d",i));
        TH1D* hKetXi_bkg = (TH1D*)f->Get(Form("xiCorrelationRapidity/KET_xi_bkg_pt%d",i));

        int nEntries = 0;
        double KetTotal = 0;
        for(int j=hKetXi->FindFirstBinAbove(0,1); j<=hKetXi->FindLastBinAbove(0,1); j++)
        {
            double nKet = hKetXi->GetBinContent(j);
            double Ket = nKet*(hKetXi->GetBinCenter(j));
            nEntries+=nKet;
            KetTotal += Ket;
        }
        for(int j=hKetXi_bkg->FindFirstBinAbove(0,1); j<=hKetXi_bkg->FindLastBinAbove(0,1); j++)
        {
            double nKet_bkg = hKetXi_bkg->GetBinContent(j);
            double Ket_bkg = nKet_bkg*(hKetXi_bkg->GetBinCenter(j));
            nEntries += nKet_bkg;
            KetTotal += Ket_bkg;
        }
        AvgKetXi.push_back(KetTotal/nEntries);

            dPhiPeak[i] = new TH1D( Form( "dPhiPeak%d",i ), "#Xi - h^{#pm} ", 31, -( 0.5 -
                        1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI  );
            TH1D *dPhiHad = new TH1D( "dPhiHad", "h^{#pm}- h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );
            //Pull 2D Histograms
            TH2D *hbackgroundPeak = (TH2D*) f->Get( Form( "xiCorrelationRapidity/BackgroundPeak_pt%d",i ) );
            TH2D *hsignalPeak     = (TH2D*) f->Get( Form( "xiCorrelationRapidity/SignalPeak_pt%d",i ) );
            TH2D *hBackgroundHad  = (TH2D*) fhad->Get( "xiCorrelationRapidity/BackgroundHad" );
            TH2D *hSignalHad      = (TH2D*) fhad->Get( "xiCorrelationRapidity/SignalHad" );

            TH1::SetDefaultSumw2(  );
            //Project Phi

            // For projecting both shoulders
            TH1D* hbPhiTotPeak = hbackgroundPeak->ProjectionY( "PhiBkgTotPeak", 1, binlow );
            TH1D* hbPhiOthPeak = hbackgroundPeak->ProjectionY( "PhiBkgOthPeak", binhigh, -1 );
            TH1D* hsPhiTotPeak = hsignalPeak->ProjectionY( "PhiSigTotPeak", 1, binlow );
            TH1D* hsPhiOthPeak = hsignalPeak->ProjectionY( "PhiSigOthPeak", binhigh, -1 );
            TH1D* hbHadPhiTot = hBackgroundHad->ProjectionY( "PhiBkgHadTot", 1, binlow );
            TH1D* hbHadPhiOth = hBackgroundHad->ProjectionY( "PhiBkgHadOth", binhigh, -1 );
            TH1D* hsHadPhiTot = hSignalHad->ProjectionY( "PhiSigHadTot", 1, binlow );
            TH1D* hsHadPhiOth = hSignalHad->ProjectionY( "PhiSigHadOth", binhigh, -1 );

            hbPhiTotPeak->Add( hbPhiOthPeak );
            hsPhiTotPeak->Add( hsPhiOthPeak );

            hbHadPhiTot->Add( hbHadPhiOth );
            hsHadPhiTot->Add( hsHadPhiOth );


            //Divide
            dPhiPeak[i]->Divide( hsPhiTotPeak, hbPhiTotPeak );
            dPhiHad->Divide( hsHadPhiTot, hbHadPhiTot );

            //Clone histograms for display without fit functions
            dPhiFourierPeak[i] = ( TH1D* )dPhiPeak[i]->Clone(  );
            TH1D* dPhiHadFourier = ( TH1D* )dPhiHad->Clone(  );

            FourierFitXi[i] = new TF1( Form( "FourierFitXi%d",i ), FourierHad, -1.5, 5, 4 );
            FourierFitXi[i]->SetNpx( 250 );
            FourierFitXi[i]->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *FourierPeak = new TCanvas( "FourierPeak", "Fourier Peak", 800,800 );
            FourierPeak->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            //dPhiFourierPeak->Fit( "FourierFitXi","","",0,PI );
            dPhiFourierPeak[i]->Fit( Form( "FourierFitXi%d",i ) );
            dPhiFourierPeak[i]->SetStats( kFALSE );
            v2values_peak.push_back( FourierFitXi[i]->GetParameter( 2 ) );
            v2errors_peak.push_back( FourierFitXi[i]->GetParError( 2 ) );
            cout << "---------------------------------" << endl;
            cout << "Peak V2 for xi-h is " << FourierFitXi[i]->GetParameter( 2 ) << endl;
            cout << "---------------------------------" << endl;

            TF1 *FourierFitHad = new TF1( "FourierFitHad", FourierHad, -1.5, 5, 4 );
            FourierFitHad->SetNpx( 250 );
            FourierFitHad->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *FourierHadron = new TCanvas( "FourierHadron", "Fourier Hadron", 800,800 );
            FourierHadron->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            //dPhiHadFourier->Fit( "FourierFitHad","","",0,PI );
            dPhiHadFourier->Fit( "FourierFitHad");
            dPhiHadFourier->SetStats( kFALSE );
            v2value_h.push_back(FourierFitHad->GetParameter(2));
            v2error_h.push_back(FourierFitHad->GetParError(2));

            double maxBinContent = dPhiFourierPeak[i]->GetBinContent( dPhiFourierPeak[i]->GetMaximumBin(  ) );
            double minBinContent = dPhiFourierPeak[i]->GetBinContent( dPhiFourierPeak[i]->GetMinimumBin(  ) );
            double minRange = minBinContent - 0.005*minBinContent;
            double maxRange = minRange + 2*( maxBinContent - minBinContent );

            //ZYAM FITS
            TCanvas *ZYAMFitPeak = new TCanvas( "ZYAMFitPeak", "ZYAM Fit Peak", 800,800 );
            ZYAMFitPeak->cd(  );

            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiPeak[i]->SetMarkerStyle( 21 );
            dPhiPeak[i]->SetMarkerColor( 4 );
            dPhiPeak[i]->SetTitleOffset( 2, "Y" );
            dPhiPeak[i]->SetTitle( "Peak" );
            dPhiPeak[i]->GetYaxis(  )->SetRangeUser( minRange , maxRange );
            dPhiPeak[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiPeak[i]->GetYaxis(  )->CenterTitle( true );
            dPhiPeak[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} " );
            dPhiPeak[i]->SetTitleOffset( 1.5, "X" );
            dPhiPeak[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiPeak[i]->GetXaxis(  )->CenterTitle( true );
            dPhiPeak[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
            //dPhiPeak[i]->Draw( "E1" );
            dPhiPeak[i]->Fit( "pol2","","", 0.4,2.4 );
            dPhiPeak[i]->SetStats( !publish );

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

            maxBinContent = dPhiHadFourier->GetBinContent( dPhiHadFourier->GetMaximumBin(  ) );
            minBinContent = dPhiHadFourier->GetBinContent( dPhiHadFourier->GetMinimumBin(  ) );
            minRange = minBinContent - 0.005*minBinContent;
            maxRange = minRange + 2*( maxBinContent - minBinContent );

            TCanvas *c3 = new TCanvas( "c3", "", 800,800 );
            c3->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiHad->SetMarkerStyle( 34 );
            dPhiHad->SetMarkerSize( 1.5 );
            dPhiHad->SetTitleOffset( 2, "Y" );
            dPhiHad->SetTitle( "" );
            dPhiHad->GetYaxis(  )->SetRangeUser( minRange , maxRange );
            dPhiHad->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiHad->GetYaxis(  )->CenterTitle( true );
            dPhiHad->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi}" );
            dPhiHad->SetTitleOffset( 1.5, "X" );
            dPhiHad->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiHad->GetXaxis(  )->CenterTitle( true );
            dPhiHad->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
            //dPhiHad->Draw( "hist E1" );
            dPhiHad->Fit( "pol2","","", 0.4,2 );
            dPhiHad->SetStats( !publish );

            if( publish )
            {
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.74, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }

            TF1 *dPhiFitPeak = dPhiPeak[i]->GetFunction( "pol2" );
            TF1 *dPhiHadFit = dPhiHad->GetFunction( "pol2" );

            double dPhiFitMinPeak = dPhiFitPeak->GetMinimum(  );
            double dPhiHadFitMin = dPhiHadFit->GetMinimum(  );

            for( int j = 1; j < 32; j++ )
            {
                dPhiFourierPeak[i]->AddBinContent( j, -dPhiFitMinPeak );
            }

            for( int j = 1; j < 32; j++ )
            {
                dPhiHadFourier->AddBinContent( j, -dPhiHadFitMin );
            }

            FourierPeak->cd(  );
            dPhiFourierPeak[i]->SetMarkerStyle( 21 );
            dPhiFourierPeak[i]->SetMarkerColor( 4 );
            //dPhiFourierPeak[i]->Draw( "E1" );
            dPhiFourierPeak[i]->SetStats( kFALSE );
            dPhiFourierPeak[i]->Fit( Form( "FourierFitXi%d",i ) );
            os << "Peak " << PtBin[i] << "_Pt_" << PtBin[i+1];
            dPhiFourierPeak[i]->SetTitle( os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            dPhiFourierPeak[i]->SetTitleOffset( 2, "Y" );
            dPhiFourierPeak[i]->GetYaxis(  )->CenterTitle( true );
            dPhiFourierPeak[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiFourierPeak[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
            dPhiFourierPeak[i]->GetYaxis(  )->SetRangeUser( -0.0004, 0.008 );
            dPhiFourierPeak[i]->SetTitleOffset( 1.5, "X" );
            dPhiFourierPeak[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiFourierPeak[i]->GetXaxis(  )->CenterTitle( true );
            dPhiFourierPeak[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );


            if( publish )
            {
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }

            FourierHadron->cd(  );
            dPhiHadFourier->SetMarkerStyle( 34 );
            dPhiHadFourier->SetMarkerSize( 1.5 );
            dPhiHadFourier->Draw( "E1" );
            dPhiHadFourier->SetStats( kFALSE );
            dPhiHadFourier->Fit( "FourierFitHad" );
            dPhiHadFourier->SetTitle( "" );
            dPhiHadFourier->SetTitleOffset( 2, "Y" );
            dPhiHadFourier->GetYaxis(  )->CenterTitle( true );
            dPhiHadFourier->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiHadFourier->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
            dPhiHadFourier->GetYaxis(  )->SetRangeUser( -0.0004 , 0.008);
            dPhiHadFourier->SetTitleOffset( 1.5, "X" );
            dPhiHadFourier->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiHadFourier->GetXaxis(  )->CenterTitle( true );
            dPhiHadFourier->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );

            if( publish )
            {
                ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 5.02 TeV" );
                ltx3->DrawLatex( 0.2, 0.74, "L_{#lower[-0.25]{int}} = 35 nb^{-1}" );
                ltx3->DrawLatex( 0.2, 0.67,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
                ltx3->DrawLatex( 0.2, 0.60, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
                ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
                ltx2->DrawLatex( 0.7, 0.74, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
            }
            //side
            dPhiSide[i] = new TH1D( Form( "dPhiSide%d",i ), "#Xi - h^{#pm} ", 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI  );
            TH2D *hbackgroundSide = (TH2D*) f->Get( Form( "xiCorrelationRapidity/BackgroundSide_pt%d",i ) );
            TH2D *hsignalSide     = (TH2D*) f->Get( Form( "xiCorrelationRapidity/SignalSide_pt%d",i ) );

            TH1::SetDefaultSumw2(  );

            //Project Phi
            TH1D* hbPhiTotSide = hbackgroundSide->ProjectionY( "PhiBkgTot", 0, 10 );
            TH1D* hbPhiOthSide = hbackgroundSide->ProjectionY( "PhiBkgOthPeak", 23, -1 );
            TH1D* hsPhiTotSide = hsignalSide->ProjectionY( "PhiSigTot", 0, 10 );
            TH1D* hsPhiOthSide = hsignalSide->ProjectionY( "PhiSigOthPeak", 23, -1 );

            hbPhiTotSide->Add( hbPhiOthSide );
            hsPhiTotSide->Add( hsPhiOthSide );

            hbHadPhiTot->Add( hbHadPhiOth );
            hsHadPhiTot->Add( hsHadPhiOth );

            //Divide
            dPhiSide[i]->Divide( hsPhiTotSide, hbPhiTotSide );

            //Clone histograms for display without fit functions
            dPhiFourierSide[i] = ( TH1D* )dPhiSide[i]->Clone(  );

            FourierFitXi[i] = new TF1( Form( "FourierFitXi%d",i ) , FourierHad, -1.5, 5, 4 );
            FourierFitXi[i]->SetNpx( 250 );
            FourierFitXi[i]->SetParNames( "Scale", "V_{1}", "V_{2}", "V_{3}" );

            TCanvas *FourierSide = new TCanvas( "FourierSide", "Fourier Side", 800,800 );
            FourierSide->cd(  );
            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiFourierSide[i]->Fit( Form( "FourierFitXi%d",i ) );
            dPhiFourierSide[i]->SetStats( kFALSE );
            v2values_side.push_back( FourierFitXi[i]->GetParameter( 2 ) );
            v2errors_side.push_back( FourierFitXi[i]->GetParError( 2 ) );
            cout << "---------------------------------" << endl;
            cout << "Side V2 for xi-h is " << FourierFitXi[i]->GetParameter( 2 ) << endl;
            cout << "---------------------------------" << endl;

            maxBinContent = dPhiFourierSide[i]->GetBinContent( dPhiFourierSide[i]->GetMaximumBin(  ) );
            minBinContent = dPhiFourierSide[i]->GetBinContent( dPhiFourierSide[i]->GetMinimumBin(  ) );
            minRange = minBinContent - 0.005*minBinContent;
            maxRange = minRange + 2*( maxBinContent - minBinContent );


            //ZYAM FITS
            TCanvas *ZYAMFitSide = new TCanvas( "ZYAMFitSide", "ZYAM Fit Side", 800,800 );
            ZYAMFitSide->cd(  );

            gPad->SetTickx(  );
            gPad->SetTicky(  );
            dPhiSide[i]->SetMarkerStyle( 21 );
            dPhiSide[i]->SetMarkerColor( 4 );
            dPhiSide[i]->SetTitleOffset( 2, "Y" );
            dPhiSide[i]->SetTitle( "Sideband" );
            dPhiSide[i]->GetYaxis(  )->SetRangeUser( minRange , maxRange );
            dPhiSide[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiSide[i]->GetYaxis(  )->CenterTitle( true );
            dPhiSide[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} " );
            dPhiSide[i]->SetTitleOffset( 1.5, "X" );
            dPhiSide[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiSide[i]->GetXaxis(  )->CenterTitle( true );
            dPhiSide[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );
            //dPhiSide[i]->Draw( "E1" );
            dPhiSide[i]->Fit( "pol2","","", 0.4,2.4 );
            dPhiSide[i]->SetStats( !publish );

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

            TF1 *dPhiFitSide = dPhiSide[i]->GetFunction( "pol2" );

            double dPhiFitMinSide = dPhiFitSide->GetMinimum(  );

            for( int j = 1; j < 32; j++ )
            {
                dPhiFourierSide[i]->AddBinContent( j, -dPhiFitMinSide );
            }


            FourierSide->cd(  );
            dPhiFourierSide[i]->SetMarkerStyle( 21 );
            dPhiFourierSide[i]->SetMarkerColor( 4 );
            //dPhiFourierSide[i]->Draw( "E1" );
            dPhiFourierSide[i]->SetStats( kFALSE );
            dPhiFourierSide[i]->Fit( Form( "FourierFitXi%d",i ) );
            os << "SideBand " << PtBin[i] << "_Pt_" << PtBin[i+1];
            dPhiFourierSide[i]->SetTitle( os.str(  ).c_str(  ) );
            os.str( std::string(  ) );
            dPhiFourierSide[i]->SetTitleOffset( 2, "Y" );
            dPhiFourierSide[i]->GetYaxis(  )->CenterTitle( true );
            dPhiFourierSide[i]->GetYaxis(  )->SetTitleSize( 0.03 );
            dPhiFourierSide[i]->GetYaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{dN^{pair}}{d#Delta#phi} - C_{#lower[-0.3]{ZYAM}}" );
            dPhiFourierSide[i]->GetYaxis(  )->SetRangeUser( -0.0004, 0.008 );
            dPhiFourierSide[i]->SetTitleOffset( 1.5, "X" );
            dPhiFourierSide[i]->GetXaxis(  )->SetTitleSize( 0.035 );
            dPhiFourierSide[i]->GetXaxis(  )->CenterTitle( true );
            dPhiFourierSide[i]->GetXaxis(  )->SetTitle( "#Delta#phi (radians)" );

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

    OutputVnValues(2,"Peak","Xi",v2values_peak,PtBin,Xiv2PeakName);
    OutputVnErrors(2,"Peak","Xi",v2errors_peak,PtBin,Xiv2PeakName);

    OutputVnValues(2,"Side","Xi",v2values_side,PtBin,Xiv2SideName);
    OutputVnErrors(2,"Side","Xi",v2errors_side,PtBin,Xiv2SideName);

    //Calculate Flow
    vnCalculate(2,"Xi",v2values_peak,v2errors_peak,v2values_side,v2errors_side,v2value_h,v2error_h,fsig_xi,Xiv2CalculatorName);

    std::map<std::string, std::vector<double> > results_xi = vnCalculateMap(2,"Xi",v2values_peak,v2errors_peak,v2values_side,v2errors_side,v2value_h,v2error_h,fsig_xi);

    std::vector<double> AvgX_xi = AvgX(f,branchname_xi,branchname_xi_bkg,numPtBins);

    std::vector<double> AvgX_ket_xi = AvgX(f,branchname_ket_xi,branchname_ket_xi_bkg,numPtBins);

    vnGraph(results_xi,AvgX_xi,AvgX_ket_xi,"Xi",graphName);

    std::ofstream XiKET;
    XiKET.open("XiAvgKET.txt");
    XiKET << "Avg Ket values\n";
    for(unsigned i=0; i<AvgKetXi.size(); i++)
        XiKET << AvgKetXi[i] << "\n";

    XiKET << "Avg Ket/nq values\n";
    for(unsigned i=0; i<AvgKetXi.size(); i++)
        XiKET << AvgKetXi[i]/3 << "\n";

    TCanvas* FourierPeakComp = new TCanvas( "FourierPeakComp", "Fourier Peak Composite", 1200,1000 );
    FourierPeakComp->Divide( 3,3 );
    for( int i=0; i<numPtBins; i++ ){
        FourierPeakComp->cd( i+1 );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        dPhiFourierPeak[i]->Draw( "E1" );
    }
    TCanvas* FourierSideComp = new TCanvas( "FourierSideComp", "Fourier Side Composite", 1200,1000 );
    FourierSideComp->Divide( 3,3 );
    for( int i=0; i<numPtBins; i++ ){
        FourierSideComp->cd( i+1 );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        dPhiFourierSide[i]->Draw( "E1" );
    }

    //Output Publication plots
    //1D correlation functions
    if( Peak ){
        TCanvas* PubFourier = new TCanvas( "PubFourier", "Pub", 800,800 );
        PubFourier->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        TH1D* dPhiFourierPeakCopy = ( TH1D* )dPhiFourierPeak[4]->Clone(  );
        dPhiFourierPeakCopy->SetTitle( "Peak" );
        dPhiFourierPeakCopy->Draw( "E1" );

        TLatex *ltx3 = new TLatex(  );
        ltx3->SetTextSize( 0.035 );
        ltx3->SetNDC( kTRUE );
        ltx3->SetTextFont( 42 );
        if( publish )
        {
            ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}" );
            ltx3->DrawLatex( 0.2, 0.74,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            ltx3->DrawLatex( 0.2, 0.67, "0.3 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
            ltx3->DrawLatex( 0.2, 0.60, "2.8 < p_{T}^{#Xi} < 3.6 GeV" );
            ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
            ltx3->SetTextSize( 0.045 );
            ltx3->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
        }
    }
    else{
        TCanvas* PubFourier = new TCanvas( "PubFourier", "Pub", 800,800 );
        PubFourier->cd(  );
        gPad->SetTickx(  );
        gPad->SetTicky(  );
        TH1D* dPhiFourierSideCopy = ( TH1D* )dPhiFourierSide[4]->Clone(  );
        dPhiFourierSideCopy->SetTitle( "SideBand" );
        dPhiFourierSideCopy->Draw( "E1" );


        TLatex *ltx3 = new TLatex(  );
        ltx3->SetTextSize( 0.035 );
        ltx3->SetNDC( kTRUE );
        ltx3->SetTextFont( 42 );
        if( publish )
        {
            ltx3->DrawLatex( 0.2, 0.82, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}" );
            ltx3->DrawLatex( 0.2, 0.74,"185  #leq  N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
            ltx3->DrawLatex( 0.2, 0.67, "0.3 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
            ltx3->DrawLatex( 0.2, 0.60, "2.8 < p_{T}^{#Xi} < 3.6 GeV" );
            ltx3->DrawLatex( 0.2, 0.53, "Long range (|#Delta#eta| > 2)" );
            ltx3->SetTextSize( 0.045 );
            ltx3->DrawLatex( 0.7, 0.74, "#Xi#kern[-0.3]{#lower[0.2]{{}^{#pm}}}- h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}" );
        }
    }

    //2D Correlation function 1-3 GeV associated
    //TCanvas* TwoDCorrelation = new TCanvas( "TwoDCorrelation", "", 1000, 1000 );
    //TwoDCorrelation->SetLeftMargin( 0.2 );

    //TH2D* Signal = ( TH2D* )f->Get( "xiCorrelationRapidity/SignalXiHad" );
    //TH2D* Background = ( TH2D* )f->Get( "xiCorrelationRapidity/BackgroundXiHad" );

    //TGaxis::SetMaxDigits( 1 );
    
    //TH2D* Correlation = ( TH2D* )Signal->Clone(  );
    //Correlation->Divide( Background );
    //Correlation->GetXaxis(  )->SetRangeUser( -4.0,4.0 );
    //Correlation->GetYaxis(  )->SetRangeUser( -PI/2.0,4.5 );
    //Correlation->GetXaxis(  )->SetTitle( "#Delta#eta" );
    //Correlation->GetXaxis(  )->SetTitleOffset( 1.4 );
    //Correlation->GetXaxis(  )->CenterTitle( true );
    //Correlation->GetYaxis(  )->SetTitle( "#Delta#phi (radians)" );
    //Correlation->GetYaxis(  )->SetTitleOffset( 1.4 );
    //Correlation->GetYaxis(  )->CenterTitle( true );
    //Correlation->GetZaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{d^{2}N^{pair}}{d#Delta#eta d#Delta#phi}" );
    //Correlation->GetZaxis(  )->SetTitleOffset( 2.3 );
    //Correlation->GetZaxis(  )->CenterTitle( true );
    //Correlation->GetXaxis(  )->SetNdivisions( 405 );
    //Correlation->GetYaxis(  )->SetNdivisions( 405 );
    //Correlation->GetZaxis(  )->SetNdivisions( 4 );
    //Correlation->SetTitle( "" );
    //Correlation->SetStats( kFALSE );
    //Correlation->Scale( 10 );

    //const Int_t NRGBs = 5;
    //const Int_t NCont = 20;

    //Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    //Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    //Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    //Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    //TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    
    //TwoDCorrelation->cd(  );
    //Correlation->Draw( "SURF1 FB " );

    //TLatex *ltx0 = new TLatex(  );
    //ltx0->SetTextSize( 0.031 );
    //ltx0->SetNDC( kTRUE );
    //ltx0->SetTextFont( 42 );

    //ltx0->DrawLatex( 0.05, 0.95, "CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV, L_{#lower[-0.25]{int}} = 62 nb^{-1}" );
    //ltx0->DrawLatex( 0.05, 0.88, "185 #leq N_{#lower[-0.3]{trk}}#kern[-0.47]{#lower[0.1]{{}^{offline}}}< 220" );
    //ltx0->DrawLatex( 0.05, 0.81, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{assoc}}} < 3 GeV" );
    //ltx0->DrawLatex( 0.05, 0.75, "1 < p_{T}#kern[-0.3]{#lower[0.1]{{}^{trig}}} < 3 GeV" );
    //ltx0->SetTextSize( 0.04 );
    //ltx0->DrawLatex( 0.85, 0.88, "#Xi#kern[-0.3]{#lower[0.02]{{}^{#pm}}}- h^{#pm}" );

    //TCanvas* SigAndBkg = new TCanvas( "SigAndBkg", "", 1600, 800 );
    //SigAndBkg->Divide( 2,1 );
    //SigAndBkg->SetLeftMargin( 0.2 );

    //SigAndBkg->cd( 1 );
    //Signal->GetXaxis(  )->SetRangeUser( -4.0,4.0 );
    //Signal->GetYaxis(  )->SetRangeUser( -PI/2.0,4.5 );
    //Signal->GetXaxis(  )->SetTitle( "#Delta#eta" );
    //Signal->GetXaxis(  )->SetTitleOffset( 1.4 );
    //Signal->GetXaxis(  )->CenterTitle( true );
    //Signal->GetYaxis(  )->SetTitle( "#Delta#phi (radians)" );
    //Signal->GetYaxis(  )->SetTitleOffset( 1.4 );
    //Signal->GetYaxis(  )->CenterTitle( true );
    //Signal->GetZaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{d^{2}N^{pair}}{d#Delta#eta d#Delta#phi}" );
    //Signal->GetZaxis(  )->SetTitleOffset( 2.5 );
    //Signal->GetZaxis(  )->CenterTitle( true );
    //Signal->GetXaxis(  )->SetNdivisions( 405 );
    //Signal->GetYaxis(  )->SetNdivisions( 405 );
    //Signal->GetZaxis(  )->SetNdivisions( 4 );
    //Signal->SetTitle( "" );
    //Signal->SetStats( kFALSE );
    //Signal->Draw( "Surf1 FB" );

    //SigAndBkg->cd( 2 );
    //Background->GetXaxis(  )->SetRangeUser( -4.0,4.0 );
    //Background->GetYaxis(  )->SetRangeUser( -PI/2.0,4.5 );
    //Background->GetXaxis(  )->SetTitle( "#Delta#eta" );
    //Background->GetXaxis(  )->SetTitleOffset( 1.4 );
    //Background->GetXaxis(  )->CenterTitle( true );
    //Background->GetYaxis(  )->SetTitle( "#Delta#phi (radians)" );
    //Background->GetYaxis(  )->SetTitleOffset( 1.4 );
    //Background->GetYaxis(  )->CenterTitle( true );
    //Background->GetZaxis(  )->SetTitle( "#frac{1}{N_{#lower[-0.3]{trig}}} #frac{d^{2}N^{pair}}{d#Delta#eta d#Delta#phi}" );
    //Background->GetZaxis(  )->SetTitleOffset( 2.5 );
    //Background->GetZaxis(  )->CenterTitle( true );
    //Background->GetXaxis(  )->SetNdivisions( 405 );
    //Background->GetYaxis(  )->SetNdivisions( 405 );
    //Background->GetZaxis(  )->SetNdivisions( 4 );
    //Background->SetTitle( "" );
    //Background->SetStats( kFALSE );
    //Background->Draw( "Surf1 FB" );

    

}

