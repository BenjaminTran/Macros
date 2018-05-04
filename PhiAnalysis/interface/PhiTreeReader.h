#ifndef PHIANALSIS__PHITREEREADER_H
#define PHIANALSIS__PHITREEREADER_H
//Includes
#include <TStyle.h>
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TMathText.h"
#include "TImage.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THnSparse.h"
#include "TTree.h"
#include <TString.h>
#include "TStyle.h"
#include "TString.h"
#include "TGaxis.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
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
#include "RooBifurGauss.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "TString.h"
#include "TGaxis.h"
#include "TClass.h"
#include "TDataType.h"
#include "TBenchmark.h"
//#include "TROOT.h"

//#include <vector>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <sys/stat.h> //For file existance checking
#include <chrono>

//typedef TTreeReaderValue<auto> AnyReaderValue;

class PhiTreeReader{
    private :

        long int event_start;
        long int event_stop;
        int check_interval;

    public :
        static constexpr double kaonMass = 0.493677;

        struct kaon{
            double p;
            double pt;
            double px;
            double py;
            double pz;
            double dedx;
            double energy;
            int charge;
            int nhits;

            kaon(double p_, double pt_, double px_, double py_, double pz_, double dedx_, double energy_, int charge_, int nhits_) :
                p(p_), pt(pt_), px(px_), py(py_), pz(pz_), dedx(dedx_), energy(energy_), charge(charge_), nhits(nhits_)  {}
        };

        PhiTreeReader(long int event_start_, long int event_stop_, int check_interval_, std::string filename_);
        ~PhiTreeReader();


        std::string filename;
        static bool CheckValue(ROOT::Internal::TTreeReaderValueBase& value);

        static TTreeReaderValue<float> SetValues(TTreeReader& reader, std::string branchName );

        static TTreeReaderValue<int> SetValuesInt(TTreeReader& reader, std::string branchName );

        static void SparseVarSetUp(std::vector<float> &varBins_, int nBins_, float minBin_, float maxBin_);

        static void CombinatorialMass(std::vector<kaon> PKp, std::vector<kaon> PKm, TH1D* h_mass_);

        static std::vector<double> CombinatorialMassSparse(std::vector<kaon> PKp, std::vector<kaon> PKm);

        bool ReadTreeV2();
        bool ReadTreeV3();
        static void ReadSparse();

        //std::map<std::string, TTreeReaderValue<auto> > SetTTreeReaderValue(TFile* file, std::string path_to_tree, TTreeReader &reader);

        //Member variables

        TH2D* h_DeDx;
        TH1D* h_mass;
        TH1I* h_nEvt;

        std::vector<PhiTreeReader::kaon> Pkp;
        std::vector<PhiTreeReader::kaon> Pkm;

        double cut_pt;
        double cut_eta;
        double cut_rapidity;
        double cut_nhits;
        double cut_DCAz;
        double cut_DCAxy;
        double cut_mult;

};

#endif
