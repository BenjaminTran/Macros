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

    TFile *f_mb = TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/HadCorr/Combine_MB0_corr_ref.root");
    TFile *f_hm = TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationRapidityTotal_08_20_2017.root");

    TFile *f_perisub = TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/AllCorrelation/PeripheralSubtractionMB.root");
    TFile *f3 = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityCorrectMultB_09_19_17.root");

    TLatex* tex = new TLatex();
    tex->SetNDC();
    
}
