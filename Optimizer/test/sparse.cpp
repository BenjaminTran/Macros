#include <TLatex.h>
#include <TStyle.h>
#include "TF1.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TLegend.h"
#include "TMathText.h"
#include "TImage.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include <TString.h>
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
#include "TGaxis.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TROOT.h"

//#include <vector>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <sys/stat.h> //For file existance checking


void sparse()
{
    const int nDim = 3;
    int nBins[nDim] = {20,4,4};
    double mins[nDim] = {0,0,0};
    double maxs[nDim] = {1,1,1};
    THnSparseF* h = new THnSparseF("sparseHist","sparseHist",nDim,nBins,mins,maxs);
    const int nVarBins = 20;
    double varBins[nVarBins+1] = {0., 0.01, 0.1, 0.11, 0.2, 0.21, 0.30, 0.31, 0.40, 0.41, 0.5, 0.51, 0.60, 0.61, 0.70, 0.71, 0.80, 0.81, 0.90, 0.91, 1.0};
    h->GetAxis(0)->Set(nVarBins, varBins);
    double val1[nDim] = {0.01,0.1,1.5};
    double val2[nDim] = {0.3,0.9,0.1};
    double val3[nDim] = {0.9,0.1,0.6};
    h->Fill(val1);
    h->Fill(val1);
    h->Fill(val2);
    h->Fill(val3);
    h->PrintEntries();

    TFile f("test.root","RECREATE");
    TCanvas* c1 = new TCanvas("c1","", 600,600);
    h->Write();
    f.Close();
    //h->GetAxis(1)->SetRange(0,5);
    h->GetAxis(0)->SetRange(2,20);
    h->GetAxis(1)->SetRange(1,3);
    h->Projection(0)->Draw();
    //TH1D* i = (TH1D*)h->Projection(2);
    //i->GetXaxis()->SetRange(0,5);
    //i->Draw();
}

void readSparse()
{
    TFile* f = new TFile("test.root");
    THnSparseF* hh = (THnSparseF*)f->Get("sparseHist");
    hh->PrintEntries();
}
