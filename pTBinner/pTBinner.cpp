#include "TH1D.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGaxis.h"
#include <iostream>

void pTBinner(  )
{
    TH1::SetDefaultSumw2(  );

    TFile* f = new TFile( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/8TeV/XiAnalysisNoPtCut8TeVPD1246Incomplete.root" );

    TH1D* pTSpectrum = ( TH1D* )f->Get( "xiCorrelation/pT_Xi" );

    //double pTintegration[]

    double pT1 = 0;
    double pT2 = 0;
    double pT3 = 0;
    double pT4 = 0;

    pT1  = pTSpectrum->Integral( 7,16 );
    pT2 = pTSpectrum->Integral( 17,19 );
    pT3  = pTSpectrum->Integral( 20,21 );
    pT4 = pTSpectrum->Integral( 22,24 );

    cout << "1 to 1.6: " << pT1 << endl;
    cout << "1.7 to 1.9: " << pT2 << endl;
    cout << "2 to 2.2: " << pT3 << endl;
    cout << "2.3 to 2.4: " << pT4 << endl;

    TCanvas* c1 = new TCanvas( "pT", "", 800,800 );
    c1->cd(  );
    gPad->SetTickx(  );
    gPad->SetTicky(  );
    pTSpectrum->Draw(  );
}
