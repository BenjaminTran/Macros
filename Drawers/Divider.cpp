#include "TH1D.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGaxis.h"

#define PI 3.1416
void Divider(  )
{
    TH1::SetDefaultSumw2(  );

    //TFile* f = new TFile( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/results/NoPtCutPeakAndSide/XiAnalysisSeparated.root" );
    TFile* f = new TFile( "/volumes/MacHD/Users/blt1/research/CascadeV2pPb/RootFiles/Flow/Thesis/XiAnalysisCorrelationPtCut8TeVPD1_4_ForFinal.root");

    TH1D* InvMass = ( TH1D* )f->Get( "xiCorrelation/InvMassXi" );
    TH2D* SignalPeak = ( TH2D* )f->Get( "xiCorrelation/SignalPeak" );
    TH2D* BackgroundPeak = ( TH2D* )f->Get( "xiCorrelation/BackgroundPeak" );
    TH2D* CorrelationRoot = ( TH2D* )f->Get( "xiCorrelation/CorrelationPeak" );
    TH2D* CorrelationPeak = new TH2D( "CorrelationPeak", "Correlation Peak", 33, -4.95, 4.95, 31, -( 0.5 - 1.0/32 )*PI, ( 1.5 - 1.0/32 )*PI );



    TCanvas* c1 = new TCanvas( "Inv", "",800,800 );
    c1->cd(  );
    InvMass->Draw(  );

    TCanvas* c2 = new TCanvas( "Corr", "",800,800 );
    c2->cd(  );
    CorrelationPeak->Add( SignalPeak );
    CorrelationPeak->Divide( BackgroundPeak );
    CorrelationPeak->Draw( "Surf1" );

    TCanvas* c3 = new TCanvas( "Root", "",800,800 );
    c3->cd(  );
    CorrelationRoot->Draw( "Surf1" );

}
