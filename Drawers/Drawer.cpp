#include "TH1D.h"
#include "TPad.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGaxis.h"

void Drawer()
{
    TH1::SetDefaultSumw2(  );
    gStyle->SetTitleFontSize(0.04);

    TFile* f1 = new TFile("/volumes/MacHD/Users/blt1/research/TestRootFiles/MCEff.root");

    TH3D* XiMassPtRap = (TH3D*)f1->Get("MassPtRapidity/XiMassPtRap");
    TH3D* LaMassPtRap = (TH3D*)f1->Get("MassPtRapidity/LaMassPtRap");
    TH3D* KsMassPtRap = (TH3D*)f1->Get("MassPtRapidity/KsMassPtRap");

    TH2D* XiMassPt = ( TH2D* )MCEff->Project3D( "yx" );

    TH2D* XiMassPt = ( TH2D* )MCEff->Project3D( "yx" );

    TCanvas* c1 = new TCanvas( "c1","",600,600 );
    XiMassPt->Draw();

}
