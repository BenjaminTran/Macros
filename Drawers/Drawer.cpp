#include <iostream>
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

    //TFile* f1 = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/MC/V0/MCMassPtTotal_08_23_2017.root");

    //TH3D* XiMassPtRap = (TH3D*)f1->Get("MassPtRapidity/XiMassPtRap");
    //TH3D* LaMassPtRap = (TH3D*)f1->Get("MassPtRapidity/LaMassPtRap");
    //TH3D* KsMassPtRap = (TH3D*)f1->Get("MassPtRapidity/KsMassPtRap");

    //TH1D* MassLa = LaMassPtRap->ProjectionX( "MassLa","", )

    //TCanvas* c1 = new TCanvas( "c1","",600,600 );

    //XiMassPt->Draw();

    TFile *f = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityTotal_08_21_2017.root");
    //================================================================================
    //KET Calculations
    //================================================================================
    int i = 13;
    TH1D* hKetKs = (TH1D*)f->Get("v0CorrelationRapidity/KETkshort_pt10");
    TH1D* hKetKs_bkg = (TH1D*)f->Get("v0CorrelationRapidity/KETkshort_bkg_pt10");

    cout << "First" <<  hKetKs->FindFirstBinAbove( 0 , 1) << endl;
    cout << "Last" <<  hKetKs->FindLastBinAbove( 0, 1 ) << endl;

    //int nEntries = 0;
    //double KetTotal = 0;
    //for(int j=hKetKs; j<20000; j++)
    //{
        //double nKet = hKetKs->GetBinContent(i);
        //double nKet_bkg = hKetKs_bkg->GetBinContent(i);
        //double Ket = nKet*(hKetKs->GetBinCenter(i));
        //double Ket_bkg = nKet_bkg*(hKetKs_bkg->GetBinCenter(i));
        //double TotalKet = Ket + Ket_bkg;
        //if(TotalKet == 0) continue;
        //KetTotal += TotalKet;
        //nEntries++;
    //}
    //AvgKetKs.push_back(KetTotal/nEntries);

}
