#include "TH2.h"

void histInv(){
    TFile* f = TFile::Open("/Users/btran/research/RootFiles/Phi/BDT/PhiBDT_30M_rapLT1_DataBkg_270618.root");
    TH2D* h = (TH2D*)f->Get("BDTApp/dedx");

    TH2D* hinvert = new TH2D("dedx","Dedx v momentum",150,0,10,200,0,5);

    for(int x=1; x<201; ++x)
    {
        for(int y=1; y<151; ++y)
        {
            double content = h->GetBinContent(x,y);
            hinvert->SetBinContent(y,x,content);
        }
    }
    hinvert->Draw();
}
