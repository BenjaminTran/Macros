#include "TFile.h"
#include "TMath.h"
#include "TGraph.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TH2.h"
#include "TVirtualFitter.h"
#include "interface/GetGraphFromFile.C"
#include "interface/MITStyle.C"

void Effgraph()
{
    std::vector<double> Effks0 = {0.00236573 ,0.0182213 ,0.0307441 ,0.0437233 ,0.0571268 ,0.0762893 ,0.102232 ,0.126588 ,0.146577 ,0.161527 ,0.175623 ,0.187691 ,0.196148 ,0.204111 ,0.20618 ,0.206229 ,0.202978 ,0.200138 ,0.202294 ,0.185897 ,0.174794 ,0.148988};
    std::vector<double> Effks1 = {0.00146942 ,0.0144038 ,0.0262371 ,0.0378985 ,0.0495121 ,0.067597 ,0.0904721 ,0.113617 ,0.135351 ,0.153586 ,0.168245 ,0.183168 ,0.192943 ,0.197841 ,0.218533 ,0.201523 ,0.205986 ,0.199435 ,0.199559 ,0.187679 ,0.17213 ,0.150061};
    std::vector<double> Effks2 = {0.00133345 ,0.0139583 ,0.0256052 ,0.0370206 ,0.048376 ,0.0629577 ,0.0850432 ,0.106269 ,0.127525 ,0.145293 ,0.160316 ,0.17429 ,0.18159 ,0.190282 ,0.19488 ,0.191657 ,0.196551 ,0.190946 ,0.184583 ,0.175483 ,0.154487 ,0.145738};
    std::vector<double> Effks3 = {0.00219414 ,0.0177678 ,0.0304313 ,0.0434582 ,0.0674569 ,0.0740474 ,0.0982712 ,0.120329 ,0.138868 ,0.155286 ,0.169243 ,0.180445 ,0.19088 ,0.196175 ,0.197966 ,0.202496 ,0.197688 ,0.199294 ,0.191349 ,0.183621 ,0.168921 ,0.148934};

    std::vector<double> Effla0 = {0.000525492 ,0.0094348 ,0.0231253 ,0.030661 ,0.0400555 ,0.0497609 ,0.0604912 ,0.0716253 ,0.0815167 ,0.0877368 ,0.0907619 ,0.0923013 ,0.0913848 ,0.0921291 ,0.0919683 ,0.0869069 ,0.0825778 ,0.0725699};
    std::vector<double> Effla1 = {0.000391787 ,0.00804342 ,0.0226002 ,0.0305117 ,0.0367709 ,0.0402254 ,0.0459731 ,0.0531688 ,0.0614725 ,0.0740074 ,0.0830869 ,0.0882421 ,0.0914201 ,0.086292 ,0.0868601 ,0.0818931 ,0.0756722 ,0.0744172};
    std::vector<double> Effla2 = {0.000305291 ,0.00743217 ,0.0220162 ,0.0300331 ,0.0362827 ,0.0385964 ,0.0431911 ,0.0483456 ,0.0558968 ,0.066205 ,0.0766199 ,0.0803229 ,0.0826549 ,0.0817624 ,0.0812538 ,0.0749496 ,0.0705298 ,0.067485};
    std::vector<double> Effla3 = {0.000521257 ,0.00924198 ,0.0231543 ,0.0308568 ,0.0378431 ,0.0439883 ,0.0532437 ,0.062991 ,0.0716101 ,0.080151 ,0.0855211 ,0.086855 ,0.0867216 ,0.0848629 ,0.0859271 ,0.0825444 ,0.0789478 ,0.0643987};

    //std::vector<double> Ptks = {0.5,0.6,0.7,0.8,0.9,1.1,1.3,1.5,1.7,1.9,2.2,2.5,2.8,3.1,3.4,3.7,4.0,4.4,5.0,6.0,7.0,10.0};
    std::vector<double> Ptks = {0.25,0.55,0.65,0.75,0.85,1.,1.2,1.4,1.6,1.8,2.15,2.35,2.65,2.95,3.25,3.55,3.85,4.2,4.7,5.5,6.5,8.5};
    //std::vector<double> Ptla = {0.8,1.1,1.4,1.7,2.0,2.3,2.6,2.9,3.2,3.6,4.0,4.4,4.8,5.2,5.6,6.4,8.0,10.0};
    std::vector<double> Ptla = {0.4,0.9,1.25,1.55,1.85,2.15,2.45,2.75,3.05,3.4,3.8,4.2,4.6,5.0,5.4,6.0,7.2,9.0};

    TCanvas* can[4];
    TCanvas* c1 = new TCanvas("c1","",1200,1200);
    c1->Divide(2,2);
    TH1F* frame[4];
    TGraphErrors* Effks[4];
    TGraphErrors* Effla[4];

    Effks[0] = new TGraphErrors(22,&Ptks[0],&Effks0[0],0,0);
    Effks[1] = new TGraphErrors(22,&Ptks[0],&Effks1[0],0,0);
    Effks[2] = new TGraphErrors(22,&Ptks[0],&Effks2[0],0,0);
    Effks[3] = new TGraphErrors(22,&Ptks[0],&Effks3[0],0,0);

    Effla[0] = new TGraphErrors(18,&Ptla[0],&Effla0[0],0,0);
    Effla[1] = new TGraphErrors(18,&Ptla[0],&Effla1[0],0,0);
    Effla[2] = new TGraphErrors(18,&Ptla[0],&Effla2[0],0,0);
    Effla[3] = new TGraphErrors(18,&Ptla[0],&Effla3[0],0,0);

    for(int i=0; i<4; i++)
    {
        Effks[i]->SetMarkerColor(kRed);
        Effks[i]->SetMarkerStyle(20);
        Effks[i]->SetMarkerSize(1.5);
        Effks[i]->SetLineColor(kRed);

        Effla[i]->SetMarkerColor(kBlue-4);
        Effla[i]->SetMarkerStyle(22);
        Effla[i]->SetMarkerSize(1.5);
        Effla[i]->SetLineColor(kBlue-4);

        //can[i] = new TCanvas(Form("Eff_%d",i),Form("Rapidity Bin %d",i+1),800,800);
        c1->cd(i+1);
        frame[i] = c1->DrawFrame(0,0,10.5,0.4);
        gPad->SetTickx();
        gPad->SetTicky();
        frame[i]->GetXaxis()->CenterTitle(1);
        frame[i]->GetYaxis()->CenterTitle(1);
        frame[i]->GetXaxis()->SetTitleSize(0.05);
        frame[i]->GetXaxis()->SetTitle("p_{T} (GeV)");
        frame[i]->GetYaxis()->SetTitle("Efficiency");
        frame[i]->GetYaxis()->SetTitleSize(0.05);
        frame[i]->SetTitleOffset(1.0,"Y");
        frame[i]->SetTitleOffset(1.0,"X");

        Effks[i]->Draw("P");
        Effla[i]->Draw("P");

        TLegend* leg = new TLegend(0.15,0.55,0.27,0.75);
        leg->SetFillColor(10);
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.05);
        //leg->AddEntry(ha_v2, "h#kern[-0.3]{#lower[0.2]{{}^{#pm}}}", "P");
        leg->AddEntry(Effks[i], "K_{S}^{0}", "P");
        leg->AddEntry(Effla[i], "#Lambda", "P");
        leg->Draw();
    }
    frame[0]->SetTitle("-1.0 < y < -0.5");
    frame[1]->SetTitle("-0.5 < y < -0.1");
    frame[2]->SetTitle("-0.1 < y < 0.3");
    frame[3]->SetTitle("0.3 < y < 1.0");
    c1->Print("EffPlot.pdf");
}
