#include <TLatex.h>

void fitTut()
{
   TH1F *background;

   gStyle->SetOptStat(0); //To get rid of stat box
      gStyle->SetOptFit(0);
   //gStyle->setOptFit(1111); //For chi-squared and probability

   TCanvas *c1 = new TCanvas("c1","Lifetime of a Muon",10,10,700,500);
   c1->SetFillColor(33);
   c1->SetFrameFillColor(41);
   c1->SetGrid();

   Double_t fitfunc(Double_t *x, Double_t *par)
   {
      Double_t xx = x[0];
      Double_t br = -par[2]; //background
      Double_t sr = par[0]*TMath::Exp(par[1]*xx); //signal
      return sr + br;
   }

   TFile *f = new TFile("muonDecayDataSameBinW.root");
   TH1F *result = (TH1F*)f->Get("muonSpectra");
   
   TF1 *ftot = new TF1("ftot",fitfunc,5,10.4,3);
   ftot->SetNpx(30);
   ftot->SetParameters(4,-0.4,2);
   ftot->SetParNames("Normalization","decay","background");
   
   result->Fit("ftot","","",1.4,10.4);
   result->GetXaxis()->SetTitle("Time (#mus)");
   result->GetYaxis()->SetTitle("Decay Count");
   
}

 /* 

Fitting with a user defined exponential
 https://root.cern.ch/phpBB3/viewtopic.php?t=12210


 */
