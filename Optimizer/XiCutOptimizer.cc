//Includes
#include <TStyle.h>
#include "TF1.h"
#include "TVirtualFitter.h"
#include "TMathText.h"
#include "TImage.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TTree.h"
#include <TString.h>
#include "TStyle.h"
#include "TString.h"
#include "TGaxis.h"
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
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "TString.h"
#include "TGaxis.h"
//#include "TROOT.h"

//#include <vector>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <sys/stat.h> //For file existance checking

bool CheckValue(ROOT::Internal::TTreeReaderValueBase& value) {
   if (value.GetSetupStatus() < 0) {
      std::cerr << "Error " << value.GetSetupStatus()
                << "setting up reader for " << value.GetBranchName() << '\n';
      return false;
   }
   return true;
}

bool XiCutOptimizer(std::string name)
{
    //Initializers
    bool Cut = false;
    TH1::SetDefaultSumw2();

    //Cut parameters to be varied. The commented elements are Hong's cuts to make a plot of them remember to change the numparam value as well
    int numCuts = 6;
    std::vector<float> om_om3dipsig     = {2.6 , 2.8 , 3.0 , 3.2 , 3.3  , 3.4 , 3.5};//3.0 };//  , 2.5};
    std::vector<float> om_omKaon3dipsig = {3.5 , 3.6 , 3.7 , 3.8 , 4.0  , 4.2 , 4.4};//4.0 };// , 5.0};
    std::vector<float> om_vtrkpi3dipsig = {2.5 , 2.6 , 2.7 , 2.8 , 3.0  , 3.2 , 3.4};//3.0 };//  , 4.0};
    std::vector<float> om_vtrkp3dipsig  = {1.5 , 1.6 , 1.7 , 1.8 , 2.0  , 2.2 , 2.4};//2.0 };//  , 3.0};
    std::vector<float> om_omflightsig   = {1.5 , 1.6 , 1.7 , 1.8 , 2.0  , 2.2 , 2.4};//2.0 };//  , 3.0};
    std::vector<float> om_distancesig   = {5.0 , 6.0 , 7.0 , 8.0 , 10.0 , 12. , 14.};//10.0};//  , 12.0};

    double misIDMass = 0.015;
    double rapidity = 1.0;
    double etacut = 2.4;
    int multHigh_ = 250;
    int numparam = om_om3dipsig.size();
    const int nDim = 9;
    // mass, pt, eta, then cuts
    //These are for eta cut
    //std::vector<int> nBins = {150,9,10,static_cast<int>(om_om3dipsig.size()),static_cast<int>(om_omKaon3dipsig.size()),static_cast<int>(om_vtrkpi3dipsig.size()),static_cast<int>(om_vtrkp3dipsig.size()),static_cast<int>(om_omflightsig.size()),static_cast<int>(om_distancesig.size())};
    std::vector<int> nBins = {150,9,10,static_cast<int>(om_om3dipsig.size()-1),static_cast<int>(om_omKaon3dipsig.size()-1),static_cast<int>(om_vtrkpi3dipsig.size()-1),static_cast<int>(om_vtrkp3dipsig.size()-1),static_cast<int>(om_omflightsig.size()-1),static_cast<int>(om_distancesig.size()-1)};
    //int nBins[nDim] = {150,9,10,(int)om_om3dipsig.size()+1,(int)om_omKaon3dipsig.size()+1,(int)om_vtrkpi3dipsig.size()+1,(int)om_vtrkp3dipsig.size()+1,(int)om_omflightsig.size()+1,(int)om_distancesig.size()+1};
    //int nBins[nDim] = {150,9,10,7,7,7,7,7,7};
    //double minbins[nDim] = {1.6  , 1.0 , -2.4 , 0                   , 3.5   , 2.5   , 1.5   , 1.5   , 5.0};
    double minbins[nDim] = {1.6  , 1.0 , -2.4 , 0                   , 3.5   , 2.5   , 1.5   , 1.5   , 5.0};
    double maxbins[nDim] = {1.75 , 10  , 2.4  , om_om3dipsig.back() , om_omKaon3dipsig.back() , om_vtrkpi3dipsig.back() , om_vtrkp3dipsig.back() , om_omflightsig.back() , om_distancesig.back()};
    //double maxbins[nDim] = {1.75 , 10  , 2.4  , 3.5 , 99999 , 99999 , 99999 , 99999 , 99999};

//These are for rap cut
    //int nBins[nDim] = {150,9,12,om_om3dipsig.size()+1,om_omKaon3dipsig.size()+1,om_vtrkpi3dipsig.size()+1,om_vtrkp3dipsig.size()+1,om_omflightsig.size()+1,om_distancesig.size()+1};
//These are for eta cut
    //float minbins[nDim] = {1.6  , 1.0 , -5 , 0                   , 3.5   , 2.5   , 1.5   , 1.5   , 5.0};
    //float maxbins[nDim] = {1.75 , 10  , 5  , om_om3dipsig.back() , 99999 , 99999 , 99999 , 99999 , 99999};


    //Now make variable sized bins for everything besides mass
    std::vector<double> VvarBins_pt            = {1.0,1.5,1.8,2.2,2.8,3.6,4.6,6.0,7.2,10.0};
    std::vector<double> VvarBins_eta           = {-2.4,-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,2.0,2.4};
    //std::vector<double> VvarBins_om3dipsig     = {0   , 2.59999 , 2.79999 , 2.99999 , 3.199999  , 3.299999 , 3.399999 , 3.499999};
    //std::vector<double> VvarBins_omKaon3dipsig = {3.5 , 3.60001 , 3.70001 , 3.80001 , 4.000001  , 4.200001 , 4.400001 , 99999};
    //std::vector<double> VvarBins_vtrkpi3dipsig = {2.5 , 2.60001 , 2.70001 , 2.80001 , 3.000001  , 3.200001 , 3.400001 , 99999};
    //std::vector<double> VvarBins_vtrkp3dipsig  = {1.5 , 1.60001 , 1.70001 , 1.80001 , 2.000001  , 2.200001 , 2.400001 , 99999};
    //std::vector<double> VvarBins_omflightsig   = {1.5 , 1.60001 , 1.70001 , 1.80001 , 2.000001  , 2.200001 , 2.400001 , 99999};
    //std::vector<double> VvarBins_distancesig   = {5.0 , 6.00001 , 7.00001 , 8.00001 , 10.00001 , 12.00001 , 14.00001 , 99999};

    std::vector<double> VvarBins_om3dipsig     = {2.6 , 2.8 , 3.0 , 3.2  , 3.3 , 3.4 , 3.5};
    std::vector<double> VvarBins_omKaon3dipsig = {3.5 , 3.6 , 3.7 , 3.8 , 4.0  , 4.2 , 4.4};
    std::vector<double> VvarBins_vtrkpi3dipsig = {2.5 , 2.6 , 2.7 , 2.8 , 3.0  , 3.2 , 3.4};
    std::vector<double> VvarBins_vtrkp3dipsig  = {1.5 , 1.6 , 1.7 , 1.8 , 2.0  , 2.2 , 2.4};
    std::vector<double> VvarBins_omflightsig   = {1.5 , 1.6 , 1.7 , 1.8 , 2.0  , 2.2 , 2.4};
    std::vector<double> VvarBins_distancesig   = {5.0 , 6.0 , 7.0 , 8.0 , 10.0 , 12. , 14.};



    //for(int i=0; i<nBins[3]; i++){
        //if(i==0) VvarBins_om3dipsig.push_back(0);
        //else{
            //VvarBins_om3dipsig.push_back(om_om3dipsig[i-1]);
        //}
    //}
    //for(int i=0; i<nBins[4]; i++){
        //if(i==nBins[4]-1) VvarBins_omKaon3dipsig.push_back(99999);
        //else{
            //VvarBins_omKaon3dipsig.push_back(om_omKaon3dipsig[i]);
        //}
    //}
    //for(int i=0; i<nBins[5]; i++){
        //if(i==nBins[5]-1) VvarBins_vtrkpi3dipsig.push_back(99999);
        //else{
            //VvarBins_vtrkpi3dipsig.push_back(om_vtrkpi3dipsig[i]);
        //}
    //}
    //for(int i=0; i<nBins[6]; i++){
        //if(i==nBins[6]-1) VvarBins_vtrkp3dipsig.push_back(99999);
        //else{
            //VvarBins_vtrkp3dipsig.push_back(om_vtrkp3dipsig[i]);
        //}
    //}
    //for(int i=0; i<nBins[7]; i++){
        //if(i==nBins[7]-1) VvarBins_omflightsig.push_back(99999);
        //else{
            //VvarBins_omflightsig.push_back(om_omflightsig[i]);
        //}
    //}
    //for(int i=0; i<nBins[8]; i++){
        //if(i==nBins[8]-1) VvarBins_distancesig.push_back(99999);
        //else{
            //VvarBins_distancesig.push_back(om_distancesig[i]);
        //}
    //}

    THnSparseF* h = new THnSparseF("hSparse","Cuts",nDim,&nBins[0],minbins,maxbins);

    cout << "1" << endl;
    h->GetAxis(1)->Set(nBins[1],&VvarBins_pt[0]);
    cout << "2" << endl;
    h->GetAxis(2)->Set(nBins[2],&VvarBins_eta[0]);
    cout << "3" << endl;
    h->GetAxis(3)->Set(nBins[3],&VvarBins_om3dipsig[0]);
    cout << "4" << endl;
    h->GetAxis(4)->Set(nBins[4],&VvarBins_omKaon3dipsig[0]);
    cout << "5" << endl;
    h->GetAxis(5)->Set(nBins[5],&VvarBins_vtrkpi3dipsig[0]);
    cout << "6" << endl;
    h->GetAxis(6)->Set(nBins[6],&VvarBins_vtrkp3dipsig[0]);
    cout << "7" << endl;
    h->GetAxis(7)->Set(nBins[7],&VvarBins_omflightsig[0]);
    cout << "8" << endl;
    h->GetAxis(8)->Set(nBins[8],&VvarBins_distancesig[0]);



    //Containers

    //Hist Containers
    //TH2D* hom_om3dipsig     [numparam];
    //TH2D* hom_omKaon3dipsig   [numparam];
    //TH2D* hom_vtrkpi3dipsig [numparam];
    //TH2D* hom_vtrkp3dipsig  [numparam];
    //TH2D* hom_omflightsig   [numparam];
    //TH2D* hom_distancesig   [numparam];

    //Tree setup
    //TFile* f1=TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/V0CasTree_09_14_17.root");
    TFile* f1=TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/V0CasTreePbPb_09_25_17.root");
    //TFile* f1= new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/PbPbTreeSmall.root");

    TTreeReader reader("OmTreeProducerRapidityPbPb/OmTree",f1);

    TTreeReaderArray<float> om3dipsig(reader    ,"om3dipsig.om3dipsig");
    TTreeReaderArray<float> omKaon3dipsig(reader,"omKaon3dipsig.omKaon3dipsig");
    TTreeReaderArray<float> vtrkpi3dipsig(reader,"vtrkpi3dipsig.vtrkpi3dipsig");
    TTreeReaderArray<float> vtrkp3dipsig(reader ,"vtrkp3dipsig.vtrkp3dipsigpt");
    TTreeReaderArray<float> omflightsig(reader  ,"omflightsig.omflightsig");
    TTreeReaderArray<float> distancesig(reader  ,"distancesig.distancesig");
    TTreeReaderArray<float> mass(reader         ,"mass.mass");
    TTreeReaderArray<float> pt(reader           ,"pt.pt");
    TTreeReaderArray<float> eta(reader          ,"eta.eta");
    TTreeReaderArray<float> rap(reader          ,"rapidity.rapidity");
    //TTreeReaderArray<int> nTrkAcc(reader        ,"nTrkAcc.nTrkAcc");
    //TTreeReaderArray<float> misIDMassLapi(reader,"misIDMassLapi.misIDMassLapi");
    //TTreeReaderArray<float> misIDMasspiLa(reader,"misIDMasspiLa.misIDMasspiLa");

    //Intialize Histograms
    TH2D* hom_NoCut = NULL;
    if(!Cut) hom_NoCut = new TH2D("hom_NoCut"    ,"NoCut"  ,150,1.25,1.40,150,0,15);
    TH2D* hom_defaultcut = new TH2D("hom_Default","Default",150,1.60,1.75,400,0,40);
    //for(int j=0; j<numparam; j++)
    //{
        //hom_om3dipsig[j]     = new TH2D(Form("hom_om3dipsig_%.1f"    ,om_om3dipsig[j])    ,Form("hom_om3dipsig_%.1f"    ,om_om3dipsig[j])    ,150,1.25,1.40,150,0,15);
        //hom_omKaon3dipsig[j] = new TH2D(Form("hom_omKaon3dipsig_%.1f",om_omKaon3dipsig[j]),Form("hom_omKaon3dipsig_%.1f",om_omKaon3dipsig[j]),150,1.25,1.40,150,0,15);
        //hom_vtrkpi3dipsig[j] = new TH2D(Form("hom_vtrkpi3dipsig_%.1f",om_vtrkpi3dipsig[j]),Form("hom_vtrkpi3dipsig_%.1f",om_vtrkpi3dipsig[j]),150,1.25,1.40,150,0,15);
        //hom_vtrkp3dipsig[j]  = new TH2D(Form("hom_vtrkp3dipsig_%.1f" ,om_vtrkp3dipsig[j]) ,Form("hom_vtrkp3dipsig_%.1f" ,om_vtrkp3dipsig[j]) ,150,1.25,1.40,150,0,15);
        //hom_omflightsig[j]   = new TH2D(Form("hom_omflightsig_%.1f"  ,om_omflightsig[j])  ,Form("hom_omflightsig_%.1f"  ,om_omflightsig[j])  ,150,1.25,1.40,150,0,15);
        //hom_distancesig[j]   = new TH2D(Form("hom_distancesig_%.1f"  ,om_distancesig[j])  ,Form("hom_distancesig_%.1f"  ,om_distancesig[j])  ,150,1.25,1.40,150,0,15);
    //}



    int j=0;
    while(reader.Next())
    //while(j<100)
    {
        if(!CheckValue(om3dipsig))     return false;
        if(!CheckValue(omKaon3dipsig)) return false;
        if(!CheckValue(vtrkpi3dipsig)) return false;
        if(!CheckValue(vtrkp3dipsig))  return false;
        if(!CheckValue(omflightsig))   return false;
        if(!CheckValue(distancesig))   return false;
        if(!CheckValue(mass))          return false;
        if(!CheckValue(pt))            return false;
        if(!CheckValue(eta))           return false;

        for(int i=0; i<mass.GetSize(); i++)
        {
            double value[nDim] = {mass[i],pt[i],eta[i],om3dipsig[i],omKaon3dipsig[i],vtrkpi3dipsig[i],vtrkp3dipsig[i],omflightsig[i],distancesig[i]};
            h->Fill(value);
        }
        j++;
    }
f1->Close();
    TCanvas* c1 = new TCanvas("c1","",600,600);
    TH1D* h2 = new TH1D("1","",100,0,10);
    h2->Fill(1);
    h2->Draw();

    //Output file creation
    std::string filetype = ".root";
    name += filetype;
    //struct stat buffer;
    //if(stat(name.c_str(), &buffer) == 0)
    //{
        //cout << "File with this name already eomsts, please select a different name" << endl;
        //return false;
    //}
    TFile out(name.c_str(),"RECREATE");
    cout << name.c_str() << " created!" << endl;

    //Write histograms to root file
    if(Cut)
    {
        //for(int j=0; j<numparam; j++)
            //hom_om3dipsig[j]->Write();
        //for(int j=0; j<numparam; j++)
            //hom_omKaon3dipsig[j]->Write();
        //for(int j=0; j<numparam; j++)
            //hom_vtrkpi3dipsig[j]->Write();
        //for(int j=0; j<numparam; j++)
            //hom_vtrkp3dipsig[j]->Write();
        //for(int j=0; j<numparam; j++)
            //hom_omflightsig[j]->Write();
        //for(int j=0; j<numparam; j++)
            //hom_distancesig[j]->Write();
    }
    else
    {
        //hom_defaultcut->Write();
        h->Write();
        out.Close();
        //c2->cd(3);
        //h->Projection(4)->Draw();
        //c2->cd(4);
        //h->Projection(3)->Draw();
    }
    return true;
}

void readSparse()
{
    std::vector<double> VvarBins_pt            = {1.0,1.5,1.8,2.2,2.8,3.6,4.6,6.0,7.2,10.0};
    std::vector<double> VvarBins_eta           = {-2.4,-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,2.0,2.4};

    std::vector<double> VvarBins_om3dipsig     = {0   , 2.60000 , 2.80000 , 3.00000 , 3.200000  , 3.300000 , 3.400000 , 3.500000};
    std::vector<double> VvarBins_omKaon3dipsig = {3.5 , 3.60000 , 3.70000 , 3.80000 , 4.000000  , 4.200000 , 4.400000 , 99999};
    std::vector<double> VvarBins_vtrkpi3dipsig = {2.5 , 2.60000 , 2.70000 , 2.80000 , 3.000000  , 3.200000 , 3.400000 , 99999};
    std::vector<double> VvarBins_vtrkp3dipsig  = {1.5 , 1.60000 , 1.70000 , 1.80000 , 2.000000  , 2.200000 , 2.400000 , 99999};
    std::vector<double> VvarBins_omflightsig   = {1.5 , 1.60000 , 1.70000 , 1.80000 , 2.000000  , 2.200000 , 2.400000 , 99999};
    std::vector<double> VvarBins_distancesig   = {5.0 , 6.00000 , 7.00000 , 8.00000 , 10.00000 , 12.00000 , 14.00000 , 99999};

    //std::vector<double> VvarBins_om3dipsig     = {2.6 , 2.8 , 3.0 , 3.2  , 3.3 , 3.4 , 3.5};
    //std::vector<double> VvarBins_omKaon3dipsig = {3.5 , 3.6 , 3.7 , 3.8 , 4.0  , 4.2 , 4.4};
    //std::vector<double> VvarBins_vtrkpi3dipsig = {2.5 , 2.6 , 2.7 , 2.8 , 3.0  , 3.2 , 3.4};
    //std::vector<double> VvarBins_vtrkp3dipsig  = {1.5 , 1.6 , 1.7 , 1.8 , 2.0  , 2.2 , 2.4};
    //std::vector<double> VvarBins_omflightsig   = {1.5 , 1.6 , 1.7 , 1.8 , 2.0  , 2.2 , 2.4};
    //std::vector<double> VvarBins_distancesig   = {5.0 , 6.0 , 7.0 , 8.0 , 10.0 , 12. , 14.};

    std::ostringstream os;
    std::ostringstream osYield;
    using namespace RooFit;
    TH1::SetDefaultSumw2();
    //TGaxis::SetMaxDigits(3);
    RooMsgService::instance().setStreamStatus(0,kFALSE);
    RooMsgService::instance().setStreamStatus(1,kFALSE);
    gStyle->SetOptStat(1111111);
    TFile* f = new TFile("OmegaFineBin.root");
    //TFile* f = new TFile("SmallBins.root");
    THnSparseF* h = (THnSparseF*)f->Get("hSparse");
    //c2->Divide(1,2);
    //c2->cd(1);
    //h1->Draw();
    //c2->cd(2);
    //h->Projection(1)->Draw();
    //for(int a=1; a<VvarBins_pt.size(); a++)
    TCanvas* c2 = new TCanvas("c2","",600,600);
    std::map<int,double> SigSig;
    for(int a=5; a<6; a++)
    {
        h->GetAxis(1)->SetRange(a,a);
        //for(int b=1; b<VvarBins_om3dipsig.size(); b++)
        for(int b=2; b<5; b++)
        {
            h->GetAxis(3)->SetRange(1,b);
            //for(int c=1; c<VvarBins_omKaon3dipsig.size(); c++)
            for(int c=3; c<7; c++)
            {
                h->GetAxis(4)->SetRange(c,7);
                //for(int d=1; d<VvarBins_vtrkpi3dipsig.size(); d++)
                for(int d=3; d<7; d++)
                {
                    h->GetAxis(5)->SetRange(d,7);
                    //for(int e=1; e<VvarBins_vtrkp3dipsig.size();e++)
                    for(int e=3; e<7;e++)
                    {
                        h->GetAxis(6)->SetRange(e,7);
                        //for(int f=1; f<VvarBins_omflightsig.size(); f++)
                        for(int f=3; f<7; f++)
                        {
                            h->GetAxis(7)->SetRange(f,7);
                            //for(int i=1; i<VvarBins_distancesig.size(); i++)
                            for(int i=3; i<7; i++)
                            {
                                struct stat buffer;
                                std::ostringstream name;
                                int mask = b*1e5 + c*1e4 + d*1e3 + e*1e2 + f*1e1 + i;
                                name << "OmegaDistributions/Omega_" << mask << ".pdf";
                                if(stat(name.str().c_str(),&buffer) == 0) continue;
                                h->GetAxis(8)->SetRange(i,7);
                                TH1D* h1 = (TH1D*)h->Projection(0);
                                h1->GetYaxis()->SetRangeUser(0,1000);
                                gStyle->SetOptTitle(kFALSE);

                                TLatex* tex = new TLatex();
                                tex->SetNDC();
                                tex->SetTextFont(42);
                                tex->SetTextSize(0.025);
                                //tex->SetTextAlign(10);

                                RooRealVar x("x","mass",1.60,1.75);
                                //RooRealVar x("x","mass",1.62,1.73);
                                RooPlot* xframe_ = x.frame(150);
                                xframe_->GetXaxis()->SetTitle("#Lambda K Invariant mass (GeV)");
                                xframe_->GetYaxis()->SetTitle("Candidates / 0.002 GeV");
                                xframe_->GetXaxis()->CenterTitle(1);
                                xframe_->GetYaxis()->CenterTitle(1);
                                xframe_->GetXaxis()->SetTickSize(0.02);
                                xframe_->GetYaxis()->SetTickSize(0.02);
                                //xframe_->GetXaxis()->SetNdivisions(407);
                                //xframe_->GetYaxis()->SetNdivisions(410);
                                //xframe_->GetXaxis()->SetTitleSize(0.06);
                                //xframe_->GetYaxis()->SetTitleSize(0.06);
                                //xframe_->GetYaxis()->SetTitleOffset(1.05);
                                //xframe_->GetXaxis()->SetTitleOffset(0.5);
                                //xframe_->GetXaxis()->SetLabelSize(xframe_->GetXaxis()->GetLabelSize()*2.0);
                                //xframe_->GetYaxis()->SetLabelSize(0.1);
                                //xframe_->GetXaxis()->SetLabelSize(0.1);
                                //xframe_->GetYaxis()->SetLabelSize(xframe_->GetYaxis()->GetLabelSize()*2.0);
                                RooDataHist data("data","dataset",x,h1);
                                data.plotOn(xframe_,Name("data"));
                                RooRealVar mean("mean","mean",1.67,1.6,1.75);
                                RooRealVar sigma1("sigma1","sigma1",0.004,0.001,0.04);
                                RooRealVar sigma2("sigma2","sigma2",0.005,0.001,0.04);
                                RooRealVar sig1("sig1","signal1",5000,0,10000);
                                RooRealVar sig2("sig2","signal2",5000,0,10000);
                                RooRealVar qsig("qsig","qsig",5000,0,100000);
                                RooRealVar alpha("alpha","alpha",0.9,-1,10);
                                RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
                                RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
                                //RooRealVar a("a","a",0,-100000,100000);
                                //RooRealVar b("b","b",0,-100000,100000);
                                //RooRealVar cp("cp","cp",0,-100000,100000);
                                //RooRealVar d("d","d",0,-100000,100000);
                                //RooPolynomial background("poly","poly",x,RooArgList(a,b,cp,d));
                                //RooPolynomial background("poly","poly",x,RooArgList(a,b,cp));
                                //RooRealVar qsig("polysig","polysig",10,0,1000000000);
                                RooGenericPdf background("background", "x - (1.115683 + 0.493677)^alpha", RooArgList(x,alpha));
                                RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,background),RooArgList(sig1,sig2,qsig));

                                x.setRange("cut",1.60,1.75);
                                //x.setRange("cut",1.62,1.73);

                                RooFitResult* r_xi = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));
                                //RooChi2Var chi2_xiVar("chi2_xi","chi2",sum,data);

                                double covQual = r_xi->covQual();
                                double mean_xi = mean.getVal();

                                double gaus1F_xi = sig1.getVal();
                                double gaus2F_xi = sig2.getVal();
                                double qsig_xi   = qsig.getVal();

                                //set ranges for individual gaussian yield determination
                                x.setRange("g1", mean.getVal() - 2*sigma1.getVal(), mean.getVal() + 2*sigma1.getVal());
                                x.setRange("g2", mean.getVal() - 2*sigma2.getVal(), mean.getVal() + 2*sigma2.getVal());

                                RooAbsReal* Intgaus1_yield_xi = gaus1.createIntegral(x,x,"g1");
                                RooAbsReal* Intgaus2_yield_xi = gaus2.createIntegral(x,x,"g2");

                                double gaus1_yield_xi = gaus1F_xi*Intgaus1_yield_xi->getVal();
                                double gaus2_yield_xi = gaus2F_xi*Intgaus2_yield_xi->getVal();
                                double gausTot_yield_xi = gaus1_yield_xi + gaus2_yield_xi;

                                cout << "Yield1: " << gaus1_yield_xi << endl;
                                cout << "Yield2: " << gaus2_yield_xi << endl;

                                double rms_gaus1_sig_xi = gaus1_yield_xi/gausTot_yield_xi;
                                double rms_gaus2_sig_xi = gaus2_yield_xi/gausTot_yield_xi;
                                double rms_true_xi = TMath::Sqrt(rms_gaus1_sig_xi*sigma1.getVal()*sigma1.getVal() + rms_gaus2_sig_xi*sigma2.getVal()*sigma2.getVal());

                                x.setRange("peak", mean.getVal() - 2*rms_true_xi, mean.getVal() + 2*rms_true_xi);
                                RooAbsReal* Intgaus1_xi      = gaus1.createIntegral(x, x,  "peak");
                                RooAbsReal* Intgaus2_xi      = gaus2.createIntegral(x, x, "peak");
                                RooAbsReal* Intbackground_xi = background.createIntegral(x, x, "peak");

                                double Intgaus1E_xi      = gaus1F_xi*Intgaus1_xi->getVal();
                                double Intgaus2E_xi      = gaus2F_xi*Intgaus2_xi->getVal();
                                double IntbackgroundE_xi = qsig_xi*Intbackground_xi->getVal();
                                double totsig_xi         = Intgaus1E_xi + Intgaus2E_xi + IntbackgroundE_xi;
                                double Yield_xi          = Intgaus1E_xi + Intgaus2E_xi;


                                double Fsig_xi = Yield_xi/totsig_xi;

                                double significance = Yield_xi/sqrt(totsig_xi);

                                cout << "Yield (xi): " << Yield_xi << endl;
                                cout << "Fsig (xi): " << Fsig_xi << endl;
                                cout << "std (xi): "  << rms_true_xi  << endl;
                                cout << "mass (xi): " << mean_xi << endl;

                                cout << "covQual (xi)" << covQual << endl;
                                cout << "Signal Sig (xi)" << significance << endl;


                                sum.plotOn(xframe_,Name("sum"),NormRange("cut"),LineWidth(1),LineColor(kBlue));
                                sum.plotOn(xframe_,Components(background),NormRange("cut"),LineStyle(kDashed),LineWidth(1),LineColor(kBlue));
                                c2->cd();
                                gPad->SetBottomMargin(0.15); //gives more space for titles
                                gPad->SetLeftMargin(0.15);
                                gPad->SetTickx(  );
                                gPad->SetTicky(  );
                                xframe_->Draw();
                                c2->Update();
                                double chi2_xi = xframe_->chiSquare("sum","data",4);

                                TLine* t1 = new TLine(mean.getVal() - 2*rms_true_xi, 0, mean.getVal() - 2*rms_true_xi, gPad->GetUymax());
                                TLine* t2 = new TLine(mean.getVal() + 2*rms_true_xi, 0, mean.getVal() + 2*rms_true_xi, gPad->GetUymax());
                                t1->SetLineStyle(2);
                                t1->SetLineColor(kGreen);
                                t2->SetLineStyle(2);
                                t2->SetLineColor(kGreen);
                                t1->Draw("same");
                                t2->Draw("same");

                                double xpos = 0.7;
                                double ypos = 0.85;
                                double increment = 0.07;
                                os << "Mean: " << std::setprecision(5) << mean_xi << " GeV" << std::setprecision(6);
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());
                                os << "#sigma :" << std::setprecision(2) << rms_true_xi << " GeV" << std::setprecision(6);
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());
                                os << "CovQual: " << covQual;
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());
                                //os << "#chi^{2}/ndf: " << std::setprecision(3) << chi2_xi << std::setprecision(6);
                                //tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                //os.str(std::string());
                                osYield << "Yield: " << std::setprecision(2) << Yield_xi;
                                tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
                                osYield.str(std::string());
                                os.str(std::string());
                                os << "Sig: " << significance;
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());

                                xpos = 0.20;
                                ypos = 0.85;
                                os << "#Omega DCA < " << VvarBins_om3dipsig[b];
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());
                                os << "Kaon DCA > " << VvarBins_omKaon3dipsig[c-1];
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());
                                os << "#pi DCA > " << VvarBins_vtrkpi3dipsig[d-1];
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());
                                os << "proton DCA > " << VvarBins_vtrkp3dipsig[e-1];
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());
                                os << "#Omega DecayL > " << VvarBins_omflightsig[f-1];
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());
                                os << "#Lambda DecayL > " << VvarBins_distancesig[i-1];
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());

                                SigSig[mask] = significance;

                                c2->Print(Form("OmegaDistributions/Omega_%d.pdf",mask));
                            }
                        }
                    }
                }
            }
        }
    }
    double value = 0;
    int pairvalue = 0;
    double maxVal = 0;
    int maxpair = 0;
    for(std::map<int,double>::iterator it=SigSig.begin(); it!=SigSig.end(); ++it){
        value = it->second;
        if(value > maxVal){ 
            maxVal = value;
            maxpair = it->first;
        }
    }
    cout << "Max Sig: " << maxVal << endl;
    cout << "Max pair: " << maxpair << endl;
}

/*
                for(int i=0;i<mass.GetSize();i++)
                {
                    //if(nTrkAcc[i]                 > multHigh_)           continue;
                    if(om3dipsig[i]               > om_om3dipsig[0])     continue;
                    //if(std::fabs(eta[i])          > etacut)              continue;
                    if(std::fabs(rap[i])          > rapidity)              continue;
                    if(omKaon3dipsig[i]           < om_omKaon3dipsig[0]) continue;
                    if(vtrkpi3dipsig[i]           < om_vtrkpi3dipsig[0]) continue;
                    if(vtrkp3dipsig[i]            < om_vtrkp3dipsig[0])  continue;
                    if(omflightsig[i]             < om_omflightsig[0])   continue;
                    if(distancesig[i]             < om_distancesig[0])   continue;
                    //if(std::abs(misIDMasspiLa[i]) < misIDMass)           continue;
                    //if(std::abs(misIDMassLapi[i]) < misIDMass)           continue;

                    hom_defaultcut->Fill(mass[i],pt[i]);
                }
        //}
        */
