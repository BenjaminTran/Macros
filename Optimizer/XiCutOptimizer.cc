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
#include "RooBifurGauss.h"
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

void WriteVector(std::vector<float> cutVector, std::ofstream &myfile, std::string name){
    myfile << name << ": ";
    for(std::vector<float>::iterator it=cutVector.begin(); it!=cutVector.end(); it++){
        myfile << *it << ", ";
    }
    myfile << "\n";
}

RooRealVar MakeMeanVar(std::string name1, std::string name2, double init, double minRange, double maxRange){
    RooRealVar Var1(name1.c_str(),name2.c_str(),init,minRange,maxRange);
    return Var1;
}

TTreeReaderValue<float> SetValues(TTreeReader& reader, std::string branchName ){
    TTreeReaderValue<float> ReaderValue(reader, branchName.c_str());
    return ReaderValue;
}

bool SparseCreator(std::string name = "Test", std::string PID = "Omega")
{
    //Initializers
    TH1::SetDefaultSumw2();
    std::ofstream myfile;
    std::string txtname = name + ".txt";
    myfile.open(txtname.c_str());

    //Cut parameters to be varied. The commented elements are Hong's cuts to make a plot of them remember to change the numparam value as well
    std::vector<float> VvarBins_eta           = {-2.4,-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,2.0,2.4};

    std::vector<float> VvarBins_pt;
    std::vector<float> Vcas3dipsig;
    std::vector<float> VBat3dipsig;
    std::vector<float> Vvtrkpi3dipsig;
    std::vector<float> Vvtrkp3dipsig;
    std::vector<float> Vflightsig;
    std::vector<float> Vdistancesig;

    std::vector<float> vdz1;
    std::vector<float> vdz2;
    std::vector<float> vdxy1;
    std::vector<float> vdxy2;
    std::vector<float> vdls;
    std::vector<float> vagl;

    if(PID == "Omega")
    {
        VvarBins_pt.insert(VvarBins_pt.begin(),       {1.0,1.5,1.8,2.2,2.8,3.6,4.6,6.0,7.2,10.0});

        //v2 && v4 w/ corrected tree
        Vcas3dipsig.insert(Vcas3dipsig.begin(),       {0.5 , 1.0  , 1.5  , 1.8 , 2.0 , 2.5 , 2.8 , 3.0 , 3.5 });
        VBat3dipsig.insert(VBat3dipsig.begin(),       {3.0 , 3.5  , 4.0  , 5.0 , 6.0 , 7.0 , 8.0 , 9.0 , 10.0});
        Vvtrkpi3dipsig.insert(Vvtrkpi3dipsig.begin(), {3.0 , 3.5  , 4.0 , 4.5  , 5.0 , 6.0 , 7.0 , 9.0 , 10.0});
        Vvtrkp3dipsig.insert(Vvtrkp3dipsig.begin(),   {2.0 , 2.5  , 3.0 , 3.5  , 4.0 , 5.0 , 7.0 , 9.5 , 10.0});
        Vflightsig.insert(Vflightsig.begin(),         {2.0 , 2.5  , 3.0 , 3.5  , 4.0 , 5.0 , 7.0 , 9.5 , 10.0});
        Vdistancesig.insert(Vdistancesig.begin(),     {10. , 10.5 , 11  , 11.5 , 12. , 13  , 16  , 18  , 20  });

        //v3
        //Vcas3dipsig.insert(Vcas3dipsig.begin(),       {1.5 , 1.6, 1.7 , 1.8 , 1.9 , 2.0});
        //VBat3dipsig.insert(VBat3dipsig.begin(),       {3.5 , 3.6, 3.8 , 4.0 , 4.2 , 4.4});
        //Vvtrkpi3dipsig.insert(Vvtrkpi3dipsig.begin(), {3.0 , 3.1 , 3.2 , 3.3 , 3.4 , 3.5});
        //Vvtrkp3dipsig.insert(Vvtrkp3dipsig.begin(),   {2.0 , 2.1 , 2.2 , 2.3 , 2.4 , 2.5});
        //Vflightsig.insert(Vflightsig.begin(),         {2.0 , 2.1 , 2.2 , 2.3 , 2.4 , 2.5});
        //Vdistancesig.insert(Vdistancesig.begin(),     {10. , 10.5 , 11  , 11.5 , 12. , 13});
    }
    else if(PID == "Xi")
    {
        VvarBins_pt.insert(VvarBins_pt.begin(),       {1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.2,10.0});

        //v2
        //Vcas3dipsig.insert(Vcas3dipsig.begin(),       {0.5 , 1.0 , 1.5 , 2.0 , 2.25, 2.5 , 3.0 , 3.5 , 4.0});
        //VBat3dipsig.insert(VBat3dipsig.begin(),       {4.5 , 5.0 , 5.5 , 6.0 , 6.5 , 7.0 , 7.5 , 8.0 , 9.0});
        //Vvtrkpi3dipsig.insert(Vvtrkpi3dipsig.begin(), {3.5 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 6.5 , 7.0 , 8.0});
        //Vvtrkp3dipsig.insert(Vvtrkp3dipsig.begin(),   {2.5 , 3.0 , 3.5 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 7.0});
        //Vflightsig.insert(Vflightsig.begin(),         {2.5 , 3.0 , 3.5 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 7.0});
        //Vdistancesig.insert(Vdistancesig.begin(),     {10. , 11. , 12. , 13. , 13.5 , 14 , 14.5 , 15 , 16 });

        //v3
        //Vcas3dipsig.insert(Vcas3dipsig.begin(),       {1.5 , 2.25, 2.5 , 2.8 , 3.5 , 4.0});
        //VBat3dipsig.insert(VBat3dipsig.begin(),       {4.8 , 5.0 , 5.2 , 6.0 , 6.5 , 7.0 , 7.5 , 8.0 , 9.0});
        //Vvtrkpi3dipsig.insert(Vvtrkpi3dipsig.begin(), {3.8 , 4.0 , 4.2 , 5.0 , 5.5 , 6.0 , 6.5 , 7.0 , 8.0});
        //Vvtrkp3dipsig.insert(Vvtrkp3dipsig.begin(),   {2.8 , 3.0 , 3.2 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 7.0});
        //Vflightsig.insert(Vflightsig.begin(),         {2.5 , 3.0 , 3.5 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 7.0});
        //Vdistancesig.insert(Vdistancesig.begin(),     {10. , 11.5 , 12. , 12.5 , 13.5 , 14 , 14.5 , 15 , 16 });

        //v4
        Vcas3dipsig.insert(Vcas3dipsig.begin(),       {1.0 ,1.5 , 2.25, 2.5 , 2.8 , 3.5 , 4.0});
        VBat3dipsig.insert(VBat3dipsig.begin(),       {4.8 , 5.0 , 5.2 , 6.0 , 6.5 , 7.0 , 7.5 , 8.0 , 9.0});
        Vvtrkpi3dipsig.insert(Vvtrkpi3dipsig.begin(), {3.8 , 4.0 , 4.2 , 5.0 , 5.5 , 6.0 , 6.5 , 7.0 , 8.0});
        Vvtrkp3dipsig.insert(Vvtrkp3dipsig.begin(),   {2.8 , 3.0 , 3.2 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 7.0});
        Vflightsig.insert(Vflightsig.begin(),         {2.5 , 3.0 , 3.5 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 7.0});
        Vdistancesig.insert(Vdistancesig.begin(),     {10. , 11.5 , 12. , 12.5 , 13.5 , 14 , 14.5 , 15 , 16 });
    }
    else if (PID == "Kshort")
    {
        VvarBins_pt.insert(VvarBins_pt.begin(),       {0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0});
        vdz1.insert(vdz1.begin()  , {1.0    , 1.2   , 1.5   , 1.7   , 2.0    , 2.5    , 3.0});
        vdz2.insert(vdz2.begin()  , {1.0    , 1.2   , 1.5   , 1.7   , 2.0    , 2.5    , 3.0});
        vdxy1.insert(vdxy1.begin() , {1.0    , 1.2   , 1.5   , 1.7   , 2.0    , 2.5    , 3.0});
        vdxy2.insert(vdxy2.begin() , {1.0    , 1.2   , 1.5   , 1.7   , 2.0    , 2.5    , 3.0});
        vdls.insert(vdls.begin()  , {4      , 5     , 6     , 7     , 8      , 9      , 10});
        vagl.insert(vagl.begin()  , {0.9965 , 0.997 , 0.998 , 0.999 , 0.9995 , 0.9999 , 0.99999});
    }
    else if(PID == "Lambda")
    {
        VvarBins_pt.insert(VvarBins_pt.begin(),       {0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0});
        vdz1.insert(vdz1.begin()   , {1.0 , 1.2 , 1.3 , 1.4 , 1.5});
        vdz2.insert(vdz2.begin()   , {1.0 , 1.2 , 1.3 , 1.4 , 1.5});
        vdxy1.insert(vdxy1.begin() , {1.0 , 1.2 , 1.3 , 1.4 , 1.5});
        vdxy2.insert(vdxy2.begin() , {1.0 , 1.2 , 1.3 , 1.4 , 1.5});
        vdls.insert(vdls.begin()   , {4.5, 4.8 , 5, 5.2, 5.4});
        vagl.insert(vagl.begin()   , {0.9996, 0.9997, 0.9998, 0.9999, 0.99992});
    }

    double misIDMass = 0.015;
    double rapidity = 1.0;
    double etacut = 2.4;
    int multHigh_ = 250;
    int numparam = Vcas3dipsig.size();
    const int nDim = 9;

    std::vector<int> nBins;
    std::vector<double> minbins;
    std::vector<double> maxbins;
    if(PID == "Omega" || PID == "Xi")
    {
        nBins.insert(nBins.begin(),{150,9,10,static_cast<int>(Vcas3dipsig.size()-1),static_cast<int>(VBat3dipsig.size()-1),static_cast<int>(Vvtrkpi3dipsig.size()-1),static_cast<int>(Vvtrkp3dipsig.size()-1),static_cast<int>(Vflightsig.size()-1),static_cast<int>(Vdistancesig.size()-1)});
    }
    else if(PID == "Kshort")
    {
        nBins.insert(nBins.begin(),{270,14,10,static_cast<int>(vdz1.size()-1),static_cast<int>(vdz2.size()-1),static_cast<int>(vdxy1.size()-1),static_cast<int>(vdxy2.size()-1),static_cast<int>(vdls.size()-1),static_cast<int>(vagl.size()-1)}); //14 pt bins because we are including up to 10 geV here but we are not in the analysis
    }
    else if(PID == "Lambda")
    {
        nBins.insert(nBins.begin(),{160,11,10,static_cast<int>(vdz1.size()-1),static_cast<int>(vdz2.size()-1),static_cast<int>(vdxy1.size()-1),static_cast<int>(vdxy2.size()-1),static_cast<int>(vdls.size()-1),static_cast<int>(vagl.size()-1)}); //14 pt bins because we are including up to 10 geV here but we are not in the analysis
    }

    //Omega
    if(PID == "Omega")
    {
        minbins.insert(minbins.begin(),{1.60,1.0,-2.4,0,VBat3dipsig.front(),Vvtrkpi3dipsig.front(),Vvtrkp3dipsig.front(),Vflightsig.front(),Vdistancesig.front()});
        maxbins.insert(maxbins.begin(),{1.75,10 ,2.4 ,Vcas3dipsig.back(),VBat3dipsig.back(),Vvtrkpi3dipsig.back(),Vvtrkp3dipsig.back(),Vflightsig.back(),Vdistancesig.back()});
    }
    else if(PID == "Xi")
    {
        minbins.insert(minbins.begin(),{1.25,1.0,-2.4,0,VBat3dipsig.front(),Vvtrkpi3dipsig.front(),Vvtrkp3dipsig.front(),Vflightsig.front(),Vdistancesig.front()});
        maxbins.insert(maxbins.begin(),{1.40,10,2.4,Vcas3dipsig.back(),VBat3dipsig.back(),Vvtrkpi3dipsig.back(),Vvtrkp3dipsig.back(),Vflightsig.back(),Vdistancesig.back()});
    }
    else if(PID == "Kshort")
    {
        minbins.insert(minbins.begin(),{0.43,0.2,-2.4,vdz1.front(),vdz2.front(),vdxy1.front(),vdxy2.front(),vdls.front(),vagl.front()});
        maxbins.insert(maxbins.begin(),{0.565,10,2.4,vdz1.back(),vdz2.back(),vdxy1.back(),vdxy2.back(),vdls.back(),vagl.back()});
    }
    else if(PID == "Lambda")
    {
        minbins.insert(minbins.begin(),{1.08,0.8,-2.4,vdz1.front(),vdz2.front(),vdxy1.front(),vdxy2.front(),vdls.front(),vagl.front()});
        maxbins.insert(maxbins.begin(),{1.16,10,2.4,vdz1.back(),vdz2.back(),vdxy1.back(),vdxy2.back(),vdls.back(),vagl.back()});
    }

    //Xi

//These are for rap cut
    //int nBins[nDim] = {150,9,12,Vcas3dipsig.size()+1,VBat3dipsig.size()+1,Vvtrkpi3dipsig.size()+1,Vvtrkp3dipsig.size()+1,Vflightsig.size()+1,Vdistancesig.size()+1};
//These are for eta cut
    //float minbins[nDim] = {1.6  , 1.0 , -5 , 0                   , 3.5   , 2.5   , 1.5   , 1.5   , 5.0};
    //float maxbins[nDim] = {1.75 , 10  , 5  , Vcas3dipsig.back() , 99999 , 99999 , 99999 , 99999 , 99999};


    //Now make variable sized bins for everything besides mass
    // Insert the maximum to cut parameters
    //Vcas3dipsig.insert(Vcas3dipsig.begin(),0);
    //VBat3dipsig.push_back(99999);
    //Vvtrkpi3dipsig.push_back(99999);
    //Vvtrkp3dipsig.push_back(99999);
    //Vflightsig.push_back(99999);
    //Vdistancesig.push_back(99999);


    WriteVector(VvarBins_pt,myfile,"Pt");
    WriteVector(VvarBins_eta,myfile,"Eta");

    if(PID == "Omega")
    {
        WriteVector(Vcas3dipsig,myfile,"Om DCA");
        WriteVector(VBat3dipsig,myfile,"Bat DCA");
        WriteVector(Vvtrkpi3dipsig,myfile,"Pi DCA");
        WriteVector(Vvtrkp3dipsig,myfile,"Pr DCA");
        WriteVector(Vflightsig,myfile,"Om DLS");
        WriteVector(Vdistancesig,myfile,"La DLS");
    }
    else if(PID == "Xi")
    {
        WriteVector(Vcas3dipsig,myfile,"Xi DCA");
        WriteVector(VBat3dipsig,myfile,"Bat DCA");
        WriteVector(Vvtrkpi3dipsig,myfile,"Pi DCA");
        WriteVector(Vvtrkp3dipsig,myfile,"Pr DCA");
        WriteVector(Vflightsig,myfile,"Xi DLS");
        WriteVector(Vdistancesig,myfile,"La DLS");
    }
    else if(PID == "Kshort")
    {
        WriteVector(vdz1,myfile,"Dau1 dz");
        WriteVector(vdz2,myfile,"Dau2 dz");
        WriteVector(vdxy1,myfile,"Dau1 dxy");
        WriteVector(vdxy2,myfile,"Dau2 dxy");
        WriteVector(vdls,myfile,"Ks DLS");
        WriteVector(vagl,myfile,"Ks Agl");
    }
    else if(PID == "Lambda")
    {
        WriteVector(vdz1,myfile,"Dau1 dz");
        WriteVector(vdz2,myfile,"Dau2 dz");
        WriteVector(vdxy1,myfile,"Dau1 dxy");
        WriteVector(vdxy2,myfile,"Dau2 dxy");
        WriteVector(vdls,myfile,"La DLS");
        WriteVector(vagl,myfile,"La Agl");
    }

    THnSparseF* h = new THnSparseF("hSparse","Cuts",nDim,&nBins[0],&minbins[0],&maxbins[0]);

    if(PID == "Omega" || PID == "Xi")
    {
        h->GetAxis(0)->SetName("Mass");
        h->GetAxis(1)->SetName("Pt");
        h->GetAxis(2)->SetName("Eta");
        h->GetAxis(3)->SetName("CasDCA");
        h->GetAxis(4)->SetName("BatDCA");
        h->GetAxis(5)->SetName("LamPiDCA");
        h->GetAxis(6)->SetName("LamProDCA");
        h->GetAxis(7)->SetName("CasDLS");
        h->GetAxis(8)->SetName("LamDLS");
        cout << "1" << endl;
        h->GetAxis(1)->Set(nBins[1],&VvarBins_pt[0]);
        cout << "2" << endl;
        h->GetAxis(2)->Set(nBins[2],&VvarBins_eta[0]);
        cout << "3" << endl;
        h->GetAxis(3)->Set(nBins[3],&Vcas3dipsig[0]);
        cout << "4" << endl;
        h->GetAxis(4)->Set(nBins[4],&VBat3dipsig[0]);
        cout << "5" << endl;
        h->GetAxis(5)->Set(nBins[5],&Vvtrkpi3dipsig[0]);
        cout << "6" << endl;
        h->GetAxis(6)->Set(nBins[6],&Vvtrkp3dipsig[0]);
        cout << "7" << endl;
        h->GetAxis(7)->Set(nBins[7],&Vflightsig[0]);
        cout << "8" << endl;
        h->GetAxis(8)->Set(nBins[8],&Vdistancesig[0]);
    }
    else
    {
        h->GetAxis(0)->SetName("Mass");
        h->GetAxis(1)->SetName("Pt");
        h->GetAxis(2)->SetName("Eta");
        h->GetAxis(3)->SetName("dz1");
        h->GetAxis(4)->SetName("dz2");
        h->GetAxis(5)->SetName("dxy1");
        h->GetAxis(6)->SetName("dxy2");
        h->GetAxis(7)->SetName("dls");
        h->GetAxis(8)->SetName("agl");
        cout << "1" << endl;
        h->GetAxis(1)->Set(nBins[1],&VvarBins_pt[0]);
        cout << "2" << endl;
        h->GetAxis(2)->Set(nBins[2],&VvarBins_eta[0]);
        cout << "3" << endl;
        h->GetAxis(3)->Set(nBins[3],&vdz1[0]);
        cout << "4" << endl;
        h->GetAxis(4)->Set(nBins[4],&vdz2[0]);
        cout << "5" << endl;
        h->GetAxis(5)->Set(nBins[5],&vdxy1[0]);
        cout << "6" << endl;
        h->GetAxis(6)->Set(nBins[6],&vdxy2[0]);
        cout << "7" << endl;
        h->GetAxis(7)->Set(nBins[7],&vdls[0]);
        cout << "8" << endl;
        h->GetAxis(8)->Set(nBins[8],&vagl[0]);
    }




    //Containers


    //Tree setup
    //TFile* f1=TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/V0CasTree_09_14_17.root");
    TFile* f1=TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/PbPbTreeCorrected.root");

    TTreeReader reader;
    if(PID == "Omega")
    {
        reader.SetTree("OmTreeProducerRapidityPbPb/OmTree",f1);
    }
    else if(PID == "Xi")
    {
        reader.SetTree("XiTreeProducerRapidityPbPb/XiTree",f1);
    }
    else if(PID == "Kshort")
    {
        reader.SetTree("KsTreeProducerRapidityPbPb/KsTree",f1);
    }
    else if(PID == "Lambda")
    {
        reader.SetTree("LaTreeProducerRapidityPbPb/LaTree",f1);
    }

    TTreeReaderValue<float> rap;
    TTreeReaderValue<float> eta;
    TTreeReaderValue<float> pt;
    TTreeReaderValue<float> mass;

    TTreeReaderValue<float> cas3dipsig;
    TTreeReaderValue<float> distancesig;
    TTreeReaderValue<float> flightsig;
    TTreeReaderValue<float> vtrkp3dipsig;
    TTreeReaderValue<float> vtrkpi3dipsig;
    TTreeReaderValue<float> bat3dipsig;

    TTreeReaderValue<float> dz1;
    TTreeReaderValue<float> dz2;
    TTreeReaderValue<float> dxy1;
    TTreeReaderValue<float> dxy2;
    TTreeReaderValue<float> dls;
    TTreeReaderValue<float> agl;

    if(PID == "Omega")
    {
        cas3dipsig    = SetValues(reader,"om3dipsig.om3dipsig");
        bat3dipsig    = SetValues(reader,"omKaon3dipsig.omKaon3dipsig");
        vtrkpi3dipsig = SetValues(reader,"vtrkpi3dipsig.vtrkpi3dipsig");
        vtrkp3dipsig  = SetValues(reader,"vtrkp3dipsig.vtrkp3dipsigpt");
        flightsig     = SetValues(reader,"omflightsig.omflightsig");
        distancesig   = SetValues(reader,"distancesig.distancesig");
        mass          = SetValues(reader,"mass.mass");
        pt            = SetValues(reader,"pt.pt");
        eta           = SetValues(reader,"eta.eta");
        rap           = SetValues(reader,"rapidity.rapidity");
    }
    else if(PID == "Xi")
    {
        cas3dipsig    = SetValues(reader,"xi3dipsig.xi3dipsig");
        bat3dipsig    = SetValues(reader,"xipi3dipsig.xipi3dipsig");
        vtrkpi3dipsig = SetValues(reader,"vtrkpi3dipsig.vtrkpi3dipsig");
        vtrkp3dipsig  = SetValues(reader,"vtrkp3dipsig.vtrkp3dipsigpt");
        flightsig     = SetValues(reader,"xiflightsig.xiflightsig");
        distancesig   = SetValues(reader,"distancesig.distancesig");
        mass          = SetValues(reader,"mass.mass");
        pt            = SetValues(reader,"pt.pt");
        eta           = SetValues(reader,"eta.eta");
        rap           = SetValues(reader,"rapidity.rapidity");
    }
    else if(PID == "Kshort")
    {
        rap  = SetValues(reader,"rapidity.rapidity");
        mass = SetValues(reader,"mass.mass");
        pt   = SetValues(reader,"pt.pt");
        eta  = SetValues(reader,"eta.eta");
        dz1  = SetValues(reader,"dzSig1.dzSig1");
        dz2  = SetValues(reader,"dzSig2.dzSig2");
        dxy1 = SetValues(reader,"dxySig1.dxySig1");
        dxy2 = SetValues(reader,"dxySig2.dxySig2");
        dls  = SetValues(reader,"decayLSig.decayLSig");
        agl  = SetValues(reader,"cosTheta.cosTheta");
    }
    else if(PID == "Lambda")
    {
        rap  = SetValues(reader,"rapidity.rapidity");
        mass = SetValues(reader,"mass.mass");
        pt   = SetValues(reader,"pt.pt");
        eta  = SetValues(reader,"eta.eta");
        dz1  = SetValues(reader,"dzSig1.dzSig1");
        dz2  = SetValues(reader,"dzSig2.dzSig2");
        dxy1 = SetValues(reader,"dxySig1.dxySig1");
        dxy2 = SetValues(reader,"dxySig2.dxySig2");
        dls  = SetValues(reader,"decayLSig.decayLSig");
        agl  = SetValues(reader,"cosTheta.cosTheta");
    }


    //reader.SetEntriesRange(0,5000);
    while(reader.Next())
    {
        if(PID == "Omega" || PID == "Xi")
        {
            if(!CheckValue(cas3dipsig))     return false;
            if(!CheckValue(bat3dipsig)) return false;
            if(!CheckValue(vtrkpi3dipsig)) return false;
            if(!CheckValue(vtrkp3dipsig))  return false;
            if(!CheckValue(flightsig))   return false;
            if(!CheckValue(distancesig))   return false;
            if(!CheckValue(mass))          return false;
            if(!CheckValue(pt))            return false;
            if(!CheckValue(eta))           return false;
            if(!CheckValue(rap)) return false;

            if(fabs(*rap) > 1.0) continue;

            double value[nDim] = {*mass,*pt,*eta,*cas3dipsig,*bat3dipsig,*vtrkpi3dipsig,*vtrkp3dipsig,*flightsig,*distancesig};
            h->Fill(value);
        }
        else
        {
            if(!CheckValue(mass))          return false;
            if(!CheckValue(pt))            return false;
            if(!CheckValue(eta))           return false;
            if(!CheckValue(dz1)) return false;
            if(!CheckValue(dz2)) return false;
            if(!CheckValue(dxy1)) return false;
            if(!CheckValue(dxy2)) return false;
            if(!CheckValue(dls)) return false;
            if(!CheckValue(agl)) return false;
            if(!CheckValue(rap)) return false;

            if(fabs(*rap) > 1.0) continue;

            double value[nDim] = {*mass,*pt,*eta,fabs(*dz1),fabs(*dz2),fabs(*dxy1),fabs(*dxy2),*dls,*agl};
            h->Fill(value);
        }
    }

    f1->Close();

    //Output file creation
    std::string filetype = ".root";
    name += filetype;
    TFile out(name.c_str(),"RECREATE");
    cout << name.c_str() << " created!" << endl;

    h->Write();
    out.Close();
    return true;
}


/*
bool LaSparseCreator(std::string name)
{
    //Initializers
    TH1::SetDefaultSumw2();
    std::ofstream myfile;
    std::string txtname = name + ".txt";
    myfile.open(txtname.c_str());

    //Cut parameters to be varied. The commented elements are Hong's cuts to make a plot of them remember to change the numparam value as well
    std::vector<float> VvarBins_pt            = {0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0};
    std::vector<float> VvarBins_eta           = {-2.4,-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,2.0,2.4};

    std::vector<float> vdca = {0.9   , 1.0   , 1.2    , 1.3    , 1.4    , 1.5    , 1.6     , 1.7};
    std::vector<float> vdls = {4     , 5     , 6      , 7      , 8      , 9      , 10      , 11};
    std::vector<float> vagl = {0.998 , 0.999 , 0.9992 , 0.9995 , 0.9998 , 0.9999 , 0.99995 , 0.99999};

    int numparam = vdca.size();
    const int nDim = 6;
    // mass, pt, eta, then cuts
    //These are for eta cut
    //Lambda
    std::vector<int> nBins = {160,13,10,static_cast<int>(Vcas3dipsig.size()),static_cast<int>(VBat3dipsig.size()),static_cast<int>(Vvtrkpi3dipsig.size()),static_cast<int>(Vvtrkp3dipsig.size()),static_cast<int>(Vflightsig.size()),static_cast<int>(Vdistancesig.size())};
    //Omega
    //double minbins[nDim] = {1.6  , 1.0 , -2.4 , 0                   , 3.5   , 2.5   , 1.5   , 1.5   , 5.0};
    //double maxbins[nDim] = {1.75 , 10  , 2.4  , Vcas3dipsig.back() , VBat3dipsig.back() , Vvtrkpi3dipsig.back() , Vvtrkp3dipsig.back() , Vflightsig.back() , Vdistancesig.back()};
    //Xi
    double minbins[nDim] = {1.08 , 1.0 , -2.4 , 0                  , VBat3dipsig.front() , Vvtrkpi3dipsig.front() , Vvtrkp3dipsig.front() , Vflightsig.front() , Vdistancesig.front()};
    double maxbins[nDim] = {1.16 , 10  , 2.4  , Vcas3dipsig.back() , 99999               , 99999                  , 99999                 , 99999              , 99999};

//These are for rap cut
    //int nBins[nDim] = {150,9,12,Vcas3dipsig.size()+1,VBat3dipsig.size()+1,Vvtrkpi3dipsig.size()+1,Vvtrkp3dipsig.size()+1,Vflightsig.size()+1,Vdistancesig.size()+1};
//These are for eta cut
    //float minbins[nDim] = {1.6  , 1.0 , -5 , 0                   , 3.5   , 2.5   , 1.5   , 1.5   , 5.0};
    //float maxbins[nDim] = {1.75 , 10  , 5  , Vcas3dipsig.back() , 99999 , 99999 , 99999 , 99999 , 99999};


    //Now make variable sized bins for everything besides mass
    // Insert the maximum to cut parameters
    Vcas3dipsig.insert(Vcas3dipsig.begin(),0);
    VBat3dipsig.push_back(99999);
    Vvtrkpi3dipsig.push_back(99999);
    Vvtrkp3dipsig.push_back(99999);
    Vflightsig.push_back(99999);
    Vdistancesig.push_back(99999);


    WriteVector(VvarBins_pt,myfile,"Pt");
    WriteVector(VvarBins_eta,myfile,"Eta");

    //Omega
    //WriteVector(Vcas3dipsig,myfile,"Om DCA");
    //WriteVector(VBat3dipsig,myfile,"Ka DCA");
    //WriteVector(Vvtrkpi3dipsig,myfile,"Pi DCA");
    //WriteVector(Vvtrkp3dipsig,myfile,"Pr DCA");
    //WriteVector(Vflightsig,myfile,"Om DLS");
    //WriteVector(Vdistancesig,myfile,"La DLS");

    //Xi
    WriteVector(Vcas3dipsig,myfile,"Xi DCA");
    WriteVector(VBat3dipsig,myfile,"Bat DCA");
    WriteVector(Vvtrkpi3dipsig,myfile,"Pi DCA");
    WriteVector(Vvtrkp3dipsig,myfile,"Pr DCA");
    WriteVector(Vflightsig,myfile,"Xi DLS");
    WriteVector(Vdistancesig,myfile,"La DLS");

    THnSparseF* h = new THnSparseF("hSparse","Cuts",nDim,&nBins[0],minbins,maxbins);

    cout << "1" << endl;
    h->GetAxis(1)->Set(nBins[1],&VvarBins_pt[0]);
    cout << "2" << endl;
    h->GetAxis(2)->Set(nBins[2],&VvarBins_eta[0]);
    cout << "3" << endl;
    h->GetAxis(3)->Set(nBins[3],&Vcas3dipsig[0]);
    cout << "4" << endl;
    h->GetAxis(4)->Set(nBins[4],&VBat3dipsig[0]);
    cout << "5" << endl;
    h->GetAxis(5)->Set(nBins[5],&Vvtrkpi3dipsig[0]);
    cout << "6" << endl;
    h->GetAxis(6)->Set(nBins[6],&Vvtrkp3dipsig[0]);
    cout << "7" << endl;
    h->GetAxis(7)->Set(nBins[7],&Vflightsig[0]);
    cout << "8" << endl;
    h->GetAxis(8)->Set(nBins[8],&Vdistancesig[0]);



    //Containers


    //Tree setup
    //TFile* f1=TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/V0CasTree_09_14_17.root");
    TFile* f1=TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/V0CasTreePbPb_09_25_17.root");
    //TFile* f1= new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/PbPbTreeSmall.root");

    //TTreeReader reader("OmTreeProducerRapidityPbPb/OmTree",f1);
    TTreeReader reader("XiTreeProducerRapidityPbPb/XiTree",f1);

    //Omega
    //TTreeReaderArray<float> cas3dipsig(reader    ,"om3dipsig.om3dipsig");
    //TTreeReaderArray<float> bat3dipsig(reader,"omKaon3dipsig.omKaon3dipsig");
    //TTreeReaderArray<float> vtrkpi3dipsig(reader,"vtrkpi3dipsig.vtrkpi3dipsig");
    //TTreeReaderArray<float> vtrkp3dipsig(reader ,"vtrkp3dipsig.vtrkp3dipsigpt");
    //TTreeReaderArray<float> flightsig(reader  ,"omflightsig.omflightsig");
    //TTreeReaderArray<float> distancesig(reader  ,"distancesig.distancesig");
    //TTreeReaderArray<float> mass(reader         ,"mass.mass");
    //TTreeReaderArray<float> pt(reader           ,"pt.pt");
    //TTreeReaderArray<float> eta(reader          ,"eta.eta");
    //TTreeReaderArray<float> rap(reader          ,"rapidity.rapidity");

    //Xi
    TTreeReaderArray<float> cas3dipsig(reader    ,"xi3dipsig.xi3dipsig");
    TTreeReaderArray<float> bat3dipsig(reader,"xipi3dipsig.xipi3dipsig");
    TTreeReaderArray<float> vtrkpi3dipsig(reader,"vtrkpi3dipsig.vtrkpi3dipsig");
    TTreeReaderArray<float> vtrkp3dipsig(reader ,"vtrkp3dipsig.vtrkp3dipsigpt");
    TTreeReaderArray<float> flightsig(reader  ,"xiflightsig.xiflightsig");
    TTreeReaderArray<float> distancesig(reader  ,"distancesig.distancesig");
    TTreeReaderArray<float> mass(reader         ,"mass.mass");
    TTreeReaderArray<float> pt(reader           ,"pt.pt");
    TTreeReaderArray<float> eta(reader          ,"eta.eta");
    TTreeReaderArray<float> rap(reader          ,"rapidity.rapidity");


    //int j=0;
    while(reader.Next())
    //while(j<10000)
    {
        //j++;
        //reader.Next();
        if(!CheckValue(cas3dipsig))     return false;
        if(!CheckValue(bat3dipsig)) return false;
        if(!CheckValue(vtrkpi3dipsig)) return false;
        if(!CheckValue(vtrkp3dipsig))  return false;
        if(!CheckValue(flightsig))   return false;
        if(!CheckValue(distancesig))   return false;
        if(!CheckValue(mass))          return false;
        if(!CheckValue(pt))            return false;
        if(!CheckValue(eta))           return false;

        for(int i=0; i<mass.GetSize(); i++)
        {
            double value[nDim] = {mass[i],pt[i],eta[i],cas3dipsig[i],bat3dipsig[i],vtrkpi3dipsig[i],vtrkp3dipsig[i],flightsig[i],distancesig[i]};
            h->Fill(value);
        }
    }

    f1->Close();

    //Output file creation
    std::string filetype = ".root";
    name += filetype;
    TFile out(name.c_str(),"RECREATE");
    cout << name.c_str() << " created!" << endl;

    h->Write();
    out.Close();
    return true;
}

*/

/*
bool KsSparseCreator(std::string name)
{
    //Initializers
    TH1::SetDefaultSumw2();
    std::ofstream myfile;
    std::string txtname = name + ".txt";
    myfile.open(txtname.c_str());

    //Cut parameters to be varied. The commented elements are Hong's cuts to make a plot of them remember to change the numparam value as well
    //std::vector<float> VvarBins_pt            = {0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0};
    //std::vector<float> VvarBins_eta           = {-2.4,-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,2.0,2.4};

    std::vector<float> vdz1  = {1.0    , 1.2   , 1.5   , 1.7    , 2.0};
    std::vector<float> vdz2  = {1.0    , 1.2   , 1.5   , 1.7    , 2.0};
    std::vector<float> vdxy1 = {1.0    , 1.2   , 1.5   , 1.7    , 2.0};
    std::vector<float> vdxy2 = {1.0    , 1.2   , 1.5   , 1.7    , 2.0};
    std::vector<float> vdls  = {4      , 6     , 7     , 9      , 10};
    std::vector<float> vagl  = {0.9965 , 0.997 , 0.999 , 0.9995 , 0.9999};

    const int nDim = 9;
    // mass, pt, eta, then cuts
    //These are for eta cut
    //Lambda
    std::vector<int> nBins = {270,14,10,(int)vdz1.size(),(int)vdz2.size(),(int)vdxy1.size(),(int)vdxy2.size(),(int)vdls.size(),(int)vagl.size()};
    //Omega
    //double minbins[nDim] = {1.6  , 1.0 , -2.4 , 0                   , 3.5   , 2.5   , 1.5   , 1.5   , 5.0};
    //double maxbins[nDim] = {1.75 , 10  , 2.4  , Vcas3dipsig.back() , VBat3dipsig.back() , Vvtrkpi3dipsig.back() , Vvtrkp3dipsig.back() , Vflightsig.back() , Vdistancesig.back()};
    //Xi
    double minbins[nDim] = {0.43  , 1.0 , -2.4 , vdz1.front() , vdz2.front() , vdxy1.front() , vdxy2.front() , vdls.front() , vagl.front()};
    double maxbins[nDim] = {0.565 , 10  , 2.4  , 99999        , 99999        , 99999         , 99999         , 99999        , 1};

//These are for rap cut
    //int nBins[nDim] = {150,9,12,Vcas3dipsig.size()+1,VBat3dipsig.size()+1,Vvtrkpi3dipsig.size()+1,Vvtrkp3dipsig.size()+1,Vflightsig.size()+1,Vdistancesig.size()+1};
//These are for eta cut
    //float minbins[nDim] = {1.6  , 1.0 , -5 , 0                   , 3.5   , 2.5   , 1.5   , 1.5   , 5.0};
    //float maxbins[nDim] = {1.75 , 10  , 5  , Vcas3dipsig.back() , 99999 , 99999 , 99999 , 99999 , 99999};


    //Now make variable sized bins for everything besides mass
    // Insert the maximum to cut parameters
    vdz1.push_back(99999);
    vdz2.push_back(99999);
    vdxy1.push_back(99999);
    vdxy2.push_back(99999);
    vdls.push_back(99999);
    vagl.push_back(1);


    WriteVector(VvarBins_pt,myfile,"Pt");
    WriteVector(VvarBins_eta,myfile,"Eta");

    //Omega
    //WriteVector(Vcas3dipsig,myfile,"Om DCA");
    //WriteVector(VBat3dipsig,myfile,"Ka DCA");
    //WriteVector(Vvtrkpi3dipsig,myfile,"Pi DCA");
    //WriteVector(Vvtrkp3dipsig,myfile,"Pr DCA");
    //WriteVector(Vflightsig,myfile,"Om DLS");
    //WriteVector(Vdistancesig,myfile,"La DLS");

    //Xi
    WriteVector(vdz1,myfile,"Daughter1 dz");
    WriteVector(vdz2,myfile,"Daughter2 dz");
    WriteVector(vdxy1,myfile,"Daughter1 dxy");
    WriteVector(vdxy2,myfile,"Daughter2 dxy");
    WriteVector(vdls,myfile,"DLS");
    WriteVector(vagl,myfile,"agl");

    THnSparseF* h = new THnSparseF("hSparse","Cuts",nDim,&nBins[0],minbins,maxbins);

    cout << "1" << endl;
    h->GetAxis(1)->Set(nBins[1],&VvarBins_pt[0]);
    cout << "2" << endl;
    h->GetAxis(2)->Set(nBins[2],&VvarBins_eta[0]);
    cout << "3" << endl;
    h->GetAxis(3)->Set(nBins[3],&vdz1[0]);
    cout << "4" << endl;
    h->GetAxis(4)->Set(nBins[4],&vdz2[0]);
    cout << "5" << endl;
    h->GetAxis(5)->Set(nBins[5],&vdxy1[0]);
    cout << "6" << endl;
    h->GetAxis(6)->Set(nBins[6],&vdxy2[0]);
    cout << "7" << endl;
    h->GetAxis(7)->Set(nBins[7],&vdls[0]);
    cout << "8" << endl;
    h->GetAxis(8)->Set(nBins[8],&vagl[0]);



    //Containers


    //Tree setup
    //TFile* f1=TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/V0CasTree_09_14_17.root");
    TFile* f1=TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/V0CasTreePbPb_09_25_17.root");
    //TFile* f1= new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/PbPbTreeSmall.root");

    //TTreeReader reader("OmTreeProducerRapidityPbPb/OmTree",f1);
    TTreeReader reader("XiTreeProducerRapidityPbPb/XiTree",f1);

    //Kshort
    TTreeReaderArray<float> dz1(reader ,"dzSig1.dzSig1");
    TTreeReaderArray<float> dz2(reader ,"dzSig2.dzSig2");
    TTreeReaderArray<float> dxy1(reader ,"dxySig1.dxySig1");
    TTreeReaderArray<float> dxy2(reader ,"dxySig2.dxySig2");
    TTreeReaderArray<float> distancesig(reader  ,"decayLSig.decayLSig");
    TTreeReaderArray<float> agl(reader  ,"cosTheta.cosTheta");
    TTreeReaderArray<float> mass(reader         ,"mass.mass");
    TTreeReaderArray<float> pt(reader           ,"pt.pt");
    TTreeReaderArray<float> eta(reader          ,"eta.eta");


    //int j=0;
    while(reader.Next())
    //while(j<10000)
    {
        //j++;
        //reader.Next();
        if(!CheckValue(dz1))     return false;
        if(!CheckValue(dz2)) return false;
        if(!CheckValue(dxy1)) return false;
        if(!CheckValue(dxy2))  return false;
        if(!CheckValue(distancesig))   return false;
        if(!CheckValue(agl))   return false;
        if(!CheckValue(mass))          return false;
        if(!CheckValue(pt))            return false;
        if(!CheckValue(eta))           return false;

        for(int i=0; i<mass.GetSize(); i++)
        {
            double value[nDim] = {mass[i],pt[i],eta[i],fabs(dz1[i]),fabs(dz2[i]),fabs(dxy1[i]),fabs(dxy2[i]),distancesig[i],agl[i]};
            h->Fill(value);
        }
    }

    f1->Close();

    //Output file creation
    std::string filetype = ".root";
    name += filetype;
    TFile out(name.c_str(),"RECREATE");
    cout << name.c_str() << " created!" << endl;

    h->Write();
    out.Close();
    return true;
}
*/

void readSparse(std::string name, std::string PID, int version)
{
    TFile* f = new TFile(name.c_str());

    std::vector<float> VvarBins_eta           = {-2.4,-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,2.0,2.4};

    std::vector<float> VvarBins_pt;
    std::vector<float> Vcas3dipsig;
    std::vector<float> VBat3dipsig;
    std::vector<float> Vvtrkpi3dipsig;
    std::vector<float> Vvtrkp3dipsig;
    std::vector<float> Vflightsig;
    std::vector<float> Vdistancesig;

    std::vector<float> vdz1;
    std::vector<float> vdz2;
    std::vector<float> vdxy1;
    std::vector<float> vdxy2;
    std::vector<float> vdls;
    std::vector<float> vagl;

    if(PID == "Omega") //831111
    {
        VvarBins_pt.insert(VvarBins_pt.begin(),       {1.0,1.5,1.8,2.2,2.8,3.6,4.6,6.0,7.2,10.0});

        if(version == 2)
        {
            Vcas3dipsig.insert(Vcas3dipsig.begin(),       {0.5 , 1.0  , 1.5  , 1.8 , 2.0 , 2.5 , 2.8 , 3.0 , 3.5 });
            VBat3dipsig.insert(VBat3dipsig.begin(),       {3.0 , 3.5  , 4.0  , 5.0 , 6.0 , 7.0 , 8.0 , 9.0 , 10.0});
            Vvtrkpi3dipsig.insert(Vvtrkpi3dipsig.begin(), {3.0 , 3.5  , 4.0 , 4.5  , 5.0 , 6.0 , 7.0 , 9.0 , 10.0});
            Vvtrkp3dipsig.insert(Vvtrkp3dipsig.begin(),   {2.0 , 2.5  , 3.0 , 3.5  , 4.0 , 5.0 , 7.0 , 9.5 , 10.0});
            Vflightsig.insert(Vflightsig.begin(),         {2.0 , 2.5  , 3.0 , 3.5  , 4.0 , 5.0 , 7.0 , 9.5 , 10.0});
            Vdistancesig.insert(Vdistancesig.begin(),     {10. , 10.5 , 11  , 11.5 , 12. , 13  , 16  , 18  , 20  });
        }

        if(version == 3) //441111
        {
            Vcas3dipsig.insert(Vcas3dipsig.begin(),       {1.5 , 1.6, 1.7 , 1.8 , 1.9 , 2.0});
            VBat3dipsig.insert(VBat3dipsig.begin(),       {3.5 , 3.6, 3.8 , 4.0 , 4.2 , 4.4});
            Vvtrkpi3dipsig.insert(Vvtrkpi3dipsig.begin(), {3.0 , 3.1 , 3.2 , 3.3 , 3.4 , 3.5});
            Vvtrkp3dipsig.insert(Vvtrkp3dipsig.begin(),   {2.0 , 2.1 , 2.2 , 2.3 , 2.4 , 2.5});
            Vflightsig.insert(Vflightsig.begin(),         {2.0 , 2.1 , 2.2 , 2.3 , 2.4 , 2.5});
            Vdistancesig.insert(Vdistancesig.begin(),     {10. , 10.5 , 11  , 11.5 , 12. , 13});
        }
        if(version == 4) // Corrected tree
        {
            Vcas3dipsig.insert(Vcas3dipsig.begin(),       {0.5 , 1.0  , 1.5  , 1.8 , 2.0 , 2.5 , 2.8 , 3.0 , 3.5 });
            VBat3dipsig.insert(VBat3dipsig.begin(),       {3.0 , 3.5  , 4.0  , 5.0 , 6.0 , 7.0 , 8.0 , 9.0 , 10.0});
            Vvtrkpi3dipsig.insert(Vvtrkpi3dipsig.begin(), {3.0 , 3.5  , 4.0 , 4.5  , 5.0 , 6.0 , 7.0 , 9.0 , 10.0});
            Vvtrkp3dipsig.insert(Vvtrkp3dipsig.begin(),   {2.0 , 2.5  , 3.0 , 3.5  , 4.0 , 5.0 , 7.0 , 9.5 , 10.0});
            Vflightsig.insert(Vflightsig.begin(),         {2.0 , 2.5  , 3.0 , 3.5  , 4.0 , 5.0 , 7.0 , 9.5 , 10.0});
            Vdistancesig.insert(Vdistancesig.begin(),     {10. , 10.5 , 11  , 11.5 , 12. , 13  , 16  , 18  , 20  });
        }
    }
    else if(PID == "Xi")
    {
        VvarBins_pt.insert(VvarBins_pt.begin(),       {1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.2,10.0});

        if(version == 2)
        {
            Vcas3dipsig.insert(Vcas3dipsig.begin(),       {0.5 , 1.0 , 1.5 , 2.0 , 2.25, 2.5 , 3.0 , 3.5 , 4.0});
            VBat3dipsig.insert(VBat3dipsig.begin(),       {4.5 , 5.0 , 5.5 , 6.0 , 6.5 , 7.0 , 7.5 , 8.0 , 9.0});
            Vvtrkpi3dipsig.insert(Vvtrkpi3dipsig.begin(), {3.5 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 6.5 , 7.0 , 8.0});
            Vvtrkp3dipsig.insert(Vvtrkp3dipsig.begin(),   {2.5 , 3.0 , 3.5 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 7.0});
            Vflightsig.insert(Vflightsig.begin(),         {2.5 , 3.0 , 3.5 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 7.0});
            Vdistancesig.insert(Vdistancesig.begin(),     {10. , 11. , 12. , 13. , 13.5 , 14 , 14.5 , 15 , 16 });
        }

        if(version == 3)
        {
            Vcas3dipsig.insert(Vcas3dipsig.begin(),       {1.5 , 2.25, 2.5 , 2.8 , 3.5 , 4.0});
            VBat3dipsig.insert(VBat3dipsig.begin(),       {4.8 , 5.0 , 5.2 , 6.0 , 6.5 , 7.0 , 7.5 , 8.0 , 9.0});
            Vvtrkpi3dipsig.insert(Vvtrkpi3dipsig.begin(), {3.8 , 4.0 , 4.2 , 5.0 , 5.5 , 6.0 , 6.5 , 7.0 , 8.0});
            Vvtrkp3dipsig.insert(Vvtrkp3dipsig.begin(),   {2.8 , 3.0 , 3.2 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 7.0});
            Vflightsig.insert(Vflightsig.begin(),         {2.5 , 3.0 , 3.5 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 7.0});
            Vdistancesig.insert(Vdistancesig.begin(),     {10. , 11.5 , 12. , 12.5 , 13.5 , 14 , 14.5 , 15 , 16 });
        }
        if(version == 4)
        {
            Vcas3dipsig.insert(Vcas3dipsig.begin(),       {1.0 ,1.5 , 2.25, 2.5 , 2.8 , 3.5 , 4.0});
            VBat3dipsig.insert(VBat3dipsig.begin(),       {4.8 , 5.0 , 5.2 , 6.0 , 6.5 , 7.0 , 7.5 , 8.0 , 9.0});
            Vvtrkpi3dipsig.insert(Vvtrkpi3dipsig.begin(), {3.8 , 4.0 , 4.2 , 5.0 , 5.5 , 6.0 , 6.5 , 7.0 , 8.0});
            Vvtrkp3dipsig.insert(Vvtrkp3dipsig.begin(),   {2.8 , 3.0 , 3.2 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 7.0});
            Vflightsig.insert(Vflightsig.begin(),         {2.5 , 3.0 , 3.5 , 4.0 , 4.5 , 5.0 , 5.5 , 6.0 , 7.0});
            Vdistancesig.insert(Vdistancesig.begin(),     {10. , 11.5 , 12. , 12.5 , 13.5 , 14 , 14.5 , 15 , 16 });
        }
    }
    else if (PID == "Kshort")
    {
        VvarBins_pt.insert(VvarBins_pt.begin(),{0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0});
        vdz1.insert(vdz1.begin()  , {1.0    , 1.2   , 1.5   , 1.7   , 2.0    , 2.5    , 3.0});
        vdz2.insert(vdz2.begin()  , {1.0    , 1.2   , 1.5   , 1.7   , 2.0    , 2.5    , 3.0});
        vdxy1.insert(vdxy1.begin() , {1.0    , 1.2   , 1.5   , 1.7   , 2.0    , 2.5    , 3.0});
        vdxy2.insert(vdxy2.begin() , {1.0    , 1.2   , 1.5   , 1.7   , 2.0    , 2.5    , 3.0});
        vdls.insert(vdls.begin()  , {4      , 5     , 6     , 7     , 8      , 9      , 10});
        vagl.insert(vagl.begin()  , {0.9965 , 0.997 , 0.998 , 0.999 , 0.9995 , 0.9999 , 0.99999});
    }
    else if(PID == "Lambda")
    {
        VvarBins_pt.insert(VvarBins_pt.begin(),       {0.8,1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.0,8.5,10.0});
        if(version == 1)
        {
            vdz1.insert(vdz1.begin()   , {1.0 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.8 , 1.9 , 2.0});
            vdz2.insert(vdz2.begin()   , {1.0 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.8 , 1.9 , 2.0});
            vdxy1.insert(vdxy1.begin() , {1.0 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.8 , 1.9 , 2.0});
            vdxy2.insert(vdxy2.begin() , {1.0 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.8 , 1.9 , 2.0});
            vdls.insert(vdls.begin()   , {5   , 6   , 7   , 8   , 9   , 10  , 11  , 12  , 13});
            vagl.insert(vagl.begin()   , {0.997 , 0.998 , 0.999 , 0.9995 , 0.9999 , 0.99995 , 0.99999});
        }

        if(version == 2)
        {
            vdz1.insert(vdz1.begin()   , {1.0 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.8 , 1.9 , 2.0});
            vdz2.insert(vdz2.begin()   , {1.0 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.8 , 1.9 , 2.0});
            vdxy1.insert(vdxy1.begin() , {1.0 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.8 , 1.9 , 2.0});
            vdxy2.insert(vdxy2.begin() , {1.0 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.8 , 1.9 , 2.0});
            vdls.insert(vdls.begin()   , {4   , 4.5  , 4.8 , 5   , 5.2 , 5.4 , 5.8 , 6.0 , 6.5});
            vagl.insert(vagl.begin()   , {0.999 , 0.9992 , 0.9993 , 0.9994, 0.9995 , 0.9996, 0.9997, 0.9998 ,0.9999});
        }

        if(version == 3)
        {
            vdz1.insert(vdz1.begin()   , {1.0 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.8 , 1.9 , 2.0});
            vdz2.insert(vdz2.begin()   , {1.0 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.8 , 1.9 , 2.0});
            vdxy1.insert(vdxy1.begin() , {1.0 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.8 , 1.9 , 2.0});
            vdxy2.insert(vdxy2.begin() , {1.0 , 1.2 , 1.3 , 1.4 , 1.5 , 1.6 , 1.8 , 1.9 , 2.0});
            vdls.insert(vdls.begin()   , {2.5, 3.0, 3.25, 3.5, 3.8, 4 , 4.5  , 4.8 , 5});
            vagl.insert(vagl.begin()   , {0.999 , 0.9992 , 0.9993 , 0.9994, 0.9995 , 0.9996, 0.9997, 0.9998 ,0.9999});
        }

        if(version == 5) //Fixed Tree, skipped 4 for some reason
        {
            vdz1.insert(vdz1.begin()   , {1.0 , 1.2 , 1.3 , 1.4 , 1.5});
            vdz2.insert(vdz2.begin()   , {1.0 , 1.2 , 1.3 , 1.4 , 1.5});
            vdxy1.insert(vdxy1.begin() , {1.0 , 1.2 , 1.3 , 1.4 , 1.5});
            vdxy2.insert(vdxy2.begin() , {1.0 , 1.2 , 1.3 , 1.4 , 1.5});
            vdls.insert(vdls.begin()   , {4.5, 4.8 , 5, 5.2, 5.4});
            vagl.insert(vagl.begin()   , {0.9996, 0.9997, 0.9998, 0.9999, 0.99992});
        }

    }

    std::ostringstream os;
    std::ostringstream osYield;
    using namespace RooFit;
    TH1::SetDefaultSumw2();
    //TGaxis::SetMaxDigits(3);
    RooMsgService::instance().setStreamStatus(0,kFALSE);
    RooMsgService::instance().setStreamStatus(1,kFALSE);
    gStyle->SetOptStat(1111111);
    THnSparseF* h = (THnSparseF*)f->Get("hSparse");
    //c2->Divide(1,2);
    //c2->cd(1);
    //h1->Draw();
    //c2->cd(2);
    //h->Projection(1)->Draw();
    //for(int a=1; a<VvarBins_pt.size(); a++)
    TCanvas* c2 = new TCanvas("c2","",600,600);
    std::map<int,double> SigSig;
    int pTBin = 0;
    ofstream File;
    //cout << "3: " << h->Projection(4)->GetBinContent(h->GetAxis(4)->GetNbins()+1) << endl;
    //cout << "4: " << h->GetAxis(4)->GetNbins() << endl;
    //cout << "5: " << h->GetAxis(5)->GetNbins() << endl;
    //cout << "6: " << h->GetAxis(6)->GetNbins() << endl;
    //cout << "7: " << h->GetAxis(7)->GetNbins() << endl;
    //cout << "8: " << h->GetAxis(8)->GetNbins() << endl;
    //return;
    //for(int a=1; a<VvarBins_pt.size(); a++)
    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.025);
    //tex->SetTextAlign(10);
    if(PID == "Omega" || PID == "Xi")
    {
        //for(int a=1; a<VvarBins_pt.size()-1; a++)
            for(int a=1; a<2; a++)
        {
            //if(a == 3) continue;
            pTBin = a;
            std::ostringstream oss;
            oss << "MaskAndSig_" << pTBin << "_" << PID << "_v" << version <<  ".txt";
            File.open(oss.str().c_str(),ios_base::app);
            oss.str(std::string());
            File << "\n";
            h->GetAxis(1)->SetRange(a,a);
            //for(int b=1; b<Vcas3dipsig.size(); b++)
            //for(int b=5; b<6; b++)
            for(int b=1; b<2; b++)
            {
                h->GetAxis(3)->SetRange(0,b);
                //for(int c=1; c<h->GetAxis(4)->GetNbins()+1; c++)
                //for(int c=3; c<4; c++)
                for(int c=2; c<3; c++)
                {
                    h->GetAxis(4)->SetRange(c,h->GetAxis(4)->GetNbins()+1);
                    //for(int d=1; d<h->GetAxis(5)->GetNbins()+1; d++)
                    //for(int d=1; d<2; d++)
                    for(int d=2; d<3; d++)
                    {
                        h->GetAxis(5)->SetRange(d,h->GetAxis(5)->GetNbins()+1);
                        //for(int e=1; e<h->GetAxis(6)->GetNbins()+1; e++)
                        //for(int e=1; e<2;e++)
                        for(int e=2; e<3;e++)
                        {
                            h->GetAxis(6)->SetRange(e,h->GetAxis(6)->GetNbins()+1);
                            //for(int f=1; f<h->GetAxis(7)->GetNbins()+1; f++)
                            //for(int f=1; f<2; f++)
                            for(int f=2; f<3; f++)
                            {
                                h->GetAxis(7)->SetRange(f,h->GetAxis(7)->GetNbins()+1);
                                //for(int i=1; i<h->GetAxis(8)->GetNbins()+1; i++)
                                //for(int i=1; i<2; i++)
                                for(int i=3; i<4; i++)
                                {
                                    h->GetAxis(8)->SetRange(i,h->GetAxis(8)->GetNbins()+1);
                                    struct stat buffer;
                                    std::ostringstream name;
                                    std::string filename;
                                    int mask = b*1e5 + c*1e4 + d*1e3 + e*1e2 + f*1e1 + i;
                                    //std::string mask = b*1e5 + c*1e4 + d*1e3 + e*1e2 + f*1e1 + i;
                                    if(PID == "Omega") name << "Distributions/OmegaDistributions/Omega_" << mask << "_" << a << "_v" << version << ".pdf";
                                    if(PID == "Xi") name << "Distributions/XiDistributions/Xi_" << mask << "_" << a << "_v" << version << ".pdf";
                                    //name << "OmegaDistributions/Omega_" << Vcas3dipsig[b-1] << "_" << VBat3dipsig[c-1] << "_" << Vvtrkpi3dipsig[d-1] << "_" << Vvtrkp3dipsig[e-1] << "_" << Vflightsig[f-1] << "_" << Vdistancesig[i-1] << "_" << a << "_v" << version << "2.pdf";
                                    //filename = name.str();
                                    //if(stat(name.str().c_str(),&buffer) == 0) continue;
                                    name.str(std::string());
                                    TH1D* h1 = (TH1D*)h->Projection(0);
                                    //h1->GetYaxis()->SetRangeUser(0,1000);
                                    gStyle->SetOptTitle(kFALSE);

                                    double s1;
                                    double s2;

                                    RooRealVar x("x","mass",0,1);
                                    //RooRealVar mean("mean","mean",1.67,1.6,1.75);//Omega
                                    RooRealVar mean("mean","mean",1.32,1.29,1.33);//Xi
                                    if(PID == "Omega")
                                    {
                                        x.setRange(1.60,1.75);
                                        //mean.setVal(1.67);
                                        //mean.setRange(1.60,1.75); //Does this need to be changed to a narrower range?
                                        s1=0.004;
                                        s2=0.005;
                                    }
                                    else if(PID == "Xi")
                                    {
                                        x.setRange(1.25,1.40);
                                        //mean.setVal(1.32);
                                        //mean.setRange(1.29,1.33);
                                        s1=0.003;
                                        s2=0.003;
                                    }
                                    RooPlot* xframe_ = x.frame(150);
                                    if(PID == "Omega") 
                                    {
                                        xframe_->GetXaxis()->SetTitle("#Lambda K Invariant mass (GeV)");
                                        //xframe_->GetYaxis()->SetRangeUser(0,600);
                                    }
                                    else if(PID == "Xi") xframe_->GetXaxis()->SetTitle("#Lambda #pi Invariant mass (GeV)");
                                    xframe_->GetYaxis()->SetTitle("Candidates / 0.001 GeV");
                                    RooDataHist data("data","dataset",x,h1);
                                    data.plotOn(xframe_,Name("data"));
                                    RooRealVar sigma1("sigma1","sigma1",s1,0.001,0.04);
                                    RooRealVar sigma2("sigma2","sigma2",s2,0.001,0.04);
                                    RooRealVar sig1("sig1","signal1",10,-100,10000000);
                                    RooRealVar sig2("sig2","signal2",10,-100,10000000);
                                    RooRealVar qsig("qsig","qsig",10,0,1000000);
                                    RooRealVar alpha("alpha","alpha",0.5,0,2);
                                    RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
                                    RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
                                    //Omega
                                    RooRealVar ap("ap","ap",-0.1,-1,1);
                                    RooRealVar bp("bp","bp",-0.1,-1,1);
                                    RooRealVar cp("cp","cp",-0.1,-1,1);
                                    RooRealVar dp("dp","dp",-0.1,-1,1);
                                    //RooRealVar ep("dp","dp",-0.1,-1,1);
                                    //RooRealVar fp("dp","dp",-0.1,-1,1);
                                    //Xi
                                    //RooRealVar ap("ap","ap",0,-5,10);
                                    //RooRealVar bp("bp","bp",0.1,-5,10);
                                    //RooRealVar cp("cp","cp",-0.1,-5,10);
                                    //RooRealVar dp("dp","dp",0.1,-5,10);
                                    //RooRealVar ap("ap","ap",0,-100000,100000);
                                    //RooRealVar bp("bp","bp",0,-100000,100000);
                                    //RooRealVar cp("cp","cp",0,-100000,100000);
                                    //RooRealVar dp("dp","dp",0,-100000,100000);
                                    RooChebychev background("background","background",x,RooArgList(ap,bp,cp,dp));
                                    //RooPolynomial background("background","background",x,RooArgList(ap,bp,cp));//,dp));
                                    //RooRealVar qsig("polysig","polysig",10,0,1000000000);
                                    //RooGenericPdf background("background", "x - (1.115683 + 0.493677)^alpha", RooArgList(x,alpha));
                                    //RooGenericPdf background("background", "x - (1.115683 + 0.13957018)^alpha", RooArgList(x,alpha));
                                    RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,background),RooArgList(sig1,sig2,qsig));

                                    if(PID == "Omega") x.setRange("cut",1.60,1.75);
                                    if(PID == "Xi") x.setRange("cut",1.25,1.40);

                                    //RooFitResult* r_xi = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));
                                    //mean.setConstant(kTRUE);
                                    RooFitResult* r_xi = sum.fitTo(data,Save(),Range("cut"),PrintLevel(-1));
                                    //mean.setConstant(kFALSE);
                                    //r_xi = sum.fitTo(data,Save(),Range("cut"));
                                    //RooFitResult* r_xi = sum.fitTo(data,Save(),Range("cut"));

                                    int count = 0;
                                    while(r_xi->covQual() < 3 ){
                                        s1+=0.004;
                                        s2+=0.003;
                                        //r_xi = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));
                                        r_xi = sum.fitTo(data,Save(),Range("cut"));
                                        //if(r_xi->covQual() == 3) break;
                                        //s2+=0.001;
                                        //r_xi = sum.fitTo(data,Save(),Range("cut"));
                                        count++;
                                        if(count == 1) break;
                                    }

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

                                    //cout << "Yield1: " << gaus1_yield_xi << endl;
                                    //cout << "Yield2: " << gaus2_yield_xi << endl;

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

                                    //cout << "Yield (xi): " << Yield_xi << endl;
                                    //cout << "Fsig (xi): " << Fsig_xi << endl;
                                    //cout << "std (xi): "  << rms_true_xi  << endl;
                                    //cout << "mass (xi): " << mean_xi << endl;

                                    //cout << "covQual (xi)" << covQual << endl;
                                    //cout << "Signal Sig (xi)" << significance << endl;


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
                                    os << "Fsig: " << Fsig_xi;
                                    tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                    os.str(std::string());
                                    os << "Mean: " << std::setprecision(5) << mean_xi << " GeV" << std::setprecision(6);
                                    tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                    os.str(std::string());
                                    os << "#sigma :" << std::setprecision(2) << rms_true_xi << " GeV" << std::setprecision(6);
                                    tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                    os.str(std::string());
                                    os << "CovQual: " << covQual;
                                    tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                    os.str(std::string());
                                    osYield << "S: " << std::setprecision(2) << Yield_xi;
                                    tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
                                    osYield.str(std::string());
                                    osYield << "B: " << std::setprecision(2) << IntbackgroundE_xi;
                                    tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
                                    osYield.str(std::string());
                                    osYield << "S-B: " << std::setprecision(2) << Yield_xi - IntbackgroundE_xi;
                                    tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
                                    osYield.str(std::string());
                                    os.str(std::string());
                                    os << "Sig: " << significance;
                                    tex->DrawLatex(xpos,0.85,os.str().c_str());
                                    os.str(std::string());

                                    xpos = 0.20;
                                    ypos = 0.85;
                                    os << VvarBins_pt[a-1] << " < Pt < " << VvarBins_pt[a];
                                    tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                    os.str(std::string());
                                    if(PID == "Omega") os << "#Omega DCA < " << Vcas3dipsig[b-1];
                                    if(PID == "Xi") os << "#Xi DCA < " << Vcas3dipsig[b];
                                    tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                    os.str(std::string());
                                    os << "Bat DCA > " << VBat3dipsig[c-1];
                                    tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                    os.str(std::string());
                                    os << "#pi DCA > " << Vvtrkpi3dipsig[d-1];
                                    tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                    os.str(std::string());
                                    os << "proton DCA > " << Vvtrkp3dipsig[e-1];
                                    tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                    os.str(std::string());
                                    if(PID == "Omega") os << "#Omega DecayL > " << Vflightsig[f-1];
                                    if(PID == "Xi") os << "#Xi DecayL > " << Vflightsig[f-1];
                                    tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                    os.str(std::string());
                                    os << "#Lambda DecayL > " << Vdistancesig[i-1];
                                    tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                    os.str(std::string());


                                    if(covQual < 3)
                                    {
                                        //name << "OmegaDistributions/Omega_" << Vcas3dipsig[b] << "_" << VBat3dipsig[c] << "_" << Vvtrkpi3dipsig[d] << "_" << Vvtrkp3dipsig[e] << "_" << Vflightsig[f] << "_" << Vdistancesig[i] << "_" << a << "_v" << version << "REDO.pdf";
                                        if(PID == "Omega") name << "Distributions/OmegaDistributions/Omega_" << mask << "_" << a << "_REDO_v" << version << ".pdf";
                                        if(PID == "Xi") name << "Distributions/XiDistributions/Xi_" << mask << "_" << a << "_REDO_v" << version << ".pdf";
                                    }
                                    else
                                    {
                                        //SigSig[mask] = significance;
                                        File << mask << "\t" << significance << "\n";
                                        //name << "OmegaDistributions/Omega_" << Vcas3dipsig[b] << "_" << VBat3dipsig[c] << "_" << Vvtrkpi3dipsig[d] << "_" << Vvtrkp3dipsig[e] << "_" << Vflightsig[f] << "_" << Vdistancesig[i] << "_" << a << "_v" << version << "REDO.pdf";
                                        if(PID == "Omega") name << "Distributions/OmegaDistributions/Omega_" << mask << "_" << a << "_REDO_v" << version << ".pdf";
                                        if(PID == "Xi") name << "Distributions/XiDistributions/Xi_" << mask << "_" << a << "_REDO_v" << version << ".pdf";
                                        if( remove(name.str().c_str()) == 0) cout << "File deleted";
                                        name.str(std::string());
                                        if(PID == "Omega") name << "Distributions/OmegaDistributions/Omega_" << mask << "_" << a << "_v" << version << ".pdf";
                                        if(PID == "Xi") name << "Distributions/XiDistributions/Xi_" << mask << "_" << a << "_v" << version << ".pdf";
                                    }

                                    c2->Print(name.str().c_str());

                                    delete t1;
                                    delete t2;
                                    delete h1;
                                }
                            }
                        }
                    }
                }
            }
            File.close();
            //double value = 0;
            //int pairvalue = 0;
            //double maxVal = 0;
            //int maxpair = 0;
            //for(std::map<int,double>::iterator it=SigSig.begin(); it!=SigSig.end(); ++it){
            //value = it->second;
            //if(value > maxVal){ 
            //maxVal = value;
            //maxpair = it->first;
            //}
            //}
            //cout << "Max Sig: " << maxVal << endl;
            //cout << "Max Key: " << maxpair << endl;
            //std::ofstream myfile;
            //std::ostringstream filename;
            //filename << PID << "_AllMaxKeys_" << pTBin << "_v" << version << ".txt";
            //myfile.open(filename.str().c_str());
            //myfile << VvarBins_pt[a-1] << " < Pt < " << VvarBins_pt[a] << "\n\n";
            //myfile << "Max Sig: " << maxVal << "\n";
            //myfile << "Max Key: " << maxpair << "\n";
            //myfile << "Other Keys with the same Significance:" << "\n";
            //for(std::map<int,double>::iterator it=SigSig.begin(); it!=SigSig.end(); ++it){
            //double values = it->second;
            //int key = 0;
            //if(values == maxVal){ 
            //key = it->first;
            //myfile << key << "\n";
            //}
            //}
        }
    }
    else
    {
        for(int a=1; a<VvarBins_pt.size(); a++)
            //for(int a=1; a<2; a++)
        {
            //if(a == 3) continue;
            pTBin = a;
            std::ostringstream oss;
            oss << "MaskAndSig_" << pTBin << "_" << PID << "_v" << version <<  ".txt";
            File.open(oss.str().c_str(),ios_base::app);
            oss.str(std::string());
            File << "\n";
            h->GetAxis(1)->SetRange(a,a);
            //for(int b=1; b<h->GetAxis(3)->GetNbins()+1; b++)
            //for(int b=3; b<7; b++)
            for(int b=1; b<2; b++)
            {
                h->GetAxis(3)->SetRange(b,h->GetAxis(3)->GetNbins()+1);
                h->GetAxis(4)->SetRange(b,h->GetAxis(4)->GetNbins()+1);
                h->GetAxis(5)->SetRange(b,h->GetAxis(5)->GetNbins()+1);
                h->GetAxis(6)->SetRange(b,h->GetAxis(6)->GetNbins()+1);
                //for(int d=5; d<h->GetAxis(7)->GetNbins()+2; d++)
                for(int d=2; d<3; d++)
                //for(int d=vdls.size(); d<vdls.size()+1; d++)
                {
                    h->GetAxis(7)->SetRange(d,h->GetAxis(7)->GetNbins()+1);
                    //for(int e=h->GetAxis(8)->GetNbins()-2; e<h->GetAxis(8)->GetNbins()+2; e++)
                    for(int e=4; e<5;e++)
                    //for(int e=6; e<10;e++)
                    {
                        h->GetAxis(8)->SetRange(e,h->GetAxis(8)->GetNbins()+1);
                        struct stat buffer;
                        std::ostringstream name;
                        std::string filename;
                        int mask = b*1e2 + d*1e1 + e;
                        if(PID == "Kshort") name << "Distributions/KshortDistributions/Kshort_" << mask << "_" << a << "_v" << version << ".pdf";
                        if(PID == "Lambda") name << "Distributions/LambdaDistributions/Lambda_" << mask << "_" << a << "_v" << version << ".pdf";
                        //filename = name.str();
                        //if(stat(name.str().c_str(),&buffer) == 0) continue;
                        name.str(std::string());
                        TH1D* h1 = (TH1D*)h->Projection(0);
                        //h1->GetYaxis()->SetRangeUser(0,1000);
                        gStyle->SetOptTitle(kFALSE);

                        double s1;
                        double s2;

                        RooRealVar x("x","mass",0,1);
                        //RooRealVar mean;
                        RooPlot* xframe_;
                        if(PID == "Kshort")
                        {
                            x.setRange(0.43,0.565);
                            xframe_ = x.frame(270);
                            //mean.setVal(0.50);
                            //mean.setRange(0.49,0.51); //Does this need to be changed to a narrower range?
                            s1=0.01;
                            s2=0.01;
                        }
                        else if(PID == "Lambda")
                        {
                            x.setRange(1.08,1.16);
                            xframe_ = x.frame(160);
                            //mean = MakeMeanVar("mean","mean",1.115,1.11,1.12);
                            s1=0.005;
                            s2=0.005;
                        }
                        if(PID == "Kshort") xframe_->GetXaxis()->SetTitle("#pi^{+}#pi^{-} Invariant mass (GeV)");
                        else if(PID == "Lambda") xframe_->GetXaxis()->SetTitle("#pi p Invariant mass (GeV)");
                        xframe_->GetYaxis()->SetTitle("Candidates / 0.0005 GeV");
                        RooDataHist data("data","dataset",x,h1);
                        data.plotOn(xframe_,Name("data"));
                        //RooRealVar mean("mean","mean",1.115,1.11,1.12);
                        RooRealVar mean("mean","mean",0.5,0.49,0.51);
                        RooRealVar sigma1("sigma1","sigma1",s1,0.001,0.01);
                        RooRealVar sigma2("sigma2","sigma2",s2,0.0001,0.01);
                        RooRealVar sig1("sig1","signal1",1e5,0,100000000);
                        RooRealVar sig2("sig2","signal2",1e5,0,100000000);
                        RooRealVar qsig("qsig","qsig",1e6,0,100000000);
                        RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
                        RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
                        RooRealVar ap("ap","ap",0,-5,5);
                        RooRealVar bp("bp","bp",0.1,-5,5);
                        RooRealVar cp("cp","cp",-0.1,-5,5);
                        RooRealVar dp("dp","dp",0.1,-5,5);
                        //RooRealVar ap("a","a",10,-100000,100000);
                        //RooRealVar bp("b","b",10,-100000,100000);
                        //RooRealVar cp("cp","cp",10,-100000,100000);
                        //RooRealVar dp("d","d",10,-100000,100000);
                        RooRealVar alpha("alpha","alpha",0.5,0,10);
                        RooChebychev background("background","background",x,RooArgList(ap,bp,cp,dp));
                        //RooPolynomial background("background","background",x,RooArgList(ap,bp,cp,dp));
                        RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,background),RooArgList(sig1,sig2,qsig));

                        if(PID == "Kshort") x.setRange("cut",0.44,0.56);
                        if(PID == "Lambda")
                        {
                            if(a >= 2 && a < 9)
                            {
                                x.setRange("cut",1.0805,1.1595);
                            }
                            else if( a == 1)
                            {
                                x.setRange("cut",1.0905,1.1495);
                            }
                            else if(a == 9)
                            {
                                x.setRange("cut",1.095,1.15);
                            }
                            else
                            {
                                x.setRange("cut",1.095, 1.155);
                            }
                        }

                        //RooFitResult* r_xi = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));
                        //mean.setConstant(kTRUE);
                        RooFitResult* r_xi = sum.fitTo(data,Save(),Range("cut"),PrintLevel(-1));
                        //mean.setConstant(kFALSE);
                        //r_xi = sum.fitTo(data,Save(),Range("cut"));
                        //RooFitResult* r_xi = sum.fitTo(data,Save(),Range("cut"));

                        int count = 0;
                        while(r_xi->covQual() < 3 ){
                            s1+=0.004;
                            s2+=0.003;
                            r_xi = sum.fitTo(data,Save(),Range("cut"));
                            //if(r_xi->covQual() == 3) break;
                            //s2+=0.001;
                            //r_xi = sum.fitTo(data,Save(),Range("cut"));
                            count++;
                            if(count == 1) break;
                        }

                        //x.setRange("cut",1.62,1.73);

                        //RooFitResult* r_xi = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));
                        //RooFitResult* r_xi = sum.fitTo(data,Save(),Range("cut"));
                        //RooChi2Var chi2_xiVar("chi2_xi","chi2",sum,data);

                        double covQual = r_xi->covQual();
                        double edm = r_xi->edm();
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

                        //cout << "Yield1: " << gaus1_yield_xi << endl;
                        //cout << "Yield2: " << gaus2_yield_xi << endl;

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

                        //cout << "Yield (xi): " << Yield_xi << endl;
                        //cout << "Fsig (xi): " << Fsig_xi << endl;
                        //cout << "std (xi): "  << rms_true_xi  << endl;
                        //cout << "mass (xi): " << mean_xi << endl;

                        //cout << "covQual (xi)" << covQual << endl;
                        //cout << "Signal Sig (xi)" << significance << endl;


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
                        os << "Fsig: " << Fsig_xi;
                        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                        os.str(std::string());
                        os << "Count: " << count;
                        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                        os.str(std::string());
                        os << "Mean: " << std::setprecision(5) << mean_xi << " GeV" << std::setprecision(6);
                        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                        os.str(std::string());
                        os << "#sigma :" << std::setprecision(2) << rms_true_xi << " GeV" << std::setprecision(6);
                        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                        os.str(std::string());
                        os << "CovQual: " << covQual;
                        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                        os.str(std::string());
                        osYield << "Yield: " << std::setprecision(2) << Yield_xi;
                        tex->DrawLatex(xpos,ypos-=increment,osYield.str().c_str());
                        osYield.str(std::string());
                        os.str(std::string());
                        os << "Sig: " << significance;
                        tex->DrawLatex(xpos,0.85,os.str().c_str());
                        os.str(std::string());
                        os << "EDM: " << edm;
                        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                        os.str(std::string());

                        xpos = 0.20;
                        ypos = 0.85;
                        os << VvarBins_pt[a-1] << " < Pt < " << VvarBins_pt[a];
                        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                        os.str(std::string());
                        os << "dz1 > " << vdz1[b-1];
                        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                        os.str(std::string());
                        os << "dz2 > " << vdz2[b-1];
                        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                        os.str(std::string());
                        os << "dxy1 > " << vdxy1[b-1];
                        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                        os.str(std::string());
                        os << "dxy2 > " << vdxy2[b-1];
                        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                        os.str(std::string());
                        os << "DecayL > " << vdls[d-1];
                        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                        os.str(std::string());
                        os << "Cos#theta > " << vagl[e-1];
                        tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                        os.str(std::string());


                        if(covQual < 3)
                        {
                            //name << "OmegaDistributions/Omega_" << Vcas3dipsig[b] << "_" << VBat3dipsig[c] << "_" << Vvtrkpi3dipsig[d] << "_" << Vvtrkp3dipsig[e] << "_" << Vflightsig[f] << "_" << Vdistancesig[i] << "_" << a << "_v" << version << "REDO.pdf";
                            if(PID == "Kshort") name << "Distributions/KshortDistributions/Kshort_" << mask << "_" << a << "_REDO_v" << version << ".pdf";
                            if(PID == "Lambda") name << "Distributions/LambdaDistributions/Lambda_" << mask << "_" << a << "_REDO_v" << version << ".pdf";
                        }
                        else
                        {
                            //SigSig[mask] = significance;
                            File << mask << "\t" << significance << "\n";
                            //name << "OmegaDistributions/Omega_" << Vcas3dipsig[b] << "_" << VBat3dipsig[c] << "_" << Vvtrkpi3dipsig[d] << "_" << Vvtrkp3dipsig[e] << "_" << Vflightsig[f] << "_" << Vdistancesig[i] << "_" << a << "_v" << version << "REDO.pdf";
                            if(PID == "Kshort") name << "Distributions/KshortDistributions/Kshort_" << mask << "_" << a << "_REDO_v" << version << ".pdf";
                            if(PID == "Lambda") name << "Distributions/LambdaDistributions/Lambda_" << mask << "_" << a << "_REDO_v" << version << ".pdf";
                            if( remove(name.str().c_str()) == 0) cout << "File deleted";
                            name.str(std::string());
                            if(PID == "Kshort") name << "Distributions/KshortDistributions/Kshort_" << mask << "_" << a << "_v" << version << ".pdf";
                            if(PID == "Lambda") name << "Distributions/LambdaDistributions/Lambda_" << mask << "_" << a << "_v" << version << ".pdf";
                        }

                        c2->Print(name.str().c_str());

                        delete t1;
                        delete t2;
                        delete h1;
                    }
                }
            }
            File.close();
            //double value = 0;
            //int pairvalue = 0;
            //double maxVal = 0;
            //int maxpair = 0;
            //for(std::map<int,double>::iterator it=SigSig.begin(); it!=SigSig.end(); ++it){
            //value = it->second;
            //if(value > maxVal){ 
            //maxVal = value;
            //maxpair = it->first;
            //}
            //}
            //cout << "Max Sig: " << maxVal << endl;
            //cout << "Max Key: " << maxpair << endl;
            //std::ofstream myfile;
            //std::ostringstream filename;
            //filename << PID << "_AllMaxKeys_" << pTBin << "_v" << version << ".txt";
            //myfile.open(filename.str().c_str());
            //myfile << VvarBins_pt[a-1] << " < Pt < " << VvarBins_pt[a] << "\n\n";
            //myfile << "Max Sig: " << maxVal << "\n";
            //myfile << "Max Key: " << maxpair << "\n";
            //myfile << "Other Keys with the same Significance:" << "\n";
            //for(std::map<int,double>::iterator it=SigSig.begin(); it!=SigSig.end(); ++it){
            //double values = it->second;
            //int key = 0;
            //if(values == maxVal){ 
            //key = it->first;
            //myfile << key << "\n";
            //}
            //}
        }
    }
}

