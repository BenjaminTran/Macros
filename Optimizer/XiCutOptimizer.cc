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

void WriteVector(std::vector<float> cutVector, std::ofstream &myfile, std::string name){
    myfile << name << ": ";
    for(std::vector<float>::iterator it=cutVector.begin(); it!=cutVector.end(); it++){
        myfile << *it << ", ";
    }
    myfile << "\n";
}

TTreeReaderValue<float> SetValues(TTreeReader& reader, std::string branchName ){
    TTreeReaderValue<float> ReaderValue(reader, branchName.c_str());
    return ReaderValue;
}

bool SparseCreatorOmXi(std::string name = "Test", std::string PID = "Omega")
{
    //Initializers
    bool Cut = false;
    TH1::SetDefaultSumw2();
    std::ofstream myfile;
    std::string txtname = name + ".txt";
    myfile.open(txtname.c_str());

    //Cut parameters to be varied. The commented elements are Hong's cuts to make a plot of them remember to change the numparam value as well
    std::vector<float> VvarBins_pt            = {1.0,1.5,1.8,2.2,2.8,3.6,4.6,6.0,7.2,10.0};
    std::vector<float> VvarBins_eta           = {-2.4,-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,2.0,2.4};

    std::vector<float> Vcas3dipsig;
    std::vector<float> VBat3dipsig;
    std::vector<float> Vvtrkpi3dipsig;
    std::vector<float> Vvtrkp3dipsig;
    std::vector<float> Vflightsig;
    std::vector<float> Vdistancesig;

    if(PID == "Omega")
    {
        Vcas3dipsig.insert(Vcas3dipsig.begin(),{2.0, 3.0 , 3.4  , 3.6 , 3.8  , 4.0 , 4.2 , 4.5 , 5.0});
        VBat3dipsig.insert(VBat3dipsig.begin(),{3.0 , 3.4  , 3.8 , 4.0  , 4.2 , 4.4 , 4.8 , 5.0 , 5.2});
        Vvtrkpi3dipsig.insert(Vvtrkpi3dipsig.begin(),{3.0 , 3.2  , 3.5 , 3.8  , 4.0 , 4.5 , 5.0 , 5.5 , 6.0});
        Vvtrkp3dipsig.insert(Vvtrkp3dipsig.begin(),{2.0 , 2.2  , 2.5 , 2.8  , 4.0 , 4.5 , 5.0 , 5.5 , 6.0});
        Vflightsig.insert(Vflightsig.begin(),{1.5 , 2.0  , 2.5 , 3.0  , 3.5 , 4.0 , 4.5 , 5.0 , 5.5});
        Vdistancesig.insert(Vdistancesig.begin(),{10. , 10.5 , 11  , 11.5 , 12. , 13  , 14  , 15  , 16  });
    }
    else
    {
        Vcas3dipsig.insert(Vcas3dipsig.begin(),{2.0, 3.0 , 3.4  , 3.6 , 3.8  , 4.0 , 4.2 , 4.5 , 5.0});
        VBat3dipsig.insert(VBat3dipsig.begin(),{3.0 , 3.4  , 3.8 , 4.0  , 4.2 , 4.4 , 4.8 , 5.0 , 5.2});
        Vvtrkpi3dipsig.insert(Vvtrkpi3dipsig.begin(),{3.0 , 3.2  , 3.5 , 3.8  , 4.0 , 4.5 , 5.0 , 5.5 , 6.0});
        Vvtrkp3dipsig.insert(Vvtrkp3dipsig.begin(),{2.0 , 2.2  , 2.5 , 2.8  , 4.0 , 4.5 , 5.0 , 5.5 , 6.0});
        Vflightsig.insert(Vflightsig.begin(),{1.5 , 2.0  , 2.5 , 3.0  , 3.5 , 4.0 , 4.5 , 5.0 , 5.5});
        Vdistancesig.insert(Vdistancesig.begin(),{10. , 10.5 , 11  , 11.5 , 12. , 13  , 14  , 15  , 16  });
    }

    //std::vector<float> Vcas3dipsig    = {2.4 , 2.5 , 2.7 , 3.0 , 3.2  , 3.4 , 3.6};//3.0 };// , 2.5};
    //std::vector<float> VBat3dipsig    = {4.0 , 4.2 , 4.4 , 4.6 , 4.8 , 5.0 , 5.2};//4.0 };// , 5.0};
    //std::vector<float> Vvtrkpi3dipsig = {3.0 , 3.2 , 3.4  , 3.6 , 3.8 , 4.0 ,4.2};//3.0 };// , 4.0};
    //std::vector<float> Vvtrkp3dipsig  = {2.0 , 2.2 , 2.4  , 2.6 , 2.8 , 3.0, 3.2};//2.0 };// , 3.0};
    //std::vector<float> Vflightsig     = {2.0 , 2.2 , 2.4  , 2.6 , 2.8 , 3.0 ,3.2};//2.0 };// , 3.0};
    //std::vector<float> Vdistancesig   = {5. , 6. , 7.0 , 8.0 , 10.0 , 12.0 ,14.0};//10.0};// , 12.0};



    double misIDMass = 0.015;
    double rapidity = 1.0;
    double etacut = 2.4;
    int multHigh_ = 250;
    int numparam = Vcas3dipsig.size();
    const int nDim = 9;
    // mass, pt, eta, then cuts
    //These are for eta cut
    std::vector<int> nBins = {150,9,10,static_cast<int>(Vcas3dipsig.size()),static_cast<int>(VBat3dipsig.size()),static_cast<int>(Vvtrkpi3dipsig.size()),static_cast<int>(Vvtrkp3dipsig.size()),static_cast<int>(Vflightsig.size()),static_cast<int>(Vdistancesig.size())};
    std::vector<double> minbins;
    std::vector<double> maxbins;

    //Omega
    if(PID == "Omega")
    {
        minbins.insert(minbins.begin(),{1.6 ,1.0,-2.4,0                 ,3.5               ,2.5                  ,1.5                 ,1.5              ,5.0});
        maxbins.insert(maxbins.begin(),{1.75,10 ,2.4 ,Vcas3dipsig.back(),VBat3dipsig.back(),Vvtrkpi3dipsig.back(),Vvtrkp3dipsig.back(),Vflightsig.back(),Vdistancesig.back()});
    }
    else
    {
        minbins.insert(minbins.begin(),{1.25 , 1.0 , -2.4 , 0                  , VBat3dipsig.front() , Vvtrkpi3dipsig.front() , Vvtrkp3dipsig.front() , Vflightsig.front() , Vdistancesig.front()});
        maxbins.insert(maxbins.begin(),{1.40 , 10  , 2.4  , Vcas3dipsig.back() , 99999               , 99999                  , 99999                 , 99999              , 99999});
    }
    //Xi

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

    if(PID == "Omega")
    {
        WriteVector(Vcas3dipsig,myfile,"Om DCA");
        WriteVector(VBat3dipsig,myfile,"Bat DCA");
        WriteVector(Vvtrkpi3dipsig,myfile,"Pi DCA");
        WriteVector(Vvtrkp3dipsig,myfile,"Pr DCA");
        WriteVector(Vflightsig,myfile,"Om DLS");
        WriteVector(Vdistancesig,myfile,"La DLS");
    }
    else
    {
        WriteVector(Vcas3dipsig,myfile,"Xi DCA");
        WriteVector(VBat3dipsig,myfile,"Bat DCA");
        WriteVector(Vvtrkpi3dipsig,myfile,"Pi DCA");
        WriteVector(Vvtrkp3dipsig,myfile,"Pr DCA");
        WriteVector(Vflightsig,myfile,"Xi DLS");
        WriteVector(Vdistancesig,myfile,"La DLS");
    }

    THnSparseF* h = new THnSparseF("hSparse","Cuts",nDim,&nBins[0],&minbins[0],&maxbins[0]);

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



    //Containers


    //Tree setup
    //TFile* f1=TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/V0CasTree_09_14_17.root");
    TFile* f1=TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/V0CasTreePbPb_Production_10_02_17.root");

    TTreeReader reader("OmTreeProducerRapidityPbPb/OmTree",f1);
    //TTreeReader reader("XiTreeProducerRapidityPbPb/XiTree",f1);

    //Omega
    TTreeReaderValue<float> cas3dipsig;
    TTreeReaderValue<float> rap;
    TTreeReaderValue<float> eta;
    TTreeReaderValue<float> pt;
    TTreeReaderValue<float> mass;
    TTreeReaderValue<float> distancesig;
    TTreeReaderValue<float> flightsig;
    TTreeReaderValue<float> vtrkp3dipsig;
    TTreeReaderValue<float> vtrkpi3dipsig;
    TTreeReaderValue<float> bat3dipsig;

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
    else
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


    reader.SetEntriesRange(0,10);
    while(reader.Next())
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

        double value[nDim] = {*mass,*pt,*eta,*cas3dipsig,*bat3dipsig,*vtrkpi3dipsig,*vtrkp3dipsig,*flightsig,*distancesig};
        h->Fill(value);
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

void readSparseOmega(std::string name)
{
    TFile* f = new TFile(name.c_str());
    bool DCAExt = false;
    int version = 1;
    //TFile* f = new TFile("SmallBins.root");
    std::vector<double> VvarBins_pt            = {1.0,1.5,1.8,2.2,2.8,3.6,4.6,6.0,7.2,10.0};
    std::vector<double> VvarBins_eta           = {-2.4,-2.0,-1.5,-1.0,-0.5,0,0.5,1.0,1.5,2.0,2.4};


    std::vector<float> Vcas3dipsig    = {2.0, 3.0 , 3.4  , 3.6 , 3.8  , 4.0 , 4.2 , 4.5 , 5.0};
    std::vector<float> VBat3dipsig    = {3.0 , 3.4  , 3.8 , 4.0  , 4.2 , 4.4 , 4.8 , 5.0 , 5.2};
    std::vector<float> Vvtrkpi3dipsig = {3.0 , 3.2  , 3.5 , 3.8  , 4.0 , 4.5 , 5.0 , 5.5 , 6.0};
    std::vector<float> Vvtrkp3dipsig  = {2.0 , 2.2  , 2.5 , 2.8  , 4.0 , 4.5 , 5.0 , 5.5 , 6.0};
    std::vector<float> Vflightsig     = {1.5 , 2.0  , 2.5 , 3.0  , 3.5 , 4.0 , 4.5 , 5.0 , 5.5};
    std::vector<float> Vdistancesig   = {10. , 10.5 , 11  , 11.5 , 12. , 13  , 14  , 15  , 16  };
    Vcas3dipsig.insert(Vcas3dipsig.begin(),0);
    VBat3dipsig.push_back(99999);
    Vvtrkpi3dipsig.push_back(99999);
    Vvtrkp3dipsig.push_back(99999);
    Vflightsig.push_back(99999);
    Vdistancesig.push_back(99999);

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
    for(int a=2; a<3; a++)
    {
        pTBin = a;
        h->GetAxis(1)->SetRange(a,a);
        //for(int b=1; b<Vcas3dipsig.size(); b++)
        for(int b=3; b<4; b++)
        {
            h->GetAxis(3)->SetRange(1,b+1);
            //for(int c=1; c<VBat3dipsig.size(); c++)
            for(int c=1; c<2; c++)
            {
                h->GetAxis(4)->SetRange(c,10);
                //for(int d=1; d<Vvtrkpi3dipsig.size(); d++)
                for(int d=1; d<2; d++)
                {
                    h->GetAxis(5)->SetRange(d,10);
                    //for(int e=1; e<Vvtrkp3dipsig.size();e++)
                    for(int e=1; e<2;e++)
                    {
                        h->GetAxis(6)->SetRange(e,10);
                        //for(int f=1; f<Vflightsig.size(); f++)
                        for(int f=1; f<2; f++)
                        {
                            h->GetAxis(7)->SetRange(f,10);
                            //for(int i=1; i<Vdistancesig.size(); i++)
                            for(int i=1; i<2; i++)
                            {
                                struct stat buffer;
                                std::ostringstream name;
                                std::string filename;
                                int mask = b*1e5 + c*1e4 + d*1e3 + e*1e2 + f*1e1 + i;
                                //std::string mask = b*1e5 + c*1e4 + d*1e3 + e*1e2 + f*1e1 + i;
                                if(DCAExt)
                                {
                                    name << "OmegaDistributions/Omega_" << mask << "_" << a << ".pdf";
                                }
                                else
                                {
                                    name << "OmegaDistributions/Omega_" << mask << "_" << a << "_v" << version << ".pdf";
                                    //name << "OmegaDistributions/Omega_" << Vcas3dipsig[b-1] << "_" << VBat3dipsig[c-1] << "_" << Vvtrkpi3dipsig[d-1] << "_" << Vvtrkp3dipsig[e-1] << "_" << Vflightsig[f-1] << "_" << Vdistancesig[i-1] << "_" << a << "_v" << version << "2.pdf";
                                    //filename = name.str();
                                }
                                //if(stat(name.str().c_str(),&buffer) == 0) continue;
                                name.str(std::string());
                                h->GetAxis(8)->SetRange(i,8);
                                TH1D* h1 = (TH1D*)h->Projection(0);
                                h1->GetYaxis()->SetRangeUser(0,1000);
                                gStyle->SetOptTitle(kFALSE);

                                TLatex* tex = new TLatex();
                                tex->SetNDC();
                                tex->SetTextFont(42);
                                tex->SetTextSize(0.025);
                                //tex->SetTextAlign(10);

                                double s1=0.001;
                                double s2=0.001;

                                RooRealVar x("x","mass",1.60,1.75);
                                RooPlot* xframe_ = x.frame(150);
                                xframe_->GetXaxis()->SetTitle("#Lambda K Invariant mass (GeV)");
                                xframe_->GetYaxis()->SetTitle("Candidates / 0.002 GeV");
                                RooDataHist data("data","dataset",x,h1);
                                data.plotOn(xframe_,Name("data"));
                                RooRealVar mean("mean","mean",1.67,1.6,1.75);
                                RooRealVar sigma1("sigma1","sigma1",s1,0.001,0.04);
                                RooRealVar sigma2("sigma2","sigma2",s2,0.001,0.04);
                                RooRealVar sig1("sig1","signal1",2000,0,10000);
                                RooRealVar sig2("sig2","signal2",2000,0,10000);
                                RooRealVar qsig("qsig","qsig",2000,0,100000);
                                RooRealVar alpha("alpha","alpha",0.9,-1,9);
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

                                RooFitResult* r_xi = sum.fitTo(data,Save(),Range("cut"));
                                x.setRange("cut",1.60,1.75);

                                int count = 0;
                                while(r_xi->covQual() < 3 ){
                                    s1+=0.001;
                                    s2+=0.001;
                                    r_xi = sum.fitTo(data,Save(),Range("cut"));
                                    //if(r_xi->covQual() == 3) break;
                                    //s2+=0.001;
                                    //r_xi = sum.fitTo(data,Save(),Range("cut"));
                                    count++;
                                    if(count == 10) break;
                                }

                                //x.setRange("cut",1.62,1.73);

                                //RooFitResult* r_xi = sum.fitTo(data,Save(),Minos(kTRUE),Range("cut"));
                                //RooFitResult* r_xi = sum.fitTo(data,Save(),Range("cut"));
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
                                os << "#Omega DCA < " << Vcas3dipsig[b];
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());
                                os << "Kaon DCA > " << VBat3dipsig[c-1];
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());
                                os << "#pi DCA > " << Vvtrkpi3dipsig[d-1];
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());
                                os << "proton DCA > " << Vvtrkp3dipsig[e-1];
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());
                                os << "#Omega DecayL > " << Vflightsig[f-1];
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());
                                os << "#Lambda DecayL > " << Vdistancesig[i-1];
                                tex->DrawLatex(xpos,ypos-=increment,os.str().c_str());
                                os.str(std::string());


                                if(DCAExt)
                                {
                                    if(covQual < 3)
                                        name << "OmegaDistributions/Omega_DCAExt_" << mask << "_" << a << "_REDO" << ".pdf";
                                    else
                                    {
                                        SigSig[mask] = significance;
                                        name << "OmegaDistributions/Omega_DCAExt_" << mask << "_" << a << ".pdf";
                                    }
                                }
                                else{
                                    if(covQual < 3)
                                    {
                                    //name << "OmegaDistributions/Omega_" << Vcas3dipsig[b] << "_" << VBat3dipsig[c] << "_" << Vvtrkpi3dipsig[d] << "_" << Vvtrkp3dipsig[e] << "_" << Vflightsig[f] << "_" << Vdistancesig[i] << "_" << a << "_v" << version << "REDO.pdf";
                                        name << "OmegaDistributions/Omega_" << mask << "_" << a << "_REDO_v" << version << ".pdf";
                                    }
                                    else
                                    {
                                        SigSig[mask] = significance;
                                        //name << "OmegaDistributions/Omega_" << Vcas3dipsig[b] << "_" << VBat3dipsig[c] << "_" << Vvtrkpi3dipsig[d] << "_" << Vvtrkp3dipsig[e] << "_" << Vflightsig[f] << "_" << Vdistancesig[i] << "_" << a << "_v" << version << "REDO.pdf";
                                        name << "OmegaDistributions/Omega_" << mask << "_" << a << "_REDO_v" << version << ".pdf";
                                        if( remove(name.str().c_str()) == 0) cout << "File deleted";
                                        name.str(std::string());
                                        name << "OmegaDistributions/Omega_" << mask << "_" << a << "_v" << version << ".pdf";
                                    }
                                }

                                c2->Print(name.str().c_str());
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
    cout << "Max Key: " << maxpair << endl;
    std::ofstream myfile;
    std::ostringstream filename;
    if(DCAExt)
    {
        filename << "AllMaxKeys_DCAExt_" << pTBin << ".txt";
    }
    else{
        filename << "AllMaxKeys_" << pTBin << "_v" << version << ".txt";
    }
    myfile.open(filename.str().c_str());
    myfile << "Max Sig: " << maxVal << "\n";
    myfile << "Max Key: " << maxpair << "\n";
    myfile << "Other Keys with the same Significance:" << "\n";
    for(std::map<int,double>::iterator it=SigSig.begin(); it!=SigSig.end(); ++it){
        double values = it->second;
        int key = 0;
        if(values == maxVal){ 
            key = it->first;
            myfile << key << "\n";
        }
    }
}

