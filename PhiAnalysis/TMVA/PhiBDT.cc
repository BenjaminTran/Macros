#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

void PhiBDT()
{
    TMVA::Tools::Instance();
    TString outfileName("PhiBDT.root");

    //TFile* input = TFile::Open("/Users/btran/research/RootFiles/Phi/TMVA_SignalBackground_v1.root");
    TFile* input = TFile::Open("/Users/btran/research/RootFiles/Phi/TMVA_SignalBackground_v2.root");
    //TFile* input = TFile::Open("/Users/btran/research/RootFiles/Phi/PhiGenMatch_v4_1.root");
    //TFile* input = TFile::Open("/Users/btran/research/RootFiles/Phi/PhiGenMatch_v4_30K_throwaway.root");
    TFile* outputFile = TFile::Open(outfileName,"RECREATE");

    TTree* signalTree = (TTree*)input->Get("PhiGenMatch/SignalTree");
    TTree* backgroundTree = (TTree*)input->Get("PhiGenMatch/BackgroundTree");

    //TMVA::Factory *factory = new TMVA::Factory("Phi_BDT", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    TMVA::Factory *factory = new TMVA::Factory("Phi_BDT", outputFile, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification" );

    TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");

    dataloader->AddVariable("fabs(pt_1/ptError_1)", "Relative Pt Error 1", 'F');
    dataloader->AddVariable("fabs(dz_1/dzError_1)", "DCA z1", 'F');
    dataloader->AddVariable("fabs(dxy_1/dxyError_1)", "DCA xy1", 'F');
    dataloader->AddVariable("fabs(rapidity_1)", 'F');
    dataloader->AddVariable("nhits_1", 'I');
    dataloader->AddVariable("dedx_1", 'F');

    dataloader->AddVariable("fabs(pt_2/ptError_2)", "Relative Pt Error 2", 'F');
    dataloader->AddVariable("fabs(dz_2/dzError_2)", "DCA z2", 'F');
    dataloader->AddVariable("fabs(dxy_2/dxyError_2)", "DCA xy2", 'F');
    dataloader->AddVariable("fabs(rapidity_2)", 'F');
    dataloader->AddVariable("nhits_2", 'I');
    dataloader->AddVariable("dedx_2", 'F');


    dataloader->AddSpectator("mass", 'F');
    //dataloader->AddSpectator("chi2norm_1", 'F');
    //dataloader->AddSpectator("dedx_1", 'F');
    //dataloader->AddSpectator("eta_1", 'F');
    //dataloader->AddSpectator("phi_1", 'F');
    //dataloader->AddSpectator("vz_1", 'F');
    //dataloader->AddSpectator("vzFlip_1", 'F');
    //dataloader->AddSpectator("chi2norm_2", 'F');
    //dataloader->AddSpectator("dedx_2", 'F');
    //dataloader->AddSpectator("eta_2", 'F');
    //dataloader->AddSpectator("phi_2", 'F');
    //dataloader->AddSpectator("vz_2", 'F');
    //dataloader->AddSpectator("vzFlip_2", 'F');

    double signalWeight = 1.0;
    double backgroundWeight = 1.0;

    dataloader->AddSignalTree(signalTree, signalWeight);
    dataloader->AddBackgroundTree(backgroundTree, backgroundWeight);

    //dataloader->SetBackgroundWeightExpression("weight");

    TCut preCut = "";

    dataloader->PrepareTrainingAndTestTree(preCut, "nTrain_Signal=150000:nTrain_Background=200000:nTest_Signal=150000:nTest_Background=200000:SplitMode=Random:NormMode=NumEvents:V" );

    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
            "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );


    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    outputFile->Close();

    TMVA::TMVAGui( outfileName );
}
