//Includes
#include <TLatex.h>
#include <TStyle.h>
#include "TF1.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TLegend.h"
#include "TMathText.h"
#include "TImage.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"
#include <TString.h>
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TString.h"
#include "TGaxis.h"

#include <vector>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <sys/stat.h> //For file existance checking

void XiCutOptimizer(std::string name)
{
    //Initializers
    bool Cut = false;
    const int CONTSIZE = 10000;
    TH1::SetDefaultSumw2();

    //Cut parameters to be varied. The commented elements are Hong's cuts to make a plot of them remember to change the numparam value as well
    //float xi_xi3dipsig[]     = {8.0  , 8.5  , 9.0  , 9.5  , 10.0 , 10.5 , 11.0 , 11.5, 12.0, 12.5, 13.0, 13.5};//, 2.5};
    int Oparamindex          = 12; // For deciding which parameter that is not being varied to use
    std::vector<double> om_om3dipsig     = {3.0 };//3.0 };//, 2.5};
    std::vector<double> om_omKaon3dipsig = {4.0 };//4.0 };//, 5.0};
    std::vector<double> om_vtrkpi3dipsig = {3.0 };//3.0 };//, 4.0};
    std::vector<double> om_vtrkp3dipsig  = {2.0 };//2.0 };//, 3.0};
    std::vector<double> om_omflightsig   = {2.0 };//2.0 };//, 3.0};
    std::vector<double> om_distancesig   = {10.0};//10.0};//, 12.0};
    double misIDMass = 0.015;
    double rapidity = 1.0;
    int multHigh_ = 250;

    int numparam             = om_om3dipsig.size();

    //Containers
    //Tree values
    float tom_om3dipsig     [CONTSIZE];
    float tom_omKaon3dipsig   [CONTSIZE];
    float tom_vtrkpi3dipsig [CONTSIZE];
    float tom_vtrkp3dipsig  [CONTSIZE];
    float tom_omflightsig   [CONTSIZE];
    float tom_distancesig   [CONTSIZE];
    float tom_mass          [CONTSIZE];
    float tom_pt            [CONTSIZE];
    float tom_eta [CONTSIZE];
    float tom_rap [CONTSIZE];
    float tom_nTrkAcc [CONTSIZE];
    float tom_misIDMassLapi [CONTSIZE];
    float tom_misIDMasspiLa [CONTSIZE];
    int   tom_n;

    //Hist Containers
    TH2D* hom_om3dipsig     [numparam];
    TH2D* hom_omKaon3dipsig   [numparam];
    TH2D* hom_vtrkpi3dipsig [numparam];
    TH2D* hom_vtrkp3dipsig  [numparam];
    TH2D* hom_omflightsig   [numparam];
    TH2D* hom_distancesig   [numparam];

    //Tree setup
    TFile* f1=TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/V0CasTree_09_14_17.root");
    TTree* OmTree = (TTree*)f1->Get("OmTreeProducer/OmTree");

    OmTree->SetBranchAddress("nCand",              &tom_n);
    OmTree->SetBranchAddress("om3dipsig",      &tom_om3dipsig);
    OmTree->SetBranchAddress("omKaon3dipsig",    &tom_omKaon3dipsig);
    OmTree->SetBranchAddress("vtrkpi3dipsig",  &tom_vtrkpi3dipsig);
    OmTree->SetBranchAddress("vtrkp3dipsig",   &tom_vtrkp3dipsig);
    OmTree->SetBranchAddress("omflightsig",    &tom_omflightsig);
    OmTree->SetBranchAddress("distancesig",    &tom_distancesig);
    OmTree->SetBranchAddress("mass",           &tom_mass);
    OmTree->SetBranchAddress("pt",             &tom_pt);
    OmTree->SetBranchAddress("eta", &tom_eta);
    OmTree->SetBranchAddress("rapidity", &tom_rap);
    OmTree->SetBranchAddress("nTrkAcc", &tom_nTrkAcc);
    OmTree->SetBranchAddress("misIDMassLapi", &tom_misIDMassLapi);
    OmTree->SetBranchAddress("misIDMasspiLa", &tom_misIDMasspiLa);


    //Intialize Histograms
    TH2D* hom_NoCut = NULL;
    if(!Cut) hom_NoCut = new TH2D("hom_NoCut","NoCut",150,1.25,1.40,150,0,15);
    TH2D* hom_defaultcut = new TH2D("hom_Default","Default",150,1.60,1.75,400,0,40);
    for(int j=0; j<numparam; j++)
    {
        hom_om3dipsig[j]     = new TH2D(Form("hom_om3dipsig_%.1f",om_om3dipsig[j]),Form("hom_om3dipsig_%.1f",om_om3dipsig[j]),150,1.25,1.40,150,0,15);
        hom_omKaon3dipsig[j] = new TH2D(Form("hom_omKaon3dipsig_%.1f",om_omKaon3dipsig[j]),Form("hom_omKaon3dipsig_%.1f",om_omKaon3dipsig[j]),150,1.25,1.40,150,0,15);
        hom_vtrkpi3dipsig[j] = new TH2D(Form("hom_vtrkpi3dipsig_%.1f",om_vtrkpi3dipsig[j]),Form("hom_vtrkpi3dipsig_%.1f",om_vtrkpi3dipsig[j]),150,1.25,1.40,150,0,15);
        hom_vtrkp3dipsig[j]  = new TH2D(Form("hom_vtrkp3dipsig_%.1f",om_vtrkp3dipsig[j]),Form("hom_vtrkp3dipsig_%.1f",om_vtrkp3dipsig[j]),150,1.25,1.40,150,0,15);
        hom_omflightsig[j]   = new TH2D(Form("hom_omflightsig_%.1f",om_omflightsig[j]),Form("hom_omflightsig_%.1f",om_omflightsig[j]),150,1.25,1.40,150,0,15);
        hom_distancesig[j]   = new TH2D(Form("hom_distancesig_%.1f",om_distancesig[j]),Form("hom_distancesig_%.1f",om_distancesig[j]),150,1.25,1.40,150,0,15);
    }

    //Output file creation
    std::string filetype = ".root";
    name += filetype;
    struct stat buffer;
    if(stat(name.c_str(), &buffer) == 0)
    {
        cout << "File with this name already eomsts, please select a different name" << endl;
        return;
    }
    TFile out(name.c_str(),"RECREATE");
    cout << name.c_str() << " created!" << endl;

    //Loop over Tree
    int om_nEvent = OmTree->GetEntries();

    for(int i=0; i<om_nEvent; i++)
    {
        if(Cut)
        {
            OmTree->GetEntry(i);
            for( int j=0; j<numparam; j++)
            {
                //om3dipsig
                for(int k=0; k<tom_n; k++)
                {
                    if(tom_om3dipsig[k]     > om_om3dipsig[j])     continue;
                    if(tom_omKaon3dipsig[k]   < om_omKaon3dipsig[Oparamindex])   continue;
                    if(tom_vtrkpi3dipsig[k] < om_vtrkpi3dipsig[Oparamindex]) continue;
                    if(tom_vtrkp3dipsig[k]  < om_vtrkp3dipsig[Oparamindex])  continue;
                    if(tom_omflightsig[k]   < om_omflightsig[Oparamindex])   continue;
                    if(tom_distancesig[k]   < om_distancesig[Oparamindex])   continue;

                    hom_om3dipsig[j]->Fill(tom_mass[k],tom_pt[k]);
                }

                //ompi3dipsig
                for(int k=0; k<tom_n; k++)
                {
                    if(tom_om3dipsig[k] > om_om3dipsig[Oparamindex])         continue;
                    if(tom_omKaon3dipsig[k] < om_omKaon3dipsig[j])     continue;
                    if(tom_vtrkpi3dipsig[k] < om_vtrkpi3dipsig[Oparamindex]) continue;
                    if(tom_vtrkp3dipsig[k] < om_vtrkp3dipsig[Oparamindex])   continue;
                    if(tom_omflightsig[k] < om_omflightsig[Oparamindex])     continue;
                    if(tom_distancesig[k] < om_distancesig[Oparamindex])     continue;

                    hom_omKaon3dipsig[j]->Fill(tom_mass[k],tom_pt[k]);
                }

                //vtrkpi3dipsig
                for(int k=0; k<tom_n; k++)
                {
                    if(tom_om3dipsig[k] > om_om3dipsig[Oparamindex])         continue;
                    if(tom_omKaon3dipsig[k] < om_omKaon3dipsig[Oparamindex])     continue;
                    if(tom_vtrkpi3dipsig[k] < om_vtrkpi3dipsig[j]) continue;
                    if(tom_vtrkp3dipsig[k] < om_vtrkp3dipsig[Oparamindex])   continue;
                    if(tom_omflightsig[k] < om_omflightsig[Oparamindex])     continue;
                    if(tom_distancesig[k] < om_distancesig[Oparamindex])     continue;

                    hom_vtrkpi3dipsig[j]->Fill(tom_mass[k],tom_pt[k]);
                }

                //vtrkp3dipsig
                for(int k=0; k<tom_n; k++)
                {
                    if(tom_om3dipsig[k] > om_om3dipsig[Oparamindex])         continue;
                    if(tom_omKaon3dipsig[k] < om_omKaon3dipsig[Oparamindex])     continue;
                    if(tom_vtrkpi3dipsig[k] < om_vtrkpi3dipsig[Oparamindex]) continue;
                    if(tom_vtrkp3dipsig[k] < om_vtrkp3dipsig[j])   continue;
                    if(tom_omflightsig[k] < om_omflightsig[Oparamindex])     continue;
                    if(tom_distancesig[k] < om_distancesig[Oparamindex])     continue;

                    hom_vtrkp3dipsig[j]->Fill(tom_mass[k],tom_pt[k]);
                }

                //omflightsig
                for(int k=0; k<tom_n; k++)
                {
                    if(tom_om3dipsig[k] > om_om3dipsig[Oparamindex])         continue;
                    if(tom_omKaon3dipsig[k] < om_omKaon3dipsig[Oparamindex])     continue;
                    if(tom_vtrkpi3dipsig[k] < om_vtrkpi3dipsig[Oparamindex]) continue;
                    if(tom_vtrkp3dipsig[k] < om_vtrkp3dipsig[Oparamindex])   continue;
                    if(tom_omflightsig[k] < om_omflightsig[j])     continue;
                    if(tom_distancesig[k] < om_distancesig[Oparamindex])     continue;

                    hom_omflightsig[j]->Fill(tom_mass[k],tom_pt[k]);
                }

                //distancesig
                for(int k=0; k<tom_n; k++)
                {
                    if(tom_om3dipsig[k] > om_om3dipsig[Oparamindex])         continue;
                    if(tom_omKaon3dipsig[k] < om_omKaon3dipsig[Oparamindex])     continue;
                    if(tom_vtrkpi3dipsig[k] < om_vtrkpi3dipsig[Oparamindex]) continue;
                    if(tom_vtrkp3dipsig[k] < om_vtrkp3dipsig[Oparamindex])   continue;
                    if(tom_omflightsig[k] < om_omflightsig[Oparamindex])     continue;
                    if(tom_distancesig[k] < om_distancesig[j])     continue;

                    hom_distancesig[j]->Fill(tom_mass[k],tom_pt[k]);
                }
            }
        }
        else
        {
            OmTree->GetEntry(i);
            for(int k=0; k<tom_n; k++)
            {
                //hom_NoCut->Fill(tom_mass[k],tom_pt[k]);
                if(tom_nTrkAcc[k] > multHigh_) continue;
                if(tom_om3dipsig[k] > om_om3dipsig[0])         continue;
                if(tom_omKaon3dipsig[k] < om_omKaon3dipsig[0])     continue;
                if(tom_vtrkpi3dipsig[k] < om_vtrkpi3dipsig[0]) continue;
                if(tom_vtrkp3dipsig[k] < om_vtrkp3dipsig[0])   continue;
                if(tom_omflightsig[k] < om_omflightsig[0])     continue;
                if(tom_distancesig[k] < om_distancesig[0])     continue;
                if(std::abs(tom_rap[k]) > rapidity) continue;
                if(std::abs(tom_misIDMasspiLa[k]) < misIDMass) continue;
                if(std::abs(tom_misIDMassLapi[k]) < misIDMass) continue;

                hom_defaultcut->Fill(tom_mass[k],tom_pt[k]);
            }
        }
    }

    //Write histograms to root file
    if(Cut)
    {
        for(int j=0; j<numparam; j++)
        {
            hom_om3dipsig[j]->Write();
        }
        for(int j=0; j<numparam; j++)
        {
            hom_omKaon3dipsig[j]->Write();
        }
        for(int j=0; j<numparam; j++)
        {
            hom_vtrkpi3dipsig[j]->Write();
        }
        for(int j=0; j<numparam; j++)
        {
            hom_vtrkp3dipsig[j]->Write();
        }
        for(int j=0; j<numparam; j++)
        {
            hom_omflightsig[j]->Write();
        }
        for(int j=0; j<numparam; j++)
        {
            hom_distancesig[j]->Write();
        }
    }
    else
        //hom_NoCut->Write();
        hom_defaultcut->Write();
}

