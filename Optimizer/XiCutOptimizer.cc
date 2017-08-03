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
    bool Cut = true;
    const int CONTSIZE = 10000;
    TH1::SetDefaultSumw2();

    //Cut parameters to be varied. The commented elements are Hong's cuts to make a plot of them remember to change the numparam value as well
    //float xi_xi3dipsig[]     = {8.0  , 8.5  , 9.0  , 9.5  , 10.0 , 10.5 , 11.0 , 11.5, 12.0, 12.5, 13.0, 13.5};//, 2.5};
    int numparam             = 13;
    int Oparamindex          = 12; // For deciding which parameter that is not being varied to use
    float xi_xi3dipsig[]     = {4.0  , 4.3  , 4.6  , 4.9  , 5.2  , 5.5  , 5.8  , 6.1 , 6.4 , 6.7 , 7.0, 7.3, 7.6};//, 2.5};
    float xi_xipi3dipsig[]   = {5.0  , 4.9  , 4.8  , 4.7  , 4.6  , 4.5  , 4.4  , 4.3 , 4.2 , 4.1 , 4.0, 3.1, 3.3};//, 5.0};
    float xi_vtrkpi3dipsig[] = {4.5  , 4.6  , 4.7  , 4.8  , 4.9  , 5.0  , 5.1  , 5.2 , 5.3 , 5.4 , 5.5, 3.0, 3.3};//, 4.0};
    float xi_vtrkp3dipsig[]  = {2.5  , 2.6  , 2.7  , 2.8  , 2.9  , 3.0  , 3.1  , 3.2 , 3.3 , 3.4 , 3.5, 2.5, 2.5};//, 3.0};
    float xi_xiflightsig[]   = {2.5  , 2.6  , 2.7  , 2.8  , 2.9  , 3.0  , 3.1  , 3.2 , 3.3 , 3.4 , 3.5, 2.5, 2.5};//, 3.0};
    float xi_distancesig[]   = {12.0 , 11.8 , 11.6 , 11.4 , 11.0 , 10.5 , 10.0 , 9.5 , 9.0 , 8.5 , 8.0, 8.0, 8.5};//, 12.0};

    //Containers
    //Tree values
    float txi_xi3dipsig     [CONTSIZE];
    float txi_xipi3dipsig   [CONTSIZE];
    float txi_vtrkpi3dipsig [CONTSIZE];
    float txi_vtrkp3dipsig  [CONTSIZE];
    float txi_xiflightsig   [CONTSIZE];
    float txi_distancesig   [CONTSIZE];
    float txi_mass          [CONTSIZE];
    float txi_pt            [CONTSIZE];
    int   txi_n;

    //Hist Containers
    TH2D* hxi_xi3dipsig     [numparam];
    TH2D* hxi_xipi3dipsig   [numparam];
    TH2D* hxi_vtrkpi3dipsig [numparam];
    TH2D* hxi_vtrkp3dipsig  [numparam];
    TH2D* hxi_xiflightsig   [numparam];
    TH2D* hxi_distancesig   [numparam];

    //Tree setup
    TFile* f1=TFile::Open("/volumes/MacHD/Users/blt1/research/CascadeV2pPb/RootFiles/Flow/CasCutLoose/XiTree/XiOmTTree.root");
    TTree* XiTree = (TTree*)f1->Get("xiTree/XiTree");

    XiTree->SetBranchAddress("n",              &txi_n);
    XiTree->SetBranchAddress("xi3dipsig",      &txi_xi3dipsig);
    XiTree->SetBranchAddress("xipi3dipsig",    &txi_xipi3dipsig);
    XiTree->SetBranchAddress("vtrkpi3dipsig",  &txi_vtrkpi3dipsig);
    XiTree->SetBranchAddress("vtrkp3dipsig",   &txi_vtrkp3dipsig);
    XiTree->SetBranchAddress("xiflightsig",    &txi_xiflightsig);
    XiTree->SetBranchAddress("distancesig",    &txi_distancesig);
    XiTree->SetBranchAddress("mass",           &txi_mass);
    XiTree->SetBranchAddress("pt",             &txi_pt);

    //Intialize Histograms
    TH2D* hxi_NoCut = NULL;
    if(!Cut) hxi_NoCut = new TH2D("hxi_NoCut","NoCut",150,1.25,1.40,150,0,15);
    for(int j=0; j<numparam; j++)
    {
        hxi_xi3dipsig[j]     = new TH2D(Form("hxi_xi3dipsig_%.1f",xi_xi3dipsig[j]),Form("hxi_xi3dipsig_%.1f",xi_xi3dipsig[j]),150,1.25,1.40,150,0,15);
        hxi_xipi3dipsig[j]   = new TH2D(Form("hxi_xipi3dipsig_%.1f",xi_xipi3dipsig[j]),Form("hxi_xipi3dipsig_%.1f",xi_xipi3dipsig[j]),150,1.25,1.40,150,0,15);
        hxi_vtrkpi3dipsig[j] = new TH2D(Form("hxi_vtrkpi3dipsig_%.1f",xi_vtrkpi3dipsig[j]),Form("hxi_vtrkpi3dipsig_%.1f",xi_vtrkpi3dipsig[j]),150,1.25,1.40,150,0,15);
        hxi_vtrkp3dipsig[j]  = new TH2D(Form("hxi_vtrkp3dipsig_%.1f",xi_vtrkp3dipsig[j]),Form("hxi_vtrkp3dipsig_%.1f",xi_vtrkp3dipsig[j]),150,1.25,1.40,150,0,15);
        hxi_xiflightsig[j]   = new TH2D(Form("hxi_xiflightsig_%.1f",xi_xiflightsig[j]),Form("hxi_xiflightsig_%.1f",xi_xiflightsig[j]),150,1.25,1.40,150,0,15);
        hxi_distancesig[j]   = new TH2D(Form("hxi_distancesig_%.1f",xi_distancesig[j]),Form("hxi_distancesig_%.1f",xi_distancesig[j]),150,1.25,1.40,150,0,15);
    }

    //Output file creation
    std::string filetype = ".root";
    name += filetype;
    struct stat buffer;
    if(stat(name.c_str(), &buffer) == 0)
    {
        cout << "File with this name already exists, please select a different name" << endl;
        return;
    }
    TFile out(name.c_str(),"RECREATE");
    cout << name.c_str() << " created!" << endl;

    //Loop over Tree
    int xi_nEvent = XiTree->GetEntries();

    for(int i=0; i<xi_nEvent; i++)
    {
        if(Cut)
        {
            cout << "Event number: " << i << endl;
            XiTree->GetEntry(i);
            for( int j=0; j<numparam; j++)
            {
                //xi3dipsig
                for(int k=0; k<txi_n; k++)
                {
                    if(txi_xi3dipsig[k]     > xi_xi3dipsig[j])     continue;
                    if(txi_xipi3dipsig[k]   < xi_xipi3dipsig[Oparamindex])   continue;
                    if(txi_vtrkpi3dipsig[k] < xi_vtrkpi3dipsig[Oparamindex]) continue;
                    if(txi_vtrkp3dipsig[k]  < xi_vtrkp3dipsig[Oparamindex])  continue;
                    if(txi_xiflightsig[k]   < xi_xiflightsig[Oparamindex])   continue;
                    if(txi_distancesig[k]   < xi_distancesig[Oparamindex])   continue;

                    hxi_xi3dipsig[j]->Fill(txi_mass[k],txi_pt[k]);
                }

                //xipi3dipsig
                for(int k=0; k<txi_n; k++)
                {
                    if(txi_xi3dipsig[k] > xi_xi3dipsig[Oparamindex])         continue;
                    if(txi_xipi3dipsig[k] < xi_xipi3dipsig[j])     continue;
                    if(txi_vtrkpi3dipsig[k] < xi_vtrkpi3dipsig[Oparamindex]) continue;
                    if(txi_vtrkp3dipsig[k] < xi_vtrkp3dipsig[Oparamindex])   continue;
                    if(txi_xiflightsig[k] < xi_xiflightsig[Oparamindex])     continue;
                    if(txi_distancesig[k] < xi_distancesig[Oparamindex])     continue;

                    hxi_xipi3dipsig[j]->Fill(txi_mass[k],txi_pt[k]);
                }

                //vtrkpi3dipsig
                for(int k=0; k<txi_n; k++)
                {
                    if(txi_xi3dipsig[k] > xi_xi3dipsig[Oparamindex])         continue;
                    if(txi_xipi3dipsig[k] < xi_xipi3dipsig[Oparamindex])     continue;
                    if(txi_vtrkpi3dipsig[k] < xi_vtrkpi3dipsig[j]) continue;
                    if(txi_vtrkp3dipsig[k] < xi_vtrkp3dipsig[Oparamindex])   continue;
                    if(txi_xiflightsig[k] < xi_xiflightsig[Oparamindex])     continue;
                    if(txi_distancesig[k] < xi_distancesig[Oparamindex])     continue;

                    hxi_vtrkpi3dipsig[j]->Fill(txi_mass[k],txi_pt[k]);
                }

                //vtrkp3dipsig
                for(int k=0; k<txi_n; k++)
                {
                    if(txi_xi3dipsig[k] > xi_xi3dipsig[Oparamindex])         continue;
                    if(txi_xipi3dipsig[k] < xi_xipi3dipsig[Oparamindex])     continue;
                    if(txi_vtrkpi3dipsig[k] < xi_vtrkpi3dipsig[Oparamindex]) continue;
                    if(txi_vtrkp3dipsig[k] < xi_vtrkp3dipsig[j])   continue;
                    if(txi_xiflightsig[k] < xi_xiflightsig[Oparamindex])     continue;
                    if(txi_distancesig[k] < xi_distancesig[Oparamindex])     continue;

                    hxi_vtrkp3dipsig[j]->Fill(txi_mass[k],txi_pt[k]);
                }

                //xiflightsig
                for(int k=0; k<txi_n; k++)
                {
                    if(txi_xi3dipsig[k] > xi_xi3dipsig[Oparamindex])         continue;
                    if(txi_xipi3dipsig[k] < xi_xipi3dipsig[Oparamindex])     continue;
                    if(txi_vtrkpi3dipsig[k] < xi_vtrkpi3dipsig[Oparamindex]) continue;
                    if(txi_vtrkp3dipsig[k] < xi_vtrkp3dipsig[Oparamindex])   continue;
                    if(txi_xiflightsig[k] < xi_xiflightsig[j])     continue;
                    if(txi_distancesig[k] < xi_distancesig[Oparamindex])     continue;

                    hxi_xiflightsig[j]->Fill(txi_mass[k],txi_pt[k]);
                }

                //distancesig
                for(int k=0; k<txi_n; k++)
                {
                    if(txi_xi3dipsig[k] > xi_xi3dipsig[Oparamindex])         continue;
                    if(txi_xipi3dipsig[k] < xi_xipi3dipsig[Oparamindex])     continue;
                    if(txi_vtrkpi3dipsig[k] < xi_vtrkpi3dipsig[Oparamindex]) continue;
                    if(txi_vtrkp3dipsig[k] < xi_vtrkp3dipsig[Oparamindex])   continue;
                    if(txi_xiflightsig[k] < xi_xiflightsig[Oparamindex])     continue;
                    if(txi_distancesig[k] < xi_distancesig[j])     continue;

                    hxi_distancesig[j]->Fill(txi_mass[k],txi_pt[k]);
                }
            }
        }
        else
        {
            cout << "Event number: " << i << endl;
            XiTree->GetEntry(i);
            for(int k=0; k<txi_n; k++)
            {
                hxi_NoCut->Fill(txi_mass[k],txi_pt[k]);
            }
        }
    }

    //Write histograms to root file
    if(Cut)
    {
        for(int j=0; j<numparam; j++)
        {
            hxi_xi3dipsig[j]->Write();
        }
        for(int j=0; j<numparam; j++)
        {
            hxi_xipi3dipsig[j]->Write();
        }
        for(int j=0; j<numparam; j++)
        {
            hxi_vtrkpi3dipsig[j]->Write();
        }
        for(int j=0; j<numparam; j++)
        {
            hxi_vtrkp3dipsig[j]->Write();
        }
        for(int j=0; j<numparam; j++)
        {
            hxi_xiflightsig[j]->Write();
        }
        for(int j=0; j<numparam; j++)
        {
            hxi_distancesig[j]->Write();
        }
    }
    else
        hxi_NoCut->Write();
}

