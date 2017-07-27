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

void XiCutOptimizer(std::string name)
{
    /* Output to file example
    ofstream myfile;
    myfile.open ("Log/Test.txt", std::ios_base::app);
    for(int i = 0; i < 10; i++)
    {
        TString s = Form("Test_%d\n",i);
        myfile << s;
    }
    */

    //Initializer constants
    const int LOOPSIZE = 3;
    const int CONTSIZE = 6;

    //Containers

    //Cut parameters to be varied. First element is Hong's cut value
    int numparam             = 1;
    float xi_xi3dipsig[]     = {2.5};
    float xi_xipi3dipsig[]   = {5};
    float xi_vtrkpi3dipsig[] = {4};
    float xi_vtrkp3dipsig[]  = {3};
    float xi_xiflightsig[]   = {3};
    float xi_distancesig[]   = {12};

    //Tree values
    float txi_xi3dipsig     [CONTSIZE];
    float txi_xipi3dipsig   [CONTSIZE];
    float txi_vtrkpi3dipsig [CONTSIZE];
    float txi_vtrkp3dipsig  [CONTSIZE];
    float txi_xiflightsig   [CONTSIZE];
    float txi_distancesig   [CONTSIZE];
    float txi_mass          [CONTSIZE];
    float txi_pt            [CONTSIZE];
    //int   txi_n;

    //Hist Containers
    TH1D hxi_xi3dipsig     [numparam];
    TH1D hxi_xipi3dipsig   [numparam];
    TH1D hxi_vtrkpi3dipsig [numparam];
    TH1D hxi_vtrkp3dipsig  [numparam];
    TH1D hxi_xiflightsig   [numparam];
    TH1D hxi_distancesig   [numparam];

    //Output TFile setup
    //std::string name = "Test.root";
    TFile Out(name.c_str(),"RECREATE");

    //Tree setup
    TFile* f1=TFile::Open("/volumes/MacHD/Users/blt1/research/CascadeV2pPb/RootFiles/Flow/CasCutLoose/XiTree/XiOmTTree.root");

    TTree* XiTree = (TTree*)f1->Get("xiTree/XiTree");
    //
    //XiTree->SetBranchAddress("n"             , &txi_n);
    XiTree->SetBranchAddress("xi3dipsig",      &txi_xi3dipsig);
    XiTree->SetBranchAddress("xipi3dipsig",    &txi_xipi3dipsig);
    XiTree->SetBranchAddress("vtrkpi3dipsig",  &txi_vtrkpi3dipsig);
    XiTree->SetBranchAddress("vtrkp3dipsig",   &txi_vtrkp3dipsig);
    XiTree->SetBranchAddress("xiflightsig",    &txi_xiflightsig);
    XiTree->SetBranchAddress("distancesig",    &txi_distancesig);
    XiTree->SetBranchAddress("mass",           &txi_mass);
    XiTree->SetBranchAddress("pt",             &txi_pt);

    //Intialize Histograms
    for(int j=0; j<numparam; j++)
    {
        hxi_xi3dipsig[j]     = new TH1D(Form("hxi_xi3dipsig_%f",xi_xi3dipsig[j]),Form("hxi_xi3dipsig_%f",xi_xi3dipsig[j]),150,1.25,1.40);
        hxi_xipi3dipsig[j]   = new TH1D(Form("hxi_xipi3dipsig_%f",xi_xipi3dipsig[j]),Form("hxi_xipi3dipsig_%f",xi_xipi3dipsig[j]),150,1.25,1.40);
        hxi_vtrkpi3dipsig[j] = new TH1D(Form("hxi_vtrkpi3dipsig_%f",xi_vtrkpi3dipsig[j]),Form("hxi_vtrkpi3dipsig_%f",xi_vtrkpi3dipsig[j]),150,1.25,1.40);
        hxi_vtrkp3dipsig[j]  = new TH1D(Form("hxi_vtrkp3dipsig_%f",xi_vtrkp3dipsig[j]),Form("hxi_vtrkp3dipsig_%f",xi_vtrkp3dipsig[j]),150,1.25,1.40);
        hxi_xiflightsig[j]   = new TH1D(Form("hxi_xiflightsig_%f",xi_xiflightsig[j]),Form("hxi_xiflightsig_%f",xi_xiflightsig[j]),150,1.25,1.40);
        hxi_distancesig[j]   = new TH1D(Form("hxi_distancesig_%f",xi_distancesig[j]),Form("hxi_distancesig_%f",xi_distancesig[j]),150,1.25,1.40);
    }

    //Loop over Tree
    int xi_nEvent = XiTree->GetEntries();

    //for(int i=0; i<xi_nEvent; i++)
    for(int i=0; i<100; i++)
    {
        cout << i << endl;
        XiTree->GetEntry(i);
        for( int j=0; j<numparam; j++)
        {
            //xi3dipsig
            for(int k=0; k<LOOPSIZE; k++)
            {
                if(txi_xi3dipsig[k] < xi_xi3dipsig[j])
                {
                    hxi_xi3dipsig[j]->Fill(txi_mass[k]);
                }
            }
            Out.Write();
            //if(txi_xi3dipsig[LOOPSIZE+1] != 0 || txi_xi3dipsig[LOOPSIZE+2] != 0 || txi_xi3dipsig[LOOPSIZE+1] != 0) cout << "Loop size too small for xi3dipsig for j = " << j <<  endl;

            //xipi3dipsig
            for(int k=0; k<LOOPSIZE; k++)
            {
                if(txi_xipi3dipsig[k] < xi_xipi3dipsig[j])
                {
                    hxi_xipi3dipsig[j]->Fill(txi_mass[k]);
                }
            }
            Out.Write();
            //if(txi_xipi3dipsig[LOOPSIZE+1] != 0 || txi_xipi3dipsig[LOOPSIZE+2] != 0 || txi_xipi3dipsig[LOOPSIZE+1] != 0) cout << "Loop size too small for xipi3dipsig for j = " << j <<  endl;

            //vtrkpi3dipsig
            for(int k=0; k<LOOPSIZE; k++)
            {
                if(txi_vtrkpi3dipsig[k] < xi_vtrkpi3dipsig[j])
                {
                    hxi_vtrkpi3dipsig[j]->Fill(txi_mass[k]);
                }
            }
            Out.Write();
            //if(txi_vtrkpi3dipsig[LOOPSIZE+1] != 0 || txi_vtrkpi3dipsig[LOOPSIZE+2] != 0 || txi_vtrkpi3dipsig[LOOPSIZE+1] != 0) cout << "Loop size too small for vtrkpi3dipsig for j = " << j <<  endl;

            //vtrkp3dipsig
            for(int k=0; k<LOOPSIZE; k++)
            {
                if(txi_vtrkp3dipsig[k] < xi_vtrkp3dipsig[j])
                {
                    hxi_vtrkp3dipsig[j]->Fill(txi_mass[k]);
                }
            }
            Out.Write();
            //if(txi_vtrkp3dipsig[LOOPSIZE+1] != 0 || txi_vtrkp3dipsig[LOOPSIZE+2] != 0 || txi_vtrkp3dipsig[LOOPSIZE+1] != 0) cout << "Loop size too small for vtrkp3dipsig for j = " << j <<  endl;

            //xiflightsig
            for(int k=0; k<LOOPSIZE; k++)
            {
                if(txi_xiflightsig[k] < xi_xiflightsig[j])
                {
                    hxi_xiflightsig[j]->Fill(txi_mass[k]);
                }
            }
            Out.Write();
            //if(txi_xiflightsig[LOOPSIZE+1] != 0 || txi_xiflightsig[LOOPSIZE+2] != 0 || txi_xiflightsig[LOOPSIZE+1] != 0) cout << "Loop size too small for xiflightsig for j = " << j <<  endl;

            //distancesig
            for(int k=0; k<LOOPSIZE; k++)
            {
                if(txi_distancesig[k] < xi_distancesig[j])
                {
                    hxi_distancesig[j]->Fill(txi_mass[k]);
                }
            }
            Out.Write();
            //if(txi_distancesig[LOOPSIZE+1] != 0 || txi_distancesig[LOOPSIZE+2] != 0 || txi_distancesig[LOOPSIZE+1] != 0) cout << "Loop size too small for distancesig for j = " << j <<  endl;
        }
    }

    //Out.Write();
}
