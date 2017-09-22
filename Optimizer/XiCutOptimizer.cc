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
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TROOT.h"

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

    //Hist Containers
    TH2D* hom_om3dipsig     [numparam];
    TH2D* hom_omKaon3dipsig   [numparam];
    TH2D* hom_vtrkpi3dipsig [numparam];
    TH2D* hom_vtrkp3dipsig  [numparam];
    TH2D* hom_omflightsig   [numparam];
    TH2D* hom_distancesig   [numparam];

    //Tree setup
    TFile* f1=TFile::Open("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/All/V0CasTree_09_14_17.root");

    TTreeReader reader("OmTreeProducer/OmTree",f1);

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
    TTreeReaderArray<int> nTrkAcc(reader      ,"nTrkAcc.nTrkAcc");
    TTreeReaderArray<float> misIDMassLapi(reader,"misIDMassLapi.misIDMassLapi");
    TTreeReaderArray<float> misIDMasspiLa(reader,"misIDMasspiLa.misIDMasspiLa");

    //Intialize Histograms
    TH2D* hom_NoCut = NULL;
    if(!Cut) hom_NoCut = new TH2D("hom_NoCut"    ,"NoCut"  ,150,1.25,1.40,150,0,15);
    TH2D* hom_defaultcut = new TH2D("hom_Default","Default",150,1.60,1.75,400,0,40);
    for(int j=0; j<numparam; j++)
    {
        hom_om3dipsig[j]     = new TH2D(Form("hom_om3dipsig_%.1f"    ,om_om3dipsig[j])    ,Form("hom_om3dipsig_%.1f"    ,om_om3dipsig[j])    ,150,1.25,1.40,150,0,15);
        hom_omKaon3dipsig[j] = new TH2D(Form("hom_omKaon3dipsig_%.1f",om_omKaon3dipsig[j]),Form("hom_omKaon3dipsig_%.1f",om_omKaon3dipsig[j]),150,1.25,1.40,150,0,15);
        hom_vtrkpi3dipsig[j] = new TH2D(Form("hom_vtrkpi3dipsig_%.1f",om_vtrkpi3dipsig[j]),Form("hom_vtrkpi3dipsig_%.1f",om_vtrkpi3dipsig[j]),150,1.25,1.40,150,0,15);
        hom_vtrkp3dipsig[j]  = new TH2D(Form("hom_vtrkp3dipsig_%.1f" ,om_vtrkp3dipsig[j]) ,Form("hom_vtrkp3dipsig_%.1f" ,om_vtrkp3dipsig[j]) ,150,1.25,1.40,150,0,15);
        hom_omflightsig[j]   = new TH2D(Form("hom_omflightsig_%.1f"  ,om_omflightsig[j])  ,Form("hom_omflightsig_%.1f"  ,om_omflightsig[j])  ,150,1.25,1.40,150,0,15);
        hom_distancesig[j]   = new TH2D(Form("hom_distancesig_%.1f"  ,om_distancesig[j])  ,Form("hom_distancesig_%.1f"  ,om_distancesig[j])  ,150,1.25,1.40,150,0,15);
    }

    //Output file creation
    std::string filetype = ".root";
    name += filetype;
    struct stat buffer;
    if(stat(name.c_str(), &buffer) == 0)
    {
        cout << "File with this name already eomsts, please select a different name" << endl;
        return false;
    }
    TFile out(name.c_str(),"RECREATE");
    cout << name.c_str() << " created!" << endl;

    while(reader.Next())
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
        if(!CheckValue(rap))           return false;
        if(!CheckValue(nTrkAcc))       return false;
        if(!CheckValue(misIDMassLapi)) return false;
        if(!CheckValue(misIDMasspiLa)) return false;
        /*
        if(Cut)
        {
            OmTree->GetEntry(i);
            for( int j=0; j<numparam; j++)
            {
                //om3dipsig
                for(int k=0; k<tom_n; k++)
                {
                    if(om3dipsig[k]     > om_om3dipsig[j])     continue;
                    if(omKaon3dipsig[k]   < om_omKaon3dipsig[Oparamindex])   continue;
                    if(vtrkpi3dipsig[k] < om_vtrkpi3dipsig[Oparamindex]) continue;
                    if(vtrkp3dipsig[k]  < om_vtrkp3dipsig[Oparamindex])  continue;
                    if(omflightsig[k]   < om_omflightsig[Oparamindex])   continue;
                    if(distancesig[k]   < om_distancesig[Oparamindex])   continue;

                    hom_om3dipsig[j]->Fill(mass[k],pt[k]);
                }

                //ompi3dipsig
                for(int k=0; k<tom_n; k++)
                {
                    if(om3dipsig[k] > om_om3dipsig[Oparamindex])         continue;
                    if(omKaon3dipsig[k] < om_omKaon3dipsig[j])     continue;
                    if(vtrkpi3dipsig[k] < om_vtrkpi3dipsig[Oparamindex]) continue;
                    if(vtrkp3dipsig[k] < om_vtrkp3dipsig[Oparamindex])   continue;
                    if(omflightsig[k] < om_omflightsig[Oparamindex])     continue;
                    if(distancesig[k] < om_distancesig[Oparamindex])     continue;

                    hom_omKaon3dipsig[j]->Fill(mass[k],pt[k]);
                }

                //vtrkpi3dipsig
                for(int k=0; k<tom_n; k++)
                {
                    if(om3dipsig[k] > om_om3dipsig[Oparamindex])         continue;
                    if(omKaon3dipsig[k] < om_omKaon3dipsig[Oparamindex])     continue;
                    if(vtrkpi3dipsig[k] < om_vtrkpi3dipsig[j]) continue;
                    if(vtrkp3dipsig[k] < om_vtrkp3dipsig[Oparamindex])   continue;
                    if(omflightsig[k] < om_omflightsig[Oparamindex])     continue;
                    if(distancesig[k] < om_distancesig[Oparamindex])     continue;

                    hom_vtrkpi3dipsig[j]->Fill(mass[k],pt[k]);
                }

                //vtrkp3dipsig
                for(int k=0; k<tom_n; k++)
                {
                    if(om3dipsig[k] > om_om3dipsig[Oparamindex])         continue;
                    if(omKaon3dipsig[k] < om_omKaon3dipsig[Oparamindex])     continue;
                    if(vtrkpi3dipsig[k] < om_vtrkpi3dipsig[Oparamindex]) continue;
                    if(vtrkp3dipsig[k] < om_vtrkp3dipsig[j])   continue;
                    if(omflightsig[k] < om_omflightsig[Oparamindex])     continue;
                    if(distancesig[k] < om_distancesig[Oparamindex])     continue;

                    hom_vtrkp3dipsig[j]->Fill(mass[k],pt[k]);
                }

                //omflightsig
                for(int k=0; k<tom_n; k++)
                {
                    if(om3dipsig[k] > om_om3dipsig[Oparamindex])         continue;
                    if(omKaon3dipsig[k] < om_omKaon3dipsig[Oparamindex])     continue;
                    if(vtrkpi3dipsig[k] < om_vtrkpi3dipsig[Oparamindex]) continue;
                    if(vtrkp3dipsig[k] < om_vtrkp3dipsig[Oparamindex])   continue;
                    if(omflightsig[k] < om_omflightsig[j])     continue;
                    if(distancesig[k] < om_distancesig[Oparamindex])     continue;

                    hom_omflightsig[j]->Fill(mass[k],pt[k]);
                }

                //distancesig
                for(int k=0; k<tom_n; k++)
                {
                    if(om3dipsig[k] > om_om3dipsig[Oparamindex])         continue;
                    if(omKaon3dipsig[k] < om_omKaon3dipsig[Oparamindex])     continue;
                    if(vtrkpi3dipsig[k] < om_vtrkpi3dipsig[Oparamindex]) continue;
                    if(vtrkp3dipsig[k] < om_vtrkp3dipsig[Oparamindex])   continue;
                    if(omflightsig[k] < om_omflightsig[Oparamindex])     continue;
                    if(distancesig[k] < om_distancesig[j])     continue;

                    hom_distancesig[j]->Fill(mass[k],pt[k]);
                }
            }
        }
        else
        {
            */
                //hom_NoCut->Fill(mass[k],pt[k]);
                for(int i=0;i<mass.GetSize();i++)
                {
                    if(nTrkAcc[i]                 > multHigh_)           continue;
                    if(om3dipsig[i]               > om_om3dipsig[0])     continue;
                    if(std::fabs(rap[i])          > rapidity)            continue;
                    if(omKaon3dipsig[i]           < om_omKaon3dipsig[0]) continue;
                    if(vtrkpi3dipsig[i]           < om_vtrkpi3dipsig[0]) continue;
                    if(vtrkp3dipsig[i]            < om_vtrkp3dipsig[0])  continue;
                    if(omflightsig[i]             < om_omflightsig[0])   continue;
                    if(distancesig[i]             < om_distancesig[0])   continue;
                    if(std::abs(misIDMasspiLa[i]) < misIDMass)           continue;
                    if(std::abs(misIDMassLapi[i]) < misIDMass)           continue;

                    hom_defaultcut->Fill(mass[i],pt[i]);
                }
        //}
    }

    //Write histograms to root file
    if(Cut)
    {
        for(int j=0; j<numparam; j++)
            hom_om3dipsig[j]->Write();
        for(int j=0; j<numparam; j++)
            hom_omKaon3dipsig[j]->Write();
        for(int j=0; j<numparam; j++)
            hom_vtrkpi3dipsig[j]->Write();
        for(int j=0; j<numparam; j++)
            hom_vtrkp3dipsig[j]->Write();
        for(int j=0; j<numparam; j++)
            hom_omflightsig[j]->Write();
        for(int j=0; j<numparam; j++)
            hom_distancesig[j]->Write();
    }
    else
    {
        //hom_NoCut->Write();
        hom_defaultcut->Write();
    }
    return true;
}

