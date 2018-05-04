#include "interface/PhiTreeReader.h"
#include "TString.h"
#include "TKey.h"
#include "TRegexp.h"
#include <sys/stat.h>
#include <stdlib.h>
#include <stdio.h>
//#include <boost/filesystem.hpp>

std::vector<TFile*> InitializeFiles(std::vector<std::string> files_)
{
    std::vector<TFile*> files_container;
    for(unsigned i=0; i<files_.size(); ++i)
    {
        TFile* tmp = new TFile(files_[i].c_str());
        files_container.push_back(tmp);
    }
    return files_container;
}

void PhiMassPtReader()
{
    std::map<int,TH2D*> h_Mass_Pt;

    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.05);

    std::vector<std::string> file_names = {
        //"PhiTreeProof/outputFiles/Mass_pt_loose_10_10_24_10_5.root",
        //"PhiTreeProof/outputFiles/Mass_pt_tight_10_10_24_10_5.root"
        //"PhiTreeProof/outputFiles/Phi_MassPt_v1.root"
        //"PhiTreeProof/outputFiles/Phi_MassPt_v2.root"
        "PhiTreeProof/outputFiles/Phi_MassPt_v3.root"
    };

    std::vector<std::string> output_names = {
        "Phi_MassPt_noDedx_PhiDauRapCut",
        "Phi_MassPt_wDedx_PhiDauRapCut",
        "Phi_MassPt_DCA3p0_PhiDauRapCut"
    };

    std::vector<TString> reg_list = {
        "mass"
    };

    std::vector<double> pt_bins = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0};
    //std::vector<double> pt_bins = {1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.2,10.0};

    std::vector<TFile*> file_list = InitializeFiles(file_names);

    std::vector<TH1D*> hist_container;
    std::vector<TCanvas*> canvas_container;

    std::vector<TRegexp> re;

    for(unsigned j=0; j<file_list.size(); ++j)
    {
        h_Mass_Pt.clear();
        for(unsigned i=0; i<reg_list.size(); ++i)
        {
            TRegexp tmp_re(reg_list[i],kTRUE);
            re.push_back(tmp_re);
        }

        /*
         * Get all of the histograms with the names given in reg_list and put them into the corresponding histogram vector
         */
        for(const auto&& object : *file_list[j]->GetListOfKeys())
        {
            TKey* key = (TKey*)object;
            TString st = key->GetName();
            for(TRegexp tmp_Regexp : re)
            {
                if(st.Index(tmp_Regexp) == kNPOS) continue;
                std::string obj_name = key->GetName();
                short cycle = key->GetCycle();
                ostringstream name;
                name << obj_name + ";" << cycle;
                if(obj_name == "mass")
                    h_Mass_Pt[cycle] = (TH2D*)file_list[j]->Get(name.str().c_str());
            }
        }

        /*
         * Loop through the histograms in the file
         */
        for(auto const& mass_pt : h_Mass_Pt)
        {
            for(unsigned i=0; i<pt_bins.size()-1; ++i)
            {
                int pt_bin_min = -10;
                int pt_bin_max = -10;
                pt_bin_min = mass_pt.second->GetYaxis()->FindBin(pt_bins[i]);
                pt_bin_max = mass_pt.second->GetYaxis()->FindBin(pt_bins[i+1]) - 1;
                ostringstream hist_name;
                hist_name << "Mass_" << mass_pt.first << "_%d";
                TH1D* h_mass = mass_pt.second->ProjectionX(Form(hist_name.str().c_str(),i),pt_bin_min,pt_bin_max,"eo");

                hist_container.emplace_back(h_mass);
            }

            for(unsigned i=0; i<hist_container.size(); ++i)
            {
                ostringstream canvas_name;
                canvas_name << "Pt_bin_" << mass_pt.first << "_%d";
                TCanvas *c = new TCanvas(Form(canvas_name.str().c_str(),i),Form("Pt_bin_%d",i),600,600);
                hist_container[i]->Draw();
                std::ostringstream os;
                os << pt_bins[i] << " - " << pt_bins[i+1];
                tex->DrawLatex(0.2,0.8,os.str().c_str());
                canvas_container.emplace_back(c);
            }


            /*
             * Output the file containing the mass pt binned histograms
             */
            //if(!boost::filesystem::is_directory("massPeaks/" + output_names[mass_pt.first - 1]))
            //{
                //if(boost::filesystem::create_directory("massPeaks/" + output_names[mass_pt.first - 1]))
                    //cout << "massPeaks/" + output_names[mass_pt.first - 1] + " created!" << endl;
            //}

            mkdir(("massPeaks/" + output_names[mass_pt.first - 1]).c_str(),0775);
            TFile* out = new TFile(("massPeaks/" + output_names[mass_pt.first - 1] + "/" + output_names[mass_pt.first - 1] + ".root").c_str(),"RECREATE");
            for(unsigned i=0; i<hist_container.size(); ++i)
            {
                canvas_container[i]->Write();
                hist_container[i]->Write();

                std::string name_pdf = "massPeaks/" + output_names[mass_pt.first - 1] + "/" + (std::string)canvas_container[i]->GetName() + ".pdf";
                //std::string name_png = "massPeaks/Images/" + (std::string)canvas_container[i]->GetName() + ".png";
                canvas_container[i]->Print(name_pdf.c_str());
                //canvas_container[i]->Print(name_png.c_str());
            }
            hist_container.clear();
            canvas_container.clear();
            out->Close();
            delete out;
        }
    }

}
