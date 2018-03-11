#include "interface/PhiTreeReader.h"

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
    TH2D* h_Mass_Pt;

    TLatex* tex = new TLatex();
    tex->SetNDC();
    tex->SetTextFont(42);
    tex->SetTextSize(0.05);

    std::vector<std::string> file_names = {
        "PhiTreeProof/outputFiles/Mass_pt_loose_10_10_24_10_5.root",
        "PhiTreeProof/outputFiles/Mass_pt_tight_10_10_24_10_5.root"
    };

    std::vector<double> pt_bins = {0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0};
    //std::vector<double> pt_bins = {1.0,1.4,1.8,2.2,2.8,3.6,4.6,6.0,7.2,10.0};

    std::vector<TFile*> file_list = InitializeFiles(file_names);

    std::vector<TH1D*> hist_container;
    std::vector<TCanvas*> canvas_container;

    file_list[0]->GetObject("mass",h_Mass_Pt);

    for(unsigned i=0; i<pt_bins.size()-1; ++i)
    {
        int pt_bin_min = -10;
        int pt_bin_max = -10;
        pt_bin_min = h_Mass_Pt->GetYaxis()->FindBin(pt_bins[i]);
        pt_bin_max = h_Mass_Pt->GetYaxis()->FindBin(pt_bins[i+1]) - 1;
        TH1D* h_mass = h_Mass_Pt->ProjectionX(Form("Mass_%d",i),pt_bin_min,pt_bin_max,"eo");

        hist_container.push_back(h_mass);
    }

    for(unsigned i=0; i<hist_container.size(); ++i)
    {
        TCanvas *c = new TCanvas(Form("Pt_bin_%d",i),Form("Pt_bin_%d",i),600,600);
        hist_container[i]->Draw();
        std::ostringstream os;
        os << pt_bins[i] << " - " << pt_bins[i+1];
        tex->DrawLatex(0.2,0.8,os.str().c_str());
        canvas_container.push_back(c);
        //delete c;
    }

    TFile* out = new TFile("Mass_pt_loose.root","RECREATE");
    for(unsigned i=0; i<hist_container.size(); ++i)
    {
        canvas_container[i]->Write();
        hist_container[i]->Write();
        
        std::string name_pdf = "massPeaks/" + (std::string)canvas_container[i]->GetName() + ".pdf";
        //std::string name_png = "massPeaks/Images/" + (std::string)canvas_container[i]->GetName() + ".png";
        canvas_container[i]->Print(name_pdf.c_str());
        //canvas_container[i]->Print(name_png.c_str());
    }

    out->Close();
}
