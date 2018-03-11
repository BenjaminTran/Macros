#include "../interface/PhiTreeReader.h"

#include "TH3D.h"
#include "TSelector.h"

class PhiTreeProof : public TSelector {
    public:
        struct kaon_proof{
            double p;
            double pt;
            double px;
            double py;
            double pz;
            double dedx;
            double energy;
            int charge;
            int nhits;

            kaon_proof(double p_, double pt_, double px_, double py_, double pz_, double dedx_, double energy_, int charge_, int nhits_) :
                p(p_), pt(pt_), px(px_), py(py_), pz(pz_), dedx(dedx_), energy(energy_), charge(charge_), nhits(nhits_)  {}
        };

        TTreeReader reader;
        TTreeReaderArray<float> momentum;
        TTreeReaderArray<float> px;
        TTreeReaderArray<float> py;
        TTreeReaderArray<float> pz;
        TTreeReaderArray<float> pt;
        TTreeReaderArray<float> ptError;
        TTreeReaderArray<float> energy;
        TTreeReaderArray<float> dedx;
        TTreeReaderArray<float> dz;
        TTreeReaderArray<float> dzError;
        TTreeReaderArray<float> dxy;
        TTreeReaderArray<float> dxyError;
        TTreeReaderArray<float> eta;
        TTreeReaderArray<float> rapidity;
        TTreeReaderArray<float> phi;
        TTreeReaderArray<float> vx;
        TTreeReaderArray<float> vy;
        TTreeReaderArray<float> vz;
        TTreeReaderArray<float> chi2norm;
        TTreeReaderArray<float> nhits;
        TTreeReaderArray<float> charge;

        TH1D* h_nEvt;
        TH1D* h_mass;
        TH2D* h_dedx;
        TH2D* h_mass_pt;

        TF1 *f_Dedx_bot;
        TF1 *f_Dedx_top;

        TFile* out;

        std::vector<kaon_proof> Pkp;
        std::vector<kaon_proof> Pkm;

        std::vector<TH1D*> hist_containers;

        double mass_low;
        double mass_high;
        double mass_numBins;

        double cut_dcaz;
        double cut_dcaxy;
        double cut_eta;
        double cut_rapidity;
        double cut_nhits;
        std::string cut_dedx;

        PhiTreeProof(TTree * = 0) :
            h_nEvt(0)
            ,h_mass(0)
            ,h_dedx(0)
            ,out(0)
            ,momentum(reader,"momentum")
            ,px(reader,"px")
            ,py(reader,"py")
            ,pz(reader,"pz")
            ,pt(reader,"pt")
            ,ptError(reader,"ptError")
            ,energy(reader,"energy")
            ,dedx(reader,"dedx")
            ,dz(reader,"dz")
            ,dzError(reader,"dzError")
            ,dxy(reader,"dxy")
            ,dxyError(reader,"dxyError")
            ,eta(reader,"eta")
            ,rapidity(reader,"rapidity")
            ,phi(reader,"phi")
            ,vx(reader,"vx")
            ,vy(reader,"vy")
            ,vz(reader,"vz")
            ,chi2norm(reader,"chi2norm")
            ,nhits(reader,"nhits")
            ,charge(reader,"charge")
            ,f_Dedx_bot(0)
            ,f_Dedx_top(0)
            {
                mass_low = 1.0;
                mass_high = 1.05;
                mass_numBins = 50;
                cut_dcaz = 1.0;
                cut_dcaxy = 1.0;
                cut_eta = 2.4;
                cut_rapidity = 1.0;
                cut_nhits = 5;
                cut_dedx = "tight";
            }

        virtual ~PhiTreeProof() {}
        virtual void Init(TTree *tree);
        virtual void SlaveBegin(TTree *tree);
        virtual bool Process(Long64_t entry);
        virtual void Terminate();
        virtual int Version() const {return 2;}

        static void CombinatorialMass(std::vector<kaon_proof> PKp, std::vector<kaon_proof> PKm, TH1D* h_mass_);
        static void CombinatorialMass2D(std::vector<kaon_proof> PKp, std::vector<kaon_proof> PKm, TH2D* h_mass_);

        ClassDef(PhiTreeProof,0);

};
