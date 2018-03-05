#include "../interface/PhiTreeReader.h"

#include "TSelector.h"

class PhiTreeProof11 : public TSelector {
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

        TF1 f_Dedx_bot;
        TF1 f_Dedx_top;

        TFile* out;

        std::vector<kaon_proof> Pkp;
        std::vector<kaon_proof> Pkm;

        PhiTreeProof11(TTree * = 0) :
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
            //tight
            ,f_Dedx_bot("f_Dedx_bot","0.55*(TMath::Power(1.15/x,2) - 1.7*TMath::Power(0.6/x,1)) + 2.7",0.01,100)
            ,f_Dedx_top("f_Dedx_top","0.55*(TMath::Power(1.62/x,2) - 2.95*TMath::Power(0.6/x,1)) + 3.3",0.01,100) {}
            //loose
            //,f_Dedx_bot("f_Dedx_bot","0.55*(TMath::Power(1.15/x,2) - 1.7*TMath::Power(0.6/x,1)) + 2.5",0.01,5)
            //,f_Dedx_top("f_Dedx_top","0.55*(TMath::Power(1.62/x,2) - 2*TMath::Power(0.6/x,1)) + 3.6",0.01,5) {}

        virtual ~PhiTreeProof11() {}
        virtual void Init(TTree *tree);
        virtual void SlaveBegin(TTree *tree);
        virtual bool Process(Long64_t entry);
        virtual void Terminate();
        virtual int Version() const {return 2;}

        static void CombinatorialMass(std::vector<kaon_proof> PKp, std::vector<kaon_proof> PKm, TH1D* h_mass_);

        ClassDef(PhiTreeProof11,0);

};
