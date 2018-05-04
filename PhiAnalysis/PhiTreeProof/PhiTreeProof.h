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
        //TH1D* h_mass;
        //TH2D* h_dedx_0;
        //TH2D* h_dedx_1;
        //TH2D* h_dedx_2;
        //TH2D* h_mass_pt_0;
        //TH2D* h_mass_pt_1;
        //TH2D* h_mass_pt_2;
        //
        TH2D **h_dedx;
        TH2D **h_mass_pt;

        int num_analyze = 3;

        std::vector<TH2D*> v_hdedx;
        std::vector<TH2D*> v_hmasspt;

        //TF1 *f_Dedx_bot;
        //TF1 *f_Dedx_top;

        TFile* out;

        std::vector<kaon_proof> Pkp;
        std::vector<kaon_proof> Pkm;

        std::vector<TH1D*> hist_containers;

        double mass_low;
        double mass_high;
        double mass_numBins;

        bool doDedx;

        std::vector<double> cut_dcaz;
        std::vector<double> cut_dcaxy;
        std::vector<double> cut_eta;
        std::vector<double> cut_rapidity;
        std::vector<int> cut_nhits;
        std::vector<std::string> cut_dedx;
        std::map<std::string,TF1*> f_Dedx_bot;
        std::map<std::string,TF1*> f_Dedx_top;

        PhiTreeProof(TTree * = 0) :
            //momentum(reader,"momentum")
            //,px(reader,"px")
            //,py(reader,"py")
            //,pz(reader,"pz")
            //,pt(reader,"pt")
            //,ptError(reader,"ptError")
            //,energy(reader,"energy")
            //,dedx(reader,"dedx")
            //,dz(reader,"dz")
            //,dzError(reader,"dzError")
            //,dxy(reader,"dxy")
            //,dxyError(reader,"dxyError")
            //,eta(reader,"eta")
            //,rapidity(reader,"rapidity")
            //,phi(reader,"phi")
            //,vx(reader,"vx")
            //,vy(reader,"vy")
            //,vz(reader,"vz")
            //,chi2norm(reader,"chi2norm")
            //,nhits(reader,"nhits")
            //,charge(reader,"charge")

            /*
             * For TMVA trees
             */
            momentum(reader,"momentum_1")
            ,px(reader,"px_1")
            ,py(reader,"py_1")
            ,pz(reader,"pz_1")
            ,pt(reader,"pt_1")
            ,ptError(reader,"ptError_1")
            ,energy(reader,"energy_1")
            ,dedx(reader,"dedx_1")
            ,dz(reader,"dz_1")
            ,dzError(reader,"dzError_1")
            ,dxy(reader,"dxy_1")
            ,dxyError(reader,"dxyError_1")
            ,eta(reader,"eta_1")
            ,rapidity(reader,"rapidity_1")
            ,phi(reader,"phi_1")
            ,vx(reader,"vx_1")
            ,vy(reader,"vy_1")
            ,vz(reader,"vz_1")
            ,chi2norm(reader,"chi2norm_1")
            ,nhits(reader,"nhits_1")
            ,charge(reader,"charge_1")
            {
                mass_low = 0.0;
                mass_high = 100.05;
                mass_numBins = 50;
                h_mass_pt = 0;
                h_dedx = 0;
            }

        virtual ~PhiTreeProof() {}
        virtual void Init(TTree *tree);
        virtual void SlaveBegin(TTree *tree);
        virtual bool Process(Long64_t entry);
        virtual void Terminate();
        virtual int Version() const {return 2;}

        static void CombinatorialMass(std::vector<kaon_proof> PKp, std::vector<kaon_proof> PKm, TH1D* h_mass_);
        static void CombinatorialMass2D(std::vector<kaon_proof> PKp, std::vector<kaon_proof> PKm, TH2D* h_mass_, bool doRapCut_ = true);
        void Analyze(double cut_dcaz_, double cut_dcaxy_, double cut_eta_, double cut_rapidity_, int cut_nhits_, std::string cut_dedx_, TH2D* h_dedx_, TH2D* h_mass_pt_);

        template <typename T, size_t N>
            T* begin(T(&arr)[N]) {return &arr[0];}
        template <typename T, size_t N>
            T* end(T(&arr)[N]) {return &arr[0] + N;}

        ClassDef(PhiTreeProof,0);

};
