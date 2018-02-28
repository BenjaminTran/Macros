#include "interface/PhiTreeReader.h"

#include "TSelector.h"

class PhiTreeProof : public TSelector {
    public:
        TTreeReader reader;
        //TTreeReaderArray<float> momentum;
        //TTreeReaderArray<float> px;
        //TTreeReaderArray<float> py;
        //TTreeReaderArray<float> pz;
        TTreeReaderArray<float> pt;
        //TTreeReaderArray<float> ptError;
        //TTreeReaderArray<float> energy;
        //TTreeReaderArray<float> dedx;
        //TTreeReaderArray<float> dz;
        //TTreeReaderArray<float> dzError;
        //TTreeReaderArray<float> dxy;
        //TTreeReaderArray<float> dxyError;
        //TTreeReaderArray<float> eta;
        //TTreeReaderArray<float> rapidity;
        //TTreeReaderArray<float> phi;
        //TTreeReaderArray<float> vx;
        //TTreeReaderArray<float> vy;
        //TTreeReaderArray<float> vz;
        //TTreeReaderArray<float> chi2norm;
        //TTreeReaderArray<float> nhits;
        //TTreeReaderArray<float> charge;
        //
        TH1D* h_pt;

        PhiTreeProof(TTree * = 0) : pt(reader,"pt") { }
        virtual ~PhiTreeProof() {}
        virtual void Init(TTree *tree);
        virtual void SlaveBegin(TTree *tree);
        virtual bool Process(Long64_t entry);
        virtual void Terminate();
        virtual int Version() const {return 2;}

        ClassDef(PhiTreeProof,0);
};
