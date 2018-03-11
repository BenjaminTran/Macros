#include "PhiTreeProofSparse.h"

void PhiTreeProofSparse::SparseVarSetUp(std::vector<float> &varBins_, int nBins_, float minBin_, float maxBin_)
{
    double interval = (maxBin_ - minBin_)/nBins_;

    for(int i=0; i<nBins_; ++i)
    {
        varBins_.emplace_back(minBin_ + i*interval);
    }

    varBins_.emplace_back(maxBin_);
}

//std::vector<double> PhiTreeProofSparse::CombinatorialMassSparse(std::vector<kaon_proof> PKp, std::vector<kaon_proof> PKm)
//{
    //for(kaon_proof Pkp : PKp)
    //{
        //for(kaon_proof Pkm : PKm)
        //{
            //TVector3 dau1p(Pkp.px,Pkp.py,Pkp.pz);
            //TVector3 dau2p(Pkm.px,Pkm.py,Pkm.pz);
            //TVector3 dauPsum(dau1p + dau2p);
            //double mass = sqrt(TMath::Power(Pkp.energy+Pkm.energy,2) - dauPsum.Mag2());
            //Mass.push_back(mass);
        //}
    //}
    //return Mass;
//}

void PhiTreeProofSparse::Init(TTree *tree)
{
    //cout << "Init" << endl;
    reader.SetTree(tree);
}

void PhiTreeProofSparse::SlaveBegin(TTree *tree)
{
    //cout << "Begin" << endl;
    std::vector<std::string> AxisLabels =
    {
        "Mass",
        "Pt1",
        "Pt2",
        "Eta1",
        "Eta2",
        "Rapidity1",
        "Rapidity2",
        "Nhits1",
        "Nhits2",
        "DCAz1",
        "DCAz2",
        "DCAxy1",
        "DCAxy2",
        //"mult"
    };

    const int nDim = AxisLabels.size();
    std::vector<int> nBins      = {80,80,200,200,24,24,30,30,35,35,20,20,20,20};/* ,800};*/
    std::vector<double> minBins = {0.98,0.98,0,0,0,0,0,0,0,0,0,0,0,0};/*  ,0};*/
    std::vector<double> maxBins = {1.06,1.06,20,20,2.4,2.4,3.0,3.0,35,35,10.0,10.0,10.0,10.0};/* ,800};*/

    std::vector<std::vector<float> > VarBinContainers(nDim);

    //for(int i=0; i<nDim; ++i)
    //{
        //std::vector<float> tmp;
        //VarBinContainers.push_back(tmp);
    //}

    h = new THnSparseF("hsparse", "Cuts", nDim, &nBins[0], &minBins[0], &maxBins[0]);

    for(int i=0; i<VarBinContainers.size(); ++i)
    {
        SparseVarSetUp(VarBinContainers[i], nBins[i], minBins[i], maxBins[i]);
        h->GetAxis(i)->SetName(AxisLabels[i].c_str());
        h->GetAxis(i)->Set(nBins[i],&(VarBinContainers[i])[0]);
    }

    h_nEvt = new TH1D("nEvt","Events",10,0,10);
    //cout << "Ready to initialize sparse" << endl;

    GetOutputList()->Add(h_nEvt);
    GetOutputList()->Add(h);
    TH1::SetDefaultSumw2();
    //cout << "Initialization done" << endl;
}

bool PhiTreeProofSparse::Process(Long64_t entry)
{
    h_nEvt->Fill(1);
    reader.SetEntry(entry);

        double momentum_ = -999;
        double dedx_     = -999;
        double pt_       = -999;
        double px_       = -999;
        double py_       = -999;
        double pz_       = -999;
        double eta_      = -999;
        double rapidity_ = -999;
        double energy_   = -999;
        double charge_   = -999;
        double nhits_    = -999;
        double mult_     = 0;
        double DCA_z     = -999;
        double DCA_xy    = -999;

        for(int i=0; i<pt.GetSize(); ++i)
        {
            momentum_ = momentum[i];
            dedx_     = dedx[i];
            pt_       = pt[i];
            px_       = px[i];
            py_       = py[i];
            pz_       = pz[i];
            eta_      = fabs(eta[i]);
            rapidity_ = fabs(rapidity[i]);
            energy_   = energy[i];
            charge_   = charge[i];
            nhits_    = nhits[i];

            double dzvtx = dz[i];
            double dxyvtx = dxy[i];
            double dzerror = dzError[i];
            double dxyerror = dxyError[i];

            DCA_z = fabs(dzvtx/dzerror);
            DCA_xy = fabs(dxyvtx/dxyerror);

            if(dedx_ <= f_Dedx_bot.Eval(momentum_) || dedx_ >= f_Dedx_top.Eval(momentum_)) continue;
            if(fabs(ptError[i]/pt[i]) > 0.10) continue;

            kaon_proof pk(momentum_, pt_, px_, py_, pz_, dedx_, energy_, eta_, rapidity_, DCA_z, DCA_xy, charge_, nhits_);

            //positive kaons
            if(pk.charge == 1)
                Pkp.push_back(pk);

            //negative kaons
            if(pk.charge == -1)
                Pkm.push_back(pk);
        }

        for(kaon_proof pkp : Pkp)
        {
            for(kaon_proof pkm : Pkm)
            {
                TVector3 dau1p(pkp.px,pkp.py,pkp.pz);
                TVector3 dau2p(pkm.px,pkm.py,pkm.pz);
                TVector3 dauPsum(dau1p + dau2p);
                double mass = sqrt(TMath::Power(pkp.energy+pkm.energy,2) - dauPsum.Mag2());
                std::vector<double> value = {mass,pkp.pt,pkm.pt,pkp.eta,pkm.eta,pkp.rapidity,pkm.rapidity,pkp.nhits,pkm.nhits,pkp.DCAz,pkm.DCAz,pkp.DCAxy,pkm.DCAxy};
                h->Fill(&value[0]);
            }
        }

        Pkp.clear();
        Pkm.clear();

    return true;
}

void PhiTreeProofSparse::Terminate()
{

    out = new TFile("Sparse_v1.root","RECREATE");
    h_nEvt->Write();
    h->Write();
    out->Close();
}
