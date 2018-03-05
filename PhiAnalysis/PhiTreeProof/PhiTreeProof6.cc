#include "PhiTreeProof6.h"

void PhiTreeProof6::CombinatorialMass(std::vector<kaon_proof> PKp, std::vector<kaon_proof> PKm, TH1D* h_mass_)
{
    for(kaon_proof Pkp : PKp)
    {
        for(kaon_proof Pkm : PKm)
        {
            TVector3 dau1p(Pkp.px,Pkp.py,Pkp.pz);
            TVector3 dau2p(Pkm.px,Pkm.py,Pkm.pz);
            TVector3 dauPsum(dau1p + dau2p);
            double mass = sqrt(TMath::Power(Pkp.energy+Pkm.energy,2) - dauPsum.Mag2());
            h_mass_->Fill(mass);
        }
    }
}

void PhiTreeProof6::Init(TTree *tree)
{
    reader.SetTree(tree);
}

void PhiTreeProof6::SlaveBegin(TTree *tree)
{
    h_nEvt = new TH1D("nEvt","Events",10,0,10);
    h_mass = new TH1D("mass","Combinatorial Mass",50,1,1.05);
    h_dedx = new TH2D("dedx","DeDx",200,0,5,250,0,15);
    GetOutputList()->Add(h_nEvt);
    //GetOutputList()->Add(out);
    GetOutputList()->Add(h_mass);
    GetOutputList()->Add(h_dedx);
    TH1::SetDefaultSumw2();
}

bool PhiTreeProof6::Process(Long64_t entry)
{
    h_nEvt->Fill(1);
    reader.SetEntry(entry);

    for(int i=0; i<pt.GetSize(); ++i)
    {
        if(fabs(ptError[i]/pt[i])>0.10) continue;
        if(fabs(dz[i]/dzError[i]) > 1.0) continue;
        if(fabs(dxy[i]/dxyError[i]) > 1.0) continue;
        if(fabs(eta[i]) > 2.4) continue;
        if(fabs(rapidity[i]) > 1.0) continue;
        if(nhits[i] < 10) continue;

        if(dedx[i] >= f_Dedx_bot.Eval(momentum[i]) && dedx[i] <= f_Dedx_top.Eval(momentum[i]))
        {
            h_dedx->Fill(momentum[i],dedx[i]);
            kaon_proof pk(momentum[i], pt[i], px[i], py[i], pz[i], dedx[i], energy[i], charge[i], nhits[i]);

            //positive kaons
            if(pk.charge == 1)
                Pkp.push_back(pk);

            //negative kaons
            if(pk.charge == -1)
                Pkm.push_back(pk);
        }
    }

    CombinatorialMass(Pkp, Pkm, h_mass);
    Pkp.clear();
    Pkm.clear();

    return true;
}

void PhiTreeProof6::Terminate()
{
    TCanvas* c = new TCanvas("Output","Output",1200,900);
    c->Divide(2,2);
    c->cd(1);
    h_mass->SetMarkerStyle(20);
    h_mass->SetMarkerSize(1.2);
    h_mass->Draw();
    c->cd(2);
    h_nEvt->Draw();
    c->cd(3);
    h_dedx->Draw("COL SCAT");

    out = new TFile("outputFiles/Histograms_tight_midRap_nhits10.root","RECREATE");
    c->Write("OutputCanvas");
    h_mass->Write();
    h_nEvt->Write();
    h_dedx->Write();
    //h_nEvt->Write("numEvents");
    //h_dedx->Write("DeDx");
    //h_mass->Write("mass");
    //out->Close();
}
