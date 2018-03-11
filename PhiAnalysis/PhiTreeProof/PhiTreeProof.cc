#include "PhiTreeProof.h"

void PhiTreeProof::CombinatorialMass(std::vector<kaon_proof> PKp, std::vector<kaon_proof> PKm, TH1D* h_mass_)
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

void PhiTreeProof::CombinatorialMass2D(std::vector<kaon_proof> PKp, std::vector<kaon_proof> PKm, TH2D* h_mass_)
{
    for(kaon_proof Pkp : PKp)
    {
        for(kaon_proof Pkm : PKm)
        {
            TVector3 dau1p(Pkp.px,Pkp.py,Pkp.pz);
            TVector3 dau2p(Pkm.px,Pkm.py,Pkm.pz);
            TVector3 dauPsum(dau1p + dau2p);
            double mass = sqrt(TMath::Power(Pkp.energy+Pkm.energy,2) - dauPsum.Mag2());
            h_mass_->Fill(mass,dauPsum.Perp());
        }
    }
}

void PhiTreeProof::Init(TTree *tree)
{
    reader.SetTree(tree);
}

void PhiTreeProof::SlaveBegin(TTree *tree)
{
    if(cut_dedx == "tight")
    {
        f_Dedx_bot = new TF1("f_Dedx_bot","0.55*(TMath::Power(1.15/x,2) - 1.7*TMath::Power(0.6/x,1)) + 2.7",0.01,100);
        f_Dedx_top = new TF1("f_Dedx_top","0.55*(TMath::Power(1.62/x,2) - 2.95*TMath::Power(0.6/x,1)) + 3.3",0.01,100);
    }
    else if(cut_dedx == "loose")
    {
        f_Dedx_bot = new TF1("f_Dedx_bot","0.55*(TMath::Power(1.15/x,2) - 1.7*TMath::Power(0.6/x,1)) + 2.5",0.01,100);
        f_Dedx_top = new TF1("f_Dedx_top","0.55*(TMath::Power(1.62/x,2) - 2*TMath::Power(0.6/x,1)) + 3.6",0.01,100);
    }
    else
        cout << "unknown dedx function selected" << endl;

    h_nEvt = new TH1D("nEvt","Events",10,0,10);
    h_mass_pt = new TH2D("mass","Combinatorial Mass",mass_numBins,mass_low,mass_high,200,0,20);
    h_dedx = new TH2D("dedx","DeDx",200,0,5,250,0,15);
    GetOutputList()->Add(h_nEvt);
    GetOutputList()->Add(h_mass_pt);
    GetOutputList()->Add(h_dedx);
    TH1::SetDefaultSumw2();
}

bool PhiTreeProof::Process(Long64_t entry)
{
    h_nEvt->Fill(1);
    reader.SetEntry(entry);

    for(int i=0; i<pt.GetSize(); ++i)
    {
        if(fabs(ptError[i]/pt[i])>0.10) continue;
        if(fabs(dz[i]/dzError[i]) > cut_dcaz) continue;
        if(fabs(dxy[i]/dxyError[i]) > cut_dcaxy) continue;
        if(fabs(eta[i]) > cut_eta) continue;
        if(fabs(rapidity[i]) > cut_rapidity) continue;
        if(nhits[i] < cut_nhits) continue;

        if(dedx[i] >= f_Dedx_bot->Eval(momentum[i]) && dedx[i] <= f_Dedx_top->Eval(momentum[i]))
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

    CombinatorialMass2D(Pkp, Pkm, h_mass_pt);
    Pkp.clear();
    Pkm.clear();

    return true;
}

void PhiTreeProof::Terminate()
{
    TCanvas* c = new TCanvas("Output","Output",1200,900);
    c->Divide(2,1);
    c->cd(1);
    h_nEvt->Draw();
    c->cd(2);
    h_dedx->Draw("COL SCAT");

/*
 * Masking the file name. Multiply cut values by 10 unless already an integer. Construct name in this order:
 * Mass_pt_dedx_dcaz_dcaxy_eta_rap_nhits.root
 */
    std::ostringstream os;
    os << "outputFiles/Mass_pt_" + cut_dedx + "_" << cut_dcaz*10 << "_" << cut_dcaxy*10 << "_" << cut_eta*10 << "_" << cut_rapidity*10 << "_" << cut_nhits << ".root";
    out = new TFile(os.str().c_str(),"RECREATE");
    h_mass_pt->Write();
    h_nEvt->Write();
    h_dedx->Write();
    //h_nEvt->Write("numEvents");
    //h_dedx->Write("DeDx");
    //h_mass->Write("mass");
    //out->Close();
}
