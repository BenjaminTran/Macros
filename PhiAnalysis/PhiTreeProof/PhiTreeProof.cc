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

void PhiTreeProof::CombinatorialMass2D(std::vector<kaon_proof> PKp, std::vector<kaon_proof> PKm, TH2D* h_mass_, bool doRapCut_)
{
    for(kaon_proof Pkp : PKp)
    {
        for(kaon_proof Pkm : PKm)
        {
            TVector3 dau1p(Pkp.px,Pkp.py,Pkp.pz);
            TVector3 dau2p(Pkm.px,Pkm.py,Pkm.pz);
            TVector3 dauPsum(dau1p + dau2p);
            TLorentzVector phiLV(dauPsum,Pkp.energy + Pkm.energy);
            double rap = phiLV.Rapidity();
            if(doRapCut_)
            {
                if(fabs(rap) > 1.0)
                    continue;
            }
            double mass = sqrt(TMath::Power(Pkp.energy+Pkm.energy,2) - dauPsum.Mag2());
            h_mass_->Fill(mass,dauPsum.Perp());
        }
    }
    cout << "Combinatorial" << endl;
}

void PhiTreeProof::Analyze(double cut_dcaz_, double cut_dcaxy_, double cut_eta_, double cut_rapidity_, int cut_nhits_, std::string cut_dedx_, TH2D* h_dedx_, TH2D* h_mass_pt_)
{
    if(cut_dedx_ == "none")
        doDedx = false;
    else
        doDedx = true;

    for(unsigned i=0; i<pt.GetSize(); ++i)
    {
        if(fabs(ptError[i]/pt[i])>0.10) continue;
        if(fabs(dz[i]/dzError[i]) > cut_dcaz_) continue;
        if(fabs(dxy[i]/dxyError[i]) > cut_dcaxy_) continue;
        if(fabs(eta[i]) > cut_eta_) continue;
        //if(fabs(rapidity[i]) > cut_rapidity_) continue;
        if(nhits[i] < cut_nhits_) continue;

        if(doDedx)
        {
            if(dedx[i] >= f_Dedx_bot[cut_dedx_]->Eval(momentum[i]) && dedx[i] <= f_Dedx_top[cut_dedx_]->Eval(momentum[i]))
            {
                h_dedx_->Fill(momentum[i],dedx[i]);
                kaon_proof pk(momentum[i], pt[i], px[i], py[i], pz[i], dedx[i], energy[i], charge[i], nhits[i]);

                //positive kaons
                if(pk.charge == 1)
                    Pkp.push_back(pk);

                //negative kaons
                if(pk.charge == -1)
                    Pkm.push_back(pk);
            }
        }
        else
        {
            kaon_proof pk(momentum[i], pt[i], px[i], py[i], pz[i], dedx[i], energy[i], charge[i], nhits[i]);

            //positive kaons
            if(pk.charge == 1)
                Pkp.push_back(pk);

            //negative kaons
            if(pk.charge == -1)
                Pkm.push_back(pk);
        }
    }

    //CombinatorialMass2D(Pkp, Pkm, h_mass_pt_,false);
    Pkp.clear();
    Pkm.clear();

}

void PhiTreeProof::Init(TTree *tree)
{
    reader.SetTree(tree);
}

void PhiTreeProof::SlaveBegin(TTree *tree)
{
    double Acut_dcaz     [] = {1.0    , 1.0     , 3.0};
    double Acut_dcaxy    [] = {1.0    , 1.0     , 3.0};
    double Acut_eta      [] = {2.4    , 2.4     , 2.4};
    double Acut_rapidity [] = {1.0    , 1.0     , 1.0};
    int Acut_nhits       [] = {5      , 5       , 5};
    std::string Acut_dedx[] = {"none" , "tight" , "tight"};

    cut_dcaz.assign(begin(Acut_dcaz), end(Acut_dcaz));
    cut_dcaxy.assign(begin(Acut_dcaxy), end(Acut_dcaxy));
    cut_eta.assign(begin(Acut_eta), end(Acut_eta));
    cut_rapidity.assign(begin(Acut_rapidity), end(Acut_rapidity));
    cut_nhits.assign(begin(Acut_nhits), end(Acut_nhits));
    cut_dedx.assign(begin(Acut_dedx), end(Acut_dedx));

    f_Dedx_bot["tight"] = new TF1("f_Dedx_bot_tight","0.55*(TMath::Power(1.15/x,2) - 1.7*TMath::Power(0.6/x,1)) + 2.7",0.01,100);
    f_Dedx_bot["loose"] = new TF1("f_Dedx_bot_loose","0.55*(TMath::Power(1.15/x,2) - 1.7*TMath::Power(0.6/x,1)) + 2.5",0.01,100);
    f_Dedx_top["tight"] = new TF1("f_Dedx_top_tight","0.55*(TMath::Power(1.62/x,2) - 2.95*TMath::Power(0.6/x,1)) + 3.3",0.01,100);
    f_Dedx_top["loose"] = new TF1("f_Dedx_top_loose","0.55*(TMath::Power(1.62/x,2) - 2*TMath::Power(0.6/x,1)) + 3.6",0.01,100);

    //Initialize vector histograms

    num_analyze = cut_dcaz.size();
    for(unsigned i=0; i<cut_dcaz.size(); ++i)
    {
        //TH2D* tmp_mass_pt = new TH2D(Form("mass_pt_%d",i),Form("Combinatorial Mass v. Pt - %d",i),mass_numBins,mass_low,mass_high,200,0,20);
        //TH2D* tmp_dedx = new TH2D(Form("dedx_%d",i),Form("DeDx - %d",i),200,0,5,250,0,15);

        v_hmasspt.emplace_back( new TH2D(Form("mass_pt_%d",i),Form("Combinatorial Mass v. Pt - %d",i),mass_numBins,mass_low,mass_high,200,0,20));
        v_hdedx.emplace_back( new TH2D(Form("dedx_%d",i),Form("DeDx - %d",i),200,0,5,250,0,15));
    }

    //h_mass_pt = new TH2D*[3];
    //h_dedx = new TH2D*[3];
    //for(unsigned i=0; i<cut_dcaz.size(); ++i)
    //{
        //h_mass_pt[i] = new TH2D(Form("mass_pt_%d",i),Form("Combinatorial Mass v. Pt - %d",i),mass_numBins,mass_low,mass_high,200,0,20);
        //h_dedx[i] = new TH2D(Form("dedx_%d",i),Form("DeDx - %d",i),200,0,5,250,0,15);
    //}

    h_nEvt = new TH1D("nEvt","Events",10,0,10);
    GetOutputList()->Add(h_nEvt);

    for(unsigned i=0; i<v_hmasspt.size(); ++i)
    {
        GetOutputList()->Add(v_hmasspt[i]);
        GetOutputList()->Add(v_hdedx[i]);
    }
    TH1::SetDefaultSumw2();
}

bool PhiTreeProof::Process(Long64_t entry)
{
    h_nEvt->Fill(1);
    reader.SetEntry(entry);

    //1 -- No dedx cut
    PhiTreeProof::Analyze(cut_dcaz[0],cut_dcaxy[0],cut_eta[0],cut_rapidity[0],cut_nhits[0],cut_dedx[0],v_hdedx[0],v_hmasspt[0]);

    //2 -- With dedx cut to compare with (1)
    PhiTreeProof::Analyze(cut_dcaz[1],cut_dcaxy[1],cut_eta[1],cut_rapidity[1],cut_nhits[1],cut_dedx[1],v_hdedx[1],v_hmasspt[1]);

    //3 -- dca 3.0 to compare with dca 1.0. Compare with (2)
    PhiTreeProof::Analyze(cut_dcaz[2],cut_dcaxy[2],cut_eta[2],cut_rapidity[2],cut_nhits[2],cut_dedx[2],v_hdedx[2],v_hmasspt[2]);

    return true;
}

void PhiTreeProof::Terminate()
{
    TCanvas* c = new TCanvas("Output","Output",1200,900);
    h_nEvt->Draw();

/*
 * Masking the file name. Multiply cut values by 10 unless already an integer. Construct name in this order:
 * Mass_pt_dedx_dcaz_dcaxy_eta_rap_nhits.root
 */
    std::ostringstream os;
    //os << "outputFiles/Mass_pt_" + cut_dedx + "_" << cut_dcaz*10 << "_" << cut_dcaxy*10 << "_" << cut_eta*10 << "_" << cut_rapidity*10 << "_" << cut_nhits << ".root";
    os << "Phi_MassPt_v3.root";
    out = new TFile(os.str().c_str(),"RECREATE");
    h_nEvt->Write();
    //v_hmasspt->at(0)->Write();
    TH2D* c_mass[10];
    TH2D* c_dedx[10];
    for(unsigned i=0; i<num_analyze; ++i)
    {
        c_mass[i] = dynamic_cast<TH2D*> (GetOutputList()->FindObject(Form("mass_pt_%d",i)));
        c_dedx[i] = dynamic_cast<TH2D*> (GetOutputList()->FindObject(Form("dedx_%d",i)));

        c_mass[i]->Write("mass");
        c_dedx[i]->Write("dedx");
    }

    out->Close();
}
