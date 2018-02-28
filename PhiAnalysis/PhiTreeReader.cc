#include "interface/PhiTreeReader.h"

const double kaonMass = 0.493677;

PhiTreeReader::PhiTreeReader(long int event_start_, long int event_stop_ , int check_interval_, std::string filename_)
{
    h_DeDx = new TH2D("DeDx","DeDx",200,0,5,250,0,15);
    h_mass = new TH1D("Mass","Mass",40,1,1.04);
    h_nEvt = new TH1I("nEvt","nEvt",5,0,5);
    event_start = event_start_;
    event_stop = event_stop_;
    check_interval = check_interval_;
    filename = filename_;
}

PhiTreeReader::~PhiTreeReader() {}

bool PhiTreeReader::CheckValue(ROOT::Internal::TTreeReaderValueBase& value) {
    if (value.GetSetupStatus() < 0) {
        std::cerr << "Error " << value.GetSetupStatus()
            << "setting up reader for " << value.GetBranchName() << '\n';
        return false;
    }
    return true;
}

TTreeReaderValue<float> PhiTreeReader::SetValues(TTreeReader& reader, std::string branchName ){
    TTreeReaderValue<float> ReaderValue(reader, branchName.c_str());
    return ReaderValue;
}

TTreeReaderValue<int> PhiTreeReader::SetValuesInt(TTreeReader& reader, std::string branchName ){
    TTreeReaderValue<int> ReaderValue(reader, branchName.c_str());
    return ReaderValue;
}

void PhiTreeReader::CombinatorialMass(std::vector<kaon> PKp, std::vector<kaon> PKm, TH1D* h_mass_)
{
    for(kaon Pkp : PKp)
    {
        for(kaon Pkm : PKm)
        {
            TVector3 dau1p(Pkp.px,Pkp.py,Pkp.pz);
            TVector3 dau2p(Pkm.px,Pkm.py,Pkm.pz);
            TVector3 dauPsum(dau1p + dau2p);
            double mass = sqrt(TMath::Power(Pkp.energy+Pkm.energy,2) - dauPsum.Mag2());
            h_mass_->Fill(mass);
        }
    }
}

std::vector<double> PhiTreeReader::CombinatorialMassSparse(std::vector<kaon> PKp, std::vector<kaon> PKm)
{
    std::vector<double> Mass;
    for(kaon Pkp : PKp)
    {
        for(kaon Pkm : PKm)
        {
            TVector3 dau1p(Pkp.px,Pkp.py,Pkp.pz);
            TVector3 dau2p(Pkm.px,Pkm.py,Pkm.pz);
            TVector3 dauPsum(dau1p + dau2p);
            double mass = sqrt(TMath::Power(Pkp.energy+Pkm.energy,2) - dauPsum.Mag2());
            Mass.push_back(mass);
        }
    }
    return Mass;
}

void PhiTreeReader::SparseVarSetUp(std::vector<float> &varBins_, int nBins_, float minBin_, float maxBin_)
{
    double interval = (maxBin_ - minBin_)/nBins_;

    for(int i=0; i<nBins_; ++i)
    {
        varBins_.emplace_back(minBin_ + i*interval);
    }

    varBins_.emplace_back(maxBin_);
}


bool PhiTreeReader::ReadTreeV2()
{
    ROOT::EnableImplicitMT();
    TH1::SetDefaultSumw2();

    gBenchmark->Start("ReadTree");

    std::chrono::time_point<std::chrono::system_clock> start;
    start = std::chrono::system_clock::now();
    long int secondsElapsed = 0;

    std::string fileName = "/Volumes/MacHD/Users/blt1/research/Macros/PhiAnalysis/TreeReader/rootFiles/Test.root";

    TFile* f1 = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/Phi/Phi_TrackTree_PD1_v2_021718.root");
    std::string path_to_tree_ = "PhiTree/TrackTree";
    //TFile* f1 = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/Phi/Phi_Tree_v1_99.root");

    TTreeReader reader;
    //std::map<std::string,AnyReaderValue> BranchMap = PhiTreeReader::SetTTreeReaderValue(f1,path_to_tree_,reader);
    reader.SetTree("PhiTree/TrackTree",f1);

    TTreeReaderValue<float> momentum(reader,"momentum");
    TTreeReaderValue<float> px(reader,"px");
    TTreeReaderValue<float> py(reader,"py");
    TTreeReaderValue<float> pz(reader,"pz");
    TTreeReaderValue<float> pt(reader,"pt");
    TTreeReaderValue<float> ptError(reader,"ptError");
    TTreeReaderValue<float> energy(reader,"energy");
    TTreeReaderValue<float> dedx(reader,"dedx");
    TTreeReaderValue<float> dz(reader,"dz");
    TTreeReaderValue<float> dzError(reader,"dzError");
    TTreeReaderValue<float> dxy(reader,"dxy");
    TTreeReaderValue<float> dxyError(reader,"dxyError");
    TTreeReaderValue<float> eta(reader,"eta");
    TTreeReaderValue<float> phi(reader,"phi");
    TTreeReaderValue<float> vx(reader,"vx");
    TTreeReaderValue<float> vy(reader,"vy");
    TTreeReaderValue<float> vz(reader,"vz");
    TTreeReaderValue<float> chi2norm(reader,"chi2norm");
    TTreeReaderValue<int> nhits(reader,"nhits");
    TTreeReaderValue<int> charge(reader,"charge");
    TTreeReaderValue<int> mult(reader,"mult");
    TTreeReaderValue<int> multRaw(reader,"multRaw");

    int counter = -1;
    int progress = 0;
    long int event_num = 0;
    //TF1* f_Dedx_bot = new TF1("f_Dedx_bot","0.55*(TMath::Power(1.15/x,2) - 1.7*TMath::Power(0.6/x,1)) + 2.7",0.01,5);
    TF1* f_Dedx_bot = new TF1("f_Dedx_bot","0.55*(TMath::Power(1.15/x,2) - 1.7*TMath::Power(0.6/x,1)) + 2.5",0.01,5);
    //TF1* f_Dedx_bot2 = new TF1("f_Dedx","0.5*(TMath::Power(0.5/x,4) - 1*TMath::Power(0.50/x,2)) + 2.7",0.01,5);
    //TF1* f_Dedx_top = new TF1("f_Dedx_top","0.55*(TMath::Power(1.62/x,2) - 2.95*TMath::Power(0.6/x,1)) + 3.3",0.01,5); //magenta
    TF1* f_Dedx_top = new TF1("f_Dedx_top","0.55*(TMath::Power(1.62/x,2) - 2*TMath::Power(0.6/x,1)) + 3.6",0.01,5); //magenta
    //TF1* f_Dedx_top = new TF1("f_Dedx_top","x",0.1,5); //magenta

    //TF1* f_Dedx_bot = new TF1("f_Dedx_bot","0.55*(TMath::Power(1.62/x,2) - 2*TMath::Power(0.6/x,1)) + 3.6",0.1,5);
    //TF1* f_Dedx_top = new TF1("f_Dedx_top","0.5*(TMath::Power(0.5/x,4) - 1*TMath::Power(0.50/x,2)) + 2.7",0.1,5); //magenta
    //
    int mult_tmp = 0;
    int mult_check = 0;
    bool check = false;

    //long int range_start = 1e8+1;
    //long int range_end = 3e8;
    //reader.SetEntriesRange(range_start,range_end);
    while(reader.Next())
    {
        counter++;
        //progress++;
        if(event_num % check_interval == 0  && event_num != 0 && check)
        {
            //int secondsElapsed = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start).count();
            cout << (double)event_num/event_stop*100 << " \% completed" << endl;
            //cout << secondsElapsed << " seconds" << endl;
            check = false;
        }
        if(!CheckValue(momentum)) return false;
        if(!CheckValue(px)) return false;
        if(!CheckValue(py)) return false;
        if(!CheckValue(pz)) return false;
        if(!CheckValue(pt)) return false;
        if(!CheckValue(ptError)) return false;
        if(!CheckValue(energy)) return false;
        if(!CheckValue(dedx)) return false;
        if(!CheckValue(dz)) return false;
        if(!CheckValue(dzError)) return false;
        if(!CheckValue(dxy)) return false;
        if(!CheckValue(dxyError)) return false;
        if(!CheckValue(eta)) return false;
        if(!CheckValue(phi)) return false;
        if(!CheckValue(vx)) return false;
        if(!CheckValue(vy)) return false;
        if(!CheckValue(vz)) return false;
        if(!CheckValue(nhits)) return false;
        if(!CheckValue(charge)) return false;
        if(!CheckValue(mult)) return false;
        if(!CheckValue(multRaw)) return false;

        mult_check = *multRaw;

        if(counter == 0)
        {
            mult_tmp = mult_check;
        }

        /*If the multiplicity of the current particle does not equal the
         * multiplicty information of the previous particle then we are on a new event
         * and so start combinatorial mass.
         */
        if(mult_tmp != mult_check)
        {
            //cout << "Begin combinatorial mass" << endl;
            //Set counter to -1 because the increment at the end of the while loop
            //will bring it back to 0
            counter = -1;
            h_nEvt->Fill(1);
            event_num++;
            check = true;
            CombinatorialMass(Pkp, Pkm, h_mass);
            Pkp.clear();
            Pkm.clear();
        }
        if(event_num == event_stop)
        {
            break;
        }

        double dzvtx = *dz;
        double dxyvtx = *dxy;
        double dzerror = *dzError;
        double dxyerror = *dxyError;

        //if(*momentum > 1) continue;
        if(fabs(*ptError/ *pt)>0.10) continue;
        if(fabs(*dz/ *dzError) > 1.0) continue;
        if(fabs(*dxy/ *dxyError) > 1.0) continue;
        if(fabs(*eta) > 2.4) continue;
        if(*nhits < 11) continue;
        //if(*pt < 0.7) continue;
        //if(*chi2norm > 1) continue;
        //if(*energy < 1.2) continue;

        double val_momentum = *momentum;
        double val_dedx = *dedx;
        double val_pt = *pt;

        if(val_dedx >= f_Dedx_bot->Eval(val_momentum) && val_dedx <= f_Dedx_top->Eval(val_momentum))
        {
            h_DeDx->Fill(val_momentum,val_dedx);
            PhiTreeReader::kaon pk(*momentum, *pt, *px, *py, *pz, *dedx, *energy, *charge, *nhits);

            //positive kaons
            if(pk.charge == 1)
                Pkp.push_back(pk);

            //negative kaons
            if(pk.charge == -1)
                Pkm.push_back(pk);
        }

    }

    TCanvas* c_DeDx = new TCanvas("c_DeDx","c_DeDx",600,600);
    c_DeDx->cd();
    //f_Dedx_bot2->SetLineColor(3);
    h_DeDx->Draw("SCAT");
    f_Dedx_bot->Draw("SAME");
    //f_Dedx_bot2->Draw("SAME");
    f_Dedx_top->Draw("SAME");

    TCanvas* c_Mass = new TCanvas("mass","mass",600,600);
    c_Mass->cd();
    h_mass->Draw();

    f1->Close();

    TFile* out = new TFile(fileName.c_str(),"RECREATE");
    h_nEvt->Write();
    h_mass->Write();
    h_DeDx->Write();

    out->Close();

    //c_Mass->Print("Image/CombMassLooseBand_nhits17_6e7.pdf");
    //c_Mass->Print("Image/CombMassLooseBand_nhits17_6e7.png");
    //c_DeDx->Print("Image/DedxLooseBand_nhits17_6e7.pdf");
    //c_DeDx->Print("Image/DedxLooseBand_nhits17_6e7.png");

    gBenchmark->Show("ReadTree");
    return true;
}

bool PhiTreeReader::ReadTreeV3()
{
    //ROOT::EnableImplicitMT();
    TH1::SetDefaultSumw2();

    gBenchmark->Start("ReadTreeV3");
    gBenchmark->Start("Total");

    std::chrono::time_point<std::chrono::system_clock> start;
    start = std::chrono::system_clock::now();
    long int secondsElapsed = 0;

    TFile* f1 = TFile::Open("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/Trees/Phi/Phi_TrackTree_PD1_v3_files1-80_022518_JL21.root");
    std::string path_to_tree_ = "PhiTree/TrackTree";

    TTreeReader reader;
    reader.SetTree("PhiTree/TrackTree",f1);

    TTreeReaderArray<float> momentum(reader,"momentum");
    TTreeReaderArray<float> px(reader,"px");
    TTreeReaderArray<float> py(reader,"py");
    TTreeReaderArray<float> pz(reader,"pz");
    TTreeReaderArray<float> pt(reader,"pt");
    TTreeReaderArray<float> ptError(reader,"ptError");
    TTreeReaderArray<float> energy(reader,"energy");
    TTreeReaderArray<float> dedx(reader,"dedx");
    TTreeReaderArray<float> dz(reader,"dz");
    TTreeReaderArray<float> dzError(reader,"dzError");
    TTreeReaderArray<float> dxy(reader,"dxy");
    TTreeReaderArray<float> dxyError(reader,"dxyError");
    TTreeReaderArray<float> eta(reader,"eta");
    TTreeReaderArray<float> rapidity(reader,"rapidity");
    TTreeReaderArray<float> phi(reader,"phi");
    TTreeReaderArray<float> vx(reader,"vx");
    TTreeReaderArray<float> vy(reader,"vy");
    TTreeReaderArray<float> vz(reader,"vz");
    TTreeReaderArray<float> chi2norm(reader,"chi2norm");
    TTreeReaderArray<float> nhits(reader,"nhits");
    TTreeReaderArray<float> charge(reader,"charge");
    //TTreeReaderArray<int> mult(reader,"mult");
    //TTreeReaderArray<int> multRaw(reader,"multRaw");

    long int event_num = 0;
    TF1* f_Dedx_bot = new TF1("f_Dedx_bot","0.55*(TMath::Power(1.15/x,2) - 1.7*TMath::Power(0.6/x,1)) + 2.7",0.01,100);
    //TF1* f_Dedx_bot = new TF1("f_Dedx_bot","0.55*(TMath::Power(1.15/x,2) - 1.7*TMath::Power(0.6/x,1)) + 2.5",0.01,5);
    //TF1* f_Dedx_bot = new TF1("f_Dedx","0.5*(TMath::Power(0.5/x,4) - 1*TMath::Power(0.50/x,2)) + 2.7",0.01,5);
    TF1* f_Dedx_top = new TF1("f_Dedx_top","0.55*(TMath::Power(1.62/x,2) - 2.95*TMath::Power(0.6/x,1)) + 3.3",0.01,100); //magenta
    //TF1* f_Dedx_top = new TF1("f_Dedx_top","0.55*(TMath::Power(1.62/x,2) - 2*TMath::Power(0.6/x,1)) + 3.6",0.01,5); //magenta

    //
    bool check     = false;

    std::vector<std::string> AxisLabels =
    {
        "Mass",
        "Pt",
        "Eta",
        "Rapidity",
        "Nhits",
        "DCAz",
        "DCAxy",
        "mult",
    };

    const int nDim = AxisLabels.size();
    std::vector<int> nBins      = {80  ,200,24 ,30 ,35,20 ,20 ,800};
    std::vector<double> minBins = {0.98,0  ,0  ,0  ,0 ,0  ,0  ,0};
    std::vector<double> maxBins = {1.06,20 ,2.4,3.0,35,10.0,10.0,800};

    std::vector<std::vector<float> > VarBinContainers;

    for(int i=0; i<nDim; ++i)
    {
        std::vector<float> tmp;
        VarBinContainers.push_back(tmp);
    }

    THnSparseF* h = new THnSparseF("hsparse", "Cuts", nDim, &nBins[0], &minBins[0], &maxBins[0]);

    for(int i=0; i<VarBinContainers.size(); ++i)
    {
        SparseVarSetUp(VarBinContainers[i], nBins[i], minBins[i], maxBins[i]);
        h->GetAxis(i)->SetName(AxisLabels[i].c_str());
        h->GetAxis(i)->Set(nBins[i],&(VarBinContainers[i])[0]);
    }

    reader.SetEntriesRange(event_start,event_stop);
    while(reader.Next())
    {
        check = true;
        if(event_num == event_stop)
        {
            break;
        }

        if(event_num % check_interval == 0  && event_num != 0 && check)
        {
            cout << (double)event_num/event_stop*100 << " \% completed" << endl;
            check = false;
            gBenchmark->Show("ReadTreeV3");
        }
        if(!CheckValue(momentum)) return false;
        if(!CheckValue(px)) return false;
        if(!CheckValue(py)) return false;
        if(!CheckValue(pz)) return false;
        if(!CheckValue(pt)) return false;
        if(!CheckValue(ptError)) return false;
        if(!CheckValue(energy)) return false;
        if(!CheckValue(dedx)) return false;
        if(!CheckValue(dz)) return false;
        if(!CheckValue(dzError)) return false;
        if(!CheckValue(dxy)) return false;
        if(!CheckValue(dxyError)) return false;
        if(!CheckValue(eta)) return false;
        if(!CheckValue(rapidity)) return false;
        if(!CheckValue(phi)) return false;
        if(!CheckValue(vx)) return false;
        if(!CheckValue(vy)) return false;
        if(!CheckValue(vz)) return false;
        if(!CheckValue(nhits)) return false;
        if(!CheckValue(charge)) return false;
        //if(!CheckValue(mult)) return false;
        //if(!CheckValue(multRaw)) return false;

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

            if(dedx_ <= f_Dedx_bot->Eval(momentum_) || dedx_ >= f_Dedx_top->Eval(momentum_)) continue;
            if(fabs(ptError[i]/pt[i]) > 0.10) continue;

            h_DeDx->Fill(momentum_,dedx_);

            kaon pk(momentum_, pt_, px_, py_, pz_, dedx_, energy_, charge_, nhits_);

            //positive kaons
            if(pk.charge == 1)
                Pkp.push_back(pk);

            //negative kaons
            if(pk.charge == -1)
                Pkm.push_back(pk);

            double dzvtx = dz[i];
            double dxyvtx = dxy[i];
            double dzerror = dzError[i];
            double dxyerror = dxyError[i];

            DCA_z = fabs(dzvtx/dzerror);
            DCA_xy = fabs(dxyvtx/dxyerror);
            mult_++;
        }

        std::vector<double> MassCombinations = CombinatorialMassSparse(Pkp,Pkm);

        for(unsigned i=0; i<MassCombinations.size(); ++i)
        {
            //double value[nDim] = {MassCombinations[i],pt_,eta_,rapidity_,nhits_,DCA_z,DCA_xy,multRaw_,mult_};
            std::vector<double> value = {MassCombinations[i],pt_,eta_,rapidity_,nhits_,DCA_z,DCA_xy,mult_};
            h->Fill(&value[0]);
        }
        event_num++;
        Pkp.clear();
        Pkm.clear();
        h_nEvt->Fill(1);
    }

    TCanvas* c_DeDx = new TCanvas("c_DeDx","c_DeDx",600,600);
    c_DeDx->cd();
    //f_Dedx_bot2->SetLineColor(3);
    h_DeDx->Draw("SCAT");
    f_Dedx_bot->Draw("SAME");
    //f_Dedx_bot2->Draw("SAME");
    f_Dedx_top->Draw("SAME");

    f1->Close();

    TFile* out = new TFile(filename.c_str(),"RECREATE");
    h->Write();
    h_nEvt->Write();
    h_DeDx->Write();

    out->Close();

    //c_Mass->Print("Image/CombMassLooseBand_nhits17_6e7.pdf");
    //c_Mass->Print("Image/CombMassLooseBand_nhits17_6e7.png");
    //c_DeDx->Print("Image/DedxLooseBand_nhits17_6e7.pdf");
    //c_DeDx->Print("Image/DedxLooseBand_nhits17_6e7.png");

    gBenchmark->Show("Total");
    return true;
}

void PhiTreeReader::ReadSparse()
{
    TH1::SetDefaultSumw2();
    TFile *f = new TFile("TreeReader/rootFiles/Sparse_v1.root");
    THnSparseF *h;
    gStyle->SetOptStat(1111111);

    f->GetObject("hsparse",h);

    f->Close();

    TFile *out = new TFile("MassPeaks.root","RECREATE");

    std::vector<TH1D*> hist_containers;
    std::vector<float> VptBins  = {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100};

    h->GetAxis(2)->SetRange(0,24); //eta
    //h->GetAxis(3)->SetRange(0,10); //rapidity
    h->GetAxis(4)->SetRange(12,-1); //nhits
    h->GetAxis(5)->SetRange(0,2); //DCAz
    h->GetAxis(6)->SetRange(0,2); //DCAxy
    //h->GetAxis(7)->SetRange(); //multiplicity
    
    TH1D* h_test = h->Projection(0);

    TCanvas* c = new TCanvas("can","can",600,600);
    //h_test->Draw();

    for(int i=0; i<VptBins.size()-1; ++i)
    {
        h->GetAxis(1)->SetRange(VptBins[i],VptBins[i+1]);
        hist_containers.push_back(h->Projection(0));
    }

    for(int i=0; i<hist_containers.size(); ++i)
    {
        TCanvas* c = new TCanvas(Form("PtBin_%d",i),Form("PtBin_%d",i),600,600);
        hist_containers[i]->SetMarkerStyle(20);
        hist_containers[i]->SetMarkerSize(1.0);
        hist_containers[i]->Sumw2();
        hist_containers[i]->Draw("P E");
        c->Write(Form("PtBin_%d",i));
    }

    out->Close();
}

int main()
{
    PhiTreeReader p(1000000,2000000,10000,"/Volumes/MacHD/Users/blt1/research/Macros/PhiAnalysis/TreeReader/rootFiles/Sparse_v1_1-2M.root");
    //PhiTreeReader p(100000,200000,10000,"/Volumes/MacHD/Users/blt1/research/Macros/PhiAnalysis/TreeReader/rootFiles/Test.root");
    p.ReadTreeV3();

    //TFile* f = new TFile(p.filename.c_str());
}
