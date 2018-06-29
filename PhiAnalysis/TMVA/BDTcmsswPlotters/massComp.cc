void massComp()
{
    TFile* f_data = TFile::Open("/Users/btran/research/PhiReconstruction/TMVA/BDTcmsswPlotters/BDTPROOF/PhiBDT_30M_rapLT1_pt0.5-1.root");
    TFile* f_MC = TFile::Open("/Users/btran/research/Macros/PhiAnalysis/TMVA/BDTcmsswPlotters/BDTPROOF/OutputFiles/PhiBDT_30M_rapLT1_pt0.5-1.root");

    TH1D* h_massData = (TH1D*)f_data->Get("mass_0.25");
    TH1D* h_MCData = (TH1D*)f_MC->Get("mass_0.25");
    h_massData->SetMarkerStyle(20);
    h_massData->SetFillColor(kGreen);
    h_MCData->SetMarkerStyle(20);
    h_MCData->SetFillColor(kGray);
    THStack* hs_mass = new THStack("hs_mass","Stacked Mass Spec");
    hs_mass->Add(h_massData);
    hs_mass->Add(h_MCData);

    TCanvas* c = new TCanvas("c","c",800,800);

    TLegend* leg = new TLegend(0.88,0.88,0.75,0.75);
    leg->AddEntry(h_massData,"Data","f");
    leg->AddEntry(h_MCData,"MC","f");

    hs_mass->Draw("nostack");
    leg->Draw();
    hs_mass->GetHistogram()->GetYaxis()->SetRangeUser(50000,130000);

    h_massData->Add(h_MCData,-1);
    TCanvas* c_diff = new TCanvas("c_diff","c_diff",800,800);
    h_massData->Draw();
}
