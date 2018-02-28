void perisubrefalltestV0true(){
TH1::SetDefaultSumw2();

    TFile *_file0 = TFile::Open("V0v2_vspt_sub1020_NassFit_highpt_pt033_obs_3term_7TeVEff.root");
    TFile *_file1 = TFile::Open("V0v2_vspt_sub1020_NassFit_highpt_pt033_bkg_3terms_7TeVEff.root");
    
    double sigfrks[12] = {0.963688,0.957693,0.956367,0.958637,0.960676,0.962354,0.96037,0.955859,0.948897,0.935449,0.921266};
    double sigfrla[12] = {0.239066,0.470361,0.721976,0.802508,0.859667,0.896981,0.921889,0.935101,0.930692,0.895919,0.825003};
    
    //read vn values
    TGraphErrors* v2obs_ks = _file0->Get("kshortv2_obs_GplusPP");
    TGraphErrors* v2obs_ks_sub = _file0->Get("kshortv2sub_obs_GplusPP");
    TGraphErrors* v2obs_la = _file0->Get("lambdav2_obs_GplusPP");
    TGraphErrors* v2obs_la_sub = _file0->Get("lambdav2sub_obs_GplusPP");
    
    TGraphErrors* v2bkg_ks = _file1->Get("kshortv2_bkg_GplusPP");
    TGraphErrors* v2bkg_ks_sub = _file1->Get("kshortv2sub_bkg_GplusPP");
    TGraphErrors* v2bkg_la = _file1->Get("lambdav2_bkg_GplusPP");
    TGraphErrors* v2bkg_la_sub = _file1->Get("lambdav2sub_bkg_GplusPP");
    
    TGraphErrors* v2obs_ks_KET = _file0->Get("kshortv2_obs_GplusPP_KET");
    TGraphErrors* v2obs_ks_sub_KET = _file0->Get("kshortv2sub_obs_GplusPP_KET");
    TGraphErrors* v2obs_la_KET = _file0->Get("lambdav2_obs_GplusPP_KET");
    TGraphErrors* v2obs_la_sub_KET = _file0->Get("lambdav2sub_obs_GplusPP_KET");
    
    TGraphErrors* v2bkg_ks_KET = _file1->Get("kshortv2_bkg_GplusPP_KET");
    TGraphErrors* v2bkg_ks_sub_KET = _file1->Get("kshortv2sub_bkg_GplusPP_KET");
    TGraphErrors* v2bkg_la_KET = _file1->Get("lambdav2_bkg_GplusPP_KET");
    TGraphErrors* v2bkg_la_sub_KET = _file1->Get("lambdav2sub_bkg_GplusPP_KET");
    
    double* v2ks = v2obs_ks->GetY();
    double* v2eks = v2obs_ks->GetEY();
    double* v2ks_bkg = v2bkg_ks->GetY();
    double* v2eks_bkg = v2bkg_ks->GetEY();
    
    double* v2kssub = v2obs_ks_sub->GetY();
    double* v2ekssub = v2obs_ks_sub->GetEY();
    double* v2ks_bkgsub = v2bkg_ks_sub->GetY();
    double* v2eks_bkgsub = v2bkg_ks_sub->GetEY();
    
    double* v2la = v2obs_la->GetY();
    double* v2ela = v2obs_la->GetEY();
    double* v2la_bkg = v2bkg_la->GetY();
    double* v2ela_bkg = v2bkg_la->GetEY();
    
    double* v2lasub = v2obs_la_sub->GetY();
    double* v2elasub = v2obs_la_sub->GetEY();
    double* v2la_bkgsub = v2bkg_la_sub->GetY();
    double* v2ela_bkgsub = v2bkg_la_sub->GetEY();
    
    double* ptks = v2obs_ks->GetX();
    double* ptla = v2obs_la->GetX();

    double* KETks = v2obs_ks_KET->GetX();
    double* KETla = v2obs_la_KET->GetX();

    //bkg subtraction
    double bkgfrks[13];
    double v2tks[13];
    double v2teks[13];
    double v2tkssub[13];
    double v2tekssub[13];

    double bkgfrla[13];
    double v2tla[13];
    double v2tela[13];
    double v2tlasub[13];
    double v2telasub[13];
    
    for(int i=0;i<9;i++)
    {
        bkgfrks[i] = 1 - sigfrks[i];
        v2tks[i] = (v2ks[i] - v2ks_bkg[i]*bkgfrks[i])/sigfrks[i];
        v2teks[i] = sqrt((v2eks_bkg[i]*bkgfrks[i])**2 + v2eks[i]**2)/sigfrks[i];
        v2tkssub[i] = (v2kssub[i] - v2ks_bkgsub[i]*bkgfrks[i])/sigfrks[i];
        v2tekssub[i] = sqrt((v2eks_bkgsub[i]*bkgfrks[i])**2 + v2ekssub[i]**2)/sigfrks[i];
        
        bkgfrla[i] = 1 - sigfrla[i];
        v2tla[i] = (v2la[i] - v2la_bkg[i]*bkgfrla[i])/sigfrla[i];
        v2tela[i] = sqrt((v2ela_bkg[i]*bkgfrla[i])**2 + v2ela[i]**2)/sigfrla[i];
        v2tlasub[i] = (v2lasub[i] - v2la_bkgsub[i]*bkgfrla[i])/sigfrla[i];
        v2telasub[i] = sqrt((v2ela_bkgsub[i]*bkgfrla[i])**2 + v2elasub[i]**2)/sigfrla[i];
        if(i<2) v2tla[i] = -999;
        if(i<2) v2tlasub[i] = -999;
    }
    
    TGraphErrors *ksv2tg = new TGraphErrors(9,ptks,v2tks,0,v2teks);
    TGraphErrors *lav2tg = new TGraphErrors(9,ptla,v2tla,0,v2tela);
    TGraphErrors *ksv2tgsub = new TGraphErrors(9,ptks,v2tkssub,0,v2tekssub);
    TGraphErrors *lav2tgsub = new TGraphErrors(9,ptla,v2tlasub,0,v2telasub);
    
    TGraphErrors *ksv2tg_KET = new TGraphErrors(9,KETks,v2tks,0,v2teks);
    TGraphErrors *lav2tg_KET = new TGraphErrors(9,KETla,v2tla,0,v2tela);
    TGraphErrors *ksv2tgsub_KET = new TGraphErrors(9,KETks,v2tkssub,0,v2tekssub);
    TGraphErrors *lav2tgsub_KET = new TGraphErrors(9,KETla,v2tlasub,0,v2telasub);
    
    ksv2tg_KET->SetName("kshortv2true_KET");
    lav2tg_KET->SetName("lambdav2true_KET");
    ksv2tgsub_KET->SetName("kshortv2truesub_KET");
    lav2tgsub_KET->SetName("lambdav2truesub_KET");

    ksv2tg->SetName("kshortv2true");
    lav2tg->SetName("lambdav2true");
    ksv2tgsub->SetName("kshortv2truesub");
    lav2tgsub->SetName("lambdav2truesub");

    TFile ofile("V0v2_vspt_sub1020_NassFit_highpt_pt033_true_3term_7TeVEff.root","RECREATE");
    
    ksv2tg.Write();
    lav2tg.Write();
    ksv2tgsub.Write();
    lav2tgsub.Write();
    
    ksv2tg_KET.Write();
    lav2tg_KET.Write();
    ksv2tgsub_KET.Write();
    lav2tgsub_KET.Write();
}
