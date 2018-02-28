#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

  //Example of script showing how to copy all objects (including directories)
  //from a source file.
  //For each input file, a new directory is created in the current directory
  //with the name of the source file.
  //After execution of:
  // root > .x copyFiles.C
  //the file result.root will contain 4 subdirectories:
  // "tot100.root", "hsimple.root", "hs1.root","hs2.root"

void CopyDir(TDirectory *source) {
   //copy all objects and subdirs of directory source as a subdir of the current directory
   source->ls();
   TDirectory *savdir = gDirectory;
   TDirectory *adir;
   //TDirectory *adir = savdir->mkdir(source->GetName());
   if(savdir->GetDirectory("v0CasCorrelationRapidity")) adir = savdir->GetDirectory("v0CasCorrelationRapidity");
   else adir = savdir->mkdir("v0CasCorrelationRapidity");
   adir->cd();
   //loop on all entries of this directory
   TKey *key;
   TIter nextkey(source->GetListOfKeys());
   while ((key = (TKey*)nextkey())) {
      const char *classname = key->GetClassName();
      TClass *cl = gROOT->GetClass(classname);
      if (!cl) continue;
      if (cl->InheritsFrom(TDirectory::Class())) {
         source->cd(key->GetName());
         TDirectory *subdir = gDirectory;
         adir->cd();
         CopyDir(subdir);
         adir->cd();
      } else if (cl->InheritsFrom(TTree::Class())) {
         TTree *T = (TTree*)source->Get(key->GetName());
         adir->cd();
         TTree *newT = T->CloneTree(-1,"fast");
         newT->Write();
      } else {
         source->cd();
         TObject *obj = key->ReadObj();
         adir->cd();
         obj->Write();
         delete obj;
     }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}
void CopyFile(const char *fname) {
   //Copy all objects and subdirs of file fname as a subdir of the current directory
   TDirectory *target = gDirectory;
   TFile *f = TFile::Open(fname);
   if (!f || f->IsZombie()) {
      printf("Cannot copy file: %s\n",fname);
      target->cd();
      return;
   }
   target->cd();
   CopyDir(f);
   delete f;
   target->cd();
}
void copyFiles() {
   //prepare files to be copied
   if(gSystem->AccessPathName("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityCorrectMultB_09_19_17.root")) {
      gSystem->CopyFile("hsimple.root", "/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityCorrectMultB_09_19_17.root");
      gSystem->CopyFile("hsimple.root", "/Bolumes/MacHD/Users/blt1/research/RootFiles/Flow/OmCorr/OmCorrelationRapidityTotal_09_24_17.root");
      gSystem->CopyFile("hsimple.root", "hs2.root");
   }
   //main function copying 4 files as subdirectories of a new file
   TFile *f = new TFile("result.root","recreate");
   CopyFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityCorrectMultB_09_19_17.root");
   CopyFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/OmCorr/OmCorrelationRapidityTotal_09_24_17.root");
   CopyFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationRapidityTotal_08_20_2017.root");
   f->ls();
   delete f;
}

void Rename()
{
    TFile *f_V0 = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/V0Corr/V0CorrelationRapidityCorrectMultB_09_19_17.root");
    //TFile *f_Xi = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/XiCorr/XiCorrelationRapidityTotal_08_20_2017.root" );
    TFile *f_Xi = new TFile("XiCorrelation.root" );
    //TFile *f_Om = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/OmCorr/OmCorrelationRapidityTotal_09_24_17.root" );
    TFile *f_Om = new TFile("OmCorrelation.root" );

    //TFile* f = new TFile("XiCorrelation.root","RECREATE");
    //TFile* f = new TFile("OmCorrelation.root","RECREATE");
    TFile* f = new TFile("V0CasCorrelationRapidityMerged.root","RECREATE");
    CopyDir(f_V0->GetDirectory("v0CorrelationRapidity"));
    CopyDir(f_Xi->GetDirectory("v0CasCorrelationRapidity"));
    //CopyDir(f_Xi->GetDirectory("xiCorrelationRapidity"));
    CopyDir(f_Om->GetDirectory("v0CasCorrelationRapidity"));
    //CopyDir(f_Om->GetDirectory("omCorrelationRapidity"));
    //
    f->Close();

}

void HistRename()
{
    TFile* f = new TFile("OmegaMB_Total_0_35_11_28_17.root","UPDATE");
    std::vector<double> PtBin_om = {1.5, 2.2, 2.8, 3.6, 5.0, 8.0};//, 20.0}; //pPb
    int numPtBins_om = PtBin_om.size() - 1;

    std::string fn = "v0CasCorrelationRapidityPeriSub";

    TDirectory* currentDir = gDirectory;
    //TDirectory *adir = currentDir->mkdir("v0CasCorrelationRapidity");
    //TDirectory *adir = f->GetDirectory("v0CasCorrelationRapidity");
    TDirectory *adir = f->GetDirectory(fn.c_str());
    adir->cd();

    for(int i=1; i<9; i++){
        cout << i << endl;
            TH2D *hbackgroundPeak = (TH2D*) f->Get( Form((fn + "/BackgroundOmPeak_pt%d").c_str(),i ) );
            TH2D *hsignalPeak     = (TH2D*) f->Get( Form((fn + "/SignalOmPeak_pt%d").c_str(),i ) );
            TH2D *hbackgroundSide = (TH2D*) f->Get( Form((fn + "/BackgroundOmSide_pt%d").c_str(),i ) );
            TH2D *hsignalSide     = (TH2D*) f->Get( Form((fn + "/SignalOmSide_pt%d").c_str(),i ) );
            TH1D *hKET = (TH1D*) f     -> Get(Form((fn+  "/KET_om_pt%d" ).c_str(),i));
            TH1D *hKET_bkg = (TH1D*) f -> Get(Form((fn+  "/KET_om_bkg_pt%d" ).c_str(),i));
            TH1D *hMass = (TH1D*) f    -> Get(Form((fn+  "/Mass_om_pt%d" ).c_str(),i));
            TH1D *hPt = (TH1D*) f      -> Get(Form((fn+  "/Pt_om_pt%d" ).c_str(),i));
            TH1D *hPt_bkg = (TH1D*) f  -> Get(Form((fn+  "/Pt_om_bkg_pt%d" ).c_str(),i));
            TH1D *hEta = (TH1D*) f     -> Get(Form((fn+  "/Eta_om_pt%d" ).c_str(),i));
            TH1D *hEta_bkg = (TH1D*) f -> Get(Form((fn+  "/Eta_om_bkg_pt%d" ).c_str(),i));
            TH1D* hrap = (TH1D*) f     -> Get(Form((fn+  "/rap_om_pt%d" ).c_str(),i));
            TH1D* hrap_bkg = (TH1D*) f -> Get(Form((fn+  "/rap_om_bkg_pt%d" ).c_str(),i));
            TH1D* hmult = (TH1D*) f -> Get(Form((fn+     "/mult_om_pt%d" ).c_str(),i));
            TH1D* hmult_bkg = (TH1D*) f -> Get(Form((fn+ "/mult_om_bkg_pt%d" ).c_str(),i));

            hbackgroundPeak->SetName(Form("BackgroundOmPeak_pt%d",i-1 ));
            hsignalPeak->SetName(Form("SignalOmPeak_pt%d",i-1 ));
            hbackgroundSide->SetName(Form("BackgroundOmSide_pt%d",i-1 ));
            hsignalSide->SetName(Form("SignalOmSide_pt%d",i-1 ));
            hKET->SetName      (Form("KET_om_pt%d",i-1) );
            hKET_bkg->SetName  (Form( "KET_om_bkg_pt%d" ,i-1) );
            hMass->SetName     (Form( "Mass_om_pt%d" ,i-1) );
            hPt->SetName       (Form( "Pt_om_pt%d" ,i-1) );
            hPt_bkg->SetName   (Form( "Pt_om_bkg_pt%d" ,i-1) );
            hEta->SetName      (Form( "Eta_om_pt%d" ,i-1) );
            hEta_bkg->SetName  (Form( "Eta_om_bkg_pt%d" ,i-1) );
            hrap->SetName      (Form( "rap_om_pt%d" ,i-1) );
            hrap_bkg->SetName  (Form( "rap_om_bkg_pt%d" ,i-1) );
            hmult->SetName     (Form( "mult_om_pt%d" ,i-1) );
            hmult_bkg->SetName (Form( "mult_om_bkg_pt%d" ,i-1) );


            //hbackgroundPeak->Write(Form("BackgroundXiPeak_pt%d",i ));
            //hsignalPeak->Write(Form("SignalXiPeak_pt%d",i ));
            //hbackgroundSide->Write(Form("BackgroundXiSide_pt%d",i ));
            //hsignalSide->Write(Form("SignalXiSide_pt%d",i ));
            hbackgroundPeak->Write(Form("BackgroundOmPeak_pt%d",i-1 ));
            hsignalPeak->Write(Form("SignalOmPeak_pt%d",i-1 ));
            hbackgroundSide->Write(Form("BackgroundOmSide_pt%d",i-1 ));
            hsignalSide->Write(Form("SignalOmSide_pt%d",i-1 ));
            hKET->Write      (Form("KET_om_pt%d",i-1) );
            hKET_bkg->Write  (Form( "KET_om_bkg_pt%d" ,i-1) );
            hMass->Write     (Form( "Mass_om_pt%d" ,i-1) );
            hPt->Write       (Form( "Pt_om_pt%d" ,i-1) );
            hPt_bkg->Write   (Form( "Pt_om_bkg_pt%d" ,i-1) );
            hEta->Write      (Form( "Eta_om_pt%d" ,i-1) );
            hEta_bkg->Write  (Form( "Eta_om_bkg_pt%d" ,i-1) );
            hrap->Write      (Form( "rap_om_pt%d" ,i-1) );
            hrap_bkg->Write  (Form( "rap_om_bkg_pt%d" ,i-1) );
            hmult->Write     (Form( "mult_om_pt%d" ,i-1) );
            hmult_bkg->Write (Form( "mult_om_bkg_pt%d" ,i-1) );

            gDirectory->Delete(Form("SignalOmPeak_pt%d;1",i));
            gDirectory->Delete(Form("BackgroundOmPeak_pt%d;1",i));
            gDirectory->Delete(Form("SignalOmSide_pt%d;1",i));
            gDirectory->Delete(Form("BackgroundOmSide_pt%d;1",i));
            gDirectory->Delete(Form("rap_om_Lorentz_pt%d;1",i));
            gDirectory->Delete(Form("KET_om_pt%d;1" ,i));
            gDirectory->Delete(Form("KET_om_bkg_pt%d;1" ,i));
            gDirectory->Delete(Form("Mass_om_pt%d;1" ,i));
            gDirectory->Delete(Form("Pt_om_pt%d;1" ,i));
            gDirectory->Delete(Form("Pt_om_bkg_pt%d;1" ,i));
            gDirectory->Delete(Form("Eta_om_pt%d;1" ,i));
            gDirectory->Delete(Form("Eta_om_bkg_pt%d;1" ,i));
            gDirectory->Delete(Form("rap_om_pt%d;1" ,i));
            gDirectory->Delete(Form("rap_om_bkg_pt%d;1" ,i));
            gDirectory->Delete(Form("mult_om_pt%d;1" ,i));
            gDirectory->Delete(Form("mult_om_bkg_pt%d;1" ,i));
    }

}

void MergeBins()
{
    TFile* f = new TFile("OmegaMB_Total_0_35_11_28_17.root","UPDATE");

    std::string fn = "v0CasCorrelationRapidityPeriSub";

    TDirectory* currentDir = gDirectory;
    TDirectory *adir = f->GetDirectory(fn.c_str());
    adir->cd();

    TH2D* hBackgroundPeak0 = (TH2D*)f->Get((fn + "/BackgroundOmPeak_pt0").c_str());
    TH2D* hBackgroundPeak1 = (TH2D*)f->Get((fn + "/BackgroundOmPeak_pt1").c_str());
    TH2D* hBackgroundSide0 = (TH2D*)f->Get((fn + "/BackgroundOmSide_pt0").c_str());
    TH2D* hBackgroundSide1 = (TH2D*)f->Get((fn + "/BackgroundOmSide_pt1").c_str());

    TH2D* hSignalPeak0 = (TH2D*)f->Get((fn + "/SignalOmPeak_pt0").c_str());
    TH2D* hSignalPeak1 = (TH2D*)f->Get((fn + "/SignalOmPeak_pt1").c_str());
    TH2D* hSignalSide0 = (TH2D*)f->Get((fn + "/SignalOmSide_pt0").c_str());
    TH2D* hSignalSide1 = (TH2D*)f->Get((fn + "/SignalOmSide_pt1").c_str());

    TH1D* hMult0 = (TH1D*)f->Get((fn + "/mult_om_pt0").c_str());
    TH1D* hMult1 = (TH1D*)f->Get((fn + "/mult_om_pt1").c_str());
    TH1D* hMult_bkg0 = (TH1D*)f->Get((fn + "/mult_om_bkg_pt0").c_str());
    TH1D* hMult_bkg1 = (TH1D*)f->Get((fn + "/mult_om_bkg_pt1").c_str());

    TH1D* hPt0 = (TH1D*)f->Get((fn + "/Pt_om_pt0").c_str());
    TH1D* hPt1 = (TH1D*)f->Get((fn + "/Pt_om_pt1").c_str());
    TH1D* hPt_bkg0 = (TH1D*)f->Get((fn + "/Pt_om_bkg_pt0").c_str());
    TH1D* hPt_bkg1 = (TH1D*)f->Get((fn + "/Pt_om_bkg_pt1").c_str());

    TH1D* hKET0 = (TH1D*)f->Get((fn + "/KET_om_pt0").c_str());
    TH1D* hKET1 = (TH1D*)f->Get((fn + "/KET_om_pt1").c_str());
    TH1D* hKET_bkg0 = (TH1D*)f->Get((fn + "/KET_om_bkg_pt0").c_str());
    TH1D* hKET_bkg1 = (TH1D*)f->Get((fn + "/KET_om_bkg_pt1").c_str());

    TH1D* Mass0 = (TH1D*)f->Get((fn + "/Mass_om_pt0").c_str());
    TH1D* Mass1 = (TH1D*)f->Get((fn + "/Mass_om_pt1").c_str());

    TH1D* Eta0 = (TH1D*)f->Get((fn + "/Eta_om_pt0").c_str());
    TH1D* Eta1 = (TH1D*)f->Get((fn + "/Eta_om_pt1").c_str());
    TH1D* Eta_bkg0 = (TH1D*)f->Get((fn + "/Eta_om_bkg_pt0").c_str());
    TH1D* Eta_bkg1 = (TH1D*)f->Get((fn + "/Eta_om_bkg_pt1").c_str());

    TH1D* Rap0 = (TH1D*)f->Get((fn + "/rap_om_pt0").c_str());
    TH1D* Rap1 = (TH1D*)f->Get((fn + "/rap_om_pt1").c_str());
    TH1D* Rap_bkg0 = (TH1D*)f->Get((fn + "/rap_om_bkg_pt0").c_str());
    TH1D* Rap_bkg1 = (TH1D*)f->Get((fn + "/rap_om_bkg_pt1").c_str());

    hMult0->GetXaxis()->SetRange(2,250);
    hMult1->GetXaxis()->SetRange(2,250);
    hMult_bkg0->GetXaxis()->SetRange(2,250);
    hMult_bkg1->GetXaxis()->SetRange(2,250);

    double Mult_f0 = hMult0->GetMean(1);
    double Mult_f1 = hMult1->GetMean(1);
    double Mult_bkg_f0 = hMult_bkg0->GetMean(1);
    double Mult_bkg_f1 = hMult_bkg1->GetMean(1);

    hBackgroundPeak0->Scale(Mult_f0);
    hBackgroundPeak1->Scale(Mult_f1);
    hBackgroundSide0->Scale(Mult_bkg_f0);
    hBackgroundSide1->Scale(Mult_bkg_f1);

    hSignalPeak0->Scale(Mult_f0);
    hSignalPeak1->Scale(Mult_f1);
    hSignalSide0->Scale(Mult_bkg_f0);
    hSignalSide1->Scale(Mult_bkg_f1);

    Mass0->Add(Mass1);
    Eta0->Add(Eta1);
    Eta_bkg0->Add(Eta_bkg1);
    Rap0->Add(Rap1);
    Rap_bkg0->Add(Rap_bkg1);
    hBackgroundPeak0 -> Add(hBackgroundPeak1);
    hBackgroundSide0 -> Add(hBackgroundSide1);
    hSignalPeak0     -> Add(hSignalPeak1);
    hSignalSide0     -> Add(hSignalSide1);
    hMult0           -> Add(hMult1);
    hMult_bkg0       -> Add(hMult_bkg1);
    hPt0             -> Add(hPt1);
    hPt_bkg0         -> Add(hPt_bkg1);
    hKET0            -> Add(hKET1);
    hKET_bkg0        -> Add(hKET_bkg1);

    hMult0->GetXaxis()->SetRange(2,250);
    hMult_bkg0->GetXaxis()->SetRange(2,250);

    double Mult_f_combined = hMult0->GetMean(1);
    double Mult_f_bkg_combined = hMult_bkg0->GetMean(1);

    hBackgroundPeak0->Scale(1.0/Mult_f_combined);
    hBackgroundSide0->Scale(1.0/Mult_f_bkg_combined);
    hSignalPeak0->Scale(1.0/Mult_f_combined);
    hSignalSide0->Scale(1.0/Mult_f_bkg_combined);

    Mass0->Write("Mass_om_pt1",TObject::kOverwrite);
    Eta0->Write("Eta_om_pt1",TObject::kOverwrite);
    Eta_bkg0->Write("Eta_om_bkg_pt1",TObject::kOverwrite);
    Rap0->Write("rap_om_pt1",TObject::kOverwrite);
    Rap_bkg0->Write("rap_om_bkg_pt1",TObject::kOverwrite);
    hBackgroundPeak0->Write("BackgroundOmPeak_pt1",TObject::kOverwrite);
    hBackgroundSide0->Write("BackgroundOmSide_pt1",TObject::kOverwrite);
    hSignalPeak0    ->Write("SignalOmPeak_pt1",TObject::kOverwrite);
    hSignalSide0    ->Write("SignalOmSide_pt1",TObject::kOverwrite);
    hMult0          ->Write("mult_om_pt1",TObject::kOverwrite);
    hMult_bkg0      ->Write("mult_om_bkg_pt1",TObject::kOverwrite);
    hPt0            ->Write("Pt_om_pt1",TObject::kOverwrite);
    hPt_bkg0        ->Write("Pt_om_bkg_pt1",TObject::kOverwrite);
    hKET0           ->Write("KET_om_pt1",TObject::kOverwrite);
    hKET_bkg0       ->Write("KET_om_bkg_pt1",TObject::kOverwrite);

    gDirectory->Delete("BackgroundOmPeak_pt0;1");
    gDirectory->Delete("BackgroundOmSide_pt0;1");
    gDirectory->Delete("SignalOmPeak_pt0;1");
    gDirectory->Delete("SignalOmSide_pt0;1");
    gDirectory->Delete("mult_om_pt0;1");
    gDirectory->Delete("mult_om_bkg_pt0;1");
    gDirectory->Delete("Pt_om_pt0;1");
    gDirectory->Delete("Pt_om_bkg_pt0;1");
    gDirectory->Delete("KET_om_pt0;1");
    gDirectory->Delete("KET_om_bkg_pt0;1");
    gDirectory->Delete("Eta_om_pt0;1");
    gDirectory->Delete("Eta_om_bkg_pt0;1");
    gDirectory->Delete("rap_om_pt0;1");
    gDirectory->Delete("rap_om_bkg_pt0;1");
    gDirectory->Delete("Mass_om_pt0;1");


}
