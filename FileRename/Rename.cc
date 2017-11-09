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
    TFile *f_Xi = new TFile("xiCorrelationRapidity.root" );
    TFile *f_Om = new TFile("/volumes/MacHD/Users/blt1/research/RootFiles/Flow/OmCorr/OmCorrelationRapidityTotal_09_24_17.root" );

    //TFile* f = new TFile("XiCorrelation.root","UPDATE");
    TFile* f = new TFile("OmCorrelation.root","UPDATE");
    std::vector<double> PtBin_xi = {1.0, 1.4, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10.0};//, 20.0};
    std::vector<double> PtBin_om = {1.0, 1.5, 1.8, 2.2, 2.8, 3.6, 4.6, 6.0, 7.2, 10.0};//, 20.0}; //pPb
    int numPtBins_xi = PtBin_xi.size() - 1;
    int numPtBins_om = PtBin_om.size() - 1;

    std::string fn = "v0CasCorrelationRapidity";

    TDirectory* currentDir = gDirectory;
    //TDirectory *adir = currentDir->mkdir("v0CasCorrelationRapidity");
    //TDirectory *adir = f->GetDirectory("v0CasCorrelationRapidity");
    TDirectory *adir = f->GetDirectory(fn.c_str());
    adir->cd();

    TH2D *hMassPt = (TH2D*) f->Get((fn+"/MassPt").c_str());

    hMassPt->Write("XiMassPt");
    gDirectory->Delete("MassPt;1");
    gDirectory->Delete("BackgroundPeak_pt0;1");
    gDirectory->Delete("BackgroundSide_pt0;1");
    gDirectory->Delete("SignalPeak_pt0;1");
    gDirectory->Delete("SignalSide_pt0;1");
    gDirectory->Delete(Form("rap_xi_Lorentz_pt%d;1",0));
    gDirectory->Delete(Form("KET_xi_pt%d;1" ,0));
    gDirectory->Delete(Form("KET_xi_bkg_pt%d;1" ,0));
    gDirectory->Delete(Form("Mass_xi_pt%d;1" ,0));
    gDirectory->Delete(Form("Pt_xi_pt%d;1" ,0));
    gDirectory->Delete(Form("Pt_xi_bkg_pt%d;1" ,0));
    gDirectory->Delete(Form("Eta_xi_pt%d;1" ,0));
    gDirectory->Delete(Form("Eta_xi_bkg_pt%d;1" ,0));
    gDirectory->Delete(Form("rap_xi_pt%d;1" ,0));
    gDirectory->Delete(Form("rap_xi_bkg_pt%d;1" ,0));
    gDirectory->Delete(Form("mult_xi_pt%d;1" ,0));
    gDirectory->Delete(Form("mult_xi_bkg_pt%d;1" ,0));
    for(int i=1; i<9; i++){
        cout << i << endl;
            TH2D *hbackgroundPeak = (TH2D*) f->Get( Form((fn + "/BackgroundPeak_pt%d").c_str(),i ) );
            TH2D *hsignalPeak     = (TH2D*) f->Get( Form((fn + "/SignalPeak_pt%d").c_str(),i ) );
            TH2D *hbackgroundSide = (TH2D*) f->Get( Form((fn + "/BackgroundSide_pt%d").c_str(),i ) );
            TH2D *hsignalSide     = (TH2D*) f->Get( Form((fn + "/SignalSide_pt%d").c_str(),i ) );
            TH1D *hKET = (TH1D*) f     -> Get(Form((fn+  "/KET_xi_pt%d" ).c_str(),i));
            TH1D *hKET_bkg = (TH1D*) f -> Get(Form((fn+  "/KET_xi_bkg_pt%d" ).c_str(),i));
            TH1D *hMass = (TH1D*) f    -> Get(Form((fn+  "/Mass_xi_pt%d" ).c_str(),i));
            TH1D *hPt = (TH1D*) f      -> Get(Form((fn+  "/Pt_xi_pt%d" ).c_str(),i));
            TH1D *hPt_bkg = (TH1D*) f  -> Get(Form((fn+  "/Pt_xi_bkg_pt%d" ).c_str(),i));
            TH1D *hEta = (TH1D*) f     -> Get(Form((fn+  "/Eta_xi_pt%d" ).c_str(),i));
            TH1D *hEta_bkg = (TH1D*) f -> Get(Form((fn+  "/Eta_xi_bkg_pt%d" ).c_str(),i));
            TH1D* hrap = (TH1D*) f     -> Get(Form((fn+  "/rap_xi_pt%d" ).c_str(),i));
            TH1D* hrap_bkg = (TH1D*) f -> Get(Form((fn+  "/rap_xi_bkg_pt%d" ).c_str(),i));
            TH1D* hmult = (TH1D*) f -> Get(Form((fn+     "/mult_xi_pt%d" ).c_str(),i));
            TH1D* hmult_bkg = (TH1D*) f -> Get(Form((fn+ "/mult_xi_bkg_pt%d" ).c_str(),i));

            hbackgroundPeak->SetName(Form("BackgroundXiPeak_pt%d",i-1 ));
            hsignalPeak->SetName(Form("SignalXiPeak_pt%d",i-1 ));
            hbackgroundSide->SetName(Form("BackgroundXiSide_pt%d",i-1 ));
            hsignalSide->SetName(Form("SignalXiSide_pt%d",i-1 ));
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

            gDirectory->Delete(Form("SignalPeak_pt%d;1",i));
            gDirectory->Delete(Form("BackgroundPeak_pt%d;1",i));
            gDirectory->Delete(Form("SignalSide_pt%d;1",i));
            gDirectory->Delete(Form("BackgroundSide_pt%d;1",i));
            gDirectory->Delete(Form("rap_xi_Lorentz_pt%d;1",i));
            gDirectory->Delete(Form("KET_xi_pt%d;1" ,i));
            gDirectory->Delete(Form("KET_xi_bkg_pt%d;1" ,i));
            gDirectory->Delete(Form("Mass_xi_pt%d;1" ,i));
            gDirectory->Delete(Form("Pt_xi_pt%d;1" ,i));
            gDirectory->Delete(Form("Pt_xi_bkg_pt%d;1" ,i));
            gDirectory->Delete(Form("Eta_xi_pt%d;1" ,i));
            gDirectory->Delete(Form("Eta_xi_bkg_pt%d;1" ,i));
            gDirectory->Delete(Form("rap_xi_pt%d;1" ,i));
            gDirectory->Delete(Form("rap_xi_bkg_pt%d;1" ,i));
            gDirectory->Delete(Form("mult_xi_pt%d;1" ,i));
            gDirectory->Delete(Form("mult_xi_bkg_pt%d;1" ,i));
    }

}
