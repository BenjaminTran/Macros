#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdio>
#include <TROOT.h>
#include <TStyle.h>

#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TAxis.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "TString.h"

#include <vector>

void V0MassFit()
{
    gStyle->SetMarkerSize(0.8);

    RooMsgService::instance().setStreamStatus(0,kFALSE);
    RooMsgService::instance().setStreamStatus(1,kFALSE);

    TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/CascadeV2pPb/RootFiles/Flow/Ksla/kslaMassPtJL1.root");
    //TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/CascadeV2pPb/RootFiles/Flow/Ksla/kslaMassPtJL2.root");

    TH1D* massks;
    TH1D* massla;
	TH2D* MassKs;
	TH2D* MassLa;
	bool lambda;

    MassKs = (TH2D*)file->Get("KslaMassPt/KsMassPt");
    MassLa = (TH2D*)file->Get("KslaMassPt/LaMassPt");
	int pTksLength = 26; // the number of bins to be fitted is half of this number
	double pks[] = {3,4, 5,6, 7,8, 9,10, 11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,90, 91,120};
	double pla[] = {0,0, 0,0, 0,0, 9,10, 11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60, 61,90, 91,120};

	std::vector<double> mass_ks;
	std::vector<double> std_ks;
	std::vector<double> fsig_ks;
	std::vector<double> covQual_ks;

	std::vector<double> mass_la;
	std::vector<double> std_la;
	std::vector<double> fsig_la;
	std::vector<double> covQual_la;



    std::ostringstream os;

	int pkscounter = 0; //for correct bin counting
	int placounter = 0;
	int hbincounter =1;
	for( int i=0; i<pTksLength; i++ )
    {
        if( pla[i] ==0 ) lambda = false;
        else lambda = true;
        massks = ( TH1D* )MassKs->ProjectionX( "massks", pks[i],pks[i+1] );
        massla = ( TH1D* )MassLa->ProjectionX( "massla", pla[i],pla[i+1] );

        using namespace RooFit;
        gStyle->SetOptTitle(kFALSE);

        TCanvas* cc1 = new TCanvas("cc1","cc1",1200,450);
        cc1->Divide(2,1);

        TLatex* tex = new TLatex();
        tex->SetNDC();
        tex->SetTextFont(42);

        //char label_energy[200]={"CMS pp #sqrt{s} = 13 TeV"};
        //char label_n[200]={"105 #leq N^{offline}_{trk} < 150"};
        //char label_pid[2][200]={"K^{0}_{S}","#Lambda/#bar{#Lambda}"};
        //char label_mean[2][200]={"Mean: 0.4976 GeV","Mean: 1.1159 GeV"};
        //char label_sigma[2][200]={"Average #sigma: 0.0067 GeV","Average #sigma: 0.0031 GeV"};
        //char label_cms[2][200]={"Preliminary","L = 0.7 pb^{-1}"};

        tex->SetTextSize(tex->GetTextSize()*0.95);

        //kshort
        RooRealVar x("x","mass",0.43,0.565);
        RooDataHist data("data","dataset",x,massks);
        RooPlot* xframe = x.frame(270);
        xframe->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");
        xframe->GetYaxis()->SetTitle("Candidates / 0.0005 GeV");
        xframe->GetXaxis()->CenterTitle(1);
        xframe->GetYaxis()->CenterTitle(1);
        xframe->GetXaxis()->SetTitleSize(xframe->GetXaxis()->GetTitleSize()*1.4);
        xframe->GetXaxis()->SetTitleOffset(1);
        xframe->GetYaxis()->SetTitleSize(xframe->GetYaxis()->GetTitleSize()*1.3);
        //xframe->GetYaxis()->SetTitleOffset(1);
        data.plotOn(xframe,Name("data"));
        RooRealVar mean("mean","mean",0.50,0.49,0.51);
        RooRealVar sigma1("sigma1","sigma1",0.01,0.001,0.04);
        RooRealVar sigma2("sigma2","sigma2",0.01,0.001,0.04);
        RooRealVar sig1("sig1","signal1",10,0,1000000000);
        RooRealVar sig2("sig2","signal2",10,0,1000000000);
        RooRealVar a("a","a",0,-100000,100000);
        RooRealVar b("b","b",0,-100000,100000);
        RooRealVar cp("cp","cp",0,-100000,100000);
        RooRealVar d("d","d",0,-100000,100000);
        RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
        RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
        RooPolynomial poly("poly","poly",x,RooArgList(a,b,cp,d));
        RooRealVar polysig("polysig","polysig",10,0,1000000000);
        RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));

        x.setRange("cut",0.44,0.56);

        RooFitResult* r_ks = sum.fitTo(data,Save(  ),Minos( kTRUE ),Range("cut"));
        //sum.fitTo(data,Range("cut"));
        //sum.fitTo(data,Range("cut"));
        //sum.fitTo(data,Range("cut"));
        //sum.fitTo(data,Range("cut"));
        //sum.fitTo(data,Range("cut"));

        double mean_ks = mean.getVal(  );
        double rms_ks = TMath::Sqrt( 0.5*sigma1.getVal(  )*sigma1.getVal(  ) + 0.5*sigma2.getVal(  )*sigma2.getVal(  ) );

        x.setRange( "peak", mean.getVal(  ) - 2*rms_ks, mean.getVal(  ) + 2*rms_ks);

        double gaus1F_ks = sig1.getVal(  );
        double gaus2F_ks = sig2.getVal(  );
        double polyF_ks = poly.getVal(  );

        RooAbsReal* Intgaus1_ks = gaus1.createIntegral( x, x,  "peak" );
        RooAbsReal* Intgaus2_ks = gaus2.createIntegral( x, x, "peak" );
        RooAbsReal* Intpoly_ks = poly.createIntegral( x, x, "peak" );

        double Intgaus1E_ks = gaus1F_ks*Intgaus1_ks->getVal(  );
        double Intgaus2E_ks = gaus2F_ks*Intgaus2_ks->getVal(  );
        double IntpolyE_ks = polyF_ks*Intpoly_ks->getVal(  );
        double totsig_ks = Intgaus1E_ks + Intgaus2E_ks + IntpolyE_ks;
        double sig_ks = Intgaus1E_ks + Intgaus2E_ks;

        double Fsig_ks = sig_ks/totsig_ks;

        mass_ks.push_back( mean_ks );
        std_ks.push_back( rms_ks );
        fsig_ks.push_back( Fsig_ks );
        covQual_ks.push_back( r_ks->covQual(  ) );

        cout << "adjusted poly integral ( ks ) " << IntpolyE_ks << endl;

        cout << "Norm Int Tot peak ( ks ) " << totsig_ks << endl;
        cout << "Norm Int poly peak ( ks )" << IntpolyE_ks << endl;


        cout << "Fsig ( ks ): " << Fsig_ks << endl;
        cout << "std ( ks ): " << rms_ks << endl;

        cout << "covQual ( ks )" << r_ks->covQual(  ) << endl;

        sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(0.5),LineColor(kBlue));
        sum.plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(0.5),LineColor(kBlue));
        cc1->cd(1);
        xframe->Draw();
        //tex->DrawLatex(0.59,0.87,label_mean[0]);
        //tex->DrawLatex(0.59,0.81,label_sigma[0]);


        if( lambda )
        {
            //lambda
            RooRealVar x("x","mass",1.08,1.155);
            RooDataHist data("data","dataset",x,massla);
            RooPlot* xframe = x.frame(160);
            xframe->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");
            xframe->GetYaxis()->SetTitle("Candidates / 0.0005GeV");
            xframe->GetXaxis()->CenterTitle(1);
            xframe->GetYaxis()->CenterTitle(1);
            xframe->GetXaxis()->SetTitleSize(xframe->GetXaxis()->GetTitleSize()*1.4);
            xframe->GetYaxis()->SetTitleSize(xframe->GetYaxis()->GetTitleSize()*1.3);
            xframe->GetXaxis()->SetTitleOffset(1);
            data.plotOn(xframe,Name("data"));
            RooRealVar mean("mean","mean",1.115,1.11,1.12);
            RooRealVar sigma1("sigma1","sigma1",0.005,0.001,0.01);
            RooRealVar sigma2("sigma2","sigma2",0.005,0.001,0.01);
            RooRealVar sig1("sig1","signal1",10,0,10000000);
            RooRealVar sig2("sig2","signal2",10,0,10000000);
            RooRealVar a("a","a",0,-100000,100000);
            RooRealVar b("b","b",0,-100000,100000);
            RooRealVar cp("cp","cp",0,-100000,100000);
            RooRealVar d("d","d",0,-100000,100000);
            RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
            RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
            RooPolynomial poly("poly","poly",x,RooArgList(a,b,cp,d));
            RooRealVar polysig("polysig","polysig",10,0,10000000);
            RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));

            x.setRange("cut",1.09,1.14);

            RooFitResult* r_la = sum.fitTo(data,Save(  ),Minos( kTRUE ),Range("cut"));
            //sum.fitTo(data,Range("cut"));
            //sum.fitTo(data,Range("cut"));
            //sum.fitTo(data,Range("cut"));
            //sum.fitTo(data,Range("cut"));
            //sum.fitTo(data,Range("cut"));

            double mean_la = mean.getVal(  );
            double rms_la = TMath::Sqrt( 0.5*sigma1.getVal(  )*sigma1.getVal(  ) + 0.5*sigma2.getVal(  )*sigma2.getVal(  ) );

            x.setRange( "peak", mean.getVal(  ) - 2*rms_la, mean.getVal(  ) + 2*rms_la);

            double gaus1F_la = sig1.getVal(  );
            double gaus2F_la = sig2.getVal(  );
            double polyF_la = poly.getVal(  );

            RooAbsReal* Intgaus1_la = gaus1.createIntegral( x, x,  "peak" );
            RooAbsReal* Intgaus2_la = gaus2.createIntegral( x, x, "peak" );
            RooAbsReal* Intpoly_la = poly.createIntegral( x, x, "peak" );

            double Intgaus1E_la = gaus1F_la*Intgaus1_la->getVal(  );
            double Intgaus2E_la = gaus2F_la*Intgaus2_la->getVal(  );
            double IntpolyE_la = polyF_la*Intpoly_la->getVal(  );
            double totsig_la = Intgaus1E_la + Intgaus2E_la + IntpolyE_la;
            double sig_la = Intgaus1E_la + Intgaus2E_la;

            double Fsig_la = sig_la/totsig_la;

            mass_la.push_back( mean_la );
            std_la.push_back( rms_la );
            fsig_la.push_back( Fsig_la );
            covQual_la.push_back( r_la->covQual(  ) );

            cout << "adjusted poly integral ( la ) " << IntpolyE_la << endl;
            cout << "Norm Int Tot peak ( la ) " << totsig_la << endl;
            cout << "Norm Int poly peak ( la )" << IntpolyE_la << endl;
            cout << "Fsig ( la ): " << Fsig_la << endl;
            cout << "std ( la ): " << rms_la << endl;
            cout << "covQual ( la )" << r_la->covQual(  ) << endl;


            sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(0.5),LineColor(kRed));
            sum.plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(0.5),LineColor(kRed));
            cc1->cd(2);
            xframe->Draw();
            //tex->DrawLatex(0.59,0.87,label_mean[1]);
            //tex->DrawLatex(0.59,0.81,label_sigma[1]);
        }

        cc1->cd(1);
        os << "Pt Bin: " << (pks[i]-1)/10 << " - " << pks[i+1]/10;
        tex->DrawLatex(0.22,0.8,os.str(  ).c_str(  ) );
        //tex->DrawLatex(0.22,0.81,label_n);
        //tex->DrawLatex(0.22,0.75,label_pid[0]);
        //tex->DrawLatex(0.22,0.69,"1 < p_{T} < 3 GeV/c");
        //tex->DrawLatex(0.22,0.63,"Preliminary");
        //tex->DrawLatex(0.15,0.58,label_cms[1]);
        os.str( std::string(  ) );

        cc1->cd(2);
        os << "Pt Bin: " << (pla[i]-1)/10 << " - " << pla[i+1]/10;
        tex->DrawLatex(0.22,0.8,os.str(  ).c_str(  ) );
        //tex->DrawLatex(0.22,0.81,label_n);
        //tex->DrawLatex(0.22,0.75,label_pid[1]);
        //tex->DrawLatex(0.22,0.69,"1 < p_{T} < 3 GeV/c");
        //tex->DrawLatex(0.22,0.63,"Preliminary");
        //tex->DrawLatex(0.15,0.58,label_cms[1]);
        os.str( std::string(  ) );

        //cc1->Print("massfit.pdf");
        //cc1->Print("massfit.gif");
        cc1->Print(Form( "massfit_%d.pdf",hbincounter ));
        i++; //to access correct bins
        hbincounter++;
    }
	pkscounter = 0;
	for( unsigned i=0; i<mass_ks.size(  ); i++ )
	{
        cout <<  "====================" << endl;
        cout << "Pt Bin: " << (pks[pkscounter]-1)/10 << " - " << pks[pkscounter+1]/10 << endl;
        cout <<  "====================" << endl;
        cout << "Mass_ks: " << mass_ks[i] << endl;
        cout << "Fsig_ks: " << fsig_ks[i] << endl;
        cout << "std_ks: " << std_ks[i] << endl;
        cout << "covQual_ks " << covQual_ks[i] << endl;
        pkscounter+=2;
	}

	placounter=6;
	for( unsigned i=0; i<mass_la.size(  ); i++ )
	{
        cout <<  "====================" << endl;
        cout << "Pt Bin: " << (pla[placounter]-1)/10 << " - " << pla[placounter+1]/10 << endl;
        cout <<  "====================" << endl;
        cout << "Mass_la: " << mass_la[i] << endl;
        cout << "Fsig_la: " << fsig_la[i] << endl;
        cout << "std_la: " << std_la[i] << endl;
        cout << "covQual_la " << covQual_la[i] << endl;
        placounter+=2;
    }
}
