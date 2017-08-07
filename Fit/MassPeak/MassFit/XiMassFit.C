#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdio>

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
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "TString.h"

#include <vector>

void XiMassFit()
{
    using namespace RooFit;
    const double massla = 1.115683;
    const double masspi = 0.13957018;
    gStyle->SetMarkerSize(0.8);

    RooMsgService::instance().setStreamStatus(0,kFALSE);
    RooMsgService::instance().setStreamStatus(1,kFALSE);

    //TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/CascadeV2pPb/RootFiles/Flow/CasCutLoose/CasCutLooseJL40.root");
    TFile* file = new TFile("/Volumes/MacHD/Users/blt1/research/RootFiles/Flow/Thesis/XiAnalysisCorrelationPtCut8TeVPD1_4_ForFinal.root");

    TH1D* massxi;

	TH2D* MassXi;

	int pTxiLength = 14; // the number of bins to be fitted is half of this number
	double pxi[] = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60};

	std::vector<double> mass_xi;
	std::vector<double> std_xi;
	std::vector<double> fsig_xi;
	std::vector<double> covQual_xi;

    //MassXi = (TH2D*)file->Get("XiMassPt/MassPt");
    MassXi = (TH2D*)file->Get("xiCorrelation/MassPt");


    std::ostringstream os;

	int pxicounter = 0; //for correct bin counting
	int hbincounter =1; //histogram bin counting
	for( int i=0; i<pTxiLength; i++ )
	{
			massxi = ( TH1D* )MassXi->ProjectionX( "massxi", pxi[i],pxi[i+1] );

			gStyle->SetOptTitle(kFALSE);

			TCanvas* cc1 = new TCanvas("cc1","cc1",600,450);

			TLatex* tex = new TLatex();
			tex->SetNDC();
			tex->SetTextFont(42);

			tex->SetTextSize(tex->GetTextSize()*0.95);

			//kshort
			RooRealVar x("x","mass",1.26,1.4);
			RooDataHist data("data","dataset",x,massxi);
			RooRealVar mean("mean","mean",1.32,1.29,1.33);
			RooRealVar sigma1("sigma1","sigma1",0.01,0.001,0.04);
			RooRealVar sigma2("sigma2","sigma2",0.01,0.001,0.04);
			RooRealVar sig1("sig1","signal1",10,0,1000000000);
			RooRealVar sig2("sig2","signal2",10,0,1000000000);
			RooGaussian gaus1("gaus1","gaus1",x,mean,sigma1);
			RooGaussian gaus2("gaus2","gaus2",x,mean,sigma2);
            RooRealVar qsig( "qsig","qsig",10,0,1000000000 );
            RooRealVar alpha( "alpha","alpha",1,0,10 );
            RooGenericPdf background( "background", "x - ( 1.115683 + 0.13957018 )^alpha", RooArgList(x,alpha ) );
            /*
            RooRealVar a( "a", "a", 1.0,-1,1 );
            RooRealVar b( "b", "b", 0.1,-1,1 );
            RooRealVar c( "c", "c", -0.1,-1,1 );
            RooChebychev background( "background", "background",x, RooArgList( a,b,c ) );
            */
			RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,background),RooArgList(sig2,sig2,qsig));

			x.setRange("cut",1.28,1.38);

			RooFitResult* r_xi = sum.fitTo(data,Save(  ),Minos( kTRUE ),Range("cut"));
			//sum.fitTo(data,Range("cut"));
			//sum.fitTo(data,Range("cut"));
			//sum.fitTo(data,Range("cut"));
			//sum.fitTo(data,Range("cut"));
			//sum.fitTo(data,Range("cut"));

			double mean_xi = mean.getVal(  );
			double rms_xi = TMath::Sqrt( 0.5*sigma1.getVal(  )*sigma1.getVal(  ) + 0.5*sigma2.getVal(  )*sigma2.getVal(  ) );

			x.setRange( "peak", mean.getVal(  ) - 2*rms_xi, mean.getVal(  ) + 2*rms_xi);

			double gaus1F_xi = sig1.getVal(  );
			double gaus2F_xi = sig2.getVal(  );
			double qsig_xi = qsig.getVal(  );

			RooAbsReal* Intgaus1_xi = gaus1.createIntegral( x, x,  "peak" );
			RooAbsReal* Intgaus2_xi = gaus2.createIntegral( x, x, "peak" );
			RooAbsReal* Intbackground_xi = background.createIntegral( x, x, "peak" );

			double Intgaus1E_xi = gaus1F_xi*Intgaus1_xi->getVal(  );
			double Intgaus2E_xi = gaus2F_xi*Intgaus2_xi->getVal(  );
			double IntbackgroundE_xi = qsig_xi*Intbackground_xi->getVal(  );
			double totsig_xi = Intgaus1E_xi + Intgaus2E_xi + IntbackgroundE_xi;
			double sig_xi = Intgaus1E_xi + Intgaus2E_xi;

			double Fsig_xi = sig_xi/totsig_xi;

			mass_xi.push_back( mean_xi );
			std_xi.push_back( rms_xi );
			fsig_xi.push_back( Fsig_xi );
			covQual_xi.push_back( r_xi->covQual(  ) );

			cout << "adjusted background integral ( xi ) " << IntbackgroundE_xi << endl;

			cout << "Norm Int Tot peak ( xi ) " << totsig_xi << endl;
			cout << "Norm Int background peak ( xi )" << IntbackgroundE_xi << endl;


			cout << "Fsig ( xi ): " << Fsig_xi << endl;
			cout << "std ( xi ): " << rms_xi << endl;
            cout << "mass ( xi ): " << mean_xi << endl;

			cout << "covQual ( xi )" << r_xi->covQual(  ) << endl;

			RooPlot* xframe = x.frame(270);
			xframe->GetXaxis()->SetTitle("invariant mass (GeV/c^{2})");
			xframe->GetYaxis()->SetTitle("Candidates / 0.001 GeV");
			xframe->GetXaxis()->CenterTitle(1);
			xframe->GetYaxis()->CenterTitle(1);
			xframe->GetXaxis()->SetTitleSize(xframe->GetXaxis()->GetTitleSize()*1.4);
			xframe->GetXaxis()->SetTitleOffset(1);
			xframe->GetYaxis()->SetTitleSize(xframe->GetYaxis()->GetTitleSize()*1.3);
			//xframe->GetYaxis()->SetTitleOffset(1);
			data.plotOn(xframe,Name("data"));
			sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(1),LineColor(kBlue));
			sum.plotOn(xframe,Components(background),NormRange("cut"),LineStyle(kDashed),LineWidth(1),LineColor(kBlue));
			cc1->cd(1);
			xframe->Draw();
			//tex->DrawLatex(0.59,0.87,label_mean[0]);
			//tex->DrawLatex(0.59,0.81,label_sigma[0]);



			cc1->cd(1);
			os << "Pt Bin: " << (pxi[i]-1)/10 << " - " << pxi[i+1]/10;
			tex->DrawLatex(0.22,0.8,os.str(  ).c_str(  ) );
			//tex->DrawLatex(0.22,0.81,label_n);
			//tex->DrawLatex(0.22,0.75,label_pid[0]);
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
	pxicounter = 0;
	for( int i=0; i<mass_xi.size(  ); i++ )
	{
			cout <<  "====================" << endl;
			cout << "Pt Bin: " << (pxi[pxicounter]-1)/10 << " - " << pxi[pxicounter+1]/10 << endl;
			cout <<  "====================" << endl;
			cout << "Mass_xi: " << mass_xi[i] << endl;
			cout << "Fsig_xi: " << fsig_xi[i] << endl;
			cout << "std_xi: " << std_xi[i] << endl;
			cout << "covQual_xi " << covQual_xi[i] << endl;
			pxicounter+=2;
	}
}
