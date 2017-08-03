//Includes
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cstdio>
#include <TROOT.h>
#include <TStyle.h>
#include <vector>
#include <sys/stat.h>

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
#include "TString.h"
#include "TVector.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "TString.h"

#include "MITStyle.C"

std::vector<std::vector<float> > transpose(const std::vector<std::vector<float> > data)
{
    std::vector<std::vector<float> > result(data[0].size(), std::vector<float>(data.size()));
    for(std::vector<float>::size_type i=0; i < data[0].size(); i++)
        for(std::vector<float>::size_type j=0; j < data.size(); j++)
        {
            result[i][j] = data[j][i];
        }
    return result;
}

void SigCalcXi(std::string name, int cutNum)// uses cut parameter name
{
    //Initializers
    bool NoCut = false;
    if(name == "NoCut") NoCut = true;
    using namespace RooFit;
    gStyle->SetMarkerSize(0.8);
    TFile* file = NULL;

    RooMsgService::instance().setStreamStatus(0,kFALSE);
    RooMsgService::instance().setStreamStatus(1,kFALSE);

    std::string rootfile = "XiTreeCut" + std::to_string(cutNum) + ".root";
    if(NoCut)
        file = new TFile("XiTreeNoCut.root");
    else
        file = new TFile(rootfile.c_str());

	int pTxiLength = 14; // the number of bins to be fitted is half of this number
	double pxi[] = {11,14, 15,18, 19,22, 23,28, 29,36, 37,46, 47,60};

    //Containers and variables
	std::vector<double> mass_xi;
	std::vector<double> std_xi;
	std::vector<double> fsig_xi;
	std::vector<double> covQual_xi;
    std::vector< std::vector<float> > significance_xi2; // first index is the parameter value and the second index is the significance at a given indexed pTBin
    std::ostringstream os;
    std::ostringstream fileos;
    std::ostringstream imageos;

    TH1D* massxi;
	TH2D* MassXi;
    std::string NoCutExt = "NoCut";
    std::string LogPath = "Log/";
    std::string PdfPath = "Plots/Peaks/";
    std::string GraphPath = "Plots/Graphs/";
    std::string filename = LogPath + name + std::to_string(cutNum) + ".txt";
    if(NoCut)
    {
        filename = LogPath + name + "_" + NoCutExt + ".txt";
    }

    //Cut Parameters For Hong's cuts uncomment last entries, increase numparam and oparamnum accordingly
    //float xi_xi3dipsig[]   = {8.0  , 8.5  , 9.0  , 9.5  , 10.0 , 10.5 , 11.0 , 11.5, 12.0, 12.5, 13.0, 13.5};//, 2.5};
    int numparam             = 13;
    int Oparamnum            = 12;
    //float xi_xi3dipsig[]     = {4.0  , 4.3  , 4.6  , 4.9  , 5.2  , 5.5  , 5.8  , 6.1 , 6.4 , 6.7 , 7.0, 7.3, 7.6};//, 2.5};
    //float xi_xipi3dipsig[]   = {5.0  , 4.9  , 4.8  , 4.7  , 4.6  , 4.5  , 4.4  , 4.3 , 4.2 , 4.1 , 4.0, 3.1, 3.3};//, 5.0};
    //float xi_vtrkpi3dipsig[] = {4.5  , 4.6  , 4.7  , 4.8  , 4.9  , 5.0  , 5.1  , 5.2 , 5.3 , 5.4 , 5.5, 3.0, 3.3};//, 4.0};
    //float xi_vtrkp3dipsig[]  = {2.5  , 2.6  , 2.7  , 2.8  , 2.9  , 3.0  , 3.1  , 3.2 , 3.3 , 3.4 , 3.5, 2.5, 2.5};//, 3.0};
    //float xi_xiflightsig[]   = {2.5  , 2.6  , 2.7  , 2.8  , 2.9  , 3.0  , 3.1  , 3.2 , 3.3 , 3.4 , 3.5, 2.5, 2.5};//, 3.0};
    //float xi_distancesig[]   = {12.0 , 11.8 , 11.6 , 11.4 , 11.0 , 10.5 , 10.0 , 9.5 , 9.0 , 8.5 , 8.0, 8.0, 8.5};//, 12.0};

    float xi_xi3dipsig[]     = {2.5  , 2.6  , 2.7  , 2.8  , 2.9  , 3.0  , 3.1  , 3.2 , 3.3 , 3.4 , 3.5, 3.5, 4.0};//, 2.5};
    float xi_xipi3dipsig[]   = {5.0  , 4.9  , 4.8  , 4.7  , 4.6  , 4.5  , 4.4  , 4.3 , 4.2 , 4.1 , 4.0, 3.1, 3.3};//, 5.0};
    float xi_vtrkpi3dipsig[] = {4.5  , 4.6  , 4.7  , 4.8  , 4.9  , 5.0  , 5.1  , 5.2 , 5.3 , 5.4 , 5.5, 3.0, 3.3};//, 4.0};
    float xi_vtrkp3dipsig[]  = {2.5  , 2.6  , 2.7  , 2.8  , 2.9  , 3.0  , 3.1  , 3.2 , 3.3 , 3.4 , 3.5, 2.5, 2.5};//, 3.0};
    float xi_xiflightsig[]   = {2.5  , 2.6  , 2.7  , 2.8  , 2.9  , 3.0  , 3.1  , 3.2 , 3.3 , 3.4 , 3.5, 2.5, 2.5};//, 3.0};
    float xi_distancesig[]   = {12.0 , 11.8 , 11.6 , 11.4 , 11.0 , 10.5 , 10.0 , 9.5 , 9.0 , 8.5 , 8.0, 8.0, 8.5};//, 12.0};

    // Open output txt file
    ofstream myfile;
    struct stat buffer;
    if(stat(filename.c_str(), &buffer) == 0)
    {
        cout << "File with this name already exists, will append beginning with line OVERWRITE" << endl;
        myfile.open(filename.c_str(), std::ios_base::app);
        myfile << "OVERWRITE OVERWRITE OVERWRITE\n";
    }
    else
    {
        myfile.open(filename.c_str());
    }

    int pxicounter = 0; //for correct bin counting
    for(int hlooper=0; hlooper<numparam; hlooper++)
    {
        //hlooper=numparam-1;
        if(name == "NoCut")
        {
            hlooper = numparam-1;
            fileos << "hxi_" << name;
        }
        if(name == "xi3dipsig") fileos << "hxi_" << name << "_" << std::fixed << std::setprecision(1) << xi_xi3dipsig[hlooper];
        if(name == "xipi3dipsig") fileos << "hxi_" << name << "_" << std::fixed << std::setprecision(1) << xi_xipi3dipsig[hlooper];
        if(name == "vtrkpi3dipsig") fileos << "hxi_" << name << "_" << std::fixed << std::setprecision(1) << xi_vtrkpi3dipsig[hlooper];
        if(name == "vtrkp3dipsig") fileos << "hxi_" << name << "_" << std::fixed << std::setprecision(1) << xi_vtrkp3dipsig[hlooper];
        if(name == "xiflightsig") fileos << "hxi_" << name << "_" << std::fixed << std::setprecision(1) << xi_xiflightsig[hlooper];
        if(name == "distancesig") fileos << "hxi_" << name << "_" << std::fixed << std::setprecision(1) << xi_distancesig[hlooper];
        MassXi = (TH2D*)file->Get(fileos.str().c_str());
        cout << fileos.str() << endl;

        std::vector<float> significance_xi;

        for( int i=0; i<pTxiLength; i++ )
        {
            massxi = ( TH1D* )MassXi->ProjectionX( "massxi", pxi[i],pxi[i+1] );

            gStyle->SetOptTitle(kFALSE);

            TCanvas* cc1 = new TCanvas("cc1","cc1",600,450);

            TLatex* tex = new TLatex();
            tex->SetNDC();
            tex->SetTextFont(42);

            tex->SetTextSize(tex->GetTextSize()*0.95);

            //Fit
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

			//RooRealVar a("a","a",0,-100000,100000);
			RooRealVar a("a","a",0,-100000,100000);
			RooRealVar b("b","b",0,-100000,100000);
            RooRealVar cp("cp","cp",0,-100000,100000);
            RooRealVar d("d","d",0,-100000,100000);
			RooRealVar polysig("polysig","polysig",10,0,1000000000);
            RooPolynomial poly("poly","poly",x,RooArgList(a,b,cp));
            RooPolynomial poly4("poly","poly",x,RooArgList(a,b,cp,d));
            RooAddPdf sum("sum","sum",RooArgList(gaus1,gaus2,background),RooArgList(sig1,sig2,qsig));
            RooAddPdf sumNoCut("sumNoCut","sumNoCut",RooArgList(gaus1,gaus2,poly),RooArgList(sig1,sig2,polysig));
            RooAddPdf sumP4("sumNoCut","sumNoCut",RooArgList(gaus1,gaus2,poly4),RooArgList(sig1,sig2,polysig));

            x.setRange("cut",1.295,1.35);
            x.setRange("tight",1.3,1.35);

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
            RooFitResult* r_xi = NULL;
            if(NoCut)
            {
                r_xi = sumNoCut.fitTo(data,Save(  ),Minos( kTRUE ),Range("tight"));
                sumNoCut.plotOn(xframe,Name("sumNoCut"),NormRange("tight"),LineWidth(1),LineColor(kBlue));
                sumNoCut.plotOn(xframe,Components(poly),NormRange("tight"),LineStyle(kDashed),LineWidth(1),LineColor(kBlue));
            }
            else
            {
                //if(i==0)
                //{
                    //r_xi = sumNoCut.fitTo(data,Save(  ),Minos( kTRUE ),Range("cut"));
                    //sumNoCut.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(1),LineColor(kBlue));
                    //sumNoCut.plotOn(xframe,Components(poly),NormRange("cut"),LineStyle(kDashed),LineWidth(1),LineColor(kBlue));
                //}
                //else
                //{
                    //r_xi = sum.fitTo(data,Save(  ),Minos( kTRUE ),Range("cut"));
                    //sum.plotOn(xframe,Name("sum"),NormRange("cut"),LineWidth(1),LineColor(kBlue));
                    //sum.plotOn(xframe,Components(background),NormRange("cut"),LineStyle(kDashed),LineWidth(1),LineColor(kBlue));
                //}
                r_xi = sumP4.fitTo(data,Save(  ),Minos( kTRUE ),Range("cut"));
                sumP4.plotOn(xframe,Name("sumNoCut"),NormRange("cut"),LineWidth(1),LineColor(kBlue));
                sumP4.plotOn(xframe,Components(poly4),NormRange("cut"),LineStyle(kDashed),LineWidth(1),LineColor(kBlue));
            }
            cc1->cd(1);
            xframe->Draw();
            //tex->DrawLatex(0.59,0.87,label_mean[0]);
            //tex->DrawLatex(0.59,0.81,label_sigma[0]);

            //Calculation
            double Intgaus1E_xi      = -999;
            double Intgaus2E_xi      = -999;
            double IntbackgroundE_xi = -999;
            double totsig_xi         = -999;
            double sig_xi            = -999;
            double Fsig_xi           = -999;
            double Significance      = -999;

            RooAbsReal* Intbackground_xi = NULL;

            double mean_xi = mean.getVal(  );
            double rms_xi  = TMath::Sqrt( 0.5*sigma1.getVal(  )*sigma1.getVal(  ) + 0.5*sigma2.getVal(  )*sigma2.getVal(  ) );

            double gaus1F_xi  = sig1.getVal(  );
            double gaus2F_xi  = sig2.getVal(  );
            double qsig_xi    = qsig.getVal(  );
            double polysig_xi = polysig.getVal();

            x.setRange("int",mean_xi - 2*rms_xi,mean_xi + 2*rms_xi);

            RooAbsReal* Intgaus1_xi = gaus1.createIntegral( x, x,  "int" );
            RooAbsReal* Intgaus2_xi = gaus2.createIntegral( x, x, "int" );

            if(NoCut)
            {
                RooAbsReal* Intbackground_xi = poly.createIntegral( x, x, "int" );

                Intgaus1E_xi      = gaus1F_xi*Intgaus1_xi->getVal(  );
                Intgaus2E_xi      = gaus2F_xi*Intgaus2_xi->getVal(  );
                IntbackgroundE_xi = polysig_xi*Intbackground_xi->getVal(  );
                totsig_xi         = Intgaus1E_xi + Intgaus2E_xi + IntbackgroundE_xi;
                sig_xi            = Intgaus1E_xi + Intgaus2E_xi;

                Fsig_xi           = sig_xi/totsig_xi;
                Significance      = sig_xi/(TMath::Sqrt(totsig_xi));
            }
            else
            {
                RooAbsReal* Intbackground_xi = background.createIntegral( x, x, "int" );

                Intgaus1E_xi      = gaus1F_xi*Intgaus1_xi->getVal(  );
                Intgaus2E_xi      = gaus2F_xi*Intgaus2_xi->getVal(  );
                //IntbackgroundE_xi = qsig_xi*Intbackground_xi->getVal(  );
                IntbackgroundE_xi = polysig_xi*Intbackground_xi->getVal(  );
                totsig_xi         = Intgaus1E_xi + Intgaus2E_xi + IntbackgroundE_xi;
                sig_xi            = Intgaus1E_xi + Intgaus2E_xi;

                Fsig_xi           = sig_xi/totsig_xi;
                Significance      = sig_xi/(TMath::Sqrt(totsig_xi));
            }

            mass_xi.push_back( mean_xi );
            std_xi.push_back( rms_xi );
            fsig_xi.push_back( Fsig_xi );
            covQual_xi.push_back( r_xi->covQual(  ) );
            significance_xi.push_back(Significance);

            cout << "adjusted background integral ( xi ) " << IntbackgroundE_xi << endl;

            cout << "Norm Int Tot peak ( xi ) " << totsig_xi << endl;
            cout << "Norm Int background peak ( xi )" << IntbackgroundE_xi << endl;

            cout << "Fsig ( xi ): " << Fsig_xi << endl;
            cout << "std ( xi ): " << rms_xi << endl;
            cout << "mass ( xi ): " << mean_xi << endl;

            cout << "covQual ( xi )" << r_xi->covQual(  ) << endl;

            // Make Files
            cc1->cd(1);
            os << "Pt Bin: " << (pxi[i]-1)/10 << " - " << pxi[i+1]/10;
            tex->DrawLatex(0.22,0.8,os.str(  ).c_str(  ) );
            os.str(std::string());

            tex->SetTextSize(tex->GetTextSize()*0.80);

            os << "Significance: " << Significance;
            tex->DrawLatex(0.15,0.75,os.str().c_str());
            os.str(std::string());

            if(!NoCut)
            {
                if(name == "xi3dipsig")
                {
                    os <<"(2.5) Xi3dipsig < " << xi_xi3dipsig[hlooper];
                    tex->DrawLatex(0.63,0.85,os.str().c_str());
                    os.str(std::string());
                    os << "(5.0) Xipi3dipsig > " << xi_xipi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.81,os.str().c_str());
                    os.str(std::string());
                    os << "(4.0) vtrkpi3dipsig > " << xi_vtrkpi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.77,os.str().c_str());
                    os.str(std::string());
                    os << "(3.0) vtrkp3dipsig > " << xi_vtrkp3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.73,os.str().c_str());
                    os.str(std::string());
                    os << "(3.0) xiflightsig > " << xi_xiflightsig[Oparamnum];
                    tex->DrawLatex(0.63,0.69,os.str().c_str());
                    os.str(std::string());
                    os << "(12) distancesig > " << xi_distancesig[Oparamnum];
                    tex->DrawLatex(0.63,0.65,os.str().c_str());
                    os.str(std::string());
                }
                else if(name == "xipi3dipsig")
                {
                    os <<"(2.5) Xi3dipsig < " << xi_xi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.85,os.str().c_str());
                    os.str(std::string());
                    os << "(5.0) Xipi3dipsig > " << xi_xipi3dipsig[hlooper];
                    tex->DrawLatex(0.63,0.81,os.str().c_str());
                    os.str(std::string());
                    os << "(4.0) vtrkpi3dipsig > " << xi_vtrkpi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.77,os.str().c_str());
                    os.str(std::string());
                    os << "(3.0) vtrkp3dipsig > " << xi_vtrkp3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.73,os.str().c_str());
                    os.str(std::string());
                    os << "(3.0) xiflightsig > " << xi_xiflightsig[Oparamnum];
                    tex->DrawLatex(0.63,0.69,os.str().c_str());
                    os.str(std::string());
                    os << "(12) distancesig > " << xi_distancesig[Oparamnum];
                    tex->DrawLatex(0.63,0.65,os.str().c_str());
                    os.str(std::string());
                }
                else if(name == "vtrkpi3dipsig")
                {
                    os <<"(2.5) Xi3dipsig < " << xi_xi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.85,os.str().c_str());
                    os.str(std::string());
                    os << "(5.0) Xipi3dipsig > " << xi_xipi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.81,os.str().c_str());
                    os.str(std::string());
                    os << "(4.0) vtrkpi3dipsig > " << xi_vtrkpi3dipsig[hlooper];
                    tex->DrawLatex(0.63,0.77,os.str().c_str());
                    os.str(std::string());
                    os << "(3.0) vtrkp3dipsig > " << xi_vtrkp3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.73,os.str().c_str());
                    os.str(std::string());
                    os << "(3.0) xiflightsig > " << xi_xiflightsig[Oparamnum];
                    tex->DrawLatex(0.63,0.69,os.str().c_str());
                    os.str(std::string());
                    os << "(12) distancesig > " << xi_distancesig[Oparamnum];
                    tex->DrawLatex(0.63,0.65,os.str().c_str());
                    os.str(std::string());
                }
                else if(name == "vtrkp3dipsig")
                {
                    os <<"(2.5) Xi3dipsig < " << xi_xi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.85,os.str().c_str());
                    os.str(std::string());
                    os << "(5.0) Xipi3dipsig > " << xi_xipi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.81,os.str().c_str());
                    os.str(std::string());
                    os << "(4.0) vtrkpi3dipsig > " << xi_vtrkpi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.77,os.str().c_str());
                    os.str(std::string());
                    os << "(3.0) vtrkp3dipsig > " << xi_vtrkp3dipsig[hlooper];
                    tex->DrawLatex(0.63,0.73,os.str().c_str());
                    os.str(std::string());
                    os << "(3.0) xiflightsig > " << xi_xiflightsig[Oparamnum];
                    tex->DrawLatex(0.63,0.69,os.str().c_str());
                    os.str(std::string());
                    os << "(12) distancesig > " << xi_distancesig[Oparamnum];
                    tex->DrawLatex(0.63,0.65,os.str().c_str());
                    os.str(std::string());
                }
                else if(name == "xiflightsig")
                {
                    os <<"(2.5) Xi3dipsig < " << xi_xi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.85,os.str().c_str());
                    os.str(std::string());
                    os << "(5.0) Xipi3dipsig > " << xi_xipi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.81,os.str().c_str());
                    os.str(std::string());
                    os << "(4.0) vtrkpi3dipsig > " << xi_vtrkpi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.77,os.str().c_str());
                    os.str(std::string());
                    os << "(3.0) vtrkp3dipsig > " << xi_vtrkp3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.73,os.str().c_str());
                    os.str(std::string());
                    os << "(3.0) xiflightsig > " << xi_xiflightsig[hlooper];
                    tex->DrawLatex(0.63,0.69,os.str().c_str());
                    os.str(std::string());
                    os << "(12) distancesig > " << xi_distancesig[Oparamnum];
                    tex->DrawLatex(0.63,0.65,os.str().c_str());
                    os.str(std::string());
                }
                else if(name == "distancesig")
                {
                    os <<"(2.5) Xi3dipsig < " << xi_xi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.85,os.str().c_str());
                    os.str(std::string());
                    os << "(5.0) Xipi3dipsig > " << xi_xipi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.81,os.str().c_str());
                    os.str(std::string());
                    os << "(4.0) vtrkpi3dipsig > " << xi_vtrkpi3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.77,os.str().c_str());
                    os.str(std::string());
                    os << "(3.0) vtrkp3dipsig > " << xi_vtrkp3dipsig[Oparamnum];
                    tex->DrawLatex(0.63,0.73,os.str().c_str());
                    os.str(std::string());
                    os << "(3.0) xiflightsig > " << xi_xiflightsig[Oparamnum];
                    tex->DrawLatex(0.63,0.69,os.str().c_str());
                    os.str(std::string());
                    os << "(12) distancesig > " << xi_distancesig[hlooper];
                    tex->DrawLatex(0.63,0.65,os.str().c_str());
                    os.str(std::string());
                }
            }

            //Make PDF
            if(i==0)
            {
                if(NoCut) imageos << NoCutExt << ".pdf(";
                if(name == "xi3dipsig") imageos << PdfPath << name << "_" << xi_xi3dipsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf(";
                if(name == "xipi3dipsig") imageos << PdfPath << name << "_" << xi_xipi3dipsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf(";
                if(name == "vtrkpi3dipsig") imageos << PdfPath << name << "_" << xi_vtrkpi3dipsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf(";
                if(name == "vtrkp3dipsig") imageos << PdfPath << name << "_" << xi_vtrkp3dipsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf(";
                if(name == "xiflightsig") imageos << PdfPath << name << "_" << xi_xiflightsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf("; //changeline
                if(name == "distancesig") imageos << PdfPath << name << "_" << xi_distancesig[hlooper]*10 << "_cutNum" << cutNum << ".pdf("; //changeline
                cc1->Print(imageos.str().c_str(),"pdf");
            }
            else if(i != (pTxiLength-2))
            {
                if(NoCut) imageos << NoCutExt << ".pdf";
                if(name == "xi3dipsig") imageos << PdfPath << name << "_" << xi_xi3dipsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf";
                if(name == "xipi3dipsig") imageos << PdfPath << name << "_" << xi_xipi3dipsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf";
                if(name == "vtrkpi3dipsig") imageos << PdfPath << name << "_" << xi_vtrkpi3dipsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf";
                if(name == "vtrkp3dipsig") imageos << PdfPath << name << "_" << xi_vtrkp3dipsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf";
                if(name == "xiflightsig") imageos << PdfPath << name << "_" << xi_xiflightsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf"; //changeline
                if(name == "distancesig") imageos << PdfPath << name << "_" << xi_distancesig[hlooper]*10 << "_cutNum" << cutNum << ".pdf"; //changeline
                cc1->Print(imageos.str().c_str(),"pdf");
            }
            else
            {
                if(NoCut) imageos << NoCutExt << ".pdf)";
                if(name == "xi3dipsig") imageos << PdfPath << name << "_" << xi_xi3dipsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf)";
                if(name == "xipi3dipsig") imageos << PdfPath << name << "_" << xi_xipi3dipsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf)";
                if(name == "vtrkpi3dipsig") imageos << PdfPath << name << "_" << xi_vtrkpi3dipsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf)";
                if(name == "vtrkp3dipsig") imageos << PdfPath << name << "_" << xi_vtrkp3dipsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf)";
                if(name == "xiflightsig") imageos << PdfPath << name << "_" << xi_xiflightsig[hlooper]*10 << "_cutNum" << cutNum << ".pdf)"; //changeline
                if(name == "distancesig") imageos << PdfPath << name << "_" << xi_distancesig[hlooper]*10 << "_cutNum" << cutNum << ".pdf)"; //changeline
                cc1->Print(imageos.str().c_str(),"pdf");
            }

            i++; //to access correct bins
            imageos.str(std::string());
        }
        significance_xi2.push_back(significance_xi);
        fileos.str(std::string());
    }
    pxicounter = 0;
    int counter = 0;
    for( unsigned i=0; i<mass_xi.size(  ); i++ )
    {
        int index = counter/14;
        cout <<  "====================" << endl;
        cout << name << endl;
        cout << "Pt Bin: " << (pxi[pxicounter]-1)/10 << " - " << pxi[pxicounter+1]/10 << endl;
        cout <<  "====================" << endl;
        cout << "Mass_xi: " << mass_xi[i] << endl;
        cout << "Fsig_xi: " << fsig_xi[i] << endl;
        cout << "std_xi: " << std_xi[i] << endl;
        cout << "covQual_xi " << covQual_xi[i] << endl;
        cout << "Significance: " << significance_xi2[index][pxicounter/2] << endl;

        os.str(std::string());
        if(name == "xi3dipsig") os << name << ": " << xi_xi3dipsig[index];
        if(name == "xipi3dipsig") os << name << ": " << xi_xipi3dipsig[index];
        if(name == "vtrkpi3dipsig") os << name << ": " << xi_vtrkpi3dipsig[index];
        if(name == "vtrkp3dipsig") os << name << ": " << xi_vtrkp3dipsig[index];
        if(name == "xiflightsig") os << name << ": " << xi_xiflightsig[index];
        if(name == "distancesig") os << name << ": " << xi_distancesig[index];
        myfile << os.str() + "\n";
        os.str(std::string());
        os << "Pt Bin: " << (pxi[pxicounter]-1)/10 << " - " << pxi[pxicounter+1]/10;
        myfile << os.str() + "\n";
        os.str(std::string());
        os << "Significance: " << significance_xi2[index][pxicounter/2];
        myfile << os.str() + "\n";
        pxicounter+=2;
        counter+=2;
        if(pxicounter==14)
        {
            pxicounter=0;
            myfile << "\n";
        }
    }

    //Plot Significance
    const int numBins = pTxiLength/2;

    std::vector<std::vector<float> > transignificance_xi2 = transpose(significance_xi2);

    float* a0 = &transignificance_xi2[0][0];
    //float* a1 = &transignificance_xi2[0][1];
    //float* a2 = &transignificance_xi2[0][2];
    //float* a3 = &transignificance_xi2[0][3];
    //float* a4 = &transignificance_xi2[0][4];
    //float* a5 = &transignificance_xi2[0][5];
    //float* a6 = &transignificance_xi2[0][6];

    TGraph* distancesigplot0 = new TGraph(numparam,xi_distancesig, a0);
    //TGraph* distancesigplot1 = new TGraph(numparam,xi_distancesig, a1);
    //TGraph* distancesigplot2 = new TGraph(numparam,xi_distancesig, a2);
    //TGraph* distancesigplot3 = new TGraph(numparam,xi_distancesig, a3);
    //TGraph* distancesigplot4 = new TGraph(numparam,xi_distancesig, a4);
    //TGraph* distancesigplot5 = new TGraph(numparam,xi_distancesig, a5);
    //TGraph* distancesigplot6 = new TGraph(numparam,xi_distancesig, a6);

    //TCanvas* c1 = new TCanvas("c1","distanceSig",800,800);
    //c1->cd();

    MITStyle();
    TCanvas* distanceSigPlot = MakeCanvas("distanceSigPlot","Plot");
    distanceSigPlot->cd();
    distanceSigPlot->SetLeftMargin(0.12);

    TH1F* frame = distanceSigPlot->DrawFrame( 0,-0.01,15,30);
    gPad->SetTickx(  );
    gPad->SetTicky(  );
    frame->GetXaxis(  )->CenterTitle( 1 );
    frame->GetYaxis(  )->CenterTitle( 1 );
    frame->GetXaxis(  )->SetTitleSize( 0.05 );
    frame->GetXaxis(  )->SetTitle( "Parameter Value" );
    frame->GetYaxis(  )->SetTitle( "Significance" );
    frame->GetYaxis(  )->SetTitleSize( 0.05 );
    frame->SetTitleOffset( 1.2,"Y" );
    frame->SetTitleOffset( 1.2,"X" );

    distancesigplot0->SetMarkerStyle(33);
    distancesigplot0->SetMarkerSize(1.5);
    //distancesigplot1->SetMarkerStyle(33);
    //distancesigplot1->SetMarkerSize(1.5);
    //distancesigplot2->SetMarkerStyle(33);
    //distancesigplot2->SetMarkerSize(1.5);
    //distancesigplot3->SetMarkerStyle(33);
    //distancesigplot3->SetMarkerSize(1.5);
    //distancesigplot4->SetMarkerStyle(33);
    //distancesigplot4->SetMarkerSize(1.5);
    //distancesigplot5->SetMarkerStyle(33);
    //distancesigplot5->SetMarkerSize(1.5);
    //distancesigplot6->SetMarkerStyle(33);
    //distancesigplot6->SetMarkerSize(1.5);

    distancesigplot0->Draw("PSAME");
    //distancesigplot1->Draw("PSAME");
    //distancesigplot2->Draw("PSAME");
    //distancesigplot3->Draw("PSAME");
    //distancesigplot4->Draw("PSAME");
    //distancesigplot5->Draw("PSAME");
    //distancesigplot6->Draw("PSAME");
}
