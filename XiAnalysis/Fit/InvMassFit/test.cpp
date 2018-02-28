#include <TLatex.h>
#include <TStyle.h>
#include "TF1.h"
#include "TMath.h"
#include "TVirtualFitter.h"
#include "TLegend.h"
#include "TMathText.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TPad.h"
#include "TFile.h"
#include <TString.h>
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGaxis.h"

#include <vector>
#include <array>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <cstring>
#include <cstdio>

#define piMass 0.13957018
#define lambdaMass 1.115683
#define PI 3.1416

void test(  )
{
    int numFiles = 0;
    int FN[] = {1,2};
    for( const int &FNs : FN ){
       numFiles++; 
    }
    std::cout << numFiles << std::endl;
}
