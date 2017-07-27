#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <fstream>
#include <Riostream.h>
#include <sstream>
#include <fstream>
#include <vector>
using namespace std;

void TextToVector(TString filename, vector<vector < float > >& Data)
{
  ifstream f(filename.Data());
  string s;
  while(getline(f,s))
  {
    istringstream iss(s);
    vector<float> Line;
    float value;
    while (iss >> value) Line.push_back(value);
    Data.push_back(Line);
  }
}

void TextToArray(string filename, float &array[])
{
  ifstream f(filename.Data());
  string s;
  int i=0;
  while(getline(f,s))
  {
    istringstream iss(s);
    float value;
    while (iss >> value) {array[i] = value; i++;}
  }
}

void GetNtupleFromFile(TString txtFileName, TNtuple* ntuple, const int ndim)
{
  ifstream f(txtFileName.Data());
  string s;
  while(getline(f,s))
  {
    istringstream iss(s);
    float* data = new float[ndim];
    for(int i=0;i<ndim;i++)
      iss >> data[i];
    ntuple->Fill(data);
    delete data;
  }
}

TGraph* GetGraphFromFile(TString txtFileName, Color_t markerColor=1, Style_t markerStyle=20, Size_t markerSize=1, Style_t lineStyle=1)
{
  Float_t x_array[400],y_array[400];
  Char_t buffer[256];
  Float_t x,y;
  Int_t nlines = 0;
  ifstream infile(txtFileName.Data());

  if (!infile.is_open()) {
    cout << "Error opening file. Exiting." << endl;
  } else {
    while (!infile.eof()) {
      infile.getline(buffer,100);
      sscanf(buffer,"%f %f\n",&x,&y);
      x_array[nlines] = x;
      y_array[nlines] = y;
      nlines++;
    }
  }
  TGraph *graph =
    new TGraph(nlines,x_array,y_array);
  txtFileName.Remove(txtFileName.Index(".txt"),4);
  graph->SetName(txtFileName.Data());
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerColor(markerColor);
  graph->SetLineStyle(lineStyle);
  graph->SetLineColor(markerColor);
  graph->SetMarkerSize(markerSize);
  graph->SetLineWidth(3);
  return graph;
}

TGraphErrors* GetGraphWithSymmErrorsFromFile(TString txtFileName, Color_t markerColor=1, Style_t markerStyle=20, Size_t markerSize=1, Style_t lineStyle=1)
{
  Float_t x_array[400],ex_array[400],y_array[400],ey_array[400];
  Char_t buffer[2560];
  Float_t x,y,ex,ey;
  Int_t nlines = 0;
  ifstream infile(txtFileName.Data());

  if (!infile.is_open()) {
    cout << "Error opening file. Exiting." << endl;
  } else {
    while (!infile.eof()) {
      infile.getline(buffer,100);
      sscanf(buffer,"%.6f %.6f %.6f %.6f\n",&x,&ex,&y,&ey);
      x_array[nlines] = x;
      ex_array[nlines] = ex;
      y_array[nlines] = y;
      ey_array[nlines] = ey;
      nlines++;
    }
  }
  TGraphErrors *graph =
    new TGraphErrors(nlines,x_array,y_array,ex_array,ey_array);
  txtFileName.Remove(txtFileName.Index(".txt"),4);
  graph->SetName(txtFileName.Data());
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerColor(markerColor);
  graph->SetLineStyle(lineStyle);
  graph->SetLineColor(markerColor);
  graph->SetMarkerSize(markerSize);
  graph->SetLineWidth(1);
  return graph;
}

TGraphAsymmErrors* GetGraphWithErrorsFromFile(TString txtFileName, Color_t markerColor=1, Style_t markerStyle=20, Size_t markerSize=1, Style_t lineStyle=1)
{

  Float_t x_array[400],ex_array[400],y_array[400],eyh_array[400],eyl_array[400];
  Char_t buffer[2560];
  Float_t x,y,eyl,eyh;
  Int_t nlines = 0;
  ifstream infile(txtFileName.Data());

  if (!infile.is_open()) {
    cout << "Error opening file. Exiting." << endl;
  } else {
    while (!infile.eof()) {
      infile.getline(buffer,2000);
      sscanf(buffer,"%f %f %f %f\n",&x,&y,&eyh,&eyl);
      x_array[nlines] = x;
      ex_array[nlines] = 0.0;
      y_array[nlines] = y;
      eyl_array[nlines] = eyl;
      eyh_array[nlines] = eyh;
      nlines++;
    }
  }

  TGraphAsymmErrors *graph = 
    new TGraphAsymmErrors(nlines,x_array,y_array,ex_array,ex_array,eyl_array,eyh_array);
  txtFileName.Remove(txtFileName.Index(".dat"),4);
  graph->SetName(txtFileName.Data()); 
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerColor(markerColor);
  graph->SetLineStyle(lineStyle);
  graph->SetLineColor(markerColor);
  graph->SetMarkerSize(markerSize);
  graph->SetLineWidth(3);  
  return graph;
}

TGraphErrors* GetGraphWithSymmYErrorsFromFile(TString txtFileName, Color_t markerColor=1, Style_t markerStyle=20, Size_t markerSize=1, Style_t lineStyle=1, bool IsNoErr=0)
{
  Float_t x_array[400],ex_array[400],y_array[400],ey_array[400];
  Char_t buffer[2048];
  Float_t x,y,ex,ey;
  Int_t nlines = 0;
  ifstream infile(txtFileName.Data());

  if (!infile.is_open()) {
    cout << "Error opening file. Exiting." << endl;
  } else {
    while (!infile.eof()) {
      infile.getline(buffer,2048);
      sscanf(buffer,"%f %f %f\n",&x,&y,&ey);
      x_array[nlines] = x;
      ex_array[nlines] = 0;
      y_array[nlines] = y;
      ey_array[nlines] = ey;
      if(IsNoErr) ey_array[nlines]=0;
      nlines++;
    }
  }
  TGraphErrors *graph =
    new TGraphErrors(nlines-1,x_array,y_array,ex_array,ey_array);
  txtFileName.Remove(txtFileName.Index(".txt"),4);
  graph->SetName(txtFileName.Data());
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerColor(markerColor);
  graph->SetLineStyle(lineStyle);
  graph->SetLineColor(markerColor);
  graph->SetMarkerSize(markerSize);
  graph->SetLineWidth(1);
  return graph;
}

TGraph* GetBandWithSymmYErrorsFromFile(TString txtFileName, Color_t fillColor=1, Style_t lineStyle=1)
{
  Float_t x_array[400],ex_array[400],y_array[400],ey_array[400];
  Char_t buffer[2048];
  Float_t x,y,ex,ey;
  Int_t nlines = 0;
  ifstream infile(txtFileName.Data());

  if (!infile.is_open()) {
    cout << "Error opening file. Exiting." << endl;
  } else {
    while (!infile.eof()) {
      infile.getline(buffer,2048);
      sscanf(buffer,"%f %f %f\n",&x,&y,&ey);
      x_array[nlines] = x;
      y_array[nlines] = y+ey;
      ey_array[nlines] = ey;
      nlines++;
    }
  }

  for(int i=nlines;i<2*nlines;i++)
  {
      x_array[i] = x_array[2*nlines-i-1]; 
      y_array[i] = y_array[2*nlines-i-1]-2*ey_array[2*nlines-i-1];
  }   

  TGraph *graph =
    new TGraph(2*nlines,x_array,y_array);
  txtFileName.Remove(txtFileName.Index(".txt"),4);
  graph->SetName(txtFileName.Data());
  graph->SetFillColor(fillColor);
  graph->SetLineWidth(3);
  graph->SetLineStyle(lineStyle);
  graph->SetLineColor(fillColor);
  return graph;
}

TH1D* GetHistogramFromFile(TString txtFileName, Color_t markerColor=1, Style_t markerStyle=20, Size_t markerSize=1, Style_t lineStyle=1)
{
  Float_t x_array[400],ex_array[400],y_array[400],ey_array[400];
  Char_t buffer[1024];
  Float_t x,y,ex,ey;
  Int_t nlines = 0;
  ifstream infile(txtFileName.Data());

  if (!infile.is_open()) {
    cout << "Error opening file. Exiting." << endl;
  } else {
    while (!infile.eof()) {
      infile.getline(buffer,100);
      sscanf(buffer,"%f %f %f\n",&x,&y,&ey);
      x_array[nlines] = x;
      ex_array[nlines] = 0;
      y_array[nlines] = y;
      ey_array[nlines] = ey;
      nlines++;
    }
  }

  double binwidth = x_array[1]-x_array[0];

  TH1D* h = new TH1D("h","h",nlines,x_array[0]-binwidth/2.0,x_array[nlines-1]+binwidth/2.0);
  for(int i=1;i<=h->GetNbinsX();i++)
  {
    h->SetBinContent(i,y_array[i-1]);
    h->SetBinError(i,ey_array[i-1]);
  }
  txtFileName.Remove(txtFileName.Index(".dat"),4);
  h->SetName(txtFileName.Data());
  h->SetMarkerStyle(markerStyle);
  h->SetMarkerColor(markerColor);
  h->SetLineStyle(lineStyle);
//  h->SetLineColor(markerColor);
  h->SetMarkerSize(markerSize);
  h->SetLineWidth(3);
  return h;
}

void GraphScaleShift(TGraphErrors* gr, double scale=1.0, double shift=0.0)
{
  double x,y,yerr;

  for(int i=0;i<gr->GetN();i++)
  {
    gr->GetPoint(i,x,y);
    yerr=gr->GetErrorY(i);
    gr->SetPoint(i,x,y/scale+shift);
    gr->SetPointError(i,0.0,yerr/scale);
  }
}

void GraphScale2DShift(TGraphErrors* gr, double scale_X=1.0, double scale_Y=1.0, double shift=0.0)
{
  double x,y,yerr;

  for(int i=0;i<gr->GetN();i++)
  {
    gr->GetPoint(i,x,y);
    yerr=gr->GetErrorY(i);
    gr->SetPoint(i,x/scale_X,y/scale_Y+shift);
    gr->SetPointError(i,0.0,yerr/scale_Y);
  }
}

void GraphScale2DShiftNoErr(TGraph* gr, double scale_X=1.0, double scale_Y=1.0, double shift=0.0)
{
  double x,y,yerr;

  for(int i=0;i<gr->GetN();i++)
  {
    gr->GetPoint(i,x,y);
    gr->SetPoint(i,x/scale_X,y/scale_Y+shift);
  }
}
void GraphChangeXaxis(TGraphErrors* gr, double *newX)
{
  double x,y,yerr;
  for(int i=0;i<gr->GetN();i++)
  {
    gr->GetPoint(i,x,y);
    yerr=gr->GetErrorY(i);
    gr->SetPoint(i,newX[i],y);
    gr->SetPointError(i,0.0,yerr);
  }
}

void HistogramScaleShift(TH1D* h, double scale=1.0, double shift=0.0)
{
  for(int i=1;i<=h->GetNbinsX();i++)
  {
    h->SetBinContent(i,h->GetBinContent(i)/scale+shift);
    h->SetBinError(i,h->GetBinError(i)/scale);
  }
}

TGraphErrors* GraphDivide(TGraphErrors* gr, TGraphErrors* gr1, int offset = 0, int offset1 = 0)
{
  TGraphErrors* gr_ratio = new TGraphErrors(gr1->GetN());
  gr_ratio->SetMarkerStyle(gr1->GetMarkerStyle());
  gr_ratio->SetMarkerColor(gr1->GetMarkerColor());
  gr_ratio->SetMarkerSize(gr1->GetMarkerSize());
  gr_ratio->SetLineStyle(gr1->GetLineStyle());
  gr_ratio->SetLineColor(gr1->GetLineColor());
  gr_ratio->SetLineWidth(gr1->GetLineWidth());

  double x,y,yerr;
  double x1,y1,yerr1;

  for(int i=0;i<gr1->GetN();i++)
  {
    gr->GetPoint(i+offset,x,y);
    gr1->GetPoint(i+offset1,x1,y1);
    yerr=gr->GetErrorY(i+offset);
    yerr1=gr1->GetErrorY(i+offset1);
    gr_ratio->SetPoint(i,x,fabs(y1/y));
    gr_ratio->SetPointError(i,0.0,fabs(sqrt((yerr/y)**2+(yerr1/y1)**2)*y1/y));
  }

  return gr_ratio;
}

TGraphErrors* GraphMultiply(TGraphErrors* gr, double array[])
{
  const int npoints = gr->GetN();
  TGraphErrors* gr_new = new TGraphErrors(npoints);

  double x,y,yerr;

  for(int i=0;i<npoints;i++)
  {
    gr->GetPoint(i,x,y);
    yerr=gr->GetErrorY(i);
    gr_new->SetPoint(i,x,y*array[i]);
    gr_new->SetPointError(i,0.0,yerr*array[i]);
  }

  return gr_new;
}

