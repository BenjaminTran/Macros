void c1_ks()
{
//=========Macro generated from canvas: c1_ks/Plot_ks
//=========  (Fri Aug 25 03:28:08 2017) by ROOT version6.10/02
   TCanvas *c1_ks = new TCanvas("c1_ks", "Plot_ks",10,142,550,500);
   c1_ks->Range(-2.926829,-0.1874324,21.46341,0.5287838);
   c1_ks->SetFillColor(0);
   c1_ks->SetBorderMode(0);
   c1_ks->SetBorderSize(10);
   c1_ks->SetTickx(1);
   c1_ks->SetTicky(1);
   c1_ks->SetLeftMargin(0.12);
   c1_ks->SetRightMargin(0.06);
   c1_ks->SetBottomMargin(0.15);
   c1_ks->SetFrameFillStyle(0);
   c1_ks->SetFrameLineStyle(0);
   c1_ks->SetFrameBorderMode(0);
   c1_ks->SetFrameBorderSize(10);
   c1_ks->SetFrameFillStyle(0);
   c1_ks->SetFrameLineStyle(0);
   c1_ks->SetFrameBorderMode(0);
   c1_ks->SetFrameBorderSize(10);
   
   TH1F *hframe__1 = new TH1F("hframe__1","K_{S}^{0} Reconstruction Cuts",1000,0,20);
   hframe__1->SetMinimum(-0.08);
   hframe__1->SetMaximum(0.45);
   hframe__1->SetDirectory(0);
   hframe__1->SetStats(0);
   hframe__1->SetLineStyle(0);
   hframe__1->SetLineWidth(3);
   hframe__1->SetMarkerStyle(20);
   hframe__1->SetMarkerSize(1.2);
   hframe__1->GetXaxis()->SetTitle("p_{T} (GeV)");
   hframe__1->GetXaxis()->CenterTitle(true);
   hframe__1->GetXaxis()->SetNdivisions(505);
   hframe__1->GetXaxis()->SetLabelFont(42);
   hframe__1->GetXaxis()->SetLabelOffset(0.006);
   hframe__1->GetXaxis()->SetTitleSize(0.05);
   hframe__1->GetXaxis()->SetTitleOffset(1.2);
   hframe__1->GetXaxis()->SetTitleFont(42);
   hframe__1->GetYaxis()->SetTitle("v_{2}^{sig}");
   hframe__1->GetYaxis()->CenterTitle(true);
   hframe__1->GetYaxis()->SetNdivisions(505);
   hframe__1->GetYaxis()->SetLabelFont(42);
   hframe__1->GetYaxis()->SetLabelOffset(0.006);
   hframe__1->GetYaxis()->SetTitleSize(0.05);
   hframe__1->GetYaxis()->SetTitleOffset(1.1);
   hframe__1->GetYaxis()->SetTitleFont(42);
   hframe__1->GetZaxis()->SetNdivisions(505);
   hframe__1->GetZaxis()->SetLabelFont(42);
   hframe__1->GetZaxis()->SetLabelOffset(0.006);
   hframe__1->GetZaxis()->SetTitleSize(0.05);
   hframe__1->GetZaxis()->SetTitleOffset(1.5);
   hframe__1->GetZaxis()->SetTitleFont(42);
   hframe__1->Draw(" ");
   
   Double_t Graph0_fx1001[16] = {
   0.3666,
   0.5309,
   0.711,
   0.9046,
   1.202,
   1.591,
   1.986,
   2.465,
   3.136,
   4.008,
   5.142,
   6.431,
   7.619,
   9.142,
   11.64,
   16.86};
   Double_t Graph0_fy1001[16] = {
   0.0135543,
   0.0293041,
   0.0433881,
   0.0594498,
   0.082389,
   0.106184,
   0.124363,
   0.137749,
   0.146475,
   0.145182,
   0.135266,
   0.125505,
   0.124708,
   0.11911,
   0.133803,
   0.0390893};
   Double_t Graph0_fex1001[16] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph0_fey1001[16] = {
   0.00217719,
   0.000522063,
   0.000305164,
   0.000248705,
   0.000161989,
   0.000173887,
   0.000208535,
   0.000225897,
   0.000298544,
   0.000448324,
   0.000705022,
   0.00146462,
   0.00192261,
   0.00315165,
   0.0042099,
   0.0308448};
   TGraphErrors *gre = new TGraphErrors(16,Graph0_fx1001,Graph0_fy1001,Graph0_fex1001,Graph0_fey1001);
   gre->SetName("Graph0");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ff0000");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph1001 = new TH1F("Graph_Graph1001","Graph",100,0,18.50934);
   Graph_Graph1001->SetMinimum(0.00742005);
   Graph_Graph1001->SetMaximum(0.1606264);
   Graph_Graph1001->SetDirectory(0);
   Graph_Graph1001->SetStats(0);
   Graph_Graph1001->SetLineStyle(0);
   Graph_Graph1001->SetLineWidth(3);
   Graph_Graph1001->SetMarkerStyle(20);
   Graph_Graph1001->SetMarkerSize(1.2);
   Graph_Graph1001->GetXaxis()->SetNdivisions(505);
   Graph_Graph1001->GetXaxis()->SetLabelFont(42);
   Graph_Graph1001->GetXaxis()->SetLabelOffset(0.006);
   Graph_Graph1001->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph1001->GetXaxis()->SetTitleFont(42);
   Graph_Graph1001->GetYaxis()->SetNdivisions(505);
   Graph_Graph1001->GetYaxis()->SetLabelFont(42);
   Graph_Graph1001->GetYaxis()->SetLabelOffset(0.006);
   Graph_Graph1001->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph1001->GetYaxis()->SetTitleFont(42);
   Graph_Graph1001->GetZaxis()->SetNdivisions(505);
   Graph_Graph1001->GetZaxis()->SetLabelFont(42);
   Graph_Graph1001->GetZaxis()->SetLabelOffset(0.006);
   Graph_Graph1001->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph1001->GetZaxis()->SetTitleOffset(1.5);
   Graph_Graph1001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1001);
   
   gre->Draw("p");
   
   Double_t Graph1_fx1002[16] = {
   0.3666,
   0.5309,
   0.711,
   0.9046,
   1.202,
   1.591,
   1.986,
   2.465,
   3.136,
   4.008,
   5.142,
   6.431,
   7.619,
   9.142,
   11.64,
   16.86};
   Double_t Graph1_fy1002[16] = {
   0.00784946,
   0.0303067,
   0.0438613,
   0.0596843,
   0.0831109,
   0.107091,
   0.124913,
   0.137904,
   0.147068,
   0.146607,
   0.134644,
   0.122733,
   0.126291,
   0.109094,
   0.131269,
   0.077138};
   Double_t Graph1_fex1002[16] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph1_fey1002[16] = {
   0.00646683,
   0.00128322,
   0.000643845,
   0.000482569,
   0.000295253,
   0.000304724,
   0.000358735,
   0.000384387,
   0.000504543,
   0.000755624,
   0.00118735,
   0.00246939,
   0.00324637,
   0.00534846,
   0.00762112,
   0.07372};
   gre = new TGraphErrors(16,Graph1_fx1002,Graph1_fy1002,Graph1_fex1002,Graph1_fey1002);
   gre->SetName("Graph1");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(26);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph1002 = new TH1F("Graph_Graph1002","Graph",100,0,18.50934);
   Graph_Graph1002->SetMinimum(0.001244367);
   Graph_Graph1002->SetMaximum(0.1658055);
   Graph_Graph1002->SetDirectory(0);
   Graph_Graph1002->SetStats(0);
   Graph_Graph1002->SetLineStyle(0);
   Graph_Graph1002->SetLineWidth(3);
   Graph_Graph1002->SetMarkerStyle(20);
   Graph_Graph1002->SetMarkerSize(1.2);
   Graph_Graph1002->GetXaxis()->SetNdivisions(505);
   Graph_Graph1002->GetXaxis()->SetLabelFont(42);
   Graph_Graph1002->GetXaxis()->SetLabelOffset(0.006);
   Graph_Graph1002->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph1002->GetXaxis()->SetTitleFont(42);
   Graph_Graph1002->GetYaxis()->SetNdivisions(505);
   Graph_Graph1002->GetYaxis()->SetLabelFont(42);
   Graph_Graph1002->GetYaxis()->SetLabelOffset(0.006);
   Graph_Graph1002->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph1002->GetYaxis()->SetTitleFont(42);
   Graph_Graph1002->GetZaxis()->SetNdivisions(505);
   Graph_Graph1002->GetZaxis()->SetLabelFont(42);
   Graph_Graph1002->GetZaxis()->SetLabelOffset(0.006);
   Graph_Graph1002->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph1002->GetZaxis()->SetTitleOffset(1.5);
   Graph_Graph1002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1002);
   
   gre->Draw("p");
   
   Double_t Graph2_fx1003[16] = {
   0.3666,
   0.5309,
   0.711,
   0.9046,
   1.202,
   1.591,
   1.986,
   2.465,
   3.136,
   4.008,
   5.142,
   6.431,
   7.619,
   9.142,
   11.64,
   16.86};
   Double_t Graph2_fy1003[16] = {
   0.0103546,
   0.0302139,
   0.0433535,
   0.0596093,
   0.082594,
   0.106931,
   0.124934,
   0.138023,
   0.147092,
   0.146583,
   0.135002,
   0.122919,
   0.127094,
   0.110311,
   0.134414,
   -0.0317351};
   Double_t Graph2_fex1003[16] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t Graph2_fey1003[16] = {
   0.00362787,
   0.000865776,
   0.000504115,
   0.000409552,
   0.000265851,
   0.000284811,
   0.00034122,
   0.000369488,
   0.000488002,
   0.000731881,
   0.00114944,
   0.00238476,
   0.00312644,
   0.00512804,
   0.00682114,
   0.0492771};
   gre = new TGraphErrors(16,Graph2_fx1003,Graph2_fy1003,Graph2_fex1003,Graph2_fey1003);
   gre->SetName("Graph2");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#ff0000");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(25);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph1003 = new TH1F("Graph_Graph1003","Graph",100,0,18.50934);
   Graph_Graph1003->SetMinimum(-0.1038714);
   Graph_Graph1003->SetMaximum(0.1704392);
   Graph_Graph1003->SetDirectory(0);
   Graph_Graph1003->SetStats(0);
   Graph_Graph1003->SetLineStyle(0);
   Graph_Graph1003->SetLineWidth(3);
   Graph_Graph1003->SetMarkerStyle(20);
   Graph_Graph1003->SetMarkerSize(1.2);
   Graph_Graph1003->GetXaxis()->SetNdivisions(505);
   Graph_Graph1003->GetXaxis()->SetLabelFont(42);
   Graph_Graph1003->GetXaxis()->SetLabelOffset(0.006);
   Graph_Graph1003->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph1003->GetXaxis()->SetTitleFont(42);
   Graph_Graph1003->GetYaxis()->SetNdivisions(505);
   Graph_Graph1003->GetYaxis()->SetLabelFont(42);
   Graph_Graph1003->GetYaxis()->SetLabelOffset(0.006);
   Graph_Graph1003->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph1003->GetYaxis()->SetTitleFont(42);
   Graph_Graph1003->GetZaxis()->SetNdivisions(505);
   Graph_Graph1003->GetZaxis()->SetLabelFont(42);
   Graph_Graph1003->GetZaxis()->SetLabelOffset(0.006);
   Graph_Graph1003->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph1003->GetZaxis()->SetTitleOffset(1.5);
   Graph_Graph1003->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1003);
   
   gre->Draw("p");
   
   TLegend *leg = new TLegend(0.15,0.55,0.27,0.75,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.03);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(10);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("Graph0","Standard reconstruction","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","Tight reconstruction","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(26);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","Loose reconstruction","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(25);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   leg->Draw();
   TLatex *   tex = new TLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
tex->SetNDC();
   tex->SetTextSize(0.045);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.4,0.32,"185 #leq N_{trk}^{offline} < 250");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   TLine *line = new TLine(0,0,6,0);
   line->SetLineStyle(2);
   line->Draw();
   
   TPaveText *pt = new TPaveText(0.379562,0.9326316,0.6770073,0.9789474,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   TText *pt_LaTex = pt->AddText("K_{S}^{0} Reconstruction Cuts");
   pt->Draw();
   c1_ks->Modified();
   c1_ks->cd();
   c1_ks->SetSelected(c1_ks);
}
