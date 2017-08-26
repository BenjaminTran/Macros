void c1_la()
{
//=========Macro generated from canvas: c1_la/Plot_la
//=========  (Fri Aug 25 03:28:13 2017) by ROOT version6.10/02
   TCanvas *c1_la = new TCanvas("c1_la", "Plot_la",519,331,550,500);
   c1_la->Range(-2.926829,-0.302027,21.46341,0.7114865);
   c1_la->SetFillColor(0);
   c1_la->SetBorderMode(0);
   c1_la->SetBorderSize(10);
   c1_la->SetTickx(1);
   c1_la->SetTicky(1);
   c1_la->SetLeftMargin(0.12);
   c1_la->SetRightMargin(0.06);
   c1_la->SetBottomMargin(0.15);
   c1_la->SetFrameFillStyle(0);
   c1_la->SetFrameLineStyle(0);
   c1_la->SetFrameBorderMode(0);
   c1_la->SetFrameBorderSize(10);
   c1_la->SetFrameFillStyle(0);
   c1_la->SetFrameLineStyle(0);
   c1_la->SetFrameBorderMode(0);
   c1_la->SetFrameBorderSize(10);
   
   TH1F *hframe__2 = new TH1F("hframe__2","#Lambda Reconstruction Cuts",1000,0,20);
   hframe__2->SetMinimum(-0.15);
   hframe__2->SetMaximum(0.6);
   hframe__2->SetDirectory(0);
   hframe__2->SetStats(0);
   hframe__2->SetLineStyle(0);
   hframe__2->SetLineWidth(3);
   hframe__2->SetMarkerStyle(20);
   hframe__2->SetMarkerSize(1.2);
   hframe__2->GetXaxis()->SetTitle("p_{T} (GeV)");
   hframe__2->GetXaxis()->CenterTitle(true);
   hframe__2->GetXaxis()->SetNdivisions(505);
   hframe__2->GetXaxis()->SetLabelFont(42);
   hframe__2->GetXaxis()->SetLabelOffset(0.006);
   hframe__2->GetXaxis()->SetTitleSize(0.05);
   hframe__2->GetXaxis()->SetTitleOffset(1.2);
   hframe__2->GetXaxis()->SetTitleFont(42);
   hframe__2->GetYaxis()->SetTitle("v_{2}^{sig}");
   hframe__2->GetYaxis()->CenterTitle(true);
   hframe__2->GetYaxis()->SetNdivisions(505);
   hframe__2->GetYaxis()->SetLabelFont(42);
   hframe__2->GetYaxis()->SetLabelOffset(0.006);
   hframe__2->GetYaxis()->SetTitleSize(0.05);
   hframe__2->GetYaxis()->SetTitleOffset(1.1);
   hframe__2->GetYaxis()->SetTitleFont(42);
   hframe__2->GetZaxis()->SetNdivisions(505);
   hframe__2->GetZaxis()->SetLabelFont(42);
   hframe__2->GetZaxis()->SetLabelOffset(0.006);
   hframe__2->GetZaxis()->SetTitleSize(0.05);
   hframe__2->GetZaxis()->SetTitleOffset(1.5);
   hframe__2->GetZaxis()->SetTitleFont(42);
   hframe__2->Draw(" ");
   
   Double_t Graph0_fx1004[13] = {
   0.9252,
   1.224,
   1.603,
   1.995,
   2.485,
   3.156,
   4.007,
   5.116,
   6.414,
   7.588,
   9.117,
   11.56,
   16.89};
   Double_t Graph0_fy1004[13] = {
   0.0345669,
   0.0510091,
   0.0771498,
   0.106513,
   0.138504,
   0.171975,
   0.19385,
   0.202724,
   0.19716,
   0.178541,
   0.168714,
   0.148544,
   0.122694};
   Double_t Graph0_fex1004[13] = {
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
   Double_t Graph0_fey1004[13] = {
   0.00106898,
   0.000437383,
   0.000387893,
   0.000396151,
   0.000360204,
   0.000398977,
   0.000561998,
   0.000940498,
   0.00223511,
   0.0032202,
   0.0063431,
   0.00742138,
   0.0622323};
   TGraphErrors *gre = new TGraphErrors(13,Graph0_fx1004,Graph0_fy1004,Graph0_fex1004,Graph0_fey1004);
   gre->SetName("Graph0");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#3333ff");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#3333ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph1004 = new TH1F("Graph_Graph1004","Graph",100,0,18.48648);
   Graph_Graph1004->SetMinimum(0.01648126);
   Graph_Graph1004->SetMaximum(0.2206812);
   Graph_Graph1004->SetDirectory(0);
   Graph_Graph1004->SetStats(0);
   Graph_Graph1004->SetLineStyle(0);
   Graph_Graph1004->SetLineWidth(3);
   Graph_Graph1004->SetMarkerStyle(20);
   Graph_Graph1004->SetMarkerSize(1.2);
   Graph_Graph1004->GetXaxis()->SetNdivisions(505);
   Graph_Graph1004->GetXaxis()->SetLabelFont(42);
   Graph_Graph1004->GetXaxis()->SetLabelOffset(0.006);
   Graph_Graph1004->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph1004->GetXaxis()->SetTitleFont(42);
   Graph_Graph1004->GetYaxis()->SetNdivisions(505);
   Graph_Graph1004->GetYaxis()->SetLabelFont(42);
   Graph_Graph1004->GetYaxis()->SetLabelOffset(0.006);
   Graph_Graph1004->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph1004->GetYaxis()->SetTitleFont(42);
   Graph_Graph1004->GetZaxis()->SetNdivisions(505);
   Graph_Graph1004->GetZaxis()->SetLabelFont(42);
   Graph_Graph1004->GetZaxis()->SetLabelOffset(0.006);
   Graph_Graph1004->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph1004->GetZaxis()->SetTitleOffset(1.5);
   Graph_Graph1004->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1004);
   
   gre->Draw("p");
   
   Double_t Graph1_fx1005[13] = {
   0.9252,
   1.224,
   1.603,
   1.995,
   2.485,
   3.156,
   4.007,
   5.116,
   6.414,
   7.588,
   9.117,
   11.56,
   16.89};
   Double_t Graph1_fy1005[13] = {
   0.0312339,
   0.0484565,
   0.0762154,
   0.106016,
   0.138673,
   0.17217,
   0.192093,
   0.201094,
   0.193343,
   0.172166,
   0.178001,
   0.150506,
   -0.120121};
   Double_t Graph1_fex1005[13] = {
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
   Double_t Graph1_fey1005[13] = {
   0.00224476,
   0.000827416,
   0.000706016,
   0.000709656,
   0.000637924,
   0.000699588,
   0.000977282,
   0.00162535,
   0.003845,
   0.00555049,
   0.0111427,
   0.0131359,
   0.135443};
   gre = new TGraphErrors(13,Graph1_fx1005,Graph1_fy1005,Graph1_fex1005,Graph1_fey1005);
   gre->SetName("Graph1");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#3333ff");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#3333ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(26);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph1005 = new TH1F("Graph_Graph1005","Graph",100,0,18.48648);
   Graph_Graph1005->SetMinimum(-0.3013923);
   Graph_Graph1005->SetMaximum(0.2485477);
   Graph_Graph1005->SetDirectory(0);
   Graph_Graph1005->SetStats(0);
   Graph_Graph1005->SetLineStyle(0);
   Graph_Graph1005->SetLineWidth(3);
   Graph_Graph1005->SetMarkerStyle(20);
   Graph_Graph1005->SetMarkerSize(1.2);
   Graph_Graph1005->GetXaxis()->SetNdivisions(505);
   Graph_Graph1005->GetXaxis()->SetLabelFont(42);
   Graph_Graph1005->GetXaxis()->SetLabelOffset(0.006);
   Graph_Graph1005->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph1005->GetXaxis()->SetTitleFont(42);
   Graph_Graph1005->GetYaxis()->SetNdivisions(505);
   Graph_Graph1005->GetYaxis()->SetLabelFont(42);
   Graph_Graph1005->GetYaxis()->SetLabelOffset(0.006);
   Graph_Graph1005->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph1005->GetYaxis()->SetTitleFont(42);
   Graph_Graph1005->GetZaxis()->SetNdivisions(505);
   Graph_Graph1005->GetZaxis()->SetLabelFont(42);
   Graph_Graph1005->GetZaxis()->SetLabelOffset(0.006);
   Graph_Graph1005->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph1005->GetZaxis()->SetTitleOffset(1.5);
   Graph_Graph1005->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1005);
   
   gre->Draw("p");
   
   Double_t Graph2_fx1006[13] = {
   0.9252,
   1.224,
   1.603,
   1.995,
   2.485,
   3.156,
   4.007,
   5.116,
   6.414,
   7.588,
   9.117,
   11.56,
   16.89};
   Double_t Graph2_fy1006[13] = {
   0.0364876,
   0.0503007,
   0.0766704,
   0.106931,
   0.138954,
   0.172568,
   0.19213,
   0.201382,
   0.193664,
   0.174807,
   0.178977,
   0.14746,
   -0.0336306};
   Double_t Graph2_fex1006[13] = {
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
   Double_t Graph2_fey1006[13] = {
   0.00177145,
   0.000720517,
   0.000636639,
   0.000648237,
   0.000588564,
   0.0006516,
   0.000916771,
   0.00153283,
   0.00362961,
   0.00522624,
   0.0102896,
   0.0119938,
   0.0987539};
   gre = new TGraphErrors(13,Graph2_fx1006,Graph2_fy1006,Graph2_fex1006,Graph2_fey1006);
   gre->SetName("Graph2");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#3333ff");
   gre->SetLineColor(ci);

   ci = TColor::GetColor("#3333ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(25);
   gre->SetMarkerSize(1.5);
   
   TH1F *Graph_Graph1006 = new TH1F("Graph_Graph1006","Graph",100,0,18.48648);
   Graph_Graph1006->SetMinimum(-0.1659144);
   Graph_Graph1006->SetMaximum(0.2364448);
   Graph_Graph1006->SetDirectory(0);
   Graph_Graph1006->SetStats(0);
   Graph_Graph1006->SetLineStyle(0);
   Graph_Graph1006->SetLineWidth(3);
   Graph_Graph1006->SetMarkerStyle(20);
   Graph_Graph1006->SetMarkerSize(1.2);
   Graph_Graph1006->GetXaxis()->SetNdivisions(505);
   Graph_Graph1006->GetXaxis()->SetLabelFont(42);
   Graph_Graph1006->GetXaxis()->SetLabelOffset(0.006);
   Graph_Graph1006->GetXaxis()->SetTitleOffset(1.4);
   Graph_Graph1006->GetXaxis()->SetTitleFont(42);
   Graph_Graph1006->GetYaxis()->SetNdivisions(505);
   Graph_Graph1006->GetYaxis()->SetLabelFont(42);
   Graph_Graph1006->GetYaxis()->SetLabelOffset(0.006);
   Graph_Graph1006->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph1006->GetYaxis()->SetTitleFont(42);
   Graph_Graph1006->GetZaxis()->SetNdivisions(505);
   Graph_Graph1006->GetZaxis()->SetLabelFont(42);
   Graph_Graph1006->GetZaxis()->SetLabelOffset(0.006);
   Graph_Graph1006->GetZaxis()->SetTitleSize(0.05);
   Graph_Graph1006->GetZaxis()->SetTitleOffset(1.5);
   Graph_Graph1006->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1006);
   
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

   ci = TColor::GetColor("#3333ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","Tight reconstruction","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#3333ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(26);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","Loose reconstruction","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#3333ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(25);
   entry->SetMarkerSize(1.5);
   entry->SetTextFont(42);
   leg->Draw();
   TLatex *   tex = new TLatex(0.15,0.8,"CMS pPb #sqrt{S_{#lower[-0.3]{NN}}} = 8.16 TeV");
tex->SetNDC();
   tex->SetTextSize(0.04);
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
   
   TPaveText *pt = new TPaveText(0.350365,0.9326316,0.6478102,0.9789474,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   TText *pt_LaTex = pt->AddText("#Lambda Reconstruction Cuts");
   pt->Draw();
   c1_la->Modified();
   c1_la->cd();
   c1_la->SetSelected(c1_la);
}
