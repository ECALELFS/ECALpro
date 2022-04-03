void listFileForTimeVariation_EE()
{
//=========Macro generated from canvas: cEE/
//=========  (Fri Mar 25 16:06:44 2022) by ROOT version 6.22/09
   TCanvas *cEE = new TCanvas("cEE", "",0,0,700,700);
   gStyle->SetOptStat(0);
   cEE->Range(0,0,1,1);
   cEE->SetFillColor(0);
   cEE->SetBorderMode(0);
   cEE->SetBorderSize(2);
   cEE->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: cEE_1
   TPad *cEE_1 = new TPad("cEE_1", "cEE_1",0.01,0.01,0.49,0.99);
   cEE_1->Draw();
   cEE_1->cd();
   cEE_1->Range(0,0,1,1);
   cEE_1->SetFillColor(0);
   cEE_1->SetBorderMode(0);
   cEE_1->SetBorderSize(2);
   cEE_1->SetFrameBorderMode(0);
   cEE_1->Modified();
   cEE->cd();
  
// ------------>Primitives in pad: cEE_2
   TPad *cEE_2 = new TPad("cEE_2", "cEE_2",0.51,0.01,0.99,0.99);
   cEE_2->Draw();
   cEE_2->cd();
   cEE_2->Range(0,0,1,1);
   cEE_2->SetFillColor(0);
   cEE_2->SetBorderMode(0);
   cEE_2->SetBorderSize(2);
   cEE_2->SetFrameBorderMode(0);
   cEE_2->Modified();
   cEE->cd();
  
// ------------>Primitives in pad: pad1
   TPad *pad1 = new TPad("pad1", "Mean VS time",0,0,0.7,1);
   pad1->Draw();
   pad1->cd();
   pad1->Range(8.466667,0.76875,9.966667,1.08125);
   pad1->SetFillColor(0);
   pad1->SetBorderMode(0);
   pad1->SetBorderSize(2);
   pad1->SetFrameBorderMode(0);
   pad1->SetFrameBorderMode(0);
   
   Double_t _fx1003[1] = {
   8.716667};
   Double_t _fy1003[1] = {
   0.9069175};
   Double_t _fex1003[1] = {
   0};
   Double_t _fey1003[1] = {
   0.009233573};
   TGraphErrors *gre = new TGraphErrors(1,_fx1003,_fy1003,_fex1003,_fey1003);
   gre->SetName("");
   gre->SetTitle("");
   gre->SetFillStyle(1000);
   gre->SetLineColor(8);
   gre->SetMarkerColor(8);
   
   TH1F *Graph_1003 = new TH1F("Graph_1003","",100,8.616667,9.816667);
   Graph_1003->SetMinimum(0.8);
   Graph_1003->SetMaximum(1.05);
   Graph_1003->SetDirectory(0);
   Graph_1003->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_1003->SetLineColor(ci);
   Graph_1003->GetXaxis()->SetTitle("Time");
   Graph_1003->GetXaxis()->SetLabelFont(42);
   Graph_1003->GetXaxis()->SetTitleOffset(1);
   Graph_1003->GetXaxis()->SetTitleFont(42);
   Graph_1003->GetYaxis()->SetTitle("Normalized #pi^{0} mass");
   Graph_1003->GetYaxis()->SetLabelFont(42);
   Graph_1003->GetYaxis()->SetTitleFont(42);
   Graph_1003->GetZaxis()->SetLabelFont(42);
   Graph_1003->GetZaxis()->SetTitleOffset(1);
   Graph_1003->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_1003);
   
   gre->Draw("ap");
   TLatex *   tex = new TLatex(8.716667,0.8729167,"With Light Monitoring Correction");
   tex->SetTextColor(8);
   tex->SetTextSize(0.025);
   tex->SetLineWidth(2);
   tex->Draw();
   
   Double_t _fx1004[1] = {
   8.716667};
   Double_t _fy1004[1] = {
   0.9069175};
   Double_t _fex1004[1] = {
   0};
   Double_t _fey1004[1] = {
   0.009233573};
   gre = new TGraphErrors(1,_fx1004,_fy1004,_fex1004,_fey1004);
   gre->SetName("");
   gre->SetTitle("");
   gre->SetFillStyle(1000);
   gre->SetLineColor(2);
   gre->SetMarkerColor(2);
   
   TH1F *Graph_1004 = new TH1F("Graph_1004","",100,8.616667,9.816667);
   Graph_1004->SetMinimum(0.8);
   Graph_1004->SetMaximum(1.05);
   Graph_1004->SetDirectory(0);
   Graph_1004->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_1004->SetLineColor(ci);
   Graph_1004->GetXaxis()->SetTitle("Time");
   Graph_1004->GetXaxis()->SetLabelFont(42);
   Graph_1004->GetXaxis()->SetTitleOffset(1);
   Graph_1004->GetXaxis()->SetTitleFont(42);
   Graph_1004->GetYaxis()->SetTitle("Normalized #pi^{0} mass");
   Graph_1004->GetYaxis()->SetLabelFont(42);
   Graph_1004->GetYaxis()->SetTitleFont(42);
   Graph_1004->GetZaxis()->SetLabelFont(42);
   Graph_1004->GetZaxis()->SetTitleOffset(1);
   Graph_1004->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_1004);
   
   gre->Draw("p");
      tex = new TLatex(8.716667,0.8729167,"Without Light Monitoring Correction");
   tex->SetTextColor(2);
   tex->SetTextSize(0.025);
   tex->SetLineWidth(2);
   tex->Draw();
   pad1->Modified();
   cEE->cd();
  
// ------------>Primitives in pad: pad2
   TPad *pad2 = new TPad("pad2", "Projection",0.65,0,1,1);
   pad2->Draw();
   pad2->cd();
   pad2->Range(-0.13125,0.7654167,1.18125,1.08625);
   pad2->SetFillColor(0);
   pad2->SetBorderMode(0);
   pad2->SetBorderSize(2);
   pad2->SetFrameBorderMode(0);
   pad2->SetFrameBorderMode(0);
   
   TH1F *hEE__3 = new TH1F("hEE__3","",120,0,1.1);
   hEE__3->SetBinContent(99,1);
   hEE__3->SetEntries(1);
   hEE__3->SetDirectory(0);
   hEE__3->SetStats(0);
   hEE__3->SetFillColor(8);

   ci = TColor::GetColor("#000099");
   hEE__3->SetLineColor(ci);
   hEE__3->GetXaxis()->SetRange(88,115);
   hEE__3->GetXaxis()->SetLabelFont(42);
   hEE__3->GetXaxis()->SetTitleOffset(1);
   hEE__3->GetXaxis()->SetTitleFont(42);
   hEE__3->GetYaxis()->SetLabelFont(42);
   hEE__3->GetYaxis()->SetTitleFont(42);
   hEE__3->GetZaxis()->SetLabelFont(42);
   hEE__3->GetZaxis()->SetTitleOffset(1);
   hEE__3->GetZaxis()->SetTitleFont(42);
   hEE__3->Draw("HBAR");
      tex = new TLatex(0.2,0.8929167,"Mean = 0.90");
   tex->SetTextColor(8);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.8879167,"RMS = 0.00");
   tex->SetTextColor(8);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TH1F *hEE__4 = new TH1F("hEE__4","",120,0,1.1);
   hEE__4->SetBinContent(99,1);
   hEE__4->SetEntries(1);
   hEE__4->SetStats(0);
   hEE__4->SetFillColor(2);

   ci = TColor::GetColor("#000099");
   hEE__4->SetLineColor(ci);
   hEE__4->GetXaxis()->SetRange(88,115);
   hEE__4->GetXaxis()->SetLabelFont(42);
   hEE__4->GetXaxis()->SetTitleOffset(1);
   hEE__4->GetXaxis()->SetTitleFont(42);
   hEE__4->GetYaxis()->SetLabelFont(42);
   hEE__4->GetYaxis()->SetTitleFont(42);
   hEE__4->GetZaxis()->SetLabelFont(42);
   hEE__4->GetZaxis()->SetTitleOffset(1);
   hEE__4->GetZaxis()->SetTitleFont(42);
   hEE__4->Draw("HBARsame");
      tex = new TLatex(0.2,0.8929167,"Mean = 0.90");
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.8879167,"RMS = 0.00");
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw();
   pad2->Modified();
   cEE->cd();
   cEE->Modified();
   cEE->cd();
   cEE->SetSelected(cEE);
}
