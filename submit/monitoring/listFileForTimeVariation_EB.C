void listFileForTimeVariation_EB()
{
//=========Macro generated from canvas: cEB/
//=========  (Fri Mar 25 16:06:44 2022) by ROOT version 6.22/09
   TCanvas *cEB = new TCanvas("cEB", "",0,0,700,700);
   gStyle->SetOptStat(0);
   cEB->Range(0,0,1,1);
   cEB->SetFillColor(0);
   cEB->SetBorderMode(0);
   cEB->SetBorderSize(2);
   cEB->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: cEB_1
   TPad *cEB_1 = new TPad("cEB_1", "cEB_1",0.01,0.01,0.49,0.99);
   cEB_1->Draw();
   cEB_1->cd();
   cEB_1->Range(0,0,1,1);
   cEB_1->SetFillColor(0);
   cEB_1->SetBorderMode(0);
   cEB_1->SetBorderSize(2);
   cEB_1->SetFrameBorderMode(0);
   cEB_1->Modified();
   cEB->cd();
  
// ------------>Primitives in pad: cEB_2
   TPad *cEB_2 = new TPad("cEB_2", "cEB_2",0.51,0.01,0.99,0.99);
   cEB_2->Draw();
   cEB_2->cd();
   cEB_2->Range(0,0,1,1);
   cEB_2->SetFillColor(0);
   cEB_2->SetBorderMode(0);
   cEB_2->SetBorderSize(2);
   cEB_2->SetFrameBorderMode(0);
   cEB_2->Modified();
   cEB->cd();
  
// ------------>Primitives in pad: pad1
   TPad *pad1 = new TPad("pad1", "Mean VS time",0,0,0.7,1);
   pad1->Draw();
   pad1->cd();
   pad1->Range(8.708334,0.76875,8.758333,1.08125);
   pad1->SetFillColor(0);
   pad1->SetBorderMode(0);
   pad1->SetBorderSize(2);
   pad1->SetFrameBorderMode(0);
   pad1->SetFrameBorderMode(0);
   
   Double_t _fx1001[2] = {
   8.716667,
   8.75};
   Double_t _fy1001[2] = {
   0.8994643,
   0.8983679};
   Double_t _fex1001[2] = {
   0,
   0};
   Double_t _fey1001[2] = {
   0.0008957378,
   0.0009030353};
   TGraphErrors *gre = new TGraphErrors(2,_fx1001,_fy1001,_fex1001,_fey1001);
   gre->SetName("");
   gre->SetTitle("");
   gre->SetFillStyle(1000);
   gre->SetLineColor(8);
   gre->SetMarkerColor(8);
   
   TH1F *Graph_1001 = new TH1F("Graph_1001","",100,8.713334,8.753333);
   Graph_1001->SetMinimum(0.8);
   Graph_1001->SetMaximum(1.05);
   Graph_1001->SetDirectory(0);
   Graph_1001->SetStats(0);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   Graph_1001->SetLineColor(ci);
   Graph_1001->GetXaxis()->SetTitle("Time");
   Graph_1001->GetXaxis()->SetLabelFont(42);
   Graph_1001->GetXaxis()->SetTitleOffset(1);
   Graph_1001->GetXaxis()->SetTitleFont(42);
   Graph_1001->GetYaxis()->SetTitle("Normalized #pi^{0} mass");
   Graph_1001->GetYaxis()->SetLabelFont(42);
   Graph_1001->GetYaxis()->SetTitleFont(42);
   Graph_1001->GetZaxis()->SetLabelFont(42);
   Graph_1001->GetZaxis()->SetTitleOffset(1);
   Graph_1001->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_1001);
   
   gre->Draw("ap");
   TLatex *   tex = new TLatex(8.716667,0.8729167,"With Light Monitoring Correction");
   tex->SetTextColor(8);
   tex->SetTextSize(0.025);
   tex->SetLineWidth(2);
   tex->Draw();
   
   Double_t _fx1002[2] = {
   8.716667,
   8.75};
   Double_t _fy1002[2] = {
   0.8994643,
   0.8983679};
   Double_t _fex1002[2] = {
   0,
   0};
   Double_t _fey1002[2] = {
   0.0008957378,
   0.0009030353};
   gre = new TGraphErrors(2,_fx1002,_fy1002,_fex1002,_fey1002);
   gre->SetName("");
   gre->SetTitle("");
   gre->SetFillStyle(1000);
   gre->SetLineColor(2);
   gre->SetMarkerColor(2);
   
   TH1F *Graph_1002 = new TH1F("Graph_1002","",100,8.713334,8.753333);
   Graph_1002->SetMinimum(0.8);
   Graph_1002->SetMaximum(1.05);
   Graph_1002->SetDirectory(0);
   Graph_1002->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_1002->SetLineColor(ci);
   Graph_1002->GetXaxis()->SetTitle("Time");
   Graph_1002->GetXaxis()->SetLabelFont(42);
   Graph_1002->GetXaxis()->SetTitleOffset(1);
   Graph_1002->GetXaxis()->SetTitleFont(42);
   Graph_1002->GetYaxis()->SetTitle("Normalized #pi^{0} mass");
   Graph_1002->GetYaxis()->SetLabelFont(42);
   Graph_1002->GetYaxis()->SetTitleFont(42);
   Graph_1002->GetZaxis()->SetLabelFont(42);
   Graph_1002->GetZaxis()->SetTitleOffset(1);
   Graph_1002->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_1002);
   
   gre->Draw("p");
      tex = new TLatex(8.716667,0.8729167,"Without Light Monitoring Correction");
   tex->SetTextColor(2);
   tex->SetTextSize(0.025);
   tex->SetLineWidth(2);
   tex->Draw();
   pad1->Modified();
   cEB->cd();
  
// ------------>Primitives in pad: pad2
   TPad *pad2 = new TPad("pad2", "Projection",0.65,0,1,1);
   pad2->Draw();
   pad2->cd();
   pad2->Range(-0.2625,0.7654167,2.3625,1.08625);
   pad2->SetFillColor(0);
   pad2->SetBorderMode(0);
   pad2->SetBorderSize(2);
   pad2->SetFrameBorderMode(0);
   pad2->SetFrameBorderMode(0);
   
   TH1F *hEB__1 = new TH1F("hEB__1","",120,0,1.1);
   hEB__1->SetBinContent(99,2);
   hEB__1->SetEntries(2);
   hEB__1->SetDirectory(0);
   hEB__1->SetStats(0);
   hEB__1->SetFillColor(8);

   ci = TColor::GetColor("#000099");
   hEB__1->SetLineColor(ci);
   hEB__1->GetXaxis()->SetRange(88,115);
   hEB__1->GetXaxis()->SetLabelFont(42);
   hEB__1->GetXaxis()->SetTitleOffset(1);
   hEB__1->GetXaxis()->SetTitleFont(42);
   hEB__1->GetYaxis()->SetLabelFont(42);
   hEB__1->GetYaxis()->SetTitleFont(42);
   hEB__1->GetZaxis()->SetLabelFont(42);
   hEB__1->GetZaxis()->SetTitleOffset(1);
   hEB__1->GetZaxis()->SetTitleFont(42);
   hEB__1->Draw("HBAR");
      tex = new TLatex(0.2,0.8929167,"Mean = 0.90");
   tex->SetTextColor(8);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.8879167,"RMS = 0.00");
   tex->SetTextColor(8);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TH1F *hEB__2 = new TH1F("hEB__2","",120,0,1.1);
   hEB__2->SetBinContent(99,2);
   hEB__2->SetEntries(2);
   hEB__2->SetStats(0);
   hEB__2->SetFillColor(2);

   ci = TColor::GetColor("#000099");
   hEB__2->SetLineColor(ci);
   hEB__2->GetXaxis()->SetRange(88,115);
   hEB__2->GetXaxis()->SetLabelFont(42);
   hEB__2->GetXaxis()->SetTitleOffset(1);
   hEB__2->GetXaxis()->SetTitleFont(42);
   hEB__2->GetYaxis()->SetLabelFont(42);
   hEB__2->GetYaxis()->SetTitleFont(42);
   hEB__2->GetZaxis()->SetLabelFont(42);
   hEB__2->GetZaxis()->SetTitleOffset(1);
   hEB__2->GetZaxis()->SetTitleFont(42);
   hEB__2->Draw("HBARsame");
      tex = new TLatex(0.2,0.8929167,"Mean = 0.90");
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.2,0.8879167,"RMS = 0.00");
   tex->SetTextColor(2);
   tex->SetLineWidth(2);
   tex->Draw();
   pad2->Modified();
   cEB->cd();
   cEB->Modified();
   cEB->cd();
   cEB->SetSelected(cEB);
}
