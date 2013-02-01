#include "TFile.h"
#include "GBRForest.h"
#include <vector>
#include <string>
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"

void analyzeforest() {
  
  TFile *file = new TFile("fgbrtraintestsigmax200-nosig-100.root","READ");
  
  const std::vector<std::string> *varlist = (std::vector<std::string>*)file->Get("varlist");
  const GBRForest *forest = (GBRForest*)file->Get("gbrtrain");
  
  int numvars = varlist->size();
  int numtrees = forest->Trees().size();
  
  for (int ivar = 0; ivar<numvars; ++ivar) {
    printf("%i: %s\n",ivar, varlist->at(ivar).c_str());
  }
  
  TH2D *hvardist = new TH2D("hvardist","",numvars,-0.5,numvars-0.5,numtrees,-0.5,numtrees-0.5);
  TH1D *htreedistall = new TH1D("htreedistall","",numtrees,-0.5,numtrees-0.5);
  TH1D *hvardistall = new TH1D("hvardistall","",numvars,-0.5,numvars-0.5);
  
  for (int itree = 0; itree<numtrees; ++itree) {
    const GBRTree &tree = forest->Trees().at(itree);
    for (int inode=0; inode<tree.CutIndices().size(); ++inode) {
      if (tree.LeftIndices()[inode]==tree.RightIndices()[inode]) continue;
      hvardist->Fill(tree.CutIndices().at(inode),itree);
      hvardistall->Fill(tree.CutIndices().at(inode));
      htreedistall->Fill(itree);
    }
  }
  
  new TCanvas;
  hvardist->GetXaxis()->SetTitle("Variable");
  hvardist->GetYaxis()->SetTitle("Tree");
  hvardist->Draw("COLZ");
  
  new TCanvas;
  hvardistall->GetXaxis()->SetTitle("Variable");
  hvardistall->GetYaxis()->SetTitle("Number of DT Splits");
  hvardistall->Draw("HIST");
  
  new TCanvas;
  htreedistall->GetXaxis()->SetTitle("Tree");
  htreedistall->GetYaxis()->SetTitle("Number of DT Splits");  
  htreedistall->Draw("HIST");
  
  
}