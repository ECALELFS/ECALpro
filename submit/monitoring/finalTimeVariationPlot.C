#include <stdlib.h>     /* atoi */
#include <algorithm>
#include <functional> 

void drawCanvas(int ifile, TCanvas *canv, vector<double>& meanMass, vector<double>& time, vector<double>& xErr, vector<double>& meanUnc, TH1F *h, int color, string label, bool normalize=true, vector<string> iovs={}, vector<string> iov_times={}){
    
    // normalize to first point
    if(normalize)
    {
        // make sure we use the first point in time
        auto ref = meanMass[std::distance(time.begin(), std::min_element(time.begin(), time.end()))];
        // scale all mass values by the first point        
        for(auto& m : meanMass)
        {
            m /= ref;
            h->Fill(m);  
        }
    }
    else
    {
        for(auto& m : meanMass)
            h->Fill(m);  
    }

    canv->cd();
    TGraphErrors *g = new TGraphErrors(meanMass.size(), &time[0], &meanMass[0], &xErr[0], &meanUnc[0]);
        
    g->SetMaximum(1.05);
    g->SetMinimum(0.80);
    g->SetMarkerSize(0.7);
    g->SetMarkerColor(color+1);
    g->SetLineColor(color+1);
    g->GetXaxis()->SetTimeOffset(0);
    g->GetXaxis()->SetTimeDisplay(1);
    g->GetXaxis()->SetTimeFormat("#splitline{%d/%m}{%H:%M}");
    g->GetXaxis()->SetTitle("Time(day/month-h:m)");
    g->GetXaxis()->SetTitleOffset(1.5);
    g->GetXaxis()->SetLabelOffset(0.02);
    g->GetXaxis()->SetNdivisions(507);
    g->GetYaxis()->SetTitle("Normalized #pi^{0} mass");
    g->SetName("");
    g->SetTitle("");
    if(ifile==0) g->Draw("AP");    
    else g->Draw("Psame");

    double mean = h->GetMean();
    double rms = h->GetRMS();
     
    // rescale summary histogram to be displayed in the same pad as the graph   
    h->Rebin(5);
    h->Scale(1./h->GetMaximum());
    auto xwidth = g->GetXaxis()->GetXmax()-g->GetXaxis()->GetXmin();
    for(int i=1; i<=h->GetNbinsX(); ++i)
        h->SetBinContent(i, h->GetBinContent(i)*(xwidth-1e5)+g->GetXaxis()->GetXmin());

    h->SetFillStyle(3015);
    h->SetLineWidth(0);
    h->SetLineColor(color);
    h->SetFillColor(color);

    h->Draw("HBARsame");

    TLatex *tex = new TLatex(0.1*xwidth+g->GetXaxis()->GetXmin(),0.88,Form("Mean = %0.2f",mean));
    tex->SetTextColor(color);
    tex->Draw();
    canv->Modified();
    canv->Update();

    tex = new TLatex(0.1*xwidth+g->GetXaxis()->GetXmin(),0.86,Form("RMS = %0.2f",rms));
    tex->SetTextColor(color);
    tex->Draw();
    canv->SetGrid();
    canv->Modified();
    canv->Update();

    tex = new TLatex(0.1*xwidth+g->GetXaxis()->GetXmin(),0.9,Form("%s",label.c_str()));
    tex->SetTextSize(0.025);
    tex->SetTextColor(color);
    tex->Draw();
    canv->SetGrid();
    canv->Modified();
    canv->Update();

    TLine liov;
    liov.SetLineStyle(7);
    liov.SetLineWidth(2);
    TText tiov;
    tiov.SetTextAngle(90);
    tiov.SetTextSize(0.03);
    for(unsigned int i=0; i<iovs.size(); ++i)
    {
        TDatime dt(iov_times[i].c_str());
        liov.DrawLine(dt.Convert(), 0.8, dt.Convert(), 1.05);        
        tiov.DrawText(dt.Convert(), 0.8, ("PS-"+iovs[i]).c_str());
    }        

    //c->Print(Form("%s.png",plotName.c_str()));
}


void finalTimeVariationPlot(string fName, string prefix="", bool usePDGmass=false, vector<string> iovs={}, vector<string> iov_times={}){

    double pdg_pi0Mass = 134.9770;
    
///read in the input text files
    vector<string> inputTextFiles, labels;
    vector<int> color;
    ifstream infile;
    infile.open(fName.c_str()); 
    string line;
    
    if(!infile.is_open()){
        cout<<"Error!!! Could not open file "<<fName.c_str()<<endl;
    }
    
    while(getline(infile, line))
    {
        std::stringstream ss(line);
        
        int fstr = line.find("#",0);
        if(fstr!=string::npos){
            continue;
        }
        
        std::string tmp_fName, tmp_label, tmp_color;
        std::getline(ss,tmp_fName,',');    std::cout<<tmp_fName<<endl;
        std::getline(ss,tmp_label,','); std::cout<<tmp_label<<endl;
        std::getline(ss,tmp_color,','); std::cout<<tmp_color<<endl;
        inputTextFiles.push_back(tmp_fName);
        labels.push_back(tmp_label);
        color.push_back(atoi(tmp_color.c_str()));
    }

    ///to be plotted on the same canva
    gStyle->SetOptStat(0);
    TCanvas *cEB = new TCanvas("cEB","",700,700);
    TCanvas *cEE = new TCanvas("cEE","",700,700);

    for(int ifile=0; ifile<inputTextFiles.size(); ifile++){
        vector<double> timeEB, timeEE;
        vector<string> regionEB, regionEE;
        vector<double> meanMassEB, meanUncEB, meanMassEE, meanUncEE;
        vector<double> xErrEB, xErrEE;
        ifstream infile;
        
        TH1F *hEB = new TH1F("hEB","",500,0.8,1.05);
        hEB->SetLineColor(color[ifile]);
        
        TH1F *hEE = new TH1F("hEE","",500,0.8,1.05);
        hEE->SetLineColor(color[ifile]);

        infile.open(Form("%s",inputTextFiles[ifile].c_str())); 
        string line;
        
        if(!infile.is_open()){
            cout<<"Error!!! Could not open file "<<inputTextFiles[ifile].c_str()<<endl;
        }
        
        if(infile.is_open()){
            while ( getline (infile,line) ){                
                double tmp_time;
                string tmp_region;
                double tmp_meanMass, tmp_meanUnc;
                infile >> tmp_time >> tmp_region >> tmp_meanMass >> tmp_meanUnc;
                                
                tmp_meanMass = tmp_meanMass/pdg_pi0Mass;
                tmp_meanUnc = tmp_meanUnc/pdg_pi0Mass;
                
                int fstr = line.find("#",0);
                if(fstr!=string::npos){
                    continue;
                }

                if(tmp_region.compare("EB")==0){
                    timeEB.push_back(tmp_time); 
                    regionEB.push_back(tmp_region); 
                    meanMassEB.push_back(tmp_meanMass);
                    meanUncEB.push_back(tmp_meanUnc);
                    xErrEB.push_back(0);
                }
                
                
                if(tmp_region.compare("EE")==0){
                    timeEE.push_back(tmp_time); 
                    regionEE.push_back(tmp_region); 
                    meanMassEE.push_back(tmp_meanMass);
                    meanUncEE.push_back(tmp_meanUnc);
                    xErrEE.push_back(0);
                }
                
                
            }//while ( getline (infile,line) )
            infile.close();
        }//if(infile.is_open())
        
        if(meanMassEB.size()>0){
            drawCanvas(ifile, cEB, meanMassEB, timeEB, xErrEB, meanUncEB, hEB, color[ifile], labels[ifile], !usePDGmass, iovs, iov_times);
        }

        if(meanMassEE.size()>0){
            drawCanvas(ifile,  cEE, meanMassEE, timeEE, xErrEE, meanUncEE, hEE, color[ifile], labels[ifile], !usePDGmass, iovs, iov_times);
        }


    }//for(int ifile=0; ifile<inputTextFiles.size(); ifile++)
    
    cEB->Print(Form("%spi0stability_EB.png",prefix.c_str()));
    cEE->Print(Form("%spi0stability_EE.png",prefix.c_str()));
    cEB->Print(Form("%spi0stability_EB.root",prefix.c_str()));
    cEE->Print(Form("%spi0stability_EE.root",prefix.c_str()));

}


