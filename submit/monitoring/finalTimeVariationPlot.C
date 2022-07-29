#include <stdlib.h>     /* atoi */
#include <algorithm>
#include <functional> 

void drawCanvas(int ifile, TPad *p1, TPad *p2, vector<float>& meanMass, vector<float>& time, vector<float>& xErr, vector<float>& meanUnc, TH1F *h, int color, string label, bool normalize=true){
    
    // normalize to first point
    if(normalize)
    {
        // make sure we use the first point in time
        auto ref = meanMass[std::distance(time.begin(), std::min_element(time.begin(), time.end()))];
        // scale all mass values by the first point        
        for(auto& m : meanMass)
            m /= ref;
    }

    p1->cd();
    TGraphErrors *g = new TGraphErrors(meanMass.size(), &time[0], &meanMass[0], &xErr[0], &meanUnc[0]);

    g->SetMaximum(1.05);
    g->SetMinimum(0.80);
    g->SetMarkerSize(0.7);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->GetXaxis()->SetTitle("Time");
    g->GetYaxis()->SetTitle("Normalized #pi^{0} mass");
    g->SetName("");
    g->SetTitle("");
    if(ifile==0) g->Draw("AP");    
    else g->Draw("Psame");
        
    p2->cd();
    h->GetXaxis()->SetRangeUser(0.80,1.05);
    //h->SetMinimum(0.85);
    if(ifile==0) h->Draw("HBAR");
    else h->Draw("HBARsame");

    double mean = h->GetMean();
    double rms = h->GetRMS();
        
    TLatex *tex = new TLatex(0.2,mean-0.1,Form("Mean = %0.2f",mean));
    tex->SetTextColor(color);
    tex->Draw();
    p2->Modified();
    p2->Update();

    tex = new TLatex(0.2,mean-0.15,Form("RMS = %0.2f",rms));
    tex->SetTextColor(color);
    tex->Draw();
    p2->SetGrid();
    p2->Modified();
    p2->Update();

    ///on pad1
    p1->cd();
    tex = new TLatex(time[0],mean-0.03,Form("%s",label.c_str()));
    tex->SetTextSize(0.025);
    tex->SetTextColor(color);
    tex->Draw();
    p1->SetGrid();
    p1->Modified();
    p1->Update();

    //c->Print(Form("%s.png",plotName.c_str()));
}


void finalTimeVariationPlot(string fName, string prefix="", bool usePDGmass=false){

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
    cEB->Divide(2);
    TPad *pad1EB = new TPad("pad1", "Mean VS time",0.0,0.0,0.7,1.0,21);
    TPad *pad2EB = new TPad("pad2", "Projection",0.65,0.0,1.0,1.0,22);

    pad1EB->SetFillColor(0);
    pad2EB->SetFillColor(0);


    pad1EB->Draw();
    pad2EB->Draw();

    ///EE
    TCanvas *cEE = new TCanvas("cEE","",700,700);
    cEE->Divide(2);
    TPad *pad1EE = new TPad("pad1", "Mean VS time",0.0,0.0,0.7,1.0,21);
    TPad *pad2EE = new TPad("pad2", "Projection",0.65,0.0,1.0,1.0,22);

    pad1EE->SetFillColor(0);
    pad2EE->SetFillColor(0);


    pad1EE->Draw();
    pad2EE->Draw();



    for(int ifile=0; ifile<inputTextFiles.size(); ifile++){
        vector<float> yearEB, monthEB, dayEB, timeEB;
        vector<float> yearEE, monthEE, dayEE, timeEE;
        vector<string> regionEB, regionEE;
        vector<float> meanMassEB, meanUncEB, meanMassEE, meanUncEE;
        vector<float> xErrEB, xErrEE;
        ifstream infile;
        
        TH1F *hEB = new TH1F("hEB","",120,0,1.1);
        hEB->SetFillColor(color[ifile]);
        
        TH1F *hEE = new TH1F("hEE","",120,0,1.1);
        hEE->SetFillColor(color[ifile]);

        infile.open(Form("%s",inputTextFiles[ifile].c_str())); 
        string line;
        
        if(!infile.is_open()){
            cout<<"Error!!! Could not open file "<<inputTextFiles[ifile].c_str()<<endl;
        }
        
        if(infile.is_open()){
            while ( getline (infile,line) ){
                float tmp_year, tmp_month, tmp_day, tmp_time;
                string tmp_region;
                float tmp_meanMass, tmp_meanUnc;
                infile >> tmp_year >> tmp_month >> tmp_day >> tmp_time >> tmp_region >> tmp_meanMass >> tmp_meanUnc;
                                
                tmp_meanMass = tmp_meanMass/pdg_pi0Mass;
                tmp_meanUnc = tmp_meanUnc/pdg_pi0Mass;
                
                int fstr = line.find("#",0);
                if(fstr!=string::npos){
                    continue;
                }

                if(tmp_region.compare("EB")==0){
                    yearEB.push_back(tmp_year); 
                    monthEB.push_back(tmp_month); 
                    dayEB.push_back(tmp_day); 
                    timeEB.push_back(tmp_time); 
                    regionEB.push_back(tmp_region); 
                    meanMassEB.push_back(tmp_meanMass);
                    meanUncEB.push_back(tmp_meanUnc);
                    xErrEB.push_back(0);
                    hEB->Fill(tmp_meanMass);
                }
                
                
                if(tmp_region.compare("EE")==0){
                    yearEE.push_back(tmp_year); 
                    monthEE.push_back(tmp_month); 
                    dayEE.push_back(tmp_day); 
                    timeEE.push_back(tmp_time); 
                    regionEE.push_back(tmp_region); 
                    meanMassEE.push_back(tmp_meanMass);
                    meanUncEE.push_back(tmp_meanUnc);
                    xErrEE.push_back(0);
                    hEE->Fill(tmp_meanMass);
                }
                
                
            }//while ( getline (infile,line) )
            infile.close();
        }//if(infile.is_open())
        
        if(meanMassEB.size()>0){
            drawCanvas(ifile, pad1EB, pad2EB, meanMassEB, timeEB, xErrEB, meanUncEB, hEB, color[ifile], labels[ifile], !usePDGmass);
        }

        if(meanMassEE.size()>0){
            drawCanvas(ifile,  pad1EE, pad2EE, meanMassEE, timeEE, xErrEE, meanUncEE, hEE, color[ifile], labels[ifile], !usePDGmass);
        }


    }//for(int ifile=0; ifile<inputTextFiles.size(); ifile++)
    
    cEB->Print(Form("%spi0stability_EB.png",prefix.c_str()));
    cEE->Print(Form("%spi0stability_EE.png",prefix.c_str()));
    cEB->Print(Form("%spi0stability_EB.root",prefix.c_str()));
    cEE->Print(Form("%spi0stability_EE.root",prefix.c_str()));

}


