void fitMass(string fName, bool isLaserCalib, bool isAppend, string prefix=""){
    
    ///extracted from FitEpsilonPlot.cc
    bool useFit_RooMinuit_ = true;
    const int ngaus = 1;
    int niter = 1;

    static double upper_bound_pi0mass_EB = 0.15;
    static double upper_bound_pi0mass_EE = 0.16;
    
    static float fitRange_low_pi0 = 0.080; // value used in the fit function to define the fit range
    //original from Rome's code but had to reduce since the drop in 12X is before 0.212 or 0.222
    //static float fitRange_high_pi0 = 0.212; // value used in the fit function to define the fit range
    //static float fitRange_high_pi0_ext = 0.222;

    static float fitRange_high_pi0 = 0.2; // value used in the fit function to define the fit range
    static float fitRange_high_pi0_ext = 0.2;

    ofstream  outfile;
    if(isAppend){
        
        if(isLaserCalib) outfile.open(Form("pi0_fitMassInfo_withCalib.txt"),std::ofstream::app);
        else outfile.open(Form("pi0_fitMassInfo_withoutCalib.txt"),std::ofstream::app);
    }
    else{
        outfile.open(Form("pi0_fitMassInfo_%s.txt",fName.c_str()));
        outfile << "Year \t Month \t Day \t Time \t Region \t MeanMass \t MeanUnc\n";
    }
    
    
    

    TFile *fin = TFile::Open(Form("%s.root",fName.c_str()));
    TIter next(fin->GetListOfKeys());
    TKey *key;
    vector<string> var_vec;

    while ((key = (TKey*)next())) {
        
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (cl->InheritsFrom("TDirectory")) {
            continue;
        }
        
        if (!cl->InheritsFrom("TH1F")) continue;
        TH1F *htmp = (TH1F*)key->ReadObj();
        var_vec.push_back(htmp->GetName());
    }


    for(int ivar=0; ivar<var_vec.size(); ivar++){

        bool isEB = true;
        string var = var_vec[ivar];
        
        ////get the year, month, day, time and region ///name has to be like this: pi0_mass_2018_8_17_8.716667_EE.png
        size_t pos = 0;
        vector<std::string> token;
        
        string clone_var = var;
        string delimiter = "_";
        while ((pos = clone_var.find(delimiter)) != std::string::npos) {
            token.push_back(clone_var.substr(0, pos));
            clone_var.erase(0, pos + delimiter.length());
        }
        token.push_back(clone_var);

        /////// needed for writing a txt file
        string year = token[2];
        string month = token[3];
        string day = token[4];
        string time = token[5];
        string region = token[6];
        cout<<"Year : month : day : time : region "<<year<<" "<<month<<" "<<day<<" "<<time<<" "<<region<<endl;


        ///Canvas
        TCanvas* canvas = new TCanvas(Form("%s_c",var.c_str()),"",700,700);
        canvas->cd();
        canvas->SetTickx(1);
        canvas->SetTicky(1);
        canvas->cd();
        canvas->SetRightMargin(0.06);
        canvas->SetLeftMargin(0.15);
        
        
        
        TH1F *h = (TH1F*)fin->Get(var.c_str());
        
        int fstr = var.find("EE",0);
        if(fstr!=string::npos)
        {
            isEB = false;
        }

        
        Double_t upMassBoundaryEB = upper_bound_pi0mass_EB; 
        Double_t upMassBoundaryEE = upper_bound_pi0mass_EE; 
        Double_t upMassBoundary = isEB ? upMassBoundaryEB : upMassBoundaryEE;
        // need a patch for some crystals in EB that might have the peak around 160 MeV due to high laser corrections.
        // depending on the year, the containment corrections might also increase a bit the peak position
        // the problem is that the peak is expected to be below 150 MeV when defining the signal model, so we have to catch these exception
        Double_t xValMaxHisto = h->GetXaxis()->GetBinCenter(h->GetMaximumBin()+1); // use value just above maximum, it will be used to set the mean of the gaussian
        // check the maximum is within xlo and xhi (the histogram range is larger than the fit range)
        Double_t maxMassForGaussianMean = 0.0; //upper_bound_pi0mass_EB;
        // first check if peak is in the fit range (for EB it will be, in EE the background rises up and the maximum might not coincide wth peak)
        double xhi = fitRange_high_pi0;
        double xlo = fitRange_low_pi0;
        
        if (xValMaxHisto < xhi) {
            
            if (xValMaxHisto > upMassBoundary) {
                maxMassForGaussianMean = xValMaxHisto;
                xhi = fitRange_high_pi0_ext; //xhi + 0.012; // increase a bit the fit range
            } else {
                maxMassForGaussianMean = upMassBoundary;
            }
            
        } else {
            
            // need to loop on bins in the fit window
            Double_t ymaxHisto = 0.0;
            Int_t binYmaxHisto = -1;
            for (Int_t ibin = h->GetXaxis()->FindFixBin(xlo); ibin <= h->GetXaxis()->FindFixBin(xhi); ibin++) {
                if (h->GetBinContent(ibin) > ymaxHisto) {
                    ymaxHisto = h->GetBinContent(ibin);
                    binYmaxHisto = ibin;
                }
            }
            // check if maximum was found and it was not the last-1 bin
            // in that case, use the next bin to get max value for mass, just to avoid biases (that's why we asked last-1)
            if (binYmaxHisto > 0 && binYmaxHisto < (h->GetXaxis()->FindFixBin(xhi)-1)) {
                maxMassForGaussianMean = h->GetXaxis()->GetBinCenter(binYmaxHisto+1);
                if (maxMassForGaussianMean > upMassBoundary) xhi = fitRange_high_pi0_ext;  //xhi + 0.012; // increase a bit the fit range
                
            } else {
                maxMassForGaussianMean = upMassBoundary; // if all this mess didn't work, just use the value we would have used in the beginning
            }
            
        }
        
        RooRealVar x("x","#gamma#gamma invariant mass",xlo, xhi, "GeV/c^2");
    
        RooDataHist dh("dh","#gamma#gamma invariant mass",RooArgList(x),h);
        
        RooRealVar mean("mean","#pi^{0} peak position", 0.13,  0.105, maxMassForGaussianMean,"GeV/c^{2}");
        RooRealVar sigma("sigma","#pi^{0} core #sigma",0.011, 0.005 ,0.015,"GeV/c^{2}");
        
        
        if(!isEB)  {
            mean.setRange( 0.1, maxMassForGaussianMean);
            mean.setVal(0.13);
            sigma.setRange(0.005,0.020);
        }
        if(isEB){
            mean.setRange(0.105, maxMassForGaussianMean);
            sigma.setRange(0.003,0.030);	  
        }
        
        //RooRealVar Nsig("Nsig","#pi^{0} yield",1000.,0.,1.e7);
        RooRealVar Nsig("Nsig","#pi^{0} yield",h->Integral()*0.15,0.,h->Integral()*10.0);
        //Nsig.setVal( h->GetSum()*0.1);
        
        RooGaussian gaus("gaus","Core Gaussian",x, mean,sigma);
        
        RooRealVar sigmaTail("sigmaTail","#pi^{0} tail #sigma",0.040, 0.020,0.065,"GeV/c^{2}");
        
        RooGaussian gaus2("gaus2","Tail Gaussian",x, mean,sigmaTail);
        
        RooRealVar fcore("fcore","f_{core}",0.9,0.,1.);
        RooAddPdf  signal("signal","signal model",RooArgList(gaus,gaus2),fcore);
        
        RooRealVar p0("p0","p0", 1000.,-1.e5,1.e5);
        RooRealVar p1("p1","p1", -3000.,-1.e5,1.e5);
        RooRealVar p2("p2","p2", 10000.,-1.e5,1.e5);
        RooRealVar p3("p3","p3", -10000.,-1.e5,1.e5);
        RooRealVar p4("p4","p4",-4000.,-1.e5,1.e5);
        RooRealVar p5("p5","p5", 5.,-1.e5,1.e5);
        RooRealVar p6("p6","p6", 6.,-1.e5,1.e5);
        
        RooRealVar cb0("cb0","cb0", 0.2, -1.,1.);
        RooRealVar cb1("cb1","cb1",-0.1, -1.,1.);
        RooRealVar cb2("cb2","cb2", 0.1,  -1.,1.);
        RooRealVar cb3("cb3","cb3",-0.1, -0.5,0.5);
        RooRealVar cb4("cb4","cb4", 0.1, -1.,1.);
        RooRealVar cb5("cb5","cb5", 0.1, -1.,1.);
        RooRealVar cb6("cb6","cb6", 0.3, -1.,1.);


        //RooChebychev bkg("bkg","bkg model", x, RooArgList(cb0,cb1,cb2) );
        //RooChebychev bkg("bkg","bkg model", x, RooArgList(cb0,cb1,cb2,cb3) );
        
        RooArgList cbpars(cb0,cb1,cb2);
        
        // try to use a second order polynomial, if the fit is bad add other terms
        // if you start with many terms, the fit creates strange curvy shapes trying to fit the statistical fluctuations
        // 2nd order means a curve with no change of concavity
        
        if(niter==1){
            cbpars.add(cb3);
        }
        if(niter==2){
            cb3.setRange(-1,1.);
            cb4.setRange(-0.3,0.3);
            cbpars.add( cb3);
            cbpars.add( cb4 );     
        }
        if(niter==3){
            cb3.setRange(-1,1.);
            cb4.setRange(-1,1);
            cb5.setRange(-0.5, 0.5);
            cbpars.add( cb3);
            cbpars.add( cb4 );
            cbpars.add( cb5 );
        }
        
        RooChebychev bkg("bkg","bkg model", x, cbpars );
        
        //RooPolynomial bkg("bkg","background model",x,RooArgList(p0,p1,p2,p3,p4,p5,p6) );
        //RooPolynomial bkg("bkg","background model",x,RooArgList(p0,p1,p2,p3) );
        
        //RooRealVar Nbkg("Nbkg","background yield",1.e3,0.,1.e8);
        RooRealVar Nbkg("Nbkg","background yield",h->Integral()*0.85,0.,h->Integral()*10.0);
        //Nbkg.setVal( h->GetSum()*0.8 );
        
        RooAbsPdf* model=0;
        
        RooAddPdf model1("model","sig+bkg",RooArgList(gaus,bkg),RooArgList(Nsig,Nbkg));
        RooAddPdf model2("model","sig+bkg",RooArgList(signal,bkg),RooArgList(Nsig,Nbkg));
        
        if(ngaus==1)      model = &model1;
        else if(ngaus==2) model = &model2;
        
        
        RooNLLVar nll("nll","log likelihood var",*model,dh, RooFit::Extended(true));
        //RooAbsReal * nll = model->createNLL(dh); //suggetsed way, taht should be the same
        
        RooFitResult* res = nullptr;
        RooMinuit m(nll);
        RooMinimizer mfit(nll);
        
        if (useFit_RooMinuit_) {
            
            // // original fit
            // // obsolete: see here --> https://root-forum.cern.ch/t/roominuit-and-roominimizer-difference/18230/8
            // // better to use RooMinimizer, but please read caveat below
            m.setVerbose(kFALSE);
            //m.setVerbose(kTRUE);
            m.migrad();
            m.hesse();  // sometimes it fails, caution
            res = m.save() ;
            
    } else {
            
            // alternative fit (results are pretty much the same)
      // IMPORTANT, READ CAREFULLY: sometimes this method fails.
      // This happens because at the boundaries of the fit range the pdf goea slightly below 0 (so it is negative). The fitter tries to cope wth it and should tipically
      // manage to converge. However, I noticed that after few attemps (even though the default number of attemps should be several hundreds or thousands of times) 
      // the job crashes, and this seems to be a feature of cmssw, not of RooFit
      // The reason why the pdf gets negative could be due to the fact that, regardless the chosen fit range given by xlo and xhi, the actual fit range goes from the 
      // lower edge of the leftmost bin containing xlo to the upper edge of the rightmost one containing xhi, but then the fit tries to "pass" across the bin centers
      // Therefore, for a sharply rising (or falling) distribution, the pdf can become negative
      // The consequence is that there are large areas in the calibration map of related 2D plots that are white (because the fit there was not done succesfully)
      // The previous method using RooMinuit seems to be more robust, so I suggest we should use that one even though it is said to be obsolete
            mfit.setVerbose(kFALSE);
            mfit.setPrintLevel(-1);
            mfit.setStrategy(2);  // 0,1,2:  MINUIT strategies for dealing most efficiently with fast FCNs (0), expensive FCNs (2) and 'intermediate' FCNs (1)
            //cout << "FIT_EPSILON: Minimize" << endl;
            mfit.minimize("Minuit2","minimize");
            //cout << "FIT_EPSILON: Minimize hesse " << endl;
            mfit.minimize("Minuit2","hesse");
            //cout<<"FIT_EPSILON: Estimate minos errors for all parameters"<<endl;
            mfit.minos(RooArgSet(Nsig,Nbkg,mean));
            res = mfit.save() ;
            
        }

        RooChi2Var chi2("chi2","chi2 var",*model,dh, true);
        // use only bins in fit range for ndof (dh is made with var x that already has the restricted range, but h is the full histogram)
    //int ndof = h->GetNbinsX() - res->floatParsFinal().getSize();
        int ndof = h->FindFixBin(xhi) - h->FindFixBin(xlo) +1 - res->floatParsFinal().getSize(); 

        //compute S/B and chi2
        x.setRange("sobRange",mean.getVal()-3.*sigma.getVal(), mean.getVal()+3.*sigma.getVal());
        RooAbsReal* integralSig = gaus.createIntegral(x,RooFit::NormSet(x),RooFit::Range("sobRange"));
        
        RooAbsReal* integralBkg = bkg.createIntegral(x,RooFit::NormSet(x),RooFit::Range("sobRange"));
        
        float normSig = integralSig->getVal();
        float normBkg = integralBkg->getVal();
        
        
        double S = normSig*Nsig.getVal();
        double Serr = normSig*Nsig.getError();
        
        double B = normBkg*Nbkg.getVal();
        double Berr = normBkg*Nbkg.getError();
        
        double SoB =  S/B;
        double SoBerr =  SoB*sqrt( pow(Serr/S,2) + 
                                          pow(Berr/B,2) ) ;
        double dof = ndof;
        double nFitParam = res->floatParsFinal().getSize();
        
        
        RooPlot*  xframe = x.frame(h->GetNbinsX());
        xframe->SetName(Form("%s_rp",var.c_str()));
        xframe->SetTitle(h->GetTitle());
        dh.plotOn(xframe, RooFit::Name("data"));
        model->plotOn(xframe,RooFit::Components(bkg),RooFit::LineStyle(kDashed), RooFit::LineColor(kRed), RooFit::Name("bkgOnly"));
        model->plotOn(xframe,RooFit::Components(gaus),RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen+1), RooFit::Name("sigOnly"));
        model->plotOn(xframe, RooFit::Name("model"));
        
        // TMAth::Prob() uses Chi2, not reduced Chi2, while xframe->chiSquare() returns the reduced Chi2
        double chi2Fit = xframe->chiSquare("model","data",nFitParam) * dof;
        double probchi2 = TMath::Prob(chi2Fit, ndof);
        
        xframe->Draw();

        cout << "FIT_EPSILON: Nsig: " << Nsig.getVal() 
             << " nsig 3sig: " << normSig*Nsig.getVal()
             << " nbkg 3sig: " << normBkg*Nbkg.getVal()
             << " S/B: " << SoB << " +/- " << SoBerr
             << " chi2: " << chi2Fit
             << " chi2 reduced: " << chi2Fit / dof
             << " DOF: " << dof
             << " N(fit.param.): " << nFitParam
             << " prob(chi2): " << probchi2
             << endl;
        
        TLatex lat;
        std::string line = "";
        lat.SetNDC();
        lat.SetTextSize(0.040);
        lat.SetTextColor(1);
        
        float xmin(0.2), yhi(0.80), ypass(0.05);
        if(!isEB) yhi=0.5;
        line = Form("Yield: %.0f #pm %.0f", Nsig.getVal(), Nsig.getError() );
        lat.DrawLatex(xmin,yhi, line.c_str());
        
        line = Form("m_{#gamma#gamma}: %.2f #pm %.2f", mean.getVal()*1000., mean.getError()*1000. );
        lat.DrawLatex(xmin,yhi-ypass, line.c_str());
        
        line = Form("#sigma: %.2f #pm %.2f (%.2f%s)", sigma.getVal()*1000., sigma.getError()*1000., sigma.getVal()*100./mean.getVal(), "%" );
        lat.DrawLatex(xmin,yhi-2.*ypass, line.c_str());
        
        //sprintf(line,"S/B(3#sigma): %.2f #pm %.2f", SoB, SoBerr );
        line = Form("S/B(3#sigma): %.2f", SoB );
        lat.DrawLatex(xmin,yhi-3.*ypass, line.c_str());
        
        line = Form("#Chi^{2}: %.2f (%d dof)", chi2Fit, dof );
        lat.DrawLatex(xmin,yhi-4.*ypass, line.c_str());
        
        line = Form("B param. %d", cbpars.getSize() );
        lat.DrawLatex(xmin,yhi-5.*ypass, line.c_str());
        
        canvas->RedrawAxis("sameaxis");
        canvas->Print(Form("%s/%s.png", prefix.c_str(), var.c_str()));
        canvas->Print(Form("%s/%s.root", prefix.c_str(), var.c_str()));

        
        outfile << year <<" \t "<< month <<" \t "<< day <<" \t "<< time <<" \t "<< region <<" \t "<< mean.getVal()*1000.  <<" \t "<<  mean.getError()*1000.  <<" \n ";
    }//for(int ivar=0; ivar<var_vec.size(); ivar++)
    
    outfile.close();
}
