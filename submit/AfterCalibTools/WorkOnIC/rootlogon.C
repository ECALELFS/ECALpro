// This is the file rootlogon.C
{
   TStyle *myStyle  = new TStyle("MyStyle","My Root Styles");

   // from ROOT plain style
   myStyle->SetCanvasBorderMode(0);
   myStyle->SetPadBorderMode(0);
   myStyle->SetPadColor(0);
   myStyle->SetCanvasColor(0);
   myStyle->SetTitleColor(1);
   myStyle->SetStatColor(0);

   myStyle->SetLabelSize(0.03,"xyz"); // size of axis values

   // default canvas positioning
   myStyle->SetCanvasDefX(900);
   myStyle->SetCanvasDefY(20);
   myStyle->SetCanvasDefH(550);
   myStyle->SetCanvasDefW(540);

   myStyle->SetPadBottomMargin(0.15); // It was 0.1
   myStyle->SetPadTopMargin(0.15);
   myStyle->SetPadLeftMargin(0.15);
   myStyle->SetPadRightMargin(0.15);

   myStyle->SetPadTickX(1);
   myStyle->SetPadTickY(1);

   myStyle->SetFrameBorderMode(0);

   // US letter
   myStyle->SetPaperSize(20, 24);
   myStyle->SetOptStat(0);

   gROOT->SetStyle("MyStyle"); //comment to unset this style

   bool foundIt=true;
   // see if CMSSW has been setup
   char *cmsbase=gSystem->Getenv("CMSSW_BASE");
   if (cmsbase==NULL) {
     cout << " CMSSW environment has not been setup -- "
	  << " FWLite libraries will not be loaded\n" << endl;
     foundIt=false;
   } else {
     cout << " CMSSW environment has been setup \n" << endl;

     char *search=gSystem->Getenv("LD_LIBRARY_PATH");
     string cms_path = search;
     
     TString FWLiteLib = "libFWCoreFWLite.so";
     const char* foundlib =gSystem->Which(search, FWLiteLib, 0);
     
     if (! foundlib) {
       FWLiteLib = "libPhysicsToolsFWLite.so";
       foundlib =gSystem->Which(search, FWLiteLib, 0);
       if (! foundlib) {
	 cout << "Could not find any FWLite libraries to load " << endl;       
	 foundIt=false;
       }
     }
   }
   if (foundIt){
     //cout << "Loading: " << FWLiteLib << endl;
     //gSystem->Load(FWLiteLib);
     //AutoLibraryLoader::enable();
     gSystem->Load("libFWCoreFWLite.so");
     FWLiteEnabler::enable()

   }
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
}

/*{
 gSystem->Load("libFWCoreFWLite.so");
 AutoLibraryLoader::enable();
 gSystem->Load("libDataFormatsFWLite.so");
 gROOT->SetStyle ("Plain");
 gSystem->Load("libRooFit") ;
 using namespace RooFit ;
 cout << "loaded" << endl;
}*/
