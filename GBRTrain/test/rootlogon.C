// --> >// $Id: rootlogon.C,v 1.1 2012/10/24 13:24:23 lpernie Exp $
{
 {
  TString libstr(Form("%s/lib/%s/%s",
                      gSystem->Getenv("CMSSW_BASE"),
                      gSystem->Getenv("SCRAM_ARCH"),
                      "libCondFormatsEgammaObjects.so"));

  gSystem->Load(libstr);
 }
 {
  TString libstr(Form("%s/lib/%s/%s",
                      gSystem->Getenv("CMSSW_BASE"),
                      gSystem->Getenv("SCRAM_ARCH"),
                      "libMitEdmGBRTrain.so"));

  gSystem->Load(libstr);
 }
  gSystem->AddIncludePath("-I$ROOFITSYS/include");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/CalibCode/MitEdm/CondFormats/EgammaObjects/interface");
  gSystem->AddIncludePath("-I$CMSSW_BASE/src/CalibCode/MitEdm/GBRTrain/interface");

  gInterpreter->AddIncludePath((TString(":")+TString(gSystem->Getenv("CMSSW_BASE"))+
				TString("/src/CalibCode/MitEdm/CondFormats/EgammaObjects/interface")).Data());
  gInterpreter->AddIncludePath((TString(":")+TString(gSystem->Getenv("CMSSW_BASE"))+
                                TString("/src/CalibCode/MitEdm/GBRTrain/interface")).Data());                                
  gInterpreter->AddIncludePath((TString(":")+TString(gSystem->Getenv("ROOFITSYS"))+
				TString("/include")).Data());

  gSystem->Load("libFWCoreFWLite.so");
  AutoLibraryLoader::enable();

}

