#!/usr/bin/env python

from shutil import copyfile
import re, sys, os, os.path, subprocess, json, ROOT, copy, math
import numpy as np
from array import array

from CMS_lumi import *

ROOT.gROOT.SetBatch(True)

#########################################################################


def addStringToEnd(name, matchToAdd, notAddIfEndswithMatch=False):
    if notAddIfEndswithMatch and name.endswith(matchToAdd):
        return name
    elif not name.endswith(matchToAdd):
        return name + matchToAdd

#########################################################################

def createPlotDirAndCopyPhp(outdir):
    if outdir != "./":
        if not os.path.exists(outdir):
            os.system("mkdir -p "+outdir)
            if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/e/emanuele/public/index.php "+outdir)
    

#########################################################################

def getAxisRangeFromUser(axisNameTmp="", 
                         separator="::", 
                         rangeSeparator=","
                         ):
  
    setXAxisRangeFromUser = False;
    fields = axisNameTmp.split(separator)
    axisName = fields[0]
    
    if len(fields) > 1:
        setXAxisRangeFromUser = True;
        xmin = float(fields[1].split(rangeSeparator)[0])
        xmax = float(fields[1].split(rangeSeparator)[1])
    else:
        xmin = 0
        xmax = 0
        
    return axisName,setXAxisRangeFromUser,xmin,xmax


#########################################################################

def adjustSettings_CMS_lumi():

    ## dummy function to be called before using any other fucntion calling CMS_lumi
    ## for some reason, the settings of the very first plot are screwed up.
    ## To fix this issue, it is enough to call it to a dummy plot
    dummy = ROOT.TH1D("dummy","",10,0,10)
    for i in range(1,1+dummy.GetNbinsX()):
        dummy.SetBinContent(i,i)
    dummy.GetXaxis().SetTitle("x axis")
    dummy.GetYaxis().SetTitle("y axis")
    cdummy = ROOT.TCanvas("cdummy","",600,600)
    dummy.Draw("HE")
    CMS_lumi(cdummy,"",True,False)
    setTDRStyle()        
    ## no need to save the canvas    


#########################################################################


# function to draw 2D histograms, can also plot profile along X on top
def drawCorrelationPlot(h2D_tmp,
                        labelXtmp="xaxis", labelYtmp="yaxis", labelZtmp="zaxis",
                        canvasName="default", plotLabel="", outdir="./",
                        rebinFactorY=0,
                        rebinFactorX=0,
                        smoothPlot=True,
                        drawProfileX=True,
                        scaleToUnitArea=True,
                        draw_both0_noLog1_onlyLog2=0,
                        leftMargin=0.16,
                        rightMargin=0.20,
                        nContours=51,
                        palette=55,
                        canvasSize="700,625",
                        passCanvas=None,
                        bottomMargin=0.1,
                        plotError=False,
                        lumi=None,
                        drawOption = "colz",
                        profileLegendColor=1):


    # if h2D.GetName() == "scaleFactor_origBinPt":
    #     print "="*20
    #     print "Check: hist %s: Z axis title = %s" % (h2D.GetName(),labelZtmp)
    #     print "="*20

    ROOT.TH1.SetDefaultSumw2()
    adjustSettings_CMS_lumi()

    if (rebinFactorX): 
        if isinstance(rebinFactorX, int): h2D_tmp.RebinY(rebinFactorX)
        else:                             h2D_tmp.RebinY(len(rebinFactorX)-1,"",array('d',rebinFactorX)) # case in which rebinFactorX is a list of bin edges

    if (rebinFactorY): 
        if isinstance(rebinFactorY, int): h2D_tmp.RebinY(rebinFactorY)
        else:                             h2D_tmp.RebinY(len(rebinFactorY)-1,"",array('d',rebinFactorY)) # case in which rebinFactorX is a list of bin edges

    if plotError:
        herr = h2D_tmp.Clone(h2D_tmp.GetName()+"_err")
        herr.Reset("ICESM")
        for i in range(1,herr.GetNbinsX()+1):
            for j in range(1,herr.GetNbinsY()+1):
                herr.SetBinContent(i,j,h2D_tmp.GetBinError(i,j))
        h2D = herr
    else:
        h2D = h2D_tmp

    ROOT.TColor.CreateGradientColorTable(3,
                                         array ("d", [0.00, 0.50, 1.00]),
                                         ##array ("d", [1.00, 1.00, 0.00]),        
                                         ##array ("d", [0.70, 1.00, 0.34]),        
                                         ##array ("d", [0.00, 1.00, 0.82]),        
                                         array ("d", [0.00, 1.00, 1.00]),
                                         array ("d", [0.34, 1.00, 0.65]),
                                         array ("d", [0.82, 1.00, 0.00]),
                                         255,  0.95)

    if palette > 0: ROOT.gStyle.SetPalette(palette)  # 55:raibow palette ; 57: kBird (blue to yellow, default) ; 107 kVisibleSpectrum ; 77 kDarkRainBow 
    ROOT.gStyle.SetNumberContours(nContours) # default is 20 

    labelX,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    labelY,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)
    labelZ,setZAxisRangeFromUser,zmin,zmax = getAxisRangeFromUser(labelZtmp)
    
    cw,ch = canvasSize.split(',')
    #canvas = ROOT.TCanvas("canvas",h2D.GetTitle() if plotLabel == "ForceTitle" else "",700,625)    
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.SetBottomMargin(bottomMargin)
    canvas.cd()

    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)
    # normalize to 1
    if (scaleToUnitArea): h2D.Scale(1./h2D.Integral())

    h2DGraph = 0

    h2DPlot = 0
    if (not smoothPlot): h2DPlot = h2D
    else:
        h2DGraph = ROOT.TGraph2D()
        h2DGraph.SetNpx(300)
        h2DGraph.SetNpy(300)
        nPoint = 0
        for iBinX in range (1,1+h2D.GetNbinsX()):
            for iBinY in range(1,1+h2D.GetNbinsY()):
                h2DGraph.SetPoint(nPoint,h2D.GetXaxis().GetBinCenter(iBinX),h2D.GetYaxis().GetBinCenter(iBinY),h2D.GetBinContent(iBinX,iBinY))
                nPoint += 1
            

        h2DPlot = h2DGraph.GetHistogram()

    if plotLabel == "ForceTitle":
        h2DPlot.SetTitle(h2D_tmp.GetTitle())
  
    h2DPlot.GetXaxis().SetTitle(labelX)
    h2DPlot.GetYaxis().SetTitle(labelY)
    h2DPlot.GetXaxis().SetTitleSize(0.05)
    h2DPlot.GetXaxis().SetLabelSize(0.04)
    h2DPlot.GetXaxis().SetTitleOffset(0.95) # 1.1 goes outside sometimes, maybe depends on root version or canvas width
    h2DPlot.GetYaxis().SetTitleSize(0.05)
    h2DPlot.GetYaxis().SetLabelSize(0.04)
    h2DPlot.GetYaxis().SetTitleOffset(1.1)
    h2DPlot.GetZaxis().SetTitleSize(0.05)
    h2DPlot.GetZaxis().SetLabelSize(0.04)
    h2DPlot.GetZaxis().SetTitleOffset(1.2)

    h2DPlot.GetZaxis().SetTitle(labelZ) 
    h2DPlot.Draw(drawOption)

    if (setXAxisRangeFromUser): h2DPlot.GetXaxis().SetRangeUser(xmin,xmax)
    if (setYAxisRangeFromUser): h2DPlot.GetYaxis().SetRangeUser(ymin,ymax)
    if (setZAxisRangeFromUser): h2DPlot.GetZaxis().SetRangeUser(zmin,zmax)


    # if h2D.GetName() == "scaleFactor_origBinPt":
    #     print "="*20
    #     print "Check: hist %s: Z axis title = %s" % (h2DPlot.GetName(),h2DPlot.GetZaxis().GetTitle())
    #     print "="*20

    # # attempt to make Z axis title farther depending on how many digits are printed
    # maxZaxisVal = h2DPlot.GetBinContent(h2DPlot.GetMaximumBin())
    # if (setZAxisRangeFromUser): maxZaxisVal = zmax

    # if maxZaxisVal >= 1.0:
    #     rootYear = int(str(ROOT.gROOT.GetVersionDate())[:4])        
    #     if (rootYear > 2016):
    #         h2DPlot.GetZaxis().SetMaxDigits(3)
    #     else:
    #         print "Warning in drawCorrelationPlot: TAxis::SetMaxDigits() not implemented for ROOT versions before 2017 (rough estimate)"
    #         print "Will not exit, but instruction will be neglected"
    #     if maxZaxisVal > 9999.:
    #         h2DPlot.GetZaxis().SetTitleOffset(h2DPlot.GetZaxis().GetTitleOffset()+0.15)
    #         print "Changing title offset by 0.15"
    # else:
    #     i = 1
    #     tryNext = True
    #     while (tryNext and i < 6):
    #         tmpVal = maxZaxisVal * pow(10,i)
    #         if tmpVal >= 1.0: tryNext = False 
    #         else: i += 1
    #     if i > 1:            
    #         print "Max Z axis < 1, will try to adjust distance of Z axis title to Z axis"
    #         print "i = %d: will move Z axis offset by 0.45" % i
    #         # for numbers like 0.025 or with more 0 after ., make increase distance between Z axis title and the Z axis
    #         h2DPlot.GetZaxis().SetTitleOffset(h2DPlot.GetZaxis().GetTitleOffset()+0.45)

    h2DPlot.GetZaxis().SetTitleOffset(h2DPlot.GetZaxis().GetTitleOffset()+0.4)


    h2DProfile = 0
    if drawProfileX:
        h2DProfile = h2D.ProfileX("%s_pfx" %h2D.GetName())
        h2DProfile.SetMarkerColor(ROOT.kBlack)
        h2DProfile.SetMarkerStyle(20)
        h2DProfile.SetMarkerSize(1)
        h2DProfile.Draw("EPsame")
        
    # not yet implemented
    setTDRStyle()
    if not plotLabel == "ForceTitle": 
        if lumi != None: CMS_lumi(canvas,lumi,True,False)
        else:            CMS_lumi(canvas,"",True,False)
    #setTDRStyle()
    #print ">>>>>>>>>>>>>> check <<<<<<<<<<<<<<<<<<<"

    if plotLabel == "ForceTitle":
        ROOT.gStyle.SetOptTitle(1)        

    #h2DPlot.GetZaxis().SetMaxDigits(1)  #for N>99, should use scientific notation, I'd like to make it work only with negative exponential but haven't succeeded yet
    # canvas.Modified()
    # canvas.Update()

    leg = ROOT.TLegend(0.39,0.75,0.89,0.95)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(62)
    if plotLabel not in ["", "ForceTitle"]: leg.AddEntry(0,plotLabel,"")
    if drawProfileX: 
        l1 = leg.AddEntry(0,"Correlation = %.2f" % h2DPlot.GetCorrelationFactor(),"")
        l1.SetTextColor(profileLegendColor)
        #leg.AddEntry(h2DProfile,"Mean along Y","P")
    leg.Draw("same")

    if (draw_both0_noLog1_onlyLog2 == 0 or draw_both0_noLog1_onlyLog2 == 1):
        for ext in ['png', 'pdf']:
            canvas.SaveAs('{od}/{cn}.{ext}'.format(od=outdir, cn=canvasName, ext=ext))
        
    if (draw_both0_noLog1_onlyLog2 == 0 or draw_both0_noLog1_onlyLog2 == 2):
        canvas.SetLogz()
        for ext in ['png', 'pdf']:
            canvas.SaveAs('{od}/{cn}_logZ.{ext}'.format(od=outdir, cn=canvasName, ext=ext))
        canvas.SetLogz(0)


##########################################################



eospath = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero_Run2/mciprian/"
folder = "pi0CC_2018_EoverEtrue_foldSM_nFit10_onlyEB_fixGamma2EoverEtrue"
niter = 0
outdir = "/afs/cern.ch/user/m/mciprian//www/pi0calib/CC_EoverEtrue_2018/" + folder
histograms = ["pi0MassVsIetaEB",
              "pi0MassVsETEB",
              "photonDeltaRVsIetaEB"]

inputfile = eospath + folder + "/iter_" + str(niter) + "/" + folder + "_epsilonPlots.root"

createPlotDirAndCopyPhp(outdir)

tf = ROOT.TFile.Open(inputfile)  
if not tf:
    print "Error when opening file %s. Exit" % inputfile
    quit()

for h in histograms:
    hist = tf.Get(h)
    if not hist:
        print "Error when reading histogram %s from file. Exit" % h
        quit()
    drawCorrelationPlot(hist, hist.GetXaxis().GetTitle(), hist.GetYaxis().GetTitle(), "Events",
                        hist.GetName(), "ForceTitle", outdir, 0, 0, False, True, False, 1,
                        canvasSize = "1000,800",
                        bottomMargin = 0.12,
                        profileLegendColor=2)
