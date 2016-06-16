#! /usr/bin/env python 
import ROOT as rt

"""
Style options mostly from CMS's tdrStyle.C
"""
def customROOTstyle() :
    rt.gROOT.SetBatch(True)
    rt.gStyle.SetOptTitle(True)
    rt.gStyle.SetOptStat(False)
    rt.gStyle.SetPadTopMargin(0.06);
    rt.gStyle.SetPadBottomMargin(0.13);
    rt.gStyle.SetPadLeftMargin(0.13);
    rt.gStyle.SetPadRightMargin(.18)
    rt.gStyle.SetLabelColor(1, "XYZ");
    rt.gStyle.SetLabelFont(42, "XYZ");
    rt.gStyle.SetLabelOffset(0.007, "XYZ");
    rt.gStyle.SetLabelSize(0.05, "XYZ");
    rt.gStyle.SetTitleSize(0.05, "XYZ");
    rt.gStyle.SetTitleOffset(1.0, "X");
    rt.gStyle.SetTitleOffset(1.0, "Y");
    rt.gStyle.SetTitleOffset(1.0, "Z");
    rt.gStyle.SetAxisColor(1, "XYZ");
    rt.gStyle.SetStripDecimals(True);
    rt.gStyle.SetTickLength(0.03, "XYZ");
    rt.gStyle.SetNdivisions(510, "XYZ");
    rt.gStyle.SetPadTickX(0);
    rt.gStyle.SetPadTickY(0);
    rt.gStyle.SetMarkerStyle(20);
    rt.gStyle.SetHistLineColor(1);
    rt.gStyle.SetHistLineStyle(1);
    rt.gStyle.SetHistLineWidth(1);
    rt.gStyle.SetFrameBorderMode(0);
    rt.gStyle.SetFrameBorderSize(1);
    rt.gStyle.SetFrameFillColor(0);
    rt.gStyle.SetFrameFillStyle(0);
    rt.gStyle.SetFrameLineColor(1);
    rt.gStyle.SetFrameLineStyle(1);
    rt.gStyle.SetFrameLineWidth(1);
    rt.gStyle.SetPalette(55);
    rt.gStyle.SetNumberContours(100);

import numpy as np
def customPalette(zeropoint = 0.5):
    Number = 3;
    Red    = np.array([0  ,  100,  110],dtype=float)/255.
    Green  = np.array([0  ,  255,  0], dtype=float)/255.
    Blue   = np.array([99 ,  100,  2], dtype=float)/255.
    Length = np.array([0.0,  zeropoint, 1.0], dtype=float)
    nb=100;
    rt.TColor.CreateGradientColorTable(Number,Length,Red,Green,Blue,nb)

def doSpam(text,x1,y1,x2,y2,align=12,fill=False,textSize=0.033,_noDelete={}):
    cmsprel = rt.TPaveText(x1,y1,x2,y2,"NDC");
    cmsprel.SetTextSize(textSize);
    cmsprel.SetFillColor(0);
    cmsprel.SetFillStyle(1001 if fill else 0);
    cmsprel.SetLineStyle(2);
    cmsprel.SetLineColor(0);
    cmsprel.SetTextAlign(align);
    cmsprel.SetTextFont(42);
    cmsprel.AddText(text);
    cmsprel.Draw("same");
    _noDelete[text] = cmsprel; ## so it doesn't get deleted by PyROOT                                                                                                             
    return cmsprel

def doTinyCmsPrelim(textLeft="_default_",textRight="_default_",hasExpo=False,textSize=0.033,lumi=None, xoffs=0):
    global options
    if textLeft  == "_default_": textLeft  = "CMS Preliminary"
    if textRight == "_default_": textRight = "(13 TeV)"
    if lumi      == None       : lumi      = 0
    if   lumi > 3.54e+1: lumitext = "%.0f fb^{-1}" % lumi
    elif lumi > 3.54e+0: lumitext = "%.1f fb^{-1}" % lumi
    elif lumi > 3.54e-1: lumitext = "%.2f fb^{-1}" % lumi
    elif lumi > 3.54e-2: lumitext = "%.0f pb^{-1}" % (lumi*1000)
    elif lumi > 3.54e-3: lumitext = "%.1f pb^{-1}" % (lumi*1000)
    else               : lumitext = ""
    textLeft = textLeft.replace("%(lumi)",lumitext)
    textRight = textRight.replace("%(lumi)",lumitext)
    if textLeft not in ['', None]:
        doSpam(textLeft, (.28 if hasExpo else .15)+xoffs, .955, .65+xoffs, .995, align=12, textSize=textSize)
    if textRight not in ['', None]:
        doSpam(textRight,.75+xoffs, .955, .90+xoffs, .995, align=32, textSize=textSize)

