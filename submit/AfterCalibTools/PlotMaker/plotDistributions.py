#!/bin/env python

import ROOT, os, sys, re, array, math
import time

# hist1/2 can have format  "hname::legeEntry", where hname is the actual name inside the file, the rest is the text for the legend

ROOT.gROOT.SetBatch(True)

def createPlotDirAndCopyPhp(outdir):
    if outdir != "./":
        if not os.path.exists(outdir):
            os.system("mkdir -p "+outdir)
            if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/e/emanuele/public/index.php "+outdir)


if __name__ == "__main__":
            
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] file1 hist1 hist2')
    parser.add_option('-o','--outdir',      dest='outdir',      default='', type='string', help='output directory to save things')
    parser.add_option('-c','--cname', dest='cname', default='', type='string', help='Name of canvas (without extension, .png and .pdf added automatically)')
    #parser.add_option('-f','--outfilename', dest='outfilename', default='', type='string', help='Name of output file to save results')
    #parser.add_option('-n','--outhistname', dest='outhistname', default='', type='string', help='Name of output histogram saved in output file. If FILE is used, take same name as the output file, removing the extension')
    parser.add_option('-x','--xAxisTitle',  dest='xAxisTitle',  default='', type='string', help='X axis title. If not given, use the one from hist1')
    parser.add_option('-y','--yAxisTitle',  dest='yAxisTitle',  default='arbitrary units', type='string', help='Y axis title. If not given, use the one from hist1. if "arbitrary units", normalize histograms to unit area')
    parser.add_option('-l','--lumi', dest='lumi', default=None, type='float', help='Value to be used for luminosity')
    parser.add_option('-e','--energy', dest='energy', default='13', type='int', help='Center of mass energy')
    parser.add_option(     '--scale-hist2', dest='scaleHist2', default='0', type='int', help='If non-zero, scale histogram 2 to adjust y axis maximum. If > 0, multiply integral, if < 0 divide integral ')
    parser.add_option('-t','--legendTitle',  dest='legendTitle',  default='', type='string', help='Title for legend (e.g. can specify eta region)')
    (options, args) = parser.parse_args()

    if len(sys.argv) < 3:
        parser.print_usage()
        quit()

    f1 = args[0]
    hist1 = args[1]
    hist2 = args[2]
    if hist2 == "PTPI0": hist2 = hist1

    leg1 = "h1"
    leg2 = "h2"
    if "::" in hist1:
        hist1,leg1 = hist1.split("::")
    if "::" in hist2:
        hist2,leg2 = hist2.split("::")

    ROOT.TH1.SetDefaultSumw2()

    print ""

    if options.outdir:
        outname = options.outdir
        if not outname.endswith("/"): outname = outname + "/"
        createPlotDirAndCopyPhp(outname)
    else:
        print "Error: you should specify an output folder using option -o <name>. Exit"
        quit()
    # if not options.outfilename:
    #     print "Error: you should specify an output file name using option -f <name>. Exit"
    #     quit()
    # if not options.outhistname:
    #     print "Error: you should specify an output histogram name using option -n <name>. "
    #     print "If FILE is used, take same name as the output file, removing the extension"
    #     print "Exit"
    #     quit()

    # if options.outhistname == "FILE":
    #     options.outhistname = options.outfilename.split('.')[0]

    lumi = options.lumi

    # file 1
    tf = ROOT.TFile.Open(f1)        
    h1 =   tf.Get(hist1)
    h2 =   tf.Get(hist2)
    if (h1 == 0):
        print "Error: could not retrieve %s from input file %s. Exit" % (h1,f1)
        quit()
    else:
        h1.SetDirectory(0)
    if (h2 == 0):
        print "Error: could not retrieve %s from input file %s. Exit" % (h2,f1)
        quit()
    else:
        h2.SetDirectory(0)
    tf.Close()

    leftmargin = 0.16
    rightmargin = 0.06
    topmargin = 0.06
    canvas = ROOT.TCanvas("canvas","",800,800)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.SetGridx(1)
    canvas.SetGridy(1)
    canvas.cd()
    canvas.SetLeftMargin(leftmargin)
    canvas.SetRightMargin(rightmargin)
    canvas.SetTopMargin(topmargin)
    canvas.SetBottomMargin(0.12)
    canvas.cd()

    h1.SetTitle("")
    h2.SetTitle("")

    h1.SetLineWidth(2)
    h1.SetLineColor(ROOT.kOrange+2)
    h1.SetFillColor(ROOT.kOrange+1)
    h1.SetFillStyle(3001)
    h2.SetLineWidth(2)
    h2.SetLineColor(ROOT.kBlue+2)
    h2.SetFillColor(ROOT.kAzure+1)
    h2.SetFillStyle(3244)

    if options.xAxisTitle: h1.GetXaxis().SetTitle(options.xAxisTitle)
    h1.GetXaxis().SetTitleOffset(1.1)
    h1.GetXaxis().SetTitleSize(0.05)
    h1.GetXaxis().SetLabelSize(0.04)
    if options.yAxisTitle: h1.GetYaxis().SetTitle(options.yAxisTitle)
    h1.GetYaxis().SetTitleOffset(1.44) 
    h1.GetYaxis().SetTitleSize(0.05)
    h1.GetYaxis().SetLabelSize(0.04)
    #h1.GetYaxis().SetTickSize(0.01)
    
    h1.SetStats(0)
    h2.SetStats(0)

    if options.yAxisTitle == "arbitrary units":
        h1.Scale(1./h1.Integral())
        if hist1 != hist2: h2.Scale(1./h2.Integral())
    leg2ext = ""
    if hist1 != hist2 and options.scaleHist2:
        if options.scaleHist2 > 0: 
            h2.Scale(options.scaleHist2)
            leg2ext = " x " + str(options.scaleHist2)
        else:
            h2.Scale(1./(abs(options.scaleHist2)))
            leg2ext = " / " + str(abs(options.scaleHist2))            

    maxy1 = h1.GetBinContent(h1.GetMaximumBin())
    #maxy2 = h2.GetBinContent(h2.GetMaximumBin())
    #maxy = max(max1,max2)
    h1.GetYaxis().SetRangeUser(0., 1.15*maxy1)

    h1.Draw("HIST")
    if hist1 != hist2: h2.Draw("HIST SAME")

    canvas.RedrawAxis("sameaxis")

    leg = ROOT.TLegend(0.69-rightmargin,(0.75 if options.legendTitle else 0.8) -topmargin,1-rightmargin,1-topmargin)
    leg.SetFillColor(0)
    #leg.SetFillStyle(0)
    #leg.SetBorderSize(0)
    if options.legendTitle:
        #leg.SetHeader("#bf{%s}" % options.legendTitle)
        leg.AddEntry(None,"","")
    leg.AddEntry(h1,leg1,"LF")
    if hist1 != hist2: leg.AddEntry(h2,leg2+leg2ext,"LF")
    leg.Draw("same")

    latCMS = ROOT.TLatex()
    latCMS.SetNDC();
    latCMS.SetTextFont(42)
    latCMS.SetTextSize(0.04)
    if options.legendTitle:
        latCMS.DrawLatex(0.70 -rightmargin, 1-topmargin - 0.05, "#bf{%s}" % options.legendTitle)
    latCMS.DrawLatex(leftmargin, 0.95, '#bf{CMS} #it{Preliminary}')
    if lumi != None: latCMS.DrawLatex(0.68, 0.95, '%s fb^{-1} (%s TeV)' % (lumi, str(options.energy)))
    else:            latCMS.DrawLatex(0.75, 0.95, '(%s TeV)' % str(options.energy))

    for ext in [".png", ".pdf", ".C"]:
        canvas.SaveAs(outname + options.cname + ext)
