#!/bin/env python

import ROOT, os, sys, re, array, math
import time

ROOT.gROOT.SetBatch(True)

isPi0 = False
inputdir = "/eos/cms/store/group/dpg_ecal/alca_ecalcalib/piZero_Run2/mciprian/"
hname = "pi0pt"
particle = "#pi^{0}" if isPi0 else "#eta^{0}"
xTitle = particle + " p_{T} [GeV]"
years = ["2016", "2017", "2018"]
dirs = {"2016" : "AlCaP0_2016_ULrereco_from0" if isPi0 else "AlCaEta_2016_ULrereco",
        "2017" : "AlCaP0_AllRun2017_condor_fixEBm16" if isPi0 else "AlCaEta_2017_ULrereco_all2017data",
        "2018" : "AlCaP0_2018_ULrereco_1every2" if isPi0 else "AlCaEta_2018_ULrereco_all2018data"
}

detIds = ["1EB", "2EB", "1EE", "2EE"]  # region 3EE missing for older files
legTitles = ["Barrel: |#eta| < 1.0", 
            "Barrel: |#eta| > 1.0",
            "Endcap: |#eta| < 1.8",
            "Endcap: 1.8 < |#eta| < 2.0"]


colors = [ROOT.kBlack, ROOT.kRed+2, ROOT.kBlue]

def createPlotDirAndCopyPhp(outdir):
    if outdir != "./":
        if not os.path.exists(outdir):
            os.system("mkdir -p "+outdir)
            if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/e/emanuele/public/index.php "+outdir)


if __name__ == "__main__":
            
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('-o','--outdir',      dest='outdir',      default='', type='string', help='output directory to save things')
    parser.add_option('-x','--xAxisTitle',  dest='xAxisTitle',  default='', type='string', help='X axis title. If not given, use the one from hist1')
    parser.add_option('-y','--yAxisTitle',  dest='yAxisTitle',  default='arbitrary units', type='string', help='Y axis title. If not given, use the one from hist1. if "arbitrary units", normalize histograms to unit area')
    parser.add_option('-e','--energy', dest='energy', default='13', type='int', help='Center of mass energy')
    (options, args) = parser.parse_args()

    #if len(sys.argv) < 3:
    #    parser.print_usage()
    #    quit()

    ROOT.TH1.SetDefaultSumw2()

    print ""

    if options.outdir:
        outname = options.outdir
        if not outname.endswith("/"): outname = outname + "/"
        createPlotDirAndCopyPhp(outname)
    else:
        print "Error: you should specify an output folder using option -o <name>. Exit"
        quit()

    for jdet,detId in enumerate(detIds):

        print ""
        print "-"*30
        print detId
        
        files = []
        for y in years:
            files.append("{ind}/{d}/iter_0/{d}_epsilonPlots.root".format(ind=inputdir,d=dirs[y]))

        histos = {}

        for i,y in enumerate(years):
            n = dirs[y]
            tf = ROOT.TFile.Open(files[i])
            h = tf.Get(hname + "_afterCuts_region" + detId)
            h.SetDirectory(0)
            histos[n] = h.Clone(n)
            histos[n].SetDirectory(0)
            tf.Close()
            print ">>> File: %s" % files[i]

        leftmargin = 0.16
        rightmargin = 0.06
        topmargin = 0.06
        canvas = ROOT.TCanvas("canvas","",1200,900)
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

        # get max and min, set colors, etc
        miny = 0.0 # 10000.0
        maxy = -10000.0
        for i,y in enumerate(years):
            n = dirs[y]
            histos[n].SetLineColor(colors[i])
            #histos[n].SetFillColor(colors[i])
            histos[n].SetLineWidth(2)
            histos[n].SetTitle("")
            histos[n].SetStats(0)
            if options.yAxisTitle == "arbitrary units":
                histos[n].Scale(1./histos[n].Integral())
            #miny = min(miny,histos[n].GetBinContent(histos[n].GetMinimumBin())) # this does not exclude 0
            # the following is needed to exclude 0
            # for ib in range(1,1+histos[n].GetNbinsX()):
            #     bc = histos[n].GetBinContent(ib)
            #     if bc > 0.0:
            #         if bc < miny: miny = bc
            maxy = max(maxy,histos[n].GetBinContent(histos[n].GetMaximumBin()))

        leg = ROOT.TLegend(0.6,0.6,0.9,0.9)
        leg.SetFillColor(0)
        #leg.SetFillStyle(0)
        #leg.SetBorderSize(0)
        leg.SetHeader(legTitles[jdet])
        for i,y in enumerate(years):
            n = dirs[y]
            if i == 0:
                histos[n].Draw("HE")
                histos[n].GetYaxis().SetRangeUser(0.9*miny,1.2*maxy)
                histos[n].GetYaxis().SetTitle(options.yAxisTitle)
                histos[n].GetYaxis().SetTitleSize(0.05)
                histos[n].GetYaxis().SetLabelSize(0.04)
                histos[n].GetXaxis().SetTitle(xTitle)
                histos[n].GetXaxis().SetTitleSize(0.05)
                histos[n].GetXaxis().SetLabelSize(0.04)
            else:
                histos[n].Draw("HESAME")
            leg.AddEntry(histos[n],y,"LF")
        leg.Draw("same")

        canvas.RedrawAxis("sameaxis")

        latCMS = ROOT.TLatex()
        latCMS.SetNDC();
        latCMS.SetTextFont(42)
        latCMS.SetTextSize(0.04)
        latCMS.DrawLatex(leftmargin, 0.95, '#bf{CMS} #it{Preliminary}')
        latCMS.DrawLatex(0.82, 0.95, '(%s TeV)' % str(options.energy))

        cname = hname + "_comparisonRun2_region" + detId + ("_pi0" if isPi0 else "_eta0")
        for ext in [".png", ".pdf"]:
            canvas.SaveAs(outname + cname + ext)

        print "-"*30
        print ""
