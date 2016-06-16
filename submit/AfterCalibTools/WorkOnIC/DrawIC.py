#!/usr/bin/env python
import ROOT as rt
import subprocess, time, sys, os, math
rt.gROOT.SetBatch(True)

from PlotUtils import customROOTstyle, customPalette

class IC:
    def __init__(self,infos):
        self.val = float(infos[0])
        self.staterr = float(infos[1])
        self.systerr = float(infos[2])
        self.chi2 = float(infos[3])
        self.toterr = math.sqrt(self.staterr**2 + self.systerr**2)

class EcalDetId:
    def __init__(self,xyz):
        self.x = int(xyz[0])
        self.y = int(xyz[1])
        self.z = int(xyz[2])

    def subdet(self):
        if self.z == 0: return 'EcalBarrel'
        elif self.z == 1: return 'EcalEndcapPlus'
        elif self.z == -1: return 'EcalEndcapMinus'
        else: return 'UNKNOWN'

class EtaRings:
    def __init__(self,mapfile):
        self.EndcapEtaRings={}
        for line in open(mapfile,'r'):
            vals = line.split(' ')
            key = "_".join(vals[:2])
            self.EndcapEtaRings[key] = int(vals[3])

    def etaring(self,detid):
        if(detid.subdet()=='EcalBarrel'): return abs(detid.x)
        else: 
            key = "_".join([str(detid.x+1),str(detid.y+1)])
            return self.EndcapEtaRings[key] if key in self.EndcapEtaRings else -1

class ICplotter:
    def loadICs(self, txt):
        data = {}
        for line in open(txt,'r'):
            if line.startswith('#'): continue
            vals = line.split()
            detid = EcalDetId(vals[:3])
            if float(vals[4]) < 998: ic = IC([vals[3],vals[5],vals[6],vals[7]])
            else: ic = IC([vals[3],999.,999.,999.])
            data[detid] = ic
        return data

    def __init__(self,icfile,name):
        self.data = self.loadICs(icfile)
        self.name = name

    def plotIC2D(self,partition,zhwidth=0.07,errwidth=0.005):
        #rt.gStyle.SetOptStat(0)
        customROOTstyle()
        plots2D = []
        profiles = {}

        if partition=='EcalBarrel': 
            h = rt.TProfile2D(('%s_%s_ic_2d' % (self.name,partition)), '',360,1,360,170,-85,85)
            h.GetXaxis().SetTitle('i#phi')
            h.GetYaxis().SetTitle('i#eta')
        else: 
            h = rt.TProfile2D(('%s_%s_ic_2d' % (self.name,partition)), '',100,1,100,100,1,100)
            h.GetXaxis().SetTitle('ix')
            h.GetYaxis().SetTitle('iy')

        if partition=='EcalBarrel': 
            hsterr = rt.TProfile(str(h.GetName()).replace('ic_2d','icsterr_1d'),'',85,0,84)
        else: 
            hsterr = rt.TProfile(str(h.GetName()).replace('ic_2d','icsterr_1d'),'',38,0,37)
            
        hsterr.GetXaxis().SetTitle('#eta ring')
        hsterr.GetYaxis().SetRangeUser(0,errwidth)
        hsterr.SetMarkerStyle(rt.kFullCircle)
        hsterr.SetMarkerColor(rt.kRed)
        hsterr.SetLineColor(rt.kRed)
        hsyerr = hsterr.Clone(str(h.GetName()).replace('ic_2d','icsyerr_1d'))
        hsyerr.SetMarkerStyle(rt.kFullTriangleUp)
        hsyerr.SetMarkerColor(rt.kGreen+2)
        hsyerr.SetLineColor(rt.kGreen+2)
        htoterr = hsterr.Clone(str(h.GetName()).replace('ic_2d','ictoterr_1d'))
        htoterr.SetMarkerStyle(rt.kFullSquare)
        htoterr.SetMarkerColor(rt.kBlack)
        htoterr.SetLineColor(rt.kBlack)
        htoterr.SetMarkerSize(1.5)
        profiles['errors'] = [hsterr,hsyerr,htoterr]

        ering = EtaRings('InputFile/Endc_x_y_ring.txt')

        zmin=1-zhwidth; zmax=1+zhwidth
        for k,v in self.data.iteritems():
            if k.subdet() != partition: continue
            if(v.staterr < 999): 
                h.Fill(k.y,k.x,max(zmin,min(zmax,v.val)))
                print "x=%d y=%d z=%d etaring=%d" % (k.x,k.y,k.z,ering.etaring(k))
                hsterr.Fill(ering.etaring(k),v.staterr)
                hsyerr.Fill(ering.etaring(k),v.systerr)
                htoterr.Fill(ering.etaring(k),v.toterr)
            
        h.GetZaxis().SetRangeUser(zmin,zmax)
        
        plots2D.append(h)

        leg = rt.TLegend(0.2,0.7,0.5,0.85)
        leg.SetFillColor(0)
        leg.SetShadowColor(0)
        leg.SetLineColor(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.03)

        for k,p in profiles.iteritems():
            canv = rt.TCanvas("c","",1200,1200)
            p[2].Draw("p")
            p[1].Draw("p same")
            p[0].Draw("p same")
            leg.AddEntry(p[0], 'Stat. error', 'LP')
            leg.AddEntry(p[1], 'Syst. error', 'LP')
            leg.AddEntry(p[2], 'Tot. error', 'LP')
            leg.Draw()
            canv.SaveAs('%s.pdf' % p[2].GetName())
            canv.SaveAs('%s.png' % p[2].GetName())
            

        for p in plots2D:
            xsize = 1200
            ysize = int(xsize*170/360+0.1*xsize) if 'EcalBarrel' in h.GetName() else int(xsize*0.9)
            canv = rt.TCanvas("c","",xsize,ysize)
            if p.GetDimension()==2: p.Draw("colz")
            else: p.Draw()
            canv.SaveAs('%s.pdf' % p.GetName())
            canv.SaveAs('%s.png' % p.GetName())
             
if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] tag1.txt tag2.txt")
    parser.add_option("-n","--name", dest="name",  type="string", default='', help="prefix to give to the figures")

    (options, args) = parser.parse_args()
    if len(args) < 1: raise RuntimeError, 'Expecting at least the tag txt file'

    inputfile = args[0]
    name = inputfile.split('.')[0] if options.name == '' else options.name

    icp = ICplotter(inputfile,name)
    icp.plotIC2D('EcalBarrel')
    icp.plotIC2D('EcalEndcapPlus',0.2,0.03)
    icp.plotIC2D('EcalEndcapMinus',0.2,0.03)

