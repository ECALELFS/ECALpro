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

    def __hash__(self):
        return hash((self.x, self.y, self.z))

    def __eq__(self, other):
        return (self.x, self.y, self.z) == (other.x, other.y, other.z)

    def __ne__(self, other):
        return not(self == other)

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
        print "loading data from file ",txt
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

    def plotIC2D(self,partition,zhwidth=0.07,errwidth=0.005,outdirname=''):
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
            if outdirname == '':
                canv.SaveAs('%s.pdf' % p[2].GetName())
                canv.SaveAs('%s.png' % p[2].GetName())
            else:
                canv.SaveAs('%s/%s.pdf' % (outdirname, p[2].GetName()))
                canv.SaveAs('%s/%s.png' % (outdirname, p[2].GetName()))

        for p in plots2D:
            xsize = 1200
            ysize = int(xsize*170/360+0.1*xsize) if 'EcalBarrel' in h.GetName() else int(xsize*0.9)
            canv = rt.TCanvas("c","",xsize,ysize)
            if p.GetDimension()==2: p.Draw("colz")
            else: p.Draw()
            if outdirname == '':
                canv.SaveAs('%s.pdf' % p.GetName())
                canv.SaveAs('%s.png' % p.GetName())
            else:
                canv.SaveAs('%s/%s.pdf' % (outdirname, p.GetName()))
                canv.SaveAs('%s/%s.png' % (outdirname, p.GetName()))


    def compareIC2D(self,data2,partition,zwidth=0.07,outdirname=''):
        #rt.gStyle.SetOptStat(0)
        customROOTstyle()
        plots = []

        if partition=='EcalBarrel': 
            h = rt.TProfile2D(('%s_%s_icratio_2d' % (self.name,partition)), '',360,1,360,170,-85,85)
            h.GetXaxis().SetTitle('i#phi')
            h.GetYaxis().SetTitle('i#eta')
        else: 
            h = rt.TProfile2D(('%s_%s_icratio_2d' % (self.name,partition)), '',100,1,100,100,1,100)
            h.GetXaxis().SetTitle('ix')
            h.GetYaxis().SetTitle('iy')

        zmin=1-zwidth; zmax=1+zwidth
        if partition=='EcalBarrel': 
            h1d = rt.TH1D(str(h.GetName()).replace('icratio_2d','icratio_1d'),'',200,zmin,zmax)
        else: 
            h1d = rt.TH1D(str(h.GetName()).replace('icratio_2d','icratio_1d'),'',200,zmin,zmax)
            
        h1d.GetXaxis().SetTitle('IC ratio')
        h1d.SetLineColor(rt.kRed)

        for k,v in self.data.iteritems():
            if k.subdet() != partition: continue
            if k not in data2: continue
            icref = data2[k]
            if(v.staterr < 999): 
                h.Fill(k.y,k.x,max(zmin,min(zmax,v.val/icref.val)))
                h1d.Fill(v.val/icref.val)
            
        h.GetZaxis().SetRangeUser(zmin,zmax)
        
        plots.append(h)

        for p in plots:
            xsize = ysize = 1200
            if p.GetDimension()==2: ysize = int(xsize*170/360+0.1*xsize) if 'EcalBarrel' in h.GetName() else int(xsize*0.9)
            canv = rt.TCanvas("c","",xsize,ysize)
            if p.GetDimension()==2: p.Draw("colz")
            else: p.Draw()
            if outdirname == '':
                canv.SaveAs('%s.pdf' % p.GetName())
                canv.SaveAs('%s.png' % p.GetName())
            else:
                canv.SaveAs('%s/%s.pdf' % (outdirname, p.GetName()))
                canv.SaveAs('%s/%s.png' % (outdirname, p.GetName()))
                

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] tag1.txt tag2.txt")
    parser.add_option("-n","--name", dest="name",  type="string", default='', help="prefix to give to the figures")
    parser.add_option("--max-EB", dest="max_EB",  type="float", default=0.05, help="width of the z-axis for IC maps for barrel")
    parser.add_option("--max-EE", dest="max_EE",  type="float", default=0.5, help="width of the z-axis for IC maps for endcaps")
    parser.add_option("--max-err-EB", dest="max_err_EB",  type="float", default=0.005, help="width of the y-axis for IC 1D error profile for barrel")
    parser.add_option("--max-err-EE", dest="max_err_EE",  type="float", default=0.2, help="width of the y-axis for IC 1D error profile for endcaps")
    parser.add_option("--noEB", dest="exclude_EB", action="store_true", default=False, help="ignore barrel (useful when you only produced IC for endcap")
    parser.add_option("--noEE", dest="exclude_EE", action="store_true", default=False, help="ignore endcap (useful when you only produced IC for barrel")
    parser.add_option("-o","--output-dir", dest="output_dir",  type="string", default='', help="output directory where plots are stored")

    (options, args) = parser.parse_args()
    if len(args) < 1: raise RuntimeError, 'Expecting at least the tag txt file'

    inputfile = args[0]
    name = inputfile.split('.')[0] if options.name == '' else options.name
    if options.name == '':
        name = name.split('/')[-1]  # get file name removing path (take last element after splitting on "/")

    if options.output_dir != '':
        print "Creating local folder to store output --> " + options.output_dir 
        folderCreation = subprocess.Popen(['mkdir -p ' + options.output_dir], stdout=subprocess.PIPE, shell=True);
        folderCreation.communicate()

    icp = ICplotter(inputfile,name)
    if not options.exclude_EB:
        icp.plotIC2D('EcalBarrel',options.max_EB,options.max_err_EB,options.output_dir)
    if not options.exclude_EE:
        icp.plotIC2D('EcalEndcapMinus',options.max_EE,options.max_err_EE,options.output_dir)
        icp.plotIC2D('EcalEndcapPlus',options.max_EE,options.max_err_EE,options.output_dir)

    if len(args)==2:    
        comparefile = args[1]
        data2 = icp.loadICs(comparefile)
        if not options.exclude_EB:
            icp.compareIC2D(data2,'EcalBarrel',options.max_EB,options.output_dir)
        if not options.exclude_EE:
            icp.compareIC2D(data2,'EcalEndcapMinus',options.max_EE,options.output_dir)
            icp.compareIC2D(data2,'EcalEndcapPlus',options.max_EE,options.output_dir)

