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

    def plotIC2D(self,partition,zhwidth=0.07,errwidth=0.005,outdirname='', norm_etaring=False):
        #rt.gStyle.SetOptStat(0)
        customROOTstyle()
        plots2D = []
        profiles = {}

        if partition=='EcalBarrel': 
            h = rt.TProfile2D(('%s_%s_ic_2d' % (self.name,partition)), '',360,1,360,171,-85.5,85.5)
            h.GetXaxis().SetTitle('i#phi')
            h.GetYaxis().SetTitle('i#eta')
        else: 
            h = rt.TProfile2D(('%s_%s_ic_2d' % (self.name,partition)), '',100,1,100,100,1,100)
            h.GetXaxis().SetTitle('ix')
            h.GetYaxis().SetTitle('iy')

        if partition=='EcalBarrel': 
            hsterr = rt.TProfile(str(h.GetName()).replace('ic_2d','icsterr_1d'),'',85,0.5,85.5)
        else: 
            hsterr = rt.TProfile(str(h.GetName()).replace('ic_2d','icsterr_1d'),'',38,0.5,38.5)
            
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
                # for EB, the file has ieta in x, but in the histogram ieta is in the y axis
                if k.subdet() == "EcalBarrel":
                    h.Fill(k.y,k.x,max(zmin,min(zmax,v.val)))
                else:
                    h.Fill(k.x,k.y,max(zmin,min(zmax,v.val)))
                hsterr.Fill(ering.etaring(k),v.staterr)
                hsyerr.Fill(ering.etaring(k),v.systerr)
                htoterr.Fill(ering.etaring(k),v.toterr)
            
        h.GetZaxis().SetRangeUser(zmin,zmax)
        
        hnorm1 = h.Clone(str(h.GetName()).replace('ic_2d','ic_2d_norm1etaring'))
    
        if norm_etaring:
            if partition=='EcalBarrel':
                #hnorm1 = h.Clone(str(h.GetName()).replace('ic_2d','ic_2d_norm1etaring'))
                for ieta in range (1,172): # range excludes last value, so we have 171 values, but ieta = 0 (bin 86 for TH1) doesn't exist, so we have 170 eta rings in EB
                    if ieta == 86: 
                        continue
                    ICsum_etaring = 0.0
                    xtalsInEtaRing = 0.0
                    # if xtal is dead, bin content is 0, do not count it in the average
                    for iphi in range (1,361):
                        if h.GetBinContent(iphi,ieta) > 0.00001:
                            ICsum_etaring += h.GetBinContent(iphi,ieta)
                            xtalsInEtaRing += 1.0
                        #print "iphi, ieta, ICsum_etaring = %s %s %s" % (str(iphi), str(ieta), str(ICsum_etaring))
                    ICsum_etaring = ICsum_etaring / xtalsInEtaRing  # now it is the average in the eta-ring

                    average = 0.0
                    for iphi in range (1,h.GetNbinsX()+1):
                        if h.GetBinContent(iphi,ieta) != 1. and h.GetBinContent(iphi,ieta) > 0.00001:
                            # if the calibration costant is exactly 1., it means the fit failed or there was some other problem such that the IC was set to 1.0
                            # in this case, do not modify it dividing by the average (note that, since its IC is 1, it can still be used to compute the average)
                            hnorm1.SetBinContent(iphi,ieta,h.GetBinContent(iphi,ieta)/ICsum_etaring)
                        #average += h.GetBinContent(iphi,ieta)/(xtalsInEtaRing * ICsum_etaring)
                        average += hnorm1.GetBinContent(iphi,ieta)/xtalsInEtaRing  # can sum even if bin is 0 or 1 for the average
                    print "EB etaring %s: Nxtal %s   average %.4f" % (str(ieta-86),str(xtalsInEtaRing),average)

                plots2D.append(hnorm1)
            else:
                #hnorm1 = h.Clone(str(h.GetName()).replace('ic_2d','ic_2d_norm1etaring'))
                f_ietaring = rt.TFile("/afs/cern.ch/user/m/mciprian/public/ECALproTools/EE_xyzToEtaRing/eerings_modified.root")
                h_ietaring = f_ietaring.Get("hEEm") if "Minus" in partition else f_ietaring.Get("hEEp")              
                for ietaring in range (0,39):
                    ICsum_etaring = 0.0
                    xtalsInEtaRing = 0.0
                    for ix in range (1,101):
                        for iy in range (1,101):
                            if h_ietaring.GetBinContent(ix,iy) == ietaring:
                                if h.GetBinContent(ix,iy) > 0.00001:
                                    ICsum_etaring += h.GetBinContent(ix,iy)
                                    xtalsInEtaRing += 1.0
                    ICsum_etaring = ICsum_etaring / xtalsInEtaRing
                                
                    average = 0.0
                    for ix in range (1,101):
                        for iy in range (1,101):
                            if h_ietaring.GetBinContent(ix,iy) == ietaring:
                                if h.GetBinContent(ix,iy) != 1. and h.GetBinContent(ix,iy) > 0.00001:
                                    # if the calibration costant is exactly 1., it means the fit failed or there was some other problem such that the IC was set to 1.0
                                    # in this case, do not modify it dividing by the average (note that, since its IC is 1, it can still be used to compute the average)
                                    hnorm1.SetBinContent(ix,iy,h.GetBinContent(ix,iy)/ICsum_etaring)
                                #average += h.GetBinContent(ix,iy)/(xtalsInEtaRing * ICsum_etaring)                       
                                average += hnorm1.GetBinContent(ix,iy)/xtalsInEtaRing  # can sum even if bin is 0 or 1 for the average
                    print "EE etaring %s: Nxtal %s   average %.4f" % (str(ietaring),str(xtalsInEtaRing),average)

                plots2D.append(hnorm1)
                f_ietaring.Close()
        else:
            plots2D.append(h)

        leg = rt.TLegend(0.2,0.7,0.5,0.85)
        leg.SetFillColor(0)
        leg.SetShadowColor(0)
        leg.SetLineColor(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.03)

        for k,p in profiles.iteritems():
            canv = rt.TCanvas("c","",1200,1200)
            canv.cd();
            canv.SetTickx(1);
            canv.SetTicky(1);
            canv.cd();
            #canv.SetBottomMargin(0.1);
            canv.SetRightMargin(0.06);
            canv.SetLeftMargin(0.18);
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
            
            ## Save map in root file
            outfilename = "%s.root" % p.GetName()

            if outdirname == '':
                canv.SaveAs('%s.pdf' % p.GetName())
                canv.SaveAs('%s.png' % p.GetName())
            else:
                outfilename = outdirname + "/" + outfilename; 
                canv.SaveAs('%s/%s.pdf' % (outdirname, p.GetName()))
                canv.SaveAs('%s/%s.png' % (outdirname, p.GetName()))

            outfile = rt.TFile(outfilename,"RECREATE")
            p.Write()
            outfile.Close()


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
    parser.add_option("--normalize-etaring", dest="normalize_etaring", action="store_true", default=False, help="Normalize IC 2D map to 1 for each eta-ring")

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
        if os.path.exists("/afs/cern.ch"): 
            os.system("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+os.path.dirname(options.output_dir))

    icp = ICplotter(inputfile,name)
    if not options.exclude_EB:
        icp.plotIC2D('EcalBarrel',options.max_EB,options.max_err_EB,options.output_dir, options.normalize_etaring)
    if not options.exclude_EE:
        icp.plotIC2D('EcalEndcapMinus',options.max_EE,options.max_err_EE,options.output_dir,options.normalize_etaring)
        icp.plotIC2D('EcalEndcapPlus',options.max_EE,options.max_err_EE,options.output_dir,options.normalize_etaring)

    if len(args)==2:    
        comparefile = args[1]
        data2 = icp.loadICs(comparefile)
        if not options.exclude_EB:
            icp.compareIC2D(data2,'EcalBarrel',options.max_EB,options.output_dir)
        if not options.exclude_EE:
            icp.compareIC2D(data2,'EcalEndcapMinus',options.max_EE,options.output_dir)
            icp.compareIC2D(data2,'EcalEndcapPlus',options.max_EE,options.output_dir)

