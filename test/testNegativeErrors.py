import ROOT
import sys,os
import numpy

print "-> inserting in path cwd"
sys.path.insert(0,os.getcwd())
print "-> inserting in path cwd/python"
sys.path.insert(0,os.getcwd()+'/python')

from library import *

loadBootStrap()

print "-> Creating Matrixes"

reco, truth, resp = ConstructFromTree()
data= ConstructData(reco)

print "-> construct BootStrap"
nReg=5
b = ROOT.BootStrap()
b.SetUnfoldType(ROOT.BootStrap.kBayes) ## BootStrap
b.SetRegParam(nReg) ##BootStrap

b.SetNToys(1000)
b.SetSeed(328956)
b.SetUMatrix(reco,truth,resp)
b.SetData( data.Clone('bootstrap_data') )
b.SetToyType(ROOT.BootStrap.kBootstrap)
# try also kBootStrap, kMatrix -> smear the matrix
#b.SetSumW2(); ## MATRIX

# 
b.negCorr = ROOT.BootStrap.kNegNone
#kNegNone, kNegZero, kNegZeroProp, kNegMoveProp, kNegReplProp

#error = ROOT.BootStrap.kMin
#error = ROOT.BootStrap.kMedian
error = ROOT.BootStrap.kRms
print "-> running BootStrap I"
b.run()
gNone = b.result(error,.68)

print "-> running BootStrap II"
bZ = ROOT.BootStrap(b);
bZ.negCorr = ROOT.BootStrap.kNegZero
bZ.run()

gZero = bZero.result(error,.68) 
gZero = ROOT.utils.Shift( gZero, -0.3, True)

print "-> running BootStrap III"
bRepl = ROOT.BootStrap(b)
bRepl.negCorr= ROOT.BootStrap.kNegReplProp
bRepl.run()

gRepl = bRepl.result(error,.68)
gRepl = ROOT.utils.Shift(gRepl,+.3,True)



print "-> plotting"
c= ROOT.TCanvas("c1","c1",600,800)
p1 = ROOT.TPad("pad1","pad1", 0,.2,1,1)
p2 = ROOT.TPad("pad2","pad2", 0,0,1,0.2)
p1.Draw()
p2.Draw()


p1.SetTopMargin(0.02)
p1.SetBottomMargin(0.08)
p1.SetLeftMargin(0.08)
p1.SetRightMargin(0.02)

p2.SetTopMargin(0.02)
p2.SetBottomMargin(0.35)
p2.SetLeftMargin(0.08)
p2.SetRightMargin(0.02)

p1.cd()

truth.SetLineColor(ROOT.kRed)

gNone.SetLineColor(ROOT.kGreen+2)
gNone.SetMarkerStyle(20)
gNone.SetMarkerSize(0.8)
gNone.SetMarkerColor(ROOT.kGreen+2)

gZero.SetLineColor(ROOT.kBlue+2)
gZero.SetMarkerStyle(24)
gZero.SetMarkerSize(0.8)
gZero.SetMarkerColor(ROOT.kBlue+2)

gResp.SetLineColor(ROOT.kMagenta+2)
gResp.SetMarkerStyle(29)
gResp.SetMarkerColor(ROOT.kMagenta+2)

truth.GetYaxis().SetLabelFont(43)
truth.GetYaxis().SetLabelSize(26)
truth.GetXaxis().SetLabelFont(43)
truth.GetXaxis().SetLabelSize(26)

truth.Draw("HIST")

gNone.Draw("PE SAME")
gZero.Draw("PE SAME")
gRepl.Draw("PE SAME")

truth.Draw("AXIS SAME")
truth.Draw("AXIS X+ Y+SAME")

l = ROOT.TLegend(0.6,.6,.98,.98)
l.SetFillStyle(0)
l.SetBorderSize(0)
l.AddEntry(truth,"truth","L")
l.AddEntry(gNone,"kNegNone ","PE")
l.AddEntry(gZero,"kNegZero ","PE")
l.AddEntry(gRepl,"kNegRepl ","PE")

l.Draw()

p2.cd()

truth_r = ROOT.utils.Ratio(truth,truth)
truth_r.Draw("HIST")
truth_r.GetYaxis().SetRangeUser(0.5,1.5)

gNone_r = ROOT.utils.Ratio(gNone,truth)
gNone_r.Draw("PE SAME")

gZero_r = ROOT.utils.Ratio(gZero,truth)
gZero_r.Draw("PE SAME")

gRepl_r = ROOT.utils.Ratio(gRepl,truth)
gRepl_r.Draw("PE SAME")

truth_r.Draw("AXIS SAME")
truth_r.Draw("AXIS X+ Y+ SAME")

raw_input("ok?")
