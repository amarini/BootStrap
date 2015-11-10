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

reco, truth, resp = ConstructFromTree(100000,100)
#data= ConstructData(reco)
data= reco.Clone("reco_clone")

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
#b.SetToyType(ROOT.BootStrap.kToy)
# try also kBootStrap, kMatrix -> smear the matrix
#b.SetSumW2(); ## MATRIX

# 
b.negCorr = ROOT.BootStrap.kNegNone
#kNegNone, kNegZero, kNegZeroProp, kNegMoveProp, kNegReplProp

error = ROOT.BootStrap.kMin
#error = ROOT.BootStrap.kMedian
#error = ROOT.BootStrap.kRms
print "-> running BootStrap I"
b.run()
gNone = b.result(error,.68)

print "-> running BootStrap II"
bZ = ROOT.BootStrap(b);
bZ.negCorr = ROOT.BootStrap.kNegZero
bZ.run()

gZero = bZ.result(error,.68) 
gZero = ROOT.utils.Shift( gZero, -0.3, True)

#print "-> running BootStrap III"
#bRepl = ROOT.BootStrap(b)
#bRepl.negCorr= ROOT.BootStrap.kNegReplProp
#bRepl.run()
#
#gRepl = bRepl.result(error,.68)
#gRepl = ROOT.utils.Shift(gRepl,+.3,True)

print "-> running BootStrap IV"
bScale = ROOT.BootStrap(b)
bScale.negCorr= ROOT.BootStrap.kNegScaleProp
bScale.run()

gScale = bScale.result(error,.68)
gScale = ROOT.utils.Shift(gScale,+.3,True)


print "-> plotting"

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

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

#gRepl.SetLineColor(ROOT.kMagenta+2)
#gRepl.SetMarkerStyle(29)
#gRepl.SetMarkerColor(ROOT.kMagenta+2)

gScale.SetLineColor(ROOT.kMagenta+2)
gScale.SetMarkerStyle(29)
gScale.SetMarkerColor(ROOT.kMagenta+2)

truth.GetYaxis().SetLabelFont(43)
truth.GetYaxis().SetLabelSize(26)
truth.GetXaxis().SetLabelFont(43)
truth.GetXaxis().SetLabelSize(26)

truth.Draw("HIST")

gNone.Draw("PE SAME")
gZero.Draw("PE SAME")
#gRepl.Draw("PE SAME")
gScale.Draw("PE SAME")

truth.Draw("AXIS SAME")
truth.Draw("AXIS X+ Y+SAME")

l = ROOT.TLegend(0.6,.6,.98,.98)
l.SetFillStyle(0)
l.SetBorderSize(0)
l.AddEntry(truth,"truth","L")
l.AddEntry(gNone,"kNegNone ","PE")
l.AddEntry(gZero,"kNegZero ","PE")
#l.AddEntry(gRepl,"kNegRepl ","PE")
l.AddEntry(gScale,"kNegScale ","PE")

l.Draw()

## truth function
fPt = ROOT.TF1("f1","[0]*1e+6*TMath::Power(x,-3)*TMath::Exp(-50./x)",0,200)
fPt.SetParameter(0,1);
I=fPt.Integral(0,200);
## up canvas
fPt.SetParameter(0, 5000.*.4 * 2.0 / I ) ; ## the total area is 5000 ## last number is the bin width
#print "--- fPt Integral is",fPt.Integral(0,200)
#print "--- truth integral is", truth.Integral()

fPt2=fPt.Clone("myclone")
fPt2.Draw(" L  SAME");

p2.cd()
p2.SetGridy()
fPt.SetParameter(0, 5000.*.4  / I ) ; ## the total area is 5000 ## last number is the bin width, not here

truth_r = ROOT.utils.Ratio(truth,fPt)
truth_r.Draw("HIST")
truth_r.GetYaxis().SetRangeUser(0.45,1.55)
truth_r.GetYaxis().SetNdivisions(504);

gNone_r = ROOT.utils.Ratio(gNone,fPt)
gNone_r.Draw("PE SAME")

gZero_r = ROOT.utils.Ratio(gZero,fPt)
gZero_r.Draw("PE SAME")

#gRepl_r = ROOT.utils.Ratio(gRepl,fPt)
#gRepl_r.Draw("PE SAME")

gScale_r = ROOT.utils.Ratio(gScale,fPt)
gScale_r.Draw("PE SAME")

truth_r.Draw("AXIS SAME")
truth_r.Draw("AXIS X+ Y+ SAME")


ROOT.utils.ChangePalette(0)
c2 =ROOT.TCanvas("c2","c2",800,800)
c2.Divide(2,2)
c2.cd(1)
title = ROOT.TLatex()
title.SetTextFont(62)
title.SetTextSize(0.05)
title.SetTextAlign(22)
title.SetNDC()

b.GetUMatrixResp().Draw("COLZ")
title.DrawLatex(.5,.94,"None")

c2.cd(2)

bZ.GetUMatrixResp().Draw("COLZ")
title.DrawLatex(.5,.94,"Zero")

c2.cd(3)

#bRepl.GetUMatrixResp().Draw("COLZ")
#title.DrawLatex(.5,.94,"Repl")
bScale.GetUMatrixResp().Draw("COLZ")
title.DrawLatex(.5,.94,"Scale")

c2.cd(4)

data.Draw("PE")
reco.Draw("HIST SAME")

raw_input("ok?")

c.SaveAs("plotNeg/canv1.pdf")
c2.SaveAs("plotNeg/canv2.pdf")
