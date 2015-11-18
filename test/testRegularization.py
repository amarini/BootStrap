import ROOT
import sys,os
import numpy
print "-> inserting in path cwd"
sys.path.insert(0,os.getcwd())
print "-> inserting in path cwd/python"
sys.path.insert(0,os.getcwd()+'/python')

from library import *


loadBootStrap()


h=ROOT.TH1D("Njets","Njets",4,0-.5,4-.5)
h.SetBinContent(1,20.2)
h.SetBinError(1,9.2)

h.SetBinContent(2,4.3)
h.SetBinError(2,6.4)

h.SetBinContent(3,3.9)
h.SetBinError(3,3.8)

h.SetBinContent(4,2.5)
h.SetBinError(4,2.2)

mc = h.Clone("mc")
mc.SetBinContent(1,17.4)
mc.SetBinContent(2,8.7)
mc.SetBinContent(3,3.1)
mc.SetBinContent(4,1.1)
for i in range(0,4): mc.SetBinError(i+1,0)

corr =  ROOT.TH2D("corr","corr",4,0-.5,4-.5,4,0-.5,4-.5)

#diag
corr.SetBinContent(1,1,1)
corr.SetBinContent(2,2,1)
corr.SetBinContent(3,3,1)
corr.SetBinContent(4,4,1)

## off -diag
corr.SetBinContent(1,2,-.18)
corr.SetBinContent(2,1,-.18)
corr.SetBinContent(2,3,-.08)
corr.SetBinContent(3,2,-.08)
corr.SetBinContent(3,4,-.11)
corr.SetBinContent(4,3,-.11)

## 
corr.SetBinContent(1,3,-0.05)
corr.SetBinContent(3,1,-0.05)
corr.SetBinContent(2,4,-0.24)
corr.SetBinContent(4,2,-0.24)

corr.SetBinContent(1,4,0.10)
corr.SetBinContent(4,1,0.10)

h.Divide(mc)
cov = corr.Clone("cov")
for i in range(0,4):
   for j in range(0,4):
	   c = corr.GetBinContent(i+1,j+1)
	   e1 = h.GetBinError(i+1)
	   e2 = h.GetBinError(j+1)
	   cov.SetBinContent(i+1,j+1, e1*e2*c)

r = ROOT.Regularize()
r.useOverFlow=False

reg=[0.0,0.01,0.10,1.0,10.,100.,1000]
colors=[ROOT.kRed,ROOT.kRed+2,ROOT.kGreen+2,ROOT.kGreen,ROOT.kBlue,ROOT.kBlue+2,ROOT.kOrange]
hist=[]

for idx,delta in enumerate(reg):
	h2 = r.applyRegularization(h,cov, delta)
	h2.Multiply(mc)
	h2.SetName("Reg%f"%delta)
	h2.SetTitle("Reg%f"%delta)
	color = colors[idx % len(colors)]
	h2.SetLineColor(color)
	h2.SetMarkerColor(color)
	h2.SetMarkerStyle(20)
	h2.SetMarkerSize(0.8)
	hist.append(h2)

canv=ROOT.TCanvas("c","c",800,800)

for idx,h2 in enumerate(hist):
	if idx==0:
		h2.Draw("PE ")
		h2.GetYaxis().SetRangeUser(-5,30)
	else: h2.Draw("PE SAME")

h.Multiply(mc)
h.Draw("P E2 SAME")
h.SetMarkerStyle(29)
h.SetMarkerColor(ROOT.kBlack)
h.SetLineColor(ROOT.kBlack)
h.SetFillStyle(3004)
h.SetFillColor(ROOT.kBlack)

canv.BuildLegend()

raw_input("ok?")
