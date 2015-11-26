import ROOT
import sys,os
import numpy
print "-> inserting in path cwd"
sys.path.insert(0,os.getcwd())
print "-> inserting in path cwd/python"
sys.path.insert(0,os.getcwd()+'/python')

from library import *


loadBootStrap()


## h=ROOT.TH1D("Njets","Njets",4,0-.5,4-.5)
## h.SetBinContent(1,20.2)
## h.SetBinError(1,9.2)
## 
## h.SetBinContent(2,4.3)
## h.SetBinError(2,6.4)
## 
## h.SetBinContent(3,3.9)
## h.SetBinError(3,3.8)
## 
## h.SetBinContent(4,2.5)
## h.SetBinError(4,2.2)
## 
## mc = h.Clone("mc")
## mc.SetBinContent(1,17.4)
## mc.SetBinContent(2,8.7)
## mc.SetBinContent(3,3.1)
## mc.SetBinContent(4,1.1)
## 
## for i in range(0,4): mc.SetBinError(i+1,0)
## 
## corr =  ROOT.TH2D("corr","corr",4,0-.5,4-.5,4,0-.5,4-.5)
## 
## #diag
## corr.SetBinContent(1,1,1)
## corr.SetBinContent(2,2,1)
## corr.SetBinContent(3,3,1)
## corr.SetBinContent(4,4,1)
## 
## ## off -diag
## corr.SetBinContent(1,2,-.18)
## corr.SetBinContent(2,1,-.18)
## corr.SetBinContent(2,3,-.08)
## corr.SetBinContent(3,2,-.08)
## corr.SetBinContent(3,4,-.11)
## corr.SetBinContent(4,3,-.11)
## 
## ## 
## corr.SetBinContent(1,3,-0.05)
## corr.SetBinContent(3,1,-0.05)
## corr.SetBinContent(2,4,-0.24)
## corr.SetBinContent(4,2,-0.24)
## 
## corr.SetBinContent(1,4,0.10)
## corr.SetBinContent(4,1,0.10)

if True: # pT 
	from array import array
	ptBins=array('f',[0,15,26,43,72,125,200,250])
	h= ROOT.TH1D("pToMscaled","pToMscaled",len(ptBins)-1,ptBins)
	xsec= [9.0,2.0,3.4,6.2,4.6,2.6,0.7]
	error= [6.3,5.2,4.7,3.6,2.55,1.0,0.5]
	xsecmc=[7.5,5.9,6.1,5.2,3.4,1.4,0.6]
	mc = h.Clone("mc")
	for i in range(0, len(ptBins)-1 ) :
		h.SetBinContent(i+1, xsec[i] ) 
		h.SetBinError(i+1, error[i])
		mc.SetBinContent(i+1,xsecmc[i])
		mc.SetBinError(i+1,0)
	corr =  ROOT.TH2D("corr","corr", len(ptBins) -1 , ptBins,len(ptBins)-1,ptBins)
	c= { (0,1): -0.11, (0,2):0.07, (1,2) : -.12 , (0,3):0.04, (2,3):-0.04, (1,4):-0.07,
			(2,4):0.05, (0,5):-0.05, (1,5):0.03, (2,5):-0.04,(3,5):0.03,(4,5):0.01,
			(5,6):0.05, (4,6):-0.03,(3,6):0.04,(2,6):0.04} 
	for i in range(0, len(ptBins)-1 ) :
	    for j in range(0, len(ptBins)-1 ) :
		    if i==j: corr.SetBinContent(i+1,j+1,1.0)
		    elif (i,j) in c : corr.SetBinContent(i+1,j+1,c[ (i,j) ] )
		    elif (j,i) in c : corr.SetBinContent(i+1,j+1,c[ (j,i) ] )
		    else: corr.SetBinContent(i+1,j+1,0.)


h.Divide(mc)
cov = corr.Clone("cov")
for i in range(0,h.GetNbinsX()):
   for j in range(0,h.GetNbinsX()):
	   c = corr.GetBinContent(i+1,j+1)
	   e1 = h.GetBinError(i+1)
	   e2 = h.GetBinError(j+1)
	   cov.SetBinContent(i+1,j+1, e1*e2*c)

r = ROOT.Regularize()
r.useOverFlow=False

reg=[0.0,0.01,0.10,1.0,10.,100.,1000]
#reg=[0.0]
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
ROOT.gStyle.SetOptStat(0)

h.Multiply(mc)
h.Draw("P E2")
h.SetMarkerStyle(29)
h.SetMarkerColor(ROOT.kBlack)
h.SetLineColor(ROOT.kBlack)
h.SetFillStyle(3004)
h.SetFillColor(ROOT.kBlack)
h.GetYaxis().SetRangeUser(-5,30)

garbage=[]
for idx,h2 in enumerate(hist):
	#f = ( (idx%2)*(-1) ) * (idx/2 +1) * 0.1
	f= -0.4 + 0.1*idx
	g = ROOT.utils.Shift(h2,f,True)
	garbage.append(g)
	g.Draw("PE SAME")
	#h2.Draw("PE SAME")


shift=0.05
canv.BuildLegend(0.5 + shift , 0.67 + shift, .88 + shift,.88 +shift)

raw_input("ok?")

canv.SaveAs("plot/Reg_" +h.GetName()+ ".pdf")
canv.SaveAs("plot/Reg_" +h.GetName()+ ".png")
