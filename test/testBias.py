import ROOT
import sys,os
import numpy

print "-> inserting in path cwd"
sys.path.insert(0,os.getcwd())
print "-> inserting in path cwd/python"
sys.path.insert(0,os.getcwd()+'/python')

from library import *


print "-> Loading Library"
loadBootStrap()
print "-> Creating Matrixes"

bkg  = ConstructBackground()
gen  = ConstructTruth()
smear= ConstructSmear(2) ## 1 ~ diagonal, 2 ~ off-diagonal, 3 low eff

resp = ConstructResponse(gen,smear)
reco = ConstructReco(bkg,resp)
data = ConstructData(reco)

print "-> construct Unfolding"
nReg = 5
print "\t\tBayes, nReg=",nReg
R = ROOT.RooUnfoldResponse(reco,gen,resp)
u = ROOT.RooUnfoldBayes(R,data,nReg)
u.SetVerbose(-1)
u.SetNToys(1000)
h_bayes = u.Hreco(ROOT.RooUnfold.kCovToy)

canvas=[]
garbage=[]

ROOT.gStyle.SetOptStat(0)

for n in range(-4,5):
  canv=ROOT.TCanvas("BiasStudy_n%d"%(n),"C",800,800)
  canv.Divide(2,2)
  canvas.append(canv)
  for cidx,delta in enumerate([0,0.05,0.10,.25]):
	canv.cd( cidx+1 ) 
	gen2=gen.Clone("gen")
	mumod=gen.Clone("mu_n%d_delta%.2f"%(n,delta))
	mumod.SetLineColor(ROOT.kRed)
	garbage.append(mumod)
	for i in range(0,gen2.GetNbinsX()):
		iBin=i+1
		x = ( float(i)+.5) / float(gen2.GetNbinsX()) * ROOT.TMath.Pi() # in [0,pi]
		if n>=0: f = 1. + delta * ROOT.TMath.Cos(n*x)
		else : f = 1. + delta * ROOT.TMath.Sin(-n*x)
		mumod.SetBinContent(iBin,f)
		mumod.SetBinError(iBin,0)
		c = gen.GetBinContent(iBin)
		gen2.SetBinContent(iBin, c*f)
		gen2.SetBinError(iBin, 0 )
	resp2 = ConstructResponse(gen2,smear)
	reco2 = ConstructReco(bkg,resp2)
	u = ROOT.RooUnfoldBayes(R, reco2 ,nReg) ## unfold the reco2, no smearing with nominal R
	u.SetVerbose(-1)
	u.SetNToys(1000)
	h_bayes = u.Hreco(ROOT.RooUnfold.kCovToy)
	h_bayes.SetName("unfold_toy_n%d_delta%.2f"%(n,delta))
	h_bayes.Divide(gen2)
	h_bayes.Draw("HIST")
	mumod.Draw("HIST SAME")
	h_bayes.GetYaxis().SetRangeUser(0.8,1.2)
	garbage.append(h_bayes)
  #canvas.cd(0)
  #latex=ROOT.TLatex()
  #latex.SetNDC()
  #latex.SetTextFont(62)
  #latex.SetTextSize(0.04)

raw_input("ok?")	
for c in canvas:
	c.SaveAs("plot/"+c.GetName() + ".pdf" )
	c.SaveAs("plot/"+c.GetName() + ".png" )
