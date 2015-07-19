import ROOT
import sys,os
import numpy

print "-> Loading Library"
ROOT.gSystem.Load("bin/libBootStrap.so") # it will know where it is RooUnfold.so and load it

print "-> Creating Matrixes"

r=ROOT.TRandom3(23149)
N=30

NBkg=1.e3
NSig=1e7


## Construct generator and bkg spectra
def ConstructBackground():
	''' Generate a Background distribution accordingly to an arbitrary fuction'''
	bkg  = ROOT.TH1D("bkg" ,"bkg" ,N,0.,1.)
	for i in range(0,N):
		iBin= i+1
		bkg . SetBinContent(iBin, NBkg*  abs(ROOT.TMath.Cos( float(iBin)/N * 3.14)) * ROOT.TMath.Power(iBin,-4) ) 
	return bkg

def ConstructTruth(type=1):
	''' Generate a Truth distribution accordingly to an arbitrary function'''
	gen  = ROOT.TH1D("gen%d"%type ,"gen" ,N,0.,1.)
	for i in range(0,N):
		iBin= i+1
		if type == 1:
			gen . SetBinContent(iBin, NSig*ROOT.TMath.Power(iBin * 20. / N,-5) ) 
		elif type == 2:
			gen . SetBinContent(iBin, NSig*ROOT.TMath.Power(iBin * 20./ N,-4.9) + NSig/10000. *ROOT.TMath.Power(iBin,-3) ) 
	return gen

## construct probability smears
def ConstructSmear(type=1):
	''' Generate a probability smear arbitrary. Sum Gen <=1'''
	smear= ROOT.TH2D("smear","smear",N,0.,1.,N,0.,1.)
	for i in range(0,N):
	   for j in range(0,N):
		iBin= i+1
		jBin= j+1
		if type==1: ## pretty normal
			if (iBin == jBin):smear.SetBinContent(iBin,jBin, 0.7) 
			if (abs(iBin-jBin) == 1):smear.SetBinContent(iBin,jBin, 0.08) 
			if (abs(iBin-jBin) == 2):smear.SetBinContent(iBin,jBin, 0.03) 
		elif type ==2: #pretty off -diagonal
			if (iBin == jBin):smear.SetBinContent(iBin,jBin, 0.5) 
			if (abs(iBin-jBin) == 1):smear.SetBinContent(iBin,jBin, 0.15) 
			if (abs(iBin-jBin) == 2):smear.SetBinContent(iBin,jBin, 0.08) 
			if (abs(iBin-jBin) == 3):smear.SetBinContent(iBin,jBin, 0.02) 
	return smear

## construct response matrix
def ConstructResponse(gen,smear):
	''' Construct the response matrix folding truth and smear'''
	resp = ROOT.TH2D("resp."+gen.GetName()+"."+smear.GetName(),"resp",N,0.,1.,N,0.,1.)
	for i in range(0,N):
	   for j in range(0,N):
		iBin= i+1
		jBin= j+1
		resp.SetBinContent(iBin,jBin, smear.GetBinContent(iBin,jBin) * gen.GetBinContent(jBin) ) 
	return resp

## construct reco
def ConstructReco(bkg,resp):
	''' Construct the Reconstructed distribution adding to the projection of the response the bkg distribution'''
	reco = ROOT.TH1D("reco."+bkg.GetName()+"."+resp.GetName(),"reco",N,0.,1.)
	for i in range(0,N):
	   S=0
	   for j in range(0,N):
		iBin= i+1
		jBin= j+1
		S += resp.GetBinContent(iBin,jBin)			
	   reco.SetBinContent(iBin, S + bkg.GetBinContent(iBin) ) 
	return reco

def ConstructData(reco):
	''' Construct toy data smearing the reco distribution with poisson.'''
	data = reco.Clone("data."+reco.GetName())
	for i in range(0,N):
		iBin= i+1
		c= r.Poisson(data.GetBinContent(iBin) )
		if c<0 : print "* ERROR: negative poisson!"
		data.SetBinContent(iBin, c ) 
	return data

ROOT.gROOT.ProcessLine (\
		"struct RGB{ \
		float r;\
		float g;\
		float b;\
		void SetRGB(int color) { gROOT->GetColor(color)->GetRGB(r,g,b); } \
		};" )
from ROOT import RGB

## def rgb(color):
## 	r= float(0.)
## 	g= float(0.)
## 	b= float(0.)
## 
## 	ROOT.gROOT.GetColor(color).GetRGB(r,g,b)
## 	return (r,g,b)

class color():
	def __init__(self):
		self.r = 0
		self.g = 0
		self.b = 0
		self.rgb = RGB(0)
	def __init__(self,color):
		self.rgb = RGB()
		self.rgb.SetRGB(color)
		self.r = self.rgb.r
		self.g = self.rgb.g
		self.b = self.rgb.b

####################

bkg  = ConstructBackground()
gen  = ConstructTruth()
smear= ConstructSmear(2) ## 1 ~ diagonal, 2 ~ off-diagonal

resp = ConstructResponse(gen,smear)
reco = ConstructReco(bkg,resp)
data = ConstructData(reco)

## different model from the truth
gen2 = ConstructTruth(2)
resp2= ConstructResponse(gen2,smear)
reco2= ConstructReco(bkg,resp2)
data2= ConstructData(reco2)

print "-> Using no bias in the matrix"
data2=data

print "-> construct Unfolding"
nReg = 5
print "\t\tBayes, nReg=",nReg

R = ROOT.RooUnfoldResponse(reco,gen,resp)
u = ROOT.RooUnfoldBayes(R,data2,nReg)
u.SetVerbose(-1)
u.SetNToys(1000)
h_bayes = u.Hreco(ROOT.RooUnfold.kCovToy)

#uInv = ROOT.RooUnfoldInvert(R,reco)
#closure = uInv.Hreco()

print "-> construct BootStrap"
b = ROOT.BootStrap()
b.SetNToys(1000)
b.SetSeed(328956)
b.SetUnfoldType(ROOT.BootStrap.kBayes)
b.SetRegParam(nReg)
b.SetUMatrix(reco,gen,resp)
b.SetData(data2.Clone('bootstrap_data'))


print "-> running BootStrap"
## ## kStd/kMin/kMedian/kMean
#g_bootstrap = b.result(ROOT.BootStrapBase.kMedian,.68)
#g_bootstrap = b.result(ROOT.BootStrapBase.kMedian,.58) ### ALMOST 68%: 1 / 10. / .58
b.run()
#g_bootstrap = b.result(ROOT.BootStrapBase.kStd,.68)
g_bootstrap = b.result(ROOT.BootStrap.kMedian,.68)

g_bootstrap = ROOT.utils.Shift( g_bootstrap, 0.3, True)

print "-> plotting"

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

c= ROOT.TCanvas()
c.SetTopMargin(0.02)
c.SetBottomMargin(0.08)
c.SetLeftMargin(0.08)
c.SetRightMargin(0.02)

gen.SetLineColor(ROOT.kRed)
h_bayes.SetLineColor(ROOT.kGreen+2)
h_bayes.SetFillColor(ROOT.kGreen)

g_bootstrap.SetLineColor(ROOT.kBlue+2)
g_bootstrap.SetMarkerColor(ROOT.kBlue+2)
g_bootstrap.SetMarkerStyle(20)
g_bootstrap.SetMarkerSize(0.8)

bkg.SetLineColor(ROOT.kGray+1)
bkg.SetLineStyle(2)
bkg.SetFillColor(ROOT.kGray)

#closure.SetLineColor(ROOT.kMagenta)
#closure.SetLineStyle(7)

reco.SetLineColor(ROOT.kBlack)

data.SetMarkerStyle(24)
data.SetMarkerColor(ROOT.kMagenta)

gen.Draw("AXIS")
bkg.Draw("HIST SAME")

h_bayes.Draw("PE2 SAME")
g_bootstrap.Draw("PE SAME")
gen.Draw("HIST SAME")
#closure.Draw("HIST SAME")

reco.Draw("HIST SAME")
data.Draw("P SAME")
gen.Draw("AXIS SAME")
gen.Draw("AXIS X+ Y+ SAME")

l = ROOT.TLegend(0.6,.6,.98,.98)
l.SetFillStyle(0)
l.SetBorderSize(0)
l.AddEntry(gen,"truth","L")
l.AddEntry(h_bayes,"bayes","LF")
l.AddEntry(g_bootstrap,"bootstrap","PE")
l.AddEntry(reco,"reco","L")
l.AddEntry(bkg,"bkg","LF")
l.AddEntry(data,"data","P")
l.Draw()

print "-> Coverage Test"
tot=0
bootstrap= 0
bayes=0
for i in range(1,h_bayes.GetNbinsX()+1):
	g = gen.GetBinContent(i)

	c1 = h_bayes.GetBinContent(i)
	e1 = h_bayes.GetBinError(i)

	tot += 1
	if abs(g- c1) < e1: bayes+=1
	
	c2=g_bootstrap.GetY()[i-1]
	e2h=g_bootstrap.GetEYhigh()[i-1]
	e2l=g_bootstrap.GetEYlow()[i-1]

	if (g >= c2 and g < c2+e2h) or ( g<=c2 and g> c2-e2l) : bootstrap+=1

print " TOT      = ",tot
print " Bayes    = ",bayes," | ", float(bayes)/tot*100, "%"
print " Bootstrap= ",bootstrap," | ", float(bootstrap)/tot*100, "%"

c.SetLogy()


ROOT.utils.ChangePalette(1)
c2 = ROOT.TCanvas("c2","c2",600,10,800,600)

c2.Divide(2)
c2.cd(1)

cov = u.Ereco(ROOT.RooUnfold.kCovToy)
err = u.ErecoV(ROOT.RooUnfold.kCovToy)

corr_bayes=ROOT.TH2D("corr_bayes","corr_bayes",N,0-.5,N-.5,N,-.5,N-.5)

for i in range(0, N ):
   for j in range(0, N):
	corr_bayes.SetBinContent(i+1,j+1, cov(i,j) / (err(i)*err(j) ) )

ltx=ROOT.TLatex()
ltx.SetNDC()
ltx.SetTextAlign(22)
corr_bayes.Draw("COLZ")
corr_bayes.GetZaxis().SetRangeUser(-1,1)
ltx.DrawLatex(.5,.92,"Bayes Correlation Matrix")

c2.cd(2)
corr = b.correlation()
corr.Draw("COLZ")
corr.GetZaxis().SetRangeUser(-1,1)
ltx.DrawLatex(.5,.92,"Bootstrap Correlation Matrix")

### DEBUG 
SuperDebug=False
if SuperDebug:
	''' For these debug private/prodected members needs to be public'''
	c_debug=ROOT.TCanvas("BootStrap","BootStrap",800,800)
	c_debug.Divide(2,2)
	
	c_debug.cd(1).SetLogy()
	
	fold = b.Fold(gen)
	fold.SetLineColor(ROOT.kMagenta)
	fold.SetLineStyle(7)
	reco.Draw("HIST")
	fold.Draw("HIST SAME")
	
	c_debug.cd(2).SetLogy()
	
	b.f_resp_bkg_.Draw("HIST")
	bkg.Draw("HIST SAME")
	
	c_debug.cd(3)
	b.f_resp_smear_.SetLineColor( ROOT.kRed)
	b.f_resp_smear_.Draw("BOX")
	
	smear.SetFillColor(ROOT.kBlue-4)
	smear.Draw("BOX SAME")
	
	c_debug.cd(4).SetLogy()
	gen.Draw("HIST")
	b1=b.bootStrap()
	b1.SetName("toy_1")
	b1.Draw("HIST SAME")
	
	b2=b.bootStrap()
	b2.SetName("toy_2")
	b2.Draw("HIST SAME")
	
	b1.SetLineColor(ROOT.kBlack)
	b2.SetLineColor(ROOT.kBlue)

raw_input("Looks ok ? ")


