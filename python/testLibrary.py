import ROOT
import sys,os

print "-> Loading Library"
ROOT.gSystem.Load("bin/libBootStrap.so") # it will know where it is RooUnfold.so and load it

print "-> Creating Matrixes"

r=ROOT.TRandom3(23149)
N=20

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
			gen . SetBinContent(iBin, NSig*ROOT.TMath.Power(iBin,-5) ) 
		elif type == 2:
			gen . SetBinContent(iBin, NSig*ROOT.TMath.Power(iBin,-4.9) + NSig/10000. *ROOT.TMath.Power(iBin,-3) ) 
	return gen

## construct probability smears
def ConstructSmear():
	''' Generate a probability smear arbitrary. Sum Gen <=1'''
	smear= ROOT.TH2D("smear","smear",N,0.,1.,N,0.,1.)
	for i in range(0,N):
	   for j in range(0,N):
		iBin= i+1
		jBin= j+1
		if (iBin == jBin):smear.SetBinContent(iBin,jBin, 0.7) 
		if (abs(iBin-jBin) == 1):smear.SetBinContent(iBin,jBin, 0.08) 
		if (abs(iBin-jBin) == 2):smear.SetBinContent(iBin,jBin, 0.03) 
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
		data.SetBinContent(iBin, r.Poisson(data.GetBinContent(iBin) ) ) 
	return data

bkg  = ConstructBackground()
gen  = ConstructTruth()
smear= ConstructSmear()

resp = ConstructResponse(gen,smear)
reco = ConstructReco(bkg,resp)
data = ConstructData(reco)

## different model from the truth
gen2 = ConstructTruth(2)
resp2= ConstructResponse(gen2,smear)
reco2= ConstructReco(bkg,resp2)
data2= ConstructData(reco2)

print "-> construct Unfolding"
nReg = 1
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
b.SetUnfoldType(ROOT.BootStrap.kBayes)
b.SetRegParam(nReg)
b.SetUMatrix(reco,gen,resp)
b.SetData(data2)

print "-> running BootStrap"
b.run()
## kStd/kMin/kMedian/kMean
g_bootstrap = b.result(ROOT.BootStrapBase.kMin);

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

gen.Draw("HIST")
bkg.Draw("HIST SAME")

h_bayes.Draw("PE2 SAME")
g_bootstrap.Draw("PE SAME")
gen.Draw("HIST SAME")
#closure.Draw("HIST SAME")

reco.Draw("HIST SAME")
data.Draw("P SAME")

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


### DEBUG 
SuperDebug=False
if SuperDebug:
	''' For these debug private/prodected members needs to be public'''
	c2=ROOT.TCanvas("BootStrap","BootStrap",800,800)
	c2.Divide(2,2)
	
	c2.cd(1).SetLogy()
	
	fold = b.Fold(gen)
	fold.SetLineColor(ROOT.kMagenta)
	fold.SetLineStyle(7)
	reco.Draw("HIST")
	fold.Draw("HIST SAME")
	
	c2.cd(2).SetLogy()
	
	b.f_resp_bkg_.Draw("HIST")
	bkg.Draw("HIST SAME")
	
	c2.cd(3)
	b.f_resp_smear_.SetLineColor( ROOT.kRed)
	b.f_resp_smear_.Draw("BOX")
	
	smear.SetFillColor(ROOT.kBlue-4)
	smear.Draw("BOX SAME")
	
	c2.cd(4).SetLogy()
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


