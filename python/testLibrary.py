import ROOT
import sys,os


print "-> Loading RooUnfold"

ROOT.gSystem.Load("libRooUnfold.so")

print "-> Loading Library"
ROOT.gSystem.Load("bin/libBootStrap.so")

print "-> Creating Matrixes"

r=ROOT.TRandom3(23148)
N=20
resp = ROOT.TH2D("resp","resp",N,0.,1.,N,0.,1.)
reco = ROOT.TH1D("reco","reco",N,0.,1.)
gen  = ROOT.TH1D("gen" ,"gen" ,N,0.,1.)
bkg  = ROOT.TH1D("bkg" ,"bkg" ,N,0.,1.)
smear= ROOT.TH2D("smear","smear",N,0.,1.,N,0.,1.)


## Construct generator and bkg spectra
def ConstructBackround():
	for i in range(0,N):
		iBin= i+1
		bkg . SetBinContent(iBin, 10.*  abs(ROOT.TMath.Cos( float(iBin)/N * 3.14)) * ROOT.TMath.Power(iBin,-4) ) 

def ConstructTruth():
	for i in range(0,N):
		iBin= i+1
		gen . SetBinContent(iBin, 1000.*ROOT.TMath.Power(iBin,-5) ) 

## construct probability smears
def ConstructSmear():
	for i in range(0,N):
	   for j in range(0,N):
		iBin= i+1
		jBin= j+1
		if (iBin == jBin):smear.SetBinContent(iBin,jBin, 0.7) 
		if (abs(iBin-jBin) == 1):smear.SetBinContent(iBin,jBin, 0.08) 
		if (abs(iBin-jBin) == 2):smear.SetBinContent(iBin,jBin, 0.03) 

## construct response matrix
def ConstructResponse():
	for i in range(0,N):
	   for j in range(0,N):
		iBin= i+1
		jBin= j+1
		resp.SetBinContent(iBin,jBin, smear.GetBinContent(iBin,jBin) * gen.GetBinContent(jBin) ) 

## construct reco
def ConstructReco():
	for i in range(0,N):
	   S=0
	   for j in range(0,N):
		iBin= i+1
		jBin= j+1
		S += resp.GetBinContent(iBin,jBin)			
	   reco.SetBinContent(iBin, S + bkg.GetBinContent(iBin) ) 

def ConstructData():
	data = reco.Clone("data")
	for i in range(0,N):
		iBin= i+1
		data.SetBinContent(iBin, r.Poisson(data.GetBinContent(iBin) ) ) 

ConstructBackground()
ConstructTruth()
ConstructSmear()

ConstructResponse()
ConstructReco()
ConstructData()

print "-> construct Unfolding"
nReg = 5
R = ROOT.RooUnfoldResponse(reco,gen,resp)
u = ROOT.RooUnfoldBayes(R,data,nReg)
u.SetNToys(1000)
h_bayes = u.Hreco(ROOT.RooUnfold.kCovToy)

print "-> construct BootStrap"
b = ROOT.BootStrap()
b.SetNToys(1000)
b.SetUnfoldType(ROOT.BootStrap.kBayes)
b.SetRegParam(nReg)
b.SetUMatrix(reco,gen,resp)
b.SetData(data)

print "-> running BootStrap"
b.run()
g_bootstrap = b.result();

print "-> plotting"

c= ROOT.TCanvas()
gen.SetLineColor(ROOT.kRed)
h_bayes.SetLineColor(ROOT.kGreen+2)
h_bayes.SetFillColor(ROOT.kGreen)

g_bootstrap.SetLineColor(ROOT.kBlue+2)
g_bootstrap.SetMarkerColor(ROOT.kBlue+2)


gen.Draw("HIST")
h_bayes.Draw("PE2 SAME")
g_bootstrap.Draw("PE SAME")
gen.Draw("HIST SAME")

reco.SetLineColor(ROOT.kBlack)
reco.Draw("HIST SAME")
bkg.SetLineColor(ROOT.kGray)
bkg.SetLineStyle(2)
bkg.Draw("HIST SAME")

gen.Print("ALL")
print "BAYES"
h_bayes.Print("ALL")
print "BOOTSTRAP"
g_bootstrap.Print("V")

l = ROOT.TLegend(0.6,.6,.89,.89)
l.AddEntry(gen,"truth","L")
l.AddEntry(h_bayes,"bayes","LF")
l.AddEntry(g_bootstrap,"bayes","PE")
l.Draw()

c2=ROOT.TCanvas()
c2.Divide(2)
c2.cd(1)
smear.Draw("BOX")
c2.cd(2)
resp.Draw("BOX")

raw_input("Looks ok ? ")


