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
for i in range(0,N):
	iBin= i+1
	gen . SetBinContent(iBin, 1000.*ROOT.TMath.Power(iBin,-5) ) 
	bkg . SetBinContent(iBin, 10.*  ROOT.TMath.Cos( float(iBin)/N * 3.14) ) 

## construct probability smears
for i in range(0,N):
   for j in range(0,N):
	iBin= i+1
	jBin= j+1
	if (iBin== jBin):resp.SetBinContent(iBin,jBin, 0.7) 
	if (abs(iBin-jBin) == 1):resp.SetBinContent(iBin,jBin, 0.2) 
	if (abs(iBin-jBin) == 2):resp.SetBinContent(iBin,jBin, 0.1) 

## construct response matrix
for i in range(0,N):
   for j in range(0,N):
	iBin= i+1
	jBin= j+1
	resp.SetBinContent(iBin,jBin, smear.GetBinContent(iBin,jBin) * gen.GetBinContent(jBin) ) 

## construct reco
for i in range(0,N):
   S=0
   for j in range(0,N):
	iBin= i+1
	jBin= j+1
	S += resp.GetBinContent(iBin,jBin)			
   reco.SetBinContent(iBin, S + bkg.GetBinContent(iBin) ) 

data = reco.Clone("data")
for i in range(0,N):
	iBin= i+1
	data.SetBinContent(iBin, r.Poisson(data.GetBinContent(iBin) ) ) 


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

l = ROOT.TLegend(0.6,.6,.89,.89)
l.AddEntry(gen,"truth","L")
l.AddEntry(h_bayes,"bayes","LF")
l.AddEntry(g_bootstrap,"bayes","PE")
l.Draw()

raw_input("Looks ok ? ")


