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

## ROOT.gROOT.ProcessLine (\
## 		"struct RGB{ \
## 		float r;\
## 		float g;\
## 		float b;\
## 		void SetRGB(int color) { gROOT->GetColor(color)->GetRGB(r,g,b); }; \
## 		};" )
## from ROOT import RGB

## def rgb(color):
## 	r= float(0.)
## 	g= float(0.)
## 	b= float(0.)
## 
## 	ROOT.gROOT.GetColor(color).GetRGB(r,g,b)
## 	return (r,g,b)

### class color():
### 	def __init__(self):
### 		self.r = 0
### 		self.g = 0
### 		self.b = 0
### 		self.rgb = RGB(0)
### 	def __init__(self,color):
### 		self.rgb = RGB()
### 		self.rgb.SetRGB(color)
### 		self.r = self.rgb.r
### 		self.g = self.rgb.g
### 		self.b = self.rgb.b

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
#b = ROOT.BootStrapMatrix()
b = ROOT.BootStrap()
b.SetUnfoldType(ROOT.BootStrap.kBayes) ## BootStrap
#b.SetUnfoldType(ROOT.BootStrap.kInv) ## BootStrap
b.SetRegParam(nReg) ##BootStrap

b.SetNToys(1000)
b.SetSeed(328956)
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

g_bs2 = b.result(ROOT.BootStrap.kStd,.68)
g_bs2 = ROOT.utils.Shift( g_bs2, -.3,True)

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
gen.SetLineColor(ROOT.kRed)
h_bayes.SetLineColor(ROOT.kGreen+2)
h_bayes.SetMarkerColor(ROOT.kGreen+2)
h_bayes.SetMarkerStyle(29)
h_bayes.SetFillColor(ROOT.kGreen-9)

g_bootstrap.SetLineColor(ROOT.kBlue+2)
g_bootstrap.SetMarkerColor(ROOT.kBlue+2)
g_bootstrap.SetMarkerStyle(20)
g_bootstrap.SetMarkerSize(0.8)

g_bs2.SetLineColor(ROOT.kViolet+1)
g_bs2.SetMarkerColor(ROOT.kViolet+1)
g_bs2.SetMarkerStyle(24)
g_bs2.SetMarkerSize(0.8)

bkg.SetLineColor(ROOT.kGray+1)
bkg.SetLineStyle(2)
bkg.SetFillColor(ROOT.kGray)

#closure.SetLineColor(ROOT.kMagenta)
#closure.SetLineStyle(7)

reco.SetLineColor(ROOT.kBlack)

data.SetMarkerStyle(24)
data.SetMarkerColor(ROOT.kMagenta)

gen.GetYaxis().SetLabelFont(43)
gen.GetYaxis().SetLabelSize(26)
gen.GetXaxis().SetLabelFont(43)
gen.GetXaxis().SetLabelSize(26)

gen.Draw("AXIS")
bkg.Draw("HIST SAME")

h_bayes.Draw("PE2 SAME")
g_bootstrap.Draw("PE SAME")
g_bs2.Draw("PE SAME")
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
l.AddEntry(g_bs2,"bootstrap 2","PE")
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

p1.SetLogy()

p2.cd()

gen_r = ROOT.utils.Ratio(gen,gen)

gen_r.Draw("AXIS")
gen_r.GetYaxis().SetNdivisions(204)
gen_r.GetYaxis().SetRangeUser(-.1,2.3)

h_bayes_r = ROOT.utils.Ratio(h_bayes,gen)
h_bayes_r.Draw("PE2 SAME")
g_bs_r = ROOT.utils.Ratio(g_bootstrap, gen)
g_bs_r.Draw("PE SAME")
g_bs2_r = ROOT.utils.Ratio(g_bs2, gen)
g_bs2_r.Draw("PE SAME")

gen_r.Draw("HIST SAME") 
gen_r.Draw("AXIS SAME")
gen_r.Draw("AXIS X+ Y+ SAME")

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
corr_bayes.Draw("AXIS")
ex1 = ROOT.TExec("ex2","utils::ChangePalette(1);");
ex1.Draw("SAME")
corr_bayes.Draw("COLZ SAME")
corr_bayes.GetZaxis().SetRangeUser(-1,1)
ltx.DrawLatex(.5,.92,"Bayes Correlation Matrix")

c2.cd(2)
corr = b.correlation()
corr.Draw("AXIS")
ex1.Draw("SAME")
corr.Draw("COLZ SAME")
corr.GetZaxis().SetRangeUser(-1,1)
ltx.DrawLatex(.5,.92,"Bootstrap Correlation Matrix")

### SHOW SOME TOY DISTRIBUTION
c3 = ROOT.TCanvas("c3","c3",1200,10,800,600)
c3.Divide(2)
c3.cd(1)
bin=4
toyDistr=b.GetToyDistribution(bin)
toyDistr.Draw("HIST")
y	= g_bootstrap.GetY()[bin-1]
ylow	= y-g_bootstrap.GetEYlow()[bin-1]
yhigh	= y+g_bootstrap.GetEYhigh()[bin-1]

up = toyDistr.GetMaximum()*.3
g_median=ROOT.TGraph()
g_median.SetName("g_median")
g_median.SetPoint(0,y,0)
g_median.SetPoint(1,y,up)

g_low=ROOT.TGraph()
g_low.SetName("g_low")
g_low.SetPoint(0,ylow,0)
g_low.SetPoint(1,ylow,up)

g_high=ROOT.TGraph()
g_high.SetName("g_high")
g_high.SetPoint(0,yhigh,0)
g_high.SetPoint(1,yhigh,up)

g_median.SetLineColor(ROOT.kRed)
g_low.SetLineColor(ROOT.kBlue)
g_high.SetLineColor(ROOT.kBlue)

g_median.SetLineWidth(3)
g_low.SetLineWidth(2)
g_high.SetLineWidth(2)
g_median.Draw("L SAME")
g_low.Draw("L SAME")
g_high.Draw("L SAME")

ltx.DrawLatex(.5,.92,"Toy Distribution for bin %d"%bin)

ROOT.utils.ChangePalette(2)
c3.cd(2)
toyDistr2=b.GetToyDistribution(bin,bin+1)
toyDistr2.Draw("AXIS")
ex2 = ROOT.TExec("ex2","utils::ChangePalette(2);");
ex2.Draw()
toyDistr2.Draw("COLZ SAME")
ltx.DrawLatex(.5,.92,"Toy Correlation for bins %d-%d"%(bin,bin+1))

c4 =ROOT.TCanvas("c4","c4",400,400)
##ex3 = ROOT.TExec("ex2","utils::ChangePalette(3);");
smear.Draw("AXIS")
###ex3.Draw()
smear.Draw("BOX SAME")
smear.GetYaxis().SetTitle("truth")
smear.GetXaxis().SetTitle("measured")
smear.GetYaxis().SetTitleOffset(1.2)
smear.Draw("AXIS SAME")
smear.Draw("AXIS X+ Y+ SAME")

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

c.Update()
c2.Update()
raw_input("Looks ok ? ")


