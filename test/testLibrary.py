import ROOT
import sys,os
import numpy

print "-> inserting in path cwd"
sys.path.insert(0,os.getcwd())
print "-> inserting in path cwd/python"
sys.path.insert(0,os.getcwd()+'/python')

from library import *


loadBootStrap()
#print "-> Loading Library"
#ROOT.gSystem.Load("bin/libBootStrap.so") # it will know where it is RooUnfold.so and load it

print "-> Creating Matrixes"

####################
#NBkg *= 10
bkg  = ConstructBackground()
gen  = ConstructTruth(1)
smear= ConstructSmear(2) ## 1 ~ diagonal, 2 ~ off-diagonal, 3 low eff

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
nReg = 15
print "\t\tBayes, nReg=",nReg

R = ROOT.RooUnfoldResponse(reco,gen,resp)
#u = ROOT.RooUnfoldBayes(R,data2,nReg)
u = ROOT.RooUnfoldSvd(R,data2,nReg)
u.SetVerbose(-1)
u.SetNToys(1000)
#u.IncludeSystematics(2)  ### only MATRIX
#h_bayes = u.Hreco(ROOT.RooUnfold.kCovToy)
h_svd = u.Hreco(ROOT.RooUnfold.kCovToy)

#uInv = ROOT.RooUnfoldInvert(R,reco)
#closure = uInv.Hreco()

print "-> construct BootStrap"
#b = ROOT.BootStrapMatrix()
b = ROOT.BootStrap()
#b.SetUnfoldType(ROOT.BootStrap.kBayes) ## BootStrap
#b.SetUnfoldType(ROOT.BootStrap.kInv) ## BootStrap
b.SetUnfoldType(ROOT.BootStrap.kTUnfoldDensity) ## Density
#b.SetRegParam(nReg) ##BootStrap
b.SetRegParam(-1) ##BootStrap

b.SetNToys(1000)
b.SetSeed(328956)
b.SetUMatrix(reco,gen,resp)
b.SetData( data2.Clone('bootstrap_data') )

## As in RooUnfold for my understanding
b.SetToyType(ROOT.BootStrap.kToy)
#b.SetToyType(ROOT.BootStrap.kIterBias)


#b.SetSumW2(); ## MATRIX
#b.SetToyType(ROOT.BootStrap.kMatrix) ## MATRIX

print "-> running BootStrap"
## ## kStd/kMin/kMedian/kMean

b.run()
#g_bootstrap = b.result(ROOT.BootStrap.kMin,.68)
g_bootstrap = b.result(ROOT.BootStrap.kMedian,.68)
g_bootstrap = ROOT.utils.Shift( g_bootstrap, 0.3, True)


b2 = ROOT.BootStrap(b);

b2.SetToyType(ROOT.BootStrap.kBootstrap)
b2.run()

g_bs2 = b2.result(ROOT.BootStrap.kMedian,.68) ## 
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
h_svd.SetLineColor(ROOT.kGreen+2)
h_svd.SetMarkerColor(ROOT.kGreen+2)
h_svd.SetMarkerStyle(29)
h_svd.SetFillColor(ROOT.kGreen-9)

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

h_svd.Draw("PE2 SAME")
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
#l.AddEntry(h_bayes,"bayes","LF")
l.AddEntry(h_svd,"svd","LF")
l.AddEntry(g_bootstrap,"Toys","PE")
l.AddEntry(g_bs2,"bootstrap ","PE")
l.AddEntry(reco,"reco","L")
l.AddEntry(bkg,"bkg","LF")
l.AddEntry(data,"data","P")
l.Draw()

print "-> Coverage Test"
tot=0
bootstrap= 0
bayes=0
for i in range(1,h_svd.GetNbinsX()+1):
	g = gen.GetBinContent(i)

	c1 = h_svd.GetBinContent(i)
	e1 = h_svd.GetBinError(i)

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

h_bayes_r = ROOT.utils.Ratio(h_svd,gen)
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
	if err(i)*err(j) >0:
		corr_bayes.SetBinContent(i+1,j+1, cov(i,j) / (err(i)*err(j) ) )
	else:
		corr_bayes.SetBinContent(i+1,j+1, 0 ) 

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


