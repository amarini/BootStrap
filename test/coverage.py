#!/bin/env python

import ROOT
import sys,os
import numpy

print "-> inserting in path cwd"
sys.path.insert(0,os.getcwd())
print "-> inserting in path cwd/python"
sys.path.insert(0,os.getcwd()+'/python')

from library import *

import threading
import time
import random

class bsThread(threading.Thread):
	def __init__(self, b):
		## use the copy Constructor
		threading.Thread.__init__(self)
		self.b=ROOT.BootStrap(b)	
		self.chi2=-999
		self.chi2_u=-999
		self.ndf=-999
		self.b.SetSeed( random.randint(1, 10000000) )

	def run(self):
        	self.b.run()
        	g_bs = self.b.result(ROOT.BootStrap.kMedian,.68)
        	corr = self.b.correlation();
        	self.chi2   = ROOT.STAT.Chi2( g_bs, gen, corr ) ;
        	self.chi2_u = ROOT.STAT.Chi2( g_bs, gen ) ;
		self.ndf= g_bs.GetN()
        	corr.Delete()
        	g_bs.Delete()


#library for coverage studies

loadBootStrap()

print "-> Creating Matrixes"

####################
Ntest=1000

bkg  = ConstructBackground()
gen  = ConstructTruth()
smear= ConstructSmear(2) ## 1 ~ diagonal, 2 ~ off-diagonal

resp = ConstructResponse(gen,smear)
reco = ConstructReco(bkg,resp)

nReg = 5
b = ROOT.BootStrap()
b.SetUnfoldType(ROOT.BootStrap.kBayes) ## BootStrap
b.SetToyType(ROOT.BootStrap.kIterBias) ## ToyType
b.SetRegParam(nReg) ##BootStrap
b.SetNToys(1000)
b.SetSeed(328956)
b.SetUMatrix(reco,gen,resp)

b.SetVerbose(0);

h1  =ROOT.TH1D("chi2","chi2;#chi^{2} (%d ndof);Events"%reco.GetNbinsX(),120,0,60)
h1_u=ROOT.TH1D("chi2_uncorr","chi2;#chi^{2};Events",120,0,60)

h2  =ROOT.TH1D("prob","prob;p-value;Events",50,0,1.)
h2_u=ROOT.TH1D("prob_uncorr","prob;p-value;Events",50,0,1.)

h3  =ROOT.TH1D("prob2","prob;p-value (2) ;Events",50,0,1.)
h3_u=ROOT.TH1D("prob2_uncorr","prob;p-value (2) ;Events",50,0,1.)

b.info()
threads=[]


# coverage bin/per/bin, uncorrelated
nIn_u  = 0;
nTot_u = 0;

## coverage corr
nIn  = 0 
nTot = 0 

for i in range(0,Ntest):
	print "\r Preparing " + str(i+1)+"/"+str(Ntest),
	sys.stdout.flush()
	data = ConstructData(reco)
	b.SetData( data.Clone('bootstrap_data') )
	##t = bsThread( b )
	##threads.append(t)
	## here I should have the thread -- and this should be collected later
	b.run()
	g_bs = b.result(ROOT.BootStrap.kMedian,.68)
	corr = b.correlation();
	chi2=ROOT.STAT.Chi2( g_bs, gen, corr ) ;
	chi2_u=ROOT.STAT.Chi2( g_bs, gen ) ;
	ndf = g_bs.GetN()
	## coverage 1
	for i in range(0, ndf): 
		nTot_u += 1
		g  =gen.GetBinContent(i+1)
		c2 =g_bs.GetY()[i]
		e2h=g_bs.GetEYhigh()[i]
		e2l=g_bs.GetEYlow()[i]
		if (g >= c2 and g < c2+e2h) or ( g<=c2 and g> c2-e2l): nIn_u +=1
	## coverage 2

	nTot += 1
	if chi2 < 1: nIn += 1

	corr.Delete()
	g_bs.Delete()
	data.Delete()
	h1.Fill(chi2);
	h1_u.Fill(chi2_u)
	h2.Fill(ROOT.TMath.Prob(chi2,ndf )); ## NDF chi2 -> uniform p-value
	h2_u.Fill(ROOT.TMath.Prob(chi2_u,ndf )); ## NDF chi2 -> uniform p-value

	h3.Fill(ROOT.TMath.Prob(chi2,ndf/2 )); ## try to find some effective ndf
	h3_u.Fill(ROOT.TMath.Prob(chi2_u,ndf/2 )); ## 

## print "DONE"
## print "threads = ", len(threads)
## for i,t in enumerate(threads):
## 
## 	print "\r Starting " + str(i+1)+"/"+str(Ntest),
## 	sys.stdout.flush()
## 
## 	while threading.activeCount() > 5: 
## 		print "\n sleep (", threading.activeCount()," processes) ... ",
## 		sys.stdout.flush()
## 		time.sleep(3)
## 		print "wake up"
## 
## 	t.start()
## 
## print "DONE"
## for i,t in enumerate(threads):
## 	print "\r Processing " + str(i+1)+"/"+str(Ntest),
## 	sys.stdout.flush()
## 	t.join()
## 
## 	chi2   = t.chi2
## 	chi2_u = t.chi2_u
## 	ndf     = t.ndf
## 
## 	h1.Fill(chi2);
## 	h1_u.Fill(chi2_u)
## 	h2.Fill(ROOT.TMath.Prob(chi2,ndf )); ## NDF chi2 -> uniform p-value
## 	h2_u.Fill(ROOT.TMath.Prob(chi2_u,ndf )); ## NDF chi2 -> uniform p-value
## print "DONE"

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
c=ROOT.TCanvas()
c.Divide(2)

h1_u.SetLineColor(ROOT.kRed)
h2_u.SetLineColor(ROOT.kRed)

h3  .SetLineColor(ROOT.kCyan)
h3  .SetLineStyle(7)
h3_u.SetLineColor(ROOT.kOrange)
h3_u.SetLineStyle(2)

p1=c.cd(1)
p1.SetTopMargin(0.05)
p1.SetRightMargin(0.05)
h1.Draw("HIST")
h1_u.Draw("HIST SAME")

p2= c.cd(2)
p2.SetTopMargin(0.05)
p2.SetRightMargin(0.05)

h2.Draw("HIST")
h2_u.Draw("HIST SAME")

h3.Draw("HIST SAME")
h3_u.Draw("HIST SAME")

l = ROOT.TLegend(.30,.65,.80,.90)
l.AddEntry(h1  ,"cov. matrix")
l.AddEntry(h1_u,"diagonal only")
l.AddEntry(h3  ,"diagonal only (ndf/2)")
l.AddEntry(h3_u,"diagonal only (ndf/2)")
l.Draw()

print " --------- COVERAGE ------------"
print " Uncorrelated:", nIn_u, "/",nTot_u, "=", float(nIn_u) / nTot_u
print "   Correlated:", nIn, "/",nTot, "=", float(nIn) / nTot

raw_input("Looks ok ? ")



##
