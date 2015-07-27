#!/bin/env python

import ROOT
import sys,os
import numpy

print "-> inserting in path cwd"
sys.path.insert(0,os.getcwd())
print "-> inserting in path cwd/python"
sys.path.insert(0,os.getcwd()+'/python')

from library import *

from threading import Thread
import time

class bsThread(Thread):
	def __init__(self, b):
		## use the copy Constructor
		self.b=ROOT.BootStrap(b)	
		self.chi2=-999
	def run(self):
        	self.b.run()
        	g_bs = b.result(ROOT.BootStrap.kMedian,.68)
        	corr = b.correlation();
        	self.chi2=ROOT.STAT.Chi2( g_bs, gen, corr ) ;
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
b.SetRegParam(nReg) ##BootStrap
b.SetNToys(1000)
b.SetSeed(328956)
b.SetUMatrix(reco,gen,resp)

b.SetVerbose(0);
h=ROOT.TH1D("chi2","chi2",1000,0,20)

b.info()
threads=[]
for i in range(0,Ntest):
	print "\r Doing test " + str(i)+"/"+str(Ntest),
	data = ConstructData(reco)
	b.SetData( data.Clone('bootstrap_data') )
	#t = bsThread( b )
	#threads.append(t)
	#t.start()
	## here I should have the thread -- and this should be collected later
	b.run()
	g_bs = b.result(ROOT.BootStrap.kMedian,.68)
	corr = b.correlation();
	chi2=ROOT.STAT.Chi2( g_bs, gen, corr ) ;
	corr.Delete()
	g_bs.Delete()
	data.Delete()
	h.Fill(chi2);

c=ROOT.TCanvas()

h.Draw("HIST")

raw_input("Looks ok ? ")



##
