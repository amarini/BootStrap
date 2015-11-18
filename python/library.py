#!/bin/env python

import ROOT
import sys,os
import numpy


def loadBootStrap():
	if hasattr(ROOT,'BootStrap') : return
	print "-> Loading Library"
	ROOT.gSystem.Load("bin/libBootStrap.so") # it will know where it is RooUnfold.so and load it


N=30
NBkg=1.e3
NSig=1e7
r=ROOT.TRandom3(23149)


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
		gen . SetBinError(iBin,  gen.GetBinContent(iBin) *0.01) 
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
		if type==3: ## pretty normal, low eff
			if (iBin == jBin):smear.SetBinContent(iBin,jBin, 0.7*.1) 
			if (abs(iBin-jBin) == 1):smear.SetBinContent(iBin,jBin, 0.01*.1) 
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
		resp . SetBinError(iBin,jBin,  resp.GetBinContent(iBin,jBin) *0.01) 
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
	   reco . SetBinContent(iBin, S + bkg.GetBinContent(iBin) ) 
	   reco . SetBinError(iBin,  reco.GetBinContent(iBin) *0.01) 
	return reco

def ConstructData(reco):
	''' Construct toy data smearing the reco distribution with poisson.'''
	data = reco.Clone("data."+reco.GetName())
	for i in range(0,N):
		iBin= i+1
		c= r.Poisson(data.GetBinContent(iBin) )
		if c<0 : print "* ERROR: negative poisson!"
		data.SetBinContent(iBin, c ) 
	   	data . SetBinError(iBin,  ROOT.TMath.Sqrt( c )) 
	return data

def ConstructFromTree(N=10000,Nbins=100):
	''' Construct a matrix from a tree. This is used to test positive and negative weights. return reco/truth/resp. TODO'''
	## TODO
	reco=ROOT.TH1D("reco","reco",Nbins,0,200.)
	truth=ROOT.TH1D("reco","reco",Nbins,0,200.)
	resp = ROOT.TH2D("resp","reps",Nbins,0,200.,Nbins,0,200.)
	#gen . SetBinContent(iBin, NSig*ROOT.TMath.Power(iBin * 20./ N,-4.9) + NSig/10000. *ROOT.TMath.Power(iBin,-3) )
	
	res = 2 ## 2 GEV ?
	eff = 0.95 #
	totAccepted=0.
	tot=0.
	## number makes this function reasonable
	fPt = ROOT.TF1("f1","1e+6*TMath::Power(x,-3)*TMath::Exp(-50./x)",0,200);
	for i in range(0,N):
		if i%100 == 0:
			print "\r Doing entry: "+str(i)+"/" + str(N) + " : %.1f%%"%(float(i)/N * 100.),
			sys.stdout.flush()
		pt = fPt.GetRandom()	

		if ROOT.gRandom.Uniform(1) < .97:
			ptReco = ROOT.gRandom.Gaus(pt,res)
		else:
			ptReco = ROOT.gRandom.Gaus(pt,res*10) ## non gaus tail

		if ptReco<0:ptReco=0

		accept = ( ROOT.gRandom.Uniform(1)  < eff)

		# 1/3 of the weights are negative -> Flat, we can add a pt dependence
		w = 5000./float(N)  ## expected events : 5000. * (1-.3*2)
		if ROOT.gRandom.Uniform(1)<.40:
			w *= -1

		tot += w
		truth.Fill(pt,w)
		if accept:
			totAccepted += w
			resp.Fill(ptReco,pt,w)
			reco.Fill(ptReco,w)

	print "\n* Done"
	print "Tot=",tot, "eff=",totAccepted/tot
	return (reco,truth,resp)
