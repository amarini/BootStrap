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
