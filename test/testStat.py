import ROOT
import sys,os
import numpy

print "-> inserting in path cwd"
sys.path.insert(0,os.getcwd())
print "-> inserting in path cwd/python"
sys.path.insert(0,os.getcwd()+'/python')

from library import *


loadBootStrap()

#myf1 = ROOT.TF1("myfunc","TMath::Exp(-x/.01)",0,1)
myf1 = ROOT.TF1("myfunc"," (x<.01) ? 1000:TMath::Exp(-x/.01)",0,1)
h=ROOT.TH1D("h","h",100,0,1)
h.FillRandom("myfunc",20000)

h2=ROOT.STAT.Normalize(h)

h2.Draw("PE")
raw_input("ok?")
