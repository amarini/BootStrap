import ROOT 


## as in HIG-14-016 paper
#njets from amcatnlo
njets= ROOT.TH1D("njets","njets", 4,0-.5,4-.5)
njets.SetBinContent(1,17.4)
njets.SetBinContent(2,8.7)
njets.SetBinContent(3,3.1)
njets.SetBinContent(4,1.1)


#https://indico.cern.ch/conferenceDisplay.py?confId=186478
resp=ROOT.TH2D("resp","resp",4,0-.5,4-.5,4,0-.5,4-.5)

# resp  ( TRUTH, MEASURED )
resp.SetBinContent( 1, 1 , .962 )  ## 0 jets
resp.SetBinContent( 2, 2 , .793 ) 
resp.SetBinContent( 3, 3 , .715 ) 
resp.SetBinContent( 4, 4 , .708 ) 

resp.SetBinContent(1,2 ,  .166 ) 
resp.SetBinContent(2,1 ,  .037 ) 

resp.SetBinContent(1,3 ,  .035 ) 
resp.SetBinContent(3,1 ,  .01 ) 
resp.SetBinContent(2,3 ,  .218 ) 
resp.SetBinContent(3,2 ,  .04 ) 

resp.SetBinContent(1,4 ,  .0011 ) 
resp.SetBinContent(4,1 ,  .0 ) 
resp.SetBinContent(2,4 ,  .053 ) 
resp.SetBinContent(4,2 ,  .001 ) 
resp.SetBinContent(3,4 ,  .228 ) 
resp.SetBinContent(4,3 ,  .034 ) 


## normalize in truth
for iTruth in range(1,njets.GetNbinsX()+1):
   sum=0.
   for jMeas in range(1,njets.GetNbinsX()+1):
	sum += resp.GetBinContent(iTruth,jMeas)
   for jMeas in range(1,njets.GetNbinsX()+1):
	c = resp.GetBinContent(iTruth,jMeas)
	resp.SetBinContent(iTruth,jMeas,c/sum)

### PRINT MATRIX
print "----- MATRIX SMEAR ----- "
for iTruth in range(1,njets.GetNbinsX()+1):
   print "Truth=%3d\t"%iTruth,
   for jMeas in range(1,njets.GetNbinsX()+1):
	   print "%5.2f " % (resp.GetBinContent(iTruth,jMeas) * 100.),
   print 
print "------------------------ "


# relative error on the reco-yields. This is enlarged from Poisson by bkg-fit
error=ROOT.TH1D("error","error",4,0-.5,4-.5)
error.SetBinContent(1,.30)
error.SetBinContent(2,.40)
error.SetBinContent(3,.50)
error.SetBinContent(4,.60)

pulls=[]
pulls.append(None)
for iBin in range(1,njets.GetNbinsX()+1): pulls.append(None)

for iBin in range(1,njets.GetNbinsX()+1):
	pulls[iBin]=ROOT.TH1D("pulls_Bin%d"%iBin,"pull Bin %d"%iBin,2000,-5,5)

respfold=resp.Clone("respfold") # fold spectra and prob
respfold.Reset("ACE")

for iTruth in range(1,njets.GetNbinsX()+1):
   for jMeas in range(1,njets.GetNbinsX()+1):
	#c = respfold.GetBinContent(iTruth,jMeas) 
	s = resp.GetBinContent(iTruth,jMeas)
	t = njets.GetBinContent(iTruth)
	cont = t*s
	respfold.SetBinContent(iTruth,jMeas, cont)


truth = respfold.ProjectionX()
measured = respfold.ProjectionY()

print "---- CHECK ---"
for iTruth in range(1,njets.GetNbinsX()+1):
	print "iTruth=",iTruth, " T= ",truth.GetBinContent(iTruth)," == MC=", njets.GetBinContent(iTruth)
print "--------------"

binbybin= truth.Clone("binbybin")
binbybin.Divide(measured)

r=ROOT.TRandom3()

NToys=1000
for iToy in range(0,NToys):
	## 30pb
	for jMeas in range(1,njets.GetNbinsX()+1):
		y = measured.GetBinContent(jMeas)
		sigma = error.GetBinContent(jMeas)*y
		smeared = r.Gaus(y,sigma)
		corr = binbybin.GetBinContent(jMeas) * smeared
		corrError = sigma * corr 
		p = (truth.GetBinContent(jMeas) - corr) / corrError
		# subtract, eff is 1 by constr.  --- ML ---
		#o = y - respfold.GetBinContent(jMeas,jMeas)
		#p = (truth.GetBinContent(jMeas) - ( smeared - o)) /sigma  ## ML
		pulls[jMeas].Fill(p)

g= ROOT.TF1("gaus","TMath::Gaus(x,0,1,1)",-5,5)
canvas=[]
for jMeas in range(1,njets.GetNbinsX()+1):
	c=ROOT.TCanvas("c%d"%jMeas,"c")
	pulls[jMeas].Scale(1./ pulls[jMeas].Integral("width") ) 
	pulls[jMeas].Draw("HIST")
	g.Draw("L SAME")
	canvas.append(c)
raw_input("ok?")
for c in canvas:
	c.SaveAs( "plot/BbB_" + c.GetName() +".pdf")



