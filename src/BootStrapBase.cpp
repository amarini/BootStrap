#include "interface/BootStrapBase.hpp"
#include "TRandom3.h"
#include <iostream>

#define VERBOSE 0

// --- Constructor 
BootStrapBase::BootStrapBase()
{
	Ntoys_ = 100;
	SumW2_ = false;
	data_ = NULL;
	seed_ = 123456;
	unf_ = NULL;
	fold_ = NULL;
}

BootStrapBase::BootStrapBase(int ntoys) : BootStrapBase()
{
	SetNToys(ntoys);
}

// --- Destructor
BootStrapBase::~BootStrapBase(){
	destroyPointer(data_);	
	clearBootstrap();
}

// ---
void BootStrapBase::SetNToys(int ntoys)
{
	Ntoys_ = ntoys;
}

void BootStrapBase::SetData(TH1D* data)
{
	setPointer(data,data_);	
}

void BootStrapBase::SetSumW2(bool sumw2)
{
	SumW2_= sumw2;
}


TH1D* BootStrapBase::releaseData(){
	return releasePointer(data_);	
}


TH1D* BootStrapBase::bootStrap(){
	if (VERBOSE >0 ) cout<<"[BootStrapBase]::[bootStrap]::[DEBUG] checking if pointers are not null "<<endl;
	// create one bootstrap instance
	if (r_ == NULL) r_ = new TRandom3(seed_);
	// unfold data
	if(unf_==NULL) unf_ = Unfold(data_);
	//  fold back into an observation
	if(fold_==NULL) fold_ = Fold(unf_);
	//
	TH1D * toy = (TH1D*)fold_ -> Clone("toy_");
	toy->Reset("ACE");
	// -- smear
	if (VERBOSE >0 ) cout<<"[BootStrapBase]::[bootStrap]::[DEBUG] smear "<<endl;
	for(int i=0;i<= toy->GetNbinsX() +1 ; ++i)
	{
		double c = toy->GetBinContent(i);
		double e = toy->GetBinError(i);
		double c2;

		if( not SumW2_ and c<=0 ) continue;

		if (SumW2_) c2=r_->Gaus(c,e);
		else c2 = r_->Poisson(c); 

		toy->SetBinContent(i,c2);
	}
	// -- unfold toy
	if (VERBOSE >0 ) cout<<"[BootStrapBase]::[bootStrap]::[DEBUG] unfold "<<endl;
	TH1D * toy2 = Unfold(toy);
	toy2->SetName("toy2_");

	destroyPointer(toy);	

	return toy2;
}

void BootStrapBase::run(){
	if (VERBOSE >0 ) cout<<"[BootStrapBase]::[run]::[DEBUG] clearBootstrap "<<endl;
	clearBootstrap();
	if (VERBOSE >0 ) cout<<"[BootStrapBase]::[run]::[DEBUG] destroyPointers "<<endl;
	destroyPointer(unf_);
	destroyPointer(fold_);
	for(int iToy=0;iToy<Ntoys_;++iToy)
	{
		if (VERBOSE >0 ) cout<<"[BootStrapBase]::[run]::[DEBUG] running Toy "<< iToy <<endl;
		TH1D *toy = bootStrap();	
		toy->SetName( Form("toy_%d",iToy) ) ;
		bootstrap_ . push_back( toy );
	}
}

#include "interface/stat.hpp"
TGraphAsymmErrors* BootStrapBase::result( ResultType type,float Q)
{
	if (bootstrap_.size() != size_t(Ntoys_ ) ) run() ;

	TGraphAsymmErrors *R = new TGraphAsymmErrors();

	switch (type)
		{
		case kStd:
			R->SetName(Form("bootstrap_Ntoys%d_kStd",Ntoys_));
			break;
		case kMin:
			R->SetName(Form("bootstrap_Ntoys%d_kMin",Ntoys_));
			break;
		case kMean:
			R->SetName(Form("bootstrap_Ntoys%d_kMean",Ntoys_));
			break;
		case kMedian:
			R->SetName(Form("bootstrap_Ntoys%d_kMedian",Ntoys_));
			break;
		}
	
	// visible bins
	for(int iBin=1 ;iBin<= unf_->GetNbinsX() ;++iBin)
	{
		// get all the values for this bin
		vector<float> values;
		for(int i=0;i<Ntoys_;++i)
			values.push_back( bootstrap_[i]->GetBinContent(iBin) );
		float mean;
		pair<float,float> err;
		switch (type)
		{
		case kStd:
			mean=unf_->GetBinContent(iBin);
			STAT::ConfidenceIntervalAround(values, mean, err,Q);
			break;
		case kMin:
			mean=unf_->GetBinContent(iBin);
			STAT::ConfidenceInterval(values, err,Q);
			break;
		case kMean:
			mean=STAT::mean(values);
			STAT::ConfidenceInterval(values, err,Q);
			break;
		case kMedian:
			mean=STAT::median(values);
			STAT::ConfidenceInterval(values, err,Q);
			break;
		}

		// error are in absolute values, will scale to the mean point
		//
		if ( mean< err.first) cout<<"Error mean point is outside the error bands"<<endl;
		if ( mean> err.second) cout<<"Error mean point is outside the error bands"<<endl;
		float elow = mean-err.first;
		float ehigh = err.second -mean;

		int npoint= R->GetN();
		R->SetPoint( npoint, unf_->GetBinCenter( iBin), mean);
		R->SetPointError(npoint, 0,0, elow,ehigh);

	}
	return R;	
}

// GetCorrelation matrix
TH2D* BootStrapBase::correlation(){
	if (bootstrap_.size() != size_t(Ntoys_ ) ) run() ;

	int nBins= unf_->GetNbinsX();
	TH2D*R = new TH2D("corr","corr",nBins,0-.5,nBins-.5,nBins,0-.5,nBins-.5);

	for(int iBin=1; iBin<= nBins;++iBin)
	for(int jBin=iBin; jBin<= nBins;++jBin)
	{
		if(iBin==jBin) { R->SetBinContent(iBin,jBin,1);	continue; }
		vector<float> valuesI;
		vector<float> valuesJ;
		for(int i=0;i<Ntoys_;++i)
		{
			valuesI.push_back( bootstrap_[i]->GetBinContent(iBin) );
			valuesJ.push_back( bootstrap_[i]->GetBinContent(jBin) );
		}
		float corr = STAT::corrPearson( valuesI,valuesJ);
		R->SetBinContent(iBin,jBin,corr);
		R->SetBinContent(jBin,iBin,corr);
	}
	return R;

}


