#include "interface/BootStrapBase.hpp"
#include "interface/stat.hpp"
#include "TRandom3.h"
#include <iostream>
#include <ctime>

#define VERBOSE 0

template<>
void BootStrapBase::destroyPointer<TObject*>(TObject* &ptr){
	if (ptr != NULL) ptr->Delete();
	ptr=NULL;
}

// --- Constructor 
BootStrapBase::BootStrapBase()
{
	Ntoys_ = 100;
	SumW2_ = false;
	data_ = NULL;
	seed_ = 12345678;
	unf_ = NULL;
	fold_ = NULL;
	r_ = NULL;
	verbose_=1;
	type_= kBootstrap;
	// iterative bias
	bias_ = NULL;
	Nib_ = 5;
	Ntoysib_ = 0;
}

// copy constructor
BootStrapBase::BootStrapBase( BootStrapBase &x)
{
	if(VERBOSE>1)cout<<"[BootStrapBase]::[BootStrapBase]::[copy constructor]"<<endl;
	Ntoys_ = x.Ntoys_;
	SumW2_ = x.SumW2_;
	data_ = (x.data_) ? (TH1D*)x.data_->Clone(Form("%s_%u",x.data_ ->GetName(), (unsigned)time(NULL))) : NULL;
	seed_ = (unsigned)time(NULL);
	unf_ = NULL;
	fold_ = NULL;
	r_ = NULL;
	verbose_= x.verbose_;
	type_= x.type_;
	// iterative bias
	bias_ = NULL;
	Nib_ = x.Nib_;
	Ntoysib_ = x.Ntoysib_;
}

// --- Destructor
BootStrapBase::~BootStrapBase(){
	destroyPointer(data_);	
	destroyPointer(r_);	
	destroyPointer(unf_);	
	destroyPointer(fold_);	
	destroyPointer(bias_);
	clearBootstrap();
	clearIbtoys();
}

void BootStrapBase::clearIbtoys(){
	for( TH1D*ptr : ibtoys_ ) destroyPointer( ptr) ; 
	ibtoys_.clear(); 
}

// ---
void BootStrapBase::SetNToys(int ntoys)
{
	Ntoys_ = ntoys;
}

void BootStrapBase::SetData(TH1D* data)
{
	if(VERBOSE>0 )cout<<"[BootStrapBase]::[SetData]::[1]"<<endl;
	setPointer(data,data_);	
}

void BootStrapBase::SetSumW2(bool sumw2)
{
	SumW2_= sumw2;
}


TH1D* BootStrapBase::releaseData(){
	return releasePointer(data_);	
}

void BootStrapBase::Smear(TH1*toy)
{
	if (r_ == NULL) r_ = new TRandom3(seed_);
	for(int i=0;i<= toy->GetNbinsX() +1 ; ++i)
	{
		double c = toy->GetBinContent(i); //  fold_
		double e = toy->GetBinError(i);
		if( not SumW2_ and c<=0 ) continue;
		double c2;
		if (SumW2_) { c2=r_->Gaus(c,e);}
		else { c2 = r_->Poisson(c); }

		toy->SetBinContent(i,c2);
		toy->SetBinError(i,0);
	
	}
	return;
}

void BootStrapBase::Smear(TH2*toy)
{
	if (r_ == NULL) r_ = new TRandom3(seed_);
	for(int j=0;j<= toy->GetNbinsY() +1 ; ++j)
	for(int i=0;i<= toy->GetNbinsX() +1 ; ++i)
	{
		double c = toy->GetBinContent(i,j); //  fold_
		double e = toy->GetBinError(i,j);
		if( not SumW2_ and c<=0 ) continue;
		double c2;
		if (SumW2_) { c2=r_->Gaus(c,e);}
		else { c2 = r_->Poisson(c); }

		toy->SetBinContent(i,j,c2);
		toy->SetBinError(i,j,0);
	
	}
	return;
}

TH1D* BootStrapBase::iterativeBias(bool toy){
	if (VERBOSE >0 ) cout<<"[BootStrapBase]::[iterativeBias]::[DEBUG] iterative Bias"<<endl;

	destroyPointer(bias_);
	TH1D *data;

	// toys are smeared
	if (toy ) 
		{
		if( unf_ == NULL ) { 
				iterativeBias(0);
				destroyPointer(fold_);
				}

		data = (TH1D*)data_->Clone( Form("%s_ib_toy",data_->GetName() ));
		Smear(data);
		}
	else 	{
		if( unf_ == NULL ) {
				unf_ = Unfold(data_);
				destroyPointer(fold_);
				}
		data = (TH1D*)data_->Clone( Form("%s_ib",data_->GetName() ));
		}
	TH1D * saveUnfold = (TH1D*)unf_->Clone("unfolded_data_for_ib");

	if (VERBOSE >0 ) cout<<"[BootStrapBase]::[iterativeBias]::[DEBUG] data are"<<data_<<" " <<data_->GetName()<<endl;

	if (r_ == NULL) r_ = new TRandom3(seed_);
	if(fold_==NULL) fold_ = Fold(unf_);

	if (Ntoysib_<=0 ) Ntoysib_ = Ntoys_ / 10; // I want only the mean point out of this toys, boost t

	int Ntot= Nib_ * Ntoysib_ ;
	int iRun=1;

	if (VERBOSE >0 ) cout<<"[BootStrapBase]::[iterativeBias]::[DEBUG] Running IB iterations"<<endl;
	for( int iIB = 0 ;iIB < Nib_ ;++iIB)
	{
		if (VERBOSE >1 ) cout<<"[BootStrapBase]::[iterativeBias]::[DEBUG] IB "<<iIB<<endl;
		// clear the array of toys
		clearIbtoys();
		// sampling of the bias
		for(int iToy = 0;iToy<Ntoysib_; ++iToy)
		{
			if (verbose_>0 or VERBOSE >0) {
				cout<<"\r * " << iRun++ <<"/" << Ntot << flush;
				}
			TH1D *toy = bootStrap(); // will apply a smearing
			ibtoys_ . push_back( toy );
		}

		// compute bias
		destroyPointer(bias_);
		bias_ = (TH1D*)unf_->Clone(Form("bias_%d",iIB));
		bias_->Reset("ACE");
		
		if (VERBOSE >1 ) cout<<"[BootStrapBase]::[iterativeBias]::[DEBUG] Compute Bias "<<iIB<<endl;

		for(int iBin=0 ;iBin<= bias_->GetNbinsX()+1 ;++iBin)
		{
			 vector<float> values; 	
			for(int iToy=0;iToy<Ntoysib_;++iToy) values.push_back( ibtoys_[iToy]->GetBinContent(iBin) );
			bias_->SetBinContent(iBin, - STAT::median(values) + unf_->GetBinContent(iBin)); // OR mean ?!?
			//bias_->SetBinContent(iBin, - STAT::mean(values) + unf_->GetBinContent(iBin)); // OR mean ?!?
		}
		if(verbose_>1 or VERBOSE >0){
			// display bias
			cout <<endl<<"--- BIAS iteration "<<iIB<<" ---"<<endl;
			for(int iBin=0;iBin<bias_->GetNbinsX()+1 ;++iBin)
				cout << bias_->GetBinContent(iBin)<<" ";
			cout <<endl<<"--------------------------------"<<endl;
		}

		// destroy unf_ pointer
		destroyPointer(unf_);
		destroyPointer(fold_);
		// update it
		if (VERBOSE >1 ) cout<<"[BootStrapBase]::[iterativeBias]::[DEBUG] Update unfolding "<<endl;
		unf_ = (TH1D*)saveUnfold->Clone(Form("biasCorrected_%d",iIB));
		for(int iBin=0;iBin<=bias_->GetNbinsX()+1; ++iBin)
			{
			double c = unf_->GetBinContent(iBin);
			double b = bias_->GetBinContent(iBin);
			unf_ ->SetBinContent(iBin, c  - b) ;
			}
		fold_= Fold(unf_);

	}
	
	if (VERBOSE >1 ) cout<<"[BootStrapBase]::[iterativeBias]::[DEBUG] Compute Estimate "<<endl;
	// construct the estimate from the ibtoys
	TH1D *res= (TH1D*) unf_->Clone( Form("%s_ibres",data->GetName()) ) ;
	res->Reset("ACE");
	for(int iBin=0;iBin<=unf_->GetNbinsX()+1 ;++iBin)
	{
		vector<float> values; 	
		for(int iToy=0;iToy<Ntoysib_;++iToy) values.push_back( ibtoys_[iToy]->GetBinContent(iBin) );
		res->SetBinContent(iBin, STAT::median(values) ); //  OR mean ?!?
	}

	if (verbose_>0 )  cout <<"\r * DONE      "<<endl;

	if (toy)
		{
		swapPointers(saveUnfold,unf_ ); // put back the old result in unf_ , local data will be destroyed
		}
	else 
		{
		destroyPointer(unf_);
		swapPointers(res,unf_ ) ; // put result in unfolding and return null
		}

	destroyPointer(saveUnfold);
	destroyPointer(data);

	if (VERBOSE >0 ) cout<<"[BootStrapBase]::[directToy]::[DEBUG] DONE "<<endl;
	return res;
}

TH1D* BootStrapBase::directToy(){
	// smear data directly, and unfold the toy
	if (VERBOSE >0 ) cout<<"[BootStrapBase]::[directToy]::[DEBUG] checking if pointers are not null "<<endl;
	if (r_ == NULL) r_ = new TRandom3(seed_);
	// unfold data
	if(unf_==NULL) unf_ = Unfold(data_);
	TH1D * toy1 = (TH1D*)data_ -> Clone("toy_");
	
	Smear(toy1);	

	TH1D * toy2 = Unfold(toy1);
	destroyPointer(toy1);	

	return toy2;
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
	// -- smear
	if (VERBOSE >0 ) cout<<"[BootStrapBase]::[bootStrap]::[DEBUG] smear fold_="<< fold_<<" != NULL"<<" r"<<r_ << "toy="<<toy<<endl;

	//
	Smear(toy);
	// -- unfold toy
	if (VERBOSE >0 ) cout<<"[BootStrapBase]::[bootStrap]::[DEBUG] unfold "<<endl;
	TH1D * toy2 = Unfold(toy);
	toy2->SetName("toy2_");

	destroyPointer(toy);	

	return toy2;
}

void BootStrapBase::run(){

	if(verbose_>0)
	{
		info();
	}


	if (VERBOSE >0 ) cout<<"[BootStrapBase]::[run]::[DEBUG] clearBootstrap "<<endl;
	clearBootstrap();
	if (VERBOSE>0) cout<<"[BootStrapBase]::[run]::[DEBUG] correct"<<endl;
	correct();
	if (VERBOSE >0 ) cout<<"[BootStrapBase]::[run]::[DEBUG] destroyPointers "<<endl;
	destroyPointer(unf_);
	destroyPointer(fold_);

	if (type_ == kIterBias) 
		{ // the toy generation is a bit more complicate
		iterativeBias(0); // return NULL, and put result in unf_
		}

	int iRun=1;

	for(int iToy=0;iToy<Ntoys_;++iToy)
	{
		if (VERBOSE >0 ) cout<<"[BootStrapBase]::[run]::[DEBUG] running Toy "<< iToy <<endl;
		if (verbose_>0 ) {
				cout << "\r * " << iRun++ << "/" << Ntoys_ << flush;
				}
		TH1D*toy = NULL;
		switch (type_) 
		{
		case kBootstrap: toy = bootStrap();break;
		case kToy: toy = directToy();break;
		case kIterBias: toy = iterativeBias(1); 
				if (verbose_ >0 ) cout <<endl;
				break; // should not arrive here. It's there to complete the switch over enum.
		case kMatrix: toy = matrixSmear(); break;
		}

		toy->SetName( Form("toy_%d",iToy) ) ;
		bootstrap_ . push_back( toy );
	}
	if (verbose_>0 ) cout <<"\r * DONE                         "<<endl;
}



TGraphAsymmErrors* BootStrapBase::result( ResultType type,float Q)
{
	if (bootstrap_.size() == size_t(0) )
		{ 
			cout <<"[BootStrapBase]::[result]::[WARNING] No Bootstrap toy. Re-running bootstrap"<<endl;
			run() ;
		}
	if (bootstrap_.size() != size_t(Ntoys_ ) )
		{
			cout <<"[BootStrapBase]::[result]::[WARNING] Avalaible toys are:"<< bootstrap_.size()<<" instead of "<<Ntoys_<<endl;
		}

	TGraphAsymmErrors *R = new TGraphAsymmErrors();

	switch (type)
		{
		case kStd:
			R->SetName(Form("bootstrap_Ntoys%d_kStd",Ntoys_));
			break;
		case kRms:
			R->SetName(Form("bootstrap_Ntoys%d_kRms",Ntoys_));
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
			{
			mean=unf_->GetBinContent(iBin);
			STAT::ConfidenceIntervalAround(values, mean, err,Q);
			break;
			}
		case kRms:
			{
			mean=unf_->GetBinContent(iBin);
			float rms = STAT::rms(values);
			err.first =  mean - (rms/2.0) ;
			err.second = mean + (rms/2.0) ;
			break;
			}
		case kMin:
			{
			mean=unf_->GetBinContent(iBin);
			STAT::ConfidenceInterval(values, err,Q);
			break;
			}
		case kMean:
			{
			mean=STAT::mean(values);
			STAT::ConfidenceInterval(values, err,Q);
			break;
			}
		case kMedian:
			{
			mean=STAT::median(values);
			STAT::ConfidenceInterval(values, err,Q);
			break;
			}
		}

		// error are in absolute values, will scale to the mean point
		//
		if ( mean< err.first) cout<<"Error mean point is outside the error bands 1:"<<mean<<" < "<<err.first <<endl;
		if ( mean> err.second) cout<<"Error mean point is outside the error bands 2"<<mean<<" > "<<err.second<< endl;

		if( mean<err.first or  mean> err.second )
			{
			cout<<"[STAT]::[INFO]: Type= "<<type<<endl;
			cout<<"[STAT]::[INFO]: Mean="<< STAT::mean(values) <<endl;
			cout<<"[STAT]::[INFO]: Median="<< STAT::median(values) <<endl;
			cout<<"[STAT]::[INFO]: RMS="<< STAT::rms(values) <<endl;
			cout<<"[STAT]::[INFO]: Err="<< err.first<<" "<<err.second<<endl;
			cout<<"[STAT]::[INFO]: Values=";
			for(auto&v : values) cout<<v<<" ";
			cout<<endl;
			}

		float elow = mean-err.first;
		float ehigh = err.second -mean;

		int npoint= R->GetN();
		R->SetPoint( npoint, unf_->GetBinCenter( iBin ), mean);
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

void BootStrapBase::info(){
	cout <<"------- BOOTSTRAP BASE ------- "<<endl;
	cout <<"Seed = "<< seed_<<endl;
	cout <<"Ntoys = "<<Ntoys_<<endl;
	cout <<"SumW2 = "<<SumW2_<<endl;
	cout <<"Toy = ";
		switch(type_)
		{
		case kBootstrap: cout<<"kBootstrap"; break;
		case kToy: cout<<"kToy"; break;
		case kIterBias: cout<<"kIterBias"; break;
		case kMatrix: cout<<"kMatrix"; break;
		}
		cout<<endl;
		//<<type_<<" | kBootstrap 0 ; kToys 1 ; kIterBias 2"<<endl;
	if (type_ == kIterBias ) cout<<"N Iter Bias = "<<Nib_<<endl;
	cout <<"------------------------------ "<<endl;
}

TH1F* BootStrapBase::GetToyDistribution(int bin)
{
	if ( bootstrap_.size() == 0 )return NULL;
	if ( bin< 1 or bin > bootstrap_[0]->GetNbinsX() ) return NULL;
	vector<float> values;
	for(size_t i=0;i< bootstrap_.size() ;++i)
		values.push_back( bootstrap_[i]->GetBinContent(bin) );
	float low= STAT::min(values);
	float high= STAT::max(values);
	TH1F * h = new TH1F(Form("toy_distribution_bin%d",bin),Form("toy_bin%d",bin),100,low,high);
	//float sigma = STAT::rms(values);
	STAT::Fill(values,h);
	return h;
}

TH2F* BootStrapBase::GetToyDistribution(int bin1,int bin2)
{
	if ( bootstrap_.size() == 0 )return NULL;
	if ( bin1< 1 or bin1 > bootstrap_[0]->GetNbinsX() ) return NULL;
	if ( bin2< 1 or bin2 > bootstrap_[0]->GetNbinsX() ) return NULL;
	vector<float> valuesX;
	vector<float> valuesY;
	for(size_t i=0;i< bootstrap_.size() ;++i)
		{
		valuesX.push_back( bootstrap_[i]->GetBinContent(bin1) );
		valuesY.push_back( bootstrap_[i]->GetBinContent(bin2) );
		}
	float Xlow= STAT::min(valuesX);
	float Xhigh= STAT::max(valuesX);
	float Ylow= STAT::min(valuesY);
	float Yhigh= STAT::max(valuesY);
	TH2F * h = new TH2F(Form("toy_distribution_bin%d_bin%d",bin1,bin2),Form("toy_bin%d_%d",bin1,bin2),100,Xlow,Xhigh,100,Ylow,Yhigh);
	//float sigma = STAT::rms(values);
	STAT::Fill(valuesX,valuesY,h);
	return h;
}

