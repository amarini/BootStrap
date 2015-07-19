#include "interface/BootStrap.hpp"
#include "TMath.h"
//#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldInvert.h"

#include <iostream>

#define VERBOSE 0

BootStrap::BootStrap(): BootStrapBase() {

	u_reco_ = NULL;
	u_truth_= NULL;
	u_resp_ = NULL;

	f_reco_ = NULL;
	f_truth_= NULL;
	f_resp_ = NULL;

	f_resp_bkg_ = NULL;
	f_resp_smear_ = NULL;
}

BootStrap::~BootStrap(){

	destroyPointer( u_reco_  );
	destroyPointer( u_truth_ );
	destroyPointer( u_resp_  );
                                
	destroyPointer( f_reco_  );
	destroyPointer( f_truth_ );
	destroyPointer( f_resp_  );

}

void BootStrap::SetUMatrix(TH1D* reco, TH1D* truth, TH2D* resp)
{
	destroyPointer( u_reco_  );
	destroyPointer( u_truth_ );
	destroyPointer( u_resp_  );

	setPointer(reco,u_reco_ );
	setPointer(truth,u_truth_ );
	setPointer(resp,u_resp_ );

	if ( f_reco_ == NULL and 
			f_truth_ == NULL and 
			f_resp_ == NULL )
	{
		f_resp_bkg_ = NULL;
		f_resp_smear_ = NULL;
	}


}
void BootStrap::SetFMatrix(TH1D* reco, TH1D* truth, TH2D* resp)
{
	destroyPointer( f_reco_  );
	destroyPointer( f_truth_ );
	destroyPointer( f_resp_  );

	setPointer(reco ,f_reco_ );
	setPointer(truth,f_truth_ );
	setPointer(resp ,f_resp_ );

	f_resp_bkg_ = NULL;
	f_resp_smear_ = NULL;
}

TH1D* BootStrap::Unfold(TH1D* h)
{
	RooUnfoldResponse RU_Resp(u_reco_,u_truth_,u_resp_,"resp","resp"); // for unfolding
//	
	switch (unfType_)
	{
	case kBayes:
		{
		RooUnfoldBayes u( &RU_Resp, h, regParam_);
		u.SetVerbose(-1);
		return (TH1D*)u.Hreco( RooUnfold::kNoError);
		break;
		}
	case kSvd:
		{
		RooUnfoldSvd u( &RU_Resp, h, regParam_)	;
		u.SetVerbose(-1);
		return (TH1D*)u.Hreco( RooUnfold::kNoError);
		break;
		}
	case kInv:
		{
		RooUnfoldInvert u( &RU_Resp, h)	;
		u.SetVerbose(-1);
		return (TH1D*)u.Hreco( RooUnfold::kNoError);
		break;
		}
	}
	return (TH1D*)NULL;// problem !!
}

void BootStrap::ConstructProjections(TH1D*reco,TH1D*truth,TH2D*resp){
	if (VERBOSE >0 ) cout<<"[BootStrap]::[Fold]::[DEBUG] Folding distribution "<<endl;
	//get projections
	//"measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
	//construct bkg
	if (VERBOSE >0 ) cout<<"[BootStrap]::[Fold]::[DEBUG] Construct bkg "<<endl;
	if (f_resp_bkg_ == NULL) 
	{
		f_resp_bkg_ = resp->ProjectionX();
		for(int i=0;i<=f_resp_bkg_->GetNbinsX()+1 ;++i)
		{
			float r = reco->GetBinContent(i);
			float genreco = f_resp_bkg_->GetBinContent(i);
			f_resp_bkg_ -> SetBinContent(i, r - genreco);
		}
	}

	if (VERBOSE >0 ) cout<<"[BootStrap]::[Fold]::[DEBUG] Construct smear "<<endl;
	if( f_resp_smear_ == NULL) 
	{
		f_resp_smear_ = (TH2D*)resp->Clone("reps_smear");
		f_resp_smear_ -> Reset("ACE");
		for(int i=0;i<= resp->GetNbinsX()+1 ;++i)
		{
			for(int j=0;j<= resp->GetNbinsY()+1 ;++j)
			{
			float c = resp->GetBinContent(i,j);
			float g = truth->GetBinContent(j);
			if(g==0) continue;
			f_resp_smear_->SetBinContent(i,j, c/g) ;
			}
		}
	}

	return;
}

TH1D* BootStrap::Fold(TH1D* h)
{
	if (VERBOSE >0 ) cout<<"[BootStrap]::[Fold]::[DEBUG] Folding distribution "<<endl;
	// folding should do the opposite:
	TH1D *reco  = f_reco_;
	TH1D *truth = f_truth_;
	TH2D *resp  = f_resp_;

	if( truth == NULL and reco == NULL and resp == NULL)
		{
		if (VERBOSE >0 ) cout<<"[BootStrap]::[Fold]::[DEBUG] Using Unfold matrixes "<<endl;
		reco  = u_reco_;
		truth = u_truth_;
		resp  = u_resp_;
		}

	// construct the matrix f_resp_smear and f_resp_bkg_ if they are NULL
	ConstructProjections(reco,truth,resp);

	// ------------------------------
	if (VERBOSE >0 ) cout<<"[BootStrap]::[Fold]::[DEBUG] Folding distribution: 1. Appl smearings "<<endl;
	TH1D *fold = (TH1D*)h->Clone("fold");
	fold->Reset("ACE");
	// 1. apply eff and smearings.
	for(int i=0;i<= fold->GetNbinsX() +1  ;++i)
	{
		float S=0;
		for(int j=0;j<= f_resp_smear_->GetNbinsY() +1  ;++j)
		{
			float c = h -> GetBinContent(j);	
			float s = f_resp_smear_ -> GetBinContent(i,j);
			S += s*c;
		}
		fold -> SetBinContent(i,S);
	}

	if (VERBOSE >0 ) cout<<"[BootStrap]::[Fold]::[DEBUG] Folding distribution: 2. Appl bkg "<<endl;
	// 2. add bkg
	for( int i=0;i<= fold->GetNbinsX()+1 ;++i)
	{
		float c= fold->GetBinContent(i);	
		float e= fold->GetBinError(i);
		float bkg = f_resp_bkg_ -> GetBinContent(i);
		fold -> SetBinContent(i, c + bkg);
		fold -> SetBinError(i, e + TMath::Sqrt(bkg));
	}
	if (VERBOSE >0 ) cout<<"[BootStrap]::[Fold]::[DEBUG] Folding distribution: Return; "<<endl;
	return (TH1D*)fold;
}

void BootStrap::info(){
	BootStrapBase::info();
	cout <<"RegParam = "<<regParam_	<<endl;
	cout <<"unfType = "<<unfType_<<"| kBayes="<<kBayes<<", kSvd="<<kSvd<<", kInv="<<kInv<<endl;
	cout <<"------------------------------ "<<endl;
}
