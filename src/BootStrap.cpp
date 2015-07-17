#include "interface/BootStrap.hpp"
#include "TMath.h"
//#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldInvert.h"

BootStrap::BootStrap(): BootStrapBase() {

	u_reco_ = NULL;
	u_truth_= NULL;
	u_resp_ = NULL;

	f_reco_ = NULL;
	f_truth_= NULL;
	f_resp_ = NULL;

	f_resp_bkg_ = NULL;
	f_resp_eff_ = NULL;
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
		f_resp_eff_ = NULL;
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
	f_resp_eff_ = NULL;

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
		return (TH1D*)u.Hreco( RooUnfold::kNoError);
		break;
		}
	case kSvd:
		{
		RooUnfoldSvd u( &RU_Resp, h, regParam_)	;
		return (TH1D*)u.Hreco( RooUnfold::kNoError);
		break;
		}
	case kInv:
		{
		RooUnfoldInvert u( &RU_Resp, h)	;
		return (TH1D*)u.Hreco( RooUnfold::kNoError);
		break;
		}
	}
	return (TH1D*)NULL;// problem !!
}

TH1D* BootStrap::Fold(TH1D* h)
{
	// folding should do the opposite:
	TH1D *reco  = f_reco_;
	TH1D *truth = f_truth_;
	TH2D *resp  = f_resp_;

	if( truth == NULL and reco == NULL and resp == NULL)
		{
		reco  = u_reco_;
		truth = u_truth_;
		resp  = u_resp_;
		}

	//get projections
	//"measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
	//construct bkg
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

	// construct efficiency
	if (f_resp_eff_ == NULL) 
	{
		f_resp_eff_ = resp->ProjectionY();
		for(int i=0;i<=f_resp_eff_->GetNbinsX()+1 ;++i)
		{
			float gen = truth->GetBinContent(i);
			float genreco = f_resp_eff_->GetBinContent(i);
			float eff =0;
			if(gen != 0) eff=genreco/gen;
			f_resp_eff_ -> SetBinContent(i,eff);
		}
	}

	if( f_resp_smear_ == NULL) 
	{
	
	}

	TH1D *fold = (TH1D*)h->Clone("fold");
	// 1. apply eff.
	for( int i=0;i<= fold->GetNbinsX()+1 ;++i)
	{
		float c= fold->GetBinContent(i);	
		float e= fold->GetBinError(i);
		float eff = f_resp_eff_ -> GetBinContent(i);
		fold -> SetBinContent(i, c*eff);
		fold -> SetBinError(i, e*eff);
	}
	// 2. apply smearings.
	// 3. add bkg
	for( int i=0;i<= fold->GetNbinsX()+1 ;++i)
	{
		float c= fold->GetBinContent(i);	
		float e= fold->GetBinError(i);
		float bkg = f_resp_bkg_ -> GetBinContent(i);
		fold -> SetBinContent(i, c + bkg);
		fold -> SetBinError(i, e + TMath::Sqrt(bkg));
	}
	return (TH1D*)NULL;//TODO
}
