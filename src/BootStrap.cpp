#include "interface/BootStrap.hpp"
#include "TMath.h"
//#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldInvert.h"

#include <iostream>

#define VERBOSE 0

BootStrap::BootStrap(): BootStrapBase() {

}

BootStrap::~BootStrap(){

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


void BootStrap::info(){
	BootStrapMatrix::info();
	cout <<"---------- BOOTSTRAP ----------- "<<endl;
	cout <<"RegParam = "<<regParam_	<<endl;
	cout <<"unfType = "<<unfType_<<"| kBayes="<<kBayes<<", kSvd="<<kSvd<<", kInv="<<kInv<<endl;
	cout <<"------------------------------ "<<endl;
}
