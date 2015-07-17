#include "interface/BootStrap.hpp"
//#include "RooUnfoldResponse.h"

BootStrap::BootStrap(): BootStrapBase() {
	u_reco_ = NULL;
	u_truth_= NULL;
	u_resp_ = NULL;

	f_reco_ = NULL;
	f_truth_= NULL;
	f_resp_ = NULL;
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

}
void BootStrap::SetFMatrix(TH1D* reco, TH1D* truth, TH2D* resp)
{
	destroyPointer( f_reco_  );
	destroyPointer( f_truth_ );
	destroyPointer( f_resp_  );

	setPointer(reco ,f_reco_ );
	setPointer(truth,f_truth_ );
	setPointer(resp ,f_resp_ );

}

TH1D* BootStrap::Unfold(TH1D* h)
{
	RooUnfoldResponse RU_Resp(u_reco_,u_truth_,u_resp_,"resp","resp"); // for unfolding
//	RooUnfoldBayes( &RU_Resp, h, regParam_)	
	return (TH1D*)NULL;//TODO
}

TH1D* BootStrap::Fold(TH1D* h)
{

	return (TH1D*)NULL;//TODO
}
