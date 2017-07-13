#include "interface/BootStrap.hpp"
#include "TMath.h"
//#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldInvert.h"

#include "TUnfoldDensity.h"

#include <iostream>

#define VERBOSE 0

BootStrap::BootStrap(): BootStrapMatrix() {
	unfType_ = kInv;
	regParam_ = 3;
}

BootStrap::BootStrap( BootStrap &x) : BootStrapMatrix( x )
{
	if(VERBOSE>1)cout<<"[BootStrap]::[BootStrap]::[copy constructor]"<<endl;
	unfType_ = x.unfType_;
	regParam_ = x.regParam_;
}

BootStrap::~BootStrap(){
}

TH1D* BootStrap::UnfoldLikelihood(TH1D*h){
	// keep track of current settings
	UnfoldType mytype=unfType_;
	// set inversion
	unfType_=kInv;
	//get result
	TH1D * r= Unfold(h);
	//set back current settings
	unfType_=mytype;
	// return result
	return r;
}

TH1D* BootStrap::Unfold(TH1D* h)
{
	RooUnfoldResponse RU_Resp(u_reco_,u_truth_,u_resp_,"resp","resp"); // for unfolding
	RU_Resp.UseOverflow(); // generally better to use overflow
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
    case kTUnfoldDensity:
        {
            // make sure there are projections
            ConstructProjections(u_reco_,u_truth_,u_resp_,u_resp_bkg_,u_resp_smear_);
            //
            TUnfold::ERegMode regMode =TUnfold::kRegModeCurvature;
            //TUnfold::ERegMode regMode =TUnfold::kRegModeNone;
            // preserve the area (constrain)
            //TUnfold::EConstraint constraintMode=TUnfold::kEConstraintArea;
            TUnfold::EConstraint constraintMode=TUnfold::kEConstraintNone;
            // ?
            TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeBinWidth;

            TH2D* m= (TH2D*)u_resp_smear_->Clone("tmp");
            for(int gen=1;gen<=m->GetNbinsY();++gen)
            {
                double Sum=0;
                for(int reco=1;reco<=m->GetNbinsX();++reco)
                {
                    Sum+= m->GetBinContent(reco,gen);
                }
                m->SetBinContent(0,gen,1-Sum);
            }

            TUnfoldDensity unfold(m,TUnfold::kHistMapOutputVert, regMode, constraintMode ,densityFlags);

            TH1D *h2=(TH1D*)h->Clone("htmp");

            // check errors
            for(int i =0 ; i<= h2->GetNbinsX()+1 ;++i) 
                if (h2->GetBinError(i)==0) h2->SetBinError(i,TMath::Sqrt(h2->GetBinContent(i)));
            for(int i =0 ; i<= h2->GetNbinsX()+1 ;++i) 
                if (h2->GetBinContent(i)==0) h2->SetBinError(i,1);


            cout <<"UNFOLDING H=";
            for(int i =0 ; i<= h2->GetNbinsX()+1 ;++i) 
                cout << i << ":("<<h2->GetBinContent(i)<<","<<h2->GetBinError(i)<<") ";
            cout <<endl;

            unfold.SetInput(h2);
            // subtract background, normalized to data luminosity
            //   //  with 10% scale error each
            Double_t scale_bgr=1.0;
            Double_t dscale_bgr=0.;
            unfold.SubtractBackground(u_resp_bkg_,"background1",scale_bgr,dscale_bgr);

            int iBest=regParam_;
            if (regParam_ <0)
            {
                Int_t nScan=30;
                TSpline *logTauX,*logTauY;
                TGraph *lCurve;
                // this method scans the parameter tau and finds the kink in the L curve
                //   // finally, the unfolding is done for the best choice of tau
                iBest=unfold.ScanLcurve(nScan,0.,0.,&lCurve,&logTauX,&logTauY);
            }
            //unfold.DoUnfold(0.);

            TH1D*r= (TH1D*)unfold.GetOutput("PT(unfold,stat+bgrerr)");

            delete h2;
            delete m;

            return r;

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
	cout <<"Overflow = Yes"<<endl;
	cout <<"------------------------------ "<<endl;
}
