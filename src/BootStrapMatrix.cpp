#include "interface/BootStrapMatrix.hpp"
#include <iostream>
#include "TMath.h"

using namespace std;

#define VERBOSE 0

BootStrapMatrix::BootStrapMatrix(): BootStrapBase() {
	u_reco_ = NULL;
	u_truth_= NULL;
	u_resp_ = NULL;

	f_reco_ = NULL;
	f_truth_= NULL;
	f_resp_ = NULL;

	f_resp_bkg_ = NULL;
	f_resp_smear_ = NULL;
}

BootStrapMatrix::~BootStrapMatrix(){

	destroyPointer( u_reco_  );
	destroyPointer( u_truth_ );
	destroyPointer( u_resp_  );
                                
	destroyPointer( f_reco_  );
	destroyPointer( f_truth_ );
	destroyPointer( f_resp_  );

}

void BootStrapMatrix::SetUMatrix(TH1D* reco, TH1D* truth, TH2D* resp)
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
void BootStrapMatrix::SetFMatrix(TH1D* reco, TH1D* truth, TH2D* resp)
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

void BootStrapMatrix::ConstructProjections(TH1D*reco,TH1D*truth,TH2D*resp){
	if (VERBOSE >0 ) cout<<"[BootStrapMatrix]::[Fold]::[DEBUG] Folding distribution "<<endl;
	//get projections
	//"measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
	//construct bkg
	if (VERBOSE >0 ) cout<<"[BootStrapMatrix]::[Fold]::[DEBUG] Construct bkg "<<endl;
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

	if (VERBOSE >0 ) cout<<"[BootStrapMatrix]::[Fold]::[DEBUG] Construct smear "<<endl;
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

TH1D* BootStrapMatrix::Fold(TH1D* h)
{
	if (VERBOSE >0 ) cout<<"[BootStrapMatrix]::[Fold]::[DEBUG] Folding distribution "<<endl;
	// folding should do the opposite:
	TH1D *reco  = f_reco_;
	TH1D *truth = f_truth_;
	TH2D *resp  = f_resp_;

	if( truth == NULL and reco == NULL and resp == NULL)
		{
		if (VERBOSE >0 ) cout<<"[BootStrapMatrix]::[Fold]::[DEBUG] Using Unfold matrixes "<<endl;
		reco  = u_reco_;
		truth = u_truth_;
		resp  = u_resp_;
		}

	// construct the matrix f_resp_smear and f_resp_bkg_ if they are NULL
	ConstructProjections(reco,truth,resp);

	// ------------------------------
	if (VERBOSE >0 ) cout<<"[BootStrapMatrix]::[Fold]::[DEBUG] Folding distribution: 1. Appl smearings "<<endl;
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

	if (VERBOSE >0 ) cout<<"[BootStrapMatrix]::[Fold]::[DEBUG] Folding distribution: 2. Appl bkg "<<endl;
	// 2. add bkg
	for( int i=0;i<= fold->GetNbinsX()+1 ;++i)
	{
		float c= fold->GetBinContent(i);	
		float e= fold->GetBinError(i);
		float bkg = f_resp_bkg_ -> GetBinContent(i);
		fold -> SetBinContent(i, c + bkg);
		fold -> SetBinError(i, e + TMath::Sqrt(bkg));
	}
	if (VERBOSE >0 ) cout<<"[BootStrapMatrix]::[Fold]::[DEBUG] Folding distribution: Return; "<<endl;
	return (TH1D*)fold;
}

void BootStrapMatrix::info(){
	BootStrapBase::info();
}

void BootStrapMatrix::ConstructMatrixes(TH1D* data){
	// Construct the TMatrix used for unfolding with ML
	ConstructProjections(u_reco_,u_truth_,u_resp_);

	int nreco= u_reco_->GetNbinsX();
	int ntruth= u_truth_->GetNbinsX();
	y.ResizeTo( nreco ) ;
	l.ResizeTo( ntruth ) ;
	K.ResizeTo( ntruth, nreco );

	K=getMatrix(  u_resp_ );
	S.ResizeTo( nreco,nreco); // no Overflow
	for(int i=1 ;i<= nreco ;++i)
		{
		S(i,i) =  data->GetBinError(i);
		}

	y= getVector(data);
}

TMatrixD BootStrapMatrix::getMatrix(TH2*h, bool useOverFlow)
{
	TMatrixD r;
	int n=h->GetNbinsX();
	int m=h->GetNbinsY();
	if (useOverFlow)n+=2;
	if (useOverFlow)m+=2;
	r.ResizeTo(n,m);
	r.Zero(); //make sure is set to 0
	if(!useOverFlow) 
	{
		for(int iBin=1;iBin<=h->GetNbinsX();iBin++)
		for(int jBin=1;jBin<=h->GetNbinsY();jBin++)
			r(iBin-1,jBin-1)=h->GetBinContent(iBin,jBin);
	
	}
	else{
		for(int iBin=0;iBin<=h->GetNbinsX()+1;iBin++)
		for(int jBin=0;jBin<=h->GetNbinsY()+1;jBin++)
			r(iBin,jBin)=h->GetBinContent(iBin,jBin);
	}
	return r;
}

TVectorD BootStrapMatrix::getVector(TH1*h,bool useOverFlow)
{
	TVectorD r;
	int n=h->GetNbinsX();
	if (useOverFlow)n+=2;
	r.ResizeTo(n);
	int i=0;
	if(useOverFlow) { r(i) = h->GetBinContent(i);i++;}
	for(int iBin=1;iBin<=h->GetNbinsX();iBin++)
		{
			r(i)=h->GetBinContent(iBin);i++;
		}
	if(useOverFlow) { 
			r(i) = h->GetBinContent(h->GetNbinsX()+1);i++;

	}
	return r;
}

TH1D* BootStrapMatrix::Unfold(TH1D*data)
{
	cout<<"[BootStrapMatrix]::[Unfold]::[WARNING] procedure not checked !"<<endl;

	ConstructMatrixes(data);
	// implement the pseudo inverse

	TMatrixD Kt(K.GetNcols(),K.GetNrows()); Kt.Transpose(K);
	TMatrixD A=Kt*S*K; 
	// make sure A is >0
		float epsilon=0;
		while ( A.Determinant() == 0)
		{
			for(int iDiag=0;iDiag<A.GetNrows();iDiag++)
				A(iDiag,iDiag)+=0.0001;
			epsilon += 0.0001;
		}
		if(epsilon >0 ) cout<<"[BootStrapMatrix]::[Unfold]::[INFO] used epsilon="<<epsilon<<" for inversion purpose"<<endl;
	A.Invert();
	//l=A*Kt*S*y;
	TMatrixD B=A*Kt*S;
	l=B*y;

	// Analytical error propagation	
	TMatrixD Bt(B.GetNcols(),B.GetNrows()); Bt.Transpose(B);
	TMatrixD cov; cov.ResizeTo(l.GetNrows(),l.GetNrows());
	cov = (B*S*Bt);  //covariance matrix after linear transformation
	TH1D *result = (TH1D*)data->Clone(Form("%s_analyticalunf",data->GetName()));
	result->Reset("ACE");
	for(int i=1;i<=result->GetNbinsX() ;++i)
		{
		result->SetBinContent(i, l(i-1) );
		result->SetBinError(i, TMath::Sqrt(TMath::Max( double(cov(i-1,i-1)),double(0.))) ) ;
		}

	return result;
}


