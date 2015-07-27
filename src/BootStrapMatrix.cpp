#include "interface/BootStrapMatrix.hpp"
#include <iostream>
#include <ctime>
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

	u_resp_bkg_ = NULL;
	u_resp_smear_ = NULL;

	matrixConstructed_=false;
}

BootStrapMatrix::BootStrapMatrix( BootStrapMatrix &x) :
	BootStrapBase( x )
{
	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[BootStrapMatrix]::[copy constructor]"<<endl;
	u_reco_ = (x.u_reco_)? (TH1D*) x.u_reco_->Clone( Form("%s_%u", x.u_reco_->GetName(),(unsigned)time(NULL)) ) : NULL;
	u_truth_= (x.u_truth_)? (TH1D*) x.u_truth_->Clone( Form("%s_%u",x.u_truth_->GetName(),(unsigned)time(NULL)) ) : NULL;
	u_resp_ = (x.u_resp_)? (TH2D*) x.u_resp_->Clone( Form("%s_%u", x.u_resp_->GetName(), (unsigned)time(NULL)) ) : NULL;

	f_reco_ = (x.f_reco_)? (TH1D*) x.f_reco_->Clone( Form("%s_%u", x.f_reco_->GetName(),(unsigned)time(NULL)) ) : NULL;
	f_truth_= (x.f_truth_)? (TH1D*) x.f_truth_->Clone( Form("%s_%u",x.f_truth_->GetName(),(unsigned)time(NULL)) ) : NULL;
	f_resp_ = (x.f_resp_)? (TH2D*) x.f_resp_->Clone( Form("%s_%u", x.f_resp_->GetName(), (unsigned)time(NULL)) ) : NULL;

	f_resp_bkg_ = NULL;
	f_resp_smear_ = NULL;

	u_resp_bkg_ = NULL;
	u_resp_smear_ = NULL;

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
	matrixConstructed_=false;

	destroyPointer( u_reco_  );
	destroyPointer( u_truth_ );
	destroyPointer( u_resp_  );

	setPointer(reco,u_reco_ );
	setPointer(truth,u_truth_ );
	setPointer(resp,u_resp_ );

	destroyPointer( u_resp_bkg_ );
	destroyPointer( u_resp_smear_ );
	
	// destroy them in any-case to preserve consistency
	destroyPointer (f_resp_bkg_ ); 
	destroyPointer (f_resp_smear_ );


}
void BootStrapMatrix::SetFMatrix(TH1D* reco, TH1D* truth, TH2D* resp)
{
	destroyPointer( f_reco_  );
	destroyPointer( f_truth_ );
	destroyPointer( f_resp_  );

	setPointer(reco ,f_reco_ );
	setPointer(truth,f_truth_ );
	setPointer(resp ,f_resp_ );

	destroyPointer( f_resp_bkg_ );
	destroyPointer( f_resp_smear_ );

}

void BootStrapMatrix::ConstructProjections(TH1D*reco,TH1D*truth,TH2D*resp){
	// default construct projections for folding
	ConstructProjections(reco,truth,resp,f_resp_bkg_,f_resp_smear_);
}

void BootStrapMatrix::ConstructProjections(TH1D*reco,TH1D*truth,TH2D*resp, TH1D* &bkg, TH2D* &smear){
	if (VERBOSE >0 ) cout<<"[BootStrapMatrix]::[ConstructProjections]::[DEBUG] Folding distribution "<<endl;
	//get projections
	//"measured" and "truth" give the projections of "response" onto the X-axis and Y-axis respectively,
	//construct bkg
	if (VERBOSE >0 ) cout<<"[BootStrapMatrix]::[ConstructProjections]::[DEBUG] Construct bkg "<<endl;
	if (bkg == NULL) 
	{
		bkg = resp->ProjectionX();
		for(int i=0;i<=bkg->GetNbinsX()+1 ;++i)
		{
			float r = reco->GetBinContent(i);
			float genreco = bkg->GetBinContent(i);
			bkg -> SetBinContent(i, r - genreco);
		}
	}

	if (VERBOSE >0 ) cout<<"[BootStrapMatrix]::[ConstructProjections]::[DEBUG] Construct smear "<<endl;
	if( smear == NULL) 
	{
		smear = (TH2D*)resp->Clone(Form("%s_reps_smear",resp->GetName()));
		smear -> Reset("ACE");
		for(int i=0;i<= resp->GetNbinsX()+1 ;++i)
		{
			for(int j=0;j<= resp->GetNbinsY()+1 ;++j)
			{
			float c = resp->GetBinContent(i,j);
			float g = truth->GetBinContent(j);
			if(g==0) continue;
			smear->SetBinContent(i,j, c/g) ;
			}
		}
	}


	if (VERBOSE >0 ) cout<<"[BootStrapMatrix]::[ConstructProjections]::[DEBUG] DONE "<<endl;
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
	cout <<"------- BOOTSTRAP MATRIX ------- "<<endl;
	cout <<"-------------------------------- "<<endl;
}

void BootStrapMatrix::ConstructMatrixes(TH1D* data){
	// Construct the TMatrix used for unfolding with ML
	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[ConstructMatrixes]::[2] Constructing Matrixes"<<endl;

	//matrix depend on the u_matrixes and on the data
	//if (matrixConstructed_) return;
	//matrixConstructed_ = true;

	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[ConstructMatrixes]::[2] Matrixes not constructed"<<endl;

	ConstructProjections(u_reco_,u_truth_,u_resp_, u_resp_bkg_, u_resp_smear_);

	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[ConstructMatrixes]::[2] Resize"<<endl;
	int nreco= u_reco_->GetNbinsX();
	int ntruth= u_truth_->GetNbinsX();
	y.ResizeTo( nreco ) ;
	l.ResizeTo( ntruth ) ;
	K.ResizeTo( ntruth, nreco );

	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[ConstructMatrixes]::[2] Converting"<<endl;
	//K=getMatrix(  u_resp_ );
	//K=getMatrix(  u_resp_smear_ );
	K=getMatrix(  u_resp_ );// STT

	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[ConstructMatrixes]::[2] Converting S"<<endl;
	// scale to truth STT
	S.ResizeTo( nreco,nreco); // no Overflow
	for(int i=1 ;i<= nreco ;++i)
		{
		if(VERBOSE>1)cout<<"[BootStrapMatrix]::[ConstructMatrixes]::[2] bin "<< i<<endl;
		//S(i-1,i-1) =  data->GetBinContent(i);// sqrt  SUMW2
		if (u_truth_->GetBinContent(i) == 0 ) S(i-1,i-1) = 1;
		else S(i-1,i-1) = data->GetBinContent(i) / (u_truth_->GetBinContent(i) * u_truth_->GetBinContent(i)) ;//STT -- scale error, this is the covariance matrix
		}

	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[ConstructMatrixes]::[2] Converting y"<<endl;
	y= getVector(data); // y should be vector subtracted
	for(int i=0;i<y.GetNrows() ;++i) // no overflow
	{
		y(i) -= u_resp_bkg_->GetBinContent(i+1);
	}

	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[ConstructMatrixes]::[2] Additional matrixes. K."<<endl;
	Kt.ResizeTo( K.GetNcols(),K.GetNrows() );
	Kt.Transpose(K);
	//TMatrixD Kt(K.GetNcols(),K.GetNrows()); Kt.Transpose(K);
	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[ConstructMatrixes]::[2] Additional matrixes. A."<<endl;
	A.ResizeTo( K.GetNcols(), K.GetNcols() ) ;
	A=Kt*S*K; 

	// make sure A is >0
		float epsilon=0;
		while ( A.Determinant() == 0)
		{
			for(int iDiag=0;iDiag<A.GetNrows();iDiag++)
				A(iDiag,iDiag)+=0.0001;
			epsilon += 0.0001;
		}
		if(epsilon >0 ) cout<<"[BootStrapMatrix]::[Unfold]::[INFO] used epsilon="<<epsilon<<" for inversion purpose"<<endl;
	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[Unfold]::[2] Inverting A"<<endl;
	A.Invert();

	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[ConstructMatrixes]::[2] Additional matrixes. B."<<endl;
	B.ResizeTo(A.GetNrows(),S.GetNcols());
	B=A*Kt*S;
	Bt.ResizeTo( B.GetNcols(),B.GetNrows() );	
	Bt.Transpose( B );
	cov.ResizeTo(B.GetNrows(),Bt.GetNcols());
	cov = (B*S*Bt); // wrt scaled

	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[ConstructMatrixes]::[2] DONE"<<endl;
}

TMatrixD BootStrapMatrix::getMatrix(TH2*h, bool useOverFlow)
{
	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[getMatrix]::[2] START"<<endl;
	if (h==NULL) cout<<"[BootStrapMatrix]::[getMatrix]::[ERROR] h is NULL"<<endl;
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
	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[getMatrix]::[2] END"<<endl;
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

void BootStrapMatrix::printMatrix(TMatrixD&m, string name)
{
	cout << "---------- M "<<name<<" ---------"<<endl;
	for(int i=0;i<m.GetNrows() ;++i)
	{
		for(int j=0;j<m.GetNcols();++j)
			cout<< m(i,j) <<" ";
		cout <<endl;
	}
	cout << "-----------E "<<name<<" ---------"<<endl;
}

TH1D* BootStrapMatrix::Unfold(TH1D*data)
{
	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[Unfold]::[2] Constructing Matrixes"<<endl;
	ConstructMatrixes(data); // get y = (data sub) , K (smear matrix) , 
	// implement the pseudo inverse

	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[Unfold]::[2] Work with Matrixes"<<endl;


	if (VERBOSE >1)
		{
		printMatrix( K,"K");
		printMatrix( S,"S");
		printMatrix( A,"A");
		}
	//l=A*Kt*S*y;
	//TMatrixD B=A*Kt*S;
	l=B*y;

	// STT
	for(int i=0 ;i<l.GetNrows() ;++i)
	{
		l(i) *= u_truth_->GetBinContent(i+1);
	}

	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[Unfold]::[2] Compute Analytical covariance matrix"<<endl;
	// Analytical error propagation	
	//TMatrixD Bt(B.GetNcols(),B.GetNrows()); Bt.Transpose(B);
	//TMatrixD cov; cov.ResizeTo(l.GetNrows(),l.GetNrows());
	//cov = (B*S*Bt);  //covariance matrix after linear transformation

	// COPYING INFO
	TH1D *result = (TH1D*)u_truth_->Clone(Form("%s_analyticalunf",data->GetName()));
	result->Reset("ACE");
	for(int i=1;i<=result->GetNbinsX() ;++i)
		{
		result->SetBinContent(i, l(i-1) );
		double tr = u_truth_->GetBinContent(i); //STT
		result->SetBinError(i, TMath::Sqrt(TMath::Max( double(cov(i-1,i-1) * tr * tr),double(0.))) ) ; // * truth[i-1,j-1]
		}
	if(VERBOSE>1)cout<<"[BootStrapMatrix]::[Unfold]::[2] Done"<<endl;

	return result;
}


