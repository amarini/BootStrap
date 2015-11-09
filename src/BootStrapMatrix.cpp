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
	
	negCorr= x.negCorr;
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
			// construct back the original distribution
			// error are the difference
			float re = reco->GetBinError(i);
			float gre = bkg->GetBinError(i);
			bkg -> SetBinError(i, TMath::Sqrt(re*re - gre * gre ));
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
	cout <<" Negative Correction is: ";
	switch(negCorr)
	{
		case kNegNone: cout<<"None"<<endl; break;
		case kNegZero: cout<<"Zero"<<endl; break;
		case kNegZeroProp: cout<<"Zero & Propagate"<<endl; break;
		case kNegMoveProp: cout<<"Move & Propagate"<<endl; break;
		case kNegReplProp: cout<<"Replace & Propagate"<<endl; break;
		default: cout <<" ??? "<<endl;
	}
	cout <<endl;
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

//#define VERBOSE 2
TH1D* BootStrapMatrix::matrixSmear(){
	if (VERBOSE>1) cout <<"[BootStrapMatrix]::[matrixSmear]::[2] Start"<<endl;

	static int matrixWarning = 0;
	if (not SumW2_ and matrixWarning == 0 )
		{
		++matrixWarning ;
		cout<<"[BootStrapMatrix]::[matrixSmear]::[WARNING] SumW2 is not active.... check it ! "<<endl;
		}

	// Preliminary check
	// unfold data -- will used in some results, like median, rms...
	if(unf_==NULL) unf_ = Unfold(data_);
	// no need of folded distribution


	if (VERBOSE>1) cout <<"[BootStrapMatrix]::[matrixSmear]::[2] Construct projections"<<endl;
	// Produce a Toy with the smearing from the matrix
	ConstructProjections(u_reco_,u_truth_,u_resp_, u_resp_bkg_, u_resp_smear_ ) ;

	if (VERBOSE>1) cout <<"[BootStrapMatrix]::[matrixSmear]::[2] Construct extra projections"<<u_resp_ << "|"<<u_resp_bkg_<<endl;
	// clone matrix at this stage
	TH1D * _reco_  = NULL ;
	TH1D * _truth_ = NULL ;

	TH2D * _resp_  = (TH2D*)u_resp_ -> Clone(Form("%s_smearmatrix",u_resp_->GetName() ));
	TH1D * _bkg_   = (TH1D*)u_resp_bkg_ ->Clone(Form("%s_smearmatrix",u_resp_bkg_->GetName()) ); 

	if (VERBOSE>1) cout <<"[BootStrapMatrix]::[matrixSmear]::[2] Construct extra projections 2"<<endl;
	TH1D * _gen_notreco_ = u_resp_ -> ProjectionY();

		for(int i=0;i<=_gen_notreco_->GetNbinsX()+1 ;++i)
		{
			if (VERBOSE>1) cout <<"[BootStrapMatrix]::[matrixSmear]::[2] GenNotReco bin"<<i<<endl;
			float g = u_truth_->GetBinContent(i);
			float genreco = _gen_notreco_->GetBinContent(i);
			_gen_notreco_ -> SetBinContent(i, g - genreco);
			// construct back the original distribution
			// error are the difference
			float ge = u_truth_->GetBinError(i);
			float gre = _gen_notreco_->GetBinError(i);
			_gen_notreco_ -> SetBinError(i, TMath::Sqrt(ge*ge -gre * gre));
		}

	
	if (VERBOSE>1) cout <<"[BootStrapMatrix]::[matrixSmear]::[2] Compute smallest error"<<endl;
	// smear the matrixes, start with background. 
	// find the smallest error, 
	bool isBkgEmpty = true;
	bool isEffEmpty = true;
	float smallestError= 1;
	for( int i=0;i<= _bkg_->GetNbinsX() +1 ; ++i)
		if (_bkg_->GetBinContent(i) > 1.e-10 ) {
			isBkgEmpty = false;
			float e = _bkg_->GetBinError(i);
			if (e<smallestError ) smallestError = e;
		}
	for( int i=0;i<= _gen_notreco_->GetNbinsX() +1 ; ++i)
		if (_gen_notreco_->GetBinContent(i) > 1.e-10 ) {
			isEffEmpty= false;
			float e = _gen_notreco_->GetBinError(i);
			if (e<smallestError ) smallestError = e;
		}
	for(int i=0;i<= _resp_->GetNbinsX() +1 ; i++)
	for(int j=0;j<= _resp_->GetNbinsY() +1 ; j++)
	{
		if (_resp_->GetBinContent(i,j) > 1.e-10 ) {
			float e = _resp_->GetBinError(i,j);
			if (e<smallestError ) smallestError = e;

		}
			
	}

	if (VERBOSE>1) cout <<"[BootStrapMatrix]::[matrixSmear]::[2] Set smallest error"<<smallestError<<endl;
	// and set 0 bins to that error, may lead to an over estimation, as sumw2 set to 1 error with content 0.
	for( int i=0;i<= _bkg_->GetNbinsX() +1 and isBkgEmpty; ++i)
		{
		if (_bkg_->GetBinContent(i) <= 1.e-10 ) 
			_bkg_->SetBinError(i,smallestError);
		}

	for( int i=0;i<= _gen_notreco_->GetNbinsX() +1 ; ++i)
		if (_gen_notreco_->GetBinContent(i) < 1.e-10 ) {
			_gen_notreco_->SetBinError(i,smallestError);
		}

	for(int i=0;i<= _resp_->GetNbinsX() +1 ; i++)
	for(int j=0;j<= _resp_->GetNbinsY() +1 ; j++)
	{
		if (_resp_->GetBinContent(i,j) <= 1.e-10 ) {
			_resp_->SetBinError(i,j,smallestError);

		}
			
	}


	if (VERBOSE>1) cout <<"[BootStrapMatrix]::[matrixSmear]::[2] Smear "<<_bkg_<<"|"<<_gen_notreco_ <<"|"<<_resp_<<endl;
	if (not isBkgEmpty) Smear( _bkg_ );
	if (not isEffEmpty) Smear( _gen_notreco_ );

	Smear( _resp_ ) ;
	
	if (VERBOSE>1) cout <<"[BootStrapMatrix]::[matrixSmear]::[2] Compute projections"<<endl;

	_truth_= _resp_ -> ProjectionY(); _truth_->SetName( Form("%s_smearmatrix", u_truth_->GetName() ));
	_truth_->Add(_gen_notreco_);

	_reco_=_resp_ -> ProjectionX() ; _reco_->SetName( Form("%s_smearmatrix", u_reco_->GetName() ));
	_reco_->Add( _bkg_ ) ;


	if (VERBOSE>1) cout <<"[BootStrapMatrix]::[matrixSmear]::[2] Swap Pointers"<<endl;
	// swap pointers, unfold the data
	swapPointers( _reco_, u_reco_ );	
	swapPointers( _truth_, u_truth_ );	
	swapPointers( _resp_, u_resp_ );	

	if (VERBOSE>1) cout <<"[BootStrapMatrix]::[matrixSmear]::[2] Unfold:"<<u_reco_<<"|"<<u_truth_<<"|"<<u_resp_<<"|"<<data_<<endl;
	TH1D* toy = Unfold( data_ );
	// get back the unsmeared matrixes
	//
	swapPointers( _reco_, u_reco_ );	
	swapPointers( _truth_, u_truth_ );	
	swapPointers( _resp_, u_resp_ );	

	destroyPointer ( _reco_);
	destroyPointer ( _truth_);
	destroyPointer ( _resp_);
	// additional pointers
	destroyPointer ( _bkg_);
	destroyPointer ( _gen_notreco_);
	
	if (VERBOSE>1) cout <<"[BootStrapMatrix]::[matrixSmear]::[2] DONE"<<endl;
	return toy;
}


int BootStrapMatrix::CorrectNegative(TH1D* reco, TH1D*truth,TH2D*resp)
{
	// This is call after set matrixes, in order to correct negative weights as in amc@nlo with different solutions.
	//
	if ( negCorr == kNegNone )  return 0; // avoid loops

	int R=0;
	for(int t = 0 ; t<= truth->GetNbinsY()+1 ;++t)
	for(int r = 0 ; r<= reco->GetNbinsX()+1; ++r)
	{
	float c  = resp -> GetBinContent(r,t);	
	float rc = reco -> GetBinContent(r);
	float tc = truth -> GetBinContent(t);

	int shift=0;
  	if (r > t) shift = -1; // move r towards t
  	if (r < t) shift =  1;

	if (c < 0)
	{
		R=1;
		switch (negCorr) {
			case kNegNone: {R=0; break;} // is the only case where the negative error is not corrected.
			case kNegZero: { resp->SetBinContent(r,t,0); break; } 
			case kNegZeroProp:{
					  resp->SetBinContent(r,t,0);
					  truth->SetBinContent(t, tc - c ); // we are removing a content of c, so we need to remove c, or add |c|
					  reco ->SetBinContent(r, rc - c );
					  break;
					  }
			case kNegMoveProp:{
					  if (shift==0){cout <<"[ERROR] could not absorb negative correction in the diagonal element, abort. Inconsistent result."<<endl; return 0;}
					  resp->SetBinContent(r,t,0); // zero actual component
					  reco ->SetBinContent(r, rc - c );
					  //truth->SetBinContent(t,tc -c);
					  float c2= resp->GetBinContent(r+shift,t);
					  resp->SetBinContent(r + shift,t, c2 + c); // move towards the center
					  float rc2 = reco->GetBinContent(r+shift);
					  reco ->SetBinContent(r + shift, rc2 + c );
					  break;
				       	  }
			case kNegReplProp:{
					  resp->SetBinContent(r,t, -c); // zero actual component
					  reco ->SetBinContent(r, rc - 2*c );
					  //truth->SetBinContent(t,tc -c);
					  float c2= resp->GetBinContent(r+shift,t);
					  resp->SetBinContent(r + shift,t, c2 + 2*c); // move towards the center
					  float rc2 = reco->GetBinContent(r+shift);
					  reco ->SetBinContent(r + shift, rc2 + 2*c );
					  break;
					  }
		}
	} // if content is negative
	} // loop over all the bins

	// now check Reco and truth histograms. This is bad!
	for(int t = 0 ; t<= truth->GetNbinsY()+1 ;++t)
	{
		float tc = truth -> GetBinContent(t);
		if (tc <0 ) 
		{
		switch (negCorr) {
			case kNegNone: { R=0; break; }
			case kNegZero: { truth->SetBinContent(t,0); break; }
			// propagation methods will not work  print out general error and exit
			default : {
				  cout<<"ERROR negative integral in truth bin "<<t<<". Setting to 0."<<endl;
				  truth->SetBinContent(t,0); 
				  }
		} //end switch
		} // end if neg 
	} //end for over truth bin

	for(int r = 0 ; r<= reco->GetNbinsX()+1; ++r)
	{
	float rc = reco -> GetBinContent(r);
		if (rc <0){ //this is bad  -> bkg is expected to be positive
		switch (negCorr) {
			case kNegNone: { R=0; break; }
			case kNegZero: { reco->SetBinContent(r,0); break; }
			// propagation methods will not work  print out general error and exit
			default : {
				  cout<<"ERROR negative integral in reco bin "<<r<<". Setting to 0."<<endl;
				  reco->SetBinContent(r,0); 
				  }
		} // end switch
		} //end if neg
	} //end for over reco bins

	if(R>0) return CorrectNegative(reco,truth,resp); // make sure there is nothing left to correct.
	return 0;
}
