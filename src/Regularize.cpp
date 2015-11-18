#include "interface/Regularize.hpp"
#include "TMath.h"
#include <iostream>
#include <assert.h>
using namespace std;

TH1D * Regularize::applyRegularization(TH1D* inv, TH2D* cov, float delta)
{
	int nBins=inv->GetNbinsX();
	float xmin = inv->GetBinLowEdge(1);
	float xmax = inv->GetBinLowEdge(nBins+1);

	if (useOverFlow)
		nBins+=2;
	
	TMatrixD K;
	K.ResizeTo(nBins,nBins);
	K.Zero();
	for(int i=0;i< nBins ; ++i)
		K(i,i) = -2;
	for(int i=0;i< nBins -1; ++i)
		 K(i,i+1) = 1;
	for(int i=0;i< nBins -1; ++i)
		 K(i+1,i) = 1;
	K(0,0) = -1;
	K(nBins-1,nBins-1)=-1;

	TVectorD h = getVector(inv);
	TMatrixD S = getMatrix(cov);
	cout <<" ORIGINAL: "<<endl;
	print(h);
	print(S);
	cout <<" <--> <--> <--> <--> "<<endl;

	TMatrixD Sinv(S); Sinv.Invert();
	// -------------------------------------
	TMatrixD A(nBins,nBins);
	A = K*K;
	A *= delta;
	A += Sinv;

	if (delta == 0 )
	{
	cout <<"---D=0---"<<endl;
	 print (Sinv);	
	 print (A);
 	cout <<"----------"<<endl;	 
	}

	// --- Invert numerically A
	TMatrixD Id(nBins,nBins);
	TMatrixD U(nBins,nBins);
	Id.Zero();
	U.Zero();
	for(int i=0;i<Id.GetNrows();++i)
	{
		Id(i,i) = 1;
		U(i,i) = 1;
	}
	/// ----
	float sum=0;
	float epsilon=0.00001;
	Id*=epsilon;
	for(int i=0; i < 100 ;++i)
	{
		if (A.Determinant() != 0) break;
		A += Id;
		sum+=epsilon;
		
	}
	if (A.Determinant() == 0 ) cout<<"Determinant is 0"<<endl;
	else if( sum>0) cout<<"Added "<<sum<<" to the diagonal to make it invertible. (As SVD)"<<endl;

	A.Invert();

	//TMatrixD B(A * delta * Kt *K  );
	TMatrixD B(nBins,nBins);
	B = A*K*K; B*= delta; 
	B = U - B;
	TMatrixD Bt(nBins,nBins) ; Bt.Transpose(B);
	TVectorD v(nBins); v= B *h;
	TMatrixD err (nBins,nBins); err =  B * S * Bt ;

	cout<<"..... RESULTS .........."<<endl;
	print (v);
	print (err);
	cout <<"........................"<<endl;
	
	TH1D *res = (TH1D*)inv->Clone("regularized");

	int i=0;
	for(int iBin=0;iBin<res->GetNbinsX() + 2; ++iBin)
	{
		if ( not useOverFlow and (iBin==0 or iBin==res->GetNbinsX() + 1)) continue;
		res->SetBinContent( iBin, v(i) );
		res->SetBinError( iBin, TMath::Sqrt(err(i,i)) );
		i++;
	}
	return res;
}

// -------------------- TOOLS ---------------------
TMatrixD Regularize::getMatrix(TH2 *h){
	TMatrixD r;
	int n=h->GetNbinsX();
	int m=h->GetNbinsY();
	if (useOverFlow)n+=2;
	if (useOverFlow)m+=2;
	r.ResizeTo(n,m);
	r.Zero(); //make sure is set to 0
	if(not useOverFlow) 
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
TMatrixD Regularize::getCovMatrix(TH1 *h){
	TMatrixD r;
	int n=h->GetNbinsX();
	if (useOverFlow)n+=2;
	r.ResizeTo(n,n); 
	r.Zero(); // make sure is set to 0
	if(!useOverFlow) 
	{
		for(int iBin=1;iBin<=h->GetNbinsX();iBin++)
			r(iBin-1,iBin-1)=TMath::Power(h->GetBinError(iBin),2);
	
	}
	else{
		for(int iBin=0;iBin<=h->GetNbinsX()+1;iBin++)
			r(iBin,iBin)=TMath::Power(h->GetBinError(iBin),2);
	}
	return r;
}

TVectorD Regularize::getVector(TH1 *h){
	TVectorD r;
	int n=h->GetNbinsX();
	if (useOverFlow)n+=2;
	r.ResizeTo(n);
	int i=0;
	if(useOverFlow) { r(i) = h->GetBinContent(i);i++;}
	for(int iBin=1;iBin<=h->GetNbinsX();iBin++)
		{
			assert(i<r.GetNrows());
			r(i)=h->GetBinContent(iBin);i++;
		}
	if(useOverFlow) { 
			assert(i<r.GetNrows());
			r(i) = h->GetBinContent(h->GetNbinsX()+1);i++;

	}
	return r;
}


TH2D*	 Regularize::getCovMatrixH(TH1 *h)
{
	int nBins=h->GetNbinsX();
	Double_t Bins[nBins+2];
	for(int iBin=1;iBin<=nBins+1;iBin++)
		Bins[iBin-1]=h->GetBinLowEdge(iBin); //offset by 1
	TH2D* R=new TH2D(Form("cov_%s",h->GetName()),h->GetName(),nBins,Bins,nBins,Bins);
	for(int iBin=0;iBin<nBins+1;iBin++)
	for(int jBin=0;jBin<nBins+1;jBin++)
		{
		R->SetBinContent(iBin+1,jBin+1,0);
		if( iBin == jBin) 
			R->SetBinContent(iBin+1,jBin+1,pow(h->GetBinError(iBin+1),2));
		}
	return R;
}

TVectorD Regularize::integrateRow(TMatrixD &a)
{
	TVectorD r;
	r.ResizeTo(a.GetNcols());
	r.Zero();
	for (int i=0;i<a.GetNrows();i++)
	for (int j=0;j<a.GetNcols();j++)
	{
		r(j)+=a(i,j);
	}
	return r;
}
TVectorD Regularize::integrateCol(TMatrixD &a)
{
	TVectorD r;
	r.ResizeTo(a.GetNrows());
	r.Zero();
	for (int i=0;i<a.GetNrows();i++)
	for (int j=0;j<a.GetNcols();j++)
	{
		r(i)+=a(i,j);
	}
	return r;
}

void Regularize::print(TVectorD &v){
	cout<<"----- VECTOR -----"<<endl;
	for(int i=0;i <v.GetNrows() ;++i)
		cout <<"* "<<i<<": "<<v(i)<<endl;
	cout<<"-------------------"<<endl;
}
void Regularize::print(TMatrixD &v)
{
	cout<<"----- MATRIX -----"<<endl;
	for(int i=0;i <v.GetNrows() ;++i)
	{
	   for(int j=0;j <v.GetNcols() ;++j)
		cout <<" "<<v(i,j);
	   cout << endl;
	}
	cout<<"-------------------"<<endl;

}
