#ifndef BOOTSTRAP_MATRIX_H
#define BOOTSTRAP_MATRIX_H

/* Original Author: Andrea C. Marini
 * Date: 16 Jul 2015
 * No Warranty
 */

#include "interface/BootStrapBase.hpp"
#include "RooUnfoldResponse.h"

#include "TVectorD.h"
#include "TMatrixD.h"

// --- This class implements unfolding with matrix a la RooUnfold.
// -- inversion with pesudo inverse
class BootStrapMatrix : public BootStrapBase
{
protected: // RooUnfold will need access to these elements
	// --- matrix -- 
	// in principle can be different ... 
	// needs studies (related to model dependence)
	TH1D * u_reco_;
	TH1D * u_truth_;
	TH2D * u_resp_;
	// ------
	TH1D * f_reco_;
	TH1D * f_truth_;
	TH2D * f_resp_;

	// --- projections
	//
	TH1D * f_resp_bkg_;
	TH2D * f_resp_smear_;

	// --- projections2
	TH1D * u_resp_bkg_;
	TH2D * u_resp_smear_;

	virtual TH1D* matrixSmear(); 
	void ConstructProjections(TH1D*reco,TH1D*truth,TH2D*resp);
	void ConstructProjections(TH1D*reco,TH1D*truth,TH2D*resp, TH1D* &bkg, TH2D* &smear);

private:
	// matrix used to perform the inversion
	TMatrixD K,S; //y = K * x + b
	TMatrixD Kt, A, B, Bt, cov;
		// solution is (Kt*S*K+eI)^-1 Kt S
	TVectorD y,b,l;

	// avoid to reconstruct matrixes each time
	bool matrixConstructed_;
	void ConstructMatrixes(TH1D*data);
	TMatrixD getMatrix(TH2*h, bool useOverFlow=true);
	TVectorD getVector(TH1*h, bool useOverFlow=true);
	void printMatrix(TMatrixD&,string name="");

	// Correct the response matrix if negative. See corrections below;
	// return 0; iff nothing has been corrected
	int CorrectNegative(TH1D*reco,TH1D* truth, TH2D* resp);
public:
	// constructor 
	BootStrapMatrix();
	//copy constructor
	BootStrapMatrix( BootStrapMatrix &x );
	~BootStrapMatrix();

	// destroy also data_scaled_	
	// -- void SetData( TH1D* data); 

	// implemetn the Fold
	virtual TH1D* Unfold(TH1D*);
	virtual TH1D* Fold(TH1D*);

	// SetMatrixes: Unfolding and or Folding. If Folding is NULL or unset, will use the same matrix.
	void SetUMatrix(TH1D* reco, TH1D* truth, TH2D* resp);
	void SetFMatrix(TH1D* reco, TH1D* truth, TH2D* resp);

	// Get Matrixes, plot /debug after neg corrections
	TH2D* GetUMatrixResp(){return u_resp_;}
	TH1D* GetUMatrixReco(){return u_reco_;}
	TH1D* GetUMatrixTruth(){return u_truth_;}

	enum NegativeCorrections{ 
		kNegNone = 0,  // Do Nothing and pray for the best
		kNegZero,  // Zero all the negative element.
		kNegZeroProp, //Zero all the negative element and propagate the differences.
		kNegMoveProp, // Move the 0 bins in the closest diagonal and propagate
		kNegReplProp, // Replace the negative with a positive and counterbalance in the diagonal
		};
	NegativeCorrections negCorr=kNegNone;
	void info() ; // override

	inline void correct() { // override
		BootStrapBase::correct();	

		if ( u_reco_ && u_truth_ && u_resp_) CorrectNegative(u_reco_,u_truth_,u_resp_);
		if ( f_reco_ && f_truth_ && f_resp_) CorrectNegative(f_reco_,f_truth_,f_resp_);
	};

};
#endif
