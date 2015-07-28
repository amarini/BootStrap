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
	TMatrixD getMatrix(TH2*h, bool useOverFlow=false);
	TVectorD getVector(TH1*h, bool useOverFlow=false);
	void printMatrix(TMatrixD&,string name="");


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

	virtual void info();

};
#endif
