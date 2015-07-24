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
class BootStrapMatrix : virtual public BootStrapBase
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
	
	void ConstructProjections(TH1D*reco,TH1D*truth,TH2D*resp);
private:
	TMatrixD K,S; //y = K * x + b
		// solution is (Kt*S*K+eI)^-1 Kt S
	TVectorD y,b,l;
	void ConstructMatrixes(TH1D*data);
	TMatrixD getMatrix(TH2*h, bool useOverFlow=false);
	TVectorD getVector(TH1*h, bool useOverFlow=false);
public:
	// constructor 
	BootStrapMatrix();
	~BootStrapMatrix();
	
	// implemetn the Fold
	virtual TH1D* Unfold(TH1D*);
	virtual TH1D* Fold(TH1D*);

	// SetMatrixes: Unfolding and or Folding. If Folding is NULL or unset, will use the same matrix.
	void SetUMatrix(TH1D* reco, TH1D* truth, TH2D* resp);
	void SetFMatrix(TH1D* reco, TH1D* truth, TH2D* resp);

	virtual void info();

};
#endif