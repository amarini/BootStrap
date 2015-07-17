#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

/* Original Author: Andrea C. Marini
 * Date: 16 Jul 2015
 * No Warranty
 */

#include "interface/BootStrapBase.hpp"
#include "RooUnfoldResponse.h"

// --- this class implements the rooUnfold unfolding
class BootStrap : virtual public BootStrapBase
{
public:
	enum UnfoldType	{kBayes=0,kSvd=1,kInv=2};

private:
	// --- matrix -- 
	// in principle can be different ... needs studies (related to model dependence)
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
	TH1D * f_resp_eff_;
	TH1D * f_resp_smear_;

	int regParam_;
	UnfoldType unfType_;	

public:
	BootStrap() ;
	BootStrap(int ntoys) : BootStrapBase(){ SetNToys(ntoys);};
	~BootStrap() ;

	//
	TH1D* Unfold(TH1D*);
	TH1D* Fold(TH1D*);

	// SetMatrixes: Unfolding and or Folding. If Folding is NULL or unset, will use the same matrix.
	void SetUMatrix(TH1D* reco, TH1D* truth, TH2D* resp);
	void SetFMatrix(TH1D* reco, TH1D* truth, TH2D* resp);

	inline void SetUnfoldType( UnfoldType t){ unfType_ = t;}
};

#endif
