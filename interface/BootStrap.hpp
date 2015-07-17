#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

#include "interface/BootStrapBase.hpp"
#include "RooUnfoldResponse.h"

// --- this class implements the rooUnfold unfolding
class BootStrap : virtual public BootStrapBase
{
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

public:
	enum UnfoldType	{kBayes=0,kSVD};
private:
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
};

#endif
