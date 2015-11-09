#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

/* Original Author: Andrea C. Marini
 * Date: 16 Jul 2015
 * No Warranty
 */

#include "interface/BootStrapBase.hpp"
#include "interface/BootStrapMatrix.hpp"
#include "RooUnfoldResponse.h"

// --- this class implements the rooUnfold unfolding
class BootStrap :  public BootStrapMatrix
{
public:
	enum UnfoldType	{kBayes=0,kSvd=1,kInv=2};

private:

	int regParam_;
	UnfoldType unfType_;	

public:
	BootStrap() ;
	BootStrap( BootStrap &);
	~BootStrap() ;

	//
	virtual TH1D* Unfold(TH1D*);
	// in case we need the inversion unfolding
	virtual TH1D* UnfoldLikelihood(TH1D*h) ;


	inline void SetUnfoldType( UnfoldType t){ unfType_ = t;}
	inline void SetRegParam( int reg) { regParam_ = reg;}

	virtual void info();

};

#endif
