#ifndef BOOTSTRAP_BASE_H
#define BOOTSTRAP_BASE_H

/* Original Author: Andrea C. Marini
 * Date: 16 Jul 2015
 * No Warranty
 */

#include "TH1D.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TGraphAsymmErrors.h"

#include <vector>
using namespace std;

class BootStrapBase
{
private:
	// number of toys iterations
	int Ntoys_;
	// use gauss vs Poisson distributions
	bool SumW2_;
	// Data
	TH1D * data_;
	// BootStrap-ed distributions
	vector<TH1D*> bootstrap_;
	//
	TRandom *r_;
	long seed_;
	
	//  folded and unfolded histograms 
	//  saved for saving time
	TH1D* unf_;
	TH1D* fold_;

	int confCounter_;
	int confSigma_ ;    // how far the toys are accepted
	int confSigmaGen_ ; // how far the toys are generated -- only for SumW2
protected:
	// ------------- templates goes in the h file
	template<class T> void setPointer( T orig,T &target);

	template<class T> T releasePointer(T &target);

	template<class T> void destroyPointer(T &ptr);


	void clearBootstrap(){ for( TH1D*ptr : bootstrap_ ) destroyPointer( ptr) ; bootstrap_.clear(); }
	TH1D* bootStrap();

	TH1D* confidence();

public:
	// --- Constructor 
	BootStrapBase();
	BootStrapBase(int ntoys);
	// --- Destructor
	~BootStrapBase();

	// --- Virtual  members that do the folding / Unfolding
	virtual TH1D* Unfold(TH1D*) = 0;
	virtual TH1D* Fold(TH1D*) = 0 ;

	// member functions
	// set number of toys
	void SetNToys(int ntoys=100);
	// set if use poisson or gaus for toys
	void SetSumW2(bool sumw2=true);

	// set and take ownership of the data histogram
	void SetData( TH1D* data); 
	// return and give back ownership
	TH1D* releaseData(); 

	enum RunType{kBootstrap=0,kConfidence=1};
	// run BootStrap
	void run(RunType type=kBootstrap);	

	//--
	inline void SetConfSigma(int sigma){confSigma_ =sigma;};
	inline void SetConfSigmaGen(int sigma){confSigmaGen_ =sigma;};
	inline void SetSeed(long seed){seed_ = seed;};


	enum ResultType { kStd=0, kMin=1, kMedian=2 , kMean = 3 };
	// get results -- these are recomputed each time
	// type =kStd : 
	// 	mean points are the unfold ones. errors are asymetric in order to cover Q/2 each.
	// type =kMin :
	// 	mean points are the unfold ones. errors are the minimum interval to cover Q
	// type = kMedian :
	// 	mean points are the median of the toys and errors covers the smallest possible interval
	// type = kMean:
	// 	points are the mean of toys and errors covers the smallest possible interval
	TGraphAsymmErrors *result(ResultType type=kStd,float Q=0.68);
	TH2D *correlation();
};

// -- TEMPLATE DEFINITIONS
template<> void BootStrapBase::destroyPointer<TObject*> (TObject* &ptr);

template<class T>
void BootStrapBase::setPointer( T orig,T &target){
	destroyPointer(target);
	target=orig;
}

template<class T>
T BootStrapBase::releasePointer(T &target)
{
	T r=target;
	target=NULL;
	return r;
}

template<class T>
void BootStrapBase::destroyPointer(T &ptr){
	if (ptr != NULL) delete ptr;
	ptr=NULL;
}


#endif
