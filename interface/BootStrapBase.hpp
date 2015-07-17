#ifndef BOOTSTRAP_BASE_H
#define BOOTSTRAP_BASE_H

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

protected:
	// ------------- templates goes in the h file
	template<class T>
	void setPointer( T orig,T &target){
		destroyPointer(target);
		target=orig;
	}

	template<class T>
	T releasePointer(T &target)
	{
		T r=target;
		target=NULL;
		return r;
	}

	template<class T>
	void destroyPointer(T &ptr){
		if (ptr != NULL) delete ptr;
		ptr=NULL;
	}

	void clearBootstrap(){ for( TH1D*ptr : bootstrap_ ) destroyPointer( ptr) ; bootstrap_.clear(); }
	TH1D* bootStrap();

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

	// run BootStrap
	void run();	

	enum ResultType { kStd=0, kMin=1, kMedian=2 , kMean = 3 };
	// get results -- these are recomputed each time
	// type =0 : 
	// 	mean points are the unfold ones. errors are asymetric in order to cover Q/2 each.
	// type =1 :
	// 	mean points are the unfold ones. errors are the minimum interval to cover Q
	// type =2 :
	// 	errors are the median of the toys and errors covers the small possible interval
	// type = 3:
	// 	points are the mean of toys and errors covers the small possible interval
	TGraphAsymmErrors *result(ResultType type=kStd,float Q=0.68);
	TH2D *correlation();
};


#endif
