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
public:
	enum ToyType{ kBootstrap = 0 , kToy = 1 , kIterBias=2 , kMatrix=3};
	// kMatrix will smear the matrix
private:
	// number of toys iterations
	int Ntoys_;
	// BootStrap-ed distributions
	vector<TH1D*> bootstrap_;
	//
	

	// save bias for iterative bias corrections
	TH1D *bias_;

	int verbose_;
	ToyType type_;
protected:
	// ------------- templates goes in the h file
	template<class T> void setPointer( T orig,T &target);

	template<class T> T releasePointer(T &target);

	template<class T> void destroyPointer(T &ptr);

	template<class T> inline void swapPointers( T* &x, T* &y ){ T* z=x; x=y; y=z ; }


	void clearBootstrap(){ for( TH1D*ptr : bootstrap_ ) destroyPointer( ptr) ; bootstrap_.clear(); }
	TH1D* bootStrap();
	TH1D* directToy();

	virtual TH1D* matrixSmear()=0; // from the smearing of the matrix

	// random number generation	
	TRandom *r_;
	long seed_;

	// Data
	TH1D * data_;

	//  folded and unfolded histograms 
	//  saved for saving time
	TH1D* unf_;
	TH1D* fold_;

	// use gauss vs Poisson distributions
	bool SumW2_;

	int Nib_; // number of iterative Bias corrections
	int Ntoysib_;
	TH1D* iterativeBias(bool toy=false); // compute the iterative bias corrections
	vector<TH1D*> ibtoys_;
	void clearIbtoys();
	
	void Smear(TH1*); // smear using poisson or SumW2
	void Smear(TH2*); // smear using poisson or SumW2

public:
	// --- Constructor 
	BootStrapBase();
	// -- Copy Constructor
	BootStrapBase(BootStrapBase &);
	// --- Destructor
	virtual ~BootStrapBase();

	// --- Virtual  members that do the folding / Unfolding
	virtual TH1D* Unfold(TH1D*) = 0;
	virtual TH1D* Fold(TH1D*) = 0 ;

	// member functions
	// set number of toys
	void SetNToys(int ntoys=100);
	// set if use poisson or gaus for toys
	void SetSumW2(bool sumw2=true);

	// set and take ownership of the data histogram
	virtual void SetData( TH1D* data); 
	// return and give back ownership
	TH1D* releaseData(); 

	// run BootStrap
	void run();	

	// call correction before running
	virtual inline void correct(){};

	//--
	inline void SetSeed(long seed){seed_ = seed;};
	//--
	inline void SetToyType(ToyType t){type_=t;};
	inline ToyType GetToyType() const { return type_;}
	//--
	inline void SetNiterBias(int n){Nib_ = n;}



	enum ResultType { kStd=0, kMin=1, kMedian=2 , kMean = 3, kRms=4 };
	// get results -- these are recomputed each time
	// type =kStd : 
	// 	mean points are the unfold ones. errors are asymetric in order to cover Q/2 each.
	// type = kRms:
	// 	mean point are the unfold ones. errors is the rms. symmetric
	// type =kMin :
	// 	mean points are the unfold ones. errors are the minimum interval to cover Q
	// type = kMedian :
	// 	mean points are the median of the toys and errors covers the smallest possible interval
	// type = kMean:
	// 	points are the mean of toys and errors covers the smallest possible interval
	TGraphAsymmErrors *result(ResultType type=kStd,float Q=0.68);
	TH2D *correlation();
	virtual void info(); // print info
	inline void SetVerbose(int v){verbose_=v;}
	
	// Get Toy Distribution for bin 
	TH1F* GetToyDistribution(int bin);
	TH2F* GetToyDistribution(int bin1,int bin2);
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
