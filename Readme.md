# BootStrap

### Table Of Contents
- [Installation](#user-content-installation)
- [User Guide](#user-content-user-guide)
- [Example](#user-content-example)
- [Code Organization](#user-content-code-organization)

## Installation
This package needs the RooUnfold package, and a C++11 compiler.
It has been tested with ROOT6.
* Edit the Makefile in order to point the RooUnfold libraries and includes
```
make -j 16
```

## User Guide

* Construct a BootStrap method. Only BootStrap based on RooUnfold is doable at the moment.

The following functions are available:
* SetNToys(): Set the number of BootStrap Toys that will be throw.
* SetSeed(): Set the seed of the Random Number Generator.
* SetUnfoldType(): BootStrap::kBayes, BootStrap::kSvd, BootStrap::kInv
* SetRegParam(): Set the regularization parameter (only for svd, and bayes)
* SetUMatrix(reco,gen,resp): set the response matrix histograms.
* SetFMatrix(): if done, set a different matrix for the folding procedure.
* SetData(): Set Data to be unfolded
* SetSumW2(): Generate toys accordingly to Poisson or to Gaus. (since unfolding is run w/o errors, gaus will have err=sqrt(content) ).
* run( type ): Run the boostrap. Type can be:
    * BootStrapBase::kBootstrap : for the base Bootstrap
    * BootStrapBase::kConfidence: for the search of the confidence interval. in this case the additional parameters:
        * SetConfSigma(): for the number of sigma a toy is accepted
        * SetConfSigmaGen(): (iff SumW2 on) for the generation of toys in the truth space.

* results ( type ): get the result as TGraphAsymmError
     * BootStrapBase::kStd: mean points are the unfold ones. errors are asymetric in order to cover Q/2 each.
     * kMin: mean points are the unfold ones. errors are the minimum interval to cover Q
     * kMedian: mean points are the median of the toys and errors covers the smallest possible interval
     * kMean: mean  points are the mean of toys and errors covers the smallest possible interval

* correlation() : return the correlation matrix. elements are the pearson correlatino coefficient of the toys.

## Example
* Open python (PyROOT) and import libraries:
```
python
import ROOT
ROOT.gSystem.Load("bin/libBootStrap.so")
```
* Get Response matrix/data and feed them to the BootStrap:
```
reco = f.Get(...) / data = f.Get() ...
b= ROOT.BootStrap()
b.SetNToys(1000)
b.SetUnfoldType(ROOT.BootStrap.kBayes)
b.SetRegParam(5)
b.SetUMatrix(reg,gen,resp)
b.SetData(data)
``` 
* Run the boostrap
```
b.run() ## default is ROOT.BootStrapBase.kBootstrap
```
* Get the results
```
g_bootstrap = b.result(ROOT.BootStrapBase.kMedian,.68)
```

## Code Organization
Two classes has been implemnted so far:

**STAT**
* implements some statistical tools that will be used
* it''s a class (and not a namespace) because I had some problem in linking it in ROOT/PyROOT

**BootStrapBase**
* pure virtual (missing the implementation of the Fold and Unfold methods)
* implements the bootstrapping techinques
* run() call different methos accordingly to the production of the toys: bootStrap(), confidence().


**BootStrap**
* implements the Folding and Unfolding methods using RooUnfold.
* the Folding method assumes that matrix performs a smearing, and sums a bkg contribution. No NULL element is supported for the moment.