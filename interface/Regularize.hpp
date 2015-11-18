#include "TH1D.h"
#include "TH2D.h"

#include "TMatrixD.h"
#include "TVectorD.h"

class Regularize
{
	/* This macros is designed to regularized a distribution
	 * It solves the problem F = || ( h-v ) S^{-1} (h-v) || + \sqrt{delta} || K * v || ^2
	 * solution is:
	 *       A =delta*K*K + S^-1
	 *       B = 1 - A^{-1}*K*K*delta
	 *       v =  B h
	 * if regularization is doven wrt v/mc, should be enough to rescale h, and with respect mc, and rescale the solution (v) after
	 * and returns the histogram
	 */
	private:

		// tools
		// this macros need to know about useOverFlow
		TMatrixD getMatrix(TH2 *h);
		TVectorD getVector(TH1 *h);
		TH2D*	 getCovMatrixH(TH1 *h);
		TMatrixD getCovMatrix(TH1 *h);
		// tools
		TVectorD integrateRow(TMatrixD &a);
		TVectorD integrateCol(TMatrixD &a);
		void print(TVectorD &v);
		void print(TMatrixD &v);
	public:
		Regularize(){};
		~Regularize(){};
		bool useOverFlow = false;
		TH1D* applyRegularization(TH1D* inv, TH2D *cov, float delta=0.0);

};
