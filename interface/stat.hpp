#ifndef STAT_H
#define STAT_H

#include <vector>
#include <map>
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
using namespace std;

class STAT{

public:
	/* Statistical tools build up on vectors and maps
	 * Class can be imported in a ROOT dictionary easier than namespaces
	 */

	// return the max of the given vector
	static float max(vector<float> &a);
	// return the min of the given vector
	static float min(vector<float> &a);
	// return the mean of a vector
	static float mean(vector<float> &a );
	// return the weighted mean of a vector
	static float mean(vector<float> &a ,vector<float>&e_a);
	// return the median of a vector
	static float median(vector<float> &a);
	// return the rms of a vector
	static float rms(vector<float> &a);
	// return the Pearson correlation coefficient
	static float corrPearson(vector<float> &a, vector<float> &b);
	static float corrSpearman(vector<float> &a ,vector<float> &b);
	//return the linear regression coefficients (LSE) of the points
	static pair<float,float> regression(vector<float>&a,vector<float>&b);
	//error 0 = <0,0> 1=<1,1> 2=<0,1>
	//return the linear regression coefficient of the points LSE
	//and the covariance matrix
	static pair<float,float> regression(vector<float>&a,vector<float>&b, vector<float>&e_b,vector<float> &e2);


	//return in r, the smallest interval to contain the fraction Q of the elements in vector	
	static float ConfidenceInterval(std::vector<float> &v,std::pair<float,float>&r,float Q);

	// --- half confidence below the given value, half above
	static float ConfidenceIntervalAround(std::vector<float> &v,float value, std::pair<float,float>&r, float Q);

	// --- chi2, (x-y)^2/e^2 iff elow is empty, error are on b
	//           (x-y)^2/elow^2 iff x<y  or ehigh iff x>y
	//           if corr is not empty use also correlation matrix
	static float Chi2( vector<float> &a, vector<float> &b, vector<float> &ehigh, vector<float> &elow ,map<pair<int,int>,float> &corr);
	//  ---------------- ROOT -------------
	// draw fit errors on the histogram h (already existnig ).
	// e2 is the covariance matrix, as above.  (regression)
	static void drawFitError(TH1*h,pair<float,float> &R,vector<float> &e2,float sigma=1);

	//Fill Histo
	static void Fill(std::vector<float> &v, TH1*h);
	static void Fill(std::vector<float> &a,std::vector<float> &b, TH2*h);
	// Get Density Histogram, Gauss with radius R
	static TH1F* GetDensity(std::vector<float> &v, float R) ;

	static float Chi2(TGraphAsymmErrors *g, TH1 *h, TH2* corr=NULL);
};

#endif
