#ifndef STAT_H
#define STAT_H

#include <vector>
#include <map>
#include "TH1F.h"
using namespace std;

class STAT{

public:
	/* Statistical tools build up on vectors and maps
	 * Class can be imported in a ROOT dictionary easier than namespaces
	 */

	static float mean(vector<float> &a );
	static float mean(vector<float> &a ,vector<float>&e_a);
	static float median(vector<float> &a);
	static float rms(vector<float> &a);
	static float corrPearson(vector<float> &a, vector<float> &b);
	static float corrSpearman(vector<float> &a ,vector<float> &b);
	static pair<float,float> regression(vector<float>&a,vector<float>&b);
	//error 0 = <0,0> 1=<1,1> 2=<0,1>
	static pair<float,float> regression(vector<float>&a,vector<float>&b, vector<float>&e_b,vector<float> &e2);
	
	static void drawFitError(TH1*h,pair<float,float> &R,vector<float> &e2,float sigma=1);
	
	static float ConfidenceInterval(std::vector<float> &v,std::pair<float,float>&r,float Q);

	// --- half confidence below the given value
	static float ConfidenceIntervalAround(std::vector<float> &v,float value, std::pair<float,float>&r, float Q);
};

#endif
