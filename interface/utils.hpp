#ifndef UTILS_H
#define UTILS_H

#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TH1D.h"


class utils
{
public:
	static void ChangePalette(int type=0);
	static TGraphAsymmErrors* Shift(TH1D* h,float dx,bool fraction=false);
	static TGraphAsymmErrors* Shift(TGraphAsymmErrors* g,float dx,bool fraction=false);

	class RGB
	{
	public:
		float r;
		float g;
		float b;
		void SetRGB(int color); 
		RGB(int color){ SetRGB(color) ; } 
		RGB(){ SetRGB(0);}
	};

	static TGraphAsymmErrors* Ratio(TGraphAsymmErrors* g, TH1*base);
	static TH1* Ratio(TH1* h, TH1*base);
	static TGraphAsymmErrors* Ratio(TGraphAsymmErrors* g, TF1*base);
	static TH1* Ratio(TH1* h, TF1*base);

	static void CloneStyle(TObject *base, TObject* target );
};

#endif
