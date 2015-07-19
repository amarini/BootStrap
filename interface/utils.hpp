#ifndef UTILS_H
#define UTILS_H

#include "TGraphAsymmErrors.h"



class utils
{
public:
	static void ChangePalette(int type=0);
	//static TGraphAsymmErrors* Shift(TH1D* h,float dx,bool fraction=false);
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
};

#endif
