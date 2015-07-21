#include "interface/utils.hpp"
#include <iostream>

#include "TROOT.h"
#include "TStyle.h"
#include "TColor.h"
#include "TH1F.h"
#include "TH1D.h"

using namespace std;


void utils::RGB::SetRGB(int color) {
	gROOT->GetColor(color)->GetRGB(r,g,b); 
}


void utils::ChangePalette(int type)
{
	RGB white (kWhite);
	RGB darkRed (kRed+2);
	RGB red  (kRed);
	RGB black  (kBlack);
	RGB orange (kOrange);
	RGB yellow (kYellow);
	RGB darkBlue(kBlue+2);
	RGB blue (kBlue);
	RGB cyan (kCyan);
	switch (type)
		{
		case 0:{ // darkRed -> White -> Blue
	        	double r[] = { darkRed.r, red.r, orange.r,  yellow.r,      white.r,        cyan.r, blue.r, darkBlue.r, black.r };
	        	double g[] = { darkRed.g, red.g, orange.g,  yellow.g,      white.g,        cyan.g, blue.g, darkBlue.g, black.g };
	        	double b[] = { darkRed.b, red.b, orange.b,  yellow.b,      white.b,        cyan.b, blue.b, darkBlue.b, black.b };
	        	double stops[] = {0.,0.05,0.10, .30, .5, .7,.9,.95, 1.0 };
	        	Int_t FI = TColor::CreateGradientColorTable(9, stops, r, g, b, 255);
	        	gStyle->SetNumberContours(99);
			};break; // 
		case 1:{ // darkBlue -> White -> darkRed
	        	double r[] = { black.r, darkBlue.r, blue.r, cyan.r,  white.r, yellow.r, orange.r, red.r, darkRed.r };
	        	double g[] = { black.g, darkBlue.g, blue.g, cyan.g,  white.g, yellow.g, orange.g, red.g, darkRed.g };
	        	double b[] = { black.b, darkBlue.b, blue.b, cyan.b,  white.b, yellow.b, orange.b, red.b, darkRed.b };
	        	double stops[] = {0.,0.05, 0.10, .30, .5, .7, .9, .95, 1.0 };
	
	        	Int_t FI = TColor::CreateGradientColorTable(9, stops, r, g, b, 255);
	        	gStyle->SetNumberContours(99);
			};break; // 
	
		case 2:{
				// GAMMAGAMMA
			  double stops[] = {0.00, 1.00};
			  double red  [] = {1.00, 48./255.};
			  double green[] = {1.00, 48./255.};
			  double blue [] = {1.00, 131./255.};
	        	  Int_t FI = TColor::CreateGradientColorTable(2, stops, red, green, blue, 255);
	        	  gStyle->SetNumberContours(255);
			}; break; //gg
	
		}

}


// TGraphAsymmErrors* utils::Shift(TH1D* h,float dx,bool fraction=false);
TGraphAsymmErrors* utils::Shift(TGraphAsymmErrors* g,float dx,bool fraction)
{
	if (fraction == true ){
		cout <<"warning: fraction for tGraph works only for constant binnings"<<endl;
		dx *=  (g->GetX()[1] - g->GetX()[0]);
	}

	TGraphAsymmErrors *g_new = new TGraphAsymmErrors();
	g_new->SetName(Form("%s_shift",g->GetName() ) );

	for(int i=0;i < g->GetN() ;++i)
	{
		double x = g->GetX()[i];
		double y = g->GetY()[i];
		double eyl = g->GetEYlow()[i];
		double eyh = g->GetEYhigh()[i];
		
		// this macro is for graphical purpose if ex is meaningless inside the bins
		g_new->SetPoint(i, x + dx, y);
		g_new ->SetPointError(i, 0, 0, eyl,eyh ) ;
	}
	return g_new;
}

TGraphAsymmErrors * utils::Ratio(TGraphAsymmErrors *g, TH1*base)
{
	TGraphAsymmErrors *ratio = (TGraphAsymmErrors*)g->Clone(Form("%s_ratio",g->GetName()) ) ;
	if (g->GetN() != base->GetNbinsX() ) 
	{
		cout<<"Error: ratio bins should match"<<endl;
	}
	for(int i=0;i< g->GetN();++i)
	{
		double x= g->GetX()[i];	
		double y= g->GetY()[i];	
		double eyl= g->GetEYlow()[i];	
		double eyh= g->GetEYhigh()[i];	
		double b = base->GetBinContent(i+1); // only content. shift by one.

		ratio->SetPoint(i, x , y/b);
		ratio->SetPointError(i, 0,0 , eyl/b,eyh/b);
	}
	return ratio;
}

TH1* utils::Ratio(TH1* h, TH1*base){
	TH1* ratio= (TH1*) h->Clone(Form("%s_ratio",h->GetName()));
	TH1* base2= (TH1*) base->Clone(Form("%s_base",base->GetName()));
	// fix error propagation
	for(int i=0;i<base2->GetNbinsX()+1 ;++i)
		base2->SetBinError(i,0);
	
	ratio->Divide(base2);
	base2->Delete();
	
	return ratio;
}
