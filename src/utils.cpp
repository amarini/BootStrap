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
		case 3:{ //  White -> darkRed
	        	double r[] = { white.r, yellow.r, orange.r, red.r, darkRed.r ,black.r};
	        	double g[] = { white.g, yellow.g, orange.g, red.g, darkRed.g ,black.g};
	        	double b[] = { white.b, yellow.b, orange.b, red.b, darkRed.b ,black.b};
	        	double stops[] = {0.,0.020, 0.45, .7, .8, .95, 1.0 };
	
	        	Int_t FI = TColor::CreateGradientColorTable(7, stops, r, g, b, 255);
	        	gStyle->SetNumberContours(99);
			};break; // 
	
		}

}


TGraphAsymmErrors* utils::Shift(TH1D* h,float dx,bool fraction){
	double N= h->GetNbinsX();
	TGraphAsymmErrors *g = new TGraphAsymmErrors();
	g->SetName(Form("%s_shift",h->GetName()));
	g->SetTitle(h->GetTitle());
	for(int i=1;i <= N ;++i)	
	{
		double x=h->GetBinCenter(i);
		double w=h->GetBinWidth(i);
		double y=h->GetBinContent(i);
		double e=h->GetBinError(i);
		
		double x1=x+dx;
		if (fraction) x1 = x + dx*w ;
		int n=g->GetN();
		g->SetPoint(n, x1,y);
		g->SetPointError(n,0,0,e,e);
	}
	CloneStyle(h,g);
	return g;
}


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
	CloneStyle(g,g_new);
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

TGraphAsymmErrors * utils::Ratio(TGraphAsymmErrors *g, TF1*base)
{
	cout<<"Ratio TGraph/TF1 warning... assume constant binning"<<endl;
	TGraphAsymmErrors *ratio = (TGraphAsymmErrors*)g->Clone(Form("%s_ratio",g->GetName()) ) ;

	double width= g->GetX()[1] - g->GetX()[0];
	for(int i=0;i< g->GetN();++i)
	{
		double x= g->GetX()[i];	
		double y= g->GetY()[i];	
		double eyl= g->GetEYlow()[i];	
		double eyh= g->GetEYhigh()[i];	
		double b = base->Integral(x-width/2.,x+width/2.); // only content. shift by one.

		ratio->SetPoint(i, x , y/b);
		ratio->SetPointError(i, 0,0 , eyl/b,eyh/b);
	}
	return ratio;
}

TH1* utils::Ratio(TH1* h, TF1*base){
	TH1* ratio= (TH1*) h->Clone(Form("%s_ratio",h->GetName()));

	for(int i=0;i<h->GetNbinsX()+1 ;++i)
	{
		double xl = h->GetBinLowEdge(i);
		double xh = h->GetBinLowEdge(i+1);
		double b = base->Integral(xl,xh);
		double c = h->GetBinContent(i);
		double e = h->GetBinError(i);
		ratio->SetBinContent(i, c/b);
		ratio->SetBinError(i, e/b);
	}
	
	return ratio;
}


void utils::CloneStyle(TObject *base, TObject *target)
{
	if(   base->InheritsFrom("TAttLine") and 
	      target->InheritsFrom("TAttLine")) {
		dynamic_cast<TAttLine*>(target)->SetLineColor( dynamic_cast<TAttLine*>(base)->GetLineColor() );
		dynamic_cast<TAttLine*>(target)->SetLineStyle( dynamic_cast<TAttLine*>(base)->GetLineStyle() );
		dynamic_cast<TAttLine*>(target)->SetLineWidth( dynamic_cast<TAttLine*>(base)->GetLineWidth() );
	}	
	if(   base->InheritsFrom("TAttMarker") and 
	      target->InheritsFrom("TAttMarker")) {
		dynamic_cast<TAttMarker*>(target)->SetMarkerColor( dynamic_cast<TAttMarker*>(base)->GetMarkerColor() );
		dynamic_cast<TAttMarker*>(target)->SetMarkerStyle( dynamic_cast<TAttMarker*>(base)->GetMarkerStyle() );
		dynamic_cast<TAttMarker*>(target)->SetMarkerSize( dynamic_cast<TAttMarker*>(base)->GetMarkerSize() );
	}	
	if(   base->InheritsFrom("TAttFill") and 
	      target->InheritsFrom("TAttFill")) {
		dynamic_cast<TAttFill*>(target)->SetFillColor( dynamic_cast<TAttFill*>(base)->GetFillColor() );
		dynamic_cast<TAttFill*>(target)->SetFillStyle( dynamic_cast<TAttFill*>(base)->GetFillStyle() );
	}	
	if(   base->InheritsFrom("TAttText") and 
	      target->InheritsFrom("TAttText")) {
		dynamic_cast<TAttText*>(target)->SetTextColor( dynamic_cast<TAttText*>(base)->GetTextColor() );
		dynamic_cast<TAttText*>(target)->SetTextAngle( dynamic_cast<TAttText*>(base)->GetTextAngle() );
		dynamic_cast<TAttText*>(target)->SetTextAlign( dynamic_cast<TAttText*>(base)->GetTextAlign() );
		dynamic_cast<TAttText*>(target)->SetTextFont( dynamic_cast<TAttText*>(base)->GetTextFont() );
		dynamic_cast<TAttText*>(target)->SetTextSize( dynamic_cast<TAttText*>(base)->GetTextSize() );
	}	
}
