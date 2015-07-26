#include <vector>
#include <algorithm>
#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "interface/stat.hpp"

using namespace std;

/* Macro to analyze the points*/

float STAT::mean(vector<float> &a )
	{
	float S=0;
	for(int i=0;i<int(a.size());i++) S+=a[i];
	S/=a.size();
	return S;
	}
float STAT::median(vector<float> &a)
	{
	if( not is_sorted(a.begin(),a.end() ) ) sort(a.begin(),a.end() );
	float M=0;
	//must be sorted
	int nEntries=a.size();
	if( (nEntries%2)==0)  //even
		{
		M=(a[nEntries/2]+a[nEntries/2+1] )/2.0;
		}
	else {//odd
		M=a[nEntries/2+1];
		}
	return M;
	}
float STAT::rms(vector<float> &a)
	{
	float S=0;
	float m=mean(a);
	for(int i=0;i<int(a.size());i++) S+=(a[i]-m)*(a[i]-m);
	S/=(a.size()-1);
	return sqrt(S);
	}

float STAT::corrPearson(vector<float> &a, vector<float> &b)
	{
	if(a.size() != b.size()) { printf("SIZE=%d %d\n",int(a.size()),int(b.size()));return -2.0;}
	float sa=rms(a);
	float sb=rms(b);
	float ma=mean(a);
	float mb=mean(b);
	float S=0;
	for(int i=0;i<int(a.size());i++) S+= (a[i]-ma)*(b[i]-mb);
	S/=(sa*sb)*(a.size()-1);
	//printf("Going to return S=%f\n",S);
	return S;
	}
float STAT::corrSpearman(vector<float> &a ,vector<float> &b)
	{
	if(a.size() != b.size()) return -2.0;
	if(a.size()==0) return -2.0;
	vector<pair<float,float> > v;
	for(int i=0;i<int(a.size() );i++)
		v.push_back(pair<float,float>(a[i],b[i]));
	if( not is_sorted(v.begin(),v.end() ) ) sort(v.begin(),v.end());
	//togli i valori uguali
	vector<pair<float,float> > v2;//duplicate removal
	float current=v[0].first;
	float mean=0; int n=0;
	for(int i=0;i<int(v.size());i++)
		{
		if(v[i].first==current){mean+=v[i].second;n++;}
		else { v2.push_back(pair<float,float>(current,mean/n));
			mean=0;
			n=0;
			current=v[i].first;
			 }
		}
	if(n>0)v2.push_back(pair<float,float>(current,mean/n));	
	
	vector<float> z1; vector<float> z2; //copy v2 in z1;z2
	for(int i=0;i<int(v2.size());i++)
			{z1.push_back(v2[i].first);z2.push_back(v2[i].second);}
	return corrPearson(z1,z2);
	}
pair<float,float> STAT::regression(vector<float>&a,vector<float>&b)
	{
	if(a.size() != b.size()) {printf("ERROR\n");return pair<float,float>(-99,-99);}
	// y= mx+q
	float Sxx=0,Sxy=0;
	float ma=mean(a),mb=mean(b);
	for(int i=0;i<int(a.size());i++) {Sxx+=(a[i]-ma)*(a[i]-ma);Sxy=(a[i]-ma)*(b[i]-mb);}
	Sxx/=a.size();
	Sxy/=a.size();
	float m=Sxy/Sxx;
	float q=mb-(m*ma);
	pair<float,float> R(q,m);
	return R;
	}


float STAT::mean(vector<float> &a ,vector<float> &e_a)
	{
	if(a.size() != e_a.size() ){fprintf(stderr,"vector not have same dim in mean calculation\n");return -99;}

	float S=0;
	float InvN=0;
	for(int i=0;i<int(a.size());i++) S+=a[i]/(e_a[i]*e_a[i]);
	for(int i=0;i<int(e_a.size());i++) InvN+=1/(e_a[i]*e_a[i]);
	S*=InvN;
	return S;
	}

pair<float,float> STAT::regression(vector<float>&a,vector<float>&b,vector<float>&e_b,vector<float> &e2)
{
	vector<int> remove;
	e2.clear();
	unsigned int N=a.size();
	if( N != b.size()) {printf("ERROR\n");return pair<float,float>(-99,-99);}
	if( N != e_b.size()) {printf("ERROR\n");return pair<float,float>(-99,-99);}
	for(int i=0;i< N;i++) if (e_b[i]==0) {e_b[i]=1; fprintf(stderr,"Error: bin %d has no error.  will be removed.\n",i); remove.push_back(i);}
	for(int i=remove.size()-1;i>=0;i--) { 
					a.erase(a.begin()+remove[i]);
					b.erase(b.begin()+remove[i]);
					e_b.erase(e_b.begin()+remove[i]);
					}
	// y= mx+q
	float Sxx=0,Sxy=0,InvN=0,Sx=0,Sy=0;
	for(int i=0;i<int(a.size());i++) {
			InvN +=      1.          / (e_b[i]*e_b[i]);
			Sx   += (a[i])           / (e_b[i]*e_b[i]);
			Sy   += (b[i])           / (e_b[i]*e_b[i]);
			Sxx  += (a[i]*a[i])      / (e_b[i]*e_b[i]);
			Sxy  += (a[i]*b[i])      / (e_b[i]*e_b[i]);
			}
	float m=(InvN*Sxy - Sx*Sy )/(InvN*Sxx-Sx*Sx);
	float q=(Sy*Sxx-Sx*Sxy)/(InvN*Sxx-Sx*Sx);
	pair<float,float> R(q,m);

	//Error Computation [0]+[1]*x 0=q 1=m
	//e2[pair<int,int>(0,0)] = Sxx/(InvN*Sxx-Sx*Sx);
	//e2[pair<int,int>(1,1)] = InvN/(InvN*Sxx-Sx*Sx);
	//e2[pair<int,int>(1,0)]=e2[pair<int,int>(0,1)]=-Sx/(InvN*Sxx-Sx*Sx);
	e2.push_back( Sxx/(InvN*Sxx-Sx*Sx)   );
	e2.push_back( InvN/(InvN*Sxx-Sx*Sx)  );
	e2.push_back( -Sx/(InvN*Sxx-Sx*Sx)   );
	
	printf("Sx=%f Sy=%f Sxx=%f Sxy=%f InvN=%f\n",Sx,Sy,Sxx,Sxy,InvN);
	return R;

}

void STAT::drawFitError(TH1*h,pair<float,float> &R,vector<float> &e2,float sigma)
{
//|e f|^-1 = |a c|
//|f g|      |c b|
	float e=e2.at(0);
	float g=e2.at(1);
	float f=e2.at(2);
	
	float a=g/(-f*f+e*g);
	float c=-f/(-f*f+e*g);
	float b=e/(-f*f+e*g);
	
	float qm=R.first;
	float mm=R.second;

	for(int iBin=1;iBin<=h->GetNbinsX();iBin++){
		float x=h->GetBinCenter(iBin);
		float t=(c-a*x)/(b-c*x);
		float qt=sqrt(sigma*sigma/(a+2*c*t+b*t*t));
		float mt=t*qt;
	
		float ext1= (mm+mt)*x+(qm+qt);	
		float ext2= (mm-mt)*x+(qm-qt);	
		
		h->SetBinContent(iBin,mm*x+qm);
		h->SetBinError(iBin,fabs(mt*x+qt));
		}
	return ;
}

float STAT::ConfidenceInterval(std::vector<float> &v,std::pair<float,float>&r,float Q){
	if( not is_sorted(v.begin(),v.end() ) ) sort(v.begin(),v.end());
	int n=int(v.size()); 
	int m=ceil(n*Q);
	//Look for m consecutive bin such that the distance covered is minima
	vector<float> d;
	int min=0;
	for(int i=0;i<n-m;i++)
		{
		d.push_back(v[i+m]-v[i]);
		if(d[i]<d[min]) min=i;
		}
	r.first=v[min];
	r.second=v[min+m];
	return (r.second-r.first)/2.;
}


float STAT::ConfidenceIntervalAround(
		std::vector<float> &v,float value, 
		std::pair<float,float>&r, float Q)
{
	if( not is_sorted(v.begin(),v.end() ) )sort(v.begin(),v.end() );	
	//find position of value
	size_t vpos =0;
	for(size_t i=0;i<v.size();++i) if (v[i]<= value) vpos=i;

	int n=int(v.size()); 
	int m=ceil(n*Q/2.);
	
	size_t pos_low = (vpos > m)?vpos - m:0;
	size_t pos_high =(vpos+m <n)? vpos + m : n-1;
	
	r.first = v[pos_low];
	r.second = v[pos_high];
	return  (r.second-r.first)/2.;
}

float STAT::max(vector<float> &a)
{
	if (a.size() ==0 )
	{
		cout <<"[STAT]::[max]::[ERROR] a.size = 0 "<<endl;
		return -999;
	}
	float m=a[0];
	for(size_t i=0;i<a.size();++i) m = (m<a[i]) ? a[i] :m; 
	return m;
}

float STAT::min(vector<float> &a)
{
	if (a.size() ==0 )
	{
		cout <<"[STAT]::[min]::[ERROR] a.size = 0 "<<endl;
		return -999;
	}
	float m=a[0];
	for(size_t i=0;i<a.size();++i) m = (m>a[i]) ? a[i] :m; 
	return m;
}

// --------------- ROOT -----------
#include "TMatrixD.h"
#include "TVectorD.h"

float STAT::Chi2( vector<float> &a, vector<float> &b, vector<float> &ehigh, vector<float> &elow ,map<pair<int,int>,float> &corr)
{
	bool symm = elow.empty();
	bool useCorr = corr.empty();

	if (a.size() != b.size()) { cout <<"[STAT]::[Chi2]::[ERROR] a.size != b.size()"<<endl; return -999;}
	if (b.size() != ehigh.size()) { cout <<"[STAT]::[Chi2]::[ERROR] b.size != e.size()"<<endl; return -999.;}

	vector<float> e ;	
	if (symm) for( auto x : ehigh ) e.push_back(x);
	else
		{
		if (elow.size() != ehigh.size() ) {cout<<"[STAT]::[Chi2]::[ERROR] elow.size() != ehigh.size()"<<endl; return -999;}
		for( size_t i = 0 ;i< a.size() ;++i)
			{
			if (a[i]< b[i] ) e.push_back( elow[i] );
			else e.push_back(ehigh[i]);
			}
		} // else
	float chi2 =0;	
	if ( useCorr ) {
		for( size_t i = 0 ;i< a.size() ;++i)
			chi2 += (a[i] - b[i]) * (a[i] - b[i]) / (e[i]*e[i]);

		}
	else {
		TVectorD x;
		TMatrixD S;
		x.ResizeTo( a.size() );
		S.ResizeTo(a.size(),a.size());

		for( size_t i = 0 ;i< a.size() ;++i)
			x(i) = a[i]-b[i];

		for( size_t i = 0 ;i< a.size() ;++i)
		for( size_t j = 0 ;j< a.size() ;++j)
			{
			if (i==j and corr[pair<int,int>(i,j)] != 1)
				cout <<"[STAT]::[Chi2]::[WARNING] diagonal element are not one"<<endl;
			S(i,j) = corr[ pair<int,int>(i,j) ] * e[i]*e[j];
			}
		if (S.Determinant() == 0 ) cout<<"[STAT]::[Chi2]::[ERROR] Covariance matrix is singular"<<endl;
		S.Invert();

		chi2=x*(S*x);
		}
	return chi2;
}

void STAT::Fill(std::vector<float> &v, TH1*h)
{
	h->Reset("ACE");
	for(size_t i=0;i<v.size() ;++i)
		h->Fill(v[i]);
	return;
}

void STAT::Fill(std::vector<float> &a, std::vector<float> &b, TH2*h)
{
	h->Reset("ACE");
	if (a.size() != b.size())
		{
		cout<<"[STAT]::[Fill]::[ERROR] a and b vectors should have same size"<<endl;
		return;
		}
	for(size_t i=0;i<a.size() ;++i)
		h->Fill(a[i],b[i]);
	return;
}


#include "TMath.h"
TH1F *STAT::GetDensity(std::vector<float> &v, float R )
{
	if( not is_sorted(v.begin(),v.end() )) sort(v.begin(), v.end() );
	float low= v[0];
	float high = v[ v.size() -1 ] ;
	float diff = high -low;
	// 100 bins, in a range a bit (10%) extended
	TH1F *h = new TH1F("density","density",1000, low-.1*diff, high+.1*diff);

	for( int i=1 ;i < h->GetNbinsX() +1 ;++i) 
	{
		float mean = h->GetBinCenter(i);	
		float S=0;
		for(int j=0;j<v.size();++j)
			S+=TMath::Gaus( v[j], mean, R , kTRUE);
		S/=float(v.size());
		h->SetBinContent(i,S);
			
	}
	return h;
}


float STAT::Chi2(TGraphAsymmErrors *g, TH1 *h, TH2* corr)
{
	vector<float> x;
	vector<float> y;
	vector<float> eylow;
	vector<float> eyhigh;
	map<pair<int,int>,float> c;

	for(int i=0;i<g->GetN();i++)
	{
		y.push_back( g->GetY()[i] );	
		eylow.push_back( g->GetEYlow()[i] ) ;
		eyhigh.push_back( g->GetEYhigh()[i] ) ;
		x.push_back( h->GetBinContent(i+1) );
	}

	for(int i=0;i<g->GetN();i++)
	for(int j=0;j<g->GetN();j++)
	{
		c[pair<int,int>(i,j)] = corr->GetBinContent(i+1,j+1);
	}

	//float STAT::Chi2( vector<float> &a, vector<float> &b, vector<float> &ehigh, vector<float> &elow ,map<pair<int,int>,float> &corr)
	return Chi2(x,y,eyhigh,eylow, c);
}
