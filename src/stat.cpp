#include <vector>
#include <algorithm>

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
	size_t pos_high =(vpos+m <n)? vpos + m : n;
	
	r.first = v[pos_low];
	r.second = v[pos_high];
	return  (r.second-r.first)/2.;
}
