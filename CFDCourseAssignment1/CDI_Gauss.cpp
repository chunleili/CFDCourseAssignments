// Central Differencing Implicit. Solving tridiagonal matrix by Thomas	algorithm
#include<iostream>
#include<cmath>
#define N 100
using namespace std;
void Thomas(double a, double b, double c, double (&x)[N+1], double (&f)[N+1]);

int i;
double  maxT=1.0,t, C=0.7, dt=0.01, u[N+1]={0}, u1[N+1]={0}, a=-C/2, b=1.0, c=C/2;

int main()
{
	for(i=10; i<=20; i++)
	{
		u[i]=100;
	}
/*
	cout<<"\nafter initial u=\n";
	for(i=0; i<=100; i++)
	{
		cout<<u[i]<<"\n";
	} 
	cout<<"\nend output initial u=\n";
*/

//B.C.
	u1[0]=u1[N+1]=0;
	Thomas(a,b,c,u1,u);
//print
	for(i=0; i<=100; i++)
	{
		//cout<<u[i]<<"\n";
	} 
return 0;
}

void Thomas(double a, double b, double c, double (&x)[N+1], double (&f)[N+1])
{	
	double b_[N],f_[N];

	for( t=0.0; t<=maxT;  t+=dt)
	{	
		b_[1]=b;
		f_[1]=f[1]-a*x[0];
		f_[N]=f[N]-c*x[N+1];	
		
		for( i=2;  i<=N;  i++)
		{
			b_[i] =  b    - c 	*  a/  b_[i-1]  ;
			f_[i] =  f[i] - f_[i-1] *  a/  b_[i-1]  ;
		}
		
		x[N]=f_[N] / b_[N]; 
		for(i=N-1; i>=1; i--)
		{
			x[i] = ( f_[i]  -  c  *  x[i+1] )  /  b_[i];
		}
		//if(t>0.1 && t<0.11) cout<<"??t=??"<<t<<endl;
		//cout<<"///////////////t="<<t<<"////////////////////";
		f[1]   = f_[1]  + a * x[0];		
		f[N]   = f_[N]  + c * x[N+1];

		for(i=0; i<=N; i++)
		{
			f[i]=x[i];
			//cout<<"f["<<i<<"]="<<f[i]<<endl;			
		}		
	
	if (fabs(t-0.5)<dt/2){
	cout<<"\ntime="<<t+dt<< "\tu=\n";
	for(i=0; i<=100; i++)
	{
		cout<<u[i]<<"\n";
	} 
	cout<<"\nend output u=\n";
	}
	}
}
