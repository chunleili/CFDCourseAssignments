// Central Differencing Implicit. Solving tridiagonal matrix by Thomas	algorithm
#include<iostream>
#include<cmath>
#define N 100
using namespace std;
void Thomas( double (&x)[N+1], double (&f)[N+1]);

int i;
double  maxT=5,t, C=0.7, dt=0.01, u[N+1]={0}, u1[N+1]={0}, a=-C/2, b=1.0, c=C/2;

int main()
{
	for(i=10; i<=20; i++)
	{
		u[i]=100;
	}

//B.C.
	u1[0]=u1[N+1]=0;
	Thomas(u1,u);
//print
	for(i=0; i<=100; i++)
	{
		//cout<<u[i]<<"\n";
	} 
return 0;
}

void Thomas( double (&x)[N+1], double (&f)[N+1])
{	
	double b_[N],f_[N], P[N], Q[N];

	P[1]=-c/b;
	Q[1]=f[1]-a*x[0];

	for( i=2;  i<=N;  i++)
	{
		P[i] =  		  - c / ( a * P[i-1] + b )  ;
		Q[i] =  ( f[i] - a * Q[i-1] ) / ( a * P[i-1] + b )  ;
	}
	for( t=0.0; t<=maxT;  t+=dt)
	{	
		
		x[N]=- c/b * x[N+1] + ( f[N] - c * x[N+1] ) / b ; 
		for(i=N-1; i>=1; i--)
		{
			x[i] = P[i] * x[i+1] + Q[i] ;
		}

		for(i=0; i<=N; i++)
		{
			f[i]=x[i];		
		}		
	
		if (fabs(t-1)<dt/2)
		{
			cout<<"\n\\\\\\try3\\\\\\\ntime="<<t+dt<< "\tu=\n";
			for(i=0; i<=100; i++)
			{
				cout<<u[i]<<"\n";
			} 
			cout<<"\nend output u=\n";
		}
	}
}
