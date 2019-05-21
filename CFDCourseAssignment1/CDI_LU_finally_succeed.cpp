/*---------------------------------------------------------------------------*\
Description
    Implementing the Central Differencing Implicit Discretization for transport equation, i.e. :
	ddt(u) + a * ddx(u) == 0;
    after discretizing, it becomes :
	-C/2 * u1[i-1] + u1[i] + C/2 * u1[i+1] == u[i];
    where C = a*dt/dx is the Courant number;
    u1 denotes velocity in next time step, and u denotes velocity in current time step.
    
    Essentially, it is a tridiagonal matrix equation, i.e. :
	A x == f
    To solve this tridiagonal equation, we use Thomas Algorithm,
    which is the main object of this cpp file.

\*---------------------------------------------------------------------------*/

#include<iostream>
#define N 100

using namespace std;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
int main()
{
//define variables and matrix. maxT denotes maximum running time, a.k.a stop time
	const double dx=0.01, dt=0.01, alpha=0.7,  C=alpha*dt/dx,  maxT=0.7/alpha;

//x[] denotes solution vector, in this case standing for u[] in next time step, which is to solve.
//f[] denotes constant vector, in this case standing for u[] in current time step, which is already known.
//beta[] is the auxilary vector, storing the corresponding coefficients after LU.
	double	t=0, x[N+2]={0}, f[N+2]={0}, beta[N],  y[N];
	int i=0;

//I.C.	
	for(i=10; i<=20; i++)
	{
		f[i]=100;
	}
//B.C.	
	x[0]=10;
	
	double a=-C/2, b=1.0, c=C/2;

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//Thomas Algorithm Begin

//LU, once is enough 

	beta[1]=c/b;
	for(i=2;i<=N-1;i++)
	{
		beta[i]=c / ( b - a * beta[i-1] );
	}

//In every time step, solve Ly=f, Ux=y
	for(t=0; t<=maxT; t+=dt)
	{	
	//B.C. in the right side: zero gradient
		x[N+1]=f[N];
	//treatment of maxtix, change it to tridiagonal
		f[1] = f[1] - a * x[0];
		f[N] = f[N] - c * x[N+1];

	//chasing, to solve Ly=f
		y[1] = f[1] / b;
		for(i=2;i<=N;i++)
		{
			y[i]=( f[i] - a * y[i-1] ) / (b - a * beta[i-1] );
		}

	//passing, to solve Ux=y
		x[N]=y[N] ;
		for(i=N-1;i>=1;i--)
		{
			x[i] = y[i] - beta[i] * x[i+1];
		}

	//before moving to next time step, swap x and f
		for(i=0;i<=N+1;i++)
		{
			f[i]=x[i];
		}
	// change back f[]
		f[1] = f[1] - a * x[0];
		f[N] = f[N] - c * x[N+1];	
	}

//print
	for(i=0; i<=N; i++)
	{
		cout<<x[i]<<"\n";
	}
return 0;
}
// ************************************************************************* //	
			
		 

