#include<iostream>
#define I 101

using namespace std;

double x=0.0; 	
double const dx=0.01, dt=0.01, alpha=0.7, C=alpha*dt/dx ;
int    i=0,n=0;
int    const maxI=1.0/dx;
double const maxN=0.7/alpha/dt;
double	u[I]={0.0}, u1[I]={0.0}, u0[I]={0.0};

void update();
void first2Steps(int flag);
void Explicit();
void Implicit();
void FogJump();
void CrankNicolson();
void MacCormark();
void LaxWendroff();
int main()
{
	//u1: next time step u; u0 former timestep
// I.C. @n=0, x=[0.1,0.2] : u=100; @anywhere else, u=0;
	for(i=10, x=0.1; x<=0.2; i++, x+=dx)
	{
		u[i]=100;
	}
	cout<<"maxI="<<maxI<<endl;
	cout<<"maxN="<<maxN<<endl;
//B.C. @x=0 : u1=0
	u1[0]=10;	

//	u[459]=0;
//Explicit 
	//Explicit();
	//LaxWendroff();
	//FogJump();
	MacCormark();
//Implicit
	//Implicit();
	//CrankNicolson();
	
//print
	for(i=0; i<=maxI; i++)
		cout<<u[i]<<"\n";
return 0;
}	

void update()
{
	for(i=0; i<=maxI; i++)
	{
		u[i]=u1[i];
	}
}

inline void first2Steps(int flag)
{
	for(i=1;i<=flag;i++)
		u1[i]=u[i] - C * ( u[i+1] - u[i]);
}


	
// Euler Explicit & Central Diff Explicit
void Explicit()
{
	for( n=0; n<=maxN;  n++)
	{
		//change the flag value to accomodate different schemes
		int flag=1;	//flag==0 means no need to solve the first 2 space step with Euler Explicit
		if(flag>=1)	//flag==1 means need to solve the first step, flag==2 means need to solve the first and second step
			first2Steps(flag);
		for(i=flag; i<=maxI-1; i++) 
		{
			//u1[i]=u[i] - C * ( u[i+1] - u[i]);//Euler Forewards Explicit //instable!!
			u1[i]=u[i] - C *   ( u[i]   - u[i-1]);//Euler Backwards Explicit //flag=1
			//u1[i]=u[i] - C/2 * ( u[i+1] - u[i-1]);//Central Differencing Explicit//falg=1// instable!! 
			//u1[i] = u[i] - C/2 * ( 3*u[i+1] -4*u[i] + u[i-1] );//Triple Point Euler Central Explicit //flag=1//instable!! 
			//u1[i] = u[i] - C/2 * ( 3*u[i]   -4*u[i-1] + u[i-2] );//Triple Point Euler Backwards Explicit // flag=2//instable!!
		}
		u1[maxI]=u1[maxI-1]; //right BC: zero gradient
		update();
	}
}


//FogJump Explicit 
void FogJump()
{
	n=1; 
	for(i=1;i<=maxI;i++)			//must solved by Euler Explict before Fogjump, to get value when n=1 (t=dt).
		u1[i]=u[i] - C * ( u[i+1] - u[i]);
	for(n=2; n<=maxN;  n++)
	{
		int flag=1;
		first2Steps(flag);
		for(i=flag; i<=maxI-1; i++)
		{
			u1[i] = u0[i] - C*( u[i+1] - u[i-1] );	//FogJump instable!!
		}
		u1[maxI]=u1[maxI-1];
		update();
	}
}

//Lax-Wendroff Explicit
void LaxWendroff()
{
	int flag=1;
	first2Steps(flag);
	for( n=0; n<=maxN;  n++)
	{
		for(i=flag; i<=maxI-1; i++)
		{
			u1[i] = u[i] - C/2 * ( u[i+1] - u[i-1] ) + C*C/2 * (u[i+1] - 2*u[i] + u[i-1]);	//LaxWendroff
		}
		u1[maxI]=u1[maxI-1];
		update();
	}
}


//Euler Implicit
void Implicit()
{
	int flag=1;
	if(flag>=1) first2Steps(flag);
	for(n=0; n<=maxN;  n++)
	{		
		for(i=flag; i<=maxI-1; i++)
		{
			//u1[i]=(  u[i] - C * u1[i+1]  ) / (1-C); //Euler Forewards Implicit //instable!!//flag=0
			u1[i]=(  u[i] + C * u1[i-1]  ) / (1+C); //Euler Backwards Implicit	//flag=1	
		}
		u1[maxI]=u1[maxI-1];
		update();
	}
}


//Crank-Nicolson
void CrankNicolson()
{
	int flag=1;
	
	for(n=0;n<=maxN;n++)
	{	
		first2Steps(flag);
		for(i=flag;i<=maxI-1;i++)
		{
			u1[i] = ( C/2 * u1[i-1] + ( 1 - C/2 ) * u[i] + C/2 * u[i-1] )  / ( 1 + C/2 );
		}
		u1[maxI]=u1[maxI-1];
		update();
	}
}

//MacCormark
void MacCormark()
{
	int flag=1;
	double dudt[I], dudt1_, dudtav, u1_[I];
	for(n=0;n<=maxN;n++)
	{	
		first2Steps(flag);
		//predict
		for(i=flag;i<=maxI-1;i++)
		{

			dudt[i] = -alpha  * (u[i+1] - u[i]) /dx;
			u1_[i]  = u[i]+ dudt[i] * dt;
		}
		//revise	
		for(i=flag;i<=maxI-1;i++)
		{
			dudt1_ = -alpha * (u1_[i] - u1_[i-1] ) /dx;
			dudtav  = 0.5* (dudt1_ + dudt[i] );
			u1[i]   = u[i] + dudtav * dt;
		}
		u1[maxI]=u1[maxI-1];
		update();
	}
}
	
