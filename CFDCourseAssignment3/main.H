#ifndef main_H
#define main_H
#include<cmath>
#include<cstdlib>
#include<iostream>
using namespace std;

//						define the MACRO func								 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define is0(den) if(fabs(den)<1e-5) cout<<#den<<" is 0! cannot be denominator!"<<endl<<endl
#define isNegative(base) if((base)<0) cout<<#base<<" is NEGATIVE! cannot sqrt()!"<<endl<<endl
#define whereNegative(base, step, I) \
if((base)<0) \
cout<<"sqrt erro!!"<<"\nbase= "<<#base<<"\nstep= "<<#step<<"\nI= "<<#I<<endl<<endl

//						declear the constants								 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	const int maxSpace= 100;
	const double physicalSpace = 2;
	const double dx = physicalSpace / maxSpace;

	const int maxTime = 100000;
	const double stopTime = 0.1;
	const double CFL=0.4;

	const double GAMMA = 1.4;

//						define the types									 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //	
typedef double Field[maxSpace+1][3];
typedef double ScalarField[maxSpace+1];
typedef double Vector[3];

//						declear the inline funcs							 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
inline double WToP(Vector W)
{
    return (GAMMA-1) * ( W[2] - 0.5* W[1]*W[1]/W[0] );
}

//						define the func										 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void zeroGradBC(Field W);
void WToF(Vector W, Vector F);

#endif
