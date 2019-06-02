#ifndef main_H
#define main_H

#include<cmath>
#include<iostream>
using namespace std;

#define CONV scalarJSTConv
//  #define CONV MacCormackConv
//  #define CONV RoeConv
//#define TIME_DIS EulerFTime
#define TIME_DIS RungeKuttaTime

//						declear the constants								 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	const int maxSpace= 1000;
	const double physicalSpace = 2;
	const double dx = physicalSpace / maxSpace;

	const int maxTime = 100000;
	const double stopTime = 0.1;
	const double CFL=1.5;

	const double GAMMA = 1.4;

//						declear the inline funcs							 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
inline double p(const double W[])
{
	double u = W[1]/W[0];
    return (GAMMA-1) * W[0]* ( W[2]/W[0] - 0.5* u * u );
}

inline  double max(const double a, const double b)
{
    return a>b?a:b;
}

//计算声速c
inline  double c(const double W[3])
{
    return sqrt( GAMMA * p(W)/W[0] );
}

inline  double lambda(const double W[3])
{
    return fabs(W[1]/W[0]) + c(W);
}

inline double WToRho(const double vec[3])
{
	return vec[0];
}

inline double WToVelocity(const double vec[3])
{
	return vec[1]/vec[0];
}

//						define the func										 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void CONV(const double W[][3], const double dtLocal, double R[][3], const int I);
void TIME_DIS(double W[][3], double dtLocal, double R[][3], const int I);
void zeroGradBC(double W[][3]);
void init(double W[][3]);
double localTime(const double W[][3],const int I);
void WToF(double const W[3], double F[3]);
void print( const double W[][3]);
void printW(const double W[][3]);
double minOverScalarField(double scalarField[]);

#endif
