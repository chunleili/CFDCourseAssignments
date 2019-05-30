#ifndef main_H
#define main_H

#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdio>
#include<iomanip>

using namespace std;


//						declear the constants								 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	const int maxSpace= 1000;
	const double physicalSpace = 2;
	const double dx = physicalSpace / maxSpace;

	const int maxTime = 100000;
	const double stopTime = 0.1;

	const double GAMMA = 1.4;

//						declear the inline funcs							 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
inline double p(const double W[])
{
	double u = W[1]/W[0];
    return (GAMMA-1) * W[0]* ( W[2]/W[0] - 0.5* u * u );
}

inline  double max(double a, double b)
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


//						define the func										 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void WToF(double const W[3], double F[3]);
void MacCormackConv( double W[][3], const double dt);
void scalarJSTConv(double W[][3], const double dt, double R[][3]);
void RoeConv(double W[][3], const double dt, double R[][3]);
void init(double W[][3]);
void print( const double W[][3]);
void printW(const double W[][3]);
double localTime(const double W[][3]);
void EulerFTime(double W[][3], const double dt, const double R[][3]);
void zeroGradBC(double W[][3]);

#endif
