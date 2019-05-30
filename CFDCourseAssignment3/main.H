#ifndef main_H
#define main_H

#include<iostream>
#include<fstream>
#include<cmath>
#include<cstdio>

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
inline double p(const double U[])
{
    return (GAMMA-1) * ( U[2] - 0.5*U[1]*U[1] / U[0] );
}

inline  double max(double a, double b)
{
    return a>b?a:b;
}
//计算声速c
inline  double c(const double U[3])
{
    return sqrt( GAMMA * p(U)/U[0] );
}

inline  double lambda(const double U[3])
{
    return fabs(U[1]/U[0]) + c(U);
}


//						define the func										 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void UToF(double const U[3], double F[3]);
void MacCormackConv( double U[][3], const double dt);
void scalarJSTConv(double U[][3], const double dt, double R[][3]);
void RoeConv(double U[][3], const double dt, double R[][3]);
void init(double U[][3]);
void print( const double U[][3]);
void printU(const double U[][3]);
double localTime(const double U[][3]);
void EulerFTime(double U[][3], const double dt, const double R[][3]);

#endif
