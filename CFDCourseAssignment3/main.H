#ifndef main_H
#define main_H
#include<cmath>
#include<iostream>
using namespace std;
//						declear the constants								 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	const int maxSpace= 1000;
	const double physicalSpace = 2;
	const double dx = physicalSpace / maxSpace;

	const int maxTime = 100000;
	const double stopTime = 0.1;
	const double CFL=0.4;

	const double GAMMA = 1.4;

//						declear the inline funcs							 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
inline double calPressure(const double W[])
{
	double u = W[1]/W[0];
    return (GAMMA-1) * W[0] *( W[2]/W[0] - 0.5* u * u );
}

inline  double max(const double a, const double b)
{
    return a>b?a:b;
}

inline  double soundVelocity(const double W[3])
{
    return sqrt( GAMMA * calPressure(W)/W[0] );
}

inline  double lambda(const double W[3])
{
    return fabs(W[1]/W[0]) + soundVelocity(W);
}

inline double WToRho(const double vector[3])
{
	return vector[0];
}

inline double WToVelocity(const double vector[3])
{
	return vector[1]/vector[0];
}

//						define the func										 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void scalarJSTConv(const double W[][3],  double R[][3]);
void RungeKuttaTime(double W[][3], const double dt, double R[][3]);
void EulerFTime(double W[][3], const double dt, double R[][3]);
void zeroGradBC(double W[][3]);
void init(double W[][3]);
double localTime(const double W[][3],const int I);
void WToF(double const W[3], double F[3]);
void print( const double W[][3]);
void printW(const double W[][3]);
void printF(const double W[][3]);
double calDtGlobal(double W[][3]);
#endif
