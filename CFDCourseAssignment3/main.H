#ifndef main_H
#define main_H
#include<cmath>
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

//						define the types								 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //	
typedef double Field[maxSpace+1][3];
typedef double ScalarField[maxSpace+1];
typedef double Vector[3];

//						declear the inline funcs							 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
inline double calPressure(const double W[])
{
	double u = W[1]/W[0];
    return (GAMMA-1) * W[0] *( W[2]/W[0] - 0.5* u * u );
}

inline  double soundVelocity(const double W[3])
{
    return sqrt( GAMMA * calPressure(W)/W[0] );
}

inline  double lambda(const double W[3])
{
    return fabs(W[1]/W[0]) + soundVelocity(W);
}


//						define the func										 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void zeroGradBC(Field W);
void WToF(Vector W, Vector F);

#endif
