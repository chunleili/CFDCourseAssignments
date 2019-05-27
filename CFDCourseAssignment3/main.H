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
	const double PhysicalTime = 0.4;
	const double dt = PhysicalTime / maxTime;

	const double GAMMA = 1.4;

//						declear the inline funcs							 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
inline double p(const double U[])
{
    return (GAMMA-1) * ( U[2] - 0.5*U[1]*U[1] / U[0] );
}



//						define the func										 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void UToF(double const U[3], double F[3]);
void MacCormackSolver( double U[][3]);
void scalarJSTSolver(double U[][3]);
void init(double U[][3]);
void print( const double U[][3]);
void printU(const double U[][3]);

#endif
