#ifndef main_H
#define main_H

#include<iostream>
#include<fstream>
#include<cmath>
#include<string>

using namespace std;


//						declear the constants								 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	const int maxSpace= 100;
	const double physicalSpace = 2;
	const double dx = physicalSpace / maxSpace;
	
	const int maxTime = 100;
	const double PhysicalTime = 0.4;
	const double dt = PhysicalTime / maxTime;

	const double GAMMA = 1.4;

//						define the func										 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void UToF(double U[][maxSpace+2], double F[][maxSpace+2]);
void MacCormackSolver( double U[3][maxSpace+2], double F[3][maxSpace+2]);
void init(double U[][maxSpace+2], double F[][maxSpace+2]);


#endif
