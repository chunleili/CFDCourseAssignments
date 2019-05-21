#ifndef main_H
#define main_H

#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
using namespace std;

void setField();
double iterate();

void printHead(ofstream &);
void printNodes(ofstream &);

//						define the inline func								 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline double max(const double a, const double b)
{
	return (a>=b)?a:b;
}

inline double sq(const double a)
{
	return pow(a, 2);
}

//						declear the constants								 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
	const int Mx=300, My=300;

	const double pitL=1.5, pitDepth=0.4, totalL=3.0, longL=1.25;
	const double shortL=totalL-longL-pitL;
	const double R=( sq(pitL/2)+ sq(pitDepth) ) / (2*pitDepth);///R=0.903125
	const double arcL=2*R*asin( pitL/(2*R) ), sqR=sq(R);///arcL=1.76997
	const double Oa=shortL+pitL/2, Ob=R-pitDepth;/// Oa=1, Ob=0.503125
	///Oa & Ob stands for the x & y distance(Oa && Ob>0)  of  
	///bottom Circile Center to the left bottom corner of the plane

	const double totalCurveL=arcL+longL+shortL;
	const double dKsi=totalCurveL/Mx, dEta=totalCurveL/My;///dKsi=dEta=0.0327

	const int NHShort=round(shortL/dKsi);///8
	const int NHLong =round(longL/dKsi);///38
	const int NHArc  =Mx-NHShort-NHLong;///NHArc=54

	const int NVShort=round(shortL/dEta);
	const int NVLong =round(longL/dEta);
	const int NVArc  =My-NVShort-NVLong;
	
//						declear the globals								 	 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //	
struct Cnode
{
	double x,y;
};

extern struct Cnode node[Mx+1][My+1];

#endif
