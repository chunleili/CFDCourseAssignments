/*---------------------------------------------------------------------------*\
Description
	This File is to set the initial fields, including boundaries.
	It should be comipled secondly.

\*---------------------------------------------------------------------------*/

#include"main.H"
					  
void setInternal();
void setBoundaries();
void setLeft(const int N1=NVLong, const int N2=NVArc);
void setRight(const int N1=NVShort, const int N2=NVArc);
void setBottom(const int N1=NHShort, const int N2=NHArc);
void setTop(const int N1=NHLong, const int N2=NHArc);
static void printBoundaries(ofstream &fout, string flag);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
const double theta0=(M_PI-arcL/R)/2;
double dx, dy;
const double dThetaH=(arcL/R)/NHArc;///NHArc=NVArc=54
const double dThetaV=(arcL/R)/NVArc;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void setField()
{
	setBoundaries();

	ofstream foutBoundary("data/boundary.dat");
	printBoundaries(foutBoundary,"all");

	setInternal();

	ofstream foutInit("data/InitMesh.dat");
	printNodes(foutInit);

	cout<<"\nBoundary is written in \"boundary.dat\""<<endl;
	cout<<"Initial mesh is written in \"initMesh.dat\""<<endl;

}				  	

//						Set the Boundaries 								 	 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void setBoundaries()
{
	setLeft();
	setRight();
	setBottom();
	setTop();
}

void setLeft(const int N1, const int N2)
{
		const int i=0, M=My;
		const double Ox=-Ob, Oy=totalL-Oa, L1=longL, L2=pitL;
		double d1=L1/N1, d2=arcL/N2, d3= (totalL-L1-L2)/(M-N1-N2);
		int j;
		for(j=0;j<N1;j++)
		{
			node[i][j].x=0;
			node[i][j].y=d1*j;
		}
		int n;
		for(j=N1;j<N1+N2;j++)
		{
			n=j-N1;
			dx=R*sin(theta0+n*dThetaV);
			dy=R*cos(theta0+n*dThetaV);
			dy*=-1;///it depends on where is Circle center relative to
					///the first point of arc. Simply put, if first point is
					///on the left, down side of the Center , dx*=-1,dy*=-1;
			node[i][j].x=Ox+dx;
			node[i][j].y=Oy+dy;
		}

		for(j=N1+N2;j<=M;j++)
		{
			n=j-N1-N2;
			node[i][j].x=0;
			node[i][j].y=L1+L2+d3*n;
		}
}

void setRight(const int N1, const int N2)
{
		const int i=Mx, M=My;
		const double Ox=totalL+Ob, Oy=Oa, L1=shortL, L2=pitL;
		double d1=L1/N1,d3=(totalL-L1-L2)/(M-N1-N2);
		int j;
		for(j=0;j<N1;j++)
		{
			node[i][j].x=totalL;
			node[i][j].y=d1*j;
		}
		int n;
		for(j=N1;j<N1+N2;j++)
		{
			n=j-N1;
			dx=R*sin(theta0+n*dThetaV);
			dy=R*cos(theta0+n*dThetaV);
			dx*=-1;
			dy*=-1;
			node[i][j].x=Ox+dx;
			node[i][j].y=Oy+dy;
		}
		for(j=N1+N2;j<=M;j++)
		{
			n=j-N1-N2;
			node[i][j].x=totalL;
			node[i][j].y=L1+L2+d3*n;
		}
}


void setBottom(const int N1, const int N2)
{
		const int j=0, M=Mx;
		const double Ox=Oa, Oy=-Ob, L1=shortL, L2=pitL;
		double d1=L1/N1, d3= (totalL-L1-L2)/(M-N1-N2);
		int i;
		for(i=0;i<N1;i++)
		{
			node[i][j].x=dKsi*i;
			node[i][j].y=0;
		}
		int n;
		for(i=N1;i<N1+N2;i++)
		{
			n=i-N1;
			dx=R*cos(theta0+n*dThetaH);///sin/Cos is opposite to Left&Right
			dy=R*sin(theta0+n*dThetaH);
			dx*=-1;
			node[i][j].x=Ox+dx;
			node[i][j].y=Oy+dy;
		}
		for(i=N1+N2;i<=M;i++)
		{
			n=i-N1-N2;
			node[i][j].x=L1+L2+d3*n;
			node[i][j].y=0;
		}
}

void setTop(const int N1, const int N2)
{
		const int j=My, M=Mx;
		const double Ox=longL+pitL/2, Oy=totalL+Ob, L1=longL, L2=pitL;
		double d1=L1/N1, d3= (totalL-L1-L2)/(M-N1-N2);
		int i;
		for(i=0;i<N1;i++)
		{
			node[i][j].x=d1*i;
			node[i][j].y=totalL;
		}
		int n;
		for(i=N1;i<N1+N2;i++)
		{
			n=i-N1;
			dx=R*cos(theta0+n*dThetaH);///sin/Cos is opposite to Left&Right
			dy=R*sin(theta0+n*dThetaH);
			dx*=-1;
			dy*=-1;
			node[i][j].x=Ox+dx;
			node[i][j].y=Oy+dy;
		}
		for(i=N1+N2;i<=M;i++)
		{
			n=i-N1-N2;
			node[i][j].x=L1+L2+d3*n;
			node[i][j].y=totalL;
		}
}

//						Set the internal fields							 	 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void setInternal()
{
		int i,j;
		for(i=1;i<Mx;i++)
			for(j=1;j<My;j++)
			{
				node[i][j].x=pitDepth+dKsi*(totalL-2*pitDepth)/totalCurveL*i;
				node[i][j].y=pitDepth+dEta*(totalL-2*pitDepth)/totalCurveL*j;
			}
}

//						Print the Boundaries 							 	 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
static void printBoundaries(ofstream &fout, string flag)
{	
	printHead(fout);

//	fout<<"flag= "<<flag<<endl;
	if(!(flag.compare("left")*flag.compare("all")))
	{
	//	fout<<"begin print left"<<endl;
		const int i=0;
		for(int j=0;j<=My;j++)
			fout<<node[i][j].x<<'\t'<<node[i][j].y<<endl;
	//	fout<<"\nend print left\n\n"<<endl;
	}

	if(!(flag.compare("right")*flag.compare("all")))
	{
	//	fout<<"begin print right"<<endl;
		const int i=Mx;
		for(int j=0;j<=My;j++)
			fout<<node[i][j].x<<'\t'<<node[i][j].y<<endl;
	//	fout<<"\nend print right\n\n"<<endl;
	}

	if(!(flag.compare("right")*flag.compare("all")))
	{
	//	fout<<"begin print bottom"<<endl;
		const int j=0;
		for(int i=0;i<=Mx;i++)
		{
			fout<<node[i][j].x<<'\t'<<node[i][j].y<<endl;
		}
	//	fout<<"\nend print bottom\n\n"<<endl;
	}

	if(!(flag.compare("right")*flag.compare("all")))
	{
	//	fout<<"begin print top"<<endl;
		const int j=My;
		for(int i=0;i<=Mx;i++)
		{
			fout<<node[i][j].x<<'\t'<<node[i][j].y<<endl;
		}
	//	fout<<"\nend print top\n\n"<<endl;
	}
	
}