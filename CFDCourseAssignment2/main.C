/*---------------------------------------------------------------------------*\
Description
	This File is the main file, which aims to control the loop and gives final
	results.

\*---------------------------------------------------------------------------*/

#include"main.H"
#include<cstdlib>
#include<ctime>
const int LOOP_LIMIT=1e4;
const double RESIDUAL_LIMIT=1e-5;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
///definition of the global arries here
struct Cnode node[Mx+1][My+1];

int main()
{	
	cout<<"\nSetting the initial field...\n";
	setField();

	cout<<"\nInitial field set over\n";

	cout<<"\nResidual limit= "<<RESIDUAL_LIMIT<<endl;
	cout<<"Loop times limit = "<<LOOP_LIMIT<<endl;
	cout<<"Number of horizentoal nodes(Mx+1)= "<<Mx+1<<endl;
	cout<<"Number of vertical nodes(My+1)= "<<My+1<<endl;
	cout<<"\nBegin iterate ? (y/n)\n";
//	char c;
//	cin>>c;
//	if(c=='n')
//	{
//		exit(0);
//	}

	ofstream foutResidual("data/residual.dat");
	ofstream foutFinalMesh("data/finalMesh.dat");
	int loopTimes=0;
	double  residual;

	while(loopTimes<LOOP_LIMIT)
	{
		residual=iterate();
		foutResidual<<residual<<endl;
		
		loopTimes++;

		cout<<"loopTimes= "<<loopTimes
			<<"\tResidual= "<<residual<<endl;

		if(residual<RESIDUAL_LIMIT)
		{	
			cout<<"\nConverged!"<<endl
			<<"\nMx= "<<Mx<<"\tMy= "<<My<<endl
			<<"Final residual ="<<residual<<endl
			<<"Loop times = "<<loopTimes<<endl
			<<"Final field has been written into \"finalMesh.dat\" \n"
			<<"History of residual has been written into \"residual.plt\"";
			printNodes(foutFinalMesh);
			break;
		}
	}

	if(loopTimes>=LOOP_LIMIT)
	{
		cout<<"\nFail to converge! in "<<loopTimes<<" times\n"
		<<"Final residual ="<<residual<<"\n\n\n";
	}
	cout<<"\nWall time = "<<(double)clock()/CLOCKS_PER_SEC<<" s"<<endl;
	cout<<"CPU time= "<<(double)clock()<<" s"<<endl;
return 0; 
}
