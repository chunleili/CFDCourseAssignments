/*---------------------------------------------------------------------------*\
Description
	This File is to define functions to print out data.
	It should be compiled firstly.

\*---------------------------------------------------------------------------*/

#include"main.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void printHead(ofstream &fout)
{
	fout
	<<"Title=\"Mesh\""<<endl
	<<"Variables=\"x\",\"y\""<<endl
	<<"Zone i="<<Mx+1<<", j="<<My+1<<", f=point"<<endl;
}
void printNodes(ofstream &fout)
{
	printHead(fout);

	for(int j=0;j<=My;j++)	
	{
		for(int i=0;i<=Mx;i++)
		{
			fout<<node[i][j].x<<'\t'<<node[i][j].y<<endl;
		//	fout<<"\t\t["<<i<<"]"<<"["<<j<<"]"<<endl;
		}
	}
}
