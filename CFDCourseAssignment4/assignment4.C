/***************************include & namespace  **************************/
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<utility>
#include<cmath>
#include<cstdlib>
using namespace std;
/***************************define the consts ********************************/

const double stop=1.0;
const int stopStep=1000;
const int maxI=400, maxJ=100;

/***************************define the type ********************************/
typedef struct xy
{
    public:
    double x=0;
    double y=0;
}xy;

typedef xy MeshPoint[maxI+1][maxJ+1];
typedef double Field[maxI+1][maxJ+1][5];

/***************************declare the funcs  **************************/
void genMesh();
void init1();
void init2();
double localTimeStepping();
void solve();
void print(Field aField, string filename);

/********************************main()**************************************/
int main()
{  
    genMesh();              //生成网格
    
    Field Q, F;
    
    cout<<"please input case Number: 1 for inlet 1.5Ma; 2 for inlet 1.8Ma"<<endl;
    int caseNo;
    cin>>caseNo;
    if(caseNo==1)       init1(); //设置初场,init1表示入口1.5Ma, init2()表示入口1.8Ma
    else if(caseNo==2)  init2();
    else {
        cout<<"wrong number! exit!"<<endl;
        exit(0);
    }

    int tStep=0;
    for (double t = 0, dt=0.001; t < stop && tStep<=stopStep; t += dt, tStep++)
    {
        dt = localTimeStepping(); //当地时间步法，注意取全域最大当地时间步
        solve();                  //求解
    }

    //print(Q,"Q.txt");                //打印 
    return 0;
}

/********************************funcs define**************************************/
void genMesh()
{
    const int block1=100, block2=100, block3=maxI-block1-block2;
    const double dx=1.0/500;

    MeshPoint mesh;

    double dy=0.3/maxJ;
    for(int i=0; i<block1; i++)
    {
        for(int j=0;j<=maxJ; j++)
        {
            mesh[i][j].x=i*dx;
            mesh[i][j].y=j*dy;
        }
    }
    
    for(int i=block1; i<block1+block2; i++)
    {
        double h=0.25*i*dx-0.05;
        dy=(0.3-h)/maxJ;
        for(int j=0;j<=maxJ; j++)
        {
            mesh[i][j].x=i*dx;
            mesh[i][j].y=h+j*dy;
        }
    }

    dy=(0.3-0.05)/maxJ;
    for(int i=block1+block2; i<=maxI; i++)
    {
        for(int j=0;j<=maxJ; j++)
        {
            mesh[i][j].x=i*dx;
            mesh[i][j].y=0.05+j*dy;
        }
    }

    fstream fout("mesh.dat");
    fout
	<<"Title=\"Mesh\""<<endl
	<<"Variables=\"x\",\"y\""<<endl
	<<"Zone i="<<maxI+1<<", j="<<maxJ+1<<", f=point"<<endl;
    for(int j=0; j<=maxJ; j++)
    {
        for(int i=0; i<=maxI; i++)
        {
            fout<< mesh[i][j].x<<" "<<mesh[i][j].y<<endl;
        }
    }

}

void init1()
{

}

void init2()
{

}

double localTimeStepping()
{
    double globalTime=0;
    return globalTime;
}

void solve()
{

}

void print(Field aField, string filename)
{
    fstream fout(filename);
    for(int j=0; j<=maxJ; j++)
    {
        for(int i=0; i<=maxI; i++)
        {
            for(int k: aField[i][j])
                fout<< aField[i][j][k]<<" ";
            fout<<endl;
        }
    }
}

/*******************other utility funcs*****************/
