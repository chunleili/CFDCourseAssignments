/***************************include & namespace  **************************/
#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<ctime>
using namespace std;
/***************************define the consts ********************************/

const double STOP_TIME=1.0;
const int STOP_STEP=1000;
const int maxI=400, maxJ=100;
const int RESIDUAL_LIMIT=1e-3;

/***************************define the type ********************************/
typedef struct XY
{
    public:
    double x=0;
    double y=0;
}XY;

typedef XY MeshPoint[maxI+1][maxJ+1];        //用于存储网格点坐标
typedef double Field[maxI+1][maxJ+1][4];     //向量场,用于定义Q对象
typedef double Tensor[maxI+1][maxJ+1][4][2]; //张量场,用于定义F对象

/***************************declare the funcs  **************************/
void genMesh();
void init1();
void init2();
double localTimeStepping();
void solve();
void print(Field aField, fstream &);

//utility functions
inline void runInfo(double t, double dt, int step, double residual);
inline double calMa(double T, double magU);
inline double calSoundSpeed(double T);

/********************************main()**************************************/
int main()
{  
    genMesh();              //生成网格
    
    Field Q;
    Tensor F;

    cout<<"please input case Number: 1 for inlet 1.8Ma; 2 for 1.5Ma"<<endl;
    int caseNo; cin>>caseNo;
    switch(caseNo)
    {
        case (1): init1(); break; //设置初场,init1表示入口1.5Ma, init2()表示入口1.8Ma
        case (2): init2(); break;
        default: cout<<"wrong number! exit!"<<endl; return 1;
    }
    
    int step=0; double residual=0.0; double t=0.0;
    for (double dt=0.001; t <= STOP_TIME && step<=STOP_STEP; t+=dt, step++)
    {
        dt = localTimeStepping(); //当地时间步法，注意取全域最大当地时间步
        solve();                  //求解

        runInfo(t,dt,step, residual);
        if(residual<RESIDUAL_LIMIT)
        {
            cout<<"\nConverged!"<<endl
		    	<<"\nmaxI= "<<maxI<<"\tmaxJ= "<<maxJ<<endl
		    	<<"Final residual ="<<residual<<endl
		    	<<"Loop times = "<<step<<endl
		    	<<"Final field has been written into \"result.dat\" \n"
		    	<<"History of residual has been written into \"residual.plt\"";
            fstream fout("result.txt");
    	    print(Q, fout);
		    break;
        }
    }

    if(t >= STOP_TIME || step>=STOP_STEP)
	{
		cout<<"\nFail to converge! in "<<step<<" times\n"
		<<"Final residual ="<<residual<<"\n\n\n";
	}
    cout<<"\nWall time = "<<(double)clock()/CLOCKS_PER_SEC<<" s"<<endl;
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
    cout<<"mesh is written in \"mesh.dat\""<<endl;

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

void print(Field aField, fstream &fout)
{
    
    for(int j=0; j<=maxJ; j++)
    {
        for(int i=0; i<=maxI; i++)
        {
            for(int k: aField[i][j])
                fout<< aField[i][j][k]<<" ";
            fout<<endl;
        }
    }
    cout<<"\nWall time = "<<(double)clock()/CLOCKS_PER_SEC<<" s"<<endl;
}

/*******************other utility funcs*****************/
static inline void runInfo(double t, double dt, int step)
{
    cout<<"\tstep= "<<step
        <<"\tt="    <<t 
		<<"\tdt= "  <<dt<<endl;
}

inline double calMa(double T, double magU)
{
    return magU/(20.045*sqrt(T));
}

inline double calSoundSpeed(double T)
{
    return 20.045*sqrt(T);
}