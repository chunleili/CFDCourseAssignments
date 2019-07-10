/***************************include & funcs declare **************************/
#include<iostream>
#include<fstream>
#include<vector>
#include<utility>
#include<cmath>
#include<cstdlib>
using namespace std;

pair<int,int> genMesh();
void init1();
void init2();
double localTimeStepping();
void solve();
void print();

/***************************define the consts ********************************/

const double stop=1.0;
const int stopStep=1000;

/***************************define the type ********************************/
typedef struct xy
{
    public:
    double x;
    double y;
}xy;

/********************************main()**************************************/
int main()
{  
    pair<int,int> maxIJ=genMesh();              //生成网格
    double maxI=maxIJ.first, maxJ=maxIJ.second;
    
    vector<vector<double> > Q, F;
    
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

    print();                //打印 
    return 0;
}

/********************************funcs define**************************************/
pair<int,int> genMesh()
{

    const int maxI=400, maxJ=100;
    const int block1=100, block2=100, block3=maxI-block1-block2;

    double mesh[maxI][maxJ];
    fstream fout("mesh.txt");
    //for(auto &i: mesh)
        

    return make_pair(maxI,maxJ);
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

void print()
{

}