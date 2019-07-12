/***************************include & namespace  **************************/
#include<iostream>
#include<cmath>
using namespace std;
/***************************define the consts ********************************/

const int STOP_STEP=100;
const int maxI=400, maxJ=100;
const double RESIDUAL_LIMIT=1e-3;

/***************************define the type ********************************/
typedef struct XY
{
    public:
    double x;
    double y;
}XY;

typedef XY MeshPoint[maxI+1][maxJ+1];        //用于存储网格点坐标
typedef double Field[maxI+1][maxJ+1][4];     //向量场,用于定义Q对象
typedef double Tensor[maxI+1][maxJ+1][4][2]; //张量场,用于定义F对象

/***************************declare the funcs  **************************/
void genMesh();
void init1();
void init2();
double localTime();
void solve();
void print(Field aField, fstream &);

//utility functions
inline double calMa(double T, double magU);
inline double calSoundSpeed(double T);

inline double calMa(double T, double magU)
{
    return magU/(20.045*sqrt(T));
}

inline double calSoundSpeed(double T)
{
    return 20.045*sqrt(T);
}