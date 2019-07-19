/***************************include & namespace  **************************/
#include<iostream>
#include<cmath>
using namespace std;

/***************************MACRO             **************************/
#define forAll(codes)\
{\
    for(unsigned i=0; i<=maxI; i++)\
        for(unsigned j=0; j<=maxJ; j++) \
            {\
                codes\
            }\
}
#define forEach(codes)\
{\
    for(unsigned i=0; i<=maxI; i++)\
        for(unsigned j=0; j<=maxJ; j++) \
            for(unsigned k=0; k<=3; k++) \
                {\
                    codes\
                }\
}
/***************************define the consts ********************************/

const int STOP_STEP=100;
const int maxI=400, maxJ=100;
const double RESIDUAL_LIMIT=1e-3;
const double GAMMA=1.4;
const double CFL=0.69;
/***************************define the type ********************************/
typedef struct XY
{
    public:
    double x;
    double y;
}XY;

typedef struct AERO
{
    public:
    double rho,u,v,VV,p,T,c,Ma;  
}AERO;

typedef struct RESIDUAL
{
    public:
    double u,v,VV,p,E;  
}RESIDUAL;

typedef XY     MeshPoint[maxI+1][maxJ+1];        //用于存储网格点坐标
typedef double Field[maxI+1][maxJ+1][4];         //向量场,用于定义Q对象
typedef double ScalarField[maxI+1][maxJ+1];      //标量场,用于p,rho等场对象
typedef XY     Tensor[maxI+1][maxJ+1][4];        //张量场,用于定义F对象
typedef double Vector[4];                        //表示某一单元格的参数
typedef unsigned const Index;                    //用于传递编号,只读
/***************************declare the funcs  **************************/
double LTS(FlowField &Q);
AERO convert(Vector vec);
double safeSqrt(double xx);
/***************************declare the class  **************************/
class FlowField;
class Mesh
{
    MeshPoint mesh;
    public:
    Mesh();
    double getVolume(Index I, Index J);
    XY     getArea(Index I, Index J);
    friend void FlowField::solve();
    void print();
};

class FlowField
{
    Field Q;
    public:
    FlowField(int caseNo);
    void init1();
    void init2();
    void BC1();
    void BC2();
    void solve();
    friend double LTS();
};