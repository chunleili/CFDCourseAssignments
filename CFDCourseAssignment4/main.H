/***************************include & namespace  **************************/
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<fstream>
#include<iomanip>
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

#define SQ(a) (a*a)
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

typedef struct FLUX
{
    public:
    double f1,f2,f3,f4;  
}FLUX;


typedef XY     MeshPoint[maxI+1][maxJ+1];        //用于存储网格点坐标
typedef double Field[maxI+1][maxJ+1][4];         //向量场,用于定义Q对象
typedef double ScalarField[maxI+1][maxJ+1];      //标量场,用于p,rho等场对象
typedef XY     Tensor[maxI+1][maxJ+1][4];        //张量场,用于定义F对象
typedef double Vector[4];                        //表示某一单元格的参数
typedef unsigned const Index;                    //用于传递编号,只读
/***************************declare the funcs  **************************/
AERO AeroConvert(Vector vec);
double safeSqrt(double xx);
void toFlux(Vector Q, Vector F);
/***************************define the class  **************************/
class Mesh
{
    MeshPoint mesh;
    XY N1,N4;
    public:
    Mesh();
    double getVolume(Index I, Index J);
    XY     getArea(Index I, Index J);
    XY     getN1(Index I, Index J);
    XY     getN4(Index I, Index J);
    void print();
};

class FlowField
{
    Mesh mesh;
    Field Q;
    Tensor Fc;

    void init1();
    void init2();
    void BC1();
    void BC2();

    public:
    FlowField(int caseNo);
    void solve();
    double LTS();
    double dt;
};

class Roe
{
    Index I, J;
    double rhoL, rhoR, uL, uR, vL, vR, HL, HR, pL, pR;   
    double lambda1,lambda2,lambda3;
    double Fc1, Fc4;
    Vector F1, F234, F5;

    void setFlux1();
    void setFlux4();
	void MUSCL(Field const U,  Vector UR, Vector UL);
	void MUSCL(ScalarField const U,  double & UR, double & UL);
    void Harten();
    double calVcv1(double u_, double v_);
    double calVcv4(double u_, double v_);
    double calF1toF5(double Vcv);

    public:
    Roe(Index II, Index JJ);
    XY getFlux() const;
};

