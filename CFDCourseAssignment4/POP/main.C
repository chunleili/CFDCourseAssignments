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
#define caseNo (1)
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


typedef XY     MeshPoint[maxI+1][maxJ+1];        //用于存储网格点坐标
typedef double Field[maxI+1][maxJ+1][4];         //向量场,用于定义Q对象
typedef double ScalarField[maxI+1][maxJ+1];      //标量场,用于p,rho等场对象
typedef XY     Tensor[maxI+1][maxJ+1][4];        //张量场,用于定义F对象
typedef double Vector[4];                        //表示某一单元格的参数
typedef unsigned const Index;                    //用于传递编号,只读

/***************************declare the utility funcs  **************************/
AERO AeroConvert(Vector vec);
double safeSqrt(double xx);
void toFlux(Vector Q, Vector F);
/***************************define the global variable &funcs **************************/
MeshPoint mesh;
double volume, S1, S4, x1,x2,x3,x4, y1,y2,y3,y4; 
void genMesh();
void printMesh();
void cellGeometry();
XY   getN1();
XY   getN4();
XY   N1,N4;

Field Q;
Tensor Fc;
void init1();
void init2();
void BC1();
void BC2();
void solve();
double LTS();
double dt;

unsigned I, J;
double rhoL, rhoR, uL, uR, vL, vR, HL, HR, pL, pR;   
double lambda1,lambda2,lambda3;
double Fc1, Fc2, Fc3, Fc4;
Vector delF1, delF234, delF5;

void setFlux1();
void setFlux4();
void MUSCL(Field const U,  Vector UR, Vector UL);
void MUSCL(ScalarField const U,  double & UR, double & UL);
void Harten();
double calVcv1(double u_, double v_);
double calVcv4(double u_, double v_);
double calF1toF5(double Vcv);


/********************************Mesh**************************************/
void genMesh()
{
    const int block1=100, block2=100;
    const double dx=1.0/500;

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
}


void printMesh()
{
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

void cellGeometry()
{
    //从左下开始逆时针编号,左下点代表1,右下2,右上3,左上4
    //左下代表本单元格坐标
    x1=mesh[I  ][J  ].x;
    x2=mesh[I+1][J  ].x;
    x3=mesh[I+1][J+1].x;
    x4=mesh[I  ][J+1].x;

    y1=mesh[I  ][J  ].y;
    y2=mesh[I+1][J  ].y;
    y3=mesh[I+1][J+1].y;
    y4=mesh[I  ][J+1].y;

    S1=safeSqrt( (x1-x2)*(x1-x2)+(y2-y1)*(y2-y1) );//下侧面积S1
    S4=safeSqrt( (x4-x1)*(x4-x1)+(y1-y4)*(y1-y4) );//左侧面积S4

    N1.x=(y2-y1)/S1;
    N1.y=(x1-x2)/S1;

    N4.x=(y1-y4)/S4;
    N4.y=(x4-x1)/S4;

    volume=0.5*( (x1-x3)*(y2-y4)+ (x4-x2)*(y1-y3) );//计算单元体体积, 厚度取为1
}


//LTS=LocalTimeStepping,当地时间步法,返回全局的时间步
double LTS()
{
    double min=1e10;
    double dtLocal;
	AERO xxx;
    double lambda;
    
    forAll(
    	xxx=AeroConvert(Q[i][j]);
        lambda=(xxx.Vcv + xxx.c)*1;//咱先这样,一会再改

        dtLocal=CFL*1/lambda;

        if(dtLocal<min)
            min=dtLocal;
    );
    return min;
}

void init1()
{
    //初始全部给1.8Ma对应的速度, 压力给大气压101325, 静温300
    //声速340.68
    //速度u=613, v=0, 密度rho=p/RT=1.17
    forAll(
        Q[i][j][0]=1.17;
        Q[i][j][1]=717.21; //1.17*613
        Q[i][j][2]=0;
        Q[i][j][1]=215713; // 287/0.4*300+613;
    );
}

void init2()
{
    //初始全部给1.5Ma对应的速度, 压力给大气压101325, 静温300, 
    //速度u=511.02, v=0, 密度rho=p/RT=1.17
    forAll(
        Q[i][j][0]=1.17;
        Q[i][j][1]=597.89; //1.17*511.02
        Q[i][j][2]=0;
        Q[i][j][1]=215713; // 287/0.4*300+613;
    );
}

void BC1()
{

}

void BC2()
{

}

//三阶显式RungeKutta法
void solve()
{
    Tensor Fc;   //对流通量, 右侧x和上侧y的, 向外为正
    const double alpha[3]={0.1918, 0.4929, 1.0};
    //先定义Q0,用于保存原始的Q
    Field Q0;
    forEach(
            Q0[i][j][k] = Q[i][j][k];
    );
    
    //后面每一步都先计算残差, 后根据RK公式更新W
    Vector R;//R for Residual, 残差
    //double dV;
    //XY dS;
    for(unsigned a=0;a<=2;a++)   //a代表荣格库塔法的每一步
        for (unsigned I = 1; I <= maxI-2; I++)//I,J为单元编号
            for(unsigned J = 1; J <= maxJ-2; J++)
            {  
                //利用Roe格式计算通量
                Roe(I, J);
                //计算残差
                for (unsigned k=0; k<3; k++)
                {   
                    Fc[I][J][k]=roe.getFlux();
                    //R代表离开单元格的通量的矢量和, =右侧+左侧+上侧+下侧
                    //Fc[I][J].x/y分别代表左侧与下侧通量, 由于与坐标系正方向反向, 故取负号
                    //右侧与上侧通量分别由临近单元格左侧与下侧取负号得来, 取两次负号变为正
                    R[k]=-Fc[I+1][J][k].x + Fc[I][J][k].x - Fc[I][J+1][k].y + Fc[k][I][J].y;
                     
                    //利用荣格库塔法计更新流场
                    Q[I][J][k] = Q0[I][J][k] - alpha[a] * dt / dV * R[k];
                }
            }
}


/**********************Roe********************************/
Roe::Roe(Index II, Index JJ):I(II), J(JJ)
{}

//Roe格式计算对流通量
void setFlux1()
{  
    //先定义流场各变量
    ScalarField rho, u, v,  p, H;

    rho[I][J]=Q[I][J][0];
    u[I][J]=Q[I][J][1]/rho[I][J];
    v[I][J]=Q[I][J][1]/rho[I][J];
    p[I][J]=(GAMMA-1)* (Q[I][J][2] - rho[I][J]* (SQ(u[I][J])+SQ(v[I][J])) *0.5);
    H[I][J]=Q[I][J][3]/rho[I][J]+p[I][J]/rho[I][J];

    //然后利用MUSCL分裂这些变量, 并且求出Roe平均量
    Vector QL, QR;

    MUSCL1(rho, I, J, rhoR, rhoL);
    MUSCL1(p,   I, J, pR,   pL  );    //注意不要越界~
    MUSCL1(u,   I, J, uR,   uL  );
    MUSCL1(v,   I, J, vR,   vL  );
    MUSCL1(H,   I, J, HR,   HL  );
    MUSCL1(Q,   I, J, QR,   QL  );

    //计算Roe平均量
    const double denoLR=safeSqrt(rhoL)+ safeSqrt(rhoR);
    const double L=safeSqrt(rhoL)/ denoLR;//定义两个系数
    const double R=safeSqrt(rhoR)/ denoLR;
    
    const double rho_=safeSqrt(rhoL*rhoR);
    const double u_  =uL*L+uR*R;
    const double v_  =vL*L+vR*R;
    const double H_  =HL*L+HR*R;
    const double q_2 =u_*u_+v_*v_;
    const double c_  =safeSqrt((GAMMA-1)*(H_-q_2/2));

    XY N1=mesh.
    double Vcv_=N1.x * u_ +N1.y* v_;
    //Vcv代表contravirant velocity 逆变速度

    //计算各个del值
    const double delP  =pR-pL;
    const double delRho=rhoR-rhoL;
    const double delu  =uR-uL;

    lambda1=Herten(fabs(Vcv_-c_));
    lambda2=Herten(fabs(Vcv_   ));
    lambda3=Herten(fabs(Vcv_+c_));   

    double delVcv  =Vcv1_R-vcv1_L;
    //求出Roe矩阵相关值
        //delu,delv,delw代表三个分量, delVcv代表大写Vcv的delta
        //lambda123分别是是Vcv-c, Vcv, Vcv+c
    const double delP  =pR-pL;
    const double delRho=rhoR-rhoL;
    const double delu  =uR-uL;
    const double delVcv=VcvR-VcvL;
    XY N1=mesh.getN1(I,J);
    F1[0]  = lambda1 * (delP-rho_*c_*delVcv)/(2*c_*c_) * 1;
    F1[1]  = lambda1 * (delP-rho_*c_*delVcv)/(2*c_*c_) * (u_-c_*N1.x);
    F1[2]  = lambda1 * (delP-rho_*c_*delVcv)/(2*c_*c_) * (u_-c_*N1.y);
	F1[3]  = lambda1 * (delP-rho_*c_*delVcv)/(2*c_*c_) * (H_-c_*Vcv_);
  
    F234[0]= lambda2 * ( (delRho-delP/(c_*c_))*1  + 0); 
    F234[1]= lambda2 * ( (delRho-delP/(c_*c_))*u_ 
                   + rho_ *(delu   -  delVcv)    );
    F234[2]= lambda2 * ( (delRho-delP/(c_*c_))*u_ 
                   + rho_ *(delu   -  delVcv)    );                
    F234[3]= lambda2 * ( (delRho-delP/(c_*c_))*u_
                              +rho_ *(u_*delu - Vcv_*delVcv) );
     
    F5[0]  = lambda3 * (delP+rho_*c_*delVcv)/(2*c_*c_) *1;
    F5[1]  = lambda3 * (delP+rho_*c_*delVcv)/(2*c_*c_) *(u_+c_);
    F5[2]  = lambda3 * (delP+rho_*c_*delVcv)/(2*c_*c_) *(u_+c_);    
    F5[3]  = lambda3 * (delP+rho_*c_*delVcv)/(2*c_*c_) *(H_+c_*Vcv_);

    Vector FcR, FcL;
    toFlux(QR, FcR);
    toFlux(QL, FcL);


    //最后算出单元IJ的对流通量Fc
    for(unsigned k=0; k<3; k++)
    {
        Fc[I][J][k].x=0.5*( FcR[k] - FcL[k] - (F1[k]+F234[k]+F5[k]) );
    }
}

//带限制器的MUSCL插值,前两个参数是输入,后两个输出, 系数k^为1/3
void Roe::MUSCL1(Field const U, Vector UR, Vector UL)
{
    if(I+2>maxI || I-1<0) {cout<<"\n\nout of bound!!!\n\n";}
    const double dV=mesh.getVolume(I,J);
    const double epsilon= pow(dV,1.0/3); //限制器参数epsilon与几何尺寸相关
    
    double aR=U[I+2][J]- U[I+1][J],   bR=U[I+1][J] - U[I  ][J];
    double aL=U[I+1][J]- U[I  ][J],   bL=U[I  ][J] - U[I-1][J];
    double a, b;
    for(unsigned k=0;k<3;k++)
    {
        a=aR; b=bR;
        double deltaR=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
        a=aL; b=bL;
        double deltaL=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);

        UR[k]=U[I+1][J][k]-0.5*deltaR;
        UL[k]=U[I  ][J][k]+0.5*deltaL;
    }
}
//重载用于标量的带限制器的MUSCL插值函数 ,前两个参数是输入,后两个输出, 系数k^为1/3
void Roe::MUSCL1(ScalarField const U, double & UR, double & UL)
{
    if(I+2>maxJ || I-1<0) {cout<<"\n\nout of bound!!!\n\n";}
    const double dV=mesh.getVolume(I,J);
    const double epsilon= pow(dV,1.0/3); //限制器参数epsilon与几何尺寸相关

    double aR=U[I+2][J]- U[I+1][J],   bR=U[I+1][J] - U[I  ][J];
    double aL=U[I+1][J]- U[I  ][J],   bL=U[I  ][J] - U[I-1][J];

    double a=aR, b=bR;
    double deltaR=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    a=aL, b=bL;
    double deltaL=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    
    UR=U[I][J+1]-0.5*deltaR;
    UL=U[I][J]  +0.5*deltaL;
}





/***********************flux4****************************/
/***********************flux4****************************/
/***********************flux4****************************/
/***********************flux4****************************/

//重载用于标量的带限制器的MUSCL插值函数 ,前两个参数是输入,后两个输出, 系数k^为1/3
void Roe::MUSCL4(ScalarField const U double & UR, double & UL)
{
    if(J+2>maxJ || J-1<0) {cout<<"\n\nout of bound!!!\n\n";}
    const double epsilon= dx; //限制器参数epsilon与几何尺寸相关

    double aR=U[I][J+2]-U[I][J+1], bR=U[I][J+1]-U[I][J];
    double aL=U[I][J+1]-U[I][J],   bL=U[I][J]  -U[I][J-1];

    double a=aR, b=bR;
    double deltaR=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    a=aL, b=bL;
    double deltaL=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    
    UR=U[I][J+1]-0.5*deltaR;
    UL=U[I][J]  +0.5*deltaL;
}

double Herten(double lambda)
{
    double delta=0.1*c_; //熵修正, Harten's entropy correction
    if(fabs(lambda)<=delta)
         lambda=(lambda*lambda+ delta*delta) /(2*delta);
    return lambda;
}





int main()
{  
	int caseNo=1;
	cout<<"Case 1 for inlet 1.8Ma; 2 for 1.5Ma\n Case: "<<caseNo<<endl;
	Mesh mesh;
    FlowField field1(1);
	
    double residual=1.0;
    for (int step=0;  step<=STOP_STEP; step++)
    {    
        field1.solve(mesh);                  //求解
    }

    return 0;
}
