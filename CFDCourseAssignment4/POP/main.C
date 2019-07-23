/***************************include & namespace  **************************/
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<cstdio>
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
#define check(val) if(val<0) cout<<"\n*Where: "<<#val<<" is negative! "<<val<<endl
#define SQ(a) ((a)*(a))

#define caseNo (1)  //case1 1.8Ma, case2 1.5Ma
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
typedef double Field[maxI+1][maxJ+1][4];         //场,用于定义Q,Fc1等对象
typedef double ScalarField[maxI+1][maxJ+1];      //标量场,用于p,rho大小等场对象
typedef XY     VectorField[maxI+1][maxJ+1];      //向量场,用于定义面的单位法量 
typedef double Vector[4];                        //表示某一单元格的参数
typedef unsigned const Index;                    //用于传递编号,只读

/***************************declare the utility funcs  **************************/
void aeroConvert(Index i, Index j);
double safeSqrt(double xx);
void toFlux(Vector Q, Vector F);
void print();
/***************************define the global variable &funcs **************************/
MeshPoint mesh;
ScalarField  volume, S1,S2,S3,S4; //面积, 逆时针顺序, 依次为下右上左
VectorField  N1,N2,N3,N4;
void genMesh();
void printMesh();
void cellGeometry(); 


Field Q, Fc1,Fc2,Fc3,Fc4;
Field R; //R for Residual, 残差
//对流通量, 逆时针顺序依次为下右上左, 向外为正, 只储存其大小, 方向由N1~N4给定
void init1();
void init2();
void BC1();
void BC2();
void solve();
double LTS();
double dt;

unsigned I, J, step;
ScalarField rho, u, v,  p, H, Vcv2, Vcv3; //Vcv代表contravirant velocity 逆变速度
//double lambda1,lambda2,lambda3;

void roeFlux3();
void roeFlux2();
void MUSCL3(Field const U,  Vector UR, Vector UL);
void MUSCL3(ScalarField const U,  double & UR, double & UL);
void MUSCL2(Field const U,  Vector UR, Vector UL);
void MUSCL2(ScalarField const U,  double & UR, double & UL);
double Harten(double lambda);


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
    cout<<"\nmesh is written in \"mesh.dat\""<<endl;
}

void cellGeometry()
{
    double x1,x2,x3,x4, y1,y2,y3,y4;
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

    S1[I][J]=safeSqrt( SQ(x1-x2)+SQ(y1-y2) );//下侧面积S1
    S2[I][J]=safeSqrt( SQ(x2-x3)+SQ(y2-y3) );//右侧面积S2
    S3[I][J]=safeSqrt( SQ(x3-x4)+SQ(y3-y4) );//上侧面积S3
    S4[I][J]=safeSqrt( SQ(x4-x1)+SQ(y4-y1) );//左侧面积S4
    check(SQ(x1-x2)+SQ(y1-y2));
    check(SQ(x2-x3)+SQ(y2-y3));
    check(SQ(x3-x4)+SQ(y3-y4));
    check(SQ(x4-x1)+SQ(y4-y1));

    N1[I][J].x=(y2-y1)/S1[I][J];
    N1[I][J].y=(x1-x2)/S1[I][J];

    N2[I][J].x=(y3-y2)/S2[I][J];
    N2[I][J].y=(x2-x3)/S2[I][J];

    N3[I][J].x=(y4-y3)/S3[I][J];
    N3[I][J].y=(x3-x4)/S3[I][J];

    N4[I][J].x=(y1-y4)/S4[I][J];
    N4[I][J].y=(x4-x1)/S4[I][J];

    volume[I][J]=0.5*( (x1-x3)*(y2-y4)+ (x4-x2)*(y1-y3) );//计算单元体体积, 厚度取为1
}

/**********************LTS********************************/
//LTS=LocalTimeStepping,当地时间步法,返回当地时间步
double LTS()
{
    double dtLocal;
    double lambda11, lambda44, S11, S44;
    XY N11, N44;

    for(unsigned i=0; i<=maxI; i++)
        for(unsigned j=0; j<=maxJ; j++) 
        {
            N11.x=0.5*(N1[I][J].x-N3[I][J].x);
            N11.y=0.5*(N1[I][J].y-N3[I][J].y);

            N44.x=0.5*(N4[I][J].x-N4[I][J].x);
            N44.y=0.5*(N4[I][J].y-N4[I][J].y);

            S11=0.5*(S1[I][J]+S3[I][J]);
            S44=0.5*(S2[I][J]+S4[I][J]);

            lambda11=u[I][J]*N11.x+v[I][J]*N11.y;
            lambda44=u[I][J]*N44.x+v[I][J]*N44.y;

            dtLocal=CFL*volume[I][J]/(lambda11*S11+lambda44*S44);
        }
  
    return dtLocal;
}

/**********************init &BC********************************/
void init1()
{
    //初始全部给1.8Ma对应的速度, 压力给大气压101325, 静温300
    //声速c=sqrt(1.4*287*300)=347.1887
    //速度u=624.9397, v=0, 密度rho=p/RT=1.176829
    forAll(

            Q[i][j][0] = 1.176829;  //101325/(287*300)
            Q[i][j][1] = 735.4473; //1.17*624
            Q[i][j][2] = 0;
            Q[i][j][3] = 483117.6; // 101325/0.4+0.5*1.176829*625*625;
            aeroConvert(i, j);
        
    );
}

void init2()
{
    //初始全部给1.5Ma对应的速度, 压力给大气压101325, 静温300, 
    //速度u=1.5*347.2=520.783, v=0, 密度rho=p/RT=1.176829
    forAll(
        Q[i][j][0]=1.176829;
        Q[i][j][1]=612.873; //1.176829*520.783
        Q[i][j][2]=0;
        Q[i][j][3]=412899.3; //  101325/0.4+0.5*1.176829*520.783^2;
        aeroConvert(i,j);
    );
}

void BC1()
{
    //左右: 入口维持1.8Ma不用管, 出口用0梯度递推出来;
    for(unsigned j=1; j<maxJ; j++)
        for(unsigned k=0; k<4; k++)
        {
            Q[maxI][j][k]=Q[maxI-1][j][k]=Q[maxI-2][j][k];
        }
    //上下是壁面, 流向不变, 用0梯度外推, 垂向直接给0
    //密度,静压不变, u不变, v变为0, E随之变化
    for(unsigned i=1; i<maxI; i++)
    {
        Q[i][maxJ][0]=Q[i][maxJ-1][0]=Q[i][maxJ-2][0];
        Q[i][maxJ][1]=Q[i][maxJ-1][1]=0;
        Q[i][maxJ][2]=Q[i][maxJ-1][2]=Q[i][maxJ-2][2];
        Q[i][maxJ][3]=Q[i][maxJ-1][3]=Q[i][maxJ-2][2]-0.5*SQ(Q[i][maxJ-2][1])/Q[i][maxJ-2][0];

        Q[i][0][0]=Q[i][1][0];
        Q[i][0][1]=0;
        Q[i][0][2]=Q[i][1][2];
        Q[i][0][3]=Q[i][1][3]-0.5*SQ(Q[i][1][1])/Q[i][1][0];
    }
        
}

void BC2()
{

}

/**********************solve********************************/
//三阶显式RungeKutta法
void solve()
{
    const double alpha[3]={0.1918, 0.4929, 1.0};
    //先定义Q0,用于保存原始的Q
    Field Q0;
    forEach(
            Q0[i][j][k] = Q[i][j][k];
    );
    
    //后面每一步都先计算残差, 后根据RK公式更新W
    dt=1;//dt先给定一个值;
    double dtLocal;
    for(unsigned a=0;a<=2;a++)   //a代表荣格库塔法的每一步
        for (I = 1; I <= maxI-2; I++)//I,J为单元编号, 只在此处变动!!
            for(J = 1; J <= maxJ-2; J++)
            {   
                cout<<"\n\n******************************************************\n";
                cout<<"step= "<<step<<" a= "<<a<<"  I= "<<I<<" J= "<<J<<endl;
                cellGeometry();
                aeroConvert(I,J);
                //利用Roe格式计算通量
                roeFlux3();
                roeFlux2();

                dtLocal=LTS();
                if(dtLocal<dt)
                    dt=dtLocal;//定常问题不需要dt, 但是荣格库塔法需要

                //计算残差
                for (unsigned k=0; k<4; k++)
                {   
                    //R代表离开单元格的通量的矢量和, =右侧+左侧+上侧+下侧
                    //左侧与下侧通量分别由临近单元格右侧与上侧取负号得来
                    R[I][J][k]=Fc1[I+1][J][k] + Fc2[I][J][k] + Fc3[I][J+1][k] + Fc4[I][J][k];
                     
                    //利用荣格库塔法计更新流场
                    Q[I][J][k] = Q0[I][J][k] - alpha[a] * dt / volume[I][J] * R[I][J][k];
                }
            }
}


/**********************Roe********************************/
//Roe格式计算对流通量
/***********************flux3****************************/
/***********************flux3****************************/
void roeFlux3()
{   
    double nx=N3[I][J].x, ny=N3[I][J].y;
    
    //利用MUSCL分裂流场变量, 并且求出Roe平均量 
    double rhoR, rhoL, pR, pL, uR, uL, vR, vL, HR, HL;
    double VcvL, VcvR;
    Vcv3[I][J] =nx * u[I][J] + ny * v[I][J];

    MUSCL3(rho,  rhoR, rhoL);    //带限制器的三点MUSCL插值
    MUSCL3(p,    pR,   pL  );    //注意不要越界~
    MUSCL3(u,    uR,   uL  );
    MUSCL3(v,    vR,   vL  );
    MUSCL3(H,    HR,   HL  );
    MUSCL3(Vcv3, VcvL, VcvR);   //为了方便代码重用, 分裂Vcv3后记为VcvL和VcvR

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

    const double Vcv_=nx * u_ + ny * v_;
    check(rhoL*1.0);
    check(rhoR);
    check((GAMMA-1)*(H_-q_2/2));
    //计算lambda
    const double lambda1=Harten(fabs(Vcv_-c_));
    const double lambda2=Harten(fabs(Vcv_   ));
    const double lambda3=Harten(fabs(Vcv_+c_));   

    //求出Roe矩阵相关值
        //delu,delv,delw代表三个分量, delVcv代表大写Vcv的delta
        //lambda123分别是是Vcv-c, Vcv, Vcv+c
    const double delP   =pR-pL;
    const double delRho =rhoR-rhoL;
    const double delu   =uR-uL;
    const double delv   =vR-vL;
    const double delVcv =VcvR-VcvL;


    Vector delF1, delF234, delF5;

    const double coeff1=lambda1 * (delP-rho_*c_*delVcv)/(2*SQ(c_));//定义系数以减少计算量
    delF1[0]  = coeff1 * 1;
    delF1[1]  = coeff1 * (u_-c_*nx);
    delF1[2]  = coeff1 * (u_-c_*ny);
	delF1[3]  = coeff1 * (H_-c_*Vcv_);

    const double coeff2=lambda2 *( delRho-delP/SQ(c_) );//定义系数以减少计算量 
    const double coeff3=lambda2 * rho_;
    delF234[0]= coeff2; 
    delF234[1]= coeff2 * u_  + coeff3*(delu - delVcv*nx);
    delF234[2]= coeff2 * v_  + coeff3*(delv - delVcv*ny);                
    delF234[3]= coeff2 * q_2 + coeff3*(u_*delu + v_*delv - Vcv_*delVcv);
    
    const double coeff5=lambda3 * (delP+rho_*c_*delVcv)/(2*SQ(c_));//定义系数以减少计算量
    delF5[0]  = coeff5 *1;
    delF5[1]  = coeff5 *(u_+c_*nx);
    delF5[2]  = coeff5 *(v_+c_*ny);    
    delF5[3]  = coeff5 *(H_+c_*Vcv_);

    //分裂Fc
    Vector QL, QR, Fc3R, Fc3L;
    MUSCL3(Q,  QR,   QL  );
    toFlux(QR, Fc3R);
    toFlux(QL, Fc3L);

    //最后算出单元IJ的对流通量Fc1 与Fc3
    for(unsigned k=0; k<4; k++)
    {
        Fc3[I][J][k]=0.5*( Fc3R[k] + Fc3L[k] - (delF1[k]+delF234[k]+delF5[k]) );
        Fc1[I][J][k]=-Fc3[I][J-1][k];
    }
}

//N1方向的(下), 带限制器的三点MUSCL插值,第一个参数是输入,后两个输出, 系数k^为1/3
void MUSCL3(Field const U, Vector UR, Vector UL)
{
    if(J+2>maxJ || J-1<0) {cout<<"\n\nout of bound!!!\n\n";}
    const double epsilon= pow(volume[I][J],1.0/3); //限制器参数epsilon与几何尺寸相关

    double aR=U[I][J+2]-U[I][J+1], bR=U[I][J+1]-U[I][J];
    double aL=U[I][J+1]-U[I][J],   bL=U[I][J]  -U[I][J-1];
    double a, b;
    a = aR; b = bR;
    double deltaR=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    a = aL; b = bL;
    double deltaL=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);

    for(unsigned k=0; k<4; k++)
    {
        UR[k] = U[I][J + 1][k] - 0.5 * deltaR;
        UL[k] = U[I][J    ][k] + 0.5 * deltaL;
    }
}

//N1方向(下), 重载用于标量的带限制器的三点MUSCL插值函数 ,第一个参数是输入,后两个输出, 系数k^为1/3
void MUSCL3(ScalarField const U, double & UR, double & UL)
{
    if(J+2>maxJ || J-1<0) {cout<<"\n\nout of bound!!!\n\n";}
    const double epsilon= pow(volume[I][J],1.0/3); //限制器参数epsilon与几何尺寸相关

    double aR=U[I][J+2]-U[I][J+1], bR=U[I][J+1]-U[I][J];
    double aL=U[I][J+1]-U[I][J],   bL=U[I][J]  -U[I][J-1];

    double a=aR, b=bR;
    double deltaR=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    a=aL, b=bL;
    double deltaL=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    
    UR=U[I][J+1]-0.5*deltaR;
    UL=U[I][J]  +0.5*deltaL;
}


/***********************flux2****************************/
/***********************flux2****************************/
void roeFlux2()
{   
    double nx=N2[I][J].x, ny=N2[I][J].y;
    
    //利用MUSCL分裂流场变量, 并且求出Roe平均量 
    double rhoR, rhoL, pR, pL, uR, uL, vR, vL, HR, HL;
    double VcvL, VcvR;
    Vcv2[I][J] =nx * u[I][J] + ny * v[I][J];

    MUSCL2(rho,  rhoR, rhoL);    //带限制器的三点MUSCL插值
    MUSCL2(p,    pR,   pL  );    //注意不要越界~
    MUSCL2(u,    uR,   uL  );
    MUSCL2(v,    vR,   vL  );
    MUSCL2(H,    HR,   HL  );
    MUSCL2(Vcv2, VcvL, VcvR);   //为了方便代码重用, 分裂Vcv4后记为VcvL和VcvR


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

    const double Vcv_=nx * u_ + ny * v_;

    check(rhoL);
    check(rhoR);
    check((GAMMA-1)*(H_-q_2/2));

   //计算lambda
    const double lambda1=Harten(fabs(Vcv_-c_));
    const double lambda2=Harten(fabs(Vcv_   ));
    const double lambda3=Harten(fabs(Vcv_+c_));   

    //求出Roe矩阵相关值
        //delu,delv,delw代表三个分量, delVcv代表大写Vcv的delta
        //lambda123分别是是Vcv-c, Vcv, Vcv+c
    const double delP   =pR-pL;
    const double delRho =rhoR-rhoL;
    const double delu   =uR-uL;
    const double delv   =vR-vL;
    const double delVcv =VcvR-VcvL;


    Vector delF1, delF234, delF5;

    const double coeff1=lambda1 * (delP-rho_*c_*delVcv)/(2*SQ(c_));//定义系数以减少计算量
    delF1[0]  = coeff1 * 1;
    delF1[1]  = coeff1 * (u_-c_*nx);
    delF1[2]  = coeff1 * (u_-c_*ny);
	delF1[3]  = coeff1 * (H_-c_*Vcv_);

    const double coeff2=lambda2 *( delRho-delP/SQ(c_) );//定义系数以减少计算量 
    const double coeff3=lambda2 * rho_;
    delF234[0]= coeff2; 
    delF234[1]= coeff2 * u_  + coeff3*(delu - delVcv*nx);
    delF234[2]= coeff2 * v_  + coeff3*(delv - delVcv*ny);                
    delF234[3]= coeff2 * q_2 + coeff3*(u_*delu + v_*delv - Vcv_*delVcv);
    
    const double coeff5=lambda3 * (delP+rho_*c_*delVcv)/(2*SQ(c_));//定义系数以减少计算量
    delF5[0]  = coeff5 *1;
    delF5[1]  = coeff5 *(u_+c_*nx);
    delF5[2]  = coeff5 *(v_+c_*ny);    
    delF5[3]  = coeff5 *(H_+c_*Vcv_);

    //分裂Fc
    Vector QL, QR, Fc2R, Fc2L;
    MUSCL2(Q,  QR,   QL  );
    toFlux(QR, Fc2R);
    toFlux(QL, Fc2L);

    //最后算出单元IJ的对流通量Fc4 与Fc2
    for(unsigned k=0; k<4; k++)
    {
        Fc2[I][J][k]=0.5*( Fc2R[k] + Fc2L[k] - (delF1[k]+delF234[k]+delF5[k]) );
        Fc4[I][J][k]=-Fc4[I-1][J][k];
    }
}

//N4方向(左), 重载用于标量的带限制器的MUSCL插值函数 ,前两个参数是输入,后两个输出, 系数k^为1/3
void MUSCL2(Field const U, Vector UR, Vector UL)
{
    if(I+2>maxI || I-1<0) {cout<<"\n\nout of bound!!!\n\n";}
    
    const double epsilon= pow(volume[I][J], 1.0/3); //限制器参数epsilon与几何尺寸相关
    
    double aR=U[I+2][J]- U[I+1][J],   bR=U[I+1][J] - U[I  ][J];
    double aL=U[I+1][J]- U[I  ][J],   bL=U[I  ][J] - U[I-1][J];
    double a, b;
    a=aR; b=bR;
    double deltaR=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    a=aL; b=bL;
    double deltaL=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    for(unsigned k=0;k<4;k++)
    {
        UR[k]=U[I+1][J][k]-0.5*deltaR;
        UL[k]=U[I  ][J][k]+0.5*deltaL;
    }
}

//N4方向(左), 重载用于标量的带限制器的三点MUSCL插值函数 ,第一个参数是输入,后两个输出, 系数k^为1/3
void MUSCL2(ScalarField const U, double & UR, double & UL)
{
    if(I+2>maxI || I-1<0) {cout<<"\n\nout of bound!!!\n\n";}
    
    const double epsilon= pow(volume[I][J], 1.0/3); //限制器参数epsilon与几何尺寸相关
    
    double aR=U[I+2][J]- U[I+1][J],   bR=U[I+1][J] - U[I  ][J];
    double aL=U[I+1][J]- U[I  ][J],   bL=U[I  ][J] - U[I-1][J];
    double a, b;
    a=aR; b=bR;
    double deltaR=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    a=aL; b=bL;
    double deltaL=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);

    UR=U[I+1][J]-0.5*deltaR;
    UL=U[I  ][J]+0.5*deltaL;

}

double Harten(double lambda)
{
    double c=safeSqrt(GAMMA*p[I][J]/rho[I][J]);
    check(GAMMA*p[I][J]/rho[I][J]);

    double delta=0.1*c; //熵修正, Harten's entropy correction
    if(fabs(lambda)<=delta)
         lambda=(lambda*lambda+ delta*delta) /(2*delta);
    return lambda;

}


/***************************utility  **************************/
//Q与Fc之间的转化
void toFlux(Vector Q, Vector F)
{
    F[0] = Q[1];
    F[1] = Q[0] * u[I][J] * u[I][J] + p[I][J];
    F[2] = (Q[2] + p[I][J]) * u[I][J];
}

//气动参数转换
void aeroConvert(Index i, Index j)
{
    rho[i][j]=Q[i][j][0];
    u[i][j]=Q[i][j][1]/rho[i][j];
    v[i][j]=Q[i][j][2]/rho[i][j];
    p[i][j]=(GAMMA-1)* (Q[i][j][3] - rho[i][j]* (SQ(u[i][j])+SQ(v[i][j])) *0.5);
    H[i][j]=Q[i][j][3]+p[i][j];
}

double safeSqrt(double xx)
{
    if(xx<0) 
    {
        cout<<"\n##Erro sqrt! Value "<<xx<<" is negative!\n\n";
        getchar();
    }
    return sqrt(xx);
}

void print()
{
    FILE *fp1, *fp2, *fp3;
    fp1=fopen("pressure.dat", "w");
    fprintf(fp1,"x       y       pressure\n");
    for(unsigned i=0; i<=maxI; i++)
        for(unsigned j=0; j<=maxJ; j++)
        {
            fprintf(fp1, "%.5f %.5f %.0f\n", mesh[i][j].x, mesh[i][j].y, p[I][J] );
        } 
    fclose(fp1);

    fp2=fopen("Ma.dat", "w");
    fprintf(fp2,"x       y       Ma\n");
    for(unsigned i=0; i<=maxI; i++)
        for(unsigned j=0; j<=maxJ; j++)
        {
            double Ma=safeSqrt( ( SQ(u[i][j])+SQ(v[i][j]) ) / (GAMMA*p[i][j]/rho[i][j]) );
            check(( SQ(u[i][j])+SQ(v[i][j]) ) / (GAMMA*p[i][j]/rho[i][j]) );
            fprintf(fp2, "%.5f %.5f %.5f\n", mesh[i][j].x, mesh[i][j].y, Ma );
        } 
    fclose(fp2);

    fp3=fopen("residual.dat", "w");
    fprintf(fp3,"step\tRho\tmagU\tE\n");
    double magU, residualRho=1e5, residualMagU=1e5, residualE=1e5;
    for (unsigned ss = 0; ss <= step; ss++)
    {
        for (unsigned i = 0; i <= maxI; i++)
            for (unsigned j = 0; j <= maxJ; j++)
            {
                magU = safeSqrt((SQ(R[i][j][1]) + SQ(R[i][j][2])) / SQ(R[i][j][0]));
                check((SQ(R[i][j][1]) + SQ(R[i][j][2])) / SQ(R[i][j][0]) );
                if (residualRho > fabs(R[i][j][0]))
                    residualRho = fabs(R[i][j][0]);
                if (residualMagU > magU)
                    residualMagU = magU;
                if (residualE > fabs(R[i][j][3]))
                    residualE = fabs(R[i][j][3]);
            }
        fprintf(fp3, "%d %7f %7f %7f\n", ss, residualRho, residualMagU, residualE);
    }

    fclose(fp2);
}       

/***************************main  **************************/
int main()
{  
	cout<<"\nCase 1 for inlet 1.8Ma; 2 for 1.5Ma\n\n Current Case: "<<caseNo<<endl;
    genMesh();//生成网格
    printMesh();//打印网格

	switch(caseNo)//初始化流场
    {
        case 1:init1(); break;
        case 2:init2(); break;
    }
    print();//打印结果
    cout<<"\nInitialzation done.\n";

    for (step=0;  step<=1; step++)
    {   
        solve();                  //求解
        switch(caseNo)            
        {            
            case 1: BC1(); break;
            case 2: BC2(); break;
        }
    }

    print();//打印结果
    return 0;
}
