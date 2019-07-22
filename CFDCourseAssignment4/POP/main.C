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
typedef double Field[maxI+1][maxJ+1][4];         //场,用于定义Q,Fc1等对象
typedef double ScalarField[maxI+1][maxJ+1];      //标量场,用于p,rho大小等场对象
typedef XY     VectorField[maxI+1][maxJ+1];      //向量场,用于定义面的单位法量 
typedef double Vector[4];                        //表示某一单元格的参数
typedef unsigned const Index;                    //用于传递编号,只读

/***************************declare the utility funcs  **************************/
void aeroConvert();
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

unsigned I, J;
ScalarField rho, u, v,  p, H, Vcv1, Vcv4; //Vcv代表contravirant velocity 逆变速度
double lambda1,lambda2,lambda3;

void roeFlux1();
void roeFlux4();
void MUSCL1(Field const U,  Vector UR, Vector UL);
void MUSCL1(ScalarField const U,  double & UR, double & UL);
void MUSCL4(Field const U,  Vector UR, Vector UL);
void MUSCL4(ScalarField const U,  double & UR, double & UL);
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
    cout<<"mesh is written in \"mesh.dat\""<<endl;
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
//LTS=LocalTimeStepping,当地时间步法,返回全局的时间步
double LTS()
{
    double min=1e10;
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

            if(dtLocal<min)
                min=dtLocal;
        }
  
    return min;
}

/**********************init &BC********************************/
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

/**********************solve********************************/
//三阶显式RungeKutta法
void solve()
{
    dt=LTS();//定常问题不需要dt, 但是荣格库塔法需要

    const double alpha[3]={0.1918, 0.4929, 1.0};
    //先定义Q0,用于保存原始的Q
    Field Q0;
    forEach(
            Q0[i][j][k] = Q[i][j][k];
    );
    
    //后面每一步都先计算残差, 后根据RK公式更新W
    
    for(unsigned a=0;a<=2;a++)   //a代表荣格库塔法的每一步
        for (I = 1; I <= maxI-2; I++)//I,J为单元编号, 只在此处变动!!
            for(J = 1; J <= maxJ-2; J++)
            {  
                cellGeometry();
                aeroConvert();
                //利用Roe格式计算通量
                roeFlux1();
                roeFlux4();
                //计算残差
                for (unsigned k=0; k<3; k++)
                {   
                    //R代表离开单元格的通量的矢量和, =右侧+左侧+上侧+下侧
                    //Fc1,Fc2, 由于与坐标系正方向反向, 故取负号
                    //右侧与上侧通量分别由临近单元格左侧与下侧取负号得来
                    R[I][J][k]=-Fc1[I+1][J][k] - Fc2[I][J][k] + Fc3[I][J+1][k] + Fc4[I][J][k];
                     
                    //利用荣格库塔法计更新流场
                    Q[I][J][k] = Q0[I][J][k] - alpha[a] * dt / volume[I][J] * R[I][J][k];
                }
            }
}


/**********************Roe********************************/
//Roe格式计算对流通量
/***********************flux1****************************/
/***********************flux1****************************/
void roeFlux1()
{   
    //利用MUSCL分裂流场变量, 并且求出Roe平均量 
    double rhoR, rhoL, pR, pL, uR, uL, vR, vL, HR, HL;
    double Vcv1L, Vcv1R;

    MUSCL1(rho, rhoR, rhoL);    //带限制器的三点MUSCL插值
    MUSCL1(p,   pR,   pL  );    //注意不要越界~
    MUSCL1(u,   uR,   uL  );
    MUSCL1(v,   vR,   vL  );
    MUSCL1(H,   HR,   HL  );
    

    Vcv1[I][J] =N1[I][J].x * u[I][J] +N1[I][J].y * v[I][J];
    MUSCL1(Vcv1, Vcv1L, Vcv1R);

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
    
    const double Vcv1_=N1[I][J].x * u_ +N1[I][J].y * v_;
    

    //计算lambda
    lambda1=Harten(fabs(Vcv1_-c_));
    lambda2=Harten(fabs(Vcv1_   ));
    lambda3=Harten(fabs(Vcv1_+c_));   

    //求出Roe矩阵相关值
        //delu,delv,delw代表三个分量, delVcv代表大写Vcv的delta
        //lambda123分别是是Vcv-c, Vcv, Vcv+c
    const double delP   =pR-pL;
    const double delRho =rhoR-rhoL;
    const double delu   =uR-uL;
    const double delVcv1=Vcv1R-Vcv1L;

    Vector delF1, delF234, delF5;
    delF1[0]  = lambda1 * (delP-rho_*c_*delVcv1)/(2*c_*c_) * 1;
    delF1[1]  = lambda1 * (delP-rho_*c_*delVcv1)/(2*c_*c_) * (u_-c_*N1[I][J].x);
    delF1[2]  = lambda1 * (delP-rho_*c_*delVcv1)/(2*c_*c_) * (u_-c_*N1[I][J].y);
	delF1[3]  = lambda1 * (delP-rho_*c_*delVcv1)/(2*c_*c_) * (H_-c_*Vcv1_);
  
    delF234[0]= lambda2 * ( (delRho-delP/(c_*c_))*1  + 0); 
    delF234[1]= lambda2 * ( (delRho-delP/(c_*c_))*u_ 
                   + rho_ *(delu   -  delVcv1)    );
    delF234[2]= lambda2 * ( (delRho-delP/(c_*c_))*u_ 
                   + rho_ *(delu   -  delVcv1)    );                
    delF234[3]= lambda2 * ( (delRho-delP/(c_*c_))*u_
                              +rho_ *(u_*delu - Vcv1_*delVcv1) );
     
    delF5[0]  = lambda3 * (delP+rho_*c_*delVcv1)/(2*c_*c_) *1;
    delF5[1]  = lambda3 * (delP+rho_*c_*delVcv1)/(2*c_*c_) *(u_+c_);
    delF5[2]  = lambda3 * (delP+rho_*c_*delVcv1)/(2*c_*c_) *(u_+c_);    
    delF5[3]  = lambda3 * (delP+rho_*c_*delVcv1)/(2*c_*c_) *(H_+c_*Vcv1_);

    //分裂Fc
    Vector QL, QR, Fc1R, Fc1L;
    MUSCL1(Q,  QR,   QL  );
    toFlux(QR, Fc1R);
    toFlux(QL, Fc1L);

    //最后算出单元IJ的对流通量Fc1 与Fc3
    for(unsigned k=0; k<3; k++)
    {
        Fc1[I][J][k]=0.5*( Fc1R[k] - Fc1L[k] - (delF1[k]+delF234[k]+delF5[k]) );
        Fc3[I][J][k]=-Fc1[I][J+1][k];
    }
}

//N1方向的(下), 带限制器的三点MUSCL插值,第一个参数是输入,后两个输出, 系数k^为1/3
void MUSCL1(Field const U, Vector UR, Vector UL)
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

    for(unsigned k=0; k<3; k++)
    {
        UR[k] = U[I][J + 1][k] - 0.5 * deltaR;
        UL[k] = U[I][J    ][k] + 0.5 * deltaL;
    }
}

//N1方向(下), 重载用于标量的带限制器的三点MUSCL插值函数 ,第一个参数是输入,后两个输出, 系数k^为1/3
void MUSCL1(ScalarField const U, double & UR, double & UL)
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


/***********************flux4****************************/
/***********************flux4****************************/
void roeFlux4()
{   
    //利用MUSCL分裂流场变量, 并且求出Roe平均量 
    double rhoR, rhoL, pR, pL, uR, uL, vR, vL, HR, HL;
    double Vcv4L, Vcv4R;
    Vcv4[I][J] =N4[I][J].x * u[I][J] +N4[I][J].y * v[I][J];

    MUSCL4(rho, rhoR, rhoL);    //带限制器的三点MUSCL插值
    MUSCL4(p,   pR,   pL  );    //注意不要越界~
    MUSCL4(u,   uR,   uL  );
    MUSCL4(v,   vR,   vL  );
    MUSCL4(H,   HR,   HL  );
    MUSCL4(Vcv4, Vcv4L, Vcv4R);

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

    const double Vcv4_=N4[I][J].x * u_ +N4[I][J].y * v_;

   //计算lambda
    lambda1=Harten(fabs(Vcv4_-c_));
    lambda2=Harten(fabs(Vcv4_   ));
    lambda3=Harten(fabs(Vcv4_+c_));   

    //求出Roe矩阵相关值
        //delu,delv,delw代表三个分量, delVcv代表大写Vcv的delta
        //lambda123分别是是Vcv-c, Vcv, Vcv+c
    const double delP   =pR-pL;
    const double delRho =rhoR-rhoL;
    const double delu   =uR-uL;
    const double delVcv4=Vcv4R-Vcv4L;


    Vector delF1, delF234, delF5;
    delF1[0]  = lambda1 * (delP-rho_*c_*delVcv4)/(2*c_*c_) * 1;
    delF1[1]  = lambda1 * (delP-rho_*c_*delVcv4)/(2*c_*c_) * (u_-c_*N1[I][J].x);
    delF1[2]  = lambda1 * (delP-rho_*c_*delVcv4)/(2*c_*c_) * (u_-c_*N1[I][J].y);
	delF1[3]  = lambda1 * (delP-rho_*c_*delVcv4)/(2*c_*c_) * (H_-c_*Vcv4_);
  
    delF234[0]= lambda2 * ( (delRho-delP/(c_*c_))*1  + 0); 
    delF234[1]= lambda2 * ( (delRho-delP/(c_*c_))*u_ 
                   + rho_ *(delu   -  delVcv4)    );
    delF234[2]= lambda2 * ( (delRho-delP/(c_*c_))*u_ 
                   + rho_ *(delu   -  delVcv4)    );                
    delF234[3]= lambda2 * ( (delRho-delP/(c_*c_))*u_
                              +rho_ *(u_*delu - Vcv4_*delVcv4) );
     
    delF5[0]  = lambda3 * (delP+rho_*c_*delVcv4)/(2*c_*c_) *1;
    delF5[1]  = lambda3 * (delP+rho_*c_*delVcv4)/(2*c_*c_) *(u_+c_);
    delF5[2]  = lambda3 * (delP+rho_*c_*delVcv4)/(2*c_*c_) *(u_+c_);    
    delF5[3]  = lambda3 * (delP+rho_*c_*delVcv4)/(2*c_*c_) *(H_+c_*Vcv4_);

    //分裂Fc
    Vector QL, QR, Fc4R, Fc4L;
    MUSCL1(Q,  QR,   QL  );
    toFlux(QR, Fc4R);
    toFlux(QL, Fc4L);

    //最后算出单元IJ的对流通量Fc4 与Fc2
    for(unsigned k=0; k<3; k++)
    {
        Fc1[I][J][k]=0.5*( Fc4R[k] - Fc4L[k] - (delF1[k]+delF234[k]+delF5[k]) );
        Fc3[I][J][k]=-Fc1[I][J+1][k];
    }
}

//N4方向(左), 重载用于标量的带限制器的MUSCL插值函数 ,前两个参数是输入,后两个输出, 系数k^为1/3
void MUSCL4(Field const U, Vector UR, Vector UL)
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
    for(unsigned k=0;k<3;k++)
    {
        UR[k]=U[I+1][J][k]-0.5*deltaR;
        UL[k]=U[I  ][J][k]+0.5*deltaL;
    }
}

//N4方向(左), 重载用于标量的带限制器的三点MUSCL插值函数 ,第一个参数是输入,后两个输出, 系数k^为1/3
void MUSCL4(ScalarField const U, double & UR, double & UL)
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
    double c=sqrt(GAMMA*p[I][J]/rho[I][J]);
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
void aeroConvert()
{
    rho[I][J]=Q[I][J][0];
    u[I][J]=Q[I][J][1]/rho[I][J];
    v[I][J]=Q[I][J][1]/rho[I][J];
    p[I][J]=(GAMMA-1)* (Q[I][J][2] - rho[I][J]* (SQ(u[I][J])+SQ(v[I][J])) *0.5);
    H[I][J]=Q[I][J][3]/rho[I][J]+p[I][J]/rho[I][J];
}

void print()
{
    FILE *fp;
    fopen("pressure.dat", "w");
    for(unsigned i=0; i<=maxI; i++)
        for(unsigned j=0; j<=maxJ; j++)
        {
            fprintf(fp, "%.5f %.5f %.5f\n", mesh[i][j].x, mesh[i][j].y, p[I][J] );
        } 
    fopen("Ma.dat", "w");
    for(unsigned i=0; i<=maxI; i++)
        for(unsigned j=0; j<=maxJ; j++)
        {
            double Ma=sqrt(SQ(u[i][j])+SQ(v[i][j]))/sqrt(GAMMA*p[i][j]/rho[i][j]);
            fprintf(fp, "%.5f %.5f %.5f\n", mesh[i][j].x, mesh[i][j].y, Ma );
        } 
    fclose(fp);
}       

/***************************main  **************************/
int main()
{  
	cout<<"Case 1 for inlet 1.8Ma; 2 for 1.5Ma\n\n Current Case: "<<caseNo<<endl;
    genMesh();
	
    for (int step=0;  step<=STOP_STEP; step++)
    {   
        solve();                  //求解
    }

    print();
    return 0;
}
