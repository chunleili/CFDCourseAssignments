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
    for(unsigned i=cellBegin; i<=cellIEnd; i++)\
        for(unsigned j=cellBegin; j<=cellJEnd; j++) \
            {\
                codes\
            }\
}
//用以循环所有实际的单元格

#define check(val) if(val<0) cout<<"\n*Check: "<<#val" = "<<val<<" is negative!\n"

#define IJcheck(val,i,j) if(val<0) cout<<"\n&IJCheck: "<<#val\
<<" ("<<i<<","<<j<<") = "<<val<<" is negative!\n"

#define SQ(a) ((a)*(a))

#define caseNo (1)  //case1 1.8Ma, case2 1.5Ma
/***************************define the consts ********************************/
const int maxI=250, maxJ=50;
const int block1=(int)(maxI*0.2+0.1), block2=(int)(maxI*0.2+0.1);//注意随着maxI更改
const int cellBegin=1, cellIEnd=maxI, cellJEnd=maxJ; //0是左下虚网格,实际网格从下标1开始,到下标maxI/J结束
const int STOP_STEP=100;

const double RESIDUAL_LIMIT=1e-3;
const double GAMMA=1.4;
const double Rg=287;
const double Cp=1004.5;//1.4/0.4*287
const double Cv=717.5;//1004.5-287
const double CFL=0.7;
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

//在最左和最下侧分别铺设一层虚网格,最上和最右侧分别铺设两层虚网格
//0代表左下虚网格, maxI/J+1和maxI/J+2代表右上虚网格, 实际网格1~max, 共max个
typedef XY     MeshPoint[maxI+3][maxJ+3];        //用于存储网格点坐标
typedef double Field[maxI+3][maxJ+3][4];         //场,用于定义Q,Fc1等对象
typedef double ScalarField[maxI+3][maxJ+3];      //标量场,用于p,rho大小等场对象
typedef XY     VectorField[maxI+3][maxJ+3];      //向量场,用于定义面的单位法量 
typedef double Vector[4];                        //表示某一单元格的参数
typedef unsigned const Index;                    //用于传递编号,只读

/***************************declare the utility funcs  **************************/
void aeroConvert(Index i, Index j);
double safeSqrt(double xx);
void toFlux3();
void toFlux2();
void print();
/*********************declear the global variable &funcs *****************/
MeshPoint mesh;
ScalarField  volume, S1,S2,S3,S4; //面积, 逆时针顺序, 依次为下右上左
VectorField  N1,N2,N3,N4;       //面法向单位矢量
void genMesh();
void printMesh();
void cellGeometry(); 

Field Q, Fc1,Fc2,Fc3,Fc4;
Field R; //R for Residual, 残差
//对流通量, 逆时针顺序依次为下右上左, 向外为正, 只储存其大小, 方向由N1~N4给定
void init1();
void solve();
double LTS();

unsigned I, J, step;
ScalarField rho, u, v,  p, H, Vcv3, Vcv2; 
//Vcv代表contravirant velocity, 实际上就是面法向速度 
//double lambda1,lambda2,lambda3;

void roeFlux3();//roe格式求解通量,注意通量要比单元格数量多1
void roeFlux2();
void MUSCL3(ScalarField const U,  double & UR, double & UL);
void MUSCL2(ScalarField const U,  double & UR, double & UL);
double Harten(double lambda);

double residualRho, residualU, residualV, residualE;
/********************************Mesh**************************************/
//虚网格的几何参数和临近实际网格一致, 但是不需要生成网格点, 为了和单元格编号一一对应,也从1开始
//左下点代表的编号和本单元格编号一致
void genMesh()//生成网格点,注意点要比单元格数量分别多一层
{
    const double dx=1.0/maxI;
    double dy=0.3/maxJ;//后面会变
    for(int i=cellBegin; i<=block1; i++)
    {
        for(int j=cellBegin; j<=maxJ+1; j++)
        {
            mesh[i][j].x=(i-1)*dx;//x=(i-1)*dx
            mesh[i][j].y=(j-1)*dy;//y=(j-1)*dy下同
        }
    }
    
    for(int i=block1+1; i<=block1+block2; i++)
    {
        double h=0.25*(i-1)*dx-0.05;
        dy=(0.3-h)/maxJ;
        for(int j=cellBegin;j<=maxJ+1; j++)
        {
            mesh[i][j].x=(i-1)*dx;
            mesh[i][j].y=h+(j-1)*dy;
        }
    }

    dy=(0.3-0.05)/maxJ;
    for(int i=block1+block2+1; i<=maxI+1; i++)
    {
        for(int j=cellBegin;j<=maxJ+1; j++)
        {
            mesh[i][j].x=(i-1)*dx;
            mesh[i][j].y=0.05+(j-1)*dy;
        }
    }
}

void printMesh()
{
    ofstream fout("mesh.dat");
    fout
	<<"Title=\"Mesh\""<<endl
	<<"Variables=\"x\",\"y\""<<endl
	<<"Zone i="<<maxI+1<<", j="<<maxJ+1<<", f=point"<<endl;//注意点要比网格数多一个!!
    for(int j=cellBegin; j<=maxJ+1; j++)
    {
        for(int i=cellBegin; i<=maxI+1; i++)
        {
            fout<< mesh[i][j].x<<" "<<mesh[i][j].y<<endl;
        }
    }
    cout<<"\nmesh is written in \"mesh.dat\""<<endl;

    FILE* fpG;
    fpG=fopen("cellGeometry.txt", "w");
    fprintf(fpG,"i   j   x    y    volume   N2.x   N2.y    N3.x    N3.y   S1   S2   S3   S4\n");
    forAll( 
            fprintf(fpG, "%-3d  %-3d  %.3f %.3f %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e \n",
            i, j, mesh[i][j].x, mesh[i][j].y,
            volume[i][j], N2[i][j].x, N2[i][j].y, N3[i][j].x, N3[i][j].y,
            S1[i][j], S2[i][j], S3[i][j], S4[i][j] );
    );
    fclose(fpG);
}

void cellGeometry()
{
    double x1,x2,x3,x4, y1,y2,y3,y4;
    //从左下开始逆时针编号,左下点代表1,右下2,右上3,左上4
    //左下点代表本单元格坐标,此处ij代表点的编号
    for(unsigned i=cellBegin; i<=cellIEnd; i++)//注意范围
        for(unsigned j=cellBegin; j<=cellJEnd; j++) 
        {
            x1 = mesh[i][j].x;
            x2 = mesh[i + 1][j].x;
            x3 = mesh[i + 1][j + 1].x;
            x4 = mesh[i][j + 1].x;

            y1 = mesh[i][j].y;
            y2 = mesh[i + 1][j].y;
            y3 = mesh[i + 1][j + 1].y;
            y4 = mesh[i][j + 1].y;

            S1[i][j] = sqrt(SQ(x1 - x2) + SQ(y1 - y2)); //下侧面积S1
            S2[i][j] = sqrt(SQ(x2 - x3) + SQ(y2 - y3)); //右侧面积S2
            S3[i][j] = sqrt(SQ(x3 - x4) + SQ(y3 - y4)); //上侧面积S3
            S4[i][j] = sqrt(SQ(x4 - x1) + SQ(y4 - y1)); //左侧面积S4

            //面法向单位矢量, 实际上N1.x=sin(theta), N1.y=cos(theta), 再带上方向
            N1[i][j].x = (y2 - y1) / S1[i][j];
            N1[i][j].y = (x1 - x2) / S1[i][j];

            N2[i][j].x = (y3 - y2) / S2[i][j];
            N2[i][j].y = (x2 - x3) / S2[i][j];

            N3[i][j].x = (y4 - y3) / S3[i][j];
            N3[i][j].y = (x3 - x4) / S3[i][j];

            N4[i][j].x = (y1 - y4) / S4[i][j];
            N4[i][j].y = (x4 - x1) / S4[i][j];
            
            //计算单元体体积, 厚度取为1
            volume[i][j] = 0.5 * ((x1 - x3) * (y2 - y4) + (x4 - x2) * (y1 - y3));
            
        }    

        //虚网格也需要计算几何参数,直接赋予相邻网格的参数即可
}

/**********************LTS********************************/
//LTS=LocalTimeStepping,当地时间步法,返回当地时间步
double LTS()
{
    double dtLocal;
    double lambdaI, lambdaJ, SI, SJ;
    XY NI, NJ;
    double c=sqrt(GAMMA*p[I][J]/rho[I][J]);
    NI.x = 0.5 * (N3[I][J].x - N1[I][J].x);
    NI.y = 0.5 * (N3[I][J].y - N1[I][J].y);

    NJ.x = 0.5 * (N2[I][J].x - N4[I][J].x);
    NJ.y = 0.5 * (N2[I][J].y - N4[I][J].y);

    SI = 0.5 * (S1[I][J] + S3[I][J]);
    SJ = 0.5 * (S2[I][J] + S4[I][J]);

    lambdaI = (u[I][J] * NI.x + v[I][J] * NI.y)+ c;
    lambdaJ = (u[I][J] * NJ.x + v[I][J] * NJ.y)+ c;

    dtLocal = CFL * volume[I][J] / (lambdaI * SI + lambdaJ * SJ);

    return dtLocal;
}

/**********************init & BC********************************/
void init1()
{
    //初始全部给1.8Ma, 攻角为0, 压力给大气压101325, 静温300K, 其余推出如下:
    //声速c=sqrt(1.4*287*300)=347.1887
    //速度u=624.9397, v=0, 密度rho=p/RT=1.176829
    //单位体积总能E=p/(GAMMA-1)+0.5*rho*(u^2+v^2)
    for(unsigned i=cellBegin-1; i<=cellIEnd+2; i++)
        for(unsigned j=cellBegin-1; j<=cellJEnd+2; j++)
        {
            p[i][j]=101325;
            u[i][j]=624.9397;
            v[i][j]=0;
            rho[i][j]=1.176829;
            H[i][j]=496662.5;//Cp*300+0.5*625*625

            Q[i][j][0] = 1.176829; //101325/(287*300)
            Q[i][j][1] = 735.4473; //1.1768*625
            Q[i][j][2] = 0;
            Q[i][j][3] = 483117.6; // 101325/0.4+0.5*1.176829*625*625;
        }
}


//在最左和最下侧分别铺设一层虚网格,最上和最右侧分别铺设两层虚网格
//0代表左下虚网格, maxI/J+1和maxI/J+2代表右上虚网格, 实际网格1~max, 共max个
//上下是壁面, 壁面法向速度为0,壁面切向不变 壁面无穿透边界
void BCup()
{
    double nx,ny,p2,p3,pw;
    const unsigned j=cellJEnd;
    for(unsigned i=cellBegin; i<=cellJEnd; i++)
    {
        nx=N3[i][j].x;
        ny=N3[i][j].y;

        p2=p[i][j];
        p3=p[i][j-1];
        pw=0.5*(3*p2-p3);//壁面的压力用两点外推

        Fc3[i][j][0]=0;       
        Fc3[i][j][1]=pw*nx;
        Fc3[i][j][2]=pw*ny;
        Fc3[i][j][3]=0;
        //虚网格的值靠外推
        for(unsigned k=0; k<=3; k++)
        {
            Q[i][j+1][k]=2*Q[i][j][k]-  Q[i][j-1][k];
            Q[i][j+2][k]=3*Q[i][j][k]-2*Q[i][j-1][k];
        }
        aeroConvert(i,j+2);
        aeroConvert(i,j+1);

        Vcv2[i][j+1] = 2*Vcv2[i][j]-Vcv2[i-1][j];
        Vcv3[i][j+1] = 2*Vcv3[i][j]-Vcv3[i-1][j];

        Vcv2[i][j+2] = 3*Vcv2[i][j]-2*Vcv2[i][j-1];
        Vcv3[i][j+2] = 3*Vcv3[i][j]-2*Vcv3[i][j-1];
    }
}

//下壁面,只有一层虚网格
void BCdown()
{
    double nx, ny, p2, p3, pw;
    unsigned j= cellBegin;
    for (unsigned i = cellBegin; i <= cellJEnd; i++)
    {
        nx = N1[i][j].x;
        ny = N1[i][j].y;

        p2 = p[i][j];
        p3 = p[i][j + 1];
        pw = 0.5 * (3 * p2 - p3); //壁面的压力用两点外推

        Fc1[i][j][0] = 0;
        Fc1[i][j][1] = pw * nx;
        Fc1[i][j][2] = pw * ny;
        Fc1[i][j][3] = 0;

        //虚网格的值靠外推
        for (unsigned k = 0; k <= 3; k++)
        {
            Q[i][j - 1][k] = 2 * Q[i][j][k] - Q[i][j + 1][k];
        }
        aeroConvert(i, j - 1);
        Vcv2[i][j-1] = 2*Vcv2[i][j]-Vcv2[i+1][j];
        Vcv3[i][j-1] = 2*Vcv3[i][j]-Vcv3[i+1][j];

        if (rho[i][j - 1] < 0)
        {
            cout << "\nhere! rho<0 i=" << i << " j=" << j - 1 << " rho=" << rho[i][j - 1] << endl;
            printf("\nQ[i][0]=%f  %f  %f  %f\n", Q[i][0][0], Q[i][0][1], Q[i][0][2], Q[i][0][3]);
            printf("\nQ[i][1]=%f  %f  %f  %f\n", Q[i][1][0], Q[i][1][1], Q[i][1][2], Q[i][1][3]);
            printf("\nQ[i][2]=%f  %f  %f  %f\n", Q[i][2][0], Q[i][2][1], Q[i][2][2], Q[i][2][3]);
            exit(0);
        }
    }
}

//入口维持1.8Ma, 出口用0梯度外推出来;
void BCright()
{
    //出口虚网格
    for(unsigned j=cellBegin; j<=cellJEnd; j++)
    {
        for (unsigned k = 0; k < 4; k++)
        {
            Q[cellIEnd+2][j][k] = Q[cellIEnd+1][j][k] = Q[cellIEnd][j][k];
        }
        aeroConvert(cellIEnd+2,j);
        aeroConvert(cellIEnd+1,j);
        Vcv2[cellIEnd+2][j] = Vcv2[cellIEnd+1][j]=Vcv2[cellIEnd][j];
        Vcv3[cellIEnd+2][j] = Vcv3[cellIEnd+1][j]=Vcv3[cellIEnd][j];
    }
}

void BCleft()
{
    unsigned i=cellBegin-1;
    for(unsigned j=cellBegin; j<=cellJEnd; j++)
    {
        p[i][j] = 101325;
        u[i][j] = 624.9397;
        v[i][j] = 0;
        rho[i][j] = 1.176829;
        H[i][j] = 496662.5; //Cp*300+0.5*625*625

        Q[i][j][0] = 1.176829; //101325/(287*300)
        Q[i][j][1] = 735.4473; //1.1768*625
        Q[i][j][2] = 0;
        Q[i][j][3] = 483117.6; // 101325/0.4+0.5*1.176829*625*625;
    }
}


/**********************solve********************************/
//三阶显式RungeKutta法
void solve()
{
    //const double alpha[3]={0.1481, 0.4, 1.0};
    //因为含有激波,所以即使是二阶迎风也要用一阶的系数

    //先定义Q0,用于保存原始的Q
    //Field Q0;
    //forAll(
    //       for (unsigned k = 0; k <= 3; k++)
    //            Q0[i][j][k] = Q[i][j][k];
    //);

    //后面每一步都先计算残差, 后根据RK公式更新W
    double dt;
    double rRho,rU,rV,rE;
    //for(unsigned a=0;a<=2;a++)   //a代表荣格库塔法的每一步
    //{
        BCleft();
        BCdown();
        BCup();
        BCright();
        for (I = cellBegin; I <= cellIEnd; I++)//I,J为单元编号, 只在此处变动!!
        {
            for(J = cellBegin; J <= cellJEnd; J++)
            {   
                aeroConvert(I,J);
                Vcv3[I][J] = N3[I][J].x * u[I][J] + N3[I][J].y * v[I][J];
                Vcv2[I][J] = N2[I][J].x * u[I][J] + N2[I][J].y * v[I][J];

                //利用Roe格式计算通量
                roeFlux3();
                roeFlux2();

                dt=LTS(); //当地时间步法求该单元格的时间步             

                //计算残差
                for (unsigned k=0; k<4; k++)
                {   
                    //R代表离开单元格的通量的矢量和, =右侧+左侧+上侧+下侧
                    //左侧与下侧通量分别由临近单元格右侧与上侧取负号得来
                    Fc1[I][J][k]=-Fc3[I][J-1][k];
                    Fc4[I][J][k]=-Fc2[I-1][J][k];

                    R[I][J][k]=Fc1[I][J][k]*S1[I][J] + Fc2[I][J][k]*S2[I][J] 
                                + Fc3[I][J][k]*S3[I][J] + Fc4[I][J][k]*S4[I][J];
                     
                    //利用荣格库塔法计更新流场
                    //Q[I][J][k] = Q0[I][J][k] - alpha[a] * dt / volume[I][J] * R[I][J][k];
                    Q[I][J][k] = Q[I][J][k] -  dt / volume[I][J] * R[I][J][k];
                }
                
                //记录残差
                rRho=fabs(R[I][J][0]);
                rU  =fabs(R[I][J][1]);
                rV  =fabs(R[I][J][2]);
                rE  =fabs(R[I][J][3]);

                if (residualRho < rRho)
                    residualRho = rRho;
                if (residualU < rU)
                    residualU = rU;
                if (residualV < rV)
                    residualV = rV;                
                if (residualE < rE)
                    residualE = rE;
            }
        }
        //}
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
/*
    MUSCL3(rho,  rhoR, rhoL);    //带限制器的三点MUSCL插值
    MUSCL3(p,    pR,   pL  );    //注意不要越界~
    MUSCL3(u,    uR,   uL  );
    MUSCL3(v,    vR,   vL  );
    MUSCL3(H,    HR,   HL  );
    MUSCL3(Vcv3, VcvL, VcvR);   //为了方便代码重用, 分裂Vcv3后记为VcvL和VcvR
*/
    rhoR=rho[I][J];     rhoL=rho[I][J-1];
    pR  =p  [I][J];     pL  =p  [I][J-1];
    uR  =u  [I][J];     uL  =u  [I][J-1];
    vR  =v  [I][J];     vL  =v  [I][J-1];
    HR  =H  [I][J];     HL  =H  [I][J-1];

    VcvR=uR[I][J]*nx+vR[I][J]*ny;

    N3[I][0].x=N3[I][1].x;
    N3[I][0].y=N3[I][1].y;
    VcvL=uL[I][J-1]*N3[I][J-1].x+vL*N3[I][J-1].y;//虚网格的nx ny与本网格相同

    //计算Roe平均量
    const double denoLR=sqrt(rhoL)+ sqrt(rhoR);
    const double L=sqrt(rhoL)/ denoLR;//定义两个系数
    const double R=sqrt(rhoR)/ denoLR;
    
    const double rho_=sqrt(rhoL*rhoR);
    const double u_  =uL*L+uR*R;
    const double v_  =vL*L+vR*R;
    const double H_  =HL*L+HR*R;
    const double q_2 =u_*u_+v_*v_;
    const double c_  =sqrt((GAMMA-1)*(H_-q_2/2));

    const double Vcv_=nx * u_ + ny * v_;
    check(rhoL*1.0);
    check(rhoR*1.0);
    check((GAMMA-1)*(H_-q_2/2)*1.0);
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
    Vector  FR, FL;
    FR[0] = rhoR*VcvR;
    FR[1] = rhoR*VcvR*uR + nx*pR;
    FR[2] = rhoR*VcvR*vR + ny*pR;
    FR[3] = rhoR*VcvR*HR;

    FL[0] = rhoL*VcvL;
    FL[1] = rhoL*VcvL*uL + nx*pL;
    FL[2] = rhoL*VcvL*vL + ny*pL;
    FL[3] = rhoL*VcvL*HL;

    //最后算出单元IJ的对流通量Fc3
    for(unsigned k=0; k<4; k++)
    {
        Fc3[I][J][k]=0.5*( FR[k] + FL[k] - (delF1[k]+delF234[k]+delF5[k]) );
    }
}


//N3方向, 用于标量的带限制器的三点MUSCL插值函数 ,第一个参数是输入,后两个输出, 系数k^为1/3
void MUSCL3(ScalarField const U, double & UR, double & UL)
{
    const double epsilon= sqrt(volume[I][J]); //限制器参数epsilon与几何尺寸相关

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

/*
    MUSCL2(rho,  rhoR, rhoL);    //带限制器的三点MUSCL插值
    MUSCL2(p,    pR,   pL  );    //注意不要越界~
    MUSCL2(u,    uR,   uL  );
    MUSCL2(v,    vR,   vL  );
    MUSCL2(H,    HR,   HL  );
    MUSCL2(Vcv2, VcvL, VcvR);   //为了方便代码重用, 分裂Vcv后记为VcvL和VcvR
*/
    rhoR=rho[I  ][J];
    rhoL=rho[I-1][J];
    pR  =p  [I  ][J];
    pL  =p  [I-1][J];
    uR  =u  [I  ][J];
    uL  =u  [I-1][J];
    vR  =v  [I  ][J];
    vL  =v  [I-1][J];
    HR  =H  [I  ][J];
    HL  =H  [I-1][J];

    VcvR=uR[I][J]*nx+vR[I][J]*ny;
    N3[I][0].x=N3[I][1].x;
    N3[I][0].y=N3[I][1].y;
    VcvL=uL[I][J-1]*N3[I][J-1].x+vL*N3[I][J-1].y;//虚网格的nx ny与本网格相同

    //计算Roe平均量
    const double denoLR=sqrt(rhoL)+ sqrt(rhoR);
    const double L=sqrt(rhoL)/ denoLR;//定义两个系数
    const double R=sqrt(rhoR)/ denoLR;
    
    const double rho_=sqrt(rhoL*rhoR);
    const double u_  =uL*L+uR*R;
    const double v_  =vL*L+vR*R;
    const double H_  =HL*L+HR*R;
    const double q_2 =u_*u_+v_*v_;
    const double c_  =sqrt((GAMMA-1)*(H_-q_2/2));

    const double Vcv_=nx * u_ + ny * v_;

    IJcheck(rhoL, I, J);
    IJcheck(rhoR, I, J);
    IJcheck((GAMMA-1)*(H_-q_2/2), I, J);

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
    Vector  FR, FL;

    FR[0] = rhoR*VcvR;
    FR[1] = rhoR*VcvR*uR + nx*pR;
    FR[2] = rhoR*VcvR*vR + ny*pR;
    FR[3] = rhoR*VcvR*HR;

    FL[0] = rhoL*VcvL;
    FL[1] = rhoL*VcvL*uL + nx*pL;
    FL[2] = rhoL*VcvL*vL + ny*pL;
    FL[3] = rhoL*VcvL*HL;

    //最后算出单元IJ的对流通量Fc2
    for(unsigned k=0; k<4; k++)
    {
        Fc2[I][J][k]=0.5*( FR[k] + FL[k] - (delF1[k]+delF234[k]+delF5[k]) );
    }
}

//N2方向, 用于标量的带限制器的MUSCL插值函数 ,前两个参数是输入,后两个输出, 系数k^为1/3
void MUSCL2(ScalarField const U, double & UR, double & UL)
{
    const double epsilon= sqrt(volume[I][J]); //限制器参数epsilon与几何尺寸相关
    
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
    IJcheck(GAMMA*p[I][J]/rho[I][J], I, J);

    double delta=0.1*c; //熵修正, Harten's entropy correction
    if(fabs(lambda)<=delta)
         lambda=(lambda*lambda+ delta*delta) /(2*delta);
    return lambda;

}


/***************************utility  **************************/


//将Q转化为rho u v p H Vcv
void aeroConvert(Index i, Index j)
{
        rho[i][j] = Q[i][j][0];
        u[i][j] = Q[i][j][1] / rho[i][j];
        v[i][j] = Q[i][j][2] / rho[i][j];
        p[i][j] = (GAMMA - 1) * (Q[i][j][3] - rho[i][j] * (SQ(u[i][j]) + SQ(v[i][j])) * 0.5);
        H[i][j] = (Q[i][j][3] + p[i][j]) / rho[i][j];

        IJcheck(rho[i][j] * 1.0, i, j);
}

double safeSqrt(double xx)
{
    if(xx<0) 
    {
        cout<<"\n##Erro sqrt! Value "<<xx<<" is negative!\n";
        //getchar();
    }
    return sqrt(xx);
}

void print()
{
    FILE *fpField;
    double xx,yy,uu,vv,rrho,pp,TT,MMa;
    fpField=fopen("result.dat", "w");
    fprintf(fpField, "Title=\"Field\"\nVariables=\"x\",\"y\",\"rho\",\"u\",\"v\",\"p\",\"T\",\"Ma\"\nZone T=\"zone1\" i=%d,j=%d,f=point\n", maxI, maxJ);
    forAll(
        xx=mesh[i][j].x;
        yy=mesh[i][j].y;
        rrho=rho[i][j];
        pp=p[i][j];
        TT=pp/(rrho*287);
        uu=u[i][j];
        vv=v[i][j];
        MMa=sqrt( ( SQ(u[i][j])+SQ(v[i][j]) ) / (GAMMA*p[i][j]/rho[i][j]) );
        fprintf(fpField, "%.3f %.3f %.4e %.4e %.4e %.4e %.4e %.4e \n", xx, yy, rrho, uu,vv, pp, TT, MMa );
    ); 
    fclose(fpField);
}   



/***************************main  **************************/
int main()
{  
	cout<<"\nCase 1 for inlet 1.8Ma; 2 for 1.5Ma\n\nCurrent Case: "<<caseNo<<endl;
    genMesh();//生成网格
    cellGeometry();//计算各个单元格面积,方向矢量等
    printMesh();//打印网格

    //初始化流场
    init1();
    //else               init2(); 

    cout<<"\nInitialzation done.\n\n"; 

    FILE *fpR;
    fpR=fopen("residual.dat", "w");
    fprintf(fpR,"iter  continuity x-velocity y-velocity Energy\n");
    for (step=1;  step<=20; step++)
    {   
        cout<<"step= "<<step<<endl;
        //printf("\nbefore solve: Q[i][1]=%f  %f  %f  %f\n", Q[I][1][0], Q[I][1][1],Q[I][1][2],Q[I][1][3]);
        solve();                  //求解
        //printf("\nafter solve: Q[i][1]=%f  %f  %f  %f\n", Q[I][1][0], Q[I][1][1],Q[I][1][2],Q[I][1][3]);

        fprintf(fpR, "%-5d %.4e %.4e %.4e %.4e\n", step, residualRho, residualU, residualV, residualE);
        printf("%-5d %.4e %.4e %.4e %.4e\n", step, residualRho, residualU, residualV, residualE);
        print();
    }
    fclose(fpR);


    print();//打印结果
    return 0;
}
