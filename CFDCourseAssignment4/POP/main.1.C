/***************************include & namespace  **************************/
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<cstdio>
#include<fstream>
#include<iomanip>
using namespace std;

/***************************MACRO             **************************/
#define IJcheck(val,i,j) if(val<0) cout<<"\n&IJCheck: "<<#val\
<<" ("<<i<<","<<j<<") = "<<val<<" is negative!\n"

#define SQ(a) ((a)*(a))

#define caseNo (1)  //case1 1.8Ma, case2 1.5Ma
/***************************define the consts ********************************/
const int maxI=250, maxJ=50;
const int block1=(int)(maxI*0.2+0.1), block2=(int)(maxI*0.2+0.1);//注意随着maxI更改
const int cellBegin=1, cellIEnd=maxI, cellJEnd=maxJ; //0是左下虚网格,实际网格从下标1开始,到下标maxI/J结束

const double GAMMA=1.4;
const double Rg=287;
const double Cp=1004.5;//1.4/0.4*287
const double Cv=717.5;//1004.5-287
const double CFL=0.7;

const double p0=101325;
const double u0=624.9397;
const double v0=0;
const double T0=300;
const double rho0=1.176829268;
const double H0=496662.5;
const double E0=483117.6;
/***************************define the type ********************************/
typedef struct XY
{
    public:
    double x;
    double y;
}XY;

//在最左和最下侧分别铺设一层虚网格,最上和最右侧分别铺设两层虚网格
//0代表左下虚网格, maxI/j+1和maxI/j+2代表右上虚网格, 实际网格1~max, 共max个
typedef XY     MeshPoint[maxI+3][maxJ+3];        //用于存储网格点坐标
typedef double Field[maxI+3][maxJ+3][4];         //场,用于定义Q,Fc等对象
typedef double ScalarField[maxI+3][maxJ+3];      //标量场,用于p,rho大小等场对象
typedef XY     VectorField[maxI+3][maxJ+3];      //向量场,用于定义面的单位法量 
typedef double Vector[4];                        //表示某一单元格的参数
typedef unsigned const Index;                    //用于传递编号,只读

/***************************declare the utility funcs  **************************/
void aeroConvert(Index i, Index j);
void print();
void extrapolation(Index newI, Index newJ, Index oldI, Index oldJ);
/*********************declear the global variable &funcs *****************/
MeshPoint mesh;
ScalarField  volume, S1,S2,S3,S4; //面积, 逆时针顺序, 依次为下右上左
VectorField  N1,N2,N3,N4;       //面法向单位矢量
double nx,ny;

void genMesh();
void printMesh();
void cellGeometry(); 

Field Fc1,Fc2,Fc3,Fc4;//右侧和上侧的对流通量的大小, 向外为正, 方向由N1~N4给定
Field Q;//Q用于存储守恒变量
Field Residual; //单元格内的残差

void init1();
void iteration();
double LTS();

unsigned i, j,k, step;
ScalarField rho, u, v,  p, H;
double c; 

void roeFlux();//roe格式求解通量,注意通量要比单元格数量多1
Vector maxR;

FILE* fpDe;

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
}



/**********************LTS********************************/
//LTS=LocalTimeStepping,当地时间步法,返回当地时间步
double LTS()
{
    
    double dtLocal;
    double lambdaI, lambdaJ, SI, SJ;
    XY NI, NJ;
    NI.x = 0.5 * (N3[i][j].x - N1[i][j].x);
    NI.y = 0.5 * (N3[i][j].y - N1[i][j].y);

    NJ.x = 0.5 * (N2[i][j].x - N4[i][j].x);
    NJ.y = 0.5 * (N2[i][j].y - N4[i][j].y);

    SI = 0.5 * (S1[i][j] + S3[i][j]);
    SJ = 0.5 * (S2[i][j] + S4[i][j]);

    lambdaI = (u[i][j] * NI.x + v[i][j] * NI.y)+ c;
    lambdaJ = (u[i][j] * NJ.x + v[i][j] * NJ.y)+ c;

    dtLocal = CFL * volume[i][j] / (lambdaI * SI + lambdaJ * SJ);

    return dtLocal;
    

}

/**********************init & BC********************************/
void init1()
{
    //初始全部给1.8Ma, 攻角为0, 压力给大气压101325, 静温300K, 其余推出如下:
    //声速c=sqrt(1.4*287*300)=347.1887
    //速度u=624.9397, v=0, 密度rho=p/RT=1.176829
    //单位体积总能E=p/(GAMMA-1)+0.5*rho*(u^2+v^2)
    for(unsigned i=cellBegin-1; i<=cellIEnd+2; i++)//包括虚网格!
    {
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
}

//在最左和最下侧分别铺设一层虚网格,最上和最右侧分别铺设两层虚网格
//0代表左下虚网格, maxI/j+1和maxI/j+2代表右上虚网格, 实际网格1~max, 共max个
//上下是壁面, 壁面法向速度为0,壁面切向不变 壁面无穿透边界
void BC()
{
    double pw;
	//上边
    const unsigned j=cellJEnd;
    for(unsigned i=cellBegin; i<=cellIEnd; i++)
    {
        //pw=0.5*(3*p[i][j]-p[i][j-1]);//壁面的压力用两点外推
        pw=p[i][j];
        nx=N3[i][j].x; ny=N3[i][j].y;

        Fc1[i][j+1][0]=0;       
        Fc1[i][j+1][1]=pw*nx;
        Fc1[i][j+1][2]=pw*ny;
        Fc1[i][j+1][3]=0;
	}
	//下边
    for (unsigned i = cellBegin; i <= cellIEnd; i++)
    {
        nx=N1[i][1].x; ny=N1[i][1].y;
        //pw = 0.5 * (3 * p[i][j] - p[i][j+1]); //壁面的压力用两点外推
        pw=p[i][1];
        
        Fc1[i][1][0] = 0;     //注意,虚网格nx ny使用相邻网格的值
        Fc1[i][1][1] = pw*nx;
        Fc1[i][1][2] = pw*ny;
        Fc1[i][1][3] = 0;

      rho[i][0] = rho[i][1];
        u[i][0] =   u[i][1];
        v[i][0] =   v[i][1];
        p[i][0] =   p[i][1];
        H[i][0] =   H[i][1];
    }
	//左边
    for(unsigned j=cellBegin; j<=cellJEnd; j++)
    {
        //这些流动变量值是虚网格,仅仅用于插值的模板
        p[0][j] = p0;
        u[0][j] = u0;
        v[0][j] = v0;
      rho[0][j] = rho0;
        H[0][j] = H0; 

        //Fc4[1][j]才是真正的边界
        const double nx=1, ny=0;
        const double Vcv0=u0;
        Fc4[1][j][0]=rho0*Vcv0;
        Fc4[1][j][1]=rho0*u0*Vcv0+p0*nx;
        Fc4[1][j][2]=rho0*v0*Vcv0+p0*ny;
        Fc4[1][j][3]=rho0*H0*Vcv0;
    }
    //右边//虚网格给予流动变量,然后根据roe格式算出虚网格的Fc4[maxI+1][j],从而赋予需要的Fc2[maxI][j]
    for(unsigned j=cellBegin; j<=cellJEnd; j++)
    {
        extrapolation(maxI+1,j,maxI,j);
    }
      

}

/**********************iteration********************************/
void iteration()
{
    double dt;
    double T;
    for (i = cellBegin; i <= cellIEnd; i++) //i,J为单元编号, 只在此处变动!!
    {
        for (j = cellBegin; j <= cellJEnd; j++)
        {

            c=sqrt(GAMMA*p[i][j]/rho[i][j]);

            dt = LTS(); //当地时间步法求该单元格的时间步

            //计算残差
            for (unsigned k = 0; k < 4; k++)
            {
                //R代表离开单元格的通量的矢量和, =右侧+左侧+上侧+下侧
                //左侧与下侧通量分别由临近单元格右侧与上侧取负号得来


                Q[i][j][k] = Q[i][j][k] - dt / volume[i][j] * Residual[i][j][k];
            }
            aeroConvert(i,j);

            T=(H[i][j]-0.5*(SQ(u[i][j])+SQ(v[i][j])) )/Cp;

            IJcheck(T, i, j);
            IJcheck(rho[i][j] * 1.0, i,j);
            IJcheck(p[i][j]*1.0,i,j);
            if(T<0) exit(3);
            if(p[i][j]<0) exit(1);
            if(rho[i][j]<0) exit(2);

        }
    }
}

/**********************Roe********************************/
//Roe格式计算对流通量
Vector FR, FL;
Vector AARoe;
double rhoR, rhoL, pR, pL, uR, uL, vR, vL, HR, HL;
double VcvL, VcvR;

void calARoe()
{
    VcvR=uR*nx+vR*ny,VcvL=uL*nx+vL*ny;
    double delP = pR - pL;
    double delRho = rhoR - rhoL;
    double delu = uR - uL;
    double delv = vR - vL;
    double delVcv = VcvR - VcvL;

    double denoLR = sqrt(rhoL) + sqrt(rhoR);
    double Lcof = sqrt(rhoL) / denoLR; //定义两个系数
    double Rcof = sqrt(rhoR) / denoLR;

    double rho_ = sqrt(rhoL * rhoR);
    double u_ = uL * Lcof + uR * Rcof;
    double v_ = vL * Lcof + vR * Rcof;
    double H_ = HL * Lcof + HR * Rcof;
    double q_2 = u_ * u_ + v_ * v_;
    double c_ = sqrt((GAMMA - 1) * (H_ - q_2 / 2.0));

    double Vcv_ = nx * u_ + ny * v_;

    if (H_ - q_2 / 2 < 0)
    {
        printf("\nH_-q_2<0!!\n"); //exit(1);
    }

    //计算lambda
    double lambda1 = fabs(Vcv_ - c_);
    double lambda2 = fabs(Vcv_);
    double lambda3 = fabs(Vcv_ + c_);

    //熵修正 Harten's entropy correction
    double delta = 0.05 * c;
    if (fabs(lambda1) <= delta)
        lambda1 = (lambda1 * lambda1 + delta * delta) / (2 * delta);
    if (fabs(lambda2) <= delta)
        lambda2 = (lambda2 * lambda2 + delta * delta) / (2 * delta);
    if (fabs(lambda3) <= delta)
        lambda3 = (lambda3 * lambda3 + delta * delta) / (2 * delta);

    //求出Roe矩阵相关值
    //lambda123分别是是Vcv-c, Vcv, Vcv+c

    double delF1[4], delF234[4], delF5[4];
    const double coeff1 = lambda1 * (delP - rho_ * c_ * delVcv) / (2 * c_ * c_); //定义系数以减少计算量
    delF1[0] = coeff1 * 1;
    delF1[1] = coeff1 * (u_ - c_ * nx);
    delF1[2] = coeff1 * (v_ - c_ * ny);
    delF1[3] = coeff1 * (H_ - c_ * Vcv_);

    const double coeff2 = lambda2 * (delRho - delP / (c_ * c_));
    const double coeff3 = lambda2 * rho_;
    delF234[0] = coeff2;
    delF234[1] = coeff2 * u_ + coeff3 * (delu - delVcv * nx);
    delF234[2] = coeff2 * v_ + coeff3 * (delv - delVcv * ny);
    delF234[3] = coeff2 * q_2 / 2.0 + coeff3 * (u_ * delu + v_ * delv - Vcv_ * delVcv);

    const double coeff5 = lambda3 * (delP + rho_ * c_ * delVcv) / (2 * c_ * c_);
    delF5[0] = coeff5 * 1;
    delF5[1] = coeff5 * (u_ + c_ * nx);
    delF5[2] = coeff5 * (v_ + c_ * ny);
    delF5[3] = coeff5 * (H_ + c_ * Vcv_);

    for (unsigned k = 0; k <= 3; k++)
        AARoe[k] = (delF1[k] + delF234[k] + delF5[k]);

    //分裂Fc
    FR[0] = rhoR*VcvR;
    FR[1] = rhoR*VcvR*uR + nx*pR;
    FR[2] = rhoR*VcvR*vR + ny*pR;
    FR[3] = rhoR*VcvR*HR;

    FL[0] = rhoL*VcvL;
    FL[1] = rhoL*VcvL*uL + nx*pL;
    FL[2] = rhoL*VcvL*vL + ny*pL;
    FL[3] = rhoL*VcvL*HL;
}


void splitI()
{
	pL = p[i - 1][j], pR = p[i][j];
	rhoL = rho[i - 1][j], rhoR = rho[i][j];
	uL = u[i - 1][j], uR = u[i][j];
	vL = v[i - 1][j], vR = v[i][j];
	HL = H[i - 1][j], HR = H[i][j];
}

void splitJ()
{
	pL = p[i][j - 1], pR = p[i][j];
	rhoL = rho[i][j - 1], rhoR = rho[i][j];
	uL = u[i][j - 1], uR = u[i][j];
	vL = v[i][j - 1], vR = v[i][j];
	HL = H[i][j - 1], HR = H[i][j];
}

void roeFlux()
{
    for (unsigned k = 0; k <= 3; k++)
        maxR[k] = 0;

//xin
    for( i=1;i<=maxI+1;i++)   //注意范围,计算最右侧虚网格
        for( j=1;j<=maxJ;j++)
		{   
			c = sqrt(GAMMA*p[i][j]/rho[i][j]);
            nx = -N4[i][j].x, ny = -N4[i][j].y;
            splitI();
			calARoe();
			for (unsigned k = 0; k <= 3; k++)
				Fc4[i][j][k] = S4[i][j] * (FL[k] + FR[k] - AARoe[k]) / 2;
            nx = -N1[i][j].x, ny = -N1[i][j].y;
            splitJ();
			calARoe();
			for (unsigned k = 0; k <= 3; k++)
				Fc1[i][j][k] = S1[i][j] * (FL[k] + FR[k] - AARoe[k]) / 2;	
		}

    for ( i = 1; i <= maxI ; i++)    //注意范围
        for ( j = 1; j <= maxJ-1 ; j++)
            for ( unsigned k = 0; k < 4; k++)
            {   
                Fc2[i][j][k] = Fc4[i + 1][j][k];
                Fc3[i][j][k] = Fc1[i][j + 1][k];
                Residual[i][j][k] = -Fc1[i][j][k] + Fc2[i][j][k] - Fc4[i][j][k] + Fc3[i][j][k];
                if (Residual[i][j][k] > maxR[k])
                    maxR[k] = Residual[i][j][0];
            }
}
  
 
/***************************utility  **************************/
//将Q转化为rho u v p H 
void aeroConvert(Index i, Index j)
{
    rho[i][j] = Q[i][j][0];
    u[i][j] = Q[i][j][1] / rho[i][j];
    v[i][j] = Q[i][j][2] / rho[i][j];
    p[i][j] = (GAMMA - 1) * (Q[i][j][3] - 0.5*rho[i][j] * ( SQ(u[i][j]) + SQ(v[i][j]) ) );
    H[i][j] = (Q[i][j][3] + p[i][j]) / rho[i][j];  
}

void extrapolation(Index newI, Index newJ, Index oldI, Index oldJ)
{
  rho[newI][newJ] =rho[oldI][oldJ];
    u[newI][newJ] =  u[oldI][oldJ];
    v[newI][newJ] =  v[oldI][oldJ];
    p[newI][newJ] =  p[oldI][oldJ];
    H[newI][newJ] =  H[oldI][oldJ];
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
    for(unsigned i=cellBegin; i<=cellIEnd; i++)
    {
        for(unsigned j=cellBegin; j<=cellJEnd; j++)
        {
            fprintf(fpG, "%-3d  %-3d  %.3f %.3f %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e %.4e \n",
            i, j, mesh[i][j].x, mesh[i][j].y,
            volume[i][j], N2[i][j].x, N2[i][j].y, N3[i][j].x, N3[i][j].y,
            S1[i][j], S2[i][j], S3[i][j], S4[i][j] );
        } 
    }
    fclose(fpG);
}

void print()
{
    FILE *fpField;
    double T, Ma;
    fpField=fopen("field.dat", "w");
    fprintf(fpField, "Title=\"Field\"\nVariables=\"x\",\"y\",\"rho\",\
    \"u\",\"v\",\"p\",\"T\",\"Ma\"\nZone T=\"zone1\" i=%d,j=%d,f=point\n", maxI+1, maxJ+1);
	for (j = 1; j <= maxJ+1; j++)
	{
		for (i = 1; i <= maxI+1; i++)
		{
            T = p[i][j] / (rho[i][j] * 287.06);
            Ma = sqrt((SQ(u[i][j]) + SQ(v[i][j])) / (GAMMA * p[i][j] / rho[i][j]));
            fprintf(fpField, "%.3f %.3f %.4e %.4e %.4e %.4e %.4e %.4e \n",
                    mesh[i][j].x, mesh[i][j].y, rho[i][j], u[i][j], v[i][j],  p[i][j],
                     T, Ma);
        }
    }
    fclose(fpField);
}   



/***************************main  **************************/
int main()
{  
    fpDe=fopen("debug.txt","w");

	cout<<"\nCase 1 for inlet 1.8Ma; 2 for 1.5Ma\n\nCurrent Case: "<<caseNo<<endl;
    genMesh();//生成网格
    cellGeometry();//计算各个单元格面积,方向矢量等
    printMesh();//打印网格

    //初始化流场
    init1();

    cout<<"\nInitialzation done.\n\n"; 

    FILE *fpR;
    fpR=fopen("residual.dat", "w");
    fprintf(fpR,"iter  continuity x-velocity y-velocity Energy\n");

    for (step=1;  step<=2000; step++)
    {   
        BC();
        roeFlux();
        iteration();                  

        fprintf(fpR, "%-5d %.4e %.4e %.4e %.4e\n",
         step, maxR[0], maxR[1], maxR[2], maxR[3]);
        printf("%-5d %.4e %.4e %.4e %.4e\n",
         step, maxR[0], maxR[1], maxR[2], maxR[3]);
    }
    fclose(fpR);


    print();//打印结果
    return 0;
}



/*
            fprintf(fpDe,"T[%2d][%2d]=%.3e\t", i,j, T);
            fprintf(fpDe,"c[%2d][%2d]=%.3e\t", i,j, c);
            fprintf(fpDe,"rho[%2d][%2d]=%.3e\t", i,j, rho[i][j]);
            fprintf(fpDe,"p[%2d][%2d]=%.3e\t", i,j, p[i][j]);
            fprintf(fpDe,"u[%2d][%2d]=%.3e\t", i,j, u[i][j]);
            fprintf(fpDe,"v[%2d][%2d]=%12.3e\t", i,j, v[i][j]);
            fprintf(fpDe,"H[%2d][%2d]=%.3e\n", i,j, H[i][j]);
*/