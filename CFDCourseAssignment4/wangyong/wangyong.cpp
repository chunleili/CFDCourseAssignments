#include <stdio.h>
#include <math.h>
#include <cstdlib>
#define Rg 287.06
#define CFL 0.8
#define SQ(a) ((a)*(a))

const unsigned maxI=250, maxJ=50;
const int cellBegin=1, cellIEnd=maxI, cellJEnd=maxJ;
const int block1=(int)(maxI*0.2+0.1), block2=(int)(maxI*0.2+0.1);

typedef double Field[maxI+1][maxJ+1][4];
typedef double ScalarField[maxI+1][maxJ+1];
typedef double Vector[4];
typedef struct XY
{
	double x,y;
}XY;
typedef XY MeshPoint[maxI+1][maxJ+1];

MeshPoint mesh;
int i, j, k;
ScalarField dyi, dxi, dyj, dxj;
ScalarField S1,S2,S3, S4, volume;
MeshPoint N1,N2,N3,N4;
ScalarField rho, p, u, v, T, Ma,H;

Field  Residual, Q;
double vi, vj, chvel, dt;
double  residualRho, residualU, residualV, residualE;
const double GAMMA = 1.4;
double dtGlobal=100;

Field Fc1,Fc2,Fc3,Fc4;

void genMesh();
void init();
void cellGeometry();
void BC();

const double p0=101325;
const double u0=624.9397;
const double v0=0;
const double T0=300;
const double rho0=1.176829268;
const double H0=496662.5;
const double E0=483117.6;
double nx,ny;
double c;


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

void init1()
{

    //初始全部给1.8Ma, 攻角为0, 压力给大气压101325, 静温300K, 其余推出如下:
    //声速c=sqrt(1.4*287*300)=347.1887
    //速度u=624.9397, v=0, 密度rho=p/RT=1.176829
    //单位体积总能E=p/(GAMMA-1)+0.5*rho*(u^2+v^2)
    for(unsigned i=0; i<cellIEnd; i++)//包括虚网格!
    {
        for(unsigned j=0; j<cellJEnd; j++)
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


void cellGeometry()
{
	double x1,x2,x3,x4,y1,y2,y3,y4;
	for(unsigned i=0; i<=maxI; i++)//注意范围
        for(unsigned j=0; j<=maxJ; j++)
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

			dyi[i][j] = y1 - y2;
			dxi[i][j] = x2 - x1;
			dyj[i][j] = y4 - y1;
			dxj[i][j] = x1 - x4;

			//dyi[i][j]=-N1[i][j].x*S1[i][j];
			//dxi[i][j]=-N1[i][j].y*S1[i][j];
			//dyj[i][j]=-N4[i][j].x*S4[i][j];
			//dxj[i][j]=-N4[i][j].y*S4[i][j];
			//面法向单位矢量, 实际上N1.x=sin(theta), N1.y=cos(theta), 再带上方向
			N1[i][j].x = (y2 - y1) / S1[i][j];
			N1[i][j].y = (x1 - x2) / S1[i][j];

			N2[i][j].x = (y3 - y2) / S2[i][j];
			N2[i][j].y = (x2 - x3) / S2[i][j];

			N3[i][j].x = (y4 - y3) / S3[i][j];
			N3[i][j].y = (x3 - x4) / S3[i][j];

			N4[i][j].x = (y1 - y4) / S4[i][j];
			N4[i][j].y = (x4 - x1) / S4[i][j];

			volume[i][j] =  0.5 * ((x1 - x3) * (y2 - y4) + (x4 - x2) * (y1 - y3));

		}
}





//在最左和最下侧分别铺设一层虚网格,最上和最右侧分别铺设两层虚网格
//0代表左下虚网格, maxI/J+1和maxI/J+2代表右上虚网格, 实际网格1~max, 共max个
//上下是壁面, 壁面法向速度为0,壁面切向不变 壁面无穿透边界
void BC()
{
    double pw;
    const unsigned j=cellJEnd;
    for(unsigned i=cellBegin; i<=cellIEnd; i++)
    {
        //pw=0.5*(3*p[i][j]-p[i][j-1]);//壁面的压力用两点外推
        pw=p[i][j];
        nx=N3[i][j].x; ny=N3[i][j].y;

        Fc1[i][j][0]=0;       
        Fc1[i][j][1]=0;
        Fc1[i][j][2]=pw;
        Fc1[i][j][3]=0;
	}
    //const unsigned j= cellBegin;
    for (unsigned i = cellBegin; i <= cellIEnd; i++)
    {
        nx=N3[i][1].x; ny=N3[i][1].y;
        //pw = 0.5 * (3 * p[i][j] - p[i][j+1]); //壁面的压力用两点外推
        pw=p[i][1];
        
        Fc1[i][0][0] = 0;     //注意,虚网格nx ny使用相邻网格的值
        Fc1[i][0][1] = pw*nx;
        Fc1[i][0][2] = pw*ny;
        Fc1[i][0][3] = 0;

        //extrapolation(i,0,i,1);
      rho[i][0] = rho[i][1];
        u[i][0] =   u[i][1];
        v[i][0] =   v[i][1];
        p[i][0] =   p[i][1];
        H[i][0] =   H[i][1];
    }

    for(unsigned j=cellBegin; j<=cellJEnd; j++)
    {
        p[0][j] = 101325;
        u[0][j] = 624.9397;
        v[0][j] = 0;
      rho[0][j] = 1.176829;
        H[0][j] = 496662.5; //Cp*300+0.5*625*625

        const double nx=1, ny=0;
        const double Vcv0=624.9397;
        Fc4[0][j][0]=rho0*Vcv0;
        Fc4[0][j][1]=rho0*u0*Vcv0+p0*nx;
        Fc4[0][j][2]=rho0*v0*Vcv0+p0*ny;
        Fc4[0][j][3]=rho0*H0*Vcv0;
  
    }
}
double rhoR,rhoL,uR,uL,vR,vL,pR,pL,HR,HL;
double AARoe[4];
double FR[4], FL[4];

void calARoe()
{

	double VcvR = uR * nx + vR * ny, VcvL = uL * nx + vL * ny;

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
	FR[0] = rhoR * VcvR;
	FR[1] = rhoR * VcvR * uR + nx * pR;
	FR[2] = rhoR * VcvR * vR + ny * pR;
	FR[3] = rhoR * VcvR * HR;

	FL[0] = rhoL * VcvL;
	FL[1] = rhoL * VcvL * uL + nx * pL;
	FL[2] = rhoL * VcvL * vL + ny * pL;
	FL[3] = rhoL * VcvL * HL;
}

void solve() //利用roe格式求解
{
	residualRho = 0;
	residualU = 0;
	residualV = 0;
	residualE = 0;
	
	//计算beta与AQ
	for (j = 1; j < maxJ+1; j++)
	{
		for (i = 1; i < maxI+1; i++)
		{   
			c = sqrt(GAMMA*p[i][j]/rho[i][j]);
			double T = p[i][j]/(287*rho[i][j]);
			if (T < 0)
			{
				printf("\n****** T<0!! T= %f\n ", T); //exit(2);
			}


			nx=-N4[i][j].x; 
			ny=-N4[i][j].y;
			
			pL = p[i - 1][j], pR = p[i][j];
			rhoL = rho[i - 1][j], rhoR = rho[i][j];
			uL = u[i - 1][j], uR = u[i][j];
			vL = v[i - 1][j], vR = v[i][j];
			HL = H[i - 1][j], HR = H[i][j];

			calARoe();


			for (unsigned k = 0; k <= 3; k++)
			{
				Fc4[i][j][k] = S4[i][j] * (FL[k] + FR[k] - AARoe[k]) / 2;
			}

			//nx=dyi[i][j] / S1[i][j];
			//ny=dxi[i][j] / S1[i][j];
			nx=-N1[i][j].x;
			ny=-N1[i][j].y;

			pL   = p[i][j-1],   pR = p[i][j];
			rhoL = rho[i][j-1],  rhoR = rho[i][j];
			uL   =  u[i][j-1],     uR =  u[i][j];
			vL   =  v[i][j-1],     vR =  v[i][j];
			HL   =    H[i][j-1],   HR =    H[i][j];
 
			calARoe();

			for (unsigned k = 0; k <= 3; k++)
			{
				Fc1[i][j][k] = S1[i][j] * (FL[k] + FR[k] - AARoe[k]) / 2;
			}

		}
	}


	for (j = 1; j < maxJ; j++)
	{
		for (i = 1; i < maxI; i++)
		{
			for (unsigned k = 0; k < 4; k++)
            {
				Fc2[i][j][k]  =Fc4[i+1][j][k];
				Fc3[i][j][k]  =Fc1[i][j+1][k];

                Residual[i][j][k] = -Fc1[i][j][k]  + Fc2[i][j][k]  - Fc4[i][j][k]  + Fc3[i][j][k] ;
            }

			if (Residual[i][j][0] > residualRho)
			{
				residualRho = Residual[i][j][0];
			}
			if (Residual[i][j][1] > residualU)
			{
				residualU = Residual[i][j][1];
			}
			if (Residual[i][j][2] > residualV)
			{
				residualV = Residual[i][j][2];
			}
			if (Residual[i][j][3] > residualE)
			{
				residualE = Residual[i][j][3];
			}
		}
	} //////////
}

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


void iteration()
{
dtGlobal=100;

	//计算当地时间步

	for (j = 1; j < maxJ; j++)
	{
		for (i = 1; i < maxI; i++)
		{

			vj = -N4[i][j].x*S4[i][j] * u[i][j]  -N4[i][j].y*S4[i][j] * v[i][j];
			vi = -N1[i][j].x*S1[i][j] * u[i][j]  -N1[i][j].y*S1[i][j] * v[i][j]; 
			c = sqrt(GAMMA * p[i][j] / rho[i][j]);
			chvel = fabs(vj) + fabs(vi) + c * (S4[i][j] + S1[i][j]);

			//dt=LTS();

			dt = CFL / chvel;
			if(dtGlobal>dt) dtGlobal=dt;

			Q[i][j][0] = Q[i][j][0] - dt * Residual[i][j][0]; //迭代求解下一步Q
			Q[i][j][1] = Q[i][j][1] - dt * Residual[i][j][1];
			Q[i][j][2] = Q[i][j][2] - dt * Residual[i][j][2];
			Q[i][j][3] = Q[i][j][3] - dt * Residual[i][j][3];

    		rho[i][j] = Q[i][j][0];
    		u[i][j] = Q[i][j][1] / rho[i][j];
    		v[i][j] = Q[i][j][2] / rho[i][j];
    		p[i][j] = (GAMMA - 1) * (Q[i][j][3] - 0.5*rho[i][j] * ( SQ(u[i][j]) + SQ(v[i][j]) ) );
    		H[i][j] = (Q[i][j][3] + p[i][j]) / rho[i][j];  

		}
	} /////////////////
}


int main()
{
	FILE *fr, *fg;

	genMesh();

	init();

	cellGeometry();

	fg = fopen("error.txt", "w");

	for (k = 0; k < 2000; k++)
	{

		printf("%d    %.15f    %.15f    %.15f    %.15f      %.5e\n", k, residualRho, residualU, residualV, residualE,  dtGlobal);
		fprintf(fg, "%d    %.15f    %.15f    %.15f    %.15f\n", k, residualRho, residualU, residualV, residualE);
		BC();
		solve();
		iteration();
	}


	fr = fopen("WYplotflow.dat", "w");
	fprintf(fr, "Title=\"CFD4\"\nVariables=\"x\",\"y\",\"rho\",\"velx\",\"vely\",\"pressure\",\"T\",\"Ma\"\nZone T=\"Zone1\" i=%d,j=%d,f=point \n",maxI+1,maxJ+1);
	for (j = 0; j < maxJ+1; j++)
	{
		for (i = 0; i < maxI+1; i++)
		{
			T[i][j] = p[i][j] / Rg / rho[i][j];
			Ma[i][j] = sqrt(u[i][j] * u[i][j] + v[i][j] * v[i][j]) / sqrt(GAMMA * Rg * T[i][j]);
			fprintf(fr, "%.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f\n", mesh[i][j].x, mesh[i][j].y, rho[i][j], u[i][j], v[i][j], p[i][j], T[i][j], Ma[i][j]);
		}
	} //////////////

	
}
