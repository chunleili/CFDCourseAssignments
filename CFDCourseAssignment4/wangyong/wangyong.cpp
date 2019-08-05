#include <stdio.h>
#include <math.h>
#include <cstdlib>
#define e 2.718281828459
#define a 2
#define gama 1.4
#define R 287.06
#define CFL 0.8
#define SQ(a) ((a)*(a))
FILE *fp, *fq, *fr, *fs, *fg;
const unsigned maxI=330, maxJ=70;

typedef double Field[maxI+1][maxJ+1][4];
typedef double ScalarField[maxI+1][maxJ+1];
typedef double Vector[4];
typedef double MeshPoint[maxI+1][maxJ+1][2];

MeshPoint nodes, nodesc1;
int i, j, k;
ScalarField dyi, dxi, dyj, dxj;
ScalarField S1,S2,S3, S4, area;
MeshPoint N1,N2,N3,N4;
ScalarField rho, pre1, u, v, T, Ma,H;

ScalarField total_pre1, total_T1;
Field  Residual, Q;
double tyj, txj, tyi, txi, tsli, tslj, vi, vj, sonic, chvel, dt;
double vnorm, vtemp, maxflux, maxflux2, maxflux3, maxflux4;
const double GAMMA = 1.4;
double dtGlobal=100;

Vector FcI, FcJ, FcIright,FcJup;
Field Fc1,Fc2,Fc3,Fc4;

void mesh_generation()
{

	for (j = 0; j < maxJ+1; j++)
	{
		for (i = 0; i < 11; i++)
		{
			nodes[i][j][0] = 0.02 * i - 1.2;
			nodes[i][j][1] = 0.01 * j;
		}
		for (i = 11; i < 61; i++)
		{
			nodes[i][j][0] = (log((i - 10) * (pow(e, a) - 1) / 50 + 1)) / a - 1;
			nodes[i][j][1] = (0.0179015 * pow(nodes[i][j][0], 6) - 0.1761616 * pow(nodes[i][j][0], 5) - 0.5135756 * pow(nodes[i][j][0], 4) - 0.02836162 * pow(nodes[i][j][0], 3) + 0.4913228 * pow(nodes[i][j][0], 2) + 0.000153865 * nodes[i][j][0] + 0.5000094) / 70 * j;
		}
		for (i = 61; i < 71; i++)
		{
			nodes[131 - i][j][0] = 0.08983227479893685 - 0.08983227479893685 * (log((i - 61) * (pow(e, a) - 1) / 10 + 1)) / a;
		}
		for (i = 61; i < 71; i++)
		{
			nodes[i][j][1] = (0.0179015 * pow(nodes[i][j][0], 6) - 0.1761616 * pow(nodes[i][j][0], 5) - 0.5135756 * pow(nodes[i][j][0], 4) - 0.02836162 * pow(nodes[i][j][0], 3) + 0.4913228 * pow(nodes[i][j][0], 2) + 0.000153865 * nodes[i][j][0] + 0.5000094) / 70 * j;
		}
		for (i = 71; i < 151; i++)
		{
			nodes[i][j][0] = (1 - 0.08983227479893685) * (log((i - 70) * (pow(e, a) - 1) / 80 + 1)) / a + 0.08983227479893685;
			nodes[i][j][1] = (0.08748864 * nodes[i][j][0] + 0.4960966) / 70 * j;
		}
		for (i = 151; i < 231; i++)
		{
			nodes[381 - i][j][0] = 2 - (log((i - 151) * (pow(e, a) - 1) / 80 + 1)) / a;
			nodes[i][j][1] = 0.583585301908299 / 70 * j;
		}
		for (i = 231; i < maxI+1; i++)
		{
			nodes[i][j][0] = 0.04 * (i - 230) + 2;
			nodes[i][j][1] = 0.583585301908299 / 70 * j;
		}
	}

}
void initialize()
{

	for (j = 0; j < maxJ+1; j++)
	{
		for (i = 0; i < maxI+1; i++)
		{
			u[i][j] = 400; //3
			v[i][j] = 0;
			pre1[i][j] = 101325;
			T[i][j] = 230;
			rho[i][j] = pre1[i][j] / T[i][j] / R;
			Ma[i][j] = sqrt(u[i][j] * u[i][j] + v[i][j] * v[i][j]) / sqrt(gama * R * T[i][j]);
			H[i][j] = (pre1[i][j] / (gama - 1) + 0.5 * rho[i][j] * (pow(v[i][j], 2) + pow(u[i][j], 2)) +
			pre1[i][j]) / rho[i][j];
		}
	}
	for (j = 1; j < maxJ; j++)
	{
		for (i = 1; i < maxI; i++)
		{
			Q[i][j][0] = rho[i][j];
			Q[i][j][1] = rho[i][j] * u[i][j];
			Q[i][j][2] = rho[i][j] * v[i][j];
			Q[i][j][3] = pre1[i][j] / (gama - 1) + 0.5 * rho[i][j] * (u[i][j] * u[i][j] + v[i][j] * v[i][j]);
		}
	}
}


void area_calculate()
{
	double x1,x2,x3,x4,y1,y2,y3,y4;
	for(unsigned i=0; i<=maxI; i++)//注意范围
        for(unsigned j=0; j<=maxJ; j++)
		{
			x1 = nodes[i][j][0];
			x2 = nodes[i + 1][j][0];
			x3 = nodes[i + 1][j + 1][0];
			x4 = nodes[i][j + 1][0];

			y1 = nodes[i][j][1];
			y2 = nodes[i + 1][j][1];
			y3 = nodes[i + 1][j + 1][1];
			y4 = nodes[i][j + 1][1];

			S1[i][j] = sqrt(SQ(x1 - x2) + SQ(y1 - y2)); //下侧面积S1
			S2[i][j] = sqrt(SQ(x2 - x3) + SQ(y2 - y3)); //右侧面积S2
			S3[i][j] = sqrt(SQ(x3 - x4) + SQ(y3 - y4)); //上侧面积S3
			S4[i][j] = sqrt(SQ(x4 - x1) + SQ(y4 - y1)); //左侧面积S4

			dyi[i][j] = nodes[i][j][1] - nodes[i + 1][j][1];
			dxi[i][j] = nodes[i + 1][j][0] - nodes[i][j][0];
			dyj[i][j] = nodes[i][j + 1][1] - nodes[i][j][1];
			dxj[i][j] = nodes[i][j][0] - nodes[i][j + 1][0];

			//面法向单位矢量, 实际上N1.x=sin(theta), N1[1]=cos(theta), 再带上方向
			N1[i][j][0] = (y2 - y1) / S1[i][j];
			N1[i][j][1] = (x1 - x2) / S1[i][j];

			N2[i][j][0] = (y3 - y2) / S2[i][j];
			N2[i][j][1] = (x2 - x3) / S2[i][j];

			N3[i][j][0] = (y4 - y3) / S3[i][j];
			N3[i][j][1] = (x3 - x4) / S3[i][j];

			N4[i][j][0] = (y1 - y4) / S4[i][j];
			N4[i][j][1] = (x4 - x1) / S4[i][j];

			area[i][j] =  0.5 * ((x1 - x3) * (y2 - y4) + (x4 - x2) * (y1 - y3));

		}
}



double pai(double Ma)
{
	double paima;
	paima = pow(1 / (1 + (gama - 1) * Ma * Ma / 2), gama / (gama - 1));
	return paima;
}

void boundary_conditions()
{
	for (j = 1; j < maxJ; j++) //进口边界条件 分区1 亚音进口，一个变量外推
	{

		total_pre1[0][j] = 250000;
		total_T1[0][j] = maxI;
		v[0][j] = 0;

		T[0][j] = T[1][j];
		//静温外推
		Ma[0][j] = sqrt((total_T1[0][j] / T[0][j] - 1) * 2 / (gama - 1)); //壁面速度为0
		pre1[0][j] = total_pre1[0][j] * pai(Ma[0][j]);
		rho[0][j] = pre1[0][j] / R / T[0][j];
		u[0][j] = Ma[0][j] * sqrt(gama * R * T[0][j]);
		
	}

	for (j = 1; j < maxJ; j++)
	{
		if (Ma[maxI-1][j] >= 1) //超音全部外推
		{
			u[maxI][j] = u[maxI-1][j];
			v[maxI][j] = v[maxI-1][j];
			pre1[maxI][j] = pre1[maxI-1][j];
			T[maxI][j] = T[maxI-1][j];
			rho[maxI][j] = pre1[maxI][j] / R / T[maxI][j];
			Ma[maxI][j] = sqrt(u[maxI][j] * u[maxI][j] + v[maxI][j] * v[maxI][j]) / sqrt(gama * R * T[maxI][j]);
		}
		else //亚音3个外推   静压给定
		{
			pre1[maxI][j] = 85419;
			u[maxI][j] = u[maxI-1][j];
			v[maxI][j] = v[maxI-1][j];
			T[maxI][j] = T[maxI-1][j];
			rho[maxI][j] = pre1[maxI][j] / R / T[maxI][j];
			Ma[maxI][j] = sqrt(u[maxI][j] * u[maxI][j] + v[maxI][j] * v[maxI][j]) / sqrt(gama * R * T[maxI][j]);
		}
	}

	for (i = 1; i < 151; i++) //喷管壁面边界条件
	{

		double nx=dyi[i][maxJ] / S1[i][maxJ];
		double ny=dxi[i][maxJ] / S1[i][maxJ];

		Fc1[i][maxJ][0]=0;
		Fc1[i][maxJ][1]=pre1[i][maxJ]*nx;
		Fc1[i][maxJ][2]=pre1[i][maxJ]*ny;
		Fc1[i][maxJ][3]=0;			
	}

	for (i = 1; i < maxI; i++) //对称边界条件
	{
		rho[i][0] = rho[i][1];
		u[i][0] = u[i][1];
		v[i][0] = 0;
		pre1[i][0] = pre1[i][1];
		T[i][0] = pre1[i][0] / rho[i][0] / R;
		Ma[i][0] = sqrt(u[i][0] * u[i][0] + v[i][0] * v[i][0]) / sqrt(gama * R * T[i][0]);
	}

	////////////////////////////

	for (i = 1; i < 81; i++) //下壁面条件
	{

		double nx=dyi[i][0] / S1[i][0];
		double ny=dxi[i][0] / S1[i][0];

		Fc1[i][0][0]=0;
		Fc1[i][0][1]=pre1[i][0]*nx;
		Fc1[i][0][2]=pre1[i][0]*ny;
		Fc1[i][0][3]=0;	

	}
}

double rhoR,rhoL,uR,uL,vR,vL,nx,ny,pR,pL,HR,HL;
double AARoe[4];
double FR[4], FL[4];
double c;

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

void roe() //利用roe格式求解
{
	maxflux = 0;
	maxflux2 = 0;
	maxflux3 = 0;
	maxflux4 = 0;
	
	//计算beta与AQ
	for (j = 1; j < maxJ+1; j++)
	{
		for (i = 1; i < maxI+1; i++)
		{   
			c = sqrt(GAMMA*pre1[i][j]/rho[i][j]);
			double T = pre1[i][j]/(287*rho[i][j]);
			if (T < 0)
			{
				printf("\n****** T<0!! T= %f\n ", T); //exit(2);
			}

			nx=dyj[i][j] / S4[i][j];
			ny=dxj[i][j] / S4[i][j];
			
			pL = pre1[i - 1][j], pR = pre1[i][j];
			rhoL = rho[i - 1][j], rhoR = rho[i][j];
			uL = u[i - 1][j], uR = u[i][j];
			vL = v[i - 1][j], vR = v[i][j];
			HL = H[i - 1][j], HR = H[i][j];

			calARoe();


			for (unsigned k = 0; k <= 3; k++)
			{
				Fc4[i][j][k] = S4[i][j] * (FL[k] + FR[k] - AARoe[k]) / 2;
			}

			nx=dyi[i][j] / S1[i][j];
			ny=dxi[i][j] / S1[i][j];

			pL   = pre1[i][j-1],   pR = pre1[i][j];
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

			if (Residual[i][j][0] > maxflux)
			{
				maxflux = Residual[i][j][0];
			}
			if (Residual[i][j][1] > maxflux2)
			{
				maxflux2 = Residual[i][j][1];
			}
			if (Residual[i][j][2] > maxflux3)
			{
				maxflux3 = Residual[i][j][2];
			}
			if (Residual[i][j][3] > maxflux4)
			{
				maxflux4 = Residual[i][j][3];
			}
		}
	} //////////
}

//LTS=LocalTimeStepping,当地时间步法,返回当地时间步
/*
double LTS()
{
    double dtLocal;
    double lambdaI, lambdaJ, SI, SJ;
    double NI[2], NJ[2];
    NI[0] = 0.5 * (N3[I][J][0] - N1[I][J][0]);
    NI[1] = 0.5 * (N3[I][J][1] - N1[I][J][1]);

    NJ[0] = 0.5 * (N2[I][J][0] - N4[I][J][0]);
    NJ[1] = 0.5 * (N2[I][J][1] - N4[I][J][1]);

    SI = 0.5 * (S1[I][J] + S3[I][J]);
    SJ = 0.5 * (S2[I][J] + S4[I][J]);

    lambdaI = (u[I][J] * NI[0] + v[I][J] * NI[1])+ c;
    lambdaJ = (u[I][J] * NJ[0] + v[I][J] * NJ[1])+ c;

    dtLocal = CFL * area[I][J] / (lambdaI * SI + lambdaJ * SJ);

    return dtLocal;
}
*/

void iteration()
{
dtGlobal=100;

	//计算当地时间步

	for (j = 1; j < maxJ; j++)
	{
		for (i = 1; i < maxI; i++)
		{
			tyj = 0.5 * (dyj[i][j] + dyj[i + 1][j]);
			txj = 0.5 * (dxj[i][j] + dxj[i + 1][j]);
			tslj = 0.5 * (S4[i][j] + S4[i + 1][j]);
			tyi = 0.5 * (dyi[i][j] + dyi[i][j + 1]);
			txi = 0.5 * (dxi[i][j] + dxi[i][j + 1]);
			tsli = 0.5 * (S1[i][j] + S1[i][j + 1]);
			vj = tyj * u[i][j] + txj * v[i][j];
			vi = tyi * u[i][j] + txi * v[i][j];
			sonic = sqrt(gama * pre1[i][j] / rho[i][j]);
			chvel = fabs(vj) + fabs(vi) + sonic * (tslj + tsli);
			dt = CFL / chvel;
			if(dtGlobal>dt) dtGlobal=dt;

			Q[i][j][0] = Q[i][j][0] - dt * Residual[i][j][0]; //迭代求解下一步Q
			Q[i][j][1] = Q[i][j][1] - dt * Residual[i][j][1];
			Q[i][j][2] = Q[i][j][2] - dt * Residual[i][j][2];
			Q[i][j][3] = Q[i][j][3] - dt * Residual[i][j][3];

    		rho[i][j] = Q[i][j][0];
    		u[i][j] = Q[i][j][1] / rho[i][j];
    		v[i][j] = Q[i][j][2] / rho[i][j];
    		pre1[i][j] = (GAMMA - 1) * (Q[i][j][3] - 0.5*rho[i][j] * ( SQ(u[i][j]) + SQ(v[i][j]) ) );
    		H[i][j] = (Q[i][j][3] + pre1[i][j]) / rho[i][j];  

		}
	} /////////////////
}


int main()
{
	mesh_generation();

	initialize();

	area_calculate();

	fg = fopen("error.txt", "w");

	for (k = 0; k < 2000; k++)
	{

		printf("%d    %.15f    %.15f    %.15f    %.15f      %.5e\n", k, maxflux, maxflux2, maxflux3, maxflux4,  dtGlobal);
		fprintf(fg, "%d    %.15f    %.15f    %.15f    %.15f\n", k, maxflux, maxflux2, maxflux3, maxflux4);
		boundary_conditions();
		roe();
		iteration();
	}

	fp = fopen("1.txt", "w");

	fprintf(fp, "x                 ,y              rho[i][j],    u[i][j],    v[i][j],      pre1[i][j],      T[i][j],      Ma[i][j]\n");
	for (i = 1; i < maxI; i++)
	{
		for (j = 1; j < maxJ; j++)
		{
			T[i][j] = pre1[i][j] / R / rho[i][j];
			Ma[i][j] = sqrt(u[i][j] * u[i][j] + v[i][j] * v[i][j]) / sqrt(gama * R * T[i][j]);
			fprintf(fp, "%.10f    %.10f    %.10f    %.10f    %.10f    %.10f    %.10f    %.10f\n", nodes[i][j][0], nodes[i][j][1], rho[i][j], u[i][j], v[i][j], pre1[i][j], T[i][j], Ma[i][j]);
		}
	}

	fr = fopen("WYplotflow.dat", "w");
	fprintf(fr, "Title=\"NOZZLE\"\nVariables=\"x\",\"y\",\"dens\",\"velx\",\"vely\",\"spre\",\"ttem\",\"mach\"\nZone T=\"NOZZLE\" i=%d,j=%d,f=point \n",maxI+1,maxJ+1);
	for (j = 0; j < maxJ+1; j++)
	{
		for (i = 0; i < maxI+1; i++)
		{
			T[i][j] = pre1[i][j] / R / rho[i][j];
			Ma[i][j] = sqrt(u[i][j] * u[i][j] + v[i][j] * v[i][j]) / sqrt(gama * R * T[i][j]);
			fprintf(fr, "%.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f\n", nodes[i][j][0], nodes[i][j][1], rho[i][j], u[i][j], v[i][j], pre1[i][j], T[i][j], Ma[i][j]);
		}
	} //////////////

	
}
