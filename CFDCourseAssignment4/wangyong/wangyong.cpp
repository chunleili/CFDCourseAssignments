#include <stdio.h>
#include <math.h>
#include <cstdlib>
#define e 2.718281828459
#define a 2
#define gama 1.4
#define R 287.06
#define cfl 0.8
#define SQ(a) ((a)*(a))
FILE *fp, *fq, *fr, *fs, *fg;
double nodes[331][71][2], nodesc1[331][71][2];
int i, j, k, m, maxi, maxj, max;
int maxi2, maxj2, max2;
double dyi[331][71], dxi[331][71], dyj[331][71], dxj[331][71];
double sli[330][71], slj[331][71], area[330][71];
double pho1[331][71], pre1[331][71], vx1[331][71], vy1[331][71], T1[331][71], ma1[331][71];

double total_pre1[331][71], total_T1[331][71];
double  H[331][71], Fjr[331][71][4], Fjl[331][71][4], Fir[331][71][4], Fil[331][71][4];
double AQi[331][71][4], AQj[331][71][4], Flux[331][71][4], Q[331][71][4];
double resm1, resave[4], imax, jmax, tyj, txj, tyi, txi, tsli, tslj, vi, vj, sonic, chvel, dt, tres1, real[331];
double vnorm, vtemp, maxflux, maxflux2, maxflux3, maxflux4;
const double GAMMA = 1.4;
double dtGlobal=100;

void mesh_generation()
{

	for (j = 0; j < 71; j++)
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
		for (i = 231; i < 331; i++)
		{
			nodes[i][j][0] = 0.04 * (i - 230) + 2;
			nodes[i][j][1] = 0.583585301908299 / 70 * j;
		}
	}

}
void initialize()
{
	for (j = 0; j < 71; j++)
	{
		for (i = 0; i < 331; i++)
		{
			vx1[i][j] = 400; //3
			vy1[i][j] = 0;
			pre1[i][j] = 101325;
			T1[i][j] = 230;
			pho1[i][j] = pre1[i][j] / T1[i][j] / R;
			ma1[i][j] = sqrt(vx1[i][j] * vx1[i][j] + vy1[i][j] * vy1[i][j]) / sqrt(gama * R * T1[i][j]);
			H[i][j] = (pre1[i][j] / (gama - 1) + 0.5 * pho1[i][j] * (pow(vy1[i][j], 2) + pow(vx1[i][j], 2)) +
			 pre1[i][j]) / pho1[i][j];
		}
	}
}
void area_caculate()
{
	for (j = 0; j < 71; j++)
	{
		for (i = 0; i < 330; i++)
		{
			dyi[i][j] = nodes[i][j][1] - nodes[i + 1][j][1];
			dxi[i][j] = nodes[i + 1][j][0] - nodes[i][j][0];
			sli[i][j] = sqrt(dyi[i][j] * dyi[i][j] + dxi[i][j] * dxi[i][j]);
		}
	}
	for (j = 0; j < 70; j++)
	{
		for (i = 0; i < 331; i++)
		{
			dyj[i][j] = nodes[i][j + 1][1] - nodes[i][j][1];
			dxj[i][j] = nodes[i][j][0] - nodes[i][j + 1][0];
			slj[i][j] = sqrt(dyj[i][j] * dyj[i][j] + dxj[i][j] * dxj[i][j]);
		}
	}
	for (j = 0; j < 70; j++)
	{
		for (i = 0; i < 330; i++)
		{
			area[i][j] = (dyj[i][j] + dyj[i + 1][j]) * dxi[i][j] / 2;
		}
	} ///////////////////////////////


}


double pai(double Ma)
{
	double paima;
	paima = pow(1 / (1 + (gama - 1) * Ma * Ma / 2), gama / (gama - 1));
	return paima;
}

void boundary_conditions()
{
	for (j = 1; j < 70; j++) //进口边界条件 分区1 亚音进口，一个变量外推
	{

		total_pre1[0][j] = 250000;
		total_T1[0][j] = 330;
		vy1[0][j] = 0;

		T1[0][j] = T1[1][j];
		//静温外推
		ma1[0][j] = sqrt((total_T1[0][j] / T1[0][j] - 1) * 2 / (gama - 1)); //壁面速度为0
		pre1[0][j] = total_pre1[0][j] * pai(ma1[0][j]);
		pho1[0][j] = pre1[0][j] / R / T1[0][j];
		vx1[0][j] = ma1[0][j] * sqrt(gama * R * T1[0][j]);
		
	/*
		pho1[0][j]=1.176829;
		vx1[]
		vy1[0][j] = 0;
		*/
	}

	for (j = 1; j < 70; j++)
	{
		if (ma1[329][j] >= 1) //超音全部外推
		{
			vx1[330][j] = vx1[329][j];
			vy1[330][j] = vy1[329][j];
			pre1[330][j] = pre1[329][j];
			T1[330][j] = T1[329][j];
			pho1[330][j] = pre1[330][j] / R / T1[330][j];
			ma1[330][j] = sqrt(vx1[330][j] * vx1[330][j] + vy1[330][j] * vy1[330][j]) / sqrt(gama * R * T1[330][j]);
		}
		else //亚音3个外推   静压给定
		{
			pre1[330][j] = 85419;
			vx1[330][j] = vx1[329][j];
			vy1[330][j] = vy1[329][j];
			T1[330][j] = T1[329][j];
			pho1[330][j] = pre1[330][j] / R / T1[330][j];
			ma1[330][j] = sqrt(vx1[330][j] * vx1[330][j] + vy1[330][j] * vy1[330][j]) / sqrt(gama * R * T1[330][j]);
		}
	}

	for (i = 1; i < 151; i++) //喷管壁面边界条件
	{

		double nx=dyi[i][70] / sli[i][70];
		double ny=dxi[i][70] / sli[i][70];

		Fjl[i][70][0]=0;
		Fjl[i][70][1]=pre1[i][70]*nx;
		Fjl[i][70][0]=pre1[i][70]*ny;
		Fjl[i][70][1]=0;			
	}

	for (i = 1; i < 330; i++) //对称边界条件
	{
		pho1[i][0] = pho1[i][1];
		vx1[i][0] = vx1[i][1];
		vy1[i][0] = 0;
		pre1[i][0] = pre1[i][1];
		T1[i][0] = pre1[i][0] / pho1[i][0] / R;
		ma1[i][0] = sqrt(vx1[i][0] * vx1[i][0] + vy1[i][0] * vy1[i][0]) / sqrt(gama * R * T1[i][0]);
	}

	////////////////////////////

	for (i = 1; i < 81; i++) //下壁面条件
	{

		double nx=dyi[i][0] / sli[i][0];
		double ny=dxi[i][0] / sli[i][0];

		Fjl[i][0][0]=0;
		Fjl[i][0][1]=pre1[i][0]*nx;
		Fjl[i][0][0]=pre1[i][0]*ny;
		Fjl[i][0][1]=0;	

	}
}

double rhoR,rhoL,uR,uL,vR,vL,nx,ny,pR,pL,HR,HL;
double AARoe[4];
double FR[4], FL[4];
double c;
double SI, SJ;

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

	//计算beta与AQ
	for (j = 1; j < 71; j++)
	{
		for (i = 1; i < 331; i++)
		{   
			nx=dyj[i][j] / slj[i][j];
			ny=dxj[i][j] / slj[i][j];
			SI=slj[i][j];
			
			pL = pre1[i - 1][j], pR = pre1[i][j];
			rhoL = pho1[i - 1][j], rhoR = pho1[i][j];
			uL = vx1[i - 1][j], uR = vx1[i][j];
			vL = vy1[i - 1][j], vR = vy1[i][j];
			HL = H[i - 1][j], HR = H[i][j];

			c = sqrt(GAMMA*pre1[i][j]/pho1[i][j]);
			double T = pre1[i][j]/(287*pho1[i][j]);
			if (T < 0)
			{
				printf("\n****** T<0!! T= %f\n ", T); //exit(2);
			}

			calARoe();

			for (unsigned k = 0; k <= 3; k++)
				AQi[i][j][k] = SI*AARoe[k];

			for (unsigned k = 0; k <= 3; k++)
			{
				Fil[i][j][k] = SI * FL[k];
				Fir[i][j][k] = SI * FR[k];
			}
		}
	} 
	
	///////////////////////


	for (j = 1; j < 71; j++)
	{
		for (i = 1; i < 331; i++)
		{   
			nx=dyi[i][j] / sli[i][j];
			ny=dxi[i][j] / sli[i][j];

			SJ=sli[i][j];
			
			pL   = pre1[i][j-1],   pR = pre1[i][j];
			rhoL = pho1[i][j-1], rhoR = pho1[i][j];
			uL   =  vx1[i][j-1],   uR =  vx1[i][j];
			vL   =  vy1[i][j-1],   vR =  vy1[i][j];
			HL   =    H[i][j-1],   HR =    H[i][j];
 
			calARoe();

			for (unsigned k = 0; k <= 3; k++)
				AQj[i][j][k] = SJ*AARoe[k];

			for (unsigned k = 0; k <= 3; k++)
			{
				Fjl[i][j][k] = SJ * FL[k];
				Fjr[i][j][k] = SJ * FR[k];
			}

		}
	}
	maxflux = 0;
	maxi = 0;
	maxj = 0;
	max = 0;
	maxflux2 = 0;
	maxflux3 = 0;
	maxflux4 = 0;
	maxi2 = 0;
	maxj2 = 0;
	max2 = 0;

	double FcI[4], FcJ[4], FcIright[4],FcJup[4];
	
	for (j = 1; j < 70; j++)
	{
		for (i = 1; i < 330; i++)
		{
			for (unsigned k = 0; k < 4; k++)
            {
                //R代表离开单元格的通量的矢量和, =右侧+左侧+上侧+下侧
                //左侧与下侧通量分别由临近单元格右侧与上侧取负号得来
				FcI[k]  =(Fil[i][j][k] + Fir[i][j][k] - AQi[i][j][k]) / 2;
				FcIright[k]=(Fil[i + 1][j][k] + Fir[i + 1][j][k] - AQi[i + 1][j][k]) / 2;
				FcJ[k]  =(Fjl[i][j][k] + Fjr[i][j][k] - AQj[i][j][k]) / 2;
				FcJup[k]=(Fjl[i][j + 1][k] + Fjr[i][j + 1][k] - AQj[i][j + 1][k]) / 2;

                Flux[i][j][k] = -FcI[k]  + FcIright[k]  - FcJ[k]  + FcJup[k] ;
            }

			if (Flux[i][j][0] > maxflux)
			{
				maxflux = Flux[i][j][0];
			}
			if (Flux[i][j][1] > maxflux2)
			{
				maxflux2 = Flux[i][j][1];
			}
			if (Flux[i][j][2] > maxflux3)
			{
				maxflux3 = Flux[i][j][2];
			}
			if (Flux[i][j][3] > maxflux4)
			{
				maxflux4 = Flux[i][j][3];
			}
		}
	} //////////
}

void iteration()
{
dtGlobal=100;
	for (j = 1; j < 70; j++)
	{
		for (i = 1; i < 330; i++)
		{
			Q[i][j][0] = pho1[i][j];
			Q[i][j][1] = pho1[i][j] * vx1[i][j];
			Q[i][j][2] = pho1[i][j] * vy1[i][j];
			Q[i][j][3] = pre1[i][j] / (gama - 1) + 0.5 * pho1[i][j] * (vx1[i][j] * vx1[i][j] + vy1[i][j] * vy1[i][j]);
		}
	}
	resm1 = 0;
	for (i = 0; i < 4; i++)
	{
		resave[i] = 0;
	}
	imax = 0;
	jmax = 0;

	//计算当地时间步

	for (j = 1; j < 70; j++)
	{
		for (i = 1; i < 330; i++)
		{
			tyj = 0.5 * (dyj[i][j] + dyj[i + 1][j]);
			txj = 0.5 * (dxj[i][j] + dxj[i + 1][j]);
			tslj = 0.5 * (slj[i][j] + slj[i + 1][j]);
			tyi = 0.5 * (dyi[i][j] + dyi[i][j + 1]);
			txi = 0.5 * (dxi[i][j] + dxi[i][j + 1]);
			tsli = 0.5 * (sli[i][j] + sli[i][j + 1]);
			vj = tyj * vx1[i][j] + txj * vy1[i][j];
			vi = tyi * vx1[i][j] + txi * vy1[i][j];
			sonic = sqrt(gama * pre1[i][j] / pho1[i][j]);
			chvel = fabs(vj) + fabs(vi) + sonic * (tslj + tsli);
			dt = cfl / chvel;
			if(dtGlobal>dt) dtGlobal=dt;

			Q[i][j][0] = Q[i][j][0] - dt * Flux[i][j][0]; //迭代求解下一步Q
			Q[i][j][1] = Q[i][j][1] - dt * Flux[i][j][1];
			Q[i][j][2] = Q[i][j][2] - dt * Flux[i][j][2];
			Q[i][j][3] = Q[i][j][3] - dt * Flux[i][j][3];

    		pho1[i][j] = Q[i][j][0];
    		vx1[i][j] = Q[i][j][1] / pho1[i][j];
    		vy1[i][j] = Q[i][j][2] / pho1[i][j];
    		pre1[i][j] = (GAMMA - 1) * (Q[i][j][3] - 0.5*pho1[i][j] * ( SQ(vx1[i][j]) + SQ(vy1[i][j]) ) );
    		H[i][j] = (Q[i][j][3] + pre1[i][j]) / pho1[i][j];  

		}
	} /////////////////
}


int main()
{
	mesh_generation();

	initialize();

	area_caculate();

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

	fprintf(fp, "x                 ,y              pho1[i][j],    vx1[i][j],    vy1[i][j],      pre1[i][j],      T1[i][j],      ma1[i][j]\n");

	for (i = 1; i < 330; i++)
	{
		for (j = 1; j < 70; j++)
		{
			T1[i][j] = pre1[i][j] / R / pho1[i][j];
			ma1[i][j] = sqrt(vx1[i][j] * vx1[i][j] + vy1[i][j] * vy1[i][j]) / sqrt(gama * R * T1[i][j]);
			fprintf(fp, "%.10f    %.10f    %.10f    %.10f    %.10f    %.10f    %.10f    %.10f\n", nodes[i][j][0], nodes[i][j][1], pho1[i][j], vx1[i][j], vy1[i][j], pre1[i][j], T1[i][j], ma1[i][j]);
		}
	}


	//record();

	fr = fopen("WYplotflow.dat", "w");
	fprintf(fr, "Title=\"NOZZLE\"\nVariables=\"x\",\"y\",\"dens\",\"velx\",\"vely\",\"spre\",\"ttem\",\"mach\"\nZone T=\"NOZZLE\" i=331,j=71,f=point \n");
	for (j = 0; j < 71; j++)
	{
		for (i = 0; i < 331; i++)
		{
			T1[i][j] = pre1[i][j] / R / pho1[i][j];
			ma1[i][j] = sqrt(vx1[i][j] * vx1[i][j] + vy1[i][j] * vy1[i][j]) / sqrt(gama * R * T1[i][j]);
			fprintf(fr, "%.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f\n", nodes[i][j][0], nodes[i][j][1], pho1[i][j], vx1[i][j], vy1[i][j], pre1[i][j], T1[i][j], ma1[i][j]);
		}
	} //////////////

	
}
