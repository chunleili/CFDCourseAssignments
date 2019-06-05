#include <stdio.h>
#include <math.h>
int i;
int j;
double u[1000];
double rou[1000];
double p[1000];
double u_k[1000];
double rou_k[1000];
double p_k[1000];
double E[1000];
double F1[999];
double F2[999];
double F3[999];
void Jameson(double rouJ[], double uJ[], double pJ[]);
int main()
{
	FILE *fp;
	fp = fopen("output.txt", "w");
	double max(double a, double b);

	double R1[1000];
	double R2[1000];
	double R3[1000];
	double W1[1000];
	double W2[1000];
	double W3[1000];

	for (i = 0; i < 499; i++)
	{
		u[i] = 0;
		rou[i] = 1;
		p[i] = 1;
	}
	for (i = 499; i < 1000; i++)
	{
		u[i] = 0;
		rou[i] = 0.125;
		p[i] = 0.1;
	}

	u_k[0] = 0;
	rou_k[0] = 1;
	p_k[0] = 1;
	u_k[99] = 0;
	rou_k[99] = 0.125;
	p_k[99] = 0.1;

	for (j = 0; j < 1000; j++) //时间
	{
		Jameson(rou,u,p);
		for (i = 1; i < 999; i++)
		{
			W1[i] = W1[i] + 0.01 * (F1[i - 1] - F1[i]);
			W2[i] = W2[i] + 0.01 * (F2[i - 1] - F2[i]);
			W3[i] = W3[i] + 0.01 * (F3[i - 1] - F3[i]);
			rou_k[i]=W1[i];
			u_k[i]=W2[i]/W1[i];
			p_k[i]=0.4*(W3[i]-0.5*rou_k[i]* u_k[i]*u_k[i]);
		}

		Jameson(rou_k,u_k,p_k);
		for (i = 1; i < 999; i++)
		{
			W1[i] = 0.5 * rou[i] + 0.5 * rou_k[i] + 0.5 * 0.01 * (F1[i - 1] - F1[i]);
			W2[i] = 0.5 * rou[i] * u[i] + 0.5 * rou_k[i] * u_k[i] + 0.5 * 0.01 * (F2[i - 1] - F2[i]);
			W3[i] = 0.5 * E[i] + 0.5 * E[i] + 0.5 * 0.01 * (F3[i - 1] - F3[i]);
			//至此，得到新时刻下的rou,u,p
			rou[i] = W1[i]; //由于密度和压力恒为正数，处理数值计算中出现的负数，可能引入误差
			u[i] = W2[i] / W1[i];
			p[i] = 0.4 * (W3[i] - 0.5 * W2[i] * W2[i] / W1[i]);
		}
	}
	for (i = 0; i < 1000; i++)
	{
		fprintf(fp, "%.5f %.5f %.5f\n", rou[i], u[i], p[i]);
	}
}

double max(double a, double b)
{
	double max = a;
	if (b > a)
	{
		max = b;
	}
	return max;
}
void Jameson(double rouJ[], double uJ[], double pJ[])
{
	double gama[1000];
	double epis_2[999];
	double epis_4[999];
	double alpha[999];
	double c[1000];
	double D1[999];
	double D2[999];
	double D3[999];
	gama[0] = 0;
	gama[999] = 0;
	D1[0] = 0;
	D1[998] = 0;
	D2[0] = 0;
	D2[998] = 0;
	D3[0] = 0;
	D3[998] = 0;
	for (i = 0; i < 1000; i++)
		{

			c[i] = sqrt(1.4*p[i] / rou[i]);
			E[i] = (p[i] / 0.4) + 0.5*rou[i] * u[i] * u[i];
		}
		for (i = 1; i < 999; i++)
		{
			gama[i] = (fabs(p[i + 1] - 2 * p[i] + p[i - 1])) / (p[i + 1] + 2 * p[i] + p[i - 1]);
		}
		for (i = 0; i < 999; i++)//i代表面上时，处于节点i和i+1之间
		{
			epis_2[i] = 1*max(gama[i], gama[i+1]);
			epis_4[i] = max(0, 0.015625 - epis_2[i]);
			alpha[i] = ((fabs(u[i]) + c[i]) + (fabs(u[i + 1]) + c[i + 1]))*0.5;
		}
		for (i = 1; i < 998; i++)
		{
			D1[i] = alpha[i] * (epis_2[i] * (rou[i + 1] - rou[i]) - epis_4[i] * (rou[i + 2] - 3 * rou[i + 1] + 3 * rou[i] - rou[i - 1]));
			D2[i] = alpha[i] * (epis_2[i] * (rou[i + 1]* u[i+1] - rou[i]* u[i]) - epis_4[i] * (rou[i + 2]* u[i + 2] - 3*rou[i + 1] * u[i + 1] + 3 * rou[i]* u[i] - rou[i - 1]* u[i - 1]));
			D3[i] = alpha[i] * (epis_2[i] * ( E[i + 1] - E[i]) - epis_4[i] * ( E[i + 2] - 3 *  E[i + 1] + 3 * E[i] -  E[i - 1]));	
		}
		for (i = 0; i < 999; i++)//i代表面上时，处于节点i和i+1之间
		{
			F1[i] = 0.5*(rou[i + 1] * u[i + 1] + rou[i] * u[i])- D1[i];
			F2[i] = 0.5*(rou[i + 1] * u[i + 1] * u[i + 1] +p[i+1]+ rou[i] * u[i]*u[i]+p[i]) - D2[i];
			F3[i] = 0.5*( E[i + 1] * u[i + 1] + p[i + 1]*u[i+1] +  E[i] * u[i] + p[i]*u[i]) - D3[i];
		}
}