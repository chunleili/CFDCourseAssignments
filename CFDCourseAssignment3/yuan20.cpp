#include <stdio.h>
#include <math.h>
int main()
{
	FILE *fp;
	fp = fopen("output.txt", "w");
	double max(double a, double b);
	int i;
	int j;
	double space[1000];
	double u[1000];
	double rou[1000];
	double p[1000];
	double u_k[1000];
	double rou_k[1000];
	double p_k[1000];
	double c[1000];
	double E[1000];
	double R1[1000];
	double R2[1000];
	double R3[1000];
	double W1[1000];
	double W2[1000];
	double W3[1000];
	double gama[1000];
	double epis_2[999];
	double epis_4[999];
	double alpha[999];
	double D1[999];
	double D2[999];
	double D3[999];
	double F1[999];
	double F2[999];
	double F3[999];
	
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
	gama[0] = 0;
	gama[999] = 0;
	D1[0] = 0;
	D1[998] = 0;
	D2[0] = 0;
	D2[998] = 0;
	D3[0] = 0;
	D3[998] = 0;
	u_k[0] = 0;
	rou_k[0] = 1;
	p_k[0] = 1;
	u_k[999] = 0;
	rou_k[999] = 0.125;
	p_k[999] = 0.1;

	for (j = 0; j <1000; j++)//时间
	{
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
		for (i = 1; i < 999; i++)
		{
			R1[i] = rou[i] + 0.1*(F1[i - 1] - F1[i]);
			R2[i] = rou[i] * u[i] + 0.1*(F2[i - 1] - F2[i]);
			R3[i] =  E[i] + 0.1*(F3[i - 1] - F3[i]);
			//由R反求暂时的u_k,p_k,rou_k，龙格库塔法里的值
			rou_k[i] =  R1[i];
			u_k[i] = R2[i] / R1[i];
			p_k[i] = 0.4*(R3[i]  - 0.5*R2[i]*R2[i] / R1[i]);
		}
		
		//将系数全更新
		for (i = 0; i < 1000; i++)
		{
			c[i] = sqrt(1.4*p_k[i] / rou_k[i]);
			E[i] = (p_k[i] / 0.4) + 0.5*rou_k[i] * u_k[i] * u_k[i];
		}
		for (i = 1; i < 999; i++)
		{
			gama[i] = (fabs(p_k[i + 1] - 2 * p_k[i] + p_k[i - 1])) / (p_k[i + 1] + 2 * p_k[i] + p_k[i - 1]);

		}
		for (i = 0; i < 999; i++)//i代表面上时，处于节点i和i+1之间
		{
			epis_2[i] = 1*max(gama[i], gama[i + 1]);
			epis_4[i] = max(0, 0.015625 - epis_2[i]);
			alpha[i] = ((fabs(u_k[i]) + c[i]) + (fabs(u_k[i + 1]) + c[i + 1]))*0.5;
		}
		for (i = 1; i < 998; i++)
		{
			D1[i] = alpha[i] * (epis_2[i] * (rou_k[i + 1] - rou_k[i]) - epis_4[i] * (rou_k[i + 2] - 3 * rou_k[i + 1] + 3 * rou_k[i] - rou_k[i - 1]));
			D2[i] = alpha[i] * (epis_2[i] * (rou_k[i + 1] * u_k[i + 1] - rou_k[i] * u_k[i]) - epis_4[i] * (rou_k[i + 2] * u_k[i + 2] - 3 * rou_k[i + 1] * u_k[i + 1] + 3 * rou_k[i] * u_k[i] - rou[i - 1] * u[i - 1]));
			D3[i] = alpha[i] * (epis_2[i] * ( E[i + 1] -  E[i]) - epis_4[i] * ( E[i + 2] - 3 *  E[i + 1] + 3 * E[i] - E[i - 1]));
		}
		for (i = 0; i < 999; i++)//i代表面上时，处于节点i和i+1之间
		{
			F1[i] = 0.5*(rou_k[i + 1] * u_k[i + 1] + rou_k[i] * u_k[i]) - D1[i];
			F2[i] = 0.5*(rou_k[i + 1] * u_k[i + 1] * u_k[i + 1] + p_k[i + 1] + rou_k[i] * u_k[i] * u_k[i] + p_k[i]) - D2[i];
			F3[i] = 0.5*( E[i + 1] * u_k[i + 1] + p_k[i + 1] * u_k[i + 1] + E[i] * u_k[i] + p_k[i] * u_k[i]) - D3[i];
		}
		for (i = 1; i < 999; i++)
		{
			W1[i] = rou_k[i] + 0.5*0.1*(F1[i - 1] - F1[i]);
			W2[i] = rou_k[i] * u_k[i] + 0.5*0.1*(F2[i - 1] - F2[i]);
			W3[i] = E[i] + 0.5*0.1*(F3[i - 1] - F3[i]);
			//至此，得到新时刻下的rou,u,p
			rou[i] = W1[i];
			u[i] = W2[i] / W1[i];
			p[i] = 0.4*(W3[i]  - 0.5*W2[i] * W2[i] / W1[i]);
		}
	}

	for (i = 0; i < 1000; i++)
	{
		space[i] = -0.4995 + 0.001*i;
		fprintf(fp, "%.5f %.5f %.5f %.5f\n", space[i],rou[i],u[i],p[i]);
	}
}
double max(double a, double b)
{
	double max = a;
	if (b>a)
	{
		max = b;
	}
	return max;
}
