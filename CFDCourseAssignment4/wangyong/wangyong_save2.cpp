#include <stdio.h>
#include <math.h>
#define e 2.718281828459
#define a 2
#define gama 1.4
#define R 287.06
#define cfl 0.8
FILE *fp, *fq, *fr, *fs, *fg;
double nodes[331][71][2], nodesc1[331][71][2], nodes2[261][81][2], nodesc2[261][81][2];
int i, j, k, m, maxi, maxj, max;
int maxi2, maxj2, max2;
double dyi_1[331][71], dxi_1[331][71], dyj_1[331][71], dxj_1[331][71];

double sli_1[330][71], slj_1[331][71], area_1[330][71];
double pho1_1[331][71], pre1_1[331][71], vx1_1[331][71], vy1_1[331][71], T1_1[331][71], ma1_1[331][71];
double pho2_1[261][81], pre2_1[261][81], vx2_1[261][81], vy2_1[261][81], T2_1[261][81], ma2_1[261][81];

double total_pre1_1[331][71], total_T1_1[331][71];
double total_pre2_1[261][81], total_T2_1[261][81];
double Uir_1[331][71], Uil_1[331][71], Ujr_1[331][71], Ujl_1[331][71], H_1[331][71], Fjr_1[331][71][4], Fjl_1[331][71][4], Fir_1[331][71][4], Fil_1[331][71][4];
double shengsu_1, pho_1_ba, vx1_1_ba, vy1_1_ba, H_1_ba, U_1_ba, lanmeta1_1, lanmeta2_1, lanmeta3_1;
double beta1_1, beta2_1, beta3_1, beta4_1, beta5_1, beta6_1, beta7_1;
double AQi_1[331][71][4], AQj_1[331][71][4], Flux_1[331][71][4], Q_1[331][71][4];
double resm1, resave[4], imax, jmax, tyj, txj, tyi, txi, tsli, tslj, vi, vj, sonic, chvel, dt, tres1, real[331];
double vnorm, vtemp, maxflux, maxflux2, maxflux3, maxflux4;
    typedef struct XY
    {
        double x,y;
    }XY;
void mesh_generation()
{
/*
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
*/

    int cellBegin=0;
    int maxI=330,maxJ=70,block1=(int)(0.2*maxI+0.1),block2=(int)(0.2*maxI+0.1);
    const double dx=1.0/maxI;
    double dy=0.3/maxJ;//后面会变

    XY mesh[331][71];
    for(int i=cellBegin; i<=block1; i++)
    {
        for(int j=cellBegin; j<=maxJ; j++)
        {
            mesh[i][j].x=(i-1)*dx;//x=(i-1)*dx
            mesh[i][j].y=(j-1)*dy;//y=(j-1)*dy下同
        }
    }
    
    for(int i=block1+1; i<=block1+block2; i++)
    {
        double h=0.25*(i-1)*dx-0.05;
        dy=(0.3-h)/maxJ;
        for(int j=cellBegin;j<=maxJ; j++)
        {
            mesh[i][j].x=(i-1)*dx;
            mesh[i][j].y=h+(j-1)*dy;
        }
    }

    dy=(0.3-0.05)/maxJ;
    for(int i=block1+block2+1; i<=maxI; i++)
    {
        for(int j=cellBegin;j<=maxJ; j++)
        {
            mesh[i][j].x=(i-1)*dx;
            mesh[i][j].y=0.05+(j-1)*dy;
        }
    }
    for(int i=cellBegin;i<=maxI;i++)
    {
        for(j=cellBegin;j<=maxJ;j++)
        {
            nodes[i][j][0]=mesh[i][j].x;
            nodes[i][j][1]=mesh[i][j].y;
        }
    }
}
void initialize()
{
    /*
	for (j = 0; j < 71; j++)
	{
		for (i = 0; i < 331; i++)
		{
			vx1_1[i][j] = 400; //3
			vy1_1[i][j] = 0;
			pre1_1[i][j] = 101325;
			T1_1[i][j] = 230;
			pho1_1[i][j] = pre1_1[i][j] / T1_1[i][j] / R;
			ma1_1[i][j] = sqrt(vx1_1[i][j] * vx1_1[i][j] + vy1_1[i][j] * vy1_1[i][j]) / sqrt(gama * R * T1_1[i][j]);
		}
	}
	for (j = 0; j < 81; j++)
	{
		for (i = 0; i < 261; i++)
		{
			vx2_1[i][j] = 150;
			vy2_1[i][j] = 0;
			pre2_1[i][j] = 90000;
			T2_1[i][j] = 288.2;
			pho2_1[i][j] = pre2_1[i][j] / T2_1[i][j] / R;
			ma2_1[i][j] = sqrt(vx2_1[i][j] * vx2_1[i][j] + vy2_1[i][j] * vy2_1[i][j]) / sqrt(gama * R * T2_1[i][j]);
		}
	}
    */
   	for (j = 0; j < 71; j++)
	{
		for (i = 0; i < 331; i++)
		{
			vx1_1[i][j] = 625; //3
			vy1_1[i][j] = 0;
			pre1_1[i][j] = 101325;
			T1_1[i][j] = 300;
			pho1_1[i][j] = pre1_1[i][j] / T1_1[i][j] / R;
			ma1_1[i][j] = sqrt(vx1_1[i][j] * vx1_1[i][j] + vy1_1[i][j] * vy1_1[i][j]) / sqrt(gama * R * T1_1[i][j]);
		}
	}
	for (j = 0; j < 81; j++)
	{
		for (i = 0; i < 261; i++)
		{
			vx2_1[i][j] = 150;
			vy2_1[i][j] = 0;
			pre2_1[i][j] = 90000;
			T2_1[i][j] = 288.2;
			pho2_1[i][j] = pre2_1[i][j] / T2_1[i][j] / R;
			ma2_1[i][j] = sqrt(vx2_1[i][j] * vx2_1[i][j] + vy2_1[i][j] * vy2_1[i][j]) / sqrt(gama * R * T2_1[i][j]);
		}
	}

}
void area_caculate()
{
	for (j = 0; j < 71; j++)
	{
		for (i = 0; i < 330; i++)
		{
			dyi_1[i][j] = nodes[i][j][1] - nodes[i + 1][j][1];
			dxi_1[i][j] = nodes[i + 1][j][0] - nodes[i][j][0];
			sli_1[i][j] = sqrt(dyi_1[i][j] * dyi_1[i][j] + dxi_1[i][j] * dxi_1[i][j]);
		}
	}
	for (j = 0; j < 70; j++)
	{
		for (i = 0; i < 331; i++)
		{
			dyj_1[i][j] = nodes[i][j + 1][1] - nodes[i][j][1];
			dxj_1[i][j] = nodes[i][j][0] - nodes[i][j + 1][0];
			slj_1[i][j] = sqrt(dyj_1[i][j] * dyj_1[i][j] + dxj_1[i][j] * dxj_1[i][j]);
		}
	}
	for (j = 0; j < 70; j++)
	{
		for (i = 0; i < 330; i++)
		{
			area_1[i][j] = (dyj_1[i][j] + dyj_1[i + 1][j]) * dxi_1[i][j] / 2;
		}
	} ///////////////////////////////


}

double tao(double Ma)
{
	double taoma;
	taoma = 1 / (1 + (gama - 1) * Ma * Ma / 2);
	return taoma;
}
double pai(double Ma)
{
	double paima;
	paima = pow(1 / (1 + (gama - 1) * Ma * Ma / 2), gama / (gama - 1));
	return paima;
}
void boundary_conditions()
{
	//	printf("%d\n",k);
	//	 fp=fopen("1.txt","w");
    /*
	for (j = 1; j < 70; j++) //进口边界条件 分区1 亚音进口，一个变量外推
	{

		total_pre1_1[0][j] = 250000;
		total_T1_1[0][j] = 330;
		vy1_1[0][j] = 0;
		if (T1_1[1][j] > 330)
		{
			T1_1[0][j] = 330;
		}
		else
		{
			T1_1[0][j] = T1_1[1][j];
		}
		//静温外推
		ma1_1[0][j] = sqrt((total_T1_1[0][j] / T1_1[0][j] - 1) * 2 / (gama - 1)); //壁面速度为0
		pre1_1[0][j] = total_pre1_1[0][j] * pai(ma1_1[0][j]);
		pho1_1[0][j] = pre1_1[0][j] / R / T1_1[0][j];
		vx1_1[0][j] = ma1_1[0][j] * sqrt(gama * R * T1_1[0][j]);
		//   fprintf(fp,"%.10f    %.10f    %.10f    %.10f    %.10f    %.10f\n",vx1_1[330][j],pho1_1[330][j],pre1_1[330][j],ma1_1[330][j],vy1_1[330][j],T1_1[330][j]);
	}

	for (j = 1; j < 70; j++)
	{
		if (ma1_1[329][j] >= 1) //超音全部外推
		{
			vx1_1[330][j] = vx1_1[329][j];
			vy1_1[330][j] = vy1_1[329][j];
			pre1_1[330][j] = pre1_1[329][j];
			T1_1[330][j] = T1_1[329][j];
			pho1_1[330][j] = pre1_1[330][j] / R / T1_1[330][j];
			ma1_1[330][j] = sqrt(vx1_1[330][j] * vx1_1[330][j] + vy1_1[330][j] * vy1_1[330][j]) / sqrt(gama * R * T1_1[330][j]);
		}
		else //亚音3个外推   静压给定
		{
			pre1_1[330][j] = 85419;
			vx1_1[330][j] = vx1_1[329][j];
			vy1_1[330][j] = vy1_1[329][j];
			T1_1[330][j] = T1_1[329][j];
			pho1_1[330][j] = pre1_1[330][j] / R / T1_1[330][j];
			ma1_1[330][j] = sqrt(vx1_1[330][j] * vx1_1[330][j] + vy1_1[330][j] * vy1_1[330][j]) / sqrt(gama * R * T1_1[330][j]);
		}
	}

	for (i = 1; i < 151; i++) //喷管壁面边界条件
	{
		//vnorm=vx1_1[i][69]*dyi_1[i][70]+vy1_1[i][69]*dxi_1[i][70];
		//	vtemp=2*vnorm/sli_1[i][70]/sli_1[i][70];
		vnorm = sqrt(vy1_1[i][69] * vy1_1[i][69] + vx1_1[i][69] * vx1_1[i][69]);
		vtemp = -(dyi_1[i][70]) / dxi_1[i][70];
		vtemp = atan(vtemp);
		pho1_1[i][70] = pho1_1[i][69];
		vx1_1[i][70] = -vx1_1[i][69] + 2 * cos(vtemp) * vnorm;
		vy1_1[i][70] = -vy1_1[i][69] + 2 * sin(vtemp) * vnorm;
		//	vx1_1[i][70]=vx1_1[i][69]-vtemp*dyi_1[i][69];    //外推
		//	vy1_1[i][70]=vy1_1[i][69]-vtemp*dxi_1[i][69];
		pre1_1[i][70] = pre1_1[i][69];
		T1_1[i][70] = pre1_1[i][70] / pho1_1[i][70] / R;
		ma1_1[i][70] = sqrt(vx1_1[i][70] * vx1_1[i][70] + vy1_1[i][70] * vy1_1[i][70]) / sqrt(gama * R * T1_1[i][70]);
	}

	for (i = 1; i < 330; i++) //对称边界条件
	{
		pho1_1[i][0] = pho1_1[i][1];
		vx1_1[i][0] = vx1_1[i][1];
		vy1_1[i][0] = 0;
		pre1_1[i][0] = pre1_1[i][1];
		T1_1[i][0] = pre1_1[i][0] / pho1_1[i][0] / R;
		ma1_1[i][0] = sqrt(vx1_1[i][0] * vx1_1[i][0] + vy1_1[i][0] * vy1_1[i][0]) / sqrt(gama * R * T1_1[i][0]);
		//	fprintf(fp,"%.10f    %.10f    %.10f    %.10f    %.10f    %.10f\n",vx1_1[i][0],pho1_1[i][0],pre1_1[i][0],ma1_1[i][0],vy1_1[i][0],T1_1[i][0]);
	}

	////////////////////////////
	for (j = 1; j < 80; j++) //进口边界条件 分区1 亚音进口，一个变量外推
	{

		total_pre2_1[0][j] = 101325;
		total_T2_1[0][j] = 288.2;
		ma2_1[0][j] = 0.5;
		//vy1_1[0][j]=0;
		pre2_1[0][j] = total_pre2_1[0][j] * pai(ma2_1[0][j]);
		T2_1[0][j] = total_T2_1[0][j] * tao(ma2_1[0][j]);
		pho2_1[0][j] = pre2_1[0][j] / T2_1[0][j] / R;
		vy2_1[0][j] = vy2_1[1][j];
		vx2_1[0][j] = sqrt(ma2_1[0][j] * sqrt(gama * R * T2_1[0][j]) * ma2_1[0][j] * sqrt(gama * R * T2_1[0][j]) - vy2_1[0][j] * vy2_1[0][j]);

		//   fprintf(fp,"%.10f    %.10f    %.10f    %.10f    %.10f    %.10f\n",vx1_1[330][j],pho1_1[330][j],pre1_1[330][j],ma1_1[330][j],vy1_1[330][j],T1_1[330][j]);
	}

	for (j = 1; j < 80; j++) //////出口边界
	{
		if (ma2_1[259][j] >= 1) //超音全部外推
		{
			vx2_1[260][j] = vx2_1[259][j];
			vy2_1[260][j] = vy2_1[259][j];
			pre2_1[260][j] = pre2_1[259][j];
			T2_1[260][j] = T2_1[259][j];
			pho2_1[260][j] = pre2_1[260][j] / R / T2_1[260][j];
			ma2_1[260][j] = sqrt(vx2_1[260][j] * vx2_1[260][j] + vy2_1[260][j] * vy2_1[260][j]) / sqrt(gama * R * T2_1[260][j]);
		}
		else //亚音3个外推   静压给定
		{
			pre2_1[260][j] = 85419;
			vx2_1[260][j] = vx2_1[259][j];
			vy2_1[260][j] = vy2_1[259][j];
			T2_1[260][j] = T2_1[259][j];
			pho2_1[260][j] = pre2_1[260][j] / R / T2_1[260][j];
			ma2_1[260][j] = sqrt(vx2_1[260][j] * vx2_1[260][j] + vy2_1[260][j] * vy2_1[260][j]) / sqrt(gama * R * T2_1[260][j]);
		}
	}

	for (i = 1; i < 260; i++)  //喷管壁面边界条件
	{						   
			 
		if (vy2_1[i][79] >= 0) //出口，给一推三静压
		{
			pre2_1[i][80] = 85419;
			vx2_1[i][80] = vx2_1[i][79];
			vy2_1[i][80] = vy2_1[i][79];
			T2_1[i][80] = T2_1[i][79];
			pho2_1[i][80] = pre2_1[i][80] / R / T2_1[i][80];
			ma2_1[i][80] = sqrt(vx2_1[i][80] * vx2_1[i][80] + vy2_1[i][80] * vy2_1[i][80]) / sqrt(gama * R * T2_1[i][80]);
		}
		else //进口 ，给三推一
		{
			pre2_1[i][80] = 85419;
			ma2_1[i][80] = 0.5;
			T2_1[i][80] = 288.2 * tao(ma2_1[i][80]);
			vy2_1[i][80] = vy2_1[i][79];
			pho2_1[i][80] = pre2_1[i][80] / R / T2_1[i][80];
			vx2_1[i][80] = sqrt(ma2_1[i][80] * sqrt(gama * R * T2_1[i][80]) * ma2_1[i][80] * sqrt(gama * R * T2_1[i][80]) - vy2_1[i][80] * vy2_1[i][80]);
		}
	}
    
	for (i = 1; i < 81; i++) //下壁面条件
	{
		vnorm = sqrt(vy2_1[i][1] * vy2_1[i][1] + vx2_1[i][1] * vx2_1[i][1]);
		vtemp = -(dyi_1[i][0]) / dxi_1[i][0];
		vtemp = atan(vtemp);
		pho2_1[i][0] = pho2_1[i][1];
		vx2_1[i][0] = -vx2_1[i][1] + 2 * cos(vtemp) * vnorm;
		vy2_1[i][0] = -vy2_1[i][1] + 2 * sin(vtemp) * vnorm;
		//	vx2_1[i][80]=vx2_1[i][79]-vtemp*dyi_1[i][79];    //外推
		//	vy2_1[i][80]=vy2_1[i][79]-vtemp*dxi_1[i][79];
		pre2_1[i][0] = pre2_1[i][1];
		T2_1[i][0] = pre2_1[i][0] / pho2_1[i][0] / R;
		ma2_1[i][0] = sqrt(vx2_1[i][0] * vx2_1[i][0] + vy2_1[i][0] * vy2_1[i][0]) / sqrt(gama * R * T2_1[i][0]);

	}

	for (i = 151; i < 330; i++)
	{
		pho1_1[i][70] = pho2_1[i - 70][1];
		pre1_1[i][70] = pre2_1[i - 70][1];
		vx1_1[i][70] = vx2_1[i - 70][1];
		vy1_1[i][70] = vy2_1[i - 70][1];
		T1_1[i][70] = pre1_1[i][70] / pho1_1[i][70] / R;
		ma1_1[i][70] = sqrt(vx1_1[i][70] * vx1_1[i][70] + vy1_1[i][70] * vy1_1[i][70]) / sqrt(gama * R * T1_1[i][70]);

		pho2_1[i - 70][0] = pho1_1[i][69];
		pre2_1[i - 70][0] = pre1_1[i][69];
		vx2_1[i - 70][0] = vx1_1[i][69];
		vy2_1[i - 70][0] = vy1_1[i][69];
		T2_1[i - 70][0] = pre2_1[i - 70][0] / pho2_1[i - 70][0] / R;
		ma2_1[i - 70][0] = sqrt(vx2_1[i - 70][0] * vx2_1[i - 70][0] + vy2_1[i - 70][0] * vy2_1[i - 70][0]) / sqrt(gama * R * T2_1[i - 70][0]);
	
	}
    */
 
    //XY N1[331][71],N3[331][71];
    double cellBegin=1,cellIEnd=330,cellJEnd=70;
    const double p0=101325;
const double u0=624.9397;
const double v0=0;
const double T0=300;
const double rho0=1.176829268;
const double H0=496662.5;
const double E0=483117.6;
    double nx,ny;

    
    double pw;
	//上边
    const unsigned j=cellJEnd;
    for(unsigned i=cellBegin; i<=cellIEnd; i++)
    {
        //pw=0.5*(3*p[i][j]-p[i][j-1]);//壁面的压力用两点外推
        pw=pre1_1[i][j];
        //nx=N3[i][j].x; ny=N3[i][j].y;
        	nx=dyi_1[i][j] / sli_1[i][j];
			ny=dxi_1[i][j] / sli_1[i][j];
        Fil_1[i][j+1][0]=0;       
        Fil_1[i][j+1][1]=pw*nx;
        Fil_1[i][j+1][2]=pw*ny;
        Fil_1[i][j+1][3]=0;
	}
	//下边
    for (unsigned i = cellBegin; i <= cellIEnd; i++)
    {
        //nx=N1[i][1].x; ny=N1[i][1].y;
        //pw = 0.5 * (3 * p[i][j] - p[i][j+1]); //壁面的压力用两点外推
        pw=pre1_1[i][1];
        nx=dyj_1[i][j] / slj_1[i][j];
		ny=dxj_1[i][j] / slj_1[i][j];
        Fil_1[i][1][0] = 0;     //注意,虚网格nx ny使用相邻网格的值
        Fil_1[i][1][1] = pw*nx;
        Fil_1[i][1][2] = pw*ny;
        Fil_1[i][1][3] = 0;

      pho1_1[i][0] = pho1_1[i][1];
        vx1_1[i][0] =   vx1_1[i][1];
        vy1_1[i][0] =   vy1_1[i][1];
        pre1_1[i][0] =   pre1_1[i][1];
        H_1[i][0] =   H_1[i][1];
    }
	//左边(右边不用管,全部推出)
    for(unsigned j=cellBegin; j<=cellJEnd; j++)
    {
        pre1_1[0][j] = 101325;
        vx1_1[0][j] = 624.9397;
        vy1_1[0][j] = 0;
      pho1_1[0][j] = 1.176829;
        H_1[0][j] = 496662.5; //Cp*300+0.5*625*625

        const double nx=1, ny=0;
        const double Vcv0=624.9397;
        Fjl_1[0][j][0]=rho0*Vcv0;
        Fjl_1[0][j][1]=rho0*u0*Vcv0+p0*nx;
        Fjl_1[0][j][2]=rho0*v0*Vcv0+p0*ny;
        Fjl_1[0][j][3]=rho0*H0*Vcv0;
    }
    
}
void roe() //利用roe格式求解
{
	for (j = 1; j < 70; j++)
	{
		for (i = 1; i < 331; i++)
		{
			Uir_1[i][j] = dyj_1[i][j] / slj_1[i][j] * vx1_1[i][j] + dxj_1[i][j] / slj_1[i][j] * vy1_1[i][j]; //沿流向流动通过边界的速度
			Uil_1[i][j] = dyj_1[i][j] / slj_1[i][j] * vx1_1[i - 1][j] + dxj_1[i][j] / slj_1[i][j] * vy1_1[i - 1][j];
		}
	}
	for (j = 1; j < 71; j++)
	{
		for (i = 1; i < 330; i++)
		{
			Ujr_1[i][j] = dyi_1[i][j] / sli_1[i][j] * vx1_1[i][j] + dxi_1[i][j] / sli_1[i][j] * vy1_1[i][j]; //沿垂直流向流动通过边界的速度
			Ujl_1[i][j] = dyi_1[i][j] / sli_1[i][j] * vx1_1[i][j - 1] + dxi_1[i][j] / sli_1[i][j] * vy1_1[i][j - 1];
		}
	}
	for (j = 0; j < 71; j++)
	{
		for (i = 0; i < 331; i++)
		{
			H_1[i][j] = (pre1_1[i][j] / (gama - 1) + 0.5 * pho1_1[i][j] * (pow(vy1_1[i][j], 2) + pow(vx1_1[i][j], 2)) + pre1_1[i][j]) / pho1_1[i][j];
		}
	} ///////////////////////


	//I方向的F求解
	for (j = 1; j < 70; j++)
	{
		for (i = 1; i < 331; i++)
		{
			Fil_1[i][j][0] = slj_1[i][j] * pho1_1[i - 1][j] * Uil_1[i][j]; //左侧对流通量
			Fil_1[i][j][1] = slj_1[i][j] * pho1_1[i - 1][j] * vx1_1[i - 1][j] * Uil_1[i][j] + dyj_1[i][j] * pre1_1[i - 1][j];
			Fil_1[i][j][2] = slj_1[i][j] * pho1_1[i - 1][j] * vy1_1[i - 1][j] * Uil_1[i][j] + dxj_1[i][j] * pre1_1[i - 1][j];
			Fil_1[i][j][3] = slj_1[i][j] * pho1_1[i - 1][j] * H_1[i - 1][j] * Uil_1[i][j];

			Fir_1[i][j][0] = slj_1[i][j] * pho1_1[i][j] * Uir_1[i][j]; //右侧对流通量
			Fir_1[i][j][1] = slj_1[i][j] * pho1_1[i][j] * vx1_1[i][j] * Uir_1[i][j] + dyj_1[i][j] * pre1_1[i][j];
			Fir_1[i][j][2] = slj_1[i][j] * pho1_1[i][j] * vy1_1[i][j] * Uir_1[i][j] + dxj_1[i][j] * pre1_1[i][j];
			Fir_1[i][j][3] = slj_1[i][j] * pho1_1[i][j] * H_1[i][j] * Uir_1[i][j];
		}
	} ///////////////



	//J方向的F求解
	for (j = 1; j < 71; j++)
	{
		for (i = 1; i < 330; i++)
		{
			Fjl_1[i][j][0] = sli_1[i][j] * pho1_1[i][j - 1] * Ujl_1[i][j]; //左侧通量
			Fjl_1[i][j][1] = sli_1[i][j] * pho1_1[i][j - 1] * vx1_1[i][j - 1] * Ujl_1[i][j] + dyi_1[i][j] * pre1_1[i][j - 1];
			Fjl_1[i][j][2] = sli_1[i][j] * pho1_1[i][j - 1] * vy1_1[i][j - 1] * Ujl_1[i][j] + dxi_1[i][j] * pre1_1[i][j - 1];
			Fjl_1[i][j][3] = sli_1[i][j] * pho1_1[i][j - 1] * H_1[i][j - 1] * Ujl_1[i][j];

			Fjr_1[i][j][0] = sli_1[i][j] * pho1_1[i][j] * Ujr_1[i][j]; //右侧通量
			Fjr_1[i][j][1] = sli_1[i][j] * pho1_1[i][j] * vx1_1[i][j] * Ujr_1[i][j] + dyi_1[i][j] * pre1_1[i][j];
			Fjr_1[i][j][2] = sli_1[i][j] * pho1_1[i][j] * vy1_1[i][j] * Ujr_1[i][j] + dxi_1[i][j] * pre1_1[i][j];
			Fjr_1[i][j][3] = sli_1[i][j] * pho1_1[i][j] * H_1[i][j] * Ujr_1[i][j];
		}
	} ///////////



	//计算beta与AQ


	for (j = 1; j < 70; j++)
	{
		for (i = 1; i < 331; i++)
		{
			pho_1_ba = sqrt(pho1_1[i - 1][j] * pho1_1[i][j]); //Roe平均定义
			vx1_1_ba = (vx1_1[i - 1][j] * sqrt(pho1_1[i - 1][j]) + vx1_1[i][j] * sqrt(pho1_1[i][j])) / (sqrt(pho1_1[i - 1][j]) + sqrt(pho1_1[i][j]));
			vy1_1_ba = (vy1_1[i - 1][j] * sqrt(pho1_1[i - 1][j]) + vy1_1[i][j] * sqrt(pho1_1[i][j])) / (sqrt(pho1_1[i - 1][j]) + sqrt(pho1_1[i][j]));
			H_1_ba = (H_1[i - 1][j] * sqrt(pho1_1[i - 1][j]) + H_1[i][j] * sqrt(pho1_1[i][j])) / (sqrt(pho1_1[i - 1][j]) + sqrt(pho1_1[i][j]));
			U_1_ba = (Uil_1[i][j] * sqrt(pho1_1[i - 1][j]) + Uir_1[i][j] * sqrt(pho1_1[i][j])) / (sqrt(pho1_1[i - 1][j]) + sqrt(pho1_1[i][j]));
			shengsu_1 = sqrt((gama - 1) * (H_1_ba - (vx1_1_ba * vx1_1_ba + vy1_1_ba * vy1_1_ba) / 2));
			lanmeta1_1 = U_1_ba;
			lanmeta2_1 = U_1_ba - shengsu_1;
			lanmeta3_1 = U_1_ba + shengsu_1;

			if (fabs(lanmeta1_1) >= 0.1) //熵修正
			{
				lanmeta1_1 = fabs(lanmeta1_1);
			}
			else
			{
				lanmeta1_1 = (lanmeta1_1 * lanmeta1_1 + 0.01) / 0.2;
			}
			if (fabs(lanmeta2_1) >= 0.1)
			{
				lanmeta2_1 = fabs(lanmeta2_1);
			}
			else
			{
				lanmeta2_1 = (lanmeta2_1 * lanmeta2_1 + 0.01) / 0.2;
			}
			if (fabs(lanmeta3_1) >= 0.1)
			{
				lanmeta3_1 = fabs(lanmeta3_1);
			}
			else
			{
				lanmeta3_1 = (lanmeta3_1 * lanmeta3_1 + 0.01) / 0.2;
			}

			beta1_1 = slj_1[i][j] * lanmeta1_1 * (pho1_1[i][j] - pho1_1[i - 1][j] - (pre1_1[i][j] - pre1_1[i - 1][j]) / shengsu_1 / shengsu_1);
			beta2_1 = slj_1[i][j] * lanmeta3_1 * (pre1_1[i][j] - pre1_1[i - 1][j] + pho_1_ba * shengsu_1 * (Uir_1[i][j] - Uil_1[i][j])) / (2 * shengsu_1 * shengsu_1);
			beta3_1 = slj_1[i][j] * lanmeta2_1 * (pre1_1[i][j] - pre1_1[i - 1][j] - pho_1_ba * shengsu_1 * (Uir_1[i][j] - Uil_1[i][j])) / (2 * shengsu_1 * shengsu_1);
			beta4_1 = beta1_1 + beta2_1 + beta3_1;
			beta5_1 = shengsu_1 * (beta2_1 - beta3_1);
			beta6_1 = slj_1[i][j] * lanmeta1_1 * pho_1_ba * (vx1_1[i][j] - vx1_1[i - 1][j] - dyj_1[i][j] / slj_1[i][j] * (Uir_1[i][j] - Uil_1[i][j]));
			beta7_1 = slj_1[i][j] * lanmeta1_1 * pho_1_ba * (vy1_1[i][j] - vy1_1[i - 1][j] - dxj_1[i][j] / slj_1[i][j] * (Uir_1[i][j] - Uil_1[i][j]));

			AQi_1[i][j][0] = beta4_1;
			AQi_1[i][j][1] = vx1_1_ba * beta4_1 + dyj_1[i][j] / slj_1[i][j] * beta5_1 + beta6_1;
			AQi_1[i][j][2] = vy1_1_ba * beta4_1 + dxj_1[i][j] / slj_1[i][j] * beta5_1 + beta7_1;
			AQi_1[i][j][3] = H_1_ba * beta4_1 + U_1_ba * beta5_1 + vx1_1_ba * beta6_1 + vy1_1_ba * beta7_1 - shengsu_1 * shengsu_1 * beta1_1 / (gama - 1);
		}
	} ///////////////////////



	for (j = 1; j < 71; j++)
	{
		for (i = 1; i < 330; i++)
		{
			pho_1_ba = sqrt(pho1_1[i][j - 1] * pho1_1[i][j]);
			vx1_1_ba = (vx1_1[i][j - 1] * sqrt(pho1_1[i][j - 1]) + vx1_1[i][j] * sqrt(pho1_1[i][j])) / (sqrt(pho1_1[i][j - 1]) + sqrt(pho1_1[i][j]));
			vy1_1_ba = (vy1_1[i][j - 1] * sqrt(pho1_1[i][j - 1]) + vy1_1[i][j] * sqrt(pho1_1[i][j])) / (sqrt(pho1_1[i][j - 1]) + sqrt(pho1_1[i][j]));
			H_1_ba = (H_1[i][j - 1] * sqrt(pho1_1[i][j - 1]) + H_1[i][j] * sqrt(pho1_1[i][j])) / (sqrt(pho1_1[i][j - 1]) + sqrt(pho1_1[i][j]));
			U_1_ba = (Ujl_1[i][j] * sqrt(pho1_1[i][j - 1]) + Ujr_1[i][j] * sqrt(pho1_1[i][j])) / (sqrt(pho1_1[i][j - 1]) + sqrt(pho1_1[i][j]));
			shengsu_1 = sqrt((gama - 1) * (H_1_ba - (vx1_1_ba * vx1_1_ba + vy1_1_ba * vy1_1_ba) / 2));
			lanmeta1_1 = U_1_ba;
			lanmeta2_1 = U_1_ba - shengsu_1;
			lanmeta3_1 = U_1_ba + shengsu_1;
			if (fabs(lanmeta1_1) >= 0.1)
			{
				lanmeta1_1 = fabs(lanmeta1_1);
			}
			else
			{
				lanmeta1_1 = (lanmeta1_1 * lanmeta1_1 + 0.01) / 0.2;
			}
			if (fabs(lanmeta2_1) >= 0.1)
			{
				lanmeta2_1 = fabs(lanmeta2_1);
			}
			else
			{
				lanmeta2_1 = (lanmeta2_1 * lanmeta2_1 + 0.01) / 0.2;
			}
			if (fabs(lanmeta3_1) >= 0.1)
			{
				lanmeta3_1 = fabs(lanmeta3_1);
			}
			else
			{
				lanmeta3_1 = (lanmeta3_1 * lanmeta3_1 + 0.01) / 0.2;
			}

			beta1_1 = sli_1[i][j] * lanmeta1_1 * (pho1_1[i][j] - pho1_1[i][j - 1] - (pre1_1[i][j] - pre1_1[i][j - 1]) / shengsu_1 / shengsu_1);
			beta2_1 = sli_1[i][j] * lanmeta3_1 * (pre1_1[i][j] - pre1_1[i][j - 1] + pho_1_ba * shengsu_1 * (Ujr_1[i][j] - Ujl_1[i][j])) / (2 * shengsu_1 * shengsu_1);
			beta3_1 = sli_1[i][j] * lanmeta2_1 * (pre1_1[i][j] - pre1_1[i][j - 1] - pho_1_ba * shengsu_1 * (Ujr_1[i][j] - Ujl_1[i][j])) / (2 * shengsu_1 * shengsu_1);
			beta4_1 = beta1_1 + beta2_1 + beta3_1;
			beta5_1 = shengsu_1 * (beta2_1 - beta3_1);
			beta6_1 = sli_1[i][j] * lanmeta1_1 * pho_1_ba * (vx1_1[i][j] - vx1_1[i][j - 1] - dyi_1[i][j] / sli_1[i][j] * (Ujr_1[i][j] - Ujl_1[i][j]));
			beta7_1 = sli_1[i][j] * lanmeta1_1 * pho_1_ba * (vy1_1[i][j] - vy1_1[i][j - 1] - dxi_1[i][j] / sli_1[i][j] * (Ujr_1[i][j] - Ujl_1[i][j]));

			//	fprintf(fp,"%.10f    %.10f    %.10f    %.10f    %.10f    %.10f    %.10f\n",beta1_1,beta2_1,beta3_1,beta4_1,beta5_1,beta6_1,beta7_1);
			AQj_1[i][j][0] = beta4_1;
			AQj_1[i][j][1] = vx1_1_ba * beta4_1 + dyi_1[i][j] / sli_1[i][j] * beta5_1 + beta6_1;
			AQj_1[i][j][2] = vy1_1_ba * beta4_1 + dxi_1[i][j] / sli_1[i][j] * beta5_1 + beta7_1;
			AQj_1[i][j][3] = H_1_ba * beta4_1 + U_1_ba * beta5_1 + vx1_1_ba * beta6_1 + vy1_1_ba * beta7_1 - shengsu_1 * shengsu_1 * beta1_1 / (gama - 1);
		}
	} ////////////////////



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
	for (j = 1; j < 70; j++)
	{
		for (i = 1; i < 330; i++)
		{
			Flux_1[i][j][0] = -(Fil_1[i][j][0] + Fir_1[i][j][0] - AQi_1[i][j][0]) / 2 + (Fil_1[i + 1][j][0] + Fir_1[i + 1][j][0] - AQi_1[i + 1][j][0]) / 2 - (Fjl_1[i][j][0] + Fjr_1[i][j][0] - AQj_1[i][j][0]) / 2 + (Fjl_1[i][j + 1][0] + Fjr_1[i][j + 1][0] - AQj_1[i][j + 1][0]) / 2;
			Flux_1[i][j][1] = -(Fil_1[i][j][1] + Fir_1[i][j][1] - AQi_1[i][j][1]) / 2 + (Fil_1[i + 1][j][1] + Fir_1[i + 1][j][1] - AQi_1[i + 1][j][1]) / 2 - (Fjl_1[i][j][1] + Fjr_1[i][j][1] - AQj_1[i][j][1]) / 2 + (Fjl_1[i][j + 1][1] + Fjr_1[i][j + 1][1] - AQj_1[i][j + 1][1]) / 2;
			Flux_1[i][j][2] = -(Fil_1[i][j][2] + Fir_1[i][j][2] - AQi_1[i][j][2]) / 2 + (Fil_1[i + 1][j][2] + Fir_1[i + 1][j][2] - AQi_1[i + 1][j][2]) / 2 - (Fjl_1[i][j][2] + Fjr_1[i][j][2] - AQj_1[i][j][2]) / 2 + (Fjl_1[i][j + 1][2] + Fjr_1[i][j + 1][2] - AQj_1[i][j + 1][2]) / 2;
			Flux_1[i][j][3] = -(Fil_1[i][j][3] + Fir_1[i][j][3] - AQi_1[i][j][3]) / 2 + (Fil_1[i + 1][j][3] + Fir_1[i + 1][j][3] - AQi_1[i + 1][j][3]) / 2 - (Fjl_1[i][j][3] + Fjr_1[i][j][3] - AQj_1[i][j][3]) / 2 + (Fjl_1[i][j + 1][3] + Fjr_1[i][j + 1][3] - AQj_1[i][j + 1][3]) / 2;
			if (Flux_1[i][j][0] > maxflux)
			{
				maxflux = Flux_1[i][j][0];
			}
			if (Flux_1[i][j][1] > maxflux2)
			{
				maxflux2 = Flux_1[i][j][1];
			}
			if (Flux_1[i][j][2] > maxflux3)
			{
				maxflux3 = Flux_1[i][j][2];
			}
			if (Flux_1[i][j][3] > maxflux4)
			{
				maxflux4 = Flux_1[i][j][3];
			}
		}
	} //////////

	

}


void iteration()
{
	for (j = 1; j < 70; j++)
	{
		for (i = 1; i < 330; i++)
		{
			Q_1[i][j][0] = pho1_1[i][j];
			Q_1[i][j][1] = pho1_1[i][j] * vx1_1[i][j];
			Q_1[i][j][2] = pho1_1[i][j] * vy1_1[i][j];
			Q_1[i][j][3] = pre1_1[i][j] / (gama - 1) + 0.5 * pho1_1[i][j] * (vx1_1[i][j] * vx1_1[i][j] + vy1_1[i][j] * vy1_1[i][j]);
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
	//	fp=fopen("1.txt","w");

	for (j = 1; j < 70; j++)
	{
		for (i = 1; i < 330; i++)
		{
			tyj = 0.5 * (dyj_1[i][j] + dyj_1[i + 1][j]);
			txj = 0.5 * (dxj_1[i][j] + dxj_1[i + 1][j]);
			tslj = 0.5 * (slj_1[i][j] + slj_1[i + 1][j]);
			tyi = 0.5 * (dyi_1[i][j] + dyi_1[i][j + 1]);
			txi = 0.5 * (dxi_1[i][j] + dxi_1[i][j + 1]);
			tsli = 0.5 * (sli_1[i][j] + sli_1[i][j + 1]);
			vj = tyj * vx1_1[i][j] + txj * vy1_1[i][j];
			vi = tyi * vx1_1[i][j] + txi * vy1_1[i][j];
			sonic = sqrt(gama * pre1_1[i][j] / pho1_1[i][j]);
			chvel = fabs(vj) + fabs(vi) + sonic * (tslj + tsli);
			dt = cfl / chvel;
			//	fprintf(fp,"%.5f\n",dt);
			Q_1[i][j][0] = Q_1[i][j][0] - dt * Flux_1[i][j][0]; //迭代求解下一步Q
			Q_1[i][j][1] = Q_1[i][j][1] - dt * Flux_1[i][j][1];
			Q_1[i][j][2] = Q_1[i][j][2] - dt * Flux_1[i][j][2];
			Q_1[i][j][3] = Q_1[i][j][3] - dt * Flux_1[i][j][3];
			pho1_1[i][j] = Q_1[i][j][0];
			vx1_1[i][j] = Q_1[i][j][1] / Q_1[i][j][0];
			vy1_1[i][j] = Q_1[i][j][2] / Q_1[i][j][0];
			pre1_1[i][j] = (gama - 1) * (Q_1[i][j][3] - 0.5 * pho1_1[i][j] * (vx1_1[i][j] * vx1_1[i][j] + vy1_1[i][j] * vy1_1[i][j]));
			T1_1[i][j] = pre1_1[i][j] / pho1_1[i][j] / R;
			ma1_1[i][j] = sqrt(vx1_1[i][j] * vx1_1[i][j] + vy1_1[i][j] * vy1_1[i][j]) / sqrt(gama * R * T1_1[i][j]);
		}
	} /////////////////

}
void record()
{
	for (j = 1; j < 70; j++)
	{
		for (i = 1; i < 330; i++)
		{
			nodesc1[i][j][0] = 0.25 * (nodes[i][j][0] + nodes[i + 1][j][0] + nodes[i + 1][j + 1][0] + nodes[i][j + 1][0]);
			nodesc1[i][j][1] = 0.25 * (nodes[i][j][1] + nodes[i + 1][j][1] + nodes[i + 1][j + 1][1] + nodes[i][j + 1][1]);
		}
	}
	for (j = 1; j < 70; j++)
	{
		nodesc1[0][j][0] = 0.5 * (nodes[1][j][0] + nodes[1][j + 1][0]);
		nodesc1[0][j][1] = 0.5 * (nodes[1][j][1] + nodes[1][j + 1][1]);
		nodesc1[330][j][0] = 0.5 * (nodes[330][j][0] + nodes[330][j + 1][0]);
		nodesc1[330][j][1] = 0.5 * (nodes[330][j][1] + nodes[330][j + 1][1]);
	}
	for (i = 1; i < 330; i++)
	{
		nodesc1[i][0][0] = 0.5 * (nodes[i][0][0] + nodes[i + 1][0][0]);
		nodesc1[i][0][1] = 0.5 * (nodes[i][0][1] + nodes[i + 1][0][1]);
		nodesc1[i][70][0] = 0.5 * (nodes[i][70][0] + nodes[i + 1][70][0]);
		nodesc1[i][70][1] = 0.5 * (nodes[i][70][1] + nodes[i + 1][70][1]);
	}
	nodesc1[0][0][0] = nodes[1][1][0];
	nodesc1[0][0][1] = nodes[1][0][1];
	nodesc1[330][0][0] = nodes[330][1][0];
	nodesc1[330][0][1] = nodes[330][0][1];
	nodesc1[330][70][0] = nodes[330][70][0];
	nodesc1[330][70][1] = nodes[330][70][1];
	nodesc1[0][70][0] = nodes[1][70][0];
	nodesc1[0][70][1] = nodes[1][70][1];
	for (j = 1; j < 70; j++)
	{
		pho1_1[0][j] = 0.5 * (pho1_1[0][j] + pho1_1[1][j]);
		vx1_1[0][j] = 0.5 * (vx1_1[0][j] + vx1_1[1][j]);
		vy1_1[0][j] = 0.5 * (vy1_1[0][j] + vy1_1[1][j]);
		pre1_1[0][j] = 0.5 * (pre1_1[0][j] + pre1_1[1][j]);
		pho1_1[330][j] = 0.5 * (pho1_1[330][j] + pho1_1[259][j]);
		vx1_1[330][j] = 0.5 * (vx1_1[330][j] + vx1_1[259][j]);
		vy1_1[330][j] = 0.5 * (vy1_1[330][j] + vy1_1[259][j]);
		pre1_1[330][j] = 0.5 * (pre1_1[330][j] + pre1_1[259][j]);
	}
	for (i = 1; i < 330; i++)
	{
		pho1_1[i][0] = 0.5 * (pho1_1[i][0] + pho1_1[i][1]);
		vx1_1[i][0] = 0.5 * (vx1_1[i][0] + vx1_1[i][1]);
		vy1_1[i][0] = 0.5 * (vy1_1[i][0] + vy1_1[i][1]);
		pre1_1[i][0] = 0.5 * (pre1_1[i][0] + pre1_1[i][1]);
		pho1_1[i][70] = 0.5 * (pho1_1[i][70] + pho1_1[i][69]);
		vx1_1[i][70] = 0.5 * (vx1_1[i][70] + vx1_1[i][69]);
		vy1_1[i][70] = 0.5 * (vy1_1[i][70] + vy1_1[i][69]);
		pre1_1[i][70] = 0.5 * (pre1_1[i][70] + pre1_1[i][69]);
	}
	pho1_1[0][0] = pho1_1[1][0];
	vx1_1[0][0] = vx1_1[1][0];
	vy1_1[0][0] = vy1_1[1][0];
	pre1_1[0][0] = pre1_1[1][0];

	pho1_1[0][70] = pho1_1[1][70];
	vx1_1[0][70] = vx1_1[1][70];
	vy1_1[0][70] = vy1_1[1][70];
	pre1_1[0][70] = pre1_1[1][70];

	pho1_1[330][0] = pho1_1[259][0];
	vx1_1[330][0] = vx1_1[259][0];
	vy1_1[330][0] = vy1_1[259][0];
	pre1_1[330][0] = pre1_1[259][0];

	pho1_1[330][70] = pho1_1[259][70];
	vx1_1[330][70] = vx1_1[259][70];
	vy1_1[330][70] = vy1_1[259][70];
	pre1_1[330][70] = pre1_1[259][70];
	////////////////////////////////////////

	for (j = 1; j < 80; j++)
	{
		for (i = 1; i < 260; i++)
		{
			nodesc2[i][j][0] = 0.25 * (nodes2[i][j][0] + nodes2[i + 1][j][0] + nodes2[i + 1][j + 1][0] + nodes2[i][j + 1][0]);
			nodesc2[i][j][1] = 0.25 * (nodes2[i][j][1] + nodes2[i + 1][j][1] + nodes2[i + 1][j + 1][1] + nodes2[i][j + 1][1]);
		}
	}
	for (j = 1; j < 80; j++)
	{
		nodesc2[0][j][0] = 0.5 * (nodes2[1][j][0] + nodes2[1][j + 1][0]);
		nodesc2[0][j][1] = 0.5 * (nodes2[1][j][1] + nodes2[1][j + 1][1]);
		nodesc2[260][j][0] = 0.5 * (nodes2[260][j][0] + nodes2[260][j + 1][0]);
		nodesc2[260][j][1] = 0.5 * (nodes2[260][j][1] + nodes2[260][j + 1][1]);
	}
	for (i = 1; i < 260; i++)
	{
		nodesc2[i][0][0] = 0.5 * (nodes2[i][0][0] + nodes2[i + 1][0][0]);
		nodesc2[i][0][1] = 0.5 * (nodes2[i][0][1] + nodes2[i + 1][0][1]); /////////////
		nodesc2[i][80][0] = 0.5 * (nodes2[i][80][0] + nodes2[i + 1][80][0]);
		nodesc2[i][80][1] = 0.5 * (nodes2[i][80][1] + nodes2[i + 1][80][1]);
	}
	nodesc2[0][0][0] = nodes2[1][1][0];
	nodesc2[0][0][1] = nodes2[1][0][1];
	nodesc2[260][0][0] = nodes2[260][1][0];
	nodesc2[260][0][1] = nodes2[260][0][1];
	nodesc2[260][80][0] = nodes2[260][80][0];
	nodesc2[260][80][1] = nodes2[260][80][1];
	nodesc2[0][80][0] = nodes2[1][80][0];
	nodesc2[0][80][1] = nodes2[1][80][1];
	for (j = 1; j < 80; j++)
	{
		pho2_1[0][j] = 0.5 * (pho2_1[0][j] + pho2_1[1][j]);
		vx2_1[0][j] = 0.5 * (vx2_1[0][j] + vx2_1[1][j]);
		vy2_1[0][j] = 0.5 * (vy2_1[0][j] + vy2_1[1][j]);
		pre2_1[0][j] = 0.5 * (pre2_1[0][j] + pre2_1[1][j]);
		pho2_1[260][j] = 0.5 * (pho2_1[260][j] + pho2_1[259][j]);
		vx2_1[260][j] = 0.5 * (vx2_1[260][j] + vx2_1[259][j]);
		vy2_1[260][j] = 0.5 * (vy2_1[260][j] + vy2_1[259][j]);
		pre2_1[260][j] = 0.5 * (pre2_1[260][j] + pre2_1[259][j]);
	}
	for (i = 1; i < 260; i++)
	{
		pho2_1[i][0] = 0.5 * (pho2_1[i][0] + pho2_1[i][1]);
		vx2_1[i][0] = 0.5 * (vx2_1[i][0] + vx2_1[i][1]);
		vy2_1[i][0] = 0.5 * (vy2_1[i][0] + vy2_1[i][1]);
		pre2_1[i][0] = 0.5 * (pre2_1[i][0] + pre2_1[i][1]);
		pho2_1[i][80] = 0.5 * (pho2_1[i][80] + pho2_1[i][79]);
		vx2_1[i][80] = 0.5 * (vx2_1[i][80] + vx2_1[i][79]);
		vy2_1[i][80] = 0.5 * (vy2_1[i][80] + vy2_1[i][79]);
		pre2_1[i][80] = 0.5 * (pre2_1[i][80] + pre2_1[i][79]);
	}
	pho2_1[0][0] = pho2_1[1][0];
	vx2_1[0][0] = vx2_1[1][0];
	vy2_1[0][0] = vy2_1[1][0];
	pre2_1[0][0] = pre2_1[1][0];

	pho2_1[0][80] = pho2_1[1][80];
	vx2_1[0][80] = vx2_1[1][80];
	vy2_1[0][80] = vy2_1[1][80];
	pre2_1[0][80] = pre2_1[1][80];

	pho2_1[260][0] = pho2_1[259][0];
	vx2_1[260][0] = vx2_1[259][0];
	vy2_1[260][0] = vy2_1[259][0];
	pre2_1[260][0] = pre2_1[259][0];

	pho2_1[260][80] = pho2_1[259][80];
	vx2_1[260][80] = vx2_1[259][80];
	vy2_1[260][80] = vy2_1[259][80];
	pre2_1[260][80] = pre2_1[259][80];
}

int main()
{
	mesh_generation();

	initialize();

	area_caculate();

	fg = fopen("error.txt", "w");

	for (k = 0; k < 2000; k++)
	{

		printf("%d    %.15f    %.15f    %.15f    %.15f    %.15f\n", k, maxflux, maxflux2, maxflux3, maxflux4, T2_1[3][50]);
		fprintf(fg, "%d    %.15f    %.15f    %.15f    %.15f    %.15f\n", k, maxflux, maxflux2, maxflux3, maxflux4, T2_1[3][50]);
		boundary_conditions();
		roe();
		iteration();
	}

	fp = fopen("1.txt", "w");

	fprintf(fp, "x                 ,y              pho1_1[i][j],    vx1_1[i][j],    vy1_1[i][j],      pre1_1[i][j],      T1_1[i][j],      ma1_1[i][j]\n");

	for (i = 1; i < 330; i++)
	{
		for (j = 1; j < 70; j++)
		{
			T1_1[i][j] = pre1_1[i][j] / R / pho1_1[i][j];
			ma1_1[i][j] = sqrt(vx1_1[i][j] * vx1_1[i][j] + vy1_1[i][j] * vy1_1[i][j]) / sqrt(gama * R * T1_1[i][j]);
			fprintf(fp, "%.10f    %.10f    %.10f    %.10f    %.10f    %.10f    %.10f    %.10f\n", nodes[i][j][0], nodes[i][j][1], pho1_1[i][j], vx1_1[i][j], vy1_1[i][j], pre1_1[i][j], T1_1[i][j], ma1_1[i][j]);
		}
	}
	fq = fopen("2.txt", "w");
	fprintf(fq, "x                 ,y              pho1_1[i][j],    vx1_1[i][j],    vy1_1[i][j],      pre1_1[i][j],      T1_1[i][j],      ma1_1[i][j]\n");

	for (i = 1; i < 260; i++)
	{
		for (j = 1; j < 80; j++)
		{
			fprintf(fq, "%.10f    %.10f    %.10f    %.10f    %.10f    %.10f    %.10f    %.10f\n", nodes2[i][j][0], nodes2[i][j][1], pho2_1[i][j], vx2_1[i][j], vy2_1[i][j], pre2_1[i][j], T2_1[i][j], ma2_1[i][j]);
		}
	}

	record();

	fr = fopen("WYplotflow.dat", "w");
	fprintf(fr, "Title=\"NOZZLE\"\nVariables=\"x\",\"y\",\"dens\",\"velx\",\"vely\",\"spre\",\"ttem\",\"mach\"\nZone T=\"NOZZLE\" i=331,j=71,f=point \n");
	for (j = 0; j < 71; j++)
	{
		for (i = 0; i < 331; i++)
		{
			T1_1[i][j] = pre1_1[i][j] / R / pho1_1[i][j];
			ma1_1[i][j] = sqrt(vx1_1[i][j] * vx1_1[i][j] + vy1_1[i][j] * vy1_1[i][j]) / sqrt(gama * R * T1_1[i][j]);
			fprintf(fr, "%.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f\n", nodesc1[i][j][0], nodesc1[i][j][1], pho1_1[i][j], vx1_1[i][j], vy1_1[i][j], pre1_1[i][j], T1_1[i][j], ma1_1[i][j]);
		}
	} //////////////
	fs = fopen("WYplotflow2.dat", "w");
	fprintf(fs, "Title=\"NOZZLE\"\nVariables=\"x\",\"y\",\"dens\",\"velx\",\"vely\",\"spre\",\"ttem\",\"mach\"\nZone T=\"NOZZLE\" i=261,j=81,f=point \n");
	for (j = 0; j < 81; j++)
	{
		for (i = 0; i < 261; i++)
		{
			T2_1[i][j] = pre2_1[i][j] / R / pho2_1[i][j];
			ma2_1[i][j] = sqrt(vx2_1[i][j] * vx2_1[i][j] + vy2_1[i][j] * vy2_1[i][j]) / sqrt(gama * R * T2_1[i][j]);
			fprintf(fs, "%.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f\n", nodesc2[i][j][0], nodesc2[i][j][1], pho2_1[i][j], vx2_1[i][j], vy2_1[i][j], pre2_1[i][j], T2_1[i][j], ma2_1[i][j]);
		}
	}
	printf("nodesc[2][00][1]=%f\n", nodesc1[2][0][1]);
}
