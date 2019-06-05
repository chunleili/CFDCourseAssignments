#include "main.H"
//			     0方向梯度边界条件, 边界上的点直接等于内部的点                      //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void init(double W[][3])
{
    cout<<"\nInitializing the fields..."<<endl;
    double u1=0.0, rho1=1.0, p1=1.0;
    double u2=0.0, rho2=0.125, p2=0.1;

    for(int i=0; i< maxSpace/2; i++)
    {
        W[i][0]=rho1;
        W[i][1]=u1*rho1;
        W[i][2]=p1/(GAMMA-1) + 0.5*rho1*u1*u1;
    }

    for(int i=maxSpace/2; i<= maxSpace+2; i++)
    {
        W[i][0]=rho2;
        W[i][1]=u2*rho2;
        W[i][2]=p2/(GAMMA-1) + 0.5*rho2*u2*u2;
    }
}

//			     0方向梯度边界条件, 边界上的点直接等于内部的点                      //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void zeroGradBC(double W[][3])
{
    for (int k = 0; k < 3; k++)
    {
        W[0][k] = W[1][k];
        W[maxSpace][k] = W[maxSpace - 2][k];
        W[maxSpace][k] = W[maxSpace - 1][k];
    }
}

//			     打印                      //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include<fstream>
#include<iomanip>
void print(const double W[][3])
{

    ofstream foutAll("data/result.dat");
    foutAll.setf(ios::left);
    foutAll.precision(4);
    cout<<"\nprint the rho, u, p, in the data/result.dat, first column is distance"<<endl;
    foutAll<<"x"<<'\t'<<"rho"<<'\t'<<"u"<<'\t'<<"p"<<endl;
    for (int i = 0; i <= maxSpace+2; i++)
    {
        foutAll<<i*dx<<'\t'<<W[i][0]<<'\t'<<W[i][1]/W[i][0]<<'\t'<<calPressure(W[i])<<endl;
    } 
}

//			     W转化为F                      //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void WToF(double const W[3], double F[3])
{
    double p0=calPressure(W);
    double u=W[1]/W[0];
    F[0] = W[1];
    F[1] = W[0] * u * u + p0;
    F[2] = (W[2] + p0) * u;
}

//			     计算全局dt                     //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double calDtGlobal(double W[][3])
{
    double min=1e10;
    double dtLocal[maxSpace];
    for (int i = 1; i <=maxSpace ; i++)
    {
        dtLocal[i]=localTime(W,i);//注意此i的生存期
        if(dtLocal[i]<min)
        {
            min=dtLocal[i];
        }
    }
    return min;
}