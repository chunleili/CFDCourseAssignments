#include"main.H"

void init(double U[][3])
{
    double u1=0.0, rho1=1.0, p1=1.0;
    double u2=0.0, rho2=0.125, p2=0.1;

    for(int i=0; i<= maxSpace/2; i++)
    {
        U[i][0]=u1;
        U[i][1]=u1*rho1;
        U[i][2]=p1/(GAMMA-1) + 0.5*rho1*u1*u1;
    }

    for(int i=maxSpace/2+1; i<= maxSpace; i++)
    {
        U[i][0]=u2;
        U[i][1]=u2*rho2;
        U[i][2]=p2/(GAMMA-1) + 0.5*rho2*u2*u2;
    }

}