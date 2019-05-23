#include"main.H"

void init(double U[][maxSpace+1], double F[][maxSpace+1])
{
    double u1=0.0, rho1=1.0, p1=1.0;
    double u2=0.0, rho2=0.125, p2=0.1;

    for(int i=0; i<= maxSpace/2; i++)
    {
        U[0][i]=u1;
        U[1][i]=u1*rho1;
        U[2][i]=p1/(GAMMA-1) + 0.5*rho1*u1*u1;
    }

    for(int i=maxSpace/2+1; i<= maxSpace; i++)
    {
        U[0][i]=u2;
        U[1][i]=u2*rho2;
        U[2][i]=p2/(GAMMA-1) + 0.5*rho2*u2*u2;
    }

    UToF(U,F);
}