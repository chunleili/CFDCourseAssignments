#include"main.H"
extern double rho[maxSpace+1], u[maxSpace+1], p[maxSpace+1], E[maxSpace+1];
void init(double W[][3])
{
    cout<<"\nInitializing the fields..."<<endl;
    double u1=0.0, rho1=1.0, p1=1.0;
    double u2=0.0, rho2=0.125, p2=0.1;

    for(int i=0; i< maxSpace/2; i++)
    {
        u[i]=u1;
        rho[i]=rho1;
        p[i]=p1;

        W[i][0]=rho1;
        W[i][1]=u1*rho1;
        W[i][2]=p1/(GAMMA-1) + 0.5*rho1*u1*u1;
    }

    for(int i=maxSpace/2; i<= maxSpace+2; i++)
    {
        u[i]=u2;
        rho[i]=rho2;
        p[i]=p2;

        W[i][0]=rho2;
        W[i][1]=u2*rho2;
        W[i][2]=p2/(GAMMA-1) + 0.5*rho2*u2*u2;
    }

}