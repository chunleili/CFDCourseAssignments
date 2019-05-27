#include"main.H"

double CFL(const double U[][3])
{
    double maxVel=1e-100, Lambda;
    const double Sf=0.8;
    for (int i = 0; i <=maxSpace; i++)
    {
        Lambda=lambda(U[i]);
        if (Lambda>maxVel)
             maxVel=Lambda;
    }
    return Sf*dx/maxVel;
}