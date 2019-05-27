#include"main.H"

double localTime(const double U[][3])
{
    double maxVel=1e-100, Lambda;
    const double CFLnumber=0.5;
    for (int i = 0; i <=maxSpace; i++)
    {
        Lambda=lambda(U[i]);
        if (Lambda>maxVel)
             maxVel=Lambda;
    }
    return CFLnumber*dx/maxVel;
}