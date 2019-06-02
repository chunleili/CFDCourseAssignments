#include"main.H"
void EulerFTime(double W[][3], const double dt, const double R[][3])
{
    CONV(W,dt,R);
    for (int I = 0; I <= maxSpace; I++)
    {
        for (int k = 0; k < 3; k++)
        {
            W[I][k] = dt * (-1/dx)*R[I][k] + W[I][k];
        }
    }
}