#include"main.H"
void EulerFTime(double W[][3], const double dt,  double R[][3], const int I)
{
    CONV(W,dt,R, I);

    for (int k = 0; k < 3; k++)
    {
        W[I][k] = dt * (-1 / dx) * R[I][k] + W[I][k];
    }
}