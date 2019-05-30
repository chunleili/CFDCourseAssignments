#include"main.H"
void EulerFTime(double U[][3], const double dt, const double R[3])
{
    for (int i = 0; i <= maxSpace; i++)
    {
        for (int k = 0; k < 3; k++)
        {
            U[I][k] = dt * R[I][k] + U[I][k];
        }
    }
}