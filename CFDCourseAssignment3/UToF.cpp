#include"main.H"

void UToF(double U[][maxSpace+2], double F[][maxSpace+2])
{

    double u, p;

    for (int i = 0; i < maxSpace + 2; i++)
    {
        u = U[1][i] / U[0][i];
        p = (GAMMA - 1) * (U[2][i] - 0.5 * U[1][i] * U[1][i] / U[0][i]);

        F[0][i] = U[1][i];
        F[1][i] = U[0][i]* u * u + p;
        F[2][i] = (U[2][i] + p) * u;
    }
}