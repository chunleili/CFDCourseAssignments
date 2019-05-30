#include"main.H"

void WToF(double const W[3], double F[3])
{

    double p=p(W);
    double u=W[1]/W[0];
    F[0] = W[1];
    F[1] = W[0] * u * u + p;
    F[2] = (W[2] + p) * u;
}