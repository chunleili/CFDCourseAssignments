#include"main.H"

double localTime(const double W[][3], const int I)
{
    return CFL*dx/lambda(W[I]);
}