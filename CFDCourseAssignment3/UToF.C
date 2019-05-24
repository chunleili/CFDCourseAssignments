#include"main.H"

void UToF(double U[3], double F[3])
{

    double u;

//    for (int i = 0; i <= maxSpace ; i++)
    {
        u = U[1] / U[0];
        
        F[0] = U[1];
        F[1] = U[0]* u * u + p(U) ;
        F[2] = (U[2] + p(U) ) * u;
    }
}