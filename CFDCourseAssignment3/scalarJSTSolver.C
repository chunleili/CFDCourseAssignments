#include"main.H"

inline static double max(double a, double b)
{
    return a>b?a:b;
}

inline static double c(const double U[3])
{
    return sqrt( GAMMA * p(U)/U[0] );
}

inline static double lambda(const double U[3])
{
    return fabs(U[1]/U[0] + c(U)) * dx;
}


void scalarJSTSolver(double U[][3])
{
    const double k2=0.5, k4=0.01;
    double D, Lambda, Ep2, Ep4, Y, Y1, U_half[3], R[3], F[3];
    double U1[maxSpace+1][3];
    for (int I = 1; I <= maxSpace - 2; I++)
    {

        Lambda = 0.5 * (lambda(U[I]) + lambda(U[I + 1]));

        Y = fabs(p(U[I + 1]) - 2 * p(U[I]) + p(U[I - 1])) 
             / (p(U[I + 1]) + 2 * p(U[I]) + p(U[I - 1]));

        Y1 = fabs(p(U[I + 2]) - 2 * p(U[I + 1]) + p(U[I])) 
             / (p(U[I + 2]) + 2 * p(U[I + 1]) + p(U[I]));

        Ep2 = k2 * max(Y, Y1);
        Ep4 = max(0, k4 - Ep2);

        for (int k = 0; k < 3; k++)
        {
            D = Lambda * (Ep2 * (U[I + 1][k] - U[I][k]) 
                        - Ep4 * (U[I + 2][k] - 3 * U[I + 1][k] + 3 * U[I][k] - U[I - 1][k])
                        );

            U_half[k] = 0.5 * (U[I][k] + U[I + 1][k]);

            UToF(U_half, F);

            R[k] = F[k] * dx - D;

            U1[I][k] = dt * R[k] + U[I][k];
        }
    }

    for (int I = 0; I <= maxSpace; I++)
    {
        for (int k = 0; k < 3; k++)
        {
            U[I][k]=U1[I][k];
        }
    }

}