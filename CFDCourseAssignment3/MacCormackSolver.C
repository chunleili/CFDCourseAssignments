#include"main.H"

void MacCormackSolver( double U[3][maxSpace+1], double F[3][maxSpace+1])
{
    const double r = dt/dx;
    const double eta = 0.25;

    double theta;
    int i, k;

    double U_new[3][maxSpace+1];
    double U_half[3][maxSpace+1];
    double F_half[3][maxSpace+1];


    //开关
    for ( i = 1; i < maxSpace; i++)
    {
        theta = fabs( ( fabs(U[0][i + 1] -U[0][i]) - fabs(U[0][i] -U[0][i-1])        )   
                 /    ( fabs(U[0][i + 1] -U[0][i]) + fabs(U[0][i] -U[0][i-1]) +1e-100) ) ;

    }

    //人工粘性
    for (k = 0; k < 3; k++)
    {
        for (i = 0; i <= maxSpace; i++)
        {
            U[k][i] = U[k][i] + 0.5 * eta * theta * (U[k][i + 1] - 2 * U[k][i] + U[k][i - 1]);
        }
    }

    UToF(U, F);

    //预报
    for (k = 0; k < 3; k++)
    {
        for (i = 1; i <= maxSpace; i++)
        {
            U_half[k][i] = U[k][i] - r * (F[k][i] - F[k][i - 1]);
        }
    }

    UToF(U_half, F_half);
    
    //校正
    for (k = 0; k < 3; k++)
    {
        for (i = 1; i <= maxSpace; i++)
        {
            U_new[k][i] = 0.5 * (U[k][i] + U_half[k][i]) - 0.5 * r * (F_half[k][i + 1] - F_half[k][i]);
        }
    }
}