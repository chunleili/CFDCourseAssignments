#include"main.H"

void MacCormackSolver( double U[maxSpace+1][3], double F[maxSpace+1][3])
{
    const double r = dt/dx;
    const double eta = 0.25;

    double theta;
    int i, k;

    double U_new[maxSpace+1][3];
    double U_half[maxSpace+1][3];
    double F_half[maxSpace+1][3];


    //开关
    for ( i = 1; i < maxSpace; i++)
    {
        theta = fabs( ( fabs(U[i + 1][0]  -U[i][0]) - fabs(U[i][0] -U[i-1][0])        )  
                 /    ( fabs(U[i + 1][0]  -U[i][0]) + fabs(U[i][0] -U[i-1][0]) +1e-100) ) ;

    }

    //人工粘性
    for (k = 0; k < 3; k++)
    {
        for (i = 0; i <= maxSpace; i++)
        {
            U[i][k] = U[i][k] + 0.5 * eta * theta * (U[i + 1][k] - 2 * U[i][k] + U[i - 1][k]);
        }
    }

    UToF(U, F);

    //预报
    for (k = 0; k < 3; k++)
    {
        for (i = 1; i <= maxSpace; i++)
        {
            U_half[i][k] = U[i][k] - r * (F[i][k] - F[i - 1][k]);
        }
    }

    UToF(U_half, F_half);
    
    //校正
    for (k = 0; k < 3; k++)
    {
        for (i = 1; i <= maxSpace; i++)
        {
            U_new[i][k] = 0.5 * (U[i][k] + U_half[i][k]) - 0.5 * r * (F_half[i [k]+ 1] - F_half[i][k]);
        }
    }
}