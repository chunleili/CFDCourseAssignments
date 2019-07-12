#include"main.H"
void MacCormackConv( double U[][3], const double dt)
{
    const double r = dt/dx;
    const double eta = 0.25;

    double theta;
    int i, k;

    double Uf[maxSpace+1][3];
    double Ef[maxSpace+1][3];


    //开关
    for ( i = 1; i < maxSpace; i++)
    {
        theta = fabs( ( fabs(U[i + 1][0]  -U[i][0]) - fabs(U[i][0] -U[i-1][0])        )  
                 /    ( fabs(U[i + 1][0]  -U[i][0]) + fabs(U[i][0] -U[i-1][0]) +1e-100) ) ;
        //人工粘性
        for (k = 0; k < 3; k++)
        {
            Ef[i][k] = U[i][k] + 0.5 * eta * theta * (U[i + 1][k] - 2 * U[i][k] + U[i - 1][k]);
        }
    }

    for ( i = 1; i < maxSpace; i++)
    {
        for ( k = 0; k < 3; k++)
        {
            U[i][k]=Ef[i][k];
        }
        
    }
    
    for(i=0;i<=maxSpace;i++)
    {
        WToF(U[i], Ef[i]);
    }

    //预报
    for (i = 1; i < maxSpace; i++)
    {
        for (k = 0; k < 3; k++)
        {
            Uf[i][k] = U[i][k] - r * (Ef[i+1][k] - Ef[i][k]);
        }
    }

    //校正
    for (k = 0; k < 3; k++)
    {
        for (i = 1; i <= maxSpace; i++)
        {
            U[i][k] = 0.5 * (U[i][k] + Uf[i][k]) - 0.5 * r * (Ef[i][k] - Ef[i][k]);
        }
    }

    //边界条件
    for (int k  = 0; k < 3; k++)
    {
        //左边界
        U[0][k]=U[1][k];
        //右边界
        U[maxSpace][k]=U[maxSpace-1][k];
    }
}