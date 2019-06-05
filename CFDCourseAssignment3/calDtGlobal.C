#include"main.H"
double calDtGlobal(double W[][3])
{
    double min=1e10;
    double dtLocal[maxSpace];
    for (int i = 1; i <=maxSpace ; i++)
    {
        double u = W[i][1]/W[i][0];
        double pre= (GAMMA-1) * W[i][0] *( W[i][2]/W[i][0] - 0.5* u * u );
        double LLam=fabs(W[i][1]/W[i][0]) + ( GAMMA * pre/W[i][0] );

        dtLocal[i]=CFL*dx/LLam;//注意此i的生存期
        
        if(dtLocal[i]<min)
        {
            min=dtLocal[i];
        }
    }
    return min;
}