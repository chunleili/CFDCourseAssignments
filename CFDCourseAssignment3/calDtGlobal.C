#include"main.H"
double calDtGlobal(double W[][3])
{
    double min=1e10;
    double dtLocal;
    for (int i = 0; i <=maxSpace ; i++)
    {
        dtLocal=localTime(W,i);//注意此i的生存期
        if(dtLocal<min)
        {
            min=dtLocal;
        }
    }
    return min;
}