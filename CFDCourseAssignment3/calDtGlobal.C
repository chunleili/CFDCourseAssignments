#include"main.H"
double calDtGlobal(double W[][3])
{
    double min=1e10;
    double dtLocal;
    for (int I = 0; I <=maxSpace ; I++)
    {
        dtLocal=localTime(W,I);
        if(dtLocal<min)
        {
            min=dtLocal;
        }
    }
    return min;
}