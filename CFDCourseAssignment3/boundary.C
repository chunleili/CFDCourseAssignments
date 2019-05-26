#include"main.H"
void boundary(double U[][3])
{
    for (int k  = 0; k < 3; k++)
    {
        //左边界
        U[0][k]=U[1][k];
        //右边界
        U[maxSpace][k]=U[maxSpace-1][k];
    }
    
}