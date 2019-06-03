#include "main.H"
//			     0方向梯度边界条件, 边界上的点直接等于内部的点                      //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void zeroGradBC(double W[][3])
{
    for (int k = 0; k < 3; k++)
    {
        W[0][k] = W[1][k];
        W[maxSpace][k] = W[maxSpace - 2][k];
        W[maxSpace][k] = W[maxSpace - 1][k];
    }
}