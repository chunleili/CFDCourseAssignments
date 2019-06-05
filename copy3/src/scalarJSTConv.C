#include "main.H"

void scalarJSTConv(const double W[][3], double R[][3], const int I)
{
    const double k2 = 0.5, k4 = 1.0/64; //JST法中的常数
                                     //先定义一些系数, Lambda和Ep2, Ep4, 用于计算D[k]
                                     //其中Y与Y1是开关函数,用于切换Ep2和Ep4
    const double Lambda_f = 0.5 * (calLambda(W[I]) + calLambda(W[I + 1]));

    const double p1 = calPressure(W[I + 1]), p2 = calPressure(W[I + 2]);
    const double p0 = calPressure(W[I]), p_1 = calPressure(W[I - 1]);
    const double Y = fabs(p1 - 2 * p0 + p_1) / (p1 + 2 * p0 + p_1);
    const double Y1 = fabs(p2 - 2 * p1 + p0) / (p2 + 2 * p1 + p0);

    const double Ep2 = k2 * max(Y, Y1);
    const double Ep4 = max(0, k4 - Ep2);

    //然后计算人工黏性D, 计算face上所用的流动变量Wf(中心型简单算数平均)
    double D[3], Wf[3];
    double F_left[3];
    double static F_right[3] = {0, 1, 0}; //static变量会保存上一次调用函数时F_right的值
    //本文Euler方程中, 通量流出表面应该为负, 所以通量F应该取向内为正
    //上一空间点的右面通量恰好为这一点的左通量F_left[],取向右为正,它应该是负方向.
    for (int k = 0; k < 3; k++)
    {
        D[k] =0;// Lambda_f *(Ep2 * (W[I + 1][k] - W[I][k]) - Ep4 * (W[I + 2][k] - 3 * W[I + 1][k] + 3 * W[I][k] - W[I - 1][k]));

        Wf[k] = 0.5 * (W[I][k] + W[I + 1][k]);

        F_left[k] = -F_right[k];
    }

    //face上的流动变量Wf转换为通量F_right, 从而更新了这一步的F_right
    WToF(Wf, F_right);
    
    //最后算出残差R, 残差是通量在所有表面的和
    for (int k = 0; k < 3; k++)
    {
        R[I][k] = F_left[k] + F_right[k];
    }
}