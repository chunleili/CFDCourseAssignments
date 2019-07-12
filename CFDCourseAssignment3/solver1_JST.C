#include "main.H"

//一阶精度的三阶显式RungeKutta法
void solver1(double W[][3], const double dt, double R[][3])
{
    const double alpha1=0.1481, alpha2=0.4, alpha3=1;
    //先定义W0,用于保存原始的W
    double W0[maxSpace+3][3];
    for (int I = 1; I <= maxSpace; I++)
    {
        for (int k = 0; k < 3; k++)
        {
            W0[I][k] = W[I][k];
        }
    }
    //后面每一步都先计算残差, 后根据RK公式更新W
    scalarJSTConv(W,  R);
    for (int I = 1; I <= maxSpace; I++)
    {
        for (int k = 0; k < 3; k++)
        {
            W[I][k] = W0[I][k] - alpha1 * dt / dx * R[I][k];
        }
    }

    scalarJSTConv(W, R);
    for (int I = 1; I <= maxSpace; I++)
    {
        for (int k = 0; k < 3; k++)
        {
            W[I][k] = W0[I][k] - alpha2 * dt / dx * R[I][k];
        }
    }

    scalarJSTConv(W,  R);
    for (int I = 1; I <= maxSpace; I++)
    {
        for (int k = 0; k < 3; k++)
        {
            W[I][k] = W0[I][k] - alpha3 * dt / dx * R[I][k];
        }
    }
}



void scalarJSTConv(const double W[][3],  double R[][3])
{
    const double k2 = 0.5, k4 = 0.01;  //JST法中的常数 
    for (int I = 1; I <= maxSpace; I++)
    {
        //先定义一些系数, Lambda和Ep2, Ep4, 用于计算D[k]
        //其中Y与Y1是开关函数,用于切换Ep2和Ep4
         double Lambda_f = 0.5 * (lambda(W[I]) + lambda(W[I + 1]));

         double p1 = calPressure(W[I + 1]), p2  = calPressure(W[I + 2]);
         double p0 = calPressure(W[I]),     p_1 = calPressure(W[I - 1]);
         double Y = fabs(p1 - 2 * p0 + p_1) / (p1 + 2 * p0 + p_1);
         double Y1 = fabs(p2 - 2 * p1 + p0) / (p2 + 2 * p1 + p0);

         double Ep2 = k2 * max(Y, Y1);
         double Ep4 = max(0, k4 - Ep2);

        //然后计算人工黏性D, 计算face上所用的流动变量Wf(中心型简单算数平均)
        double D[3], Wf[3];
        for (int k = 0; k < 3; k++)
        {
            D[k] = Lambda_f *
                (  Ep2 * (W[I + 1][k] - W[I][k]) 
                 - Ep4 * (W[I + 2][k] - 3 * W[I + 1][k] + 3 * W[I][k] - W[I - 1][k])
                 );

            Wf[k] = 0.5 * (W[I][k] + W[I + 1][k]);
        }

        double F_left[3];
        double static F_right[3] = {0, 1, 0}; //static变量会保存上一次调用函数时F_right的值
        //上一空间点的右面通量恰好为这一点的左通量F_left[],取向右为正,它应该是负方向.
        //在第一个网格,Fring被初始化为0,1,0;所以经过下面的语句,Fleft也被初始化
        for (int k = 0; k < 3; k++)
        {
            F_left[k] = -F_right[k];
        }

        //face上的流动变量Wf转换为通量F, 从而更新了这一步的F_right
        WToF(Wf, F_right);

        //最后算出残差R
        for (int k = 0; k < 3; k++)
        {
            R[I][k] = F_left[k] + F_right[k] - D[k];
        }
    }
}