#include"main.H"
//一阶精度的三阶显式RungeKutta法
void RungeKuttaTime(double W[][3], double dtLocal, double R[][3], const int I)
{
    const double alpha1=0.1481, alpha2=0.4, alpha3=1;
    //先定义W0,用于保存原始的W
    double W0[3];
    for (int k = 0; k < 3; k ++)
    {
        W0[k] = W[I][k];
    }

    //后面每一步都先计算该单元格的残差, 后根据RK公式更新W
    CONV(W, dtLocal, R, I);
    for (int k = 0; k < 3; k ++)
    {
        W[I][k] = W0[k] - alpha1 * dtLocal / dx * R[I][k];
    }

    CONV(W, dtLocal, R, I);
    for (int k = 0; k < 3; k ++)
    {
        W[I][k] = W0[k] - alpha2 * dtLocal / dx * R[I][k];
    }

    CONV(W, dtLocal, R, I);
    for (int k = 0; k < 3; k ++)
    {
        W[I][k] = W0[k] - alpha3 * dtLocal / dx * R[I][k];
    }

}