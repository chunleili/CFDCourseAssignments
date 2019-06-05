#include "main.H"
#include <cmath>
//采用单元中心法, 实际存储网格从W[1] 到W[maxSpace]; W[0]与W[maxSpace+1],W[maxSpace+2]是边界虚网格
double F_right[maxSpace + 3][3] = {0, 1, 0};
double W[maxSpace + 3][3], R[maxSpace + 3][3];
double rho[maxSpace + 1], u[maxSpace + 1], p[maxSpace + 1], E[maxSpace + 1];

void JST();

int main()
{
    init(W);
    double dtGlobal;
    int timeStep = 0;
    for (double t = 0; t <= stopTime && timeStep < 10000; t += dtGlobal)
    {
        dtGlobal = calDtGlobal(W);
        double dt=dtGlobal;
        JST();
        //三阶荣格库塔法
        const double alpha1 = 0.1481, alpha2 = 0.4, alpha3 = 1;
        //先定义W0,用于保存原始的W
        double W0[maxSpace + 3][3];
        for (int I = 1; I <= maxSpace; I++)
        {
            for (int k = 0; k < 3; k++)
            {
                W0[I][k] = W[I][k];
            }
        }
        JST();
        for (int I = 1; I <= maxSpace; I++)
        {
            for (int k = 0; k < 3; k++)
            {
                W[I][k] = W0[I][k] - alpha1 * dt / dx * R[I][k];
            }
        }

        JST();
        for (int I = 1; I <= maxSpace; I++)
        {
            for (int k = 0; k < 3; k++)
            {
                W[I][k] = W0[I][k] - alpha2 * dt / dx * R[I][k];
            }
        }

        JST();
        for (int I = 1; I <= maxSpace; I++)
        {
            for (int k = 0; k < 3; k++)
            {
                W[I][k] = W0[I][k] - alpha3 * dt / dx * R[I][k];
            }
        }

        timeStep++;
        cout << "dtGlobal= " << dtGlobal << "\t time step = " << timeStep << "\t time = " << t << endl;
    }
    print(W);
    return 0;
}

void JST()
{
    const double k2 = 0.5, k4 = 0.015625; //JST法中的常数
    double c[maxSpace + 1], lam[maxSpace + 1], Y[maxSpace + 1];
    for (int I = 1; I <= maxSpace; I++)
    {
        rho[I] = W[I][0];
        u[I] = W[I][1] / W[I][0];
        E[I] = W[I][2] / W[I][0];
        p[I] = (GAMMA - 1) * rho[I] * (E[I] - 0.5 * u[I] * u[I]);

        //先定义一些系数, Lambda和Ep2, Ep4, 用于计算D[k]
        //其中Y是开关函数,用于切换Ep2和Ep4
        c[I] = sqrt(GAMMA * p[I] / rho[I]);
        lam[I] = fabs(u[I]) + c[I];
        double Lambda_f = 0.5 * (lam[I] + lam[I + 1]);
        Y[I] = fabs(p[I + 1] - 2 * p[I] + p[I - 1]) / (p[I + 1] + 2 * p[I] + p[I - 1]);
        double Ep2 = k2 * max(Y[I], Y[I + 1]);
        double Ep4 = max(0, k4 - Ep2);

        double D[3], Wf[3];
        for (int k = 0; k < 3; k++)
        {
            //然后计算人工黏性D
            D[k] = Lambda_f * (Ep2 * (W[I + 1][k] - W[I][k]) - Ep4 * (W[I + 2][k] - 3 * W[I + 1][k] + 3 * W[I][k] - W[I - 1][k]));

            //计算face上所用的流动变量Wf(中心型简单算数平均)
            Wf[k] = 0.5 * (W[I][k] + W[I + 1][k]);

            //face上的流动变量Wf转换为通量F_right, 从而更新了这一步的F_right
            F_right[I][0] = Wf[1] - D[1];
            F_right[I][1] = Wf[1] * Wf[1] / Wf[0] + p[I] - D[2];
            F_right[I][2] = (Wf[2] + p[I]) * Wf[1] / Wf[0] - D[3];

            //最后算出残差R, 残差是通量在所有表面的和
            R[I][k] = F_right[I - 1][k] - F_right[I][k];
        }
    }
}
