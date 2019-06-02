#include"main.H"

void scalarJSTConv(double W[][3], const double dt, double R[][3], const int I)
{
    const double k2=0.5, k4=0.01;
   // const double Lambda, Ep2, Ep4, Y, Y1,;
    double D[3],  Wf[3], F[3];
    if(I>=1 && I<=maxSpace-2)
    {
        //先定义一些系数, Lambda和Ep2, Ep4, 用于计算D[k]
        //其中Y与Y1是开关函数,用于切换Ep2和Ep4
        const double Lambda = 0.5 * (lambda(W[I]) + lambda(W[I + 1]));

        const double p1=p(W[I+1]), p2=p(W[I + 2]), p0=p(W[I]), p_1=p(W[I - 1]);
        const double Y = fabs(p1- 2 * p0 + p_1) 
                      / (p1+ 2 * p0 + p_1);

        const double Y1 = fabs(p2 - 2 * p1+ p0) 
                        / (p2 + 2 * p1+ p0);

        const double Ep2 = k2 * max(Y, Y1);
        const double Ep4 = max(0, k4 - Ep2);

        //然后计算人工黏性D, 计算face上所用的流动变量Wf(中心型简单算数平均)
        for (int k = 0; k < 3; k++)
        {
            D[k] = Lambda 
                    * (   Ep2 * (W[I + 1][k] - W[I][k]) 
                        - Ep4 * (W[I + 2][k] - 3 * W[I + 1][k] + 3 * W[I][k] - W[I - 1][k])
                      );

            Wf[k] = 0.5 * (W[I][k] + W[I + 1][k]);
        }
    }

    else
    {
        //外推边界网格的W, 第一个时间步内是根据初场外推的
        zeroGradBC(W);
    }

    //face上的流动变量Wf转换为通量F
//    WToF(Wf, F);
    double p0=p(W[I]);
    double u=W[I][1]/W[I][0];
    F[0] = W[I][1];
    F[1] = W[I][0] * u * u + p0;
    F[2] = (W[I][2] + p0) * u;

    //最后算出残差R
    for (int k = 0; k < 3; k++)
    {
        R[I][k] = F[k] - D[k];
    }
}