#include"main.H"

inline static double max(double a, double b)
{
    return a>b?a:b;
}
//计算声速c
inline static double c(const double U[3])
{
    return sqrt( GAMMA * p(U)/U[0] );
}

inline static double lambda(const double U[3])
{
    return fabs(U[1]/U[0]) + c(U);
}


void scalarJSTSolver(double U[][3])
{
    const double k2=0.5, k4=0.01;
    double D[3], Lambda, Ep2, Ep4, Y, Y1, Uf[3], R[3], F[3];
    double U1[maxSpace+1][3];
    for (int I = 1; I <= maxSpace - 2; I++)
    {
        //先定义一些系数, Lambda和Ep2, Ep4, 用于计算D[k]
        //其中Y与Y1是开关函数,用于切换Ep2和Ep4
        Lambda = 0.5 * (lambda(U[I]) + lambda(U[I + 1]));

        Y = fabs(p(U[I + 1]) - 2 * p(U[I]) + p(U[I - 1])) 
             / (p(U[I + 1]) + 2 * p(U[I]) + p(U[I - 1]));

        Y1 = fabs(p(U[I + 2]) - 2 * p(U[I + 1]) + p(U[I])) 
             / (p(U[I + 2]) + 2 * p(U[I + 1]) + p(U[I]));

        Ep2 = k2 * max(Y, Y1);
        Ep4 = max(0, k4 - Ep2);

        //然后计算人工黏性D, 计算face上所用的流动变量Uf(中心型简单算数平均)
        for (int k = 0; k < 3; k++)
        {
            D[k] = Lambda * (Ep2 * (U[I + 1][k] - U[I][k]) 
                        - Ep4 * (U[I + 2][k] - 3 * U[I + 1][k] + 3 * U[I][k] - U[I - 1][k])
                        );

            Uf[k] = 0.5 * (U[I][k] + U[I + 1][k]);
        }

        //face上的流动变量Uf转换为通量F
        UToF(Uf, F);
        
        //最后算出残差R, 用欧拉后插离散时间项,算出下一时间步流动变量U1
        for(int k =0; k<3; k++)
        {
            R[k] = F[k]  - D[k];

            U1[I][k] = dt * R[k] + U[I][k];
        }
    }

    //边界条件,直接取零法向梯度
    for (int k = 0; k < 3; k++)
    {
        U1[0][k] = U1[1][k];
        U1[maxSpace - 1][k] = U1[maxSpace - 2][k];
        U1[maxSpace][k] = U1[maxSpace - 1][k];
    }

    //更新下一时间步
    for (int I = 0; I <= maxSpace; I++)
    {
        for (int k = 0; k < 3; k++)
        {
            U[I][k]=U1[I][k];
        }
    }

}