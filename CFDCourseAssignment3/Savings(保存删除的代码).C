//这个文件不会被编译,仅仅用来保存不再使用但是不舍得丢掉的代码

//一阶精度的三阶显式RungeKutta法
void RungeKuttaTime(double W[][3], const double dt, double R[][3])
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

void printW(const double W[][3])
{
    ofstream fout("data/vecterW.dat");
    cout<<"\nprint W[][3]\n";
    for (unsigned i =0 ; i <=maxSpace ; i++)
    {
        fout<<W[i][0]<<"\t"<<W[i][1]<<'\t'<<W[i][2]<<endl;
    }
}

void printF(const double W[][3])
{
    ofstream fout("data/vecterF.dat");
    cout<<"\nprint F\n";
    for (unsigned i =0 ; i <=maxSpace ; i++)
    {
        double F[3];
        WToF(W[i], F);
        fout<<F[0]<<"\t"<<F[1]<<'\t'<<F[2]<<endl;
    }
}


double localTime(const double W[][3], const int I)
{
    return CFL*dx/lambda(W[I]);
}

//不带限制器的MUSCL插值,前两个参数是输入,后两个输出, 系数取1/3, 
void MUSCL(Field const U, int const I, Vector UR, Vector UL)
{
    const double kk=1.0/3;   //epsilon=1.0;

    for(int k=0;k<3;k++)
    {
        double delBUI =U[I  ][k]-U[I-1][k];
        double delFUI1=U[I+2][k]-U[I+1][k];
        double delBUI1=U[I+1][k]-U[I  ][k];
        double delFUI =U[I+1][k]-U[I  ][k];

        UR[k]=U[I+1][k]-0.25*((1+kk)* delBUI1[k] + (1-kk)* delFUI1[k]);
        UL[k]=U[I]  [k]+0.25*((1+kk)* delFUI [k] + (1-kk)* delBUI [k]);
    }
}

//重载用于标量的不带限制器的MUSCL插值函数 ,前两个参数是输入,后两个输出
void MUSCL(double const U, int const I, double & UR, double & UL)
{
    const double kk=1.0/3;   //epsilon=1.0;

    double delBUI =U[I  ]-U[I-1];
    double delFUI1=U[I+2]-U[I+1];
    double delBUI1=U[I+1]-U[I  ];
    double delFUI =U[I+1]-U[I  ];
    UR=UI1-0.25*((1+kk)* delBUI1 + (1-kk)* delFUI1);
    UL=UI +0.25*((1+kk)* delFUI  + (1-kk)* delBUI );
    
}