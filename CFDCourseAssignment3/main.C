#include "main.H"
#include <ctime>
#include <fstream>
#include <iostream>

void init(Field W);
void print( const Field W);
double LTS(Field W);
void solver3(Field W, const double dt);
//solver1是Jameson格式搭配3阶荣格库塔;
//solver2是MacCormack法
//solver3是Roe格式搭配3阶荣格库塔法

int main()
{
    Field W;
    //采用单元中心法, W[0] 与W[maxSpace]是边界;
    init(W);

    double dtGlobal;
    int step = 0;
    for (double t = 0; t <= stopTime && step<100 ; t += dtGlobal)
    {
        dtGlobal=LTS(W);
        solver3(W,dtGlobal); //当前用的是solver3
        
        step++;
        cout.width(8);
        cout << "dtGlobal= " << dtGlobal
             << "\t step = " << step
             << "\t time = " << t << endl;
    }
    print(W);
    return 0;
}

//						declear the funcs								 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void init(Field W)
{
    cout<<"\nInitializing the fields..."<<endl;
    double u1=0.0, rho1=1.0, p1=1.0;
    double u2=0.0, rho2=0.125, p2=0.1;

    for(int i=0; i< maxSpace/2; i++)
    {
        W[i][0]=rho1;
        W[i][1]=u1*rho1;
        W[i][2]=p1/(GAMMA-1) + 0.5*rho1*u1*u1;
    }

    for(int i=maxSpace/2; i<= maxSpace; i++)
    {
        W[i][0]=rho2;
        W[i][1]=u2*rho2;
        W[i][2]=p2/(GAMMA-1) + 0.5*rho2*u2*u2;
    }
}

//0方向梯度边界条件, 边界上的点直接等于内部的点 
void zeroGradBC(Field W)
{
    for (int k = 0; k < 3; k++)
    {
        W[0][k] = W[1][k];
        W[maxSpace][k] = W[maxSpace - 2][k];
        W[maxSpace][k] = W[maxSpace - 1][k];
    }
}

void WToF(Vector W, Vector F)
{
    double p0=calPressure(W);
    double u=W[1]/W[0];
    F[0] = W[1];
    F[1] = W[0] * u * u + p0;
    F[2] = (W[2] + p0) * u;
}

void print(const Field W)
{

    ofstream foutAll("result.dat");
    cout<<"\nprint the rho, u, p, in the \"result.dat\", first column is distance"<<endl;
    foutAll.setf(ios::left);
    foutAll.width(7);
    foutAll<<"x"<<'\t'<<"rho"<<'\t'<<"u"<<'\t'<<"p"<<endl;
    for (int i = 0; i <= maxSpace; i++)
    {
        foutAll.width(7);
        foutAll.precision(5);
        foutAll<<i*dx<<'\t'<<W[i][0]<<'\t'<<W[i][1]/W[i][0]<<'\t'<<calPressure(W[i])<<endl;
    } 
}


//LTS=LocalTimeStepping,当地时间步法,返回全局的时间步
double LTS(Field W)
{
    double min=1e10;
    double dtLocal;
    for (int i = 0; i <=maxSpace ; i++)
    {
        dtLocal=CFL*dx/lambda(W[i]);//注意此i的生存期
        if(dtLocal<min)
        {
            min=dtLocal;
        }
    }
    return min;
}
