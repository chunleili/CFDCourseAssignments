#include "main.H"
#include <ctime>
int main()
{
    double W[maxSpace + 3][3] , R[maxSpace + 3][3];
    //采用单元中心法, 实际存储网格从W[1] 到W[maxSpace]; W[0]与W[maxSpace+1],W[maxSpace+2]是边界虚网格
    init(W);
    double dtGlobal;
    int timeStep = 0;

    for (double t = 0; t <= stopTime && timeStep<100 ; t += dtGlobal)
    {
        dtGlobal=calDtGlobal(W);
        solver3(W, dtGlobal, R);//solver1是Jameson格式搭配3阶荣格库塔; 
                                //solver2是MacCormack法
                                //solver3是Roe格式搭配3阶荣格库塔法
        
        timeStep++;
        cout << "dtGlobal= " << dtGlobal;
        cout << "\t time step = " << timeStep;
        cout << "\t time = " << t << endl;
    }
    print(W);
    printW(W);
    printF(W);
    return 0;
}

//						declear the funcs								 //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void init(double W[][3])
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

    for(int i=maxSpace/2; i<= maxSpace+2; i++)
    {
        W[i][0]=rho2;
        W[i][1]=u2*rho2;
        W[i][2]=p2/(GAMMA-1) + 0.5*rho2*u2*u2;
    }
}

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

void WToF(double const W[3], double F[3])
{
    double p0=calPressure(W);
    double u=W[1]/W[0];
    F[0] = W[1];
    F[1] = W[0] * u * u + p0;
    F[2] = (W[2] + p0) * u;
}

void print(const double W[][3])
{

    ofstream foutAll("data/result.dat");
    foutAll.setf(ios::left);
    foutAll.precision(4);
    cout<<"\nprint the rho, u, p, in the data/result.dat, first column is distance"<<endl;
    foutAll<<"x"<<'\t'<<"rho"<<'\t'<<"u"<<'\t'<<"p"<<endl;
    for (int i = 0; i <= maxSpace+2; i++)
    {
        foutAll<<i*dx<<'\t'<<W[i][0]<<'\t'<<W[i][1]/W[i][0]<<'\t'<<calPressure(W[i])<<endl;
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

double calDtGlobal(double W[][3])
{
    double min=1e10;
    double dtLocal;
    for (int i = 0; i <=maxSpace ; i++)
    {
        dtLocal=localTime(W,i);//注意此i的生存期
        if(dtLocal<min)
        {
            min=dtLocal;
        }
    }
    return min;
}

void EulerFTime(double W[][3], const double dt, double R[][3])
{
    scalarJSTConv(W, R);

    for (int I = 1; I <= maxSpace; I++)
    {
        for (int k = 0; k < 3; k++)
        {
            W[I][k] = dt * (-1 / dx) * R[I][k] + W[I][k];
        }
    }
}