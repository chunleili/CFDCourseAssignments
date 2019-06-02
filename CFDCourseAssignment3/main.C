#include "main.H"
#include <ctime>
#include<cstdlib>
double calDtGlobal(double W[][3]);
int main()
{
    double W[maxSpace + 1][3] , R[maxSpace + 1][3] ;
    init(W);
    double dtGlobal;
    int timeStep = 0;
    for (double t = 0; t <= stopTime ; t += dtGlobal)
    {
        dtGlobal=calDtGlobal(W);
        for (int I = 1; I <= maxSpace-2; I++)
        {
            TIME_DIS(W, dtGlobal, R, I);
        }
        zeroGradBC(W);
        timeStep++;
        cout << "dtGlobal= " << dtGlobal;
        cout << "\t time step = " << timeStep;
        cout << "\t time = " << t << endl;
    }
    print(W);
    return 0;
}