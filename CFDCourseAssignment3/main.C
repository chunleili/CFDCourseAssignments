#include "main.H"
#include <ctime>

int main()
{
    double W[maxSpace + 3][3] , R[maxSpace + 3][3];
    //采用单元中心法, 实际存储网格从W[1] 到W[maxSpace]; W[0]与W[maxSpace+1],W[maxSpace+2]是边界虚网格
    init(W);
    double dtGlobal;
    int timeStep = 0;
    for (double t = 0; t <= stopTime && timeStep<1 ; t += dtGlobal)
    {
        dtGlobal=calDtGlobal(W);
        
        EulerFTime(W, dtGlobal, R);
        
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