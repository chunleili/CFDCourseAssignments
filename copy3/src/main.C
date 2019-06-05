#include "main.H"

int main()
{
    double W[maxSpace + 3][3] , R[maxSpace + 3][3];
    //采用单元中心法, 实际存储网格从W[1] 到W[maxSpace]; W[0]与W[maxSpace+1],W[maxSpace+2]是边界虚网格
    init(W);
    double dtGlobal;
    int timeStep = 0;
    for (double t = 0; t <= stopTime && timeStep<10000 ; t += dtGlobal)
    {
        dtGlobal=calDtGlobal(W);
        
        RungeKuttaTime(W, dtGlobal, R);
        
        timeStep++;
        cout << "dtGlobal= " << dtGlobal << "\t time step = " << timeStep << "\t time = " << t << endl;
    }
    print(W);
    return 0;
}
