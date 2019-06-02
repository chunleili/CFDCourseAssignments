#include "main.H"
#include <ctime>

int main()
{
    double W[maxSpace + 1][3] = {0}, R[maxSpace + 1][3] = {0};
    init(W);
    double dtGlobal = 1111;
    int timeStep = 0;
    for (double t = 0; t <= stopTime && timeStep < maxTime; t += dtGlobal)
    {
        double dtLocal;
        for (int I = 0; I <= maxSpace; I++)
        {
            dtLocal = localTime(W, I);
            TIME_DIS(W, dtLocal, R, I);
            if (dtGlobal > dtLocal)
                dtGlobal = dtLocal;
        }

        timeStep++;
        cout << "dtGlobal= " << dtGlobal;
        cout << "\t time step = " << timeStep;
        cout << "\t time = " << t << endl;
    }

    print(W);
    printW(W);
    cout << "\nFinal Wall time = " << (double)clock() / CLOCKS_PER_SEC << " s" << endl;

    return 0;
}