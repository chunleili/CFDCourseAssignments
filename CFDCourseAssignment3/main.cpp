#include"main.H"

int main()
{
    double U[3][maxSpace+2], F[3][maxSpace+2];

    init(U,F);

    for( unsigned int timeStep =0; timeStep<=maxTime; timeStep++)
    {
        MacCormackSolver(U, F);
    }
}