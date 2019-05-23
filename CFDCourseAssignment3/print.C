#include"main.H"
void print( double U[][maxSpace+1])
{
    ofstream foutRho("data/rho.dat");
    ofstream foutP("data/p.dat");
    ofstream foutU("data/u.dat");

    double p;
    cout<<"\nprint the rho, p, u, in the data/*.dat respectively"<<endl;
    for (int i = 0; i <= maxSpace; i++)
    {
        foutRho<<i*dx<<'\t'<<U[0][i]<<endl;
        foutU<<i*dx<<'\t'<<U[1][i]/U[0][i]<<endl;
        
        p=(GAMMA-1) * ( U[2][i] - 0.5*U[1][i]*U[1][i] / U[0][i] );
        foutP<<i*dx<<'\t'<<p<<endl;
    }
    
}