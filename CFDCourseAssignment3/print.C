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
        foutRho<<i*dx<<'\t'<<U[i][0]<<endl;
        foutU<<i*dx<<'\t'<<U[i][1]/U[i][0]<<endl;
        
        p=(GAMMA-1) * ( U[i][2] - 0.5*U[i][1]*U[i][1] / U[i][0] );
        foutP<<i*dx<<'\t'<<p<<endl;
    }
    
}