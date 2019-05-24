#include"main.H"
void print( double U[][3])
{
    ofstream foutRho("data/rho.dat");
    ofstream foutP("data/p.dat");
    ofstream foutU("data/u.dat");

    cout<<"\nprint the rho, p, u, in the data/*.dat respectively"<<endl;
    for (int i = 0; i <= maxSpace; i++)
    {
        foutRho<<i*dx<<'\t'<<U[i][0]<<endl;
        foutU<<i*dx<<'\t'<<U[i][1]/U[i][0]<<endl;
        foutP<<i*dx<<'\t'<<p(U[i])<<endl;
    }
    
}