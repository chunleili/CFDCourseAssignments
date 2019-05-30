#include"main.H"

void print(const double U[][3])
{

    /*
    ofstream foutRho("data/rho.dat");
    ofstream foutP("data/p.dat");
    ofstream foutU("data/u.dat");

    cout<<"\nprint the rho, p, u, in the data/.dat respectively"<<endl;
    for (int i = 0; i <= maxSpace; i++)
    {
        foutRho<<i*dx<<'\t'<<U[i][0]<<endl;
        foutU<<i*dx<<'\t'<<U[i][1]/U[i][0]<<endl;
        foutP<<i*dx<<'\t'<<p(U[i])<<endl;
    }
    */
    ofstream foutAll("data/result.dat");
    foutAll.setf(ios::fixed);
    foutAll.setf(ios::showpoint);
    foutAll.precision(3);
    cout<<"\nprint the rho, u, p, in the data/result.dat, first column is distance"<<endl;
    foutAll<<"x"<<'\t'<<"rho"<<'\t'<<"u"<<'\t'<<"p"<<endl;
    for (int i = 0; i <= maxSpace; i++)
    {
        foutAll<<i*dx<<'\t'<<U[i][0]<<'\t'<<U[i][1]/U[i][0]<<'\t'<<p(U[i])<<endl;
    } 
}

void printU(const double U[][3])
{
    ofstream fout("data/vecterU.dat");
    cout<<"\nprint U[][3]\n";
    for (unsigned i =0 ; i <=maxSpace ; i++)
    {
        fout<<U[i][0]<<"\t"<<U[i][1]<<'\t'<<U[i][2]<<endl;
    }
    
}