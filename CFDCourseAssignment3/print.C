#include"main.H"
#include<fstream>
#include<iomanip>
void print(const double W[][3])
{

    ofstream foutAll("data/result.dat");
    foutAll.setf(ios::left);
    foutAll.precision(3);
    cout<<"\nprint the rho, u, p, in the data/result.dat, first column is distance"<<endl;
    foutAll<<"x"<<'\t'<<"rho"<<'\t'<<"u"<<'\t'<<"p"<<endl;
    for (int i = 0; i <= maxSpace; i++)
    {
        foutAll<<i*dx<<'\t'<<W[i][0]<<'\t'<<W[i][1]/W[i][0]<<'\t'<<p(W[i])<<endl;
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
