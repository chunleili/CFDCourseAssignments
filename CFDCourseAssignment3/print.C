#include"main.H"
#include<fstream>
#include<iostream>
#include<cstdio>
#include<iomanip>
void print(const double W[][3])
{

    /*
    ofstream foutRho("data/rho.dat");
    ofstream foutP("data/p.dat");
    ofstream foutW("data/u.dat");

    cout<<"\nprint the rho, p, u, in the data/.dat respectively"<<endl;
    for (int i = 0; i <= maxSpace; i++)
    {
        foutRho<<i*dx<<'\t'<<W[i][0]<<endl;
        foutW<<i*dx<<'\t'<<W[i][1]/W[i][0]<<endl;
        foutP<<i*dx<<'\t'<<p(W[i])<<endl;
    }
    */
    ofstream foutAll("data/result.dat");
    //foutAll.setf(ios::fixed);
    //foutAll.setf(ios::showpoint);
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