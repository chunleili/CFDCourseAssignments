#include"main.H"
#include<cstdlib>
#include<ctime>
#define CONV scalarJSTConv
//  #define CONV MacCormackConv;
//  #define CONV RoeConv;
int main()
{
    double U[maxSpace+1][3]={0}, R[maxSpace+1][3]={0};

    cout<<"\nInitializing the fields..."<<endl;
    init(U);
    printU(U);
    print(U);

    cout<<"\n Ready iterating?(y/n)"<<endl;
    char key;
    cin>>key;
    if (key == 'n' || key == 'N')
    {
        cout<<"User terminated!"<<endl;
        exit(0);
    }

    cout<<"\nMarching in the time steps"<<endl;
    double dt=localTime(U), Rmax=0;
    int timeStep=0;
    for( double t =0; t<=stopTime && timeStep<maxTime; t+=dt)
    {
        cout<<"dt= "<<dt;
        cout<<"\t time step = "<<timeStep;
        cout<<"\t time = "<<t;
     //   MacCormackConv(U,dt,R);
        Conv(U,dt,R);
     //   RoeConv(U,dt,R);
        EulerFTime(U,dt,R);
        cout<<"\t Residual = "<<Rmax<<endl;
        dt=localTime(U);
        timeStep++;
    }

    cout<<"\nCalculate over, printing results..."<<endl;
    print(U);
    printU(U);
    cout<<"\nFinal Wall time = "<<(double)clock()/CLOCKS_PER_SEC<<" s"<<endl;
	cout<<"Final CPU time= "<<(double)clock()<<" s"<<endl;

    return 0;
}