#include"main.H"
#include<cstdlib>
#include<ctime>
#define CONV scalarJSTConv
//  #define CONV MacCormackConv
//  #define CONV RoeConv
#define TIME_DIS EulerFTime
int main()
{
    double W[maxSpace+1][3]={0}, R[maxSpace+1][3]={0};

    cout.width(15);
    cout.setf(ios::left);
    cout<<"\nInitializing the fields..."<<endl;
    init(W);
    printW(W);
    print(W);

    cout<<"\n Ready iterating?(y/n)"<<endl;
    char key;
    cin>>key;
    if (key == 'n' || key == 'N')
    {
        cout<<"terminated!"<<endl;
        exit(0);
    }

    double dt=localTime(W);
 //   double Rmax=0;
    int timeStep=0;
    for( double t =0; t<=stopTime && timeStep<maxTime; t+=dt)
    {
        cout<<"dt= "<<dt;
        cout<<"\t time step = "<<timeStep;
        cout<<"\t time = "<<t<<endl;
        CONV(W,dt,R);
        zeroGradBC(W);
        TIME_DIS(W,dt,R);
//        cout<<"\t Residual = "<<Rmax<<endl;
        dt=localTime(W);
        timeStep++;
    }

    cout<<"\nCalculate over, printing results..."<<endl;
    print(W);
    printW(W);
    cout<<"\nFinal Wall time = "<<(double)clock()/CLOCKS_PER_SEC<<" s"<<endl;

    return 0;
}