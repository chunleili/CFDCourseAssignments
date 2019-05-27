#include"main.H"
#include<cstdlib>
#include<ctime>

int main()
{
    double U[maxSpace+1][3]={0};

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
    for( int timeStep =0; timeStep<=maxTime; timeStep++)
    {
        cout<<"Current time step = "<<timeStep;
        cout<<"\t Current physical time = "<<timeStep*dt<<endl;
        MacCormackSolver(U);
    //    scalarJSTSolver(U);
    }


    cout<<"\nCalculate over, printing results..."<<endl;
    print(U);
    printU(U);
    cout<<"\nFinal Wall time = "<<(double)clock()/CLOCKS_PER_SEC<<" s"<<endl;
	cout<<"Final CPU time= "<<(double)clock()<<" s"<<endl;

    return 0;
}