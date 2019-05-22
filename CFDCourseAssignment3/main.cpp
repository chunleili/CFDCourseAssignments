#include"main.H"
#include<cstdlib>
#include<ctime>



int main()
{
    double U[3][maxSpace+1], F[3][maxSpace+1];

    cout<<"\nInitializing the fields..."<<endl;
    init(U,F);

    cout<<"\n Ready iterating?(y/n)"<<endl;
    char key;
    cin>>key;
    if (key == 'n' || key == 'N')
    {
        cout<<"User terminated!"<<endl;
        exit(0);
    }

    cout<<"\nMarching in the time steps"<<endl;
    for( unsigned int timeStep =0; timeStep<=maxTime; timeStep++)
    {
        cout<<"Current time step = "<<timeStep;
        cout<<"\t Current physical time = "<<timeStep*dt;
        MacCormackSolver(U, F);

        cout<<"\tWall time = "<<(double)clock()/CLOCKS_PER_SEC<<" s";
	//    cout<<"\tCPU time= "<<(double)clock()<<" s"<<endl;
    }


    cout<<"\nCalculate over, printing results..."<<endl;
    print(U);

    cout<<"\nFinal Wall time = "<<(double)clock()/CLOCKS_PER_SEC<<" s"<<endl;
	cout<<"Final CPU time= "<<(double)clock()<<" s"<<endl;

    return 0;
}