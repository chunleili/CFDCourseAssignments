#include"main.H"
#include<cstdlib>
#include<ctime>
#include<iostream>
#include<cstdio>
#include<iomanip>

int main()
{
    double W[maxSpace+1][3]={0}, R[maxSpace+1][3]={0};

    cout<<"\nInitializing the fields..."<<endl;
    init(W);
    printW(W);
    print(W);

/*    cout<<"\n Ready iterating?(y/n)"<<endl;
    char key;
    cin>>key;
    if (key == 'n' || key == 'N')
    {
        cout<<"terminated!"<<endl;
        exit(0);
    }*/

    double dtGlobal=1111;
    int timeStep=0;
  //  double dtLocalField[maxSpace+1];

   /* for (int I = 0; I <=maxSpace; I++)
    {
        dtLocal=localTime(W,I);
        if(dtLocal<dtGlobal)
            dtGlobal=dtLocal;
    }*/
    
    for( double t =0; t<=stopTime && timeStep<maxTime; t+=dtGlobal)
    {
        
        double dtLocal;
        for (int I = 0; I <=maxSpace; I++)
        {
      //      double dtLocal;
            dtLocal=localTime(W,I);
        //    dtLocalField[I]=dtLocal;

            TIME_DIS(W, dtLocal, R, I);

            if (dtGlobal>dtLocal)
                dtGlobal=dtLocal;
        }
        
//        dtGlobal=minOverScalarField(dtLocalField);
        timeStep++;

        cout<<"dtGlobal= "<<dtGlobal;
        cout<<"\t time step = "<<timeStep;
        cout<<"\t time = "<<t<<endl;
    }

    cout<<"\nCalculate over, printing results..."<<endl;
    print(W);
    printW(W);
    cout<<"\nFinal Wall time = "<<(double)clock()/CLOCKS_PER_SEC<<" s"<<endl;

    return 0;
}