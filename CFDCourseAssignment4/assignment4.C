/***************************include & namespace  **************************/
#include"main.H"
#include<cstdlib>
#include<ctime>
#include<fstream>
/***************************declare static funcs  **************************/
static inline void chooseCaseAndInit();
static inline void runInfo(double t, double dt, int step, double residual);
static inline void ifConverge(double residual, int step, Field Q);

/********************************main()**************************************/
int main()
{  
    genMesh();              //生成网格
    
    Field Q;
    Tensor F;

    chooseCaseAndInit();

    int step=0; double residual=0.0; double t=0.0;
    for (double dt=0.001; t <= STOP_TIME && step<=STOP_STEP; t+=dt, step++)
    {
        //dt = localTime();       //当地时间步法，注意取全域最大当地时间步
        solve();                  //求解

        runInfo(t,dt,step, residual);
        ifConverge(residual, step, Q);
    }
    return 0;
}

/********************************funcs define**************************************/
/********************************genMesh()**************************************/
void genMesh()
{
    const int block1=100, block2=100, block3=maxI-block1-block2;
    const double dx=1.0/500;

    MeshPoint mesh;

    double dy=0.3/maxJ;
    for(int i=0; i<block1; i++)
    {
        for(int j=0;j<=maxJ; j++)
        {
            mesh[i][j].x=i*dx;
            mesh[i][j].y=j*dy;
        }
    }
    
    for(int i=block1; i<block1+block2; i++)
    {
        double h=0.25*i*dx-0.05;
        dy=(0.3-h)/maxJ;
        for(int j=0;j<=maxJ; j++)
        {
            mesh[i][j].x=i*dx;
            mesh[i][j].y=h+j*dy;
        }
    }

    dy=(0.3-0.05)/maxJ;
    for(int i=block1+block2; i<=maxI; i++)
    {
        for(int j=0;j<=maxJ; j++)
        {
            mesh[i][j].x=i*dx;
            mesh[i][j].y=0.05+j*dy;
        }
    }

    fstream fout("mesh.dat");
    fout
	<<"Title=\"Mesh\""<<endl
	<<"Variables=\"x\",\"y\""<<endl
	<<"Zone i="<<maxI+1<<", j="<<maxJ+1<<", f=point"<<endl;
    for(int j=0; j<=maxJ; j++)
    {
        for(int i=0; i<=maxI; i++)
        {
            fout<< mesh[i][j].x<<" "<<mesh[i][j].y<<endl;
        }
    }
    cout<<"mesh is written in \"mesh.dat\""<<endl;

}
/********************************init1()**************************************/
void init1()
{
    
}

void init2()
{

}


/********************************other funcs**************************************/
double localTime()
{
    double globalTime=0;
    return globalTime;
}


void print(Field aField, fstream &fout)
{
    
    for(int j=0; j<=maxJ; j++)
    {
        for(int i=0; i<=maxI; i++)
        {
            for(int k: aField[i][j])
                fout<< aField[i][j][k]<<" ";
            fout<<endl;
        }
    }
    cout<<"\nWall time = "<<(double)clock()/CLOCKS_PER_SEC<<" s"<<endl;
}

/********************************static funcs**************************************/
static inline void runInfo(double t, double dt, int step, double residual)
{
    cout<<"step= "   <<step
        <<"\tt="       <<t 
		<<"\tdt= "     <<dt
        <<"\tresidual="<<residual
        <<endl;
    fstream fout("residual.dat");
    fout<<"\tt="        <<t 
        <<"\tresidual=" <<residual
        <<endl;
}

static inline void ifConverge(double residual, int step, Field Q)
{
    if(residual<RESIDUAL_LIMIT)
    {
        cout<<"\nConverged!"<<endl
	    	<<"\nmaxI= "<<maxI<<"\tmaxJ= "<<maxJ<<endl
	    	<<"Final residual ="<<residual<<endl
	    	<<"Loop times = "<<step<<endl
	    	<<"Final field has been written into \"result.dat\" \n"
	    	<<"History of residual has been written into \"residual.dat\"";
        fstream fout("result.dat");
        print(Q, fout);
	    exit(0);
        
    }
    else if(step>=STOP_STEP)
	{
		cout<<"\nFail to converge! in "<<STOP_STEP<<" times\n"
		<<"Final residual ="<<residual<<"\n\n\n";
        cout<<"\nWall time = "<<(double)clock()/CLOCKS_PER_SEC<<" s"<<endl;
	}
    cout<<"\nWall time = "<<(double)clock()/CLOCKS_PER_SEC<<" s"<<endl;

}

static inline void chooseCaseAndInit()
{
    cout<<"please input case Number: 1 for inlet 1.8Ma; 2 for 1.5Ma"<<endl;
    int caseNo; cin>>caseNo;
    switch(caseNo)
    {
        case (1): init1(); break; //设置初场,init1表示入口1.5Ma, init2()表示入口1.8Ma
        case (2): init2(); break;
        default: cout<<"wrong number! exit!"<<endl; exit(1);
    }
}