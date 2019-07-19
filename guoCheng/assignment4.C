/***************************include & namespace  **************************/
#include"main.H"
#include<cstdlib>
#include<cstdio>
#include<ctime>
#include<fstream>
#include<iomanip>
/***************************declare funcs  **************************/
static  void chooseCaseAndInit(Field Q);
static  void runInfo(int step, double residual);

static  void init1(Field Q);
static  void init2(Field Q);
static  void  print(Field Q); 
static  double LTS();
/********************************main()**************************************/
int main()
{  
    FlowField field1;

    double residual=1.0;
    for (int step=0;  step<=STOP_STEP; step++)
    {
        dt = field1.LTS();       //当地时间步法，注意取全域最大当地时间步
        solve();                  //求解
        
        runInfo(step, residual);
    }

    print(Q);
    return 0;
}

/********************************funcs define**************************************/
/********************************Mesh**************************************/
Mesh::Mesh()
{
    const int block1=100, block2=100;
    const double dx=1.0/500;

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
}
void Mesh::print()
{
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

//计算单元体体积, 厚度取为1
double Mesh::getVolume(Index I, Index J)
{
    double volume, x1,x2,x3,x4, y1,y2,y3,y4; 
    //从左下开始逆时针编号,左下点代表1,右下2,右上3,左上4
    //左下代表本单元格坐标
    x1=mesh[I  ][J  ].x;
    x2=mesh[I+1][J  ].x;
    x3=mesh[I+1][J+1].x;
    x4=mesh[I  ][J+1].x;

    y1=mesh[I  ][J  ].y;
    y2=mesh[I+1][J  ].y;
    y3=mesh[I+1][J+1].y;
    y4=mesh[I  ][J+1].y;

    volume=0.5*( (x1-x3)*(y2-y4)+ (x4-x2)*(y1-y3) );
    return volume;
}

//计算单元体右上侧面积
XY Mesh::getArea( Index I, Index J)
{
    double x1,x2,x3,x4, y1,y2,y3,y4;
    XY S; 
    //从左下开始逆时针编号,左下点代表1,右下2,右上3,左上4
    //左下代表本单元格坐标
    x1=mesh[I  ][J  ].x;
    x2=mesh[I+1][J  ].x;
    x3=mesh[I+1][J+1].x;
    x4=mesh[I  ][J+1].x;

    y1=mesh[I  ][J  ].y;
    y2=mesh[I+1][J  ].y;
    y3=mesh[I+1][J+1].y;
    y4=mesh[I  ][J+1].y;

    S.x=sqrt( (x3-x2)*(x3-x2)+(y3-y2)*(y3-y2) );//右侧面积
    S.y=sqrt( (x4-x3)*(x4-x3)+(y4-y3)*(y4-y3) );//上侧面积
    
    return S;
}


/********************************FlowField**************************************/
FlowField::FlowField(int caseNo)
{
    cout<<"please input case Number: 1 for inlet 1.8Ma; 2 for 1.5Ma"<<endl;
    int caseNo; cin>>caseNo;
    switch(caseNo)
    {
        case (1): init1(Q); break; //设置初场,init1表示入口1.5Ma, init2()表示入口1.8Ma
        case (2): init2(Q); break;
        default: cout<<"wrong number! exit!"<<endl; exit(1);
    }
}

void FlowField::init1()
{
    //初始全部给1.8Ma对应的速度, 压力给大气压101325, 静温300
    //声速340.68
    //速度u=613, v=0, 密度rho=p/RT=1.17
    forAll(
        Q[i][j][0]=1.17;
        Q[i][j][1]=717.21; //1.17*613
        Q[i][j][2]=0;
        Q[i][j][1]=215713; // 287/0.4*300+613;
    );
}

void FlowField::init2()
{
    //初始全部给1.5Ma对应的速度, 压力给大气压101325, 静温300, 
    //速度u=511.02, v=0, 密度rho=p/RT=1.17
    forAll(
        Q[i][j][0]=1.17;
        Q[i][j][1]=597.89; //1.17*511.02
        Q[i][j][2]=0;
        Q[i][j][1]=215713; // 287/0.4*300+613;
    );
}

void FlowField::BC1()
{

}

void FlowField::BC2()
{

}

/********************************utility**************************************/
//气动参数转换
AERO convert(Vector vec)
{
    AERO ff;
    double rho,u,v,VV,p,T,c,Ma;
    ff.rho=vec[0];     
    ff.u  =vec[1] / vec[0];
    ff.v  =vec[2] / vec[0];
    ff.VV =u*u+v*v;
    ff.p  =(GAMMA-1)*(vec[3]-0.5*rho*(u*u+v*v));
    ff.T  =p/(rho*287);
    ff.c  =20.04*sqrt(T);
    ff.Ma =(u*u+v*v)/c;
    return ff;
}

//Q与Fc之间的转化
void QToF(Vector W, Vector F)
{
    double   u=W[1]/W[0];
    double   p0=WToP(W);
    is0(W[0]);
    F[0] = W[1];
    F[1] = W[0] * u * u + p0;
    F[2] = (W[2] + p0) * u;
}



double safeSqrt(double xx)
{
    if(xx<0)
    {
        cout<<"\n\n\nerro: sqrt NEGATIVE!!!!\n\n\n"<<endl;
        return -1e100000;
    }
    else
        return sqrt(xx);
}

/********************************other funcs**************************************/
//LTS=LocalTimeStepping,当地时间步法,返回全局的时间步
double LTS(FlowField &field)
{
    double min=1e10;
    double dtLocal;
    double lambda, c, p, rho, u;
    forAll(
        rho=field.Q[i][j][0];
        u=field.Q[i][j][1]/rho;
        p=QToP(field.Q[i][j]);
        c= safeSqrt(GAMMA*p/rho);
        lambda=fabs(u) + c;

        dtLocal=CFL*dx/lambda;

        if(dtLocal<min)
            min=dtLocal;
    );
    return min;
}


static void print(Field Q, MeshPoint mesh)
{
    fstream fout("result.dat");
    cout << "\nprint the field in the \"result.dat\"" << endl;
    fout.setf(ios::left);
    fout << setw(7) << "x"   << '\t' << setw(7) << "y" << '\t' 
         << setw(7) << "rho" << '\t' << setw(7) << "u" << '\t' 
         << setw(7) << "v"   << '\t' << setw(7) << "p" << '\t'
         << setw(7) << "p"   << '\t' << setw(7) << "T" << '\t'
         << setw(7) << "c"   << '\t' << setw(7) << "Ma"<< endl;
    
    const double dx=1.0/500;
    double y,rho,u,v,p,T,c,Ma;
    for(unsigned j=0; j<=maxJ; j++)
    { 
        for (unsigned i = 0; i <= maxI; i++)
        {   
            rho=Q[i][j][0];     
            u=Q[i][j][1] / Q[i][j][0];
            v=Q[i][j][2] / Q[i][j][0];
            p=(GAMMA-1)*(Q[i][j][3]-0.5*rho*(u*u+v*v));
            T=p/(rho*287);
            c=1.1832*sqrt(p/Q[i][j][0]);
            Ma=(u*u+v*v)/c;
            fout.precision(5);
            fout << setw(7) << mesh[i][j].x  << '\t'
                 << setw(7) << mesh[i][j].y  << '\t' 
                 << setw(7) << rho           << '\t'
                 << setw(7) << u             << '\t'
                 << setw(7) << v             << '\t'
                 << setw(7) << p             << '\t'
                 << setw(7) << T             << '\t'
                 << setw(7) << c             << '\t'
                 << setw(7) << Ma            << endl;
        }
    }
}

static void runInfo(int step, double residual)
{
    cout.setf(ios::left);
    cout.width(15);
    cout.precision(6);
    cout<<"step= "   <<step
        <<"\tresidual="<<residual
        <<endl;
    
    fstream fout("residual.dat");
    fout.setf(ios::left);
    fout.width(15);
    fout.precision(6);
    fout<<"\tstep="     <<step 
        <<"\tresidual=" <<residual
        <<endl;
}
