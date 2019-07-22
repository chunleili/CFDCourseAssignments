//这个是一维的!待更改!
#include"main.H"

FlowField::FlowField(int caseNo)//构造函数,选一个方案
{
    switch(caseNo)
    {
        case (1): init1(); break; //设置初场,init1表示入口1.5Ma, init2()表示入口1.8Ma
        case (2): init2(); break;
        default: cout<<"wrong number! exit!"<<endl; exit(1);
    }
}

//LTS=LocalTimeStepping,当地时间步法,返回全局的时间步
double FlowField::LTS()
{
    double min=1e10;
    double dtLocal;
	AERO xxx;
    double lambda;
    
    forAll(
    	xxx=AeroConvert(Q[i][j]);
        lambda=(xxx.Vcv + xxx.c)*1;//咱先这样,一会再改

        dtLocal=CFL*1/lambda;

        if(dtLocal<min)
            min=dtLocal;
    );
    return min;
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

//三阶显式RungeKutta法
void FlowField::solve()
{
    Tensor Fc;   //对流通量, 右侧x和上侧y的, 向外为正
    const double alpha[3]={0.1918, 0.4929, 1.0};
    //先定义Q0,用于保存原始的Q
    Field Q0;
    forEach(
            Q0[i][j][k] = Q[i][j][k];
    );
    
    //后面每一步都先计算残差, 后根据RK公式更新W
    Vector R;//R for Residual, 残差
    //double dV;
    //XY dS;
    for(unsigned a=0;a<=2;a++)   //a代表荣格库塔法的每一步
        for (unsigned I = 1; I <= maxI-2; I++)//I,J为单元编号
            for(unsigned J = 1; J <= maxJ-2; J++)
            {  
                //利用Roe格式计算通量
                Roe roe(I, J);
                //计算残差
                for (unsigned k=0; k<3; k++)
                {   
                    Fc[I][J][k]=roe.getFlux();
                    //R代表离开单元格的通量的矢量和, =右侧+左侧+上侧+下侧
                    //Fc[I][J].x/y分别代表左侧与下侧通量, 由于与坐标系正方向反向, 故取负号
                    //右侧与上侧通量分别由临近单元格左侧与下侧取负号得来, 取两次负号变为正
                    R[k]=-Fc[I+1][J][k].x + Fc[I][J][k].x - Fc[I][J+1][k].y + Fc[k][I][J].y;
                     
                    //利用荣格库塔法计更新流场
                    Q[I][J][k] = Q0[I][J][k] - alpha[a] * dt / dV * R[k];
                }
            }
}