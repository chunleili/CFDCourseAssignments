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
    	xxx=convert(Q[i][j]);
        lambda=(xxx.VV + xxx.c)*1;//咱先这样,一会再改

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
void FlowField::solve(Mesh& mesh)
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
    double dV;
    XY dS;
    for(unsigned a=0;a<=2;a++)   //a代表荣格库塔法的每一步
        for (unsigned I = 1; I <= maxI-2; I++)//I,J为单元编号
            for(unsigned J = 1; J <= maxJ-2; J++)
            {  
                dV=mesh.getVolume(I,J);
                dS=mesh.getArea(I,J);
                //利用Roe格式计算通量
                Roe(I, J, dS, dV);
                //计算残差
                for (unsigned k=0; k<3; k++)
                {
                    R[k]=Fc[I][J][k].x-Fc[I][J][k].y+Fc[I-1][J][k].x-Fc[I][J-1][k].y;//右侧-左侧+上侧-下侧 
                //利用荣格库塔法计更新流场
                    Q[I][J][k] = Q0[I][J][k] - alpha[a] * dt / dV * R[k];
                }
            }
}

//Roe格式计算对流通量
void FlowField::Roe(Index I,Index J, XY dS, double dV)
{  
    //先定义流场各变量
    ScalarField rho, u, v, VV,  p, H;

    rho[I][J]=Q[I][J][0];
    u[I][J]=Q[I][J][1]/rho[I][J];
    v[I][J]=Q[I][J][1]/rho[I][J];
    VV[I][J]=u[I][J]*u[I][J]+v[I][J]*v[I][J];
    p[I][J]=(GAMMA-1)* (Q[I][J][2] - rho[I][J]*VV[I][J]*0.5);
    H[I][J]=Q[I][J][3]/rho[I][J]+p[I][J]/rho[I][J];

    //然后利用MUSCL分裂这些变量, 并且求出Roe平均量
    Vector QL, QR;
    double rhoL, rhoR, uL, uR, vL, vR, HL, HR, pL, pR;
    Vector F1, F234, F5;
    double LAMC;

    MUSCL(rho, I, J, rhoR, rhoL);
    MUSCL(p,   I, J, pR,   pL);    //注意不要越界~
    MUSCL(u,   I, J, uR,   uL);
    MUSCL(H,   I, J, HR,   HL);
    MUSCL(Q,   I, J, QR,   QL);

    //计算Roe平均量
    const double LR=safeSqrt(rhoL)+ safeSqrt(rhoR);
    const double L=safeSqrt(rhoL)/ LR;//定义两个系数
    const double R=safeSqrt(rhoR)/ LR;

    const double rho_=safeSqrt(rhoL*rhoR);
    const double u_  =uL*L+uR*R;
    const double v_  =vL*L+vR*R;
    const double H_  =HL*L+HR*R;
    const double q_2 =u_*u_+v_*v_;
    const double c_  =safeSqrt((GAMMA-1)*(H_-q_2/2));
    XY V_;
    V_.x=calV_(); V_.y=calV_();

    //计算各个del值
    const double delP  =pR-pL;
    const double delRho=rhoR-rhoL;
    const double delu  =uR-uL;
    const double delV  =delu;

    //求出Roe矩阵相关值
        //delu,delv,delw代表三个分量, delV代表大写V的delta
        //lambda123分别是是V-c, V, V+c
    F1[0]  = lambda1 * (delP-rho_*c_*delV)/(2*c_*c_) * 1;
    F1[1]  = lambda1 * (delP-rho_*c_*delV)/(2*c_*c_) * (u_-c_);
	F1[2]  = lambda1 * (delP-rho_*c_*delV)/(2*c_*c_) * (H_-c_*V_);
  
    F234[0]= lambda2 * ( (delRho-delP/(c_*c_))*1  + 0); 
    F234[1]= lambda2 * ( (delRho-delP/(c_*c_))*u_ 
                   + rho_ *(delu   -  delV)    );
    F234[2]= lambda2 * ( (delRho-delP/(c_*c_))*u_
                              +rho_ *(u_*delu - V_*delV) );
     
    F5[0]  = lambda3 * (delP+rho_*c_*delV)/(2*c_*c_)*1;
    F5[1]  = lambda3 * (delP+rho_*c_*delV)/(2*c_*c_) *(u_+c_);
    F5[0]  = lambda3 * (delP+rho_*c_*delV)/(2*c_*c_)*(H_+c_*V_);

    double delta=0.1*c_; //熵修正, Harten's entropy correction
    LAMC=lambda3;
    if(fabs(LAMC)<=delta)
         LAMC=(LAMC*LAMC+ delta*delta) /(2*delta);
    

    Vector FcR, FcL;
    WToF(WR, FcR);
    WToF(WL, FcL);

    //最后算出单元I的对流通量Fc
    for(unsigned k=0; k<3; k++)
        Fc[I][k]=0.5*( FcR[k]+ FcL[k]- (F1[k]+F234[k]+F5[k]) );
}

//带限制器的MUSCL插值,前两个参数是输入,后两个输出, 系数k^为1/3
void FlowField::MUSCL(Field const U, Index I, Index J, Vector UR, Vector UL)
{
    if(I+2>maxI || I-1<0) {cout<<"\n\nout of bound!!!\n\n";}
    if(J+2>maxJ || J-1<0) {cout<<"\n\nout of bound!!!\n\n";}
    const double epsilon= dx; //限制器参数epsilon与几何尺寸相关
    
    double aR=U[I+2]-U[I+1], bR=U[I+1]-U[I];
    double aL=U[I+1]-U[I],   bL=U[I]  -U[I-1];
    double a, b;
    for(unsigned k=0;k<3;k++)
    {
        a=aR, b=bR;
        double deltaR=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
        a=aL; b=bL;
        double deltaL=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);

        UR[k]=U[I+1][k]-0.5*deltaR;
        UL[k]=U[I]  [k]+0.5*deltaL;
    }
}

//重载用于标量的带限制器的MUSCL插值函数 ,前两个参数是输入,后两个输出, 系数k^为1/3
void FlowField::MUSCL(ScalarField const U, Index I, Index J, double & UR, double & UL)
{
    if(I+2>maxI || I-1<0) {cout<<"\n\nout of bound!!!\n\n";}
    if(J+2>maxJ || J-1<0) {cout<<"\n\nout of bound!!!\n\n";}
    const double epsilon= dx; //限制器参数epsilon与几何尺寸相关

    double aR=U[I+2]-U[I+1], bR=U[I+1]-U[I];
    double aL=U[I+1]-U[I],   bL=U[I]  -U[I-1];

    double a=aR, b=bR;
    is0(( 2*a*a +2*b*b -a*b +3*epsilon));
    double deltaR=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    a=aL, b=bL;
    is0(( 2*a*a +2*b*b -a*b +3*epsilon));
    double deltaL=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    
    UR=U[I+1]-0.5*deltaR;
    UL=U[I]  +0.5*deltaL;
}


void FlowField::calV(Index I, Index J, Mesh &mesh)
{
    double x1,x2,x3,x4, y1,y2,y3,y4;
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
    double u_nx=(y4-y1)*u_/safeSqrt((x1-x4)*(x1-x4)+(y1-y4)*(y1-y4));
    double v_ny=(x2-x1)*v_/safeSqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}