#define SOLVER3
#ifdef SOLVER3

#include"main.H"
void Roe(Field W, int const I, Field Fc);
void MUSCL(Field const U, int const I, Vector UR, Vector UL);
void MUSCL(ScalarField const U, int const I, double & UR, double & UL);

//一阶精度的三阶显式RungeKutta法
void solver3(Field W, const double dt)
{
    Field Fc;
    const double alpha[3]={0.1481, 0.4, 1.0};
    //先定义W0,用于保存原始的W
    Field W0;
    for (int I = 0; I <= maxSpace; I++)
        for (int k = 0; k < 3; k++)
            W0[I][k] = W[I][k];
    
    //后面每一步都先计算残差, 后根据RK公式更新W
    Vector R;//R for Residual, 残差
    for(int a=0;a<=2;a++)   //a代表荣格库塔法的每一步
    {
        
        for (int I = 0; I <= maxSpace; I++)
        {  
            //利用Roe格式计算通量
            Roe(W, I, Fc);
            //计算残差
            for (int k = 0; k < 3; k++)
                R[k]=Fc[I+1][k]-Fc[I][k];
    
            //利用荣格库塔法计更新流场
            for (int k = 0; k < 3; k++)
                W[I][k] = W0[I][k] - alpha[a] * dt / dx * R[k];
        }
    }
}

//Roe格式计算对流通量
void Roe(Field W, int const I, Field Fc)
{  
    //先定义流场各变量
    ScalarField u, V, rho, p, H;

    rho[I]=W[I][0];
    u[I]=W[I][1]/rho[I];
    V[I]=u[I];
    H[I]=W[I][2]/rho[I]+p[I]/rho[I];
    p[I]=(GAMMA-1)* (W[I][2] - rho[I]*fabs(V[I]*V[I])/2);
    
    //然后利用MUSCL分裂这些变量, 并且求出Roe平均量
    Vector WL, WR;
    double rhoL, rhoR, uL, uR, HL, HR, pL, pR;
    Vector F1, F234, F5;
    double LAMC;

    MUSCL(p, I, pR, pL);
    MUSCL(rho, I, rhoR, rhoL);
    MUSCL(u, I, uR, uL);
    MUSCL(H, I, HR, HL);
    MUSCL(W, I, WR, WL);

    //计算Roe平均量
    const double L=sqrt(rhoL)/ (sqrt(rhoL)+ sqrt(rhoR));//定义两个系数
    const double R=sqrt(rhoR)/ (sqrt(rhoL)+ sqrt(rhoR));

    const double rho_=sqrt(rhoL*rhoR);
    const double u_  =uL*L+uR*R;
    const double H_  =HL*L+HR*R;
    const double q_2 =u_*u_;
    const double c_  =sqrt((GAMMA-1)*(H_-q_2/2));
    const double V_  =u_;

    //计算各个del值
    const double delP  =pR-pL;
    const double delRho=rhoR-rhoL;
    const double delu  =uR-uL;
    const double delV  =delu;

    //求出Roe矩阵相关值
        //delu,delv,delw代表三个分量, delV代表大写V的delta
    F1[0]  = fabs(V_-c_) * (delP-rho_*c_*delV)/(2*c_*c_) * 1;
    F1[1]  = fabs(V_-c_) * (delP-rho_*c_*delV)/(2*c_*c_) * (u_-c_);
	F1[2]  = fabs(V_-c_) * (delP-rho_*c_*delV)/(2*c_*c_) * (H_-c_*V_);
  
    F234[0]= fabs(V_)*( (delRho-delP/(c_*c_))*1  + 0); 
    F234[1]= fabs(V_)*( (delRho-delP/(c_*c_))*u_ 
                  +rho_*(delu   -  delV)    );
    F234[2]= fabs(V_)*( (delRho-delP/(c_*c_))*u_
                             +rho_*(u_*delu - V_*delV) );
    
    F5[0]  = fabs(V_+c_)*(delP+rho_*c_*delV)/(2*c_*c_)*1;
    F5[1]  = fabs(V_+c_)*(delP+rho_*c_*delV)/(2*c_*c_)*(u_+c_);
    F5[0]  = fabs(V_+c_)*(delP+rho_*c_*delV)/(2*c_*c_)*(H_+c_*V_);


    double delta=0.1*soundVelocity(W[I]);
    LAMC=fabs(V_+c_);
    if(fabs(LAMC)<=delta)
    LAMC=(LAMC*LAMC+ delta*delta) /(2*delta);

    Vector FcR, FcL;
    WToF(WR, FcR);
    WToF(WL, FcL);

    //最后算出单元I的对流通量Fc
    for(int k=0; k<3; k++)
        Fc[I][k]=0.5*( FcR[k]+ FcL[k]- (F1[k]+F234[k]+F5[k]) );
}

//带限制器的MUSCL插值,前两个参数是输入,后两个输出, 系数k^为1/3
void MUSCL(Field const U, int const I, Vector UR, Vector UL)
{
    const double epsilon= dx; //限制器参数epsilon与几何尺寸相关
    
    double aR=U[I+2]-U[I+1], bR=U[I+1]-U[I];
    double aL=U[I+1]-U[I],   bL=U[I]  -U[I-1];
    double a, b;
    for(int k=0;k<3;k++)
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
void MUSCL(ScalarField const U, int const I, double & UR, double & UL)
{
    const double epsilon= dx; //限制器参数epsilon与几何尺寸相关

    double aR=U[I+2]-U[I+1], bR=U[I+1]-U[I];
    double aL=U[I+1]-U[I],   bL=U[I]  -U[I-1];

    double a=aR, b=bR;
    double deltaR=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    a=aL, b=bL;
    double deltaL=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    
    UR=U[I+1]-0.5*deltaR;
    UL=U[I]  +0.5*deltaL;

}

#endif