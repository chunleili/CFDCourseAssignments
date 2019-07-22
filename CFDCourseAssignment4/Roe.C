#include"main.H"
/**********************Roe********************************/
Roe::Roe(Index II, Index JJ):I(II), J(JJ)
{}

//Roe格式计算对流通量
void setFlux1()
{  
    //先定义流场各变量
    ScalarField rho, u, v,  p, H;

    rho[I][J]=Q[I][J][0];
    u[I][J]=Q[I][J][1]/rho[I][J];
    v[I][J]=Q[I][J][1]/rho[I][J];
    p[I][J]=(GAMMA-1)* (Q[I][J][2] - rho[I][J]* (SQ(u[I][J])+SQ(v[I][J])) *0.5);
    H[I][J]=Q[I][J][3]/rho[I][J]+p[I][J]/rho[I][J];

    //然后利用MUSCL分裂这些变量, 并且求出Roe平均量
    Vector QL, QR;

    MUSCL1(rho, I, J, rhoR, rhoL);
    MUSCL1(p,   I, J, pR,   pL  );    //注意不要越界~
    MUSCL1(u,   I, J, uR,   uL  );
    MUSCL1(v,   I, J, vR,   vL  );
    MUSCL1(H,   I, J, HR,   HL  );
    MUSCL1(Q,   I, J, QR,   QL  );

    //计算Roe平均量
    const double denoLR=safeSqrt(rhoL)+ safeSqrt(rhoR);
    const double L=safeSqrt(rhoL)/ denoLR;//定义两个系数
    const double R=safeSqrt(rhoR)/ denoLR;
    
    const double rho_=safeSqrt(rhoL*rhoR);
    const double u_  =uL*L+uR*R;
    const double v_  =vL*L+vR*R;
    const double H_  =HL*L+HR*R;
    const double q_2 =u_*u_+v_*v_;
    const double c_  =safeSqrt((GAMMA-1)*(H_-q_2/2));

    XY N1=mesh.
    double Vcv_=N1.x * u_ +N1.y* v_;
    //Vcv代表contravirant velocity 逆变速度

    //计算各个del值
    const double delP  =pR-pL;
    const double delRho=rhoR-rhoL;
    const double delu  =uR-uL;

    lambda1=Herten(fabs(Vcv_-c_));
    lambda2=Herten(fabs(Vcv_   ));
    lambda3=Herten(fabs(Vcv_+c_));   

    double delVcv  =Vcv1_R-vcv1_L;
    //求出Roe矩阵相关值
        //delu,delv,delw代表三个分量, delVcv代表大写Vcv的delta
        //lambda123分别是是Vcv-c, Vcv, Vcv+c
    const double delP  =pR-pL;
    const double delRho=rhoR-rhoL;
    const double delu  =uR-uL;
    const double delVcv=VcvR-VcvL;
    XY N1=mesh.getN1(I,J);
    F1[0]  = lambda1 * (delP-rho_*c_*delVcv)/(2*c_*c_) * 1;
    F1[1]  = lambda1 * (delP-rho_*c_*delVcv)/(2*c_*c_) * (u_-c_*N1.x);
    F1[2]  = lambda1 * (delP-rho_*c_*delVcv)/(2*c_*c_) * (u_-c_*N1.y);
	F1[3]  = lambda1 * (delP-rho_*c_*delVcv)/(2*c_*c_) * (H_-c_*Vcv_);
  
    F234[0]= lambda2 * ( (delRho-delP/(c_*c_))*1  + 0); 
    F234[1]= lambda2 * ( (delRho-delP/(c_*c_))*u_ 
                   + rho_ *(delu   -  delVcv)    );
    F234[2]= lambda2 * ( (delRho-delP/(c_*c_))*u_ 
                   + rho_ *(delu   -  delVcv)    );                
    F234[3]= lambda2 * ( (delRho-delP/(c_*c_))*u_
                              +rho_ *(u_*delu - Vcv_*delVcv) );
     
    F5[0]  = lambda3 * (delP+rho_*c_*delVcv)/(2*c_*c_) *1;
    F5[1]  = lambda3 * (delP+rho_*c_*delVcv)/(2*c_*c_) *(u_+c_);
    F5[2]  = lambda3 * (delP+rho_*c_*delVcv)/(2*c_*c_) *(u_+c_);    
    F5[3]  = lambda3 * (delP+rho_*c_*delVcv)/(2*c_*c_) *(H_+c_*Vcv_);

    Vector FcR, FcL;
    toFlux(QR, FcR);
    toFlux(QL, FcL);


    //最后算出单元IJ的对流通量Fc
    for(unsigned k=0; k<3; k++)
    {
        Fc[I][J][k].x=0.5*( FcR[k] - FcL[k] - (F1[k]+F234[k]+F5[k]) );
    }
}

//带限制器的MUSCL插值,前两个参数是输入,后两个输出, 系数k^为1/3
void Roe::MUSCL1(Field const U, Vector UR, Vector UL)
{
    if(I+2>maxI || I-1<0) {cout<<"\n\nout of bound!!!\n\n";}
    const double dV=mesh.getVolume(I,J);
    const double epsilon= pow(dV,1.0/3); //限制器参数epsilon与几何尺寸相关
    
    double aR=U[I+2][J]- U[I+1][J],   bR=U[I+1][J] - U[I  ][J];
    double aL=U[I+1][J]- U[I  ][J],   bL=U[I  ][J] - U[I-1][J];
    double a, b;
    for(unsigned k=0;k<3;k++)
    {
        a=aR; b=bR;
        double deltaR=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
        a=aL; b=bL;
        double deltaL=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);

        UR[k]=U[I+1][J][k]-0.5*deltaR;
        UL[k]=U[I  ][J][k]+0.5*deltaL;
    }
}
//重载用于标量的带限制器的MUSCL插值函数 ,前两个参数是输入,后两个输出, 系数k^为1/3
void Roe::MUSCL1(ScalarField const U, double & UR, double & UL)
{
    if(I+2>maxJ || I-1<0) {cout<<"\n\nout of bound!!!\n\n";}
    const double dV=mesh.getVolume(I,J);
    const double epsilon= pow(dV,1.0/3); //限制器参数epsilon与几何尺寸相关

    double aR=U[I+2][J]- U[I+1][J],   bR=U[I+1][J] - U[I  ][J];
    double aL=U[I+1][J]- U[I  ][J],   bL=U[I  ][J] - U[I-1][J];

    double a=aR, b=bR;
    double deltaR=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    a=aL, b=bL;
    double deltaL=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    
    UR=U[I][J+1]-0.5*deltaR;
    UL=U[I][J]  +0.5*deltaL;
}





/***********************flux4****************************/
/***********************flux4****************************/
/***********************flux4****************************/
/***********************flux4****************************/

//重载用于标量的带限制器的MUSCL插值函数 ,前两个参数是输入,后两个输出, 系数k^为1/3
void Roe::MUSCL4(ScalarField const U double & UR, double & UL)
{
    if(J+2>maxJ || J-1<0) {cout<<"\n\nout of bound!!!\n\n";}
    const double epsilon= dx; //限制器参数epsilon与几何尺寸相关

    double aR=U[I][J+2]-U[I][J+1], bR=U[I][J+1]-U[I][J];
    double aL=U[I][J+1]-U[I][J],   bL=U[I][J]  -U[I][J-1];

    double a=aR, b=bR;
    double deltaR=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    a=aL, b=bL;
    double deltaL=( (2*a*a+epsilon)*b+(b*b+2*epsilon)*a ) /( 2*a*a +2*b*b -a*b +3*epsilon);
    
    UR=U[I][J+1]-0.5*deltaR;
    UL=U[I][J]  +0.5*deltaL;
}

double Herten(double lambda)
{
    double delta=0.1*c_; //熵修正, Harten's entropy correction
    if(fabs(lambda)<=delta)
         lambda=(lambda*lambda+ delta*delta) /(2*delta);
    return lambda;
}

XY Roe::getFlux() const
{
    XY flux;
    flux.x=Fc1;
    flux.y=Fc4;
    return flux;
}