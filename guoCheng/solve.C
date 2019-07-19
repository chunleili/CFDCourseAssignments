//这个是一维的!待更改!
#include"main.H"
void Roe(Field Q, Index I,Index J, Tensor Fc);
void MUSCL(Field const U, Index I,Index J, Vector UR, Vector UL);
void MUSCL(ScalarField const U, Index I,Index J, double & UR, double & UL);

//三阶显式RungeKutta法
void solve(const double dt)
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
void Roe(Field Q, Index I, Index J, Tensor Fc, const XY dS, const double dV)
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

    MUSCL(rho, I, rhoR, rhoL);
    MUSCL(p, I, pR, pL);    //注意不要越界~
    MUSCL(u, I, uR, uL);
    MUSCL(H, I, HR, HL);
    MUSCL(Q, I, QR, QL);

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
    V_.x=u_; V_.y=v_;

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

    double delta=0.1*c_;
    LAMC=fabs(V_+c_);
    if(fabs(LAMC)<=delta)
         LAMC=(LAMC*LAMC+ delta*delta) /(2*delta);
    
    is0(delta);

    Vector FcR, FcL;
    WToF(WR, FcR);
    WToF(WL, FcL);

    //最后算出单元I的对流通量Fc
    for(unsigned k=0; k<3; k++)
        Fc[I][k]=0.5*( FcR[k]+ FcL[k]- (F1[k]+F234[k]+F5[k]) );
}

//带限制器的MUSCL插值,前两个参数是输入,后两个输出, 系数k^为1/3
void MUSCL(Field const U, Index I, Vector UR, Vector UL)
{
    if(I+2>maxSpace || I-1<0) {cout<<"out of bound!!!"<<endl;  }

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
void MUSCL(ScalarField const U, Index I, double & UR, double & UL)
{
    if(I+2>maxSpace || I-1<0) {cout<<"out of bound!!!"<<endl; }
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


