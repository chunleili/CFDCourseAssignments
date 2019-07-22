#include"main.H"
/***************************utility  **************************/
//Q与Fc之间的转化
void toFlux(Vector Q, Vector F)
{
    double   u=Q[1]/Q[0];
    double   p0=convert(Q).p;
    F[0] = Q[1];
    F[1] = Q[0] * u * u + p0;
    F[2] = (Q[2] + p0) * u;
}

//气动参数转换
AERO AeroConvert(Vector vec)
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