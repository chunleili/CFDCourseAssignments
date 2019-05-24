#include"main.H"

static double max(double a, double b)
{
    return a>b?a:b;
}

void scalarJSTSolver(double U[][3])
{
    const double k2=0.5, k4=0.01;
    double D, Lambda, Lambda1, Ep2, Ep4, c, Y, Y1, U_half[3], R[3], F[3];
    double U1[maxSpace+1][3];
    for (int I = 1; I <= maxSpace - 2; I++)
    {
        for (int k = 0; k < 3; k++)
        {
            U_half[k] = 0.5*( U[I][k] + U[I+1][k] );

            c=sqrt( GAMMA * p(U[I])/U[I][0] );

            Lambda = ( fabs( U[I][1]/U[I][0] ) + c )* dx;
            Lambda1= ( fabs( U[I+1][1]/U[I+1][0] ) + c )* dx;
            Lambda = 0.5 *(Lambda + Lambda1);

            Y= fabs( p(U[I+1]) -2*p(U[I]) + p(U[I-1]) )
                /  ( p(U[I+1]) +2*p(U[I]) + p(U[I-1]) );

            Y1= fabs( p(U[I+2]) -2*p(U[I+1]) + p(U[I]) )
                /   ( p(U[I+2]) +2*p(U[I+1]) + p(U[I]) );

            Ep2= k2 * max(Y, Y1);
            Ep4= max(0, k4-Ep2);

            D= Lambda * (
                     Ep2 *( U[I+1][k]- U[I][k] )
                    -Ep4 *( U[I+2][k] - 3* U[I+1][k] + 3* U[I][k] - U[I-1][k] )
                        );

            UToF(U_half, F);

            R[k] = F[k] * dx -D;

            U1[I][k]=dt*R[k] + U[I][k];
        }
    }

    for (int I = 1; I <= maxSpace - 2; I++)
    {
        for (int k = 0; k < 3; k++)
        {
            U[I][k]=U1[I][k];
        }
    }

}