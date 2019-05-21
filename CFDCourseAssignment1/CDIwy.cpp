#include <iostream>
using namespace std;
int main()
{

	int t,i,j,flag;
	double x[100],x_n1[100];
	double a[100],b[100],c[100];
   	 double l[100],m[100],u[100],y[100]; 
    //取deta t=0.001,deta x=0.001,a=0.8 
	for(i=0;i<100;i++)
	{
		x[i]=0;
	}
	for(i=10;i<20;i++)
	{
		x[i]=100;
	}
	for(i=20;i<101;i++)
	{
		x[i]=0;
	}             //设置初值 



     	for(t=0;t<=100;t++)
     	{
     		for(j=1;j<100;j++)
     		{
     			a[j]=1;
     			b[j]=-0.35;
     			c[j]=0.35;
		 }
     		b[1]=0;c[99]=0;  //系数矩阵 
     		l[1]=a[1];
     		m[1]=0;

     		for(i=2;i<100;i++)
     		{
     			m[i]=b[i];
     			u[i-1]=c[i-1]/l[i-1];
     			l[i]=a[i]-b[i]*u[i-1];
       		 }
       		 u[99]=0;
        	y[1]=x[i]/l[1];

        	for(i=2;i<100;i++)
        	{
        		y[i]=(x[i]-m[i]*y[i-1])/l[i];
        	}

     		x_n1[99]=y[99];

     		for(i=98;i>0;i--)
     		{
			 x_n1[i]=y[i]-u[i]*x_n1[i+1];
        	}

     		x_n1[0]=0;x_n1[100]=0;
		
		for(i=0;i<101;i++)
		{
			x[i]=x_n1[i];
		}
	} 
 
//print
	for(i=0; i<=100; i++)
	{
		cout<<x[i]<<"\n";
	} 
return 0;
}

