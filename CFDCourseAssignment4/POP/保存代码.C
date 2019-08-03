#ifdef 0
    FILE *fp;
    fp=fopen("something.dat", "w");
    fprintf(fp,"x    y    pressure\n");
    for(unsigned i=0; i<=maxI; i++)
        for(unsigned j=0; j<=maxJ; j++)
        {
            fprintf(fp1, "%.2f %.2f %.5e\n", mesh[i][j].x, mesh[i][j].y, p[I][J] );
        } 
    fclose(fp1);
//上壁面
        rho1=Q[i][j][0];        

            Vt=u[i][j]*ny+v[i][j]*nx;
        newU=Vt*ny;//速度u,v要重新分配
        newV=Vt*nx;

        oldRhoV2=(SQ(Q[i][j][1])+SQ(Q[i][j][2]))/rho1;
        newRhoV2=rho1*(SQ(newU)+SQ(newV)); 


        Q[i][maxJ][0]=Q[i][maxJ-1][0]=rho1;
        Q[i][maxJ][1]=Q[i][maxJ-1][1]=rho1*newU;
        Q[i][maxJ][2]=Q[i][maxJ-1][2]=rho1*newV;
        Q[i][maxJ][3]=Q[i][maxJ-1][3]=Q[i][j][3]-0.5*(newRhoV2-oldRhoV2);

//下壁面
        rho1=Q[i][j][0];

                Vt=u[i][j]*ny+v[i][j]*nx;
        newU=Vt*ny;
        newV=Vt*nx;

        oldRhoV2=(SQ(Q[i][j][1])+SQ(Q[i][j][2]))/rho1;
        newRhoV2=rho1*(SQ(newU)+SQ(newV));

        Q[i][0][0]=rho1;
        Q[i][0][1]=rho1*newU;
        Q[i][0][2]=rho1*newV;
        Q[i][0][3]=Q[i][j][3]-0.5*(newRhoV2-oldRhoV2);








void wallBC()
{
    //上下是壁面, 壁面法向速度为0,壁面切向不变 壁面无穿透边界
    //密度,静压不变, E随之变化
    double nx,ny,p2,p3,pw;
    unsigned j;
    for(unsigned i=cellBegin; i<=cellJEnd; i++)
    {
        //上壁面
        j=cellJEnd;

        nx=N3[i][j].x;
        ny=N3[i][j].y;

        p2=p[i][cellJEnd];
        p3=p[i][cellJEnd-1];
        pw=0.5*(3*p2-p3);//壁面的压力用两点外推

        Fc3[i][j][0]=0;       
        Fc3[i][j][1]=pw*nx;
        Fc3[i][j][2]=pw*ny;
        Fc3[i][j][3]=0;
        
        //虚网格的值靠外推
        for(unsigned k=0; k<=3; k++)
        {
            Q[i][cellJEnd+1][k]=2*Q[i][j][k]-  Q[i][j-1][k];
            Q[i][cellJEnd+2][k]=3*Q[i][j][k]-2*Q[i][j-1][k];
        }
        aeroConvert(i,cellJEnd+2);
        aeroConvert(i,cellJEnd+1);
        Vcv2[i][cellJEnd+2]=Vcv2[i][cellJEnd+1] = Vcv2[i][cellJEnd];
        Vcv3[i][cellJEnd+2]=Vcv3[i][cellJEnd+1] = Vcv3[i][cellJEnd];

        //下壁面,只有一层虚网格
        j=cellBegin;

        nx=N1[i][j].x;
        ny=N1[i][j].y;

        p2=p[i][j];
        p3=p[i][j+1];
        pw=0.5*(3*p2-p3);//壁面的压力用两点外推

        Fc1[i][j][0]=0;       
        Fc1[i][j][1]=pw*nx;
        Fc1[i][j][2]=pw*ny;
        Fc1[i][j][3]=0;
        
        //虚网格的值靠外推
        for(unsigned k=0; k<=3; k++)
        {
            Q[i][cellBegin-1][k]=2*Q[i][j][k]-  Q[i][j+1][k];
        }
        aeroConvert(i,cellBegin-1);
        Vcv2[i][cellBegin-1] = Vcv2[i][cellBegin];
        Vcv3[i][cellBegin-1] = Vcv3[i][cellBegin];

        if(rho[i][cellBegin-1]<0) 
        {
            cout<<"here! rho<0 i="<<i<<" j="<<cellBegin-1<<" rho="<<rho[i][cellBegin-1]<<endl;
            printf("\nQ[i][0]=%f  %f  %f  %f\n", Q[i][0][0], Q[i][0][1],Q[i][0][2],Q[i][0][3]);
            printf("\nQ[i][1]=%f  %f  %f  %f\n", Q[i][1][0], Q[i][1][1],Q[i][1][2],Q[i][1][3]);
            printf("\nQ[i][2]=%f  %f  %f  %f\n", Q[i][2][0], Q[i][2][1],Q[i][2][2],Q[i][2][3]);
        }
    }     
}

//出口超音速,出口全部外推
void BC1()
{
    //入口维持1.8Ma, 出口用0梯度外推出来;
    //const double MaInf=1.8, UInf= ; 
    for(unsigned j=cellBegin; j<=cellJEnd; j++)
    {
        //入口虚网格不用变(初始化时候已经给定)

        //出口虚网格
        for (unsigned k = 0; k < 4; k++)
        {
            Q[cellIEnd+2][j][k] = Q[cellIEnd+1][j][k] = Q[cellIEnd][j][k];
        }
        aeroConvert(cellIEnd+2,j);
        aeroConvert(cellIEnd+1,j);
        Vcv2[cellIEnd+2][j] = Vcv2[cellIEnd+1][j]=Vcv2[cellIEnd][j];
        Vcv3[cellIEnd+2][j] = Vcv3[cellIEnd+1][j]=Vcv3[cellIEnd][j];
    }
    
    wallBC();
}

//出口亚音速,给压力,外推u,v,rho
void BC2()
{        
    double pout=200000.0, rho1, u1, v1;
    for(unsigned j=cellBegin; j<=cellJEnd; j++)
    {
        //入口全部给定
        Q[0][j][0]=1.176829;
        Q[0][j][1]=612.873; //1.176829*520.783
        Q[0][j][2]=0;
        Q[0][j][3]=412899.3; //  101325/0.4+0.5*1.176829*520.783^2;
        //出口给三推一
        for(unsigned k=0; k<3; k++)//注意k从0到2, 前三个量 rho, rho*u, rho*v外推
        {
            Q[cellIEnd+2][j][k] = Q[cellIEnd+1][j][k] = Q[cellIEnd][j][k];
        }
        //E靠p给定
        rho1=Q[maxI-2][j][0];
        u1=Q[maxI-2][j][1]/rho1;
        v1=Q[maxI-2][j][2]/rho1;
        Q[maxI][j][3]=Q[maxI-1][j][3]=pout/0.4+0.5*rho1*(u1*u1+v1*v1);
    }
    
    wallBC();
}



void init2()
{
    //初始全部给1.5Ma对应的速度, 压力给大气压101325, 静温300, 
    //速度u=1.5*347.2=520.783, v=0, 密度rho=p/RT=1.176829
    forAll(
        Q[i][j][0]=1.176829;
        Q[i][j][1]=612.873; //1.176829*520.783
        Q[i][j][2]=0;
        Q[i][j][3]=412899.3; //  101325/0.4+0.5*1.176829*520.783^2;
    );
}




//Q与Fc之间的转化
/*
void toFlux2(Index i, Index j)
{
    const double H1=(Q[i][j][3]+p[i][j])/rho[i][j];
    F[i][j][0] = Q[i][j][0];
    F[i][j][1] = Q[i][j][1]*Vcv2[i][j]+N2[i][j].x*p[i][j];
    F[i][j][2] = Q[i][j][2]*Vcv2[i][j]+N2[i][j].y*p[i][j];
    F[i][j][3] = Q[i][j][0]*H1*Vcv2[i][j];
}
void toFlux3(Index i, Index j)
{
    const double H1=(Q[i][j][3]+p[i][j])/rho[i][j];
    F[i][j][0] = Q[i][j][0];
    F[i][j][1] = Q[i][j][1]*Vcv3[i][j]+N3[i][j].x*p[i][j];
    F[i][j][2] = Q[i][j][2]*Vcv3[i][j]+N3[i][j].y*p[i][j];
    F[i][j][3] = Q[i][j][0]*H1*Vcv3[i][j];
}
*/

//if (caseNo==1) BC1();
//else           BC2();

//cout<<"\n\n******************************************************\n";
//cout<<"step= "<<step<<" a= "<<a<<"  I= "<<I<<" J= "<<J<<endl;
//cout<<"step= "<<step<<"  I= "<<I<<" J= "<<J<<endl;


        //左边
        for(unsigned j=cellBegin; j<=cellJEnd; j++)
        {
            N1[cellBegin-1][j].x=N1[cellBegin][j].x;
            N1[cellBegin-1][j].y=N1[cellBegin][j].y;
        }

        //右边
        for(unsigned j=cellBegin; j<=cellJEnd; j++)
        {
            N1[cellIEnd+2][j].x=N1[cellIEnd+1][j].x=N1[cellIEnd][j].x;
            N1[cellIEnd+2][j].y=N1[cellIEnd+1][j].y=N1[cellIEnd][j].y;
        }

        //上边
        for(unsigned i=cellBegin; i<=cellIEnd ;i++)
        {
            N1[i][cellJEnd+2].x=N1[i][cellJEnd+1].x=N1[i][cellJEnd].x;
            N1[i][cellJEnd+2].y=N1[i][cellJEnd+1].y=N1[i][cellJEnd].y;
        }
#endif