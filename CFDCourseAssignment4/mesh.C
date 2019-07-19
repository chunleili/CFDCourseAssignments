#include"main.H"
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

//计算单元体左下侧面积
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

    S.x=safeSqrt( (x1-x2)*(x1-x2)+(y2-y1)*(y2-y1) );//下侧面积S1
    S.y=safeSqrt( (x4-x1)*(x4-x1)+(y1-y4)*(y1-y4) );//左侧面积S4
    
    return S;
}
