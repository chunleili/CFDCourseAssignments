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
    #endif