    FILE *fp;
    fp=fopen("something.dat", "w");
    fprintf(fp,"x    y    pressure\n");
    for(unsigned i=0; i<=maxI; i++)
        for(unsigned j=0; j<=maxJ; j++)
        {
            fprintf(fp1, "%.2f %.2f %.5e\n", mesh[i][j].x, mesh[i][j].y, p[I][J] );
        } 
    fclose(fp1);