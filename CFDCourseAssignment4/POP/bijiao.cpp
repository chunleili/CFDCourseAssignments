//yuan
for (i = 1; i <= maxI; i++)
    for (j = 1; j <= maxJ; j++)
    {
        c = sqrt(GAMMA * p[i][j] / rho[i][j]);
        splitI();
        calARoe();
        for (k = 0; k <= 3; k++)
            Fc2[i][j][k] = S2[i][j] * (FR[k] + FL[k] - AARoe[k]) / 2;
        splitJ();
        calARoe();
        for (k = 0; k <= 3; k++)
            Fc3[i][j][k] = S3[i][j] * (FR[k] + FL[k] - AARoe[k]) / 2;
    }

for (i = 1; i <= maxI - 1; i++)
    for (j = 1; j <= maxJ - 1; j++)
        for (k = 0; k <= 3; k++)
        {
            Fc4[i][j][k] = Fc2[i - 1][j][k];
            Fc1[i][j][k] = Fc3[i][j - 1][k];
            Residual[i][j][k] = -Fc1[i][j][k] + Fc2[i][j][k] - Fc3[i][j][k] + Fc4[i][j][k];
        }
    