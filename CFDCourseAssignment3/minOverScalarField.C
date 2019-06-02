#include"main.H"
double minOverScalarField(double scalarField[])
{
	double minVal;
	for (int i = 0; i <=maxSpace ; i++)
	{
		if (minVal<scalarField[i])
		{
			minVal=scalarField[i];
		}
	}
	return minVal;
}