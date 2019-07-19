#include"main.H"
int main()
{  
	int caseNo=1;
	cout<<"Case 1 for inlet 1.8Ma; 2 for 1.5Ma\n Case: "<<caseNo<<endl;
	Mesh mesh;
    FlowField field1(1);
	
    double residual=1.0;
    for (int step=0;  step<=STOP_STEP; step++)
    {    
        field1.solve(mesh);                  //求解
    }

    return 0;
}