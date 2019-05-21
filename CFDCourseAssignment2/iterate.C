/*---------------------------------------------------------------------------*\
Description
	This File is to iterate the loop once. It contains the core algorithm.

\*---------------------------------------------------------------------------*/

#include"main.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
double iterate()
{
	double alpha, beta, gamma, bw, be, bs, bn, bp, cpx, cpy;
	double r, rMax=0.0;
	int i,j;
	for(i=1;i<Mx;i++)
	{
		for(j=1;j<My;j++)
		{
		
			//					Calculate the coefficients					 //
			// 	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
			alpha=	sq
					(
					 	( node[i][j+1].x - node[i][j-1].x ) / ( 2*dEta ) 
					)
				  +	sq
				  	(
				   		( node[i][j+1].y - node[i][j-1].y ) / ( 2*dEta )
				   	);
			
			beta=	( node[i+1][j].x - node[i-1][j].x ) / ( 2*dKsi ) 
				  *	( node[i][j+1].x - node[i][j-1].x ) / ( 2*dEta )
				  +
				  	( node[i+1][j].y - node[i-1][j].y ) / ( 2*dKsi ) 
				  *	( node[i][j+1].y - node[i][j-1].y ) / ( 2*dEta );
					
			gamma=	sq
					(
						( node[i+1][j].x - node[i-1][j].x ) / ( 2*dKsi )
					)
				  +	sq
				  	(
				  	 	( node[i+1][j].y - node[i-1][j].y ) / ( 2*dKsi )
				  	);
			
			bw=	be= alpha/(dKsi*dKsi);
			
			bs=	bn=	gamma/(dEta*dEta);
			
			bp=bw+be+bs+bn;
			
			cpx= - beta
				*(
					 node[i+1][j+1].x
					-node[i+1][j-1].x
					-node[i-1][j+1].x
					+node[i-1][j-1].x
				 )
				/	(2*dKsi*dEta);
				 
			cpy= - beta
				*(
					 node[i+1][j+1].y
					-node[i+1][j-1].y
					-node[i-1][j+1].y
					+node[i-1][j-1].y
				 )
				/	(2*dKsi*dEta);

			//				Calculate the new nodes value					 //
			// 	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //	 
			double x_new=(
							  bw*node[i-1][j].x
							+ be*node[i+1][j].x
							+ bs*node[i][j-1].x
							+ bn*node[i][j+1].x
							+ cpx
						) / bp	;
						
			double y_new=(
							  bw*node[i-1][j].y
							+ be*node[i+1][j].y
							+ bs*node[i][j-1].y
							+ bn*node[i][j+1].y
							+ cpy
						) / bp	;

			//				calculate the residual							 //
			// 	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //						
			r=max
				(
					fabs( x_new - node[i][j].x ) ,
					fabs( y_new - node[i][j].y )
				);
		    if(r>rMax)
		    	rMax=r;

			//				update the nodes, move to next					 //
			// 	* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  //
			node[i][j].x=x_new;
			node[i][j].y=y_new;

		}
	}
	return rMax;
}									
