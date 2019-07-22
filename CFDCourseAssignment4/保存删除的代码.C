#ifdef 0
ifConverge(residual, step, Q);
static  void ifConverge(double residual, int step, Field Q);
ifConverge(residual, step, Q);
static  void ifConverge(double residual, int step, Field Q)
{
    if(residual<RESIDUAL_LIMIT)
    {
        cout<<"\nConverged!"<<endl
	    	<<"\nmaxI= "<<maxI<<"\tmaxJ= "<<maxJ<<endl
	    	<<"Final residual ="<<residual<<endl
	    	<<"Loop times = "<<step<<endl
	    	<<"Final field has been written into \"result.dat\" \n"
	    	<<"History of residual has been written into \"residual.dat\"";
        fstream fout("result.dat");
        print(Q, fout);
        cout<<"\nWall time = "<<(double)clock()/CLOCKS_PER_SEC<<" s"<<endl;
	    exit(0);
    }
    else if(step>=STOP_STEP)
	{
		cout<<"\nFail to converge! in "<<STOP_STEP<<" steps\n"
		<<"Final residual ="<<residual<<"\n\n\n";
        cout<<"\nWall time = "<<(double)clock()/CLOCKS_PER_SEC<<" s"<<endl;
	}
}

//utility functions
inline double calMa(double T, double magU);
inline double calSoundSpeed(double T);

inline double calMa(double T, double magU)
{
    return magU/(20.045*sqrt(T));
}

inline double calSoundSpeed(double T)
{
    return 20.045*sqrt(T);
}
#define safeSqrt(xx)\
{\
    if(xx<0) cout<<"\n\nerro sqrt! value is negative!\n\n\n";\
    sqrt(xx)\
}
#endif
