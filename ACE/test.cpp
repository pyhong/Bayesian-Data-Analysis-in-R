#include <iostream>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <fstream>
using namespace std;

int main()
{
	gsl_rng * r = gsl_rng_alloc(gsl_rng_default);
	gsl_rng_set(r, time(NULL));
	char filename[] = "gamma_number.txt";
	ofstream fout(filename);
	for(int i = 0; i < 2000; ++i)
	{
		fout<<gsl_ran_gamma(r, 500, 0.005)<<" ";
	}
	fout<<endl;
	return 0;
}
