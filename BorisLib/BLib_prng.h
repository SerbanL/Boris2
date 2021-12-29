#pragma once

#include <omp.h>
#include <vector>
#include "Funcs_Math_base.h"

//Very simple pseudo-random number generator based on LCG (linear congruential generator) with textbook values (modulus 2^32 used). Can be used in multi-threaded code.
//This is much faster than the standard random number generator in <random> (over 5 times faster in testing)

//NOTE : if you want to generate the same sequence of random numbers on different computers for the same seed, don't generate them in a parallel loop.

//linear congruential generator values obtained from "Numerical Recipes" - W.H. Press et al., Cambridge University Press, 3rd edition (2007) :
//Xn+1 = (a * Xn + c) mod m, where:
//a = 1664525
//c = 1013904223
//m = 2^32

class BorisRand {

private:

	std::vector<unsigned> prn;

	//count number of random numbers generated between calls to check_periodicity : if this divides the LCG period it could be problematic so need to adjust
	std::vector<unsigned> period;
	//set to true on first call to check_periodicity
	bool calculate_period = false;

public:

	BorisRand(unsigned seed, bool multithreaded_ = true)
	{
		int OmpThreads = omp_get_num_procs();
		prn.resize(OmpThreads);
		period.resize(OmpThreads);

		//seed all threads
		for (int idx = 0; idx < OmpThreads; idx++) {

			prn[idx] = seed * (idx + 1);
			period[idx] = 0;
		}
	}

	void check_periodicity(void)
	{
		calculate_period = true;

		for (int idx = 0; idx < period.size(); idx++) {

			//if the generation period at this point matches the LCG period then notch generation by 1 point, i.e. increase period by 1.
			if (period[idx] && (unsigned)4294967295 % period[idx] == period[idx] - 1) {

				prn[idx] = ((unsigned)1664525 * prn[idx] + (unsigned)1013904223);

				//reset period : set to 1 since a point has already been generated
				period[idx] = 1;
			}
			else period[idx] = 0;
		}
	}

	//unsigned integer value out : 0 to 2^32 - 1
	unsigned randi(void)
	{
		int tn = omp_get_thread_num();

		//LCG equation used to generate next number in sequence : the modulo operation is free since unsigned is 32 bits wide
		prn[tn] = ((unsigned)1664525 * prn[tn] + (unsigned)1013904223);

		//count number of points generated on this thread since last call to check_periodicity
		if (calculate_period) period[tn]++;

		return prn[tn];
	}

	//floating point value out in interval [0, 1]
	double rand(void)
	{
		int tn = omp_get_thread_num();

		prn[tn] = ((unsigned)1664525 * prn[tn] + (unsigned)1013904223);

		//count number of points generated on this thread since last call to check_periodicity
		if (calculate_period) period[tn]++;

		return (double)prn[tn] / (unsigned)4294967295;
	}

	//Box-Muller transform to generate Gaussian distribution from uniform distribution
	double rand_gauss(double mean, double std)
	{
		thread_local double z1;
		thread_local bool generate;
		generate = !generate;

		if (!generate) return z1 * std + mean;

		double u1, u2;
		do {
			u1 = this->rand();
			u2 = this->rand();
		} while (u1 <= 1e-12);

		double z0;
		z0 = sqrt(-2.0 * log(u1)) * cos(TWO_PI * u2);
		z1 = sqrt(-2.0 * log(u1)) * sin(TWO_PI * u2);

		return z0 * std + mean;
	}
};
