#pragma once

#include <omp.h>

//Very simple pseudo-random number generator based on LCG (linear congruential generator) with textbook values (modulus 2^32 used). Can be used in multi-threaded code.
//This is much faster than the standard random number generator in <random> (over 5 times faster in testing)

//To run without multi-threading set the flag in the constructor. This is useful if you want to generate random numbers reproducibly for the same seed on different computers.
//In multithreaded mode, computers with different number of corse will generate different

//linear congruential generator values obtained from "Numerical Recipes" - W.H. Press et al., Cambridge University Press, 3rd edition (2007) :
//Xn+1 = (a * Xn + c) mod m, where:
//a = 1664525
//c = 1013904223
//m = 2^32

class BorisRand {

private:

	std::vector<unsigned> prn;

public:

	BorisRand(unsigned seed, bool multithreaded_ = true)
	{
		int OmpThreads = omp_get_num_procs();
		prn.resize(OmpThreads);

		//seed all threads
		for (int idx = 0; idx < OmpThreads; idx++) {

			prn[idx] = seed * (idx + 1);
		}
	}

	//unsigned integer value out : 0 to 2^32 - 1
	unsigned randi(void)
	{
		int tn = omp_get_thread_num();

		//LCG equation used to generate next number in sequence : the modulo operation is free since unsigned is 32 bits wide
		prn[tn] = ((unsigned)1664525 * prn[tn] + (unsigned)1013904223);

		return prn[tn];
	}

	//floating point value out in interval [0, 1]
	double rand(void)
	{
		int tn = omp_get_thread_num();

		prn[tn] = ((unsigned)1664525 * prn[tn] + (unsigned)1013904223);

		return (double)prn[tn] / (unsigned)4294967295;
	}
};
