//General purpose math functions

#pragma once

#include <math.h>

#define TWO_PI 6.2831853071795864769252867666
#define PI 3.1415926535897932384626433833
#define MU0 1.256637061435916e-6					//Units kg m / s^2 A^2 = N / A^2
#define MUB 9.27400968e-24							//Bohr magneton, units Am^2 = J / T
#define MUB_MU0	1.165406427200340E-29				//mu0 * muB
#define ECHARGE 1.60217662e-19						//electron charge (C)
#define MUB_E 5.788381608E-05						//Bohr magneton / electron charge
#define HBAR_E 6.5821194172E-16						//hbar / e
#define BOLTZMANN 1.3806488e-23						//Boltzmann's constant, kB (m^2 kg / s^2 K)
#define GAMMA 2.212761569e5							//modulus of MU0 * gamma_e (m/As), where gamma_e is the electron gyromagnetic ratio
#define GAMMA_E	-1.76086E+11						//electron gyromagnetic ratio (note the -ve sign!). Units 1/Ts
#define GAMMA_E_MOD	1.76086E+11						//modulus of electron gyromagnetic ratio
#define GMUB_2E 5.795094E-05						//units m^2/s, this is -hbar * gamma_e / 2e = g mu_b / 2e, where gamma_e = -g e / 2m_e = -g mu_b / hbar

///////////////////////////////////////////////////////////////////////////////

template <typename RType>
RType maximum(RType param1, RType param2)
{
	return (param1 > param2 ? param1 : param2);
}

template <typename RType, typename ... PType>
RType maximum(RType param, PType ... params)
{
	RType value = maximum(params...);
	return (param > value ? param : value);
}

template <typename RType>
RType minimum(RType param1, RType param2)
{
	return (param1 < param2 ? param1 : param2);
}

template <typename RType, typename ... PType>
RType minimum(RType param, PType ... params)
{
	RType value = minimum(params...);
	return (param < value ? param : value);
}

///////////////////////////////////////////////////////////////////////////////

//value = return_value * 10^decexp, where |decexp| is a multiple of 3
inline double decexp_eng(const double value, int *decexp, const double minimum = 1e-18, const double maximum = 1e+18)
{
	double avalue = fabs(value);

	if (avalue < minimum || avalue >= maximum) { *decexp = 0; return value; }

	*decexp = 0;
	double divisor = 1;

	while (true) {

		//do not use IsGE or similar here - the epsilon is too large
		if (avalue / divisor >= 1e3) {

			*decexp += 3;
			divisor *= 1e3;
			continue;
		}

		if (avalue / divisor < 1) {

			*decexp -= 3;
			divisor *= 1e-3;
			continue;
		}

		break;
	}

	return (value / divisor);
}

//value = return_value * 10^decexp, where return_value is a real in the interval [1, 10)
inline double decexp(const double value, int *decexp, const double minimum = 1e-18, const double maximum = 1e+18)
{
	double avalue = fabs(value);

	if (avalue < minimum || avalue >= maximum) { *decexp = 0; return value; }

	*decexp = 0;
	double divisor = 1;

	while (true) {

		//do not use IsGE or similar here - the epsilon is too large
		if (avalue / divisor >= 10) {

			*decexp += 1;
			divisor *= 10;
			continue;
		}

		if (avalue / divisor < 1) {

			*decexp -= 1;
			divisor *= 0.1;
			continue;
		}

		break;
	}

	return (value / divisor);
}