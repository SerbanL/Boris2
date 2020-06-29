#pragma once

#include "Funcs_Math_base.h"
#include "Types_VAL.h"
#include "Funcs_Algorithms.h"
#include "TEquation_FSPEC.h"

//Special functions which require initialization to use.
//Initialization may be computationally expensive so need to be careful to re-use code as much as possible.

class Funcs_Special {

	//smallest value to calculate in values : store it at values[0]
	double start = 0.0;

	//number of points per unit to store in values
	int resolution = 1;

	std::vector<double> values;

public:

	///////////////////////////////////////////////////////////////

	Funcs_Special(EqComp::FUNC_ type, double start_ = 0.0, int resolution_ = 1) :
		start(start_), resolution(resolution_)
	{
		switch (type) {

		case EqComp::FUNC_CURIEWEISS:
			Initialize_CurieWeiss(0.0, 1.0);
			break;

		case EqComp::FUNC_CURIEWEISS1:
			Initialize_CurieWeiss1(DBL2(0.5), DBL2(0.5), DBL2(0.0), 1.0);
			break;

		case EqComp::FUNC_CURIEWEISS2:
			Initialize_CurieWeiss1(DBL2(0.5), DBL2(0.5), DBL2(0.0), 1.0);
			break;

		case EqComp::FUNC_LONGRELSUS:
			Initialize_CurieWeiss(0.0, 1.0);
			Initialize_LongitudinalRelSusceptibility(values, 1.0, 1.0);
			break;

		case EqComp::FUNC_LONGRELSUS1:
			Initialize_CurieWeiss1(DBL2(0.5), DBL2(0.5), DBL2(0.0), 1.0);
			Initialize_LongitudinalRelSusceptibility1(values, values, DBL2(0.5), DBL2(0.5), DBL2(1.0), 1.0);
			break;

		case EqComp::FUNC_LONGRELSUS2:
			Initialize_CurieWeiss2(DBL2(0.5), DBL2(0.5), DBL2(0.0), 1.0);
			Initialize_LongitudinalRelSusceptibility2(values, values, DBL2(0.5), DBL2(0.5), DBL2(1.0), 1.0);
			break;

		case EqComp::FUNC_ALPHA1:
			Initialize_Alpha1(0.5, 0.5);
			break;

		case EqComp::FUNC_ALPHA2:
			Initialize_Alpha2(0.5, 0.5);
			break;

		default:
			values = { 0.0 };
			break;
		}
	}

	///////////////////////////////////////////////////////////////
	// CURIE-WEISS

	//me = B(me * 3Tc / T + mu * mu0*Ha / kBT), where B(x) = coth(x) - 1 / x
	//Set Tc = 1 (normalized Curie-Weiss) and calculate between 0 and 2 (2 times Tc).
	//Let d = mu * mu0*Ha / kB
	//Thus this function calculates me(T):
	//me = B((3 * me + d) / T), for T in 0 to 2
	//Note, because we've normalised to Tc, if d is not zero (there's an applied field) we need to divide it by the actual required Tc to get the correct scaling.
	void Initialize_CurieWeiss(double d, double Tc)
	{
		if (values.size() != resolution * 2 + 1) values.resize(resolution * 2 + 1);

		double me = 1.0;
		values[0] = me;

		for (int t_idx = 1; t_idx <= resolution * 2; t_idx++) {

			double t_value = (double)t_idx / resolution;

			double c = 3 / t_value;

			//This is F(x) = B(cx + d) - x
			std::function<double(double)> F = [&](double x) -> double {

				double y = c * x + (d / Tc) / t_value;
				return ((1 + exp(-2 * y)) / (1 - exp(-2 * y))) - 1 / y - x;
			};

			//This is F'(x) = B'(cx + d) * c - 1
			std::function<double(double)> Fdiff = [&](double x) -> double {

				double y = c * x + (d / Tc) / t_value;
				double coth = (1 + exp(-2 * y)) / (1 - exp(-2 * y));
				return ((1 - coth * coth) + 1 / (y*y)) * c - 1;
			};

			//next starting value for the root finder is me - this makes the NewtonRaphson algorithm extremely efficient (typically a single iteration is required to keep within 1e-6 accuracy!)

			double me_previous = me;

			me = Root_NewtonRaphson(F, Fdiff, me, 1e-9, 100);

			//me needs to be decreasing and positive : when me is very small, numerical noise becomes a big factor (above Tc for small fields).
			if (me > me_previous || me < 0) me = me_previous;

			values[t_idx] = me;
		}
	}

	///////////////////////////////////////////////////////////////
	// CURIE-WEISS FOR 2-SUBLATTICE MODEL

	//me1 = B[ (me1*tau1 + me2*tau12) * 3rTc/T + mu1*mu0*Ha/kBT ], where B(x) = coth(x) - 1 / x
	//me2 = B[ (me2*tau2 + me1*tau21) * 3rTc/T + mu2*mu0*Ha/kBT ], where B(x) = coth(x) - 1 / x

	//rTc is the renormalized Tc value as rTc = Tc / ((tau1 + tau2 + sqrt((tau1 - tau2)^2 + 4 * tau12 * tau21)) / 2);

	//Set Tc = 1 (normalized Curie-Weiss) and calculate between 0 and 2 (2 times Tc).
	//Let d1 = mu1 * mu0*Ha / kB, d2 = mu2 * mu0*Ha / kB

	//Thus this function calculates me(T):
	//me1 = B[ (me1*tau1 + me2*tau12) * 3r/T + d1/T ] = B[ (me1*tau1 + me2*tau12) * c + d1/T ], for T in 0 to 2
	//me2 = B[ (me2*tau2 + me1*tau21) * 3r/T + d2/T ] = B[ (me2*tau2 + me1*tau21) * c + d2/T ], for T in 0 to 2

	//Note, because we've normalised to Tc, if d is not zero (there's an applied field) we need to divide it by the actual required Tc (renormalized) to get the correct scaling.
	
	//this function calculates me1 only and stores it here so we can access it separately
	void Initialize_CurieWeiss1(DBL2 tau_ii, DBL2 tau_ij, DBL2 d, double Tc)
	{
		if (values.size() != resolution * 2 + 1) values.resize(resolution * 2 + 1);

		DBL2 me = DBL2(1.0);
		values[0] = me.i;

		double Tc_renorm = (tau_ii.i + tau_ii.j + sqrt(pow(tau_ii.i - tau_ii.j, 2) + 4 * tau_ij.i * tau_ij.j)) / 2;

		for (int t_idx = 1; t_idx <= resolution * 2; t_idx++) {

			double T = (double)t_idx / resolution;

			double c = 3 / (T * Tc_renorm);
			double d1 = d.i * Tc_renorm / (T * Tc);
			double d2 = d.j * Tc_renorm / (T * Tc);

			//This is F1(x1, x2) = B((x1 * tau1 + x2 * tau12) * c + mu1 * d) - x1
			std::function<double(DBL2)> F1 = [&](DBL2 x) -> double {

				double y = c * (tau_ii.i * x.i + x.j * tau_ij.i) + d1;
				return ((1 + exp(-2 * y)) / (1 - exp(-2 * y))) - 1 / y - x.i;
			};

			//This is F2(x1, x2) = B((x2 * tau2 + x1 * tau21) * c + mu2 * d) - x2
			std::function<double(DBL2)> F2 = [&](DBL2 x) -> double {

				double y = c * (tau_ii.j * x.j + x.i * tau_ij.j) + d2;
				return ((1 + exp(-2 * y)) / (1 - exp(-2 * y))) - 1 / y - x.j;
			};

			//This is F11'(x1, x2) = B'((x1 * tau1 + x2 * tau12) * c + mu1 * d) * c * tau1 - 1
			//B'(x) = (1 - coth(x) * coth(x)) + 1 / x^2
			std::function<double(DBL2)> F11diff = [&](DBL2 x) -> double {

				double y = c * (tau_ii.i * x.i + x.j * tau_ij.i) + d1;
				double coth = (1 + exp(-2 * y)) / (1 - exp(-2 * y));
				return ((1 - coth * coth) + 1 / (y*y)) * c * tau_ii.i - 1;
			};

			//This is F22'(x1, x2) = B'((x2 * tau2 + x1 * tau21) * c + mu2 * d) * c * tau2 - 1
			std::function<double(DBL2)> F22diff = [&](DBL2 x) -> double {

				double y = c * (tau_ii.j * x.j + x.i * tau_ij.j) + d2;
				double coth = (1 + exp(-2 * y)) / (1 - exp(-2 * y));
				return ((1 - coth * coth) + 1 / (y*y)) * c * tau_ii.j - 1;
			};

			//This is F12'(x1, x2) = B'((x1 * tau1 + x2 * tau12) * c + mu1 * d) * c * tau12
			std::function<double(DBL2)> F12diff = [&](DBL2 x) -> double {

				double y = c * (tau_ii.i * x.i + x.j * tau_ij.i) + d1;
				double coth = (1 + exp(-2 * y)) / (1 - exp(-2 * y));
				return ((1 - coth * coth) + 1 / (y*y)) * c * tau_ij.i;
			};

			//This is F21'(x1, x2) = B'((x2 * tau2 + x1 * tau21) * c + mu2 * d) * c * tau21
			std::function<double(DBL2)> F21diff = [&](DBL2 x) -> double {

				double y = c * (tau_ii.j * x.j + x.i * tau_ij.j) + d2;
				double coth = (1 + exp(-2 * y)) / (1 - exp(-2 * y));
				return ((1 - coth * coth) + 1 / (y*y)) * c * tau_ij.j;
			};

			//next starting value for the root finder is me - this makes the NewtonRaphson algorithm extremely efficient (typically a single iteration is required to keep within 1e-6 accuracy!)

			DBL2 me_previous = me;

			me = Root_NewtonRaphson(F1, F2, F11diff, F12diff, F21diff, F22diff, me, 1e-9, 2);

			//me needs to be decreasing and positive : when me is very small, numerical noise becomes a big factor (above Tc for small fields).
			if (me.i > me_previous.i || me.i < 0) me.i = me_previous.i;
			if (me.j > me_previous.j || me.j < 0) me.j = me_previous.j;

			values[t_idx] = me.i;
		}
	}

	//this function calculates me2 only and stores it here so we can access it separately
	void Initialize_CurieWeiss2(DBL2 tau_ii, DBL2 tau_ij, DBL2 d, double Tc)
	{
		if (values.size() != resolution * 2 + 1) values.resize(resolution * 2 + 1);

		DBL2 me = DBL2(1.0);
		values[0] = me.j;

		double Tc_renorm = (tau_ii.i + tau_ii.j + sqrt(pow(tau_ii.i - tau_ii.j, 2) + 4 * tau_ij.i * tau_ij.j)) / 2;

		for (int t_idx = 1; t_idx <= resolution * 2; t_idx++) {

			double T = (double)t_idx / resolution;

			double c = 3 / (T * Tc_renorm);
			double d1 = d.i * Tc_renorm / (T * Tc);
			double d2 = d.j * Tc_renorm / (T * Tc);

			//This is F1(x1, x2) = B((x1 * tau1 + x2 * tau12) * c + mu1 * d) - x1
			std::function<double(DBL2)> F1 = [&](DBL2 x) -> double {

				double y = c * (tau_ii.i * x.i + x.j * tau_ij.i) + d1;
				return ((1 + exp(-2 * y)) / (1 - exp(-2 * y))) - 1 / y - x.i;
			};

			//This is F2(x1, x2) = B((x2 * tau2 + x1 * tau21) * c + mu2 * d) - x2
			std::function<double(DBL2)> F2 = [&](DBL2 x) -> double {

				double y = c * (tau_ii.j * x.j + x.i * tau_ij.j) + d2;
				return ((1 + exp(-2 * y)) / (1 - exp(-2 * y))) - 1 / y - x.j;
			};

			//This is F11'(x1, x2) = B'((x1 * tau1 + x2 * tau12) * c + mu1 * d) * c * tau1 - 1
			//B'(x) = (1 - coth(x) * coth(x)) + 1 / x^2
			std::function<double(DBL2)> F11diff = [&](DBL2 x) -> double {

				double y = c * (tau_ii.i * x.i + x.j * tau_ij.i) + d1;
				double coth = (1 + exp(-2 * y)) / (1 - exp(-2 * y));
				return ((1 - coth * coth) + 1 / (y*y)) * c * tau_ii.i - 1;
			};

			//This is F22'(x1, x2) = B'((x2 * tau2 + x1 * tau21) * c + mu2 * d) * c * tau2 - 1
			std::function<double(DBL2)> F22diff = [&](DBL2 x) -> double {

				double y = c * (tau_ii.j * x.j + x.i * tau_ij.j) + d2;
				double coth = (1 + exp(-2 * y)) / (1 - exp(-2 * y));
				return ((1 - coth * coth) + 1 / (y*y)) * c * tau_ii.j - 1;
			};

			//This is F12'(x1, x2) = B'((x1 * tau1 + x2 * tau12) * c + mu1 * d) * c * tau12
			std::function<double(DBL2)> F12diff = [&](DBL2 x) -> double {

				double y = c * (tau_ii.i * x.i + x.j * tau_ij.i) + d1;
				double coth = (1 + exp(-2 * y)) / (1 - exp(-2 * y));
				return ((1 - coth * coth) + 1 / (y*y)) * c * tau_ij.i;
			};

			//This is F21'(x1, x2) = B'((x2 * tau2 + x1 * tau21) * c + mu2 * d) * c * tau21
			std::function<double(DBL2)> F21diff = [&](DBL2 x) -> double {

				double y = c * (tau_ii.j * x.j + x.i * tau_ij.j) + d2;
				double coth = (1 + exp(-2 * y)) / (1 - exp(-2 * y));
				return ((1 - coth * coth) + 1 / (y*y)) * c * tau_ij.j;
			};

			//next starting value for the root finder is me - this makes the NewtonRaphson algorithm extremely efficient (typically a single iteration is required to keep within 1e-6 accuracy!)

			DBL2 me_previous = me;

			me = Root_NewtonRaphson(F1, F2, F11diff, F12diff, F21diff, F22diff, me, 1e-9, 2);

			//me needs to be decreasing and positive : when me is very small, numerical noise becomes a big factor (above Tc for small fields).
			if (me.i > me_previous.i || me.i < 0) me.i = me_previous.i;
			if (me.j > me_previous.j || me.j < 0) me.j = me_previous.j;

			values[t_idx] = me.j;
		}
	}

	///////////////////////////////////////////////////////////////
	// LONGITUDINAL RELATIVE SUSCEPTIBILITY FUNCTION (i.e. normalised to mu0 * Ms0)

	//given by dMe / dHa at Ha = 0, where Me = me * Ms0, Ms0 being the zero temperature saturation magnetization
	//Thus suspar = mu0 Ms0 d0 B'(c me) / (1 - cB'(c me)), where d0 = mu/kBT
	//divide by mu0 Ms0 to remove the Ms0 dependency, thus obtaining the longitudinal relative susceptibility as a scaling function with temperature
	//This is calculated to a normalized Tc value of 1. Thus when it is used for a real Tc value, as chi(T/Tc), we need to divide it here by Tc so the output chi(T/Tc) function is correct
	void Initialize_LongitudinalRelSusceptibility(std::vector<double>& me_values, double atomic_moment, double Tc)
	{
		if (values.size() != resolution * 2 + 1) values.resize(resolution * 2 + 1);

		//This is B'(x) = 1 - coth^2(x) + 1/x^2
		auto Bdiff = [](double x) -> double {

			double coth = (1 + exp(-2 * x)) / (1 - exp(-2 * x));
			return ((1 - coth * coth) + 1 / (x*x));
		};

		for (int t_idx = 1; t_idx <= resolution * 2; t_idx++) {

			double t_value = (double)t_idx / resolution;

			double c = 3 / t_value;
			double d0 = (MUB / BOLTZMANN) * atomic_moment / (t_value * Tc);

			//susrel scaling
			if (t_idx <= resolution) {

				values[t_idx] = d0 * Bdiff(c * me_values[t_idx]) / (1.0 - c * Bdiff(c * me_values[t_idx]));
			}
			else {

				//Above Tc numerical noise is a problem : susrel must be positive and not increasing so use this info
				double value = d0 * Bdiff(c * me_values[t_idx]) / (1.0 - c * Bdiff(c * me_values[t_idx]));
				
				if (value > 0 && value <= values[t_idx - 1]) values[t_idx] = value;
				else values[t_idx] = values[t_idx - 1];
			}
		}

		//don't set it to zero at T = 0 to avoid NaNs when used in simulations
		values[0] = values[1];
	}

	///////////////////////////////////////////////////////////////
	// LONGITUDINAL RELATIVE SUSCEPTIBILITY FUNCTION (i.e. normalised to mu0 * Ms0) FOR 2-SUBLATTICE MODEL

	//relative longitudinal suscpeptibility in the 2-sublattice model : component 1

	//B1 = B[ (me1*tau1 + me2*tau12) * 3r/T ] = B[ (me1*tau1 + me2*tau12) * c ], for T in 0 to 2
	//B2 = B[ (me2*tau2 + me1*tau21) * 3r/T ] = B[ (me2*tau2 + me1*tau21) * c ], for T in 0 to 2

	//rTc is the renormalized Tc value as rTc = Tc / ((tau1 + tau2 + sqrt((tau1 - tau2)^2 + 4 * tau12 * tau21)) / 2);
	//Set Tc = 1 (normalized Curie-Weiss) and calculate between 0 and 2 (2 times Tc).

	void Initialize_LongitudinalRelSusceptibility1(std::vector<double>& me1_values, std::vector<double>& me2_values, DBL2 tau_ii, DBL2 tau_ij, DBL2 mu, double Tc)
	{
		if (values.size() != resolution * 2 + 1) values.resize(resolution * 2 + 1);

		double tau1 = tau_ii.i;
		double tau2 = tau_ii.j;
		double tau12 = tau_ij.i;
		double tau21 = tau_ij.j;

		double Tc_renorm = (tau1 + tau2 + sqrt(pow(tau1 - tau2, 2) + 4 * tau12 * tau21)) / 2;

		for (int t_idx = 1; t_idx <= resolution * 2; t_idx++) {

			double T = (double)t_idx / resolution;

			double c = 3 / (T * Tc_renorm);

			double d1 = (MUB / BOLTZMANN) * mu.i * Tc_renorm / (T * Tc);
			double d2 = (MUB / BOLTZMANN) * mu.j * Tc_renorm / (T * Tc);

			//B'(x) = (1 - coth(x) * coth(x)) + 1 / x^2
			std::function<double(DBL2)> B1diff = [&](DBL2 x) -> double {

				double y = c * (tau1 * x.i + x.j * tau12);
				double coth = (1 + exp(-2 * y)) / (1 - exp(-2 * y));
				return ((1 - coth * coth) + 1 / (y*y));
			};

			std::function<double(DBL2)> B2diff = [&](DBL2 x) -> double {

				double y = c * (tau2 * x.j + x.i * tau21);
				double coth = (1 + exp(-2 * y)) / (1 - exp(-2 * y));
				return ((1 - coth * coth) + 1 / (y*y));
			};

			DBL2 me = DBL2(me1_values[t_idx], me2_values[t_idx]);

			//susrel scaling
			if (t_idx <= resolution) {

				values[t_idx] = (d1 * B1diff(me) * (1 - c * tau2*B2diff(me)) + d2 * c*tau12*B1diff(me)*B2diff(me)) / ((1 - c * tau1*B1diff(me)) * (1 - c * tau2*B2diff(me)) - c * c*tau12*tau21*B1diff(me)*B2diff(me));
			}
			else {

				//Above Tc numerical noise is a problem : susrel must be positive and not increasing so use this info
				double value = (d1 * B1diff(me) * (1 - c * tau2*B2diff(me)) + d2 * c*tau12*B1diff(me)*B2diff(me)) / ((1 - c * tau1*B1diff(me)) * (1 - c * tau2*B2diff(me)) - c * c*tau12*tau21*B1diff(me)*B2diff(me));

				if (value > 0 && value <= values[t_idx - 1]) values[t_idx] = value;
				else values[t_idx] = values[t_idx - 1];
			}
		}

		//don't set it to zero at T = 0 to avoid NaNs when used in simulations
		values[0] = values[1];
	}

	//relative longitudinal suscpeptibility in the 2-sublattice model : component 2
	void Initialize_LongitudinalRelSusceptibility2(std::vector<double>& me1_values, std::vector<double>& me2_values, DBL2 tau_ii, DBL2 tau_ij, DBL2 mu, double Tc)
	{
		if (values.size() != resolution * 2 + 1) values.resize(resolution * 2 + 1);

		double tau1 = tau_ii.i;
		double tau2 = tau_ii.j;
		double tau12 = tau_ij.i;
		double tau21 = tau_ij.j;

		double Tc_renorm = (tau1 + tau2 + sqrt(pow(tau1 - tau2, 2) + 4 * tau12 * tau21)) / 2;

		for (int t_idx = 1; t_idx <= resolution * 2; t_idx++) {

			double T = (double)t_idx / resolution;

			double c = 3 / (T * Tc_renorm);

			double d1 = (MUB / BOLTZMANN) * mu.i * Tc_renorm / (T * Tc);
			double d2 = (MUB / BOLTZMANN) * mu.j * Tc_renorm / (T * Tc);

			//B'(x) = (1 - coth(x) * coth(x)) + 1 / x^2
			std::function<double(DBL2)> B1diff = [&](DBL2 x) -> double {

				double y = c * (tau1 * x.i + x.j * tau12);
				double coth = (1 + exp(-2 * y)) / (1 - exp(-2 * y));
				return ((1 - coth * coth) + 1 / (y*y));
			};

			std::function<double(DBL2)> B2diff = [&](DBL2 x) -> double {

				double y = c * (tau2 * x.j + x.i * tau21);
				double coth = (1 + exp(-2 * y)) / (1 - exp(-2 * y));
				return ((1 - coth * coth) + 1 / (y*y));
			};

			DBL2 me = DBL2(me1_values[t_idx], me2_values[t_idx]);

			//susrel scaling
			if (t_idx <= resolution) {

				values[t_idx] = (d2 * B2diff(me) * (1 - c * tau1*B1diff(me)) + d1 * c*tau21*B1diff(me)*B2diff(me)) / ((1 - c * tau1*B1diff(me)) * (1 - c * tau2*B2diff(me)) - c * c*tau12*tau21*B1diff(me)*B2diff(me));
			}
			else {

				//Above Tc numerical noise is a problem : susrel must be positive and not increasing so use this info
				double value = (d2 * B2diff(me) * (1 - c * tau1*B1diff(me)) + d1 * c*tau21*B1diff(me)*B2diff(me)) / ((1 - c * tau1*B1diff(me)) * (1 - c * tau2*B2diff(me)) - c * c*tau12*tau21*B1diff(me)*B2diff(me));

				if (value > 0 && value <= values[t_idx - 1]) values[t_idx] = value;
				else values[t_idx] = values[t_idx - 1];
			}
		}

		//don't set it to zero at T = 0 to avoid NaNs when used in simulations
		values[0] = values[1];
	}

	///////////////////////////////////////////////////////////////
	// TRANSVERSE DAMPING SCALING FOR 2-SUBLATTICE MODEL

	//1st component, Tc normalized to 1.
	void Initialize_Alpha1(std::vector<double>& me1_values, std::vector<double>& me2_values, DBL2 tau_ii, DBL2 tau_ij)
	{
		if (values.size() != resolution * 2 + 1) values.resize(resolution * 2 + 1);

		double tau1 = tau_ii.i;
		double tau2 = tau_ii.j;
		double tau12 = tau_ij.i;
		double tau21 = tau_ij.j;

		double Tc_renorm = (tau1 + tau2 + sqrt(pow(tau1 - tau2, 2) + 4 * tau12 * tau21)) / 2;

		for (int t_idx = 0; t_idx <= resolution * 2; t_idx++) {

			double T = (double)t_idx / resolution;

			if (T < 1.0) values[t_idx] = 1.0 - T * Tc_renorm / (3 * (tau1 + tau12 * me2_values[t_idx] / me1_values[t_idx]));
			else values[t_idx] = 2 * T / 3;
		}
	}

	//version where me1 and me2 are the same
	void Initialize_Alpha1(DBL2 tau_ii, DBL2 tau_ij)
	{
		if (values.size() != resolution * 2 + 1) values.resize(resolution * 2 + 1);

		double tau1 = tau_ii.i;
		double tau2 = tau_ii.j;
		double tau12 = tau_ij.i;
		double tau21 = tau_ij.j;

		double Tc_renorm = (tau1 + tau2 + sqrt(pow(tau1 - tau2, 2) + 4 * tau12 * tau21)) / 2;

		for (int t_idx = 0; t_idx <= resolution * 2; t_idx++) {

			double T = (double)t_idx / resolution;

			if (T < 1.0) values[t_idx] = 1.0 - T * Tc_renorm / (3 * (tau1 + tau12));
			else values[t_idx] = 2 * T / 3;
		}
	}

	//2nd component, Tc normalized to 1.
	void Initialize_Alpha2(std::vector<double>& me1_values, std::vector<double>& me2_values, DBL2 tau_ii, DBL2 tau_ij)
	{
		if (values.size() != resolution * 2 + 1) values.resize(resolution * 2 + 1);

		double tau1 = tau_ii.i;
		double tau2 = tau_ii.j;
		double tau12 = tau_ij.i;
		double tau21 = tau_ij.j;

		double Tc_renorm = (tau1 + tau2 + sqrt(pow(tau1 - tau2, 2) + 4 * tau12 * tau21)) / 2;

		for (int t_idx = 0; t_idx <= resolution * 2; t_idx++) {

			double T = (double)t_idx / resolution;

			if (T < 1.0) values[t_idx] = 1.0 - T * Tc_renorm / (3 * (tau2 + tau21 * me1_values[t_idx] / me2_values[t_idx]));
			else values[t_idx] = 2 * T / 3 ;
		}
	}

	//version where me1 and me2 are the same
	void Initialize_Alpha2(DBL2 tau_ii, DBL2 tau_ij)
	{
		if (values.size() != resolution * 2 + 1) values.resize(resolution * 2 + 1);

		double tau1 = tau_ii.i;
		double tau2 = tau_ii.j;
		double tau12 = tau_ij.i;
		double tau21 = tau_ij.j;

		double Tc_renorm = (tau1 + tau2 + sqrt(pow(tau1 - tau2, 2) + 4 * tau12 * tau21)) / 2;

		for (int t_idx = 0; t_idx <= resolution * 2; t_idx++) {

			double T = (double)t_idx / resolution;

			if (T < 1.0) values[t_idx] = 1.0 - T * Tc_renorm / (3 * (tau2 + tau21));
			else values[t_idx] = 2 * T / 3;
		}
	}

	///////////////////////////////////////////////////////////////

	double evaluate(double x) const
	{
		double findex = (x - start) * resolution;
		int index = (int)floor_epsilon(findex);

		if (index + 1 < values.size() && index >= 0) {

			return values[index] * (double(index + 1) - findex) + values[index + 1] * (findex - double(index));
		}
		else {

			if (index >= 0) return values.back();
			else return values[0];
		}
	}

	///////////////////////////////////////////////////////////////

	std::vector<double>& get_data(void) { return values; }

	double get_start(void) { return start; }
	int get_resolution(void) { return resolution; }
};
