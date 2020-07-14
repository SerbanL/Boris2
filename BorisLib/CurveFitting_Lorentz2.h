#pragma once

#include "VEC.h"

template <class CurveFitting_LMA>
class CurveFitting_LorentzSA
{

protected:
	
	/////////////// LORENTZ PEAK FUNCTION : f(x) = y0 + S * (dH + A*(x-H0)) / ( 4*(x-H0)^2 + dH^2 )

	//params[0] is S, params[1] is A, params[2] is H0, params[3] is dH, params[4] is y0
	double GetLorentzSA(double x, std::vector<double> &params)
	{
		//Lorentzian peak function: f(x) = y0 + S * (dH + A*(x-H0)) / ( 4*(x-H0)^2 + dH^2 )
		//params[0] is S, params[1] is A, params[2] is H0, params[3] is dH, params[4] is y0

		double denom = 4 * (x - params[2])*(x - params[2]) + params[3] * params[3];

		if (denom) return (params[4] + params[0] * (params[3] + params[1] * (x - params[2])) / denom);
		else return 0;
	}

	//Lorentzian function derivatives wrt parameters
	double GetLorentzSA_dS(double x, std::vector<double> &params)
	{
		//Lorentzian peak function: f(x) = y0 + S * (dH + A*(x-H0)) / ( 4*(x-H0)^2 + dH^2 )
		//params[0] is S, params[1] is A, params[2] is H0, params[3] is dH, params[4] is y0

		double denom = 4 * (x - params[2])*(x - params[2]) + params[3] * params[3];

		if (denom) return ((params[3] + params[1] * (x - params[2])) / denom);
		else return 0;
	}

	double GetLorentzSA_dA(double x, std::vector<double> &params)
	{
		//Lorentzian peak function: f(x) = y0 + S * (dH + A*(x-H0)) / ( 4*(x-H0)^2 + dH^2 )
		//params[0] is S, params[1] is A, params[2] is H0, params[3] is dH, params[4] is y0

		double denom = 4 * (x - params[2])*(x - params[2]) + params[3] * params[3];

		if (denom) return (params[0]*(x - params[2]) / denom);
		else return 0;
	}

	double GetLorentzSA_dH0(double x, std::vector<double> &params)
	{
		//Lorentzian peak function: f(x) = y0 + S * (dH + A*(x-H0)) / ( 4*(x-H0)^2 + dH^2 )
		//params[0] is S, params[1] is A, params[2] is H0, params[3] is dH, params[4] is y0

		double denom = 4 * (x - params[2])*(x - params[2]) + params[3] * params[3];

		if (denom) return ((8 * params[0] * (params[3] + params[1] * (x - params[2])) * (x - params[2])) / (denom*denom) - params[0]*params[1] / denom);
		else return 0;
	}

	double GetLorentzSA_ddH(double x, std::vector<double> &params)
	{
		//Lorentzian peak function: f(x) = y0 + S * (dH + A*(x-H0)) / ( 4*(x-H0)^2 + dH^2 )
		//params[0] is S, params[1] is A, params[2] is H0, params[3] is dH, params[4] is y0

		double denom = 4 * (x - params[2])*(x - params[2]) + params[3] * params[3];

		if (denom) return (params[0] / denom - 2 * params[0] * params[3] * (params[3] + params[1] * (x - params[2])) / (denom*denom));
		else return 0;
	}

	double GetLorentzSA_dy0(double x, std::vector<double> &params) { return 1; }

	//Lorentzian fitting algorithm starting parameters setting method
	void LorentzSA_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params, int points)
	{
		//find starting values for parameters
		int idxMax = 0;
		double max = xy[0].j;

		for (int idx = 1; idx < (points > 0 && points <= xy.size() ? points : xy.size()); idx++) {

			if (max < xy[idx].j) {

				max = xy[idx].j;
				idxMax = idx;
			}
		}

		//starting y0: just take the first y point
		params[4] = xy[0].j;

		//starting H0 : value at which maximum y value occurs
		params[2] = xy[idxMax].i;

		//starting dH : just take H0/10 for simplicity
		params[3] = params[2] / 10;

		//starting S : calculate from current guess for H0 and f(H0) = max
		params[0] = (max - xy[0].j) * params[3];

		//starting A : just assume 0.0, in most cases this part is small compared to the symmetric part
		params[1] = 0.0;
	}

public:

	CurveFitting_LorentzSA(void) {}

	//For the fitting functions below you can specify the number of points to use from xy - can be lower than xy.size(). : int points
	//You can also specify if the starting values for the parameters need to be determined, or if params already has reasonable values (e.g. from a previous fit if used repeatedly) : bool find_start (default true to find initial params values)

	/////////////// Fit Lorentz Peak Function to xy data : fill in params and their standard deviations

	//params[0] is S, params[1] is H0, params[2] is dH, params[3] is y0
	int FitLorentzSA_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &stdParams, int points = -1, bool find_start = true)
	{
		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetLorentzSA);
		evalFunc.push_back(&CurveFitting_LMA::GetLorentzSA_dS);
		evalFunc.push_back(&CurveFitting_LMA::GetLorentzSA_dA);
		evalFunc.push_back(&CurveFitting_LMA::GetLorentzSA_dH0);
		evalFunc.push_back(&CurveFitting_LMA::GetLorentzSA_ddH);
		evalFunc.push_back(&CurveFitting_LMA::GetLorentzSA_dy0);

		//returned fitting parameters and std in these
		params.resize(evalFunc.size() - 1);
		stdParams.resize(evalFunc.size() - 1);

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		if (find_start) {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, stdParams, evalFunc, &CurveFitting_LMA::LorentzSA_SetFittingStart, points);
		}
		else {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, stdParams, evalFunc, points);
		}
	}

	int FitLorentzSA_LMA(std::vector<DBL2> &xy, std::vector<double> &params, int points = -1, bool find_start = true)
	{
		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetLorentzSA);
		evalFunc.push_back(&CurveFitting_LMA::GetLorentzSA_dS);
		evalFunc.push_back(&CurveFitting_LMA::GetLorentzSA_dA);
		evalFunc.push_back(&CurveFitting_LMA::GetLorentzSA_dH0);
		evalFunc.push_back(&CurveFitting_LMA::GetLorentzSA_ddH);
		evalFunc.push_back(&CurveFitting_LMA::GetLorentzSA_dy0);

		//returned fitting parameters and std in these
		params.resize(evalFunc.size() - 1);

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		if (find_start) {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, evalFunc, &CurveFitting_LMA::LorentzSA_SetFittingStart, points);
		}
		else {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, evalFunc, points);
		}
	}
};