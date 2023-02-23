#pragma once

#include "VEC.h"

template <class CurveFitting_LMA>
class CurveFitting_DW
{

protected:

	/////////////// DOMAIN WALL TANH COMPONENT FUNCTION : M(x) = A * tanh((x-x0) / (D / PI))
	//Works on magnetization component which follows a tanh function as a function of profile position x
	//|A| : amplitude
	//x0 : domain wall centre
	//D : domain wall width

	//params[0] is A, params[1] is x0, params[2] is D
	double GetDW(double x, std::vector<double> &params)
	{
		return params[0] * tanh(PI * (x - params[1]) / params[2]);
	}

	//Domain wall function derivatives wrt parameters
	double GetDW_dA(double x, std::vector<double> &params)
	{
		return tanh(PI * (x - params[1]) / params[2]);
	}

	double GetDW_dx0(double x, std::vector<double> &params)
	{
		double th = tanh(PI * (x - params[1]) / params[2]);
		return params[0] * (th*th - 1.0)  * PI / params[2];
	}

	double GetDW_dD(double x, std::vector<double> &params)
	{
		double th = tanh(PI * (x - params[1]) / params[2]);
		return params[0] * (th*th - 1.0) * PI * (x - params[1]) / (params[2]* params[2]);
	}

	//Domain wall fitting algorithm starting parameters setting method
	void DW_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params, int points)
	{
		//To estimate A find value of last point
		//To estimate x0 find mid-way between zero crossings.
		//To estimate D just set it to one fifth the length of the profile : reasonable starting value is sufficient

		//A estimation
		params[0] = xy[xy.size() - 1].y;

		for (int idx = 1; idx < (points > 0 && points <= xy.size() ? points : xy.size()); idx++) {

			if (xy[idx].j * xy[idx - 1].j <= 0) {

				//x0 estimate
				params[1] = xy[idx].x;
				break;
			}
		}

		//For D just set a reasonable starting value
		params[2] = (xy[xy.size() - 1].x - xy[0].x) / 5;
	}

public:

	CurveFitting_DW(void) {}

	//For the fitting functions below you can specify the number of points to use from xy - can be lower than xy.size(). : int points
	//You can also specify if the starting values for the parameters need to be determined, or if params already has reasonable values (e.g. from a previous fit if used repeatedly) : bool find_start (default true to find initial params values)

	/////////////// Fit Domain Wall TANH Function to xy data : fill in params and their standard deviations

	//params[0] is A, params[1] is x0, params[2] is D
	int FitDW_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &stdParams, int points = -1, bool find_start = true)
	{
		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetDW);
		evalFunc.push_back(&CurveFitting_LMA::GetDW_dA);
		evalFunc.push_back(&CurveFitting_LMA::GetDW_dx0);
		evalFunc.push_back(&CurveFitting_LMA::GetDW_dD);

		//returned fitting parameters and std in these
		params.resize(evalFunc.size() - 1);
		stdParams.resize(evalFunc.size() - 1);

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		if (find_start) {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, stdParams, evalFunc, &CurveFitting_LMA::DW_SetFittingStart, points);
		}
		else {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, stdParams, evalFunc, points);
		}
	}

	int FitDW_LMA(std::vector<DBL2> &xy, std::vector<double> &params, int points = -1, bool find_start = true)
	{
		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetDW);
		evalFunc.push_back(&CurveFitting_LMA::GetDW_dA);
		evalFunc.push_back(&CurveFitting_LMA::GetDW_dx0);
		evalFunc.push_back(&CurveFitting_LMA::GetDW_dD);

		//returned fitting parameters in this
		params.resize(evalFunc.size() - 1);

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		if (find_start) {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, evalFunc, &CurveFitting_LMA::DW_SetFittingStart, points);
		}
		else {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, evalFunc, points);
		}
	}
};
