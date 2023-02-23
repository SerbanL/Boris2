#pragma once

#include "VEC.h"

template <class CurveFitting_LMA>
class CurveFitting_SOT
{
	double SHAeff = 0.0;
	double std_SHAeff = 0.0;

	std::vector<double>* pdampinglike_x;
	std::vector<double>* pfieldlike_x;

protected:
	
	/////////////// Spin-Orbit Torque fitting function. Tsot = (SHAeff * MUB_E / dF) * [ (m ^ (m ^ p)) + flST * (m ^ p) ], where p = z ^ J
	//
	//to extract SHAeff and flST we separate the orthogonal damping-like and field-like terms and fit them separately
	//Thus first fit (m ^ (m ^ p)) * Ts = SHAeff * MUB_E / dF * |m ^ (m ^ p)|^2 -> this gives SHAeff
	//Next, using the obtained SHAeff, fit (m ^ p) * Ts = -(SHAeff * MUB_E / dF) * flST * |m ^ p|^2 -> this gives flST
	//Once we have beta we can go back and get the real P from P'

	//SHAeff fitting

	//params[0] is SHAeff
	double GetSOT_SHAeff(double x, std::vector<double> &params)
	{
		return params[0] * MUB_E * (*pdampinglike_x)[(int)x];
	}

	//params[0] is SHAeff
	double GetSOT_dSHAeff(double x, std::vector<double> &params)
	{
		return MUB_E * (*pdampinglike_x)[(int)x];
	}

	//params[0] is SHAeff
	void SOT_SHAeff_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params, int points)
	{
		//just take reasonable starting values and you'll be fine
		params[0] = 0.1;
	}

	//flST fitting

	//params[0] is flST
	double GetSOT_flST(double x, std::vector<double> &params)
	{
		return SHAeff * params[0] * MUB_E * (*pfieldlike_x)[(int)x];
	}

	//params[0] is flST
	double GetSOT_dflST(double x, std::vector<double> &params)
	{
		return SHAeff * MUB_E * (*pfieldlike_x)[(int)x];
	}

	//params[0] is flST
	void SOT_flST_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params, int points)
	{
		//just take reasonable starting values and you'll be fine
		params[0] = 0.1;
	}

public:

	CurveFitting_SOT(void) {}

	//For the fitting functions below you can specify the number of points to use from xy - can be lower than xy.size(). : int points
	//You can also specify if the starting values for the parameters need to be determined, or if params already has reasonable values (e.g. from a previous fit if used repeatedly) : bool find_start (default true to find initial params values)

	/////////////// Fit SOT Function

	//returned params[0] is SHAeff, params[1] is flST
	int FitSOT_LMA(
		std::vector<DBL2>& dampinglike_y, std::vector<DBL2>& fieldlike_y,
		std::vector<double>& dampinglike_x, std::vector<double>& fieldlike_x,
		std::vector<double> &params, std::vector<double> &stdParams, double *pRsq = nullptr)
	{
		if (dampinglike_y.size() != fieldlike_y.size()) return 0;

		pdampinglike_x = &dampinglike_x;
		pfieldlike_x = &fieldlike_x;

		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetSOT_SHAeff);
		evalFunc.push_back(&CurveFitting_LMA::GetSOT_dSHAeff);

		params.resize(1);
		stdParams.resize(1);

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		//Find SHAeff first

		fitting_algorithm.FitNonLinearFunction_LMA(dampinglike_y, params, stdParams, evalFunc, &CurveFitting_LMA::SOT_SHAeff_SetFittingStart, dampinglike_y.size());

		SHAeff = params[0];
		std_SHAeff = stdParams[0];
		
		//find flST next

		evalFunc.clear();
		evalFunc.push_back(&CurveFitting_LMA::GetSOT_flST);
		evalFunc.push_back(&CurveFitting_LMA::GetSOT_dflST);

		int iterations = fitting_algorithm.FitNonLinearFunction_LMA(fieldlike_y, params, stdParams, evalFunc, &CurveFitting_LMA::SOT_flST_SetFittingStart, fieldlike_y.size());
		
		params.insert(params.begin(), SHAeff);
		stdParams.insert(stdParams.begin(), std_SHAeff);
		
		if (pRsq) {

			std::vector<double> Y(dampinglike_y.size()), F(dampinglike_y.size());

#pragma omp parallel for
			for (int idx = 0; idx < dampinglike_y.size(); idx++) {

				Y[idx] = dampinglike_y[idx].y + fieldlike_y[idx].y;
				F[idx] = params[0] * MUB_E * (dampinglike_x[idx] + params[1] * fieldlike_x[idx]);
			}

			*pRsq = fitting_algorithm.Rsquared_Measure(Y, F, dampinglike_y.size());
		}

		return iterations;
	}

	int FitSOT_LMA(
		std::vector<DBL2>& dampinglike_y, std::vector<DBL2>& fieldlike_y,
		std::vector<double>& dampinglike_x, std::vector<double>& fieldlike_x,
		std::vector<double> &params)
	{
		if (dampinglike_y.size() != fieldlike_y.size()) return 0;

		pdampinglike_x = &dampinglike_x;
		pfieldlike_x = &fieldlike_x;

		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetSOT_SHAeff);
		evalFunc.push_back(&CurveFitting_LMA::GetSOT_dSHAeff);

		params.resize(1);

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		//Find SHAeff first

		fitting_algorithm.FitNonLinearFunction_LMA(dampinglike_y, params, evalFunc, &CurveFitting_LMA::SOT_SHAeff_SetFittingStart, dampinglike_y.size());

		SHAeff = params[0];

		//find flST next

		evalFunc.clear();
		evalFunc.push_back(&CurveFitting_LMA::GetSOT_flST);
		evalFunc.push_back(&CurveFitting_LMA::GetSOT_dflST);

		int iterations = fitting_algorithm.FitNonLinearFunction_LMA(fieldlike_y, params, evalFunc, &CurveFitting_LMA::SOT_flST_SetFittingStart, fieldlike_y.size());

		params.insert(params.begin(), SHAeff);

		return iterations;
	}
};