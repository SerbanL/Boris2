#pragma once

#include "VEC.h"

template <class CurveFitting_LMA>
class CurveFitting_STT
{
	double P = 0.0;
	double std_P = 0.0;

	std::vector<double>* padiabatic_x;
	std::vector<double>* pnonadiabatic_x;

protected:
	
	/////////////// Spin-Transfer Torque fitting function. Ts = (P / (1 + beta*beta)) * MUB_E * [ (J.del)m - beta * m ^ (J.del)m ]; Fit a self-consistently calculated torque Ts for P and beta, assuming Ts takes the form of Zhang-Li STT only.
	//
	//to extract P and beta we separate the orthogonal adiabatic and non-adiabatic terms and fit them separately
	//Thus first fit (J.del)m * Ts = P' * MUB_E * [(J.del)m]^2 -> this gives P' = P / (1 +beta*beta)
	//Next, using the obtained P', fit (m^(J.del)m) * Ts = -P' * MUB_E * beta * [m^(J.del)m]^2 -> this gives beta
	//Once we have beta we can go back and get the real P from P'
	
	//P fitting

	//params[0] is P
	double GetSTT_P(double x, std::vector<double> &params)
	{
		return params[0] * MUB_E * (*padiabatic_x)[(int)x];
	}

	//params[0] is P
	double GetSTT_dP(double x, std::vector<double> &params)
	{
		return MUB_E * (*padiabatic_x)[(int)x];
	}

	//params[0] is P
	void STT_P_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params, int points)
	{
		//just take reasonable starting values and you'll be fine
		params[0] = 0.4;
	}

	//beta fitting

	//params[0] is beta
	double GetSTT_beta(double x, std::vector<double> &params)
	{
		return -P * params[0] * MUB_E * (*pnonadiabatic_x)[(int)x];
	}

	//params[0] is beta
	double GetSTT_dbeta(double x, std::vector<double> &params)
	{
		return -P * MUB_E * (*pnonadiabatic_x)[(int)x];
	}

	//params[0] is beta
	void STT_beta_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params, int points)
	{
		//just take reasonable starting values and you'll be fine
		params[0] = 0.04;
	}

public:

	CurveFitting_STT(void) {}

	//For the fitting functions below you can specify the number of points to use from xy - can be lower than xy.size(). : int points
	//You can also specify if the starting values for the parameters need to be determined, or if params already has reasonable values (e.g. from a previous fit if used repeatedly) : bool find_start (default true to find initial params values)

	/////////////// Fit STT Function

	//returned params[0] is P, params[1] is beta
	int FitSTT_LMA(
		std::vector<DBL2>& adiabatic_y, std::vector<DBL2>& nonadiabatic_y,
		std::vector<double>& adiabatic_x, std::vector<double>& nonadiabatic_x,
		std::vector<double> &params, std::vector<double> &stdParams, double *pRsq = nullptr)
	{
		if (adiabatic_y.size() != nonadiabatic_y.size()) return 0;

		padiabatic_x = &adiabatic_x;
		pnonadiabatic_x = &nonadiabatic_x;

		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetSTT_P);
		evalFunc.push_back(&CurveFitting_LMA::GetSTT_dP);

		params.resize(1);
		stdParams.resize(1);

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		//Find P first

		fitting_algorithm.FitNonLinearFunction_LMA(adiabatic_y, params, stdParams, evalFunc, &CurveFitting_LMA::STT_P_SetFittingStart, adiabatic_y.size());

		P = params[0];
		std_P = stdParams[0];
		
		//find beta next

		evalFunc.clear();
		evalFunc.push_back(&CurveFitting_LMA::GetSTT_beta);
		evalFunc.push_back(&CurveFitting_LMA::GetSTT_dbeta);

		int iterations = fitting_algorithm.FitNonLinearFunction_LMA(nonadiabatic_y, params, stdParams, evalFunc, &CurveFitting_LMA::STT_beta_SetFittingStart, nonadiabatic_y.size());
		
		double beta = params[0];

		//Adjust 
		params.insert(params.begin(), P * (1 + beta*beta));
		stdParams.insert(stdParams.begin(), std_P * (1 + beta*beta));
		
		if (pRsq) {

			std::vector<double> Y(adiabatic_y.size()), F(adiabatic_y.size());

#pragma omp parallel for
			for (int idx = 0; idx < adiabatic_y.size(); idx++) {

				Y[idx] = adiabatic_y[idx].y + nonadiabatic_y[idx].y;
				F[idx] = P * MUB_E * (adiabatic_x[idx] - params[1] * nonadiabatic_x[idx]);
			}

			*pRsq = fitting_algorithm.Rsquared_Measure(Y, F, adiabatic_y.size());
		}

		return iterations;
	}

	//returned params[0] is P, params[1] is beta
	int FitSTT_LMA(
		std::vector<DBL2>& adiabatic_y, std::vector<DBL2>& nonadiabatic_y,
		std::vector<double>& adiabatic_x, std::vector<double>& nonadiabatic_x,
		std::vector<double> &params)
	{
		if (adiabatic_y.size() != nonadiabatic_y.size()) return 0;

		padiabatic_x = &adiabatic_x;
		pnonadiabatic_x = &nonadiabatic_x;

		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetSTT_P);
		evalFunc.push_back(&CurveFitting_LMA::GetSTT_dP);

		params.resize(1);

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		//Find P first

		fitting_algorithm.FitNonLinearFunction_LMA(adiabatic_y, params, evalFunc, &CurveFitting_LMA::STT_P_SetFittingStart, adiabatic_y.size());

		P = params[0];

		//find beta next

		evalFunc.clear();
		evalFunc.push_back(&CurveFitting_LMA::GetSTT_beta);
		evalFunc.push_back(&CurveFitting_LMA::GetSTT_dbeta);

		int iterations = fitting_algorithm.FitNonLinearFunction_LMA(nonadiabatic_y, params, evalFunc, &CurveFitting_LMA::STT_beta_SetFittingStart, nonadiabatic_y.size());

		double beta = params[0];

		//Adjust 
		params.insert(params.begin(), P * (1 + beta * beta));

		return iterations;
	}
};