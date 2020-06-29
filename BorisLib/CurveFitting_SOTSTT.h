#pragma once

#include "VEC.h"

template <class CurveFitting_LMA>
class CurveFitting_SOTSTT
{
	//1. First find P
	double P = 0.0;
	double std_P = 0.0;

	std::vector<double>* padiabatic_fit_main_x;
	std::vector<double>* padiabatic_fit_dampinglike_x;
	std::vector<double>* padiabatic_fit_fieldlike_x;

	//2. Next find beta
	double beta = 0.0;
	double std_beta = 0.0;

	std::vector<double>* pnonadiabatic_fit_main_x;
	std::vector<double>* pnonadiabatic_fit_dampinglike_x;
	std::vector<double>* pnonadiabatic_fit_fieldlike_x;

	//3. Next find SHAeff
	double SHAeff = 0.0;
	double std_SHAeff = 0.0;

	std::vector<double>* pdampinglike_fit_main_x;
	std::vector<double>* pdampinglike_fit_adiabatic_x;
	std::vector<double>* pdampinglike_fit_nonadiabatic_x;
	
	//4. Finally find flST
	double flST = 0.0;
	double std_flST = 0.0;

	std::vector<double>* pfieldlike_fit_main_x;
	std::vector<double>* pfieldlike_fit_adiabatic_x;
	std::vector<double>* pfieldlike_fit_nonadiabatic_x;

protected:
	
	/////////////// This combined STT and SOT fitting. Extract P, beta, SHAeff, flST in turn (see separate SOT and STT fitting for details).
	
	//1. P fitting

	//params[0] is P, params[1] is SHAeff, params[2] is flST
	double GetSOTSTT_P(double x, std::vector<double> &params)
	{
		return params[0] * MUB_E * (*padiabatic_fit_main_x)[(int)x] + params[1] * MUB_E * ( (*padiabatic_fit_dampinglike_x)[(int)x] + params[2] * (*padiabatic_fit_fieldlike_x)[(int)x] );
	}

	//params[0] is P, params[1] is SHAeff, params[2] is flST
	double GetSOTSTT_P_dP(double x, std::vector<double> &params)
	{
		return MUB_E * (*padiabatic_fit_main_x)[(int)x];
	}

	//params[0] is P, params[1] is SHAeff, params[2] is flST
	double GetSOTSTT_P_dSHAeff(double x, std::vector<double> &params)
	{
		return MUB_E * ((*padiabatic_fit_dampinglike_x)[(int)x] + params[2] * (*padiabatic_fit_fieldlike_x)[(int)x]);
	}

	//params[0] is P, params[1] is SHAeff, params[2] is flST
	double GetSOTSTT_P_dflST(double x, std::vector<double> &params)
	{
		return params[1] * MUB_E * (*padiabatic_fit_fieldlike_x)[(int)x];
	}

	//params[0] is P, params[1] is SHAeff, params[2] is flST
	void SOTSTT_P_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params, int points)
	{
		//just take reasonable starting values and you'll be fine
		params[0] = 0.4;
		params[1] = 0.1;
		params[2] = 0.1;
	}

	//2. beta fitting (have P)

	//params[0] is beta, params[1] is SHAeff, params[2] is flST
	double GetSOTSTT_beta(double x, std::vector<double> &params)
	{
		return -P * params[0] * MUB_E * (*pnonadiabatic_fit_main_x)[(int)x] + params[1] * MUB_E * ((*pnonadiabatic_fit_dampinglike_x)[(int)x] + params[2] * (*pnonadiabatic_fit_fieldlike_x)[(int)x]);
	}

	//params[0] is beta, params[1] is SHAeff, params[2] is flST
	double GetSOTSTT_beta_dbeta(double x, std::vector<double> &params)
	{
		return -P * MUB_E * (*pnonadiabatic_fit_main_x)[(int)x];
	}

	//params[0] is beta, params[1] is SHAeff, params[2] is flST
	double GetSOTSTT_beta_dSHAeff(double x, std::vector<double> &params)
	{
		return MUB_E * ((*pnonadiabatic_fit_dampinglike_x)[(int)x] + params[2] * (*pnonadiabatic_fit_fieldlike_x)[(int)x]);
	}

	//params[0] is beta, params[1] is SHAeff, params[2] is flST
	double GetSOTSTT_beta_dflST(double x, std::vector<double> &params)
	{
		return params[1] * MUB_E * (*pnonadiabatic_fit_fieldlike_x)[(int)x];
	}

	//params[0] is beta, params[1] is SHAeff, params[2] is flST
	void SOTSTT_beta_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params, int points)
	{
		//just take reasonable starting values and you'll be fine
		params[0] = 0.04;
		params[1] = 0.1;
		params[2] = 0.1;
	}

	//3. SHAeff fitting (have P, beta)

	//params[0] is SHAeff
	double GetSOTSTT_SHAeff(double x, std::vector<double> &params)
	{
		return P * MUB_E * ( (*pdampinglike_fit_adiabatic_x)[(int)x] - beta * (*pdampinglike_fit_nonadiabatic_x)[(int)x] ) + params[0] * MUB_E * (*pdampinglike_fit_main_x)[(int)x];
	}

	//params[0] is SHAeff
	double GetSOTSTT_SHAeff_dSHAeff(double x, std::vector<double> &params)
	{
		return MUB_E * (*pdampinglike_fit_main_x)[(int)x];
	}

	//params[0] is SHAeff
	void SOTSTT_SHAeff_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params, int points)
	{
		//just take reasonable starting values and you'll be fine
		params[0] = 0.1;
	}

	//4. SHAeff fitting (have P, beta, SHAeff)

	//params[0] is flST
	double GetSOTSTT_flST(double x, std::vector<double> &params)
	{
		return P * MUB_E * ((*pfieldlike_fit_adiabatic_x)[(int)x] - beta * (*pfieldlike_fit_nonadiabatic_x)[(int)x]) + SHAeff * params[0] * MUB_E * (*pfieldlike_fit_main_x)[(int)x];
	}

	//params[0] is flST
	double GetSOTSTT_flST_dflST(double x, std::vector<double> &params)
	{
		return SHAeff * MUB_E * (*pfieldlike_fit_main_x)[(int)x];
	}

	//params[0] is flST
	void SOTSTT_flST_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params, int points)
	{
		//just take reasonable starting values and you'll be fine
		params[0] = 0.1;
	}

public:

	CurveFitting_SOTSTT(void) {}

	//For the fitting functions below you can specify the number of points to use from xy - can be lower than xy.size(). : int points
	//You can also specify if the starting values for the parameters need to be determined, or if params already has reasonable values (e.g. from a previous fit if used repeatedly) : bool find_start (default true to find initial params values)

	/////////////// Fit STT Function

	//returned params[0] is P, params[1] is beta
	int FitSOTSTT_LMA(
		std::vector<DBL2>& adiabatic_y, std::vector<DBL2>& nonadiabatic_y,
		std::vector<DBL2>& dampinglike_y, std::vector<DBL2>& fieldlike_y,
		std::vector<double>& adiabatic_fit_main_x, std::vector<double>& adiabatic_fit_dampinglike_x, std::vector<double>& adiabatic_fit_fieldlike_x,
		std::vector<double>& nonadiabatic_fit_main_x, std::vector<double>& nonadiabatic_fit_dampinglike_x, std::vector<double>& nonadiabatic_fit_fieldlike_x,
		std::vector<double>& dampinglike_fit_main_x, std::vector<double>& dampinglike_fit_adiabatic_x, std::vector<double>& dampinglike_fit_nonadiabatic_x,
		std::vector<double>& fieldlike_fit_main_x, std::vector<double>& fieldlike_fit_adiabatic_x, std::vector<double>& fieldlike_fit_nonadiabatic_x,
		std::vector<double> &params, std::vector<double> &stdParams, double *pRsq = nullptr)
	{
		padiabatic_fit_main_x = &adiabatic_fit_main_x;
		padiabatic_fit_dampinglike_x = &adiabatic_fit_dampinglike_x;
		padiabatic_fit_fieldlike_x = &adiabatic_fit_fieldlike_x;

		pnonadiabatic_fit_main_x = &nonadiabatic_fit_main_x;
		pnonadiabatic_fit_dampinglike_x = &nonadiabatic_fit_dampinglike_x;
		pnonadiabatic_fit_fieldlike_x = &nonadiabatic_fit_fieldlike_x;

		pdampinglike_fit_main_x = &dampinglike_fit_main_x;
		pdampinglike_fit_adiabatic_x = &dampinglike_fit_adiabatic_x;
		pdampinglike_fit_nonadiabatic_x = &dampinglike_fit_nonadiabatic_x;

		pfieldlike_fit_main_x = &fieldlike_fit_main_x;
		pfieldlike_fit_adiabatic_x = &fieldlike_fit_adiabatic_x;
		pfieldlike_fit_nonadiabatic_x = &fieldlike_fit_nonadiabatic_x;

		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		//1. Find P

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_P);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_P_dP);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_P_dSHAeff);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_P_dflST);

		params.resize(evalFunc.size() - 1);
		stdParams.resize(evalFunc.size() - 1);

		fitting_algorithm.FitNonLinearFunction_LMA(adiabatic_y, params, stdParams, evalFunc, &CurveFitting_LMA::SOTSTT_P_SetFittingStart, adiabatic_y.size());

		P = params[0];
		std_P = stdParams[0];
		
		//2. Find beta (have P)

		evalFunc.clear();
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_beta);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_beta_dbeta);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_beta_dSHAeff);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_beta_dflST);

		params.resize(evalFunc.size() - 1);
		stdParams.resize(evalFunc.size() - 1);

		fitting_algorithm.FitNonLinearFunction_LMA(nonadiabatic_y, params, stdParams, evalFunc, &CurveFitting_LMA::SOTSTT_beta_SetFittingStart, nonadiabatic_y.size());
		
		beta = params[0];
		std_beta = stdParams[0];

		//3. Find SHAeff (have P, beta)

		evalFunc.clear();
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_SHAeff);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_SHAeff_dSHAeff);

		params.resize(evalFunc.size() - 1);
		stdParams.resize(evalFunc.size() - 1);

		fitting_algorithm.FitNonLinearFunction_LMA(dampinglike_y, params, stdParams, evalFunc, &CurveFitting_LMA::SOTSTT_SHAeff_SetFittingStart, dampinglike_y.size());

		SHAeff = params[0];
		std_SHAeff = stdParams[0];

		//4. Find flST (have P, beta, SHAeff)

		evalFunc.clear();
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_flST);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_flST_dflST);

		params.resize(evalFunc.size() - 1);
		stdParams.resize(evalFunc.size() - 1);

		int iterations = fitting_algorithm.FitNonLinearFunction_LMA(fieldlike_y, params, stdParams, evalFunc, &CurveFitting_LMA::SOTSTT_flST_SetFittingStart, fieldlike_y.size());

		flST = params[0];
		std_flST = stdParams[0];

		params.clear();
		stdParams.clear();

		params.push_back(P * (1 + beta * beta));
		params.push_back(beta);
		params.push_back(SHAeff);
		params.push_back(flST);

		stdParams.push_back(std_P * (1 + beta * beta));
		stdParams.push_back(std_beta);
		stdParams.push_back(std_SHAeff);
		stdParams.push_back(std_flST);
		
		if (pRsq) {

			std::vector<double> Y(adiabatic_y.size()), F(adiabatic_y.size());

#pragma omp parallel for
			for (int idx = 0; idx < adiabatic_y.size(); idx++) {

				Y[idx] = adiabatic_y[idx].y + nonadiabatic_y[idx].y + dampinglike_y[idx].y + fieldlike_y[idx].y;

				F[idx] = P * MUB_E * adiabatic_fit_main_x[idx] + SHAeff * MUB_E * (adiabatic_fit_dampinglike_x[idx] + flST * adiabatic_fit_fieldlike_x[idx]);
				F[idx] += -P * beta * MUB_E * nonadiabatic_fit_main_x[idx] + SHAeff * MUB_E * (nonadiabatic_fit_dampinglike_x[idx] + flST * nonadiabatic_fit_fieldlike_x[idx]);
				F[idx] += P * MUB_E * (dampinglike_fit_adiabatic_x[idx] - beta * dampinglike_fit_nonadiabatic_x[idx]) + SHAeff * MUB_E * dampinglike_fit_main_x[idx];
				F[idx] += P * MUB_E * (fieldlike_fit_adiabatic_x[idx] - beta * fieldlike_fit_nonadiabatic_x[idx]) + SHAeff * flST * MUB_E * fieldlike_fit_main_x[idx];
			}

			*pRsq = fitting_algorithm.Rsquared_Measure(Y, F, adiabatic_y.size());
		}
		
		return iterations;
	}

	int FitSOTSTT_LMA(
		std::vector<DBL2>& adiabatic_y, std::vector<DBL2>& nonadiabatic_y,
		std::vector<DBL2>& dampinglike_y, std::vector<DBL2>& fieldlike_y,
		std::vector<double>& adiabatic_fit_main_x, std::vector<double>& adiabatic_fit_dampinglike_x, std::vector<double>& adiabatic_fit_fieldlike_x,
		std::vector<double>& nonadiabatic_fit_main_x, std::vector<double>& nonadiabatic_fit_dampinglike_x, std::vector<double>& nonadiabatic_fit_fieldlike_x,
		std::vector<double>& dampinglike_fit_main_x, std::vector<double>& dampinglike_fit_adiabatic_x, std::vector<double>& dampinglike_fit_nonadiabatic_x,
		std::vector<double>& fieldlike_fit_main_x, std::vector<double>& fieldlike_fit_adiabatic_x, std::vector<double>& fieldlike_fit_nonadiabatic_x,
		std::vector<double> &params)
	{
		padiabatic_fit_main_x = &adiabatic_fit_main_x;
		padiabatic_fit_dampinglike_x = &adiabatic_fit_dampinglike_x;
		padiabatic_fit_fieldlike_x = &adiabatic_fit_fieldlike_x;

		pnonadiabatic_fit_main_x = &nonadiabatic_fit_main_x;
		pnonadiabatic_fit_dampinglike_x = &nonadiabatic_fit_dampinglike_x;
		pnonadiabatic_fit_fieldlike_x = &nonadiabatic_fit_fieldlike_x;

		pdampinglike_fit_main_x = &dampinglike_fit_main_x;
		pdampinglike_fit_adiabatic_x = &dampinglike_fit_adiabatic_x;
		pdampinglike_fit_nonadiabatic_x = &dampinglike_fit_nonadiabatic_x;

		pfieldlike_fit_main_x = &fieldlike_fit_main_x;
		pfieldlike_fit_adiabatic_x = &fieldlike_fit_adiabatic_x;
		pfieldlike_fit_nonadiabatic_x = &fieldlike_fit_nonadiabatic_x;

		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		//1. Find P

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_P);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_P_dP);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_P_dSHAeff);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_P_dflST);

		params.resize(evalFunc.size() - 1);

		fitting_algorithm.FitNonLinearFunction_LMA(adiabatic_y, params, evalFunc, &CurveFitting_LMA::SOTSTT_P_SetFittingStart, adiabatic_y.size());

		P = params[0];

		//2. Find beta (have P)

		evalFunc.clear();
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_beta);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_beta_dbeta);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_beta_dSHAeff);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_beta_dflST);

		params.resize(evalFunc.size() - 1);

		fitting_algorithm.FitNonLinearFunction_LMA(nonadiabatic_y, params, evalFunc, &CurveFitting_LMA::SOTSTT_beta_SetFittingStart, nonadiabatic_y.size());

		beta = params[0];

		//3. Find SHAeff (have P, beta)

		evalFunc.clear();
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_SHAeff);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_SHAeff_dSHAeff);

		params.resize(evalFunc.size() - 1);

		fitting_algorithm.FitNonLinearFunction_LMA(dampinglike_y, params, evalFunc, &CurveFitting_LMA::SOTSTT_SHAeff_SetFittingStart, dampinglike_y.size());

		SHAeff = params[0];

		//4. Find flST (have P, beta, SHAeff)

		evalFunc.clear();
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_flST);
		evalFunc.push_back(&CurveFitting_LMA::GetSOTSTT_flST_dflST);

		params.resize(evalFunc.size() - 1);

		int iterations = fitting_algorithm.FitNonLinearFunction_LMA(fieldlike_y, params, evalFunc, &CurveFitting_LMA::SOTSTT_flST_SetFittingStart, fieldlike_y.size());

		flST = params[0];

		params.clear();

		params.push_back(P * (1 + beta * beta));
		params.push_back(beta);
		params.push_back(SHAeff);
		params.push_back(flST);

		return iterations;
	}
};