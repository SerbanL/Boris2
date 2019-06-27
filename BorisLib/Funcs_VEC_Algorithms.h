#pragma once

#include <algorithm>

#include "VEC.h"

class CurveFitting {

	//function pointer to function evaluations: input x value and a vector of parameters passed through reference. Return function value.
	typedef double(CurveFitting::*f_eval)(double, std::vector<double>&);

	//function pointer to function fitting parameters start: set initial parameters for a given function at the start of the fitting algorithm. Input (x, y) and parameters passed through reference.
	typedef void(CurveFitting::*fitStart)(std::vector<DBL2>&, std::vector<double>&);

#define CALLMETHOD(function) (this->*function)

private:

	//maximum number of iterations allowed for fitting
	int maxIterations = 10000;

	//convergence threshold for acceptance of fit
	double threshold = 1e-8;

private:

	/////////////// Levenberg-Marquardt algorithm

	//Non-linear curve fitting using the Levenberg-Marquardt algorithm. Also return standard deviation of fitting parameter error.
	//Input xi and yi vectors: the data to fit
	//params and stdParams will contain the function fitting parameter values and their standard deviation after the fit
	//fill evalFunc vector of function points with the actual function evaluation method and derivatives wrt fitting parameters
	//fitFunc is a function pointer to the method which sets the initial parameters to prime the algorithm

	void IterateLMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &paramsNext, std::vector<f_eval> evalFunc, double damping);
	double ResidualSumofSquares(std::vector<DBL2> &xy, std::vector<double> &params, f_eval evalFunc);
	int FitNonLinearFunction_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &stdParams, std::vector<f_eval> evalFunc, fitStart fitFunc);
	
	/////////////// LORENTZ PEAK FUNCTION : f(x) = y0 + S * dH / ( 4*(x-H0)^2 + dH^2 )

	//params[0] is S, params[1] is H0, params[2] is dH, params[3] is y0
	double GetLorentz(double x, std::vector<double> &params);

	//Lorentzian function derivatives wrt parameters
	double GetLorentz_dS(double x, std::vector<double> &params);
	double GetLorentz_dH0(double x, std::vector<double> &params);
	double GetLorentz_ddH(double x, std::vector<double> &params);
	double GetLorentz_dy0(double x, std::vector<double> &params) { return 1; }

	//Lorentzian fitting algorithm starting parameters setting method
	void Lorentz_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params);

public:

	CurveFitting(void) {}

	CurveFitting(double threshold_) :
		threshold(threshold_)
	{}

	CurveFitting(double threshold_, int maxIterations_) :
		threshold(threshold_), maxIterations(maxIterations_)
	{}

	/////////////// Fit Lorentz Peak Function to xy data : fill in params and their standard deviations

	int FitLorentz_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &stdParams);
};

///////////////////////////////////////////

//Lorentzian peak function: f(x) = y0 + S * dH / ( 4*(x-H0)^2 + dH^2 )
//params[0] is S, params[1] is H0, params[2] is dH, params[3] is y0
inline double CurveFitting::GetLorentz(double x, std::vector<double> &params)
{
	//Lorentzian peak function: f(x) = y0 + S * dH / ( 4*(x-H0)^2 + dH^2 )
	//params[0] is S, params[1] is H0, params[2] is dH, params[3] is y0

	double denom = 4 * (x - params[1])*(x - params[1]) + params[2] * params[2];

	if (denom) return (params[3] + params[0] * params[2] / denom);
	else return 0;
}

//Lorentzian function derivatives wrt parameters
inline double CurveFitting::GetLorentz_dS(double x, std::vector<double> &params)
{
	//Lorentzian peak function: f(x) = S * dH / ( 4*(x-H0)^2 + dH^2 )
	//params[0] is S, params[1] is H0, params[2] is dH

	double denom = 4 * (x - params[1])*(x - params[1]) + params[2] * params[2];

	if (denom) return (params[2] / denom);
	else return 0;
}

inline double CurveFitting::GetLorentz_dH0(double x, std::vector<double> &params)
{
	//Lorentzian peak function: f(x) = S * dH / ( 4*(x-H0)^2 + dH^2 )
	//params[0] is S, params[1] is H0, params[2] is dH

	double denom = 4 * (x - params[1])*(x - params[1]) + params[2] * params[2];

	if (denom) return (8 * params[0] * params[2] * (x - params[1]) / (denom*denom));
	else return 0;
}

inline double CurveFitting::GetLorentz_ddH(double x, std::vector<double> &params)
{
	//Lorentzian peak function: f(x) = S * dH / ( 4*(x-H0)^2 + dH^2 )
	//params[0] is S, params[1] is H0, params[2] is dH

	double denom = 4 * (x - params[1])*(x - params[1]) + params[2] * params[2];

	if (denom) return (params[0] * (4 * (x - params[1])*(x - params[1]) - params[2] * params[2]) / (denom*denom));
	else return 0;
}

//Lorentzian fitting algorithm starting parameters setting method
inline void CurveFitting::Lorentz_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params)
{
	//find starting values for parameters
	int idxMax = 0;
	double max = xy[0].j;

	for (int idx = 1; idx < xy.size(); idx++) {

		if (max < xy[idx].j) {

			max = xy[idx].j;
			idxMax = idx;
		}
	}

	//starting y0: just take the first y point
	params[3] = xy[0].j;

	//starting H0 : value at which maximum y value occurs
	params[1] = xy[idxMax].i;

	//starting dH : just take H0/10 for simplicity
	params[2] = params[1] / 10;

	//starting S : calculate from current guess for H0 and f(H0) = max
	params[0] = (max - xy[0].j) * params[2];
}

/////////////////////////////////////////////////////////////////////////////////

//Non-linear curve fitting using the Levenberg-Marquardt algorithm. Also return standard deviation of fitting parameter error.
//Input xi and yi vectors: the data to fit
//params and stdParams will contain the function fitting parameter values and their standard deviation after the fit
//fill evalFunc vector of function points with the actual function evaluation method and derivatives wrt fitting parameters
//fitFunc is a function pointer to the method which sets the initial parameters to prime the algorithm

inline void CurveFitting::IterateLMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &paramsNext, std::vector<f_eval> evalFunc, double damping)
{
	VEC<double> resid;
	VEC<double> Jacobian, Jacobian_T, Matrix_Prod;

	//number of data points to fit
	int points = xy.size();

	//number of fitting parameters
	int fParams = evalFunc.size() - 1;

	resid.resize(SZ3(1, points, 1));
	Jacobian.resize(SZ3(fParams, points, 1));
	Jacobian_T.resize(SZ3(points, fParams, 1));
	Matrix_Prod.resize(SZ3(fParams, fParams, 1));

	//Build Jacobian (points rows, fParams columns)
#pragma omp parallel for
	for (int i = 0; i < points; i++) {
		for (int j = 0; j < fParams; j++) {

			Jacobian[i*fParams + j] = CALLMETHOD(evalFunc[j + 1])(xy[i].i, params);
		}
	}

	//Jacobian transpose
	Jacobian.transpose_xy(Jacobian_T);

	//J_T * J : output fParams x fParams matrix
	Matrix_Prod.matrix_mul(Jacobian_T, Jacobian);

	//J_T * J + alpha*diag(J_T * J)
	Matrix_Prod.matrix_muldiag((1.0 + damping));

	//(J_T * J + alpha*diag(J_T * J))^-1
	Matrix_Prod.matrix_inverse();

	//Final Jacobian working matrix: (J_T * J + alpha*diag(J_T * J))^-1 * J_T
	//Jacobian will have dimensions points x fParams x 1 : points columns and fParams rows
	Jacobian.matrix_mul(Matrix_Prod, Jacobian_T);

	//Find residuals for current iteration parameters
#pragma omp parallel for
	for (int idx = 0; idx < points; idx++) {

		resid[idx] = (xy[idx].j - CALLMETHOD(evalFunc[0])(xy[idx].i, params));
	}

	VEC<double> vecparamsNext;

	//Multiply residuals with Jacobian working matrix : output size 1 x fParams
	vecparamsNext.matrix_mul(Jacobian, resid);

	//Form next iteration values for parameters
#pragma omp parallel for
	for (int idx = 0; idx < paramsNext.size(); idx++) {

		paramsNext[idx] = params[idx] + vecparamsNext[idx];
	}
}

inline double CurveFitting::ResidualSumofSquares(std::vector<DBL2> &xy, std::vector<double> &params, f_eval evalFunc)
{
	double rSos = 0;

	//number of data points to fit
	int points = xy.size();

	//Find residuals for current iteration parameters
#pragma omp parallel for reduction(+:rSos)
	for (int idx = 0; idx < points; idx++) {

		rSos += pow(xy[idx].j - CALLMETHOD(evalFunc)(xy[idx].i, params), 2);
	}

	return rSos;
}

inline int CurveFitting::FitNonLinearFunction_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &stdParams, std::vector<f_eval> evalFunc, fitStart fitFunc)
{
	//evalFunc[0] is the main function evaluation
	//evalFunc[1] ... up to evalFunc[number of parameters] are the evaluations of the dfferentials of the function wrt to the parameters.

	//number of data points to fit
	int points = xy.size();

	//number of fitting parameters
	int fParams = evalFunc.size() - 1;

	//must have at least fParams points : we have fParams fitting parameters
	if (points < fParams) return 0;

	std::vector<double> resid;
	VEC<double> Jacobian, Jacobian_T, Matrix_Prod;

	//used to calculate fitting parameters
	std::vector<double> paramsNext, paramsNext1, paramsNext2;

	resid.resize(points);

	Jacobian.resize(SZ3(fParams, points, 1));
	Jacobian_T.resize(SZ3(points, fParams, 1));
	Matrix_Prod.resize(SZ3(fParams, fParams, 1));

	params.resize(fParams);
	paramsNext.resize(fParams);
	paramsNext1.resize(fParams);
	paramsNext2.resize(fParams);

	//find starting values for parameters
	CALLMETHOD(fitFunc)(xy, params);

	int iter = 0;

	//LMA damping parameter
	double lam0 = 0.9, v = 1.01;

	for (iter = 0; iter < maxIterations; iter++) {

		//Find initial residual sum of squares (for current parameters)
		double rSos0 = ResidualSumofSquares(xy, params, evalFunc[0]);

		while (true) {

			//Find residual sum of squares after one iteration using lambda = lambda0 damping parameter
			IterateLMA(xy, params, paramsNext1, evalFunc, lam0);
			double rSos1 = ResidualSumofSquares(xy, paramsNext1, evalFunc[0]);

			//Find residual sum of squares for one iteration using lambda = lambda0/v damping parameter
			IterateLMA(xy, params, paramsNext2, evalFunc, lam0 / v);
			double rSos2 = ResidualSumofSquares(xy, paramsNext2, evalFunc[0]);

			//both evaluations are worse: set new damping parameter and continue
			if (rSos1 > rSos0 && rSos2 > rSos0) {

				lam0 *= v;
			}

			//Evaluation with lam0/v results in improvement: accept this as new damping and accept new parameters
			if (rSos2 <= rSos0) {

				lam0 = lam0 / v;
				paramsNext = paramsNext2;
				break;
			}

			//Evaluation with lam0 results in improvement: accept new parameters
			if (rSos1 <= rSos0) {

				paramsNext = paramsNext1;
				break;
			}
		}

		//calculate change in parameters from one iteration to the next and break if below set threshold
		double change = get_distance(params, paramsNext);
		double mag = GetMagnitude(params);
		if (mag) change /= mag;

		params = paramsNext;
		if (change < threshold) break;
	}

	//Find residuals for current iteration parameters
#pragma omp parallel for
	for (int idx = 0; idx < points; idx++) {

		resid[idx] = xy[idx].j - CALLMETHOD(evalFunc[0])(xy[idx].i, params);
	}

	//calculate standard error in fit parameters
	double stdy = GetMagnitude(resid) * sqrt((double)1.0 / (points - fParams + 1));

	//Build Jacobian (points rows, fParams columns)
#pragma omp parallel for
	for (int i = 0; i < points; i++) {
		for (int j = 0; j < fParams; j++) {

			Jacobian[i*fParams + j] = CALLMETHOD(evalFunc[j + 1])(xy[i].i, params);
		}
	}

	//Jacobian transpose (to fParams rows, points columns)
	Jacobian.transpose_xy(Jacobian_T);

	//Product is now (fParams rows, fParams columns)
	Matrix_Prod.matrix_mul(Jacobian_T, Jacobian);

	//[J_T * J]^-1
	Matrix_Prod.matrix_inverse();

	Matrix_Prod.matrix_getdiagonal(stdParams);

#pragma omp parallel for
	for (int idx = 0; idx < stdParams.size(); idx++)
		stdParams[idx] = sqrt(stdParams[idx]) * stdy;

	return iter;
}

/////////////////////////////////

inline int CurveFitting::FitLorentz_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &stdParams)
{
	std::vector<f_eval> evalFunc;

	//Build vector of function evaluations to be used in LMA algorithm
	evalFunc.push_back(&CurveFitting::GetLorentz);
	evalFunc.push_back(&CurveFitting::GetLorentz_dS);
	evalFunc.push_back(&CurveFitting::GetLorentz_dH0);
	evalFunc.push_back(&CurveFitting::GetLorentz_ddH);
	evalFunc.push_back(&CurveFitting::GetLorentz_dy0);

	return FitNonLinearFunction_LMA(xy, params, stdParams, evalFunc, &CurveFitting::Lorentz_SetFittingStart);
}