#pragma once

#include <algorithm>

#include "BLib_VEC.h"

#include "CurveFitting_Lorentz.h"
#include "CurveFitting_Lorentz2.h"
#include "CurveFitting_Skyrmion.h"
#include "CurveFitting_STT.h"
#include "CurveFitting_SOT.h"
#include "CurveFitting_SOTSTT.h"
#include "CurveFitting_DW.h"

class CurveFitting_LMA :
	//CurveFitting_LMA is the generic LMA algorithm, and it inherits particular fitting functions here 
	//To use a particular fitting function just call one of the public functions in the classes below
	//These classes use CRTP to access CurveFitting_LMA data and methods, hence the template parameter
	public CurveFitting_Lorentz<CurveFitting_LMA>,
	public CurveFitting_LorentzSA<CurveFitting_LMA>,
	public CurveFitting_Skyrmion<CurveFitting_LMA>,
	public CurveFitting_STT<CurveFitting_LMA>,
	public CurveFitting_SOT<CurveFitting_LMA>,
	public CurveFitting_SOTSTT<CurveFitting_LMA>,
	public CurveFitting_DW<CurveFitting_LMA>
{
	//use friend declaration so I don't have to make the data members below public
	friend CurveFitting_Lorentz;
	friend CurveFitting_LorentzSA;
	friend CurveFitting_Skyrmion;
	friend CurveFitting_STT;
	friend CurveFitting_SOT;
	friend CurveFitting_SOTSTT;
	friend CurveFitting_DW;

private:

#define CALLMETHOD(function) (this->*function)

	//maximum number of iterations allowed for fitting
	int maxIterations = 1000;

	//convergence threshold for acceptance of fit
	double threshold = 1e-8;

	//function pointer to function evaluations: input x value and a vector of parameters passed through reference. Return function value.
	typedef double(CurveFitting_LMA::*f_eval)(double, std::vector<double>&);

	//function pointer to function fitting parameters start: set initial parameters for a given function at the start of the fitting algorithm. Input (x, y) and parameters passed through reference.
	typedef void(CurveFitting_LMA::*fitStart)(std::vector<DBL2>&, std::vector<double>&, int);

private:

	/////////////// Levenberg-Marquardt algorithm

	//Non-linear curve fitting using the Levenberg-Marquardt algorithm. Also return standard deviation of fitting parameter error.
	//Input xi and yi vectors: the data to fit
	//params and stdParams will contain the function fitting parameter values and their standard deviation after the fit
	//fill evalFunc vector of function points with the actual function evaluation method and derivatives wrt fitting parameters
	//fitFunc is a function pointer to the method which sets the initial parameters to prime the algorithm

	//if you want to use just the first n points from xy, set value in points : default -1 means use all points in xy.

	void IterateLMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &paramsNext, std::vector<f_eval> evalFunc, double damping, int points);
	double ResidualSumofSquares(std::vector<DBL2> &xy, std::vector<double> &params, f_eval evalFunc, int points);
	
	//Fitting algorithm start : parameters, std, and fitting start
	int FitNonLinearFunction_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &stdParams, std::vector<f_eval> evalFunc, fitStart fitFunc, int points);

	//Fitting algorithm start : parameters, std, no fitting start (params should already have reasonable values in)
	int FitNonLinearFunction_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &stdParams, std::vector<f_eval> evalFunc, int points);

	//Fitting algorithm start : parameters, and fitting start
	int FitNonLinearFunction_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<f_eval> evalFunc, fitStart fitFunc, int points);

	//Fitting algorithm start : parameters, no fitting start (params should already have reasonable values in)
	int FitNonLinearFunction_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<f_eval> evalFunc, int points);

	//Calculate an R^2 fitting measure between Y and F
	double Rsquared_Measure(std::vector<double>& Y, std::vector<double>& F, int points);

public:

	CurveFitting_LMA(void) {}

	CurveFitting_LMA(double threshold_) :
		threshold(threshold_)
	{}

	CurveFitting_LMA(double threshold_, int maxIterations_) :
		threshold(threshold_), maxIterations(maxIterations_)
	{}
};

/////////////////////////////////////////////////////////////////////////////////

//Non-linear curve fitting using the Levenberg-Marquardt algorithm. Also return standard deviation of fitting parameter error.
//Input xi and yi vectors: the data to fit
//params and stdParams will contain the function fitting parameter values and their standard deviation after the fit
//fill evalFunc vector of function points with the actual function evaluation method and derivatives wrt fitting parameters
//fitFunc is a function pointer to the method which sets the initial parameters to prime the algorithm

inline void CurveFitting_LMA::IterateLMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &paramsNext, std::vector<f_eval> evalFunc, double damping, int points)
{
	VEC<double> resid;
	VEC<double> Jacobian, Jacobian_T, Matrix_Prod;

	//number of data points to fit
	points = (points > 0 && points <= xy.size() ? points : xy.size());

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
	for (int idx = 0; idx < fParams; idx++) {

		paramsNext[idx] = params[idx] + vecparamsNext[idx];
	}
}

inline double CurveFitting_LMA::ResidualSumofSquares(std::vector<DBL2> &xy, std::vector<double> &params, f_eval evalFunc, int points)
{
	double rSos = 0;

	//number of data points to fit
	points = (points > 0 && points <= xy.size() ? points : xy.size());

	//Find residuals for current iteration parameters
#pragma omp parallel for reduction(+:rSos)
	for (int idx = 0; idx < points; idx++) {

		rSos += pow(xy[idx].j - CALLMETHOD(evalFunc)(xy[idx].i, params), 2);
	}

	return rSos;
}

//Fitting algorithm start : parameters, std, and fitting start
inline int CurveFitting_LMA::FitNonLinearFunction_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &stdParams, std::vector<f_eval> evalFunc, fitStart fitFunc, int points_)
{
	//evalFunc[0] is the main function evaluation
	//evalFunc[1] ... up to evalFunc[number of parameters] are the evaluations of the dfferentials of the function wrt to the parameters.

	//number of data points to fit
	int points = (points_ > 0 && points_ <= xy.size() ? points_ : xy.size());

	//number of fitting parameters
	int fParams = evalFunc.size() - 1;

	//must have at least fParams points : we have fParams fitting parameters
	if (points < fParams) return 0;

	//find starting values for parameters
	CALLMETHOD(fitFunc)(xy, params, points);

	return FitNonLinearFunction_LMA(xy, params, stdParams, evalFunc, points);
}

//Fitting algorithm start : parameters, std, no fitting start (params should already have reasonable values in)
inline int CurveFitting_LMA::FitNonLinearFunction_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &stdParams, std::vector<f_eval> evalFunc, int points_)
{
	//evalFunc[0] is the main function evaluation
	//evalFunc[1] ... up to evalFunc[number of parameters] are the evaluations of the dfferentials of the function wrt to the parameters.

	//number of data points to fit
	int points = (points_ > 0 && points_ <= xy.size() ? points_ : xy.size());

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

	paramsNext = params;
	paramsNext1 = params;
	paramsNext2 = params;

	int iter = 0;

	//LMA damping parameter
	double lam0 = 0.9, v = 1.01;

	for (iter = 0; iter < maxIterations; iter++) {

		//Find initial residual sum of squares (for current parameters)
		double rSos0 = ResidualSumofSquares(xy, params, evalFunc[0], points);

		while (true) {

			//Find residual sum of squares after one iteration using lambda = lambda0 damping parameter
			IterateLMA(xy, params, paramsNext1, evalFunc, lam0, points);
			double rSos1 = ResidualSumofSquares(xy, paramsNext1, evalFunc[0], points);

			//Find residual sum of squares for one iteration using lambda = lambda0/v damping parameter
			IterateLMA(xy, params, paramsNext2, evalFunc, lam0 / v, points);
			double rSos2 = ResidualSumofSquares(xy, paramsNext2, evalFunc[0], points);

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

	for (int idx = 0; idx < fParams; idx++) {

		stdParams[idx] = sqrt(stdParams[idx]) * stdy;
	}

	return iter;
}

//Fitting algorithm start : parameters, and fitting start
inline int CurveFitting_LMA::FitNonLinearFunction_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<f_eval> evalFunc, fitStart fitFunc, int points_)
{
	//evalFunc[0] is the main function evaluation
	//evalFunc[1] ... up to evalFunc[number of parameters] are the evaluations of the dfferentials of the function wrt to the parameters.

	//number of data points to fit
	int points = (points_ > 0 && points_ <= xy.size() ? points_ : xy.size());

	//number of fitting parameters
	int fParams = evalFunc.size() - 1;

	//must have at least fParams points : we have fParams fitting parameters
	if (points < fParams) return 0;

	//find starting values for parameters
	CALLMETHOD(fitFunc)(xy, params, points);

	return FitNonLinearFunction_LMA(xy, params, evalFunc, points);
}

//Fitting algorithm start : parameters, no fitting start (params should already have reasonable values in)
inline int CurveFitting_LMA::FitNonLinearFunction_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<f_eval> evalFunc, int points_)
{
	//evalFunc[0] is the main function evaluation
	//evalFunc[1] ... up to evalFunc[number of parameters] are the evaluations of the dfferentials of the function wrt to the parameters.

	//number of data points to fit
	int points = (points_ > 0 && points_ <= xy.size() ? points_ : xy.size());

	//number of fitting parameters
	int fParams = evalFunc.size() - 1;

	//must have at least fParams points : we have fParams fitting parameters
	if (points < fParams) return 0;

	//used to calculate fitting parameters
	std::vector<double> paramsNext, paramsNext1, paramsNext2;

	paramsNext = params;
	paramsNext1 = params;
	paramsNext2 = params;

	int iter = 0;

	//LMA damping parameter
	double lam0 = 0.9, v = 1.01;

	for (iter = 0; iter < maxIterations; iter++) {

		//Find initial residual sum of squares (for current parameters)
		double rSos0 = ResidualSumofSquares(xy, params, evalFunc[0], points);

		while (true) {

			//Find residual sum of squares after one iteration using lambda = lambda0 damping parameter
			IterateLMA(xy, params, paramsNext1, evalFunc, lam0, points);
			double rSos1 = ResidualSumofSquares(xy, paramsNext1, evalFunc[0], points);

			//Find residual sum of squares for one iteration using lambda = lambda0/v damping parameter
			IterateLMA(xy, params, paramsNext2, evalFunc, lam0 / v, points);
			double rSos2 = ResidualSumofSquares(xy, paramsNext2, evalFunc[0], points);

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

	return iter;
}

//Calculate an R^2 fitting measure between Y and F
inline double CurveFitting_LMA::Rsquared_Measure(std::vector<double>& Y, std::vector<double>& F, int points)
{
	//Find R^2 measure

	double Y_mean = 0.0;

#pragma omp parallel for reduction (+:Y_mean)
	for (int idx = 0; idx < points; idx++) {

		Y_mean += Y[idx] / points;
	}

	double Stot = 0.0, Sres = 0.0;

#pragma omp parallel for reduction (+:Stot, Sres)
	for (int idx = 0; idx < points; idx++) {

		Sres += (Y[idx] - F[idx]) * (Y[idx] - F[idx]);
		Stot += (Y[idx] - Y_mean) * (Y[idx] - Y_mean);
	}

	if (Stot) {

		return 1.0 - Sres / Stot;
	}
	else return 0.0;
}