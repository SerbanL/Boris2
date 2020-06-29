#pragma once

#include "VEC.h"

template <class CurveFitting_LMA>
class CurveFitting_Skyrmion
{

protected:

	/////////////// SKYRMION Z COMPONENT FUNCTION : Mz(x) = Ms * cos(2*arctan(sinh(R/w)/sinh((x-x0)/w))), x != x0, = -Ms for x = x0
	//Works for both skyrmion topological charges, but Ms allowed to be negative if core direction is +1.
	//R : radius
	//x0 : skyrmion center position
	//w : fitting factor, but in practice should be equal to w = PI * D / 4K, where K = Ku - mu0 Ms^2 / 2
	//Ms also used as a fitting factor, and this should come out to the expected value (same comment for w)
	//in practice we only extract R and x0, and w and Ms can be used to check they are as expected.

	//params[0] is R, params[1] is x0, params[2] is Ms, params[3] is w
	double GetSkyrmion(double x, std::vector<double> &params)
	{
		if (IsNZ(x - params[1])) return params[2] * cos(2 * atan(sinh(params[0] / params[3]) / sinh((x - params[1]) / params[3])));
		else return -params[2];
	}

	//Skyrmion function derivatives wrt parameters
	double GetSkyrmion_dR(double x, std::vector<double> &params)
	{
		if (IsNZ(x - params[1])) {

			double sinh_r_w = sinh((x - params[1]) / params[3]);
			double sinh_R_w = sinh(params[0] / params[3]);

			return -2.0 * params[2] * cosh(params[0] / params[3]) * sin(2 * atan(sinh_R_w / sinh_r_w)) / (params[3] * sinh_r_w * (1 + sinh_R_w * sinh_R_w / (sinh_r_w * sinh_r_w)));
		}
		else return 0.0;
	}

	double GetSkyrmion_dx0(double x, std::vector<double> &params)
	{
		if (IsNZ(x - params[1])) {

			double sinh_r_w = sinh((x - params[1]) / params[3]);
			double sinh_R_w = sinh(params[0] / params[3]);

			return -2.0 * (params[2] / params[3]) * sinh_R_w * cosh((x - params[1]) / params[3]) * sin(2 * atan(sinh_R_w / sinh_r_w)) / (sinh_R_w * sinh_R_w + sinh_r_w * sinh_r_w);
		}
		else return 0.0;
	}

	double GetSkyrmion_dMs(double x, std::vector<double> &params)
	{
		if (IsNZ(x - params[1])) return cos(2 * atan(sinh(params[0] / params[3]) / sinh((x - params[1]) / params[3])));
		else return -1.0;
	}

	double GetSkyrmion_dw(double x, std::vector<double> &params)
	{
		if (IsNZ(x - params[1])) {

			double sinh_r_w = sinh((x - params[1]) / params[3]);
			double sinh_R_w = sinh(params[0] / params[3]);
			double cosh_r_w = cosh((x - params[1]) / params[3]);
			double cosh_R_w = cosh(params[0] / params[3]);

			return -2.0 * params[2] * sin(2 * atan(sinh_R_w / sinh_r_w)) * ((x - params[1]) * sinh_R_w * cosh_r_w - params[0] * sinh_r_w * cosh_R_w)
				/ (params[3] * params[3] * sinh_r_w * sinh_r_w * (1 + sinh_R_w * sinh_R_w / (sinh_r_w * sinh_r_w)));
		}
		else return 0.0;
	}

	//Skyrmion fitting algorithm starting parameters setting method
	void Skyrmion_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params, int points)
	{
		//To estimate R find points where the z skyrmion profile crosses 0.
		//To estimate x0 find mid-way between zero crossings.
		//To estimate Ms just use first y value (this can be negative to allow skyrmions with core up to be fitted too)

		int num_switches = 0;

		double first_crossing = 0.0;
		double second_crossing = 0.0;

		for (int idx = 1; idx < (points > 0 && points <= xy.size() ? points : xy.size()); idx++) {

			if (xy[idx].j * xy[idx - 1].j <= 0) {

				num_switches++;

				if (num_switches == 1) first_crossing = xy[idx].x;
				else if (num_switches == 2) second_crossing = xy[idx].x;
				else break;
			}
		}

		//R estimation
		if (num_switches == 2) {

			//R
			params[0] = (second_crossing - first_crossing) / 2;
			//x0
			params[1] = (second_crossing + first_crossing) / 2;
		}

		else {

			//couldn't detect 2 switches, something is not right. default to sensible values
			//R
			params[0] = 30e-9;
			//x0
			params[1] = (xy.front().x + xy.back().x) / 2;
		}

		//Ms estimation
		params[2] = xy[0].y;

		//For w just set a reasonable starting value
		params[3] = 5e-9;
	}

	/////////////// SKYRMION Z COMPONENT FUNCTION WITH KNOWN CONSTANTS : Mz(x) = Ms * cos(2*arctan(sinh(R/w)/sinh((x-x0)/w))), x != x0, = -Ms for x = x0
	//Works for both skyrmion topological charges
	//R : radius
	//x0 : skyrmion center position
	//w : w = PI * |D| / 4|K|, where K = |Ku| - mu0 Ms^2 / 2
	//D, Ku, and Ms are given.

	//Skyrmion fitting algorithm starting parameters setting method
	//params[0] is R, params[1] is x0, params[2] is Ms adjusted for core sign, params[3] is w (params[2] and params[3] are determined from input parameters and not used as fitting constants
	void Skyrmion_FixedConstants_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params, int points)
	{
		//To estimate R find points where the z skyrmion profile crosses 0.
		//To estimate x0 find mid-way between zero crossings.
		//To estimate Ms just use first y value (this can be negative to allow skyrmions with core up to be fitted too)

		int num_switches = 0;

		double first_crossing = 0.0;
		double second_crossing = 0.0;

		for (int idx = 1; idx < (points > 0 && points <= xy.size() ? points : xy.size()); idx++) {

			if (xy[idx].j * xy[idx - 1].j <= 0) {

				num_switches++;

				if (num_switches == 1) first_crossing = xy[idx].x;
				else if (num_switches == 2) second_crossing = xy[idx].x;
				else break;
			}
		}

		//R estimation
		if (num_switches == 2) {

			//R
			params[0] = (second_crossing - first_crossing) / 2;
			//x0
			params[1] = (second_crossing + first_crossing) / 2;
		}

		else {

			//couldn't detect 2 switches, something is not right. default to sensible values
			//R
			params[0] = 30e-9;
			//x0
			params[1] = (xy.front().x + xy.back().x) / 2;
		}

		//params[2] and params[3] have already been set and are now fixed.

		//we just need to adjust params[2] for core sign.
		if (xy.front().y < 0) params[2] *= -1;
	}

	/////////////// SKYRMION X or Y COMPONENT FUNCTION : Mx(x) = Ms * sin(2*arctan(sinh(R/w)/sinh((x-x0)/w))), x != x0, = 0 for x = x0
	//Works for both skyrmion topological charges, but Ms allowed to be negative if core direction is +1.
	//R : radius
	//x0 : skyrmion center position
	//w : fitting factor, but in practice should be equal to w = PI * D / 4K, where K = Ku - mu0 Ms^2 / 2
	//Ms also used as a fitting factor, and this should come out to the expected value (same comment for w)
	//in practice we only extract R and x0, and w and Ms can be used to check they are as expected.

	//params[0] is R, params[1] is x0, params[2] is Ms, params[3] is w
	double GetSkyrmion_Longitudinal(double x, std::vector<double> &params)
	{
		if (IsNZ(x - params[1])) return params[2] * sin(2 * atan(sinh(params[0] / params[3]) / sinh((x - params[1]) / params[3])));
		else return 0.0;
	}

	//Skyrmion function derivatives wrt parameters
	double GetSkyrmion_Longitudinal_dR(double x, std::vector<double> &params)
	{
		if (IsNZ(x - params[1])) {

			double sinh_r_w = sinh((x - params[1]) / params[3]);
			double sinh_R_w = sinh(params[0] / params[3]);

			return 2.0 * params[2] * cosh(params[0] / params[3]) * cos(2 * atan(sinh_R_w / sinh_r_w)) / (params[3] * sinh_r_w * (1 + sinh_R_w * sinh_R_w / (sinh_r_w * sinh_r_w)));
		}
		else return 0.0;
	}

	double GetSkyrmion_Longitudinal_dx0(double x, std::vector<double> &params)
	{
		if (IsNZ(x - params[1])) {

			double sinh_r_w = sinh((x - params[1]) / params[3]);
			double sinh_R_w = sinh(params[0] / params[3]);

			return 2.0 * (params[2] / params[3]) * sinh_R_w * cosh((x - params[1]) / params[3]) * cos(2 * atan(sinh_R_w / sinh_r_w)) / (sinh_R_w * sinh_R_w + sinh_r_w * sinh_r_w);
		}
		else return 0.0;
	}

	double GetSkyrmion_Longitudinal_dMs(double x, std::vector<double> &params)
	{
		if (IsNZ(x - params[1])) return sin(2 * atan(sinh(params[0] / params[3]) / sinh((x - params[1]) / params[3])));
		else return 0.0;
	}

	double GetSkyrmion_Longitudinal_dw(double x, std::vector<double> &params)
	{
		if (IsNZ(x - params[1])) {

			double sinh_r_w = sinh((x - params[1]) / params[3]);
			double sinh_R_w = sinh(params[0] / params[3]);
			double cosh_r_w = cosh((x - params[1]) / params[3]);
			double cosh_R_w = cosh(params[0] / params[3]);

			return 2.0 * params[2] * cos(2 * atan(sinh_R_w / sinh_r_w)) * ((x - params[1]) * sinh_R_w * cosh_r_w - params[0] * sinh_r_w * cosh_R_w)
				/ (params[3] * params[3] * sinh_r_w * sinh_r_w * (1 + sinh_R_w * sinh_R_w / (sinh_r_w * sinh_r_w)));
		}
		else return 0.0;
	}

	//Skyrmion fitting algorithm starting parameters setting method
	void Skyrmion_Longitudinal_SetFittingStart(std::vector<DBL2> &xy, std::vector<double> &params, int points)
	{
		//To estimate R find points where the skyrmion profile has maximum, minimum
		//To estimate x0 find mid-way between the above extrema
		//To estimate Ms just use -1 * value of first extremum point

		int num_switches = 0;

		DBL2 maximum, minimum;

		for (int idx = 1; idx < (points > 0 && points <= xy.size() ? points : xy.size()); idx++) {

			if (maximum.y < xy[idx].y) maximum = xy[idx];
			if (minimum.y > xy[idx].y) minimum = xy[idx];
		}

		//R
		params[0] = abs(maximum.x - minimum.x) / 2;
		//x0
		params[1] = (maximum.x + minimum.x) / 2;

		//Ms estimation
		if (minimum.x < maximum.x) {

			params[2] = -1.0 * minimum.y;
		}
		else {

			params[2] = -1.0 * maximum.y;
		}

		//For w just set a reasonable starting value
		params[3] = 5e-9;
	}

public:

	CurveFitting_Skyrmion(void) {}

	//For the fitting functions below you can specify the number of points to use from xy - can be lower than xy.size(). : int points
	//You can also specify if the starting values for the parameters need to be determined, or if params already has reasonable values (e.g. from a previous fit if used repeatedly) : bool find_start (default true to find initial params values)

	/////////////// Fit Skyrmion Z Function to xy data : fill in params and their standard deviations

	//params[0] is R, params[1] is x0, params[2] is Ms, params[3] is w
	int FitSkyrmion_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &stdParams, int points = -1, bool find_start = true)
	{
		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_dR);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_dx0);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_dMs);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_dw);

		//returned fitting parameters and std in these
		params.resize(evalFunc.size() - 1);
		stdParams.resize(evalFunc.size() - 1);

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		if (find_start) {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, stdParams, evalFunc, &CurveFitting_LMA::Skyrmion_SetFittingStart, points);
		}
		else {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, stdParams, evalFunc, points);
		}
	}

	int FitSkyrmion_LMA(std::vector<DBL2> &xy, std::vector<double> &params, int points = -1, bool find_start = true)
	{
		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_dR);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_dx0);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_dMs);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_dw);

		//returned fitting parameters and std in these
		params.resize(evalFunc.size() - 1);

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		if (find_start) {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, evalFunc, &CurveFitting_LMA::Skyrmion_SetFittingStart, points);
		}
		else {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, evalFunc, points);
		}
	}

	//params[0] is R, params[1] is x0
	int FitSkyrmion_FixedConstants_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &stdParams, double D, double Ku, double Ms, int points = -1, bool find_start = true)
	{
		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_dR);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_dx0);

		//returned fitting parameters and std in these
		params.resize(evalFunc.size() - 1 + 2);
		stdParams.resize(evalFunc.size() - 1 + 2);

		double K = abs(abs(Ku) - MU0 * Ms*Ms / 2);
		double w = PI * abs(D) / (4 * K);

		//Ms value but not adjusted for core direction - fitting start function will adjust it if needed
		params[2] = Ms;
		params[3] = w;

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		if (find_start) {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, stdParams, evalFunc, &CurveFitting_LMA::Skyrmion_FixedConstants_SetFittingStart, points);
		}
		else {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, stdParams, evalFunc, points);
		}
	}

	int FitSkyrmion_FixedConstants_LMA(std::vector<DBL2> &xy, std::vector<double> &params, double D, double Ku, double Ms, int points = -1, bool find_start = true)
	{
		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_dR);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_dx0);

		//returned fitting parameters and std in these
		params.resize(evalFunc.size() - 1 + 2);

		double K = abs(abs(Ku) - MU0 * Ms*Ms / 2);
		double w = PI * abs(D) / (4 * K);

		//Ms value but not adjusted for core direction - fitting start function will adjust it if needed
		params[2] = Ms;
		params[3] = w;

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		if (find_start) {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, evalFunc, &CurveFitting_LMA::Skyrmion_FixedConstants_SetFittingStart, points);
		}
		else {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, evalFunc, points);
		}
	}

	/////////////// Fit Skyrmion X Function to xy data : fill in params and their standard deviations

	//params[0] is R, params[1] is x0, params[2] is Ms, params[3] is w
	int FitSkyrmion_Longitudinal_LMA(std::vector<DBL2> &xy, std::vector<double> &params, std::vector<double> &stdParams, int points = -1, bool find_start = true)
	{
		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_Longitudinal);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_Longitudinal_dR);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_Longitudinal_dx0);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_Longitudinal_dMs);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_Longitudinal_dw);

		//returned fitting parameters and std in these
		params.resize(evalFunc.size() - 1);
		stdParams.resize(evalFunc.size() - 1);

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		if (find_start) {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, stdParams, evalFunc, &CurveFitting_LMA::Skyrmion_Longitudinal_SetFittingStart, points);
		}
		else {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, stdParams, evalFunc, points);
		}
	}

	int FitSkyrmion_Longitudinal_LMA(std::vector<DBL2> &xy, std::vector<double> &params, int points = -1, bool find_start = true)
	{
		std::vector<typename CurveFitting_LMA::f_eval> evalFunc;

		//Build vector of function evaluations to be used in LMA algorithm
		//These functions are defined here but will be available in the derived class
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_Longitudinal);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_Longitudinal_dR);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_Longitudinal_dx0);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_Longitudinal_dMs);
		evalFunc.push_back(&CurveFitting_LMA::GetSkyrmion_Longitudinal_dw);

		//returned fitting parameters and std in these
		params.resize(evalFunc.size() - 1);

		//Use CRTP!
		CurveFitting_LMA& fitting_algorithm = static_cast<CurveFitting_LMA&>(*this);

		if (find_start) {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, evalFunc, &CurveFitting_LMA::Skyrmion_Longitudinal_SetFittingStart, points);
		}
		else {

			return fitting_algorithm.FitNonLinearFunction_LMA(xy, params, evalFunc, points);
		}
	}
};
