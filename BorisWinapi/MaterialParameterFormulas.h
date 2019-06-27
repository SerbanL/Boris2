#pragma once

#include "BorisLib.h"

#include "SimSharedData.h"

#include "ParametersDefs.h"

using namespace std;

class MatPFormula :
	public SimulationSharedData
{

protected:

	//selector for pre-defined temperature scaling formula (if this is MATPFORM_NONE then t_scaling is used if set)
	int formula_selector = MATPFORM_NONE;

	//coefficients
	vector<double> coeff = {0, 0};

private:

	//--------- MATPFORM_NONE

	double __MATPFORM_NONE(double T) const { return 1.0; }

	//--------- MATPFORM_LINEAR

	// y = x0 * T + 1
	double __MATPFORM_LINEAR(double T) const { return (coeff[0] * T + 1); }

	//--------- MATPFORM_PARABOLIC

	// y = x0 * T^2 + x1 * T + 1
	double __MATPFORM_PARABOLIC(double T) const { return (coeff[0] * T*T + coeff[1] * T + 1); }

	//--------- MATPFORM_INVERSELINEAR

	// y = 1 / (x0 * T + 1). coeff[0] not allowed to be negative.
	double __MATPFORM_INVERSELINEAR(double T) const { return 1 / (coeff[0] * T + 1); }

protected:

	//---------Constructors

	MatPFormula(void) {}

	//---------Assignment operator

	MatPFormula& operator=(const MatPFormula& copy_this)
	{
		formula_selector = copy_this.formula_selector;

		clear_vector(coeff);
		coeff = copy_this.coeff;

		return *this;
	}

	//---------Return function object depending on selector

	//get function object for currently set selector
	function<double(const MatPFormula&, double)> get_function(void);

	//set selector then get function object
	function<double(const MatPFormula&, double)> get_set_function(MATPFORM_ formula_selector_ = MATPFORM_NONE);

	//get string describing currently set formula (even if MATPFORM_NONE), including any set parameters; intended for console display
	string get_formula_info(void) const;

	//set coeffcients for currently selected formula and make any check for value restrictions
	void set_formula_coeffs(vector<double> coefficients);

public:

	int get_formula_selector(void) const { return formula_selector; }

	vector<double>& formula_coeff_ref(void) { return coeff; }
};