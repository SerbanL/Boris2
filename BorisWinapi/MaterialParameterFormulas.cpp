#include "stdafx.h"
#include "MaterialParameterFormulas.h"

//get function object for currently set selector
function<double(const MatPFormula&, double)> MatPFormula::get_function(void)
{
	switch (formula_selector) {

	case MATPFORM_NONE:
		return &MatPFormula::__MATPFORM_NONE;
		break;

	case MATPFORM_LINEAR:
		return &MatPFormula::__MATPFORM_LINEAR;
		break;

	case MATPFORM_PARABOLIC:
		return &MatPFormula::__MATPFORM_PARABOLIC;
		break;

	case MATPFORM_INVERSELINEAR:
		return &MatPFormula::__MATPFORM_INVERSELINEAR;
		break;
	}

	return &MatPFormula::__MATPFORM_NONE;
}

function<double(const MatPFormula&, double)> MatPFormula::get_set_function(MATPFORM_ formula_selector_)
{
	formula_selector = (int)formula_selector_;

	return get_function();
}

string MatPFormula::get_formula_info(void) const
{
	int num_params = formula_descriptor(formula_selector);
	string name = formula_descriptor.get_key_from_ID(formula_selector);

	string info = name;

	for (int idx = 0; idx < num_params; idx++)
		info += " " + ToString(coeff[idx]);

	return info;
}

void MatPFormula::set_formula_coeffs(vector<double> coefficients)
{
	//list any cases which need restrictions on values here - the default is to just copy the coefficients locally
	switch (formula_selector) {

	case MATPFORM_INVERSELINEAR:
		if (coefficients.size() && coefficients[0] >= 0) coeff[0] = coefficients[0];
		else coeff[0] = 0;
		break;

	default:
	{
		int numCoeffs = formula_descriptor(formula_selector);

		if (coefficients.size() == numCoeffs)
			copy(coefficients.begin(), coefficients.begin() + numCoeffs, coeff.begin());
	}
	break;
	}
}