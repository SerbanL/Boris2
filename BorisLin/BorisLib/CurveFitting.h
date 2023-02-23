#pragma once

#include "CurveFitting_LMA.h"

class CurveFitting :
	public CurveFitting_LMA
{

public:

	CurveFitting(void) {}

	CurveFitting(double threshold_) :
		CurveFitting_LMA(threshold_)
	{}

	CurveFitting(double threshold_, int maxIterations_) :
		CurveFitting_LMA(threshold_, maxIterations_)
	{}
};