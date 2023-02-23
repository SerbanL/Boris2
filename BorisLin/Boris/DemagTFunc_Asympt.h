#pragma once

#include "BorisLib.h"

/////////////////////////////////////////// ASYMPTOTIC EXPANSION APPROXIMATIONS FOR LONGER DISTANCES

class DemagAsymptoticDiag
{
	//various coefficients computed in constructor
	bool cubic_cell;

	double lead_weight;

	double a1, a2, a3, a4, a5, a6;
	double b1, b2, b3, b4, b5, b6, b7, b8, b9, b10;
	double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15;

public:

	DemagAsymptoticDiag(double hx, double hy, double hz);

	//Evaluation method
	double AsymptoticLdia(double x, double y, double z);
};

class DemagAsymptoticOffDiag
{
	//various coefficients computed in constructor
	bool cubic_cell;

	double lead_weight;

	double a1, a2, a3;
	double b1, b2, b3, b4, b5, b6;
	double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10;

public:

	DemagAsymptoticOffDiag(double hx, double hy, double hz);

	//Evaluation method
	double AsymptoticLodia(double x, double y, double z);
};