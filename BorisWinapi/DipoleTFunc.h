#pragma once

#include "BorisLib.h"

#include <math.h>

//NOTE : formulas for N factors here are equivalent to those obtained from Newell formulas, and also applicable for non-cubic cells

namespace DipoleTFunc {

	//Functions for calculating demag tensor for a general rectangular prism. 
	//(x, y, z) is the distance from the centre of the rectangular prism to the point where the field is obtained - DBL3 xyz
	//This is H = -N*M, N being the usual demagnetizing tensor and M the uniform magnetization of the prism.
	//(a, b, c) are the dimensions of the rectangular prism - DBL3 abc
	double Nxx(DBL3 xyz, DBL3 abc);
	double Nyy(DBL3 xyz, DBL3 abc);
	double Nzz(DBL3 xyz, DBL3 abc);

	double Nxy(DBL3 xyz, DBL3 abc);
	double Nxz(DBL3 xyz, DBL3 abc);
	double Nyz(DBL3 xyz, DBL3 abc);

	double fx(DBL3 xyz, DBL3 abc);
	double fy(DBL3 xyz, DBL3 abc);
	double fz(DBL3 xyz, DBL3 abc);

	double Fxy(DBL3 xyz, DBL3 abc);
	double Fxz(DBL3 xyz, DBL3 abc);
	double Fyz(DBL3 xyz, DBL3 abc);
}