#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "DemagTFunc_Defs.h"
#include "BorisCUDALib.h"

/////////////////////////////////////////// ASYMPTOTIC EXPANSION APPROXIMATIONS FOR LONGER DISTANCES

class DemagAsymptoticDiagCUDA
{
	//various coefficients computed in constructor
	bool cubic_cell;

	double lead_weight;

	double a1, a2, a3, a4, a5, a6;
	double b1, b2, b3, b4, b5, b6, b7, b8, b9, b10;
	double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15;

public:

	//void constructor
	__host__ void construct_cu_obj(void) {}

	__host__ void destruct_cu_obj(void) {}

	//setup values depending on cellsize
	__host__ void setup(double hx, double hy, double hz);
	__device__ void setup_aux(double hx, double hy, double hz);

	//Evaluation method
	__device__ double AsymptoticLdia(double x, double y, double z);
};

class DemagAsymptoticOffDiagCUDA
{
	//various coefficients computed in constructor
	bool cubic_cell;

	double lead_weight;

	double a1, a2, a3;
	double b1, b2, b3, b4, b5, b6;
	double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10;

public:

	//void constructor
	__host__ void construct_cu_obj(void) {}

	__host__ void destruct_cu_obj(void) {}

	//setup values depending on cellsize
	__host__ void setup(double hx, double hy, double hz);
	__device__ void setup_aux(double hx, double hy, double hz);

	//Evaluation method
	__device__ double AsymptoticLodia(double x, double y, double z);
};

#endif