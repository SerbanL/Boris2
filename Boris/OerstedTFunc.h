#pragma once

#include "BorisLib.h"

#include <math.h>

namespace OeFunc {

	double G1(double x, double y, double z);

	double I2(double x, double y, double z, double lx, double ly, double lz);

	//Oersted tensor (rank-2) coefficients
	//Kyx = -Kxy
	//Kzy = -Kyz
	//Kxz = -Kzx
	//Kxx = Kyy = Kzz = 0
	double CalcKxy(DBL3 pos, DBL3 h);
	double CalcKyz(DBL3 pos, DBL3 h);
	double CalcKxz(DBL3 pos, DBL3 h);

	void CalcOerstedTensors(VEC<DBL3>& DOe, INT3 n, INT3 N, DBL3 h);
}
