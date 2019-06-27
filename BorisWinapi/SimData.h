#pragma once

#include "BorisLib.h"

using namespace std;

class SimData {

protected: //Protected data

	//Gilbert damping
	double alpha;
	
	//Saturation magnetisation (A/m)
	double Ms;
	
	//in-plane demagnetizing factors (used for Demag_N module)
	double Nx, Ny;

	//Exchange stifness (J/m)
	double A;

	//Magneto-crystalline anisotropy K1 and K2 constants (J/m^3) and easy axes directions. For uniaxial anisotropy only ea1 is needed, for cubic ea1 and ea2 should be orthogonal.
	double K1, K2;
	DBL3 mcanis_ea1, mcanis_ea2;

public:	 //Public data

private: //Private methods

	void SetDefaultParameters(void);

public:	 //Public methods

	SimData(void);
	~SimData() {}
};
