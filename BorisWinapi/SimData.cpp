#include "stdafx.h"
#include "SimData.h"

SimData::SimData(void) {

	SetDefaultParameters();
}

void SimData::SetDefaultParameters(void) {

	//Default parameters for permalloy where applicable (Ni80Fe20)

	alpha = 0.02;

	//A/m
	Ms = 8e5;

	Nx = 0; Ny = 0;

	//J/m
	A = 1.3e-11;

	//J/m^3
	K1 = 1e4;
	K2 = 0;
	mcanis_ea1 = DBL3(1,0,0);
	mcanis_ea2 = DBL3(0,1,0);
}