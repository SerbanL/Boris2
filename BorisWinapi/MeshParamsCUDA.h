#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "MaterialParameterCUDA.h"

class MeshParams;

class MeshParamsCUDA {

private:

	MeshParams *pmeshParams;

public:

	//Relative electron gyromagnetic ratio
	cu_obj<MatPCUDA<cuReal, cuReal>> grel;

	//Gilbert damping
	cu_obj<MatPCUDA<cuReal, cuReal>> alpha;

	//Saturation magnetisation (A/m)
	cu_obj<MatPCUDA<cuReal, cuReal>> Ms;

	//in-plane demagnetizing factors (used for Demag_N module)
	cu_obj<MatPCUDA<cuReal2, cuReal>> Nxy;

	//Exchange stifness (J/m)
	cu_obj<MatPCUDA<cuReal, cuReal>> A;

	//Dzyaloshinskii-Moriya exchange constant (J/m^2)
	cu_obj<MatPCUDA<cuReal, cuReal>> D;

	//bilinear surface exchange coupling (J/m^2) : J1, bottom and top layer values
	//biquadratic surface exchange coupling (J/m^2) : J2, bottom and top layer values
	cu_obj<MatPCUDA<cuReal, cuReal>> J1;
	cu_obj<MatPCUDA<cuReal, cuReal>> J2;

	//Magneto-crystalline anisotropy K1 and K2 constants (J/m^3) and easy axes directions. For uniaxial anisotropy only ea1 is needed, for cubic ea1 and ea2 should be orthogonal.
	cu_obj<MatPCUDA<cuReal, cuReal>> K1;
	cu_obj<MatPCUDA<cuReal, cuReal>> K2;
	cu_obj<MatPCUDA<cuReal3, cuReal3>> mcanis_ea1;
	cu_obj<MatPCUDA<cuReal3, cuReal3>> mcanis_ea2;

	//longitudinal (parallel) susceptibility relative to mu0*Ms0, i.e. divided by mu0*Ms0, Ms0 is the 0K Ms value - for use with LLB equation. Units As^2/kg
	cu_obj<MatPCUDA<cuReal, cuReal>> susrel;

	//perpendicular (transverse) susceptibility relative to mu0*Ms0, i.e. divided by mu0*Ms0, Ms0 is the 0K Ms value - for use with LLB equation. Units As^2/kg
	cu_obj<MatPCUDA<cuReal, cuReal>> susprel;

	//applied field spatial variation coefficient (unitless)
	cu_obj<MatPCUDA<cuReal, cuReal>> cHA;

	//electrical conductivity (units S/m).
	//this is the value at 0K for Ni80Fe20. Temperature dependence typically scaled by 1 / (1 + alpha*(T-T0)), where alpha = 0.003, T0 = 293K with sigma = 1.7e6 S/m and 293K.
	//Using scaling 1 / (1 + alpha0 * T) on the zero-temperature conductivity gives sigma0 = sigmaT0 / (1 - alpha*T0), alpha0 = alpha / (1 - alpha*T0), so alpha0 = 0.025.
	cu_obj<MatPCUDA<cuReal, cuReal>> elecCond;

	//anisotropic magnetoresistance as a percentage (of base resistance)
	cu_obj<MatPCUDA<cuReal, cuReal>> amrPercentage;

	//spin current polarization and non-adiabaticity (for Zhang-Li STT).
	cu_obj<MatPCUDA<cuReal, cuReal>> P;
	cu_obj<MatPCUDA<cuReal, cuReal>> beta;

	//parameters for spin current solver

	//electron diffusion constant (m^2/s)
	cu_obj<MatPCUDA<cuReal, cuReal>> De;

	//diffusion spin polarisation (unitless)
	cu_obj<MatPCUDA<cuReal, cuReal>> betaD;

	//spin Hall angle (unitless)
	cu_obj<MatPCUDA<cuReal, cuReal>> SHA;

	//"inverse" spin Hall angle (unitless) -> should normally be the same as SHA but e.g. can be set to zero to turn off the inverse SHE in the spin transport equation
	cu_obj<MatPCUDA<cuReal, cuReal>> iSHA;

	//field-like spin torque coefficient (unitless)
	cu_obj<MatPCUDA<cuReal, cuReal>> flSOT;

	//spin-flip length (m)
	cu_obj<MatPCUDA<cuReal, cuReal>> l_sf;

	//spin exchange rotation length (m)
	cu_obj<MatPCUDA<cuReal, cuReal>> l_ex;

	//spin dephasing length (m)
	cu_obj<MatPCUDA<cuReal, cuReal>> l_ph;

	//interface spin-dependent conductivity (spin-up and spin-down) (S/m^2)
	cu_obj<MatPCUDA<cuReal2, cuReal>> Gi;

	//interface spin-mixing conductivity (real and imaginary parts) (S/m^2)
	cu_obj<MatPCUDA<cuReal2, cuReal>> Gmix;

	//spin accumulation torque efficiency in the bulk (unitless, varies from 0 : no torque, up to 1 : full torque)
	cu_obj<MatPCUDA<cuReal, cuReal>> ts_eff;

	//spin accumulation torque efficiency at interfaces (unitless, varies from 0 : no torque, up to 1 : full torque)
	cu_obj<MatPCUDA<cuReal, cuReal>> tsi_eff;

	//spin pumping efficiency (unitless, varies from 0 : no spin pumping, up to 1 : full strength)
	cu_obj<MatPCUDA<cuReal, cuReal>> pump_eff;

	//the mesh base temperature (K)
	cu_obj<cuReal> base_temperature;

	//Curie temperture (K)
	cu_obj<cuReal> T_Curie;

	//thermal conductivity (W/mK) - default for permalloy
	cu_obj<MatPCUDA<cuReal, cuReal>> thermCond;

	//mass density (kg/m^3) - default for permalloy
	cu_obj<MatPCUDA<cuReal, cuReal>> density;

	//specific heat capacity (J/kgK) - default for permalloy
	cu_obj<MatPCUDA<cuReal, cuReal>> shc;

public:

	MeshParamsCUDA(MeshParams *pmeshParams);
	virtual ~MeshParamsCUDA();
};

#endif
