#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "MaterialParameterCUDA.h"
#include "ParametersDefs.h"

class Atom_MeshParams;

class Atom_MeshParamsCUDA {

private:

	Atom_MeshParams *pameshParams;

protected:

public:

	//-----------SIMPLE CUBIC

	//Relative electron gyromagnetic ratio
	cu_obj<MatPCUDA<cuBReal, cuBReal>> grel;

	//Gilbert damping (atomistic: intrinsic)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> alpha;

	//atomic moment (units of muB) - default for bcc Fe
	cu_obj<MatPCUDA<cuBReal, cuBReal>> mu_s;

	//Exchange constant (units of J) - default for bcc Fe
	cu_obj<MatPCUDA<cuBReal, cuBReal>> J;

	//DMI exchange constant : (units of J)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> D;

	//Interfacial DMI symmetry axis direction, used by vector interfacial DMI module
	cu_obj<MatPCUDA<cuReal3, cuReal3>> D_dir;

	//Surface exchange coupling, used by the surfexchange module to couple two spins on different meshes at the surface (units of J)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> Js;
	//Secondary surface exchange coupling constant, used for coupling atomistic meshes to micromagnetic 2-sublattice meshes.
	cu_obj<MatPCUDA<cuBReal, cuBReal>> Js2;

	//Magneto-crystalline anisotropy constants (J) and easy axes directions. For uniaxial anisotropy only ea1 is needed.
	cu_obj<MatPCUDA<cuBReal, cuBReal>> K1;
	cu_obj<MatPCUDA<cuBReal, cuBReal>> K2;
	cu_obj<MatPCUDA<cuBReal, cuBReal>> K3;

	//Magneto-crystalline anisotropy easy axes directions
	cu_obj<MatPCUDA<cuReal3, cuReal3>> mcanis_ea1;
	cu_obj<MatPCUDA<cuReal3, cuReal3>> mcanis_ea2;
	cu_obj<MatPCUDA<cuReal3, cuReal3>> mcanis_ea3;

	//tensorial anisotropy. each term is a contribution to the anisotropy energy density as d*a^n1 b^n2 c^n3. Here a = m.mcanis_ea1, b = m.mcanis_ea2, c = m.mcanis_ea3.
	//For 2nd order we aditionally multiply by K1, 4th order K2, 6th order K3. Any other orders d coefficient contains anisotropy energy density.
	//each DBL4 stores (d, n1, n2, n3), where d != 0, n1, n2, n3 >= 0, n1+n2+n3>0. Odd order terms allowed.
	cu_obj<cuVEC<cuReal4>> Kt;

	//-----------BCC (2 per unit cell)

	//-----------FCC (4 per unit cell)

	//-----------HCP (4 per effective unit cell)

	//-----------Others

	//in-plane demagnetizing factors (used for Atom_Demag_N module)
	cu_obj<MatPCUDA<cuReal2, cuBReal>> Nxy;

	//applied field spatial variation coefficient (unitless)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> cHA;

	//Magneto-Optical field strength (A/m)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> cHmo;

	//Stochasticity efficiency parameter
	cu_obj<MatPCUDA<cuBReal, cuBReal>> s_eff;

	//electrical conductivity (units S/m).
	//this is the value at RT for Ni80Fe20.
	cu_obj<MatPCUDA<cuBReal, cuBReal>> elecCond;

	//TMR RA products for parallel and antiparallel states (Ohms m^2)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> RAtmr_p;
	cu_obj<MatPCUDA<cuBReal, cuBReal>> RAtmr_ap;

	//anisotropic magnetoresistance as a percentage (of base resistance)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> amrPercentage;

	//spin current polarization and non-adiabaticity (for Zhang-Li STT).
	cu_obj<MatPCUDA<cuBReal, cuBReal>> P;
	cu_obj<MatPCUDA<cuBReal, cuBReal>> beta;

	//parameters for spin current solver

	//electron diffusion constant (m^2/s)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> De;

	//electron carrier density (1/m^3)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> n_density;

	//diffusion spin polarization (unitless)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> betaD;

	//spin Hall angle (unitless)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> SHA;

	//field-like spin torque coefficient (unitless)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> flSOT;

	//Slonczewski macrospin torques q+, q- parameters as in PRB 72, 014446 (2005) (unitless)
	cu_obj<MatPCUDA<cuReal2, cuBReal>> STq;

	//Slonczewski macrospin torques A, B parameters as in PRB 72, 014446 (2005) (unitless)
	cu_obj<MatPCUDA<cuReal2, cuBReal>> STa;

	//Slonczewski macrospin torques spin polarization unit vector as in PRB 72, 014446 (2005) (unitless)
	cu_obj<MatPCUDA<cuReal3, cuReal3>> STp;

	//spin-flip length (m)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> l_sf;

	//spin exchange rotation length (m)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> l_ex;

	//spin dephasing length (m)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> l_ph;

	//interface spin-dependent conductivity (spin-up and spin-down) (S/m^2)
	cu_obj<MatPCUDA<cuReal2, cuBReal>> Gi;

	//interface spin-mixing conductivity (real and imaginary parts) (S/m^2)
	cu_obj<MatPCUDA<cuReal2, cuBReal>> Gmix;

	//spin accumulation torque efficiency in the bulk (unitless, varies from 0 : no torque, up to 1 : full torque)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> ts_eff;

	//spin accumulation torque efficiency at interfaces (unitless, varies from 0 : no torque, up to 1 : full torque)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> tsi_eff;

	//spin pumping efficiency (unitless, varies from 0 : no spin pumping, up to 1 : full strength)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> pump_eff;

	//charge pumping efficiency (unitless, varies from 0 : no charge pumping, up to 1 : full strength)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> cpump_eff;

	//topological Hall effect efficiency (unitless, varies from 0 : none, up to 1 : full strength)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> the_eff;

	//the mesh base temperature (K)
	cu_obj<cuBReal> base_temperature;

	//thermal conductivity (W/mK) - default for permalloy
	cu_obj<MatPCUDA<cuBReal, cuBReal>> thermCond;

	//mass density (kg/m^3) - default for permalloy
	cu_obj<MatPCUDA<cuBReal, cuBReal>> density;

	//specific heat capacity (J/kgK) - default for permalloy
	cu_obj<MatPCUDA<cuBReal, cuBReal>> shc;

	//electron specific heat capacity at room temperature used in many-temperature models (J/kgK); Note, if used you should assign a temperature dependence to it, e.g. linear with temperature for the free electron approximation; none assigned by default.
	cu_obj<MatPCUDA<cuBReal, cuBReal>> shc_e;

	//electron-lattice coupling constant (W/m^3K) used in two-temperature model.
	cu_obj<MatPCUDA<cuBReal, cuBReal>> G_e;

	//set temperature spatial variation coefficient (unitless) - used with temperature settings in a simulation schedule only, not with console command directly
	cu_obj<MatPCUDA<cuBReal, cuBReal>> cT;

	//Heat source stimulus in heat equation. Ideally used with a spatial variation. (W//m3)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> Q;

private:

public:

	Atom_MeshParamsCUDA(Atom_MeshParams *pameshParams);
	virtual ~Atom_MeshParamsCUDA();
};

#endif

