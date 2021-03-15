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

	//Gilbert damping (atomistic: intrinsic)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> alpha;

	//atomic moment (units of muB) - default for bcc Fe
	cu_obj<MatPCUDA<cuBReal, cuBReal>> mu_s;

	//Exchange constant (units of J) - default for bcc Fe
	cu_obj<MatPCUDA<cuBReal, cuBReal>> J;

	//DMI exchange constant : (units of J)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> D;

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

	//electrical conductivity (units S/m).
	//this is the value at RT for Ni80Fe20.
	cu_obj<MatPCUDA<cuBReal, cuBReal>> elecCond;

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

