#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "MaterialParameterCUDA.h"
#include "ParametersDefs.h"

class MeshParams;

class MeshParamsCUDA {

private:

	MeshParams *pmeshParams;

protected:

	//Special functions to be set in material parameters text equations when needed

	//resolution of 10000 means e.g. for Tc = 1000 the Curie-Weiss function will be available with a resolution of 0.1 K
	cu_obj<ManagedFuncs_Special_CUDA> CurieWeiss_CUDA;
	cu_obj<ManagedFuncs_Special_CUDA> LongRelSus_CUDA;

	cu_obj<ManagedFuncs_Special_CUDA> CurieWeiss1_CUDA;
	cu_obj<ManagedFuncs_Special_CUDA> CurieWeiss2_CUDA;
	cu_obj<ManagedFuncs_Special_CUDA> LongRelSus1_CUDA;
	cu_obj<ManagedFuncs_Special_CUDA> LongRelSus2_CUDA;

	cu_obj<ManagedFuncs_Special_CUDA> Alpha1_CUDA;
	cu_obj<ManagedFuncs_Special_CUDA> Alpha2_CUDA;

public:

	//A number of parameters have _AFM termination. These are used for antiferromagnetic meshes with 2-sublattice local approximation and are doubled-up for sub-lattices A, B

	//Relative electron gyromagnetic ratio
	cu_obj<MatPCUDA<cuBReal, cuBReal>> grel;
	cu_obj<MatPCUDA<cuReal2, cuBReal>> grel_AFM;

	//Gilbert damping
	cu_obj<MatPCUDA<cuBReal, cuBReal>> alpha;
	cu_obj<MatPCUDA<cuReal2, cuBReal>> alpha_AFM;

	//Saturation magnetization (A/m)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> Ms;
	cu_obj<MatPCUDA<cuReal2, cuBReal>> Ms_AFM;

	//in-plane demagnetizing factors (used for Demag_N module)
	cu_obj<MatPCUDA<cuReal2, cuBReal>> Nxy;

	//Exchange stiffness (J/m)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> A;
	cu_obj<MatPCUDA<cuReal2, cuBReal>> A_AFM;

	//Homogeneous AFM coupling between sub-lattices A and B, defined as A / a*a (J/m^3), where A is the homogeneous antiferromagnetic exchange stifness (negative), and a is the lattice constant.
	//e.g. a = 0.3nm, A = -1pJ/m gives Ah as -1e7 J/m^3 to order of magnitude.
	cu_obj<MatPCUDA<cuReal2, cuBReal>> Ah;

	//Nonhomogeneous AFM coupling between sub-lattices A and B (J/m)
	cu_obj<MatPCUDA<cuReal2, cuBReal>> Anh;

	//Dzyaloshinskii-Moriya exchange constant (J/m^2)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> D;
	cu_obj<MatPCUDA<cuReal2, cuBReal>> D_AFM;

	//Coupling between exchange integral and critical temperature (Neel or Curie temperature) for 2-sublattice model : intra-lattice term, 0.5 for ideal antiferromagnet
	//J = 3 * tau * kB * Tc
	cu_obj<MatPCUDA<cuReal2, cuBReal>> tau_ii;

	//Coupling between exchange integral and critical temperature (Neel or Curie temperature) for 2-sublattice model : inter-lattice, or cross-lattice term, 0.5 for ideal antiferromagnet.
	//J = 3 * tau * kB * Tc
	cu_obj<MatPCUDA<cuReal2, cuBReal>> tau_ij;

	//bilinear surface exchange coupling (J/m^2) : J1, bottom and top layer values
	//biquadratic surface exchange coupling (J/m^2) : J2, bottom and top layer values
	cu_obj<MatPCUDA<cuBReal, cuBReal>> J1;
	cu_obj<MatPCUDA<cuBReal, cuBReal>> J2;

	//surface exchange coupling per magnetization from a diamagnet (J/Am) - EXPERIMENTAL : Hse = neta * sus * Hext / mu0 * Ms * tF
	cu_obj<MatPCUDA<cuBReal, cuBReal>> neta_dia;

	//Magneto-crystalline anisotropy K1 and K2 constants (J/m^3) and easy axes directions. For uniaxial anisotropy only ea1 is needed, for cubic ea1 and ea2 should be orthogonal.
	cu_obj<MatPCUDA<cuBReal, cuBReal>> K1;
	cu_obj<MatPCUDA<cuBReal, cuBReal>> K2;
	cu_obj<MatPCUDA<cuBReal, cuBReal>> K3;
	cu_obj<MatPCUDA<cuReal3, cuReal3>> mcanis_ea1;
	cu_obj<MatPCUDA<cuReal3, cuReal3>> mcanis_ea2;
	cu_obj<MatPCUDA<cuReal3, cuReal3>> mcanis_ea3;

	//Anisotropy values for 2-sublattice model
	cu_obj<MatPCUDA<cuReal2, cuBReal>> K1_AFM;
	cu_obj<MatPCUDA<cuReal2, cuBReal>> K2_AFM;
	cu_obj<MatPCUDA<cuReal2, cuBReal>> K3_AFM;

	//tensorial anisotropy. each term is a contribution to the anisotropy energy density as d*a^n1 b^n2 c^n3. Here a = m.mcanis_ea1, b = m.mcanis_ea2, c = m.mcanis_ea3.
	//For 2nd order we aditionally multiply by K1, 4th order K2, 6th order K3. Any other orders d coefficient contains anisotropy energy density.
	//each DBL4 stores (d, n1, n2, n3), where d != 0, n1, n2, n3 >= 0, n1+n2+n3>0. Odd order terms allowed.
	cu_obj<cuVEC<cuReal4>> Kt, Kt2;

	//longitudinal (parallel) susceptibility relative to mu0*Ms0, i.e. divided by mu0*Ms0, Ms0 is the 0K Ms value - for use with LLB equation. Units As^2/kg
	cu_obj<MatPCUDA<cuBReal, cuBReal>> susrel;

	//longitudinal (parallel) susceptibility relative to mu0*Ms0, i.e. divided by mu0*Ms0, Ms0 is the 0K Ms value - for use with LLB equation 2-sublattice model. Units As^2/kg
	cu_obj<MatPCUDA<cuReal2, cuBReal>> susrel_AFM;

	//perpendicular (transverse) susceptibility relative to mu0*Ms0, i.e. divided by mu0*Ms0, Ms0 is the 0K Ms value - for use with LLB equation. Units As^2/kg
	cu_obj<MatPCUDA<cuBReal, cuBReal>> susprel;

	//applied field spatial variation coefficient (unitless)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> cHA;

	//Magneto-Optical field strength (A/m)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> cHmo;

	//electrical conductivity (units S/m).
	//this is the value at 0K for Ni80Fe20. Temperature dependence typically scaled by 1 / (1 + alpha*(T-T0)), where alpha = 0.003, T0 = 293K with sigma = 1.7e6 S/m and 293K.
	//Using scaling 1 / (1 + alpha0 * T) on the zero-temperature conductivity gives sigma0 = sigmaT0 / (1 - alpha*T0), alpha0 = alpha / (1 - alpha*T0), so alpha0 = 0.025.
	cu_obj<MatPCUDA<cuBReal, cuBReal>> elecCond;

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

	//"inverse" spin Hall angle (unitless) -> should normally be the same as SHA but e.g. can be set to zero to turn off the inverse SHE in the spin transport equation
	cu_obj<MatPCUDA<cuBReal, cuBReal>> iSHA;

	//field-like spin torque coefficient (unitless)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> flSOT;

	//Slonczewski macrospin torques q+, q- parameters as in PRB 72, 014446 (2005) (unitless)
	cu_obj<MatPCUDA<cuReal2, cuBReal>> STq;

	//Slonczewski macrospin torques A, B parameters as in PRB 72, 014446 (2005) (unitless)
	cu_obj<MatPCUDA<cuReal2, cuBReal>> STa;

	//Slonczewski macrospin torques spin polarization unit vector as in PRB 72, 014446 (2005) (unitless)
	cu_obj<MatPCUDA<cuReal3, cuBReal>> STp;

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

	//Curie temperture (K)
	cu_obj<cuBReal> T_Curie;

	//The atomic magnetic moment as a multiple of the Bohr magneton - default 1 ub for permalloy.
	cu_obj<MatPCUDA<cuBReal, cuBReal>> atomic_moment;

	//atomic moments for 2-sublattice model (again multiples of the Bohr magneton)
	cu_obj<MatPCUDA<cuReal2, cuBReal>> atomic_moment_AFM;

	//thermal conductivity (W/mK)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> thermCond;

	//mass density (kg/m^3)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> density;

	//Magneto-elastic coefficients (J/m^3)
	cu_obj<MatPCUDA<cuReal2, cuBReal>> MEc;

	//Young's modulus (Pa)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> Ym;

	//Poisson's ratio (unitless)
	cu_obj<MatPCUDA<cuBReal, cuBReal>> Pr;

	//specific heat capacity (J/kgK)
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

	//set pre-calculated Funcs_Special objects in material parameters
	void set_special_functions(PARAM_ paramID = PARAM_ALL);

public:

	MeshParamsCUDA(MeshParams *pmeshParams);
	virtual ~MeshParamsCUDA();

	void set_special_functions_data(void);
};

#endif
