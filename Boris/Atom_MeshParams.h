#pragma once

#include "BorisLib.h"

#include "MaterialParameter.h"
#include "ErrorHandler.h"
#include "Boris_Enums_Defs.h"

#include "MeshParamsBase.h"

#if COMPILECUDA == 1
#include "Atom_MeshParamsCUDA.h"
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////
//
//	The Atom_Mesh material parameters
//

class Atom_MeshParams :
	virtual public MeshParamsBase			//need virtual inheritance : see MeshParamsBase.h comments
{

#if COMPILECUDA == 1
	friend Atom_MeshParamsCUDA;
#endif

	friend MeshParamsBase;					//so we can access run_on_param from it

private:

protected:

public:

	//-------------------------------- LIST ALL MESH PARAMETERS HERE

	//-----------SIMPLE CUBIC

	//Relative electron gyromagnetic ratio
	MatP<double, double> grel = 1.0;

	//Gilbert damping (atomistic: intrinsic)
	MatP<double, double> alpha = 0.02;

	//atomic moment (units of muB) - default for bcc Fe
	MatP<double, double> mu_s = 2.22;

	//Exchange constant : (units of J) - default for bcc Fe
	MatP<double, double> J = 7.05e-21;

	//DMI exchange constant : (units of J)
	MatP<double, double> D = 5e-23;

	//Interfacial DMI symmetry axis direction, used by vector interfacial DMI module
	MatP<DBL3, DBL3> D_dir = DBL3(0, 0, 1);

	//Surface exchange coupling, used by the surfexchange module to couple two spins on different meshes at the surface (units of J)
	MatP<double, double> Js = -1e-21;
	//Secondary surface exchange coupling constant, used for coupling atomistic meshes to micromagnetic 2-sublattice meshes.
	MatP<double, double> Js2 = -1e-21;

	//Magneto-crystalline anisotropy constants (J) and easy axes directions. For uniaxial anisotropy only ea1 is needed.
	MatP<double, double> K1 = 5.65e-25;
	MatP<double, double> K2 = 0.0;
	MatP<double, double> K3 = 0.0;

	//Magneto-crystalline anisotropy easy axes directions
	MatP<DBL3, DBL3> mcanis_ea1 = DBL3(1, 0, 0);
	MatP<DBL3, DBL3> mcanis_ea2 = DBL3(0, 1, 0);
	MatP<DBL3, DBL3> mcanis_ea3 = DBL3(0, 0, 1);

	//tensorial anisotropy. each term is a contribution to the anisotropy energy density as d*a^n1 b^n2 c^n3. Here a = m.mcanis_ea1, b = m.mcanis_ea2, c = m.mcanis_ea3.
	//For 2nd order we aditionally multiply by K1, 4th order K2, 6th order K3. Any other orders d coefficient contains anisotropy energy density.
	//each DBL4 stores (d, n1, n2, n3), where d != 0, n1, n2, n3 >= 0, n1+n2+n3>0. Odd order terms allowed.
	std::vector<DBL4> Kt;

	//-----------BCC (2 per unit cell)

	//-----------FCC (4 per unit cell)

	//-----------HCP (4 per effective unit cell)

	//-----------Others

	//in-plane demagnetizing factors (used for Atom_Demag_N module)
	MatP<DBL2, double> Nxy = DBL2(0);

	//applied field spatial variation coefficient (unitless)
	MatP<double, double> cHA = 1.0;

	//Magneto-Optical field strength (A/m)
	MatP<double, double> cHmo = 0.0;

	//Stochasticity efficiency parameter
	MatP<double, double> s_eff = 1.0;

	//electrical conductivity (units S/m).
	//this is the value at RT for Ni80Fe20.
	MatP<double, double> elecCond = 7e6;

	//TMR RA products for parallel and antiparallel states (Ohms m^2)
	MatP<double, double> RAtmr_p = 0.075e-12;
	MatP<double, double> RAtmr_ap = 0.225e-12;

	//anisotropic magnetoresistance as a percentage (of base resistance)
	MatP<double, double> amrPercentage = 0.0;

	//spin current polarization (also the charge current spin polarization for the spin current solver) and non-adiabaticity (for Zhang-Li STT). (unitless)
	MatP<double, double> P = 0.4;
	MatP<double, double> beta = 0.04;

	//parameters for spin current solver

	//electron diffusion constant (m^2/s)
	MatP<double, double> De = 1e-2;

	//electron carrier density (1/m^3)
	MatP<double, double> n_density = 1.8e29;

	//diffusion spin polarization (unitless)
	MatP<double, double> betaD = 0.5;

	//spin Hall angle (unitless)
	MatP<double, double> SHA = 0.1;

	//field-like spin torque coefficient (unitless)
	MatP<double, double> flSOT = 0.0;

	//Slonczewski macrospin torques q+, q- parameters as in PRB 72, 014446 (2005) (unitless)
	MatP<DBL2, double> STq = DBL2(1.0, 0.0);

	//Slonczewski macrospin torques A, B parameters as in PRB 72, 014446 (2005) (unitless)
	MatP<DBL2, double> STa = DBL2(0.6, 0.4);

	//Slonczewski macrospin torques spin polarization unit vector as in PRB 72, 014446 (2005) (unitless); or SOT symmetry axis direction (e.g. z direction for HM/FM bilayer).
	MatP<DBL3, DBL3> STp = DBL3(0, 0, 1);

	//spin-flip length (m)
	MatP<double, double> l_sf = 10e-9;

	//spin exchange rotation length (m)
	MatP<double, double> l_ex = 2e-9;

	//spin dephasing length (m)
	MatP<double, double> l_ph = 4e-9;

	//interface spin-dependent conductivity (spin-up and spin-down) (S/m^2)
	MatP<DBL2, double> Gi = DBL2(1e15, 1e14);

	//interface spin-mixing conductivity (real and imaginary parts) (S/m^2)
	MatP<DBL2, double> Gmix = DBL2(1e15, 1e14);

	//spin accumulation torque efficiency in the bulk (unitless, varies from 0 : no torque, up to 1 : full torque)
	MatP<double, double> ts_eff = 1;

	//spin accumulation torque efficiency at interfaces (unitless, varies from 0 : no torque, up to 1 : full torque)
	MatP<double, double> tsi_eff = 1;

	//spin pumping efficiency (unitless, varies from 0 : no spin pumping, up to 1 : full strength)
	//disabled by default
	MatP<double, double> pump_eff = 0;

	//charge pumping efficiency (unitless, varies from 0 : no charge pumping, up to 1 : full strength)
	//disabled by default
	MatP<double, double> cpump_eff = 0;

	//topological Hall effect efficiency (unitless, varies from 0 : none, up to 1 : full strength)
	//disabled by default
	MatP<double, double> the_eff = 0;

	//thermal conductivity (W/mK) - default for permalloy
	MatP<double, double> thermCond = 46.4;

	//mass density (kg/m^3) - default for permalloy
	MatP<double, double> density = 8740;

	//specific heat capacity (J/kgK) - default for permalloy
	MatP<double, double> shc = 430;

	//electron specific heat capacity at room temperature used in many-temperature models (J/kgK); Note, if used you should assign a temperature dependence to it, e.g. linear with temperature for the free electron approximation; none assigned by default.
	MatP<double, double> shc_e = 40;

	//electron-lattice coupling constant (W/m^3K) used in two-temperature model.
	MatP<double, double> G_e = 1e18;

private:

	//-------------------------Parameter control

	//run set code (run_this) using set parameters (run_this_args) on a MatP object selected through the paramID selector (majorID in meshParams identifying the required material parameter).
	template <typename RType, typename Lambda, typename ... PType>
	RType run_on_param(PARAM_ paramID, Lambda& run_this, PType& ... run_this_args);

protected:

public:

	//------------------------CTOR/DTOR

	//inherited by Mesh implementations
	Atom_MeshParams(std::vector<PARAM_>& enabledParams);
	
	virtual ~Atom_MeshParams() {}

	//-------------------------Setters

	//copy all parameters from another Mesh
	void copy_parameters(Atom_MeshParams& copy_this);

	//-------------------------Getters

	std::string get_tensorial_anisotropy_string(void);

	//-------------------------Setters/Updaters : text equations

	//set the mesh parameter temperature equation with given user constants
	void set_meshparam_t_equation(PARAM_ paramID, std::string& equationText, vector_key<double>& userConstants);

	//update text equations for mesh parameters with user constants, mesh dimensions, base temperature
	bool update_meshparam_equations(PARAM_ paramID, vector_key<double>& userConstants, DBL3 meshDimensions);
};

//-------------------------Parameter control

template <typename RType, typename Lambda, typename ... PType>
RType Atom_MeshParams::run_on_param(PARAM_ paramID, Lambda& run_this, PType& ... run_this_args)
{
	switch (paramID) {
		
	case PARAM_GREL:
		return run_this(grel, run_this_args...);
		break;

	case PARAM_ATOM_SC_DAMPING:
		return run_this(alpha, run_this_args...);
		break;

	case PARAM_ATOM_SC_MUS:
		return run_this(mu_s, run_this_args...);
		break;

	case PARAM_ATOM_SC_J:
		return run_this(J, run_this_args...);
		break;

	case PARAM_ATOM_SC_D:
		return run_this(D, run_this_args...);
		break;

	case PARAM_DMI_DIR:
		return run_this(D_dir, run_this_args...);
		break;

	case PARAM_ATOM_JS:
		return run_this(Js, run_this_args...);
		break;

	case PARAM_ATOM_JS2:
		return run_this(Js2, run_this_args...);
		break;

	case PARAM_ATOM_SC_K1:
		return run_this(K1, run_this_args...);
		break;

	case PARAM_ATOM_SC_K2:
		return run_this(K2, run_this_args...);
		break;

	case PARAM_ATOM_SC_K3:
		return run_this(K3, run_this_args...);
		break;

	case PARAM_ATOM_EA1:
		return run_this(mcanis_ea1, run_this_args...);
		break;

	case PARAM_ATOM_EA2:
		return run_this(mcanis_ea2, run_this_args...);
		break;

	case PARAM_ATOM_EA3:
		return run_this(mcanis_ea3, run_this_args...);
		break;

	case PARAM_DEMAGXY:
		return run_this(Nxy, run_this_args...);
		break;

	case PARAM_HA:
		return run_this(cHA, run_this_args...);
		break;

	case PARAM_HMO:
		return run_this(cHmo, run_this_args...);
		break;

	case PARAM_S_EFF:
		return run_this(s_eff, run_this_args...);
		break;

	case PARAM_ELC:
		return run_this(elecCond, run_this_args...);
		break;

	case PARAM_RATMR_P:
		return run_this(RAtmr_p, run_this_args...);
		break;

	case PARAM_RATMR_AP:
		return run_this(RAtmr_ap, run_this_args...);
		break;
		
	case PARAM_AMR:
		return run_this(amrPercentage, run_this_args...);
		break;

	case PARAM_P:
		return run_this(P, run_this_args...);
		break;

	case PARAM_BETA:
		return run_this(beta, run_this_args...);
		break;

	case PARAM_DE:
		return run_this(De, run_this_args...);
		break;

	case PARAM_BETAD:
		return run_this(betaD, run_this_args...);
		break;

	case PARAM_SHA:
		return run_this(SHA, run_this_args...);
		break;

	case PARAM_FLSOT:
		return run_this(flSOT, run_this_args...);
		break;

	case PARAM_STQ:
		return run_this(STq, run_this_args...);
		break;

	case PARAM_STA:
		return run_this(STa, run_this_args...);
		break;

	case PARAM_STP:
		return run_this(STp, run_this_args...);
		break;

	case PARAM_LSF:
		return run_this(l_sf, run_this_args...);
		break;

	case PARAM_LEX:
		return run_this(l_ex, run_this_args...);
		break;

	case PARAM_LPH:
		return run_this(l_ph, run_this_args...);
		break;

	case PARAM_GI:
		return run_this(Gi, run_this_args...);
		break;

	case PARAM_GMIX:
		return run_this(Gmix, run_this_args...);
		break;

	case PARAM_TSEFF:
		return run_this(ts_eff, run_this_args...);
		break;

	case PARAM_TSIEFF:
		return run_this(tsi_eff, run_this_args...);
		break;

	case PARAM_PUMPEFF:
		return run_this(pump_eff, run_this_args...);
		break;

	case PARAM_CPUMP_EFF:
		return run_this(cpump_eff, run_this_args...);
		break;

	case PARAM_THE_EFF:
		return run_this(the_eff, run_this_args...);
		break;

	case PARAM_NDENSITY:
		return run_this(n_density, run_this_args...);
		break;
		
	case PARAM_THERMCOND:
		return run_this(thermCond, run_this_args...);
		break;

	case PARAM_DENSITY:
		return run_this(density, run_this_args...);
		break;

	case PARAM_SHC:
		return run_this(shc, run_this_args...);
		break;

	case PARAM_SHC_E:
		return run_this(shc_e, run_this_args...);
		break;

	case PARAM_G_E:
		return run_this(G_e, run_this_args...);
		break;

	case PARAM_T:
		return run_this(cT, run_this_args...);
		break;

	case PARAM_Q:
		return run_this(Q, run_this_args...);
		break;

	default:
		//this is needed to stop the "not all control paths return a value" error, but should never get here
		return run_this(alpha, run_this_args...);
		break;
	}
}