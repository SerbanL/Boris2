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

	//electrical conductivity (units S/m).
	//this is the value at RT for Ni80Fe20.
	MatP<double, double> elecCond = 7e6;

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

	case PARAM_ELC:
		return run_this(elecCond, run_this_args...);
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