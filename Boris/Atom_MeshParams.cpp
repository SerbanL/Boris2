#include "stdafx.h"
#include "Atom_MeshParams.h"

Atom_MeshParams::Atom_MeshParams(std::vector<PARAM_>& enabledParams)
{
	//see comment in MeshParamsBase
	set_meshparamsbase_implementation(this, MESHTYPE_ATOMISTIC);

	//store all simulation parameters in meshParams : allows easy handling in the console, as well as saving and loading. In the computational routines use the parmeters directly, not through simParams
	for (int idx = 0; idx < (int)enabledParams.size(); idx++) {

		//List all parameters here, but only add them if configured to
		switch (enabledParams[idx]) {
			
		case PARAM_GREL:
			meshParams.push_back("grel", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_GREL);
			break;

		case PARAM_ATOM_SC_DAMPING:
			meshParams.push_back("damping", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_ATOM_SC_DAMPING);
			break;

		case PARAM_ATOM_SC_MUS:
			meshParams.push_back("mu_s", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "muB"), PARAM_ATOM_SC_MUS);
			break;

		case PARAM_ATOM_SC_J:
			meshParams.push_back("J", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J"), PARAM_ATOM_SC_J);
			break;

		case PARAM_ATOM_SC_D:
			meshParams.push_back("D", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J"), PARAM_ATOM_SC_D);
			break;

		case PARAM_DMI_DIR:
			meshParams.push_back("Ddir", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_DMI_DIR);
			break;

		case PARAM_ATOM_JS:
			meshParams.push_back("Js", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J"), PARAM_ATOM_JS);
			break;

		case PARAM_ATOM_SC_K1:
			meshParams.push_back("K1", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J"), PARAM_ATOM_SC_K1);
			break;

		case PARAM_ATOM_SC_K2:
			meshParams.push_back("K2", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J"), PARAM_ATOM_SC_K2);
			break;

		case PARAM_ATOM_SC_K3:
			meshParams.push_back("K3", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J"), PARAM_ATOM_SC_K3);
			break;

		case PARAM_ATOM_EA1:
			meshParams.push_back("ea1", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_ATOM_EA1);
			break;

		case PARAM_ATOM_EA2:
			meshParams.push_back("ea2", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_ATOM_EA2);
			break;

		case PARAM_ATOM_EA3:
			meshParams.push_back("ea3", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_ATOM_EA3);
			break;

		case PARAM_DEMAGXY:
			meshParams.push_back("Nxy", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_DEMAGXY);
			break;

		case PARAM_HA:
			meshParams.push_back("cHa", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_HA);
			break;

		case PARAM_HMO:
			meshParams.push_back("cHmo", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_HMO);
			break;

		case PARAM_T:
			meshParams.push_back("cT", MeshParamDescriptor(PARAMTYPE_THERMAL), PARAM_T);
			break;

		case PARAM_Q:
			meshParams.push_back("Q", MeshParamDescriptor(PARAMTYPE_THERMAL, "W/m3"), PARAM_Q);
			break;

		case PARAM_ELC:
			meshParams.push_back("elC", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "S/m"), PARAM_ELC);
			break;
			
		case PARAM_AMR:
			meshParams.push_back("amr", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "%"), PARAM_AMR);
			break;

		case PARAM_P:
			meshParams.push_back("P", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_P);
			break;

		case PARAM_BETA:
			meshParams.push_back("beta", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_BETA);
			break;

		case PARAM_DE:
			meshParams.push_back("De", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "m2/s"), PARAM_DE);
			break;

		case PARAM_BETAD:
			meshParams.push_back("betaD", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_BETAD);
			break;

		case PARAM_SHA:
			meshParams.push_back("SHA", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_SHA);
			break;

		case PARAM_FLSOT:
			meshParams.push_back("flST", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_FLSOT);
			break;

		case PARAM_STQ:
			meshParams.push_back("STq", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_STQ);
			break;

		case PARAM_STA:
			meshParams.push_back("STa", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_STA);
			break;

		case PARAM_STP:
			meshParams.push_back("STp", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_STP);
			break;

		case PARAM_LSF:
			meshParams.push_back("l_sf", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "m"), PARAM_LSF);
			break;

		case PARAM_LEX:
			meshParams.push_back("l_J", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "m"), PARAM_LEX);
			break;

		case PARAM_LPH:
			meshParams.push_back("l_phi", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "m"), PARAM_LPH);
			break;

		case PARAM_GI:
			meshParams.push_back("Gi", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "S/m2"), PARAM_GI);
			break;

		case PARAM_GMIX:
			meshParams.push_back("Gmix", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "S/m2"), PARAM_GMIX);
			break;

		case PARAM_TSEFF:
			meshParams.push_back("ts_eff", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_TSEFF);
			break;

		case PARAM_TSIEFF:
			meshParams.push_back("tsi_eff", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_TSIEFF);
			break;

		case PARAM_PUMPEFF:
			meshParams.push_back("pump_eff", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_PUMPEFF);
			break;

		case PARAM_CPUMP_EFF:
			meshParams.push_back("cpump_eff", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_CPUMP_EFF);
			break;

		case PARAM_THE_EFF:
			meshParams.push_back("the_eff", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_THE_EFF);
			break;

		case PARAM_NDENSITY:
			meshParams.push_back("n", MeshParamDescriptor(PARAMTYPE_ELECTRIC, "m-3"), PARAM_NDENSITY);
			break;
			
		case PARAM_THERMCOND:
			meshParams.push_back("thermK", MeshParamDescriptor(PARAMTYPE_THERMAL, "W/mK"), PARAM_THERMCOND);
			break;

		case PARAM_DENSITY:
			meshParams.push_back("density", MeshParamDescriptor(PARAMTYPE_MECHANICAL, "kg/m3"), PARAM_DENSITY);
			break;

		case PARAM_SHC:
			meshParams.push_back("shc", MeshParamDescriptor(PARAMTYPE_THERMAL, "J/kgK"), PARAM_SHC);
			break;

		case PARAM_SHC_E:
			meshParams.push_back("shc_e", MeshParamDescriptor(PARAMTYPE_THERMAL, "J/kgK"), PARAM_SHC_E);
			break;

		case PARAM_G_E:
			meshParams.push_back("G_e", MeshParamDescriptor(PARAMTYPE_THERMAL, "W/m3K"), PARAM_G_E);
			break;
		}
	}
}

//-------------------------Getters

std::string Atom_MeshParams::get_tensorial_anisotropy_string(void)
{
	std::string Ktstring;

	for (int idx = 0; idx < Kt.size(); idx++) {

		Ktstring += ToString(Kt[idx].i) + "x" + ToString(Kt[idx].j) + "y" + ToString(Kt[idx].k) + "z" + ToString(Kt[idx].l);
		if (idx != Kt.size() - 1) Ktstring += " ";
	}

	return Ktstring;
}

//-------------------------Parameter control

//copy all parameters from another Mesh
void Atom_MeshParams::copy_parameters(Atom_MeshParams& copy_this)
{
	auto copy_param = [&](PARAM_ paramID) {
		
		auto copy_from = [&](auto& MatP_object, Atom_MeshParams& copy_this) {

			auto get_meshparam = [&](decltype(MatP_object)& MatP_object2) -> decltype(MatP_object)& {

				return MatP_object2;
			};

			MatP_object = copy_this.run_on_param<decltype(MatP_object)>(paramID, get_meshparam);
		};

		run_on_param<void>(paramID, copy_from, copy_this);
	};

	//copy all MatP parameters
	for (int index = 0; index < meshParams.size(); index++) {

		PARAM_ paramID = (PARAM_)meshParams.get_ID_from_index(index);

		copy_param(paramID);
	}

	//now copy special values which are not MatP
	base_temperature = copy_this.base_temperature;
}

//-------------------------Setters/Updaters : text equations

//set the mesh parameter equation with given user constants
void Atom_MeshParams::set_meshparam_t_equation(PARAM_ paramID, std::string& equationText, vector_key<double>& userConstants)
{
	auto code = [](auto& MatP_object, std::string& equationText, vector_key<double>& userConstants, double base_temperature) -> void {

		MatP_object.set_t_scaling_equation(equationText, userConstants, 0.0, base_temperature);
	};

	if (paramID == PARAM_ALL) {

		for (int index = 0; index < meshParams.size(); index++) {

			run_on_param<void>((PARAM_)meshParams.get_ID_from_index(index), code, equationText, userConstants, base_temperature);
		}
	}
	else run_on_param<void>(paramID, code, equationText, userConstants, base_temperature);
}

//update text equations for mesh parameters with user constants, mesh dimensions, base temperature
bool Atom_MeshParams::update_meshparam_equations(PARAM_ paramID, vector_key<double>& userConstants, DBL3 meshDimensions)
{
	auto code = [](auto& MatP_object, vector_key<double>& userConstants, DBL3 meshDimensions, double base_temperature) -> bool {

		return MatP_object.update_equations(userConstants, meshDimensions, 0.0, base_temperature);
	};

	return run_on_param<bool>(paramID, code, userConstants, meshDimensions, base_temperature);
}