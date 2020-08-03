#include "stdafx.h"
#include "MeshParams.h"

MeshParams::MeshParams(vector<PARAM_>& enabledParams)
{
	//see comment in MeshParamsBase
	set_meshparamsbase_implementation(this, MESHTYPE_MICROMAGNETIC);

	//store all simulation parameters in meshParams : allows easy handling in the console, as well as saving and loading. In the computational routines use the parmeters directly, not through simParams
	for (int idx = 0; idx < (int)enabledParams.size(); idx++) {

		//List all parameters here, but only add them if configured to
		switch (enabledParams[idx]) {

		case PARAM_GREL:
			meshParams.push_back("grel", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_GREL);
			break;

		case PARAM_GREL_AFM:
			meshParams.push_back("grel_AFM", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_GREL_AFM);
			break;

		case PARAM_GDAMPING:
			meshParams.push_back("damping", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_GDAMPING);
			break;

		case PARAM_GDAMPING_AFM:
			meshParams.push_back("damping_AFM", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_GDAMPING_AFM);
			break;

		case PARAM_MS:
			meshParams.push_back("Ms", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "A/m"), PARAM_MS);
			break;

		case PARAM_MS_AFM:
			meshParams.push_back("Ms_AFM", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "A/m"), PARAM_MS_AFM);
			break;

		case PARAM_DEMAGXY:
			meshParams.push_back("Nxy", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_DEMAGXY);
			break;

		case PARAM_A:
			meshParams.push_back("A", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m"), PARAM_A);
			break;

		case PARAM_A_AFM:
			meshParams.push_back("A_AFM", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m"), PARAM_A_AFM);
			break;

		case PARAM_A_AFH:
			meshParams.push_back("Ah", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m3"), PARAM_A_AFH);
			break;

		case PARAM_AFTAU:
			//add it to mesh but mark it as hidden
			meshParams.push_back("tau_ii", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "", true), PARAM_AFTAU);
			break;

		case PARAM_AFTAUCROSS:
			//add it to mesh but mark it as hidden
			meshParams.push_back("tau_ij", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "", true), PARAM_AFTAUCROSS);
			break;

		case PARAM_A_AFNH:
			meshParams.push_back("Anh", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m"), PARAM_A_AFNH);
			break;

		case PARAM_D:
			meshParams.push_back("D", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m2"), PARAM_D);
			break;

		case PARAM_D_AFM:
			meshParams.push_back("D_AFM", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m2"), PARAM_D_AFM);
			break;

		case PARAM_J1:
			meshParams.push_back("J1", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m2"), PARAM_J1);
			break;

		case PARAM_J2:
			meshParams.push_back("J2", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m2"), PARAM_J2);
			break;

		case PARAM_NETADIA:
			meshParams.push_back("neta", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/Am"), PARAM_NETADIA);
			break;

		case PARAM_K1:
			meshParams.push_back("K1", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m3"), PARAM_K1);
			break;

		case PARAM_K2:
			meshParams.push_back("K2", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m3"), PARAM_K2);
			break;

		case PARAM_K1_AFM:
			meshParams.push_back("K1_AFM", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m3"), PARAM_K1_AFM);
			break;

		case PARAM_K2_AFM:
			meshParams.push_back("K2_AFM", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "J/m3"), PARAM_K2_AFM);
			break;

		case PARAM_TC:
			//add it to mesh but mark it as hidden
			meshParams.push_back("Tc", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "K", true), PARAM_TC);
			break;

		case PARAM_MUB:
			//add it to mesh but mark it as hidden
			meshParams.push_back("muB", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "uB", true), PARAM_MUB);
			break;

		case PARAM_MUB_AFM:
			//add it to mesh but mark it as hidden
			meshParams.push_back("muB_AFM", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "uB", true), PARAM_MUB_AFM);
			break;

		case PARAM_EA1:
			meshParams.push_back("ea1", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_EA1);
			break;

		case PARAM_EA2:
			meshParams.push_back("ea2", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_EA2);
			break;

		case PARAM_SUSREL:
			meshParams.push_back("susrel", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "As2/kg"), PARAM_SUSREL);
			break;

		case PARAM_SUSREL_AFM:
			meshParams.push_back("susrel_AFM", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "As2/kg"), PARAM_SUSREL_AFM);
			break;

		case PARAM_SUSPREL:
			meshParams.push_back("susprel", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "As2/kg"), PARAM_SUSPREL);
			break;

		case PARAM_HA:
			meshParams.push_back("cHa", MeshParamDescriptor(PARAMTYPE_MAGNETIC), PARAM_HA);
			break;

		case PARAM_HMO:
			meshParams.push_back("Hmo", MeshParamDescriptor(PARAMTYPE_MAGNETIC, "A/m"), PARAM_HMO);
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

		case PARAM_ISHA:
			meshParams.push_back("iSHA", MeshParamDescriptor(PARAMTYPE_ELECTRIC), PARAM_ISHA);
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

		case PARAM_MECOEFF:
			meshParams.push_back("MEc", MeshParamDescriptor(PARAMTYPE_MECHANICAL, "J/m3"), PARAM_MECOEFF);
			break;

		case PARAM_YOUNGSMOD:
			meshParams.push_back("Ym", MeshParamDescriptor(PARAMTYPE_MECHANICAL, "Pa"), PARAM_YOUNGSMOD);
			break;

		case PARAM_POISSONRATIO:
			meshParams.push_back("Pr", MeshParamDescriptor(PARAMTYPE_MECHANICAL), PARAM_POISSONRATIO);
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

	/// Special Functions

	//resolution of 10000 means e.g. for Tc = 1000 the Curie-Weiss function will be available with a resolution of 0.1 K
	pCurieWeiss = shared_ptr<Funcs_Special>(new Funcs_Special(EqComp::FUNC_CURIEWEISS, 0.0, 10000));
	pLongRelSus = shared_ptr<Funcs_Special>(new Funcs_Special(EqComp::FUNC_LONGRELSUS, 0.0, 10000));
	pLongRelSus->Initialize_LongitudinalRelSusceptibility(pCurieWeiss->get_data(), atomic_moment, T_Curie_material);

	pCurieWeiss1 = shared_ptr<Funcs_Special>(new Funcs_Special(EqComp::FUNC_CURIEWEISS1, 0.0, 10000));
	pCurieWeiss2 = shared_ptr<Funcs_Special>(new Funcs_Special(EqComp::FUNC_CURIEWEISS2, 0.0, 10000));
	pLongRelSus1 = shared_ptr<Funcs_Special>(new Funcs_Special(EqComp::FUNC_LONGRELSUS1, 0.0, 10000));
	pLongRelSus2 = shared_ptr<Funcs_Special>(new Funcs_Special(EqComp::FUNC_LONGRELSUS2, 0.0, 10000));
	pLongRelSus1->Initialize_LongitudinalRelSusceptibility1(pCurieWeiss1->get_data(), pCurieWeiss2->get_data(), tau_ii, tau_ij, atomic_moment_AFM, T_Curie_material);
	pLongRelSus2->Initialize_LongitudinalRelSusceptibility2(pCurieWeiss1->get_data(), pCurieWeiss2->get_data(), tau_ii, tau_ij, atomic_moment_AFM, T_Curie_material);

	pAlpha1 = shared_ptr<Funcs_Special>(new Funcs_Special(EqComp::FUNC_ALPHA1, 0.0, 10000));
	pAlpha2 = shared_ptr<Funcs_Special>(new Funcs_Special(EqComp::FUNC_ALPHA2, 0.0, 10000));

	//make sure special functions are set by default for all material parameters text equations
	set_special_functions();
}

//-------------------------Parameter control

//set pre-calculated Funcs_Special objects in material parameters
void MeshParams::set_special_functions(PARAM_ paramID)
{
	auto set_param_special_functions = [&](PARAM_ update_paramID) {

		auto code = [&](auto& MatP_object) -> void {

			MatP_object.set_t_scaling_special_functions(pCurieWeiss, pLongRelSus, pCurieWeiss1, pCurieWeiss2, pLongRelSus1, pLongRelSus2, pAlpha1, pAlpha2);
		};

		run_on_param<void>(update_paramID, code);
	};

	if (paramID == PARAM_ALL) {

		for (int index = 0; index < meshParams.size(); index++) {

			set_param_special_functions((PARAM_)meshParams.get_ID_from_index(index));
		}
	}
	else set_param_special_functions(paramID);
}

//copy all parameters from another Mesh
void MeshParams::copy_parameters(MeshParams& copy_this)
{
	auto copy_param = [&](PARAM_ paramID) {
		
		auto copy_from = [&](auto& MatP_object, MeshParams& copy_this) {

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
	T_Curie = copy_this.T_Curie;
}

//-------------------------Setters/Updaters : text equations

//set the mesh parameter equation with given user constants
void MeshParams::set_meshparam_t_equation(PARAM_ paramID, string& equationText, vector_key<double>& userConstants)
{
	auto code = [](auto& MatP_object, string& equationText, vector_key<double>& userConstants, double T_Curie, double base_temperature) -> void {

		MatP_object.set_t_scaling_equation(equationText, userConstants, T_Curie, base_temperature);
	};

	if (paramID == PARAM_ALL) {

		for (int index = 0; index < meshParams.size(); index++) {

			run_on_param<void>((PARAM_)meshParams.get_ID_from_index(index), code, equationText, userConstants, T_Curie, base_temperature);
		}
	}
	else run_on_param<void>(paramID, code, equationText, userConstants, T_Curie, base_temperature);
}

//update text equations for mesh parameters with user constants, mesh dimensions, Curie temperature, base temperature
bool MeshParams::update_meshparam_equations(PARAM_ paramID, vector_key<double>& userConstants, DBL3 meshDimensions)
{
	auto code = [](auto& MatP_object, vector_key<double>& userConstants, DBL3 meshDimensions, double T_Curie, double base_temperature) -> bool {

		return MatP_object.update_equations(userConstants, meshDimensions, T_Curie, base_temperature);
	};

	return run_on_param<bool>(paramID, code, userConstants, meshDimensions, T_Curie, base_temperature);
}