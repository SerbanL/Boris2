#include "stdafx.h"
#include "SimSharedData.h"
#include "Mesh.h"
#include "Modules.h"
#include "MeshParams.h"
#include "SimulationData.h"

vector_key_lut<int> SimulationSharedData::formula_descriptor;
vector_key_lut<string> SimulationSharedData::vargenerator_descriptor;

vector_lut<string> SimulationSharedData::meshDisplay_unit;
vector_lut<string> SimulationSharedData::displayHandles;

vector_lut< vector<MOD_> > SimulationSharedData::modules_for_meshtype;

vector_lut< vector<MESHDISPLAY_> > SimulationSharedData::meshAllowedDisplay;

vector_lut< vector<PARAM_> > SimulationSharedData::params_for_meshtype;

exclusions<MOD_> SimulationSharedData::exclusiveModules;

exclusions<MOD_> SimulationSharedData::superMeshExclusiveModules;

exclusions<MOD_> SimulationSharedData::superMeshCompanionModules;

vector_lut<DatumConfig> SimulationSharedData::saveDataList;

bool SimulationSharedData::shape_change_individual = false;

bool SimulationSharedData::static_transport_solver = false;

bool SimulationSharedData::cudaEnabled = false;

size_t SimulationSharedData::gpuMemFree_MB = 0;
size_t SimulationSharedData::gpuMemTotal_MB = 0;
size_t SimulationSharedData::cpuMemFree_MB = 0;
size_t SimulationSharedData::cpuMemTotal_MB = 0;

SimulationSharedData::SimulationSharedData(bool called_from_Simulation)
{ 
	//only make the statically shared objects once when called from Simulation object - other objects only inherit these but do not make any changes
	if (called_from_Simulation) {

		//--------------

		formula_descriptor.push_back("none", 0, MATPFORM_NONE);
		formula_descriptor.push_back("linear", 1, MATPFORM_LINEAR);
		formula_descriptor.push_back("parabolic", 2, MATPFORM_PARABOLIC);
		formula_descriptor.push_back("ilinear", 1, MATPFORM_INVERSELINEAR);

		//--------------

		vargenerator_descriptor.push_back("custom", "0, 1; pngfile", MATPVAR_CUSTOM);
		//DBL2 range, int seed
		vargenerator_descriptor.push_back("random", "0.9, 1.1; 1", MATPVAR_RANDOM);
		//DBL2 range, double spacing, int seed
		vargenerator_descriptor.push_back("jagged", "0.9, 1.1; 20nm; 1", MATPVAR_JAGGED);
		//DBL2 range, DBL2 diameter_range, double spacing, int seed
		vargenerator_descriptor.push_back("defects", "0.9, 1.1; 20nm, 50nm; 100nm; 1", MATPVAR_DEFECTS);
		//DBL2 range, DBL2 length_range, DBL2 orientation_range, double spacing, int seed
		vargenerator_descriptor.push_back("faults", "0.9, 1.1; 50nm, 100nm; -20, 20; 100nm; 1", MATPVAR_FAULTS);
		//DBL2 range, double spacing, int seed
		vargenerator_descriptor.push_back("vor2D", "0.9, 1.1; 40nm; 1", MATPVAR_VORONOI2D);
		//DBL2 range, double spacing, int seed
		vargenerator_descriptor.push_back("vor3D", "0.9, 1.1; 40nm; 1", MATPVAR_VORONOI3D);
		//DBL2 range, double spacing, int seed
		vargenerator_descriptor.push_back("vorbnd2D", "0.9, 1.1; 40nm; 1", MATPVAR_VORONOIBND2D);
		//DBL2 range, double spacing, int seed
		vargenerator_descriptor.push_back("vorbnd3D", "0.9, 1.1; 40nm; 1", MATPVAR_VORONOIBND3D);
		//DBL2 theta (degrees), DBL2 phi (degrees), double spacing, int seed
		vargenerator_descriptor.push_back("vorrot2D", "70, 110; -45, 45; 40nm; 1", MATPVAR_VORONOIROT2D);
		//DBL2 theta (degrees), DBL2 phi (degrees), double spacing, int seed
		vargenerator_descriptor.push_back("vorrot3D", "70, 110; -45, 45; 40nm; 1", MATPVAR_VORONOIROT3D);

		//--------------

		meshDisplay_unit.push_back("", MESHDISPLAY_NONE);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_SM_DEMAG);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_SM_OERSTED);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_SM_STRAYH);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_MAGNETIZATION);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_MAGNETIZATION2);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_MAGNETIZATION12);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_EFFECTIVEFIELD);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_EFFECTIVEFIELD2);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_EFFECTIVEFIELD12);
		meshDisplay_unit.push_back("A/m3", MESHDISPLAY_CURRDENSITY);
		meshDisplay_unit.push_back("V", MESHDISPLAY_VOLTAGE);
		meshDisplay_unit.push_back("S/m", MESHDISPLAY_ELCOND);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_SACCUM);
		meshDisplay_unit.push_back("A/s", MESHDISPLAY_JSX);
		meshDisplay_unit.push_back("A/s", MESHDISPLAY_JSY);
		meshDisplay_unit.push_back("A/s", MESHDISPLAY_JSZ);
		meshDisplay_unit.push_back("A/ms", MESHDISPLAY_TS);
		meshDisplay_unit.push_back("A/ms", MESHDISPLAY_TSI);
		meshDisplay_unit.push_back("K", MESHDISPLAY_TEMPERATURE);
		meshDisplay_unit.push_back("", MESHDISPLAY_PARAMVAR);
		meshDisplay_unit.push_back("", MESHDISPLAY_ROUGHNESS);
		meshDisplay_unit.push_back("", MESHDISPLAY_CUSTOM_VEC);
		meshDisplay_unit.push_back("", MESHDISPLAY_CUSTOM_SCA);

		displayHandles.push_back("Nothing", MESHDISPLAY_NONE);
		displayHandles.push_back("Hdemag", MESHDISPLAY_SM_DEMAG);
		displayHandles.push_back("HOe", MESHDISPLAY_SM_OERSTED);
		displayHandles.push_back("Hstray", MESHDISPLAY_SM_STRAYH);
		displayHandles.push_back("M", MESHDISPLAY_MAGNETIZATION);
		displayHandles.push_back("M2", MESHDISPLAY_MAGNETIZATION2);
		displayHandles.push_back("M12", MESHDISPLAY_MAGNETIZATION12);
		displayHandles.push_back("Heff", MESHDISPLAY_EFFECTIVEFIELD);
		displayHandles.push_back("Heff2", MESHDISPLAY_EFFECTIVEFIELD2);
		displayHandles.push_back("Heff12", MESHDISPLAY_EFFECTIVEFIELD12);
		displayHandles.push_back("Jc", MESHDISPLAY_CURRDENSITY);
		displayHandles.push_back("V", MESHDISPLAY_VOLTAGE);
		displayHandles.push_back("elC", MESHDISPLAY_ELCOND);
		displayHandles.push_back("S", MESHDISPLAY_SACCUM);
		displayHandles.push_back("Jsx", MESHDISPLAY_JSX);
		displayHandles.push_back("Jsy", MESHDISPLAY_JSY);
		displayHandles.push_back("Jsz", MESHDISPLAY_JSZ);
		displayHandles.push_back("Ts", MESHDISPLAY_TS);
		displayHandles.push_back("Tsi", MESHDISPLAY_TSI);
		displayHandles.push_back("Temp", MESHDISPLAY_TEMPERATURE);
		displayHandles.push_back("ParamVar", MESHDISPLAY_PARAMVAR);
		displayHandles.push_back("Roughness", MESHDISPLAY_ROUGHNESS);
		displayHandles.push_back("Cust_V", MESHDISPLAY_CUSTOM_VEC);
		displayHandles.push_back("Cust_S", MESHDISPLAY_CUSTOM_SCA);

		//--------------

		//specify forbidden module combinations - each set is an exclusive modules set, i.e. only one of each can be active at any one time
		//also specify non-exclusive modules (i.e. entries in exclusiveModules with only one entry per set)
		exclusiveModules.storeset(MOD_DEMAG_N, MOD_DEMAG, MOD_SDEMAG_DEMAG);
		exclusiveModules.storeset(MOD_EXCHANGE6NGBR, MOD_DMEXCHANGE, MOD_IDMEXCHANGE);
		exclusiveModules.storeset(MOD_SURFEXCHANGE);
		exclusiveModules.storeset(MOD_ANIUNI, MOD_ANICUBI);
		exclusiveModules.storeset(MOD_ZEEMAN);
		exclusiveModules.storeset(MOD_TRANSPORT);
		exclusiveModules.storeset(MOD_HEAT);
		exclusiveModules.storeset(MOD_SOTFIELD);
		exclusiveModules.storeset(MOD_ROUGHNESS);

		//--------------

		//for some supermesh modules, specify a number of modules which run on individual meshes, which should not run if the supermesh version is active
		superMeshExclusiveModules.storeset(MODS_SDEMAG, MOD_DEMAG_N, MOD_DEMAG);

		//this is the opposite of above: if a module in a superMeshCompanionModules set is active, then all the other ones must be active too
		superMeshCompanionModules.storeset(MODS_STRANSPORT, MOD_TRANSPORT);
		superMeshCompanionModules.storeset(MODS_SHEAT, MOD_HEAT);

		//---------------

		//assign possible modules for each mesh type
		modules_for_meshtype.push_back(make_vector(MODS_SDEMAG, MODS_STRAYFIELD, MODS_STRANSPORT, MODS_OERSTED, MODS_SHEAT), MESH_SUPERMESH);

		modules_for_meshtype.push_back(make_vector(
			MOD_DEMAG_N, MOD_DEMAG, MOD_SDEMAG_DEMAG,
			MOD_EXCHANGE6NGBR, MOD_DMEXCHANGE, MOD_IDMEXCHANGE, MOD_SURFEXCHANGE, 
			MOD_ZEEMAN,
			MOD_ANIUNI, MOD_ANICUBI,
			MOD_TRANSPORT,
			MOD_HEAT, 
			MOD_SOTFIELD, MOD_ROUGHNESS), MESH_FERROMAGNETIC);

		modules_for_meshtype.push_back(make_vector(
			MOD_DEMAG_N, MOD_DEMAG, MOD_SDEMAG_DEMAG,
			MOD_EXCHANGE6NGBR, MOD_DMEXCHANGE, MOD_IDMEXCHANGE, 
			MOD_ZEEMAN, 
			MOD_ANIUNI, MOD_ANICUBI, 
			MOD_ROUGHNESS), MESH_ANTIFERROMAGNETIC);

		modules_for_meshtype.push_back(make_vector(MOD_TRANSPORT, MOD_HEAT), MESH_DIPOLE);

		modules_for_meshtype.push_back(make_vector(MOD_TRANSPORT, MOD_HEAT), MESH_METAL);

		modules_for_meshtype.push_back(make_vector(MOD_HEAT), MESH_INSULATOR);

		//----------------

		meshAllowedDisplay.push_back(make_vector(MESHDISPLAY_NONE, MESHDISPLAY_SM_DEMAG, MESHDISPLAY_SM_OERSTED, MESHDISPLAY_SM_STRAYH), MESH_SUPERMESH);

		meshAllowedDisplay.push_back(make_vector(
			MESHDISPLAY_NONE, MESHDISPLAY_MAGNETIZATION, MESHDISPLAY_EFFECTIVEFIELD, 
			MESHDISPLAY_CURRDENSITY, MESHDISPLAY_VOLTAGE, MESHDISPLAY_ELCOND, MESHDISPLAY_SACCUM, MESHDISPLAY_JSX, MESHDISPLAY_JSY, MESHDISPLAY_JSZ, MESHDISPLAY_TS, MESHDISPLAY_TSI,
			MESHDISPLAY_TEMPERATURE, MESHDISPLAY_PARAMVAR, MESHDISPLAY_ROUGHNESS, MESHDISPLAY_CUSTOM_VEC, MESHDISPLAY_CUSTOM_SCA), MESH_FERROMAGNETIC);

		meshAllowedDisplay.push_back(make_vector(
			MESHDISPLAY_NONE, MESHDISPLAY_MAGNETIZATION, MESHDISPLAY_MAGNETIZATION2, MESHDISPLAY_MAGNETIZATION12, MESHDISPLAY_EFFECTIVEFIELD, MESHDISPLAY_EFFECTIVEFIELD2, MESHDISPLAY_EFFECTIVEFIELD12,
			MESHDISPLAY_PARAMVAR, MESHDISPLAY_ROUGHNESS), MESH_ANTIFERROMAGNETIC);

		meshAllowedDisplay.push_back(make_vector(
			MESHDISPLAY_NONE, MESHDISPLAY_MAGNETIZATION, 
			MESHDISPLAY_CURRDENSITY, MESHDISPLAY_VOLTAGE, MESHDISPLAY_ELCOND, MESHDISPLAY_SACCUM, MESHDISPLAY_JSX, MESHDISPLAY_JSY, MESHDISPLAY_JSZ,
			MESHDISPLAY_TEMPERATURE, MESHDISPLAY_PARAMVAR), MESH_DIPOLE);

		meshAllowedDisplay.push_back(make_vector(
			MESHDISPLAY_NONE, 
			MESHDISPLAY_CURRDENSITY, MESHDISPLAY_VOLTAGE, MESHDISPLAY_ELCOND, MESHDISPLAY_SACCUM, MESHDISPLAY_JSX, MESHDISPLAY_JSY, MESHDISPLAY_JSZ,
			MESHDISPLAY_TEMPERATURE, MESHDISPLAY_PARAMVAR), MESH_METAL);

		meshAllowedDisplay.push_back(make_vector(MESHDISPLAY_NONE, MESHDISPLAY_TEMPERATURE, MESHDISPLAY_PARAMVAR), MESH_INSULATOR);

		//----------------

		//assign possible parameters for each mesh type
		params_for_meshtype.push_back(make_vector(PARAM_GREL, PARAM_GDAMPING, PARAM_MS, PARAM_DEMAGXY, PARAM_A, PARAM_D, PARAM_J1, PARAM_J2, PARAM_K1, PARAM_K2, PARAM_EA1, PARAM_EA2, PARAM_TC, PARAM_MUB, PARAM_SUSREL, PARAM_SUSPREL, PARAM_HA,
			PARAM_ELC, PARAM_AMR, PARAM_P, PARAM_BETA, PARAM_DE, PARAM_NDENSITY, PARAM_SHA, PARAM_FLSOT, PARAM_BETAD, PARAM_LSF, PARAM_LEX, PARAM_LPH, PARAM_GI, PARAM_GMIX, PARAM_PUMPEFF, PARAM_CPUMP_EFF, PARAM_THE_EFF, PARAM_TSEFF, PARAM_TSIEFF,
			PARAM_THERMCOND, PARAM_DENSITY, PARAM_SHC, PARAM_T, PARAM_Q), MESH_FERROMAGNETIC);

		params_for_meshtype.push_back(make_vector(PARAM_GREL_AFM, PARAM_GDAMPING_AFM, PARAM_MS_AFM, PARAM_DEMAGXY, PARAM_A_AFM, PARAM_A12, PARAM_D_AFM, PARAM_K1, PARAM_K2, PARAM_EA1, PARAM_EA2, PARAM_HA), MESH_ANTIFERROMAGNETIC);

		params_for_meshtype.push_back(make_vector(PARAM_MS, PARAM_TC,
			PARAM_ELC, PARAM_AMR, PARAM_P, PARAM_DE, PARAM_NDENSITY, PARAM_BETAD, PARAM_LSF, PARAM_LEX, PARAM_LPH, PARAM_GI, PARAM_GMIX,
			PARAM_THERMCOND, PARAM_DENSITY, PARAM_SHC, PARAM_T, PARAM_Q), MESH_DIPOLE);

		params_for_meshtype.push_back(make_vector(
			PARAM_ELC, PARAM_DE, PARAM_NDENSITY, PARAM_SHA, PARAM_ISHA, PARAM_LSF, PARAM_GI, PARAM_GMIX,
			PARAM_THERMCOND, PARAM_DENSITY, PARAM_SHC, PARAM_T, PARAM_Q), MESH_METAL);

		params_for_meshtype.push_back(make_vector(PARAM_THERMCOND, PARAM_DENSITY, PARAM_SHC), MESH_INSULATOR);

		//there's also an entry for MESH_SUPERMESH : this includes all possible material parameters which could be saved in a materials data base
		params_for_meshtype.push_back(make_vector(PARAM_GREL, PARAM_GDAMPING, PARAM_MS, PARAM_A, PARAM_D, PARAM_J1, PARAM_J2, PARAM_K1, PARAM_K2, PARAM_EA1, PARAM_EA2, PARAM_TC, PARAM_MUB,
			PARAM_ELC, PARAM_AMR, PARAM_P, PARAM_BETA, PARAM_DE, PARAM_BETAD, PARAM_SHA, PARAM_ISHA, PARAM_FLSOT, PARAM_LSF, PARAM_LEX, PARAM_LPH, PARAM_GI, PARAM_GMIX,
			PARAM_THERMCOND, PARAM_DENSITY, PARAM_SHC), MESH_SUPERMESH);
	}
}