#include "stdafx.h"
#include "SimSharedData.h"
#include "Mesh.h"
#include "Modules.h"
#include "MeshParams.h"
#include "SimulationData.h"

vector_lut<std::string> SimulationSharedData::temperature_dependence_type;

vector_key_lut<std::string> SimulationSharedData::vargenerator_descriptor;

vector_lut<std::string> SimulationSharedData::meshDisplay_unit;
vector_lut<std::string> SimulationSharedData::displayHandles;

vector_lut< std::vector<MOD_> > SimulationSharedData::modules_for_meshtype;
vector_lut< std::vector<MOD_> > SimulationSharedData::displaymodules_for_meshtype;

vector_lut< std::vector<MESHDISPLAY_> > SimulationSharedData::meshAllowedDisplay;

vector_lut< std::vector<PARAM_> > SimulationSharedData::params_for_meshtype;

vector_lut<std::pair<bool, bool>> SimulationSharedData::params_enabled_props;

exclusions<MOD_> SimulationSharedData::exclusiveModules;

exclusions<MOD_> SimulationSharedData::superMeshExclusiveModules;

exclusions<MOD_> SimulationSharedData::superMeshCompanionModules;

vector_lut<DatumConfig> SimulationSharedData::saveDataList;

std::string SimulationSharedData::directory;

DBL2 SimulationSharedData::displayTransparency = DBL2(0.7, 1.0);
DBL2 SimulationSharedData::displayThresholds = DBL2();
int SimulationSharedData::displayThresholdTrigger = (int)VEC3REP_Z;

bool SimulationSharedData::shape_change_individual = false;

bool SimulationSharedData::static_transport_solver = false;

bool SimulationSharedData::disabled_transport_solver = false;

bool SimulationSharedData::cudaEnabled = false;

int SimulationSharedData::cudaDeviceSelect = 0;

INT2 SimulationSharedData::stage_step = INT2();

size_t SimulationSharedData::gpuMemFree_MB = 0;
size_t SimulationSharedData::gpuMemTotal_MB = 0;
size_t SimulationSharedData::cpuMemFree_MB = 0;
size_t SimulationSharedData::cpuMemTotal_MB = 0;

//collect available cuda device major versions here, indexed from 0
std::vector<std::pair<int, std::string>> SimulationSharedData::cudaDeviceVersions;

vector_key<double> SimulationSharedData::userConstants;

SimulationSharedData::SimulationSharedData(bool called_from_Simulation)
{ 
	//only make the statically shared objects once when called from Simulation object - other objects only inherit these but do not make any changes
	if (called_from_Simulation) {

		//--------------

		temperature_dependence_type.push_back("none", MATPTDEP_NONE);
		temperature_dependence_type.push_back("array", MATPTDEP_ARRAY);
		temperature_dependence_type.push_back("equation", MATPTDEP_EQUATION);

		//--------------

		vargenerator_descriptor.push_back("none", "", MATPVAR_NONE);
		vargenerator_descriptor.push_back("mask", "0, 1; pngfile", MATPVAR_MASK);
		vargenerator_descriptor.push_back("shape", "shape", MATPVAR_SHAPE);
		vargenerator_descriptor.push_back("ovf2", "ovf2file", MATPVAR_OVF2);
		vargenerator_descriptor.push_back("equation", "1", MATPVAR_EQUATION);
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
		//DBL2 x sides ratios, DBL2 y sides ratios, DBL2 z sides ratios, DBL3 inner mimimum, outer maximum, polynomial exponent
		vargenerator_descriptor.push_back("abl_pol", "0.25, 0.25; 0, 0; 0, 0; 0, 1, 1", MATPVAR_ABLPOL);
		//DBL2 x sides ratios, DBL2 y sides ratios, DBL2 z sides ratios, DBL3 inner mimimum, outer maximum, tanh sigma in nm
		vargenerator_descriptor.push_back("abl_tanh", "0.25, 0.25; 0, 0; 0, 0; 0, 1, 10.0", MATPVAR_ABLTANH);
		//DBL2 x sides ratios, DBL2 y sides ratios, DBL2 z sides ratios, DBL3 inner mimimum, outer maximum, exp sigma in nm
		vargenerator_descriptor.push_back("abl_exp", "0.25, 0.25; 0, 0; 0, 0; 0, 1, 10.0", MATPVAR_ABLEXP);

		//--------------

		meshDisplay_unit.push_back("", MESHDISPLAY_NONE);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_SM_DEMAG);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_SM_OERSTED);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_SM_STRAYH);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_MAGNETIZATION);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_MAGNETIZATION2);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_MAGNETIZATION12);
		meshDisplay_unit.push_back("J/T", MESHDISPLAY_MOMENT);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_EFFECTIVEFIELD);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_EFFECTIVEFIELD2);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_EFFECTIVEFIELD12);
		meshDisplay_unit.push_back("J/m3", MESHDISPLAY_ENERGY);
		meshDisplay_unit.push_back("J/m3", MESHDISPLAY_ENERGY2);
		meshDisplay_unit.push_back("A/m2", MESHDISPLAY_CURRDENSITY);
		meshDisplay_unit.push_back("V", MESHDISPLAY_VOLTAGE);
		meshDisplay_unit.push_back("S/m", MESHDISPLAY_ELCOND);
		meshDisplay_unit.push_back("A/m", MESHDISPLAY_SACCUM);
		meshDisplay_unit.push_back("A/s", MESHDISPLAY_JSX);
		meshDisplay_unit.push_back("A/s", MESHDISPLAY_JSY);
		meshDisplay_unit.push_back("A/s", MESHDISPLAY_JSZ);
		meshDisplay_unit.push_back("A/ms", MESHDISPLAY_TS);
		meshDisplay_unit.push_back("A/ms", MESHDISPLAY_TSI);
		meshDisplay_unit.push_back("K", MESHDISPLAY_TEMPERATURE);
		meshDisplay_unit.push_back("m", MESHDISPLAY_UDISP);
		meshDisplay_unit.push_back("", MESHDISPLAY_STRAINDIAG);
		meshDisplay_unit.push_back("", MESHDISPLAY_STRAINODIAG);
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
		displayHandles.push_back("mu", MESHDISPLAY_MOMENT);
		displayHandles.push_back("Heff", MESHDISPLAY_EFFECTIVEFIELD);
		displayHandles.push_back("Heff2", MESHDISPLAY_EFFECTIVEFIELD2);
		displayHandles.push_back("Heff12", MESHDISPLAY_EFFECTIVEFIELD12);
		displayHandles.push_back("Ed", MESHDISPLAY_ENERGY);
		displayHandles.push_back("Ed2", MESHDISPLAY_ENERGY2);
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
		displayHandles.push_back("u", MESHDISPLAY_UDISP);
		displayHandles.push_back("S_d", MESHDISPLAY_STRAINDIAG);
		displayHandles.push_back("S_od", MESHDISPLAY_STRAINODIAG);
		displayHandles.push_back("ParamVar", MESHDISPLAY_PARAMVAR);
		displayHandles.push_back("Roughness", MESHDISPLAY_ROUGHNESS);
		displayHandles.push_back("Cust_V", MESHDISPLAY_CUSTOM_VEC);
		displayHandles.push_back("Cust_S", MESHDISPLAY_CUSTOM_SCA);

		//--------------

		//specify forbidden module combinations - each set is an exclusive modules set, i.e. only one of each can be active at any one time
		//also specify non-exclusive modules (i.e. entries in exclusiveModules with only one entry per set)
		exclusiveModules.storeset(MOD_DEMAG_N, MOD_DEMAG, MOD_SDEMAG_DEMAG, MOD_ATOM_DIPOLEDIPOLE);
		exclusiveModules.storeset(MOD_EXCHANGE, MOD_DMEXCHANGE, MOD_IDMEXCHANGE);
		exclusiveModules.storeset(MOD_SURFEXCHANGE);
		exclusiveModules.storeset(MOD_ANIUNI, MOD_ANICUBI, MOD_ANIBI, MOD_ANITENS);
		exclusiveModules.storeset(MOD_MELASTIC);
		exclusiveModules.storeset(MOD_ZEEMAN);
		exclusiveModules.storeset(MOD_MOPTICAL);
		exclusiveModules.storeset(MOD_TRANSPORT);
		exclusiveModules.storeset(MOD_HEAT);
		exclusiveModules.storeset(MOD_SOTFIELD);
		exclusiveModules.storeset(MOD_STFIELD);
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

		//FERROMAGNETIC
		modules_for_meshtype.push_back(make_vector(
			MOD_DEMAG_N, MOD_DEMAG, MOD_SDEMAG_DEMAG,
			MOD_EXCHANGE, MOD_DMEXCHANGE, MOD_IDMEXCHANGE, MOD_SURFEXCHANGE, 
			MOD_ZEEMAN, MOD_MOPTICAL, MOD_MELASTIC, MOD_ROUGHNESS,
			MOD_ANIUNI, MOD_ANICUBI, MOD_ANIBI, MOD_ANITENS, 
			MOD_TRANSPORT,
			MOD_HEAT, 
			MOD_SOTFIELD, MOD_STFIELD), MESH_FERROMAGNETIC);

		displaymodules_for_meshtype.push_back(make_vector(
			MOD_DEMAG_N, MOD_DEMAG, MOD_SDEMAG_DEMAG,
			MOD_EXCHANGE, MOD_DMEXCHANGE, MOD_IDMEXCHANGE, MOD_SURFEXCHANGE,
			MOD_ZEEMAN, MOD_MOPTICAL, MOD_MELASTIC, MOD_ROUGHNESS,
			MOD_ANIUNI, MOD_ANICUBI, MOD_ANIBI, MOD_ANITENS,
			MOD_SOTFIELD, MOD_STFIELD), MESH_FERROMAGNETIC);

		//ANTIFERROMAGNETIC
		modules_for_meshtype.push_back(make_vector(
			MOD_DEMAG_N, MOD_DEMAG, MOD_SDEMAG_DEMAG,
			MOD_EXCHANGE, MOD_DMEXCHANGE, MOD_IDMEXCHANGE, MOD_SURFEXCHANGE,
			MOD_ZEEMAN, MOD_MOPTICAL,
			MOD_ANIUNI, MOD_ANICUBI, MOD_ANIBI, MOD_ANITENS,
			MOD_TRANSPORT,
			MOD_HEAT,
			MOD_SOTFIELD, MOD_ROUGHNESS), MESH_ANTIFERROMAGNETIC);

		displaymodules_for_meshtype.push_back(make_vector(
			MOD_DEMAG_N, MOD_DEMAG, MOD_SDEMAG_DEMAG,
			MOD_EXCHANGE, MOD_DMEXCHANGE, MOD_IDMEXCHANGE, MOD_SURFEXCHANGE,
			MOD_ZEEMAN, MOD_MOPTICAL,
			MOD_ANIUNI, MOD_ANICUBI, MOD_ANIBI, MOD_ANITENS,
			MOD_SOTFIELD, MOD_ROUGHNESS), MESH_ANTIFERROMAGNETIC);

		//DIAMAGNETIC
		modules_for_meshtype.push_back(make_vector(
			MOD_DEMAG, MOD_SDEMAG_DEMAG,
			MOD_SURFEXCHANGE,
			MOD_ZEEMAN, MOD_MOPTICAL,
			MOD_TRANSPORT,
			MOD_HEAT), MESH_DIAMAGNETIC);

		displaymodules_for_meshtype.push_back(make_vector(
			MOD_DEMAG, MOD_SDEMAG_DEMAG,
			MOD_SURFEXCHANGE,
			MOD_ZEEMAN, MOD_MOPTICAL), MESH_DIAMAGNETIC);

		//DIPOLE
		modules_for_meshtype.push_back(make_vector(MOD_TRANSPORT, MOD_HEAT), MESH_DIPOLE);
		displaymodules_for_meshtype.push_back({}, MESH_DIPOLE);

		//METAL
		modules_for_meshtype.push_back(make_vector(MOD_TRANSPORT, MOD_HEAT), MESH_METAL);
		displaymodules_for_meshtype.push_back({}, MESH_METAL);

		//INSULATOR
		modules_for_meshtype.push_back(make_vector(MOD_HEAT), MESH_INSULATOR);
		displaymodules_for_meshtype.push_back({}, MESH_INSULATOR);

		//ATOMISTIC SIMPLE CUBIC
		modules_for_meshtype.push_back(make_vector(
			MOD_DEMAG_N, MOD_DEMAG, MOD_ATOM_DIPOLEDIPOLE,
			MOD_EXCHANGE, MOD_DMEXCHANGE, MOD_IDMEXCHANGE,
			MOD_ZEEMAN, MOD_MOPTICAL,
			MOD_ANIUNI, MOD_ANICUBI, MOD_ANIBI, MOD_ANITENS,
			MOD_HEAT), MESH_ATOM_CUBIC);

		displaymodules_for_meshtype.push_back(make_vector(
			MOD_DEMAG_N, MOD_DEMAG, MOD_ATOM_DIPOLEDIPOLE,
			MOD_EXCHANGE, MOD_DMEXCHANGE, MOD_IDMEXCHANGE,
			MOD_ZEEMAN, MOD_MOPTICAL,
			MOD_ANIUNI, MOD_ANICUBI, MOD_ANIBI, MOD_ANITENS), MESH_ATOM_CUBIC);

		//----------------

		meshAllowedDisplay.push_back(make_vector(MESHDISPLAY_NONE, MESHDISPLAY_SM_DEMAG, MESHDISPLAY_SM_OERSTED, MESHDISPLAY_SM_STRAYH), MESH_SUPERMESH);

		meshAllowedDisplay.push_back(make_vector(
			MESHDISPLAY_NONE, MESHDISPLAY_MAGNETIZATION, MESHDISPLAY_EFFECTIVEFIELD, MESHDISPLAY_ENERGY,
			MESHDISPLAY_CURRDENSITY, MESHDISPLAY_VOLTAGE, MESHDISPLAY_ELCOND, MESHDISPLAY_SACCUM, MESHDISPLAY_JSX, MESHDISPLAY_JSY, MESHDISPLAY_JSZ, MESHDISPLAY_TS, MESHDISPLAY_TSI,
			MESHDISPLAY_TEMPERATURE, MESHDISPLAY_UDISP, MESHDISPLAY_STRAINDIAG, MESHDISPLAY_STRAINODIAG, MESHDISPLAY_PARAMVAR, MESHDISPLAY_ROUGHNESS, MESHDISPLAY_CUSTOM_VEC, MESHDISPLAY_CUSTOM_SCA), MESH_FERROMAGNETIC);

		meshAllowedDisplay.push_back(make_vector(
			MESHDISPLAY_NONE, MESHDISPLAY_MAGNETIZATION, MESHDISPLAY_MAGNETIZATION2, MESHDISPLAY_MAGNETIZATION12, MESHDISPLAY_EFFECTIVEFIELD, MESHDISPLAY_EFFECTIVEFIELD2, MESHDISPLAY_EFFECTIVEFIELD12, MESHDISPLAY_ENERGY, MESHDISPLAY_ENERGY2,
			MESHDISPLAY_CURRDENSITY, MESHDISPLAY_VOLTAGE, MESHDISPLAY_ELCOND,
			MESHDISPLAY_TEMPERATURE,
			MESHDISPLAY_PARAMVAR, MESHDISPLAY_ROUGHNESS, MESHDISPLAY_CUSTOM_VEC, MESHDISPLAY_CUSTOM_SCA), MESH_ANTIFERROMAGNETIC);

		meshAllowedDisplay.push_back(make_vector(
			MESHDISPLAY_NONE, MESHDISPLAY_MAGNETIZATION, MESHDISPLAY_EFFECTIVEFIELD, MESHDISPLAY_ENERGY,
			MESHDISPLAY_CURRDENSITY, MESHDISPLAY_VOLTAGE, MESHDISPLAY_ELCOND, MESHDISPLAY_SACCUM, MESHDISPLAY_JSX, MESHDISPLAY_JSY, MESHDISPLAY_JSZ,
			MESHDISPLAY_TEMPERATURE, MESHDISPLAY_PARAMVAR), MESH_DIAMAGNETIC);

		meshAllowedDisplay.push_back(make_vector(
			MESHDISPLAY_NONE, MESHDISPLAY_MAGNETIZATION, 
			MESHDISPLAY_CURRDENSITY, MESHDISPLAY_VOLTAGE, MESHDISPLAY_ELCOND, MESHDISPLAY_SACCUM, MESHDISPLAY_JSX, MESHDISPLAY_JSY, MESHDISPLAY_JSZ,
			MESHDISPLAY_TEMPERATURE, MESHDISPLAY_PARAMVAR), MESH_DIPOLE);

		meshAllowedDisplay.push_back(make_vector(
			MESHDISPLAY_NONE, 
			MESHDISPLAY_CURRDENSITY, MESHDISPLAY_VOLTAGE, MESHDISPLAY_ELCOND, MESHDISPLAY_SACCUM, MESHDISPLAY_JSX, MESHDISPLAY_JSY, MESHDISPLAY_JSZ,
			MESHDISPLAY_TEMPERATURE, MESHDISPLAY_PARAMVAR), MESH_METAL);

		meshAllowedDisplay.push_back(make_vector(MESHDISPLAY_NONE, MESHDISPLAY_TEMPERATURE, MESHDISPLAY_PARAMVAR), MESH_INSULATOR);

		meshAllowedDisplay.push_back(make_vector(
			MESHDISPLAY_NONE, MESHDISPLAY_MOMENT, MESHDISPLAY_EFFECTIVEFIELD, MESHDISPLAY_ENERGY,
			MESHDISPLAY_TEMPERATURE, MESHDISPLAY_PARAMVAR, MESHDISPLAY_CUSTOM_VEC, MESHDISPLAY_CUSTOM_SCA), MESH_ATOM_CUBIC);

		//----------------

		//assign possible parameters for each mesh type
		params_for_meshtype.push_back(make_vector(
			PARAM_GREL, PARAM_GDAMPING, PARAM_MS, PARAM_DEMAGXY, 
			PARAM_A, PARAM_D, PARAM_J1, PARAM_J2, 
			PARAM_K1, PARAM_K2, PARAM_K3, PARAM_EA1, PARAM_EA2, PARAM_EA3,
			PARAM_TC, PARAM_MUB, PARAM_SUSREL,
			PARAM_HA, PARAM_HMO,
			PARAM_ELC, PARAM_AMR, PARAM_P, PARAM_BETA, PARAM_DE, PARAM_NDENSITY, PARAM_SHA, PARAM_FLSOT, PARAM_STQ, PARAM_STA, PARAM_STP, PARAM_BETAD, PARAM_LSF, PARAM_LEX, PARAM_LPH, PARAM_GI, PARAM_GMIX, PARAM_PUMPEFF, PARAM_CPUMP_EFF, PARAM_THE_EFF, PARAM_TSEFF, PARAM_TSIEFF,
			PARAM_THERMCOND, PARAM_DENSITY, PARAM_MECOEFF, PARAM_YOUNGSMOD, PARAM_POISSONRATIO, 
			PARAM_SHC, PARAM_SHC_E, PARAM_G_E, PARAM_T, PARAM_Q), MESH_FERROMAGNETIC);

		params_for_meshtype.push_back(make_vector(
			PARAM_GREL_AFM, PARAM_GDAMPING_AFM, PARAM_MS_AFM, PARAM_DEMAGXY, 
			PARAM_A_AFM, PARAM_A_AFH, PARAM_A_AFNH, PARAM_D_AFM, PARAM_AFTAU, PARAM_AFTAUCROSS, PARAM_J1, PARAM_J2, 
			PARAM_K1_AFM, PARAM_K2_AFM, PARAM_K3_AFM, PARAM_EA1, PARAM_EA2, PARAM_EA3,
			PARAM_TC, PARAM_MUB_AFM, PARAM_SUSREL_AFM,
			PARAM_HA, PARAM_HMO,
			PARAM_ELC, PARAM_P, PARAM_BETA, PARAM_SHA, PARAM_FLSOT,
			PARAM_THERMCOND, PARAM_DENSITY, PARAM_SHC, 
			PARAM_SHC_E, PARAM_G_E, PARAM_T, PARAM_Q), MESH_ANTIFERROMAGNETIC);

		params_for_meshtype.push_back(make_vector(PARAM_NETADIA, PARAM_SUSREL, PARAM_HA, PARAM_HMO,
			PARAM_ELC, PARAM_DE, PARAM_NDENSITY, PARAM_SHA, PARAM_ISHA, PARAM_LSF, PARAM_GI, PARAM_GMIX,
			PARAM_THERMCOND, PARAM_DENSITY, PARAM_SHC, PARAM_SHC_E, PARAM_G_E, PARAM_T, PARAM_Q), MESH_DIAMAGNETIC);

		params_for_meshtype.push_back(make_vector(PARAM_MS, PARAM_TC,
			PARAM_ELC, PARAM_AMR, PARAM_P, PARAM_DE, PARAM_NDENSITY, PARAM_BETAD, PARAM_LSF, PARAM_LEX, PARAM_LPH, PARAM_GI, PARAM_GMIX,
			PARAM_THERMCOND, PARAM_DENSITY, PARAM_SHC, PARAM_SHC_E, PARAM_G_E, PARAM_T, PARAM_Q), MESH_DIPOLE);

		params_for_meshtype.push_back(make_vector(
			PARAM_ELC, PARAM_DE, PARAM_NDENSITY, PARAM_SHA, PARAM_ISHA, PARAM_LSF, PARAM_GI, PARAM_GMIX,
			PARAM_THERMCOND, PARAM_DENSITY, PARAM_SHC, PARAM_SHC_E, PARAM_G_E, PARAM_T, PARAM_Q), MESH_METAL);

		params_for_meshtype.push_back(make_vector(PARAM_THERMCOND, PARAM_DENSITY, PARAM_SHC), MESH_INSULATOR);

		params_for_meshtype.push_back(make_vector(
			PARAM_ATOM_SC_DAMPING, PARAM_ATOM_SC_MUS, PARAM_DEMAGXY,
			PARAM_ATOM_SC_J, PARAM_ATOM_SC_D,
			PARAM_ATOM_SC_K1, PARAM_ATOM_SC_K2, PARAM_ATOM_SC_K3, PARAM_ATOM_EA1, PARAM_ATOM_EA2, PARAM_ATOM_EA3,
			PARAM_HA, PARAM_HMO,
			PARAM_ELC,
			PARAM_THERMCOND, PARAM_DENSITY,
			PARAM_SHC, PARAM_SHC_E, PARAM_G_E, PARAM_T, PARAM_Q), MESH_ATOM_CUBIC);

		//there's also an entry for MESH_SUPERMESH : this includes all possible material parameters which could be saved in a materials data base
		params_for_meshtype.push_back(make_vector(PARAM_GREL, PARAM_GDAMPING, PARAM_MS, PARAM_A, PARAM_D, PARAM_J1, PARAM_J2, PARAM_K1, PARAM_K2, PARAM_EA1, PARAM_EA2, PARAM_TC, PARAM_MUB,
			PARAM_ELC, PARAM_AMR, PARAM_P, PARAM_BETA, PARAM_DE, PARAM_BETAD, PARAM_SHA, PARAM_ISHA, PARAM_FLSOT, PARAM_LSF, PARAM_LEX, PARAM_LPH, PARAM_GI, PARAM_GMIX,
			PARAM_THERMCOND, PARAM_DENSITY, PARAM_SHC), MESH_SUPERMESH);

		//----------------

		//entries from PARAM_, specifying if temperature dependence (first) and spatial variation (second) are enabled.
		params_enabled_props.push_back({ true, true }, PARAM_GREL);
		params_enabled_props.push_back({ true, true }, PARAM_GDAMPING);
		params_enabled_props.push_back({ true, true }, PARAM_MS);
		params_enabled_props.push_back({ false, false }, PARAM_DEMAGXY);
		params_enabled_props.push_back({ true, true }, PARAM_A);
		params_enabled_props.push_back({ true, true }, PARAM_D);
		params_enabled_props.push_back({ true, true }, PARAM_J1);
		params_enabled_props.push_back({ true, true }, PARAM_J2);
		params_enabled_props.push_back({ true, true }, PARAM_K1);
		params_enabled_props.push_back({ true, true }, PARAM_K2);
		params_enabled_props.push_back({ true, true }, PARAM_K3);
		params_enabled_props.push_back({ false, true }, PARAM_EA1);
		params_enabled_props.push_back({ false, true }, PARAM_EA2);
		params_enabled_props.push_back({ false, true }, PARAM_EA3);
		params_enabled_props.push_back({ true, false }, PARAM_SUSREL);
		params_enabled_props.push_back({ false, false }, PARAM_SUSPREL);
		params_enabled_props.push_back({ true, true }, PARAM_ELC);
		params_enabled_props.push_back({ true, true }, PARAM_AMR);
		params_enabled_props.push_back({ true, true }, PARAM_P);
		params_enabled_props.push_back({ true, true }, PARAM_BETA);
		params_enabled_props.push_back({ true, true }, PARAM_DE);
		params_enabled_props.push_back({ true, true }, PARAM_BETAD);
		params_enabled_props.push_back({ true, true }, PARAM_SHA);
		params_enabled_props.push_back({ true, true }, PARAM_ISHA);
		params_enabled_props.push_back({ true, true }, PARAM_LSF);
		params_enabled_props.push_back({ true, true }, PARAM_LEX);
		params_enabled_props.push_back({ true, true }, PARAM_LPH);
		params_enabled_props.push_back({ true, true }, PARAM_GI);
		params_enabled_props.push_back({ true, true }, PARAM_GMIX);
		params_enabled_props.push_back({ true, true }, PARAM_TSEFF);
		params_enabled_props.push_back({ true, true }, PARAM_TSIEFF);
		params_enabled_props.push_back({ true, true }, PARAM_PUMPEFF);
		params_enabled_props.push_back({ true, true }, PARAM_THERMCOND);
		params_enabled_props.push_back({ true, true }, PARAM_DENSITY);
		params_enabled_props.push_back({ true, true }, PARAM_SHC);
		params_enabled_props.push_back({ true, true }, PARAM_FLSOT);
		params_enabled_props.push_back({ true, false }, PARAM_STQ);
		params_enabled_props.push_back({ true, false }, PARAM_STA);
		params_enabled_props.push_back({ true, false }, PARAM_STP);
		params_enabled_props.push_back({ true, true }, PARAM_HA);
		params_enabled_props.push_back({ false, true }, PARAM_TC);
		params_enabled_props.push_back({ false, true }, PARAM_MUB);
		params_enabled_props.push_back({ false, true }, PARAM_T);
		params_enabled_props.push_back({ false, true }, PARAM_Q);
		params_enabled_props.push_back({ true, true }, PARAM_GREL_AFM);
		params_enabled_props.push_back({ true, true }, PARAM_GDAMPING_AFM);
		params_enabled_props.push_back({ true, true }, PARAM_MS_AFM);
		params_enabled_props.push_back({ true, true }, PARAM_A_AFM);
		params_enabled_props.push_back({ true, true }, PARAM_A_AFNH);
		params_enabled_props.push_back({ true, true }, PARAM_D_AFM);
		params_enabled_props.push_back({ true, true }, PARAM_CPUMP_EFF);
		params_enabled_props.push_back({ true, true }, PARAM_THE_EFF);
		params_enabled_props.push_back({ false, true }, PARAM_NDENSITY);
		params_enabled_props.push_back({ true, true }, PARAM_MECOEFF);
		params_enabled_props.push_back({ true, true }, PARAM_YOUNGSMOD);
		params_enabled_props.push_back({ true, true }, PARAM_POISSONRATIO);
		params_enabled_props.push_back({ true, true }, PARAM_NETADIA);
		params_enabled_props.push_back({ true, true }, PARAM_SHC_E);
		params_enabled_props.push_back({ true, true }, PARAM_G_E);
		params_enabled_props.push_back({ true, true }, PARAM_A_AFH);
		params_enabled_props.push_back({ true, true }, PARAM_SUSREL_AFM);
		params_enabled_props.push_back({ false, true }, PARAM_AFTAU);
		params_enabled_props.push_back({ false, true }, PARAM_AFTAUCROSS);
		params_enabled_props.push_back({ false, true }, PARAM_MUB_AFM);
		params_enabled_props.push_back({ true, true }, PARAM_K1_AFM);
		params_enabled_props.push_back({ true, true }, PARAM_K2_AFM);
		params_enabled_props.push_back({ true, true }, PARAM_K3_AFM);
		params_enabled_props.push_back({ true, true }, PARAM_HMO);
		params_enabled_props.push_back({ false, true }, PARAM_ATOM_SC_DAMPING);
		params_enabled_props.push_back({ false, true }, PARAM_ATOM_SC_MUS);
		params_enabled_props.push_back({ false, true }, PARAM_ATOM_SC_J);
		params_enabled_props.push_back({ false, true }, PARAM_ATOM_SC_D);
		params_enabled_props.push_back({ false, true }, PARAM_ATOM_SC_K1);
		params_enabled_props.push_back({ false, true }, PARAM_ATOM_SC_K2);
		params_enabled_props.push_back({ false, true }, PARAM_ATOM_SC_K3);
		params_enabled_props.push_back({ false, true }, PARAM_ATOM_EA1);
		params_enabled_props.push_back({ false, true }, PARAM_ATOM_EA2);
		params_enabled_props.push_back({ false, true }, PARAM_ATOM_EA3);
	}
}