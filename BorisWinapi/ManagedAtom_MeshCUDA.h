#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "Atom_MeshParamsCUDA.h"

#include "ManagedAtom_DiffEq_CommonCUDA.h"

class Atom_MeshCUDA;

//This holds pointers to managed objects in Atom_MeshCUDA (and inherited Atom_MeshParamsCUDA) : set and forget. They are available for use in cuda kernels by passing a cu_obj-managed object ManagedAtom_MeshCUDA

class ManagedAtom_MeshCUDA {

public:

	//Material Parameters

	//-----------SIMPLE CUBIC

	//Gilbert damping (atomistic: intrinsic)
	MatPCUDA<cuBReal, cuBReal>* palpha;

	//atomic moment (units of muB) - default for bcc Fe
	MatPCUDA<cuBReal, cuBReal>* pmu_s;

	//Exchange constant : (units of J) - default for bcc Fe
	MatPCUDA<cuBReal, cuBReal>* pJ;

	//DMI exchange constant : (units of J)
	MatPCUDA<cuBReal, cuBReal>* pD;

	//Magneto-crystalline anisotropy constants (J) and easy axes directions. For uniaxial anisotropy only ea1 is needed.
	MatPCUDA<cuBReal, cuBReal>* pK;

	//Magneto-crystalline anisotropy easy axes directions
	MatPCUDA<cuReal3, cuReal3>* pmcanis_ea1;
	MatPCUDA<cuReal3, cuReal3>* pmcanis_ea2;

	//-----------BCC (2 per unit cell)

	//-----------FCC (4 per unit cell)

	//-----------HCP (4 per effective unit cell)

	//-----------Others

	//applied field spatial variation coefficient (unitless)
	MatPCUDA<cuBReal, cuBReal>* pcHA;

	//Magneto-Optical field strength (A/m)
	MatPCUDA<cuBReal, cuBReal>* pcHmo;

	//electrical conductivity (units S/m).
	//this is the value at RT for Ni80Fe20.
	MatPCUDA<cuBReal, cuBReal>* pelecCond;

	//thermal conductivity (W/mK) - default for permalloy
	MatPCUDA<cuBReal, cuBReal>* pthermCond;

	//the mesh base temperature (K)
	cuBReal* pbase_temperature;

	//mass density (kg/m^3) - default for permalloy
	MatPCUDA<cuBReal, cuBReal>* pdensity;

	//specific heat capacity (J/kgK) - default for permalloy
	MatPCUDA<cuBReal, cuBReal>* pshc;

	//electron specific heat capacity at room temperature used in many-temperature models (J/kgK); Note, if used you should assign a temperature dependence to it, e.g. linear with temperature for the free electron approximation; none assigned by default.
	MatPCUDA<cuBReal, cuBReal>* pshc_e;

	//electron-lattice coupling constant (W/m^3K) used in two-temperature model.
	MatPCUDA<cuBReal, cuBReal>* pG_e;

	//set temperature spatial variation coefficient (unitless) - used with temperature settings in a simulation schedule only, not with console command directly
	MatPCUDA<cuBReal, cuBReal>* pcT;

	//Heat source stimulus in heat equation. Ideally used with a spatial variation. (W//m3)
	MatPCUDA<cuBReal, cuBReal>* pQ;

	//-----Magnetic properties

	//Magnetization
	cuVEC_VC<cuReal3>* pM1;

	//effective field - sum total field of all the added modules
	cuVEC<cuReal3>* pHeff1;

	//-----Electric conduction properties (Electron charge and spin Transport)

	//electrical potential - on n_e, h_e mesh
	cuVEC_VC<cuBReal>* pV;

	//electrical conductivity - on n_e, h_e mesh
	cuVEC_VC<cuBReal>* pelC;

	//electrical field - on n_e, h_e mesh
	cuVEC_VC<cuReal3>* pE;

	//spin accumulation - on n_e, h_e mesh
	cuVEC_VC<cuReal3>* pS;

	//-----Thermal conduction properties

	//temperature calculated by Heat module (primary temperature, always used for 1-temperature model; for multi-temperature models in metals this is the itinerant electron temperature)
	cuVEC_VC<cuBReal>* pTemp;

	//lattice temperature used in many-T models
	cuVEC_VC<cuBReal>* pTemp_l;

	//mechanical displacement vectors - on n_m, h_m mesh
	cuVEC_VC<cuReal3>* pu_disp;

	//strain tensor (symmetric):
	//diagonal and off-diagonal components - on n_m, h_m mesh
	//xx, yy, zz
	cuVEC_VC<cuReal3>* pstrain_diag;
	//yz, xz, xy
	cuVEC_VC<cuReal3>* pstrain_odiag;

	//-----

	//Managed cuda mesh pointer so all mesh data can be accessed in device code
	ManagedAtom_DiffEq_CommonCUDA* pcuaDiffEq;

private:

	//----------------------------------- RUNTIME PARAMETER UPDATERS (AUXILIARY) (Atom_MeshParamsControlCUDA.h)

	//UPDATER M COARSENESS - PRIVATE

	//SPATIAL DEPENDENCE ONLY - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_mcoarse_spatial(int mcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version; position not calculated
	template <typename PType, typename SType>
	__device__ void update_parameters_mcoarse_spatial(int mcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TEMPERATURE DEPENDENCE - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_mcoarse_full(int mcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	__device__ void update_parameters_mcoarse_full(int mcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value);

	//UPDATER E COARSENESS - PRIVATE

	//SPATIAL DEPENDENCE ONLY - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_ecoarse_spatial(int ecell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version; position not calculated
	template <typename PType, typename SType>
	__device__ void update_parameters_ecoarse_spatial(int ecell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TEMPERATURE DEPENDENCE - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_ecoarse_full(int ecell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	__device__ void update_parameters_ecoarse_full(int ecell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value);

	//UPDATER T COARSENESS - PRIVATE

	//SPATIAL DEPENDENCE ONLY - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_tcoarse_spatial(int tcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version; position not calculated
	template <typename PType, typename SType>
	__device__ void update_parameters_tcoarse_spatial(int tcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TEMPERATURE DEPENDENCE - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_tcoarse_full(int tcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	__device__ void update_parameters_tcoarse_full(int tcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value);

public:

	void construct_cu_obj(void) {}

	void destruct_cu_obj(void) {}

	BError set_pointers(Atom_MeshCUDA* paMeshCUDA);

	//----------------------------------- RUNTIME PARAMETER UPDATERS (Atom_MeshParamsControlCUDA.h)

	//SPATIAL DEPENDENCE ONLY - HAVE POSITION

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_spatial(const cuReal3& position, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	__device__ void update_parameters_spatial(const cuReal3& position, MatPCUDA<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TEMPERATURE DEPENDENCE - HAVE POSITION

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_full(const cuReal3& position, const cuBReal& Temperature, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	__device__ void update_parameters_full(const cuReal3& position, const cuBReal& Temperature, MatPCUDA<PType, SType>& matp, PType& matp_value);

	//UPDATER M COARSENESS - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	__device__ void update_parameters_mcoarse(int mcell_idx, MeshParam_List& ... params);

	//UPDATER E COARSENESS - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	__device__ void update_parameters_ecoarse(int ecell_idx, MeshParam_List& ... params);

	//UPDATER T COARSENESS - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	__device__ void update_parameters_tcoarse(int tcell_idx, MeshParam_List& ... params);

	//UPDATER POSITION KNOWN - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	__device__ void update_parameters_atposition(const cuReal3& position, MeshParam_List& ... params);
};

#endif