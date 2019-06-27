#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "MeshParamsCUDA.h"

class MeshCUDA;

//This holds pointers to managed objects in MeshCUDA (and inherited MeshParamsCUDA) : set and forget. They are available for use in cuda kernels by passing a cu_obj-managed object ManagedMeshCUDA

class ManagedMeshCUDA {

public:

	//Material Parameters

	//Relative electron gyromagnetic ratio
	MatPCUDA<cuReal, cuReal>* pgrel;

	//Gilbert damping
	MatPCUDA<cuReal, cuReal>* palpha;

	//Saturation magnetisation (A/m)
	MatPCUDA<cuReal, cuReal>* pMs;

	//in-plane demagnetizing factors (used for Demag_N module)
	MatPCUDA<cuReal2, cuReal>* pNxy;

	//Exchange stifness (J/m)
	MatPCUDA<cuReal, cuReal>* pA;

	//Dzyaloshinskii-Moriya exchange constant (J/m2)
	MatPCUDA<cuReal, cuReal>* pD;

	//bilinear surface exchange coupling (J/m^2) : J1, bottom and top layer values
	//biquadratic surface exchange coupling (J/m^2) : J2, bottom and top layer values
	MatPCUDA<cuReal, cuReal>* pJ1;
	MatPCUDA<cuReal, cuReal>* pJ2;

	//Magneto-crystalline anisotropy K1 and K2 constants (J/m^3) and easy axes directions. For uniaxial anisotropy only ea1 is needed, for cubic ea1 and ea2 should be orthogonal.
	MatPCUDA<cuReal, cuReal>* pK1;
	MatPCUDA<cuReal, cuReal>* pK2;
	MatPCUDA<cuReal3, cuReal3>* pmcanis_ea1;
	MatPCUDA<cuReal3, cuReal3>* pmcanis_ea2;

	//longitudinal (parallel) susceptibility relative to mu0*Ms0, i.e. divided by mu0*Ms0, Ms0 is the 0K Ms value - for use with LLB equation. Units As^2/kg
	MatPCUDA<cuReal, cuReal>* psusrel;

	//perpendicular (transverse) susceptibility relative to mu0*Ms0, i.e. divided by mu0*Ms0, Ms0 is the 0K Ms value - for use with LLB equation. Units As^2/kg
	MatPCUDA<cuReal, cuReal>* psusprel;

	//applied field spatial variation coefficient (unitless)
	MatPCUDA<cuReal, cuReal>* pcHA;

	//electrical conductivity (units S/m).
	//this is the value at 0K for Ni80Fe20. Temperature dependence typically scaled by 1 / (1 + alpha*(T-T0)), where alpha = 0.003, T0 = 293K with sigma = 1.7e6 S/m and 293K.
	//Using scaling 1 / (1 + alpha0 * T) on the zero-temperature conductivity gives sigma0 = sigmaT0 / (1 - alpha*T0), alpha0 = alpha / (1 - alpha*T0), so alpha0 = 0.025.
	MatPCUDA<cuReal, cuReal>* pelecCond;

	//anisotropic magnetoresistance as a percentage (of base resistance)
	MatPCUDA<cuReal, cuReal>* pamrPercentage;

	//spin current polarization and non-adiabaticity (for Zhang-Li STT).
	MatPCUDA<cuReal, cuReal>* pP;
	MatPCUDA<cuReal, cuReal>* pbeta;

	MatPCUDA<cuReal, cuReal>* pDe;
	MatPCUDA<cuReal, cuReal>* pbetaD;

	MatPCUDA<cuReal, cuReal>* pSHA;
	MatPCUDA<cuReal, cuReal>* piSHA;
	MatPCUDA<cuReal, cuReal>* pflSOT;

	MatPCUDA<cuReal, cuReal>* pl_sf;
	MatPCUDA<cuReal, cuReal>* pl_ex;
	MatPCUDA<cuReal, cuReal>* pl_ph;

	MatPCUDA<cuReal2, cuReal>* pGi;
	MatPCUDA<cuReal2, cuReal>* pGmix;

	MatPCUDA<cuReal, cuReal>* pts_eff;
	MatPCUDA<cuReal, cuReal>* ptsi_eff;

	MatPCUDA<cuReal, cuReal>* ppump_eff;

	//the mesh base temperature (K)
	cuReal* pbase_temperature;

	//Curie temperature (K)
	cuReal* pT_Curie;

	//thermal conductivity (W/mK) - default for permalloy
	MatPCUDA<cuReal, cuReal>* pthermCond;

	//mass density (kg/m^3) - default for permalloy
	MatPCUDA<cuReal, cuReal>* pdensity;

	//specific heat capacity (J/kgK) - default for permalloy
	MatPCUDA<cuReal, cuReal>* pshc;

	//-----Ferromagnetic properties

	//Magnetization
	cuVEC_VC<cuReal3>* pM;

	//effective field - sum total field of all the added modules
	cuVEC<cuReal3>* pHeff;

	//-----Electric conduction properties (Electron charge and spin Transport)

	//electrical potential - on n_e, h_e mesh
	cuVEC_VC<cuReal>* pV;

	//electrical conductivity - on n_e, h_e mesh
	cuVEC_VC<cuReal>* pelC;

	//electrical current density - on n_e, h_e mesh
	cuVEC_VC<cuReal3>* pJc;

	//spin accumulation - on n_e, h_e mesh
	cuVEC_VC<cuReal3>* pS;

	//-----Thermal conduction properties

	//temperature calculated by Heat module
	cuVEC_VC<cuReal>* pTemp;

private:

	//----------------------------------- RUNTIME PARAMETER UPDATERS (AUXILIARY) (MeshParamsControlCUDA.h)

	//UPDATER M CORSENESS - PRIVATE

	//SPATIAL DEPENDENCE ONLY - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_mcoarse_spatial(int mcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version; position not calculated
	template <typename PType, typename SType>
	__device__ void update_parameters_mcoarse_spatial(int mcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TIME DEPENDENCE - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_mcoarse_full(int mcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	__device__ void update_parameters_mcoarse_full(int mcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value);

	//UPDATER E CORSENESS - PRIVATE

	//SPATIAL DEPENDENCE ONLY - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_ecoarse_spatial(int ecell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version; position not calculated
	template <typename PType, typename SType>
	__device__ void update_parameters_ecoarse_spatial(int ecell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TIME DEPENDENCE - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_ecoarse_full(int ecell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	__device__ void update_parameters_ecoarse_full(int ecell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value);

	//UPDATER T CORSENESS - PRIVATE

	//SPATIAL DEPENDENCE ONLY - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_tcoarse_spatial(int tcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version; position not calculated
	template <typename PType, typename SType>
	__device__ void update_parameters_tcoarse_spatial(int tcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TIME DEPENDENCE - NO POSITION YET

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_tcoarse_full(int tcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	__device__ void update_parameters_tcoarse_full(int tcell_idx, MatPCUDA<PType, SType>& matp, PType& matp_value);

public:

	void construct_cu_obj(void) {}

	void destruct_cu_obj(void) {}

	BError set_pointers(MeshCUDA* pMeshCUDA);

	//----------------------------------- RUNTIME PARAMETER UPDATERS (MeshParamsControlCUDA.h)

	//SPATIAL DEPENDENCE ONLY - HAVE POSITION

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_spatial(const cuReal3& position, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	__device__ void update_parameters_spatial(const cuReal3& position, MatPCUDA<PType, SType>& matp, PType& matp_value);

	//SPATIAL AND TIME DEPENDENCE - HAVE POSITION

	//update parameters in the list for spatial dependence only
	template <typename PType, typename SType, typename ... MeshParam_List>
	__device__ void update_parameters_full(const cuReal3& position, const cuReal& Temperature, MatPCUDA<PType, SType>& matp, PType& matp_value, MeshParam_List& ... params);

	//update parameters in the list for spatial dependence only - single parameter version
	template <typename PType, typename SType>
	__device__ void update_parameters_full(const cuReal3& position, const cuReal& Temperature, MatPCUDA<PType, SType>& matp, PType& matp_value);

	//UPDATER M CORSENESS - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	__device__ void update_parameters_mcoarse(int mcell_idx, MeshParam_List& ... params);

	//UPDATER E CORSENESS - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	__device__ void update_parameters_ecoarse(int ecell_idx, MeshParam_List& ... params);

	//UPDATER T CORSENESS - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	__device__ void update_parameters_tcoarse(int tcell_idx, MeshParam_List& ... params);

	//UPDATER POSITION KNOWN - PUBLIC

	//Update parameter values if temperature dependent at the given cell index - M cell index; position not calculated
	template <typename ... MeshParam_List>
	__device__ void update_parameters_atposition(const cuReal3& position, MeshParam_List& ... params);
};

#endif