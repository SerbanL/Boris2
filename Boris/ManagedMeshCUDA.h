#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "MeshParamsCUDA.h"

#include "ManagedDiffEq_CommonCUDA.h"

class MeshCUDA;
class ManagedAtom_MeshCUDA;

//This holds pointers to managed objects in MeshCUDA (and inherited MeshParamsCUDA) : set and forget. They are available for use in cuda kernels by passing a cu_obj-managed object ManagedMeshCUDA

class ManagedMeshCUDA {

public:

	//Material Parameters

	//Relative electron gyromagnetic ratio
	MatPCUDA<cuBReal, cuBReal>* pgrel;
	MatPCUDA<cuReal2, cuBReal>* pgrel_AFM;

	//Gilbert damping
	MatPCUDA<cuBReal, cuBReal>* palpha;
	MatPCUDA<cuReal2, cuBReal>* palpha_AFM;

	//Saturation magnetization (A/m)
	MatPCUDA<cuBReal, cuBReal>* pMs;
	MatPCUDA<cuReal2, cuBReal>* pMs_AFM;

	//in-plane demagnetizing factors (used for Demag_N module)
	MatPCUDA<cuReal2, cuBReal>* pNxy;

	//Exchange stiffness (J/m)
	MatPCUDA<cuBReal, cuBReal>* pA;
	MatPCUDA<cuReal2, cuBReal>* pA_AFM;
	
	//Homogeneous AFM coupling between sub-lattices A and B, defined as A / a*a (J/m^3), where A is the homogeneous antiferromagnetic exchange stifness (negative), and a is the lattice constant.
	//e.g. a = 0.3nm, A = -1pJ/m gives Ah as -1e7 J/m^3 to order of magnitude.
	MatPCUDA<cuReal2, cuBReal>* pAh;

	//Nonhomogeneous AFM coupling between sub-lattices A and B (J/m)
	MatPCUDA<cuReal2, cuBReal>* pAnh;

	//Dzyaloshinskii-Moriya exchange constant (J/m2)
	MatPCUDA<cuBReal, cuBReal>* pD;
	MatPCUDA<cuReal2, cuBReal>* pD_AFM;

	//Homogeneous DMI constant for 2-sublattice models (J/m^3)
	MatPCUDA<cuBReal, cuBReal>* pDh;
	//Homogeneous DMI term orientation
	MatPCUDA<cuReal3, cuReal3>* pdh_dir;

	//Interfacial DMI symmetry axis direction, used by vector interfacial DMI module
	MatPCUDA<cuReal3, cuReal3>* pD_dir;

	//Coupling between exchange integral and critical temperature (Neel or Curie temperature) for 2-sublattice model : intra-lattice term, 0.5 for ideal antiferromagnet
	//J = 3 * tau * kB * Tc
	MatPCUDA<cuReal2, cuBReal>* ptau_ii;

	//Coupling between exchange integral and critical temperature (Neel or Curie temperature) for 2-sublattice model : inter-lattice, or cross-lattice term, 0.5 for ideal antiferromagnet.
	//J = 3 * tau * kB * Tc
	MatPCUDA<cuReal2, cuBReal>* ptau_ij;

	//bilinear surface exchange coupling (J/m^2) : J1, bottom and top layer values
	//biquadratic surface exchange coupling (J/m^2) : J2, bottom and top layer values
	MatPCUDA<cuBReal, cuBReal>* pJ1;
	MatPCUDA<cuBReal, cuBReal>* pJ2;

	//Magneto-crystalline anisotropy K1 and K2 constants (J/m^3) and easy axes directions. For uniaxial anisotropy only ea1 is needed, for cubic ea1 and ea2 should be orthogonal.
	MatPCUDA<cuBReal, cuBReal>* pK1;
	MatPCUDA<cuBReal, cuBReal>* pK2;
	MatPCUDA<cuBReal, cuBReal>* pK3;
	MatPCUDA<cuReal3, cuReal3>* pmcanis_ea1;
	MatPCUDA<cuReal3, cuReal3>* pmcanis_ea2;
	MatPCUDA<cuReal3, cuReal3>* pmcanis_ea3;

	//Anisotropy values for 2-sublattice model
	MatPCUDA<cuReal2, cuBReal>* pK1_AFM;
	MatPCUDA<cuReal2, cuBReal>* pK2_AFM;
	MatPCUDA<cuReal2, cuBReal>* pK3_AFM;

	//Tensorial anisotropy terms
	cuVEC<cuReal4> *pKt, *pKt2;

	//longitudinal (parallel) susceptibility relative to mu0*Ms0, i.e. divided by mu0*Ms0, Ms0 is the 0K Ms value - for use with LLB equation. Units As^2/kg
	MatPCUDA<cuBReal, cuBReal>* psusrel;

	//longitudinal (parallel) susceptibility relative to mu0*Ms0, i.e. divided by mu0*Ms0, Ms0 is the 0K Ms value - for use with LLB equation 2-sublattice model. Units As^2/kg
	MatPCUDA<cuReal2, cuBReal>* psusrel_AFM;

	//perpendicular (transverse) susceptibility relative to mu0*Ms0, i.e. divided by mu0*Ms0, Ms0 is the 0K Ms value - for use with LLB equation. Units As^2/kg
	MatPCUDA<cuBReal, cuBReal>* psusprel;

	//applied field spatial variation coefficient (unitless)
	MatPCUDA<cuBReal, cuBReal>* pcHA;

	//Magneto-Optical field strength (A/m)
	MatPCUDA<cuBReal, cuBReal>* pcHmo;

	//Stochasticity efficiency parameter
	MatPCUDA<cuBReal, cuBReal>* ps_eff;

	//electrical conductivity (units S/m).
	//this is the value at 0K for Ni80Fe20. Temperature dependence typically scaled by 1 / (1 + alpha*(T-T0)), where alpha = 0.003, T0 = 293K with sigma = 1.7e6 S/m and 293K.
	//Using scaling 1 / (1 + alpha0 * T) on the zero-temperature conductivity gives sigma0 = sigmaT0 / (1 - alpha*T0), alpha0 = alpha / (1 - alpha*T0), so alpha0 = 0.025.
	MatPCUDA<cuBReal, cuBReal>* pelecCond;

	//TMR RA products for parallel and antiparallel states (Ohms m^2)
	MatPCUDA<cuBReal, cuBReal>* pRAtmr_p;
	MatPCUDA<cuBReal, cuBReal>* pRAtmr_ap;

	//anisotropic magnetoresistance as a percentage (of base resistance)
	MatPCUDA<cuBReal, cuBReal>* pamrPercentage;

	//spin current polarization and non-adiabaticity (for Zhang-Li STT).
	MatPCUDA<cuBReal, cuBReal>* pP;
	MatPCUDA<cuBReal, cuBReal>* pbeta;

	MatPCUDA<cuBReal, cuBReal>* pDe;
	MatPCUDA<cuBReal, cuBReal>* pn_density;
	MatPCUDA<cuBReal, cuBReal>* pbetaD;

	MatPCUDA<cuBReal, cuBReal>* pSHA;
	MatPCUDA<cuBReal, cuBReal>* piSHA;
	MatPCUDA<cuBReal, cuBReal>* pflSOT;

	//Slonczewski macrospin torques q+, q- parameters as in PRB 72, 014446 (2005) (unitless)
	MatPCUDA<cuReal2, cuBReal>* pSTq;

	//Slonczewski macrospin torques A, B parameters as in PRB 72, 014446 (2005) (unitless)
	MatPCUDA<cuReal2, cuBReal>* pSTa;

	//Slonczewski macrospin torques spin polarization unit vector as in PRB 72, 014446 (2005) (unitless)
	MatPCUDA<cuReal3, cuReal3>* pSTp;

	MatPCUDA<cuBReal, cuBReal>* pl_sf;
	MatPCUDA<cuBReal, cuBReal>* pl_ex;
	MatPCUDA<cuBReal, cuBReal>* pl_ph;

	MatPCUDA<cuReal2, cuBReal>* pGi;
	MatPCUDA<cuReal2, cuBReal>* pGmix;

	MatPCUDA<cuBReal, cuBReal>* pts_eff;
	MatPCUDA<cuBReal, cuBReal>* ptsi_eff;

	MatPCUDA<cuBReal, cuBReal>* ppump_eff;
	MatPCUDA<cuBReal, cuBReal>* pcpump_eff;
	MatPCUDA<cuBReal, cuBReal>* pthe_eff;

	//the mesh base temperature (K)
	cuBReal* pbase_temperature;

	//Curie temperature (K)
	cuBReal* pT_Curie;

	//The atomic magnetic moment as a multiple of the Bohr magneton - default 1 ub for permalloy.
	MatPCUDA<cuBReal, cuBReal>* patomic_moment;

	//atomic moments for 2-sublattice model (again multiples of the Bohr magneton)
	MatPCUDA<cuReal2, cuBReal>* patomic_moment_AFM;

	//thermal conductivity (W/mK) - default for permalloy
	MatPCUDA<cuBReal, cuBReal>* pthermCond;

	//mass density (kg/m^3) - default for permalloy
	MatPCUDA<cuBReal, cuBReal>* pdensity;

	//Magneto-elastic coefficients (J/m^3) - default for Ni
	MatPCUDA<cuReal2, cuBReal>* pMEc;

	//Young's modulus (Pa) - default for permalloy
	MatPCUDA<cuBReal, cuBReal>* pYm ;

	//Poisson's ratio (unitless) - default for permalloy
	MatPCUDA<cuBReal, cuBReal>* pPr;

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

	//-----Ferromagnetic properties

	//Magnetization
	cuVEC_VC<cuReal3>* pM;
	cuVEC_VC<cuReal3>* pM2;

	//effective field - sum total field of all the added modules
	cuVEC<cuReal3>* pHeff;
	cuVEC<cuReal3>* pHeff2;

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
	ManagedDiffEq_CommonCUDA* pcuDiffEq;

	//----------------------------------- MONTE-CARLO DATA FROM MODULES

	//SurfExchangeCUDA
	//arrays with pointers to other meshes in surface exchange coupling with the mesh holding this module, top and bottom (set from SurfExchangeCUDA)
	ManagedMeshCUDA* pMeshFM_Top;
	ManagedMeshCUDA* pMeshFM_Bot;
	size_t pMeshFM_Top_size, pMeshFM_Bot_size;

	ManagedMeshCUDA* pMeshAFM_Top;
	ManagedMeshCUDA* pMeshAFM_Bot;
	size_t pMeshAFM_Top_size, pMeshAFM_Bot_size;

	ManagedAtom_MeshCUDA* pMeshAtom_Top;
	ManagedAtom_MeshCUDA* pMeshAtom_Bot;
	size_t pMeshAtom_Top_size, pMeshAtom_Bot_size;

	//DemagCUDA (or SDemag_DemagCUDA)
	cuVEC<cuReal3>*  pDemag_Heff;

	//RoughnessCUDA
	cuVEC<cuReal3>* pFmul_rough;
	cuVEC<cuReal3>* pFomul_rough;

	//ZeemanCUDA
	cuVEC<cuReal3>* pHavec;

	//StrayField_MeshCUDA
	cuVEC<cuReal3>* pstrayField;

private:

	//----------------------------------- RUNTIME PARAMETER UPDATERS (AUXILIARY) (MeshParamsControlCUDA.h)

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

	BError set_pointers(MeshCUDA* pMeshCUDA);

	//----------------------------------- RUNTIME PARAMETER UPDATERS (MeshParamsControlCUDA.h)

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

	//----------------------------------- MONTE-CARLO METHODS FOR ENERGY COMPUTATION

	//Energy Deltas

	////////////////////////////////////
	//
	// Ferromagnetic

	//switch function which adds all assigned energy contributions in this mesh to calculate energy change from current spin to Mnew spin : return energy change as new - old
	//If Mnew is passed in as cuReal3(), then this function returns the current spin energy only - all functions in the switch statement implement this eventuality.
	__device__ cuBReal Get_EnergyChange_FM(int spin_index, cuReal3 Mnew, int*& cuModules, int numModules, cuReal3& Ha);

	//Demag_N
	__device__ cuBReal Get_EnergyChange_FM_DemagNCUDA(int spin_index, cuReal3 Mnew);

	//Demag
	__device__ cuBReal Get_EnergyChange_FM_DemagCUDA(int spin_index, cuReal3 Mnew);

	//StrayField_Mesh
	__device__ cuBReal Get_EnergyChange_FM_StrayField_MeshCUDA(int spin_index, cuReal3 Mnew);

	//Exch_6ngbr_Neu
	__device__ cuBReal Get_EnergyChange_FM_ExchangeCUDA(int spin_index, cuReal3 Mnew);

	//DMExchangeCUDA
	__device__ cuBReal Get_EnergyChange_FM_DMExchangeCUDA(int spin_index, cuReal3 Mnew);

	//iDMExchangeCUDA
	__device__ cuBReal Get_EnergyChange_FM_iDMExchangeCUDA(int spin_index, cuReal3 Mnew);

	//viDMExchangeCUDA
	__device__ cuBReal Get_EnergyChange_FM_viDMExchangeCUDA(int spin_index, cuReal3 Mnew);

	//SurfExchangeCUDA
	__device__ cuBReal Get_EnergyChange_FM_SurfExchangeCUDA(int spin_index, cuReal3 Mnew);

	//ZeemanCUDA
	__device__ cuBReal Get_EnergyChange_FM_ZeemanCUDA(int spin_index, cuReal3 Mnew, cuReal3& Ha);

	//MOpticalCUDA
	__device__ cuBReal Get_EnergyChange_FM_MOpticalCUDA(int spin_index, cuReal3 Mnew);

	//AnisotropyCUDA
	__device__ cuBReal Get_EnergyChange_FM_AnisotropyCUDA(int spin_index, cuReal3 Mnew);

	//AnisotropyCubiCUDA
	__device__ cuBReal Get_EnergyChange_FM_AnisotropyCubiCUDA(int spin_index, cuReal3 Mnew);

	//AnisotropyBiaxialCUDA
	__device__ cuBReal Get_EnergyChange_FM_AnisotropyBiaxialCUDA(int spin_index, cuReal3 Mnew);

	//AnisotropyTensorialCUDA
	__device__ cuBReal Get_EnergyChange_FM_AnisotropyTensorialCUDA(int spin_index, cuReal3 Mnew);

	//RoughnessCUDA
	__device__ cuBReal Get_EnergyChange_FM_RoughnessCUDA(int spin_index, cuReal3 Mnew);

	//MElasticCUDA
	__device__ cuBReal Get_EnergyChange_FM_MElasticCUDA(int spin_index, cuReal3 Mnew);

	////////////////////////////////////
	//
	// Antiferromagnetic

	//switch function which adds all assigned energy contributions in this mesh to calculate energy change from current spin to Mnew spin : return energy change as new - old
	//If Mnew is passed in as cuReal3(), then this function returns the current spin energy only - all functions in the switch statement implement this eventuality.
	__device__ cuReal2 Get_EnergyChange_AFM(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B, int*& cuModules, int numModules, cuReal3& Ha);

	//Demag_N
	__device__ cuReal2 Get_EnergyChange_AFM_DemagNCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B);

	//Demag
	__device__ cuReal2 Get_EnergyChange_AFM_DemagCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B);

	//StrayField_Mesh
	__device__ cuReal2 Get_EnergyChange_AFM_StrayField_MeshCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B);

	//Exch_6ngbr_Neu
	__device__ cuReal2 Get_EnergyChange_AFM_ExchangeCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B);

	//DMExchangeCUDA
	__device__ cuReal2 Get_EnergyChange_AFM_DMExchangeCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B);

	//iDMExchangeCUDA
	__device__ cuReal2 Get_EnergyChange_AFM_iDMExchangeCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B);

	//viDMExchangeCUDA
	__device__ cuReal2 Get_EnergyChange_AFM_viDMExchangeCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B);

	//SurfExchangeCUDA
	__device__ cuReal2 Get_EnergyChange_AFM_SurfExchangeCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B);

	//ZeemanCUDA
	__device__ cuReal2 Get_EnergyChange_AFM_ZeemanCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B, cuReal3& Ha);

	//MOpticalCUDA
	__device__ cuReal2 Get_EnergyChange_AFM_MOpticalCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B);

	//AnisotropyCUDA
	__device__ cuReal2 Get_EnergyChange_AFM_AnisotropyCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B);

	//AnisotropyCubiCUDA
	__device__ cuReal2 Get_EnergyChange_AFM_AnisotropyCubiCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B);

	//AnisotropyBiaxialCUDA
	__device__ cuReal2 Get_EnergyChange_AFM_AnisotropyBiaxialCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B);

	//AnisotropyTensorialCUDA
	__device__ cuReal2 Get_EnergyChange_AFM_AnisotropyTensorialCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B);

	//RoughnessCUDA
	__device__ cuReal2 Get_EnergyChange_AFM_RoughnessCUDA(int spin_index, cuReal3 Mnew_A, cuReal3 Mnew_B);
};

#endif