#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TMR

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

#include "TransportBaseCUDA.h"
#include "Transport_Defs.h"

class MeshCUDA;
class Mesh;
class ManagedMeshCUDA;
class ManagedAtom_MeshCUDA;
class SuperMesh;

class TMR;
class STransportCUDA;

class TMRCUDA :
	public ModulesCUDA,
	public TransportBaseCUDA
{
	friend TMR;
	friend STransportCUDA;

	friend TransportCUDA_Spin_V_Funcs;
	friend TransportCUDA_Spin_S_Funcs;

private:

	//pointer to cpu version of SurfExchange
	TMR* pTMR;

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

	//pointer to cpu version of MeshCUDA
	Mesh* pMesh;

	//Contacting top and bottom mehes used to calculate TMR
	//All magnetic meshes used (Ferromagnetic and Atomistic)
	cu_arr<ManagedMeshCUDA> pMeshFM_Bot;
	cu_arr<ManagedMeshCUDA> pMeshFM_Top;

	cu_arr<ManagedAtom_MeshCUDA> pMeshAtom_Bot;
	cu_arr<ManagedAtom_MeshCUDA> pMeshAtom_Top;

	//bias dependence equations for TMR RA values
	TEquationCUDA<cuBReal, cuBReal, cuBReal> RAtmr_p_equation, RAtmr_ap_equation;

private:

	//-------------------Calculation Methods
	
	//calculate electrical conductivity based on TMR formula
	void CalculateElectricalConductivity_TMR(TMR_ TMR_type);

	//calculate electric field as the negative gradient of V
	void CalculateElectricField(bool open_potential = false);

	//Charge transport only : V

	//take a single iteration of the Transport solver in this Mesh (cannot solve fully in one go as there may be other meshes so need boundary conditions set externally). Use SOR. Return relaxation value (tends to zero).
	void IterateChargeSolver_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value);

	//Calculation Methods used by Spin Current Solver only

	//before iterating the spin solver (charge part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
	void PrimeSpinSolver_Charge(void);

	//take a single iteration of the charge transport solver (within the spin current solver) in this Mesh (cannot solve fully in one go as there may be other meshes so need boundary conditions set externally). Use SOR. Return relaxation value (tends to zero).
	//use_NNeu flag indicates if we need to use the homogeneous or non-homogeneous Neumann boundary conditions versions
	void IterateSpinSolver_Charge_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value, bool use_NNeu);

	//before iterating the spin solver (spin part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
	void PrimeSpinSolver_Spin(void);

	//solve for spin accumulation using Poisson equation for delsq_S. Use SOR.
	//use_NNeu flag indicates if we need to use the homogeneous or non-homogeneous Neumann boundary conditions versions
	void IterateSpinSolver_Spin_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value, bool use_NNeu);

	//------------------Others

	//check if dM_dt Calculation should be enabled
	bool Need_dM_dt_Calculation(void);

	//check if the delsq_V_fixed VEC is needed
	bool Need_delsq_V_fixed_Precalculation(void);

	//check if the delsq_S_fixed VEC is needed
	bool Need_delsq_S_fixed_Precalculation(void);

public:

	TMRCUDA(TMR* pTMR_);
	~TMRCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	void UpdateField(void);

	//-------------------Output Data

	//calculate resistance in given rectangle (relative to mesh)
	cuBReal GetResistance(cuRect rect);

	//-------------------Public calculation Methods

	//calculate the field resulting from a spin accumulation (spin current solver enabled) so a spin accumulation torque results when using the LLG or LLB equations
	void CalculateSAField(void);

	//Calculate the field resulting from interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set)
	void CalculateSAInterfaceField(TransportBaseCUDA* ptrans_sec, CMBNDInfoCUDA& contactCUDA, bool primary_top);

	//Calculate the interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set), accumulating result in displayVEC
	void CalculateDisplaySAInterfaceTorque(TransportBaseCUDA* ptrans_sec, CMBNDInfoCUDA& contactCUDA, bool primary_top);

	//calculate elC VEC using AMR and temperature information
	void CalculateElectricalConductivity(bool force_recalculate = false);

	//-------------------Display Calculation Methods

	//return x, y, or z component of spin current (component = 0, 1, or 2)
	cu_obj<cuVEC<cuReal3>>& GetSpinCurrent(int component);

	cu_obj<cuVEC_VC<cuReal3>>& GetChargeCurrent(void);

	//return spin torque computed from spin accumulation
	cu_obj<cuVEC<cuReal3>>& GetSpinTorque(void);
	
	//-------------------Other Settings

	BError SetBiasEquation_Parallel(const std::vector< std::vector<EqComp::FSPEC> >& fspec);
	BError SetBiasEquation_AntiParallel(const std::vector< std::vector<EqComp::FSPEC> >& fspec);
};

#else

class TMRCUDA
{
};

#endif

#endif


