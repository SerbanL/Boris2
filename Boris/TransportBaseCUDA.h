#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "BorisCUDALib.h"

class MeshBaseCUDA;
class MeshBase;
class SuperMeshCUDA;
class SuperMesh;

class STransportCUDA;

class TransportBaseCUDA
{

	friend STransportCUDA;

private:

protected:

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshBaseCUDA* pMeshBaseCUDA;
	MeshBase* pMeshBase;

	//Same for the supermesh
	SuperMeshCUDA* pSMeshCUDA;
	SuperMesh* pSMesh;

	//spin transport solver type (see Transport_Defs.h)
	int stsolve;

private:

	//-------------------Auxiliary

	//set the stsolve indicator depending on current configuration
	//void Set_STSolveType(void);

	//------------------Others

	//set fixed potential cells in this mesh for given rectangle
	//bool SetFixedPotentialCells(cuRect rectangle, cuBReal potential);

	//void ClearFixedPotentialCells(void);

	virtual void Set_Linear_PotentialDrop(cuReal3 ground_electrode_center, cuBReal ground_potential, cuReal3 electrode_center, cuBReal electrode_potential) = 0;

	//prepare displayVEC ready for calculation of display quantity
	//bool PrepareDisplayVEC(DBL3 cellsize);

	//prepare displayVEC_VC ready for calculation of display quantity - used for charge current
	//bool PrepareDisplayVEC_VC(DBL3 cellsize);

	//check if dM_dt Calculation should be enabled
	//bool Need_dM_dt_Calculation(void);

	//check if the delsq_V_fixed VEC is needed
	//bool Need_delsq_V_fixed_Precalculation(void);

	//check if the delsq_S_fixed VEC is needed
	//bool Need_delsq_S_fixed_Precalculation(void);

	//-------------------Calculation Methods

	//calculate electrical conductivity with AMR present
	//void CalculateElectricalConductivity_AMR(void);

	//calculate electrical conductivity without AMR
	//void CalculateElectricalConductivity_NoAMR(void);

	//calculate electric field as the negative gradient of V
	virtual void CalculateElectricField(void) = 0;

	//calculate total current flowing into an electrode with given box - return ground_current and net_current values
	//cuBReal CalculateElectrodeCurrent(cuBox electrode_box);

	//Charge transport only : V

	//take a single iteration of the Transport solver in this Mesh (cannot solve fully in one go as there may be other meshes so need boundary conditions set externally). Use SOR. Return relaxation value (tends to zero).
	virtual void IterateChargeSolver_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value) = 0;

	//Calculation Methods used by Spin Current Solver only

	//before iterating the spin solver (charge part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
	virtual void PrimeSpinSolver_Charge(void) = 0;

	//take a single iteration of the charge transport solver (within the spin current solver) in this Mesh (cannot solve fully in one go as there may be other meshes so need boundary conditions set externally). Use SOR. Return relaxation value (tends to zero).
	//use_NNeu flag indicates if we need to use the homogeneous or non-homogeneous Neumann boundary conditions versions
	virtual void IterateSpinSolver_Charge_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value, bool use_NNeu) = 0;

	//before iterating the spin solver (spin part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
	virtual void PrimeSpinSolver_Spin(void) = 0;

	//solve for spin accumulation using Poisson equation for delsq_S. Use SOR.
	//use_NNeu flag indicates if we need to use the homogeneous or non-homogeneous Neumann boundary conditions versions
	virtual void IterateSpinSolver_Spin_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value, bool use_NNeu) = 0;

public:

	TransportBaseCUDA() {}
	virtual ~TransportBaseCUDA() {}

	//-------------------

	int Get_STSolveType(void) { return stsolve; }

	//-------------------Public calculation Methods

	//calculate the field resulting from a spin accumulation (spin current solver enabled) so a spin accumulation torque results when using the LLG or LLB equations
	virtual void CalculateSAField(void) = 0;

	//Calculate the field resulting from interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set)
	//void CalculateSAInterfaceField(TransportCUDA* ptrans_sec, CMBNDInfoCUDA& contactCUDA, bool primary_top);

	//Calculate the interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set), accumulating result in displayVEC
	//void CalculateDisplaySAInterfaceTorque(TransportCUDA* ptrans_sec, CMBNDInfoCUDA& contactCUDA, bool primary_top);

	//calculate elC VEC using AMR and temperature information
	virtual void CalculateElectricalConductivity(bool force_recalculate = false) = 0;

	//-------------------Display Calculation Methods

	//return x, y, or z component of spin current (component = 0, 1, or 2)
	cu_obj<cuVEC<cuReal3>>& GetSpinCurrent(int component);

	cu_obj<cuVEC_VC<cuReal3>>& GetChargeCurrent(void);

	//return spin torque computed from spin accumulation
	cu_obj<cuVEC<cuReal3>>& GetSpinTorque(void);
};

#endif

#endif