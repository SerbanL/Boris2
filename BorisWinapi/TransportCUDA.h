#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_TRANSPORT

#include "BorisCUDALib.h"
#include "ModulesCUDA.h"

#include "TransportCUDA_Poisson.h"
#include "TransportCUDA_Poisson_Spin_V.h"
#include "TransportCUDA_Poisson_Spin_S.h"

#include "Transport_Defs.h"

class MeshCUDA;
class Mesh;
class SuperMeshCUDA;
class SuperMesh;

class Transport;
class STransportCUDA;

class TransportCUDA :
	public ModulesCUDA
{

	friend Transport;
	friend STransportCUDA;
	friend TransportCUDA_Spin_V_Funcs;
	friend TransportCUDA_Spin_S_Funcs;

private:

	//spin transport solver type (see Transport_Defs.h)
	int stsolve;

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshCUDA* pMeshCUDA;

	//pointer to cpu version of MeshCUDA
	Mesh* pMesh;

	//Same for the supermesgh
	SuperMeshCUDA* pSMeshCUDA;
	SuperMesh* pSMesh;

	//pointer to cpu version of this module
	Transport* pTransport;

	//TransportCUDA_V_Funcs holds the Poisson_RHS method used in solving the Poisson equation.
	//Pass the managed TransportCUDA_V_Funcs object to IteratePoisson_SOR method in V -  see IterateTransportSolver_SOR method here. IteratePoisson_SOR will then use the Poisson_RHS method defined in TransportCUDA_V_Funcs.
	cu_obj<TransportCUDA_V_Funcs> poisson_V;

	//similar for the full spin current solver
	cu_obj<TransportCUDA_Spin_V_Funcs> poisson_Spin_V;
	cu_obj<TransportCUDA_Spin_S_Funcs> poisson_Spin_S;

	//dM_dt VEC when we need to do vector calculus operations on it
	cu_obj<cuVEC_VC<cuReal3>> dM_dt;

	//for Poisson equations for V and S some values are fixed during relaxation, so pre-calculate them and store here to re-use.
	cu_obj <cuVEC<cuBReal>> delsq_V_fixed;
	cu_obj <cuVEC<cuReal3>> delsq_S_fixed;

	//used to compute spin current components and spin torque for display purposes - only calculated and memory allocated when asked to do so.
	cu_obj<cuVEC<cuReal3>> displayVEC;

	//used to compute charge current for display purposes (and in some cases we need to run vector calculus operations on it)
	cu_obj<cuVEC_VC<cuReal3>> displayVEC_VC;

private:

	//-------------------Auxiliary

	//set the stsolve indicator depending on current configuration
	void Set_STSolveType(void);

	//------------------Others

	//set fixed potential cells in this mesh for given rectangle
	bool SetFixedPotentialCells(cuRect rectangle, cuBReal potential);

	void ClearFixedPotentialCells(void);

	//prepare displayVEC ready for calculation of display quantity
	bool PrepareDisplayVEC(DBL3 cellsize);

	//prepare displayVEC_VC ready for calculation of display quantity - used for charge current
	bool PrepareDisplayVEC_VC(DBL3 cellsize);

	//check if dM_dt Calculation should be enabled
	bool Need_dM_dt_Calculation(void);

	//check if the delsq_V_fixed VEC is needed
	bool Need_delsq_V_fixed_Precalculation(void);

	//check if the delsq_S_fixed VEC is needed
	bool Need_delsq_S_fixed_Precalculation(void);

	//-------------------Calculation Methods

	//calculate electrical conductivity with AMR present
	void CalculateElectricalConductivity_AMR(void);

	//calculate electrical conductivity without AMR
	void CalculateElectricalConductivity_NoAMR(void);

	//calculate electric field as the negative gradient of V
	void CalculateElectricField(void);

	//calculate total current flowing into an electrode with given box - return ground_current and net_current values
	cuBReal CalculateElectrodeCurrent(cuBox electrode_box);

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

public:

	TransportCUDA(Mesh* pMesh_, SuperMesh* pSMesh_, Transport* pTransport_);
	~TransportCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);

	//-------------------Public calculation Methods

	//calculate the field resulting from a spin accumulation (spin current solver enabled) so a spin accumulation torque results when using the LLG or LLB equations
	void CalculateSAField(void);

	//Calculate the field resulting from interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set)
	void CalculateSAInterfaceField(TransportCUDA* ptrans_sec, CMBNDInfoCUDA& contactCUDA, bool primary_top);

	//Calculate the interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set), accumulating result in displayVEC
	void CalculateDisplaySAInterfaceTorque(TransportCUDA* ptrans_sec, CMBNDInfoCUDA& contactCUDA, bool primary_top);

	//calculate elC VEC using AMR and temperature information
	void CalculateElectricalConductivity(bool force_recalculate = false);

	//-------------------Display Calculation Methods

	//return x, y, or z component of spin current (component = 0, 1, or 2)
	cu_obj<cuVEC<cuReal3>>& GetSpinCurrent(int component);

	cu_obj<cuVEC_VC<cuReal3>>& GetChargeCurrent(void);

	//return spin torque computed from spin accumulation
	cu_obj<cuVEC<cuReal3>>& GetSpinTorque(void);
};

#else

class TransportCUDA
{
};

#endif

#endif