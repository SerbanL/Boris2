#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TRANSPORT

#include "BorisCUDALib.h"

#include "TransportCUDA_Poisson.h"
#include "TransportCUDA_Poisson_Spin_V.h"
#include "TransportCUDA_Poisson_Spin_S.h"

class SuperMeshCUDA;
class SuperMesh;

class MeshBaseCUDA;

class TransportBase;

class TransportBaseCUDA
{
	friend class TransportBase;

	friend TransportCUDA_Spin_V_Funcs;
	friend TransportCUDA_Spin_S_Funcs;

	friend class STransportCUDA;
	friend class Atom_TransportCUDA;
	friend class TransportCUDA;

private:

protected:

	//pointer to cpu version of this Transport Base
	TransportBase* pTransportBase;

	//pointer to CUDA version of mesh object holding the effective field module holding this CUDA module
	MeshBaseCUDA* pMeshBaseCUDA;

	//Same for the supermesh
	SuperMesh* pSMesh;
	SuperMeshCUDA* pSMeshCUDA;

	//spin transport solver type (see Transport_Defs.h)
	int stsolve;

	//auxiliary for reductions
	cu_obj<cuBReal> auxReal;

	//used to compute spin current components and spin torque for display purposes - only calculated and memory allocated when asked to do so.
	cu_obj<cuVEC<cuReal3>> displayVEC;

	//used to compute charge current for display purposes (and in some cases we need to run vector calculus operations on it)
	cu_obj<cuVEC_VC<cuReal3>> displayVEC_VC;

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

protected:

	//-------------------Auxiliary

	//set the stsolve indicator depending on current configuration
	void Set_STSolveType(void);

	//prepare displayVEC ready for calculation of display quantity
	bool PrepareDisplayVEC(DBL3 cellsize);

	//prepare displayVEC_VC ready for calculation of display quantity - used for charge current
	bool PrepareDisplayVEC_VC(DBL3 cellsize);

	//-------------------Calculation Methods

	//calculate electric field as the negative gradient of V
	virtual void CalculateElectricField(void) = 0;

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

	//Other Calculation Methods

	//calculate total current flowing into an electrode with given box - return ground_current and net_current values
	cuBReal CalculateElectrodeCurrent(cuBox electrode_box, cuINT3 sign);

	//------------------Others

	//set fixed potential cells in this mesh for given rectangle
	bool SetFixedPotentialCells(cuRect rectangle, cuBReal potential);

	void ClearFixedPotentialCells(void);

	void Set_Linear_PotentialDrop(cuRect contact1, cuBReal potential1, cuRect contact2, cuBReal potential2, cuReal2 degeneracy);

	//check if dM_dt Calculation should be enabled
	virtual bool Need_dM_dt_Calculation(void) = 0;

	//check if the delsq_V_fixed VEC is needed
	virtual bool Need_delsq_V_fixed_Precalculation(void) = 0;

	//check if the delsq_S_fixed VEC is needed
	virtual bool Need_delsq_S_fixed_Precalculation(void) = 0;

public:

	TransportBaseCUDA(TransportBase* pTransportBase_);
	virtual ~TransportBaseCUDA();

	//-------------------Properties

	int Get_STSolveType(void) { return stsolve; }

	bool GInterface_Enabled(void);

	bool iSHA_nonzero(void);
	bool SHA_nonzero(void);

	//-------------------Public calculation Methods

	//calculate the field resulting from a spin accumulation (spin current solver enabled) so a spin accumulation torque results when using the LLG or LLB equations
	virtual void CalculateSAField(void) = 0;

	//Calculate the field resulting from interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set)
	virtual void CalculateSAInterfaceField(TransportBaseCUDA* ptrans_sec, CMBNDInfoCUDA& contactCUDA, bool primary_top) = 0;

	//Calculate the interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set), accumulating result in displayVEC
	virtual void CalculateDisplaySAInterfaceTorque(TransportBaseCUDA* ptrans_sec, CMBNDInfoCUDA& contactCUDA, bool primary_top) = 0;

	//calculate elC VEC using AMR and temperature information
	virtual void CalculateElectricalConductivity(bool force_recalculate = false) = 0;

	//-------------------Display Calculation Methods

	//return x, y, or z component of spin current (component = 0, 1, or 2)
	virtual cu_obj<cuVEC<cuReal3>>& GetSpinCurrent(int component) = 0;

	virtual cu_obj<cuVEC_VC<cuReal3>>& GetChargeCurrent(void) = 0;

	//return spin torque computed from spin accumulation
	virtual cu_obj<cuVEC<cuReal3>>& GetSpinTorque(void) = 0;
};

#endif

#endif
