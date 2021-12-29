#pragma once

#include "BorisLib.h"
#include "Boris_Enums_Defs.h"
#include "ErrorHandler.h"

class MeshBase;
class SuperMesh;
class STransport;

#ifdef MODULE_COMPILATION_TRANSPORT

#include "Transport_Defs.h"

class TransportBase {

	friend STransport;

protected:

	//pointer to mesh object holding this effective field module
	MeshBase* pMeshBase;

	SuperMesh* pSMesh;

	//spin transport solver type (see Transport_Defs.h)
	STSOLVE_ stsolve;

protected:

	//-------------------Auxiliary

	//-------------------Calculation Methods

	//calculate electric field as the negative gradient of V
	virtual void CalculateElectricField(void) = 0;

	//Charge transport only : V

	//take a single iteration of the charge transport solver in this Mesh (cannot solve fully in one go as there may be other meshes so need boundary conditions set externally). Use SOR. 
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	virtual DBL2 IterateChargeSolver_SOR(double damping) = 0;

	//call-back method for Poisson equation to evaluate RHS
	virtual double Evaluate_ChargeSolver_delsqV_RHS(int idx) const = 0;

	//Calculation Methods used by Spin Current Solver only

	//take a single iteration of the charge transport solver (within the spin current solver) in this Mesh (cannot solve fully in one go as there may be other meshes so need boundary conditions set externally). Use SOR. 
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	virtual DBL2 IterateSpinSolver_Charge_SOR(double damping) = 0;

	//before iterating the spin solver (charge part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
	virtual void PrimeSpinSolver_Charge(void) = 0;

	//call-back method for Poisson equation to evaluate RHS
	virtual double Evaluate_SpinSolver_delsqV_RHS(int idx) const = 0;

	//solve for spin accumulation using Poisson equation for delsq_S. Use SOR. 
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	virtual DBL2 IterateSpinSolver_Spin_SOR(double damping) = 0;

	//before iterating the spin solver (spin part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
	virtual void PrimeSpinSolver_Spin(void) = 0;

	//call-back method for Poisson equation for S
	virtual DBL3 Evaluate_SpinSolver_delsqS_RHS(int idx) const = 0;

	//Other Calculation Methods

	//calculate total current flowing into an electrode with given box - return ground_current and net_current values
	virtual double CalculateElectrodeCurrent(Rect &electrode_rect) = 0;

	//------------------Others

	//set fixed potential cells in this mesh for given rectangle
	virtual bool SetFixedPotentialCells(Rect rectangle, double potential) = 0;

	virtual void ClearFixedPotentialCells(void) = 0;

	virtual void Set_Linear_PotentialDrop(DBL3 ground_electrode_center, double ground_potential, DBL3 electrode_center, double electrode_potential) = 0;

	//prepare displayVEC ready for calculation of display quantity - used for spin current
	virtual bool PrepareDisplayVEC(DBL3 cellsize) = 0;

	//prepare displayVEC_VC ready for calculation of display quantity - used for charge current
	virtual bool PrepareDisplayVEC_VC(DBL3 cellsize) = 0;

	//check if dM_dt Calculation should be enabled
	virtual bool Need_dM_dt_Calculation(void) = 0;

	//check if the delsq_V_fixed VEC is needed
	virtual bool Need_delsq_V_fixed_Precalculation(void) = 0;

	//check if the delsq_S_fixed VEC is needed
	virtual bool Need_delsq_S_fixed_Precalculation(void) = 0;

public:

	TransportBase(void)
	{
	}

	virtual ~TransportBase() {}

	//-------------------

	STSOLVE_ Get_STSolveType(void) { return stsolve; }

	//-------------------CMBND computation methods
	
	//CMBND values set based on continuity of a potential and flux

	//Charge transport only : V

	//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface
	virtual double afunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const = 0;
	virtual double afunc_V_pri(int cell1_idx, int cell2_idx, DBL3 shift) const = 0;

	virtual double bfunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const = 0;
	virtual double bfunc_V_pri(int cell1_idx, int cell2_idx) const = 0;

	//second order differential of V at cells either side of the boundary; delsq V = -grad V * grad elC / elC
	virtual double diff2_V_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const = 0;
	virtual double diff2_V_pri(int cell1_idx, DBL3 shift) const = 0;

	//Spin transport : V and S
	//CMBND for V
	virtual double afunc_st_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const = 0;
	virtual double afunc_st_V_pri(int cell1_idx, int cell2_idx, DBL3 shift) const = 0;

	virtual double bfunc_st_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const = 0;
	virtual double bfunc_st_V_pri(int cell1_idx, int cell2_idx) const = 0;

	virtual double diff2_st_V_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const = 0;
	virtual double diff2_st_V_pri(int cell1_idx, DBL3 shift) const = 0;

	//CMBND for S
	virtual DBL3 afunc_st_S_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const = 0;
	virtual DBL3 afunc_st_S_pri(int cell1_idx, int cell2_idx, DBL3 shift) const = 0;

	virtual double bfunc_st_S_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const = 0;
	virtual double bfunc_st_S_pri(int cell1_idx, int cell2_idx) const = 0;

	virtual DBL3 diff2_st_S_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const = 0;
	virtual DBL3 diff2_st_S_pri(int cell1_idx, DBL3 shift) const = 0;

	//multiply spin accumulation by these to obtain spin potential, i.e. Vs = (De / elC) * (e/muB) * S, evaluated at the boundary
	virtual double cfunc_sec(DBL3 relpos, DBL3 stencil) const = 0;
	virtual double cfunc_pri(int cell_idx) const = 0;

	//-------------------Properties

	virtual bool GInterface_Enabled(void) = 0;

	//-------------------Others

	//called by MoveMesh method in this mesh - move relevant transport quantities
	virtual void MoveMesh_Transport(double x_shift) = 0;

	//-------------------Public calculation Methods

	//calculate the field resulting from a spin accumulation (spin current solver enabled) so a spin accumulation torque results when using the LLG or LLB equations
	virtual void CalculateSAField(void) = 0;

	//calculate elC VEC using AMR and temperature information
	virtual void CalculateElectricalConductivity(bool force_recalculate = false) = 0;

	//-------------------Display Calculation / Get Methods

	//return x, y, or z component of spin current (component = 0, 1, or 2)
	virtual VEC<DBL3>& GetSpinCurrent(int component) = 0;
	//get value, assuming quantity already calculated with a call to the above method
	virtual DBL3 GetCalculatedSpinCurrentValue(const DBL3& rel_pos) = 0;

	virtual DBL3 GetAverageSpinCurrent(int component, Rect rectangle = Rect(), std::vector<MeshShape> shapes = {}) = 0;

	//calculate charge current density over the mesh : applies to both charge-only transport and spin transport solvers (if check used inside this method)
	virtual VEC_VC<DBL3>& GetChargeCurrent(void) = 0;
	//get value, assuming quantity already calculated with a call to the above method
	virtual DBL3 GetCalculatedChargeCurrentValue(const DBL3& rel_pos) = 0;

	virtual DBL3 GetAverageChargeCurrent(Rect rectangle = Rect(), std::vector<MeshShape> shapes = {}) = 0;

	//return spin torque computed from spin accumulation
	virtual VEC<DBL3>& GetSpinTorque(void) = 0;

	//get value, assuming quantity already calculated with a call to the above method
	virtual DBL3 GetCalculatedSpinTorqueValue(const DBL3& rel_pos) = 0;

	//Get average bulk spin torque in given rectangle - calculate it first
	virtual DBL3 GetAverageSpinTorque(Rect rectangle = Rect(), std::vector<MeshShape> shapes = {}) = 0;

	//Get average interfacial spin torque in given rectangle - must have been calculated already in displayVEC through the supermesh transport module
	virtual DBL3 GetAverageInterfacialSpinTorque(Rect rectangle = Rect(), std::vector<MeshShape> shapes = {}) = 0;
};

#endif