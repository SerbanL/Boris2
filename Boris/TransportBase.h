#pragma once

#include "BorisLib.h"
#include "Boris_Enums_Defs.h"
#include "ErrorHandler.h"

#ifdef MODULE_COMPILATION_TRANSPORT

class MeshBase;
class SuperMesh;
class STransport;

class TransportBaseCUDA;

#include "Transport_Defs.h"

#if COMPILECUDA == 1
#include "TransportBaseCUDA.h"
#endif

class TransportBase 
{
	friend class Atom_Transport;
	friend class Transport;
	friend class STransport;

protected:

	//pointer to mesh object holding this effective field module
	MeshBase* pMeshBase;
	SuperMesh* pSMesh;

	//pointer to STransport supermesh module; this is set from STransport in its UpdateConfiguration method
	//need it to access some STransport properties
	STransport* pSTrans = nullptr;

#if COMPILECUDA == 1
	friend TransportBaseCUDA;
	TransportBaseCUDA* pTransportBaseCUDA = nullptr;
#endif

	//spin transport solver type (see Transport_Defs.h)
	STSOLVE_ stsolve;

	//used to compute spin current components and spin torque for display purposes - only calculated and memory allocated when asked to do so.
	VEC<DBL3> displayVEC;

	//used to compute charge current for display purposes (and in some cases we need to run vector calculus operations on it)
	VEC_VC<DBL3> displayVEC_VC;

	//dM_dt VEC when we need to do vector calculus operations on it
	//used inside a Atom_Transport const function, so needs to be mutable, else the compiler is unhappy!
	mutable VEC_VC<DBL3> dM_dt;

	//for Poisson equations for V and S some values are fixed during relaxation, so pre-calculate them and store here to re-use.
	mutable VEC<double> delsq_V_fixed;
	mutable VEC<DBL3> delsq_S_fixed;

	//for meshes with thermoelectric effect enabled count net current (A) coming out of the mesh at cmbnd cells, if open potential enabled.
	double mesh_thermoelectric_net_current = 0.0;

	//does current mesh have thermoelectric effect enabled? (i.e. must have heat module and Sc not zero)
	bool is_thermoelectric_mesh = false;

	//user test equation for TAMR. If not set use default formula.
	TEquation<double, double, double, double> TAMR_conductivity_equation;

protected:

	//-------------------Auxiliary

	//set the stsolve indicator depending on current configuration
	void Set_STSolveType(void);

	//prepare displayVEC ready for calculation of display quantity - used for spin current
	bool PrepareDisplayVEC(DBL3 cellsize);

	//prepare displayVEC_VC ready for calculation of display quantity - used for charge current
	bool PrepareDisplayVEC_VC(DBL3 cellsize);

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

	//before iterating the spin solver (charge part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
	virtual void PrimeSpinSolver_Charge(void) = 0;

	//take a single iteration of the charge transport solver (within the spin current solver) in this Mesh (cannot solve fully in one go as there may be other meshes so need boundary conditions set externally). Use SOR. 
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	virtual DBL2 IterateSpinSolver_Charge_SOR(double damping) = 0;

	//call-back method for Poisson equation to evaluate RHS
	virtual double Evaluate_SpinSolver_delsqV_RHS(int idx) const = 0;

	//before iterating the spin solver (spin part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
	virtual void PrimeSpinSolver_Spin(void) = 0;

	//solve for spin accumulation using Poisson equation for delsq_S. Use SOR. 
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	virtual DBL2 IterateSpinSolver_Spin_SOR(double damping) = 0;

	//call-back method for Poisson equation for S
	virtual DBL3 Evaluate_SpinSolver_delsqS_RHS(int idx) const = 0;

	//Non-homogeneous Neumann boundary condition for V' - call-back method for Poisson equation for V
	virtual DBL3 NHNeumann_Vdiff(int idx) const = 0;

	//Non-homogeneous Neumann boundary condition for S' - call-back method for Poisson equation for S
	virtual DBL33 NHNeumann_Sdiff(int idx) const = 0;

	//Other Calculation Methods

	//calculate total current flowing into an electrode with given box - return ground_current and net_current values
	double CalculateElectrodeCurrent(Rect &electrode_rect);

	//------------------Others

	//set fixed potential cells in this mesh for given rectangle
	bool SetFixedPotentialCells(Rect rectangle, double potential);

	void ClearFixedPotentialCells(void);

	//check if dM_dt Calculation should be enabled
	virtual bool Need_dM_dt_Calculation(void) = 0;

	//check if the delsq_V_fixed VEC is needed
	virtual bool Need_delsq_V_fixed_Precalculation(void) = 0;

	//check if the delsq_S_fixed VEC is needed
	virtual bool Need_delsq_S_fixed_Precalculation(void) = 0;

public:

	TransportBase(void) :
		TAMR_conductivity_equation({ "TAMR", "d1", "d2", "d3" })
	{}

	TransportBase(MeshBase* pMeshBase_);

	virtual ~TransportBase() {}

	//-------------------Properties

	STSOLVE_ Get_STSolveType(void) { return stsolve; }

	bool GInterface_Enabled(void);

	bool iSHA_nonzero(void);
	bool SHA_nonzero(void);

	bool IsOpenPotential(void);

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

	//-------------------TAMR

	BError Set_TAMR_Conductivity_Equation(std::string equation_string);

	//-------------------Others

	//called by MoveMesh method in this mesh - move relevant transport quantities
	void MoveMesh_Transport(double x_shift);

	//-------------------Public calculation Methods

	//calculate the field resulting from a spin accumulation (spin current solver enabled) so a spin accumulation torque results when using the LLG or LLB equations
	virtual void CalculateSAField(void) = 0;

	//Calculate the field resulting from interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set)
	virtual void CalculateSAInterfaceField(TransportBase* ptrans_sec, CMBNDInfo& contact) = 0;

	//Calculate the interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set), accumulating result in displayVEC
	virtual void CalculateDisplaySAInterfaceTorque(TransportBase* ptrans_sec, CMBNDInfo& contact) = 0;

	//calculate elC VEC using AMR and temperature information
	virtual void CalculateElectricalConductivity(bool force_recalculate = false) = 0;

	//-------------------Display Calculation / Get Methods

	//return x, y, or z component of spin current (component = 0, 1, or 2)
	virtual VEC<DBL3>& GetSpinCurrent(int component) = 0;
	
	//get value, assuming quantity already calculated with a call to the above method
	DBL3 GetCalculatedSpinCurrentValue(const DBL3& rel_pos) { if (displayVEC.linear_size()) return displayVEC[rel_pos]; else return DBL3(); }

	DBL3 GetAverageSpinCurrent(int component, Rect rectangle = Rect(), std::vector<MeshShape> shapes = {});

	//calculate charge current density over the mesh : applies to both charge-only transport and spin transport solvers (if check used inside this method)
	virtual VEC_VC<DBL3>& GetChargeCurrent(void) = 0;

	//get value, assuming quantity already calculated with a call to the above method
	DBL3 GetCalculatedChargeCurrentValue(const DBL3& rel_pos) { if (displayVEC_VC.linear_size()) return displayVEC_VC[rel_pos]; else return DBL3(); }
	DBL3 GetAverageChargeCurrent(Rect rectangle = Rect(), std::vector<MeshShape> shapes = {});

	//return spin torque computed from spin accumulation
	virtual VEC<DBL3>& GetSpinTorque(void) = 0;

	//get value, assuming quantity already calculated with a call to the above method
	DBL3 GetCalculatedSpinTorqueValue(const DBL3& rel_pos) { if (displayVEC.linear_size()) return displayVEC[rel_pos]; else return DBL3(); }

	//Get average bulk spin torque in given rectangle - calculate it first
	DBL3 GetAverageSpinTorque(Rect rectangle = Rect(), std::vector<MeshShape> shapes = {});

	//Get average interfacial spin torque in given rectangle - must have been calculated already in displayVEC through the supermesh transport module
	DBL3 GetAverageInterfacialSpinTorque(Rect rectangle = Rect(), std::vector<MeshShape> shapes = {});

#if COMPILECUDA == 1
	cu_obj<cuVEC_VC<cuReal3>>& GetChargeCurrentCUDA(void);
	cu_obj<cuVEC<cuReal3>>& GetSpinCurrentCUDA(int component);
	cu_obj<cuVEC<cuReal3>>& GetSpinTorqueCUDA(void);
#endif

	//------------------Others

	void Set_Linear_PotentialDrop(Rect contact1, double potential1, Rect contact2, double potential2, DBL2 degeneracy);
};

#endif