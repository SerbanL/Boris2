#pragma once

#include "BorisLib.h"
#include "Modules.h"



class Mesh;
class SuperMesh;
class STransport;

#ifdef MODULE_COMPILATION_TRANSPORT

#include "Transport_Defs.h"

#if COMPILECUDA == 1
#include "TransportCUDA.h"
#endif

class Transport :
	public Modules,
	public ProgramState<Transport, std::tuple<>,	std::tuple<>>
{
	
	friend STransport;

#if COMPILECUDA == 1
	friend TransportCUDA;

	//dynamic cast pModuleCUDA to this after new operator, so we don't have to cast it again later
	//do not check this for nullptr, but check pModuleCUDA instead; this is only to replace repeated use of dynamic_cast
	//thus if you check pModuleCUDA first, you can safely use pTransportCUDA since this is cast from pModuleCUDA only after pModuleCUDA allocated
	TransportCUDA* pTransportCUDA = nullptr;
#endif

private:

	//spin transport solver type (see Transport_Defs.h)
	STSOLVE_ stsolve;

	//pointer to mesh object holding this effective field module
	Mesh* pMesh;

	SuperMesh* pSMesh;

	//used to compute spin current components and spin torque for display purposes - only calculated and memory allocated when asked to do so.
	VEC<DBL3> displayVEC;

	//used to compute charge current for display purposes (and in some cases we need to run vector calculus operations on it)
	VEC_VC<DBL3> displayVEC_VC;

	//dM_dt VEC when we need to do vector calculus operations on it
	//used inside a Transport const function, so needs to be mutable, else the compiler is unhappy!
	mutable VEC_VC<DBL3> dM_dt;

	//for Poisson equations for V and S some values are fixed during relaxation, so pre-calculate them and store here to re-use.
	mutable VEC<double> delsq_V_fixed;
	mutable VEC<DBL3> delsq_S_fixed;

private:

	//-------------------Auxiliary

	//set the stsolve indicator depending on current configuration
	void Set_STSolveType(void);

	//-------------------Calculation Methods

	//calculate electric field as the negative gradient of V
	void CalculateElectricField(void);

	//Charge transport only : V

	//take a single iteration of the charge transport solver in this Mesh (cannot solve fully in one go as there may be other meshes so need boundary conditions set externally). Use SOR. 
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	DBL2 IterateChargeSolver_SOR(double damping);

	//call-back method for Poisson equation to evaluate RHS
	double Evaluate_ChargeSolver_delsqV_RHS(int idx) const;

	//Calculation Methods used by Spin Current Solver only

	//take a single iteration of the charge transport solver (within the spin current solver) in this Mesh (cannot solve fully in one go as there may be other meshes so need boundary conditions set externally). Use SOR. 
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	DBL2 IterateSpinSolver_Charge_SOR(double damping);

	//before iterating the spin solver (charge part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
	void PrimeSpinSolver_Charge(void);

	//call-back method for Poisson equation to evaluate RHS
	double Evaluate_SpinSolver_delsqV_RHS(int idx) const;

	//solve for spin accumulation using Poisson equation for delsq_S. Use SOR. 
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	DBL2 IterateSpinSolver_Spin_SOR(double damping);

	//before iterating the spin solver (spin part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
	void PrimeSpinSolver_Spin(void);

	//call-back method for Poisson equation for S
	DBL3 Evaluate_SpinSolver_delsqS_RHS(int idx) const;

	//Non-homogeneous Neumann boundary condition for V' - call-back method for Poisson equation for V
	DBL3 NHNeumann_Vdiff(int idx) const;

	//Non-homogeneous Neumann boundary condition for S' - call-back method for Poisson equation for S
	DBL33 NHNeumann_Sdiff(int idx) const;

	//Other Calculation Methods

	//calculate total current flowing into an electrode with given box - return ground_current and net_current values
	double CalculateElectrodeCurrent(Rect &electrode_rect);

	//------------------Others

	//set fixed potential cells in this mesh for given rectangle
	bool SetFixedPotentialCells(Rect rectangle, double potential);

	void ClearFixedPotentialCells(void);

	//prepare displayVEC ready for calculation of display quantity - used for spin current
	bool PrepareDisplayVEC(DBL3 cellsize);

	//prepare displayVEC_VC ready for calculation of display quantity - used for charge current
	bool PrepareDisplayVEC_VC(DBL3 cellsize);

	//check if dM_dt Calculation should be enabled
	bool Need_dM_dt_Calculation(void);

	//check if the delsq_V_fixed VEC is needed
	bool Need_delsq_V_fixed_Precalculation(void);

	//check if the delsq_S_fixed VEC is needed
	bool Need_delsq_S_fixed_Precalculation(void);

public:

	Transport(Mesh *pMesh_);
	~Transport();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------CMBND computation methods

	//CMBND values set based on continuity of a potential and flux

	//Charge transport only : V

	//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface
	double afunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const;
	double afunc_V_pri(int cell1_idx, int cell2_idx, DBL3 shift) const;

	double bfunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const;
	double bfunc_V_pri(int cell1_idx, int cell2_idx) const;

	//second order differential of V at cells either side of the boundary; delsq V = -grad V * grad elC / elC
	double diff2_V_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const;
	double diff2_V_pri(int cell1_idx, DBL3 shift) const;

	//Spin transport : V and S
	//CMBND for V
	double afunc_st_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const;
	double afunc_st_V_pri(int cell1_idx, int cell2_idx, DBL3 shift) const;

	double bfunc_st_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const;
	double bfunc_st_V_pri(int cell1_idx, int cell2_idx) const;

	double diff2_st_V_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const;
	double diff2_st_V_pri(int cell1_idx, DBL3 shift) const;

	//CMBND for S
	DBL3 afunc_st_S_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const;
	DBL3 afunc_st_S_pri(int cell1_idx, int cell2_idx, DBL3 shift) const;

	double bfunc_st_S_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const;
	double bfunc_st_S_pri(int cell1_idx, int cell2_idx) const;

	DBL3 diff2_st_S_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const;
	DBL3 diff2_st_S_pri(int cell1_idx, DBL3 shift) const;

	//multiply spin accumulation by these to obtain spin potential, i.e. Vs = (De / elC) * (e/muB) * S, evaluated at the boundary
	double cfunc_sec(DBL3 relpos, DBL3 stencil) const;
	double cfunc_pri(int cell_idx) const;

	//-------------------Others

	//called by MoveMesh method in this mesh - move relevant transport quantities
	void MoveMesh_Transport(double x_shift);

	//-------------------Public calculation Methods

	//calculate the field resulting from a spin accumulation (spin current solver enabled) so a spin accumulation torque results when using the LLG or LLB equations
	void CalculateSAField(void);

	//Calculate the field resulting from interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set)
	void CalculateSAInterfaceField(Transport* ptrans_sec, CMBNDInfo& contact);

	//Calculate the interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set), accumulating result in displayVEC
	void CalculateDisplaySAInterfaceTorque(Transport* ptrans_sec, CMBNDInfo& contact);

	//calculate elC VEC using AMR and temperature information
	void CalculateElectricalConductivity(bool force_recalculate = false);

	//-------------------Display Calculation / Get Methods

	//return x, y, or z component of spin current (component = 0, 1, or 2)
	VEC<DBL3>& GetSpinCurrent(int component);
	//get value, assuming quantity already calculated with a call to the above method
	DBL3 GetCalculatedSpinCurrentValue(const DBL3& rel_pos) { if (displayVEC.linear_size()) return displayVEC[rel_pos]; else return DBL3(); }

	DBL3 GetAverageSpinCurrent(int component, Rect rectangle = Rect(), std::vector<MeshShape> shapes = {});

	//calculate charge current density over the mesh : applies to both charge-only transport and spin transport solvers (if check used inside this method)
	VEC_VC<DBL3>& GetChargeCurrent(void);
	//get value, assuming quantity already calculated with a call to the above method
	DBL3 GetCalculatedChargeCurrentValue(const DBL3& rel_pos) { if (displayVEC_VC.linear_size()) return displayVEC_VC[rel_pos]; else return DBL3(); }

	DBL3 GetAverageChargeCurrent(Rect rectangle = Rect(), std::vector<MeshShape> shapes = {});

	//return spin torque computed from spin accumulation
	VEC<DBL3>& GetSpinTorque(void);

	//get value, assuming quantity already calculated with a call to the above method
	DBL3 GetCalculatedSpinTorqueValue(const DBL3& rel_pos) { if (displayVEC.linear_size()) return displayVEC[rel_pos]; else return DBL3(); }

	//Get average bulk spin torque in given rectangle - calculate it first
	DBL3 GetAverageSpinTorque(Rect rectangle = Rect(), std::vector<MeshShape> shapes = {});

	//Get average interfacial spin torque in given rectangle - must have been calculated already in displayVEC through the supermesh transport module
	DBL3 GetAverageInterfacialSpinTorque(Rect rectangle = Rect(), std::vector<MeshShape> shapes = {});

#if COMPILECUDA == 1
	//Only use these if pModuleCUDA is enabled
	cu_obj<cuVEC<cuReal3>>& GetSpinCurrentCUDA(int component);
	cu_obj<cuVEC_VC<cuReal3>>& GetChargeCurrentCUDA(void);

	cu_obj<cuVEC<cuReal3>>& GetSpinTorqueCUDA(void);
#endif
};

#else

class Transport :
	public Modules
{
	friend STransport;

private:

	//used to compute spin current components and spin torque for display purposes - only calculated and memory allocated when asked to do so.
	VEC<DBL3> displayVEC;

	//used to compute charge current for display purposes (and in some cases we need to run vector calculus operations on it)
	VEC_VC<DBL3> displayVEC_VC;

#if COMPILECUDA == 1
	cu_obj<cuVEC<cuReal3>> displayVEC_CUDA;
	cu_obj<cuVEC_VC<cuReal3>> displayVEC_VC_CUDA;
#endif

private:

public:

	Transport(Mesh *pMesh_) {}
	~Transport() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Others

	//called by MoveMesh method in this mesh - move relevant transport quantities
	void MoveMesh_Transport(double x_shift) {}

	//-------------------Public calculation Methods

	//calculate the field resulting from a spin accumulation (spin current solver enabled) so a spin accumulation torque results when using the LLG or LLB equations
	void CalculateSAField(void) {}

	//Calculate the field resulting from interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set)
	void CalculateSAInterfaceField(Transport* ptrans_sec, CMBNDInfo& contact) {}

	//Calculate the interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set), accumulating result in displayVEC
	void CalculateDisplaySAInterfaceTorque(Transport* ptrans_sec, CMBNDInfo& contact) {}

	//calculate elC VEC using AMR and temperature information
	void CalculateElectricalConductivity(bool force_recalculate = false) {}

	//-------------------Display Calculation Methods

	//return x, y, or z component of spin current (component = 0, 1, or 2)
	VEC<DBL3>& GetSpinCurrent(int component) { return displayVEC; }
	//get value, assuming quantity already calculated with a call to the above method
	DBL3 GetCalculatedSpinCurrentValue(const DBL3& rel_pos) { return DBL3(); }

	DBL3 GetAverageSpinCurrent(int component, Rect rectangle = Rect(), std::vector<MeshShape> shapes = {}) { return DBL3(); }

	//calculate charge current density over the mesh : applies to both charge-only transport and spin transport solvers (if check used inside this method)
	VEC_VC<DBL3>& GetChargeCurrent(void) { return displayVEC_VC; }
	//get value, assuming quantity already calculated with a call to the above method
	DBL3 GetCalculatedChargeCurrentValue(const DBL3& rel_pos) { return DBL3(); }

	DBL3 GetAverageChargeCurrent(Rect rectangle = Rect(), std::vector<MeshShape> shapes = {}) { return DBL3(); }

	//return spin torque computed from spin accumulation
	VEC<DBL3>& GetSpinTorque(void) { return displayVEC; }

	//get value, assuming quantity already calculated with a call to the above method
	DBL3 GetCalculatedSpinTorqueValue(const DBL3& rel_pos) { return DBL3(); }

	//Get average bulk spin torque in given rectangle - calculate it first
	DBL3 GetAverageSpinTorque(Rect rectangle = Rect(), std::vector<MeshShape> shapes = {}) { return DBL3(); }

	//Get average interfacial spin torque in given rectangle - must have been calculated already in displayVEC through the supermesh transport module
	DBL3 GetAverageInterfacialSpinTorque(Rect rectangle = Rect(), std::vector<MeshShape> shapes = {}) { return DBL3(); }

#if COMPILECUDA == 1
	cu_obj<cuVEC<cuReal3>>& GetSpinCurrentCUDA(int component) { return displayVEC_CUDA; }
	cu_obj<cuVEC_VC<cuReal3>>& GetChargeCurrentCUDA(void) { return displayVEC_VC_CUDA;  }
	cu_obj<cuVEC<cuReal3>>& GetSpinTorqueCUDA(void) { return displayVEC_CUDA; }
#endif
};

#endif