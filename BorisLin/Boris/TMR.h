#pragma once

#include "BorisLib.h"
#include "Modules.h"
#include "TransportBase.h"

//This module is held in an insulator mesh (InsulatorMesh), which implements Mesh
class Mesh;
class Atom_Mesh;
class STransport;
class Transport;

#ifdef MODULE_COMPILATION_TMR

#if COMPILECUDA == 1
#include "TMRCUDA.h"
#endif

class TMR :
	public Modules,
	public TransportBase,
	public ProgramState<TMR, std::tuple<TEquation<double, double, double>, TEquation<double, double, double>>, std::tuple<>>
{
	friend STransport;

private:

	//pointer to mesh object holding this effective field module
	Mesh * pMesh;

	//the formula used to calculate conductivity dependence on angle (copy of flag from  InsulatorMesh)
	TMR_ TMR_type;

#if COMPILECUDA == 1
	friend TMRCUDA;

	//dynamic cast pModuleCUDA to this after new operator, so we don't have to cast it again later
	//do not check this for nullptr, but check pModuleCUDA instead; this is only to replace repeated use of dynamic_cast
	//thus if you check pModuleCUDA first, you can safely use pTMRCUDA since this is cast from pModuleCUDA only after pModuleCUDA allocated
	TMRCUDA* pTMRCUDA = nullptr;
#endif

	//Contacting top and bottom mehes used to calculate TMR
	//All magnetic meshes used (Ferromagnetic and Atomistic)
	std::vector<Mesh*> pMeshFM_Bot, pMeshFM_Top;
	std::vector<Atom_Mesh*> pMeshAtom_Bot, pMeshAtom_Top;

	//bias dependence equations for TMR RA values
	TEquation<double, double, double> RAtmr_p_equation, RAtmr_ap_equation;

private:

	//Update TEquation object with user constants values
	void UpdateTEquationUserConstants(bool makeCuda = true);

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

	//before iterating the spin solver (charge part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
	void PrimeSpinSolver_Charge(void);

	//take a single iteration of the charge transport solver (within the spin current solver) in this Mesh (cannot solve fully in one go as there may be other meshes so need boundary conditions set externally). Use SOR. 
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	DBL2 IterateSpinSolver_Charge_SOR(double damping);

	//call-back method for Poisson equation to evaluate RHS
	double Evaluate_SpinSolver_delsqV_RHS(int idx) const;

	//before iterating the spin solver (spin part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
	void PrimeSpinSolver_Spin(void);

	//solve for spin accumulation using Poisson equation for delsq_S. Use SOR. 
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	DBL2 IterateSpinSolver_Spin_SOR(double damping);

	//call-back method for Poisson equation for S
	DBL3 Evaluate_SpinSolver_delsqS_RHS(int idx) const;

	//Non-homogeneous Neumann boundary condition for V' - call-back method for Poisson equation for V
	DBL3 NHNeumann_Vdiff(int idx) const;

	//Non-homogeneous Neumann boundary condition for S' - call-back method for Poisson equation for S
	DBL33 NHNeumann_Sdiff(int idx) const;

	//------------------Others

	//check if dM_dt Calculation should be enabled
	bool Need_dM_dt_Calculation(void);

	//check if the delsq_V_fixed VEC is needed
	bool Need_delsq_V_fixed_Precalculation(void);

	//check if the delsq_S_fixed VEC is needed
	bool Need_delsq_S_fixed_Precalculation(void);

public:

	TMR(Mesh *pMesh_);
	~TMR();

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Output Data

	//calculate resistance in given rectangle (relative to mesh)
	double GetResistance(Rect rect);

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

	//-------------------Public calculation Methods

	//calculate the field resulting from a spin accumulation (spin current solver enabled) so a spin accumulation torque results when using the LLG or LLB equations
	void CalculateSAField(void);

	//Calculate the field resulting from interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set)
	void CalculateSAInterfaceField(TransportBase* ptrans_sec, CMBNDInfo& contact);

	//Calculate the interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set), accumulating result in displayVEC
	void CalculateDisplaySAInterfaceTorque(TransportBase* ptrans_sec, CMBNDInfo& contact);

	//calculate elC VEC using AMR and temperature information
	void CalculateElectricalConductivity(bool force_recalculate = false);

	//-------------------Display Calculation / Get Methods

	//return x, y, or z component of spin current (component = 0, 1, or 2)
	VEC<DBL3>& GetSpinCurrent(int component);

	//calculate charge current density over the mesh : applies to both charge-only transport and spin transport solvers (if check used inside this method)
	VEC_VC<DBL3>& GetChargeCurrent(void);

	//return spin torque computed from spin accumulation
	VEC<DBL3>& GetSpinTorque(void);

	//-------------------Other Settings

	BError SetBiasEquation_Parallel(std::string equation_string);
	BError SetBiasEquation_AntiParallel(std::string equation_string);
};

#else

class TMR :
	public Modules,
	public TransportBase
{

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

	//Update TEquation object with user constants values
	void UpdateTEquationUserConstants(bool makeCuda = true) {}

	//-------------------Calculation Methods

	//calculate electric field as the negative gradient of V
	void CalculateElectricField(void) {}

	//Charge transport only : V

	//take a single iteration of the charge transport solver in this Mesh (cannot solve fully in one go as there may be other meshes so need boundary conditions set externally). Use SOR. 
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	DBL2 IterateChargeSolver_SOR(double damping) { return DBL2(); }

	//call-back method for Poisson equation to evaluate RHS
	double Evaluate_ChargeSolver_delsqV_RHS(int idx) const { return 0.0; }

	//Calculation Methods used by Spin Current Solver only

	//before iterating the spin solver (charge part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
	void PrimeSpinSolver_Charge(void) {}

	//take a single iteration of the charge transport solver (within the spin current solver) in this Mesh (cannot solve fully in one go as there may be other meshes so need boundary conditions set externally). Use SOR. 
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	DBL2 IterateSpinSolver_Charge_SOR(double damping) { return DBL2(); }

	//call-back method for Poisson equation to evaluate RHS
	double Evaluate_SpinSolver_delsqV_RHS(int idx) const { return 0.0; }

	//before iterating the spin solver (spin part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
	void PrimeSpinSolver_Spin(void) {}

	//solve for spin accumulation using Poisson equation for delsq_S. Use SOR. 
	//Return un-normalized error (maximum change in quantity from one iteration to the next) - first - and maximum value  -second - divide them to obtain normalized error
	DBL2 IterateSpinSolver_Spin_SOR(double damping) { return DBL2(); }

	//call-back method for Poisson equation for S
	DBL3 Evaluate_SpinSolver_delsqS_RHS(int idx) const { return DBL3(); }

	//Non-homogeneous Neumann boundary condition for V' - call-back method for Poisson equation for V
	DBL3 NHNeumann_Vdiff(int idx) const { return DBL3(); }

	//Non-homogeneous Neumann boundary condition for S' - call-back method for Poisson equation for S
	DBL33 NHNeumann_Sdiff(int idx) const { return DBL33(); }

	//------------------Others

	//check if dM_dt Calculation should be enabled
	bool Need_dM_dt_Calculation(void) { return false; }

	//check if the delsq_V_fixed VEC is needed
	bool Need_delsq_V_fixed_Precalculation(void) { return false; }

	//check if the delsq_S_fixed VEC is needed
	bool Need_delsq_S_fixed_Precalculation(void) { return false; }

public:

	TMR(Mesh *pMesh_) {}
	~TMR() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Output Data

	//calculate resistance in given rectangle (relative to mesh)
	double GetResistance(Rect rect) { return 0.0; }

	//-------------------CMBND computation methods

	//CMBND values set based on continuity of a potential and flux

	//Charge transport only : V

	//For V only : V and Jc are continuous; Jc = -sigma * grad V = a + b * grad V -> a = 0 and b = -sigma taken at the interface
	double afunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const { return 0.0; }
	double afunc_V_pri(int cell1_idx, int cell2_idx, DBL3 shift) const { return 0.0; }

	double bfunc_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const { return 0.0; }
	double bfunc_V_pri(int cell1_idx, int cell2_idx) const { return 0.0; }

	//second order differential of V at cells either side of the boundary; delsq V = -grad V * grad elC / elC
	double diff2_V_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const { return 0.0; }
	double diff2_V_pri(int cell1_idx, DBL3 shift) const { return 0.0; }

	//Spin transport : V and S
	//CMBND for V
	double afunc_st_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const { return 0.0; }
	double afunc_st_V_pri(int cell1_idx, int cell2_idx, DBL3 shift) const { return 0.0; }

	double bfunc_st_V_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const { return 0.0; }
	double bfunc_st_V_pri(int cell1_idx, int cell2_idx) const { return 0.0; }

	double diff2_st_V_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const { return 0.0; }
	double diff2_st_V_pri(int cell1_idx, DBL3 shift) const { return 0.0; }

	//CMBND for S
	DBL3 afunc_st_S_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const { return DBL3(); }
	DBL3 afunc_st_S_pri(int cell1_idx, int cell2_idx, DBL3 shift) const { return DBL3(); }

	double bfunc_st_S_sec(DBL3 relpos_m1, DBL3 shift, DBL3 stencil) const { return 0.0; }
	double bfunc_st_S_pri(int cell1_idx, int cell2_idx) const { return 0.0; }

	DBL3 diff2_st_S_sec(DBL3 relpos_m1, DBL3 stencil, DBL3 shift) const { return DBL3(); }
	DBL3 diff2_st_S_pri(int cell1_idx, DBL3 shift) const { return DBL3(); }

	//multiply spin accumulation by these to obtain spin potential, i.e. Vs = (De / elC) * (e/muB) * S, evaluated at the boundary
	double cfunc_sec(DBL3 relpos, DBL3 stencil) const { return 0.0; }
	double cfunc_pri(int cell_idx) const { return 0.0; }

	//-------------------Public calculation Methods

	//calculate the field resulting from a spin accumulation (spin current solver enabled) so a spin accumulation torque results when using the LLG or LLB equations
	void CalculateSAField(void) {}

	//Calculate the field resulting from interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set)
	void CalculateSAInterfaceField(TransportBase* ptrans_sec, CMBNDInfo& contact) {}

	//Calculate the interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set), accumulating result in displayVEC
	void CalculateDisplaySAInterfaceTorque(TransportBase* ptrans_sec, CMBNDInfo& contact) {}

	//calculate elC VEC using AMR and temperature information
	void CalculateElectricalConductivity(bool force_recalculate = false) {}

	//-------------------Display Calculation Methods

	//return x, y, or z component of spin current (component = 0, 1, or 2)
	VEC<DBL3>& GetSpinCurrent(int component) { return displayVEC; }

	//calculate charge current density over the mesh : applies to both charge-only transport and spin transport solvers (if check used inside this method)
	VEC_VC<DBL3>& GetChargeCurrent(void) { return displayVEC_VC; }

	//return spin torque computed from spin accumulation
	VEC<DBL3>& GetSpinTorque(void) { return displayVEC; }

	//-------------------Other Settings

	BError SetBiasEquation_Parallel(std::string equation_string) { return BError(); }
	BError SetBiasEquation_AntiParallel(std::string equation_string) { return BError(); }
};

#endif