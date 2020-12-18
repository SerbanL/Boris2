#pragma once

#include "BorisLib.h"
#include "Modules.h"
#include "Transport.h"

#if COMPILECUDA == 1
#include "STransportCUDA.h"
#endif



class SuperMesh;

#ifdef MODULE_COMPILATION_TRANSPORT

class STransport :
	public Modules,
	public ProgramState<STransport, std::tuple<vector_lut<Rect>, std::vector<double>, int, double, double, double, double, bool, double, int, double, int, DBL2, TEquation<double>, TEquation<double>>, std::tuple<>>
{

#if COMPILECUDA == 1
	friend STransportCUDA;

	//dynamic cast pModuleCUDA to this after new operator, so we don't have to cast it again later
	//do not check this for nullptr, but check pModuleCUDA instead; this is only to replace repeated use of dynamic_cast
	//thus if you check pModuleCUDA first, you can safely use pTransportCUDA since this is cast from pModuleCUDA only after pModuleCUDA allocated
	STransportCUDA* pSTransportCUDA = nullptr;
#endif

private:

	//pointer to supermesh
	SuperMesh* pSMesh;

	//---------------------- CMBND data

	//CMBND contacts for all contacting transport meshes - these are ordered by first vector index; for each transport mesh there could be multiple contacting meshes and these are ordered by second vector index
	//CMBNDInfo describes the contact between 2 meshes, allowing calculation of values at cmbnd cells based on continuity of a potential and flux
	std::vector< std::vector<CMBNDInfo> > CMBNDcontacts;

	//list of all transport modules in transport meshes (same ordering as first vector in CMBNDcontacts)
	std::vector<Transport*> pTransport;

	//vector of pointers to all V - need this to set cmbnd flags (same ordering as first vector in CMBNDcontacts)
	std::vector<VEC_VC<double>*> pV;

	//vector of pointers to all S - need this to set cmbnd flags (same ordering as first vector in CMBNDcontacts)
	std::vector<VEC_VC<DBL3>*> pS;

	//----------------------

	//electrodes set Dirichlet boundary conditions : use Rects (metric values) to mark these (only cells intersecting with these Rects on mesh surfaces will be marked).
	vector_lut<Rect> electrode_rects;
	
	//fixed potential values on the electrodes (the Dirichlet boundary values)
	std::vector<double> electrode_potentials;

	//index of electrode designated as ground (-1 if none). Note: need a ground to define current polarity and calculate total current.
	int ground_electrode_index = -1;

	//potential (V) drop to ground
	double current = 0.0, net_current = 0.0, potential = 0.0;

	//Set potential (or current) using user equation allowing dependence on stage time (t); stage step (Ss) is introduced as user constant.
	TEquation<double> V_equation, I_equation;

	//sample resistance defined as potential / current
	double resistance = 0.0;

	bool constant_current_source = false;

	//threshold for considering laplace (or Poisson) equation solved : should actually be lower than this in host code (cuda device code struggles to relax below this value so setting this higher value instead).
	double errorMaxLaplace = 1e-6;
	
	//maximum number of iterations to take when solving the Laplace (or Poisson) equation
	int maxLaplaceIterations = 500;

	//iterations taken to converge to given errorMaxLaplace
	int iters_to_conv = 0;

	//maximum error (for convergence) and number of iterations used for spin accumulation solver
	double s_errorMax = 1e-5;
	int s_maxIterations = 200;

	//iterations taken by spin accumulation solver to converge to given s_errorMax
	int s_iters_to_conv = 0;

	//fixed SOR damping to use for V (first value) and S (second value) Poisson equations
	DBL2 SOR_damping = DBL2(1.4, 0.5);

	//after transport solver has relaxed below errorMaxLaplace, it only needs to be updated if relevant quantities change (e.g. potential, conductivity)
	//When these changes occur this flag is set to true.
	bool recalculate_transport = true;

	//after Transport solver has run following a recalculate_transport flag check, set transport_recalculated to true so other dependent modules can also update - this flag will be reset here on next UpdateField.
	bool transport_recalculated = true;

private:

	//adjust potential values on electrodes to give the specified potential drop (set an asymmetric potential drop), but do not change any constant current source settings
	void adjust_potential(double potential_);

	//set potential values using a slope between the potential values of ground and another electrode (if set)
	void initialize_potential_values(void);

	//-----Charge Transport only

	//solve for V and Jc in all meshes using SOR
	void solve_charge_transport_sor(void);

	//calculate and set values at composite media boundaries for V (charge transport only) after all other cells have been computed and set
	void set_cmbnd_charge_transport(void);

	//-----Spin and Charge Transport

	//solve for V, Jc and S in all meshes using SOR for Poisson equation and FTCS for S equation
	void solve_spin_transport_sor(void);

	//calculate and set values at composite media boundaries for V (when using spin transport solver)
	void set_cmbnd_spin_transport_V(void);

	//calculate and set values at composite media boundaries for S
	void set_cmbnd_spin_transport_S(void);

	//functions to specify boundary conditions for interface conductance approach for charge : Jc_N = Jc_F = A + B * dV
	double Afunc_V(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, Transport& trans_sec, Transport& trans_pri) const;
	double Bfunc_V(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, Transport& trans_sec, Transport& trans_pri) const;

	//functions to specify boundary conditions for interface conductance approach for spin : Js_N = A_N + B_N * dVs, Js_F = A_F + B_F * dVs
	DBL3 Afunc_N_S(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, Transport& trans_sec, Transport& trans_pri) const;
	DBL3 Afunc_F_S(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, Transport& trans_sec, Transport& trans_pri) const;

	DBL33 Bfunc_N_S(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, Transport& trans_sec, Transport& trans_pri) const;
	DBL33 Bfunc_F_S(int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 shift, DBL3 stencil, Transport& trans_sec, Transport& trans_pri) const;

	//Calculate interface spin accumulation torque (in magnetic meshes for NF interfaces with G interface conductance set)
	void CalculateSAInterfaceField(void);

	//Update TEquation object with user constants values
	void UpdateTEquationUserConstants(void);

public:

	STransport(SuperMesh *pSMesh_);
	~STransport() {}

	//-------------------Implement ProgramState method

	void RepairObjectState(void) {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; recalculate_transport = true; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage);

	BError MakeCUDAModule(void);

	double UpdateField(void);

	//-------------------Display Calculation Methods

	//return interfacial spin torque in given mesh with matching transport module
	VEC<DBL3>& GetInterfacialSpinTorque(Transport* pMeshTrans);

#if COMPILECUDA == 1
	//return interfacial spin torque in given mesh with matching transport module
	cu_obj<cuVEC<cuReal3>>& GetInterfacialSpinTorqueCUDA(Transport* pMeshTrans);
#endif

	//-------------------Getters

	//get currently set potential
	double GetPotential(void) { return potential; }
	
	//calculate new current value for the set potential. If a constant current source is used, the potential is adjusted instead in order to keep the set current value.
	DBL2 GetCurrent(void);

	//get sample resistance by calculating current first
	double GetResistance(void);

	bool UsingConstantCurrentSource(void) { return constant_current_source; }

	int GetItersToConv(void) { return iters_to_conv; }
	int GetSItersToConv(void) { return s_iters_to_conv; }

	bool Transport_Recalculated(void) { return transport_recalculated; }

	double GetConvergenceError(void) { return errorMaxLaplace; }
	int GetConvergenceTimeout(void) { return maxLaplaceIterations; }

	double GetSConvergenceError(void) { return s_errorMax; }
	int GetSConvergenceTimeout(void) { return s_maxIterations; }

	//get fixed SOR damping values (for V and S solvers)
	DBL2 GetSORDamping(void) { return SOR_damping; }

	//-------------------Setters

	void Flag_Recalculate_Transport(void) { recalculate_transport = true; }

	//set potential value, also reset any constant current source settings; by default clear the text equation setting, unless we set the value from the equation evaluation (set flag to false then)
	void SetPotential(double potential_, bool clear_equation = true);

	//set current value, also reset any constant potential source settings; by default clear the text equation setting, unless we set the value from the equation evaluation (set flag to false then)
	void SetCurrent(double current_, bool clear_equation = true);

	void SetConvergenceError(double errorMaxLaplace_, int maxLaplaceIterations_);
	void SetSConvergenceError(double s_errorMax_, int s_maxIterations_);

	//set fixed SOR damping values (for V and S solvers)
	void SetSORDamping(DBL2 _SOR_damping);

	//set text equation from std::string
	BError SetPotentialEquation(std::string equation_string, int step);
	BError SetCurrentEquation(std::string equation_string, int step);

	//-------------------Electrodes methods

	//add a new electrode
	void AddElectrode(double electrode_potential, Rect electrode_rect);
	
	//delete electrode with given index
	bool DelElectrode(int index);
	
	//delete all electrodes
	void ClearElectrodes(void);
	
	void DesignateGroundElectrode(int index);

	//set box for existing electrode with given index
	bool SetElectrodeRect(int index, Rect electrode_rect);

	//set potential on existing electrode with given index.
	//NOTE : all potential value changes (including adjustments due to constant current source) eventually call this method. Set flag to recalculate transport.
	bool SetElectrodePotential(int index, double electrode_potential);
	
	bool IsGroundElectrode(int index) { return (index == ground_electrode_index); }
	
	//get index of ground electrode, if any. If none then return -1
	int GetGroundElectrodeIndex(void) { return ground_electrode_index; }

	//get electrode index for electrode with given minor Id (major Id is 0)
	int GetElectrodeIndex(int minorId) { return electrode_rects.get_index_from_id(INT2(0, minorId)); }

	//return number of electrodes
	int GetNumberofElectrodes(void) { return (int)electrode_rects.size(); }

	//get electrode info (Rect and potential value) of electrode with given index
	std::pair<Rect, double> GetElectrodeInfo(int index);

	//get the minor id of electrode with given index
	int GetElectrodeid(int index);
};

#else

class STransport :
	public Modules
{

#if COMPILECUDA == 1
	friend STransportCUDA;
#endif

private:

private:

public:

	STransport(SuperMesh *pSMesh_) {}
	~STransport() {}

	//-------------------Abstract base class method implementations

	void Uninitialize(void) {}

	BError Initialize(void) { return BError(); }

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage) { return BError(); }
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	BError MakeCUDAModule(void) { return BError(); }

	double UpdateField(void) { return 0.0; }

	//-------------------Display Calculation Methods

	//return interfacial spin torque in given mesh with matching transport module
	VEC<DBL3>& GetInterfacialSpinTorque(Transport* pMeshTrans) { return pMeshTrans->displayVEC; }

#if COMPILECUDA == 1
	//return interfacial spin torque in given mesh with matching transport module
	cu_obj<cuVEC<cuReal3>>& GetInterfacialSpinTorqueCUDA(Transport* pMeshTrans) { return pMeshTrans->displayVEC_CUDA; }
#endif

	//-------------------Getters

	//get currently set potential
	double GetPotential(void) { return 0.0; }

	//calculate new current value for the set potential. If a constant current source is used, the potential is adjusted instead in order to keep the set current value.
	DBL2 GetCurrent(void) { return DBL2(); }

	//get sample resistance by calculating current first
	double GetResistance(void) { return 0.0; }

	bool UsingConstantCurrentSource(void) { return false; }

	int GetItersToConv(void) { return 0; }
	int GetSItersToConv(void) { return 0; }

	bool Transport_Recalculated(void) { return true; }

	double GetConvergenceError(void) { return 0.0; }
	int GetConvergenceTimeout(void) { return 0; }

	double GetSConvergenceError(void) { return 0.0; }
	int GetSConvergenceTimeout(void) { return 0; }

	//get fixed SOR damping values (for V and S solvers)
	DBL2 GetSORDamping(void) { return DBL2(); }

	//-------------------Setters

	void Flag_Recalculate_Transport(void) {}

	//set potential value, also reset any constant current source settings
	void SetPotential(double potential_, bool clear_equation) {}

	void SetCurrent(double current_, bool clear_equation) {}

	void SetConvergenceError(double errorMaxLaplace_, int maxLaplaceIterations_) {}
	void SetSConvergenceError(double s_errorMax_, int s_maxIterations_) {}

	//set fixed SOR damping values (for V and S solvers)
	void SetSORDamping(DBL2 _SOR_damping) {}

	//set text equation from std::string
	BError SetPotentialEquation(std::string equation_string, int step) { return BError(); }
	BError SetCurrentEquation(std::string equation_string, int step) { return BError(); }

	//-------------------Electrodes methods

	//add a new electrode
	void AddElectrode(double electrode_potential, Rect electrode_rect) {}

	//delete electrode with given index
	bool DelElectrode(int index) { return true; }

	//delete all electrodes
	void ClearElectrodes(void) {}

	void DesignateGroundElectrode(int index) {}

	//set box for existing electrode with given index
	bool SetElectrodeRect(int index, Rect electrode_rect) { return true; }

	//set potential on existing electrode with given index.
	//NOTE : all potential value changes (including adjustments due to constant current source) eventually call this method. Set flag to recalculate transport.
	bool SetElectrodePotential(int index, double electrode_potential) { return true; }

	bool IsGroundElectrode(int index) { return false; }

	//get index of ground electrode, if any. If none then return -1
	int GetGroundElectrodeIndex(void) { return -1; }

	//get electrode index for electrode with given minor Id (major Id is 0)
	int GetElectrodeIndex(int minorId) { return -1; }

	//return number of electrodes
	int GetNumberofElectrodes(void) { return 0; }

	//get electrode info (Rect and potential value) of electrode with given index
	std::pair<Rect, double> GetElectrodeInfo(int index) { return std::pair<Rect, double>(); }

	//get the minor id of electrode with given index
	int GetElectrodeid(int index) { return -1; }
};

#endif
