#pragma once

#include "BorisLib.h"

#include "ErrorHandler.h"

#include "Boris_Enums_Defs.h"

#include "DiffEq_Defs.h"

#if COMPILECUDA == 1
#include "DiffEqCUDA.h"
#endif

using namespace std;

//------------------------------------------------------------------------------------------------------

class Mesh;
class FMesh;
class SuperMesh;

class DifferentialEquation;
typedef DBL3 (DifferentialEquation::*Equation)(int idx);

//"static" base class for DifferentialEquation objects.
//Each DifferentialEquation object is held in its ferromagnetic mesh and has different internal storage dimensions depending on the mesh dimensions
//There are a number of properties that are in common between all the DifferentialEquation objects, since the same equation and evaluation method must be used for all the ferromagnetic meshes.
//Moreover all the DifferentialEquation objects must be coordinated by a top-level object : ODECommon also provides this by keeping track of all of them 

class ODECommon :
	public ProgramState<ODECommon,
	tuple<int, int, double, double, double, int, int, double, bool, bool, double, double>,
	tuple<>>
{
#if COMPILECUDA == 1
	friend ODECommonCUDA;
#endif
private:

	//parameters for iterations and stages control
	static int iteration, stageiteration;

	//time and stagetime
	static double time, stagetime;

	//evaluation step for higher order evaluation methods
	static int evalStep;

	//flag to indicate mesh magnetization is available for use in output saves - e.g. at the last step of an evaluation method
	static bool available;

	//some evaluation methods need to be primed (e.g. ABM)
	static bool primed;

	//final mxh value after each iteration
	static double mxh;

protected:
	
	//only renormalize M after solving for a time step for some equations. For LLB type equations must not renormalize.
	static bool renormalize;

	//does the ODE require full spin current solver? (set by SA ODE versions, which require spin accumulation to evaluate spin torques)
	static bool solve_spin_current;

	//currently set evaluation method
	static int setODE;

	//set evaluation method for given ODE
	static int evalMethod;

	//time-step
	static double dT;
	//dT_last = dT, unless an adaptive time step method is used, where dT is saved before being changed - needed to calculate dM/dt for adaptive time step methods
	static double dT_last;

	//used to alternate between past equation evaluations (e.g. for ABM)
	static bool alternator;

	//collection of all set ODE in ferromagnetic meshes : this is handled by the DifferentialEquation constructor and destructor
	//When a new ODE is made in a ferromagnetic mesh, DifferentialEquation constructor is called and it pushes a pointer in this vector, together with a unique id for it.
	//When a ODE is deleted in a ferromagnetic mesh, DifferentialEquation destructor is called, and it erases the entry n this vector using the unique odeId previously generated.
	static vector_lut<DifferentialEquation*> pODE;

	//function pointer to equation to solve
	static Equation equation;

	//use moving mesh algorithm?
	static bool moving_mesh;

	//moving mesh algorithm type: antisymmetric for DWs, symmetric for skyrmions
	static bool moving_mesh_antisymmetric;

	//normalised threshold value used to trigger a mesh shift
	static double moving_mesh_threshold;

	//current dw shift as resulting from moving mesh algorithm
	static double moving_mesh_dwshift;

#if COMPILECUDA == 1
	static ODECommonCUDA *pODECUDA;
#endif

private:

	//out of all currently set meshes find maximum mxh value and set it in mxh. This is called after all meshes have iterated, so their respective mxh values are in mxh_reduction[0]
	void Set_mxh(void);

	//get largest local truncation error from all the meshes
	double Get_lte(void);

public:

	ODECommon(bool called_from_derived = false);
	virtual ~ODECommon();

	//implement pure virtual method from ProgramState
	void RepairObjectState(void);

	//---------------------------------------- SET-UP METHODS : DiffEq_Common.cpp

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	//switch CUDA state on/off
	BError SwitchCUDAState(bool cudaState);

	BError SetODE(ODE_ setODE_, EVAL_ evalMethod_, bool set_eval_method = true);

	BError SetEvaluationMethod(EVAL_ evalMethod_);

	void Reset(void);

	void NewStage(void);

	void SetdT(double dT);

	void SetMoveMeshTrigger(bool status, int meshId = -1);
	void SetMoveMeshAntisymmetric(bool antisymmetric) { moving_mesh_antisymmetric = antisymmetric; }
	void SetMoveMeshThreshold(double threshold) { moving_mesh_threshold = threshold; }

	//---------------------------------------- ITERATE METHODS : DiffEq_Iterate.cpp (and DiffEq_IterateCUDA.cpp)

	//advance time using the set ODE solver
	void Iterate(void);

#if COMPILECUDA == 1
	void IterateCUDA(void);
#endif

	//moving mesh algorithm : shift all meshes using their MoveMesh method
	void MovingMeshAlgorithm(SuperMesh* pSMesh);

	//---------------------------------------- GET METHODS : DiffEqCommon.cpp

	int GetIteration(void) { return iteration; }
	int GetStageIteration(void) { return stageiteration; }

	double GetTime(void) { return time; }
	double GetStageTime(void) { return stagetime; }
	double GetTimeStep(void) { return dT; }
	double Get_mxh(void);

	bool IsMovingMeshSet(void) { return moving_mesh; }
	int GetId_of_MoveMeshTrigger(void);

	bool MoveMeshAntisymmetric(void) { return moving_mesh_antisymmetric; }
	double MoveMeshThreshold(void) { return moving_mesh_threshold; }

	double Get_dwshift(void) { return moving_mesh_dwshift; }

	bool TimeStepSolved(void) { return available; }

	bool SolveSpinCurrent(void) { return solve_spin_current; }

	void QueryODE(ODE_ &setODE_, EVAL_ &evalMethod_) { setODE_ = (ODE_)setODE; evalMethod_ = (EVAL_)evalMethod; }
};

//The actual differential equation specific to a certain mesh (the ferromagnetic pMesh), but with settings in ODECommon common to all meshes
class DifferentialEquation : 
	public ODECommon 
{
	friend ODECommon;

#if COMPILECUDA == 1
	friend DifferentialEquationCUDA;
#endif

private: //private data

	//if the object couldn't be created properly in the constructor an error is set here
	BError error_on_create;

	//reductions
	OmpReduction<double> mxh_reduction;
	OmpReduction<DBL3> mxh_av_reduction;
	OmpReduction<double> lte_reduction;

	//Used to save starting magnetization - all evaluation methods do this, even when not needed by the method itself, so we can calculate dM/dt when needed.
	VEC<DBL3> sM1;

	//Used for RK4 (0, 1, 2); ABM (0, 1)
	VEC<DBL3> sEval0, sEval1, sEval2;

	//Additional for use with RKF45
	VEC<DBL3> sEval3, sEval4;

	//Thermal field and torques, enabled only for the stochastic equations
	VEC<DBL3> H_Thermal, Torque_Thermal;

	//random number generator
	BorisRand prng;

	FMesh *pMesh;

	//unique odeId generated when a new entry is made in the pODE vector : used to delete it in the destructor.
	INT2 odeId;

#if COMPILECUDA == 1
	DifferentialEquationCUDA *pmeshODECUDA = nullptr;
#endif

private: //private methods

	//---------------------------------------- SOLVER METHODS : DiffEq_Evals.cpp

	//Euler evaluation of ODE
	void RunEuler(void);

	//Trapezoidal Euler evaluation of ODE
	void RunTEuler_Step0(void);
	void RunTEuler_Step1(void);

	//RK4 evaluation of ODE
	void RunRK4_Step0(void);
	void RunRK4_Step1(void);
	void RunRK4_Step2(void);
	void RunRK4_Step3(void);

	//ABM
	void RunABM_Predictor(void);
	void RunABM_Corrector(void);
	void RunABM_TEuler0(void);
	void RunABM_TEuler1(void);

	//RKF45
	void RunRKF45_Step0(void);
	void RunRKF45_Step1(void);
	void RunRKF45_Step2(void);
	void RunRKF45_Step3(void);
	void RunRKF45_Step4(void);
	void RunRKF45_Step5(void);

	//---------------------------------------- OTHERS : DiffEq_Evals.cpp

	//Restore magnetisation after a failed step for adaptive time-step methods
	void RestoreMagnetisation(void);

	//---------------------------------------- OTHER CALCULATION METHODS : DiffEq_SEquations.cpp

	//called when using stochastic equations
	void GenerateThermalField(void);
	void GenerateThermalField_and_Torque(void);

	//---------------------------------------- SET-UP METHODS  : DiffEq.cpp

	//allocate memory depending on set evaluation method - also cleans up previously allocated memory by calling CleanupMemory();
	BError AllocateMemory(void);

	//deallocate memory before re-allocating it (depending on evaluation method previously allocated memory might not be used again, so need clean-up before)
	void CleanupMemory(void);

	//---------------------------------------- EQUATIONS : DiffEq_Equations.cpp and DiffEq_SEquations.cpp

	//Landau-Lifshitz-Gilbert equation
	DBL3 LLG(int idx);

	//Landau-Lifshitz-Gilbert equation with Zhang-Li STT
	DBL3 LLGSTT(int idx);

	//Landau-Lifshitz-Bloch equation
	DBL3 LLB(int idx);

	//Landau-Lifshitz-Bloch equation with Zhang-Li STT
	DBL3 LLBSTT(int idx);

	//Stochastic Landau-Lifshitz-Gilbert equation
	DBL3 SLLG(int idx);

	//Stochastic Landau-Lifshitz-Gilbert equation with Zhang-Li STT
	DBL3 SLLGSTT(int idx);

	//Stochastic Landau-Lifshitz-Bloch equation
	DBL3 SLLB(int idx);

	//Stochastic Landau-Lifshitz-Bloch equation with Zhang-Li STT
	DBL3 SLLBSTT(int idx);

public:  //public methods

	DifferentialEquation(FMesh *pMesh);
	~DifferentialEquation();

	BError Error_On_Create(void) { return error_on_create; }

	//---------------------------------------- SET-UP METHODS : DiffEq.cpp

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage = UPDATECONFIG_GENERIC);

	//switch CUDA state on/off
	BError SwitchCUDAState(bool cudaState);

	//---------------------------------------- GETTERS

	//return dM by dT - should only be used when evaluation sequence has ended (TimeStepSolved() == true)
	DBL3 dMdt(int idx);

#if COMPILECUDA == 1
	//get stored cuda differential equation pointer (pmeshODECUDA)
	DifferentialEquationCUDA* Get_DifferentialEquationCUDA_ptr(void) { return pmeshODECUDA; }
#endif
};
