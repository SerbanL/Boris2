#pragma once

#include "BorisLib.h"

#include "ErrorHandler.h"
#include "CompileFlags.h"

#include "Boris_Enums_Defs.h"

#include "DiffEq_Defs.h"

#if COMPILECUDA == 1
#include "DiffEq_CommonCUDA.h"
#endif

using namespace std;

//------------------------------------------------------------------------------------------------------

class Mesh;
class SuperMesh;

class DifferentialEquation;
typedef DBL3(DifferentialEquation::*Equation)(int idx);

//"static" base class for DifferentialEquation objects.
//Each DifferentialEquation type object is held in its (anti)ferromagnetic mesh and has different internal storage dimensions depending on the mesh dimensions
//There are a number of properties that are in common between all the DifferentialEquation objects, since the same equation and evaluation method must be used for all the (anti)ferromagnetic meshes.
//Moreover all the DifferentialEquation objects must be coordinated by a top-level object : ODECommon also provides this by keeping track of all of them 

class ODECommon :
	public ProgramState<ODECommon,
	tuple<int, int, double, double, double, double, int, int, double, double, double, double, double, double, double, int, bool, bool, double, double>,
	tuple<>>
{
#if COMPILECUDA == 1
	friend ODECommonCUDA;
#endif
private:

	//-----------------------------------Primary data

	//parameters for iterations and stages control
	static int iteration, stageiteration;

	//time and stagetime
	static double time, stagetime;

	//to avoid calculating mxh every single iteration whether it's needed or not, only calculate it if this flag is set
	//after updating mxh value set this flag to false again. Thus when a call to Get_mxh is received, set this flag to true.
	//this means the actual mxh value returned by Get_mxh will be off by 1 iteration but it doesn't really matter
	//if a schedule with mxh stop condition is used, all it means is the schedule will run for an extra iteration
	//saved mxh values will also be off by 1 iteration but again it doesn't matter.
	static bool calculate_mxh;

	//final mxh value after each iteration
	static double mxh;

	//as for mxh but for dmdt stopping condition
protected:
	static bool calculate_dmdt;
private:

	//final dmdt value after each iteration
	static double dmdt;

	//-----------------------------------Multistep evaluation methods control

	//evaluation step for higher order evaluation methods
	static int evalStep;

	//flag to indicate mesh magnetization is available for use in output saves - e.g. at the last step of an evaluation method
	static bool available;

	//some evaluation methods need to be primed (e.g. ABM); 
	static bool primed;

	//-----------------------------------Adaptive time step control

	//parameters for adaptive time step calculation method
	//fail above this relative error : repeat with lower time step
	static double err_high_fail;

	//decrease time step above this relative error - decrease factor depends on error ratio.
	static double err_high;

	//increase time step above this relative error
	static double err_low;

	//when increasing dT multiply by this (> 1.0)
	static double dT_increase;

	//-----------------------------------Properties flags

	//does the ODE require full spin current solver? (set by SA ODE versions, which require spin accumulation to evaluate spin torques)
	static bool solve_spin_current;

	//use evaluation speedup by skipping certain effective field contributions at specific steps in the evaluation?
	//this takes on a value from EVALSPEEDUP_ enum
	static int use_evaluation_speedup;

protected:

	//-----------------------------------Collection of ferromagnetic meshes

	//collection of all set ODE in ferromagnetic meshes : this is handled by the DifferentialEquation constructor and destructor
	//When a new ODE is made in a ferromagnetic mesh, DifferentialEquation constructor is called and it pushes a pointer in this vector, together with a unique id for it.
	//When a ODE is deleted in a ferromagnetic mesh, DifferentialEquation destructor is called, and it erases the entry n this vector using the unique odeId previously generated.
	static vector_lut<DifferentialEquation*> pODE;

	//-----------------------------------Equation and Evaluation method values

	//currently set evaluation method
	static int setODE;

	//set evaluation method for given ODE
	static int evalMethod;

	//function pointer to equation to solve
	static Equation equation;

	//-----------------------------------Time step

	//time-step
	static double dT;
	//dT_last = dT, unless an adaptive time step method is used, where dT is saved before being changed - needed to calculate dM/dt for adaptive time step methods
	static double dT_last;

	//-----------------------------------Adaptive time step control

	//maximum and minimum dT values allowed
	static double dT_min;
	static double dT_max;

	//-----------------------------------Steepest Descent Solver

	//quantities used to calculate Barzilai-Borwein stepsizes across multiple meshes
	//Accumulate values in these quantities, then obtain stepsizes as:
	//step1 = delta_M_sq / delta_M_dot_delta_G
	//step2 = delta_M_dot_delta_G / delta_G_sq
	static double delta_M_sq, delta_M2_sq;
	static double delta_G_sq, delta_G2_sq;
	static double delta_M_dot_delta_G, delta_M2_dot_delta_G2;

	//-----------------------------------Evaluation method modifiers

	//only renormalize M after solving for a time step for some equations. For LLB type equations must not renormalize.
	static bool renormalize;

	//-----------------------------------Special values

	//used to alternate between past equation evaluations (e.g. for ABM)
	static bool alternator;

	//-----------------------------------Moving mesh data

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

	//out of all currently set meshes find maximum dmdt value and set it in dmdt. This is called after all meshes have iterated, so their respective mxh values are in dmdt_reduction[0]
	void Set_dmdt(void);

	//get largest local truncation error from all the meshes
	double Get_lte(void);

	//calculate time step for adaptive methods based on current error values and evaluation method settings - in DiffEq_Iterate.cpp
	//return true if current time step is good, else false if we must repeat it
	//this uses a 2 level error threshold -> above the high threshold fail, adjust step based on max_error / error ratio. Below the low error threshold increase step by a small constant factor.
	bool SetAdaptiveTimeStep(void);
	//Adaptive time step using a single error threshold -> above it fail; adjust step based on max_error / error ratio.
	bool SetAdaptiveTimeStep_SingleThreshold(void);
	bool SetAdaptiveTimeStep_SingleThreshold2(void);

#if COMPILECUDA == 1
	bool SetAdaptiveTimeStepCUDA(void);
	bool SetAdaptiveTimeStepCUDA_SingleThreshold(void);
	bool SetAdaptiveTimeStepCUDA_SingleThreshold2(void);
#endif

public:

	ODECommon(bool called_from_derived = false);
	virtual ~ODECommon();

	//implement pure virtual method from ProgramState
	void RepairObjectState(void);

	//---------------------------------------- SET-UP METHODS : DiffEq_Common.cpp

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//switch CUDA state on/off
	BError SwitchCUDAState(bool cudaState);

	BError SetODE(ODE_ setODE_, EVAL_ evalMethod_, bool set_eval_method = true);

	BError SetEvaluationMethod(EVAL_ evalMethod_);

	void Reset(void);

	void NewStage(void);

	void SetdT(double dT);

	void SetAdaptiveTimeStepCtrl(double err_high_fail, double err_high, double err_low, double dT_increase, double dT_min, double dT_max);

	void SetMoveMeshTrigger(bool status, int meshId = -1);
	void SetMoveMeshAntisymmetric(bool antisymmetric) { moving_mesh_antisymmetric = antisymmetric; }
	void SetMoveMeshThreshold(double threshold) { moving_mesh_threshold = threshold; }

	void SetEvaluationSpeedup(int status) { if (status >= EVALSPEEDUP_NONE && status < EVALSPEEDUP_NUMENTRIES) use_evaluation_speedup = status; }

	//---------------------------------------- ITERATE METHODS : DiffEq_Iterate.cpp (and DiffEq_IterateCUDA.cpp)

	//advance time using the set ODE solver
	void Iterate(void);

#if COMPILECUDA == 1
	void IterateCUDA(void);
#endif

	//moving mesh algorithm : shift all meshes using their MoveMesh method
	void MovingMeshAlgorithm(SuperMesh* pSMesh);

	//depending on the set evaluation method, some effective fields (demag field in particular) might not need to be updated at certain steps in the evaluation - keep the previous calculated field contribution
	//Modules can call this method to check if they should be updating the field or not.
	//Return values:
	//EVALSPEEDUPSTEP_SKIP : do not update field, use previous calculation if available
	//EVALSPEEDUPSTEP_COMPUTE_NO_SAVE : update field and do not save calculation for next step (since at the next step we'll have to calculate field again so no point saving it)
	//EVALSPEEDUPSTEP_COMPUTE_AND_SAVE : update field and save calculation for next time (since at the next step we'll need to re-use calculation)
	//To enable this mode you need to set use_evaluation_speedup != EVALSPEEDUP_NONE
	int Check_Step_Update(void);

	//---------------------------------------- GET METHODS : DiffEqCommon.cpp

	int GetIteration(void) { return iteration; }
	int GetStageIteration(void) { return stageiteration; }

	double GetTime(void) { return time; }
	double GetStageTime(void) { return stagetime; }
	double GetTimeStep(void) { return dT; }
	double Get_mxh(void);
	double Get_dmdt(void);

	DBL3 Get_AStepRelErrCtrl(void) { return DBL3(err_high_fail, err_high, err_low); }
	DBL3 Get_AStepdTCtrl(void) { return DBL3(dT_increase, dT_min, dT_max); }

	bool IsMovingMeshSet(void) { return moving_mesh; }
	int GetId_of_MoveMeshTrigger(void);

	bool MoveMeshAntisymmetric(void) { return moving_mesh_antisymmetric; }
	double MoveMeshThreshold(void) { return moving_mesh_threshold; }

	double Get_dwshift(void) { return moving_mesh_dwshift; }

	bool TimeStepSolved(void) { return available; }

	bool SolveSpinCurrent(void) { return solve_spin_current; }

	int EvaluationSpeedup(void) { return use_evaluation_speedup; }

	void QueryODE(ODE_ &setODE_, EVAL_ &evalMethod_) { setODE_ = (ODE_)setODE; evalMethod_ = (EVAL_)evalMethod; }
	void QueryODE(ODE_ &setODE_) { setODE_ = (ODE_)setODE; }

#if COMPILECUDA == 1
	ODECommonCUDA* Get_pODECUDA(void) { return pODECUDA; }
#endif
};
