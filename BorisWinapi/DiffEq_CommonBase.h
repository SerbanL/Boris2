#pragma once

#include "BorisLib.h"

#include "ErrorHandler.h"
#include "CompileFlags.h"

#include "Boris_Enums_Defs.h"

#include "DiffEq_Defs.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	"static" base class (SBC) for micromagnetic and atomistic ODE classes (i.e. all data members are static, so when inherited by different classes they all share these data and associated methods)
//
//	ODECommon_Base -> ODECommon		 -(SBC)-> DifferentialEquation      -(ABC implementation)-> DifferentialEquationFM (etc.)
//				   -> Atom_ODECommon -(SBC)-> Atom_DifferentialEquation -(ABC implementation)-> Atom_DifferentialEquationCubic (etc.)
//
// instances of ODECommon and Atom_ODECommon are held by the SuperMesh, so this level holds only supermesh level data
// implementations of the ABCs are held in individual meshes, so they all see the common SBC data and can make use of associated methods (and also modify these data)
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using namespace std;

class SuperMesh;

class ODECommon;
class Atom_ODECommon;

#if COMPILECUDA == 1
#include "DiffEq_CommonBaseCUDA.h"
#endif

class ODECommon_Base {

#if COMPILECUDA == 1
	friend ODECommon_BaseCUDA;
#endif

private:

	//-----------------------------------Pointers

	//hold pointers to next level up SBCs so we can coordinate them from here
	static ODECommon *podeSolver;
	static Atom_ODECommon *patom_odeSolver;

protected:

	//-----------------------------------Primary Data

	//parameters for iterations and stages control
	static int iteration, stageiteration;

	//time and stagetime
	static double time, stagetime;

	//-----------------------------------Time step

	//time-step
	static double dT;
	//dT_last = dT, unless an adaptive time step method is used, where dT is saved before being changed - needed to calculate dM/dt for adaptive time step methods
	static double dT_last;

	//time-step for stochastic field generation (must be greater or equal to dT)
	static double dTstoch;
	//last time value the stochastic field was generated : used in conjunction with dTstoch to trigger stochastic field generation.
	static double time_stoch;
	//if set to true, dT is used instead of dT (default)
	static bool link_dTstoch;

	//-----------------------------------Evaluation Method Data

	//flag to indicate if evaluation method has completed a full iteration
	static bool available;

	//set evaluation method for given ODE
	static int evalMethod;

	//evaluation step for higher order evaluation methods
	static int evalStep;

	//use evaluation speedup by skipping certain effective field contributions at specific steps in the evaluation?
	//this takes on a value from EVALSPEEDUP_ enum
	static int use_evaluation_speedup;

	//-----------------------------------mxh and dmdt

	//to avoid calculating mxh every single iteration whether it's needed or not, only calculate it if this flag is set
	//after updating mxh value set this flag to false again. Thus when a call to Get_mxh is received, set this flag to true.
	//this means the actual mxh value returned by Get_mxh will be off by 1 iteration but it doesn't really matter
	//if a schedule with mxh stop condition is used, all it means is the schedule will run for an extra iteration
	//saved mxh values will also be off by 1 iteration but again it doesn't matter.
	static bool calculate_mxh;

	//final mxh value after each iteration
	static double mxh;

	//as for mxh but for dmdt stopping condition
	static bool calculate_dmdt;

	//final dmdt value after each iteration
	static double dmdt;

	//-----------------------------------Adaptive time step control

	//local truncation error used to adjust time step
	static double lte;

	//parameters for adaptive time step calculation method
	//fail above this relative error : repeat with lower time step
	static double err_high_fail;

	//decrease time step above this relative error - decrease factor depends on error ratio.
	static double err_high;

	//increase time step above this relative error
	static double err_low;

	//when increasing dT multiply by this (> 1.0)
	static double dT_increase;

	//maximum and minimum dT values allowed
	static double dT_min;
	static double dT_max;

	//-----------------------------------Special evaluation values

	//used to alternate between past equation evaluations (e.g. for ABM)
	static bool alternator;

	//some evaluation methods need to be primed (e.g. ABM); 
	static bool primed;

	//-----------------------------------Moving mesh data

	//use moving mesh algorithm?
	static bool moving_mesh;

	//moving mesh algorithm type: antisymmetric for DWs, symmetric for skyrmions
	static bool moving_mesh_antisymmetric;

	//normalised threshold value used to trigger a mesh shift
	static double moving_mesh_threshold;

	//current dw shift as resulting from moving mesh algorithm
	static double moving_mesh_dwshift;

	//-----------------------------------Special Properties

	//is the currently set equation a SA version? (set by SA ODE versions, which require spin accumulation to evaluate spin torques)
	static bool solve_spin_current;

private:

	//----------------------------------- Runtime Iteration Helpers

	//calculate time step for adaptive methods based on current error values and evaluation method settings
	//return true if current time step is good, else false if we must repeat it
	//this uses a 2 level error threshold -> above the high threshold fail, adjust step based on max_error / error ratio. Below the low error threshold increase step by a small constant factor.
	bool SetAdaptiveTimeStep(void);

protected:
	
	//----------------------------------- Runtime Iteration Helpers

	//restore ODE solvers so iteration can be redone (call if AdvanceIteration has failed)
	virtual void Restore(void) = 0;
	
#if COMPILECUDA == 1
	virtual void RestoreCUDA(void) = 0;
#endif

	//----------------------------------- CTOR/DTOR

	ODECommon_Base() {}
	virtual ~ODECommon_Base() {}

public:

	//----------------------------------- Setup

	//once ODECommon and Atom_ODECommon have been constructed, this function should be called to set pointers here (in SuperMesh constructor).
	void set_pointers(ODECommon& odeSolver, Atom_ODECommon& atom_odeSolver);

	//----------------------------------- Important Control Methods : DiffEq_CommonBase.cpp

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//switch CUDA state on/off
	BError SwitchCUDAState(bool cudaState);

	//----------------------------------- Runtime Iteration : DiffEq_CommonBase_Iterate.cpp

	//use if both types of solvers are active (also works if just one is active): micromagnetic and atomistic
	void Iterate(void);

#if COMPILECUDA == 1
	void IterateCUDA(void);
#endif

	//depending on the set evaluation method, some effective fields (demag field in particular) might not need to be updated at certain steps in the evaluation - keep the previous calculated field contribution
	//Modules can call this method to check if they should be updating the field or not.
	//Return values:
	//EVALSPEEDUPSTEP_SKIP : do not update field, use previous calculation if available
	//EVALSPEEDUPSTEP_COMPUTE_NO_SAVE : update field and do not save calculation for next step (since at the next step we'll have to calculate field again so no point saving it)
	//EVALSPEEDUPSTEP_COMPUTE_AND_SAVE : update field and save calculation for next time (since at the next step we'll need to re-use calculation)
	//To enable this mode you need to set use_evaluation_speedup != EVALSPEEDUP_NONE
	int Check_Step_Update(void);

	//----------------------------------- Evaluation Method and Control: DiffEq_CommonBase_Control.cpp

	BError SetEvaluationMethod(EVAL_ evalMethod_);

	void Reset(void);

	void NewStage(void);

	void SetdT(double dT);

	void SetdTstoch(double dTstoch);
	void SetLink_dTstoch(bool link_dTstoch);

	void SetAdaptiveTimeStepCtrl(double err_high_fail, double err_high, double err_low, double dT_increase, double dT_min, double dT_max);

	void SetEvaluationSpeedup(int status) { if (status >= EVALSPEEDUP_NONE && status < EVALSPEEDUP_NUMENTRIES) use_evaluation_speedup = status; }

	//----------------------------------- Moving Mesh Methods : DiffEq_CommonBase_MovingMesh.cpp

	void SetMoveMeshTrigger(bool status, int meshId = -1);
	void SetMoveMeshAntisymmetric(bool antisymmetric) { moving_mesh_antisymmetric = antisymmetric; }
	void SetMoveMeshThreshold(double threshold) { moving_mesh_threshold = threshold; }

	//moving mesh algorithm : shift all meshes using their MoveMesh method
	void MovingMeshAlgorithm(SuperMesh* pSMesh);

	bool IsMovingMeshSet(void) { return moving_mesh; }
	int GetId_of_MoveMeshTrigger(void);

	bool MoveMeshAntisymmetric(void) { return moving_mesh_antisymmetric; }
	double MoveMeshThreshold(void) { return moving_mesh_threshold; }

	double Get_dwshift(void) { return moving_mesh_dwshift; }

	//----------------------------------- Get Primary Data

	int GetIteration(void) { return iteration; }
	int GetStageIteration(void) { return stageiteration; }

	double GetTime(void) { return time; }
	double GetStageTime(void) { return stagetime; }
	double GetTimeStep(void) { return dT; }

	double GetStochTimeStep(void) { return (link_dTstoch ? dT : dTstoch); }
	bool GetLink_dTstoch(void) { return link_dTstoch; }

	DBL3 Get_AStepRelErrCtrl(void) { return DBL3(err_high_fail, err_high, err_low); }
	DBL3 Get_AStepdTCtrl(void) { return DBL3(dT_increase, dT_min, dT_max); }

	//----------------------------------- Status Getters

	bool TimeStepSolved(void) { return available; }

	bool SolveSpinCurrent(void) { return solve_spin_current; }

	int EvaluationSpeedup(void) { return use_evaluation_speedup; }

	//----------------------------------- Value Getters

	double Get_mxh(void);
	double Get_dmdt(void);
};
