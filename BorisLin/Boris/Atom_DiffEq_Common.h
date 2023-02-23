#pragma once

#include "DiffEq_CommonBase.h"

class Atom_Mesh;

class Atom_DifferentialEquation;
typedef DBL3(Atom_DifferentialEquation::*Atom_Equation)(int idx);

#if COMPILECUDA == 1
#include "Atom_DiffEq_CommonCUDA.h"
#endif

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	"static" base class for DifferentialEquation objects (i.e. all data members are static, so when inherited by different classes they all share these data and associated methods)
//
//	Each DifferentialEquation type object is held in its atomistic mesh and has different internal storage dimensions depending on the mesh dimensions
//	There are a number of properties that are in common between all the DifferentialEquation objects, since the same equation and evaluation method must be used for all the atomistic meshes.
//	Moreover all the DifferentialEquation objects must be coordinated by a top-level object : ODECommon also provides this by keeping track of all of them 
//	There's another level up, which is the ODECommon_Base (also a "static" base class), and this is shared by micromagnetic and atomistic meshes
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Atom_ODECommon :
	public ProgramState<Atom_ODECommon,
	std::tuple<
	int, int, 
	double, double, 
	double, double, int, int, 
	double, double, double, bool,
	double, double, bool,
	double, double, double, double, double, double, int,
	int,
	bool, bool, double, double>,
	std::tuple<>>,
	public ODECommon_Base
{
	friend ODECommon_Base;

#if COMPILECUDA == 1
	friend Atom_ODECommonCUDA;
#endif

protected:

	//-----------------------------------Collection of atomistic magnetic meshes

	//collection of all set ODE in ferromagnetic meshes : this is handled by the DifferentialEquation constructor and destructor
	//When a new ODE is made in a ferromagnetic mesh, DifferentialEquation constructor is called and it pushes a pointer in this vector, together with a unique id for it.
	//When a ODE is deleted in a ferromagnetic mesh, DifferentialEquation destructor is called, and it erases the entry n this vector using the unique odeId previously generated.
	static vector_lut<Atom_DifferentialEquation*> pODE;

	//-----------------------------------Equation and Evaluation method values

	//currently set evaluation method
	static int setODE;

	//function pointer to equation to solve
	static Atom_Equation equation;

	//-----------------------------------Evaluation method modifiers

	//only renormalize M after solving for a time step for some equations. For LLB type equations must not renormalize.
	static bool renormalize;

	//-----------------------------------Steepest Descent Solver

	//quantities used to calculate Barzilai-Borwein stepsizes across multiple meshes
	//Accumulate values in these quantities, then obtain stepsizes as:
	//step1 = delta_m_sq / delta_m_dot_delta_G
	//step2 = delta_m_dot_delta_G / delta_G_sq
	static double delta_m_sq, delta_m2_sq;
	static double delta_G_sq, delta_G2_sq;
	static double delta_m_dot_delta_G, delta_m2_dot_delta_G2;

#if COMPILECUDA == 1
	static Atom_ODECommonCUDA *pODECUDA;
#endif

private:

	//out of all currently set meshes find maximum mxh value and set it in mxh. This is called after all meshes have iterated, so their respective mxh values are in mxh_reduction[0]
	void Set_mxh(void);

	//out of all currently set meshes find maximum dmdt value and set it in dmdt. This is called after all meshes have iterated, so their respective mxh values are in dmdt_reduction[0]
	void Set_dmdt(void);

	//set largest local truncation error from all the meshes in the static lte
	void Set_lte(void);

public:

	Atom_ODECommon(bool called_from_derived = false);
	virtual ~Atom_ODECommon();

	//implement pure virtual method from ProgramState
	void RepairObjectState(void);

	//---------------------------------------- SET-UP METHODS : DiffEq_Common.cpp

	BError SetODE(ODE_ setODE_, EVAL_ evalMethod_, bool set_eval_method = true);

	//---------------------------------------- ITERATE METHODS : Atom_DiffEq_Iterate.cpp (and Atom_DiffEq_IterateCUDA.cpp)

	//restore ODE solvers so iteration can be redone (call if AdvanceIteration has failed)
	void Restore(void);

#if COMPILECUDA == 1
	void RestoreCUDA(void);
#endif

	//---------------------------------------- GET METHODS : DiffEqCommon.cpp

	//get number of meshes which use this ODE solver
	int size(void) { return pODE.size(); }

	void QueryODE(ODE_ &setODE_, EVAL_ &evalMethod_) { setODE_ = (ODE_)setODE; evalMethod_ = (EVAL_)evalMethod; }
	void QueryODE(ODE_ &setODE_) { setODE_ = (ODE_)setODE; }

#if COMPILECUDA == 1
	Atom_ODECommonCUDA* Get_pODECUDA(void) { return pODECUDA; }
#endif
};