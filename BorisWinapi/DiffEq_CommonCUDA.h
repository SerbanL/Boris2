#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "MeshParamsCUDA.h"

#include "ManagedDiffEq_CommonCUDA.h"

//NOTES : renormalize at last step for all equations except all LLB versions

class ODECommon;

class ODECommonCUDA
{
	friend ODECommon;
	friend ManagedDiffEq_CommonCUDA;

private:

	static ODECommon *pODE;

	//ManagedDiffEq_CommonCUDA holds pointers to data in ODECommonCUDA in an object in gpu memory.
	//pass cuDiffEq to a cuda kernel then all gpu data held here in cu_obj objects can be accessed in device code.
	//Initialize ManagedDiffEq_CommonCUDA with all the pointers you need then forget about it - no book-keeping required.
	cu_obj<ManagedDiffEq_CommonCUDA> cuDiffEq;

protected:

	//-----------------------------------Time and Stage Time

	static cu_obj<cuBReal>* ptime;
	static cu_obj<cuBReal>* pstagetime;

	//-----------------------------------Time step

	//these need to be pointers, not cu_obj directly : we only want to make the cuda objects when ODECommonCUDA is made (cuda switched on), not at the start of the program - what if cuda not available on the system?
	static cu_obj<cuBReal>* pdT;
	static cu_obj<cuBReal>* pdT_last;

	//-----------------------------------Primary data

	static cu_obj<cuBReal>* pmxh;
	static cu_obj<cuReal3>* pmxh_av;
	static cu_obj<size_t>* pavpoints;

	static cu_obj<cuBReal>* pdmdt;
	static cu_obj<cuReal3>* pdmdt_av;
	static cu_obj<size_t>* pavpoints2;

	static cu_obj<cuBReal>* plte;

	//-----------------------------------Evaluation method modifiers

	static cu_obj<bool>* prenormalize;

	//-----------------------------------Properties flags

	static cu_obj<bool>* psolve_spin_current;

	//-----------------------------------Equation and Evaluation method values

	static cu_obj<int>* psetODE;

	//-----------------------------------Special values

	static cu_obj<bool>* palternator;

	//-----------------------------------Steepest Descent Solver

	//quantities used to calculate Barzilai-Borwein stepsizes across multiple meshes
	//Accumulate values in these quantities, then obtain stepsizes as:
	//step1 = delta_M_sq / delta_M_dot_delta_G
	//step2 = delta_M_dot_delta_G / delta_G_sq
	static cu_obj<cuBReal>* pdelta_M_sq;
	static cu_obj<cuBReal>* pdelta_G_sq;
	static cu_obj<cuBReal>* pdelta_M_dot_delta_G;

	static cu_obj<cuBReal>* pdelta_M2_sq;
	static cu_obj<cuBReal>* pdelta_G2_sq;
	static cu_obj<cuBReal>* pdelta_M2_dot_delta_G2;

private:

	//---------------------------------------- SET-UP METHODS : DiffEqCUDA.cpp, DiffEq_EvalsCUDA.cu

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//zero all main reduction values : mxh, dmdt, lte
	void Zero_reduction_values(void);
	void Zero_mxh_lte_values(void);
	void Zero_dmdt_lte_values(void);
	void Zero_lte_value(void);

	void mxhav_to_mxh(void);
	void dmdtav_to_dmdt(void);

	//set all cuda values here from their cpu values held in ODECommon
	void SyncODEValues(void);

	//set specific cuda values (used often)
	void Sync_time(void);
	void Sync_dT(void);
	void Sync_dT_last(void);
	void Sync_alternator(void);

	//specific to SD solver
	void Zero_SD_Solver_BB_Values(void);

	void Get_SD_Solver_BB_Values(
		double* pdelta_M_sq_cpu, double* pdelta_G_sq_cpu, double* pdelta_M_dot_delta_G_cpu,
		double* pdelta_M2_sq_cpu, double* pdelta_G2_sq_cpu, double* pdelta_M2_dot_delta_G2_cpu)
	{
		*pdelta_M_sq_cpu = pdelta_M_sq->to_cpu();
		*pdelta_G_sq_cpu = pdelta_G_sq->to_cpu();
		*pdelta_M_dot_delta_G_cpu = pdelta_M_dot_delta_G->to_cpu();

		*pdelta_M2_sq_cpu = pdelta_M2_sq->to_cpu();
		*pdelta_G2_sq_cpu = pdelta_G2_sq->to_cpu();
		*pdelta_M2_dot_delta_G2_cpu = pdelta_M2_dot_delta_G2->to_cpu();
	}

public:

	ODECommonCUDA(void) {}
	ODECommonCUDA(ODECommon *pODE_);

	virtual ~ODECommonCUDA();

	//---------------------------------------- GET METHODS

	cuBReal Get_mxh(void) { return pmxh->to_cpu(); }
	cuBReal Get_dmdt(void) { return pdmdt->to_cpu(); }
	cuBReal Get_lte(void) { return plte->to_cpu(); }

	//----------------------------------- GETTERS

	//get reference to stored managed cuda differential equation object (cuDiffEq)
	cu_obj<ManagedDiffEq_CommonCUDA>& Get_ManagedDiffEq_CommonCUDA(void) { return cuDiffEq; }

};

#endif

