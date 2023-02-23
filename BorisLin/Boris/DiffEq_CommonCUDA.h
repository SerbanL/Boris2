#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

#include "MeshParamsCUDA.h"

#include "ManagedDiffEq_CommonCUDA.h"

#include "DiffEq_CommonBaseCUDA.h"

class ODECommon;
class ODECommon_Base;

class ODECommonCUDA :
	public ODECommon_BaseCUDA
{
	friend ODECommon;
	friend ODECommon_Base;
	friend ManagedDiffEq_CommonCUDA;
	
private:

	//-----------------------------------CPU version pointer

	//pointer to CPU version
	static ODECommon *pODE;

	//-----------------------------------Managed DiffEq

	//ManagedDiffEq_CommonCUDA holds pointers to data in ODECommonCUDA in an object in gpu memory.
	//pass cuDiffEq to a cuda kernel then all gpu data held here in cu_obj objects can be accessed in device code.
	//Initialize ManagedDiffEq_CommonCUDA with all the pointers you need then forget about it - no book-keeping required.
	static cu_obj<ManagedDiffEq_CommonCUDA>* pcuDiffEq;

protected:

	//-----------------------------------Equation

	//set equation identifier in GPU memory
	static cu_obj<int>* psetODE;

	//-----------------------------------Evaluation method modifiers

	static cu_obj<bool>* prenormalize;

	//-----------------------------------Steepest Descent Solver

	//quantities used to calculate Barzilai-Borwein stepsizes across multiple meshes
	//Accumulate values in these quantities, then obtain stepsizes as:
	//step1 = delta_m_sq / delta_m_dot_delta_G
	//step2 = delta_m_dot_delta_G / delta_G_sq
	static cu_obj<cuBReal>* pdelta_M_sq;
	static cu_obj<cuBReal>* pdelta_G_sq;
	static cu_obj<cuBReal>* pdelta_M_dot_delta_G;

	static cu_obj<cuBReal>* pdelta_M2_sq;
	static cu_obj<cuBReal>* pdelta_G2_sq;
	static cu_obj<cuBReal>* pdelta_M2_dot_delta_G2;

private:

	//----------------------------------- SET-UP METHODS : DiffEq_CommonCUDA.cpp

	//Allocate memory for all static data; deletion only happens in the destructor, however allocation can also be triggered by UpdateConfiguration since the static data can be deleted by another instance which inherits same static data
	void AllocateStaticData(void);

	//----------------------------------- Auxiliary : DiffEq_CommonCUDA.cu, DiffEq_CommonCUDA.cpp

	//specific to SD solver
	void Zero_SD_Solver_BB_Values(void);

	void Get_SD_Solver_BB_Values(
		double* pdelta_M_sq_cpu, double* pdelta_G_sq_cpu, double* pdelta_M_dot_delta_G_cpu,
		double* pdelta_M2_sq_cpu, double* pdelta_G2_sq_cpu, double* pdelta_M2_dot_delta_G2_cpu);

	//----------------------------------- GPU <-> CPU sync : DiffEq_CommonCUDA.cpp

	//set all cuda values here from their cpu values held in ODECommon
	void SyncODEValues(void);

public:

	ODECommonCUDA(void) {}
	ODECommonCUDA(ODECommon *pODE_);

	virtual ~ODECommonCUDA();

	//----------------------------------- Important Control Methods

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	//----------------------------------- GETTERS

	//get reference to stored managed cuda differential equation object (cuDiffEq)
	cu_obj<ManagedDiffEq_CommonCUDA>& Get_ManagedDiffEq_CommonCUDA(void) { return *pcuDiffEq; }

};

#endif

