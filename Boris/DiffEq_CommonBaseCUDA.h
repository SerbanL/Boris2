#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

#include "ErrorHandler.h"

class ODECommon_Base;

class ODECommon_BaseCUDA
{

private:

	//-----------------------------------CPU version pointer

	//pointer to CPU version
	static ODECommon_Base *pODEBase;

protected:

	//-----------------------------------Primary Data

	static cu_obj<cuBReal>* ptime;
	static cu_obj<cuBReal>* pstagetime;

	//-----------------------------------Time step

	//these need to be pointers, not cu_obj directly : we only want to make the cuda objects when ODECommonCUDA is made (cuda switched on), not at the start of the program - what if cuda not available on the system?
	static cu_obj<cuBReal>* pdT;
	static cu_obj<cuBReal>* pdT_last;

	//-----------------------------------mxh and dmdt

	static cu_obj<cuBReal>* pmxh;
	static cu_obj<cuReal3>* pmxh_av;
	static cu_obj<size_t>* pavpoints;

	static cu_obj<cuBReal>* pdmdt;
	static cu_obj<cuReal3>* pdmdt_av;
	static cu_obj<size_t>* pavpoints2;

	//----------------------------------Adaptive time step control

	static cu_obj<cuBReal>* plte;

	//-----------------------------------Special evaluation values

	static cu_obj<bool>* palternator;

	//-----------------------------------Special Properties

	static cu_obj<bool>* psolve_spin_current;

protected:

	//----------------------------------- SET-UP METHODS : DiffEq_CommonBaseCUDA.cpp

	//Allocate memory for all static data; deletion only happens in the destructor, however allocation can also be triggered by UpdateConfiguration since the static data can be deleted by another instance which inherits same static data
	void AllocateStaticData(void);

	//----------------------------------- Auxiliary : DiffEq_CommonBaseCUDA.cu
	
	//zero all main reduction values : mxh, dmdt, lte
	void Zero_reduction_values(void);
	void Zero_mxh_lte_values(void);
	void Zero_dmdt_lte_values(void);
	void Zero_lte_value(void);

	void mxhav_to_mxh(void);
	void dmdtav_to_dmdt(void);

	//----------------------------------- GPU <-> CPU sync : DiffEq_CommonBaseCUDA.cpp

	//set all cuda values here from their cpu values held in ODECommon
	void SyncODEValues(void);

	//set specific cuda values (used often)
	void Sync_time(void);
	void Sync_dT(void);
	void Sync_dT_last(void);
	void Sync_alternator(void);

public:

	ODECommon_BaseCUDA(void) {}
	ODECommon_BaseCUDA(ODECommon_Base *pODEBase_);

	virtual ~ODECommon_BaseCUDA();

	//---------------------------------------- GET METHODS

	cuBReal Get_mxh(void) { return pmxh->to_cpu(); }
	cuBReal Get_dmdt(void) { return pdmdt->to_cpu(); }
	cuBReal Get_lte(void) { return plte->to_cpu(); }
};

#endif
