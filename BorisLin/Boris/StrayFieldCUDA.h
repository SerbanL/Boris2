#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_STRAYFIELD

#include "ModulesCUDA.h"

#include "StrayField_BaseCUDA.h"

class SuperMesh;
class StrayField;

class StrayFieldCUDA :
	public StrayField_BaseCUDA,
	public ModulesCUDA
{

private:

	//the cpu version module holding this CUDA version
	StrayField* pStrayField;

public:

	StrayFieldCUDA(SuperMesh* pSMesh_, StrayField* pStrayField_);
	~StrayFieldCUDA();

	//-------------------Abstract base class method implementations

	void Uninitialize(void) { initialized = false; }

	BError Initialize(void);

	BError UpdateConfiguration(UPDATECONFIG_ cfgMessage);
	void UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage) {}

	void UpdateField(void);

	//-------------------Getters

	cu_obj<cuVEC<cuReal3>>& GetStrayField(void) { return strayField; }
};

#else

class StrayFieldCUDA
{
};

#endif

#endif




