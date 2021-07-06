#pragma once

#include "Boris_Enums_Defs.h"

#ifdef MODULE_COMPILATION_TRANSPORT

//spin transport solver type : 

//0. no spin transport (just charge transport)

//1. normal metal (applicable to metal meshes)

//2. ferromagnetic (applicable to ferromagnetic meshes, and also dipole meshes)

//3. antiferromagnetic (to be implemented - will be applicable to antiferromagnetic meshes, currently these are set to STSOLVE_NONE)

//4. tunneling (to be implemented eventually - will be applicable to insulator meshes)

enum STSOLVE_ {

	STSOLVE_NONE,
	STSOLVE_NORMALMETAL, 
	STSOLVE_FERROMAGNETIC, 
	STSOLVE_ANTIFERROMAGNETIC,
	STSOLVE_TUNNELING
};

#endif
