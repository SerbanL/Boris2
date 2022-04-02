#pragma once

//spin transport solver type : 

//0. no spin transport (just charge transport)

//1. normal metal (applicable to metal meshes)

//2. ferromagnetic (applicable to ferromagnetic meshes, and also dipole meshes)

//3. antiferromagnetic (to be implemented - will be applicable to antiferromagnetic meshes, currently these are set to STSOLVE_NONE)

//4. tunneling (to be implemented eventually - will be applicable to insulator meshes)

//5. ferromagnetic in atomistic mesh

enum STSOLVE_ {

	STSOLVE_NONE,
	STSOLVE_NORMALMETAL, 
	STSOLVE_FERROMAGNETIC, 
	STSOLVE_ANTIFERROMAGNETIC,
	STSOLVE_TUNNELING,
	STSOLVE_FERROMAGNETIC_ATOM
};

enum TMR_ {

	TMR_NONE = -1,

	//simple (1-cos) dependence of resistance
	TMR_COS = 0,
	//cos dependence of conductance
	TMR_SLONCZEWSKI = 1,

	//number of options in this enum
	TMR_NUMOPTIONS
};