#pragma once

//Includes ProgramState, allowing saving to/from a file of current program state.

#include "ProgramState.h"

//CIRCULAR INCLUSION CHECK : PASSED 

/*

#include "Types_VAL.h"
-> Funcs_Conv
	-> Introspection_base
	-> Types_Conversion
	-> Introspection_base
	-> Funcs_Math_base
-> Funcs_Strings
	-> Funcs_Vectors

#include "Introspection.h"
-> Introspection_base
-> Types_VAL
	-> Funcs_Conv
		-> Introspection_base
		-> Types_Conversion
		-> Introspection_base
		-> Funcs_Math_base
	-> Funcs_Strings
		-> Funcs_Vectors
	-> Funcs_Aux_base

#include "Funcs_Vectors.h"
-> Funcs_Conv
	-> Types_Conversion
		-> Introspection_base
		-> Funcs_Math_base

*/