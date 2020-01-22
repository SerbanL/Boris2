#pragma once

#include "cuBRealComplex.h"
#include "cuFuncs_Aux.h"
#include "cuFuncs_Math.h"

#include "TEquation_FSPEC.h"

template <typename ... BVarType>
class ManagedFunctionCUDA {

	typedef cuBReal(ManagedFunctionCUDA::*pFunc)(BVarType...);

public:

	int varlevel;

	//some functions need a parameter
	cuBReal param;

	//for FUNC_POWER_EXPCONST or FUNC_POWER_BASECONST functions this is the base or exponent value.
	//cannot set it in param since that could be used with FUNC_POWER_EXPCONST_PMUL or FUNC_POWER_BASECONST_PMUL
	cuBReal base_or_exponent;

	//further functions for unary functions and binary operators
	ManagedFunctionCUDA<BVarType...>* pFunc1;
	ManagedFunctionCUDA<BVarType...>* pFunc2;

	pFunc func;

public:

	//BASIS FUNCTIONS

	__device__ cuBReal __F_bvar(int varlevelsearch, cuBReal bvar)
	{
		if (varlevel == varlevelsearch) {

			return bvar;
		}
		else return 0.0;
	}

	template <typename ... __BVarType>
	__device__ cuBReal __F_bvar(int varlevelsearch, cuBReal bvar, __BVarType... bvars)
	{
		if (varlevel == varlevelsearch) {

			return bvar;
		}
		else {

			return __F_bvar(++varlevelsearch, bvars...);
		}
	}

	//variable index 0
	template <typename ... __BVarType>
	__device__ cuBReal get_0(const cuBReal& bvar0, __BVarType... bvars) const { return bvar0; }
	__device__ cuBReal get_0(const cuBReal& bvar0) const { return bvar0; }
	__device__ cuBReal get_0(...) const { return 0.0; }

	//variable index 1
	template <typename ... __BVarType>
	__device__ cuBReal get_1(const cuBReal& bvar0, const cuBReal& bvar1, __BVarType... bvars) const { return bvar1; }
	__device__ cuBReal get_1(const cuBReal& bvar0, const cuBReal& bvar1) const { return bvar1; }
	__device__ cuBReal get_1(...) const { return 0.0; }

	//variable index 2
	template <typename ... __BVarType>
	__device__ cuBReal get_2(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, __BVarType... bvars) const { return bvar2; }
	__device__ cuBReal get_2(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2) const { return bvar2; }
	__device__ cuBReal get_2(...) const { return 0.0; }

	//variable index 3
	template <typename ... __BVarType>
	__device__ cuBReal get_3(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, const cuBReal& bvar3, __BVarType... bvars) const { return bvar3; }
	__device__ cuBReal get_3(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, const cuBReal& bvar3) const { return bvar3; }
	__device__ cuBReal get_3(...) const { return 0.0; }

	//variable index 4
	template <typename ... __BVarType>
	__device__ cuBReal get_4(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, const cuBReal& bvar3, const cuBReal& bvar4, __BVarType... bvars) const { return bvar4; }
	__device__ cuBReal get_4(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, const cuBReal& bvar3, const cuBReal& bvar4) const { return bvar4; }
	__device__ cuBReal get_4(...) const { return 0.0; }

	//variable index 5
	template <typename ... __BVarType>
	__device__ cuBReal get_5(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, const cuBReal& bvar3, const cuBReal& bvar4, const cuBReal& bvar5, __BVarType... bvars) const { return bvar5; }
	__device__ cuBReal get_5(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, const cuBReal& bvar3, const cuBReal& bvar4, const cuBReal& bvar5) const { return bvar5; }
	__device__ cuBReal get_5(...) const { return 0.0; }

	//variable index 6
	template <typename ... __BVarType>
	__device__ cuBReal get_6(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, const cuBReal& bvar3, const cuBReal& bvar4, const cuBReal& bvar5, const cuBReal& bvar6, __BVarType... bvars) const { return bvar6; }
	__device__ cuBReal get_6(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, const cuBReal& bvar3, const cuBReal& bvar4, const cuBReal& bvar5, const cuBReal& bvar6) const { return bvar6; }
	__device__ cuBReal get_6(...) const { return 0.0; }

	//variable index 7
	template <typename ... __BVarType>
	__device__ cuBReal get_7(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, const cuBReal& bvar3, const cuBReal& bvar4, const cuBReal& bvar5, const cuBReal& bvar6, const cuBReal& bvar7, __BVarType... bvars) const { return bvar7; }
	__device__ cuBReal get_7(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, const cuBReal& bvar3, const cuBReal& bvar4, const cuBReal& bvar5, const cuBReal& bvar6, const cuBReal& bvar7) const { return bvar7; }
	__device__ cuBReal get_7(...) const { return 0.0; }

	//variable index 8
	template <typename ... __BVarType>
	__device__ cuBReal get_8(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, const cuBReal& bvar3, const cuBReal& bvar4, const cuBReal& bvar5, const cuBReal& bvar6, const cuBReal& bvar7, const cuBReal& bvar8, __BVarType... bvars) const { return bvar8; }
	__device__ cuBReal get_8(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, const cuBReal& bvar3, const cuBReal& bvar4, const cuBReal& bvar5, const cuBReal& bvar6, const cuBReal& bvar7, const cuBReal& bvar8) const { return bvar8; }
	__device__ cuBReal get_8(...) const { return 0.0; }

	//variable index 9
	template <typename ... __BVarType>
	__device__ cuBReal get_9(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, const cuBReal& bvar3, const cuBReal& bvar4, const cuBReal& bvar5, const cuBReal& bvar6, const cuBReal& bvar7, const cuBReal& bvar8, const cuBReal& bvar9, __BVarType... bvars) const { return bvar9; }
	__device__ cuBReal get_9(const cuBReal& bvar0, const cuBReal& bvar1, const cuBReal& bvar2, const cuBReal& bvar3, const cuBReal& bvar4, const cuBReal& bvar5, const cuBReal& bvar6, const cuBReal& bvar7, const cuBReal& bvar8, const cuBReal& bvar9) const { return bvar9; }
	__device__ cuBReal get_9(...) const { return 0.0; }

	//user variable : search and return the correct one depending on varlevel value
	__device__ cuBReal F_bvar(BVarType... bvars)
	{
		//hard-coded up to 10 variables; anything above this (really!?) you'll need function recursion (slower) to get them.
		switch (varlevel) {
			
		case 0:
			return get_0(bvars...);

		case 1:
			return get_1(bvars...);

		case 2:
			return get_2(bvars...);

		case 3:
			return get_3(bvars...);

		case 4:
			return get_4(bvars...);

		case 5:
			return get_5(bvars...);

		case 6:
			return get_6(bvars...);

		case 7:
			return get_7(bvars...);

		case 8:
			return get_8(bvars...);

		case 9:
			return get_9(bvars...);
			
		default:
			return __F_bvar(10, bvars...);
		}
	}

	__device__ cuBReal F_bvar_pmul(BVarType... bvars)
	{
		//hard-coded up to 10 variables; anything above this (really!?) you'll need function recursion (slower) to get them.
		switch (varlevel) {
			
		case 0:
			return param * get_0(bvars...);

		case 1:
			return param * get_1(bvars...);

		case 2:
			return param * get_2(bvars...);

		case 3:
			return param * get_3(bvars...);

		case 4:
			return param * get_4(bvars...);

		case 5:
			return param * get_5(bvars...);

		case 6:
			return param * get_6(bvars...);

		case 7:
			return param * get_7(bvars...);

		case 8:
			return param * get_8(bvars...);

		case 9:
			return param * get_9(bvars...);
			
		default:
			return param * __F_bvar(10, bvars...);
		}
	}

	//constant
	__device__ cuBReal F_const(BVarType... bvars) { return param; }

	//UNARY FUNCTIONS

	__device__ cuBReal F_sin(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return sin((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_sin_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * sin((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_cos(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return cos((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_cos_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * cos((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_tan(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return tan((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_tan_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * tan((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_sinh(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return sinh((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_sinh_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * sinh((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_cosh(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return cosh((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_cosh_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * cosh((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_tanh(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return tanh((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_tanh_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * tanh((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	//square root
	__device__ cuBReal F_sqrt(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return sqrt((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_sqrt_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * sqrt((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_exp(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return exp((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_exp_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * exp((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_asin(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return asin((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_asin_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * asin((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_acos(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return acos((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_acos_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * acos((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_atan(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return atan((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_atan_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * atan((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_asinh(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return asinh((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_asinh_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * asinh((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_acosh(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return acosh((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_acosh_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * acosh((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_atanh(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return atanh((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_atanh_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * atanh((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	//natual logarithm
	__device__ cuBReal F_ln(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return log((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_ln_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * log((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	//logarithm base 10
	__device__ cuBReal F_log(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return log10((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_log_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * log10((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	//modulus
	__device__ cuBReal F_abs(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return fabs((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_abs_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * fabs((cuFunc1.*(cuFunc1.func))(bvars...));
	}

	//step function: step(t) = 0 for t < 0, = 1 to t>= 0
	__device__ cuBReal F_step(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		if ((cuFunc1.*(cuFunc1.func))(bvars...) < 0) return 0.0;
		else return 1.0;
	}

	__device__ cuBReal F_step_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		if ((cuFunc1.*(cuFunc1.func))(bvars...) < 0) return 0.0;
		else return param;
	}

	//square wave period 2*PI range -1 to 1 : swav(0+) = +1, swav(0-) = -1
	__device__ cuBReal F_swav(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		cuBReal value = (cuFunc1.*(cuFunc1.func))(bvars...);

		return -2 * (((int)floor(fabs(value) / (cuBReal)PI) + (cu_get_sign(value) < 0)) % 2) + 1;
	}

	__device__ cuBReal F_swav_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		cuBReal value = (cuFunc1.*(cuFunc1.func))(bvars...);

		return param * (-2 * (((int)floor(fabs(value) / (cuBReal)PI) + (cu_get_sign(value) < 0)) % 2) + 1);
	}

	//triangular wave period 2*PI range -1 to 1 : twav(0 to PI) = +1 to -1, twav(-PI to 0) = -1 to +1
	__device__ cuBReal F_twav(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		cuBReal value = (cuFunc1.*(cuFunc1.func))(bvars...);

		if ((int)floor(fabs(value) / (cuBReal)PI) % 2) return (2 * fmod(fabs(value), (cuBReal)PI) / (cuBReal)PI - 1);
		else return (1 - 2 * fmod(fabs(value), (cuBReal)PI) / (cuBReal)PI);
	}

	__device__ cuBReal F_twav_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		cuBReal value = (cuFunc1.*(cuFunc1.func))(bvars...);

		if ((int)floor(fabs(value) / (cuBReal)PI) % 2) return param * (2 * fmod(fabs(value), (cuBReal)PI) / (cuBReal)PI - 1);
		else return param * (1 - 2 * fmod(fabs(value), (cuBReal)PI) / (cuBReal)PI);
	}

	//power function : constant exponent
	__device__ cuBReal F_pow_expconst(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return pow((cuFunc1.*(cuFunc1.func))(bvars...), base_or_exponent);
	}

	__device__ cuBReal F_pow_expconst_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * pow((cuFunc1.*(cuFunc1.func))(bvars...), base_or_exponent);
	}

	//power function : constant base
	__device__ cuBReal F_pow_baseconst(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return pow(base_or_exponent, (cuFunc1.*(cuFunc1.func))(bvars...));
	}

	__device__ cuBReal F_pow_baseconst_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		return param * pow(base_or_exponent, (cuFunc1.*(cuFunc1.func))(bvars...));
	}

	//BINARY FUNCTIONS
	__device__ cuBReal F_add(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		ManagedFunctionCUDA<BVarType...>& cuFunc2 = *pFunc2;

		return (cuFunc1.*(cuFunc1.func))(bvars...) + (cuFunc2.*(cuFunc2.func))(bvars...);
	}

	__device__ cuBReal F_add_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		ManagedFunctionCUDA<BVarType...>& cuFunc2 = *pFunc2;

		return param * ((cuFunc1.*(cuFunc1.func))(bvars...) + (cuFunc2.*(cuFunc2.func))(bvars...));
	}

	__device__ cuBReal F_sub(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		ManagedFunctionCUDA<BVarType...>& cuFunc2 = *pFunc2;

		return (cuFunc1.*(cuFunc1.func))(bvars...) - (cuFunc2.*(cuFunc2.func))(bvars...);
	}

	__device__ cuBReal F_sub_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		ManagedFunctionCUDA<BVarType...>& cuFunc2 = *pFunc2;

		return param * ((cuFunc1.*(cuFunc1.func))(bvars...) - (cuFunc2.*(cuFunc2.func))(bvars...));
	}

	__device__ cuBReal F_mul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		ManagedFunctionCUDA<BVarType...>& cuFunc2 = *pFunc2;

		return (cuFunc1.*(cuFunc1.func))(bvars...) * (cuFunc2.*(cuFunc2.func))(bvars...);
	}

	__device__ cuBReal F_mul_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		ManagedFunctionCUDA<BVarType...>& cuFunc2 = *pFunc2;

		return param * (cuFunc1.*(cuFunc1.func))(bvars...) * (cuFunc2.*(cuFunc2.func))(bvars...);
	}

	__device__ cuBReal F_div(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		ManagedFunctionCUDA<BVarType...>& cuFunc2 = *pFunc2;

		return (cuFunc1.*(cuFunc1.func))(bvars...) / (cuFunc2.*(cuFunc2.func))(bvars...);
	}

	__device__ cuBReal F_div_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		ManagedFunctionCUDA<BVarType...>& cuFunc2 = *pFunc2;

		return param * ((cuFunc1.*(cuFunc1.func))(bvars...) / (cuFunc2.*(cuFunc2.func))(bvars...));
	}

	__device__ cuBReal F_pow(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		ManagedFunctionCUDA<BVarType...>& cuFunc2 = *pFunc2;

		return pow((cuFunc1.*(cuFunc1.func))(bvars...), (cuFunc2.*(cuFunc2.func))(bvars...));
	}

	__device__ cuBReal F_pow_pmul(BVarType... bvars)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc1 = *pFunc1;
		ManagedFunctionCUDA<BVarType...>& cuFunc2 = *pFunc2;

		return param * pow((cuFunc1.*(cuFunc1.func))(bvars...), (cuFunc2.*(cuFunc2.func))(bvars...));
	}

public:

	////////////////////////////////////////////////

	__host__ void construct_cu_obj(void)
	{
		nullgpuptr(pFunc1);
		nullgpuptr(pFunc2);
		set_gpu_value(param, (cuBReal)1.0);
		set_gpu_value(varlevel, (int)0);
		set_gpu_value(base_or_exponent, (cuBReal)0);
	}

	__host__ void destruct_cu_obj(void)
	{
	}

	////////////////////////////////////////////////

	//BASIS FUNCTIONS
	//Setup this object as a basis function
	__host__ void Function_Basis(EqComp::FSPEC fspec);

	//UNARY FUNCTIONS
	//Setup this object as a unary function on cuFunc1
	__host__ void Function_Unary(EqComp::FSPEC fspec, ManagedFunctionCUDA<BVarType...>& cuFunc1);

	//BINARY FUNCTIONS
	//Setup this object as a binary operator on cuFunc1 and cuFunc2
	__host__ void Operator_Binary(EqComp::FSPEC fspec, ManagedFunctionCUDA<BVarType...>& cuFunc1, ManagedFunctionCUDA<BVarType...>& cuFunc2);

	////////////////////////////////////////////////

	__device__ cuBReal evaluate(BVarType ... bvar)
	{
		ManagedFunctionCUDA<BVarType...>& cuFunc = *this;

		return (cuFunc.*(cuFunc.func))(bvar...);
	}
};

