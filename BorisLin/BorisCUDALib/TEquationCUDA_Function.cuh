#include "TEquationCUDA.h"

/////////////////////////////////////////////////////////////////////////////////////////////

//BASIS FUNCTIONS
//Setup this object as a basis function

template <typename ... BVarType>
__global__ void set_ManagedFunctionCUDA_Basis(ManagedFunctionCUDA<BVarType...>& cuFunc, int type)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx == 0) {

		switch (type) {

		case EqComp::FUNC_BVAR:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_bvar;
			break;

		case EqComp::FUNC_BVAR_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_bvar_pmul;
			break;

		case EqComp::FUNC_CONST:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_const;
			break;
		}
	}
}

template <typename ... BVarType>
__host__ void ManagedFunctionCUDA<BVarType...>::Function_Basis(EqComp::FSPEC fspec)
{
	set_gpu_value(param, (cuBReal)fspec.param);
	set_gpu_value(varlevel, fspec.varlevel);

	set_ManagedFunctionCUDA_Basis<BVarType...> << <1, 1 >> > (*this, fspec.type);
}

/////////////////////////////////////////////////////////////////////////////////////////////

//UNARY FUNCTIONS
//Setup this object as a unary function on cuFunc1

template <typename ... BVarType>
__global__ void set_ManagedFunctionCUDA_Unary(ManagedFunctionCUDA<BVarType...>& cuFunc, ManagedFunctionCUDA<BVarType...>& cuFunc1, int type)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx == 0) {

		cuFunc.pFunc1 = &cuFunc1;

		switch (type) {

		case EqComp::FUNC_SIN:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_sin;
			break;

		case EqComp::FUNC_SIN_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_sin_pmul;
			break;

		case EqComp::FUNC_SINC:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_sinc;
			break;

		case EqComp::FUNC_SINC_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_sinc_pmul;
			break;

		case EqComp::FUNC_COS:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_cos;
			break;

		case EqComp::FUNC_COS_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_cos_pmul;
			break;

		case EqComp::FUNC_TAN:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_tan;
			break;

		case EqComp::FUNC_TAN_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_tan_pmul;
			break;

		case EqComp::FUNC_SINH:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_sinh;
			break;

		case EqComp::FUNC_SINH_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_sinh_pmul;
			break;

		case EqComp::FUNC_COSH:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_cosh;
			break;

		case EqComp::FUNC_COSH_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_cosh_pmul;
			break;

		case EqComp::FUNC_TANH:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_tanh;
			break;

		case EqComp::FUNC_TANH_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_tanh_pmul;
			break;

		case EqComp::FUNC_SQRT:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_sqrt;
			break;

		case EqComp::FUNC_SQRT_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_sqrt_pmul;
			break;

		case EqComp::FUNC_EXP:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_exp;
			break;

		case EqComp::FUNC_EXP_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_exp_pmul;
			break;

		case EqComp::FUNC_ASIN:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_asin;
			break;

		case EqComp::FUNC_ASIN_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_asin_pmul;
			break;

		case EqComp::FUNC_ACOS:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_acos;
			break;

		case EqComp::FUNC_ACOS_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_acos_pmul;
			break;

		case EqComp::FUNC_ATAN:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_atan;
			break;

		case EqComp::FUNC_ATAN_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_atan_pmul;
			break;

		case EqComp::FUNC_ASINH:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_asinh;
			break;

		case EqComp::FUNC_ASINH_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_asinh_pmul;
			break;

		case EqComp::FUNC_ACOSH:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_acosh;
			break;

		case EqComp::FUNC_ACOSH_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_acosh_pmul;
			break;

		case EqComp::FUNC_ATANH:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_atanh;
			break;

		case EqComp::FUNC_ATANH_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_atanh_pmul;
			break;

		case EqComp::FUNC_LN:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_ln;
			break;

		case EqComp::FUNC_LN_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_ln_pmul;
			break;

		case EqComp::FUNC_LOG:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_log;
			break;

		case EqComp::FUNC_LOG_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_log_pmul;
			break;

		case EqComp::FUNC_ABS:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_abs;
			break;

		case EqComp::FUNC_ABS_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_abs_pmul;
			break;

		case EqComp::FUNC_SGN:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_sgn;
			break;

		case EqComp::FUNC_SGN_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_sgn_pmul;
			break;

		case EqComp::FUNC_CEIL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_ceil;
			break;

		case EqComp::FUNC_CEIL_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_ceil_pmul;
			break;

		case EqComp::FUNC_FLOOR:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_floor;
			break;

		case EqComp::FUNC_FLOOR_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_floor_pmul;
			break;

		case EqComp::FUNC_ROUND:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_round;
			break;

		case EqComp::FUNC_ROUND_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_round_pmul;
			break;

		case EqComp::FUNC_STEP:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_step;
			break;

		case EqComp::FUNC_STEP_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_step_pmul;
			break;

		case EqComp::FUNC_SWAV:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_swav;
			break;

		case EqComp::FUNC_SWAV_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_swav_pmul;
			break;

		case EqComp::FUNC_TWAV:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_twav;
			break;

		case EqComp::FUNC_TWAV_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_twav_pmul;
			break;

		case EqComp::FUNC_POWER_EXPCONST:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_pow_expconst;
			break;

		case EqComp::FUNC_POWER_EXPCONST_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_pow_expconst_pmul;
			break;

		case EqComp::FUNC_POWER_BASECONST:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_pow_baseconst;
			break;

		case EqComp::FUNC_POWER_BASECONST_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_pow_baseconst_pmul;
			break;

		case EqComp::FUNC_CURIEWEISS:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_CurieWeiss;
			break;

		case EqComp::FUNC_CURIEWEISS1:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_CurieWeiss1;
			break;

		case EqComp::FUNC_CURIEWEISS2:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_CurieWeiss2;
			break;

		case EqComp::FUNC_LONGRELSUS:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_LongRelSus;
			break;

		case EqComp::FUNC_LONGRELSUS1:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_LongRelSus1;
			break;

		case EqComp::FUNC_LONGRELSUS2:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_LongRelSus2;
			break;

		case EqComp::FUNC_ALPHA1:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_Alpha1;
			break;

		case EqComp::FUNC_ALPHA2:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_Alpha2;
			break;
		}
	}
}

template <typename ... BVarType>
__host__ void ManagedFunctionCUDA<BVarType...>::Function_Unary(EqComp::FSPEC fspec, ManagedFunctionCUDA<BVarType...>& cuFunc1)
{
	set_gpu_value(param, (cuBReal)fspec.param);
	set_gpu_value(base_or_exponent, (cuBReal)fspec.base_or_exponent);

	set_ManagedFunctionCUDA_Unary<BVarType...> <<<1, 1>>> (*this, cuFunc1, fspec.type);
}

/////////////////////////////////////////////////////////////////////////////////////////////

//BINARY FUNCTIONS
//Setup this object as a binary operator on cuFunc1 and cuFunc2

template <typename ... BVarType>
__global__ void set_ManagedFunctionCUDA_Binary(ManagedFunctionCUDA<BVarType...>& cuFunc, ManagedFunctionCUDA<BVarType...>& cuFunc1, ManagedFunctionCUDA<BVarType...>& cuFunc2, int type)
{
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx == 0) {

		cuFunc.pFunc1 = &cuFunc1;
		cuFunc.pFunc2 = &cuFunc2;

		switch (type) {

		case EqComp::FUNC_ADD:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_add;
			break;

		case EqComp::FUNC_ADD_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_add_pmul;
			break;

		case EqComp::FUNC_SUB:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_sub;
			break;

		case EqComp::FUNC_SUB_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_sub_pmul;
			break;

		case EqComp::FUNC_MUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_mul;
			break;

		case EqComp::FUNC_MUL_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_mul_pmul;
			break;

		case EqComp::FUNC_DIV:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_div;
			break;

		case EqComp::FUNC_DIV_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_div_pmul;
			break;

		case EqComp::FUNC_POW:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_pow;
			break;

		case EqComp::FUNC_POW_PMUL:
			cuFunc.func = &ManagedFunctionCUDA<BVarType...>::F_pow_pmul;
			break;
		}
	}
}

template <typename ... BVarType>
__host__ void ManagedFunctionCUDA<BVarType...>::Operator_Binary(EqComp::FSPEC fspec, ManagedFunctionCUDA<BVarType...>& cuFunc1, ManagedFunctionCUDA<BVarType...>& cuFunc2)
{
	set_gpu_value(param, (cuBReal)fspec.param);

	set_ManagedFunctionCUDA_Binary<BVarType...> << <1, 1 >> > (*this, cuFunc1, cuFunc2, fspec.type);
}

