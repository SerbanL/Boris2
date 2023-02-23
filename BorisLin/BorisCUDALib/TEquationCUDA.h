#pragma once

#include "TEquationCUDA_Component.h"

template <typename ... BVarType>
class TEquationCUDA
{

private:

	TEquationCUDA_Component<BVarType...> eq_component_1;
	TEquationCUDA_Component<BVarType...> eq_component_2;
	TEquationCUDA_Component<BVarType...> eq_component_3;

	//Special Functions

	cu_obj<ManagedFuncs_Special_CUDA>* pCurieWeiss = nullptr;
	cu_obj<ManagedFuncs_Special_CUDA>* pLongRelSus = nullptr;
	cu_obj<ManagedFuncs_Special_CUDA>* pCurieWeiss1 = nullptr;
	cu_obj<ManagedFuncs_Special_CUDA>* pCurieWeiss2 = nullptr;
	cu_obj<ManagedFuncs_Special_CUDA>* pLongRelSus1 = nullptr;
	cu_obj<ManagedFuncs_Special_CUDA>* pLongRelSus2 = nullptr;
	cu_obj<ManagedFuncs_Special_CUDA>* pAlpha1 = nullptr;
	cu_obj<ManagedFuncs_Special_CUDA>* pAlpha2 = nullptr;

private:

	void clear_Funcs(void)
	{
		eq_component_1.clear_Funcs();
		eq_component_2.clear_Funcs();
		eq_component_3.clear_Funcs();
	}

	void Set_StoredSpecialFunctions(void)
	{
		auto set = [&](EqComp::FUNC_ type, cu_obj<ManagedFuncs_Special_CUDA>* pSpecialFunction) -> void {

			if (pSpecialFunction) {

				eq_component_1.Set_SpecialFunction(type, pSpecialFunction);
				eq_component_2.Set_SpecialFunction(type, pSpecialFunction);
				eq_component_3.Set_SpecialFunction(type, pSpecialFunction);
			}
		};

		set(EqComp::FUNC_CURIEWEISS, pCurieWeiss);
		set(EqComp::FUNC_LONGRELSUS, pLongRelSus);
		set(EqComp::FUNC_CURIEWEISS1, pCurieWeiss1);
		set(EqComp::FUNC_CURIEWEISS2, pCurieWeiss2);
		set(EqComp::FUNC_LONGRELSUS1, pLongRelSus1);
		set(EqComp::FUNC_LONGRELSUS2, pLongRelSus2);
		set(EqComp::FUNC_ALPHA1, pAlpha1);
		set(EqComp::FUNC_ALPHA2, pAlpha2);
	}

public:

	/////////////////////////////////////////////////////////
	//
	// CONSTRUCTOR

	TEquationCUDA(void) {}

	~TEquationCUDA()
	{
		clear_Funcs();
	}

	/////////////////////////////////////////////////////////
	//
	// MAKE EQUATION

	bool make_scalar(const std::vector< std::vector<EqComp::FSPEC> >& fspec)
	{
		bool success = eq_component_1.make(fspec);

		if (success) Set_StoredSpecialFunctions();

		return success;
	}

	bool make_dual(const std::vector<std::vector< std::vector<EqComp::FSPEC> >>& fspec)
	{
		if (fspec.size() < 2) return false;

		bool success = true;

		success &= eq_component_1.make(fspec[0]);
		success &= eq_component_2.make(fspec[1]);

		if (success) Set_StoredSpecialFunctions();

		return success;
	}

	bool make_vector(const std::vector<std::vector< std::vector<EqComp::FSPEC> >>& fspec)
	{
		if (fspec.size() < 3) return false;

		bool success = true;

		success &= eq_component_1.make(fspec[0]);
		success &= eq_component_2.make(fspec[1]);
		success &= eq_component_3.make(fspec[2]);

		if (success) Set_StoredSpecialFunctions();

		return success;
	}

	/////////////////////////////////////////////////////////
	//
	// DESTROY EQUATION

	void clear(void) { clear_Funcs(); }

	/////////////////////////////////////////////////////////
	//
	// EQUATION TYPE CHECKING

	bool is_set(void) const { return eq_component_1.is_set(); }
	bool is_set_dual(void) const { return eq_component_1.is_set() && eq_component_2.is_set(); }
	bool is_set_vector(void) const { return eq_component_1.is_set() && eq_component_2.is_set() && eq_component_3.is_set(); }

	/////////////////////////////////////////////////////////
	//
	// LOGICAL EQUATION GETTERS

	__host__ operator ManagedFunctionCUDA<BVarType...>&() { return *eq_component_1.Funcs.back().back(); }

	__host__ ManagedFunctionCUDA<BVarType...>& get_x(void) { return *eq_component_1.Funcs.back().back(); }
	__host__ ManagedFunctionCUDA<BVarType...>& get_y(void) { return *eq_component_2.Funcs.back().back(); }
	__host__ ManagedFunctionCUDA<BVarType...>& get_z(void) { return *eq_component_3.Funcs.back().back(); }

	__host__ cu_obj<ManagedFunctionCUDA<BVarType...>>* get_pcu_obj_x(void)
	{
		if (eq_component_1.Funcs.size()) return eq_component_1.Funcs.back().back();
		else return nullptr;
	}
	__host__ cu_obj<ManagedFunctionCUDA<BVarType...>>* get_pcu_obj_y(void)
	{
		if (eq_component_2.Funcs.size()) return eq_component_2.Funcs.back().back();
		else return nullptr;
	}

	__host__ cu_obj<ManagedFunctionCUDA<BVarType...>>* get_pcu_obj_z(void)
	{
		if (eq_component_3.Funcs.size()) return eq_component_3.Funcs.back().back();
		else return nullptr;
	}

	/////////////////////////////////////////////////////////
	//
	// SET SPECIAL FUNCTIONS

	void Set_SpecialFunction(EqComp::FUNC_ type, cu_obj<ManagedFuncs_Special_CUDA>* pSpecialFunc)
	{
		switch (type) {

		case EqComp::FUNC_CURIEWEISS:
			pCurieWeiss = pSpecialFunc;
			break;

		case EqComp::FUNC_CURIEWEISS1:
			pCurieWeiss1 = pSpecialFunc;
			break;

		case EqComp::FUNC_CURIEWEISS2:
			pCurieWeiss2 = pSpecialFunc;
			break;

		case EqComp::FUNC_LONGRELSUS:
			pLongRelSus = pSpecialFunc;
			break;

		case EqComp::FUNC_LONGRELSUS1:
			pLongRelSus1 = pSpecialFunc;
			break;

		case EqComp::FUNC_LONGRELSUS2:
			pLongRelSus2 = pSpecialFunc;
			break;

		case EqComp::FUNC_ALPHA1:
			pAlpha1 = pSpecialFunc;
			break;

		case EqComp::FUNC_ALPHA2:
			pAlpha2 = pSpecialFunc;
			break;
		}

		eq_component_1.Set_SpecialFunction(type, pSpecialFunc);
		eq_component_2.Set_SpecialFunction(type, pSpecialFunc);
		eq_component_3.Set_SpecialFunction(type, pSpecialFunc);
	}
};