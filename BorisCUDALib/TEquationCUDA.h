#pragma once

#include "TEquationCUDA_Component.h"

template <typename ... BVarType>
class TEquationCUDA
{

private:

	TEquationCUDA_Component<BVarType...> eq_component_1;
	TEquationCUDA_Component<BVarType...> eq_component_2;
	TEquationCUDA_Component<BVarType...> eq_component_3;

private:

	void clear_Funcs(void)
	{
		eq_component_1.clear_Funcs();
		eq_component_2.clear_Funcs();
		eq_component_3.clear_Funcs();
	}

public:

	TEquationCUDA(void) {}

	~TEquationCUDA()
	{
		clear_Funcs();
	}

	bool make_scalar(const std::vector< std::vector<EqComp::FSPEC> >& fspec)
	{
		return eq_component_1.make(fspec);
	}

	bool make_vector(const std::vector<std::vector< std::vector<EqComp::FSPEC> >>& fspec)
	{
		if (fspec.size() != 3) return false;

		bool success = true;

		success &= eq_component_1.make(fspec[0]);
		success &= eq_component_2.make(fspec[1]);
		success &= eq_component_3.make(fspec[2]);

		return success;
	}

	bool is_set(void) { return eq_component_1.is_set(); }

	void clear(void) { clear_Funcs(); }

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
};















