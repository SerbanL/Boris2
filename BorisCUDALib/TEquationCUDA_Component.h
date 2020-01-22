#pragma once

#include "cuObject.h"

#include "TEquationCUDA_Function.h"

template <typename ... BVarType>
class TEquationCUDA_Component
{
public:

	std::vector< std::vector< cu_obj<ManagedFunctionCUDA<BVarType...>>* > > Funcs;

public:

	void clear_Funcs(void)
	{
		for (int idx = 0; idx < Funcs.size(); idx++) {

			for (int idx_tree = 0; idx_tree < Funcs[idx].size(); idx_tree++) {

				delete Funcs[idx][idx_tree];
			}

			Funcs[idx].clear();
		}

		Funcs.clear();
	}

	bool make(const std::vector< std::vector<EqComp::FSPEC> >& fspec)
	{
		clear_Funcs();

		for (int idx = 0; idx < fspec.size(); idx++) {

			//new branch
			Funcs.push_back(std::vector< cu_obj<ManagedFunctionCUDA<BVarType...>>* >());

			if (fspec[idx][0].is_binary_operator()) {

				int bin_idx1 = fspec[idx][0].bin_idx1;
				int bin_idx2 = fspec[idx][0].bin_idx2;

				if (bin_idx1 >= Funcs.size() || bin_idx2 >= Funcs.size()) return false;

				Funcs.back().push_back(new cu_obj<ManagedFunctionCUDA<BVarType...>>);
				(*(Funcs.back().back()))()->Operator_Binary(fspec[idx][0], *Funcs[bin_idx1].back(), *Funcs[bin_idx2].back());
			}
			else {

				Funcs.back().push_back(new cu_obj<ManagedFunctionCUDA<BVarType...>>);
				(*(Funcs.back().back()))()->Function_Basis(fspec[idx][0]);
			}

			for (int idx_tree = 1; idx_tree < fspec[idx].size(); idx_tree++) {

				Funcs.back().push_back(new cu_obj<ManagedFunctionCUDA<BVarType...>>);
				(*(Funcs.back().back()))()->Function_Unary(fspec[idx][idx_tree], *Funcs.back()[idx_tree - 1]);
			}
		}

		return true;
	}

	bool is_set(void) { return Funcs.size(); }

};