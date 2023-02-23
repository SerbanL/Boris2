#include "stdafx.h"
#include "MeshBase.h"

//----------------------------------- DISPLAY-ASSOCIATED GET/SET METHODS AUXILIARY

//average into profile_storage_dbl / profile_storage_dbl3
void MeshBase::average_mesh_profile(std::vector<double>& line_profile)
{
	if (profile_storage_dbl.size() != line_profile.size()) {

		profile_storage_dbl.resize(line_profile.size());
		num_profile_averages = 0;
	}

	if (num_profile_averages) {

#pragma omp parallel for
		for (int idx = 0; idx < line_profile.size(); idx++) {

			profile_storage_dbl[idx] = (profile_storage_dbl[idx] * num_profile_averages + line_profile[idx]) / (num_profile_averages + 1);
		}
	}
	else {

#pragma omp parallel for
		for (int idx = 0; idx < line_profile.size(); idx++) {

			profile_storage_dbl[idx] = line_profile[idx];
		}
	}

	//if num_profile_averages == 0, it's possible profile_storage_sca is not initialized and currently storing a nan

	num_profile_averages++;
}

void MeshBase::average_mesh_profile(std::vector<DBL3>& line_profile)
{
	if (profile_storage_dbl3.size() != line_profile.size()) {

		profile_storage_dbl3.resize(line_profile.size());
		num_profile_averages = 0;
	}

	if (num_profile_averages) {

#pragma omp parallel for
		for (int idx = 0; idx < line_profile.size(); idx++) {

			profile_storage_dbl3[idx] = (profile_storage_dbl3[idx] * num_profile_averages + line_profile[idx]) / (num_profile_averages + 1);
		}
	}
	else {

#pragma omp parallel for
		for (int idx = 0; idx < line_profile.size(); idx++) {

			profile_storage_dbl3[idx] = line_profile[idx];
		}
	}

	num_profile_averages++;
}