#include "stdafx.h"
#include "Mesh_Ferromagnetic.h"
#include "SuperMesh.h"

#ifdef MESH_COMPILATION_FERROMAGNETIC

//compute magnitude histogram data
//extract histogram between magnitudes min and max with given number of bins. if min max not given (set them to zero) then determine them first. 
//output probabilities in histogram_p, corresponding to values set in histogram_x min, min + bin, ..., max, where bin = (max - min) / (num_bins - 1)
//if macrocell_dims is not INT3(1) then first average in macrocells containing given number of individual mesh cells, then obtain histogram
bool FMesh::Get_Histogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, INT3 macrocell_dims)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		return pMeshCUDA->M()->get_mag_histogram(histogram_x, histogram_p, num_bins, min, max, M.get_nonempty_cells(), macrocell_dims);
	}
#endif

	return M.get_mag_histogram(histogram_x, histogram_p, num_bins, min, max, M.get_nonempty_cells(), macrocell_dims);
}

//angular deviation histogram computed from ndir unit vector direction. If ndir not given (DBL3()), then angular deviation computed from average magnetization direction
bool FMesh::Get_AngHistogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, INT3 macrocell_dims, DBL3 ndir)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		return pMeshCUDA->M()->get_ang_histogram(histogram_x, histogram_p, num_bins, min, max, M.get_nonempty_cells(), macrocell_dims, ndir);
	}
#endif

	return M.get_ang_histogram(histogram_x, histogram_p, num_bins, min, max, M.get_nonempty_cells(), macrocell_dims, ndir);
}

//As for Get_Histogram, but use thermal averaging in each macrocell
bool FMesh::Get_ThAvHistogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, INT3 macrocell_dims)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		return reinterpret_cast<FMeshCUDA*>(pMeshCUDA)->Get_ThAvHistogram(histogram_x, histogram_p, num_bins, min, max, macrocell_dims);
	}
#endif

	//First do thermal cell-wise pre-averaging

	//allocate required memory for auxVEC
	SZ3 num_av_cells = round((DBL3)n / macrocell_dims);
	auxVEC.resize(num_av_cells);

	//cell-wise pre-averaging
	for (int k_seg = 0; k_seg < num_av_cells.k; k_seg++) {
#pragma omp parallel for
		for (int j_seg = 0; j_seg < num_av_cells.j; j_seg++) {
			for (int i_seg = 0; i_seg < num_av_cells.i; i_seg++) {

				DBL3 th_average = DBL3();
				double Z_cell = 0.0;

				for (int k = 0; k < macrocell_dims.k; k++) {
					for (int j = 0; j < macrocell_dims.j; j++) {
						for (int i = 0; i < macrocell_dims.i; i++) {

							INT3 ijk = INT3(i_seg * macrocell_dims.i + i, j_seg * macrocell_dims.j + j, k_seg * macrocell_dims.k + k);
							int idx = ijk.i + ijk.j * M.n.x + ijk.k * M.n.x*M.n.y;

							if (idx < M.linear_size() && M.is_not_empty(idx)) {

								double Ei = 0.0;
								for (int mod_idx = 0; mod_idx < pMod.size(); mod_idx++) {

									Ei += pMod[mod_idx]->Get_Energy(idx);
								}

								double Temperature;
								if (Temp.linear_size()) Temperature = Temp[M.cellidx_to_position(idx)];
								else Temperature = base_temperature;

								//Longitudinal term contribution
								double Ms_val = Ms;
								double susrel_val = susrel;
								update_parameters_mcoarse(idx, Ms, Ms_val, susrel, susrel_val);

								double Ms0 = Ms.get0();
								double m = M[idx].norm() / Ms0;

								if (Temperature <= T_Curie) {

									double me = Ms_val / Ms0;
									double diff = m * m - me * me;

									Ei += h.dim() * (Ms0 / (8 * susrel_val * me*me)) * diff * diff;
								}
								else {

									double r = 3 * T_Curie / (10 * (Temperature - T_Curie));
									double m_sq = m * m;
									Ei += h.dim() * (Ms0 / (2 * susrel_val)) * m_sq * (1 + r * m_sq);
								}
								
								double w = m*m*exp(-Ei / (BOLTZMANN * Temperature));
								Z_cell += w;

								th_average += w * M[idx];
							}
						}
					}
				}

				if (Z_cell) auxVEC[INT3(i_seg, j_seg, k_seg)] = th_average / Z_cell;
				else auxVEC[INT3(i_seg, j_seg, k_seg)] = DBL3();
			}
		}
	}

	//get histogram from auxVEC
	return auxVEC.get_mag_histogram(histogram_x, histogram_p, num_bins, min, max);
}

//As for Get_AngHistogram, but use thermal averaging in each macrocell
bool FMesh::Get_ThAvAngHistogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, INT3 macrocell_dims, DBL3 ndir)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		return reinterpret_cast<FMeshCUDA*>(pMeshCUDA)->Get_ThAvAngHistogram(histogram_x, histogram_p, num_bins, min, max, macrocell_dims, ndir);
	}
#endif

	//First do thermal cell-wise pre-averaging

	//allocate required memory for auxVEC
	SZ3 num_av_cells = round((DBL3)n / macrocell_dims);
	auxVEC.resize(num_av_cells);

	//cell-wise pre-averaging
	for (int k_seg = 0; k_seg < num_av_cells.k; k_seg++) {
#pragma omp parallel for
		for (int j_seg = 0; j_seg < num_av_cells.j; j_seg++) {
			for (int i_seg = 0; i_seg < num_av_cells.i; i_seg++) {

				DBL3 th_average = DBL3();
				double Z_cell = 0.0;

				for (int k = 0; k < macrocell_dims.k; k++) {
					for (int j = 0; j < macrocell_dims.j; j++) {
						for (int i = 0; i < macrocell_dims.i; i++) {

							INT3 ijk = INT3(i_seg * macrocell_dims.i + i, j_seg * macrocell_dims.j + j, k_seg * macrocell_dims.k + k);
							int idx = ijk.i + ijk.j * M.n.x + ijk.k * M.n.x*M.n.y;

							if (idx < M.linear_size() && M.is_not_empty(idx)) {

								double Ei = 0.0;
								for (int mod_idx = 0; mod_idx < pMod.size(); mod_idx++) {

									Ei += pMod[mod_idx]->Get_Energy(idx);
								}

								double Temperature;
								if (Temp.linear_size()) Temperature = Temp[M.cellidx_to_position(idx)];
								else Temperature = base_temperature;

								//Longitudinal term contribution
								double Ms_val = Ms;
								double susrel_val = susrel;
								update_parameters_mcoarse(idx, Ms, Ms_val, susrel, susrel_val);

								double Ms0 = Ms.get0();
								double m = M[idx].norm() / Ms0;

								if (Temperature <= T_Curie) {

									double me = Ms_val / Ms0;
									double diff = m * m - me * me;

									Ei += h.dim() * (Ms0 / (8 * susrel_val * me*me)) * diff * diff;
								}
								else {

									double r = 3 * T_Curie / (10 * (Temperature - T_Curie));
									double m_sq = m * m;
									Ei += h.dim() * (Ms0 / (2 * susrel_val)) * m_sq * (1 + r * m_sq);
								}

								double w = m*m*exp(-Ei / (BOLTZMANN * Temperature));
								Z_cell += w;

								th_average += w * M[idx];
							}
						}
					}
				}

				if (Z_cell) auxVEC[INT3(i_seg, j_seg, k_seg)] = th_average / Z_cell;
				else auxVEC[INT3(i_seg, j_seg, k_seg)] = DBL3();
			}
		}
	}

	if (ndir.IsNull()) ndir = GetThermodynamicAverageMagnetization(Rect()).normalized();

	//get histogram from auxVEC
	return auxVEC.get_ang_histogram(histogram_x, histogram_p, num_bins, min, max, num_av_cells.dim(), INT3(1), ndir);
}

DBL3 FMesh::GetThermodynamicAverageMagnetization(Rect rectangle)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		return reinterpret_cast<FMeshCUDA*>(pMeshCUDA)->GetThermodynamicAverageMagnetization(rectangle);
	}
#endif

	double Z = 0.0;
	double Mthav_x = 0.0, Mthav_y = 0.0, Mthav_z = 0.0;

	if (rectangle.IsNull()) rectangle = meshRect;

#pragma omp parallel for reduction(+:Z, Mthav_x, Mthav_y, Mthav_z)
	for (int idx = 0; idx < M.linear_size(); idx++) {

		if (M.is_not_empty(idx) && rectangle.contains(M.cellidx_to_position(idx))) {

			double Ei = 0.0;
			for (int mod_idx = 0; mod_idx < pMod.size(); mod_idx++) {

				Ei += pMod[mod_idx]->Get_Energy(idx);
			}

			double Temperature;
			if (Temp.linear_size()) Temperature = Temp[M.cellidx_to_position(idx)];
			else Temperature = base_temperature;

			//Longitudinal term contribution
			double Ms_val = Ms;
			double susrel_val = susrel;
			update_parameters_mcoarse(idx, Ms, Ms_val, susrel, susrel_val);

			double Ms0 = Ms.get0();
			double m = M[idx].norm() / Ms0;

			if (Temperature <= T_Curie) {

				double me = Ms_val / Ms0;
				double diff = m * m - me * me;

				Ei += h.dim() * (Ms0 / (8 * susrel_val * me*me)) * diff * diff;
			}
			else {

				double r = 3 * T_Curie / (10 * (Temperature - T_Curie));
				double m_sq = m * m;
				Ei += h.dim() * (Ms0 / (2 * susrel_val)) * m_sq * (1 + r * m_sq);
			}

			double w = m*m*exp(-Ei / (BOLTZMANN * Temperature));
			Z += w;

			Mthav_x += w * M[idx].x;
			Mthav_y += w * M[idx].y;
			Mthav_z += w * M[idx].z;
		}
	}

	return DBL3(Mthav_x, Mthav_y, Mthav_z) / Z;
}

#endif