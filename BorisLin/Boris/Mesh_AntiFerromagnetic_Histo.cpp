#include "stdafx.h"
#include "Mesh_AntiFerromagnetic.h"
#include "SuperMesh.h"

#ifdef MESH_COMPILATION_ANTIFERROMAGNETIC

//compute magnitude histogram data
//extract histogram between magnitudes min and max with given number of bins. if min max not given (set them to zero) then determine them first. 
//output probabilities in histogram_p, corresponding to values set in histogram_x min, min + bin, ..., max, where bin = (max - min) / (num_bins - 1)
//if macrocell_dims is not INT3(1) then first average in macrocells containing given number of individual mesh cells, then obtain histogram
bool AFMesh::Get_Histogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, INT3 macrocell_dims)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		return pMeshCUDA->M()->get_mag_histogram(histogram_x, histogram_p, num_bins, min, max, M.get_nonempty_cells(), macrocell_dims);
	}
#endif

	return M.get_mag_histogram(histogram_x, histogram_p, num_bins, min, max, M.get_nonempty_cells(), macrocell_dims);
}

//angular deviation histogram computed from ndir unit vector direction. If ndir not given (DBL3()), then angular deviation computed from average magnetization direction
bool AFMesh::Get_AngHistogram(std::vector<double>& histogram_x, std::vector<double>& histogram_p, int num_bins, double& min, double& max, INT3 macrocell_dims, DBL3 ndir)
{
#if COMPILECUDA == 1
	if (pMeshCUDA) {

		return pMeshCUDA->M()->get_ang_histogram(histogram_x, histogram_p, num_bins, min, max, M.get_nonempty_cells(), macrocell_dims, ndir);
	}
#endif

	return M.get_ang_histogram(histogram_x, histogram_p, num_bins, min, max, M.get_nonempty_cells(), macrocell_dims, ndir);
}

#endif