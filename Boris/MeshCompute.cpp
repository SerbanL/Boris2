#include "stdafx.h"
#include "Mesh.h"

#include "Exchange.h"
#include "DMExchange.h"
#include "iDMExchange.h"

//compute topological charge density spatial dependence and have it available to display in Cust_S
//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
void Mesh::Compute_TopoChargeDensity(void)
{
	if (M.linear_size()) {

		displayVEC_SCA.resize(h, meshRect);

#if COMPILECUDA == 1
		if (pMeshCUDA) {

			//compute topological charge density spatial variation in MeshCUDA::aux_vec_sca
			pMeshCUDA->Compute_TopoChargeDensity();

			//now transfer it to displayVEC_SCA ready to display
			pMeshCUDA->copy_aux_vec_sca(displayVEC_SCA);
			return;
		}
#endif

#pragma omp parallel for
		for (int idx = 0; idx < M.linear_size(); idx++) {

			if (M.is_not_empty(idx)) {

				double M_mag = M[idx].norm();

				DBL33 M_grad = M.grad_neu(idx);

				DBL3 dm_dx = M_grad.x / M_mag;
				DBL3 dm_dy = M_grad.y / M_mag;

				//divide by number of z cells (this is intended for 2D layers)
				displayVEC_SCA[idx] = (M[idx] / M_mag) * (dm_dx ^ dm_dy) * M.h.x * M.h.y / (4 * PI * M.n.z);
			}
			else displayVEC_SCA[idx] = 0.0;
		}
	}
}