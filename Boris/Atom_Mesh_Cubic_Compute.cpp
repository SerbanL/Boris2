#include "stdafx.h"
#include "Atom_Mesh_Cubic.h"

#ifdef MESH_COMPILATION_ATOM_CUBIC

//compute topological charge density spatial dependence and have it available to display in Cust_S
//Use formula Qdensity = m.(dm/dx x dm/dy) / 4PI
void Atom_Mesh_Cubic::Compute_TopoChargeDensity(void)
{
	if (M1.linear_size()) {

		displayVEC_SCA.resize(h, meshRect);

#if COMPILECUDA == 1
		if (paMeshCUDA) {

			//compute topological charge density spatial variation in MeshCUDA::aux_vec_sca
			paMeshCUDA->Compute_TopoChargeDensity();

			//now transfer it to displayVEC_SCA ready to display
			paMeshCUDA->copy_aux_vec_sca(displayVEC_SCA);
			return;
		}
#endif

#pragma omp parallel for
		for (int idx = 0; idx < M1.linear_size(); idx++) {

			if (M1.is_not_empty(idx)) {

				double Mnorm = M1[idx].norm();

				DBL33 M_grad = M1.grad_neu(idx);

				DBL3 dm_dx = M_grad.x / Mnorm;
				DBL3 dm_dy = M_grad.y / Mnorm;

				//divide by number of z cells (this is intended for 2D layers)
				displayVEC_SCA[idx] = (M1[idx] / Mnorm) * (dm_dx ^ dm_dy) * M1.h.x * M1.h.y / (4 * PI * M1.n.z);
			}
			else displayVEC_SCA[idx] = 0.0;
		}
	}
}

//return phase transition temperature (K) based on formula Tc = J*e*z/3kB
double Atom_Mesh_Cubic::Show_Transition_Temperature(void)
{
	return J * spinwave_factor * coordination_number / (3 * BOLTZMANN);
}

//return saturation magnetization (A/m) based on formula Ms = mu_s*n/a^3
double Atom_Mesh_Cubic::Show_Ms(void)
{
	return mu_s * (MUB / h.dim()) * atoms_per_cell;
}

//return exchange stiffness (J/m) based on formula A = J*n/2a
double Atom_Mesh_Cubic::Show_A(void)
{
	//formula applicable for a cubic cell
	double h_av = (h.x + h.y + h.z) / 3;

	return J * atoms_per_cell / (2 * h_av);
}

//return uniaxial anisotropy constant (J/m^3) based on formula K = k*n/a^3
double Atom_Mesh_Cubic::Show_Ku(void)
{
	return K * atoms_per_cell / h.dim();
}

#endif