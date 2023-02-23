#include "stdafx.h"
#include "MElastic.h"

#ifdef MODULE_COMPILATION_MELASTIC

#include "Mesh_Ferromagnetic.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "MElasticCUDA.h"
#endif

//------------------- STRAIN GENERATION without SOLVER

void MElastic::SetUniformStress(DBL3 Tsig_xyz)
{
	Tsig = Tsig_xyz;

	//From uniform stress set corresponding displacement vectors

	//Young's modulus of permalloy (Pa)
	double Ym = pMesh->Ym;
	//Poisson's ratio
	double Pr = pMesh->Pr;

	//Strain vector for uniform stress
	DBL3 strain = DBL33(DBL3(1.0, -Pr, -Pr), DBL3(-Pr, 1.0, -Pr), DBL3(-Pr, -Pr, 1.0)) * Tsig / Ym;

#pragma omp parallel for
	for (int j = 0; j < pMesh->n_m.j; j++) {
		for (int k = 0; k < pMesh->n_m.k; k++) {
			for (int i = 0; i < pMesh->n_m.i; i++) {

				//make origin the 0, 0, 0 displacement point
				//ignore any type of translational displacement
				//ignore any shear displacement
				DBL3 u = DBL3(
					(i + 0.5) * pMesh->h_m.x * strain.x,
					(j + 0.5) * pMesh->h_m.y * strain.y,
					(k + 0.5) * pMesh->h_m.z * strain.z);

				pMesh->u_disp[INT3(i, j, k)] = u;
			}
		}
	}

	Calculate_Strain();

	//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pModuleCUDA) dynamic_cast<MElasticCUDA*>(pModuleCUDA)->copy_VECs_to_GPU();
#endif
}

//Displacement : from this calculate strain tensor
BError MElastic::Load_Displacement_OVF2(std::string fileName)
{
	BError error(CLASS_STR(MElastic));

	OVF2 ovf2;
	VEC<DBL3> data;
	error = ovf2.Read_OVF2_VEC(fileName, data);
	if (error) return error;

	//copy displacement to current mesh size, then calculate corresponding strain
	pMesh->u_disp.copy_values(data);
	Calculate_Strain();

	//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pModuleCUDA) return dynamic_cast<MElasticCUDA*>(pModuleCUDA)->copy_VECs_to_GPU();
#endif

	return error;
}

//Tensor : displacement is not calculated, as we only need the strain tensor to obtain the effective fields at runtime
//For a strain tensor for cubic crystals, we only have 6 elements : diagonal and off-diagonal (symmetric).
//These are specified using 2 separate OVF2 files containing vector data:
//one for the xx, yy, zz elements (diagonal)
//the other for the yz, xz, xy elements (off-diagonal, in this order)
BError MElastic::Load_Strain_OVF2(std::string fileName_Diag, std::string fileName_ODiag)
{
	BError error(CLASS_STR(MElastic));

	OVF2 ovf2;
	VEC<DBL3> data;

	error = ovf2.Read_OVF2_VEC(fileName_Diag, data);
	if (error) return error;
	pMesh->strain_diag.copy_values(data);

	error = ovf2.Read_OVF2_VEC(fileName_ODiag, data);
	if (error) return error;
	pMesh->strain_odiag.copy_values(data);

	//-------------------------- CUDA mirroring

#if COMPILECUDA == 1
	if (pModuleCUDA) return dynamic_cast<MElasticCUDA*>(pModuleCUDA)->copy_VECs_to_GPU();
#endif

	return error;
}

//----------------------------------------------- Auxiliary

//Run-time auxiliary to set strain directly from user supplied text formulas
void MElastic::Set_Strain_From_Formula(void)
{
	double time = pSMesh->GetStageTime();

	for (int k = 0; k < pMesh->n_m.k; k++) {
#pragma omp parallel for
		for (int j = 0; j < pMesh->n_m.j; j++) {
			for (int i = 0; i < pMesh->n_m.i; i++) {

				int idx = i + j * pMesh->n_m.i + k * pMesh->n_m.i * pMesh->n_m.j;

				if (pMesh->u_disp.is_empty(idx)) continue;

				DBL3 pos = DBL3(i + 0.5, j + 0.5, k + 0.5) & pMesh->h_m;

				if (Sd_equation.is_set_vector()) pMesh->strain_diag[idx] = Sd_equation.evaluate_vector(pos.x, pos.y, pos.z, time);
				else pMesh->strain_diag[idx] = DBL3();

				if (Sod_equation.is_set_vector()) pMesh->strain_odiag[idx] = Sod_equation.evaluate_vector(pos.x, pos.y, pos.z, time);
				else pMesh->strain_odiag[idx] = DBL3();
			}
		}
	}
}

//----------------------------------------------- Computational Helpers

//compute strain form mechanical displacement
void MElastic::Calculate_Strain(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n_m.dim(); idx++) {

		if (pMesh->u_disp.is_not_empty(idx)) {

			//get all 9 first-order differentials of u
			DBL33 grad_u = pMesh->u_disp.grad_sided(idx);

			//diagonal components
			pMesh->strain_diag[idx] = DBL3(grad_u.x.x, grad_u.y.y, grad_u.z.z);

			//off-diagonal components (yz, xz, xy)
			pMesh->strain_odiag[idx] = 0.5 * DBL3(grad_u.y.z + grad_u.z.y, grad_u.x.z + grad_u.z.x, grad_u.x.y + grad_u.y.x);
		}
		else {

			pMesh->strain_diag[idx] = DBL3();
			pMesh->strain_odiag[idx] = DBL3();
		}
	}
}

#endif
