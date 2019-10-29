#include "stdafx.h"
#include "Roughness.h"

#ifdef MODULE_ROUGHNESS

#include "Mesh.h"
#include "SuperMesh.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

Roughness::Roughness(Mesh *pMesh_) :
	Modules(),
	Convolution<RoughnessKernel<Roughness>>(pMesh_->GetMeshSize(), pMesh_->GetMeshCellsize(), this),
	ProgramStateNames(this, {VINFO(refine), VINFO(Mshape_fine)}, {})
{
	pMesh = pMesh_;

	Uninitialize();

	error_on_create = Convolution_Error_on_Create();

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);
}

void Roughness::RepairObjectState(void) {}

//set shape of M from Mshape_fine : any coarse cells which have at least one fine cell non-empty must be set non-empty
void Roughness::Set_Coarse_M_Inclusion(void)
{
	pMesh->M.resize(pMesh->h, pMesh->meshRect, Mshape_fine);

	if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

		pMesh->M2.resize(pMesh->h, pMesh->meshRect, Mshape_fine);
	}

	//it's possible mesh was reset and previously empty cells are now not empty : must set a value in these as they will now have zero value
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->M.linear_size(); idx++) {

		if (pMesh->M.is_not_empty(idx) && pMesh->M[idx] == DBL3(0.0)) {

			pMesh->M[idx] = DBL3(-pMesh->Ms.get0(), 0, 0);

			if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) pMesh->M2[idx] = pMesh->M[idx] * -1.0;
		}
	}

#if COMPILECUDA == 1
	if (pMesh->pMeshCUDA) pMesh->pMeshCUDA->M()->set_from_cpuvec(pMesh->M);
#endif
}

//set values in Mshape_fine so it display nicely with color coding along z axis, indicating roughness
void Roughness::Calculate_Mshape_fine(void)
{
	SZ3 n = Mshape_fine.n;

#pragma omp parallel for
	for (int j = 0; j < n.y; j++) {
		for (int i = 0; i < n.x; i++) {

			if (n.z >= 2) {

				double value = 1.0;
				int cells_thick = 0;
				int half = n.z / 2;

				//go through top half : count number of non-empty cells in top half along z direction
				for (int k = half; k < n.z; k++) {

					int idx = i + j * n.x + k * n.x*n.y;

					if (Mshape_fine.is_not_empty(idx)) cells_thick++;
					else break;
				}

				//value ranges from 0 to 1
				value = (double)cells_thick / (n.z - half);

				//set value in top half
				for (int k = half; k < n.z; k++) {

					int idx = i + j * n.x + k * n.x*n.y;

					if (Mshape_fine.is_not_empty(idx)) Mshape_fine[idx] = value;
				}

				value = 1.0;
				cells_thick = 0;

				//go through bottom half : count number of non-empty cells in top half along z direction
				for (int k = half - 1; k >= 0; k--) {

					int idx = i + j * n.x + k * n.x*n.y;

					if (Mshape_fine.is_not_empty(idx)) cells_thick++;
					else break;
				}

				//value ranges from 0 to 1
				value = (double)cells_thick / half;

				//set value in top half
				for (int k = half - 1; k >= 0; k--) {

					int idx = i + j * n.x + k * n.x*n.y;

					if (Mshape_fine.is_not_empty(idx)) Mshape_fine[idx] = value;
				}
			}
			else Mshape_fine[i + j * n.x] = 1.0;
		}
	}
}

//copy roughness from another Roughness module (e.g. from a different mesh - the values, including shape, and discretisation are copied)
BError Roughness::copy_roughness(Roughness* pRough_copy_this)
{
	BError error(__FUNCTION__);

	Uninitialize();

	//refinement must be the same
	refine = pRough_copy_this->get_refine();

	//discretisation may have changed so must update (also takes care of CUDA version)
	UpdateConfiguration(UPDATECONFIG_ROUGHNESS_CHANGE);

	//copy values
	Mshape_fine.copy_values(pRough_copy_this->Mshape_fine);

	//set shape of M from Mshape_fine (also takes care of CUDA version)
	Set_Coarse_M_Inclusion();

	return error;
}

//clear roughness: set fine shape to coarse shape
void Roughness::clear_roughness(void)
{
	Uninitialize();
	
	//rest shape back to full rectangle
	Mshape_fine.setrect(Mshape_fine.rect, 1);

	//now apply the coarse shape
	Mshape_fine.resize(pMesh->h / refine, pMesh->meshRect, pMesh->M);

	Calculate_Mshape_fine();
}

//compute Fmul_rough (diagonal) or Fomul_rough (off-diagonal) after the Roughness Kernel has been computed
BError Roughness::Initialize_Fmul(VEC<DBL3>& Fout, bool off_diagonal)
{
	BError error(__FUNCTION__);

	off_diagonal_convolution = off_diagonal;

	//Obtain VECs F1_fine, F2_fine which we'll use at the end to get Fmul_rough
	
	VEC<DBL3> F1_fine, F2_fine, Fconv;

	if (!F1_fine.resize(pMesh->h / refine, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
	if (!F2_fine.resize(pMesh->h / refine, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
	if (!Fconv.resize(pMesh->h / refine, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);

	//number of non-empty cells in "roughened" body : Nvr
	int Nvr = Mshape_fine.get_nonempty_cells();

	//there must be some non-empty cells in Mshape_fine
	if (Nvr == 0) return error(BERROR_INCORRECTCONFIG);

	//number of non-empty cells in "smooth" body : Nv, counted on the same fine mesh as Nvr, hence multiplying by refine.dim()
	int Nv = pMesh->M.get_nonempty_cells() * refine.dim();
	
	//Setup convolution input function for F1_fine
#pragma omp parallel for
	for (int idx = 0; idx < Mshape_fine.linear_size(); idx++) {

		if (Mshape_fine.is_not_empty(idx)) {

			Fconv[idx] = DBL3((double)Nv / (double)Nvr - 1);
		}
		else {
			
			if (pMesh->M.is_empty(Mshape_fine.cellidx_to_position(idx))) {

				Fconv[idx] = DBL3(0.0);
			}
			else Fconv[idx] = DBL3 (-1.0);
		}
	}

	//Now convolute to get F1_fine
	Convolute(Fconv, F1_fine, true);

	//Setup convolution input function for F2_fine
#pragma omp parallel for
	for (int idx = 0; idx < Mshape_fine.linear_size(); idx++) {

		if (pMesh->M.is_empty(Mshape_fine.cellidx_to_position(idx))) {

			Fconv[idx] = DBL3(0.0);
		}
		else Fconv[idx] = DBL3(-1.0);
	}
	
	//Now convolute to get F2_fine
	Convolute(Fconv, F2_fine, true);
	
	//Combine F1_fine and F2_fine into a fine version of Fmul_rough -> set result in F1_fine
#pragma omp parallel for
	for (int idx = 0; idx < Mshape_fine.linear_size(); idx++) {

		//in V - Vr the fine version of Fmul_rough takes on F2 value
		//In Vr it takes on F1 value

		if (Mshape_fine.is_empty(idx)) {

			//idx is outside of Vr, so take on F2 value
			F1_fine[idx] = F2_fine[idx];
		}
	}

	//now average F1_fine up to Fmul_rough (or Fomul_rough) - Fout - for final result
#pragma omp parallel for
	for (int idx = 0; idx < Fout.linear_size(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			DBL3 position = Fout.cellidx_to_position(idx);

			//average from fine mesh to coarse mesh
			Fout[idx] = F1_fine.average(Rect(position - pMesh->h/2, position + pMesh->h/2));
		}
		else Fout[idx] = DBL3(0);
	}

	return error;
}

BError Roughness::Initialize(void)
{
	BError error(CLASS_STR(Roughness));

	if (!initialized) {

		error = Calculate_Roughness_Kernels();

		//now that we have the fine demag kernel need to compute the Fmul_rough function for runtime UpdateField
		if (!error) error = Initialize_Fmul(Fmul_rough, false);
		if (!error) error = Initialize_Fmul(Fomul_rough, true);

		if (!error) initialized = true;
	}

	return error;
}

BError Roughness::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(Roughness));

	if (refine.x == 0 || refine.y == 0 || refine.z == 0) return error(BERROR_INCORRECTCONFIG);

	//Do we need to take into account PBC when calculating the fine roughness kernel?
	//If PBCs are being used we have to set same settings here for the roughness formulas to apply.
	//If PBCs are set, these will be found either in the individual mesh Demag module, or else in the SDemag module.
	//Must be careful though if using many PBC images and fine roughness discretisation : this will result in long kernel computation times.
	INT3 pbc_images;

	if (pMesh->IsModuleSet(MOD_DEMAG)) {

		pbc_images = pMesh->CallModuleMethod(&Demag::Get_PBC);
	}
	else if (pMesh->pSMesh->IsSuperMeshModuleSet(MODS_SDEMAG)) {

		pbc_images = pMesh->pSMesh->CallModuleMethod(&SDemag::Get_PBC);
	}
	else pbc_images = INT3();

	//only need to uninitialize if n or h have changed
	if (!CheckDimensions(pMesh->n & refine, pMesh->h / refine, pbc_images)) {

		Uninitialize();
		error = SetDimensions(pMesh->n & refine, pMesh->h / refine, true, pbc_images);
	}

	//resize Fmul_rough and Fomul_rough (output multiplicative functions) for the coarse F mesh
	if (!Fmul_rough.resize(pMesh->h, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
	if (!Fomul_rough.resize(pMesh->h, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);

	//resize Mshape_fine for the fine mesh
	if (Mshape_fine.linear_size()) {

		if (!Mshape_fine.resize(pMesh->h / refine, pMesh->meshRect)) return error(BERROR_OUTOFMEMORY_CRIT);
	}
	else {

		if (!Mshape_fine.assign(pMesh->h / refine, pMesh->meshRect, 1.0)) return error(BERROR_OUTOFMEMORY_CRIT);

		//When starting module, M may already have a shape and we want to transfer it to Mshape_fine
		Mshape_fine.resize(pMesh->h / refine, pMesh->meshRect, pMesh->M);
	}

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError Roughness::MakeCUDAModule(void)
{
	BError error(CLASS_STR(Roughness));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new RoughnessCUDA(pMesh->pMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double Roughness::UpdateField(void)
{
	double energy = 0;

	if (pMesh->GetMeshType() == MESH_FERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				DBL3 Hrough = DBL33(
					DBL3(Fmul_rough[idx].x, Fomul_rough[idx].x, Fomul_rough[idx].y),
					DBL3(Fomul_rough[idx].x, Fmul_rough[idx].y, Fomul_rough[idx].z),
					DBL3(Fomul_rough[idx].y, Fomul_rough[idx].z, Fmul_rough[idx].z)) * pMesh->M[idx];

				pMesh->Heff[idx] += Hrough;

				energy += Hrough * pMesh->M[idx];
			}
		}
	}

	else if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy)
		for (int idx = 0; idx < pMesh->n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				DBL3 Hrough = DBL33(
					DBL3(Fmul_rough[idx].x, Fomul_rough[idx].x, Fomul_rough[idx].y),
					DBL3(Fomul_rough[idx].x, Fmul_rough[idx].y, Fomul_rough[idx].z),
					DBL3(Fomul_rough[idx].y, Fomul_rough[idx].z, Fmul_rough[idx].z)) * (pMesh->M[idx] + pMesh->M2[idx]) / 2;

				pMesh->Heff[idx] += Hrough;
				pMesh->Heff2[idx] += Hrough;

				energy += Hrough * (pMesh->M[idx] + pMesh->M2[idx]) / 2;
			}
		}
	}

	//average energy density
	if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / (2 * (pMesh->M.get_nonempty_cells()));
	else energy = 0;

	this->energy = energy;

	return energy;
}

#endif
