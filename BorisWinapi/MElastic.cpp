#include "stdafx.h"
#include "MElastic.h"

#ifdef MODULE_MELASTIC

#include "Mesh_Ferromagnetic.h"
#include "MeshParamsControl.h"

#include "SuperMesh.h"

#if COMPILECUDA == 1
#include "MElasticCUDA.h"
#endif

MElastic::MElastic(Mesh *pMesh_) :
	Modules(),
	ProgramStateNames(this, { VINFO(Tsig) }, {})
{
	pMesh = pMesh_;

	pSMesh = pMesh->pSMesh;

	Tsig = DBL3(0, 0, 0);

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

MElastic::~MElastic()
{
	//free memory used for heat solver
	pMesh->u_disp.clear();
	pMesh->strain_diag.clear();
	pMesh->strain_odiag.clear();
}

BError MElastic::Initialize(void)
{
	BError error(CLASS_STR(MElastic));

	initialized = true;

	return error;
}

BError MElastic::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(MElastic));

	Uninitialize();

	bool success = true;

	if (ucfg::check_cfgflags(cfgMessage, UPDATECONFIG_MESHSHAPECHANGE, UPDATECONFIG_MESHCHANGE)) {

		//make sure the cellsize divides the mesh rectangle
		pMesh->n_m = round(pMesh->meshRect / pMesh->h_m);
		if (pMesh->n_m.x < 2) pMesh->n_m.x = 2;
		if (pMesh->n_m.y < 2) pMesh->n_m.y = 2;
		if (pMesh->n_m.z < 2) pMesh->n_m.z = 2;
		pMesh->h_m = pMesh->meshRect / pMesh->n_m;

		//make sure correct memory is assigned for heat solver quantities

		if (pMesh->M.linear_size()) {

			//in a ferromagnet the magnetization sets the shape only on initialization. If already initialized then shape already set.
			if (pMesh->u_disp.linear_size()) {

				success = pMesh->u_disp.resize(pMesh->h_m, pMesh->meshRect);
			}
			else {

				success = pMesh->u_disp.assign(pMesh->h_m, pMesh->meshRect, DBL3(), pMesh->M);
			}
		}
		else {

			//the shape is set directly
			if (pMesh->u_disp.linear_size()) {

				success = pMesh->u_disp.resize(pMesh->h_m, pMesh->meshRect);
			}
			else {

				success = pMesh->u_disp.assign(pMesh->h_m, pMesh->meshRect, DBL3());
			}
		}

		//strain tensor - set empty cells using information in u_disp
		if (pMesh->strain_diag.linear_size()) {

			success &= pMesh->strain_diag.resize(pMesh->h_m, pMesh->meshRect, pMesh->u_disp);
			success &= pMesh->strain_odiag.resize(pMesh->h_m, pMesh->meshRect, pMesh->u_disp);
		}
		else {

			success &= pMesh->strain_diag.assign(pMesh->h_m, pMesh->meshRect, DBL3(), pMesh->u_disp);
			success &= pMesh->strain_odiag.assign(pMesh->h_m, pMesh->meshRect, DBL3(), pMesh->u_disp);
		}
	}

	if (!success) return error(BERROR_OUTOFMEMORY_CRIT);

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

void MElastic::UpdateConfiguration_Values(UPDATECONFIG_ cfgMessage)
{
	if (cfgMessage == UPDATECONFIG_TEQUATION_CONSTANTS) {

		
	}
	else if (cfgMessage == UPDATECONFIG_TEQUATION_CLEAR) {

		
	}

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		pModuleCUDA->UpdateConfiguration_Values(cfgMessage);
	}
#endif
}

BError MElastic::MakeCUDAModule(void)
{
	BError error(CLASS_STR(MElastic));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new MElasticCUDA(pMesh, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double MElastic::UpdateField(void)
{
	/////////////////////////////////////////
	// Fixed set stress
	/////////////////////////////////////////

	double energy = 0;

	//Energy density formula:

	//Let e1, e2, e3 be the orthogonal cubic axes, e.g. x, y, z in the simplest case

	//diagonal terms in stress tensor :
	//let Sd = (exx, eyy, ezz) be the diagonal terms.

	//off-diagonal terms in stress tensor (remember it is symmetric so we only have 3 of these):
	//let Sod = (eyz, exz, exy)

	//Then:

	//energy density from diagonal terms:
	//emel_d = B1 * [ (m.e1)^2*(Sd.e1) + (m.e2)^2*(Sd.e2) + (m.e3)^2*(Sd.e3) ]

	//energy density from off-diagonal terms (but remember this is zero for uniform strains):
	//emel_od = 2 * B2 * [ (m.e1)*(m.e2)*(Sod.e3) + (m.e1)*(m.e3)*(Sod.e2) + (m.e2)*(m.e3)*(Sod.e1) ]

	//Obtain Hmel as usual : Hmel = -1/mu0Ms * de/dm, where exchange stiffness-related terms are ignored here.

#pragma omp parallel for reduction(+:energy)
	for (int idx = 0; idx < pMesh->n.dim(); idx++) {

		if (pMesh->M.is_not_empty(idx)) {

			double Ms = pMesh->Ms;
			DBL3 mcanis_ea1 = pMesh->mcanis_ea1;
			DBL3 mcanis_ea2 = pMesh->mcanis_ea2;
			DBL2 MEc = pMesh->MEc;
			pMesh->update_parameters_mcoarse(idx, pMesh->Ms, Ms, pMesh->MEc, MEc, pMesh->mcanis_ea1, mcanis_ea1, pMesh->mcanis_ea2, mcanis_ea2);

			//vector product of ea1 and ea2 : the third orthogonal axis
			DBL3 mcanis_ea3 = mcanis_ea1 ^ mcanis_ea2;

			DBL3 position = pMesh->M.cellidx_to_position(idx);
			//xx, yy, zz
			DBL3 Sd = pMesh->strain_diag[position];
			//yz, xz, xy
			DBL3 Sod = pMesh->strain_odiag[position];
			
			//normalised magnetisation
			//Magneto-elastic term here applicable for a cubic crystal. We use the mcanis_ea1 and mcanis_ea2 axes to fix the cubic lattice orientation, thus rotate the m, Sd and Sod vectors.
			
			DBL3 m = DBL3(pMesh->M[idx] * mcanis_ea1, pMesh->M[idx] * mcanis_ea2, pMesh->M[idx] * mcanis_ea3) / Ms;
			Sd = DBL3(Sd * mcanis_ea1, Sd * mcanis_ea2, Sd * mcanis_ea3);
			Sod = DBL3(Sod * mcanis_ea1, Sod * mcanis_ea2, Sod * mcanis_ea3);

			DBL3 Hmel_1 = (-2.0 * MEc.i / (MU0 * Ms)) * DBL3(
				m.x*Sd.x*mcanis_ea1.x + m.y*Sd.y*mcanis_ea2.x + m.z*Sd.z*mcanis_ea3.x,
				m.x*Sd.x*mcanis_ea1.y + m.y*Sd.y*mcanis_ea2.y + m.z*Sd.z*mcanis_ea3.y,
				m.x*Sd.x*mcanis_ea1.z + m.y*Sd.y*mcanis_ea2.z + m.z*Sd.z*mcanis_ea3.z);

			DBL3 Hmel_2 = (-2.0 * MEc.j / (MU0 * Ms)) * DBL3(
				Sod.z * (mcanis_ea1.x*m.y + mcanis_ea2.x*m.x) + Sod.y * (mcanis_ea1.x*m.z + mcanis_ea3.x*m.x) + Sod.x * (mcanis_ea2.x*m.z + mcanis_ea3.x*m.y),
				Sod.z * (mcanis_ea1.y*m.y + mcanis_ea2.y*m.x) + Sod.y * (mcanis_ea1.y*m.z + mcanis_ea3.y*m.x) + Sod.x * (mcanis_ea2.y*m.z + mcanis_ea3.y*m.y),
				Sod.z * (mcanis_ea1.z*m.y + mcanis_ea2.z*m.x) + Sod.y * (mcanis_ea1.z*m.z + mcanis_ea3.z*m.x) + Sod.x * (mcanis_ea2.z*m.z + mcanis_ea3.z*m.y));

			pMesh->Heff[idx] += Hmel_1 + Hmel_2;

			energy += pMesh->M[idx] * (Hmel_1 + Hmel_2);
		}
	}
	
	if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / (2 * pMesh->M.get_nonempty_cells());
	else energy = 0;

	this->energy = energy;

	return this->energy;
}

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
BError MElastic::Load_Displacement_OVF2(string fileName)
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
BError MElastic::Load_Strain_OVF2(string fileName_Diag, string fileName_ODiag)
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

//----------------------------------------------- Computational Helpers

//compute strain form mechanical displacement
void MElastic::Calculate_Strain(void)
{
#pragma omp parallel for
	for (int idx = 0; idx < pMesh->n_m.dim(); idx++) {

		//get all 9 first-order differentials of u
		DBL33 grad_u = pMesh->u_disp.grad_neu(idx);

		//diagonal components
		pMesh->strain_diag[idx] = DBL3(grad_u.x.x, grad_u.y.y, grad_u.z.z);

		//off-diagonal components (yz, xz, xy)
		pMesh->strain_odiag[idx] = 0.5 * DBL3(grad_u.y.z + grad_u.z.y, grad_u.x.z + grad_u.z.x, grad_u.x.y + grad_u.y.x);
	}
}

#endif
