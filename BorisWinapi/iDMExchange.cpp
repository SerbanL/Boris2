#include "stdafx.h"
#include "iDMExchange.h"

#ifdef MODULE_IDMEXCHANGE

#include "Mesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "iDMExchangeCUDA.h"
#endif

//For z-axis symmetry breaking Dzyaloshinskii-Moriya exchange, we have :
//
//Hdm,ex = -2D/(mu0*Ms) * (dmz/dx, dmz/dy, -dmx/dx - dmy/dy)
//
//Hdm,ex adds to the usual isotropic direct exchange term, energy calculated using the total exchange field as usual.
//
//Here D is the DM exchange constant (J/m^2)

iDMExchange::iDMExchange(Mesh *pMesh_) :
	Modules(),
	ExchangeBase(pMesh_),
	ProgramStateNames(this, {}, {})
{
	pMesh = pMesh_;

	error_on_create = UpdateConfiguration(UPDATECONFIG_FORCEUPDATE);

	//-------------------------- Is CUDA currently enabled?

	//If cuda is enabled we also need to make the cuda module version
	if (pMesh->cudaEnabled) {

		if (!error_on_create) error_on_create = SwitchCUDAState(true);
	}
}

iDMExchange::~iDMExchange()
{
}

BError iDMExchange::Initialize(void)
{
	BError error(CLASS_STR(iDMExchange));

	error = ExchangeBase::Initialize();

	if (!error) initialized = true;

	return error;
}

BError iDMExchange::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(iDMExchange));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError iDMExchange::MakeCUDAModule(void)
{
	BError error(CLASS_STR(iDMExchange));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new iDMExchangeCUDA(pMesh->pMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double iDMExchange::UpdateField(void)
{
	double energy = 0;

	INT3 n = pMesh->n;

	///////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////// FERROMAGNETIC MESH /////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	if (pMesh->GetMeshType() == MESH_FERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy) 
		for (int idx = 0; idx < n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				double Ms = pMesh->Ms;
				double A = pMesh->A;
				double D = pMesh->D;
				pMesh->update_parameters_mcoarse(idx, pMesh->A, A, pMesh->D, D, pMesh->Ms, Ms);

				double Aconst = 2 * A / (MU0 * Ms * Ms);
				double Dconst = -2 * D / (MU0 * Ms * Ms);

				DBL3 Hexch;

				if (pMesh->M.is_plane_interior(idx)) {

					//interior point : can use cheaper neu versions

					//direct exchange contribution
					Hexch = Aconst * pMesh->M.delsq_neu(idx);

					//Dzyaloshinskii-Moriya interfacial exchange contribution

					//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
					DBL33 Mdiff = pMesh->M.grad_neu(idx);

					//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
					Hexch += Dconst * DBL3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);
				}
				else {

					//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
					DBL3 bnd_dm_dx = (D / (2 * A)) * DBL3(pMesh->M[idx].z, 0, -pMesh->M[idx].x);
					DBL3 bnd_dm_dy = (D / (2 * A)) * DBL3(0, pMesh->M[idx].z, -pMesh->M[idx].y);
					DBL33 bnd_nneu = DBL33(bnd_dm_dx, bnd_dm_dy, DBL3());

					//direct exchange contribution
					//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_nneu evaluates to zero in the CMBND coupling direction.
					Hexch = Aconst * pMesh->M.delsq_nneu(idx, bnd_nneu);

					//Dzyaloshinskii-Moriya interfacial exchange contribution

					//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
					//For cmbnd cells grad_nneu does not evaluate to zero in the CMBND coupling direction, but sided differentials are used - when setting values at CMBND cells for exchange coupled meshes must correct for this.
					DBL33 Mdiff = pMesh->M.grad_nneu(idx, bnd_nneu);

					//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
					Hexch += Dconst * DBL3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);
				}

				pMesh->Heff[idx] += Hexch;

				energy += pMesh->M[idx] * Hexch;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////// ANTIFERROMAGNETIC MESH ///////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	else if (pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC) {

#pragma omp parallel for reduction(+:energy) 
		for (int idx = 0; idx < n.dim(); idx++) {

			if (pMesh->M.is_not_empty(idx)) {

				DBL2 Ms_AFM = pMesh->Ms_AFM;
				DBL2 A_AFM = pMesh->A_AFM;
				double A12 = pMesh->A12;
				DBL2 D_AFM = pMesh->D_AFM;
				pMesh->update_parameters_mcoarse(idx, pMesh->A_AFM, A_AFM, pMesh->A12, A12, pMesh->D_AFM, D_AFM, pMesh->Ms_AFM, Ms_AFM);

				DBL2 Aconst = 2 * A_AFM / (MU0 * (Ms_AFM & Ms_AFM));
				DBL2 Dconst = -2 * D_AFM / (MU0 * (Ms_AFM & Ms_AFM));

				DBL3 Hexch, Hexch2;

				if (pMesh->M.is_plane_interior(idx)) {

					//interior point : can use cheaper neu versions

					//direct exchange contribution + AFM contribution
					Hexch = Aconst.i * pMesh->M.delsq_neu(idx) + (4 * A12 / (MU0*Ms_AFM.i*Ms_AFM.j)) * pMesh->M2[idx];
					Hexch2 = Aconst.j * pMesh->M2.delsq_neu(idx) + (4 * A12 / (MU0*Ms_AFM.i*Ms_AFM.j)) * pMesh->M[idx];

					//Dzyaloshinskii-Moriya interfacial exchange contribution

					//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
					DBL33 Mdiff = pMesh->M.grad_neu(idx);

					//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
					Hexch += Dconst.i * DBL3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);

					Mdiff = pMesh->M2.grad_neu(idx);
					Hexch2 += Dconst.j * DBL3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);
				}
				else {

					//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
					DBL3 bnd_dm_dx = (D_AFM.i / (2 * A_AFM.i)) * DBL3(pMesh->M[idx].z, 0, -pMesh->M[idx].x);
					DBL3 bnd_dm_dy = (D_AFM.i / (2 * A_AFM.i)) * DBL3(0, pMesh->M[idx].z, -pMesh->M[idx].y);
					DBL33 bnd_nneu = DBL33(bnd_dm_dx, bnd_dm_dy, DBL3());

					//direct exchange contribution + AFM contribution
					//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_nneu evaluates to zero in the CMBND coupling direction.
					Hexch = Aconst.i * pMesh->M.delsq_nneu(idx, bnd_nneu) + (4 * A12 / (MU0*Ms_AFM.i*Ms_AFM.j)) * pMesh->M2[idx];

					//Dzyaloshinskii-Moriya interfacial exchange contribution

					//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
					//For cmbnd cells grad_nneu does not evaluate to zero in the CMBND coupling direction, but sided differentials are used - when setting values at CMBND cells for exchange coupled meshes must correct for this.
					DBL33 Mdiff = pMesh->M.grad_nneu(idx, bnd_nneu);

					//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
					Hexch += Dconst.i * DBL3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);

					//same thing on sub-lattice B (2)

					bnd_dm_dx = (D_AFM.j / (2 * A_AFM.j)) * DBL3(pMesh->M2[idx].z, 0, -pMesh->M2[idx].x);
					bnd_dm_dy = (D_AFM.j / (2 * A_AFM.j)) * DBL3(0, pMesh->M2[idx].z, -pMesh->M2[idx].y);
					bnd_nneu = DBL33(bnd_dm_dx, bnd_dm_dy, DBL3());

					Hexch2 = Aconst.j * pMesh->M2.delsq_nneu(idx, bnd_nneu) + (4 * A12 / (MU0*Ms_AFM.i*Ms_AFM.j)) * pMesh->M[idx];

					Mdiff = pMesh->M2.grad_nneu(idx, bnd_nneu);
					Hexch2 += Dconst.j * DBL3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);
				}

				pMesh->Heff[idx] += Hexch;
				pMesh->Heff2[idx] += Hexch2;

				energy += (pMesh->M[idx] * Hexch + pMesh->M2[idx] * Hexch2) / 2;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////// COUPLING ACROSS MULTIPLE MESHES ///////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	//if exchange coupling across multiple meshes, this is calculation method to use
	std::function<double(int, int, DBL3, DBL3, DBL3, Mesh&, Mesh&)> calculate_coupling = [&](int cell1_idx, int cell2_idx, DBL3 relpos_m1, DBL3 stencil, DBL3 hshift_primary, Mesh& Mesh_pri, Mesh& Mesh_sec) -> double {

		double energy_ = 0.0;

		double hR = hshift_primary.norm();
		double hRsq = hR * hR;

		if (Mesh_pri.GetMeshType() == MESH_ANTIFERROMAGNETIC && Mesh_sec.GetMeshType() == MESH_ANTIFERROMAGNETIC) {

			//Both meshes antiferromagnetic : both sub-lattices couple, but do not include anti AF coupling here as it has already been included in the main loop above
			//here we only compute differential operators across a boundary.

			DBL2 Ms_AFM = Mesh_pri.Ms_AFM;
			DBL2 A_AFM = Mesh_pri.A_AFM;
			DBL2 D_AFM = Mesh_pri.D_AFM;
			Mesh_pri.update_parameters_mcoarse(cell1_idx, Mesh_pri.A_AFM, A_AFM, Mesh_pri.D_AFM, D_AFM, Mesh_pri.Ms_AFM, Ms_AFM);

			DBL3 Hexch, Hexch_B;

			//values at cells -1, 1
			DBL3 M_1 = Mesh_pri.M[cell1_idx];
			DBL3 M_m1 = Mesh_sec.M.weighted_average(relpos_m1, stencil);

			DBL3 M_1_B = Mesh_pri.M2[cell1_idx];
			DBL3 M_m1_B = Mesh_sec.M2.weighted_average(relpos_m1, stencil);

			if (cell2_idx < Mesh_pri.n.dim() && Mesh_pri.M.is_not_empty(cell2_idx)) {

				//cell2_idx is valid and M is not empty there
				DBL3 M_2 = Mesh_pri.M[cell2_idx];
				DBL3 M_2_B = Mesh_pri.M2[cell2_idx];

				//set effective field value contribution at cell 1 : direct exchange coupling + AFM coupling
				Hexch = (2 * A_AFM.i / (MU0*Ms_AFM.i*Ms_AFM.i)) * (M_2 + M_m1 - 2 * M_1) / hRsq;
				Hexch_B = (2 * A_AFM.j / (MU0*Ms_AFM.j*Ms_AFM.j)) * (M_2_B + M_m1_B - 2 * M_1_B) / hRsq;

				//add iDMI contributions at CMBND cells, correcting for the sided differentials already applied here
				//the contributions are different depending on the CMBND coupling direction
				if (IsNZ(hshift_primary.x)) {

					//along x
					Hexch += (-2 * D_AFM.i / (MU0*Ms_AFM.i*Ms_AFM.i)) * DBL3(M_2.z + M_m1.z - 2 * M_1.z, 0, -M_2.x - M_m1.x + 2 * M_1.x) / (2 * hshift_primary.x);
					Hexch_B += (-2 * D_AFM.j / (MU0*Ms_AFM.j*Ms_AFM.j)) * DBL3(M_2_B.z + M_m1_B.z - 2 * M_1_B.z, 0, -M_2_B.x - M_m1_B.x + 2 * M_1_B.x) / (2 * hshift_primary.x);
				}
				else if (IsNZ(hshift_primary.y)) {

					//along y
					Hexch += (-2 * D_AFM.i / (MU0*Ms_AFM.i*Ms_AFM.i)) * DBL3(0, M_2.z + M_m1.z - 2 * M_1.z, -M_2.y - M_m1.y + 2 * M_1.y) / (2 * hshift_primary.y);
					Hexch_B += (-2 * D_AFM.j / (MU0*Ms_AFM.j*Ms_AFM.j)) * DBL3(0, M_2_B.z + M_m1_B.z - 2 * M_1_B.z, -M_2_B.y - M_m1_B.y + 2 * M_1_B.y) / (2 * hshift_primary.y);
				}
			}
			else {

				//set effective field value contribution at cell 1 : direct exchange coupling + AFM coupling
				Hexch = (2 * A_AFM.i / (MU0*Ms_AFM.i*Ms_AFM.i)) * (M_m1 - M_1) / hRsq;
				Hexch_B = (2 * A_AFM.j / (MU0*Ms_AFM.j*Ms_AFM.j)) * (M_m1_B - M_1_B) / hRsq;

				//no iDMI contribution here
			}

			Mesh_pri.Heff[cell1_idx] += Hexch;
			Mesh_pri.Heff2[cell1_idx] += Hexch_B;

			energy_ = (M_1 * Hexch + M_1_B * Hexch_B) / 2;
		}

		//FM to FM
		else if (Mesh_pri.GetMeshType() == MESH_FERROMAGNETIC && Mesh_sec.GetMeshType() == MESH_FERROMAGNETIC) {

			//both meshes are ferromagnetic
			//here we only compute differential operators across a boundary.

			double Ms, A, D;
			Ms = Mesh_pri.Ms;
			A = Mesh_pri.A;
			D = Mesh_pri.D;
			Mesh_pri.update_parameters_mcoarse(cell1_idx, Mesh_pri.A, A, Mesh_pri.D, D, Mesh_pri.Ms, Ms);

			DBL3 Hexch;

			//values at cells -1, 1
			DBL3 M_1 = Mesh_pri.M[cell1_idx];
			DBL3 M_m1 = Mesh_sec.M.weighted_average(relpos_m1, stencil);

			if (cell2_idx < Mesh_pri.n.dim() && Mesh_pri.M.is_not_empty(cell2_idx)) {

				//cell2_idx is valid and M is not empty there
				DBL3 M_2 = Mesh_pri.M[cell2_idx];

				//set effective field value contribution at cell 1 : direct exchange coupling
				Hexch = (2 * A / (MU0*Ms*Ms)) * (M_2 + M_m1 - 2 * M_1) / hRsq;
				
				//add iDMI contributions at CMBND cells, correcting for the sided differentials already applied here
				//the contributions are different depending on the CMBND coupling direction
				if (IsNZ(hshift_primary.x)) {

					//along x
					Hexch += (-2 * D / (MU0*Ms*Ms)) * DBL3(M_2.z + M_m1.z - 2 * M_1.z, 0, -M_2.x - M_m1.x + 2 * M_1.x) / (2 * hshift_primary.x);
				}
				else if (IsNZ(hshift_primary.y)) {

					//along y
					Hexch += (-2 * D / (MU0*Ms*Ms)) * DBL3(0, M_2.z + M_m1.z - 2 * M_1.z, -M_2.y - M_m1.y + 2 * M_1.y) / (2 * hshift_primary.y);
				}
			}
			else {

				//set effective field value contribution at cell 1 : direct exchange coupling
				Hexch = (2 * A / (MU0*Ms*Ms)) * (M_m1 - M_1) / hRsq;
			}

			Mesh_pri.Heff[cell1_idx] += Hexch;

			energy_ = M_1 * Hexch;
		}

		return energy_;
	};

	//if exchange coupled to other meshes calculate the exchange field at marked cmbnd cells and accumulate energy density contribution
	if (pMesh->GetMeshExchangeCoupling()) CalculateExchangeCoupling(energy, calculate_coupling);

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////// FINAL ENERGY DENSITY VALUE //////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	//average energy density, this is correct, see notes - e_ex ~= -(mu0/2) M.H as can be derived from the full expression for e_ex. It is not -mu0 M.H, which is the energy density for a dipole in a field !
	if (pMesh->M.get_nonempty_cells()) energy *= -MU0 / (2 * (pMesh->M.get_nonempty_cells()));
	else energy = 0;

	this->energy = energy;

	return this->energy;
}

#endif