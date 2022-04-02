#include "stdafx.h"
#include "iDMExchange.h"

#ifdef MODULE_COMPILATION_IDMEXCHANGE

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

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect, 
		(MOD_)pMesh->Get_ActualModule_Heff_Display() == MOD_IDMEXCHANGE || pMesh->IsOutputDataSet_withRect(DATA_E_EXCH),
		(MOD_)pMesh->Get_ActualModule_Energy_Display() == MOD_IDMEXCHANGE || pMesh->IsOutputDataSet_withRect(DATA_E_EXCH),
		pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC);
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

	SZ3 n = pMesh->n;

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

				DBL3 Hexch_A, Hexch_D;

				if (pMesh->M.is_plane_interior(idx)) {

					//interior point : can use cheaper neu versions

					//direct exchange contribution
					if (pMesh->base_temperature > 0.0 && pMesh->T_Curie > 0.0) {

						//for finite temperature simulations the magnetization length may have a spatial variation
						//this will not affect the transverse torque (mxH), but will affect the longitudinal term in the sLLB equation (m.H) and cannot be neglected when close to Tc.

						DBL33 Mg = pMesh->M.grad_neu(idx);
						DBL3 dMdx = Mg.x, dMdy = Mg.y, dMdz = Mg.z;

						double delsq_Msq = 2 * pMesh->M[idx] * (pMesh->M.dxx_neu(idx) + pMesh->M.dyy_neu(idx) + pMesh->M.dzz_neu(idx)) + 2 * (dMdx * dMdx + dMdy * dMdy + dMdz * dMdz);
						double Mnorm = pMesh->M[idx].norm();
						Hexch_A = Aconst * (pMesh->M.delsq_neu(idx) - pMesh->M[idx] * delsq_Msq / (2 * Mnorm*Mnorm));
					}
					else {

						//zero temperature simulations : magnetization length could still vary but will only affect mxH term, so not needed for 0K simulations.
						Hexch_A = Aconst * pMesh->M.delsq_neu(idx);
					}

					//Dzyaloshinskii-Moriya interfacial exchange contribution

					//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
					DBL33 Mdiff = pMesh->M.grad_neu(idx);

					//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
					Hexch_D = Dconst * DBL3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);
				}
				else {

					//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
					DBL3 bnd_dm_dx = (D / (2 * A)) * DBL3(pMesh->M[idx].z, 0, -pMesh->M[idx].x);
					DBL3 bnd_dm_dy = (D / (2 * A)) * DBL3(0, pMesh->M[idx].z, -pMesh->M[idx].y);
					DBL33 bnd_nneu = DBL33(bnd_dm_dx, bnd_dm_dy, DBL3());

					//direct exchange contribution
					//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_nneu evaluates to zero in the CMBND coupling direction.
					if (pMesh->base_temperature > 0.0 && pMesh->T_Curie > 0.0) {

						//for finite temperature simulations the magnetization length may have a spatial variation
						//this will not affect the transverse torque (mxH), but will affect the longitudinal term in the sLLB equation (m.H) and cannot be neglected when close to Tc.

						DBL33 Mg = pMesh->M.grad_nneu(idx, bnd_nneu);
						DBL3 dMdx = Mg.x, dMdy = Mg.y, dMdz = Mg.z;

						double delsq_Msq = 2 * pMesh->M[idx] * (pMesh->M.dxx_nneu(idx, bnd_nneu) + pMesh->M.dyy_nneu(idx, bnd_nneu) + pMesh->M.dzz_nneu(idx, bnd_nneu)) + 2 * (dMdx * dMdx + dMdy * dMdy + dMdz * dMdz);
						double Mnorm = pMesh->M[idx].norm();
						Hexch_A = Aconst * (pMesh->M.delsq_nneu(idx, bnd_nneu) - pMesh->M[idx] * delsq_Msq / (2 * Mnorm*Mnorm));
					}
					else {

						//zero temperature simulations : magnetization length could still vary but will only affect mxH term, so not needed for 0K simulations.
						Hexch_A = Aconst * pMesh->M.delsq_nneu(idx, bnd_nneu);
					}

					//Dzyaloshinskii-Moriya interfacial exchange contribution

					//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
					//For cmbnd cells grad_nneu does not evaluate to zero in the CMBND coupling direction, but sided differentials are used - when setting values at CMBND cells for exchange coupled meshes must correct for this.
					DBL33 Mdiff = pMesh->M.grad_nneu(idx, bnd_nneu);

					//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
					Hexch_D = Dconst * DBL3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);
				}

				pMesh->Heff[idx] += Hexch_A + Hexch_D;

				energy += pMesh->M[idx] * (Hexch_A + Hexch_D);

				//spatial dependence display of effective field and energy density
				if (Module_Heff.linear_size() && Module_energy.linear_size()) {
					
					if ((MOD_)pMesh->Get_Module_Heff_Display() == MOD_EXCHANGE) {

						//total : direct and DMI
						Module_Heff[idx] = Hexch_A + Hexch_D;
						Module_energy[idx] = -MU0 * (pMesh->M[idx] * (Hexch_A + Hexch_D)) / 2;
					}
					else {

						//just DMI
						Module_Heff[idx] = Hexch_D;
						Module_energy[idx] = -MU0 * (pMesh->M[idx] * Hexch_D) / 2;
					}
				}
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
				DBL2 Ah = pMesh->Ah;
				DBL2 Anh = pMesh->Anh;
				DBL2 D_AFM = pMesh->D_AFM;
				pMesh->update_parameters_mcoarse(idx, pMesh->A_AFM, A_AFM, pMesh->Ah, Ah, pMesh->Anh, Anh, pMesh->D_AFM, D_AFM, pMesh->Ms_AFM, Ms_AFM);

				DBL2 Aconst = 2 * A_AFM / (MU0 * (Ms_AFM & Ms_AFM));
				DBL2 Dconst = -2 * D_AFM / (MU0 * (Ms_AFM & Ms_AFM));

				DBL3 Hexch_A, Hexch_A2, Hexch_D, Hexch_D2;

				if (pMesh->M.is_plane_interior(idx)) {

					//interior point : can use cheaper neu versions

					//1. direct exchange contribution + AFM contribution
					DBL3 delsq_M_A = pMesh->M.delsq_neu(idx);
					DBL3 delsq_M_B = pMesh->M2.delsq_neu(idx);

					DBL2 M = DBL2(pMesh->M[idx].norm(), pMesh->M2[idx].norm());

					Hexch_A = Aconst.i * delsq_M_A + (-4 * Ah.i * (pMesh->M[idx] ^ (pMesh->M[idx] ^ pMesh->M2[idx])) / (M.i*M.i) + Anh.i * delsq_M_B) / (MU0*Ms_AFM.i*Ms_AFM.j);
					Hexch_A2 = Aconst.j * delsq_M_B + (-4 * Ah.j * (pMesh->M2[idx] ^ (pMesh->M2[idx] ^ pMesh->M[idx])) / (M.j*M.j) + Anh.j * delsq_M_A) / (MU0*Ms_AFM.i*Ms_AFM.j);

					//2. Dzyaloshinskii-Moriya interfacial exchange contribution

					//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
					DBL33 Mdiff_A = pMesh->M.grad_neu(idx);
					DBL33 Mdiff_B = pMesh->M2.grad_neu(idx);

					//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
					Hexch_D = Dconst.i * DBL3(Mdiff_A.x.z, Mdiff_A.y.z, -Mdiff_A.x.x - Mdiff_A.y.y);
					Hexch_D2 = Dconst.j * DBL3(Mdiff_B.x.z, Mdiff_B.y.z, -Mdiff_B.x.x - Mdiff_B.y.y);
				}
				else {

					DBL33 bndA_nneu, bndB_nneu;

					DBL2 nhconst = Anh / (2 * A_AFM);

					if (fabs(nhconst.i) != 1.0) {

						//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
						DBL3 bnd_dm_dx = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(pMesh->M[idx].z - nhconst.i * pMesh->M2[idx].z, 0, -pMesh->M[idx].x + nhconst.i * pMesh->M2[idx].x);
						DBL3 bnd_dm_dy = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(0, pMesh->M[idx].z - nhconst.i * pMesh->M2[idx].z, -pMesh->M[idx].y + nhconst.i * pMesh->M2[idx].y);
						
						bndA_nneu = DBL33(bnd_dm_dx, bnd_dm_dy, DBL3());
					}
					else {

						DBL3 bnd_dm_dx = (D_AFM.i / (4 * A_AFM.i)) * DBL3(pMesh->M[idx].z, 0, -pMesh->M[idx].x);
						DBL3 bnd_dm_dy = (D_AFM.i / (4 * A_AFM.i)) * DBL3(0, pMesh->M[idx].z, -pMesh->M[idx].y);

						bndA_nneu = DBL33(bnd_dm_dx, bnd_dm_dy, DBL3());
					}

					if (fabs(nhconst.j) != 1.0) {

						DBL3 bnd_dm_dx = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(pMesh->M2[idx].z - nhconst.j * pMesh->M[idx].z, 0, -pMesh->M2[idx].x + nhconst.j * pMesh->M[idx].x);
						DBL3 bnd_dm_dy = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(0, pMesh->M2[idx].z - nhconst.j * pMesh->M[idx].z, -pMesh->M2[idx].y + nhconst.j * pMesh->M[idx].y);

						bndB_nneu = DBL33(bnd_dm_dx, bnd_dm_dy, DBL3());
					}
					else {

						DBL3 bnd_dm_dx = (D_AFM.j / (4 * A_AFM.j)) * DBL3(pMesh->M2[idx].z, 0, -pMesh->M2[idx].x);
						DBL3 bnd_dm_dy = (D_AFM.j / (4 * A_AFM.j)) * DBL3(0, pMesh->M2[idx].z, -pMesh->M2[idx].y);

						bndB_nneu = DBL33(bnd_dm_dx, bnd_dm_dy, DBL3());
					}

					DBL3 delsq_M_A = pMesh->M.delsq_nneu(idx, bndA_nneu);
					DBL3 delsq_M_B = pMesh->M2.delsq_nneu(idx, bndB_nneu);

					DBL2 M = DBL2(pMesh->M[idx].norm(), pMesh->M2[idx].norm());

					//1. direct exchange contribution + AFM contribution

					//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_nneu evaluates to zero in the CMBND coupling direction.
					Hexch_A = Aconst.i * delsq_M_A + (-4 * Ah.i * (pMesh->M[idx] ^ (pMesh->M[idx] ^ pMesh->M2[idx])) / (M.i*M.i) + Anh.i * delsq_M_B) / (MU0*Ms_AFM.i*Ms_AFM.j);
					Hexch_A2 = Aconst.j * delsq_M_B + (-4 * Ah.j * (pMesh->M2[idx] ^ (pMesh->M2[idx] ^ pMesh->M[idx])) / (M.j*M.j) + Anh.j * delsq_M_A) / (MU0*Ms_AFM.i*Ms_AFM.j);

					//2. Dzyaloshinskii-Moriya interfacial exchange contribution

					//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
					//For cmbnd cells grad_nneu does not evaluate to zero in the CMBND coupling direction, but sided differentials are used - when setting values at CMBND cells for exchange coupled meshes must correct for this.
					DBL33 Mdiff_A = pMesh->M.grad_nneu(idx, bndA_nneu);
					DBL33 Mdiff_B = pMesh->M2.grad_nneu(idx, bndB_nneu);

					//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
					Hexch_D = Dconst.i * DBL3(Mdiff_A.x.z, Mdiff_A.y.z, -Mdiff_A.x.x - Mdiff_A.y.y);
					Hexch_D2 = Dconst.j * DBL3(Mdiff_B.x.z, Mdiff_B.y.z, -Mdiff_B.x.x - Mdiff_B.y.y);
				}

				pMesh->Heff[idx] += (Hexch_A + Hexch_D);
				pMesh->Heff2[idx] += (Hexch_A2 + Hexch_D2);

				energy += (pMesh->M[idx] * (Hexch_A + Hexch_D) + pMesh->M2[idx] * (Hexch_A2 + Hexch_D2)) / 2;

				//spatial dependence display of effective field and energy density
				if (Module_Heff.linear_size() && Module_energy.linear_size()) {

					if ((MOD_)pMesh->Get_Module_Heff_Display() == MOD_EXCHANGE) {

						//total : direct and DMI
						Module_Heff[idx] = Hexch_A + Hexch_D;
						Module_energy[idx] = -MU0 * (pMesh->M[idx] * (Hexch_A + Hexch_D)) / 2;
					}
					else {

						//just DMI
						Module_Heff[idx] = Hexch_D;
						Module_energy[idx] = -MU0 * (pMesh->M[idx] * Hexch_D) / 2;
					}
				}

				if (Module_Heff2.linear_size() && Module_energy2.linear_size()) {

					if ((MOD_)pMesh->Get_Module_Heff_Display() == MOD_EXCHANGE) {

						//total : direct and DMI
						Module_Heff2[idx] = Hexch_A2 + Hexch_D2;
						Module_energy2[idx] = -MU0 * (pMesh->M2[idx] * (Hexch_A2 + Hexch_D2)) / 2;
					}
					else {

						//just DMI
						Module_Heff2[idx] = Hexch_D2;
						Module_energy2[idx] = -MU0 * (pMesh->M2[idx] * Hexch_D2) / 2;
					}
				}
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
			DBL2 Anh = pMesh->Anh;
			Mesh_pri.update_parameters_mcoarse(cell1_idx, Mesh_pri.A_AFM, A_AFM, Mesh_pri.D_AFM, D_AFM, Mesh_pri.Ms_AFM, Ms_AFM, pMesh->Anh, Anh);

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

				DBL3 delsq_M_A = (M_2 + M_m1 - 2 * M_1) / hRsq;
				DBL3 delsq_M_B = (M_2_B + M_m1_B - 2 * M_1_B) / hRsq;

				//set effective field value contribution at cell 1 : direct exchange coupling + AFM coupling
				Hexch = (2 * A_AFM.i / (MU0*Ms_AFM.i*Ms_AFM.i)) * delsq_M_A + (Anh.i / (MU0*Ms_AFM.i*Ms_AFM.j)) * delsq_M_B;
				Hexch_B = (2 * A_AFM.j / (MU0*Ms_AFM.j*Ms_AFM.j)) * delsq_M_B + (Anh.j / (MU0*Ms_AFM.i*Ms_AFM.j)) * delsq_M_A;

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

				DBL3 delsq_M_A = (M_m1 - M_1) / hRsq;
				DBL3 delsq_M_B = (M_m1_B - M_1_B) / hRsq;

				//set effective field value contribution at cell 1 : direct exchange coupling + AFM coupling
				Hexch = (2 * A_AFM.i / (MU0*Ms_AFM.i*Ms_AFM.i)) * delsq_M_A + (Anh.i / (MU0*Ms_AFM.i*Ms_AFM.j)) * delsq_M_B;
				Hexch_B = (2 * A_AFM.j / (MU0*Ms_AFM.j*Ms_AFM.j)) * delsq_M_B + (Anh.j / (MU0*Ms_AFM.i*Ms_AFM.j)) * delsq_M_A;

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

//-------------------Energy methods

//FM Mesh
double iDMExchange::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index)) {

		double Ms = pMesh->Ms;
		double A = pMesh->A;
		double D = pMesh->D;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->A, A, pMesh->D, D, pMesh->Ms, Ms);

		double Aconst = 2 * A / (MU0 * Ms * Ms);
		double Dconst = -2 * D / (MU0 * Ms * Ms);

		auto Get_Energy = [&](void) -> double
		{
			DBL3 Hexch_A, Hexch_D;

			if (pMesh->M.is_plane_interior(spin_index)) {

				//interior point : can use cheaper neu versions

				//direct exchange contribution
				Hexch_A = Aconst * pMesh->M.delsq_neu(spin_index);

				//Dzyaloshinskii-Moriya interfacial exchange contribution

				//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
				DBL33 Mdiff = pMesh->M.grad_neu(spin_index);

				//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
				Hexch_D = Dconst * DBL3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);
			}
			else {

				//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
				DBL3 bnd_dm_dx = (D / (2 * A)) * DBL3(pMesh->M[spin_index].z, 0, -pMesh->M[spin_index].x);
				DBL3 bnd_dm_dy = (D / (2 * A)) * DBL3(0, pMesh->M[spin_index].z, -pMesh->M[spin_index].y);
				DBL33 bnd_nneu = DBL33(bnd_dm_dx, bnd_dm_dy, DBL3());

				//direct exchange contribution
				//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_nneu evaluates to zero in the CMBND coupling direction.
				Hexch_A = Aconst * pMesh->M.delsq_nneu(spin_index, bnd_nneu);

				//Dzyaloshinskii-Moriya interfacial exchange contribution

				//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
				//For cmbnd cells grad_nneu does not evaluate to zero in the CMBND coupling direction, but sided differentials are used - when setting values at CMBND cells for exchange coupled meshes must correct for this.
				DBL33 Mdiff = pMesh->M.grad_nneu(spin_index, bnd_nneu);

				//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
				Hexch_D = Dconst * DBL3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);
			}

			return pMesh->M[spin_index] * (Hexch_A + Hexch_D);
		};

		double energy_ = Get_Energy();

		if (Mnew != DBL3()) {

			//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
			Mnew.renormalize(pMesh->M[spin_index].norm());

			//new spin energy
			DBL3 Mold = pMesh->M[spin_index];
			pMesh->M[spin_index] = Mnew;
			double energynew_ = Get_Energy();
			pMesh->M[spin_index] = Mold;

			//do not divide by 2 as we are not double-counting here
			return -MU0 * pMesh->h.dim() * (energynew_ - energy_);
		}
		else return -MU0 * pMesh->h.dim() * energy_;
	}
	else return 0.0;
}

//AFM mesh
DBL2 iDMExchange::Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B)
{	
	if (pMesh->M.is_not_empty(spin_index) && pMesh->M2.is_not_empty(spin_index)) {

		DBL2 Ms_AFM = pMesh->Ms_AFM;
		DBL2 A_AFM = pMesh->A_AFM;
		DBL2 Ah = pMesh->Ah;
		DBL2 Anh = pMesh->Anh;
		DBL2 D_AFM = pMesh->D_AFM;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->A_AFM, A_AFM, pMesh->Ah, Ah, pMesh->Anh, Anh, pMesh->Ms_AFM, Ms_AFM, pMesh->D_AFM, D_AFM);

		DBL2 Aconst = 2 * A_AFM / (MU0 * (Ms_AFM & Ms_AFM));
		DBL2 Dconst = -2 * D_AFM / (MU0 * (Ms_AFM & Ms_AFM));

		auto Get_Energy = [&](void) -> DBL2
		{
			DBL3 Hexch_A, Hexch_A2, Hexch_D, Hexch_D2;

			if (pMesh->M.is_plane_interior(spin_index)) {

				//interior point : can use cheaper neu versions

				//1. direct exchange contribution + AFM contribution
				DBL3 delsq_M_A = pMesh->M.delsq_neu(spin_index);
				DBL3 delsq_M_B = pMesh->M2.delsq_neu(spin_index);

				Hexch_A = Aconst.i * delsq_M_A + (4 * Ah.i * pMesh->M2[spin_index] + Anh.i * delsq_M_B) / (MU0*Ms_AFM.i*Ms_AFM.j);
				Hexch_A2 = Aconst.j * delsq_M_B + (4 * Ah.j * pMesh->M[spin_index] + Anh.j * delsq_M_A) / (MU0*Ms_AFM.i*Ms_AFM.j);

				//2. Dzyaloshinskii-Moriya interfacial exchange contribution

				//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
				DBL33 Mdiff_A = pMesh->M.grad_neu(spin_index);
				DBL33 Mdiff_B = pMesh->M2.grad_neu(spin_index);

				//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
				Hexch_D = Dconst.i * DBL3(Mdiff_A.x.z, Mdiff_A.y.z, -Mdiff_A.x.x - Mdiff_A.y.y);
				Hexch_D2 = Dconst.j * DBL3(Mdiff_B.x.z, Mdiff_B.y.z, -Mdiff_B.x.x - Mdiff_B.y.y);
			}
			else {

				DBL33 bndA_nneu, bndB_nneu;

				DBL2 nhconst = Anh / (2 * A_AFM);

				if (fabs(nhconst.i) != 1.0) {

					//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
					DBL3 bnd_dm_dx = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(pMesh->M[spin_index].z - nhconst.i * pMesh->M2[spin_index].z, 0, -pMesh->M[spin_index].x + nhconst.i * pMesh->M2[spin_index].x);
					DBL3 bnd_dm_dy = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(0, pMesh->M[spin_index].z - nhconst.i * pMesh->M2[spin_index].z, -pMesh->M[spin_index].y + nhconst.i * pMesh->M2[spin_index].y);

					bndA_nneu = DBL33(bnd_dm_dx, bnd_dm_dy, DBL3());
				}
				else {

					DBL3 bnd_dm_dx = (D_AFM.i / (4 * A_AFM.i)) * DBL3(pMesh->M[spin_index].z, 0, -pMesh->M[spin_index].x);
					DBL3 bnd_dm_dy = (D_AFM.i / (4 * A_AFM.i)) * DBL3(0, pMesh->M[spin_index].z, -pMesh->M[spin_index].y);

					bndA_nneu = DBL33(bnd_dm_dx, bnd_dm_dy, DBL3());
				}

				if (fabs(nhconst.j) != 1.0) {

					DBL3 bnd_dm_dx = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(pMesh->M2[spin_index].z - nhconst.j * pMesh->M[spin_index].z, 0, -pMesh->M2[spin_index].x + nhconst.j * pMesh->M[spin_index].x);
					DBL3 bnd_dm_dy = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(0, pMesh->M2[spin_index].z - nhconst.j * pMesh->M[spin_index].z, -pMesh->M2[spin_index].y + nhconst.j * pMesh->M[spin_index].y);

					bndB_nneu = DBL33(bnd_dm_dx, bnd_dm_dy, DBL3());
				}
				else {

					DBL3 bnd_dm_dx = (D_AFM.j / (4 * A_AFM.j)) * DBL3(pMesh->M2[spin_index].z, 0, -pMesh->M2[spin_index].x);
					DBL3 bnd_dm_dy = (D_AFM.j / (4 * A_AFM.j)) * DBL3(0, pMesh->M2[spin_index].z, -pMesh->M2[spin_index].y);

					bndB_nneu = DBL33(bnd_dm_dx, bnd_dm_dy, DBL3());
				}

				DBL3 delsq_M_A = pMesh->M.delsq_nneu(spin_index, bndA_nneu);
				DBL3 delsq_M_B = pMesh->M2.delsq_nneu(spin_index, bndB_nneu);

				//1. direct exchange contribution + AFM contribution

				//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_nneu evaluates to zero in the CMBND coupling direction.
				Hexch_A = Aconst.i * delsq_M_A + (4 * Ah.i * pMesh->M2[spin_index] + Anh.i * delsq_M_B) / (MU0*Ms_AFM.i*Ms_AFM.j);
				Hexch_A2 = Aconst.j * delsq_M_B + (4 * Ah.j * pMesh->M[spin_index] + Anh.j * delsq_M_A) / (MU0*Ms_AFM.i*Ms_AFM.j);

				//2. Dzyaloshinskii-Moriya interfacial exchange contribution

				//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
				//For cmbnd cells grad_nneu does not evaluate to zero in the CMBND coupling direction, but sided differentials are used - when setting values at CMBND cells for exchange coupled meshes must correct for this.
				DBL33 Mdiff_A = pMesh->M.grad_nneu(spin_index, bndA_nneu);
				DBL33 Mdiff_B = pMesh->M2.grad_nneu(spin_index, bndB_nneu);

				//Hdm, ex = -2D / (mu0*Ms) * (dmz / dx, dmz / dy, -dmx / dx - dmy / dy)
				Hexch_D = Dconst.i * DBL3(Mdiff_A.x.z, Mdiff_A.y.z, -Mdiff_A.x.x - Mdiff_A.y.y);
				Hexch_D2 = Dconst.j * DBL3(Mdiff_B.x.z, Mdiff_B.y.z, -Mdiff_B.x.x - Mdiff_B.y.y);
			}

			return DBL2(pMesh->M[spin_index] * (Hexch_A + Hexch_D), pMesh->M2[spin_index] * (Hexch_A2 + Hexch_D2));
		};

		DBL2 energy_ = Get_Energy();

		if (Mnew_A != DBL3() && Mnew_B != DBL3()) {

			//NOTE : here we only need the change in energy due to spin rotation only. Thus the longitudinal part, which is dependent on spin length only, cancels out. Enforce this by making Mnew length same as old one.
			Mnew_A.renormalize(pMesh->M[spin_index].norm());
			Mnew_B.renormalize(pMesh->M2[spin_index].norm());

			DBL3 Mold_A = pMesh->M[spin_index];
			DBL3 Mold_B = pMesh->M2[spin_index];

			pMesh->M[spin_index] = Mnew_A;
			pMesh->M2[spin_index] = Mnew_B;

			DBL2 energynew_ = Get_Energy();
			
			pMesh->M[spin_index] = Mold_A;
			pMesh->M2[spin_index] = Mold_B;

			//do not divide by 2 as we are not double-counting here
			return -MU0 * pMesh->h.dim() * (energynew_ - energy_);
		}
		//If Mnew is null then this method is used to obtain current energy only, not energy change
		else return -MU0 * pMesh->h.dim() * energy_;
	}
	else return DBL2();
}

//-------------------Torque methods

DBL3 iDMExchange::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return pModuleCUDA->GetTorque(avRect);
#endif

	return CalculateTorque(pMesh->M, avRect);
}

#endif