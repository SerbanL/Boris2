#include "stdafx.h"
#include "viDMExchange.h"

#ifdef MODULE_COMPILATION_VIDMEXCHANGE

#include "Mesh.h"
#include "MeshParamsControl.h"

#if COMPILECUDA == 1
#include "viDMExchangeCUDA.h"
#endif

//For z-axis symmetry DMI, we have :
//
//Hdm,ex,z = -2D/(mu0*Ms) * (dmz/dx, dmz/dy, -dmx/dx - dmy/dy)
//
//Hdm,ex adds to the usual isotropic direct exchange term, energy calculated using the total exchange field as usual.
//
//Here D is the DM exchange constant (J/m^2)
//
//For y-axis symmetry DMI, we have:
//Hdm,ex,y = -2D/(mu0*Ms) * (dmy/dx, -dmx/dx - dmz/dz, dmy/dz)

//For x-axis symmetry DMI, we have:
//Hdm,ex,x = -2D/(mu0*Ms) * (-dmy/dy - dmz/dz, dmx/dy, dmx/dz)
//
//Then if we have an arbitrary symmetry axis unit direction as d = alpha*x + beta*y + gamma*z
//we have the DMI field:

//Hdm,ex = alpha * Hdm,ex,x + beta * Hdm,ex,y + gamma * Hdm,ex,z

viDMExchange::viDMExchange(Mesh *pMesh_) :
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

viDMExchange::~viDMExchange()
{
}

BError viDMExchange::Initialize(void)
{
	BError error(CLASS_STR(viDMExchange));

	error = ExchangeBase::Initialize();

	//Make sure display data has memory allocated (or freed) as required
	error = Update_Module_Display_VECs(
		pMesh->h, pMesh->meshRect,
		(MOD_)pMesh->Get_ActualModule_Heff_Display() == MOD_VIDMEXCHANGE || pMesh->IsOutputDataSet_withRect(DATA_E_EXCH),
		(MOD_)pMesh->Get_ActualModule_Energy_Display() == MOD_VIDMEXCHANGE || pMesh->IsOutputDataSet_withRect(DATA_E_EXCH),
		pMesh->GetMeshType() == MESH_ANTIFERROMAGNETIC);
	if (!error) initialized = true;

	return error;
}

BError viDMExchange::UpdateConfiguration(UPDATECONFIG_ cfgMessage)
{
	BError error(CLASS_STR(viDMExchange));

	Uninitialize();

	//------------------------ CUDA UpdateConfiguration if set

#if COMPILECUDA == 1
	if (pModuleCUDA) {

		if (!error) error = pModuleCUDA->UpdateConfiguration(cfgMessage);
	}
#endif

	return error;
}

BError viDMExchange::MakeCUDAModule(void)
{
	BError error(CLASS_STR(viDMExchange));

#if COMPILECUDA == 1

	if (pMesh->pMeshCUDA) {

		//Note : it is posible pMeshCUDA has not been allocated yet, but this module has been created whilst cuda is switched on. This will happen when a new mesh is being made which adds this module by default.
		//In this case, after the mesh has been fully made, it will call SwitchCUDAState on the mesh, which in turn will call this SwitchCUDAState method; then pMeshCUDA will not be nullptr and we can make the cuda module version
		pModuleCUDA = new viDMExchangeCUDA(pMesh->pMeshCUDA, this);
		error = pModuleCUDA->Error_On_Create();
	}

#endif

	return error;
}

double viDMExchange::UpdateField(void)
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
				DBL3 D_dir = pMesh->D_dir;
				pMesh->update_parameters_mcoarse(idx, pMesh->A, A, pMesh->D, D, pMesh->Ms, Ms, pMesh->D_dir, D_dir);

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

					DBL3 hexch_D_x = DBL3(-Mdiff.y.y - Mdiff.z.z, Mdiff.y.x, Mdiff.z.x);
					DBL3 hexch_D_y = DBL3(Mdiff.x.y, -Mdiff.x.x - Mdiff.z.z, Mdiff.z.y);
					DBL3 hexch_D_z = DBL3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);

					Hexch_D = Dconst * (D_dir.x * hexch_D_x + D_dir.y * hexch_D_y + D_dir.z * hexch_D_z);
				}
				else {

					//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
					DBL3 bnd_dm_dy_x = (D / (2 * A)) * DBL3(-pMesh->M[idx].y, pMesh->M[idx].x, 0);
					DBL3 bnd_dm_dz_x = (D / (2 * A)) * DBL3(-pMesh->M[idx].z, 0, pMesh->M[idx].x);
					DBL33 bnd_nneu_x = DBL33(DBL3(), bnd_dm_dy_x, bnd_dm_dz_x);

					DBL3 bnd_dm_dx_y = (D / (2 * A)) * DBL3(pMesh->M[idx].y, -pMesh->M[idx].x, 0);
					DBL3 bnd_dm_dz_y = (D / (2 * A)) * DBL3(0, -pMesh->M[idx].z, pMesh->M[idx].y);
					DBL33 bnd_nneu_y = DBL33(bnd_dm_dx_y, DBL3(), bnd_dm_dz_y);

					DBL3 bnd_dm_dx_z = (D / (2 * A)) * DBL3(pMesh->M[idx].z, 0, -pMesh->M[idx].x);
					DBL3 bnd_dm_dy_z = (D / (2 * A)) * DBL3(0, pMesh->M[idx].z, -pMesh->M[idx].y);
					DBL33 bnd_nneu_z = DBL33(bnd_dm_dx_z, bnd_dm_dy_z, DBL3());

					DBL33 bnd_nneu = D_dir.x * bnd_nneu_x + D_dir.y * bnd_nneu_y + D_dir.z * bnd_nneu_z;

					//direct exchange contribution
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
					DBL33 Mdiff = pMesh->M.grad_nneu(idx, bnd_nneu);

					DBL3 hexch_D_x = DBL3(-Mdiff.y.y - Mdiff.z.z, Mdiff.y.x, Mdiff.z.x);
					DBL3 hexch_D_y = DBL3(Mdiff.x.y, -Mdiff.x.x - Mdiff.z.z, Mdiff.z.y);
					DBL3 hexch_D_z = DBL3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);

					Hexch_D = Dconst * (D_dir.x * hexch_D_x + D_dir.y * hexch_D_y + D_dir.z * hexch_D_z);
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
				DBL3 D_dir = pMesh->D_dir;
				pMesh->update_parameters_mcoarse(idx, pMesh->A_AFM, A_AFM, pMesh->Ah, Ah, pMesh->Anh, Anh, pMesh->D_AFM, D_AFM, pMesh->Ms_AFM, Ms_AFM, pMesh->D_dir, D_dir);

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

					DBL3 hexch_D_A_x = DBL3(-Mdiff_A.y.y - Mdiff_A.z.z, Mdiff_A.y.x, Mdiff_A.z.x);
					DBL3 hexch_D_A_y = DBL3(Mdiff_A.x.y, -Mdiff_A.x.x - Mdiff_A.z.z, Mdiff_A.z.y);
					DBL3 hexch_D_A_z = DBL3(Mdiff_A.x.z, Mdiff_A.y.z, -Mdiff_A.x.x - Mdiff_A.y.y);

					Hexch_D = Dconst.i * (D_dir.x * hexch_D_A_x + D_dir.y * hexch_D_A_y + D_dir.z * hexch_D_A_z);

					DBL3 hexch_D_B_x = DBL3(-Mdiff_B.y.y - Mdiff_B.z.z, Mdiff_B.y.x, Mdiff_B.z.x);
					DBL3 hexch_D_B_y = DBL3(Mdiff_B.x.y, -Mdiff_B.x.x - Mdiff_B.z.z, Mdiff_B.z.y);
					DBL3 hexch_D_B_z = DBL3(Mdiff_B.x.z, Mdiff_B.y.z, -Mdiff_B.x.x - Mdiff_B.y.y);

					Hexch_D2 = Dconst.j * (D_dir.x * hexch_D_B_x + D_dir.y * hexch_D_B_y + D_dir.z * hexch_D_B_z);
				}
				else {

					DBL33 bndA_nneu, bndB_nneu;

					DBL2 nhconst = Anh / (2 * A_AFM);

					if (fabs(nhconst.i) != 1.0) {

						//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
						DBL3 bnd_dm_dy_x = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(-pMesh->M[idx].y + nhconst.i * pMesh->M2[idx].y, pMesh->M[idx].x - nhconst.i * pMesh->M2[idx].x, 0);
						DBL3 bnd_dm_dz_x = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(-pMesh->M[idx].z + nhconst.i * pMesh->M2[idx].z, 0, pMesh->M[idx].x - nhconst.i * pMesh->M2[idx].x);
						DBL33 bndA_nneu_x = DBL33(DBL3(), bnd_dm_dy_x, bnd_dm_dz_x);

						DBL3 bnd_dm_dx_y = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(pMesh->M[idx].y - nhconst.i * pMesh->M2[idx].y, -pMesh->M[idx].x + nhconst.i * pMesh->M2[idx].x, 0);
						DBL3 bnd_dm_dz_y = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(0, -pMesh->M[idx].z + nhconst.i * pMesh->M2[idx].z, pMesh->M[idx].y - nhconst.i * pMesh->M2[idx].y);
						DBL33 bndA_nneu_y = DBL33(bnd_dm_dx_y, DBL3(), bnd_dm_dz_y);

						DBL3 bnd_dm_dx_z = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(pMesh->M[idx].z - nhconst.i * pMesh->M2[idx].z, 0, -pMesh->M[idx].x + nhconst.i * pMesh->M2[idx].x);
						DBL3 bnd_dm_dy_z = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(0, pMesh->M[idx].z - nhconst.i * pMesh->M2[idx].z, -pMesh->M[idx].y + nhconst.i * pMesh->M2[idx].y);
						DBL33 bndA_nneu_z = DBL33(bnd_dm_dx_z, bnd_dm_dy_z, DBL3());

						bndA_nneu = D_dir.x * bndA_nneu_x + D_dir.y * bndA_nneu_y + D_dir.z * bndA_nneu_z;
					}
					else {

						DBL3 bnd_dm_dy_x = (D_AFM.i / (4 * A_AFM.i)) * DBL3(-pMesh->M[idx].y, pMesh->M[idx].x, 0);
						DBL3 bnd_dm_dz_x = (D_AFM.i / (4 * A_AFM.i)) * DBL3(-pMesh->M[idx].z, 0, pMesh->M[idx].x);
						DBL33 bndA_nneu_x = DBL33(DBL3(), bnd_dm_dy_x, bnd_dm_dz_x);

						DBL3 bnd_dm_dx_y = (D_AFM.i / (4 * A_AFM.i)) * DBL3(pMesh->M[idx].y, -pMesh->M[idx].x, 0);
						DBL3 bnd_dm_dz_y = (D_AFM.i / (4 * A_AFM.i)) * DBL3(0, -pMesh->M[idx].z, pMesh->M[idx].y);
						DBL33 bndA_nneu_y = DBL33(bnd_dm_dx_y, DBL3(), bnd_dm_dz_y);

						DBL3 bnd_dm_dx_z = (D_AFM.i / (4 * A_AFM.i)) * DBL3(pMesh->M[idx].z, 0, -pMesh->M[idx].x);
						DBL3 bnd_dm_dy_z = (D_AFM.i / (4 * A_AFM.i)) * DBL3(0, pMesh->M[idx].z, -pMesh->M[idx].y);
						DBL33 bndA_nneu_z = DBL33(bnd_dm_dx_z, bnd_dm_dy_z, DBL3());

						bndA_nneu = D_dir.x * bndA_nneu_x + D_dir.y * bndA_nneu_y + D_dir.z * bndA_nneu_z;
					}

					if (fabs(nhconst.j) != 1.0) {

						DBL3 bnd_dm_dy_x = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(-pMesh->M2[idx].y + nhconst.j * pMesh->M[idx].y, pMesh->M2[idx].x - nhconst.j * pMesh->M[idx].x, 0);
						DBL3 bnd_dm_dz_x = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(-pMesh->M2[idx].z + nhconst.j * pMesh->M[idx].z, 0, pMesh->M2[idx].x - nhconst.j * pMesh->M[idx].x);
						DBL33 bndB_nneu_x = DBL33(DBL3(), bnd_dm_dy_x, bnd_dm_dz_x);

						DBL3 bnd_dm_dx_y = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(pMesh->M2[idx].y - nhconst.j * pMesh->M[idx].y, -pMesh->M2[idx].x + nhconst.j * pMesh->M[idx].x, 0);
						DBL3 bnd_dm_dz_y = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(0, -pMesh->M2[idx].z + nhconst.j * pMesh->M[idx].z, pMesh->M2[idx].y - nhconst.j * pMesh->M[idx].y);
						DBL33 bndB_nneu_y = DBL33(bnd_dm_dx_y, DBL3(), bnd_dm_dz_y);

						DBL3 bnd_dm_dx_z = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(pMesh->M2[idx].z - nhconst.j * pMesh->M[idx].z, 0, -pMesh->M2[idx].x + nhconst.j * pMesh->M[idx].x);
						DBL3 bnd_dm_dy_z = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(0, pMesh->M2[idx].z - nhconst.j * pMesh->M[idx].z, -pMesh->M2[idx].y + nhconst.j * pMesh->M[idx].y);
						DBL33 bndB_nneu_z = DBL33(bnd_dm_dx_z, bnd_dm_dy_z, DBL3());

						bndB_nneu = D_dir.x * bndB_nneu_x + D_dir.y * bndB_nneu_y + D_dir.z * bndB_nneu_z;
					}
					else {

						DBL3 bnd_dm_dy_x = (D_AFM.j / (4 * A_AFM.j)) * DBL3(-pMesh->M2[idx].y, pMesh->M2[idx].x, 0);
						DBL3 bnd_dm_dz_x = (D_AFM.j / (4 * A_AFM.j)) * DBL3(-pMesh->M2[idx].z, 0, pMesh->M2[idx].x);
						DBL33 bndB_nneu_x = DBL33(DBL3(), bnd_dm_dy_x, bnd_dm_dz_x);

						DBL3 bnd_dm_dx_y = (D_AFM.j / (4 * A_AFM.j)) * DBL3(pMesh->M2[idx].y, -pMesh->M2[idx].x, 0);
						DBL3 bnd_dm_dz_y = (D_AFM.j / (4 * A_AFM.j)) * DBL3(0, -pMesh->M2[idx].z, pMesh->M2[idx].y);
						DBL33 bndB_nneu_y = DBL33(bnd_dm_dx_y, DBL3(), bnd_dm_dz_y);

						DBL3 bnd_dm_dx_z = (D_AFM.j / (4 * A_AFM.j)) * DBL3(pMesh->M2[idx].z, 0, -pMesh->M2[idx].x);
						DBL3 bnd_dm_dy_z = (D_AFM.j / (4 * A_AFM.j)) * DBL3(0, pMesh->M2[idx].z, -pMesh->M2[idx].y);
						DBL33 bndB_nneu_z = DBL33(bnd_dm_dx_z, bnd_dm_dy_z, DBL3());

						bndB_nneu = D_dir.x * bndB_nneu_x + D_dir.y * bndB_nneu_y + D_dir.z * bndB_nneu_z;
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

					DBL3 hexch_D_A_x = DBL3(-Mdiff_A.y.y - Mdiff_A.z.z, Mdiff_A.y.x, Mdiff_A.z.x);
					DBL3 hexch_D_A_y = DBL3(Mdiff_A.x.y, -Mdiff_A.x.x - Mdiff_A.z.z, Mdiff_A.z.y);
					DBL3 hexch_D_A_z = DBL3(Mdiff_A.x.z, Mdiff_A.y.z, -Mdiff_A.x.x - Mdiff_A.y.y);

					Hexch_D = Dconst.i * (D_dir.x * hexch_D_A_x + D_dir.y * hexch_D_A_y + D_dir.z * hexch_D_A_z);

					DBL3 hexch_D_B_x = DBL3(-Mdiff_B.y.y - Mdiff_B.z.z, Mdiff_B.y.x, Mdiff_B.z.x);
					DBL3 hexch_D_B_y = DBL3(Mdiff_B.x.y, -Mdiff_B.x.x - Mdiff_B.z.z, Mdiff_B.z.y);
					DBL3 hexch_D_B_z = DBL3(Mdiff_B.x.z, Mdiff_B.y.z, -Mdiff_B.x.x - Mdiff_B.y.y);

					Hexch_D2 = Dconst.j * (D_dir.x * hexch_D_B_x + D_dir.y * hexch_D_B_y + D_dir.z * hexch_D_B_z);
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

	//Not implemented for this module

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
double viDMExchange::Get_EnergyChange(int spin_index, DBL3 Mnew)
{
	//For CUDA there are separate device functions used by CUDA kernels.

	if (pMesh->M.is_not_empty(spin_index)) {

		double Ms = pMesh->Ms;
		double A = pMesh->A;
		double D = pMesh->D; 
		DBL3 D_dir = pMesh->D_dir;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->A, A, pMesh->D, D, pMesh->Ms, Ms, pMesh->D_dir, D_dir);

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

				DBL3 hexch_D_x = DBL3(-Mdiff.y.y - Mdiff.z.z, Mdiff.y.x, Mdiff.z.x);
				DBL3 hexch_D_y = DBL3(Mdiff.x.y, -Mdiff.x.x - Mdiff.z.z, Mdiff.z.y);
				DBL3 hexch_D_z = DBL3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);

				Hexch_D = Dconst * (D_dir.x * hexch_D_x + D_dir.y * hexch_D_y + D_dir.z * hexch_D_z);
			}
			else {

				//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
				DBL3 bnd_dm_dy_x = (D / (2 * A)) * DBL3(-pMesh->M[spin_index].y, pMesh->M[spin_index].x, 0);
				DBL3 bnd_dm_dz_x = (D / (2 * A)) * DBL3(-pMesh->M[spin_index].z, 0, pMesh->M[spin_index].x);
				DBL33 bnd_nneu_x = DBL33(DBL3(), bnd_dm_dy_x, bnd_dm_dz_x);

				DBL3 bnd_dm_dx_y = (D / (2 * A)) * DBL3(pMesh->M[spin_index].y, -pMesh->M[spin_index].x, 0);
				DBL3 bnd_dm_dz_y = (D / (2 * A)) * DBL3(0, -pMesh->M[spin_index].z, pMesh->M[spin_index].y);
				DBL33 bnd_nneu_y = DBL33(bnd_dm_dx_y, DBL3(), bnd_dm_dz_y);

				DBL3 bnd_dm_dx_z = (D / (2 * A)) * DBL3(pMesh->M[spin_index].z, 0, -pMesh->M[spin_index].x);
				DBL3 bnd_dm_dy_z = (D / (2 * A)) * DBL3(0, pMesh->M[spin_index].z, -pMesh->M[spin_index].y);
				DBL33 bnd_nneu_z = DBL33(bnd_dm_dx_z, bnd_dm_dy_z, DBL3());

				DBL33 bnd_nneu = D_dir.x * bnd_nneu_x + D_dir.y * bnd_nneu_y + D_dir.z * bnd_nneu_z;

				//direct exchange contribution
				//cells marked with cmbnd are calculated using exchange coupling to other ferromagnetic meshes - see below; the delsq_nneu evaluates to zero in the CMBND coupling direction.
				Hexch_A = Aconst * pMesh->M.delsq_nneu(spin_index, bnd_nneu);

				//Differentials of M components (we only need 4, not all 9 so this could be optimised). First index is the differential direction, second index is the M component
				DBL33 Mdiff = pMesh->M.grad_nneu(spin_index, bnd_nneu);

				DBL3 hexch_D_x = DBL3(-Mdiff.y.y - Mdiff.z.z, Mdiff.y.x, Mdiff.z.x);
				DBL3 hexch_D_y = DBL3(Mdiff.x.y, -Mdiff.x.x - Mdiff.z.z, Mdiff.z.y);
				DBL3 hexch_D_z = DBL3(Mdiff.x.z, Mdiff.y.z, -Mdiff.x.x - Mdiff.y.y);

				Hexch_D = Dconst * (D_dir.x * hexch_D_x + D_dir.y * hexch_D_y + D_dir.z * hexch_D_z);
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
DBL2 viDMExchange::Get_EnergyChange(int spin_index, DBL3 Mnew_A, DBL3 Mnew_B)
{
	if (pMesh->M.is_not_empty(spin_index) && pMesh->M2.is_not_empty(spin_index)) {

		DBL2 Ms_AFM = pMesh->Ms_AFM;
		DBL2 A_AFM = pMesh->A_AFM;
		DBL2 Ah = pMesh->Ah;
		DBL2 Anh = pMesh->Anh;
		DBL2 D_AFM = pMesh->D_AFM;
		DBL3 D_dir = pMesh->D_dir;
		pMesh->update_parameters_mcoarse(spin_index, pMesh->A_AFM, A_AFM, pMesh->Ah, Ah, pMesh->Anh, Anh, pMesh->Ms_AFM, Ms_AFM, pMesh->D_AFM, D_AFM, pMesh->D_dir, D_dir);

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

				DBL3 hexch_D_A_x = DBL3(-Mdiff_A.y.y - Mdiff_A.z.z, Mdiff_A.y.x, Mdiff_A.z.x);
				DBL3 hexch_D_A_y = DBL3(Mdiff_A.x.y, -Mdiff_A.x.x - Mdiff_A.z.z, Mdiff_A.z.y);
				DBL3 hexch_D_A_z = DBL3(Mdiff_A.x.z, Mdiff_A.y.z, -Mdiff_A.x.x - Mdiff_A.y.y);

				Hexch_D = Dconst.i * (D_dir.x * hexch_D_A_x + D_dir.y * hexch_D_A_y + D_dir.z * hexch_D_A_z);

				DBL3 hexch_D_B_x = DBL3(-Mdiff_B.y.y - Mdiff_B.z.z, Mdiff_B.y.x, Mdiff_B.z.x);
				DBL3 hexch_D_B_y = DBL3(Mdiff_B.x.y, -Mdiff_B.x.x - Mdiff_B.z.z, Mdiff_B.z.y);
				DBL3 hexch_D_B_z = DBL3(Mdiff_B.x.z, Mdiff_B.y.z, -Mdiff_B.x.x - Mdiff_B.y.y);

				Hexch_D2 = Dconst.j * (D_dir.x * hexch_D_B_x + D_dir.y * hexch_D_B_y + D_dir.z * hexch_D_B_z);
			}
			else {

				DBL33 bndA_nneu, bndB_nneu;

				DBL2 nhconst = Anh / (2 * A_AFM);

				if (fabs(nhconst.i) != 1.0) {

					//Non-homogeneous Neumann boundary conditions apply when using DMI. Required to ensure Brown's condition is fulfilled, i.e. m x h -> 0 when relaxing.
					DBL3 bnd_dm_dy_x = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(-pMesh->M[spin_index].y + nhconst.i * pMesh->M2[spin_index].y, pMesh->M[spin_index].x - nhconst.i * pMesh->M2[spin_index].x, 0);
					DBL3 bnd_dm_dz_x = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(-pMesh->M[spin_index].z + nhconst.i * pMesh->M2[spin_index].z, 0, pMesh->M[spin_index].x - nhconst.i * pMesh->M2[spin_index].x);
					DBL33 bndA_nneu_x = DBL33(DBL3(), bnd_dm_dy_x, bnd_dm_dz_x);

					DBL3 bnd_dm_dx_y = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(pMesh->M[spin_index].y - nhconst.i * pMesh->M2[spin_index].y, -pMesh->M[spin_index].x + nhconst.i * pMesh->M2[spin_index].x, 0);
					DBL3 bnd_dm_dz_y = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(0, -pMesh->M[spin_index].z + nhconst.i * pMesh->M2[spin_index].z, pMesh->M[spin_index].y - nhconst.i * pMesh->M2[spin_index].y);
					DBL33 bndA_nneu_y = DBL33(bnd_dm_dx_y, DBL3(), bnd_dm_dz_y);

					DBL3 bnd_dm_dx_z = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(pMesh->M[spin_index].z - nhconst.i * pMesh->M2[spin_index].z, 0, -pMesh->M[spin_index].x + nhconst.i * pMesh->M2[spin_index].x);
					DBL3 bnd_dm_dy_z = (D_AFM.i / (2 * A_AFM.i * (1 - nhconst.i * nhconst.i))) * DBL3(0, pMesh->M[spin_index].z - nhconst.i * pMesh->M2[spin_index].z, -pMesh->M[spin_index].y + nhconst.i * pMesh->M2[spin_index].y);
					DBL33 bndA_nneu_z = DBL33(bnd_dm_dx_z, bnd_dm_dy_z, DBL3());

					bndA_nneu = D_dir.x * bndA_nneu_x + D_dir.y * bndA_nneu_y + D_dir.z * bndA_nneu_z;
				}
				else {

					DBL3 bnd_dm_dy_x = (D_AFM.i / (4 * A_AFM.i)) * DBL3(-pMesh->M[spin_index].y, pMesh->M[spin_index].x, 0);
					DBL3 bnd_dm_dz_x = (D_AFM.i / (4 * A_AFM.i)) * DBL3(-pMesh->M[spin_index].z, 0, pMesh->M[spin_index].x);
					DBL33 bndA_nneu_x = DBL33(DBL3(), bnd_dm_dy_x, bnd_dm_dz_x);

					DBL3 bnd_dm_dx_y = (D_AFM.i / (4 * A_AFM.i)) * DBL3(pMesh->M[spin_index].y, -pMesh->M[spin_index].x, 0);
					DBL3 bnd_dm_dz_y = (D_AFM.i / (4 * A_AFM.i)) * DBL3(0, -pMesh->M[spin_index].z, pMesh->M[spin_index].y);
					DBL33 bndA_nneu_y = DBL33(bnd_dm_dx_y, DBL3(), bnd_dm_dz_y);

					DBL3 bnd_dm_dx_z = (D_AFM.i / (4 * A_AFM.i)) * DBL3(pMesh->M[spin_index].z, 0, -pMesh->M[spin_index].x);
					DBL3 bnd_dm_dy_z = (D_AFM.i / (4 * A_AFM.i)) * DBL3(0, pMesh->M[spin_index].z, -pMesh->M[spin_index].y);
					DBL33 bndA_nneu_z = DBL33(bnd_dm_dx_z, bnd_dm_dy_z, DBL3());

					bndA_nneu = D_dir.x * bndA_nneu_x + D_dir.y * bndA_nneu_y + D_dir.z * bndA_nneu_z;
				}

				if (fabs(nhconst.j) != 1.0) {

					DBL3 bnd_dm_dy_x = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(-pMesh->M2[spin_index].y + nhconst.j * pMesh->M[spin_index].y, pMesh->M2[spin_index].x - nhconst.j * pMesh->M[spin_index].x, 0);
					DBL3 bnd_dm_dz_x = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(-pMesh->M2[spin_index].z + nhconst.j * pMesh->M[spin_index].z, 0, pMesh->M2[spin_index].x - nhconst.j * pMesh->M[spin_index].x);
					DBL33 bndB_nneu_x = DBL33(DBL3(), bnd_dm_dy_x, bnd_dm_dz_x);

					DBL3 bnd_dm_dx_y = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(pMesh->M2[spin_index].y - nhconst.j * pMesh->M[spin_index].y, -pMesh->M2[spin_index].x + nhconst.j * pMesh->M[spin_index].x, 0);
					DBL3 bnd_dm_dz_y = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(0, -pMesh->M2[spin_index].z + nhconst.j * pMesh->M[spin_index].z, pMesh->M2[spin_index].y - nhconst.j * pMesh->M[spin_index].y);
					DBL33 bndB_nneu_y = DBL33(bnd_dm_dx_y, DBL3(), bnd_dm_dz_y);

					DBL3 bnd_dm_dx_z = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(pMesh->M2[spin_index].z - nhconst.j * pMesh->M[spin_index].z, 0, -pMesh->M2[spin_index].x + nhconst.j * pMesh->M[spin_index].x);
					DBL3 bnd_dm_dy_z = (D_AFM.j / (2 * A_AFM.j * (1 - nhconst.j * nhconst.j))) * DBL3(0, pMesh->M2[spin_index].z - nhconst.j * pMesh->M[spin_index].z, -pMesh->M2[spin_index].y + nhconst.j * pMesh->M[spin_index].y);
					DBL33 bndB_nneu_z = DBL33(bnd_dm_dx_z, bnd_dm_dy_z, DBL3());

					bndB_nneu = D_dir.x * bndB_nneu_x + D_dir.y * bndB_nneu_y + D_dir.z * bndB_nneu_z;
				}
				else {

					DBL3 bnd_dm_dy_x = (D_AFM.j / (4 * A_AFM.j)) * DBL3(-pMesh->M2[spin_index].y, pMesh->M2[spin_index].x, 0);
					DBL3 bnd_dm_dz_x = (D_AFM.j / (4 * A_AFM.j)) * DBL3(-pMesh->M2[spin_index].z, 0, pMesh->M2[spin_index].x);
					DBL33 bndB_nneu_x = DBL33(DBL3(), bnd_dm_dy_x, bnd_dm_dz_x);

					DBL3 bnd_dm_dx_y = (D_AFM.j / (4 * A_AFM.j)) * DBL3(pMesh->M2[spin_index].y, -pMesh->M2[spin_index].x, 0);
					DBL3 bnd_dm_dz_y = (D_AFM.j / (4 * A_AFM.j)) * DBL3(0, -pMesh->M2[spin_index].z, pMesh->M2[spin_index].y);
					DBL33 bndB_nneu_y = DBL33(bnd_dm_dx_y, DBL3(), bnd_dm_dz_y);

					DBL3 bnd_dm_dx_z = (D_AFM.j / (4 * A_AFM.j)) * DBL3(pMesh->M2[spin_index].z, 0, -pMesh->M2[spin_index].x);
					DBL3 bnd_dm_dy_z = (D_AFM.j / (4 * A_AFM.j)) * DBL3(0, pMesh->M2[spin_index].z, -pMesh->M2[spin_index].y);
					DBL33 bndB_nneu_z = DBL33(bnd_dm_dx_z, bnd_dm_dy_z, DBL3());

					bndB_nneu = D_dir.x * bndB_nneu_x + D_dir.y * bndB_nneu_y + D_dir.z * bndB_nneu_z;
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

				DBL3 hexch_D_A_x = DBL3(-Mdiff_A.y.y - Mdiff_A.z.z, Mdiff_A.y.x, Mdiff_A.z.x);
				DBL3 hexch_D_A_y = DBL3(Mdiff_A.x.y, -Mdiff_A.x.x - Mdiff_A.z.z, Mdiff_A.z.y);
				DBL3 hexch_D_A_z = DBL3(Mdiff_A.x.z, Mdiff_A.y.z, -Mdiff_A.x.x - Mdiff_A.y.y);

				Hexch_D = Dconst.i * (D_dir.x * hexch_D_A_x + D_dir.y * hexch_D_A_y + D_dir.z * hexch_D_A_z);

				DBL3 hexch_D_B_x = DBL3(-Mdiff_B.y.y - Mdiff_B.z.z, Mdiff_B.y.x, Mdiff_B.z.x);
				DBL3 hexch_D_B_y = DBL3(Mdiff_B.x.y, -Mdiff_B.x.x - Mdiff_B.z.z, Mdiff_B.z.y);
				DBL3 hexch_D_B_z = DBL3(Mdiff_B.x.z, Mdiff_B.y.z, -Mdiff_B.x.x - Mdiff_B.y.y);

				Hexch_D2 = Dconst.j * (D_dir.x * hexch_D_B_x + D_dir.y * hexch_D_B_y + D_dir.z * hexch_D_B_z);
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

DBL3 viDMExchange::GetTorque(Rect& avRect)
{
#if COMPILECUDA == 1
	if (pModuleCUDA) return pModuleCUDA->GetTorque(avRect);
#endif

	return CalculateTorque(pMesh->M, avRect);
}

#endif