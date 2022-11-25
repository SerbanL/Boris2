#include "TMRCUDA.h"

#if COMPILECUDA == 1

#ifdef MODULE_COMPILATION_TMR

#include "BorisCUDALib.cuh"

#include "MeshCUDA.h"
#include "SuperMeshCUDA.h"
#include "MeshParamsControlCUDA.h"

#include "ManagedDiffEqFMCUDA.h"

//-------------------Calculation Methods : Iterate Spin-Charge Solver

void TMRCUDA::IterateSpinSolver_Charge_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value, bool use_NNeu)
{
	pMeshCUDA->V()->IteratePoisson_SOR(pMeshCUDA->n_e.dim(), (TransportCUDA_Spin_V_Funcs&)poisson_Spin_V, damping, max_error, max_value);
}

//------------------- PRIME SPIN-CHARGE SOLVER

//before iterating the spin solver (charge part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
void TMRCUDA::PrimeSpinSolver_Charge(void)
{
}

//-------------------Calculation Methods : Iterate Spin-Spin Solver

//solve for spin accumulation using Poisson equation for delsq_S, solved using SOR algorithm
void TMRCUDA::IterateSpinSolver_Spin_SOR(cu_obj<cuBReal>& damping, cu_obj<cuBReal>& max_error, cu_obj<cuBReal>& max_value, bool use_NNeu)
{
	pMeshCUDA->S()->IteratePoisson_SOR(pMeshCUDA->n_e.dim(), (TransportCUDA_Spin_S_Funcs&)poisson_Spin_S, damping, max_error, max_value);
}

//------------------- PRIME SPIN-SPIN SOLVER

//before iterating the spin solver (spin part) we need to prime it : pre-compute values which do not change as the spin solver relaxes.
void TMRCUDA::PrimeSpinSolver_Spin(void)
{
}

//--------------------------------------------------------------- Effective field from spin accumulation

//Spin accumulation field
void TMRCUDA::CalculateSAField(void)
{
}

//--------------------------------------------------------------- Effective field from interface spin accumulation drop

//Calculate the field resulting from interface spin accumulation torque for a given contact (in magnetic meshes for NF interfaces with G interface conductance set)
void TMRCUDA::CalculateSAInterfaceField(TransportBaseCUDA* ptrans_sec, CMBNDInfoCUDA& contactCUDA, bool primary_top)
{
}

#endif

#endif