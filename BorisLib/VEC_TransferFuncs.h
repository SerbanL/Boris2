#pragma once

#include "VEC.h"
#include "VEC_MeshTransfer.h"

//--------------------------------------------MESH TRANSFER : VEC_MeshTransfer.h

//SINGLE INPUT, SINGLE OUTPUT

//set-up mesh transfers, ready to use - return false if failed (not enough memory)
//can set a multiplier value to apply to all contributions (1.0 by default, so no effect)
template <typename VType>
bool VEC<VType>::Initialize_MeshTransfer(const std::vector< VEC<VType>* >& mesh_in, const std::vector< VEC<VType>* >& mesh_out, int correction_type, double multiplier)
{ 
	return transfer.initialize_transfer(mesh_in, mesh_out, correction_type, multiplier); 
}

//MULTIPLE INPUTS, SINGLE OUTPUT
	
//mesh_in and mesh_in2 vectors must have same sizes
//All VECs in mesh_in should be non-empty
//Some VECs in mesh_in2 allowed to be non-empty (in this case single input is used), but otherwise should have exactly same dimensions as the corresponding VECs in mesh_in
//can set a multiplier value to apply to all contributions (1.0 by default, so no effect)
template <typename VType>
bool VEC<VType>::Initialize_MeshTransfer_AveragedInputs(const std::vector< VEC<VType>* >& mesh_in, const std::vector< VEC<VType>* >& mesh_in2, const std::vector< VEC<VType>* >& mesh_out, int correction_type, double multiplier)
{ 
	return transfer.initialize_transfer_averagedinputs(mesh_in, mesh_in2, mesh_out, correction_type, multiplier);
}

//mesh_in2 must be a vector of VEC<double> inputs
//can set a multiplier value to apply to all contributions (1.0 by default, so no effect)
template <typename VType>
bool VEC<VType>::Initialize_MeshTransfer_MultipliedInputs(const std::vector< VEC<VType>* >& mesh_in, const std::vector< VEC<double>* >& mesh_in2_double, const std::vector< VEC<VType>* >& mesh_out, int correction_type, double multiplier)
{
	return transfer.initialize_transfer_multipliedinputs(mesh_in, mesh_in2_double, mesh_out, correction_type, multiplier);
}

//MULTIPLE INPUT, MULTIPLE OUTPUT

//mesh_in and mesh_in2 vectors must have same sizes; same as mesh_out, mesh_out2
//All VECs in mesh_in and mesh_out should be non-empty
//Some VECs in mesh_in2 and mesh_out2 allowed to be non-empty (in this case single input/output is used), but otherwise should have exactly same dimensions as the corresponding VECs in mesh_in, mesh_out
//Also if a VEC in mesh_in2 is non-empty the corresponding VEC in mesh_out2 should also be non-empty.
//can set a multiplier value to apply to all contributions (1.0 by default, so no effect)
template <typename VType>
bool VEC<VType>::Initialize_MeshTransfer_AveragedInputs_DuplicatedOutputs(const std::vector< VEC<VType>* >& mesh_in, const std::vector< VEC<VType>* >& mesh_in2, const std::vector< VEC<VType>* >& mesh_out, const std::vector< VEC<VType>* >& mesh_out2, int correction_type, double multiplier)
{ 
	return transfer.initialize_transfer_averagedinputs_duplicatedoutputs(mesh_in, mesh_in2, mesh_out, mesh_out2, correction_type, multiplier);
}

//SINGLE INPUT, SINGLE OUTPUT

//do the actual transfer of values to and from this mesh using these
template <typename VType>
void VEC<VType>::transfer_in(void) 
{ 
	transfer.transfer_from_external_meshes(); 
}
	
template <typename VType>
void VEC<VType>::transfer_out(bool setOutput) 
{ 
	transfer.transfer_to_external_meshes(setOutput); 
}

//AVERAGED INPUTS

template <typename VType>
void VEC<VType>::transfer_in_averaged(void) 
{ 
	transfer.transfer_from_external_meshes_averaged(); 
}

//MULTIPLIED INPUTS

//only use this if mesh_in2 was initialized as a vector of VEC<double> inputs : see corresponding Initialize_MeshTransfer_MultipliedInputs
template <typename VType>
void VEC<VType>::transfer_in_multiplied(void) 
{ 
	transfer.transfer_from_external_meshes_multiplied(); 
}

//DUPLICATED OUTPUTS

template <typename VType>
void VEC<VType>::transfer_out_duplicated(bool setOutput) 
{ 
	transfer.transfer_to_external_meshes_duplicated(setOutput); 
}

//flattened in and out transfer sizes (i.e. total number of cell contributions
template <typename VType>
size_t VEC<VType>::size_transfer_in(void)
{ 
	return transfer.size_transfer_in(); 
}

template <typename VType>
size_t VEC<VType>::size_transfer_out(void)
{ 
	return transfer.size_transfer_out(); 
}

//this is used to pass transfer information to a cuVEC for copying to gpu memory : for gpu computations we use "flattened" transfers so it can be parallelized better
//return type: vector of transfers, where INT3 contains : i - input mesh index, j - input mesh cell index, k - super-mesh cell index. the double entry is the weight for the value contribution
template <typename VType>
std::vector<std::pair<INT3, double>> VEC<VType>::get_flattened_transfer_in_info(void) 
{ 
	return transfer.get_flattened_transfer_in_info(); 
}

//this is used to pass transfer information to a cuVEC for copying to gpu memory : for gpu computations we use "flattened" transfers so it can be parallelized better
//return type: vector of transfers, where INT3 contains : i - output mesh index, j - output mesh cell index, k - super-mesh cell index. the double entry is the weight for the value contribution
template <typename VType>
std::vector<std::pair<INT3, double>> VEC<VType>::get_flattened_transfer_out_info(void) 
{ 
	return transfer.get_flattened_transfer_out_info(); 
}