#pragma once

#include "VEC.h"
#include "VEC_MeshTransfer.h"

//--------------------------------------------MESH TRANSFER : VEC_MeshTransfer.h

//GET
template <typename VType>
Transfer<VType>& VEC<VType>::get_transfer(void)
{
	return transfer;
}

template <typename VType>
Transfer<VType>& VEC<VType>::get_transfer2(void)
{
	return transfer2;
}


//CLEAR
template <typename VType>
void VEC<VType>::clear_transfer(void)
{
	transfer.clear();
}

template <typename VType>
void VEC<VType>::clear_transfer2(void)
{
	transfer2.clear();
}

//SINGLE INPUT, SINGLE OUTPUT

//set-up mesh transfers, ready to use - return false if failed (not enough memory)
//can set a multiplier value to apply to all contributions (1.0 by default, so no effect)
template <typename VType>
bool VEC<VType>::Initialize_MeshTransfer(const std::vector< VEC<VType>* >& mesh_in, const std::vector< VEC<VType>* >& mesh_out, int correction_type, double multiplier)
{ 
	return transfer.initialize_transfer(mesh_in, mesh_out, correction_type, multiplier); 
}

template <typename VType>
bool VEC<VType>::Initialize_MeshTransfer2(const std::vector< VEC<VType>* >& mesh_in, const std::vector< VEC<VType>* >& mesh_out, int correction_type, double multiplier)
{
	return transfer2.initialize_transfer(mesh_in, mesh_out, correction_type, multiplier);
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

template <typename VType>
bool VEC<VType>::Initialize_MeshTransfer2_AveragedInputs(const std::vector< VEC<VType>* >& mesh_in, const std::vector< VEC<VType>* >& mesh_in2, const std::vector< VEC<VType>* >& mesh_out, int correction_type, double multiplier)
{
	return transfer2.initialize_transfer_averagedinputs(mesh_in, mesh_in2, mesh_out, correction_type, multiplier);
}

//mesh_in2 must be a vector of VEC<double> inputs
//can set a multiplier value to apply to all contributions (1.0 by default, so no effect)
template <typename VType>
bool VEC<VType>::Initialize_MeshTransfer_MultipliedInputs(const std::vector< VEC<VType>* >& mesh_in, const std::vector< VEC<double>* >& mesh_in2_double, const std::vector< VEC<VType>* >& mesh_out, int correction_type, double multiplier)
{
	return transfer.initialize_transfer_multipliedinputs(mesh_in, mesh_in2_double, mesh_out, correction_type, multiplier);
}

template <typename VType>
bool VEC<VType>::Initialize_MeshTransfer2_MultipliedInputs(const std::vector< VEC<VType>* >& mesh_in, const std::vector< VEC<double>* >& mesh_in2_double, const std::vector< VEC<VType>* >& mesh_out, int correction_type, double multiplier)
{
	return transfer2.initialize_transfer_multipliedinputs(mesh_in, mesh_in2_double, mesh_out, correction_type, multiplier);
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

template <typename VType>
bool VEC<VType>::Initialize_MeshTransfer2_AveragedInputs_DuplicatedOutputs(const std::vector< VEC<VType>* >& mesh_in, const std::vector< VEC<VType>* >& mesh_in2, const std::vector< VEC<VType>* >& mesh_out, const std::vector< VEC<VType>* >& mesh_out2, int correction_type, double multiplier)
{
	return transfer2.initialize_transfer_averagedinputs_duplicatedoutputs(mesh_in, mesh_in2, mesh_out, mesh_out2, correction_type, multiplier);
}

//SINGLE INPUT, SINGLE OUTPUT

//do the actual transfer of values to and from this mesh using these
template <typename VType>
void VEC<VType>::transfer_in(bool clear_input)
{ 
	transfer.transfer_from_external_meshes(clear_input);
}

template <typename VType>
void VEC<VType>::transfer2_in(bool clear_input)
{
	transfer2.transfer_from_external_meshes(clear_input);
}
	
template <typename VType>
void VEC<VType>::transfer_out(bool setOutput) 
{ 
	transfer.transfer_to_external_meshes(setOutput); 
}

template <typename VType>
void VEC<VType>::transfer2_out(bool setOutput)
{
	transfer2.transfer_to_external_meshes(setOutput);
}

//AVERAGED INPUTS

template <typename VType>
void VEC<VType>::transfer_in_averaged(bool clear_input)
{ 
	transfer.transfer_from_external_meshes_averaged(clear_input);
}

template <typename VType>
void VEC<VType>::transfer2_in_averaged(bool clear_input)
{
	transfer2.transfer_from_external_meshes_averaged(clear_input);
}

//MULTIPLIED INPUTS

//only use this if mesh_in2 was initialized as a vector of VEC<double> inputs : see corresponding Initialize_MeshTransfer_MultipliedInputs
template <typename VType>
void VEC<VType>::transfer_in_multiplied(bool clear_input)
{ 
	transfer.transfer_from_external_meshes_multiplied(clear_input);
}

template <typename VType>
void VEC<VType>::transfer2_in_multiplied(bool clear_input)
{
	transfer2.transfer_from_external_meshes_multiplied(clear_input);
}

//DUPLICATED OUTPUTS

template <typename VType>
void VEC<VType>::transfer_out_duplicated(bool setOutput) 
{ 
	transfer.transfer_to_external_meshes_duplicated(setOutput); 
}

template <typename VType>
void VEC<VType>::transfer2_out_duplicated(bool setOutput)
{
	transfer2.transfer_to_external_meshes_duplicated(setOutput);
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

//flattened in and out transfer sizes (i.e. total number of cell contributions
template <typename VType>
size_t VEC<VType>::size_transfer2_in(void)
{
	return transfer2.size_transfer_in();
}

template <typename VType>
size_t VEC<VType>::size_transfer2_out(void)
{
	return transfer2.size_transfer_out();
}