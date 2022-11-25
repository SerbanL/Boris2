#pragma once

#include <vector>
#include <string>

#include "Funcs_Conv.h"

///////////////////////////////////////////////////////////////////////////////
//VECTORS

//----------------------------------------

//joins rhs to lhs std::vector, also converting elements in rhs std::vector to type of elements in lhs std::vector (must be convertible)
template <typename VType, typename PType>
void JoinToVector(std::vector<VType> &lhsVector, const std::vector<PType> &rhsVector)
{
	lhsVector.reserve(lhsVector.size() + rhsVector.size());
	lhsVector.insert(lhsVector.end(), rhsVector.begin(), rhsVector.end());
}

//JointToVector, multiple vectors version
template <typename VType, typename Type, typename ... PType>
void JoinToVector(std::vector<VType> &lhsVector, const std::vector<Type> &rhsVector, const std::vector<PType>&... furtherVectors)
{
	JoinToVector(lhsVector, rhsVector);
	JoinToVector(lhsVector, furtherVectors...);
}

//----------------------------------------

//make a std::vector using parameter pack - the first parameter determines the std::vector stored type
template <typename VType> 
std::vector<VType> make_vector(const VType& value) { return { value }; }

template <typename VType, typename... PType> 
std::vector<VType> make_vector(const VType& value, const PType&... values) { return { value, (VType)values... }; }

template <typename VType> 
std::vector<VType> make_vector(void) { return std::vector<VType>(); }

template <typename... PType> 
std::vector<std::string> make_strings_vector(const PType&... values) { return { ToString(values)... }; }

//----------------------------------------

//check if std::vector contains value
template <typename VType> 
bool vector_contains(const std::vector<VType> &searchVector, const VType& value) 
{	
	if(find(searchVector.begin(), searchVector.end(), value) != searchVector.end()) return true;
	else return false;
}

//check if std::vector contains value and return index of found element (-1 if not found)
template <typename VType> 
int search_vector(const std::vector<VType> &searchVector, const VType& value)
{
	auto it = find(searchVector.begin(), searchVector.end(), value);

	if (it == searchVector.end()) return -1;
	else return (it - searchVector.begin());
}

//----------------------------------------

//return subvector as a new copy. subvector starts from startIdx and is endIdx - startIdx elements long. Leave endIdx out to return the rest of the std::vector.
template <typename VType> 
std::vector<VType> subvec(const std::vector<VType> &vec, int startIdx, int endIdx = -1) 
{
	if(endIdx < 0) {

		std::vector<VType> newVec(vec.begin() + startIdx, vec.end());
		return newVec;
	}
	else {
		
		std::vector<VType> newVec(vec.begin() + startIdx, vec.begin() + endIdx);
		return newVec;
	}
}

//return subvector as a new copy, where vec is a linear mapping of 3D memory with dimensions ni, nj, nk
//extract sub-vector starting at cell coordinates si, sj, sk, up to cell (not including) with coordinates ei, ej, ek.
template <typename VType> 
std::vector<VType> subvec(
	const std::vector<VType> &vec, 
	int si, int sj, int sk, 
	int ei, int ej, int ek, 
	int ni, int nj, int nk)
{
	std::vector<VType> newVec((ei - si)*(ej - sj)*(ek - sk));

	for (int k = sk; k < ek; k++) {
#pragma omp parallel for
		for (int j = sj; j < ej; j++) {
			for (int i = si; i < ei; i++) {
			
				newVec[(i - si) + (j - sj)*(ei - si) + (k - sk)*(ei - si)*(ej - sj)] = vec[i + j*ni + k*ni*nj];
			}
		}
	}

	return newVec;
}

//----------------------------------------

//copy values from vec_from to vec_to using strides
//vec_to is a linear mapping of 3D memory with dimensions ni, nj, nk
//copy starting at coordinates si, sj, sk in vec_to, up to ei, ej, ek in vec_to
//here vec_from must have dimensions (ei - si), (ej - sj), (ek - sk)
template <typename VType>
bool strided_copy(
	std::vector<VType> &vec_to, std::vector<VType> &vec_from,
	int si, int sj, int sk,
	int ei, int ej, int ek,
	int ni, int nj, int nk)
{
	if (vec_from.size() != (ei - si)*(ej - sj)*(ek - sk) || vec_to.size() != ni*nj*nk) return false;

	for (int k = sk; k < ek; k++) {
#pragma omp parallel for
		for (int j = sj; j < ej; j++) {
			for (int i = si; i < ei; i++) {

				vec_to[i + j*ni + k*ni*nj] = vec_from[(i - si) + (j - sj)*(ei - si) + (k - sk)*(ei - si)*(ej - sj)];
			}
		}
	}

	return true;
}

//----------------------------------------

//set value at given index in a collection of vectors (value must be convertible to each of the vectors' types)
template <typename Type, typename VType>
void set_all_vectors(int index, Type value, std::vector<VType>& vec) { vec[index] = (VType)value; }

template <typename Type, typename VType, typename ... PType>
void set_all_vectors(int index, Type value, std::vector<VType>& vec, std::vector<PType>&... further_vec)
{
	vec[index] = (VType)value;
	set_all_vectors(index, value, further_vec...);
}

//----------------------------------------

//swap elements with indexes j and k in given std::vector, or collection of vectors
template <typename VType>
void swap(int j, int k, std::vector<VType>& vec)
{
	VType sv = vec[j];
	vec[j] = vec[k];
	vec[k] = sv;
}

template <typename VType, typename ... PType>
void swap(int j, int k, std::vector<VType>& vec, std::vector<PType>&... further_vec)
{
	swap(j, k, vec);
	swap(j, k, further_vec...);
}

//----------------------------------------

//erase an entry at given index in a vector - single and multiple vectors versions
template <typename VType>
void erase_vectors_entry(int index, std::vector<VType>& vec)
{
	if (index >= 0 && index < vec.size()) vec.erase(vec.begin() + index);
}

template <typename VType, typename ... PType>
void erase_vectors_entry(int index, std::vector<VType>& vec, std::vector<PType>&... further_vec)
{
	if (index >= 0 && index < vec.size()) vec.erase(vec.begin() + index);

	erase_vectors_entry(index, further_vec...);
}

//----------------------------------------

template <typename Vector>
void clear_vector(Vector& vec, std::true_type pointers)
{
	//storing pointers : deep clear
	
	for (auto it = vec.begin(); it != vec.end(); ++it) {

		if (*it != nullptr) delete *it;
		*it = nullptr;
	}

	//free memory, including reserved memory (capacity)
	vec.clear();
	vec.shrink_to_fit();
}

template <typename Vector>
void clear_vector(Vector& vec, std::false_type pointers)
{
	//not storing pointers : just free memory, including reserved memory (capacity)
	vec.clear();
	vec.shrink_to_fit();
}

//deep clear a vector - if it contains pointers then delete them all before clearing the vector
template <typename Vector>
void clear_vector(Vector& vec)
{
	//stored element type
	using SType = typename contained_type<Vector>::type;
	
	clear_vector(vec, std::is_pointer<SType>());
}

//----------------------------------------

//convert input vector to output vector element by element
template <typename OType, typename IType>
std::vector<OType> vec_convert(const std::vector<IType> &in_vec)
{
	std::vector<OType> out_vec;
	out_vec.resize(in_vec.size());

	#pragma omp parallel for
	for (int idx = 0; idx < (int)in_vec.size(); idx++) {

		out_vec[idx] = Convert<OType, IType>(in_vec[idx]);
	}

	return out_vec;
}

//convert input vector to output vector element by element.
template <typename OType, typename IType>
bool vec_convert(const std::vector<IType> &in_vec, std::vector<OType>& out_vec)
{
	if (!malloc_vector(out_vec, in_vec.size())) return false;

	#pragma omp parallel for
	for (int idx = 0; idx < (int)in_vec.size(); idx++) {

		out_vec[idx] = Convert<OType, IType>(in_vec[idx]);
	}

	return true;
}

//----------------------------------------

//delete value repeats when next to each other - single and multiple vectors versions
template <typename VType>
void delete_repeats(std::vector<VType>& vec)
{
	if (!vec.size()) return;

	VType value = vec[0];

	int idx = 1;
	while(idx < (int)vec.size()) {

		if (vec[idx] == value) {

			vec.erase(vec.begin() + idx);
		}
		else value = vec[idx++];
	}
}

template <typename VType, typename ... PType>
void delete_repeats(std::vector<VType>& vec, std::vector<PType>&... dep_vec)
{
	if (!vec.size()) return;

	VType value = vec[0];

	int idx = 1;
	while (idx < (int)vec.size()) {

		if (vec[idx] == value) {

			erase_vectors_entry(idx, vec, dep_vec...);
		}
		else value = vec[idx++];
	}
}

//----------------------------------------

//reserve memory for vector up to available physical memory
template <typename VType>
bool mreserve_vector(std::vector<VType>& vec, size_t num_elements)
{
	if (num_elements == vec.size()) return true;

	try {

		vec.reserve(num_elements);
	}
	catch (...) {

		return false;
	}

	return true;
}

//resize vector to given number of elements : return false if failed. Only resize if the vector can fit in available physical memory
template <typename VType>
bool malloc_vector(std::vector<VType>& vec, size_t num_elements)
{
	if (num_elements == vec.size()) return true;

	if (!mreserve_vector(vec, num_elements)) return false;

	vec.resize(num_elements);

	//managed to resize : make sure there's no unused capacity
	vec.shrink_to_fit();

	return true;
}

template <typename VType>
bool malloc_vector(std::vector<VType>& vec, size_t num_elements, VType value)
{
	if (!mreserve_vector(vec, num_elements)) return false;

	vec.assign(num_elements, value);

	//managed to resize : make sure there's no unused capacity
	vec.shrink_to_fit();

	return true;
}