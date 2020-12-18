#pragma once

#include <string>
#include <sstream>
#include <cuda_runtime.h>

//Example use:
//
//std::string cuda_error_string = gpuErrchk(cudaPeekAtLastError());
//if(cuda_error_string.length()) { ... show error std::string ... }
//

#define gpuErrchk(ans) gpuAssert((ans), __FILE__, __LINE__);

inline std::string gpuAssert(cudaError_t code, const char *file, int line)
{
	if (code != cudaSuccess) {

		std::stringstream ss;

		ss << std::string(cudaGetErrorString(code)) << " : " << std::string(file) << " : " << line;
	
		return ss.str();
	}
	else return "";
}