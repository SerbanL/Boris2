#pragma once

#include "Boris_Enums_Defs.h"
#if COMPILECUDA == 1

#include "BorisCUDALib.h"

//------------------------------

//Functions for calculating demag tensor for a general rectangular prism. 
//(x, y, z) is the distance from the centre of the rectangular prism to the point where the field is obtained - DBL3 xyz
//This is H = -N*M, N being the usual demagnetizing tensor and M the uniform magnetization of the prism.
//(a, b, c) are the dimensions of the rectangular prism - DBL3 abc
__device__ inline cuBReal fx(cuReal3 xyz, cuReal3 abc)
{
	return (abc.y / 2 - xyz.y)*(abc.z / 2 - xyz.z) /
		((abc.x / 2 - xyz.x) * sqrt((abc.x / 2 - xyz.x)*(abc.x / 2 - xyz.x) + (abc.y / 2 - xyz.y)*(abc.y / 2 - xyz.y) + (abc.z / 2 - xyz.z)*(abc.z / 2 - xyz.z)));
}

__device__ inline cuBReal fy(cuReal3 xyz, cuReal3 abc)
{
	return (abc.x / 2 - xyz.x)*(abc.z / 2 - xyz.z) /
		((abc.y / 2 - xyz.y) * sqrt((abc.x / 2 - xyz.x)*(abc.x / 2 - xyz.x) + (abc.y / 2 - xyz.y)*(abc.y / 2 - xyz.y) + (abc.z / 2 - xyz.z)*(abc.z / 2 - xyz.z)));
}

__device__ inline cuBReal fz(cuReal3 xyz, cuReal3 abc)
{
	return (abc.y / 2 - xyz.y)*(abc.x / 2 - xyz.x) /
		((abc.z / 2 - xyz.z) * sqrt((abc.x / 2 - xyz.x)*(abc.x / 2 - xyz.x) + (abc.y / 2 - xyz.y)*(abc.y / 2 - xyz.y) + (abc.z / 2 - xyz.z)*(abc.z / 2 - xyz.z)));
}

__device__ inline cuBReal Fxy(cuReal3 xyz, cuReal3 abc)
{
	return (abc.z / 2 - xyz.z) + sqrt((abc.x / 2 + xyz.x)*(abc.x / 2 + xyz.x) + (abc.y / 2 - xyz.y)*(abc.y / 2 - xyz.y) + (abc.z / 2 - xyz.z)*(abc.z / 2 - xyz.z));
}

__device__ inline cuBReal Fxz(cuReal3 xyz, cuReal3 abc)
{
	return (abc.y / 2 - xyz.y) + sqrt((abc.x / 2 + xyz.x)*(abc.x / 2 + xyz.x) + (abc.y / 2 - xyz.y)*(abc.y / 2 - xyz.y) + (abc.z / 2 - xyz.z)*(abc.z / 2 - xyz.z));
}

__device__ inline cuBReal Fyz(cuReal3 xyz, cuReal3 abc)
{
	return (abc.x / 2 + xyz.x) + sqrt((abc.x / 2 + xyz.x)*(abc.x / 2 + xyz.x) + (abc.y / 2 - xyz.y)*(abc.y / 2 - xyz.y) + (abc.z / 2 - xyz.z)*(abc.z / 2 - xyz.z));
}

__device__ inline cuBReal Nxx(cuReal3 xyz, cuReal3 abc)
{
	cuBReal sum_this;

	sum_this = atan(fx(xyz & cuReal3(-1, 1, 1), abc));
	sum_this += atan(fx(xyz & cuReal3(1, 1, 1), abc));
	sum_this += atan(fx(xyz & cuReal3(-1, -1, 1), abc));
	sum_this += atan(fx(xyz & cuReal3(-1, 1, -1), abc));
	sum_this += atan(fx(xyz & cuReal3(1, -1, 1), abc));
	sum_this += atan(fx(xyz & cuReal3(-1, -1, -1), abc));
	sum_this += atan(fx(xyz & cuReal3(1, 1, -1), abc));
	sum_this += atan(fx(xyz & cuReal3(1, -1, -1), abc));

	return sum_this / (4 * (cuBReal)PI);
}

__device__ inline cuBReal Nyy(cuReal3 xyz, cuReal3 abc)
{
	cuBReal sum_this;

	sum_this = atan(fy(xyz & cuReal3(-1, 1, 1), abc));
	sum_this += atan(fy(xyz & cuReal3(1, 1, 1), abc));
	sum_this += atan(fy(xyz & cuReal3(-1, -1, 1), abc));
	sum_this += atan(fy(xyz & cuReal3(-1, 1, -1), abc));
	sum_this += atan(fy(xyz & cuReal3(1, -1, 1), abc));
	sum_this += atan(fy(xyz & cuReal3(-1, -1, -1), abc));
	sum_this += atan(fy(xyz & cuReal3(1, 1, -1), abc));
	sum_this += atan(fy(xyz & cuReal3(1, -1, -1), abc));

	return sum_this / (4 * (cuBReal)PI);
}

__device__ inline cuBReal Nzz(cuReal3 xyz, cuReal3 abc)
{
	cuBReal sum_this;

	sum_this = atan(fz(xyz & cuReal3(-1, 1, 1), abc));
	sum_this += atan(fz(xyz & cuReal3(1, 1, 1), abc));
	sum_this += atan(fz(xyz & cuReal3(-1, -1, 1), abc));
	sum_this += atan(fz(xyz & cuReal3(-1, 1, -1), abc));
	sum_this += atan(fz(xyz & cuReal3(1, -1, 1), abc));
	sum_this += atan(fz(xyz & cuReal3(-1, -1, -1), abc));
	sum_this += atan(fz(xyz & cuReal3(1, 1, -1), abc));
	sum_this += atan(fz(xyz & cuReal3(1, -1, -1), abc));

	return sum_this / (4 * (cuBReal)PI);
}

__device__ inline cuBReal Nxy(cuReal3 xyz, cuReal3 abc)
{
	cuBReal sum_this = 0.0;

	cuBReal val = Fxy(xyz, abc);
	if (val) sum_this = log(val);

	val = Fxy(xyz, abc & cuReal3(-1, -1, 1));
	if (val) sum_this += log(val);

	val = Fxy(xyz, abc & cuReal3(1, -1, -1));
	if (val) sum_this += log(val);

	val = Fxy(xyz, abc & cuReal3(-1, 1, -1));
	if (val) sum_this += log(val);

	val = Fxy(xyz, abc & cuReal3(1, -1, 1));
	if (val) sum_this -= log(val);

	val = Fxy(xyz, abc & cuReal3(-1, 1, 1));
	if (val) sum_this -= log(val);

	val = Fxy(xyz, abc & cuReal3(1, 1, -1));
	if (val) sum_this -= log(val);

	val = Fxy(xyz, abc & cuReal3(-1, -1, -1));
	if (val) sum_this -= log(val);

	return sum_this / (4 * (cuBReal)PI);
}

__device__ inline cuBReal Nxz(cuReal3 xyz, cuReal3 abc)
{
	cuBReal sum_this = 0.0;

	cuBReal val = Fxz(xyz, abc);
	if (val) sum_this = log(val);

	val = Fxz(xyz, abc & cuReal3(-1, -1, 1));
	if (val) sum_this += log(val);

	val = Fxz(xyz, abc & cuReal3(1, -1, -1));
	if (val) sum_this += log(val);

	val = Fxz(xyz, abc & cuReal3(-1, 1, -1));
	if (val) sum_this += log(val);

	val = Fxz(xyz, abc & cuReal3(1, -1, 1));
	if (val) sum_this -= log(val);

	val = Fxz(xyz, abc & cuReal3(-1, 1, 1));
	if (val) sum_this -= log(val);

	val = Fxz(xyz, abc & cuReal3(1, 1, -1));
	if (val) sum_this -= log(val);

	val = Fxz(xyz, abc & cuReal3(-1, -1, -1));
	if (val) sum_this -= log(val);

	return sum_this / (4 * (cuBReal)PI);
}

__device__ inline cuBReal Nyz(cuReal3 xyz, cuReal3 abc)
{
	cuBReal sum_this = 0.0;

	cuBReal val = Fyz(xyz, abc);
	if (val) sum_this = log(val);

	val = Fyz(xyz, abc & cuReal3(-1, -1, 1));
	if (val) sum_this += log(val);
	
	val = Fyz(xyz, abc & cuReal3(1, -1, -1));
	if (val) sum_this += log(val);
	
	val = Fyz(xyz, abc & cuReal3(-1, 1, -1));
	if (val) sum_this += log(val);
	
	val = Fyz(xyz, abc & cuReal3(1, -1, 1));
	if (val) sum_this -= log(val);
	
	val = Fyz(xyz, abc & cuReal3(-1, 1, 1));
	if (val) sum_this -= log(val);
	
	val = Fyz(xyz, abc & cuReal3(1, 1, -1));
	if (val) sum_this -= log(val);
	
	val = Fyz(xyz, abc & cuReal3(-1, -1, -1));
	if (val) sum_this -= log(val);

	return sum_this / (4 * (cuBReal)PI);
}

#endif