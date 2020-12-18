#pragma once

#include "Types_VAL.h"

/////////////////////////////////////////////////////////////////////////////////////////

enum MSHAPE_ {

	//error or no shape
	MSHAPE_NONE = -1,
	
	//2D shapes
	MSHAPE_DISK, MSHAPE_RECT, MSHAPE_TRIANGLE, 
	
	//3D shapes
	MSHAPE_ELLIPSOID, MSHAPE_PYRAMID, MSHAPE_TETRAHEDRON, MSHAPE_CONE, MSHAPE_TORUS
};

enum MSHAPEMETHOD_ {

	//error
	MSHAPEMETHOD_NONE = -1,

	//add into mesh, replacing any existing parts
	MSHAPEMETHOD_ADD,

	//subtract from mesh, deleting any existing parts
	MSHAPEMETHOD_SUB, 

	//only replace existing parts
	MSHAPEMETHOD_AND,

	//add into mesh, but do not replace existing parts
	MSHAPEMETHOD_XOR,

	MSHAPEMETHOD_NUMMETHODS
};

/////////////////////////////////////////////////////////////////////////////////////////

//Definition of a mesh shape - convenient grouping of common parameters
struct MeshShape {

	//Shape identifier - used by the actual shape generator
	MSHAPE_ id = MSHAPE_NONE;

	//Dimensions (in mesh units)
	DBL3 dimensions = DBL3();

	//Centre position (relative coordinates in mesh)
	DBL3 centre_pos = DBL3(); 

	//Rotation of shape using psi (around y), theta (around x) and phi (around z) in degrees. 0, 0, 0 means no rotation
	DBL3 rotation = DBL3(); 

	//Number of repetitions of unit shape as a x y z array. Default value of 1 1 1: no repetitions.
	INT3 repetitions = INT3(1, 1, 1); 

	//Displacement in mesh units if generating an array
	DBL3 displacements = DBL3(); 

	//method used by shape generator to combine with other shapes
	MSHAPEMETHOD_ method = MSHAPEMETHOD_NONE;

	/////////////////////////////////////////////////////////////////////////////////////////

	//void
	MeshShape(void) {}

	/////////////////////////////////////////////////////////////////////////////////////////

	//full ctor
	MeshShape(MSHAPE_ id, DBL3 dimensions, DBL3 centre_pos, DBL3 rotation, INT3 repetitions, DBL3 displacements, MSHAPEMETHOD_ method)
	{
		this->id = id;
		this->dimensions = dimensions;
		this->centre_pos = centre_pos;
		this->rotation = rotation;
		this->repetitions = repetitions;
		this->displacements = displacements;
		this->method = method;
	}

	/////////////////////////////////////////////////////////////////////////////////////////

	//set shape only : by id
	MeshShape(MSHAPE_ id)
	{
		this->id = id;
	}

	//set shape only : by name
	MSHAPE_ set_shape(std::string name)
	{
		if (name == "disk") id = MSHAPE_DISK;
		else if (name == "rect") id = MSHAPE_RECT;
		else if (name == "triangle") id = MSHAPE_TRIANGLE;
		else if (name == "ellipsoid") id = MSHAPE_ELLIPSOID;
		else if (name == "pyramid") id = MSHAPE_PYRAMID;
		else if (name == "tetrahedron") id = MSHAPE_TETRAHEDRON;
		else if (name == "cone") id = MSHAPE_CONE;
		else if (name == "torus") id = MSHAPE_TORUS;
		else id = MSHAPE_NONE;

		return id;
	}

	/////////////////////////////////////////////////////////////////////////////////////////

	//set shape method only : by id
	MeshShape(MSHAPEMETHOD_ method)
	{
		this->method = method;
	}

	//set shape method only : by name
	MSHAPEMETHOD_ set_method(std::string name)
	{
		if (name == "add") method = MSHAPEMETHOD_ADD;
		else if (name == "sub") method = MSHAPEMETHOD_SUB;
		else if (name == "and") method = MSHAPEMETHOD_AND;
		else if (name == "xor") method = MSHAPEMETHOD_XOR;
		else method = MSHAPEMETHOD_NONE;

		return method;
	}

	/////////////////////////////////////////////////////////////////////////////////////////

	//get shape name from set id
	std::string name(void)
	{
		switch (id) {

		case MSHAPE_DISK:
			return "disk";
			break;

		case MSHAPE_RECT:
			return "rect";
			break;

		case MSHAPE_TRIANGLE:
			return "triangle";
			break;

		case MSHAPE_ELLIPSOID:
			return "ellipsoid";
			break;

		case MSHAPE_PYRAMID:
			return "pyramid";
			break;

		case MSHAPE_TETRAHEDRON:
			return "tetrahedron";
			break;

		case MSHAPE_CONE:
			return "cone";
			break;

		case MSHAPE_TORUS:
			return "torus";
			break;

		default:
			return "";
			break;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////

	//get shape method name from set id
	std::string method_name(void)
	{
		switch (method) {

		case MSHAPEMETHOD_ADD:
			return "add";
			break;

		case MSHAPEMETHOD_SUB:
			return "sub";
			break;

		case MSHAPEMETHOD_AND:
			return "and";
			break;

		case MSHAPEMETHOD_XOR:
			return "xor";
			break;

		default:
			return "";
			break;
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////
};