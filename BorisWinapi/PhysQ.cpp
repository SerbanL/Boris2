#include "stdafx.h"
#include "PhysQ.h"

//------------------------------------------- Constructors

//disabled mesh display
PhysQ::PhysQ(Rect rect_, DBL3 h_, unsigned displayedType_)
{
	rect = rect_;
	h = h_;
	n = SZ3();
	displayedType = displayedType_;
}

PhysQ::PhysQ(VEC<DBL3> *pVEC_, unsigned displayedType_, VEC3REP_ vec3rep_)
{
	pVEC = pVEC_;
	vtype = VECTYPE_VEC;
	vstype = VECSTYPE_DBL3;

	rect = pVEC_->rect;
	h = pVEC_->h;
	n = pVEC_->n;
	displayedType = displayedType_;

	vec3rep = vec3rep_;
}

PhysQ::PhysQ(VEC<FLT3> *pVEC_, unsigned displayedType_, VEC3REP_ vec3rep_)
{
	pVEC = pVEC_;
	vtype = VECTYPE_VEC;
	vstype = VECSTYPE_FLT3;

	rect = pVEC_->rect;
	h = pVEC_->h;
	n = pVEC_->n;
	displayedType = displayedType_;

	vec3rep = vec3rep_;
}

PhysQ::PhysQ(VEC<DBL3> *pVEC_, VEC<DBL3> *pVEC2_, unsigned displayedType_, VEC3REP_ vec3rep_)
{
	pVEC = pVEC_;
	pVEC2 = pVEC2_;
	vtype = VECTYPE_VEC;
	vstype = VECSTYPE_DBL3;

	rect = pVEC_->rect;
	h = pVEC_->h;
	n = pVEC_->n;
	displayedType = displayedType_;

	vec3rep = vec3rep_;
}

PhysQ::PhysQ(VEC<FLT3> *pVEC_, VEC<FLT3> *pVEC2_, unsigned displayedType_, VEC3REP_ vec3rep_)
{
	pVEC = pVEC_;
	pVEC2 = pVEC2_;
	vtype = VECTYPE_VEC;
	vstype = VECSTYPE_FLT3;

	rect = pVEC_->rect;
	h = pVEC_->h;
	n = pVEC_->n;
	displayedType = displayedType_;

	vec3rep = vec3rep_;
}

PhysQ::PhysQ(VEC<double> *pVEC_, unsigned displayedType_)
{
	pVEC = pVEC_;
	vtype = VECTYPE_VEC;
	vstype = VECSTYPE_DOUBLE;

	rect = pVEC_->rect;
	h = pVEC_->h;
	n = pVEC_->n;
	displayedType = displayedType_;
}

PhysQ::PhysQ(VEC<float> *pVEC_, unsigned displayedType_)
{
	pVEC = pVEC_;
	vtype = VECTYPE_VEC;
	vstype = VECSTYPE_FLOAT;

	rect = pVEC_->rect;
	h = pVEC_->h;
	n = pVEC_->n;
	displayedType = displayedType_;
}

PhysQ::PhysQ(VEC_VC<DBL3> *pVEC_, unsigned displayedType_, VEC3REP_ vec3rep_)
{
	pVEC = pVEC_;
	vtype = VECTYPE_VEC_VC;
	vstype = VECSTYPE_DBL3;

	rect = pVEC_->rect;
	h = pVEC_->h;
	n = pVEC_->n;
	displayedType = displayedType_;

	vec3rep = vec3rep_;
}

PhysQ::PhysQ(VEC_VC<FLT3> *pVEC_, unsigned displayedType_, VEC3REP_ vec3rep_)
{
	pVEC = pVEC_;
	vtype = VECTYPE_VEC_VC;
	vstype = VECSTYPE_FLT3;

	rect = pVEC_->rect;
	h = pVEC_->h; 
	n = pVEC_->n;
	displayedType = displayedType_;

	vec3rep = vec3rep_;
}

PhysQ::PhysQ(VEC_VC<DBL3> *pVEC_, VEC_VC<DBL3> *pVEC2_, unsigned displayedType_, VEC3REP_ vec3rep_)
{
	pVEC = pVEC_;
	pVEC2 = pVEC2_;
	vtype = VECTYPE_VEC_VC;
	vstype = VECSTYPE_DBL3;

	rect = pVEC_->rect;
	h = pVEC_->h;
	n = pVEC_->n;
	displayedType = displayedType_;

	vec3rep = vec3rep_;
}

PhysQ::PhysQ(VEC_VC<FLT3> *pVEC_, VEC_VC<FLT3> *pVEC2_, unsigned displayedType_, VEC3REP_ vec3rep_)
{
	pVEC = pVEC_;
	pVEC2 = pVEC2_;
	vtype = VECTYPE_VEC_VC;
	vstype = VECSTYPE_FLT3;

	rect = pVEC_->rect;
	h = pVEC_->h;
	n = pVEC_->n;
	displayedType = displayedType_;

	vec3rep = vec3rep_;
}

PhysQ::PhysQ(VEC_VC<double> *pVEC_, unsigned displayedType_)
{
	pVEC = pVEC_;
	vtype = VECTYPE_VEC_VC;
	vstype = VECSTYPE_DOUBLE;

	rect = pVEC_->rect;
	h = pVEC_->h;
	n = pVEC_->n;
	displayedType = displayedType_;
}

PhysQ::PhysQ(VEC_VC<float> *pVEC_, unsigned displayedType_)
{
	pVEC = pVEC_;
	vtype = VECTYPE_VEC_VC;
	vstype = VECSTYPE_FLOAT;

	rect = pVEC_->rect;
	h = pVEC_->h;
	n = pVEC_->n;
	displayedType = displayedType_;
}

PhysQ& PhysQ::operator=(const PhysQ& copyThis)
{
	pVEC = copyThis.pVEC;
	pVEC2 = copyThis.pVEC2;

	vtype = copyThis.vtype;
	vstype = copyThis.vstype;

	displayedType = copyThis.displayedType;
	displayTransparency = copyThis.displayTransparency;
	displayThresholds = copyThis.displayThresholds;
	displayThresholdTrigger = copyThis.displayThresholdTrigger;

	rect = copyThis.rect;
	h = copyThis.h;
	n = copyThis.n;
	meshInFocus = copyThis.meshInFocus;
	meshName = copyThis.meshName;

	vec3rep = copyThis.vec3rep;

	return *this;
}

//------------------------------------------- Access

//get value at coordinate, where we know if the stored type is vectorial or scalar (e.g. DBL3 vs double)
DBL3 PhysQ::get_vec_point(DBL3 coordinate)
{
	//cast depending on stored type : cast to VEC since VEC_VC is derived from it
	switch (vstype) {

	case VECSTYPE_NONE:
		return DBL3();

	case VECSTYPE_FLT3:
		return (DBL3)reinterpret_cast<VEC<FLT3>*>(pVEC)->weighted_average(coordinate - rect.s, h);

	case VECSTYPE_DBL3:
		return reinterpret_cast<VEC<DBL3>*>(pVEC)->weighted_average(coordinate - rect.s, h);
	}

	return DBL3();
}

double PhysQ::get_sca_point(DBL3 coordinate)
{
	//cast depending on stored type : cast to VEC since VEC_VC is derived from it
	switch (vstype) {

	case VECSTYPE_NONE:
		return 0.0;

	case VECSTYPE_FLOAT:
		return (double)reinterpret_cast<VEC<float>*>(pVEC)->weighted_average(coordinate - rect.s, h);

	case VECSTYPE_DOUBLE:
		return reinterpret_cast<VEC<double>*>(pVEC)->weighted_average(coordinate - rect.s, h);
	}

	return 0.0;
}

//------------------------------------------- Properties

SZ3 PhysQ::size(void)
{
	//cast depending on stored type : cast to VEC since VEC_VC is derived from it
	switch (vstype) {

	case VECSTYPE_NONE:
		return SZ3();

	case VECSTYPE_FLOAT:
		return reinterpret_cast<VEC<float>*>(pVEC)->size();

	case VECSTYPE_DOUBLE:
		return reinterpret_cast<VEC<double>*>(pVEC)->size();

	case VECSTYPE_FLT3:
		return reinterpret_cast<VEC<FLT3>*>(pVEC)->size();

	case VECSTYPE_DBL3:
		return reinterpret_cast<VEC<DBL3>*>(pVEC)->size();
	}

	return SZ3();
}