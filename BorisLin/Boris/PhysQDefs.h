#pragma once

//options for representing vectorial quantities in a PhysQ, used by PhysQRep when calculating a representation : 
//default is VEC3REP_FULL, which means represent it as a vectorial quantity
//VEC3REP_X : show only the x component as a scalar quantity
//VEC3REP_Y : show only the y component as a scalar quantity
//VEC3REP_Z : show only the z component as a scalar quantity
//VEC3REP_DIRECTION : show only the direction component as a scalar quantity (i.e. color coding of vector orientation only)
//VEC3REP_MAGNITUDE : show magnitude as a scalar quantity
//VEC3REP_NUMOPTIONS gives the number of options available
enum VEC3REP_ { VEC3REP_FULL = 0, VEC3REP_X, VEC3REP_Y, VEC3REP_Z, VEC3REP_DIRECTION, VEC3REP_MAGNITUDE, VEC3REP_NUMOPTIONS };

#include "CompileFlags.h"
#if GRAPHICS == 1

#include "BorisGraphics.h"

//number of logical units to display on screen along the width
#define DISPLAYUNITS	100

//multipliers for detail level adjustment
#define DETAILNOTCHUP 1.25
#define DETAILNOTCHDN 0.75

//aim to draw in maximum this amount of time (ms). If this is exceeded then try to reduce level of detail.
#define MAXDRAWTIMEALLOWED_MS	1000

//If there are too many cells to draw then start using rendering speedup tricks so drawing time remains reasonable
//Level 1: start using elements with fewer vertexes, e.g. cones instead of arrows
#define NUMDRAWNCELLS_RENDERSPEEDUP1	0
//Level 2: don't draw surrounded cells
#define NUMDRAWNCELLS_RENDERSPEEDUP2	1000000
//Level 3: only draw cells on a checkerboard pattern
#define NUMDRAWNCELLS_RENDERSPEEDUP3	2000000

//default detail level in order to fit a maximum number of cells along maximum mesh dimensions
#define DETAILDEFAULTCELLS	30

//smallest and largest values of detail level (m). Smaller means greater detail.
#define DETAILVALUEMIN	1e-10
#define DETAILVALUEMAX	1e-6

//reduce the fit by scaling the meter to logical units conversion constant - makes the mesh to view window fit looser (UNITSSCALEBACK = 1 for exact fit)
#define UNITSSCALEBACK	0.75f
#define UNITSSCALEBACK_ZOOMED 0.9f

//Refresh time for animations (in ms)
#define ANIMATIONREFRESH_MS	10.0f

//animation duration when chaning mesh focus (in ms)
#define FOCUSCHANGEDURATION_MS	500.0f

//cutoff for tanh animation curve
#define TANHCUTOFF	1.5f

//alpha blending based on object position from camera as : alpha = (distance from camera / (CAMERACUTOFFDISTANCE)) ^ ALPHABLENDPOW up to alpha = 1
#define CAMERACUTOFFDISTANCE	120.0f
#define ALPHABLENDPOW	3.0f
//maximum alpha to use when rendering mesh objects
#define ALPHAMAXIMUM	0.95f

#define HIGHLIGHTCELLSCALE	1.05f
#define HIGHLIGHTCELLCOLOR	XMFLOAT4(0.3f, 0.3f, 0.3f, 0.7f)

#define MESHFRAMECOLOR	XMFLOAT4(0.3f, 0.3f, 0.3f, 0.5f)

#endif