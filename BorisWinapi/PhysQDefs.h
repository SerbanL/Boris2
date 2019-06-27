#pragma once

#include "CompileFlags.h"
#if GRAPHICS == 1

#include "BorisGraphics.h"

//number of logical units to display on screen along the width
#define DISPLAYUNITS	100

//multipliers for detail level adjustment
#define DETAILNOTCHUP 1.25
#define DETAILNOTCHDN 0.75

//default detail level in order to fit a maximum number of cells along maximum mesh dimensions
#define DETAILDEFAULTCELLS	30

//smallest and largest values of detail level (m). Smaller means greater detail.
#define DETAILVALUEMIN	1e-9
#define DETAILVALUEMAX	1e-6

//reduce the fit by scaling the meter to logical units conversion constant - makes the mesh to view window fit looser (UNITSSCALEBACK = 1 for exact fit)
#define UNITSSCALEBACK	0.75f

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