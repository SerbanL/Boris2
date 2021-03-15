#pragma once

//stencil ratio for extracting end values (length ratio)
#define DWPOS_ENDSTENCIL	0.25

//stencil value for extracting centre position (length ratio)
#define DWPOS_STENCIL	0.2

//threshold for extracting domain wall width (only get dw width for normalised absolute y values falling in tanh profile within these bounds)
#define DWPOS_YTHRESHOLD_MAX	0.9
#define DWPOS_YTHRESHOLD_MIN	0.1

//profile must have a minimum number of points
#define DWPOS_MINPROFILEPOINTS	5
