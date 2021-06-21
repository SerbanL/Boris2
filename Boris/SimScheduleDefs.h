#pragma once

//simulation stage/step settings -> add new values at the end to keep older simulation files compatible
enum SS_ {
	SS_RELAX = 0,
	SS_HFIELDXYZ, SS_HFIELDXYZSEQ, SS_HPOLARSEQ, SS_HFMR,
	SS_V, SS_VSEQ, SS_VSIN, SS_VCOS,
	SS_I, SS_ISEQ, SS_ISIN, SS_ICOS,
	SS_T, SS_TSEQ,
	SS_Q, SS_QSEQ,
	SS_HFIELDEQUATION, SS_HFIELDEQUATIONSEQ,
	SS_VEQUATION, SS_VEQUATIONSEQ,
	SS_IEQUATION, SS_IEQUATIONSEQ,
	SS_TEQUATION, SS_TEQUATIONSEQ,
	SS_QEQUATION, SS_QEQUATIONSEQ,
	SS_HFIELDFILE, SS_VFILE, SS_IFILE, SS_TFILE, SS_QFILE,
	SS_TSIGPOLAR,

	SS_MONTECARLO
};

//simulation stage stop conditions -> add new values at the end to keep older simulation files compatible
enum STOP_ { STOP_NOSTOP = 0, STOP_ITERATIONS, STOP_MXH, STOP_TIME, STOP_DMDT, STOP_MXH_ITER, STOP_DMDT_ITER };

//data save conditions -> add new values at the end to keep older simulation files compatible
enum DSAVE_ { DSAVE_NONE = 0, DSAVE_STAGE, DSAVE_STEP, DSAVE_ITER, DSAVE_TIME };