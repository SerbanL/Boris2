#pragma once

//moving mesh algorithm - ends are fixed. This value is the ratio of wire length to set end lengths.
#define MOVEMESH_ENDRATIO	0.1
//default threshold to trigger a mesh shift when using an antisymmetric moving mesh (e.g. domain wall)
#define MOVEMESH_ANTISYMMETRIC_THRESHOLD	0.25
//default threshold to trigger a mesh shift when using a symmetric moving mesh (e.g. skyrmion)
#define MOVEMESH_SYMMETRIC_THRESHOLD	0.15

//default dT
#define EULER_DEFAULT_DT	0.3e-14

//default dT
#define TEULER_DEFAULT_DT	0.5e-13

//default dT
#define AHEUN_DEFAULT_DT	0.5e-13
//above this relative error the evaluation has failed and will be redone with a lower time step
#define AHEUN_RELERRFAIL	1e-4
//When increasing the time step limit by this
#define AHEUN_DTINCREASE	2
//maximum time step the method can reach
#define AHEUN_MAXDT	1e-12
//minimum time step the method can reach
#define AHEUN_MINDT	1e-15

//default dT
#define RK4_DEFAULT_DT	0.5e-12

//parameters for ABM adaptive time step
#define ABM_RELERRFAIL	1e-4
#define ABM_DTINCREASE	2
#define ABM_MAXDT	1e-12
#define ABM_MINDT	1e-15
#define ABM_DEFAULT_DT	0.1e-12

//fixed parameters for RK23 adaptive time step
#define RK23_RELERRFAIL	1e-4
#define RK23_DTINCREASE	2
#define RK23_MAXDT	2e-12
#define RK23_MINDT	1e-15
#define RK23_DEFAULT_DT	0.1e-12

//fixed parameters for RKF45 adaptive time step
#define RKF_RELERRFAIL	1e-4
#define RKF_DTINCREASE	2
#define RKF_MAXDT	3e-12
#define RKF_MINDT	1e-15
#define RKF_DEFAULT_DT	0.5e-12

//fixed parameters for RKCK45 adaptive time step
#define RKCK_RELERRFAIL	1e-4
#define RKCK_DTINCREASE	2
#define RKCK_MAXDT	3e-12
#define RKCK_MINDT	1e-15
#define RKCK_DEFAULT_DT	0.5e-12

//fixed parameters for RKDP54 adaptive time step
#define RKDP_RELERRFAIL	1e-4
#define RKDP_DTINCREASE	2
#define RKDP_MAXDT	3e-12
#define RKDP_MINDT	1e-15
#define RKDP_DEFAULT_DT	0.5e-12

//default dT -> for the SD solver this acts as the starting timestep and the value it resets to when needed
#define SD_DEFAULT_DT	1e-15
#define SD_MAXDT	1e-9
#define SD_MINDT	SD_DEFAULT_DT

//difficult to simulate when temperature is very close to the Curie temperature due to numerical instability, especially with stochastic equations. instead use an epsilon approach (units of Kelvin).
#define TCURIE_EPSILON	0.5

//Available equations to solve enum - to keep bsm files backward compatible add new entries at the end
enum ODE_ {

	ODE_ERROR = -1, 
	
	ODE_LLG, ODE_LLGSTT, 
	
	ODE_LLB, ODE_LLBSTT, 
	
	ODE_SLLG, ODE_SLLGSTT, 
	
	ODE_SLLB, ODE_SLLBSTT, 
	
	ODE_LLGSA, ODE_SLLGSA, 
	
	ODE_LLBSA, ODE_SLLBSA, 
	
	ODE_LLGSTATIC, ODE_LLGSTATICSA
};

//ODE evaluation methods enum - to keep bsm files backward compatible add new entries at the end
enum EVAL_ { 

	EVAL_ERROR = -1, 
	
	//Fixed time-step methods
	EVAL_EULER = 0, EVAL_TEULER = 1, EVAL_RK4 = 2, 
	
	//Adpative 2nd order
	EVAL_AHEUN = 7,
	
	//Adaptive linear multistep 2nd order
	EVAL_ABM = 3,

	//Adaptive embedded error estimator, 2nd order
	EVAL_RK23 = 5, 

	//Adaptive embedded error estimator, 4th order
	EVAL_RKF45 = 4, EVAL_RKCK45 = 8, 

	//Adaptive embedded error estimator, 5th order
	EVAL_RKDP54 = 9, EVAL_RKF56 = 10,

	//Energy minimizers
	EVAL_SD = 6

}; //Current maximum : 10

//EVALSPEEDUP_NONE : evaluate all fields every step (default)
//EVALSPEEDUP_STEP : use previously computed demag field
//EVALSPEEDUP_LINEAR : linear interpolation using 2 previously computed demag fields
//EVALSPEEDUP_QUADRATIC : quadratic (polynomial) interpolation using 3 previously computed demag fields
enum EVALSPEEDUP_ { EVALSPEEDUP_NONE = 0, EVALSPEEDUP_STEP, EVALSPEEDUP_LINEAR, EVALSPEEDUP_QUADRATIC, EVALSPEEDUP_CUBIC, EVALSPEEDUP_QUARTIC, EVALSPEEDUP_QUINTIC, EVALSPEEDUP_NUMENTRIES };