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
#define AHEUN_RELERRFAIL	1.5e-4
//above this relative error the time step will be reduced
#define AHEUN_RELERRMAX	9e-5
//below this relative error the time step will be increased
#define AHEUN_RELERRMIN	8e-5
//When increasing the time step multiply it with this
#define AHEUN_DTINCREASE	1.001
//maximum time step the method can reach
#define AHEUN_MAXDT	1e-12
//minimum time step the method can reach
#define AHEUN_MINDT	1e-15

//default dT
#define RK4_DEFAULT_DT	0.5e-12

//parameters for ABM adaptive time step
//above this relative error the evaluation has failed and will be redone with a lower time step
#define ABM_RELERRFAIL	1.5e-4
//above this relative error the time step will be reduced
#define ABM_RELERRMAX	9e-5
//below this relative error the time step will be increased
#define ABM_RELERRMIN	4e-5
//When increasing the time step multiply it with this
#define ABM_DTINCREASE	1.001
//maximum time step the method can reach
#define ABM_MAXDT	1e-12
//minimum time step the method can reach
#define ABM_MINDT	1e-15
//default dT
#define ABM_DEFAULT_DT	0.1e-12

//fixed parameters for RK23 adaptive time step
//above this relative error the evaluation has failed and will be redone with a lower time step
#define RK23_RELERRFAIL	1e-4
//above this relative error the time step will be reduced
#define RK23_RELERRMAX	5e-5
//below this relative error the time step will be increased
#define RK23_RELERRMIN	1e-5
//When increasing the time step multiply it with this
#define RK23_DTINCREASE	1.001
//maximum time step the method can reach
#define RK23_MAXDT	2e-12
//minimum time step the method can reach
#define RK23_MINDT	1e-15
//default dT
#define RK23_DEFAULT_DT	0.1e-12

//fixed parameters for RKF45 adaptive time step
//above this relative error the evaluation has failed and will be redone with a lower time step
#define RKF_RELERRFAIL	1e-4
//above this relative error the time step will be reduced
#define RKF_RELERRMAX	5e-5
//below this relative error the time step will be increased
#define RKF_RELERRMIN	1e-5
//When increasing the time step multiply it with this
#define RKF_DTINCREASE	1.001
//maximum time step the method can reach
#define RKF_MAXDT	3e-12
//minimum time step the method can reach
#define RKF_MINDT	1e-15
//default dT
#define RKF_DEFAULT_DT	0.5e-12

//fixed parameters for RKCK45 adaptive time step
//above this relative error the evaluation has failed and will be redone with a lower time step
#define RKCK_RELERRFAIL	1e-4
//above this relative error the time step will be reduced
#define RKCK_RELERRMAX	5e-5
//below this relative error the time step will be increased
#define RKCK_RELERRMIN	1e-5
//When increasing the time step multiply it with this
#define RKCK_DTINCREASE	1.001
//maximum time step the method can reach
#define RKCK_MAXDT	3e-12
//minimum time step the method can reach
#define RKCK_MINDT	1e-15
//default dT
#define RKCK_DEFAULT_DT	0.5e-12

//fixed parameters for RKDP54 adaptive time step
//above this relative error the evaluation has failed and will be redone with a lower time step
#define RKDP_RELERRFAIL	1e-4
//above this relative error the time step will be reduced
#define RKDP_RELERRMAX	5e-5
//below this relative error the time step will be increased
#define RKDP_RELERRMIN	1e-5
//When increasing the time step multiply it with this
#define RKDP_DTINCREASE	1.001
//maximum time step the method can reach
#define RKDP_MAXDT	3e-12
//minimum time step the method can reach
#define RKDP_MINDT	1e-15
//default dT
#define RKDP_DEFAULT_DT	0.5e-12

//default dT -> for the SD solver this acts as the starting timestep and the value it resets to when needed
#define SD_DEFAULT_DT	1e-15

//Available equations to solve enum - to keep bsm files backward compatible add new entries at the end
enum ODE_ { ODE_ERROR = -1, ODE_LLG, ODE_LLGSTT, ODE_LLB, ODE_LLBSTT, ODE_SLLG, ODE_SLLGSTT, ODE_SLLB, ODE_SLLBSTT, ODE_LLGSA, ODE_SLLGSA, ODE_LLBSA, ODE_SLLBSA, ODE_LLGSTATIC };

//ODE evaluation methods enum - to keep bsm files backward compatible add new entries at the end
enum EVAL_ { EVAL_ERROR = -1, EVAL_EULER, EVAL_TEULER, EVAL_RK4, EVAL_ABM, EVAL_RKF, EVAL_RK23, EVAL_SD, EVAL_AHEUN, EVAL_RKCK, EVAL_RKDP };

//EVALSPEEDUP_NONE : evaluate all fields every step (default)
//EVALSPEEDUP_ACCURATE : skip field evaluation but only if accuracy not "significantly" affected
//EVALSPEEDUP_AGGRESIVE : allow some moderate loss of accuracy by skipping more steps where applicable (in many cases this is the same as EVALSPEEDUP_ACCURATE)
//EVALSPEEDUP_EXTREME : only one field evaluation per time step, irrespective of how many steps the evaluation method uses
enum EVALSPEEDUP_ { EVALSPEEDUP_NONE = 0, EVALSPEEDUP_ACCURATE, EVALSPEEDUP_AGGRESIVE, EVALSPEEDUP_EXTREME, EVALSPEEDUP_NUMENTRIES };

//return values for Check_Step_Update method in DiffEq
//EVALSPEEDUPSTEP_SKIP : do not update field, use previous calculation if available
//EVALSPEEDUPSTEP_COMPUTE_NO_SAVE : update field and do not save calculation for next step (since at the next step we'll have to calculate field again so no point saving it)
//EVALSPEEDUPSTEP_COMPUTE_AND_SAVE : update field and save calculation for next time (since at the next step we'll need to re-use calculation)
enum EVALSPEEDUPSTEP_ { EVALSPEEDUPSTEP_SKIP = 0, EVALSPEEDUPSTEP_COMPUTE_NO_SAVE, EVALSPEEDUPSTEP_COMPUTE_AND_SAVE };