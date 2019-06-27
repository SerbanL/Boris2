#pragma once

//moving mesh algorithm - ends are fixed. This value is the ratio of wire length to set end lengths.
#define MOVEMESH_ENDRATIO	0.1
//default threshold to trigger a mesh shift when using an antisymmetric moving mesh (e.g. domain wall)
#define MOVEMESH_ANTISYMMETRIC_THRESHOLD	0.25
//default threshold to trigger a mesh shift when using a symmetric moving mesh (e.g. skyrmion)
#define MOVEMESH_SYMMETRIC_THRESHOLD	0.15

//parameters for ABM adaptive time step
//above this relative error the evaluation has failed and will be redone with a lower time step
#define ABM_RELERRFAIL	1.5e-4
//above this relative error the time step will be reduced
#define ABM_RELERRMAX	9e-5
//below this relative error the time step will be increased
#define ABM_RELERRMIN	4e-5
//When reducing the time step multiply it with this
#define ABM_DTREDUCE	0.9
//When increasing the time step multiply it with this
#define ABM_DTINCREASE	1.001
//maximum time step the method can reach
#define ABM_MAXDT	1e-12
//minimum time step the method can reach
#define ABM_MINDT	1e-15

//fixed parameters for RKF45 adaptive time step
//above this relative error the evaluation has failed and will be redone with a lower time step
#define RKF_RELERRFAIL	2e-4
//above this relative error the time step will be reduced
#define RKF_RELERRMAX	1e-4
//below this relative error the time step will be increased
#define RKF_RELERRMIN	4e-5
//When reducing the time step multiply it with this
#define RKF_DTREDUCE	0.9
//When increasing the time step multiply it with this
#define RKF_DTINCREASE	1.001
//maximum time step the method can reach
#define RKF_MAXDT	3e-12
//minimum time step the method can reach
#define RKF_MINDT	1e-15

//Available equations to solve enum
enum ODE_ { ODE_ERROR = -1, ODE_LLG, ODE_LLGSTT, ODE_LLB, ODE_LLBSTT, ODE_SLLG, ODE_SLLGSTT, ODE_SLLB, ODE_SLLBSTT, ODE_LLGSA, ODE_SLLGSA, ODE_LLBSA, ODE_SLLBSA };

//ODE evaluation methods enum
enum EVAL_ { EVAL_ERROR = -1, EVAL_EULER, EVAL_TEULER, EVAL_RK4, EVAL_ABM, EVAL_RKF };