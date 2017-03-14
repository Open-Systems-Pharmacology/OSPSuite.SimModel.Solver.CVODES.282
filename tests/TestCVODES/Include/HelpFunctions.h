#ifndef __HELP_FUNCTIONS_h
#define __HELP_FUNCTIONS_h

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <cvodes/cvodes.h>           /* prototypes for CVODES fcts. and consts. */
#include <cvodes/cvodes_dense.h>     /* prototype for CVDENSE fcts. and constants */
#include <nvector/nvector_serial.h>  /* defs. of serial NVECTOR fcts. and macros  */
#include <sundials/sundials_types.h> /* def. of type realtype */
#include <sundials/sundials_math.h>  /* definition of ABS */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* i-th vector component i=1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* (i,j)-th matrix component i,j=1..NEQ */

/* Prototypes of functions by CVODES */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int Jac(long int N, realtype t,
	N_Vector y, N_Vector fy,
	DlsMat J, void *user_data,
	N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot,
	int iS, N_Vector yS, N_Vector ySdot,
	void *user_data, N_Vector tmp1, N_Vector tmp2);

static int ewt(N_Vector y, N_Vector w, void *user_data);

/* Prototypes of private functions */

void ProcessArgs(int argc, char *argv[], booleantype *sensi, int *sensi_meth, booleantype *err_con);
void WrongArgs(char *name);
void PrintOutput(void *cvode_mem, realtype t, N_Vector u);
void PrintOutputS(N_Vector *uS, int noOfVariables, int noOfSensitivityParameters);
void PrintFinalStats(void *cvode_mem, booleantype sensi);
int check_flag(void *flagvalue, char *funcname, int opt);


#endif