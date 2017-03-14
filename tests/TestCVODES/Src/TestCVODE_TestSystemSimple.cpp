//---- solving the system:
//
// y1' = (P1+P2)*y2 + (P3-2)*y1 + (y3-2)
// y2' = y1 + P4
// y3' = 0 + 0
//
// y1(0) = 2
// y2(0) = 0
// y3(0) = 2
//
// P1 = 1
// P2 = 0
// P3 = 2
// P4 = 0
//
// System is aquivalent to:
//
// y1' = y2   y1(0) = 2
// y2' = y1   y2(0) = 0
// y3' = 0    y3(0) = 2
//
// Analytical solution is:
//
// y1 = exp(Time) + exp(-Time)
// y2 = exp(Time) - exp(-Time)
// y3 = 2

#include "HelpFunctions.h"

typedef struct {
  realtype p[4];           /* problem parameters */
} *UserData;

int main(int argc, char *argv[])
{
	const bool sensi = true;
    void *cvode_mem = NULL;
    UserData data = NULL;
    realtype t, tout;
    N_Vector y = NULL;
    int iout, flag, is; 
    N_Vector *yS = NULL;

    //User data structure
    data = (UserData) malloc(sizeof *data);
    if (check_flag((void *)data, "malloc", 2)) return(1);
    data->p[0] = 1.0; data->p[1] = 0.0; data->p[2] = 2.0; data->p[3] = 0.0;

    //Initial conditions
    y = N_VNew_Serial(3);
    if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

    Ith(y,1) = 2.0; Ith(y,2) = 0.0; Ith(y,3) = 2.0;

    //Create CVODES object
    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

    //Allocate space for CVODES
    flag = CVodeInit(cvode_mem, f, 0.0, y);
    if (check_flag(&flag, "CVodeInit", 1)) return(1);

	CVodeSetMaxOrd(cvode_mem, 5);
	CVodeSetMaxNumSteps(cvode_mem, 100000);
	CVodeSetMaxHnilWarns(cvode_mem, 10);
	CVodeSetInitStep(cvode_mem, 0.0);
	CVodeSetMaxStep(cvode_mem, 60.0);
	CVodeSetMinStep(cvode_mem, 0.0);

    //tolerance
    //N_Vector absTol_NV = N_VNew_Serial(3);
    //Ith(absTol_NV, 1) = Ith(absTol_NV, 2) = Ith(absTol_NV, 3) = 1.0e-10;
    //flag = CVodeSVtolerances(cvode_mem, 1.0e-4, absTol_NV);
	flag = CVodeSStolerances(cvode_mem, 1.0e-9, 1.0e-12);
    if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

    //Attach user data
    flag = CVodeSetUserData(cvode_mem, data);
    if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);

    //Attach linear solver
    flag = CVDense(cvode_mem, 3);
    if (check_flag(&flag, "CVDense", 1)) return(1);
//todo - commented code
    //flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
    //if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);

	if (sensi)
	{
		 yS = N_VCloneVectorArray_Serial(4, y);
		 if (check_flag((void *)yS, "N_VCloneVectorArray_Serial", 0)) return(1);
		 for (is=0;is<4;is++) 
		 N_VConst(0.0, yS[is]);

		 //TODO wenn fS implementeirt ist - ersetzen
		 //flag = CVodeSensInit1(cvode_mem, 4, CV_STAGGERED, fS, yS);
		 flag = CVodeSensInit1(cvode_mem, 4, CV_STAGGERED, NULL, yS);
		 if(check_flag(&flag, "CVodeSensInit", 1)) return(1);

		 flag = CVodeSensEEtolerances(cvode_mem);
		 if(check_flag(&flag, "CVodeSensEEtolerances", 1)) return(1);

		 flag = CVodeSetSensErrCon(cvode_mem, false);
		 if (check_flag(&flag, "CVodeSetSensErrCon", 1)) return(1);

		 flag = CVodeSetSensParams(cvode_mem, data->p, NULL, NULL);
		 if (check_flag(&flag, "CVodeSetSensParams", 1)) return(1);
	}
  
    printf("\n\n=======================================================================\n");
    printf("     T     Q       H      NST           y1           y2           y3    \n");
    printf("=======================================================================\n");

	for (iout = 1; iout <= 3; iout++)
	{
		tout = 1.0 * iout;
		flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
		if (check_flag(&flag, "CVode", 1)) break;

		PrintOutput(cvode_mem, t, y);

		if (sensi)
		{
			flag = CVodeGetSens(cvode_mem, &t, yS);
			if (check_flag(&flag, "CVodeGetSens", 1)) break;
			PrintOutputS(yS, 3, 4);
		}

		printf("-----------------------------------------------------------------------\n");
	}

//    PrintFinalStats(cvode_mem, sensi);
	
    return(0);
}

//f routine. Compute f(t,y)
static int f(double t, N_Vector y, N_Vector ydot, void *user_data)
{
	double p1, p2, p3, p4;
	double y1, y2, y3;

	UserData data = (UserData)user_data;
	p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2]; p4 = data->p[3];

	y1 = Ith(y, 1);
	y2 = Ith(y, 2);
	y3 = Ith(y, 3);

	Ith(ydot, 1) = (p1 + p2)*y2 + (p3 - 2.0)*y1 + (y3 - 2.0);
	Ith(ydot, 2) = y1 + p4;
	Ith(ydot, 3) = 0.0;

  return 0;
}

//Jacobian routine. Compute J(t,y)
static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, 
               DlsMat J, void *user_data, 
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
	throw "not implemented";

	//IJth(J, 1, 1) = 0.0;
	//IJth(J, 1, 2) = 1.0;
	//IJth(J, 1, 3) = 0.0;

	//IJth(J, 2, 1) = 1.0;
	//IJth(J, 2, 2) = 0.0;
	//IJth(J, 2, 3) = 0.0;

	//IJth(J, 3, 1) = 0.0;
	//IJth(J, 3, 2) = 0.0;
	//IJth(J, 3, 3) = 0.0;

  //realtype y1, y2, y3;
  //UserData data;
  //realtype p1, p2, p3;
 
  //y1 = Ith(y,1); y2 = Ith(y,2); y3 = Ith(y,3);
  //data = (UserData) user_data;
  //p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];
 
  //IJth(J,1,1) = -p1;  IJth(J,1,2) = p2*y3;          IJth(J,1,3) = p2*y2;
  //IJth(J,2,1) =  p1;  IJth(J,2,2) = -p2*y3-2*p3*y2; IJth(J,2,3) = -p2*y2;
  //                    IJth(J,3,2) = 2*p3*y2;

  return(0);
}
 
//fS routine. Compute sensitivity r.h.s
static int fS(int Ns, realtype t, N_Vector y, N_Vector ydot, 
              int iS, N_Vector yS, N_Vector ySdot, 
              void *user_data, N_Vector tmp1, N_Vector tmp2)
{
	throw "not implemented";
  //UserData data;
  //realtype p1, p2, p3;
  //realtype y1, y2, y3;
  //realtype s1, s2, s3;
  //realtype sd1, sd2, sd3;

  //data = (UserData) user_data;
  //p1 = data->p[0]; p2 = data->p[1]; p3 = data->p[2];

  //y1 = Ith(y,1);  y2 = Ith(y,2);  y3 = Ith(y,3);
  //s1 = Ith(yS,1); s2 = Ith(yS,2); s3 = Ith(yS,3);

  //sd1 = -p1*s1 + p2*y3*s2 + p2*y2*s3;
  //sd3 = 2*p3*y2*s2;
  //sd2 = -sd1-sd3;

  //switch (iS) {
  //case 0:
  //  sd1 += -y1;
  //  sd2 +=  y1;
  //  break;
  //case 1:
  //  sd1 +=  y2*y3;
  //  sd2 += -y2*y3;
  //  break;
  //case 2:
  //  sd2 += -y2*y2;
  //  sd3 +=  y2*y2;
  //  break;
  //}
  //
  //Ith(ySdot,1) = sd1;
  //Ith(ySdot,2) = sd2;
  //Ith(ySdot,3) = sd3;

  //return(0);
}
