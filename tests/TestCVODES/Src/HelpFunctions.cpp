#include "HelpFunctions.h"

/*
* Print current t, step count, order, stepsize, and solution.
*/

void PrintOutput(void *cvode_mem, realtype t, N_Vector u)
{
	long int nst;
	int qu, flag;
	realtype hu, *udata;

	udata = NV_DATA_S(u);

	flag = CVodeGetNumSteps(cvode_mem, &nst);
	check_flag(&flag, "CVodeGetNumSteps", 1);
	flag = CVodeGetLastOrder(cvode_mem, &qu);
	check_flag(&flag, "CVodeGetLastOrder", 1);
	flag = CVodeGetLastStep(cvode_mem, &hu);
	check_flag(&flag, "CVodeGetLastStep", 1);

#if defined(SUNDIALS_EXTENDED_PRECISION)
	printf("%8.3Le %2d  %8.3Le %5ld\n", t, qu, hu, nst);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
	printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
#else
	printf("%8.3e %2d  %8.3e %5ld\n", t, qu, hu, nst);
#endif

	printf("                  Solution       ");

#if defined(SUNDIALS_EXTENDED_PRECISION)
	printf("%12.4Le %12.4Le %12.4Le \n", udata[0], udata[1], udata[2]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
	printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
#else
	printf("%12.4e %12.4e %12.4e \n", udata[0], udata[1], udata[2]);
#endif

}

/*
* Print sensitivities.
*/

void PrintOutputS(N_Vector *uS, int noOfVariables, int noOfSensitivityParameters)
{
	realtype *sdata;

	for (int i = 0; i < noOfSensitivityParameters; i++)
	{
		sdata = NV_DATA_S(uS[i]);
		printf("                  Sensitivity %d  ", i+1);

		for (int j = 0; j < noOfVariables; j++)
		{
			printf("%12.4e ", sdata[j]);
		}
		printf("\n");
	}
}

/*
* Print some final statistics from the CVODES memory.
*/

void PrintFinalStats(void *cvode_mem, booleantype sensi)
{
	long int nst;
	long int nfe, nsetups, nni, ncfn, netf;
	long int nfSe, nfeS, nsetupsS, nniS, ncfnS, netfS;
	long int nje, nfeLS;
	int flag;

	flag = CVodeGetNumSteps(cvode_mem, &nst);
	check_flag(&flag, "CVodeGetNumSteps", 1);
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
	check_flag(&flag, "CVodeGetNumRhsEvals", 1);
	flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
	check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
	flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
	check_flag(&flag, "CVodeGetNumErrTestFails", 1);
	flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
	check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
	flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
	check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

	if (sensi) {
		flag = CVodeGetSensNumRhsEvals(cvode_mem, &nfSe);
		check_flag(&flag, "CVodeGetSensNumRhsEvals", 1);
		flag = CVodeGetNumRhsEvalsSens(cvode_mem, &nfeS);
		check_flag(&flag, "CVodeGetNumRhsEvalsSens", 1);
		flag = CVodeGetSensNumLinSolvSetups(cvode_mem, &nsetupsS);
		check_flag(&flag, "CVodeGetSensNumLinSolvSetups", 1);
		flag = CVodeGetSensNumErrTestFails(cvode_mem, &netfS);
		check_flag(&flag, "CVodeGetSensNumErrTestFails", 1);
		flag = CVodeGetSensNumNonlinSolvIters(cvode_mem, &nniS);
		check_flag(&flag, "CVodeGetSensNumNonlinSolvIters", 1);
		flag = CVodeGetSensNumNonlinSolvConvFails(cvode_mem, &ncfnS);
		check_flag(&flag, "CVodeGetSensNumNonlinSolvConvFails", 1);
	}

	flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
	check_flag(&flag, "CVDlsGetNumJacEvals", 1);
	flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
	check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

	printf("\nFinal Statistics\n\n");
	printf("nst     = %5ld\n\n", nst);
	printf("nfe     = %5ld\n", nfe);
	printf("netf    = %5ld    nsetups  = %5ld\n", netf, nsetups);
	printf("nni     = %5ld    ncfn     = %5ld\n", nni, ncfn);

	if (sensi) {
		printf("\n");
		printf("nfSe    = %5ld    nfeS     = %5ld\n", nfSe, nfeS);
		printf("netfs   = %5ld    nsetupsS = %5ld\n", netfS, nsetupsS);
		printf("nniS    = %5ld    ncfnS    = %5ld\n", nniS, ncfnS);
	}

	printf("\n");
	printf("nje    = %5ld    nfeLS     = %5ld\n", nje, nfeLS);

}

/*
* Check function return value.
*    opt == 0 means SUNDIALS function allocates memory so check if
*             returned NULL pointer
*    opt == 1 means SUNDIALS function returns a flag so check if
*             flag >= 0
*    opt == 2 means function allocates memory so check if returned
*             NULL pointer
*/

int check_flag(void *flagvalue, char *funcname, int opt)
{
	int *errflag;

	/* Check if SUNDIALS function returned NULL pointer - no memory allocated */
	if (opt == 0 && flagvalue == NULL) {
		fprintf(stderr,
			"\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
			funcname);
		return(1);
	}

	/* Check if flag < 0 */
	else if (opt == 1) {
		errflag = (int *)flagvalue;
		if (*errflag < 0) {
			fprintf(stderr,
				"\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
				funcname, *errflag);
			return(1);
		}
	}

	/* Check if function returned NULL pointer - no memory allocated */
	else if (opt == 2 && flagvalue == NULL) {
		fprintf(stderr,
			"\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
			funcname);
		return(1);
	}

	return(0);
}
