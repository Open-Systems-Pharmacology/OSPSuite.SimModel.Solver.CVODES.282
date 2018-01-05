#ifndef __SimModelSolver_CVODES282_H_
#define __SimModelSolver_CVODES282_H_

#ifdef _WINDOWS
#pragma warning(disable:4786)
#endif

#include "cvodes/cvodes.h"    /* prototypes for CVodeMalloc, CVode, and CVodeFree, */
                      /* constants OPT_SIZE, BDF, NEWTON, SV, SUCCESS,     */
                      /* NST, NFE, NSETUPS, NNI, NCFN, NETF                */
#include "cvodes/cvodes_dense.h"  /* prototype for CVDense, constant DENSE_NJE         */
#include "cvodes/cvodes_band.h"  /* prototype for CVBand         */
#include "nvector/nvector_serial.h"

#include "SimModelSolverBase/SimModelSolverBase.h"
#include "SimModelSolverBase/SimModelSolverErrorData.h"

#ifdef _WINDOWS
#define CVODES282_EXPORT __declspec(dllexport)
#endif
#ifdef linux
#define CVODES282_EXPORT 
#endif

class UserData;

class SimModelSolver_CVODES282 : public SimModelSolverBase
{
private:

	//Density matrix used for conversion between dense and band jacobian
	DlsMat _interimDenseJacobian;

	//type of linear multistep method to be used (ADAMS or BDF)
	int _lmm;

	//type of iteration used to solve non linear system that arises during each internal time step
	//legal values are FUNCTIONAL and NEWTON
	int _iter;

	int _maxOrd;

	int _mxHNil;

	//absolute tolerance
	N_Vector _absTol_NV;      
    
	//initial values of ODE system
	N_Vector _initialData;
	
	//solution vector
	N_Vector _solution;

	N_Vector * _sensitivityValues;

	//a pointer to CVODE problem memory
	void * _cvodeMem;

	//relative tolerance
	realtype _relTol_CVODE;

	//fill CVOde specific Solver options in datatypes required by CVODE
	void FillSolverOptions(void);

	//Call to Rhs function
	static int Rhs (realtype t, N_Vector y, N_Vector ydot, void * user_data);

	//Call to jacobian function for the dense linear solver
	static int CVODE_JacFn_Dense(long int N, 
		                         realtype t,
                                 N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

	//Call to jacobian function for the band linear solver
	static int CVODE_JacFn_Band(long int N, 
		                        long int mupper, long int mlower,
		                        realtype t, 
								N_Vector y, N_Vector fy, DlsMat J, void *user_data,
                                N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

	//Call to sensitivity RHS function
	static int CVODE_SensitivityRhsFunction(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, 
		                                    N_Vector yS, N_Vector ySdot, void *user_data,
		                                    N_Vector tmp1, N_Vector tmp2);

	std::string ToString (double dValue);

	void setupSensitivityProblem();

public:
	UserData * CVODES_UserData;

	CVODES282_EXPORT SimModelSolver_CVODES282(ISolverCaller * pSolverCaller, int problemSize, int numberOfSensitivityParameters);

	CVODES282_EXPORT virtual ~SimModelSolver_CVODES282 ();

	//Get vector with information of "Non-Standard" solver options
	CVODES282_EXPORT std::vector < OptionInfo > GetSolverOptionsInfo ();

	//-----------------------------------------------------------------------------------------------------
	//Solver dependent initialization routine
	//MUST be called before first call to PerformSolverStep after all solver properties are set
	//
	//Inherited class Init routine:
	// - MUST call SimModelSolverInterface::Init() first
	// - MUST set initialized status = true at the end (in case of success)
	//-----------------------------------------------------------------------------------------------------
	CVODES282_EXPORT void Init ();

	//-----------------------------------------------------------------------------------------------------
	//Get solution of ODE/DDE system at the "next" timepoint.
	// - [IN] tout: next time at which a solution is required
	// - [OUT] tret: time point reached by solver:
	//                        - returns tout, if no error occured
	//                        - returns value < tout otherwise (e.g. max. no. of internal solver steps reached; 
	//                                                                                     no convergence etc.)
	// - [OUT] y: Solution vector at time tret
	// - [OUT] yS: Parameter sensitivities at time tret. yS[i][j]=dy_i/dp_j
	//Returns:
	// - 0 if successful
	// - positive value if a recoverable error occurred (e.g. max. no. of internal solver steps reached)
	// - negative value if an unrecoverable error occurred (e.g. illegal input)
	//-----------------------------------------------------------------------------------------------------
	CVODES282_EXPORT int PerformSolverStep(double tout, double * y, double ** yS, double & tret);

	//-----------------------------------------------------------------------------------------------------
	//Reinitialize DE system (e.g. in case of bigger discontinuities)
	//New relative / absolute tolerance should be set by caller prior to ReInit (if required)
	// - [IN] t0: continue integration from this time point
	// - [IN] y0: new initial value at t0
	//Returns:
	// - 0 if successful
	// - positive value if a recoverable error occurred 
	// - negative value if an unrecoverable error occurred 
	//
	//Inherited class:
	// - MUST call SimModelSolverInterface::ReInit first
	// - MUST perform solver-dependent checks (if applies)
	//-----------------------------------------------------------------------------------------------------
	CVODES282_EXPORT int ReInit (double t0, const std::vector < double > & y0);

	//Solver dependent clean up routine (clear memory etc.)
	CVODES282_EXPORT void Terminate ();

	//Retrieve solver error message from any error code returned by PerformSolverStep or ReInit
	CVODES282_EXPORT std::string GetSolverErrMsg (int SolverRetVal);

	CVODES282_EXPORT void SetOption(const std::string & name, double value);

	CVODES282_EXPORT SimModelSolverErrorData::errNumber GetErrorNumberFromSolverReturnValue(int solverRetVal);

	DlsMat InterimDenseJacobian(void);
};

class UserData
{
public:
	SimModelSolver_CVODES282 * Solver;
	double * SensitivityParameters;
	double * ScalingFactors;
	UserData();
	virtual ~UserData();
};

#endif