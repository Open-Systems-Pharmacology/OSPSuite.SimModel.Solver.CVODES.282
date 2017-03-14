#ifndef _ExceptionHelper_H_
#define _ExceptionHelper_H_

#include "SimModelSolverBase/SimModelSolverErrorData.h"
#include <string>

ref class ExceptionHelper
{
public:
	static void ThrowExceptionFrom(SimModelSolverErrorData & ED);
	static void ThrowExceptionFrom(const std::string message);
	static void ThrowExceptionFromUnknown(std::string source);
	static void ThrowExceptionFromUnknown();
};

#endif //_ExceptionHelper_H_

