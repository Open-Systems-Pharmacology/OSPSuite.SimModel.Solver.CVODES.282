#include "SimModelSolver_CVODES_2_8_2Specs/ExceptionHelper.h"

using namespace std;

void ExceptionHelper::ThrowExceptionFrom(const string message)
{
	throw gcnew System::Exception(gcnew System::String(message.c_str()));
}

void ExceptionHelper::ThrowExceptionFrom(SimModelSolverErrorData & ED)
{
	string message = "C++ Exception in "+ED.GetSource()+":\n"+ED.GetDescription();
	ExceptionHelper::ThrowExceptionFrom(message);
}

void ExceptionHelper::ThrowExceptionFromUnknown(string source)
{
	string message = "C++ unknown Exception";
	if (source != "")
		message += " in" + source;

	ExceptionHelper::ThrowExceptionFrom(message);
}

void ExceptionHelper::ThrowExceptionFromUnknown()
{
	ExceptionHelper::ThrowExceptionFromUnknown("");
}
