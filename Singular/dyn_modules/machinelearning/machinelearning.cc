#include <stdio.h>
#include <readline/history.h>

#include "kernel/mod2.h"
#include "Singular/ipid.h"
#include "Singular/mod_lib.h"

#include "mlpredict.h"

static BOOLEAN openHelpFile(leftv result, leftv arg)
{
	char buffer[100];
	if (arg->rtyp != STRING_CMD) {
		// We really need the argument passed to be a string
		return TRUE;
	}

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
	//define something for Windows (32-bit and 64-bit, this part is common)
	#ifdef _WIN64
		//define something for Windows (64-bit only)
		sprintf(buffer,
			"cmd /c start ~/.singular/helpfiles/singular/html/%s",
			(char *)arg->data);
	#else
		//define something for Windows (32-bit only)
		sprintf(buffer,
			"cmd /c start ~/.singular/helpfiles/singular/html/%s",
			(char *)arg->data);
	#endif
#elif __APPLE__
	#include <TargetConditionals.h>
	#if TARGET_OS_MAC
		// Other kinds of Mac OS
		sprintf(buffer, "open ~/.singular/helpfiles/singular/html/%s",
				(char *)arg->data);
	#else
		Print("Unsupported Apple device:  openHelpFile is not \
				implemented for this device");
		return TRUE;
	#endif
#elif __linux__
	// linux
	sprintf(buffer, "xdg-open ~/.singular/helpfiles/singular/html/%s",
			(char *)arg->data);
#elif __unix__ // all unices not caught above
	// Unix
	Print("Unsupported Unix device:  openHelpFile is not \
			implemented for this device");
		return TRUE;
#elif defined(_POSIX_VERSION)
    // POSIX
#else
	Print("Unsupported Posix device:  openHelpFile is not \
			implemented for this device");
		return TRUE;
#endif
	Print("%s", buffer);
	PrintLn();

	if (system(buffer)) {
		Print("The System complained when running the open command");
		PrintLn();
	}
	result->rtyp=NONE; // set the result type
	return FALSE;
}

static BOOLEAN predictHelp(leftv result, leftv arg)
{
	char *buffer[5];
	int lengths[5];
	char filename[15] = "~/.history";
	int i;
	lists L = (lists)omAllocBin(slists_bin);

	buffer[0] = NULL;
	buffer[1] = NULL;
	buffer[2] = NULL;
	buffer[3] = NULL;
	buffer[4] = NULL;

	if (write_history(NULL)) {
		printf("Failed to write history\n");
		return TRUE;
	}

	if (!ml_initialise()) {
		Print("Initialise of ml failed.");
		PrintLn();
		/* Notify singular that an error occured */
		return TRUE;
	}

	ml_make_prediction(filename, buffer, lengths, _omStrDup);

	L->Init(5);

	for (i = 0; i < 5; i++) {
		//printf("prediction %d: %s\n", i, buffer[i]);
		L->m[i].rtyp = STRING_CMD;
		L->m[i].data = buffer[i];
	}

	// pass the resultant list to the res datastructure
	result->rtyp=LIST_CMD;
	result->data=(void *)L;

	ml_finalise();
	return FALSE;

	//result->rtyp=NONE; // set the result type
}

// initialisation of the module
extern "C" int SI_MOD_INIT(machinelearning)(SModulFunctions* psModulFunctions)
{
	// this is the initialization routine of the module
	// adding the routine predictHelp:
	psModulFunctions->iiAddCproc(
			(currPack->libname? currPack->libname: ""),
			"predictHelp",// for the singular interpreter
			FALSE, // should enter the global name space
			predictHelp); // the C/C++ routine
	psModulFunctions->iiAddCproc(
			(currPack->libname? currPack->libname: ""),
			"openHelpFile",// for the singular interpreter
			FALSE, // should enter the global name space
			openHelpFile); // the C/C++ routine

	return MAX_TOK;
}
