#include <readline/history.h>

#include "kernel/mod2.h"
#include "Singular/ipid.h"

#include "mlpredict.h"

BOOLEAN predictHelp(leftv result, leftv arg)
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

extern "C" int mod_init(SModulFunctions* psModulFunctions)
{
	// this is the initialization routine of the module
	// adding the routine predictHelp:
	psModulFunctions->iiAddCproc(
			(currPack->libname? currPack->libname: ""),
			"predictHelp",// for the singular interpreter
			FALSE, // should enter the global name space
			predictHelp); // the C/C++ routine
	return 1;
}