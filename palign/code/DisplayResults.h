#ifndef _DISPLAY_RESULTS_H
#define _DISPLAY_RESULTS_H

#include "stdinc.h"
#include "Params.h"

class Results {
public:
	Results();
	Results(Params *params);
	virtual ~Results();

	virtual void display(Nullset *nullset, int nullsetIndex);
	static void displaySeq(ostream &out, int *seq, int len);
private:
	Params *params; //pointer - do not deallocate
	void _displayNullset(Nullset *nullset, int nullsetIndex, Params *params);
};

#endif

