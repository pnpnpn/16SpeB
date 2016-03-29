#ifndef _STDINC_H
#define _STDINC_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <assert.h>
#include <time.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#if DEBUG
	#define DEBUG0 1 //sanity check
#else
	#define DEBUG0 0
#endif

#if VERBOSE
	#define DEBUG1 1 //verbose mode
#else
	#define DEBUG1 0
#endif

#define DEBUG2 0

#define TRUE true
#define FALSE false

#define NUMALPHAS 4
#define NUMCHARS 5 //includes extra gap characters
#define GAP_CHAR 4
#define SPAN_MIN_CAPACITY 3
#define SPAN_CAPACITY 1000
#define LOGZERO -1e100

typedef bool boolean;

using namespace std;

#endif

