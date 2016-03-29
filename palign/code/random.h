
#ifndef _RANDOM_H
#define _RANDOM_H

#include "mt19937ar.h"

inline void sRandom(unsigned long seed) { init_genrand(seed);} 
inline double Random() { return genrand_real2(); }

#endif

