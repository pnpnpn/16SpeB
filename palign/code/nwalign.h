#ifndef _NWALIGN_H
#define _NWALIGN_H

#include "stdinc.h"

typedef struct {
	int *align1;
	int *align2;
	int len;
	int capacity;

} AlignPair;

typedef struct {
	int match;
	int mismatch;
	int gapopen;
	int gapext;

	double **dpm; //capacity by capacity
	double **Ix; //capacity by capacity
	double **Iy; //capacity by capacity
	int **tb_dpm;
	int **tb_Ix;
	int **tb_Iy;
	int matrix_capacity; //seq_maxlen + 1 because 0 positions are for no alignment

} NWAlignParams;

//return aligned sequence-pair
extern void nwalign(NWAlignParams*, int *seq1, int len1, int *seq2, int len2, AlignPair* result);

extern AlignPair* constructAlignPair(int len1, int len2);
extern void nilAlignPair(AlignPair *alignPair);

extern NWAlignParams* constructNWAlignParams(int match, int mismatch, int gapopen, int gapext, int seq_maxlen);
extern void nilNWAlignParams(NWAlignParams *params);

//defined as number of identities divded by number of non-gap aligned characters
extern double computePidOverNongap(int *align1, int *align2, int len);
extern double computePidOverAlignlen(int *align1, int *align2, int len);

#endif
