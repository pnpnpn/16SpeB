#ifndef _DATASET_IO_H
#define _DATASET_IO_H

#include "stdinc.h"
#include "symbols.h"

static const int FASTA_HEADER_CAPACITY = 200;

// A C G T
// 0 1 2 3

typedef struct {
	//sequence and site index starts at 0.
	int **seqs; //seqs[numseqs][seqlen[i]] = {0, 1, 2, 3}
	//int numalphas;
	int	numseqs; //number of sequences
	int *seqlen; //length of each sequence
	int maxseqlen;
	int minseqlen;
	double avgseqlen;
	int	count[NUMALPHAS]; // residue count
	double ntFreq[NUMALPHAS]; //nucleotide frequency
	int totalCount; //sum(count) - different from total seqlen over all sequences because of GAP_CHAR
	bool **isBadPos; //sequence positions with W or N

	//reverse-complementary
	bool useRevcompl;
	int minspan;
	int maxspan;

	int **isValidSite[SPAN_CAPACITY]; //[span][numseqs][seqlen[i]]
	int *numValidSites[SPAN_CAPACITY]; //[span][numseqs]

	//headers
	char **headers;
} Dataset;

extern Dataset *openDataset(char *filename, bool useRevcompl, int minspan, int maxspan);
extern Dataset *openBackgroundData(char *filename, bool useRevcompl);
extern void nilDataset(Dataset *);

//reverse-complementary
extern int getConcatPosOfOppStrand(Dataset* data, int seqind, int concatPos, int span);
extern bool isForwardStrand(Dataset* data, int seqind, int concatPos);
extern int getConcat2DoubleStrandPos(Dataset* data, int seqind, int concatPos);
extern int getDoubleStrand2ConcatPos(Dataset *data, int seqind, int dspos);

extern bool isValidConcatPos(Dataset* data, int seqind, int concatPos, int span);
extern bool isBoundedDoubleStrandPos(Dataset* data, int seqind, int dspos);


#endif
