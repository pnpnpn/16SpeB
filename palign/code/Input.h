#ifndef _INPUT_H
#define _INPUT_H

#include "symbols.h"

#define INVALID_POS -1


class Seqset {
public:
	int **seqs;
	int numseqs;
	int *seqlen;
	int maxseqlen;
	int minseqlen;

	Seqset();
	Seqset(int **seqs, int numseqs, int *seqlen);
	Seqset(const Seqset &src); //copy constructor
	virtual ~Seqset();

	virtual void createAndCopySeqs(int **seqs, int numseqs, int *seqlen);
	virtual void removeGaps();
	static void revcompl(int *oldseq, int *newseq, int len);
	virtual int* createSingleSeq(int &seqlen); 
};

class Nullset : public Seqset {
public:
	Nullset();
	Nullset(const Seqset &src);
	virtual ~Nullset();

	virtual void randomize(); //randomize by permutation
	virtual int incrementIters();
	virtual int getNumIters();

protected:
	Nullset(const Nullset &src); //copy constructor
private:
	int numIters;
};

class Input {
public:
	Input();
	Input(string fastaFilename);
	Input(string fastaFilename, string bgFastaFilename);
	virtual ~Input();

	Seqset *seqset;
	Seqset *gappedSeqset;

	Seqset *bgSeqset;

	//headers 
	vector<string> fastaHeaders; //headers without '>' character

	//superimpose gaps from an original input to null sequence
	void superimposeGaps(int *gapseq, int gapseqLen, int *nongapseq, int nongapseqLen, string &result); 
	bool hasSeparateBgfsa();
protected:

private:
	//int argc;
	//char **argv;
	void openFastaMold(string fastaFilename);
	bool separateBgfsa;

};



#endif

