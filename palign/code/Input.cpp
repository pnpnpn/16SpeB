using namespace std;

#include "dataset.h"
#include "Input.h"
#include "random.h"

//--------------------------------------------------------------------------
// Seqset
//--------------------------------------------------------------------------

static
int _getMaxSeqlen(int *seqlen, int numseqs) {
	int max = INT_MIN;
	for(int i = 0; i < numseqs; i++) {
		if(seqlen[i] > max) {
			max = seqlen[i];
		}
	}
	return max;
}

static
int _getMinSeqlen(int *seqlen, int numseqs) {
	int min = INT_MAX;
	for(int i = 0; i < numseqs; i++) {
		if(seqlen[i] < min) {
			min = seqlen[i];
		}
	}
	return min;
}

Seqset::Seqset() {
	this->maxseqlen = 0;
	this->minseqlen = 0;
	this->numseqs = 0;
	this->seqlen = NULL;
	this->seqs = NULL;
}

Seqset::Seqset(const Seqset &src) {
	this->createAndCopySeqs(src.seqs, src.numseqs, src.seqlen);
	
	this->minseqlen = src.minseqlen;
	this->maxseqlen = src.maxseqlen;

	if(DEBUG0) {
		assert(this->minseqlen == _getMinSeqlen(this->seqlen, this->numseqs));
		assert(this->maxseqlen == _getMaxSeqlen(this->seqlen, this->numseqs));
	}

}

Seqset::Seqset(int **seqs, int numseqs, int *seqlen) {
	this->createAndCopySeqs(seqs, numseqs, seqlen);
	
	this->minseqlen = _getMinSeqlen(this->seqlen, this->numseqs);
	this->maxseqlen = _getMaxSeqlen(this->seqlen, this->numseqs);
}

void Seqset::createAndCopySeqs(int **seqs, int numseqs, int *seqlen) {
	this->numseqs = numseqs;

	//seqlen
	this->seqlen = new int[this->numseqs];
	for(int i = 0; i < this->numseqs; i++) {
		this->seqlen[i] = seqlen[i];
	}

	//seqs
	this->seqs = new int*[this->numseqs];
	for(int i = 0; i < this->numseqs; i++) {
		this->seqs[i] = new int[seqlen[i]];
	}
	for(int i = 0; i < this->numseqs; i++) {
		for(int j = 0; j < this->seqlen[i]; j++) {
			this->seqs[i][j] = seqs[i][j];
		}
	}
}
Seqset::~Seqset() {
	delete[] this->seqlen;
	for(int i = 0; i < this->numseqs; i++) {
		delete[] this->seqs[i];
	}
	delete[] this->seqs;
}

int* Seqset::createSingleSeq(int &newSeqlen) {
	newSeqlen = 0;
	for(int i = 0; i < numseqs; i++) {
		newSeqlen += seqlen[i];
	}
	int *newSeq = new int[newSeqlen];

	int count = 0;
	for(int i = 0; i < numseqs; i++) {
		for(int j = 0; j < seqlen[i]; j++) {
			newSeq[count] = seqs[i][j];
			if(DEBUG0) {
				assert(count < newSeqlen);
			}
			count++;
		}
	}

	return newSeq;
}


void Seqset::revcompl(int *oldseq, int *newseq, int len) {
	//assume newseq has at least len

	//complement
	for(int i = 0; i < len; i++) {
		newseq[i] = complChar(oldseq[i]);
	}

	//reverse
	int left = 0; 
	int right = len-1;
	while(left < right) {
		int temp = newseq[left];
		newseq[left] = newseq[right];
		newseq[right] = temp;

		left++;
		right--;
	}
}



void Seqset::removeGaps() {
	//this algorithm is in-place
	for(int i = 0; i < this->numseqs; i++) {
		int count = 0; //number of non-gaps
		for(int j = 0; j < this->seqlen[i]; j++) {
			if(this->seqs[i][j] != GAP_CHAR) {
				this->seqs[i][count++] = this->seqs[i][j];
			}
		}
		this->seqlen[i] = count;
	}
	this->minseqlen = _getMinSeqlen(this->seqlen, this->numseqs);
	this->maxseqlen = _getMaxSeqlen(this->seqlen, this->numseqs);
}


//--------------------------------------------------------------------------
// Nullset
//--------------------------------------------------------------------------
Nullset::Nullset() {
	this->numIters = 0;

	if(DEBUG0) {
		//verifies constructor was called
		assert(this->maxseqlen == 0);
		assert(this->minseqlen == 0);
		assert(this->numseqs ==0);
		assert(this->seqlen == NULL);
		assert(this->seqs == NULL);
	}
}

Nullset::Nullset(const Seqset &src) : Seqset(src) {
	this->numIters = 0;

	if(DEBUG0) {
		assert(this->minseqlen == src.minseqlen);
		assert(this->maxseqlen == src.maxseqlen);
	}

}

Nullset::Nullset(const Nullset &src) : Seqset(src){
	this->numIters = src.numIters;
	if(DEBUG0) {
		assert(this->minseqlen == src.minseqlen);
		assert(this->maxseqlen == src.maxseqlen);
	}
}


Nullset::~Nullset() {
	//automatically calls ~Seqset
}

int Nullset::getNumIters() {
	return this->numIters;
}

int Nullset::incrementIters() {
	this->numIters++;
	return (this->numIters);
}


void Nullset::randomize() {
	//randomize by shuffling within the window starting at current location
	int wndsize = 25; 
	for(int i = 0; i < this->numseqs; i++) {
		for(int j = 0; j < this->seqlen[i]; j++) {
			int w = (wndsize < this->seqlen[i] - j ? wndsize : this->seqlen[i] - j); //min()
			int r = (int) (Random() * w);
			if(DEBUG0) {
				if(j+r >= this->seqlen[i]) {
					cerr<<"Error: the position "<<j+r<<" is longer than the sequence length "<<this->seqlen[i]<<endl;
					abort();
				}
			}
			//swap nucleotide
			int temp = this->seqs[i][j];
			this->seqs[i][j] = this->seqs[i][j+r];
			this->seqs[i][j+r] = temp;
		}
	}
	
}


//--------------------------------------------------------------------------
// Input 
//--------------------------------------------------------------------------
Input::Input() {
	abort();
}

Input::Input(string fastaFilename) {
	openFastaMold(fastaFilename);
	this->bgSeqset = seqset; //reference copy
	this->separateBgfsa = false;
}

Input::Input(string fastaFilename, string bgFastaFilename) {
	openFastaMold(fastaFilename);

	Dataset *dataset = openBackgroundData((char*) bgFastaFilename.c_str(), 1); //use revcompl
	this->bgSeqset = new Seqset(dataset->seqs, dataset->numseqs, dataset->seqlen); 
	this->bgSeqset->removeGaps();
	nilDataset(dataset);

	this->separateBgfsa = true;
}

bool Input::hasSeparateBgfsa(){
	return this->separateBgfsa;
}

void Input::openFastaMold(string fastaFilename) {
	//open file using old c-code from GibbsMarkov
	Dataset *dataset = openBackgroundData((char*)fastaFilename.c_str(), 0);

	//Gapped sequence-set
	this->gappedSeqset = new Seqset(dataset->seqs, dataset->numseqs, dataset->seqlen);

	//Nongap sequence-set
	this->seqset = new Seqset(*(this->gappedSeqset));
	this->seqset->removeGaps();

	//headers
	this->fastaHeaders.resize(dataset->numseqs);
	for(int i = 0; i < dataset->numseqs; i++) {
		string str(dataset->headers[i]);
		this->fastaHeaders[i].clear();

		if(DEBUG1) {
			cerr<<str<<endl;
		}

		//remove '>' characters
		for(int k = 0; k < (int)str.length(); k++) {
			if(str[k] != '>') {
				this->fastaHeaders[i].push_back(str[k]);
			}
		}
	}

	//cleanup
	nilDataset(dataset);
}

Input::~Input() {
	delete this->gappedSeqset;
	delete this->seqset;
	if(this->separateBgfsa) {
		delete this->bgSeqset;
	}
}

void Input::superimposeGaps(int *gapseq, int gapseqLen, int *nongapseq, int nongapseqLen, string &result) {
	//sanity check
	if(DEBUG0) {
		int count = 0;
		for(int i = 0; i < gapseqLen; i++) {
			if(gapseq[i] != GAP_CHAR) {
				count++;
			}
		}
		if(count != nongapseqLen) {
			cerr<<"The number of non-gap characters in gapseq and nongapseq does not match"<<endl;
			cerr<<"gapseq: "<<count<<"; nongapseq: "<<nongapseqLen<<endl;
			abort();
		}
		count = 0;
		for(int i = 0; i < nongapseqLen; i++) {
			assert(nongapseq[i] != GAP_CHAR);
		}
	}

	//allocate space
	result.resize(gapseqLen);

	//superimpose
	int nongap_pos = 0;
	for(int i = 0; i < gapseqLen; i++) {
		if(gapseq[i] == GAP_CHAR) {
			result[i] = 'N';
		}
		else {
			result[i] = numToChar(nongapseq[nongap_pos++]);
		}
	}
	if(nongap_pos != nongapseqLen) {
		cerr << "Fatal error: nongap_pos="<<nongap_pos<<" is not equal to nongapseqLen="<<nongapseqLen<<endl;
		abort();
	}
}


