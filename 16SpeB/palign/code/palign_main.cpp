#include "stdinc.h"
#include "Input.h"
#include "nwalign.h"
#include "DisplayResults.h"
#include "random.h"

using namespace std;

enum PairType { ALL_PAIR, NEXT_PAIR, RAND_PAIR};

static
void printHelp() {
	cout << "Pairwise global alignment" << endl << endl
		<< "Usage: <program name> <seqset-FASTA> [OPTIONS]" << endl <<endl
		<< "-s <UINT>" <<endl
		<< "-quiet             Does not display alignment" <<endl
		<< "-print-fsa         Print FASTA in STDERR" <<endl
		<< endl
		<< "-all-pair          All possible pairs (n-choose-2 pairs)" <<endl
		<< "-next-pair         Every next pair (n/2 pairs)" <<endl
		<< "-rand-pair <INT>   Sample specified number of pairs "<<endl
		<<endl;
	exit(1);
}

static
void alignHelper(
		int seqind1, 
		int seqind2, 
		bool printFsa, 
		bool quietOut, 
		AlignPair *pair, 
		NWAlignParams *nwparams, 
		Input *input
		) {
	int *seq1 = input->seqset->seqs[seqind1];
	int *seq2 = input->seqset->seqs[seqind2];

	int seqlen1 = input->seqset->seqlen[seqind1];
	int seqlen2 = input->seqset->seqlen[seqind2];

	nwalign(nwparams, seq1, seqlen1, seq2, seqlen2, pair);
	double pidOverNongap = computePidOverNongap(pair->align1, pair->align2, pair->len);
	double pidOverAlignlen = computePidOverAlignlen(pair->align1, pair->align2, pair->len);

	//display
	if(printFsa) {
		cerr<<">"<<input->fastaHeaders[seqind1]
			<<"; PID1-over-non-gap="<<pidOverNongap<<"; PID1-over-alignlen="<<pidOverAlignlen<<endl;
		Results::displaySeq(cerr, pair->align1, pair->len);
		cerr<<">"<<input->fastaHeaders[seqind2]
			<<"; PID2-over-non-gap="<<pidOverNongap<<"; PID2-over-alignlen="<<pidOverAlignlen<<endl;
		Results::displaySeq(cerr, pair->align2, pair->len);
		cerr<<endl;
	}

	cout<<">"<<input->fastaHeaders[seqind1]<<endl;
	cout<<">"<<input->fastaHeaders[seqind2]<<endl;

	if(!quietOut) {
		cout<<"Original:"<<endl;
		Results::displaySeq(cout, seq1, seqlen1);
		Results::displaySeq(cout, seq2, seqlen2);
		cout<<endl;

		cout<<"Alignment:"<<endl;
		Results::displaySeq(cout, pair->align1, pair->len);
		Results::displaySeq(cout, pair->align2, pair->len);
		cout<<endl;
	}

	cout<<"PID over non-gap: "<< pidOverNongap <<endl;
	cout<<"PID over alignment-length: "<< pidOverAlignlen<<endl;
	cout<<endl;

	cout<<"==================================================================="<<endl;
	cout<<endl;


}

int main(int argc, char** argv) {
	if(DEBUG0) {
		string str = "WARNING: running under DEBUG mode\n\n";
		cout<<str;
		cerr<<str;
	}
	if(DEBUG1) {
		string str = "WARNING: running under VERBOSE mode\n\n";
		cout<<str;
		cerr<<str;
	}
	printf("Compiled on " __DATE__ " " __TIME__ "\n");
	printf("\n");
	printf("ChangeLog\n");
	printf("\n");
	for(int i = 0; i < argc; i++) {
		cout<<argv[i]<<" ";
	}
	cout<<endl<<endl;

	if(argc == 1) {
		printHelp();
	}

	//set variables from argv
	//argv[1] FASTA file
	string fastaFilename(argv[1]);
	PairType pairMode = NEXT_PAIR;
	bool quietOut = false;
	bool printFsa = false;
	int numRandPairs = 0;
    unsigned int randomSeed = (unsigned int)time(NULL);

	int i = 2;
	while(i < argc) {
		if (!strcmp(argv[i],"-all-pair")) {
			pairMode = ALL_PAIR;
		}
		else if (!strcmp(argv[i],"-next-pair")) {
			pairMode = NEXT_PAIR;
		}
		else if (!strcmp(argv[i],"-rand-pair")) {
			pairMode = RAND_PAIR;
			i++;
			int err = sscanf(argv[i], "%d", &(numRandPairs));
			if(err<1) printHelp();
		}
		else if (!strcmp(argv[i],"-s")) {
			i++;
			int err = sscanf(argv[i], "%d", &(randomSeed));
			if(err<1) printHelp();
		}
		else if (!strcmp(argv[i],"-quiet")) {
			quietOut = true;
		}
		else if (!strcmp(argv[i],"-print-fsa")) {
			printFsa = true;
		}
		else {
			printf("Unknown command: %s\n", argv[i]);
			printHelp();
		}
		i++;
	}

	clock_t startClock = clock();
    sRandom(randomSeed);

	Input *input = new Input(fastaFilename);
	int seq_maxlen = input->seqset->maxseqlen;
    cout<< "Random seed: " << randomSeed << endl;
    cout<< "Number of sequences: "<<input->seqset->numseqs <<endl;

	if(pairMode == NEXT_PAIR && input->seqset->numseqs % 2 != 0) {
		cerr<<"Error: FASTA file should have even number of sequences in -next-pair mode."<<endl;
		exit(1);
	}
 
    //matlab has 5,-4,-8, 
    //blastn has 1, -2, -5, -2
    int match = 1;
    int mismatch = -2;
    int gapopen = -5;
    int gapext = -2;
	cout<<"match "<<match<<endl;
	cout<<"mismatch "<<mismatch<<endl;
	cout<<"gapopen "<<gapopen<<endl;
	cout<<"gapext "<<gapext<<endl;
	cout<<endl;

	NWAlignParams* nwparams = constructNWAlignParams(match, mismatch, gapopen, gapext, seq_maxlen);
	AlignPair *pair = constructAlignPair(seq_maxlen, seq_maxlen);

	int pairsCount = 0;

	if(pairMode == NEXT_PAIR) {
		for(int i = 0; i < input->seqset->numseqs; i+=2) {
			int seqind1 = i;
			int seqind2 = i+1;
			alignHelper(seqind1, seqind2, printFsa, quietOut, pair, nwparams, input);
			pairsCount++;
		}
	}
	else if(pairMode == ALL_PAIR) {
		for(int i = 0; i < input->seqset->numseqs; i++) {
			for(int j = i+1; j < input->seqset->numseqs; j++) {
				int seqind1 = i;
				int seqind2 = j;
				alignHelper(seqind1, seqind2, printFsa, quietOut, pair, nwparams, input);
				pairsCount++;
			}
		}
	}
	else if(pairMode == RAND_PAIR) {
		for(int i = 0; i < numRandPairs; i++) {
			int seqind1;
			int seqind2;
			do {
				seqind1 = (int) (Random() * input->seqset->numseqs);
				seqind2 = (int) (Random() * input->seqset->numseqs);
			}while(seqind1 == seqind2);
			alignHelper(seqind1, seqind2, printFsa, quietOut, pair, nwparams, input);
			pairsCount++;
		}
	}
	else {
		cerr<<"Invalid pairMode"<<endl;
		exit(1);
	}

	cout<<"Number of pairs aligned: "<<pairsCount<<endl;
	double elapsed = ((double) ( clock() - startClock )) / CLOCKS_PER_SEC;
	printf("Total elapsed CPU time (in seconds): %.2lf\n",elapsed );

	nilAlignPair(pair);
	nilNWAlignParams(nwparams);
	delete input;
}
