#include "dataset.h"

#define print_error(str) for(fprintf(stderr,"%s\n",str); TRUE; exit(1))

static const int FASTA_TAG_STRLEN = 5;

static
FILE *openFile(char *filename)
{
	FILE *fptr;
	if((fptr = fopen(filename,"r")) == NULL) {
		fprintf(stderr,"Could not open file \"%s\"\n",filename);
		print_error("File does not exist!\n");
	}
	return fptr;
}


static
void printSeqs(Dataset *data) {
	for(int i = 0; i < data->numseqs; i++) {
		for(int j = 0; j < data->seqlen[i]; j++) {
			fprintf(stderr, "%c", numToChar(data->seqs[i][j]));
		}
		fprintf(stderr, "\n");
	}
}

//---------------------------------------------------------------------------------
// Converting to Reverse-complementary (double strand)
//---------------------------------------------------------------------------------
static
void augumentRevcompl2Dataset(Dataset *data) {
	if(!data->useRevcompl) {
		printf("Error: invalid usage of revcompl\n");
		exit(1);
	}

	//concatenate reverse complement to each strand
	for(int i = 0; i < data->numseqs; i++) {
		for(int j = 0; j < data->seqlen[i]; j++) {
			data->seqs[i][ 2*data->seqlen[i] - 1 - j ] = complChar(data->seqs[i][j]);
		}
	}

	//modify length of sequence (2 * current length)
	for(int i = 0; i < data->numseqs; i++) {
		data->seqlen[i] = data->seqlen[i] * 2;
	}

}

//---------------------------------------------------------------------------------
// Reverse-complementary helper functions
//---------------------------------------------------------------------------------

//Can be used for both reverse-complementary or single strand
//Double-strand-pos starts with position 1 or -1
int getConcat2DoubleStrandPos(Dataset* data, int seqind, int concatPos) {
	int dspos;
	if(!data->useRevcompl) {
		dspos = concatPos + 1;
	}
	else {
		//the true sequence length is a half of seqlen
		if(concatPos < data->seqlen[seqind]/2) {
			dspos = concatPos + 1;
		}
		else {
			dspos =  concatPos - data->seqlen[seqind];
		}
	}
	if(DEBUG0) {
		//check to see if numbers are within bound
		if(data->useRevcompl) {
			assert( abs(dspos) <= data->seqlen[seqind]/2);
			assert( 1 <= abs(dspos) );
		}
	}
	return dspos;
}


int getDoubleStrand2ConcatPos(Dataset *data, int seqind, int dspos) {
	int concatPos;
	if(dspos >= 1) {
		concatPos = dspos - 1;
	}
	else if(dspos <= -1){
		concatPos = data->seqlen[seqind] + dspos;
	}
	else {
		fprintf(stderr, "Error: double-strand position cannot be 0");
		exit(1);
	}
	if(DEBUG0) {
		assert(dspos == getConcat2DoubleStrandPos(data, seqind, concatPos));
	}
	if(concatPos < 0 || concatPos >= data->seqlen[seqind]) {
		fprintf(stderr, "Error: incorrect dspos input for getDoubleStrand2ConcatPos()");
		exit(1);
	}
	return concatPos;
}

//Can be used for single strand as well
bool isForwardStrand(Dataset* data, int seqind, int concatPos) {
	int ds = getConcat2DoubleStrandPos(data, seqind, concatPos);
	return ds > 0;
}


//If single strand is used, then it just returns concatPos.
//Else it returns the starting position of the k-mer on opposite strand (position span-1 of k-mer)
//Span must be greater 1. If the user wants the concatPos of the revcompl strand, set span=1.
int getConcatPosOfOppStrand(Dataset* data, int seqind, int concatPos, int span) {
	int newpos = -1;

	if(!data->useRevcompl) {
		newpos = concatPos;
	}
	else {
		//same for both forward and backward. I got this by drawing.
		newpos = data->seqlen[seqind] - concatPos - span;
	}

	if(DEBUG0) {
		assert(span > 0);
	}

	return newpos;
}


static
void determineValidSites(Dataset *data) {
	int minspan = data->minspan;
	int maxspan = data->maxspan;

	for(int s = 0; s < SPAN_CAPACITY; s++) {
		data->isValidSite[s] = NULL;
		data->numValidSites[s] = NULL;
	}

	for(int s = minspan; s <= maxspan; s++) {
		data->isValidSite[s] = (int**) malloc(sizeof(int*) * data->numseqs);
		for(int i = 0; i < data->numseqs; i++) {
			data->isValidSite[s][i] = (int*) malloc(sizeof(int) * data->seqlen[i]);
		}
		data->numValidSites[s] = (int*) calloc(data->numseqs, sizeof(int));
	}

	for(int s = minspan; s <=maxspan; s++) {
		for(int i = 0; i < data->numseqs;i++) {
			for(int j = 0; j < data->seqlen[i]; j++) {
				//Set invalid sites in between forward and backward strand if applicable
				if(!data->useRevcompl) {
					data->isValidSite[s][i][j] = ( j + s - 1 < data->seqlen[i] ? 1 : 0);
				}
				else {
					data->isValidSite[s][i][j] = (
						(j + s - 1 < data->seqlen[i]/2) ||  //forward
						(data->seqlen[i]/2 <= j && j + s - 1 < data->seqlen[i]) //backward
						? 1 : 0);
				}

				//sites that contain character N or W are invalid sites
				if( j + s - 1 < data->seqlen[i]) {
					for(int k = 0; k < s; k++) {
						if(data->isBadPos[i][j+k]) {
							data->isValidSite[s][i][j] = 0;
						}
					}
				}
			}
		}
	}

	//ensure each sequence has a valid site
	for(int s = minspan; s <=maxspan; s++) {
		for(int i = 0; i < data->numseqs;i++) {
			bool hasValidSite = false;
			for(int j = 0; j < data->seqlen[i]; j++) {
				if(data->isValidSite[s][i][j]) {
					hasValidSite = true;
				}
			}
			if(data->seqlen[i] < 1 || !hasValidSite) {
				fprintf(stderr, "ERROR: Sequence %d does not have a valid site for a motif\n", i);
				exit(1);
			}
		}
	}

	//count number of valid sites
	for(int s = minspan; s <=maxspan; s++) {
		for(int i = 0; i < data->numseqs;i++) {
			data->numValidSites[s][i] = 0;
			for(int j = 0; j < data->seqlen[i]; j++) {
				if(data->isValidSite[s][i][j]) {
					data->numValidSites[s][i]++;
				}
			}
		}
	}

	if(DEBUG0) {
		if(data->useRevcompl) {
			// reverse-complementary mode should have even sequence length
			for(int i = 0; i < data->numseqs; i++) {
				assert(data->seqlen[i] % 2 == 0);
			}
		}
	}

	if(DEBUG1) {
		for(int s = minspan; s <=maxspan; s++) {
			for(int i = 0; i < data->numseqs;i++) {
				fprintf(stderr, "Number of valid sites for span %d seq %4d: %5d\n", s, i, data->numValidSites[s][i]);
			}
		}
	}
}

//Test whether is concatPos is a valid site, i.e. it is not at the border between forward/backward strand
//Can be used for both reverse-complementary or single strand
bool isValidConcatPos(Dataset* data, int seqind, int concatPos, int span) {
	for(int m = 0; m < span; m++) {
		if(data->isBadPos[seqind][concatPos + m]) {
			return false;
		}
	}
	if(concatPos < 0) {
		return false;
	}
	if(!data->useRevcompl) {
		return(concatPos + span - 1< data->seqlen[seqind]);

	}
	else {
		return (concatPos + span - 1 < data->seqlen[seqind]/2) || //forward
			(data->seqlen[seqind]/2 <= concatPos && concatPos + span - 1 < data->seqlen[seqind]); //backward
	}
}

//This has less constraint then isValidConcatPos. This only looks at whether the dspos is within bound of the sequence
bool isBoundedDoubleStrandPos(Dataset *data, int seqind, int dspos) {
	if(data->useRevcompl) {
		return ( 1 <= abs(dspos) && abs(dspos) <= data->seqlen[seqind]/2 );
	}
	else {
		return ( 1 <= dspos && dspos <= data->seqlen[seqind]);
	}
}

//
////Can be used for both reverse-complementary or single strand
//bool isValidConcatPos(Dataset* data, int seqind, int concatPos, int span) {
//	return (!data->useRevcompl && concatPos + span - 1< data->seqlen[seqind])
//		|| (data->useRevcompl &&
//			(concatPos + span - 1 < data->seqlen[seqind]/2 || //forward
//			(data->seqlen[seqind]/2 <= concatPos && concatPos + span - 1 < data->seqlen[seqind]) //backward
//			));
//
//}


//---------------------------------------------------------------------------------
// Reading sequences from file input
//---------------------------------------------------------------------------------

//returns the length of the sequence
//This function can easily leads to buffer overflow, but I don't care about security now.
static
int readSingleSeq(FILE *fptr, Dataset *data, int seqInd) {
	int siteInd = 0;
	int num;
	char c;
	while((c = fgetc(fptr)) != '>' && c!= EOF) {
		if( c == '\n' || c== '\r') {
			continue;
		}
		if(siteInd >= data->seqlen[seqInd]) {
			print_error("ERROR: exceeds seqlen\n");
		}

		num = charToNum(c);
		data->seqs[seqInd][siteInd] = num;
		siteInd++;
	}
	ungetc(c, fptr);

	return siteInd;
}

static
void readSeqs(FILE *fptr, Dataset *data) {
	char c;
	int seqInd = 0;
	boolean inHeader = FALSE;
	while ((c = fgetc(fptr)) != EOF) {
		if( c == '>' ) {
			if(inHeader) {
				print_error("Error: There is extra > in header.\n");
			}
			inHeader = TRUE;
			if(seqInd >= data->numseqs) {
				print_error("ERROR in readSeqs: number of seqs exceed max.");
			}
		}
		else if ( inHeader && c == '\n') {
			inHeader = FALSE;
			int seqlen = readSingleSeq(fptr, data, seqInd);
			if(seqlen != data->seqlen[seqInd]) {
				print_error("ERROR: length of seqs does not match");
			}
			seqInd++;
		}
	}
	if(seqInd != data->numseqs) {
		print_error("ERROR in readSeqs: Inconsistent number of sequences.");
	}
}

static
void readHeaders(FILE *fptr, Dataset *data) {
	data->headers = (char**) malloc(sizeof(char*) * data->numseqs);
	for(int i = 0; i < data->numseqs; i++) {
		data->headers[i] = (char*) malloc(sizeof(char) * (FASTA_HEADER_CAPACITY+1));
	}

	char c;
	int seqInd = 0;
	int strLen = 0;
	boolean inHeader = FALSE;
	while ((c = fgetc(fptr)) != EOF) {
		if( c == '\r') {
			continue; //skip over these annoying characters
		}
		if( c == '>' ) {
			if(inHeader) {
				print_error("Error: There is extra > in header.\n");
			}
			inHeader = TRUE;
			if(seqInd >= data->numseqs) {
				print_error("ERROR in readSeqs: number of seqs exceed max.");
			}

			strLen = 0;
			data->headers[seqInd][strLen++] = c;
		}
		else if ( inHeader ){
			if( c == '\n') {
				inHeader = FALSE;
				data->headers[seqInd][strLen++] = '\0'; //null-termination character
				seqInd++;
			}
			else {
				if(strLen < FASTA_HEADER_CAPACITY) { //one extra space for null character
					data->headers[seqInd][strLen++] = c;
				}
			}
		}
	}
	if(seqInd != data->numseqs) {
		print_error("ERROR in readHeaders: Inconsistent number of sequences.");
	}
}

static
int getNumOfSeqs(FILE *fptr) {
	char c;
	int seqInd = 0;
	while ((c = fgetc(fptr)) != EOF) {
		if( c == '>' ) {
			seqInd++;
		}
	}
	return seqInd;
}

static
void getLenOfSeqs(FILE *fptr, Dataset *data) {
	char c;
	int seqInd = -1;
	boolean inHeader = FALSE;
	boolean inSeq = FALSE;
	while ((c = fgetc(fptr)) != EOF) {
		if( c == '>' ) {
			inHeader = TRUE;
			inSeq = FALSE;
			seqInd++;
			if(seqInd >= data->numseqs) {
				print_error("ERROR in getLenOfSeqs: number of seqs exceed max.");
			}
			data->seqlen[seqInd] = 0;
		}
		else if ( inHeader && c == '\n') {
			inHeader = FALSE;
			inSeq = TRUE;
		}
		else if (inSeq && isChar(c)) {
			data->seqlen[seqInd]++;
		}
	}
}

//should be used after augmentation of revcompl strand
static
void constructBadPos(Dataset *data) {
	data->isBadPos = (bool**) malloc(data->numseqs * sizeof(bool*));
	for(int i = 0; i < data->numseqs; i++) {
		data->isBadPos[i] = (bool*) malloc(data->seqlen[i] * sizeof(bool));
	}

	for(int i = 0; i < data->numseqs; i++) {
		for(int j = 0; j < data->seqlen[i]; j++) {
			data->isBadPos[i][j] = !isNucleotide(data->seqs[i][j]);
		}
	}
}

//

static
Dataset* _readFasta(char *filename, bool useRevcompl) {
	Dataset *data;
	FILE *fptr;

	if(DEBUG1) {
		fprintf(stderr, "FASTA open: %s\n", filename);
	}

	fptr = openFile(filename);
	data = (Dataset*) malloc(sizeof(Dataset));
	data->useRevcompl = useRevcompl;

	//set isvalid site and badpos to NULL
	for(int i = 0; i < SPAN_CAPACITY; i++) {
		data->isValidSite[i] = NULL;
	}
	data->isBadPos = NULL;

	//first character is '>', which is a FASTA format style
	if(fgetc(fptr) != '>') {
		fprintf(stderr, "Error: %s does not start with >\n", filename);
		exit(1);
	}

	//number of sequence
	rewind(fptr);
	data->numseqs = getNumOfSeqs(fptr);

	data->seqlen = (int*) malloc( data->numseqs * sizeof(int));
	data->seqs = (int**) malloc(data->numseqs * sizeof(int*));

	//len per seq
	rewind(fptr);
	getLenOfSeqs(fptr, data);

	for(int i = 0; i < data->numseqs; i++) {
		int len;

		//reverse-complementary is twice as long after concatenating
		if(data->useRevcompl) {
			len = data->seqlen[i] * 2;
		}
		else {
			len = data->seqlen[i];
		}
		data->seqs[i] = (int*) malloc(len * sizeof(int));
	}

	//read sequences (core of the function)
	rewind(fptr);
	readSeqs(fptr, data);
	rewind(fptr);
	readHeaders(fptr, data);
	fclose(fptr);

	if(DEBUG1) {
		fprintf(stderr, "FASTA closed\n");
	}

	//double-strand should be transparent as a single strand after this point
	if(data->useRevcompl) {
		//concatenate reverse complement to each strand
		//modify length of sequence (2 * current length)
		augumentRevcompl2Dataset(data);
	}

	if(DEBUG2) {
		printSeqs(data);
	}

	//Find avg/max/min of sequence length
	data->maxseqlen = 0;
	for(int i = 0; i < data->numseqs; i++) {
		if(data->maxseqlen < data->seqlen[i]) {
			data->maxseqlen = data->seqlen[i];
		}
	}
	data->minseqlen = INT_MAX;
	for(int i = 0; i < data->numseqs; i++) {
		if(data->minseqlen > data->seqlen[i]) {
			data->minseqlen = data->seqlen[i];
		}
	}
	data->avgseqlen = 0.0;
	for(int i = 0; i < data->numseqs; i++) {
		data->avgseqlen += data->seqlen[i];
	}
	data->avgseqlen /= data->numseqs;

	//determine count
	int count[NUMCHARS];
	for(int a = 0; a < NUMCHARS; a++) {
		count[a] = 0;
	}
	for(int i = 0; i < data->numseqs; i++) {
		for(int j = 0; j < data->seqlen[i]; j++) {
			count[ data->seqs[i][j] ]++;
		}
	}
	data->totalCount = 0;
	for(int a = 0; a < NUMALPHAS; a++) {
		data->count[a] = count[a];
		data->totalCount += count[a];
	}
	for(int a = 0; a < NUMALPHAS; a++) {
		data->ntFreq[a] = ((double)data->count[a]) / data->totalCount;
	}

	return data;
}
//---------------------------------------------------------------------------------
// Main functions
//---------------------------------------------------------------------------------

Dataset *openBackgroundData(char *filename, bool useRevcompl) {
	Dataset *data = _readFasta(filename, useRevcompl);

	data->minspan = 0; //hack for minspan and maxspan, i.e. minspan > maxspan
	data->maxspan = -1;

	return data;
}

Dataset *openDataset(char *filename, bool useRevcompl, int minspan, int maxspan) {
	Dataset *data = _readFasta(filename, useRevcompl);

	data->minspan = minspan;
	data->maxspan = maxspan;

	if(DEBUG1) {
		fprintf(stderr, "Constructing badPos...\n");
	}
	constructBadPos(data); //should be used after augument of revcompl

	if(DEBUG1) {
		fprintf(stderr, "Constructing validSites...\n");
	}
	determineValidSites(data);


	if(DEBUG1) {
		fprintf(stderr, "Complete reading dataset: %s\n", filename);
	}

	return data;
}


void nilDataset(Dataset *data) {
	free(data->seqlen);

	for(int i = 0; i< data->numseqs; i++) {
		free(data->seqs[i]);
		free(data->headers[i]);

		if(data->isBadPos != NULL) {
			free(data->isBadPos[i]);
		}
	}
	free(data->seqs);
	free(data->headers);
	if(data->isBadPos != NULL) {
		free(data->isBadPos);
	}

	for(int s = data->minspan; s <= data->maxspan; s++) {
		if(data->isValidSite[s] != NULL) {
			for(int i = 0; i < data->numseqs; i++) {
				free(data->isValidSite[s][i]);
			}
			free(data->isValidSite[s]);
		}
		if(data->numValidSites[s] != NULL) {
			free(data->numValidSites[s]);
		}
	}
	free(data);
}

