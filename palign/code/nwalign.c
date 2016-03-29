#include "nwalign.h"
#include "symbols.h"

enum DirectionType { DIR_ERR, DIR_M, DIR_IX, DIR_IY};

static
void _display_seq(FILE *fptr, int *seq, int len) {
	for(int i = 0; i < len; i++ ) {
		if(seq[i] == GAP_CHAR) {
			fprintf(fptr, "-");;
		}
		else {
			fprintf(fptr, "%c", numToChar(seq[i]));
		}
	}
	fprintf(fptr, "\n");
}

static
void inplace_reverse_seq(int *seq, int len) {
	int left = 0;
	int right = len-1;
	while(left < right) {
		int temp = seq[left];
		seq[left] = seq[right];
		seq[right] = temp;

		left++;
		right--;
	}
}

static
void traceback_align(NWAlignParams *params, int *seq1, int len1, int *seq2, int len2, AlignPair *result) {
	int *align1 = result->align1;
	int *align2 = result->align2;

	int align_pos = 0;
	int i = len1; //remember that tbmatrix is size len1+1 by len2+1
	int j = len2;

	//starting direction/matrix
	int dirtyp;
	if(params->Ix[len1][len2] >= params->Iy[len1][len2] && params->Ix[len1][len2] >= params->dpm[len1][len2]) {
		dirtyp = DIR_IX;
	}
	else if(params->Iy[len1][len2] >= params->Ix[len1][len2] && params->Iy[len1][len2] >= params->dpm[len1][len2]) {
		dirtyp = DIR_IY;
	}
	else {
		dirtyp = DIR_M;
	}


	//create the alignment in reverse
	while( i > 0 || j > 0 ) {
		if(DEBUG0) {
			assert(align_pos < result->capacity);
		}

		if(dirtyp == DIR_IX) {
			align1[align_pos] = seq1[i-1];
			align2[align_pos] = GAP_CHAR;
			dirtyp = params->tb_Ix[i][j];
			i--;
		}
		else if(dirtyp == DIR_IY) {
			align1[align_pos] = GAP_CHAR;
			align2[align_pos] = seq2[j-1];
			dirtyp = params->tb_Iy[i][j];
			j--;
		}
		else if(dirtyp == DIR_M) {
			align1[align_pos] = seq1[i-1];
			align2[align_pos] = seq2[j-1];
			dirtyp = params->tb_dpm[i][j];
			i--;
			j--;
		}
		else {
			if(DEBUG0) {
				fprintf(stderr, "Error: DIR_ERR found at position (%d,%d) and align_pos %d\n", i,j, align_pos);
			}
			abort();
		}

		align_pos++;
	}

	//reverse the alignment
	inplace_reverse_seq(align1, align_pos);
	inplace_reverse_seq(align2, align_pos);

	result->len = align_pos;

	if(DEBUG0) {
		assert(align_pos <= len1 + len2);

		//check if all alphabets are present without gaps
		int pos = 0;
		for(int i = 0; i < result->len; i++) {
			if(result->align1[i] != GAP_CHAR) {
				if(seq1[pos] != result->align1[i]) {
					fprintf(stderr, "seq1 does not match align1 at position (seq1, %d)\n", pos);
					_display_seq(stderr, seq1, len1);
					_display_seq(stderr, result->align1, result->len);
					abort();
				}
				pos++;
			}
		}
		pos = 0;
		for(int i = 0; i < result->len; i++) {
			if(result->align2[i] != GAP_CHAR) {
				if(seq2[pos] != result->align2[i]) {
					fprintf(stderr, "seq2 does not match align2 at position (seq2, %d)\n", pos);
					_display_seq(stderr, seq2, len2);
					_display_seq(stderr, result->align2, result->len);
					abort();
				}
				pos++;
			}
		}
	}
}

static void _print_intmatrix(double **mat, int dim1, int dim2) {
	FILE *fptr = stderr;
	for(int i = 0; i < dim1; i++) {
		for(int j = 0; j < dim2; j++) {
			fprintf(fptr, "%8.0lf", mat[i][j]);
		}
		fprintf(fptr, "\n");
	}
}


//See Durbin p. 29, equation (2.16)
void nwalign(NWAlignParams *params, int *seq1, int len1, int *seq2, int len2, AlignPair* result){
	if(DEBUG0) {
		if(len1+1 > params->matrix_capacity || len2+1 > params->matrix_capacity) {
			fprintf(stderr, "DP matrix out of range.\n");
			abort();
		}
	}

	double **dpm = params->dpm;
	double **Ix = params->Ix;
	double **Iy = params->Iy;

	int **tb_dpm = params->tb_dpm; //traceback
	int **tb_Ix = params->tb_Ix; //traceback
	int **tb_Iy = params->tb_Iy; //traceback

	int match = params->match;
	int mismatch = params->mismatch;
	int gapopen = params->gapopen;
	int gapext = params->gapext;

	if(DEBUG1) {
		cerr<<"match "<<match<<endl;
		cerr<<"mismatch "<<mismatch<<endl;
		cerr<<"gapopen "<<gapopen<<endl;
		cerr<<"gapext "<<gapext<<endl;
	}

	//initialize
	dpm[0][0] = 0; //M(i,j) is the best score up to (i,j) given that x_i is aligned to y_i
	Ix[0][0] = -INFINITY; //so that gap open must start from dpm
	Iy[0][0] = -INFINITY;
	tb_dpm[0][0] = DIR_ERR;
	tb_Ix[0][0] = DIR_ERR;
	tb_Iy[0][0] = DIR_ERR;
	for(int i = 1; i < len1+1; i++) {
		dpm[i][0] = -INFINITY;
		Ix[i][0] = gapopen + (i-1) * gapext;
		Iy[i][0] = -INFINITY;

		tb_dpm[i][0] = DIR_ERR;
		tb_Ix[i][0] = (i == 1 ? DIR_M : DIR_IX);
		tb_Iy[i][0] = DIR_ERR;
	}
	for(int j = 1; j < len2+1; j++) {
		dpm[0][j] = -INFINITY;
		Ix[0][j] = -INFINITY;
		Iy[0][j] = gapopen + (j-1) * gapext;

		tb_dpm[0][j] = DIR_ERR;
		tb_Ix[0][j] = DIR_ERR;
		tb_Iy[0][j] = (j == 1 ? DIR_M : DIR_IY);
	}

	//set DP matrix
	for(int i = 1; i < len1 + 1; i++) {
		for(int j = 1; j < len2 +1; j++) {
			double m_val, ix_val, iy_val;
			int s = (seq1[i-1] == seq2[j-1] ? match : mismatch);

			//Setting dpm 
			m_val = dpm[i-1][j-1]+s;
			ix_val = Ix[i-1][j-1]+s;
			iy_val = Iy[i-1][j-1]+s;
		
			if(ix_val >= m_val && ix_val >= iy_val) {
				dpm[i][j] = ix_val;
				tb_dpm[i][j] = DIR_IX;
			}
			else if(iy_val >= m_val && iy_val >= ix_val) {
				dpm[i][j] = iy_val;
				tb_dpm[i][j] = DIR_IY;
			}
			else {
				dpm[i][j] = m_val;
				tb_dpm[i][j] = DIR_M;
			}

			//Setting Ix
			m_val = dpm[i-1][j] + gapopen;
			ix_val = Ix[i-1][j] + gapext;

			if(ix_val >= m_val) {
				Ix[i][j] = ix_val;
				tb_Ix[i][j] = DIR_IX;
			}
			else {
				Ix[i][j] = m_val;
				tb_Ix[i][j] = DIR_M;
			}

			//Setting Iy
			m_val = dpm[i][j-1] + gapopen;
			iy_val = Iy[i][j-1] + gapext;

			if(iy_val >= m_val) {
				Iy[i][j] = iy_val;
				tb_Iy[i][j] = DIR_IY;
			}
			else {
				Iy[i][j] = m_val;
				tb_Iy[i][j] = DIR_M;
			}
		}
	}

	traceback_align(params, seq1, len1, seq2, len2, result);

	if(DEBUG1) {
		//fprintf(stderr, "DP matrix\n");
		//_print_intmatrix(dpm, len1 + 1, len2+1) ;
		//fprintf(stderr, "\n");
		//fprintf(stderr, "Traceback matrix\n");
		//_print_intmatrix(tbm, len1 + 1, len2+1) ;
		//fprintf(stderr, "\n");

		fprintf(stderr, "dpm_score=%.1lf\n", dpm[len1][len2]);
		fprintf(stderr, "Ix_score=%.1lf\n", Ix[len1][len2]);
		fprintf(stderr, "Iy_score=%.1lf\n", Iy[len1][len2]);
		_display_seq(stderr, result->align1, result->len);
		_display_seq(stderr, result->align2, result->len);
		fprintf(stderr, "\n");
	}
}

double computePidOverAlignlen(int *align1, int *align2, int len) {
	int ident = 0;
	for(int i = 0; i < len; i++) {
		if(align1[i] != GAP_CHAR && align2[i] != GAP_CHAR) {
			if(align1[i] == align2[i]) {
				ident++;
			}
		}
	}
	return ((double)ident) / len;
}


double computePidOverNongap(int *align1, int *align2, int len) {
	int ident = 0;
	int nongap = 0;
	for(int i = 0; i < len; i++) {
		if(align1[i] != GAP_CHAR && align2[i] != GAP_CHAR) {
			nongap++;
			if(align1[i] == align2[i]) {
				ident++;
			}
		}
	}
	if(nongap == 0) {
		return 0; //because ident must be 0 as well
	}
	else {
		return ((double)ident)/nongap;
	}

}

//-----------------------------------------------------------------------------------------
// Constructors/Destructors
//-----------------------------------------------------------------------------------------

AlignPair* constructAlignPair(int len1, int len2) {
	AlignPair *pair = (AlignPair*) malloc(sizeof(AlignPair));
	pair->len = 0;
	pair->capacity = len1+len2;
	pair->align1 = (int*) malloc(sizeof(int) * pair->capacity);
	pair->align2 = (int*) malloc(sizeof(int) * pair->capacity);

	if(pair->align1 <= 0 || pair->align2 <=0 ) {
		fprintf(stderr, "Out of memory at constructAlignPair()\n");
		abort();
	}

	return pair;
}

void nilAlignPair(AlignPair *pair) {
	free(pair->align1);
	free(pair->align2);
	free(pair);
}

NWAlignParams* constructNWAlignParams(int match, int mismatch, int gapopen, int gapext, int seq_maxlen) {
	NWAlignParams *params = (NWAlignParams*) malloc(sizeof(NWAlignParams));
	params->gapopen = gapopen;
	params->gapext = gapext;
	params->match = match;
	params->mismatch = mismatch;
	params->matrix_capacity = seq_maxlen + 1;

	params->dpm = (double**) malloc(sizeof(double*) * params->matrix_capacity);
	for(int i = 0; i < params->matrix_capacity; i++) {
		params->dpm[i] = (double*) malloc(sizeof(double) * params->matrix_capacity);
	}
	params->Ix= (double**) malloc(sizeof(double*) * params->matrix_capacity);
	for(int i = 0; i < params->matrix_capacity; i++) {
		params->Ix[i] = (double*) malloc(sizeof(double) * params->matrix_capacity);
	}
	params->Iy = (double**) malloc(sizeof(double*) * params->matrix_capacity);
	for(int i = 0; i < params->matrix_capacity; i++) {
		params->Iy[i] = (double*) malloc(sizeof(double) * params->matrix_capacity);
	}

	params->tb_dpm = (int**) malloc(sizeof(int*) * params->matrix_capacity);
	for(int i = 0; i < params->matrix_capacity; i++) {
		params->tb_dpm[i] = (int*) malloc(sizeof(int) * params->matrix_capacity);
	}
	params->tb_Ix = (int**) malloc(sizeof(int*) * params->matrix_capacity);
	for(int i = 0; i < params->matrix_capacity; i++) {
		params->tb_Ix[i] = (int*) malloc(sizeof(int) * params->matrix_capacity);
	}
	params->tb_Iy = (int**) malloc(sizeof(int*) * params->matrix_capacity);
	for(int i = 0; i < params->matrix_capacity; i++) {
		params->tb_Iy[i] = (int*) malloc(sizeof(int) * params->matrix_capacity);
	}

	return params;
}

void nilNWAlignParams(NWAlignParams *params) {
	for(int i = 0; i < params->matrix_capacity; i++) {
		free(params->dpm[i]);
		free(params->Ix[i]);
		free(params->Iy[i]);
		free(params->tb_dpm[i]);
		free(params->tb_Ix[i]);
		free(params->tb_Iy[i]);
	}
	free(params->dpm);
	free(params->Ix);
	free(params->Iy);
	free(params->tb_dpm);
	free(params->tb_Ix);
	free(params->tb_Iy);

	free(params);
}


