#include "symbols.h"



//only modify if the character is between 'a' and 'z'
static
char lowerToUppercase(char c) {
	if('a' <= c && c <= 'z') {
		return c + 'A' - 'a';
	}
	return c;
}

int charToNum(char alpha) {
	switch (lowerToUppercase(alpha)) {
		case '-': 
		case 'X': 
		case 'R': 
		case 'Y': 
		case 'M': 
		case 'K': 
		case 'S': 
		case 'W': 
		case 'B': 
		case 'D': 
		case 'H': 
		case 'V': 
		case 'N': return GAP_CHAR;

		case 'A': return 0; 
		case 'C': return 1; 
		case 'G': return 2; 
		case 'T': return 3; 
	}
	fprintf(stderr, "ERROR: charToNum has incorrect input: %c\n", alpha);
	exit(1);
}

char numToChar(int num) {
	switch (num) {
		case 0: return 'A'; 
		case 1: return 'C'; 
		case 2: return 'G'; 
		case 3: return 'T'; 

		case GAP_CHAR: return 'X'; 
	}
	fprintf(stderr, "ERROR: numToChar has incorrect input\n");
	exit(1);
}


//nucleotide is synonymous to alpha
bool isNucleotide(int num) {
	if(num == GAP_CHAR) {
		return false;
	}
	else if(num >=0 && num < NUMALPHAS) {
		return true;
	}
	else {
		fprintf(stderr, "ERROR: isNotNt has incorrect input: %d\n", num);
		exit(1);
	}
}

boolean isChar(char c) {
	switch(lowerToUppercase(c)) {
		case '-': 
		case 'X': 
		case 'R': 
		case 'Y': 
		case 'M': 
		case 'K': 
		case 'S': 
		case 'W': 
		case 'B': 
		case 'D': 
		case 'H': 
		case 'V': 
		case 'N': 

		case 'A': 
		case 'C':  
		case 'G':  
		case 'T': return true; 
		default: return false;
	}
}


int complChar(int num) {
	switch (num) {
		case 0: return 3; //A to T
		case 1: return 2; //C to G
		case 2: return 1; //G to C
		case 3: return 0; //T to A

		case GAP_CHAR: return GAP_CHAR;
	}
	fprintf(stderr, "ERROR: compl has incorrect input characeter %d\n", num);
	exit(1);
}

