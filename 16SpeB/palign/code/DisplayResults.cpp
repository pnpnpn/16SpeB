#include "DisplayResults.h"

using namespace std;

Results::Results() {
}

Results::Results(Params *params) {
	this->params = params;
}

Results::~Results() {
}

void Results::displaySeq(ostream &out, int *seq, int len) {
	for(int i = 0; i < len; i++ ) {
		if(seq[i] == GAP_CHAR) {
			out << '-';
		}
		else {
			out << numToChar(seq[i]);
		}
	}
	out << endl;
}


void Results::_displayNullset(Nullset *nullset, int nullsetIndex, Params *params) {
	//For display:
	//nullset index starts at 1
	//sequence index starts at 1
	string ntseq;
	Input *input = params->input;

	//create nullsetName
	ostringstream oss;
	char buffer[50];
	sprintf(buffer, "null-%07d", nullsetIndex+1);
	oss << buffer;
	string nullsetName = oss.str();
	
	//display sequences
	cout<<"#BEGIN " << nullsetName <<endl;
	for(int i = 0; i < nullset->numseqs; i++) {
		printf(">%s | seqind=%04d | header=%s\n", nullsetName.c_str(), i, input->fastaHeaders[i].c_str());
		params->input->superimposeGaps(input->gappedSeqset->seqs[i], input->gappedSeqset->seqlen[i], 
			nullset->seqs[i], nullset->seqlen[i], ntseq);
		cout << ntseq <<endl;
	}
	cout<<"#END " << nullsetName <<endl;
	cout<<endl;
}

void Results::display(Nullset *nullset, int nullsetIndex) {
	//cout <<endl;
	//cout << "########################################################"<<endl<<endl;;
	//for(int i = 0; i < params->getNumNullsets(); i++) {
	//	_displayNullset(params->getNullset(i), i, params);
	//	cout << "########################################################"<<endl<<endl;;
	//}
	_displayNullset(nullset, nullsetIndex, params);
	cout << "########################################################"<<endl<<endl;;
	
}
