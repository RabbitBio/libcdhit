#ifndef CDHIT_H
#define CDHIT_H
#include "cdhit-common.h"
#include "input_sequence.h"

class cluster{
public:
//	vector<int> parent;
	cluster() {
	}

	void cdhit_cluster(std::vector<Sequence_new>& seq, std::vector<int>& parent);
	void cdhit_cluster(std::vector<Sequence_new>& seq, std::vector<int>& parent, int neededThread);
};
#endif
