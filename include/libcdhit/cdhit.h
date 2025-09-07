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
void cluster_sequences(
		std::vector<Sequence_new>& sequences,
		std::vector<int>& parent,
		int kmer_size = 5,
		double tau = 0.36,
		int nthreads = 1
		);
void cluster_sequences_st(
		std::vector<Sequence_new>& sequences,
		std::vector<int>& parent,
		int kmer_size = 5,
		double tau = 0.36
		);

void cluster_sequences_st_less10(
		std::vector<Sequence_new>& seqs,
		std::vector<int>& parent,
		int kmer_size = 5,
		double tau = 0.36
		);

#endif
