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
class ClusterBuffer{
public:
	ClusterBuffer(int seqs_no, int kmer_size, int max_len){
		max_seq_len = max_len;
		int table_size = 1;
		for (int i = 0; i < kmer_size; ++i) 
			table_size *= MAX_UAA;

		word_table.reserve(table_size);
		word_table.resize(table_size);

		//for(int i = 0; i < table_size; i++)
		//	word_table[i].reserve(64);

		word_encodes.reserve(max_seq_len);
		word_encodes_no.reserve(max_seq_len);

		A.reserve(seqs_no);
		counts.reserve(seqs_no);
		
		visited.reserve(1<<10);
	}

	void clear(){
		for(int i = 0; i < word_table.size(); i++)
			word_table[i].clear();
		A.clear();
		//std::fill(counts.begin(), counts.end(), 0);
		counts.clear();
	}

	vector<vector<pair<int, int>>> word_table;//(table_size)
	vector<int> word_encodes;//(max_seq_len);
	vector<int> word_encodes_no;//(max_seq_len);
	std::vector<int> A;//(seqs.size());
	std::vector<int> counts;//(N, 0)
	std::vector<int> visited; //buffering
	std::vector<std::pair<int,int>> out_pairs; //buffering             
	
	
	int max_seq_len;	
};
void cluster_sequences(
		std::vector<Sequence_new>& sequences,
		std::vector<int>& parent,
		int kmer_size = 5,
		double tau = 0.36,
		int nthreads = 1
		);
void cluster_sequences_st(
		ClusterBuffer & buff,
		std::vector<Sequence_new>& sequences,
		std::vector<int>& parent,
		int kmer_size = 5,
		double tau = 0.36
		);


#endif
