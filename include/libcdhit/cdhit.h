#ifndef CDHIT_H
#define CDHIT_H
#include "cdhit-common.h"
#include "input_sequence.h"
#include "robin_hood.h"

class cluster{
public:
//	vector<int> parent;
	cluster() {
	}

	void cdhit_cluster(std::vector<Sequence_new>& seq, std::vector<int>& parent);
	void cdhit_cluster(std::vector<Sequence_new>& seq, std::vector<int>& parent, int neededThread);
};

// 按需放到合适的头文件
struct ClusterWS {
	// 稀疏倒排表：token -> postings (seq_id, count)
	robin_hood::unordered_map<int, std::vector<std::pair<int,int>>> word_table;

	// 复用缓冲
	std::vector<int> word_encodes;
	std::vector<int> word_encodes_no;
	std::vector<int> counts;                 // N 大小
	std::vector<int> visited;                // 稀疏触达列表
	std::vector<std::pair<int,int>> out_pairs;
	std::vector<int> A;                      // |k-mers| for each sequence

	// 预留容量（可根据业务微调）
	//void reserve_for(int N, int max_seq_len, size_t expected_tokens_per_seq = 0) {
	//	// 倒排表：保留哈希桶；不清 map，本轮开始前把 row 清空即可
	//	if (expected_tokens_per_seq) {
	//		// 如果能估计“每条序列的不同 token 数”，可设置更接近的 reserve，减少 rehash
	//		size_t expect = (size_t)N * expected_tokens_per_seq;
	//		if (expect > word_table.size()) word_table.reserve(expect);
	//	} else {
	//		// 经验：N * 8 个 key（视数据而定）
	//		if (word_table.bucket_count() < (size_t)std::max(16, N * 8))
	//			word_table.reserve(std::max(16, N * 8));
	//	}

	//	word_encodes.resize(max_seq_len);
	//	word_encodes_no.resize(max_seq_len);
	//	counts.resize(N, 0);                 // 注意：保持为 N 大小
	//	visited.clear();  visited.reserve(1 << 14);
	//	out_pairs.clear(); out_pairs.reserve(1 << 14);
	//	A.resize(N);
	//}

void reserve_for(int N, int max_seq_len, size_t expected_tokens_per_seq = 0) {
	// 估算一个目标容量并 reserve
	size_t target = expected_tokens_per_seq
		? (size_t)N * expected_tokens_per_seq
		: std::max<size_t>(16, (size_t)N * 8);

	// robin_hood::unordered_map 没有 bucket_count()，可以直接 reserve
	word_table.reserve(target);   // 或者先判断 capacity() 再 reserve

	word_encodes.resize(max_seq_len);
	word_encodes_no.resize(max_seq_len);
	counts.resize(N, 0);
	visited.clear();  visited.reserve(1 << 10);
	out_pairs.clear(); out_pairs.reserve(1 << 10);
	A.resize(N);
}

// 在新一轮构建前清空所有 bucket（保留容量以复用）
void clear_word_table_rows_keep_capacity() {
		for (auto& kv : word_table) kv.second.clear();
	}
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
void cluster_sequences_st_reuse(
		std::vector<Sequence_new>& seqs,
		std::vector<int>& parent,
		int kmer_size,
		double tau,
		ClusterWS& ws);

void cluster_sequences_st_less10(
		std::vector<Sequence_new>& seqs,
		std::vector<int>& parent,
		int kmer_size = 5,
		double tau = 0.36
		);

void cluster_sequence_singleThread_smallScale_cArray(
    std::vector<Sequence_new>& seqs,
    std::vector<int>& parent,
    int kmer_size,
    double tau
);

#endif
