#include "cdhit.h"

#include <zlib.h>
#include <omp.h>
#include <iostream>
#include <string>
#include <cstring>  // 添加以支持 strlen
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <sys/time.h>
#include <atomic>
#include <unordered_map>  // 替换 robin_hood 为 std::unordered_map (若无 robin_hood，可用此；否则添加 #include "robin_hood.h")
#include <unordered_set>  
#include <cassert>

using namespace std;

// 时间函数
double get_time() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
}

// DSU 结构
struct DSU {
	std::vector<int> p, sz;
	DSU(int n = 0) : p(n, -1), sz(n, 1) {
		for (int i = 0; i < n; ++i) p[i] = i;
	}
	int find(int x) { return p[x] == x ? x : p[x] = find(p[x]); }
	void unite(int a, int b) {
		a = find(a); b = find(b);
		if (a == b) return;
		if (sz[a] < sz[b]) std::swap(a, b);
		p[b] = a; sz[a] += sz[b];
	}
};

// NAA 常量初始化
//#define MAX_UCA 21
//int NAA1, NAA2, NAA3, NAA4, NAA5, NAA6, NAA7, NAA8, NAA9, NAA10, NAA11, NAA12;
//int NAAN_array[13] = {1};
//
//void InitNAA(int max) {
//	NAA1 = NAAN_array[1] = max;
//	NAA2 = NAAN_array[2] = NAA1 * NAA1;
//	NAA3 = NAAN_array[3] = NAA1 * NAA2;
//	NAA4 = NAAN_array[4] = NAA2 * NAA2;
//	NAA5 = NAAN_array[5] = NAA2 * NAA3;
//	NAA6 = NAAN_array[6] = NAA3 * NAA3;
//	NAA7 = NAAN_array[7] = NAA3 * NAA4;
//	NAA8 = NAAN_array[8] = NAA4 * NAA4;
//	NAA9 = NAAN_array[9] = NAA4 * NAA5;
//	NAA10 = NAAN_array[10] = NAA5 * NAA5;
//	NAA11 = NAAN_array[11] = NAA5 * NAA6;
//	NAA12 = NAAN_array[12] = NAA6 * NAA6;
//}

int aa_map[256];

void init_aa_map() {
	fill(begin(aa_map), end(aa_map), 0);
	aa_map['A'] = 0; aa_map['C'] = 1; aa_map['D'] = 2; aa_map['E'] = 3; aa_map['F'] = 4;
	aa_map['G'] = 5; aa_map['H'] = 6; aa_map['I'] = 7; aa_map['K'] = 8; aa_map['L'] = 9;
	aa_map['M'] = 10; aa_map['N'] = 11; aa_map['P'] = 12; aa_map['Q'] = 13; aa_map['R'] = 14;
	aa_map['S'] = 15; aa_map['T'] = 16; aa_map['V'] = 17; aa_map['W'] = 18; aa_map['Y'] = 19;
	aa_map['X'] = 20;
}

// EncodeWords 函数（修改为使用 aa_map 计算编码，无需预映射）
int EncodeWords(const Sequence_new &seq, vector<int> &word_encodes, vector<int> &word_encodes_no, int NAA) {
	const char* seqi = seq.data;
	int len = strlen(seqi);  // 使用 strlen
	int aan_no = len - NAA + 1;
	int j;
	unsigned char k, k1;
	for (j = 0; j < aan_no; j++) {
		const char* word = seqi + j;
		int encode = 0;
		for (k = 0, k1 = NAA - 1; k < NAA; k++, k1--) {
			encode += aa_map[(unsigned char)word[k]] * NAAN_array[k1];  // 修改为使用 aa_map
		}
		word_encodes[j] = encode;
	}
	std::sort(word_encodes.begin(), word_encodes.begin() + aan_no);
	for (j = 0; j < aan_no; j++) word_encodes_no[j] = 1;
	for (j = aan_no - 1; j; j--) {
		if (word_encodes[j] == word_encodes[j - 1]) {
			word_encodes_no[j - 1] += word_encodes_no[j];
			word_encodes_no[j] = 0;
		}
	}
	return 0;
}

// CountWords_SA 函数（保持原样）
int CountWords_SA(int aan_no,
		const std::vector<int>& word_encodes,
		const std::vector<int>& word_encodes_no,
		const std::vector<std::vector<std::pair<int,int>>>& word_table,
		int min_rest,
		int qid,
		std::vector<int>& counts,
		std::vector<int>& visited,
		std::vector<std::pair<int,int>>& out_pairs) {
	out_pairs.clear();
	visited.clear();
	if (aan_no <= 0) return 0;

	auto add_count = [&](int tid, int add) {
		if (counts[tid] == 0) visited.push_back(tid);
		counts[tid] += add;
	};

	for (int j0 = 0; j0 < aan_no; ++j0) {
		int bucket = word_encodes[j0];
		int qcnt = word_encodes_no[j0];
		if (qcnt == 0) continue;

		const auto& hits = word_table[bucket];
		int rest = aan_no - j0 + 1;

		for (size_t k = 0; k < hits.size(); ++k) {
			int tid = hits[k].first;
			if (tid <= qid) break;
			int tcnt = hits[k].second;

			if (min_rest > 0 && rest < min_rest && counts[tid] == 0) continue;
			int add = (qcnt < tcnt) ? qcnt : tcnt;
			add_count(tid, add);
		}
	}

	out_pairs.reserve(visited.size());
	for (int tid : visited) {
		out_pairs.emplace_back(tid, counts[tid]);
		counts[tid] = 0;
	}
	return 0;
}

// Jaccard 计算
static inline double jaccard_from_CAB(int C, int A, int B) {
	int denom = A + B - C;
	return denom > 0 ? double(C) / double(denom) : 0.0;
}

// precompute_edges_jaccard 函数（调整为 Sequence_new）
void precompute_edges_jaccard(
		const std::vector<Sequence_new>& seqs,
		std::vector<std::vector<std::pair<int,int>>>& word_table,
		int kmer_size, double tau,
		std::vector<std::vector<int>>& neigh,
		int nthreads
		) {
	const int N = (int)seqs.size();
	neigh.assign(N, {});

	std::vector<int> A(N);
	for (int i = 0; i < N; ++i) {
		int L = strlen(seqs[i].data);  // 使用 strlen
		A[i] = std::max(0, L - kmer_size + 1);
	}

	size_t max_len = 0;
	for (auto& s : seqs) {
		max_len = std::max(max_len, strlen(s.data));  // 使用 strlen
	}

	std::atomic<int> progress{0};
//	int nthreads = 1;
//#pragma omp parallel
//	{
//#pragma omp single
//		nthreads = omp_get_num_threads();
//	}

	std::vector<std::vector<std::pair<int,int>>> thread_edges(nthreads);

#pragma omp parallel num_threads(nthreads)
	{
		int tid = omp_get_thread_num();
		auto &edges = thread_edges[tid];
		edges.clear();
		edges.reserve(1<<20);

		std::vector<int> word_encodes(max_len);
		std::vector<int> word_encodes_no(max_len);

		std::vector<int> counts(seqs.size(), 0);
		std::vector<int> visited; visited.reserve(1<<14);
		std::vector<std::pair<int,int>> out_pairs;

#pragma omp for schedule(dynamic,1)
		for (int i = 0; i < N; ++i) {
			if (A[i] <= 0) { ++progress; continue; }

			EncodeWords(seqs[i], word_encodes, word_encodes_no, kmer_size);  // const &
			CountWords_SA(A[i], word_encodes, word_encodes_no, word_table, 0, i, counts, visited, out_pairs);

			for (auto &pr : out_pairs) {
				int j = pr.first;
				if (A[j] <= 0) continue;
				int C = pr.second;
				double jac = jaccard_from_CAB(C, A[i], A[j]);
				if (jac >= tau) {
					edges.emplace_back(i, j);
				}
			}

//			int p = ++progress;
//			if ((p % 1000) == 0) {
//				double percent = 100.0 * p / N;
//				std::cerr << "\rPhase A (precompute): " << p << "/" << N
//					<< " (" << percent << "%)" << std::flush;
//			}
		}
	}

	// 合并边
	std::vector<size_t> cnt(N, 0);
	for (auto &vec : thread_edges)
		for (auto &e : vec)
			++cnt[e.first];

	for (int u = 0; u < N; ++u)
		if (cnt[u]) neigh[u].reserve(neigh[u].size() + cnt[u]);

	for (auto &vec : thread_edges)
		for (auto &e : vec)
			neigh[e.first].push_back(e.second);

//	std::cerr << "\n";
}

// 主函数实现
void cluster_sequences(
		std::vector<Sequence_new>& seqs,
		std::vector<int>& parent,
		int kmer_size,
		double tau,
		int nthreads
		) {
	//std::vector<Sequence_new> seqs = sequences;  // 复制
//	if (max_num_seqs > 0 && (int)seqs.size() > max_num_seqs) {
//		seqs.resize(max_num_seqs);
//	}

	assert(nthreads>=1);
	InitNAA(MAX_UAA);
	init_aa_map();

	int N = (int)seqs.size();
	int max_seq_len = 0;
	for (const auto& s : seqs) {
		max_seq_len = std::max(max_seq_len, (int)strlen(s.data));
	}

	// 无预处理映射（data 为 const）

	// 排序序列按长度降序
	sort(seqs.begin(), seqs.end(),
			[](const Sequence_new& a, const Sequence_new& b) {
			return strlen(a.data) > strlen(b.data);
			});

	// 构建 word_table
	int table_size = 1;
	for (int i = 0; i < kmer_size; ++i) table_size *= MAX_UAA;
	vector<vector<pair<int, int>>> word_table(table_size);

//	int nthreads = 1;
//#pragma omp parallel
//	{
//#pragma omp single
//		nthreads = omp_get_num_threads();
//	}

	vector<vector<vector<pair<int,int>>>> local_tables(
			nthreads, vector<vector<pair<int,int>>>(table_size)
			);

	double t1 = get_time();
#pragma omp parallel num_threads(nthreads)
	{
		int tid = omp_get_thread_num();
		auto& local = local_tables[tid];

		vector<int> word_encodes(max_seq_len);
		vector<int> word_encodes_no(max_seq_len);

#pragma omp for schedule(dynamic,1)
		for (int seq_id = 0; seq_id < N; ++seq_id) {
			auto& s = seqs[seq_id];
			int len = strlen(s.data);
			if (len < kmer_size) continue;

			EncodeWords(s, word_encodes, word_encodes_no, kmer_size);

			int kmer_no = len - kmer_size + 1;
			for (int j = 0; j < kmer_no; ++j) {
				int bucket = word_encodes[j];
				int count = word_encodes_no[j];
				if (count > 0) {
					local[bucket].emplace_back(seq_id, count);
				}
			}
		}
	}

#pragma omp parallel for schedule(static) num_threads(nthreads)
	for (long long b = 0; b < (long long)table_size; ++b) {
		size_t add = 0;
		for (int t = 0; t < nthreads; ++t)
			add += local_tables[t][b].size();
		auto& dst = word_table[b];
		if (add) dst.reserve(dst.size() + add);

		for (int t = 0; t < nthreads; ++t) {
			auto& src = local_tables[t][b];
			if (!src.empty()) {
				dst.insert(dst.end(),
						make_move_iterator(src.begin()),
						make_move_iterator(src.end()));
				src.clear();
				src.shrink_to_fit();
			}
		}
	}

	// 排序每个 bucket 按 seq_id 降序
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
	for (size_t i = 0; i < word_table.size(); ++i) {
		auto& row = word_table[i];
		std::sort(row.begin(), row.end(),
				[](const std::pair<int,int>& a, const std::pair<int,int>& b) {
				return a.first > b.first;
				});
	}

	double t2 = get_time();
	//std::cerr << "Word table build time: " << (t2 - t1) << " seconds" << std::endl;

	// 预计算边缘
	double t3 = get_time();
	std::vector<std::vector<int>> neigh;
	precompute_edges_jaccard(seqs, word_table, kmer_size, tau, neigh, nthreads);

	size_t total_edges = 0;
	for (const auto& vec : neigh) total_edges += vec.size();
	//std::cerr << "Total edges: " << total_edges << std::endl;

	double t4 = get_time();

	// DSU 聚类
	DSU dsu(N);
	for (int u = 0; u < N; ++u) {
		for (int v : neigh[u]) dsu.unite(u, v);
	}

	// 输出到 parent: parent[i] = root of i (root points to itself)
	// parent.resize(N);
	/// i = local seq_id
	/// to global id
	for (int i = 0; i < N; ++i) {
		parent[seqs[i].seq_id] = seqs[dsu.find(i)].seq_id;
	}

	double t5 = get_time();
	//std::cerr << "Jaccard filtering time: " << (t4 - t3) << " s" << std::endl;
	//std::cerr << "DSU clustering time: " << (t5 - t4) << " s" << std::endl;
	//std::cerr << "Number of clusters: " << dsu.sz.size() - std::count(parent.begin(), parent.end(), -1) << std::endl;  // 粗略计数
	//std::unordered_set<int> unique_roots(parent.begin(), parent.end());
	//std::cerr << "Number of clusters: " << unique_roots.size() << std::endl;
}


