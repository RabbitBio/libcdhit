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
#include <unordered_map>  // 替换 robin_hood 为 std::unordered_map (若无 robin_hood，可用此；否则添加 
#include "robin_hood.h"
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


// Jaccard 计算
static inline double jaccard_from_CAB(int C, int A, int B) {
	int denom = A + B - C;
	return denom > 0 ? double(C) / double(denom) : 0.0;
}
// 仅统计 id<qid（行按 seq_id 升序存储，遇到 >=qid 早停）
// 输出 out_pairs: vector<pair<tid, C>>
int CountWords_SA(int aan_no,
		const std::vector<int>& word_encodes,
		const std::vector<int>& word_encodes_no,
		const std::vector<std::vector<std::pair<int,int>>>& word_table,
		int min_rest,
		int qid,
		// 稀疏累加器（线程私有）
		std::vector<int>& counts,
		std::vector<int>& visited,
		std::vector<std::pair<int,int>>& out_pairs)
{
	out_pairs.clear();
	visited.clear();
	if (aan_no <= 0) return 0;

	auto add_count = [&](int tid, int add){
		if (counts[tid] == 0) visited.push_back(tid);
		counts[tid] += add;
	};

	for (int j0 = 0; j0 < aan_no; ++j0) {
		int bucket = word_encodes[j0];
		int qcnt   = word_encodes_no[j0];
		if (qcnt == 0) continue;

		const auto& hits = word_table[bucket];  // 已按 seq_id 升序
		int rest = aan_no - j0 + 1;

		// 只扫 id<qid 的前缀
		for (size_t k = 0; k < hits.size(); ++k) {
			int tid  = hits[k].first;
			if (tid >= qid) break;  // 及早终止
			int tcnt = hits[k].second;

			if (min_rest > 0 && rest < min_rest && counts[tid] == 0) continue;
			int add = (qcnt < tcnt) ? qcnt : tcnt;
			add_count(tid, add);
		}
	}

	// 导出并清零
	out_pairs.reserve(visited.size());
	for (int tid : visited) {
		out_pairs.emplace_back(tid, counts[tid]);
		counts[tid] = 0;
	}
	return 0;
}

//using WordTable = std::unordered_map<int, std::vector<std::pair<int,int>>>;
int CountWords_SA(
		int aan_no,
		const std::vector<int>& word_encodes,
		const std::vector<int>& word_encodes_no,
		robin_hood::unordered_map<int, std::vector<std::pair<int,int>>> & word_table,   // 修改：由 vector<vector<...>> 改为 unordered_map
		int min_rest,
		int qid,
		// 稀疏累加器（线程私有）
		std::vector<int>& counts,
		std::vector<int>& visited,
		std::vector<std::pair<int,int>>& out_pairs)
{
	out_pairs.clear();
	visited.clear();
	if (aan_no <= 0) return 0;

	auto add_count = [&](int tid, int add){
		if (counts[tid] == 0) visited.push_back(tid);
		counts[tid] += add;
	};

	for (int j0 = 0; j0 < aan_no; ++j0) {
		const int bucket = word_encodes[j0];
		const int qcnt   = word_encodes_no[j0];
		if (qcnt == 0) continue;

		const int rest = aan_no - j0 + 1;

		// 用 find()，没有该 bucket 就跳过；不会插入空向量
		auto it = word_table.find(bucket);
		if (it == word_table.end()) continue;

		const auto& hits = it->second;  // 已保证按 seq_id 升序
		for (size_t k = 0; k < hits.size(); ++k) {
			const int tid  = hits[k].first;
			if (tid >= qid) break;      // 及早终止：只扫 id<qid
			const int tcnt = hits[k].second;

			if (min_rest > 0 && rest < min_rest && counts[tid] == 0) continue;
			const int add = (qcnt < tcnt) ? qcnt : tcnt;
			add_count(tid, add);
		}
	}

	// 导出并清零
	out_pairs.reserve(visited.size());
	for (int tid : visited) {
		out_pairs.emplace_back(tid, counts[tid]);
		counts[tid] = 0;
	}
	return 0;
}

void precompute_edges_jaccard(
		const std::vector<Sequence_new>& seqs,
		std::vector<std::vector<std::pair<int,int>>>& word_table,
		int kmer_size, double tau,
		DSU& global_dsu,// 输出：全局 DSU
		int nthreads
		)
{
	const int N = (int)seqs.size();
	// 预先计算 A[i] = |k-mers(i)| = Li - k + 1
	std::vector<int> A(N);
	for (int i = 0; i < N; ++i) {
		int L = strlen(seqs[i].data);
		A[i] = std::max(0, L - kmer_size + 1);
	}

	// 为线程私有缓存做准备
	size_t max_len = 0;
	for (auto& s : seqs) max_len = std::max(max_len, strlen(s.data));

	std::atomic<int> progress{0};
	double tA = get_time();

	assert(nthreads >= 1);

	// 每个线程一个本地 DSU
	std::vector<DSU> thread_dsu(nthreads, DSU(N));

#pragma omp parallel num_threads(nthreads)
	{
		int tid = omp_get_thread_num();

		std::vector<int> word_encodes(max_len);
		std::vector<int> word_encodes_no(max_len);

		std::vector<int> counts(N, 0);                // 线程私有
		std::vector<int> visited; 
		visited.reserve(1<<14);       // 线程私有 //TODO: reserve how many?
		std::vector<std::pair<int,int>> out_pairs;              // 线程私有


#pragma omp for schedule(dynamic,1)
		for (int i = 0; i < N; ++i) {
			if (A[i] <= 0) { ++progress; continue; }

			EncodeWords(seqs[i], word_encodes, word_encodes_no, kmer_size);
			CountWords_SA(A[i], word_encodes, word_encodes_no, word_table, 0, i, counts, visited, out_pairs);

			// 只保留 Jaccard 过阈值的边（i>j）
			for (auto &pr : out_pairs) {
				int j = pr.first;       // j < i
				if (A[j] <= 0) continue;
				int C = pr.second;
				double jac = jaccard_from_CAB(C, A[i], A[j]);
				if (jac >= tau) {
					thread_dsu[tid].unite(i, j); // 线程本地 unite
				}
			}

			int p = ++progress;
			if ((p % 1000) == 0) {
				double percent = 100.0 * p / N;
				std::cout << "\rPhase A (precompute): " << p << "/" << N
					<< " (" << percent << "%)" << std::flush;
			}
		}
	} // 并行区结束

	double tB = get_time();
	// ===== 合并阶段（单线程）=====
	// 合并所有线程的 DSU
	global_dsu = std::move(thread_dsu[0]);
	for (int t = 1; t < nthreads; ++t) {
		for (int i = 0; i < N; ++i) {
			int rg = global_dsu.find(i);
			int rt = thread_dsu[t].find(i);
			global_dsu.unite(rg, rt);
		}
	}

	std::cout << "\n";
	double tC = get_time();
	//std::cerr << "Precompute (edges) time: " << (tC - tA) << " s\n";
	//std::cerr << "Precompute (edges) time (parallel region): " << (tB - tA) << " s\n";
	//std::cerr << "Merge DSU time: " << (tC - tB) << " s\n";
}


// 主函数实现
void cluster_sequences(
		std::vector<Sequence_new>& seqs,
		std::vector<int>& parent,
		int kmer_size,
		double tau,
		int nthreads
		) {

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

	vector<vector<pair<int, int>>> word_table(table_size); //TODO: buffering

	vector<vector<vector<pair<int,int>>>> local_tables(
			nthreads, vector<vector<pair<int,int>>>(table_size)
			);//TODO: buffering

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

	local_tables.clear();
	local_tables.shrink_to_fit();

	// 排序每个 bucket 按 seq_id 升序
#pragma omp parallel for schedule(dynamic) num_threads(nthreads)
	for (size_t i = 0; i < word_table.size(); ++i) {
		auto& row = word_table[i];
		std::sort(row.begin(), row.end(),
				[](const std::pair<int,int>& a, const std::pair<int,int>& b) {
				return a.first < b.first;
				});
	}

	double t2 = get_time();
	//std::cerr << "Word table build time: " << (t2 - t1) << " seconds" << std::endl;

	// 预计算边缘
	double t3 = get_time();

	//DSU dsu(seqs.size());
	DSU dsu;
	precompute_edges_jaccard(seqs, word_table, kmer_size, tau, dsu, nthreads);

	double t4 = get_time();
	
	// 输出到 parent: parent[i] = root of i (root points to itself)
	/// i = local seq_id
	/// to global id
	for (int i = 0; i < N; ++i) {
		parent[seqs[i].seq_id] = seqs[dsu.find(i)].seq_id;
	}

	double t5 = get_time();
	//std::cerr << "Jaccard filtering time: " << (t4 - t3) << " s" << std::endl;
	//std::cerr << "DSU clustering time: " << (t5 - t4) << " s" << std::endl;
	//std::unordered_set<int> unique_roots(parent.begin(), parent.end());
	//std::cerr << "Number of clusters: " << unique_roots.size() << std::endl;
}

void cluster_sequences_st_old(
		std::vector<Sequence_new>& seqs,
		std::vector<int>& parent,
		int kmer_size,
		double tau
		) {

	InitNAA(MAX_UAA); //TODO:move out
	init_aa_map(); //TODO:move out

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

	vector<vector<pair<int, int>>> word_table(table_size); //TODO: buffering
	vector<int> word_encodes(max_seq_len);
	vector<int> word_encodes_no(max_seq_len);

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
				word_table[bucket].emplace_back(seq_id, count);
			}
		}
	}

	// 排序每个 bucket 按 seq_id 升序
	for (size_t i = 0; i < word_table.size(); ++i) {
		auto& row = word_table[i];
		std::sort(row.begin(), row.end(),
				[](const std::pair<int,int>& a, const std::pair<int,int>& b) {
				return a.first < b.first;
				});
	}


	std::vector<int> A(N);//TODO:buffering
	for (int i = 0; i < N; ++i) {
		int L = strlen(seqs[i].data);
		A[i] = std::max(0, L - kmer_size + 1);
	}

	DSU dsu(seqs.size()); //buffering
	std::vector<int> counts(N, 0);//buffering
	std::vector<int> visited; //buffering
	visited.reserve(1<<14);   //TODO: reserve how many?
	std::vector<std::pair<int,int>> out_pairs; //buffering             

	for (int i = 0; i < N; ++i) {

		EncodeWords(seqs[i], word_encodes, word_encodes_no, kmer_size);
		CountWords_SA(A[i], word_encodes, word_encodes_no, word_table, 0, i, counts, visited, out_pairs);

		// 只保留 Jaccard 过阈值的边（i>j）
		for (auto &pr : out_pairs) {
			int j = pr.first;       // j < i
			int C = pr.second;
			double jac = jaccard_from_CAB(C, A[i], A[j]);
			if (jac >= tau) {
				dsu.unite(i, j); 
			}
		}

	}

	for (int i = 0; i < N; ++i) {
		parent[seqs[i].seq_id] = seqs[dsu.find(i)].seq_id;
	}

}
void cluster_sequences_st(
		std::vector<Sequence_new>& seqs,
		std::vector<int>& parent,
		int kmer_size,
		double tau)
{
	InitNAA(MAX_UAA); // TODO: 可外移
	init_aa_map();    // TODO: 可外移

	const int N = (int)seqs.size();
	int max_seq_len = 0;
	for (auto& s : seqs) {
		s.length = strlen(s.data);
		max_seq_len = std::max(max_seq_len, s.length);
	}

	// 按长度降序
	std::sort(seqs.begin(), seqs.end(),
			[](const Sequence_new& a, const Sequence_new& b){
			//return strlen(a.data) > strlen(b.data);
			return a.length > b.length;
			});

	// —— 改为基于 hashmap 的稀疏 word table ——
	robin_hood::unordered_map<int, std::vector<std::pair<int,int>>> word_table;
	// 估计一个合理的 bucket 数量做 reserve（经验值，可按数据集调）
	// 假设平均每条序列有 ~ (len - k + 1) 个 k-mer，且有一定重复：
	// 这里保守按 N*8 预留，避免频繁 rehash（不强制）
	word_table.reserve(std::max(16, N * 8));

	std::vector<int> word_encodes(max_seq_len);
	std::vector<int> word_encodes_no(max_seq_len);

	// 构建稀疏表：只为实际出现的 bucket 建立项
	for (int seq_id = 0; seq_id < N; ++seq_id) {
		const auto& s = seqs[seq_id];
		const int len = s.length;
		if (len < kmer_size) continue;

		EncodeWords(s, word_encodes, word_encodes_no, kmer_size);
		const int kmer_no = len - kmer_size + 1;

		for (int j = 0; j < kmer_no; ++j) {
			const int bucket = word_encodes[j];
			const int count  = word_encodes_no[j];
			if (count > 0) {
				word_table[bucket].emplace_back(seq_id, count);
			}
		}
	}

	// 每个 bucket 内按 seq_id 升序，保证 CountWords_SA 的“及早终止”成立
	for (auto& kv : word_table) {
		auto& row = kv.second;
		std::sort(row.begin(), row.end(),
				[](const std::pair<int,int>& a, const std::pair<int,int>& b){
				return a.first < b.first;
				});
	}

	// 预计算每条序列的 A[i] = 有效 k-mer 数（重复保留；与你原逻辑一致）
	std::vector<int> A(N);
	for (int i = 0; i < N; ++i) {
		const int L = (int)strlen(seqs[i].data);
		A[i] = std::max(0, L - kmer_size + 1);
	}

	DSU dsu(seqs.size());
	std::vector<int> counts(N, 0);
	std::vector<int> visited;  visited.reserve(1<<14);
	std::vector<std::pair<int,int>> out_pairs;

	for (int i = 0; i < N; ++i) {
		EncodeWords(seqs[i], word_encodes, word_encodes_no, kmer_size);

		// 调用修改后的 hashmap 版本
		CountWords_SA(A[i], word_encodes, word_encodes_no,
				word_table, /*min_rest=*/0, /*qid=*/i,
				counts, visited, out_pairs);

		// 只保留 Jaccard 过阈值的边（i>j）
		for (auto &pr : out_pairs) {
			const int j = pr.first;   // j < i
			const int C = pr.second;
			const double jac = jaccard_from_CAB(C, A[i], A[j]);
			if (jac >= tau) {
				dsu.unite(i, j);
			}
		}
	}

	// 写回代表元（保持与原始 seq_id 的对应）
	for (int i = 0; i < N; ++i) {
		parent[seqs[i].seq_id] = seqs[dsu.find(i)].seq_id;
	}
}
