#include "kseq.h"
#include "CLI11.hpp"
#include <zlib.h>
#include <omp.h>
#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <sys/time.h>
#include <atomic>
#include "robin_hood.h"

using namespace std;

#define MAX_UCA 21
KSEQ_INIT(gzFile, gzread)

double get_time()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
}

struct Sequence
{
	string name;
	string comment;
	string seq;
};

struct IndexCount
{
	int index;
	int count;
};

/// This is for NAAN array.
/// FIXME: reimplement it with constexpr
int NAA1;
int NAA2;
int NAA3;
int NAA4;
int NAA5;
int NAA6;
int NAA7;
int NAA8;
int NAA9;
int NAA10;
int NAA11;
int NAA12;
int NAAN_array[13] = { 1 };

void InitNAA(int max)
{
	NAA1 = NAAN_array[1] = max;
	NAA2 = NAAN_array[2] = NAA1 * NAA1;
	NAA3 = NAAN_array[3] = NAA1 * NAA2;
	NAA4 = NAAN_array[4] = NAA2 * NAA2;
	NAA5 = NAAN_array[5] = NAA2 * NAA3;
	NAA6 = NAAN_array[6] = NAA3 * NAA3;
	NAA7 = NAAN_array[7] = NAA3 * NAA4;
	NAA8 = NAAN_array[8] = NAA4 * NAA4;
	NAA9 = NAAN_array[9] = NAA4 * NAA5;
	NAA10 = NAAN_array[10] = NAA5 * NAA5;
	NAA11 = NAAN_array[11] = NAA5 * NAA6;
	NAA12 = NAAN_array[12] = NAA6 * NAA6;
}
int aa_map[256];

void init_aa_map() {
    fill(begin(aa_map), end(aa_map), 0);
    aa_map['A'] = 0;
    aa_map['C'] = 1;
    aa_map['D'] = 2;
    aa_map['E'] = 3;
    aa_map['F'] = 4;
    aa_map['G'] = 5;
    aa_map['H'] = 6;
    aa_map['I'] = 7;
    aa_map['K'] = 8;
    aa_map['L'] = 9;
    aa_map['M'] = 10;
    aa_map['N'] = 11;
    aa_map['P'] = 12;
    aa_map['Q'] = 13;
    aa_map['R'] = 14;
    aa_map['S'] = 15;
    aa_map['T'] = 16;
    aa_map['V'] = 17;
    aa_map['W'] = 18;
    aa_map['Y'] = 19;
    aa_map['X'] = 20;
}

/// This function is imported from cdhit
/// Raw seq is map to [0-21] when reading data from disk
/// NAA mean k-mer size
/// word_encodes and word_encodes_no are pre-allocated with max seq length
int EncodeWords(Sequence &seq, vector<int> & word_encodes, vector<int> & word_encodes_no, int NAA)
{
    const char* seqi = seq.seq.c_str();
	int len = seq.seq.size();
	// check_word_encodes
	int aan_no = len - NAA + 1;
	int i, j, i0, i1;
	unsigned char k, k1;
	for (j = 0; j < aan_no; j++) {
		const char* word = seqi + j;
		int encode = 0;
		//for (k = 0, k1 = NAA - 1; k < NAA; k++, k1--) encode += aa_map[(unsigned char)word[k]] * NAAN_array[k1];
		for (k = 0, k1 = NAA - 1; k < NAA; k++, k1--) encode += word[k] * NAAN_array[k1];
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

/// This is modified from cdhit
int CountWords(int aan_no,
               vector<int>& word_encodes,
               vector<int>& word_encodes_no,
               robin_hood::unordered_map<int, int>& lookCounts,      // 输出：候选目标及公共k-mer数
               vector<vector<pair<int,int>>>& word_table,
               int min_rest, int qid)                         // 与原逻辑一致：剩余kmer阈值过滤
{
	lookCounts.clear();

	// 2) 按 query 的编码遍历倒排表并累积
	int j0 = 0;
	const int* we  = word_encodes.data();

	for (; j0 < aan_no; ++j0) {
		int bucket = word_encodes[j0];
		int qcnt   = word_encodes_no[j0];
		if (qcnt == 0) continue;

		const auto& hits = word_table[bucket];
		int rest = aan_no - j0 + 1;          // 剩余k-mer数（与原代码一致）
		for (const auto& hit : hits) {
			int target_id   = hit.first;           // target seq id
			if (target_id <= qid) break;
			int tcnt  = hit.second;          // 该 target 在此 bucket 的计数
			if (rest < min_rest && !lookCounts.count(target_id)) {
				continue; // 剪枝：首次遇到且剩余不足
			}
			int add   = std::min(qcnt, tcnt);
			lookCounts[target_id] += add;
		}
	}

	// 如果你需要与老代码完全兼容的“哨兵”结尾，可取消注释：
	// lookCounts.push_back({-1, 0});

	return 0;
}
// 仅统计 id>qid（行按 seq_id 降序存储，遇到 <=qid 早停）
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

		const auto& hits = word_table[bucket];  // 已按 seq_id 降序
		int rest = aan_no - j0 + 1;

		// 只扫 id>qid 的前缀
		for (size_t k = 0; k < hits.size(); ++k) {
			int tid  = hits[k].first;
			if (tid <= qid) break;  // 及早终止
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
// ===== helper: Jaccard =====
static inline double jaccard_from_CAB(int C, int A, int B) 
{
	int denom = A + B - C;
	//return denom > 0 ? double(C) / double(denom) : 0.0;
	return double(C) / double(denom);
}

// ===== 阶段A：并行预计算 >=阈值 的稀疏图（只存 i<j 的邻接） =====
// 依赖你现有的 EncodeWords / CountWords / word_table
// word_table[bucket] = vector<pair<seq_id, count>>
void precompute_edges_jaccard(
		const std::vector<Sequence>& seqs,
		std::vector<std::vector<std::pair<int,int>>>& word_table,
		int kmer_size, double tau,
		std::vector<std::vector<int>>& neigh  // 输出：neigh[u] 存 v(>u)
		)
{
	const int N = (int)seqs.size();
	neigh.assign(N, {});                  // 只保留 i<j 方向
	// 预先计算 A[i] = |k-mers(i)| = Li - k + 1
	std::vector<int> A(N);
	for (int i = 0; i < N; ++i) {
		int L = (int)seqs[i].seq.size();
		A[i] = std::max(0, L - kmer_size + 1);
	}

	// 为线程私有缓存做准备
	size_t max_len = 0;
	for (auto& s : seqs) max_len = std::max(max_len, s.seq.size());

	std::atomic<int> progress{0};
	double tA = get_time();
	int nthreads = 1;
#pragma omp parallel
	{
#pragma omp single
		nthreads = omp_get_num_threads();
	}
	// 每个线程一块本地边缓冲
	std::vector<std::vector<std::pair<int,int>>> thread_edges(nthreads);

#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		auto &edges = thread_edges[tid];
		edges.clear();
		edges.reserve(1<<20); // 视数据量调整

		std::vector<int> word_encodes(max_len);
		std::vector<int> word_encodes_no(max_len);
		//robin_hood::unordered_map<int,int> local;
		
		std::vector<int> counts(seqs.size(), 0);                // 线程私有
		std::vector<int> visited; visited.reserve(1<<14);       // 线程私有
		std::vector<std::pair<int,int>> out_pairs;              // 线程私有


#pragma omp for schedule(dynamic,1)
		for (int i = 0; i < N; ++i) {
			if (A[i] <= 0) { ++progress; continue; }

			EncodeWords((Sequence&)seqs[i], word_encodes, word_encodes_no, kmer_size);
			//local.clear();
			//CountWords(A[i], word_encodes, word_encodes_no, local, word_table, /*min_rest=*/0, i);
			CountWords_SA(A[i], word_encodes, word_encodes_no, word_table, /*min_rest=*/A[i]/2, i, counts, visited, out_pairs);

//			for (auto &kv : local) {
//				int j = kv.first;
//				if (j == i || A[j] <= 0) continue;
//				int u = i < j ? i : j;
//				int v = i < j ? j : i;
//				int C = kv.second;
//				double jac = jaccard_from_CAB(C, A[u], A[v]);
//				if (jac >= tau) {
//					edges.emplace_back(u, v);   // 仅写本线程缓冲，无需 critical
//				}
//			}

			// 只保留 Jaccard 过阈值的边（i<j）
			for (auto &pr : out_pairs) {
				int j = pr.first;       // j > i
				if (A[j] <= 0) continue;
				int C = pr.second;
				double jac = jaccard_from_CAB(C, A[i], A[j]);
				if (jac >= tau) {
					edges.emplace_back(i, j); // 线程私有 edges，无需加锁
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
	// 1) 统计每个 u 的边数，便于一次性 reserve
	std::vector<size_t> cnt(N, 0);
	for (auto &vec : thread_edges)
		for (auto &e : vec)
			++cnt[e.first];

	for (int u = 0; u < N; ++u)
		if (cnt[u]) neigh[u].reserve(neigh[u].size() + cnt[u]);

	// 2) 顺序合并
	for (auto &vec : thread_edges)
		for (auto &e : vec)
			neigh[e.first].push_back(e.second);

	std::cout << "\n";
	double tC = get_time();
	std::cerr << "Precompute (edges) time: " << (tC - tA) << " s\n";
	std::cerr << "Precompute (edges) time (parallel region): " << (tB - tA) << " s\n";
}

// ===== 阶段B：顺序提交（贪心增量，长度降序） =====
struct GreedyResult {
	std::vector<int>  parent;     // 冗余点的代表；代表自身为 -1
	std::vector<char> redundant;  // 是否被吸收
	int centers = 0;              // 代表数量
};

GreedyResult greedy_commit_by_order_already_sorted(
		const std::vector<Sequence>& seqs,
		const std::vector<std::vector<int>>& neigh,
		int kmer_size
		)
{
	const int N = (int)seqs.size();
	std::vector<char> redundant(N, 0);
	std::vector<int>  parent(N, -1);

	double tC0 = get_time();
	// 注意：seqs 已经按长度降序排好，直接用 i=0..N-1
	for (int i = 0; i < N; ++i) {
		if (redundant[i]) continue;   // 若已被更长代表吸收则跳过

		// i 成为代表，吸收它的邻居（仅存了 i<j 的边）
		for (int v : neigh[i]) {
			if (!redundant[v]) {
				redundant[v] = 1;
				parent[v] = i;
			}
		}

		if ((i % 1000) == 0) {
			double percent = 100.0 * (i+1) / N;
			std::cout << "\rPhase B (commit): " << (i+1) << "/" << N
				<< " (" << percent << "%)" << std::flush;
		}
	}
	std::cout << "\n";
	double tC1 = get_time();

	int centers = 0;
	for (int i = 0; i < N; ++i) if (!redundant[i]) ++centers;

	std::cerr << "Commit (greedy) time: " << (tC1 - tC0) << " s\n";
	return { std::move(parent), std::move(redundant), centers };
}

static inline uint64_t splitmix64(uint64_t x){
	x += 0x9e3779b97f4a7c15ull;
	x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
	x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
	return x ^ (x >> 31);
}
static inline double u01_from_hash(uint64_t h){
	const double k = 1.0 / (double(1ull<<53));
	return ((h>>11) + 0.5) * k; // (0,1)
}
static inline double gamma2_from_2u(double u1, double u2){
	return -std::log(u1) - std::log(u2); // Gamma(k=2,theta=1)
}
struct RandTriplet { double r, c, beta; };
static inline RandTriplet rand_triple(uint32_t feature_code, uint32_t dim, uint64_t seed){
	uint64_t h1 = splitmix64(seed ^ (uint64_t(feature_code)<<32) ^ dim);
	uint64_t h2 = splitmix64(h1 + 0x9e3779b97f4a7c15ull);
	uint64_t h3 = splitmix64(h2 + 0x9e3779b97f4a7c15ull);
	double u1 = u01_from_hash(h1), u2 = u01_from_hash(h2), u3 = u01_from_hash(h3);
	RandTriplet rt;
	rt.r    = std::max(1e-12, gamma2_from_2u(u1, u2));
	rt.c    = std::max(1e-12, gamma2_from_2u(u2, u3));
	rt.beta = std::min(1.0-1e-12, std::max(1e-12, u3));
	return rt;
}

// 输入：升序唯一的 codes/cnts（长度 nnz）
// 输出：S 维 argmin_code（每维保存赢家的 k-mer 编码值）
static inline void weighted_minhash_icws(const std::vector<int>& codes,
		const std::vector<int>& cnts,
		int nnz,
		int S, uint64_t seed,
		std::vector<uint32_t>& argmin_code)
{
	argmin_code.assign(S, 0);
	std::vector<double> best_a(S, std::numeric_limits<double>::infinity());

	for (int t = 0; t < nnz; ++t) {
		int code = codes[t], w = cnts[t];
		if (w <= 0) continue;
		double logw = std::log(double(w));
		for (int d = 0; d < S; ++d) {
			RandTriplet rt = rand_triple((uint32_t)code, (uint32_t)d, seed);
			double tt = std::floor(logw / rt.r + rt.beta);
			double y  = std::exp(rt.r * (tt - rt.beta));   // (0,w]
			double a  = rt.c / (y * std::exp(rt.r));       // 目标：最小
			if (a < best_a[d]) {
				best_a[d]      = a;
				argmin_code[d] = (uint32_t)code;           // 保存 k-mer 编码
			}
		}
	}
}

// ======= 构建每条序列的 WMH sketch =======

void build_wmh_sketches(std::vector<Sequence>& seqs,
		int kmer_size,
		int S, uint64_t seed,
		std::vector<std::vector<uint32_t>>& wmh_argmin)
{
	const int N = (int)seqs.size();
	wmh_argmin.assign(N, std::vector<uint32_t>(S, 0));

	size_t max_len = 0;
	for (auto& s : seqs) max_len = std::max(max_len, s.seq.size());

#pragma omp parallel
	{
		std::vector<int> codes(max_len), cnts(max_len);

#pragma omp for schedule(dynamic,1)
		for (int i = 0; i < N; ++i) {
			Sequence & s = seqs[i];
			if ((int)s.seq.size() < kmer_size) {
				std::fill(wmh_argmin[i].begin(), wmh_argmin[i].end(), 0u);
				continue;
			}
			// 你的 EncodeWords 应当输出去重后的 (code,count)，且 codes 升序
			EncodeWords(s, codes, cnts, kmer_size);
			int nnz = std::max(0, (int)s.seq.size() - kmer_size + 1); // 或者 EncodeWords 返回的有效长度

			weighted_minhash_icws(codes, cnts, nnz, S, seed, wmh_argmin[i]);
		}
	}
}

// ======= WMH 倒排表（band=1）：每一维一张，以 kmer code 为索引，只存 seq_id =======

// vector 结构：S 张表；每张表是 table_size 个桶，每桶一个 vector<int>（存 seq_id，最终升序）
using WMHWordTables = std::vector<std::vector<std::vector<int>>>;

void build_wmh_word_tables(const std::vector<std::vector<uint32_t>>& wmh_argmin,
		int S,
		size_t table_size,  // = MAX_UCA^kmer_size，与 EncodeWords 的编码空间一致
		WMHWordTables& tables_out)
{
	const int N = (int)wmh_argmin.size();
	tables_out.clear();
	tables_out.resize(S);
	for (int d = 0; d < S; ++d) {
		tables_out[d].assign(table_size, {});
	}

	// 填充
	for (int i = 0; i < N; ++i) {
		const auto& sk = wmh_argmin[i];
		for (int d = 0; d < S; ++d) {
			uint32_t code = sk[d];
			// 假定 code < table_size
			tables_out[d][(size_t)code].push_back(i);
		}
	}

	// 每个 posting 升序，便于查询时用二分跳到 j>i 的后缀
#pragma omp parallel for schedule(dynamic)
	for (int d = 0; d < S; ++d) {
		auto& table = tables_out[d];
		for (size_t b = 0; b < table.size(); ++b) {
			auto& lst = table[b];
			if (!lst.empty()) std::sort(lst.begin(), lst.end());
		}
	}
}

// ======= 预计算邻接：WMH-only，使用每维倒排（vector 结构），不做精算 =======

void precompute_edges_wmh_only_tables(
		const std::vector<std::vector<uint32_t>>& wmh_argmin, // [N][S]
		const WMHWordTables& tables,                          // S 个表，每表 table_size 桶
		double tau,
		std::vector<std::vector<int>>& neigh                  // 输出：neigh[i] 存 j(>i)
		){
	const int N = (int)wmh_argmin.size();
	const int S = N ? (int)wmh_argmin[0].size() : 0;
	const int need = (int)std::ceil(tau * S);

	neigh.assign(N, {});

	// 稀疏累加器：线程私有（counts[j] 是票数；visited 记录本次触及的 j）
	std::atomic<int> progress{0};

#pragma omp parallel
	{
		std::vector<int> counts(N, 0);
		std::vector<int> visited; visited.reserve(1<<14);
		std::vector<int> out;     out.reserve(1024);

		auto add_vote = [&](int j){
			if (counts[j] == 0) visited.push_back(j);
			++counts[j];
		};

#pragma omp for schedule(dynamic,1)
		for (int i = 0; i < N; ++i) {
			const auto& sk = wmh_argmin[i];
			visited.clear();
			out.clear();

			// 对每一维：从该维倒排中取 code=sk[d] 的 posting list（升序），只对 j>i 计票
			for (int d = 0; d < S; ++d) {
				uint32_t code = sk[d];
				const auto& lst = tables[d][(size_t)code];
				if (lst.empty()) continue;

				// 二分跳到 j>i 的后缀
				auto it = std::upper_bound(lst.begin(), lst.end(), i);
				for (; it != lst.end(); ++it) add_vote(*it);
			}

			// 导出达到 need 的 j（WMH-only，直接通过阈值）
			for (int j : visited) {
				if (counts[j] >= need) out.push_back(j);
				counts[j] = 0; // 清零
			}
			std::sort(out.begin(), out.end());
			neigh[i].swap(out);

			int p = ++progress;
			if ((p % 1000) == 0) {
				double percent = 100.0 * p / N;
				std::cout << "\rWMH-only precompute (band=1): " << p << "/" << N
					<< " (" << percent << "%)" << std::flush;
			}
		}
	}
	std::cout << "\n";
}

// ======= 示例：整体管线 =======
// 1) 生成 WMH sketches
// 2) 用 argmin_code 构建每维的倒排表（vector 结构，索引即 k-mer code）
// 3) WMH-only 预计算邻接（i<j）
// 后续可直接进入你已有的“顺序提交（贪心增量）”阶段

struct WMHParams {
	int S = 256;                 // sketch 维度
	uint64_t seed = 0xC0FFEE12345678ull;
	int kmer_size = 5;
	int MAX_UCA_v = 21;            // 氨基酸字母表
};

static inline size_t safe_pow_size_t(size_t a, int b){
	size_t r = 1;
	for (int i=0;i<b;++i){
		if (r > std::numeric_limits<size_t>::max()/a) {
			throw std::runtime_error("table_size overflow");
		}
		r *= a;
	}
	return r;
}

void build_graph_wmh_only_tables(
		std::vector<Sequence>& seqs,
		double tau,
		const WMHParams& P,
		std::vector<std::vector<int>>& neigh  // 输出：neigh[i] 存 j>i
		){
	const int N = (int)seqs.size();
	std::cerr << "[WMH] N="<<N<<", S="<<P.S<<", k="<<P.kmer_size<<", MAX_UCA="<<P.MAX_UCA_v<<"\n";

	double t1 = get_time();

	// 1) 生成 WMH sketches
	std::vector<std::vector<uint32_t>> wmh_argmin;
	build_wmh_sketches(seqs, P.kmer_size, P.S, P.seed, wmh_argmin);

	// 2) 构建每维倒排表（vector 结构；索引= k-mer code）
	size_t table_size = safe_pow_size_t((size_t)P.MAX_UCA_v, P.kmer_size);
	WMHWordTables tables;
	build_wmh_word_tables(wmh_argmin, P.S, table_size, tables);
	
	double t2 = get_time();

	// 3) WMH-only 预计算邻接
	precompute_edges_wmh_only_tables(wmh_argmin, tables, tau, neigh);

	double t3 = get_time();

	std::cerr << "[WMH] done. neigh built for i<j.\n";
	cerr << "Building WMinhash time: " << t2 - t1 << " seconds" << endl;
	cerr << "Compute WMH jaccard time: " << t3 - t2 << " seconds" << endl;
}

int main(int argc, char* argv[])
{
	gzFile fp1;
	kseq_t *ks1;

	CLI::App app{"build word table v.0.0.1, build word table from fasta protein files"};
	string filename = "";
	string res_file = "";
	int64_t max_num_seqs = 0;
	int kmer_size = 5;

	auto option_input = app.add_option("-i, --input", filename, "input file name, fasta or gziped fasta formats");
	auto option_max_n = app.add_option("-n, --max_num_seqs", max_num_seqs, "max number of seqs for building word table");
	auto option_k = app.add_option("-k, --kmer_size", kmer_size, "set kmer size, default 5.");
	auto option_output = app.add_option("-o, --output", res_file, "output files");

	option_input->required();

	CLI11_PARSE(app, argc, argv);

	fp1 = gzopen(filename.c_str(),"r");

	InitNAA(MAX_UCA);
	init_aa_map();
	
	WMHParams P; P.S = 100; P.kmer_size = 5; P.MAX_UCA_v = 21;

	if(NULL == fp1){
		cerr << "Fail to open file: " << filename << endl;
		return 0;
	}

	ks1 = kseq_init(fp1);

	vector<Sequence> seqs;
	int64_t pos = 0;
	int64_t number_seqs = 0;
	int max_seq_len = 0;

	//read fasta to seqs
	while(1)
	{
		int length = kseq_read(ks1);
		if(length < 0) break;

		max_seq_len = max<int>(max_seq_len, ks1->seq.l);

		Sequence seq;
		seq.name = ks1->name.s;
		seq.comment = ks1->comment.s;
		seq.seq = ks1->seq.s;
		char * seq_data = seq.seq.data();
		for(int i = 0; i < seq.seq.size(); i++)
		{
			seq_data[i] = aa_map[seq_data[i]];
		}
		seqs.emplace_back(seq);

		number_seqs++;

		if(max_num_seqs > 0 && number_seqs >= max_num_seqs) break;
	}

	cerr << "number of sequences: " << seqs.size() << endl;
	cerr << "max seq length:" << max_seq_len << endl;
	if(seqs.size() < max_num_seqs) 
		cerr << "No more seqs than required! Total: "<< seqs.size() << "\tRequired: " << max_num_seqs <<endl;
	double t3 = get_time();
	double tau = 0.36;
	std::vector<std::vector<int>> neigh;
	build_graph_wmh_only_tables(seqs, tau, P, neigh);


	// 阶段B：顺序提交（按长度降序），得到 parent/redundant/centers
	double t4 = get_time();
	auto gres = greedy_commit_by_order_already_sorted(seqs, neigh, kmer_size);
	double t5 = get_time();

	std::cerr << "Greedy (two-stage) centers: " << gres.centers << "\n";
	std::cerr << "Two-stage jaccard filtering time: " << (t4 - t3) << " s\n";
	std::cerr << "Two-stage serial commit time: " << (t5 - t4) << " s\n";

	gzclose(fp1);	

	return 0;
}
