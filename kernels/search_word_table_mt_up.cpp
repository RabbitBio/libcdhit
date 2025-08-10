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
#define MAX_UCA 21
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
// 只统计 target_id < qid 的部分
int CountWordsUp(int aan_no,
		const std::vector<int>& word_encodes,
		const std::vector<int>& word_encodes_no,
		robin_hood::unordered_map<int,int>& lookCounts,
		const std::vector<std::vector<std::pair<int,int>>>& word_table,
		int min_rest,
		int qid)
{
	lookCounts.clear();
	if (aan_no <= 0) return 0;

	for (int j0 = 0; j0 < aan_no; ++j0) {
		int bucket = word_encodes[j0];
		int qcnt   = word_encodes_no[j0];
		if (qcnt == 0) continue;

		const auto& hits = word_table[bucket]; // 行内按 seq_id 降序
		// 找到第一个 <= qid 的位置（降序比较）
		auto it = std::lower_bound(
				hits.begin(), hits.end(), qid,
				[](const std::pair<int,int>& a, int val){ return a.first > val; }
				);
		int rest = aan_no - j0 + 1;

		// 只遍历后缀：target_id <= qid
		for (; it != hits.end(); ++it) {
			int tid  = it->first;
			if (tid == qid) continue; // 排除自身
			int tcnt = it->second;
			if (min_rest > 0 && rest < min_rest && !lookCounts.count(tid)) continue;
			int add = (qcnt < tcnt) ? qcnt : tcnt;
			lookCounts[tid] += add;
		}
	}
	return 0;
}
struct ClassifyUpResult {
    std::vector<int> parent; // parent[i] = 代表ID；无代表则 -1
    int centers = 0;
};

static inline double jaccard_CAB(int C, int A, int B) {
    int denom = A + B - C;
    return denom > 0 ? double(C) / double(denom) : 0.0;
}

ClassifyUpResult classify_up_parallel(
		const std::vector<Sequence>& seqs,
		const std::vector<std::vector<std::pair<int,int>>>& word_table,
		int kmer_size, double tau)
{
	const int N = (int)seqs.size();
	// 预计算 A[i] = |k-mers(i)|
	std::vector<int> A(N);
	for (int i = 0; i < N; ++i) {
		int L = (int)seqs[i].seq.size();
		A[i] = std::max(0, L - kmer_size + 1);
	}

	size_t max_len = 0;
	for (auto& s : seqs) max_len = std::max(max_len, s.seq.size());

	std::vector<int> parent(N, -1);
	std::atomic<int> progress{0};

	double t0 = get_time();
#pragma omp parallel
	{
		std::vector<int> wenc(max_len);
		std::vector<int> wnum(max_len);
		robin_hood::unordered_map<int,int> local; // tid -> C

#pragma omp for schedule(static,128)
		for (int qid = 0; qid < N; ++qid) {
			if (A[qid] <= 0) { parent[qid] = -1; ++progress; continue; }

			// 编码
			EncodeWords((Sequence&)seqs[qid], wenc, wnum, kmer_size);

			// 只向上（id < qid）计数
			local.clear();
			CountWordsUp(A[qid], wenc, wnum, local, word_table, /*min_rest=*/0, qid);

			// 选代表：Jaccard 最大且 >= tau（也可改成“先过阈值的最早 id”）
			int best_id = -1; double best_sim = -1.0;
			for (auto& kv : local) {
				int tid = kv.first;
				if (tid >= qid || A[tid] <= 0) continue; // 只看更长
				int C = kv.second;
				double jac = jaccard_CAB(C, A[qid], A[tid]);
				if (jac >= tau && jac > best_sim) {
					best_sim = jac;
					best_id = tid;
				}
			}
			parent[qid] = best_id; // -1 表示它是中心

			int p = ++progress;
			if ((p & 1023) == 0) {
				double percent = 100.0 * p / N;
				std::cout << "\rClassify-up progress: " << p << "/" << N
					<< " (" << percent << "%)" << std::flush;
			}
		}
	}
	std::cout << "\n";
	double t1 = get_time();

	int centers = 0;
	for (int i = 0; i < N; ++i) if (parent[i] == -1) ++centers;

	std::cerr << "Classify-up time: " << (t1 - t0) << " s\n";
	return { std::move(parent), centers };
}

// 迭代版 find_root + 路径压缩（保证 parent[x] < x 或 -1）
int find_root_compress(std::vector<int>& parent, int x) {
	int r = x;
	while (parent[r] != -1) r = parent[r];   // 找根
	// 路径压缩
	while (x != r) {
		int p = parent[x];
		parent[x] = r;
		x = p;
	}
	return r;
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


	/// build word table
	/// step 1: encode 
	/// step 2: add to wordtable
	
	vector<int> word_encodes;
	word_encodes.resize(max_seq_len);	
	vector<int> word_encodes_no;
	word_encodes_no.resize(max_seq_len);	

	/// init word table
	/// vector<vector<pair<seq_id, count>>>
	int table_size = 1;
	for(int i = 0; i < kmer_size; i++) table_size *= MAX_UCA;
	vector<vector<pair<int, int>>> word_table(table_size);

	cerr << "buckets size: " << table_size << endl;

	/// sort seq by len
    sort(seqs.begin(), seqs.end(),
         [](const Sequence& a, const Sequence& b) {
             return a.seq.size() > b.seq.size();
         });

	int nthreads = 1;
    #pragma omp parallel
    {
        #pragma omp single
        nthreads = omp_get_num_threads();
    }

    vector<vector<vector<pair<int,int>>>> local_tables(
        nthreads, vector<vector<pair<int,int>>>(table_size)
    );

	double t1 = get_time();
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        auto& local = local_tables[tid];

        vector<int> word_encodes(max_seq_len);
        vector<int> word_encodes_no(max_seq_len);

        #pragma omp for schedule(dynamic,1)
        for (int seq_id = 0; seq_id < (int)seqs.size(); ++seq_id) {
            auto& s = seqs[seq_id];
            if ((int)s.seq.size() < kmer_size) continue;

            EncodeWords(s, word_encodes, word_encodes_no, kmer_size);

            int kmer_no = (int)s.seq.size() - kmer_size + 1;
            for (int j = 0; j < kmer_no; ++j) {
                int bucket = word_encodes[j];
                int count  = word_encodes_no[j];
                if (count > 0) {
                    local[bucket].emplace_back(seq_id, count); 
                }
            }
        }
    }

	#pragma omp parallel for schedule(static)
	for (long long b = 0; b < (long long)table_size; ++b) {
		// 预估容量，减少 realloc
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

	//sort each row in word table in ascending order
	#pragma omp parallel for schedule(dynamic)
	for (size_t i = 0; i < word_table.size(); ++i) {
		auto& row = word_table[i];
		std::sort(row.begin(), row.end(),
				[](const std::pair<int,int>& a, const std::pair<int,int>& b) {
				return a.first < b.first; // seq_id 升序
				});
	}
	double t2 = get_time();
 
	int element_no = 0;
	int empty_no = 0;
	for(int i = 0; i < word_table.size(); i++)
	{
		element_no += word_table[i].size();
		if(word_table[i].size() == 0) empty_no++;
	}
	
	cerr << "Elements in word table: " << element_no << endl;
	cerr << "Empty buckets in word table: " << empty_no << endl;
	cerr << "Encode and insert word table time: " << t2 - t1 << " seconds" << endl;

	//search all-vs-all
	//unordered_map<int, int> lookCounts;
	//vector<int> indexMapping(seqs.size());
	//int min_rest = 0; //FIXME: using real min_rest
	//int64_t all_hits = 0;
	//double t3 = get_time();
	//for(int i = 0; i < seqs.size(); i++)
	//{
    //    EncodeWords(seqs[i], word_encodes, word_encodes_no, kmer_size);
	//	int kmer_no = seqs[i].seq.size() - kmer_size + 1;
	//	CountWords(kmer_no,
    //           word_encodes,
    //           word_encodes_no,
    //           lookCounts,      // 输出：候选目标及公共k-mer数
    //           word_table,
	//		   min_rest);
	//
	//	all_hits += lookCounts.size();

	//	if (i % 100 == 0) {
	//		double percent = 100.0 * i / seqs.size();
	//		std::cout << "\rProgress: " << i << "/" << seqs.size()
	//			<< " (" << percent << "%)" << std::flush;
	//	}

	//}
	//cerr << endl;
	std::atomic<int> progress{0};
	long long all_hits = 0;  // 归约变量
	const int min_rest = 0;  // FIXME: 按需设置

//	double t3 = get_time();
//#pragma omp parallel
//	{
//		// 线程私有缓存
//		vector<int> word_encodes(max_seq_len);
//		vector<int> word_encodes_no(max_seq_len);
//		robin_hood::unordered_map<int,int> lookCounts;   // 线程私有，避免竞争
//		//lookCounts.reserve(1 << 12);         // 视数据量调整
//
//		long long local_hits = 0; // 每线程局部计数，最后归并
//
//#pragma omp for schedule(dynamic, 1) nowait
//		for (int i = 0; i < seqs.size(); ++i) {
//			// 编码
//			EncodeWords(seqs[i], word_encodes, word_encodes_no, kmer_size);
//			int kmer_no = (int)seqs[i].seq.size() - kmer_size + 1;
//
//			// 搜索计数（写入线程私有 lookCounts）
//			lookCounts.clear();
//			CountWords(kmer_no, word_encodes, word_encodes_no,
//					lookCounts, word_table, min_rest);
//
//			local_hits += (long long)lookCounts.size();
//
//			// 进度（每处理100条打印）
//			int p = ++progress;
//			if (p % 100 == 0) {
//				double percent = 100.0 * p / seqs.size();
//				std::cout << "\rProgress: " << p << "/" << seqs.size()
//					<< " (" << percent << "%)" << std::flush;
//			}
//		}
//
//		// 归并 all_hits
//#pragma omp atomic
//		all_hits += local_hits;
//	}
//
//	std::cerr << std::endl;
//	double t4 = get_time();
//
//	std::cerr << "All vs all kmer hits: " << all_hits << '\n';
//	std::cerr << "CountWords time: " << (t4 - t3) << " seconds" << std::endl;

// ===== 两阶段：A) 并行预计算稀疏图(按阈值存边)  B) 顺序提交（贪心） =====
	double tau = 0.5;   // Jaccard 阈值：按需设置

	// 阶段A：预计算邻接（只存 i<j）
	double t3 = get_time();
	auto cres = classify_up_parallel(seqs, word_table, kmer_size, tau);
	// 阶段B：顺序提交（按长度降序），得到 parent/redundant/centers
	double t4 = get_time();
	// 一次性并行压缩到最终中心
	auto & parent = cres.parent;
#pragma omp parallel for schedule(static, 1024)
	for (int i = 0; i < (int)parent.size(); ++i) {
		if (parent[i] != -1)
			(void)find_root_compress(parent, i);
	}

	// 压缩完成后：parent[i]==-1 的是中心；否则 parent[i] 就是最终代表 id
	int centers = 0;
	for (int i = 0; i < (int)parent.size(); ++i)
		if (parent[i] == -1) ++centers;	
	double t5 = get_time();

	std::cerr << "Centers (classify-up): " << centers << "\n";
	std::cerr << "Classify-up time: " << (t4 - t3) << " s\n";
	std::cerr << "Root compress time: " << (t5 - t4) << " s\n";

	gzclose(fp1);	

	return 0;
}
