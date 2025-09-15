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

#include <immintrin.h>

using namespace std;

// #define __AVX2__

double sort_time = 0.0;

// 时间函数
static double get_time() {
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

void cluster_sequences_st(
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


// EncodeWord for std::pair format,instead of two array
// by mgl
void EncodeWordsPair(const Sequence_new &seq,std::vector<pair<int,int>>& word_encodes,int NAA){
    const char* seqi = seq.data;
    int len = strlen(seqi);
    int aan_no = len-NAA+1;
    int j;

    unsigned char k,k1;
    for(j=0;j<aan_no;j++){
        const char* word = seqi+j;
        int encode = 0;
        for (k = 0, k1 = NAA - 1; k < NAA; k++, k1--) {
            encode += aa_map[(unsigned char)word[k]] * NAAN_array[k1];  // 修改为使用 aa_map
        }
        word_encodes[j].first = encode;
        word_encodes[j].second = 1;
    }
    std::sort(word_encodes.begin(),word_encodes.end(),[](const std::pair<int,int>&a, const std::pair<int,int>&b){return a.first<b.first;});
    // cout<<"sort finished"<<endl;
    int w = -1;
    int i = 0;
    int current_key = -1;
    for(i=0;i<aan_no;i++){
        if(current_key != word_encodes[i].first){
            current_key = word_encodes[i].first;
            w++;
            word_encodes[w].first = word_encodes[i].first;
            word_encodes[w].second = 0;
        }
        word_encodes[w].second++;
    }
    word_encodes.resize(w+1);
}


// template <class T, class U, class Compare = std::less<T>>
// void SortByKey(std::vector<T>& keys, std::vector<U>& vals, Compare comp = Compare{}) {
//     const size_t n = keys.size();
//     if (vals.size() != n) return;
//     if (n <= 1) return;

//     std::vector<std::pair<T,U>> items;
//     items.reserve(n);
//     for (size_t i = 0; i < n; ++i) {
//         items.emplace_back(std::move(keys[i]), std::move(vals[i]));
//     }
//     std::sort(items.begin(), items.end(),
//               [&](const auto& a, const auto& b){ return comp(a.first, b.first); });

//     for (size_t i = 0; i < n; ++i) {
//         keys[i] = std::move(items[i].first);
//         vals[i] = std::move(items[i].second);
//     }
// }

void CountWeightedJaccard(
    const std::vector<pair<int,int>>& seqi,
    const std::vector<pair<int,int>>& seqj,
    double& jac
){
    size_t idx_i = 0,idx_j = 0;
    const size_t len_i = seqi.size(),len_j = seqj.size();
    long long inter_val = 0,union_val = 0;
    while(idx_i<len_i && idx_j<len_j){
        if(seqi[idx_i].first==seqj[idx_j].first){
            const int a = seqi[idx_i].second;
            const int b = seqj[idx_j].second;
            inter_val+=std::min(a,b);
            union_val+=a+b;
            idx_i++,idx_j++;
        }
        else if(seqi[idx_i].first<seqj[idx_j].second){
            union_val+=seqi[idx_i].second;
            idx_i++;
        }
        else{
            union_val+=seqj[idx_j].second;
            idx_j++;
        }
    }
    while(idx_i<len_i) union_val+=seqi[idx_i++].second;
    while(idx_j<len_j) union_val+=seqj[idx_j++].second;
    union_val -= inter_val;
    jac = (union_val!=0)? static_cast<double>(inter_val) / static_cast<double>(union_val) : 0.0;
}

void CountWeightedJaccard_SoA(
    const std::vector<int>& seqi,
    const std::vector<int>& counti,
    const std::vector<int>& seqj,
    const std::vector<int>& countj,
    double& jac
){
    size_t idx_i = 0,idx_j = 0;
    const size_t len_i = seqi.size(),len_j = seqj.size();
    long long inter_val = 0,union_val = 0;
    while(idx_i<len_i && idx_j<len_j){
        if(seqi[idx_i]==seqj[idx_j]){
            const int a = counti[idx_i];
            const int b = countj[idx_j];
            inter_val+=std::min(a,b);
            union_val+=a+b;
            idx_i++,idx_j++;
        }
        else if(seqi[idx_i]<seqj[idx_j]){
            union_val+=counti[idx_i];
            idx_i++;
        }
        else{
            union_val+=countj[idx_j];
            idx_j++;
        }
    }
    while(idx_i<len_i) union_val+=counti[idx_i++];
    while(idx_j<len_j) union_val+=countj[idx_j++];
    union_val -= inter_val;
    jac = (union_val!=0)? static_cast<double>(inter_val) / static_cast<double>(union_val) : 0.0;
	// cout<<jac<<endl;
}

static inline int u32_weighted_jaccard_scalar_tail(
    const int* a, const int* wa, size_t na,
    const int* b, const int* wb, size_t nb
){
    size_t i = 0, j = 0;
    int inter_val = 0;

    while (i < na && j < nb) {
        int ka = a[i];
        int kb = b[j];

        if (ka == kb) {
            int fa = wa[i];
            int fb = wb[j];
            inter_val += (fa < fb) ? fa : fb;
            ++i; ++j;
        } else if (ka < kb) {
            ++i;
        } else {
            ++j;
        }
    }
    return inter_val;
}

// #define __AVX512F__
#ifdef __AVX512F__

alignas(64) static const int u32_rot_idx[16][16] = {
    // r = 0
    {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 },
    // r = 1
    { 15,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14 },
    // r = 2
    { 14, 15,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13 },
    // r = 3
    { 13, 14, 15,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12 },
    // r = 4
    { 12, 13, 14, 15,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11 },
    // r = 5
    { 11, 12, 13, 14, 15,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10 },
    // r = 6
    { 10, 11, 12, 13, 14, 15,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9 },
    // r = 7
    {  9, 10, 11, 12, 13, 14, 15,  0,  1,  2,  3,  4,  5,  6,  7,  8 },
    // r = 8
    {  8,  9, 10, 11, 12, 13, 14, 15,  0,  1,  2,  3,  4,  5,  6,  7 },
    // r = 9
    {  7,  8,  9, 10, 11, 12, 13, 14, 15,  0,  1,  2,  3,  4,  5,  6 },
    // r = 10
    {  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,  0,  1,  2,  3,  4,  5 },
    // r = 11
    {  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,  0,  1,  2,  3,  4 },
    // r = 12
    {  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,  0,  1,  2,  3 },
    // r = 13
    {  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,  0,  1,  2 },
    // r = 14
    {  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,  0,  1 },
    // r = 15
    {  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,  0 },
};

void u32_WeightedJaccard_vetcor_AVX512(
    const vector<int>& keysA,const vector<int>& freqA,
    const vector<int>& keysB,const vector<int>& freqB,
    double& jac
){
    const size_t n = keysA.size(),m = keysB.size();
    if (n == 0 && m == 0) { jac = 0.0; return; }
    assert(n == freqA.size() && m == freqB.size());
    int union_val = 0;
    for(int v : freqA) union_val+=v;
    for(int v : freqB) union_val+=v;
    const int *ka = keysA.data();
    const int *kb = keysB.data();
    const int *fa = freqA.data();
    const int *fb = freqB.data();
    
    size_t ia = 0, ib = 0;
    const size_t st_a = (n/16)*16;
    const size_t st_b = (m/16)*16;

    int inter_val;
    if(n<16||m<16){
        inter_val = u32_weighted_jaccard_scalar_tail(ka,fa,n,kb,fb,m);
        union_val -= inter_val;
        jac = (union_val==0)? 0.0 : (double)inter_val/(double)union_val;
        return;
    }

    alignas(64) int tmp_store[16];
    while(ia<st_a && ib<st_b){
        __m512i va_k = _mm512_loadu_si512((const __m512i*)(ka+ia));
        __m512i va_f = _mm512_loadu_si512((const __m512i*)(fa+ia));
        __m512i vb_k = _mm512_loadu_si512((const __m512i*)(kb+ib));
        __m512i vb_f = _mm512_loadu_si512((const __m512i*)(fb+ib));

        int a_max = ka[ia+15];
        int b_max = kb[ib+15];
        if(a_max<=b_max) ia+=16;
        if(b_max<=a_max) ib+=16;

        for (int r = 0; r < 16; ++r) {
            const __m512i rot_idx = _mm512_load_si512((const __m512i*)u32_rot_idx[r]);

            __m512i vb_k_rot = _mm512_permutexvar_epi32(rot_idx, vb_k);
            __mmask16 mcmp   = _mm512_cmpeq_epi32_mask(va_k, vb_k_rot);
            if (mcmp == 0) continue;

            __m512i vb_f_rot = _mm512_permutexvar_epi32(rot_idx, vb_f);
            __m512i vminf    = _mm512_min_epi32(va_f, vb_f_rot);

            __m512i vsel = _mm512_maskz_mov_epi32(mcmp, vminf);
            _mm512_store_si512((__m512i*)tmp_store, vsel);

            long long lane_sum = 0;
            for (int t = 0; t < 16; ++t) lane_sum += tmp_store[t];
            inter_val += lane_sum;
        }
    }
    if (ia < n && ib < m) {
        inter_val += u32_weighted_jaccard_scalar_tail(
            ka + ia, fa + ia, n - ia,
            kb + ib, fb + ib, m - ib
        );
    }

    union_val -= inter_val;
    jac = (union_val == 0) ? 0.0 : (double)inter_val / (double)union_val;
}

#endif

#ifdef __AVX2__

// 对应 8 种 v_b 变换下的“b-lane 映射”表：给定 a 的 lane（0..7），返回匹配到的 b 的 lane 下标。
// 变换定义与下面代码的 rot 顺序对应：
// 0: identity
// 1: within-128 right rotate by 1 (permute_ps SHUFFLE(0,3,2,1))
// 2: within-128 "between" (SHUFFLE(1,0,3,2))
// 3: within-128 left rotate by 1 (SHUFFLE(2,1,0,3))
// 4: swap 128-bit halves
// 5: swap + right rotate
// 6: swap + between
// 7: swap + left rotate
static constexpr uint8_t B_LANE_MAP[8][8] = {
    // 0: identity
    {0,1,2,3,4,5,6,7},
    // 1: right rot within each 128b half: [3,0,1,2 | 7,4,5,6]
    {3,0,1,2,7,4,5,6},
    // 2: between (1,0,3,2 | 5,4,7,6)
    {1,0,3,2,5,4,7,6},
    // 3: left rot within each half: [1,2,3,0 | 5,6,7,4]
    {1,2,3,0,5,6,7,4},
    // 4: swap halves: [4,5,6,7 | 0,1,2,3]
    {4,5,6,7,0,1,2,3},
    // 5: swap + right rot: [7,4,5,6 | 3,0,1,2]
    {7,4,5,6,3,0,1,2},
    // 6: swap + between: [5,4,7,6 | 1,0,3,2]
    {5,4,7,6,1,0,3,2},
    // 7: swap + left rot: [5,6,7,4 | 1,2,3,0]
    {5,6,7,4,1,2,3,0},
};

void u32_WeightedJaccard_vector_AVX2(
    const vector<int>& keysA,const vector<int>& freqA,
    const vector<int>& keysB,const vector<int>& freqB,
    double& jac
){
	// cout<<"use avx2"<<endl;
    const size_t n = keysA.size(),m = keysB.size();
    if(n == 0&& m==0){
        jac = 0.0;
		// cout<<jac<<endl;
        return;
    }
    assert(n==freqA.size()&&m==freqB.size());
    int union_val = 0;
    for(auto v: freqA) union_val+=v;
    for(auto v: freqB) union_val+=v;

    const int* ka = keysA.data();
    const int* kb = keysB.data();
    const int* fa = freqA.data();
    const int* fb = freqB.data();
    size_t ia = 0, ib = 0;
    const size_t st_a = (n/8)*8;
    const size_t st_b = (m/8)*8;
    
    int inter_val = 0;
    if(n<8||m<8){   
        inter_val = u32_weighted_jaccard_scalar_tail(ka,fa,n,kb,fb,m);
        union_val -= inter_val;
        jac = (union_val == 0)? 0.0 : (double)inter_val/(double) union_val;
		// cout<<jac<<endl;
        return;
    }
    
    const int32_t sh_right = _MM_SHUFFLE(0,3,2,1);
    const int32_t sh_left  = _MM_SHUFFLE(2,1,0,3);
    const int32_t sh_between = _MM_SHUFFLE(1,0,3,2);

    while(ia<st_a && ib<st_b){
        __m256i v_a = _mm256_loadu_si256((const __m256i*)(ka + ia));
        __m256i v_b = _mm256_loadu_si256((const __m256i*)(kb + ib));

        int a_max = ka[ia+7];
        int b_max = kb[ib+7];
        if(a_max<=b_max) ia+=8;
        if(b_max<=a_max) ib+=8;

        __m256 vb_ps = (__m256) v_b;
        __m256 vb_r1 = _mm256_permute_ps(vb_ps,sh_right);
        __m256 vb_bt = _mm256_permute_ps(vb_ps,sh_between);
        __m256 vb_l1 = _mm256_permute_ps(vb_ps,sh_left);

        // 128bit 半宽

        __m256 vb_sw = _mm256_permute2f128_ps(vb_ps,vb_ps,1);
        __m256 vb_sw_r1 = _mm256_permute_ps(vb_sw,sh_right);
        __m256 vb_sw_bt = _mm256_permute_ps(vb_sw,sh_between);
        __m256 vb_sw_l1 = _mm256_permute_ps(vb_sw,sh_left);

        __m256i cmp0 = _mm256_cmpeq_epi32(v_a,v_b);
        __m256i cmp1 = _mm256_cmpeq_epi32(v_a,(__m256i)vb_r1);
        __m256i cmp2 = _mm256_cmpeq_epi32(v_a,(__m256i)vb_bt);
        __m256i cmp3 = _mm256_cmpeq_epi32(v_a,(__m256i)vb_l1);
        __m256i cmp4 = _mm256_cmpeq_epi32(v_a,(__m256i)vb_sw);
        __m256i cmp5 = _mm256_cmpeq_epi32(v_a,(__m256i)vb_sw_r1);
        __m256i cmp6 = _mm256_cmpeq_epi32(v_a,(__m256i)vb_sw_bt);
        __m256i cmp7 = _mm256_cmpeq_epi32(v_a,(__m256i)vb_sw_l1);

        const __m256i cmps[8] = {cmp0,cmp1,cmp2,cmp3,cmp4,cmp5,cmp6,cmp7};

        for(int p=0;p<8;p++){
            uint32_t mask = (uint32_t)_mm256_movemask_ps((__m256)cmps[p]);
            while(mask){
#if defined(__BMI__)
                uint32_t laneA = (uint32_t)_tzcnt_u32(mask);
#else
                uint32_t laneA = (uint32_t)__builtin_ctz(mask);
#endif
                mask &= mask - 1;
                uint32_t laneB = B_LANE_MAP[p][laneA];
                size_t idxA = (ia-8) + laneA;
                size_t idxB = (ib-8) + laneB;
                int freA = fa[idxA];
                int freB = fb[idxB];
                inter_val += (freA<freB)? freA:freB;
            }
        }
    }
    if(ia<n && ib<m){
        inter_val += u32_weighted_jaccard_scalar_tail(ka+ia,fa+ia,n-ia,kb+ib,fb+ib,m-ib);

    }

    union_val -= inter_val;
    jac = (union_val==0) ? 0.0 : (double)inter_val/(double)union_val;
	// cout<<jac<<endl;
}
#endif

void LSD_sort_8bit(vector<int>& word_encodes,vector<int>& word_encodes_no){
    // double t1 = get_time();
    const size_t n = word_encodes.size();
    if(n<=1){
        word_encodes_no[0]=1;
        return;
    }
    vector<int> buff(n);
    vector<int> cnt(256),start(256),next_pos(256);
    auto pass = [&](int shift, const vector<int>& src, vector<int>& dst) {
        // 计数
        fill(cnt.begin(), cnt.end(), 0);
        for (uint32_t v : src) {
            uint32_t d = (v >> shift) & 0xFFu;
            ++cnt[d];
        }
        // 前缀（每个桶的起始写入位置）
        size_t run = 0;
        for (int d = 0; d < 256; ++d) {
            start[d] = run;
            run += cnt[d];
        }
        // 稳定分配
        next_pos = start;
        for (uint32_t v : src) {
            uint32_t d = (v >> shift) & 0xFFu;
            dst[next_pos[d]++] = v;
        }
    };
    pass(0,word_encodes,buff);
    pass(8,buff,word_encodes);
    pass(16,word_encodes,buff);
    word_encodes.swap(buff);
}

void LSD_sort_11bit(vector<int>& word_encodes,vector<int>& word_encodes_no){
    const size_t n = word_encodes.size();
    if(n<=1){
        word_encodes_no[0]=1;
        return;
    }
    vector<int> buff(n);
    vector<int> cnt(2048),start(2048),next_pos(2048);
    auto pass = [&](int shift, const vector<int>& src, vector<int>& dst) {
        // 计数
        fill(cnt.begin(), cnt.end(), 0);
        for (int v : src) {
            unsigned d = (unsigned(v) >> shift) & 0x7FFu;  // 取 11 位
            ++cnt[d];
        }
        // 前缀：每个桶的起始位置
        size_t run = 0;
        for (size_t d = 0; d < 2048; ++d) {
            start[d] = run;
            run += cnt[d];
        }
        // 稳定分配
        next_pos = start;
        for (int v : src) {
            unsigned d = (unsigned(v) >> shift) & 0x7FFu;
            dst[next_pos[d]++] = v;
        }
    };
    pass(0,word_encodes,buff);
    pass(11,buff,word_encodes);
}

void std_sort(vector<int>& word_encodes,const int len){
    // double t1 = get_time();
    sort(word_encodes.begin(),word_encodes.begin()+len);
    // double t2 = get_time();
    // sort_time+=t2-t1;
    // cerr << "Sorting time: " << t2 - t1 << " s" << endl;
}

void MergeKmerFreq(vector<int>& word_encodes,vector<int>&word_encodes_no){
    int kmer = -1;
    int w = -1;
    for(int i=0;i<word_encodes.size();i++){
        if(word_encodes[i]!=kmer){
            kmer=word_encodes[i];
            w++;
            word_encodes[w]=kmer;
        }
        word_encodes_no[w]++;
    }
    word_encodes.resize(w+1);
    word_encodes_no.resize(w+1);
}

void EncodeWordsSoA(const Sequence_new &seq, vector<int> &word_encodes, vector<int> &word_encodes_no, int NAA) {
    const char* seqi = seq.data;
    int len = strlen(seqi);
    int aan_no = len - NAA + 1;
    int j;
    unsigned char k, k1;
    for (j = 0; j < aan_no; j++) {
        const char* word = seqi + j;
        int encode = 0;
        for (k = 0, k1 = NAA - 1; k < NAA; k++, k1--) {
            encode += aa_map[(unsigned char)word[k]] * NAAN_array[k1];
        }
        word_encodes[j] = encode;
    }
    word_encodes.resize(aan_no);
    // cout<<"<"<<word_encodes.size()<<", ";
    
    LSD_sort_8bit(word_encodes,word_encodes_no);
    // std_sort(word_encodes,aan_no);

    // for(int i=0;i<10;i++){
    //     cout<<"<"<<word_encodes[i]<<"> ";
    // }
    // cout<<"\n";
    // std_sort(word_encodes,word_encodes_no);
    MergeKmerFreq(word_encodes,word_encodes_no);
    // cout<<word_encodes.size()<<"> "<<endl;
    // for(int i=0;i<5;i++){
    //     cout<<"<"<<word_encodes[i]<<", "<<word_encodes_no[i]<<"> ";
    // }
    // cout<<"\n";
}

void cluster_sequences_st_less10(
        std::vector<Sequence_new>& seqs,
        std::vector<int>& parent,
        int kmer_size,
        double tau
){
    InitNAA(MAX_UAA);
    init_aa_map();
    int N=(int)seqs.size();

    int max_seq_len = 0;
    sort(seqs.begin(), seqs.end(),
        [](const Sequence_new& a, const Sequence_new& b) {
        return strlen(a.data) > strlen(b.data);
        });
    max_seq_len = strlen(seqs[0].data);
    
    // vector<vector<pair<int,int>>> word_encodes(N);
    vector<vector<int>> word_encodes(N);
    vector<vector<int>> word_encodes_no(N);
    
    for (int seq_id = 0; seq_id < N; ++seq_id) {
        word_encodes[seq_id].resize(max_seq_len);
        word_encodes_no[seq_id].resize(max_seq_len);
        // word_encodes_no[seq_id].resize(max_seq_len);
        auto& s = seqs[seq_id];
        int len = strlen(s.data);
        if (len < kmer_size) continue;
        // EncodeWordsPair(s,word_encodes[seq_id],kmer_size);
        EncodeWordsSoA(s,word_encodes[seq_id],word_encodes_no[seq_id],kmer_size);
    }

    DSU dsu(N);
    double jac = 0.0;
    for(int seq_i=0;seq_i<N;seq_i++){
        for(int seq_j=seq_i+1;seq_j<N;seq_j++){
            if(dsu.find(seq_i)==dsu.find(seq_j)) continue;
			// cout<<"<"<<seq_i<<","<<seq_j<<">\t";
#if defined(__AVX512F__)
            // AVX-512 编译：用刚刚写的 AVX512 版本
            u32_WeightedJaccard_vetcor_AVX512(
                word_encodes[seq_i], word_encodes_no[seq_i],
                word_encodes[seq_j], word_encodes_no[seq_j],
                jac
            );
#elif defined(__AVX2__)
            // AVX2 编译：用已有的 AVX2 版本
            u32_WeightedJaccard_vector_AVX2(
                word_encodes[seq_i], word_encodes_no[seq_i],
                word_encodes[seq_j], word_encodes_no[seq_j],
                jac
            );
#else
            // 既没有 AVX-512 也没有 AVX2：用默认标量/SoA 版本
            CountWeightedJaccard_SoA(
                word_encodes[seq_i], word_encodes_no[seq_i],
                word_encodes[seq_j], word_encodes_no[seq_j],
                jac
            );
#endif

            // CountWeightedJaccard_SoA(word_encodes[seq_i],word_encodes_no[seq_i],word_encodes[seq_j],word_encodes_no[seq_j],jac);
            // cout<<jac<<endl;
            if(jac>=tau){
                dsu.unite(seq_i,seq_j);
            }
        }
    }
    
    for (int i = 0; i < N; ++i) {
        parent[seqs[i].seq_id] = seqs[dsu.find(i)].seq_id;
    }
    // cerr<<"Sorting total time: "<<sort_time<<" s\n";
}

void std_sort_cArray(int* word_encodes,const int len){
    // double t1 = get_time();
    sort(word_encodes,word_encodes+len);
    // double t2 = get_time();
    // sort_time+=t2-t1;
    // cerr << "Sorting time: " << t2 - t1 << " s" << endl;
}

void LSD_sort_8bit_cArray(int* word_encodes,const int n){
    // double t1 = get_time();
    // const size_t n = word_encodes.size();
    if(n<=1){
        return;
    }

    int* buff = new int[n];
    // 辅助数组
    int cnt[256];
    int start[256];
    int next_pos[256];

    auto pass = [&](int shift, int* src, int* dst) {
        // 计数
        memset(cnt, 0, sizeof(cnt));
        for (int i = 0; i < n; i++) {
            uint32_t v = (uint32_t)src[i];
            uint32_t d = (v >> shift) & 0xFFu;
            ++cnt[d];
        }
        // 前缀和：桶起始位置
        int run = 0;
        for (int d = 0; d < 256; d++) {
            start[d] = run;
            run += cnt[d];
        }
        // 稳定分配
        memcpy(next_pos, start, sizeof(start));
        for (int i = 0; i < n; i++) {
            uint32_t v = (uint32_t)src[i];
            uint32_t d = (v >> shift) & 0xFFu;
            dst[next_pos[d]++] = v;
        }
    };

    pass(0, word_encodes, buff);
    pass(8, buff, word_encodes);
    pass(16, word_encodes, buff);

    // 拷回结果
    for (int i = 0; i < n; i++) {
        word_encodes[i] = buff[i];
    }
    delete[] buff;
}

void MergeKmerFreq_cArray(int* word_encodes,int* word_encodes_no,int& len){
    int kmer = -1;
    int w = -1;
    for(int i=0;i<len;i++){
        if(word_encodes[i]!=kmer){
            kmer=word_encodes[i];
            w++;
            word_encodes[w]=kmer;
        }
        word_encodes_no[w]++;
    }
    len = w+1;
}

void EncodeWordsSoA_cArray(
    const Sequence_new &seq,
    int* word_encodes,
    int* word_encodes_no,
    int& len,
    int NAA
){
    const char* seqi = seq.data;
    len = strlen(seqi) - NAA + 1;
    int j;
    unsigned char k, k1;
    for (j = 0; j < len; j++) {
        const char* word = seqi + j;
        int encode = 0;
        for (k = 0, k1 = NAA - 1; k < NAA; k++, k1--) {
            encode += aa_map[(unsigned char)word[k]] * NAAN_array[k1];
        }
        word_encodes[j] = encode;
    }
    // cout<<"<"<<len<<", "<<endl;
    // std_sort_cArray(word_encodes,len);
    LSD_sort_8bit_cArray(word_encodes,len);
    // for(int i=0;i<10;i++){
    //     cout<<"<"<<word_encodes[i]<<"> ";
    // }
    // cout<<"\n";
    MergeKmerFreq_cArray(word_encodes,word_encodes_no,len);
    // cout<<len<<"> "<<endl;
}

void CountWeightedJaccard_SoA_cArray(
    const int* seqi,const int* counti,const int len_i,
    const int* seqj,const int* countj,const int len_j,
    double& jac
){
    size_t idx_i = 0,idx_j = 0;
    long long inter_val = 0,union_val = 0;
    while(idx_i<len_i && idx_j<len_j){
        if(seqi[idx_i]==seqj[idx_j]){
            const int a = counti[idx_i];
            const int b = countj[idx_j];
            inter_val+=std::min(a,b);
            union_val+=a+b;
            idx_i++,idx_j++;
        }
        else if(seqi[idx_i]<seqj[idx_j]){
            union_val+=counti[idx_i];
            idx_i++;
        }
        else{
            union_val+=countj[idx_j];
            idx_j++;
        }
    }
    while(idx_i<len_i) union_val+=counti[idx_i++];
    while(idx_j<len_j) union_val+=countj[idx_j++];
    union_val -= inter_val;
    // cout<<" "<<inter_val<<" , "<< union_val<<"\t";
    jac = (union_val!=0)? static_cast<double>(inter_val) / static_cast<double>(union_val) : 0.0;
}

#ifdef __AVX2__

void u32_WeightedJaccard_vector_AVX2_cArray(
    const int* ka,const int* fa,const int n,
    const int* kb,const int* fb,const int m,
    double& jac
){
	// cout<<"use avx2"<<endl;
    if(n == 0&& m==0){
        jac = 0.0;
		// cout<<jac<<endl;
        return;
    }

    int union_val = 0;
    // for(auto v: freqA) union_val+=v;
    // for(auto v: freqB) union_val+=v;
    for(int i=0;i<n;i++) union_val += fa[i];
    for(int i=0;i<m;i++) union_val+=fb[i];

    size_t ia = 0, ib = 0;
    const size_t st_a = (n/8)*8;
    const size_t st_b = (m/8)*8;
    
    int inter_val = 0;
    if(n<8||m<8){   
        inter_val = u32_weighted_jaccard_scalar_tail(ka,fa,n,kb,fb,m);
        union_val -= inter_val;
        jac = (union_val == 0)? 0.0 : (double)inter_val/(double) union_val;
		// cout<<jac<<endl;
        return;
    }
    
    const int32_t sh_right = _MM_SHUFFLE(0,3,2,1);
    const int32_t sh_left  = _MM_SHUFFLE(2,1,0,3);
    const int32_t sh_between = _MM_SHUFFLE(1,0,3,2);

    while(ia<st_a && ib<st_b){
        __m256i v_a = _mm256_loadu_si256((const __m256i*)(ka + ia));
        __m256i v_b = _mm256_loadu_si256((const __m256i*)(kb + ib));

        int a_max = ka[ia+7];
        int b_max = kb[ib+7];
        if(a_max<=b_max) ia+=8;
        if(b_max<=a_max) ib+=8;

        __m256 vb_ps = (__m256) v_b;
        __m256 vb_r1 = _mm256_permute_ps(vb_ps,sh_right);
        __m256 vb_bt = _mm256_permute_ps(vb_ps,sh_between);
        __m256 vb_l1 = _mm256_permute_ps(vb_ps,sh_left);

        // 128bit 半宽

        __m256 vb_sw = _mm256_permute2f128_ps(vb_ps,vb_ps,1);
        __m256 vb_sw_r1 = _mm256_permute_ps(vb_sw,sh_right);
        __m256 vb_sw_bt = _mm256_permute_ps(vb_sw,sh_between);
        __m256 vb_sw_l1 = _mm256_permute_ps(vb_sw,sh_left);

        __m256i cmp0 = _mm256_cmpeq_epi32(v_a,v_b);
        __m256i cmp1 = _mm256_cmpeq_epi32(v_a,(__m256i)vb_r1);
        __m256i cmp2 = _mm256_cmpeq_epi32(v_a,(__m256i)vb_bt);
        __m256i cmp3 = _mm256_cmpeq_epi32(v_a,(__m256i)vb_l1);
        __m256i cmp4 = _mm256_cmpeq_epi32(v_a,(__m256i)vb_sw);
        __m256i cmp5 = _mm256_cmpeq_epi32(v_a,(__m256i)vb_sw_r1);
        __m256i cmp6 = _mm256_cmpeq_epi32(v_a,(__m256i)vb_sw_bt);
        __m256i cmp7 = _mm256_cmpeq_epi32(v_a,(__m256i)vb_sw_l1);

        const __m256i cmps[8] = {cmp0,cmp1,cmp2,cmp3,cmp4,cmp5,cmp6,cmp7};

        for(int p=0;p<8;p++){
            uint32_t mask = (uint32_t)_mm256_movemask_ps((__m256)cmps[p]);
            while(mask){
#if defined(__BMI__)
                uint32_t laneA = (uint32_t)_tzcnt_u32(mask);
#else
                uint32_t laneA = (uint32_t)__builtin_ctz(mask);
#endif
                mask &= mask - 1;
                uint32_t laneB = B_LANE_MAP[p][laneA];
                size_t idxA = (ia-8) + laneA;
                size_t idxB = (ib-8) + laneB;
                int freA = fa[idxA];
                int freB = fa[idxB];
                inter_val += (freA<freB)? freA:freB;
            }
        }
    }
    if(ia<n && ib<m){
        inter_val += u32_weighted_jaccard_scalar_tail(ka+ia,fa+ia,n-ia,kb+ib,fb+ib,m-ib);

    }

    union_val -= inter_val;
    jac = (union_val==0) ? 0.0 : (double)inter_val/(double)union_val;
	// cout<<jac<<endl;
}
#endif

void cluster_sequence_singleThread_smallScale_cArray(
    std::vector<Sequence_new>& seqs,
    std::vector<int>& parent,
    int kmer_size,
    double tau
){
    InitNAA(MAX_UAA);
    init_aa_map();
    int N=(int)seqs.size();

    int max_seq_len = 0;
    sort(seqs.begin(), seqs.end(),
        [](const Sequence_new& a, const Sequence_new& b) {
        return strlen(a.data) > strlen(b.data);
        });
    max_seq_len = strlen(seqs[0].data);

    int** word_encodes = new int*[N];
    int** word_encodes_no = new int*[N];
    int* buff = new int[max_seq_len];
    int* lens = new int[N];
    
    for (int seq_id = 0; seq_id < N; ++seq_id) {
        // word_encodes_no[seq_id].resize(max_seq_len);
        word_encodes[seq_id] = new int[max_seq_len];
        word_encodes_no[seq_id] = new int[max_seq_len];
        auto& s = seqs[seq_id];
        int len = strlen(s.data);
        if (len < kmer_size) continue;
        // EncodeWordsPair(s,word_encodes[seq_id],kmer_size);
        EncodeWordsSoA_cArray(s,word_encodes[seq_id],word_encodes_no[seq_id],lens[seq_id],kmer_size);
    }

    DSU dsu(N);
    double jac = 0.0;
    for(int seq_i=0;seq_i<N;seq_i++){
        for(int seq_j=seq_i+1;seq_j<N;seq_j++){
            if(dsu.find(seq_i)==dsu.find(seq_j)) continue;
			// cout<<"<"<<seq_i<<","<<seq_j<<">\t";
#ifdef __AVX2__
			u32_WeightedJaccard_vector_AVX2_cArray(word_encodes[seq_i],word_encodes_no[seq_i],lens[seq_i],
                word_encodes[seq_j],word_encodes_no[seq_j],lens[seq_j],jac);
#else
        	CountWeightedJaccard_SoA_cArray(word_encodes[seq_i],word_encodes_no[seq_i],lens[seq_i],
                word_encodes[seq_j],word_encodes_no[seq_j],lens[seq_j],jac);
#endif
            // CountWeightedJaccard_SoA_cArray(word_encodes[seq_i],word_encodes_no[seq_i],lens[seq_i],
                // word_encodes[seq_j],word_encodes_no[seq_j],lens[seq_j],jac);
            // cout<<jac<<endl;
            if(jac>=tau){
                dsu.unite(seq_i,seq_j);
            }
        }
    }
    
    for (int i = 0; i < N; ++i) {
        parent[seqs[i].seq_id] = seqs[dsu.find(i)].seq_id;
    }

    // free memory
    for(int i=0;i<N;i++){
        delete[] word_encodes[i];
        delete[] word_encodes_no[i];
    }
    delete[] word_encodes;
    delete[] word_encodes_no;
    delete[] buff;
    // cerr<<"Sorting total time: "<<sort_time<<" s\n";
}