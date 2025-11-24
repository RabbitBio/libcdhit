#include <immintrin.h>
#include <sys/time.h>
#include <zlib.h>

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

#include "CLI11.hpp"
#include "kseq.h"

using Clock = std::chrono::steady_clock;

#define NAA1 21
#define MAX_SEQ 655360
#define MAX_AA 23
#define NAA2 (NAA1 * NAA1) // 二元对的总数
#define MAX_DIAG (MAX_SEQ << 1)
#define FAILED_FUNC 1
#define OK_FUNC 0
using namespace std;
template <class TYPE>
class Vector : public vector<TYPE> {
public:
    Vector() : vector<TYPE>() {}
    Vector(size_t size) : vector<TYPE>(size) {}
    Vector(size_t size, const TYPE& deft) : vector<TYPE>(size, deft) {}

    void Append(const TYPE& item) {
        size_t n = this->size();
        if (n + 1 >= this->capacity()) this->reserve(n + n / 5 + 1);
        this->push_back(item);
    }
    int size() const { return (int)vector<TYPE>::size(); }
};

template <class TYPE>
class NVector {
public:
    TYPE* items;
    int size;
    int capacity;

    NVector() {
        size = capacity = 0;
        items = NULL;
    }
    NVector(int n, const TYPE& v = TYPE()) {
        size = capacity = 0;
        items = NULL;
        Resize(n, v);
    }
    NVector(const NVector& other) {
        size = capacity = 0;
        items = NULL;
        if (other.items) {
            Resize(other.size);
            memcpy(items, other.items, other.size * sizeof(TYPE));
        }
    }

    ~NVector() {
        if (items) free(items);
    }

    int Size() const { return size; }
    void Clear() {
        if (items) free(items);
        size = capacity = 0;
        items = NULL;
    }

    void Resize(int n, const TYPE& value = TYPE()) {
        if (n == size && capacity > 0) return;
        int i;
        // When resize() is called, probably this is the intended size,
        // and will not be changed frequently.
        if (n != capacity) {
            capacity = n;
            items = (TYPE*)realloc(items, capacity * sizeof(TYPE));
        }
        for (i = size; i < n; i++) items[i] = value;
        size = n;
    }
    void Append(const TYPE& item) {
        if (size + 1 >= capacity) {
            capacity = size + size / 5 + 1;
            items = (TYPE*)realloc(items, capacity * sizeof(TYPE));
        }
        items[size] = item;
        size++;
    }

    TYPE& operator[](const int i) {
        // if( i <0 or i >= size ) printf( "out of range\n" );
        return items[i];
    }
    TYPE& operator[](const int i) const {
        // if( i <0 or i >= size ) printf( "out of range\n" );
        return items[i];
    }
};
typedef NVector<int64_t> VectorInt64;
typedef Vector<VectorInt64> MatrixInt64;
typedef NVector<int> VectorInt;
typedef Vector<VectorInt> MatrixInt;
enum { DP_BACK_NONE = 0, DP_BACK_LEFT_TOP = 1, DP_BACK_LEFT = 2, DP_BACK_TOP = 3 };

KSEQ_INIT(gzFile, gzread)
int aa2idx[] = {0, 2, 4, 3, 6, 13, 7, 8, 9, 20, 11, 10, 12, 2, 20, 14, 5, 1, 15, 16, 20, 19, 17, 20, 18, 6};

int BLOSUM62[] = {
    4,                                                                                     // A
    -1, 5,                                                                                 // R
    -2, 0,  6,                                                                             // N
    -2, -2, 1,  6,                                                                         // D
    0,  -3, -3, -3, 9,                                                                     // C
    -1, 1,  0,  0,  -3, 5,                                                                 // Q
    -1, 0,  0,  2,  -4, 2,  5,                                                             // E
    0,  -2, 0,  -1, -3, -2, -2, 6,                                                         // G
    -2, 0,  1,  -1, -3, 0,  0,  -2, 8,                                                     // H
    -1, -3, -3, -3, -1, -3, -3, -4, -3, 4,                                                 // I
    -1, -2, -3, -4, -1, -2, -3, -4, -3, 2,  4,                                             // L
    -1, 2,  0,  -1, -3, 1,  1,  -2, -1, -3, -2, 5,                                         // K
    -1, -1, -2, -3, -1, 0,  -2, -3, -2, 1,  2,  -1, 5,                                     // M
    -2, -3, -3, -3, -2, -3, -3, -3, -1, 0,  0,  -3, 0,  6,                                 // F
    -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7,                             // P
    1,  -1, 1,  0,  -1, 0,  0,  0,  -1, -2, -2, 0,  -1, -2, -1, 4,                         // S
    0,  -1, 0,  -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1,  5,                     // T
    -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1,  -4, -3, -2, 11,                // W
    -2, -2, -2, -3, -2, -1, -2, -3, 2,  -1, -1, -2, -1, 3,  -3, -2, -2, 2,  7,             // Y
    0,  -3, -3, -3, -1, -2, -2, -3, -3, 3,  1,  -2, 1,  -1, -2, -2, 0,  -3, -1, 4,         // V
    -2, -1, 3,  4,  -3, 0,  1,  -1, 0,  -3, -4, 0,  -3, -3, -2, 0,  -1, -4, -3, -3, 4,     // B
    -1, 0,  0,  1,  -3, 3,  4,  -2, 0,  -3, -3, 1,  -1, -3, -1, 0,  -1, -3, -2, -2, 1,  4, // Z
    0,  -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, 0,  0,  -2, -1, -1, -1, -1,
    -1 // X
    // A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
    // 0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19  2  6 20
};
int BLOSUM62_na[] = {
    2,                     // A
    -2, 2,                 // C
    -2, -2, 2,             // G
    -2, -2, -2, 2,         // T
    -2, -2, -2, 1,  2,     // U
    -2, -2, -2, -2, -2, 1, // N
    0,  0,  0,  0,  0,  0,
    1 // X
    // A  C  G  T  U  N  X
    // 0  1  2  3  3  4  5
};
double get_time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
}
struct Sequence {
    string name;
    string comment;
    string seq;
    int len;
};
// WorkingBuffer结构体定义
struct WorkingBuffer {
    vector<int> taap;        // 用于存储二元对出现次数
    vector<int> aap_begin;   // 存储每个二元对的起始位置
    vector<int> aap_list;    // 存储实际的二元对出现位置
    vector<int> diag_score;  // 存储对角线命中次数
    vector<int> diag_score2; // 存储加权命中次数
    MatrixInt64 score_mat;
    MatrixInt back_mat;

    alignas(64) int64_t avx_i_arr[8];
    alignas(64) int64_t avx_j_arr[8];
    alignas(64) int64_t avx_m_arr[8];
    alignas(64) int64_t avx_sij_arr[8];

    alignas(64) int64_t* flat_matrix;

    // 构造函数，初始化大小
    WorkingBuffer(int size) {
        taap.resize(NAA2, 0);
        aap_begin.resize(NAA2, 0);
        aap_list.resize(size);
        diag_score.resize(MAX_DIAG, 0);
        diag_score2.resize(MAX_DIAG, 0);
    }
};
////////// Class definition //////////
class ScoreMatrix { // Matrix
private:
public:
    int matrix[MAX_AA][MAX_AA];
    alignas(64) int64_t* flat_matrix;
    int gap, ext_gap;

    ScoreMatrix();
    void init();
    void set_gap(int gap1, int ext_gap1);
    void set_matrix(int* mat1);
    void set_to_na();
    void set_match(int score);
    void set_mismatch(int score);
    void update_flat_matrix();
}; // END class ScoreMatrix

ScoreMatrix mat;
/////////////////
ScoreMatrix::ScoreMatrix() {
    flat_matrix = (int64_t*)aligned_alloc(64, MAX_AA * MAX_AA * sizeof(int64_t));
    if (!flat_matrix) {
        fprintf(stderr, "Failed to allocate flat_matrix\n");
        exit(1);
    }
    init();
}

void ScoreMatrix::init() {
    set_gap(-11, -1);
    set_matrix(BLOSUM62);
}

void ScoreMatrix::set_gap(int gap1, int ext_gap1) {
    int i;
    gap = MAX_SEQ * gap1;
    ext_gap = MAX_SEQ * ext_gap1;
}

void ScoreMatrix::set_matrix(int* mat1) {
    int i, j, k;
    k = 0;
    for (i = 0; i < MAX_AA; i++)
        for (j = 0; j <= i; j++) matrix[j][i] = matrix[i][j] = MAX_SEQ * mat1[k++];
    update_flat_matrix();
}

void ScoreMatrix::update_flat_matrix() {
    for (int i = 0; i < MAX_AA; i++) {
        for (int j = 0; j < MAX_AA; j++) {
            flat_matrix[i * MAX_AA + j] = (int64_t)matrix[i][j];
        }
    }
}

void ScoreMatrix::set_to_na() {
    set_gap(-6, -1);
    set_matrix(BLOSUM62_na);
}
// Only for est
void ScoreMatrix::set_match(int score) {
    int i;
    for (i = 0; i < 5; i++) matrix[i][i] = MAX_SEQ * score;
    // matrix[3][4] = matrix[4][3] = MAX_SEQ * score;
}
// Only for est
void ScoreMatrix::set_mismatch(int score) {
    int i, j;
    for (i = 0; i < MAX_AA; i++)
        for (j = 0; j < i; j++) matrix[j][i] = matrix[i][j] = MAX_SEQ * score;
    matrix[3][4] = matrix[4][3] = MAX_SEQ;
}

// 计算二元对的倒排表
void ComputeAAP(const char* seqi, int size, WorkingBuffer& buffer) {
    int len1 = size - 1;
    int sk, j1, mm, c22;

    // 统计每个二元对出现的次数
    for (sk = 0; sk < NAA2; sk++) buffer.taap[sk] = 0;
    for (j1 = 0; j1 < len1; j1++) {
        c22 = seqi[j1] * NAA1 + seqi[j1 + 1];
        buffer.taap[c22]++;
    }

    // 计算aap_begin数组
    for (sk = 0, mm = 0; sk < NAA2; sk++) {
        buffer.aap_begin[sk] = mm;
        mm += buffer.taap[sk];
        buffer.taap[sk] = 0;
    }

    // 填充aap_list数组
    for (j1 = 0; j1 < len1; j1++) {
        c22 = seqi[j1] * NAA1 + seqi[j1 + 1];
        buffer.aap_list[buffer.aap_begin[c22] + buffer.taap[c22]++] = j1;
    }
}

// 计算最佳带宽和最佳命中数
int diag_test_aapn(const char* iseq2, int len1, int len2, WorkingBuffer& buffer, int& best_sum, int band_width,
                   int& band_left, int& band_center, int& band_right, int required_aa1) {
    int i, i1, j, k;
    int nall = len1 + len2 - 1; // 总对角线数
    vector<int>& taap = buffer.taap;
    vector<int>& aap_begin = buffer.aap_begin;
    vector<int>& aap_list = buffer.aap_list;
    vector<int>& diag_score = buffer.diag_score;
    vector<int>& diag_score2 = buffer.diag_score2;

    if (nall > MAX_DIAG) {
        cerr << "in diag_test_aapn, MAX_DIAG reached" << endl;
        return -1;
    }

    // 初始化diag_score和diag_score2
    fill(diag_score.begin(), diag_score.end(), 0);
    fill(diag_score2.begin(), diag_score2.end(), 0);

    int c22, cpx;
    int len11 = len1 - 1;
    int len22 = len2 - 1;
    i1 = len11;

    // 填充diag_score和diag_score2
    for (i = 0; i < len22; i++, i1++) {
        c22 = iseq2[i] * NAA1 + iseq2[i + 1];
        cpx = 1 + (iseq2[i] != iseq2[i + 1]);
        if ((j = taap[c22]) == 0) continue;
        int m = aap_begin[c22];
        for (k = 0; k < j; k++) {
            diag_score[i1 - aap_list[m + k]]++;
            diag_score2[i1 - aap_list[m + k]] += cpx;
        }
    }

    // 找到最佳带宽范围
    int band_b = required_aa1 - 1 >= 0 ? required_aa1 - 1 : 0;
    int band_e = nall - band_b;

    int band_m = (band_b + band_width - 1 < band_e) ? band_b + band_width - 1 : band_e;
    int best_score = 0;
    int best_score2 = 0;
    int max_diag = 0;
    int max_diag2 = 0;
    int imax_diag = 0;

    for (i = band_b; i <= band_m; i++) {
        best_score += diag_score[i];
        best_score2 += diag_score2[i];
        if (diag_score2[i] > max_diag2) {
            max_diag2 = diag_score2[i];
            max_diag = diag_score[i];
            imax_diag = i;
        }
    }

    int from = band_b;
    int end = band_m;
    int score = best_score;
    int score2 = best_score2;
    for (k = from, j = band_m + 1; j < band_e; j++, k++) {
        score -= diag_score[k];
        score += diag_score[j];
        score2 -= diag_score2[k];
        score2 += diag_score2[j];
        if (score2 > best_score2) {
            from = k + 1;
            end = j;
            best_score = score;
            best_score2 = score2;
            if (diag_score2[j] > max_diag2) {
                max_diag2 = diag_score2[j];
                max_diag = diag_score[j];
                imax_diag = j;
            }
        }
    }

    int mlen = imax_diag;
    if (imax_diag > len1) mlen = nall - imax_diag;
    int emax = int((1.0 - 0.9) * mlen) + 1; // 假设 options.cluster_thd = 0.5
    for (j = from; j < imax_diag; j++) {
        if ((imax_diag - j) > emax || diag_score[j] < 1) {
            best_score -= diag_score[j];
            from++;
        } else
            break;
    }
    for (j = end; j > imax_diag; j--) {
        if ((j - imax_diag) > emax || diag_score[j] < 1) {
            best_score -= diag_score[j];
            end--;
        } else
            break;
    }

    // 计算最终带宽
    band_left = from - len1 + 1;
    band_right = end - len1 + 1;
    band_center = imax_diag - len1 + 1;
    best_sum = best_score;

    return 0;
}

int diag_no_table(const char* seqi, const size_t sizei, const char* seqj, const size_t sizej, WorkingBuffer& buffer,
                  int& best_sum, int band_width, int& band_left, int& band_center, int& band_right, int required_aa1) {
    if (sizei < 2 || sizej < 2) {
        buffer.diag_score.clear();
        buffer.diag_score2.clear();
        return 0;
    }
    const int len1 = static_cast<int>(sizei) - 1; // 序列1二元对个数
    const int len2 = static_cast<int>(sizej) - 1; // 序列2二元对个数
    const int nall = len1 + len2 - 1;
    if (nall <= 0 || nall > MAX_DIAG) {
        buffer.diag_score.clear();
        buffer.diag_score2.clear();
        return 0;
    }

    // --- 准备输出打分区 ---
    buffer.diag_score.assign(nall, 0);
    buffer.diag_score2.assign(nall, 0);
    auto& diag_score = buffer.diag_score;
    auto& diag_score2 = buffer.diag_score2;

    // --- 构造 (code, pos) 列表，同时统计频次（用于前缀和） ---
    std::vector<std::pair<int, int>> seqi_list;
    seqi_list.reserve(len1);
    std::vector<std::pair<int, int>> seqj_list;
    seqj_list.reserve(len2);

    // 注意：把 char 转成 unsigned，避免负值
    std::vector<int> cnt_i(NAA2, 0), cnt_j(NAA2, 0);

    for (int i = 0; i < len1; ++i) {
        int a = static_cast<unsigned char>(seqi[i]);
        int b = static_cast<unsigned char>(seqi[i + 1]);
        int code = a * NAA1 + b; // 0..NAA2-1
        seqi_list.emplace_back(code, i);
        ++cnt_i[code];
    }
    for (int j = 0; j < len2; ++j) {
        int a = static_cast<unsigned char>(seqj[j]);
        int b = static_cast<unsigned char>(seqj[j + 1]);
        int code = a * NAA1 + b;
        seqj_list.emplace_back(code, j);
        ++cnt_j[code];
    }

    std::sort(seqi_list.begin(), seqi_list.end());
    std::sort(seqj_list.begin(), seqj_list.end());

    // --- 建立“前缀和跳转表”：start[c] = 代码块 c 的起始下标；并设哨兵 start[NAA2] = size ---
    std::vector<int> start_i(NAA2 + 1, 0), start_j(NAA2 + 1, 0);
    for (int c = 0; c < NAA2; ++c) start_i[c + 1] = start_i[c] + cnt_i[c];
    for (int c = 0; c < NAA2; ++c) start_j[c + 1] = start_j[c] + cnt_j[c];
    // 现在：
    //   seqi 中 code=c 的块区间是 [ start_i[c], start_i[c+1] )
    //   seqj 中 code=c 的块区间是 [ start_j[c], start_j[c+1] )

    // --- 双指针 + 快速跳转 ---
    size_t i = 0, j = 0;
    while (i < seqi_list.size() && j < seqj_list.size()) {
        int code_i = seqi_list[i].first;
        int code_j = seqj_list[j].first;

        if (code_i < code_j) {
            // 一步跳到 seqi 中 code >= code_j 的第一个位置
            size_t ni = static_cast<size_t>(start_i[code_j]);
            i = (ni > i ? ni : i); // 防御：保证单调前进
            if (i >= seqi_list.size()) break;
            continue;
        }
        if (code_i > code_j) {
            // 一步跳到 seqj 中 code >= code_i 的第一个位置
            size_t nj = static_cast<size_t>(start_j[code_i]);
            j = (nj > j ? nj : j);
            if (j >= seqj_list.size()) break;
            continue;
        }

        // code_i == code_j：处理这一个 code 的整块
        const int code = code_i;
        size_t i_end = static_cast<size_t>(start_i[code + 1]);
        size_t j_end = static_cast<size_t>(start_j[code + 1]);

        // 对两个块做笛卡尔积，把票投到对角线
        for (size_t jj = j; jj < j_end; ++jj) {
            const int jpos = seqj_list[jj].second;              // seqj 的二元对起点
            const int i1 = (len1 - 1) + jpos;                   // 对角线坐标偏移
            const int cpx = 1 + (seqj[jpos] != seqj[jpos + 1]); // 加权：相同=1，不同=2
            for (size_t ii = i; ii < i_end; ++ii) {
                const int ipos = seqi_list[ii].second; // seqi 的二元对起点
                const int d = i1 - ipos;               // 对角线下标 0..nall-1
                if ((unsigned)d < (unsigned)nall) {
                    diag_score[d] += 1;
                    diag_score2[d] += cpx;
                }
            }
        }

        // 整块处理完，一次性跳到块末尾
        i = i_end;
        j = j_end;
    }
    int band_b = required_aa1 - 1 >= 0 ? required_aa1 - 1 : 0;
    int band_e = nall - band_b;

    int band_m = (band_b + band_width - 1 < band_e) ? band_b + band_width - 1 : band_e;
    int best_score = 0;
    int best_score2 = 0;
    int max_diag = 0;
    int max_diag2 = 0;
    int imax_diag = 0;

    for (i = band_b; i <= band_m; i++) {
        best_score += diag_score[i];
        best_score2 += diag_score2[i];
        if (diag_score2[i] > max_diag2) {
            max_diag2 = diag_score2[i];
            max_diag = diag_score[i];
            imax_diag = i;
        }
    }

    int from = band_b;
    int end = band_m;
    int score = best_score;
    int score2 = best_score2;
    for (int k = from, j = band_m + 1; j < band_e; j++, k++) {
        score -= diag_score[k];
        score += diag_score[j];
        score2 -= diag_score2[k];
        score2 += diag_score2[j];
        if (score2 > best_score2) {
            from = k + 1;
            end = j;
            best_score = score;
            best_score2 = score2;
            if (diag_score2[j] > max_diag2) {
                max_diag2 = diag_score2[j];
                max_diag = diag_score[j];
                imax_diag = j;
            }
        }
    }

    int mlen = imax_diag;
    if (imax_diag > len1) mlen = nall - imax_diag;
    int emax = int((1.0 - 0.9) * mlen) + 1; // 假设 options.cluster_thd = 0.5
    for (j = from; j < imax_diag; j++) {
        if ((imax_diag - j) > emax || diag_score[j] < 1) {
            best_score -= diag_score[j];
            from++;
        } else
            break;
    }
    for (j = end; j > imax_diag; j--) {
        if ((j - imax_diag) > emax || diag_score[j] < 1) {
            best_score -= diag_score[j];
            end--;
        } else
            break;
    }

    // 计算最终带宽
    band_left = from - len1 + 1;
    band_right = end - len1 + 1;
    band_center = imax_diag - len1 + 1;
    best_sum = best_score;

    return 0;
}

// #define __AVX2__
#ifdef __AVX2__

int Init_Matrix_AVX2(char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix& mat, int& best_score, int& iden_no,
                     int& alnln, float& dist, int* alninfo, int band_left, int band_center, int band_right,
                     WorkingBuffer& buffer) {
    if ((band_right >= len2) || (band_left <= -len1) || (band_left > band_right)) return FAILED_FUNC;

    int band_width = band_right - band_left + 1;
    int band_width1 = 13;

    MatrixInt64& score_mat = buffer.score_mat;
    MatrixInt& back_mat = buffer.back_mat;

    int L = band_left, R = band_right;
    int kmin = (R < 0) ? -R : (L > 0) ? L : 0;
    int kmax = (R < len2 - len1) ? (R + 2 * len1) : (L > len2 - len1) ? (2 * len2 - L) : (len1 + len2);

    int lenY = kmax - kmin + 1;
    if (score_mat.size() <= lenY) {
        VectorInt row(band_width1, 0);
        VectorInt64 row2(band_width1, 0);
        while (score_mat.size() <= lenY) {
            score_mat.Append(row2);
            back_mat.Append(row);
        }
    }
    for (int i = 0; i <= lenY; i++) {
        if (score_mat[i].size < band_width1) score_mat[i].Resize(band_width1);
        if (back_mat[i].size < band_width1) back_mat[i].Resize(band_width1);
    }
    return OK_FUNC;
}

inline __m256i div2_s64(__m256i x) {
    __m256i sign = _mm256_srli_epi64(x, 63);
    __m256i mask = _mm256_sub_epi64(_mm256_setzero_si256(), sign);
    __m256i shifted = _mm256_srli_epi64(x, 1);
    __m256i result = _mm256_or_si256(shifted, _mm256_slli_epi64(mask, 63));
    return result;
}

inline __m128i avx2_cvtepi64_epi32(__m256i a) {
    __m256i shuffled = _mm256_shuffle_epi32(a, 0xD8); // 0b11011000 = 3,1,2,0
    __m128i low = _mm256_castsi256_si128(shuffled);
    __m128i high = _mm256_extracti128_si256(shuffled, 1);
    return _mm_unpacklo_epi64(low, high);
}

int rotation_band_align_AVX2(char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix& mat, int& best_score,
                             int& iden_no, int& alnln, float& dist, int* alninfo, int band_left, int band_center,
                             int band_right, WorkingBuffer& buffer) {
    auto t1 = Clock::now();

    int i, j, k, j1;
    int jj, kk;
    int x, y;
    int iden_no1;
    int64_t best_score1;
    iden_no = 0;

    if ((band_right >= len2) || (band_left <= -len1) || (band_left > band_right)) return FAILED_FUNC;

    int band_width = band_right - band_left + 1;
    int band_width1 = 13;

    MatrixInt64& score_mat = buffer.score_mat;
    MatrixInt& back_mat = buffer.back_mat;

    int L = band_left, R = band_right;
    int kmin = (R < 0) ? -R : (L > 0) ? L : 0;
    int kmax = (R < len2 - len1) ? (R + 2 * len1) : (L > len2 - len1) ? (2 * len2 - L) : (len1 + len2);

    int lenY = kmax - kmin + 1;

    if (score_mat.size() <= lenY) {
        VectorInt row(band_width1, 0);
        VectorInt64 row2(band_width1, 0);
        while (score_mat.size() <= lenY) {
            score_mat.Append(row2);
            back_mat.Append(row);
        }
    }
    for (int i = 0; i <= lenY; i++) {
        if (score_mat[i].size < band_width1) score_mat[i].Resize(band_width1);
        if (back_mat[i].size < band_width1) back_mat[i].Resize(band_width1);
    }

    best_score = 0;

    auto t5 = Clock::now();

    // 初始化边界
    if (L < 0) {
        int T = (R < 0) ? R : 0;
        for (int X = L; X <= T; ++X) {
            int I = -X, J = 0;
            if (I < 0 || I > len1) continue;
            int n = (I + J) - kmin;
            int m = (X - L + 1) >> 1;
            score_mat[n][m] = (int64_t)mat.ext_gap * I;
            back_mat[n][m] = DP_BACK_TOP;
        }
        back_mat[(-T) - kmin][(T - L + 1) / 2] = DP_BACK_NONE;
    }

    if (R >= 0) {
        int T = (L > 0) ? L : 0;
        for (int X = T; X <= R; ++X) {
            int I = 0, J = X;
            if (J < 0 || J > len2) continue;
            int n = (I + J) - kmin;
            int m = (X - L + 1) >> 1;
            score_mat[n][m] = (int64_t)mat.ext_gap * J;
            back_mat[n][m] = DP_BACK_LEFT;
        }
        back_mat[T - kmin][(T - L + 1) / 2] = DP_BACK_NONE;
    }

    int gap_open[2] = {mat.gap, mat.ext_gap};
    int max_diag = band_center - band_left;
    int extra_score[4] = {4, 3, 2, 1};
    auto t2 = Clock::now();

    // ============ AVX-2 向量化变量定义 ============
    const int SIMD_WIDTH = 4;

    int64_t* i_arr = buffer.avx_i_arr;
    int64_t* j_arr = buffer.avx_j_arr;
    int64_t* m_arr = buffer.avx_m_arr;
    int64_t* sij_arr = buffer.avx_sij_arr;

    __m256i vec_0 = _mm256_setzero_si256();
    __m256i vec_1 = _mm256_set1_epi64x(1);
    __m256i vec_3 = _mm256_set1_epi64x(3);
    __m256i vec_4 = _mm256_set1_epi64x(4);

    __m256i vec_len1 = _mm256_set1_epi64x(len1);
    __m256i vec_len2 = _mm256_set1_epi64x(len2);
    __m256i vec_lenY = _mm256_set1_epi64x(lenY);
    __m256i vec_R = _mm256_set1_epi64x(R);

    __m256i vec_band_center = _mm256_set1_epi64x(band_center);

    __m256i vec_ext_gap = _mm256_set1_epi64x(mat.ext_gap);
    __m256i vec_gap_open = _mm256_set1_epi64x(mat.gap);

    __m128i vec_DP_BACK_LEFT = _mm_set1_epi32(DP_BACK_LEFT);
    __m128i vec_DP_BACK_TOP = _mm_set1_epi32(DP_BACK_TOP);
    __m128i vec_DP_BACK_LEFT_TOP = _mm_set1_epi32(DP_BACK_LEFT_TOP);
    __m256i vec_offset = _mm256_setr_epi64x(0, 2, 4, 6);

    // ============ 主循环 ============
    for (int y = kmin + 1; y <= kmax; y++) {
        int offset = (abs(y + L)) % 2;
        int x_start = L + offset;

        __m256i vec_y = _mm256_set1_epi64x(y);
        __m256i vec_y_minus_kmin = _mm256_set1_epi64x(y - kmin);

        int num_elements = 12;
        int vec_iterations = num_elements / SIMD_WIDTH;

        int index_y = y - kmin;
        int index_y1 = y - 1 - kmin;
        int index_y2 = y - 2 - kmin;

        // 预取下一行数据
        if (y + 1 <= kmax) {
            _mm_prefetch((const char*)&score_mat[y + 1 - kmin][0], _MM_HINT_T0);
            _mm_prefetch((const char*)&back_mat[y + 1 - kmin][0], _MM_HINT_T0);
        }

        // ============ 向量化主循环 ============
        for (int vec_idx = 0; vec_idx < vec_iterations; vec_idx++) {
            int x_base = x_start + vec_idx * SIMD_WIDTH * 2;
            int index_x = (x_base - L + 1) >> 1;
            int index_x_L = (x_base - 1 - L + 1) >> 1;
            int index_x_R = (x_base + 1 - L + 1) >> 1;
            // cout<<"x_base:"<<x_base<<endl;

            __m256i vec_x = _mm256_add_epi64(_mm256_set1_epi64x(x_base), vec_offset);

            __m256i vec_i = _mm256_sub_epi64(vec_y, vec_x);
            vec_i = div2_s64(vec_i);
            __m256i vec_j = _mm256_add_epi64(vec_x, vec_y);
            vec_j = div2_s64(vec_j);

            __m256i vec_i_gt_0 = _mm256_cmpgt_epi64(vec_i, vec_0);
            __m256i vec_i_le_len1 =
                _mm256_or_si256(_mm256_cmpgt_epi64(vec_len1, vec_i), _mm256_cmpeq_epi64(vec_len1, vec_i));

            __m256i vec_j_gt_0 = _mm256_cmpgt_epi64(vec_j, vec_0);
            __m256i vec_j_le_len2 =
                _mm256_or_si256(_mm256_cmpgt_epi64(vec_len2, vec_j), _mm256_cmpeq_epi64(vec_len2, vec_j));

            __m256i vec_mask = _mm256_and_si256(_mm256_and_si256(vec_i_gt_0, vec_i_le_len1),
                                                _mm256_and_si256(vec_j_gt_0, vec_j_le_len2));

            int mask = _mm256_movemask_pd(_mm256_castsi256_pd(vec_mask));

            int count = _mm_popcnt_u32((unsigned int)mask);
            if (count == 0) continue;

            // printf(" (0x%02X)\n", (unsigned int)mask);
            __m256i vec_i_minus_1 = _mm256_sub_epi64(vec_i, vec_1);

            __m256i vec_j_minus_1 = _mm256_sub_epi64(vec_j, vec_1);

            _mm256_store_si256((__m256i*)i_arr, vec_i_minus_1);
            _mm256_store_si256((__m256i*)j_arr, vec_j_minus_1);

            // 获取序列值和矩阵分数
            for (int idx = 0; idx < SIMD_WIDTH; idx++) {
                int bit = (mask >> idx) & 1;
                if (!bit) {
                    sij_arr[idx] = 0;
                } else {
                    sij_arr[idx] = mat.matrix[iseq1[i_arr[idx]]][iseq2[j_arr[idx]]];
                }
            }
            __m256i vec_sij = _mm256_load_si256((__m256i*)sij_arr);

            // int extra = extra_score[abs(x-band_center)&3];
            __m256i vec_extra = _mm256_sub_epi64(vec_x, vec_band_center);

            // vec_extra = _mm256_abs_epi64(vec_extra);
            __m256i vec_extra_sign = _mm256_cmpgt_epi64(vec_0, vec_extra);
            __m256i vec_extra_neg = _mm256_sub_epi64(vec_0, vec_extra);
            vec_extra = _mm256_blendv_epi8(vec_extra, vec_extra_neg, vec_extra_sign);

            vec_extra = _mm256_and_si256(vec_extra, vec_3);
            vec_extra = _mm256_sub_epi64(vec_4, vec_extra);

            __m256i vec_sij_gt_0 = _mm256_cmpgt_epi64(vec_sij, vec_0);

            vec_extra = _mm256_and_si256(vec_extra, vec_sij_gt_0);
            vec_sij = _mm256_add_epi64(vec_sij, vec_extra);

            __m256i vec_best_score1 = _mm256_loadu_si256((__m256i*)&score_mat[index_y][index_x]);
            __m128i vec_back = _mm_loadu_si128((__m128i*)&back_mat[index_y][index_x]);

            // left-top
            if (y - 2 >= kmin) {
                __m256i vec_score_y2 = _mm256_loadu_si256((__m256i*)&score_mat[index_y2][index_x]);
                vec_best_score1 = _mm256_add_epi64(vec_score_y2, vec_sij);
                vec_back = _mm_set1_epi32(DP_BACK_LEFT_TOP);
            }

            __m256i vec_mask_gap =
                _mm256_or_si256(_mm256_cmpeq_epi64(vec_i, vec_len1), _mm256_cmpeq_epi64(vec_j, vec_len2));

            __m256i vec_gap0 = _mm256_blendv_epi8(vec_gap_open, vec_ext_gap, vec_mask_gap);
            __m256i vec_gap;
            __m256i vec_score;

            // left
            if (y - 1 >= kmin) {
                vec_gap = vec_gap0;

                __m256i vec_score_y1 = _mm256_loadu_si256((__m256i*)&score_mat[index_y1][index_x_L]);
                __m128i vec_back_y1 = _mm_loadu_si128((__m128i*)&back_mat[index_y1][index_x_L]);

                __m128i gap_flag = _mm_cmpeq_epi32(vec_back_y1, vec_DP_BACK_LEFT);
                __m256i gap_flag_64 = _mm256_cvtepi32_epi64(gap_flag);

                vec_gap = _mm256_blendv_epi8(vec_gap0, vec_ext_gap, gap_flag_64);
                vec_score = _mm256_add_epi64(vec_score_y1, vec_gap);

                __m256i modify_flag =
                    _mm256_and_si256(_mm256_cmpgt_epi64(vec_score, vec_best_score1),
                                     _mm256_set_epi64x(-1LL, -1LL, -1LL, ((vec_idx) | (offset)) ? -1LL : 0));
                vec_best_score1 = _mm256_blendv_epi8(vec_best_score1, vec_score, modify_flag);

                // __m128i modify_flag_32 = _mm256_cvtepi64_epi32(modify_flag);
                __m128i modify_flag_32 = avx2_cvtepi64_epi32(modify_flag);
                vec_back = _mm_blendv_epi8(vec_back, vec_DP_BACK_LEFT, modify_flag_32);
            }

            // top
            if (y - 1 >= kmin) {
                vec_gap = vec_gap0;

                __m256i vec_score_y1 = _mm256_loadu_si256((__m256i*)&score_mat[index_y1][index_x_R]);
                __m128i vec_back_y1 = _mm_loadu_si128((__m128i*)&back_mat[index_y1][index_x_R]);

                __m128i gap_flag = _mm_cmpeq_epi32(vec_back_y1, vec_DP_BACK_TOP);
                __m256i gap_flag_64 = _mm256_cvtepi32_epi64(gap_flag);

                vec_gap = _mm256_blendv_epi8(vec_gap0, vec_ext_gap, gap_flag_64);
                vec_score = _mm256_add_epi64(vec_score_y1, vec_gap);

                __m256i score_batter = _mm256_cmpgt_epi64(vec_score, vec_best_score1);
                __m256i vec_x_plus_1 = _mm256_add_epi64(vec_x, vec_1);
                __m256i x_in_range =
                    _mm256_or_si256(_mm256_cmpgt_epi64(vec_R, vec_x_plus_1), _mm256_cmpeq_epi64(vec_R, vec_x_plus_1));

                __m256i modify_flag = _mm256_and_si256(score_batter, x_in_range);
                vec_best_score1 = _mm256_blendv_epi8(vec_best_score1, vec_score, modify_flag);

                // __m128i modify_flag_32 = _mm256_cvtepi64_epi32(modify_flag);
                __m128i modify_flag_32 = avx2_cvtepi64_epi32(modify_flag);

                vec_back = _mm_blendv_epi8(vec_back, vec_DP_BACK_TOP, modify_flag_32);
            }

            __m256i old_score = _mm256_loadu_si256((__m256i*)&score_mat[index_y][index_x]);
            __m256i result_score = _mm256_blendv_epi8(old_score, vec_best_score1, vec_mask);
            _mm256_storeu_si256((__m256i*)&score_mat[index_y][index_x], result_score);

            __m128i old_back = _mm_loadu_si128((__m128i*)&back_mat[index_y][index_x]);
            // __m128i result_back = _mm_blendv_epi8(old_back, vec_back, _mm256_cvtepi64_epi32(vec_mask));
            __m128i result_back = _mm_blendv_epi8(old_back, vec_back, avx2_cvtepi64_epi32(vec_mask));
            _mm_storeu_si128((__m128i*)&back_mat[index_y][index_x], result_back);
        }
        // // if(y <= kmax - 20) continue;
        // for(int x=L+(abs(y+L))%2;x<=R;x+=2){
        // 	// avx<<score_mat[y-kmin][(x-L+1)>>1]<<" ";
        // 	cout << score_mat[y-kmin][(x-L+1)>>1]<<" ";
        // }
        // // avx<<endl;
        // cout << endl;
        // if(y > kmin+20) exit(0);
    }
    auto t3 = Clock::now();
    // cout<<t<<endl;

    // ============ 回溯部分（保持原样）============
    x = (R < len2 - len1) ? R : (L > len2 - len1) ? L : len2 - len1;
    y = kmax;
    i = (-x + y) >> 1;
    j = (x + y) >> 1;
    // printf("\n(%d,%d)\n", i, j);
    best_score = score_mat[y - kmin][(x - L + 1) >> 1];
    best_score1 = score_mat[y - kmin][(x - L + 1) >> 1];
    // printf("%2li(%2i) ", best_score1, iden_no1);

    int back = back_mat[y - kmin][(x - L + 1) >> 1];
    int last = back;

    int count = 0, count2 = 0, count3 = 0;
    int match, begin1, begin2, end1, end2;
    int gbegin1 = 0, gbegin2 = 0, gend1 = 0, gend2 = 0;
    int64_t score, smin = best_score1, smax = best_score1 - 1;
    int posmin, posmax, pos = 0;
    int bl, dlen = 0, dcount = 0;
    posmin = posmax = 0;
    begin1 = begin2 = end1 = end2 = 0;

    int masked = 0;
    int indels = 0;
    int max_indels = 0;

    while (back != DP_BACK_NONE) {
        switch (back) {
            case DP_BACK_TOP:
                bl = (last != back) & (j != 1) & (j != len2);
                dlen += bl;
                dcount += bl;
                score = score_mat[y - kmin][(x - L + 1) >> 1];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                i -= 1;
                x = x + 1;
                y = y - 1;
                break;
            case DP_BACK_LEFT:
                bl = (last != back) & (i != 1) & (i != len1);
                dlen += bl;
                dcount += bl;
                score = score_mat[y - kmin][(x - L + 1) >> 1];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                j -= 1;
                x = x - 1;
                y = y - 1;
                break;
            case DP_BACK_LEFT_TOP:
                if (alninfo && true) {
                    if (i == 1 || j == 1) {
                        gbegin1 = i - 1;
                        gbegin2 = j - 1;
                    } else if (i == len1 || j == len2) {
                        gend1 = i - 1;
                        gend2 = j - 1;
                    }
                }
                score = score_mat[y - kmin][(x - L + 1) >> 1];
                i -= 1;
                j -= 1;
                y -= 2;
                match = iseq1[i] == iseq2[j];
                if (score > smax) {
                    count = 0;
                    smax = score;
                    posmax = pos;
                    end1 = i;
                    end2 = j;
                }
                if (false && (iseq1[i] > 4 || iseq2[j] > 4)) {
                    masked += 1;
                } else {
                    dlen += 1;
                    dcount += !match;
                    count += match;
                    count2 += match;
                    count3 += match;
                }
                if (score < smin) {
                    int mm = match == 0;
                    count2 = 0;
                    smin = score;
                    posmin = pos - mm;
                    begin1 = i + mm;
                    begin2 = j + mm;
                }
                break;
            default:
                printf("%i\n", back);
                break;
        }
        pos += 1;
        last = back;
        back = back_mat[y - kmin][(x - L + 1) >> 1];
    }

    iden_no = true ? count3 : count - count2;
    alnln = posmin - posmax + 1 - masked;
    dist = dcount / (float)dlen;

    int umtail1 = len1 - 1 - end1;
    int umtail2 = len2 - 1 - end2;
    int umhead = begin1 < begin2 ? begin1 : begin2;
    int umtail = umtail1 < umtail2 ? umtail1 : umtail2;
    int umlen = umhead + umtail;

    auto t4 = Clock::now();

    auto d1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t5 - t1);
    auto d4 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t5);
    auto d2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2);
    auto d3 = std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3);

    cout << "Alloc Time: " << d1.count() << " ns\n";
    cout << "Init  Time: " << d4.count() << " ns\n";
    cout << "Vect  Time: " << d2.count() << " ns\n";
    cout << "Back  Time: " << d3.count() << " ns\n";

    if (umlen > 99999999) return FAILED_FUNC;
    if (umlen > len1 * 1.0) return FAILED_FUNC;
    if (umlen > len2 * 1.0) return FAILED_FUNC;

    if (alninfo) {
        alninfo[0] = begin1;
        alninfo[1] = end1;
        alninfo[2] = begin2;
        alninfo[3] = end2;
        alninfo[4] = masked;
        if (true) {
            alninfo[0] = gbegin1;
            alninfo[1] = gend1;
            alninfo[2] = gbegin2;
            alninfo[3] = gend2;
        }
    }
    return OK_FUNC;
}

#endif

// FIXME:
// tags: AVX512, Rotation(Anti-diag), Compact
// Finished 20/10/2025 mgl
#define __AVX512F__
#ifdef __AVX512F__

static inline void print_m512i_epi64(__m512i v, const char* label) {
    alignas(64) int64_t data[8];
    _mm512_store_si512(data, v);

    if (label) printf("%s: ", label);
    for (int i = 0; i < 8; i++) {
        printf("%lld ", (long long)data[i]);
    }
    printf("\n");
}

static inline void print_m256i_epi32(__m256i v, const char* label) {
    alignas(32) int32_t data[8];
    _mm256_store_si256((__m256i*)data, v);

    if (label) printf("%s: ", label);
    for (int i = 0; i < 8; i++) {
        printf("%lld ", data[i]);
    }
    printf("\n");
}

static inline void pirnt_m512i_epi64_mask(__m512i v, __mmask8 mask, const char* label) {
    alignas(64) int64_t data[8];
    _mm512_store_si512(data, v);

    if (label) printf("%s: ", label);
    for (int index = 0; index < 8; index++) {
        int bit = (mask >> index) & 1;
        if (!bit) continue;
        cout << data[index] << " ";
    }
    cout << endl;
}

int Init_Matrix_AVX512(char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix& mat, int& best_score, int& iden_no,
                       int& alnln, float& dist, int* alninfo, int band_left, int band_center, int band_right,
                       WorkingBuffer& buffer) {
    if ((band_right >= len2) || (band_left <= -len1) || (band_left > band_right)) return FAILED_FUNC;

    int band_width = band_right - band_left + 1;
    int band_width1 = 17;

    MatrixInt64& score_mat = buffer.score_mat;
    MatrixInt& back_mat = buffer.back_mat;

    int L = band_left, R = band_right;
    int kmin = (R < 0) ? -R : (L > 0) ? L : 0;
    int kmax = (R < len2 - len1) ? (R + 2 * len1) : (L > len2 - len1) ? (2 * len2 - L) : (len1 + len2);

    int lenY = kmax - kmin + 1;
    if (score_mat.size() <= lenY) {
        VectorInt row(band_width1, 0);
        VectorInt64 row2(band_width1, 0);
        while (score_mat.size() <= lenY) {
            score_mat.Append(row2);
            back_mat.Append(row);
        }
    }
    for (int i = 0; i <= lenY; i++) {
        if (score_mat[i].size < band_width1) score_mat[i].Resize(band_width1);
        if (back_mat[i].size < band_width1) back_mat[i].Resize(band_width1);
    }
    return OK_FUNC;
}

alignas(256) static const uint8_t REVERSE_MASK_LUT[256] = {
    0x00, 0x80, 0x40, 0xC0, 0x20, 0xA0, 0x60, 0xE0, 0x10, 0x90, 0x50, 0xD0, 0x30, 0xB0, 0x70, 0xF0, 0x08, 0x88, 0x48,
    0xC8, 0x28, 0xA8, 0x68, 0xE8, 0x18, 0x98, 0x58, 0xD8, 0x38, 0xB8, 0x78, 0xF8, 0x04, 0x84, 0x44, 0xC4, 0x24, 0xA4,
    0x64, 0xE4, 0x14, 0x94, 0x54, 0xD4, 0x34, 0xB4, 0x74, 0xF4, 0x0C, 0x8C, 0x4C, 0xCC, 0x2C, 0xAC, 0x6C, 0xEC, 0x1C,
    0x9C, 0x5C, 0xDC, 0x3C, 0xBC, 0x7C, 0xFC, 0x02, 0x82, 0x42, 0xC2, 0x22, 0xA2, 0x62, 0xE2, 0x12, 0x92, 0x52, 0xD2,
    0x32, 0xB2, 0x72, 0xF2, 0x0A, 0x8A, 0x4A, 0xCA, 0x2A, 0xAA, 0x6A, 0xEA, 0x1A, 0x9A, 0x5A, 0xDA, 0x3A, 0xBA, 0x7A,
    0xFA, 0x06, 0x86, 0x46, 0xC6, 0x26, 0xA6, 0x66, 0xE6, 0x16, 0x96, 0x56, 0xD6, 0x36, 0xB6, 0x76, 0xF6, 0x0E, 0x8E,
    0x4E, 0xCE, 0x2E, 0xAE, 0x6E, 0xEE, 0x1E, 0x9E, 0x5E, 0xDE, 0x3E, 0xBE, 0x7E, 0xFE, 0x01, 0x81, 0x41, 0xC1, 0x21,
    0xA1, 0x61, 0xE1, 0x11, 0x91, 0x51, 0xD1, 0x31, 0xB1, 0x71, 0xF1, 0x09, 0x89, 0x49, 0xC9, 0x29, 0xA9, 0x69, 0xE9,
    0x19, 0x99, 0x59, 0xD9, 0x39, 0xB9, 0x79, 0xF9, 0x05, 0x85, 0x45, 0xC5, 0x25, 0xA5, 0x65, 0xE5, 0x15, 0x95, 0x55,
    0xD5, 0x35, 0xB5, 0x75, 0xF5, 0x0D, 0x8D, 0x4D, 0xCD, 0x2D, 0xAD, 0x6D, 0xED, 0x1D, 0x9D, 0x5D, 0xDD, 0x3D, 0xBD,
    0x7D, 0xFD, 0x03, 0x83, 0x43, 0xC3, 0x23, 0xA3, 0x63, 0xE3, 0x13, 0x93, 0x53, 0xD3, 0x33, 0xB3, 0x73, 0xF3, 0x0B,
    0x8B, 0x4B, 0xCB, 0x2B, 0xAB, 0x6B, 0xEB, 0x1B, 0x9B, 0x5B, 0xDB, 0x3B, 0xBB, 0x7B, 0xFB, 0x07, 0x87, 0x47, 0xC7,
    0x27, 0xA7, 0x67, 0xE7, 0x17, 0x97, 0x57, 0xD7, 0x37, 0xB7, 0x77, 0xF7, 0x0F, 0x8F, 0x4F, 0xCF, 0x2F, 0xAF, 0x6F,
    0xEF, 0x1F, 0x9F, 0x5F, 0xDF, 0x3F, 0xBF, 0x7F, 0xFF};

alignas(32) static const int32_t LEFT_SHIFT_TABLE[8][8] = {
    {0, 1, 2, 3, 4, 5, 6, 7}, {7, 0, 1, 2, 3, 4, 5, 6}, {7, 7, 0, 1, 2, 3, 4, 5}, {7, 7, 7, 0, 1, 2, 3, 4},
    {7, 7, 7, 7, 0, 1, 2, 3}, {7, 7, 7, 7, 7, 0, 1, 2}, {7, 7, 7, 7, 7, 7, 0, 1}, {7, 7, 7, 7, 7, 7, 7, 0}};

alignas(32) static const int32_t RIGHT_SHIFT_TABLE[8][8] = {
    {7, 6, 5, 4, 3, 2, 1, 0}, {6, 5, 4, 3, 2, 1, 0, 7}, {5, 4, 3, 2, 1, 0, 7, 7}, {4, 3, 2, 1, 0, 7, 7, 7},
    {3, 2, 1, 0, 7, 7, 7, 7}, {2, 1, 0, 7, 7, 7, 7, 7}, {1, 0, 7, 7, 7, 7, 7, 7}, {0, 7, 7, 7, 7, 7, 7, 7}};
alignas(32) static const int32_t SHUFFLE_TABLE[8][8] = {
    {0, 1, 2, 3, 4, 5, 6, 7}, {1, 2, 3, 4, 5, 6, 7, 0}, {2, 3, 4, 5, 6, 7, 0, 1}, {3, 4, 5, 6, 7, 0, 1, 2},
    {4, 5, 6, 7, 0, 1, 2, 3}, {5, 6, 7, 0, 1, 2, 3, 4}, {6, 7, 0, 1, 2, 3, 4, 5}, {7, 0, 1, 2, 3, 4, 5, 6}};

void _mm256_load_char_array_forward(__m256i& vec_index, const char* arr, __m256i& seq_vals, __mmask8 mask) {
    __mmask8 copy_mask = mask;
    uint32_t lowbit = copy_mask & (-copy_mask);
    int pos = _mm_popcnt_u32(lowbit - 1);
    __m256i shuffle_perm = _mm256_load_si256((__m256i*)SHUFFLE_TABLE[pos]);
    vec_index = seq_vals = _mm256_permutevar8x32_epi32(vec_index, shuffle_perm);
    int base_index = _mm256_extract_epi32(vec_index, 0);
    copy_mask >>= pos;
    __mmask16 mask_16 = copy_mask;
    __m128i bytes = _mm_maskz_loadu_epi8(mask_16, arr + base_index);
    seq_vals = _mm256_cvtepu8_epi32(bytes);
    __m256i perm = _mm256_load_si256((__m256i*)LEFT_SHIFT_TABLE[pos]);
    seq_vals = _mm256_permutevar8x32_epi32(seq_vals, perm);
    seq_vals = _mm256_maskz_mov_epi32(mask, seq_vals);
}

void _mm256_load_char_array_backward(__m256i& vec_index, const char* arr, __m256i& seq_vals, __mmask8 mask) {
    __mmask8 reverse_mask = REVERSE_MASK_LUT[mask];
    __mmask8 copy_mask = reverse_mask;
    uint32_t lowbit = copy_mask & (-copy_mask);
    int pos = _mm_popcnt_u32(lowbit - 1);
    __m256i shuffle_perm = _mm256_load_si256((__m256i*)SHUFFLE_TABLE[7 - pos]);
    vec_index = seq_vals = _mm256_permutevar8x32_epi32(vec_index, shuffle_perm);
    int base_index = _mm256_extract_epi32(vec_index, 0);
    copy_mask >>= pos;
    __mmask16 mask_16 = copy_mask;
    __m128i bytes = _mm_maskz_loadu_epi8(mask_16, arr + base_index);
    seq_vals = _mm256_cvtepu8_epi32(bytes);
    __m256i perm = _mm256_load_si256((__m256i*)RIGHT_SHIFT_TABLE[pos]);
    seq_vals = _mm256_permutevar8x32_epi32(seq_vals, perm);
    seq_vals = _mm256_maskz_mov_epi32(mask, seq_vals);
}

int rotation_band_align_AVX512(char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix& mat, int& best_score,
                               int& iden_no, int& alnln, float& dist, int* alninfo, int band_left, int band_center,
                               int band_right, WorkingBuffer& buffer) {
    auto t1 = Clock::now();

    int i, j, k, j1;
    int jj, kk;
    int x, y;
    int iden_no1;
    int64_t best_score1;
    iden_no = 0;

    if ((band_right >= len2) || (band_left <= -len1) || (band_left > band_right)) return FAILED_FUNC;

    int band_width = band_right - band_left + 1;
    int band_width1 = 17;

    MatrixInt64& score_mat = buffer.score_mat;
    MatrixInt& back_mat = buffer.back_mat;

    int L = band_left, R = band_right;
    int kmin = (R < 0) ? -R : (L > 0) ? L : 0;
    int kmax = (R < len2 - len1) ? (R + 2 * len1) : (L > len2 - len1) ? (2 * len2 - L) : (len1 + len2);

    int lenY = kmax - kmin + 1;

    if (score_mat.size() <= lenY) {
        VectorInt row(band_width1, 0);
        VectorInt64 row2(band_width1, 0);
        while (score_mat.size() <= lenY) {
            score_mat.Append(row2);
            back_mat.Append(row);
        }
    }
    for (int i = 0; i <= lenY; i++) {
        if (score_mat[i].size < band_width1) score_mat[i].Resize(band_width1);
        if (back_mat[i].size < band_width1) back_mat[i].Resize(band_width1);
    }

    best_score = 0;

    auto t5 = Clock::now();

    // 初始化边界
    if (L < 0) {
        int T = (R < 0) ? R : 0;
        for (int X = L; X <= T; ++X) {
            int I = -X, J = 0;
            if (I < 0 || I > len1) continue;
            int n = (I + J) - kmin;
            int m = (X - L + 1) >> 1;
            score_mat[n][m] = (int64_t)mat.ext_gap * I;
            back_mat[n][m] = DP_BACK_TOP;
        }
        back_mat[(-T) - kmin][(T - L + 1) / 2] = DP_BACK_NONE;
    }

    if (R >= 0) {
        int T = (L > 0) ? L : 0;
        for (int X = T; X <= R; ++X) {
            int I = 0, J = X;
            if (J < 0 || J > len2) continue;
            int n = (I + J) - kmin;
            int m = (X - L + 1) >> 1;
            score_mat[n][m] = (int64_t)mat.ext_gap * J;
            back_mat[n][m] = DP_BACK_LEFT;
        }
        back_mat[T - kmin][(T - L + 1) / 2] = DP_BACK_NONE;
    }

    int gap_open[2] = {mat.gap, mat.ext_gap};
    int max_diag = band_center - band_left;
    int extra_score[4] = {4, 3, 2, 1};
    auto t2 = Clock::now();

    // ============ AVX-512 向量化变量定义 ============
    const int SIMD_WIDTH = 8;

    __m512i vec_0 = _mm512_setzero_si512();
    __m512i vec_1 = _mm512_set1_epi64(1);
    __m512i vec_3 = _mm512_set1_epi64(3);
    __m512i vec_4 = _mm512_set1_epi64(4);
    __m256i vec_DP_BACK_LEFT = _mm256_set1_epi32(DP_BACK_LEFT);
    __m256i vec_DP_BACK_TOP = _mm256_set1_epi32(DP_BACK_TOP);
    __m256i vec_DP_BACK_LEFT_TOP = _mm256_set1_epi32(DP_BACK_LEFT_TOP);
    __m512i vec_offset = _mm512_setr_epi64(0, 2, 4, 6, 8, 10, 12, 14);

    __m512i vec_len1 = _mm512_set1_epi64(len1);
    __m512i vec_len2 = _mm512_set1_epi64(len2);
    __m512i vec_lenY = _mm512_set1_epi64(lenY);
    __m512i vec_R = _mm512_set1_epi64(R);

    __m512i vec_band_center = _mm512_set1_epi64(band_center);

    __m512i vec_ext_gap = _mm512_set1_epi64(mat.ext_gap);
    __m512i vec_gap_open = _mm512_set1_epi64(mat.gap);

    // ============ 主循环 ============
    for (int y = kmin + 1; y <= kmax; y++) {
        int offset = (abs(y + L)) % 2;
        int x_start = L + offset;

        __m512i vec_y = _mm512_set1_epi64(y);
        __m512i vec_y_minus_kmin = _mm512_set1_epi64(y - kmin);

        int num_elements = 16;
        int vec_iterations = num_elements / SIMD_WIDTH;

        int index_y = y - kmin;
        int index_y1 = y - 1 - kmin;
        int index_y2 = y - 2 - kmin;

        if (y + 1 <= kmax) {
            _mm_prefetch((const char*)&score_mat[y + 1 - kmin][0], _MM_HINT_T0);
            _mm_prefetch((const char*)&back_mat[y + 1 - kmin][0], _MM_HINT_T0);
        }

        // ============ 向量化主循环 ============
        for (int vec_idx = 0; vec_idx < vec_iterations; vec_idx++) {
            int x_base = x_start + vec_idx * SIMD_WIDTH * 2;
            int index_x = (x_base - L + 1) >> 1;
            int index_x_L = (x_base - 1 - L + 1) >> 1;
            int index_x_R = (x_base + 1 - L + 1) >> 1;
            // cout<<"x_base:"<<x_base<<endl;

            __m512i vec_x = _mm512_add_epi64(_mm512_set1_epi64(x_base), vec_offset);

            __m512i vec_i = _mm512_sub_epi64(vec_y, vec_x);
            vec_i = _mm512_srai_epi64(vec_i, 1);
            __m512i vec_j = _mm512_add_epi64(vec_x, vec_y);
            vec_j = _mm512_srai_epi64(vec_j, 1);

            __mmask8 vec_i_gt_0 = _mm512_cmpgt_epi64_mask(vec_i, vec_0);
            __mmask8 vec_i_le_len1 = _mm512_cmple_epi64_mask(vec_i, vec_len1);
            __mmask8 vec_j_gt_0 = _mm512_cmpgt_epi64_mask(vec_j, vec_0);
            __mmask8 vec_j_le_len2 = _mm512_cmple_epi64_mask(vec_j, vec_len2);
            __mmask8 mask = vec_i_gt_0 & vec_i_le_len1 & vec_j_gt_0 & vec_j_le_len2;
            int count = _mm_popcnt_u32((unsigned int)mask);
            if (count == 0) continue;

#ifdef GATHER
            // printf(" (0x%02X)\n", (unsigned int)mask);
            __m512i vec_i_minus_1 = _mm512_sub_epi64(vec_i, vec_1);

            __m512i vec_j_minus_1 = _mm512_sub_epi64(vec_j, vec_1);

            vec_i_minus_1 = _mm512_mask_blend_epi64(mask, vec_0, vec_i_minus_1);
            vec_j_minus_1 = _mm512_mask_blend_epi64(mask, vec_0, vec_j_minus_1);

            __m256i vec_i_index = _mm512_cvtepi64_epi32(vec_i_minus_1);
            __m256i vec_j_index = _mm512_cvtepi64_epi32(vec_j_minus_1);

            __m256i seq1_i32 =
                _mm256_mmask_i32gather_epi32(_mm256_setzero_si256(), mask, vec_i_index, (const int*)iseq1, 1);

            __m256i seq2_i32 =
                _mm256_mmask_i32gather_epi32(_mm256_setzero_si256(), mask, vec_j_index, (const int*)iseq2, 1);
            __m256i seq1_vals = _mm256_and_si256(seq1_i32, _mm256_set1_epi32(0xFF));
            __m256i seq2_vals = _mm256_and_si256(seq2_i32, _mm256_set1_epi32(0xFF));
#else

            __m512i vec_i64_index = _mm512_sub_epi64(vec_i, vec_1);
            __m512i vec_j64_index = _mm512_sub_epi64(vec_j, vec_1);
            __m256i vec_i_index = _mm512_cvtepi64_epi32(vec_i64_index);
            __m256i vec_j_index = _mm512_cvtepi64_epi32(vec_j64_index);

            __m256i seq1_vals, seq2_vals;
            _mm256_load_char_array_backward(vec_i_index, iseq1, seq1_vals, mask);
            _mm256_load_char_array_forward(vec_j_index, iseq2, seq2_vals, mask);
#endif

            __m256i matrix_indices = _mm256_mullo_epi32(seq1_vals, _mm256_set1_epi32(MAX_AA));
            matrix_indices = _mm256_add_epi32(matrix_indices, seq2_vals);

            __m512i vec_sij =
                _mm512_mask_i32gather_epi64(_mm512_setzero_si512(), mask, matrix_indices, mat.flat_matrix, 8);

            // _mm512_store_epi64(i_arr, vec_i);
            // _mm512_store_epi64(j_arr, vec_j);
            // // double t1 =get_time();
            // for(int index=0; index<SIMD_WIDTH; index++){
            // 	int bit = (mask >> index) & 1;
            // 	if(!bit) sij_arr[index] = 0;
            // 	else sij_arr[index] = mat.matrix[iseq1[i_arr[index]-1]][iseq2[j_arr[index]-1]];
            // }
            // __m512i vec_sij = _mm512_load_si512(sij_arr);

            // double t2 = get_time();
            // t += t2 -t1;

            // int extra = extra_score[abs(x-band_center)&3];
            __m512i vec_extra = _mm512_sub_epi64(vec_x, vec_band_center);
            vec_extra = _mm512_abs_epi64(vec_extra);
            vec_extra = _mm512_and_si512(vec_extra, vec_3);
            vec_extra = _mm512_sub_epi64(vec_4, vec_extra);

            __mmask8 mask_sij = _mm512_cmpgt_epi64_mask(vec_sij, vec_0);
            mask_sij = mask_sij & mask;
            vec_extra = _mm512_mask_blend_epi64(mask_sij, vec_0, vec_extra);
            vec_sij = _mm512_add_epi64(vec_sij, vec_extra);

            __m512i vec_best_score1 = _mm512_loadu_si512(&score_mat[index_y][index_x]);
            __m256i vec_back = _mm256_loadu_si256((__m256i*)&back_mat[index_y][index_x]);

            // left-top
            if (y - 2 >= kmin) {
                __m512i vec_score_y2 = _mm512_loadu_si512(&score_mat[index_y2][index_x]);
                vec_best_score1 = _mm512_add_epi64(vec_score_y2, vec_sij);
                vec_back = _mm256_set1_epi32(DP_BACK_LEFT_TOP);
            }

            __mmask8 mask_gap = _mm512_cmpeq_epi64_mask(vec_i, vec_len1) | _mm512_cmpeq_epi64_mask(vec_j, vec_len2);
            __m512i vec_gap0 = _mm512_mask_blend_epi64(mask_gap, vec_gap_open, vec_ext_gap);
            __m512i vec_gap;
            __m512i vec_score;

            // left
            if (y - 1 >= kmin) {
                vec_gap = vec_gap0;
                __m512i vec_score_y1;
                __m256i vec_back_y1;

                vec_score_y1 = _mm512_loadu_si512(&score_mat[index_y1][index_x_L]);
                vec_back_y1 = _mm256_loadu_si256((__m256i*)&back_mat[index_y1][index_x_L]);

                __mmask8 gap_flag = _mm256_cmpeq_epi32_mask(vec_back_y1, vec_DP_BACK_LEFT);

                vec_gap = _mm512_mask_blend_epi64(gap_flag, vec_gap0, vec_ext_gap);

                vec_score = _mm512_add_epi64(vec_score_y1, vec_gap);

                __mmask8 modify_flag =
                    _mm512_cmpgt_epi64_mask(vec_score, vec_best_score1) & (0xFE | ((offset) | (vec_idx)));
                vec_best_score1 = _mm512_mask_blend_epi64(modify_flag, vec_best_score1, vec_score);
                vec_back = _mm256_mask_blend_epi32(modify_flag, vec_back, vec_DP_BACK_LEFT);
            }

            // top
            if (y - 1 >= kmin) {
                vec_gap = vec_gap0;
                __m512i vec_score_y1;
                __m256i vec_back_y1;
                vec_score_y1 = _mm512_loadu_si512(&score_mat[index_y1][index_x_R]);
                vec_back_y1 = _mm256_loadu_si256((__m256i*)&back_mat[index_y1][index_x_R]);

                __mmask8 gap_flag = _mm256_cmpeq_epi32_mask(vec_back_y1, vec_DP_BACK_TOP);
                vec_gap = _mm512_mask_blend_epi64(gap_flag, vec_gap0, vec_ext_gap);

                vec_score = _mm512_add_epi64(vec_score_y1, vec_gap);
                __mmask8 modify_flag = _mm512_cmpgt_epi64_mask(vec_score, vec_best_score1) &
                                       _mm512_cmple_epi64_mask(_mm512_add_epi64(vec_x, vec_1), vec_R);
                vec_best_score1 = _mm512_mask_blend_epi64(modify_flag, vec_best_score1, vec_score);
                vec_back = _mm256_mask_blend_epi32(modify_flag, vec_back, vec_DP_BACK_TOP);
            }
            _mm512_mask_storeu_epi64(&score_mat[index_y][index_x], mask, vec_best_score1);
            _mm256_mask_storeu_epi32((__m256i*)&back_mat[index_y][index_x], mask, vec_back);
        }
        // if (y >= kmin + 3) exit(0);
        // for(int x=L+(abs(y+L))%2;x<=R;x+=2){
        // 	// avx<<score_mat[y-kmin][(x-L+1)>>1]<<" ";
        // 	cout << score_mat[y-kmin][(x-L+1)>>1]<<" ";
        // }
        // // avx<<endl;
        // cout << endl;
        // if(y > kmin+20) exit(0);
    }
    auto t3 = Clock::now();
    // cout<<t<<endl;

    // ============ 回溯部分（保持原样）============
    x = (R < len2 - len1) ? R : (L > len2 - len1) ? L : len2 - len1;
    y = kmax;
    i = (-x + y) >> 1;
    j = (x + y) >> 1;
    // printf("\n(%d,%d)\n", i, j);
    best_score = score_mat[y - kmin][(x - L + 1) >> 1];
    best_score1 = score_mat[y - kmin][(x - L + 1) >> 1];
    printf("%2li(%2i) ", best_score1, iden_no1);

    int back = back_mat[y - kmin][(x - L + 1) >> 1];
    int last = back;

    int count = 0, count2 = 0, count3 = 0;
    int match, begin1, begin2, end1, end2;
    int gbegin1 = 0, gbegin2 = 0, gend1 = 0, gend2 = 0;
    int64_t score, smin = best_score1, smax = best_score1 - 1;
    int posmin, posmax, pos = 0;
    int bl, dlen = 0, dcount = 0;
    posmin = posmax = 0;
    begin1 = begin2 = end1 = end2 = 0;

    int masked = 0;
    int indels = 0;
    int max_indels = 0;

    while (back != DP_BACK_NONE) {
        switch (back) {
            case DP_BACK_TOP:
                bl = (last != back) & (j != 1) & (j != len2);
                dlen += bl;
                dcount += bl;
                score = score_mat[y - kmin][(x - L + 1) >> 1];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                i -= 1;
                x = x + 1;
                y = y - 1;
                break;
            case DP_BACK_LEFT:
                bl = (last != back) & (i != 1) & (i != len1);
                dlen += bl;
                dcount += bl;
                score = score_mat[y - kmin][(x - L + 1) >> 1];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                j -= 1;
                x = x - 1;
                y = y - 1;
                break;
            case DP_BACK_LEFT_TOP:
                if (alninfo && true) {
                    if (i == 1 || j == 1) {
                        gbegin1 = i - 1;
                        gbegin2 = j - 1;
                    } else if (i == len1 || j == len2) {
                        gend1 = i - 1;
                        gend2 = j - 1;
                    }
                }
                score = score_mat[y - kmin][(x - L + 1) >> 1];
                i -= 1;
                j -= 1;
                y -= 2;
                match = iseq1[i] == iseq2[j];
                if (score > smax) {
                    count = 0;
                    smax = score;
                    posmax = pos;
                    end1 = i;
                    end2 = j;
                }
                if (false && (iseq1[i] > 4 || iseq2[j] > 4)) {
                    masked += 1;
                } else {
                    dlen += 1;
                    dcount += !match;
                    count += match;
                    count2 += match;
                    count3 += match;
                }
                if (score < smin) {
                    int mm = match == 0;
                    count2 = 0;
                    smin = score;
                    posmin = pos - mm;
                    begin1 = i + mm;
                    begin2 = j + mm;
                }
                break;
            default:
                printf("%i\n", back);
                break;
        }
        pos += 1;
        last = back;
        back = back_mat[y - kmin][(x - L + 1) >> 1];
    }

    iden_no = true ? count3 : count - count2;
    alnln = posmin - posmax + 1 - masked;
    dist = dcount / (float)dlen;

    int umtail1 = len1 - 1 - end1;
    int umtail2 = len2 - 1 - end2;
    int umhead = begin1 < begin2 ? begin1 : begin2;
    int umtail = umtail1 < umtail2 ? umtail1 : umtail2;
    int umlen = umhead + umtail;

    auto t4 = Clock::now();

    auto d1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t5 - t1);
    auto d4 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t5);
    auto d2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2);
    auto d3 = std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3);

    cout << "Alloc Time: " << d1.count() << " ns\n";
    cout << "Init  Time: " << d4.count() << " ns\n";
    cout << "Vect  Time: " << d2.count() << " ns\n";
    cout << "Back  Time: " << d3.count() << " ns\n";

    if (umlen > 99999999) return FAILED_FUNC;
    if (umlen > len1 * 1.0) return FAILED_FUNC;
    if (umlen > len2 * 1.0) return FAILED_FUNC;

    if (alninfo) {
        alninfo[0] = begin1;
        alninfo[1] = end1;
        alninfo[2] = begin2;
        alninfo[3] = end2;
        alninfo[4] = masked;
        if (true) {
            alninfo[0] = gbegin1;
            alninfo[1] = gend1;
            alninfo[2] = gbegin2;
            alninfo[3] = gend2;
        }
    }
    return OK_FUNC;
}

#endif

int Init_Matrix_Compact(char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix& mat, int& best_score, int& iden_no,
                        int& alnln, float& dist, int* alninfo, int band_left, int band_center, int band_right,
                        WorkingBuffer& buffer) {
    if ((band_right >= len2) || (band_left <= -len1) || (band_left > band_right)) return FAILED_FUNC;
    int band_width = band_right - band_left + 1;
    int band_width1 = (band_width + 1) / 2 + 1;

    MatrixInt64& score_mat = buffer.score_mat;
    MatrixInt& back_mat = buffer.back_mat;

    int L = band_left, R = band_right;
    int kmin = (R < 0) ? -R : (L > 0) ? L : 0;
    int kmax = (R < len2 - len1) ? (R + 2 * len1) : (L > len2 - len1) ? (2 * len2 - L) : (len1 + len2);

    int lenY = kmax - kmin + 1;
    if (score_mat.size() <= lenY) {
        VectorInt row(band_width1, 0);
        VectorInt64 row2(band_width1, 0);
        while (score_mat.size() <= lenY) {
            score_mat.Append(row2);
            back_mat.Append(row);
        }
    }
    for (int i = 0; i <= lenY; i++) {
        if (score_mat[i].size < band_width1) score_mat[i].Resize(band_width1);
        if (back_mat[i].size < band_width1) back_mat[i].Resize(band_width1);
    }
    return OK_FUNC;
}

// FIXME:
// Modified after bug fixes,
//	completed the tight array implementation of rotation and anti-diagonal alignment
// Finished 20/10/2025 mgl
int rotation_compact_band_align(char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix& mat, int& best_score,
                                int& iden_no, int& alnln, float& dist, int* alninfo, int band_left, int band_center,
                                int band_right, WorkingBuffer& buffer) {
    std::ofstream noavx("noavx.tmp");
    int i, j, k, j1;
    int jj, kk;
    int x, y;
    int iden_no1;
    int64_t best_score1;
    iden_no = 0;
    if ((band_right >= len2) || (band_left <= -len1) || (band_left > band_right)) return FAILED_FUNC;
    int band_width = band_right - band_left + 1;
    int band_width1 = (band_width + 1) / 2 + 1;

    MatrixInt64& score_mat = buffer.score_mat;
    MatrixInt& back_mat = buffer.back_mat;

    int L = band_left, R = band_right;
    int kmin = (R < 0) ? -R : (L > 0) ? L : 0;
    int kmax = (R < len2 - len1) ? (R + 2 * len1) : (L > len2 - len1) ? (2 * len2 - L) : (len1 + len2);

    int lenY = kmax - kmin + 1;
    if (score_mat.size() <= lenY) {
        VectorInt row(band_width1, 0);
        VectorInt64 row2(band_width1, 0);
        while (score_mat.size() <= lenY) {
            score_mat.Append(row2);
            back_mat.Append(row);
        }
    }
    for (int i = 0; i <= lenY; i++) {
        if (score_mat[i].size < band_width1) score_mat[i].Resize(band_width1);
        if (back_mat[i].size < band_width1) back_mat[i].Resize(band_width1);
    }
    best_score = 0;
    if (L < 0) {
        int T = (R < 0) ? R : 0; // X = J - I = -I ∈ [L..T]
        for (int X = L; X <= T; ++X) {
            int I = -X, J = 0; // Y = I + J = -X
            if (I < 0 || I > len1) continue;
            int n = (I + J) - kmin; // = -X - Kmin
            int m = (X - L + 1) >> 1;
            score_mat[n][m] = (int64_t)mat.ext_gap * I;
            back_mat[n][m] = DP_BACK_TOP;
        }
        back_mat[(-T) - kmin][(T - L + 1) / 2] = DP_BACK_NONE;
    }

    if (R >= 0) {
        int T = (L > 0) ? L : 0; // X = J - 0 = J ∈ [T..R]
        for (int X = T; X <= R; ++X) {
            int I = 0, J = X; // Y = X
            if (J < 0 || J > len2) continue;
            int n = (I + J) - kmin; // = X - Kmin
            int m = (X - L + 1) >> 1;
            score_mat[n][m] = (int64_t)mat.ext_gap * J;
            back_mat[n][m] = DP_BACK_LEFT; // 从左延伸
        }
        back_mat[T - kmin][(T - L + 1) / 2] = DP_BACK_NONE;
    }

    int gap_open[2] = {mat.gap, mat.ext_gap};
    int max_diag = band_center - band_left;
    int extra_score[4] = {4, 3, 2, 1};

    // double t0 = get_time();

    for (int y = kmin + 1; y <= kmax; y++) {
        for (int x = L + (abs(y + L)) % 2; x <= R; x += 2) {
            // int x = (index<<1)+offset+L;
            int i = (y - x) >> 1, j = (x + y) >> 1;
            if (i < 0 || i > len1 || j < 0 || j > len2) continue;
            if (i == 0 || j == 0) continue;
            int ci = iseq1[i - 1];
            int cj = iseq2[j - 1];
            int sij = mat.matrix[ci][cj];
            int s1, k0, back;
            int extra = extra_score[abs(x - band_center) & 3];
            sij += extra * (sij > 0);

            back = back_mat[y - kmin][(x - L + 1) >> 1];
            best_score1 = score_mat[y - kmin][(x - L + 1) >> 1];

            // try (x,y-2)/(i-1,j-1) the top left conerner of the original array
            if (y - 2 >= kmin) {
                best_score1 = score_mat[y - 2 - kmin][(x - L + 1) >> 1] + sij;
                back = DP_BACK_LEFT_TOP;
            }

            int gap0 = gap_open[(i == len1) | (j == len2)];
            int gap = 0;
            int64_t score;
            // try (x-1,y-1)/(i,j-1) the left of original array
            if (x - 1 >= L && y - 1 >= kmin) {
                gap = gap0;

                if (back_mat[y - 1 - kmin][(x - 1 - L + 1) >> 1] == DP_BACK_LEFT) gap = mat.ext_gap;
                // cout<<score_mat[y-1-kmin][(x-1-L+1)>>1]+gap<<" ";
                if ((score = score_mat[y - 1 - kmin][(x - 1 - L + 1) >> 1] + gap) > best_score1) {
                    back = DP_BACK_LEFT;
                    best_score1 = score;
                }
            }
            // try (x+1,y-1)/(i-1,j) the top of original array
            if (x + 1 <= R && y - 1 >= kmin) {
                gap = gap0;
                if (back_mat[y - 1 - kmin][(x + 1 - L + 1) >> 1] == DP_BACK_TOP) gap = mat.ext_gap;
                if ((score = score_mat[y - 1 - kmin][(x + 1 - L + 1) >> 1] + gap) > best_score1) {
                    back = DP_BACK_TOP;
                    best_score1 = score;
                }
            }
            score_mat[y - kmin][(x - L + 1) >> 1] = best_score1;
            back_mat[y - kmin][(x - L + 1) >> 1] = back;
            // cout<<best_score1<<" ";
        }

        // // if(y <= kmax - 20) continue;
        // for(int x=L+(abs(y+L))%2;x<=R;x+=2){
        // 	// noavx<<score_mat[y-kmin][(x-L+1)>>1]<<" ";
        // 	cout<<score_mat[y-kmin][(x-L+1)>>1]<<" ";
        // }
        // // noavx<<endl;
        // cout<<endl;
        // if(y > kmin+20) return 0;
    }

    // double t1 = get_time();
    // cout << "Alignment Time: "<<t1-t0<<" seconds" << endl;

    x = (R < len2 - len1) ? R : (L > len2 - len1) ? L : len2 - len1;
    y = kmax;
    i = (-x + y) >> 1, j = (x + y) >> 1;
    // printf("\n(%d,%d)\n",i,j);
    best_score = score_mat[y - kmin][(x - L + 1) >> 1];
    best_score1 = score_mat[y - kmin][(x - L + 1) >> 1];
    // cout<<best_score<<endl;
    // printf("%2li(%2i) ",best_score1,iden_no1);

    int back = back_mat[y - kmin][(x - L + 1) >> 1];
    int last = back;

    // printf("\n[DEBUG INFO](%d,%d)\n",back,last);

    int count = 0, count2 = 0, count3 = 0;
    int match, begin1, begin2, end1, end2;
    int gbegin1 = 0, gbegin2 = 0, gend1 = 0, gend2 = 0;
    int64_t score, smin = best_score1, smax = best_score1 - 1;
    int posmin, posmax, pos = 0;
    int bl, dlen = 0, dcount = 0;
    posmin = posmax = 0;
    begin1 = begin2 = end1 = end2 = 0;

    int masked = 0;
    int indels = 0;
    int max_indels = 0;
    while (back != DP_BACK_NONE) {
        switch (back) {
            case DP_BACK_TOP:
                // update from (x+1,y-1)/(i-1,j)
                bl = (last != back) & (j != 1) & (j != len2);
                dlen += bl;
                dcount += bl;
                score = score_mat[y - kmin][(x - L + 1) >> 1];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                i -= 1;
                x = x + 1;
                y = y - 1;
                break;
            case DP_BACK_LEFT:
                // update from (x-1,y-1)/(i,j-1)
                bl = (last != back) & (i != 1) & (i != len1);
                dlen += bl;
                dcount += bl;
                score = score_mat[y - kmin][(x - L + 1) >> 1];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                j -= 1;
                x = x - 1;
                y = y - 1;
                break;
            case DP_BACK_LEFT_TOP:
                // update from (x,y-2)/(i-1,j-1)
                if (alninfo && true) {
                    if (i == 1 || j == 1) {
                        gbegin1 = i - 1;
                        gbegin2 = j - 1;
                    } else if (i == len1 || j == len2) {
                        gend1 = i - 1;
                        gend2 = j - 1;
                    }
                }
                score = score_mat[y - kmin][(x - L + 1) >> 1];
                i -= 1;
                j -= 1;
                y -= 2;
                match = iseq1[i] == iseq2[j];
                if (score > smax) {
                    count = 0;
                    smax = score;
                    posmax = pos;
                    end1 = i;
                    end2 = j;
                }
                if (false && (iseq1[i] > 4 || iseq2[j] > 4)) {
                    masked += 1;
                } else {
                    dlen += 1;
                    dcount += !match;
                    count += match;
                    count2 += match;
                    count3 += match;
                }
                if (score < smin) {
                    int mm = match == 0;
                    count2 = 0;
                    smin = score;
                    posmin = pos - mm;
                    begin1 = i + mm;
                    begin2 = j + mm;
                }
                break;
            default:
                printf("%i/n", back);
                break;
        }
        pos += 1;
        last = back;
        back = back_mat[y - kmin][(x - L + 1) >> 1];
    }
    // double t2 = get_time();
    // cout<<"BackTrace Time: "<<t2-t1<<" seconds"<<endl;
    iden_no = true ? count3 : count - count2;
    alnln = posmin - posmax + 1 - masked;
    // printf("\n[DEBUG INFO] posmin:%d posmax:%d masked:%d\n",posmin,posmax,masked);
    dist = dcount / (float)dlen;
    // dist = - 0.75 * log( 1.0 - dist * 4.0 / 3.0 );
    int umtail1 = len1 - 1 - end1;
    int umtail2 = len2 - 1 - end2;
    int umhead = begin1 < begin2 ? begin1 : begin2;
    int umtail = umtail1 < umtail2 ? umtail1 : umtail2;
    int umlen = umhead + umtail;
    if (umlen > 99999999) return FAILED_FUNC;
    if (umlen > len1 * 1.0) return FAILED_FUNC;
    if (umlen > len2 * 1.0) return FAILED_FUNC;
    if (alninfo) {
        alninfo[0] = begin1;
        alninfo[1] = end1;
        alninfo[2] = begin2;
        alninfo[3] = end2;
        alninfo[4] = masked;
        if (true) {
            alninfo[0] = gbegin1;
            alninfo[1] = gend1;
            alninfo[2] = gbegin2;
            alninfo[3] = gend2;
        }
    }
    return OK_FUNC;
}

// FIXME:
// An implementation of a simple rotating array calculation.
// Finished 13/10/2025 mgl

int Init_Matrix_Rotation(char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix& mat, int& best_score,
                         int& iden_no, int& alnln, float& dist, int* alninfo, int band_left, int band_center,
                         int band_right, WorkingBuffer& buffer) {
    if ((band_right >= len2) || (band_left <= -len1) || (band_left > band_right)) return FAILED_FUNC;
    int band_width = band_right - band_left + 1;
    int band_width1 = band_width + 1;
    MatrixInt64& score_mat = buffer.score_mat;
    MatrixInt& back_mat = buffer.back_mat;

    int L = band_left, R = band_right;
    int kmin = (R < 0) ? -R : (L > 0) ? L : 0;
    int kmax = (R < len2 - len1) ? (R + 2 * len1) : (L > len2 - len1) ? (2 * len2 - L) : (len1 + len2);

    int lenY = kmax - kmin + 1;
    if (score_mat.size() <= lenY) {
        VectorInt row(band_width1, 0);
        VectorInt64 row2(band_width1, 0);
        while (score_mat.size() <= lenY) {
            score_mat.Append(row2);
            back_mat.Append(row);
        }
    }
    for (int i = 0; i <= lenY; i++) {
        if (score_mat[i].size < band_width1) score_mat[i].Resize(band_width1);
        if (back_mat[i].size < band_width1) back_mat[i].Resize(band_width1);
    }
    return OK_FUNC;
}

int rotation_band_align(char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix& mat, int& best_score, int& iden_no,
                        int& alnln, float& dist, int* alninfo, int band_left, int band_center, int band_right,
                        WorkingBuffer& buffer) {
    int i, j, k, j1;
    int jj, kk;
    int x, y;
    int iden_no1;
    int64_t best_score1;
    iden_no = 0;
    if ((band_right >= len2) || (band_left <= -len1) || (band_left > band_right)) return FAILED_FUNC;
    int band_width = band_right - band_left + 1;
    int band_width1 = band_width + 1;
    MatrixInt64& score_mat = buffer.score_mat;
    MatrixInt& back_mat = buffer.back_mat;

    int L = band_left, R = band_right;
    int kmin = (R < 0) ? -R : (L > 0) ? L : 0;
    int kmax = (R < len2 - len1) ? (R + 2 * len1) : (L > len2 - len1) ? (2 * len2 - L) : (len1 + len2);

    int lenY = kmax - kmin + 1;
    if (score_mat.size() <= lenY) {
        VectorInt row(band_width1, 0);
        VectorInt64 row2(band_width1, 0);
        while (score_mat.size() <= lenY) {
            score_mat.Append(row2);
            back_mat.Append(row);
        }
    }
    for (int i = 0; i <= lenY; i++) {
        if (score_mat[i].size < band_width1) score_mat[i].Resize(band_width1);
        if (back_mat[i].size < band_width1) back_mat[i].Resize(band_width1);
    }
    best_score = 0;
    //	For coordinate Mapping,first swicth to the coordinate form of diagonal and anti-diagonal,
    //		and then offset it for the convenience of memory storage.
    //	The initial grid coordinates are (i,j), the diagonal coordinates are (x,y), and the subsequent storage
    // coordinates are (n,m) 	Then there is the following coordinate conversion form: 		x = (-i+j), y = (i+j) ||
    // i = 1/2*(-x+y), j = 1/2*(x+y) 		n = y-kmin, m = x-L  ||  y = n+kmin, x = m+L
    if (L < 0) {
        int T = (R < 0) ? R : 0; // X = J - I = -I ∈ [L..T]
        for (int X = L; X <= T; ++X) {
            int I = -X, J = 0; // Y = I + J = -X
            if (I < 0 || I > len1) continue;
            int n = (I + J) - kmin; // = -X - Kmin
            int m = X - L;
            score_mat[n][m] = (int64_t)mat.ext_gap * I;
            back_mat[n][m] = DP_BACK_TOP; // 从上延伸
        }
        // 回溯终止点
        back_mat[(-T) - kmin][T - L] = DP_BACK_NONE;
    }

    if (R >= 0) {
        int T = (L > 0) ? L : 0; // X = J - 0 = J ∈ [T..R]
        for (int X = T; X <= R; ++X) {
            int I = 0, J = X; // Y = X
            if (J < 0 || J > len2) continue;
            int n = (I + J) - kmin; // = X - Kmin
            int m = X - L;
            score_mat[n][m] = (int64_t)mat.ext_gap * J;
            back_mat[n][m] = DP_BACK_LEFT; // 从左延伸
        }
        back_mat[T - kmin][T - L] = DP_BACK_NONE;
    }

    int gap_open[2] = {mat.gap, mat.ext_gap};
    int max_diag = band_center - band_left;
    int extra_score[4] = {4, 3, 2, 1};

    for (int y = kmin + 1; y <= kmax; y++) {
        for (int x = L + (abs(y + L)) % 2; x <= R; x += 2) {
            int i = (y - x) >> 1, j = (x + y) >> 1;
            if (i < 0 || i > len1 || j < 0 || j > len2) continue;
            if (i == 0 || j == 0) continue;
            int ci = iseq1[i - 1];
            int cj = iseq2[j - 1];
            int sij = mat.matrix[ci][cj];
            int s1, k0, back;
            int extra = extra_score[abs(x - band_center) & 3];
            sij += extra * (sij > 0);

            back = back_mat[y - kmin][x - L];
            best_score1 = score_mat[y - kmin][x - L];

            // try (x,y-2)/(i-1,j-1) the top left conerner of the original array
            if (y - 2 >= kmin) {
                best_score1 = score_mat[y - 2 - kmin][x - L] + sij;
                back = DP_BACK_LEFT_TOP;
            }
            int gap0 = gap_open[(i == len1) | (j == len2)];
            int gap = 0;
            int64_t score;
            // try (x-1,y-1)/(i,j-1) the left of original array
            if (x - 1 >= L && y - 1 >= kmin) {
                gap = gap0;
                if (back_mat[y - 1 - kmin][x - 1 - L] == DP_BACK_LEFT) gap = mat.ext_gap;
                if ((score = score_mat[y - 1 - kmin][x - 1 - L] + gap) > best_score1) {
                    back = DP_BACK_LEFT;
                    best_score1 = score;
                }
            }
            // try (x+1,y-1)/(i-1,j) the top of original array
            if (x + 1 <= R && y - 1 >= kmin) {
                gap = gap0;
                if (back_mat[y - 1 - kmin][x + 1 - L] == DP_BACK_TOP) gap = mat.ext_gap;
                if ((score = score_mat[y - 1 - kmin][x + 1 - L] + gap) > best_score1) {
                    back = DP_BACK_TOP;
                    best_score1 = score;
                }
            }
            score_mat[y - kmin][x - L] = best_score1;
            back_mat[y - kmin][x - L] = back;
        }
    }
    x = (R < len2 - len1) ? R : (L > len2 - len1) ? L : len2 - len1;
    y = kmax;
    i = (-x + y) / 2, j = (x + y) / 2;
    // printf("\n(%d,%d)\n",i,j);
    best_score = score_mat[y - kmin][x - L];
    best_score1 = score_mat[y - kmin][x - L];
    // cout<<best_score<<endl;
    // printf("%2i(%2i) ",best_score1,iden_no1);

#if 1
    const char* letters = "acgtnx";
    const char* letters2 = "ACGTNX";
#else
    const char* letters = "arndcqeghilkmfpstwyvbzx";
    const char* letters2 = "ARNDCQEGHILKMFPSTWYVBZX";
#endif
    int back = back_mat[y - kmin][x - L];
    int last = back;

    // printf("\n[DEBUG INFO](%d,%d)\n",back,last);

    int count = 0, count2 = 0, count3 = 0;
    int match, begin1, begin2, end1, end2;
    int gbegin1 = 0, gbegin2 = 0, gend1 = 0, gend2 = 0;
    int64_t score, smin = best_score1, smax = best_score1 - 1;
    int posmin, posmax, pos = 0;
    int bl, dlen = 0, dcount = 0;
    posmin = posmax = 0;
    begin1 = begin2 = end1 = end2 = 0;

#ifdef PRINT
#define PRINT
    printf("%i %i\n", best_score, score_mat[i][j1]);
    printf("%i %i %i\n", band_left, band_center, band_right);
    printf("%i %i %i %i\n", i, j, j1, len2);
#endif

#ifdef MAKEALIGN
#define MAKEALIGN
    char AA[MAX_SEQ], BB[MAX_SEQ];
    int NN = 0;
    int IA, IB;
    for (IA = len1; IA > i; IA--) {
        AA[NN] = letters[iseq1[IA - 1]];
        BB[NN++] = '-';
    }
    for (IB = len2; IB > j; IB--) {
        AA[NN] = '-';
        BB[NN++] = letters[iseq2[IB - 1]];
    }
#endif

    int masked = 0;
    int indels = 0;
    int max_indels = 0;
    while (back != DP_BACK_NONE) {
        switch (back) {
            case DP_BACK_TOP:
                // update from (x+1,y-1)/(i-1,j)
#ifdef PRINT
                printf("%5i: %c %c %9i\n", pos, letters[iseq1[i]], '|', score_mat[y - 1 - kmin][x + 1 - L]);
#endif
#ifdef MAKEALIGN
                AA[NN] = letters[iseq1[i]];
                BB[NN++] = '-';
#endif
                bl = (last != back) & (j != 1) & (j != len2);
                dlen += bl;
                dcount += bl;
                score = score_mat[y - kmin][x - L];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                i -= 1;
                x = x + 1;
                y = y - 1;
                break;
            case DP_BACK_LEFT:
                // update from (x-1,y-1)/(i,j-1)
#ifdef PRINT
                printf("%5i: %c %c %9i\n", pos, '|', letters[iseq2[j]], score_mat[y - 1 - kmin][x - 1 - L]);
#endif
#ifdef MAKEALIGN
                AA[NN] = '-';
                BB[NN++] = letters[iseq2[j]];
#endif
                bl = (last != back) & (i != 1) & (i != len1);
                dlen += bl;
                dcount += bl;
                score = score_mat[y - kmin][x - L];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                j -= 1;
                x = x - 1;
                y = y - 1;
                break;
            case DP_BACK_LEFT_TOP:
                // update from (x,y-2)/(i-1,j-1)
#ifdef PRINT
                if (iseq1[i] == iseq2[j]) {
                    printf("%5i: %c %c %9i\n", pos, letters2[iseq1[i]], letters2[iseq2[j]],
                           score_mat[y - 2 - kmin][x - L]);
                } else {
                    printf("%5i: %c %c %9i\n", pos, letters[iseq1[i]], letters[iseq2[j]],
                           score_mat[y - 2 - kmin][x - L]);
                }
#endif
#ifdef MAKEALIGN
                if (iseq1[i] == iseq2[j]) {
                    AA[NN] = letters2[iseq1[i]];
                    BB[NN++] = letters2[iseq2[j]];
                } else {
                    AA[NN] = letters[iseq1[i]];
                    BB[NN++] = letters[iseq2[j]];
                }
#endif
                if (alninfo && true) {
                    if (i == 1 || j == 1) {
                        gbegin1 = i - 1;
                        gbegin2 = j - 1;
                    } else if (i == len1 || j == len2) {
                        gend1 = i - 1;
                        gend2 = j - 1;
                    }
                }
                score = score_mat[y - kmin][x - L];
                i -= 1;
                j -= 1;
                y -= 2;
                match = iseq1[i] == iseq2[j];
                if (score > smax) {
                    count = 0;
                    smax = score;
                    posmax = pos;
                    end1 = i;
                    end2 = j;
                }
                if (false && (iseq1[i] > 4 || iseq2[j] > 4)) {
                    masked += 1;
                } else {
                    dlen += 1;
                    dcount += !match;
                    count += match;
                    count2 += match;
                    count3 += match;
                }
                if (score < smin) {
                    int mm = match == 0;
                    count2 = 0;
                    smin = score;
                    posmin = pos - mm;
                    begin1 = i + mm;
                    begin2 = j + mm;
                }
                break;
            default:
                printf("%i/n", back);
                break;
        }
        pos += 1;
        last = back;
        back = back_mat[y - kmin][x - L];
    }
    iden_no = true ? count3 : count - count2;
    alnln = posmin - posmax + 1 - masked;
    // printf("\n[DEBUG INFO] posmin:%d posmax:%d masked:%d\n",posmin,posmax,masked);
    dist = dcount / (float)dlen;
    // dist = - 0.75 * log( 1.0 - dist * 4.0 / 3.0 );
    int umtail1 = len1 - 1 - end1;
    int umtail2 = len2 - 1 - end2;
    int umhead = begin1 < begin2 ? begin1 : begin2;
    int umtail = umtail1 < umtail2 ? umtail1 : umtail2;
    int umlen = umhead + umtail;
    if (umlen > 99999999) return FAILED_FUNC;
    if (umlen > len1 * 1.0) return FAILED_FUNC;
    if (umlen > len2 * 1.0) return FAILED_FUNC;
    if (alninfo) {
        alninfo[0] = begin1;
        alninfo[1] = end1;
        alninfo[2] = begin2;
        alninfo[3] = end2;
        alninfo[4] = masked;
        if (true) {
            alninfo[0] = gbegin1;
            alninfo[1] = gend1;
            alninfo[2] = gbegin2;
            alninfo[3] = gend2;
        }
    }
#ifdef PRINT
    printf("%6i %6i:  %4i %4i %4i %4i\n", alnln, iden_no, begin1, end1, begin2, end2);
    printf("%6i %6i:  %4i %4i\n", posmin, posmax, posmin - posmax, count - count2);
    printf("smin = %9i, smax = %9i\n", smin, smax);
    printf("dlen = %5i, dcount = %5i, dist = %.3f\n", dlen, dcount, dcount / (float)dlen);
#endif
#ifdef MAKEALIGN
    float identity = iden_no / (float)(options.global_identity ? (len1 - masked) : alnln);
    if (identity < options.cluster_thd) return OK_FUNC;
    while (i >= 0) {
        AA[NN] = letters[iseq1[i]];
        BB[NN++] = '-';
    }
    while (j >= 0) {
        AA[NN] = '-';
        BB[NN++] = letters[iseq2[j]];
    }
    AA[NN] = '\0';
    BB[NN] = '\0';

    for (i = 0; i < NN / 2; i++) {
        char aa = AA[i], bb = BB[i];
        AA[i] = AA[NN - i - 1];
        BB[i] = BB[NN - i - 1];
        AA[NN - i - 1] = aa;
        BB[NN - i - 1] = bb;
    }
    static int fcount = 0;
    fcount += 1;
    FILE* fout = fopen("alignments.txt", "a");
    if (fout == NULL) {
        if (fcount <= 1) printf("alignment files open failed\n");
        return OK_FUNC;
    }
    fprintf(fout, "\n\n######################################################\n");
    fprintf(fout, "# length X = %i\n", len2);
    fprintf(fout, "# length Y = %i\n", len1);
    fprintf(fout, "# best align X: %i-%i\n", begin2 + 1, end2 + 1);
    fprintf(fout, "# best align Y: %i-%i\n", begin1 + 1, end1 + 1);
    if (alninfo) {
        fprintf(fout, "# align X: %i-%i\n", alninfo[2] + 1, alninfo[3] + 1);
        fprintf(fout, "# align Y: %i-%i\n", alninfo[0] + 1, alninfo[1] + 1);
    }
    fprintf(fout, "# alignment length: %i\n", alnln);
    fprintf(fout, "# identity count: %i\n", iden_no);
    fprintf(fout, "# identity: %g\n", identity);
    fprintf(fout, "# distance: %g\n", dist);
    if (options.is454) fprintf(fout, "# max indel: %i\n", max_indels);
#if 0
	fprintf( fout, "%i %s\n", seq1->index, AA );
	fprintf( fout, "%i %s\n", seq2->index, BB );
#else
    bool printaa = true;
    IB = IA = 0;
    fprintf(fout, "\n\nX ");
    while (IA < NN) {
        if (printaa) {
            fprintf(fout, "%c", BB[IB]);
            IB += 1;
            if (IB % 75 == 0 or IB == NN) printaa = false, fprintf(fout, "\nY ");
        } else {
            fprintf(fout, "%c", AA[IA]);
            IA += 1;
            if (IA % 75 == 0) printaa = true, fprintf(fout, "\n\nX ");
        }
    }
#endif
    fclose(fout);
#endif
    return OK_FUNC;
}

int Init_Matrix(char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix& mat, int& best_score, int& iden_no,
                int& alnln, float& dist, int* alninfo, int band_left, int band_center, int band_right,
                WorkingBuffer& buffer) {
    if ((band_right >= len2) || (band_left <= -len1) || (band_left > band_right)) return FAILED_FUNC;

    // allocate mem for score_mat[len1][len2] etc
    int band_width = band_right - band_left + 1;
    int band_width1 = band_width + 1;

    // score_mat, back_mat [i][j]: i index of seqi (0 to len(seqi)-1), j index of band (0 to band_width-1)
    MatrixInt64& score_mat = buffer.score_mat;
    MatrixInt& back_mat = buffer.back_mat;

    // printf( "%i  %i\n", band_right, band_left );

    if (score_mat.size() <= len1) {
        VectorInt row(band_width1, 0);
        VectorInt64 row2(band_width1, 0);
        while (score_mat.size() <= len1) {
            score_mat.Append(row2);
            back_mat.Append(row);
        }
    }
    for (int i = 0; i <= len1; i++) {
        if (score_mat[i].Size() < band_width1) score_mat[i].Resize(band_width1);
        if (back_mat[i].Size() < band_width1) back_mat[i].Resize(band_width1);
    }
    return OK_FUNC;
}

int local_band_align(char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix& mat, int& best_score, int& iden_no,
                     int& alnln, float& dist, int* alninfo, int band_left, int band_center, int band_right,
                     WorkingBuffer& buffer) {
    int i, j, k, j1;
    int jj, kk;
    int iden_no1;
    int64_t best_score1;
    iden_no = 0;

    if ((band_right >= len2) || (band_left <= -len1) || (band_left > band_right)) return FAILED_FUNC;

    // allocate mem for score_mat[len1][len2] etc
    int band_width = band_right - band_left + 1;
    int band_width1 = band_width + 1;

    // score_mat, back_mat [i][j]: i index of seqi (0 to len(seqi)-1), j index of band (0 to band_width-1)
    MatrixInt64& score_mat = buffer.score_mat;
    MatrixInt& back_mat = buffer.back_mat;

    // printf( "%i  %i\n", band_right, band_left );

    if (score_mat.size() <= len1) {
        VectorInt row(band_width1, 0);
        VectorInt64 row2(band_width1, 0);
        while (score_mat.size() <= len1) {
            score_mat.Append(row2);
            back_mat.Append(row);
        }
    }
    for (i = 0; i <= len1; i++) {
        if (score_mat[i].Size() < band_width1) score_mat[i].Resize(band_width1);
        if (back_mat[i].Size() < band_width1) back_mat[i].Resize(band_width1);
    }

    best_score = 0;
    /* seq1 is query, seq2 is rep
                  seq2    len2 = 17       seq2    len2 = 17    seq2    len2 = 17
                  01234567890123456       01234567890123456    01234567890123456
       0          xxxxxxxxxxxxxxxxx \\\\\\XXXxxxxxxxxxxxxxx    xXXXXXXXxxxxxxxxx
       1     \\\\\Xxxxxxxxxxxxxxxxx  \\\\\Xxx\xxxxxxxxxxxxx    xx\xxxxx\xxxxxxxx
       2      \\\\X\xxxxxxxxxxxxxxx   \\\\Xxxx\xxxxxxxxxxxx    xxx\xxxxx\xxxxxxx
  seq1 3       \\\Xx\xxxxxxxxxxxxxx    \\\Xxxxx\xxxxxxxxxxx    xxxx\xxxxx\xxxxxx
  len1 4        \\Xxx\xxxxxxxxxxxxx     \\Xxxxxx\xxxxxxxxxx    xxxxx\xxxxx\xxxxx
  = 11 5         \Xxxx\xxxxxxxxxxxx      \Xxxxxxx\xxxxxxxxx    xxxxxx\xxxxx\xxxx
       6          Xxxxx\xxxxxxxxxxx       Xxxxxxxx\xxxxxxxx    xxxxxxx\xxxxx\xxx
       7          x\xxxx\xxxxxxxxxx       x\xxxxxxx\xxxxxxx    xxxxxxxx\xxxxx\xx
       8          xx\xxxx\xxxxxxxxx       xx\xxxxxxx\xxxxxx    xxxxxxxxx\xxxxx\x
       9          xxx\xxxx\xxxxxxxx       xxx\xxxxxxx\xxxxx    xxxxxxxxxx\xxxxx\
       0          xxxx\xxxx\xxxxxxx       xxxx\xxxxxxx\xxxx    xxxxxxxxxxx\xxxxx
                  band_left < 0 (-6)      band_left < 0 (-6)   band_left >=0 (1)
                  band_right < 0 (-1)     band_right >=0 (2)   band_right >=0(7)
                  band_width 6            band_width 9         band_width 7
       init score_mat, and iden_mat (place with upper 'X')
     */

    if (band_left < 0) { // set score to left border of the matrix within band
        int tband = (band_right < 0) ? band_right : 0;
        // for (k=band_left; k<tband; k++)
        for (k = band_left; k <= tband; k++) { // fixed on 2006 11 14
            i = -k;
            j1 = k - band_left;
            // penalty for leading gap opening = penalty for gap extension
            // each of the left side query hunging residues give ext_gap (-1)
            score_mat[i][j1] = mat.ext_gap * i;
            back_mat[i][j1] = DP_BACK_TOP;
        }
        back_mat[-tband][tband - band_left] = DP_BACK_NONE;
    }

    if (band_right >= 0) { // set score to top border of the matrix within band
        int tband = (band_left > 0) ? band_left : 0;
        for (j = tband; j <= band_right; j++) {
            j1 = j - band_left;
            score_mat[0][j1] = mat.ext_gap * j;
            back_mat[0][j1] = DP_BACK_LEFT;
        }
        back_mat[0][tband - band_left] = DP_BACK_NONE;
    }

    int gap_open[2] = {mat.gap, mat.ext_gap};
    int max_diag = band_center - band_left;
    int extra_score[4] = {4, 3, 2, 1};

    for (i = 1; i <= len1; i++) {   // 次数
        int J0 = 1 - band_left - i; // 从下网上
        int J1 = len2 - band_left - i;
        if (J0 < 0) J0 = 0;
        if (J1 >= band_width) J1 = band_width;
        for (j1 = J0; j1 <= J1; j1++) { // 次数
            j = j1 + i + band_left;

            int ci = iseq1[i - 1];
            int cj = iseq2[j - 1];
            int sij = mat.matrix[ci][cj]; // sij 是当前两个字符比对的得分
            // int iden_ij = (ci == cj);
            int s1, k0, back;

            /* extra score according to the distance to the best diagonal */
            int extra = extra_score[abs(j1 - max_diag) & 3]; // max distance 3   模 4
            sij += extra * (sij > 0);

            back = DP_BACK_LEFT_TOP;
            best_score1 = score_mat[i - 1][j1] + sij;
            int gap0 = gap_open[(i == len1) | (j == len2)];
            int gap = 0;
            int64_t score;

            if (j1 > 0) {
                gap = gap0;
                if (back_mat[i][j1 - 1] == DP_BACK_LEFT) gap = mat.ext_gap;
                if ((score = score_mat[i][j1 - 1] + gap) > best_score1) {
                    back = DP_BACK_LEFT;
                    best_score1 = score;
                }
            }
            if (j1 + 1 < band_width) {
                gap = gap0;
                if (back_mat[i - 1][j1 + 1] == DP_BACK_TOP) gap = mat.ext_gap;
                if ((score = score_mat[i - 1][j1 + 1] + gap) > best_score1) {
                    back = DP_BACK_TOP;
                    best_score1 = score;
                }
            }
            score_mat[i][j1] = best_score1;
            back_mat[i][j1] = back;
            // printf( "%2i(%2i) ", best_score1, iden_no1 );
        }
        // printf( "%2i(%2i) ", best_score1, iden_no1 );
        // printf( "\n" );
    }
    // printf( "%2i(%2i) ", best_score1, iden_no1 );
    i = j = 0;
    if (len2 - band_left < len1) {
        i = len2 - band_left;
        j = len2;
    } else if (len1 + band_right < len2) {
        i = len1;
        j = len1 + band_right;
    } else {
        i = len1;
        j = len2;
    }
    printf("\n(%d,%d)\n", i, j);
    j1 = j - i - band_left;
    best_score = score_mat[i][j1];
    best_score1 = score_mat[i][j1];
    // cout<<best_score<<endl;
    printf("%2i(%2i) ", best_score1, iden_no1);
#if 1
    const char* letters = "acgtnx";
    const char* letters2 = "ACGTNX";
#else
    const char* letters = "arndcqeghilkmfpstwyvbzx";
    const char* letters2 = "ARNDCQEGHILKMFPSTWYVBZX";
#endif
    int back = back_mat[i][j1];
    int last = back;

    // printf("\n[DEBUG INFO](%d,%d)\n",back,last);

    int count = 0, count2 = 0, count3 = 0;
    int match, begin1, begin2, end1, end2;
    int gbegin1 = 0, gbegin2 = 0, gend1 = 0, gend2 = 0;
    int64_t score, smin = best_score1, smax = best_score1 - 1;
    int posmin, posmax, pos = 0;
    int bl, dlen = 0, dcount = 0;
    posmin = posmax = 0;
    begin1 = begin2 = end1 = end2 = 0;

#ifdef PRINT
#define PRINT
    printf("%i %i\n", best_score, score_mat[i][j1]);
    printf("%i %i %i\n", band_left, band_center, band_right);
    printf("%i %i %i %i\n", i, j, j1, len2);
#endif

#ifdef MAKEALIGN
#define MAKEALIGN
    char AA[MAX_SEQ], BB[MAX_SEQ];
    int NN = 0;
    int IA, IB;
    for (IA = len1; IA > i; IA--) {
        AA[NN] = letters[iseq1[IA - 1]];
        BB[NN++] = '-';
    }
    for (IB = len2; IB > j; IB--) {
        AA[NN] = '-';
        BB[NN++] = letters[iseq2[IB - 1]];
    }
#endif

    int masked = 0;
    int indels = 0;
    int max_indels = 0;
    while (back != DP_BACK_NONE) {
        switch (back) {
            case DP_BACK_TOP:
#ifdef PRINT
                printf("%5i: %c %c %9i\n", pos, letters[iseq1[i - 1]], '|', score_mat[i][j1]);
#endif
#ifdef MAKEALIGN
                AA[NN] = letters[iseq1[i - 1]];
                BB[NN++] = '-';
#endif
                bl = (last != back) & (j != 1) & (j != len2);
                dlen += bl;
                dcount += bl;
                score = score_mat[i][j1];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                i -= 1;
                j1 += 1;
                break;
            case DP_BACK_LEFT:
#ifdef PRINT
                printf("%5i: %c %c %9i\n", pos, '|', letters[iseq2[j - 1]], score_mat[i][j1]);
#endif
#ifdef MAKEALIGN
                AA[NN] = '-';
                BB[NN++] = letters[iseq2[j - 1]];
#endif
                bl = (last != back) & (i != 1) & (i != len1);
                dlen += bl;
                dcount += bl;
                score = score_mat[i][j1];
                if (score < smin) {
                    count2 = 0;
                    smin = score;
                    posmin = pos - 1;
                    begin1 = i;
                    begin2 = j;
                }
                j1 -= 1;
                j -= 1;
                break;
            case DP_BACK_LEFT_TOP:
#ifdef PRINT
                if (iseq1[i - 1] == iseq2[j - 1]) {
                    printf("%5i: %c %c %9i\n", pos, letters2[iseq1[i - 1]], letters2[iseq2[j - 1]], score_mat[i][j1]);
                } else {
                    printf("%5i: %c %c %9i\n", pos, letters[iseq1[i - 1]], letters[iseq2[j - 1]], score_mat[i][j1]);
                }
#endif
#ifdef MAKEALIGN
                if (iseq1[i - 1] == iseq2[j - 1]) {
                    AA[NN] = letters2[iseq1[i - 1]];
                    BB[NN++] = letters2[iseq2[j - 1]];
                } else {
                    AA[NN] = letters[iseq1[i - 1]];
                    BB[NN++] = letters[iseq2[j - 1]];
                }
#endif
                if (alninfo && true) {
                    if (i == 1 || j == 1) {
                        gbegin1 = i - 1;
                        gbegin2 = j - 1;
                    } else if (i == len1 || j == len2) {
                        gend1 = i - 1;
                        gend2 = j - 1;
                    }
                }
                score = score_mat[i][j1];
                i -= 1;
                j -= 1;
                match = iseq1[i] == iseq2[j];
                if (score > smax) {
                    count = 0;
                    smax = score;
                    posmax = pos;
                    end1 = i;
                    end2 = j;
                }
                if (false && (iseq1[i] > 4 || iseq2[j] > 4)) {
                    masked += 1;
                } else {
                    dlen += 1;
                    dcount += !match;
                    count += match;
                    count2 += match;
                    count3 += match;
                }
                if (score < smin) {
                    int mm = match == 0;
                    count2 = 0;
                    smin = score;
                    posmin = pos - mm;
                    begin1 = i + mm;
                    begin2 = j + mm;
                }
                break;
            default:
                printf("%i\n", back);
                break;
        }

        pos += 1;
        last = back;
        back = back_mat[i][j1];
    }
    iden_no = true ? count3 : count - count2;
    alnln = posmin - posmax + 1 - masked;
    // printf("\n[DEBUG INFO] posmin:%d posmax:%d masked:%d\n",posmin,posmax,masked);
    dist = dcount / (float)dlen;
    // dist = - 0.75 * log( 1.0 - dist * 4.0 / 3.0 );
    int umtail1 = len1 - 1 - end1;
    int umtail2 = len2 - 1 - end2;
    int umhead = begin1 < begin2 ? begin1 : begin2;
    int umtail = umtail1 < umtail2 ? umtail1 : umtail2;
    int umlen = umhead + umtail;
    if (umlen > 99999999) return FAILED_FUNC;
    if (umlen > len1 * 1.0) return FAILED_FUNC;
    if (umlen > len2 * 1.0) return FAILED_FUNC;
    if (alninfo) {
        alninfo[0] = begin1;
        alninfo[1] = end1;
        alninfo[2] = begin2;
        alninfo[3] = end2;
        alninfo[4] = masked;
        if (true) {
            alninfo[0] = gbegin1;
            alninfo[1] = gend1;
            alninfo[2] = gbegin2;
            alninfo[3] = gend2;
        }
    }
#ifdef PRINT
    printf("%6i %6i:  %4i %4i %4i %4i\n", alnln, iden_no, begin1, end1, begin2, end2);
    printf("%6i %6i:  %4i %4i\n", posmin, posmax, posmin - posmax, count - count2);
    printf("smin = %9i, smax = %9i\n", smin, smax);
    printf("dlen = %5i, dcount = %5i, dist = %.3f\n", dlen, dcount, dcount / (float)dlen);
#endif
#ifdef MAKEALIGN
    float identity = iden_no / (float)(options.global_identity ? (len1 - masked) : alnln);
    if (identity < options.cluster_thd) return OK_FUNC;
    while (i--) {
        AA[NN] = letters[iseq1[i - 1]];
        BB[NN++] = '-';
    }
    while (j--) {
        AA[NN] = '-';
        BB[NN++] = letters[iseq2[j - 1]];
    }
    AA[NN] = '\0';
    BB[NN] = '\0';
    for (i = 0; i < NN / 2; i++) {
        char aa = AA[i], bb = BB[i];
        AA[i] = AA[NN - i - 1];
        BB[i] = BB[NN - i - 1];
        AA[NN - i - 1] = aa;
        BB[NN - i - 1] = bb;
    }
    static int fcount = 0;
    fcount += 1;
    FILE* fout = fopen("alignments.txt", "a");
    if (fout == NULL) {
        if (fcount <= 1) printf("alignment files open failed\n");
        return OK_FUNC;
    }
    fprintf(fout, "\n\n######################################################\n");
    fprintf(fout, "# length X = %i\n", len2);
    fprintf(fout, "# length Y = %i\n", len1);
    fprintf(fout, "# best align X: %i-%i\n", begin2 + 1, end2 + 1);
    fprintf(fout, "# best align Y: %i-%i\n", begin1 + 1, end1 + 1);
    if (alninfo) {
        fprintf(fout, "# align X: %i-%i\n", alninfo[2] + 1, alninfo[3] + 1);
        fprintf(fout, "# align Y: %i-%i\n", alninfo[0] + 1, alninfo[1] + 1);
    }
    fprintf(fout, "# alignment length: %i\n", alnln);
    fprintf(fout, "# identity count: %i\n", iden_no);
    fprintf(fout, "# identity: %g\n", identity);
    fprintf(fout, "# distance: %g\n", dist);
    if (options.is454) fprintf(fout, "# max indel: %i\n", max_indels);
#if 0
	fprintf( fout, "%i %s\n", seq1->index, AA );
	fprintf( fout, "%i %s\n", seq2->index, BB );
#else
    bool printaa = true;
    IB = IA = 0;
    fprintf(fout, "\n\nX ");
    while (IA < NN) {
        if (printaa) {
            fprintf(fout, "%c", BB[IB]);
            IB += 1;
            if (IB % 75 == 0 or IB == NN) printaa = false, fprintf(fout, "\nY ");
        } else {
            fprintf(fout, "%c", AA[IA]);
            IA += 1;
            if (IA % 75 == 0) printaa = true, fprintf(fout, "\n\nX ");
        }
    }
#endif
    fclose(fout);
#endif

    return OK_FUNC;
} // END int local_band_align

int main(int argc, char* argv[]) {
    gzFile fp1;
    kseq_t* ks1;
    int max_seq_len = 0;
    CLI::App app{"build word table v.0.0.1, build word table from fasta protein files"};
    string filename = "";
    auto option_input = app.add_option("-i, --input", filename, "input file name, fasta or gziped fasta formats");
    option_input->required();
    CLI11_PARSE(app, argc, argv);
    fp1 = gzopen(filename.c_str(), "r");
    ks1 = kseq_init(fp1);
    vector<Sequence> seqs;
    //   cerr<<"11111"<<endl;
    while (1) {
        int length = kseq_read(ks1);
        if (length < 0) break;

        max_seq_len = max<int>(max_seq_len, ks1->seq.l);

        Sequence seq;
        seq.name = ks1->name.s;
        // seq.comment = ks1->comment.s;
        seq.seq = ks1->seq.s;
        seq.len = length;
        if (seq.seq.size() < 50) continue;
        //  cerr<<"11111"<<endl;
        char* seq_data = seq.seq.data();
        for (int i = 0; i < seq.seq.size(); i++) {
            seq_data[i] = aa2idx[seq_data[i] - 'A'];
        }
        // cerr<<"11111"<<endl;
        seqs.emplace_back(seq);

        // number_seqs++;

        // if(max_num_seqs > 0 && number_seqs >= max_num_seqs) break;
    }
    // exit(0);
    // const char seq1[] = {1, 2, 1, 2, 1};  // 假设的蛋白质序列1
    // const char seq2[] = {0, 1, 2, 1, 4};  // 假设的蛋白质序列2
    int band_width = 20, required_aa1 = 1851;

    WorkingBuffer buffer(max_seq_len); // 创建工作区

    char* seq1 = seqs[1].seq.data();
    char* seq2 = seqs[0].seq.data();
    int len1 = seqs[1].len;
    int len2 = seqs[0].len;

    auto t0 = Clock::now();
    // 计算二元对的倒排表
    ComputeAAP(seq1, len1, buffer);

    auto t1 = Clock::now();

    // 计算最佳带宽和匹配结果
    int best_sum, band_left, band_center, band_right, best_score, tiden_no, alnln;
    int talign_info[5];
    float tiden_pc, distance = 0;
    diag_test_aapn(seq2, len1, len2, buffer, best_sum, band_width, band_left, band_center, band_right, required_aa1);

    auto t2 = Clock::now();

    // diag_no_table(seq1,len1,seq2,len2,buffer, best_sum, band_width, band_left, band_center, band_right,
    // required_aa1); int required_aa2 = 31740; if ( best_sum < required_aa2 ) exit(0);
    int rc = FAILED_FUNC;
    int tmp = FAILED_FUNC;

    // double t2 = get_time();

    cout << "Best sum: " << best_sum << endl;
    cout << "Band left: " << band_left << ", Band center: " << band_center << ", Band right: " << band_right << endl;

#ifdef ORIGIN
    tmp = Init_Matrix(seq1, seq2, len1, len2, mat, best_score, tiden_no, alnln, distance, talign_info, band_left,
                      band_center, band_right, buffer);
#endif

#ifdef COMPACT
    tmp = Init_Matrix_Compact(seq1, seq2, len1, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                              band_left, band_center, band_right, buffer);
#endif

#ifdef ROTATION
    tmp = Init_Matrix_Rotation(seq1, seq2, len1, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                               band_left, band_center, band_right, buffer);
#endif

#ifdef __AVX512F__
#if !defined(ORIGIN) && !defined(COMPACT)
    tmp = Init_Matrix_AVX512(seq1, seq2, len1, len2, mat, best_score, tiden_no, alnln, distance, talign_info, band_left,
                             band_center, band_right, buffer);
#endif
#endif

#ifdef __AVX2__
#if !defined(ORIGIN) && !defined(COMPACT) && !defined(__AVX512F__)

    tmp = Init_Matrix_AVX2(seq1, seq2, len1, len2, mat, best_score, tiden_no, alnln, distance, talign_info, band_left,
                           band_center, band_right, buffer);
#endif
#endif

    auto t3 = Clock::now();

#ifdef ORIGIN
    cout << "\n<Method Info> Local Band Align" << endl;
    rc = local_band_align(seq1, seq2, len1, len2, mat, best_score, tiden_no, alnln, distance, talign_info, band_left,
                          band_center, band_right, buffer);
#endif

#ifdef COMPACT
    cout << "\n<Method Info> Rotation Compact Band Align" << endl;
    rc = rotation_compact_band_align(seq1, seq2, len1, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                     band_left, band_center, band_right, buffer);
#endif

#ifdef ROTATION
    cout << "\n<Method Info> Rotation Band Align" << endl;
    rc = rotation_band_align(seq1, seq2, len1, len2, mat, best_score, tiden_no, alnln, distance, talign_info, band_left,
                             band_center, band_right, buffer);
#endif

#ifdef __AVX512F__
#if !defined(ORIGIN) && !defined(COMPACT) && !defined(ROTATION)
    cout << "\n<Method Info> AVX512 Band Align" << endl;
    for (int i = 1; i <= 1; i++)
        rc = rotation_band_align_AVX512(seq1, seq2, len1, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                        band_left, band_center, band_right, buffer);
#endif
#endif

#ifdef __AVX2__
#if !defined(ORIGIN) && !defined(COMPACT) && !defined(__AVX512F__)
    cout << "\n<Method Info> AVX2 Band Align" << endl;
    rc = rotation_band_align_AVX2(seq1, seq2, len1, len2, mat, best_score, tiden_no, alnln, distance, talign_info,
                                  band_left, band_center, band_right, buffer);
#endif
#endif

    auto t4 = Clock::now();

    cout << "best_score " << best_score << endl;
    cout << "tiden_no  " << tiden_no << endl;
    cout << "alnln  " << alnln << endl;
    cout << "distance  " << distance << endl;

    auto d0 = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0);
    auto d1 = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
    auto d2 = std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2);
    auto d3 = std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3);

    cout << "[Time Info] ComputeAPP\tTime: " << d0.count() << " ns\n";
    cout << "[Time Info] diag_test_aapn\tTime: " << d1.count() << " ns\n";
    cout << "[Time Info] Init_Matrix\tTime: " << d2.count() << " ns\n";
    cout << "[Time Info] band_align\tTime: " << d3.count() << " ns\n";

    // double t3 = get_time();
    // cerr << "\nlocal_band_align time : " << t3 - t2 << " seconds" << endl;
    return 0;
}
