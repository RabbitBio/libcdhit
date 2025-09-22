#include <iostream>
#include <vector>
#include <algorithm>
#include "kseq.h"
#include "CLI11.hpp"
#include <zlib.h>
#include <sys/time.h>
#define NAA1 21      // 假设我们有 5 种氨基酸
#define MAX_SEQ 655360
#define MAX_AA 23
#define NAA2 (NAA1 * NAA1) // 二元对的总数
#define MAX_DIAG (MAX_SEQ<<1)    
#define FAILED_FUNC 1
#define OK_FUNC 0
using namespace std;
template<class TYPE>
class Vector : public vector<TYPE>
{
	public:
		Vector() : vector<TYPE>(){}
		Vector( size_t size ) : vector<TYPE>( size ){}
		Vector( size_t size, const TYPE & deft ) : vector<TYPE>( size, deft ){}

		void Append( const TYPE & item ){
			size_t n = this->size();
			if( n + 1 >= this->capacity() ) this->reserve( n + n/5 + 1 );
			this->push_back( item );
		}
		int size()const{ return (int)vector<TYPE>::size(); }
};

template<class TYPE>
class NVector
{
	public:
		TYPE   *items;
		int     size;
		int     capacity;

		NVector(){ size = capacity = 0; items = NULL; }
		NVector( int n, const TYPE & v=TYPE() ){ 
			size = capacity = 0; items = NULL; 
			Resize( n, v );
		}
		NVector( const NVector & other ){
			size = capacity = 0; items = NULL; 
			if( other.items ){
				Resize( other.size );
				memcpy( items, other.items, other.size * sizeof(TYPE) );
			}
		}

		~NVector(){ if( items ) free( items ); }

		int  Size()const{ return size; }
		void Clear(){
			if( items ) free( items );
			size = capacity = 0; items = NULL; 
		}

		void Resize( int n, const TYPE & value=TYPE() ){
			if( n == size && capacity > 0 ) return;
			int i;
			// When resize() is called, probably this is the intended size,
			// and will not be changed frequently.
			if( n != capacity ){
				capacity = n;
				items = (TYPE*)realloc( items, capacity*sizeof(TYPE) );
			}
			for(i=size; i<n; i++ ) items[i] = value;
			size = n;
		}
		void Append( const TYPE & item ){
			if( size + 1 >= capacity ){
				capacity = size + size/5 + 1;
				items = (TYPE*)realloc( items, capacity*sizeof(TYPE) );
			}
			items[size] = item;
			size ++;
		}

		TYPE& operator[]( const int i ){
			//if( i <0 or i >= size ) printf( "out of range\n" );
			return items[i];
		}
		TYPE& operator[]( const int i )const{
			//if( i <0 or i >= size ) printf( "out of range\n" );
			return items[i];
		}
};
typedef NVector<int64_t>   VectorInt64;
typedef Vector<VectorInt64> MatrixInt64;
typedef NVector<int>      VectorInt;
typedef Vector<VectorInt> MatrixInt;
enum { DP_BACK_NONE=0, DP_BACK_LEFT_TOP=1, DP_BACK_LEFT=2, DP_BACK_TOP=3 };

KSEQ_INIT(gzFile, gzread)
int aa2idx[] = {0, 2, 4, 3, 6, 13,7, 8, 9,20,11,10,12, 2,20,14,
                5, 1,15,16,20,19,17,20,18, 6};

int BLOSUM62[] = {
  4,                                                                  // A
 -1, 5,                                                               // R
 -2, 0, 6,                                                            // N
 -2,-2, 1, 6,                                                         // D
  0,-3,-3,-3, 9,                                                      // C
 -1, 1, 0, 0,-3, 5,                                                   // Q
 -1, 0, 0, 2,-4, 2, 5,                                                // E
  0,-2, 0,-1,-3,-2,-2, 6,                                             // G
 -2, 0, 1,-1,-3, 0, 0,-2, 8,                                          // H
 -1,-3,-3,-3,-1,-3,-3,-4,-3, 4,                                       // I
 -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,                                    // L
 -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,                                 // K
 -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5,                              // M
 -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,                           // F
 -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,                        // P
  1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4,                     // S
  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,                  // T
 -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11,               // W
 -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,            // Y
  0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,         // V
 -2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4,      // B
 -1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,   // Z
  0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1 // X
//A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X
//0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19  2  6 20
};
int BLOSUM62_na[] = {
  2,                  // A
 -2, 2,               // C
 -2,-2, 2,            // G
 -2,-2,-2, 2,         // T
 -2,-2,-2, 1, 2,      // U
 -2,-2,-2,-2,-2, 1,   // N
  0, 0, 0, 0, 0, 0, 1 // X
//A  C  G  T  U  N  X
//0  1  2  3  3  4  5
};
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
    int len;
};
// WorkingBuffer结构体定义
struct WorkingBuffer {
    vector<int> taap;     // 用于存储二元对出现次数
    vector<int> aap_begin; // 存储每个二元对的起始位置
    vector<int> aap_list;  // 存储实际的二元对出现位置
    vector<int> diag_score; // 存储对角线命中次数
    vector<int> diag_score2; // 存储加权命中次数
    MatrixInt64  score_mat;
    MatrixInt    back_mat;
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
class ScoreMatrix { //Matrix
	private:

	public:
		int matrix[MAX_AA][MAX_AA];
		int gap, ext_gap;

		ScoreMatrix();
		void init();
		void set_gap(int gap1, int ext_gap1);
		void set_matrix(int *mat1);
		void set_to_na();
		void set_match( int score );
		void set_mismatch( int score );
}; // END class ScoreMatrix

ScoreMatrix mat;
/////////////////
ScoreMatrix::ScoreMatrix()
{
	init();
}

void ScoreMatrix::init() 
{
	set_gap( -11, -1 );
	set_matrix( BLOSUM62 );
}

void ScoreMatrix::set_gap(int gap1, int ext_gap1)
{
	int i;
	gap = MAX_SEQ * gap1;
	ext_gap = MAX_SEQ * ext_gap1;
}

void ScoreMatrix::set_matrix(int *mat1)
{
	int i, j, k;
	k = 0;
	for ( i=0; i<MAX_AA; i++)
		for ( j=0; j<=i; j++)
			matrix[j][i] = matrix[i][j] = MAX_SEQ * mat1[ k++ ];
}

void ScoreMatrix::set_to_na()
{
	set_gap( -6, -1 );
	set_matrix( BLOSUM62_na );
}
// Only for est
void ScoreMatrix::set_match( int score )
{
	int i;
	for ( i=0; i<5; i++) matrix[i][i] = MAX_SEQ * score;
	//matrix[3][4] = matrix[4][3] = MAX_SEQ * score;
}
// Only for est
void ScoreMatrix::set_mismatch( int score )
{
	int i, j;
	for ( i=0; i<MAX_AA; i++)
		for ( j=0; j<i; j++)
			matrix[j][i] = matrix[i][j] = MAX_SEQ * score;
	matrix[3][4] = matrix[4][3] = MAX_SEQ;
}

// 计算二元对的倒排表
void ComputeAAP(const char *seqi, int size, WorkingBuffer &buffer) {
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
int diag_test_aapn( const char *iseq2, int len1, int len2, WorkingBuffer &buffer,
                   int &best_sum, int band_width, int &band_left, int &band_center, int &band_right,
                   int required_aa1) {
    
    int i, i1, j, k;
    int nall = len1 + len2 - 1; // 总对角线数
    vector<int> &taap = buffer.taap;
    vector<int> &aap_begin = buffer.aap_begin;
    vector<int> &aap_list = buffer.aap_list;
    vector<int> &diag_score = buffer.diag_score;
    vector<int> &diag_score2 = buffer.diag_score2;

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
    int emax = int((1.0 - 0.9) * mlen) + 1;  // 假设 options.cluster_thd = 0.5
    for (j = from; j < imax_diag; j++) {
        if ((imax_diag - j) > emax || diag_score[j] < 1) {
            best_score -= diag_score[j];
            from++;
        } else break;
    }
    for (j = end; j > imax_diag; j--) {
        if ((j - imax_diag) > emax || diag_score[j] < 1) {
            best_score -= diag_score[j];
            end--;
        } else break;
    }

    // 计算最终带宽
    band_left = from - len1 + 1;
    band_right = end - len1 + 1;
    band_center = imax_diag - len1 + 1;
    best_sum = best_score;

    return 0;
}

int diag_no_table(
	const char* seqi, const size_t sizei,
	const char* seqj, const size_t sizej,
	WorkingBuffer &buffer,
    int &best_sum, int band_width, int &band_left, int &band_center, int &band_right,
    int required_aa1
){
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
    auto& diag_score  = buffer.diag_score;
    auto& diag_score2 = buffer.diag_score2;

    // --- 构造 (code, pos) 列表，同时统计频次（用于前缀和） ---
    std::vector<std::pair<int,int>> seqi_list; seqi_list.reserve(len1);
    std::vector<std::pair<int,int>> seqj_list; seqj_list.reserve(len2);

    // 注意：把 char 转成 unsigned，避免负值
    std::vector<int> cnt_i(NAA2, 0), cnt_j(NAA2, 0);

    for (int i = 0; i < len1; ++i) {
        int a = static_cast<unsigned char>(seqi[i]);
        int b = static_cast<unsigned char>(seqi[i + 1]);
        int code = a * NAA1 + b;           // 0..NAA2-1
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
            const int jpos = seqj_list[jj].second;               // seqj 的二元对起点
            const int i1   = (len1 - 1) + jpos;                  // 对角线坐标偏移
            const int cpx  = 1 + (seqj[jpos] != seqj[jpos + 1]); // 加权：相同=1，不同=2
            for (size_t ii = i; ii < i_end; ++ii) {
                const int ipos = seqi_list[ii].second;           // seqi 的二元对起点
                const int d    = i1 - ipos;                      // 对角线下标 0..nall-1
                if ((unsigned)d < (unsigned)nall) {
                    diag_score[d]  += 1;
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
    int emax = int((1.0 - 0.9) * mlen) + 1;  // 假设 options.cluster_thd = 0.5
    for (j = from; j < imax_diag; j++) {
        if ((imax_diag - j) > emax || diag_score[j] < 1) {
            best_score -= diag_score[j];
            from++;
        } else break;
    }
    for (j = end; j > imax_diag; j--) {
        if ((j - imax_diag) > emax || diag_score[j] < 1) {
            best_score -= diag_score[j];
            end--;
        } else break;
    }

    // 计算最终带宽
    band_left = from - len1 + 1;
    band_right = end - len1 + 1;
    band_center = imax_diag - len1 + 1;
    best_sum = best_score;

    return 0;
}


int local_band_align( char iseq1[], char iseq2[], int len1, int len2, ScoreMatrix &mat, 
		int &best_score, int &iden_no, int &alnln, float &dist, int *alninfo,
		int band_left, int band_center, int band_right, WorkingBuffer & buffer)
{
	int i, j, k, j1;
	int jj, kk;
	int iden_no1;
	int64_t best_score1;
	iden_no = 0;

	if ( (band_right >= len2 ) ||
			(band_left  <= -len1) ||
			(band_left  > band_right) ) return FAILED_FUNC;

	// allocate mem for score_mat[len1][len2] etc
	int band_width = band_right - band_left + 1;
	int band_width1 = band_width + 1;

    // score_mat, back_mat [i][j]: i index of seqi (0 to len(seqi)-1), j index of band (0 to band_width-1)
	MatrixInt64 & score_mat = buffer.score_mat;
	MatrixInt   & back_mat = buffer.back_mat;

	//printf( "%i  %i\n", band_right, band_left );

	if( score_mat.size() <= len1 ){
		VectorInt   row( band_width1, 0 );
		VectorInt64 row2( band_width1, 0 );
		while( score_mat.size() <= len1 ){
			score_mat.Append( row2 );
			back_mat.Append( row );
		}
	}
	for(i=0; i<=len1; i++){
		if( score_mat[i].Size() < band_width1 ) score_mat[i].Resize( band_width1 );
		if( back_mat[i].Size() < band_width1 ) back_mat[i].Resize( band_width1 );
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

	if (band_left < 0) {  //set score to left border of the matrix within band
		int tband = (band_right < 0) ? band_right : 0;
		//for (k=band_left; k<tband; k++)
		for (k=band_left; k<=tband; k++) { // fixed on 2006 11 14
			i = -k;
			j1 = k-band_left;
			// penalty for leading gap opening = penalty for gap extension
            // each of the left side query hunging residues give ext_gap (-1)
			score_mat[i][j1] =  mat.ext_gap * i;
			back_mat[i][j1] = DP_BACK_TOP;
		}
		back_mat[-tband][tband-band_left] = DP_BACK_NONE;
	}

	if (band_right >=0) { //set score to top border of the matrix within band
		int tband = (band_left > 0) ? band_left : 0;
		for (j=tband; j<=band_right; j++) {
			j1 = j-band_left;
			score_mat[0][j1] = mat.ext_gap * j;
			back_mat[0][j1] = DP_BACK_LEFT;
		}
		back_mat[0][tband-band_left] = DP_BACK_NONE;
	}

	int gap_open[2] = { mat.gap, mat.ext_gap };
	int max_diag = band_center - band_left;
	int extra_score[4] = { 4, 3, 2, 1 };


	for (i=1; i<=len1; i++) {//次数
		int J0 = 1 - band_left - i;//从下网上
		int J1 = len2 - band_left - i;
		if( J0 < 0 ) J0 = 0;
		if( J1 >= band_width ) J1 = band_width;
		for (j1=J0; j1<=J1; j1++){//次数
			j = j1+i+band_left;

			int ci = iseq1[i-1];
			int cj = iseq2[j-1];
			int sij = mat.matrix[ci][cj];//sij 是当前两个字符比对的得分
			//int iden_ij = (ci == cj);
			int s1, k0, back;

			/* extra score according to the distance to the best diagonal */
			int extra = extra_score[ abs(j1 - max_diag) & 3 ]; // max distance 3   模 4
			sij += extra * (sij>0);

			back = DP_BACK_LEFT_TOP;
			best_score1 = score_mat[i-1][j1] + sij;
			int gap0 = gap_open[ (i == len1) | (j == len2) ];
			int gap = 0;
			int64_t score;

			if( j1 > 0 ){
				gap = gap0;
				if( back_mat[i][j1-1] == DP_BACK_LEFT ) gap = mat.ext_gap;
				if( (score = score_mat[i][j1-1] + gap) > best_score1 ){
					back = DP_BACK_LEFT;
					best_score1 = score;
				}
			}
			if(j1+1<band_width){
				gap = gap0;
				if( back_mat[i-1][j1+1] == DP_BACK_TOP ) gap = mat.ext_gap;
				if( (score = score_mat[i-1][j1+1] + gap) > best_score1 ){
					back = DP_BACK_TOP;
					best_score1 = score;
				}
			}
			score_mat[i][j1] = best_score1;
			back_mat[i][j1]  = back;
			// printf( "%2i(%2i) ", best_score1, iden_no1 );

		}
        // printf( "%2i(%2i) ", best_score1, iden_no1 );
		//printf( "\n" );
	}
    // printf( "%2i(%2i) ", best_score1, iden_no1 );
	i = j = 0;
	if( len2 - band_left < len1 ){
		i = len2 - band_left;
		j = len2;
	}else if( len1 + band_right < len2 ){
		i = len1;
		j = len1 + band_right;
	}else{
		i = len1;
		j = len2;
	}
	j1 = j - i - band_left;
	best_score = score_mat[i][j1];
	best_score1 = score_mat[i][j1];
      printf( "%2i(%2i) ", best_score1, iden_no1 );
#if 1
	const char *letters = "acgtnx";
	const char *letters2 = "ACGTNX";
#else
	const char *letters = "arndcqeghilkmfpstwyvbzx";
	const char *letters2 = "ARNDCQEGHILKMFPSTWYVBZX";
#endif
	int back = back_mat[i][j1];
	int last = back;
	int count = 0, count2 = 0, count3 = 0;
	int match, begin1, begin2, end1, end2;
	int gbegin1=0, gbegin2=0, gend1=0, gend2=0;
	int64_t score, smin = best_score1, smax = best_score1 - 1;
	int posmin, posmax, pos = 0;
	int bl, dlen = 0, dcount = 0;
	posmin = posmax = 0;
	begin1 = begin2 = end1 = end2 = 0;

#ifdef PRINT
#define PRINT
	printf( "%i %i\n", best_score, score_mat[i][j1] );
	printf( "%i %i %i\n", band_left, band_center, band_right );
	printf( "%i %i %i %i\n", i, j, j1, len2 );
#endif
#ifdef MAKEALIGN
#define MAKEALIGN
	char AA[ MAX_SEQ ], BB[ MAX_SEQ ];
	int NN = 0;
	int IA, IB;
	for(IA=len1;IA>i;IA--){
		AA[NN] = letters[ iseq1[IA-1] ];
		BB[NN++] = '-';
	}
	for(IB=len2;IB>j;IB--){
		AA[NN] = '-';
		BB[NN++] = letters[ iseq2[IB-1] ];
	}
#endif

	int masked = 0;
	int indels = 0;
	int max_indels = 0;
	while( back != DP_BACK_NONE ){
		switch( back ){
		case DP_BACK_TOP  :
#ifdef PRINT
			printf( "%5i: %c %c %9i\n", pos, letters[ iseq1[i-1] ], '|', score_mat[i][j1] );
#endif
#ifdef MAKEALIGN
			AA[NN] = letters[ iseq1[i-1] ];
			BB[NN++] = '-';
#endif
			bl = (last != back) & (j != 1) & (j != len2);
			dlen += bl;
			dcount += bl;
			score = score_mat[i][j1];
			if( score < smin ){
				count2 = 0;
				smin = score;
				posmin = pos - 1;
				begin1 = i;
				begin2 = j;
			}
			i -= 1;
			j1 += 1;
			break;
		case DP_BACK_LEFT :
#ifdef PRINT
			printf( "%5i: %c %c %9i\n", pos, '|', letters[ iseq2[j-1] ], score_mat[i][j1] );
#endif
#ifdef MAKEALIGN
			AA[NN] = '-';
			BB[NN++] = letters[ iseq2[j-1] ];
#endif
			bl = (last != back) & (i != 1) & (i != len1);
			dlen += bl;
			dcount += bl;
			score = score_mat[i][j1];
			if( score < smin ){
				count2 = 0;
				smin = score;
				posmin = pos - 1;
				begin1 = i;
				begin2 = j;
			}
			j1 -= 1;
			j -= 1;
			break;
		case DP_BACK_LEFT_TOP :
#ifdef PRINT
			if( iseq1[i-1] == iseq2[j-1] ){
				printf( "%5i: %c %c %9i\n", pos, letters2[ iseq1[i-1] ], letters2[ iseq2[j-1] ], score_mat[i][j1] );
			}else{
				printf( "%5i: %c %c %9i\n", pos, letters[ iseq1[i-1] ], letters[ iseq2[j-1] ], score_mat[i][j1] );
			}
#endif
#ifdef MAKEALIGN
			if( iseq1[i-1] == iseq2[j-1] ){
				AA[NN] = letters2[ iseq1[i-1] ];
				BB[NN++] = letters2[ iseq2[j-1] ];
			}else{
				AA[NN] = letters[ iseq1[i-1] ];
				BB[NN++] = letters[ iseq2[j-1] ];
			}
#endif
			if( alninfo && true ){
				if( i == 1 || j == 1 ){
					gbegin1 = i-1;
					gbegin2 = j-1;
				}else if( i == len1 || j == len2 ){
					gend1 = i-1;
					gend2 = j-1;
				}
			}
			score = score_mat[i][j1];
			i -= 1;
			j -= 1;
			match = iseq1[i] == iseq2[j];
			if( score > smax ){
				count = 0;
				smax = score;
				posmax = pos;
				end1 = i;
				end2 = j;
			}
			if( false && (iseq1[i] > 4 || iseq2[j] > 4) ){
				masked += 1;
			}else{
				dlen += 1;
				dcount += ! match;
				count += match;
				count2 += match;
				count3 += match;
			}
			if( score < smin ){
				int mm = match == 0;
				count2 = 0;
				smin = score;
				posmin = pos - mm;
				begin1 = i + mm;
				begin2 = j + mm;
			}
			break;
		default : printf( "%i\n", back ); break;
		}

		pos += 1;
		last = back;
		back = back_mat[i][j1];
	}
	iden_no = true ? count3 : count - count2;
	alnln = posmin - posmax + 1 - masked;
	dist = dcount/(float)dlen;
	//dist = - 0.75 * log( 1.0 - dist * 4.0 / 3.0 );
	int umtail1 = len1 - 1 - end1;
	int umtail2 = len2 - 1 - end2;
	int umhead = begin1 < begin2 ? begin1 : begin2;
	int umtail = umtail1 < umtail2 ? umtail1 : umtail2;
	int umlen = umhead + umtail;
	if( umlen > 99999999 ) return FAILED_FUNC;
	if( umlen > len1 * 1.0) return FAILED_FUNC;
	if( umlen > len2 * 1.0 ) return FAILED_FUNC;
	if( alninfo ){
		alninfo[0] = begin1;
		alninfo[1] = end1;
		alninfo[2] = begin2;
		alninfo[3] = end2;
		alninfo[4] = masked;
		if( true ){
			alninfo[0] = gbegin1;
			alninfo[1] = gend1;
			alninfo[2] = gbegin2;
			alninfo[3] = gend2;
		}
	}
#ifdef PRINT
	printf( "%6i %6i:  %4i %4i %4i %4i\n", alnln, iden_no, begin1, end1, begin2, end2 );
	printf( "%6i %6i:  %4i %4i\n", posmin, posmax, posmin - posmax, count - count2 );
	printf( "smin = %9i, smax = %9i\n", smin, smax );
	printf( "dlen = %5i, dcount = %5i, dist = %.3f\n", dlen, dcount, dcount/(float)dlen );
#endif
#ifdef MAKEALIGN
	float identity = iden_no / (float)( options.global_identity ? (len1 - masked) : alnln);
	if( identity < options.cluster_thd ) return OK_FUNC;
	while(i--){
		AA[NN] = letters[ iseq1[i-1] ];
		BB[NN++] = '-';
	}
	while(j--){
		AA[NN] = '-';
		BB[NN++] = letters[ iseq2[j-1] ];
	}
	AA[NN] = '\0';
	BB[NN] = '\0';
	for(i=0; i<NN/2; i++){
		char aa = AA[i], bb = BB[i];
		AA[i] = AA[NN-i-1];
		BB[i] = BB[NN-i-1];
		AA[NN-i-1] = aa;
		BB[NN-i-1] = bb;
	}
	static int fcount = 0; 
	fcount += 1;
	FILE *fout = fopen( "alignments.txt", "a" );
	if( fout == NULL ){
		if( fcount <= 1 ) printf( "alignment files open failed\n" );
		return OK_FUNC;
	}
	fprintf( fout, "\n\n######################################################\n" );
	fprintf( fout, "# length X = %i\n", len2 );
	fprintf( fout, "# length Y = %i\n", len1 );
	fprintf( fout, "# best align X: %i-%i\n", begin2+1, end2+1 );
	fprintf( fout, "# best align Y: %i-%i\n", begin1+1, end1+1 );
	if( alninfo ){
		fprintf( fout, "# align X: %i-%i\n", alninfo[2]+1, alninfo[3]+1 );
		fprintf( fout, "# align Y: %i-%i\n", alninfo[0]+1, alninfo[1]+1 );
	}
	fprintf( fout, "# alignment length: %i\n", alnln );
	fprintf( fout, "# identity count: %i\n", iden_no );
	fprintf( fout, "# identity: %g\n", identity );
	fprintf( fout, "# distance: %g\n", dist );
	if( options.is454 ) fprintf( fout, "# max indel: %i\n", max_indels );
#if 0
	fprintf( fout, "%i %s\n", seq1->index, AA );
	fprintf( fout, "%i %s\n", seq2->index, BB );
#else
	bool printaa = true;
	IB = IA = 0;
	fprintf( fout, "\n\nX " );
	while( IA < NN ){
		if( printaa ){
			fprintf( fout, "%c", BB[IB] );
			IB += 1;
			if( IB % 75 ==0 or IB == NN ) printaa = false, fprintf( fout, "\nY " );
		}else{
			fprintf( fout, "%c", AA[IA] );
			IA += 1;
			if( IA % 75 ==0 ) printaa = true, fprintf( fout, "\n\nX " );
		}
	}
#endif
	fclose( fout );
#endif

	return OK_FUNC;
} // END int local_band_align
int main(int argc, char* argv[]){
    gzFile fp1;
	kseq_t *ks1;
    int max_seq_len = 0;
    CLI::App app{"build word table v.0.0.1, build word table from fasta protein files"};
    string filename = "";
    auto option_input = app.add_option("-i, --input", filename, "input file name, fasta or gziped fasta formats");
    option_input->required();
    CLI11_PARSE(app, argc, argv);
    fp1 = gzopen(filename.c_str(),"r");
    ks1 = kseq_init(fp1);
    vector<Sequence> seqs;
    //   cerr<<"11111"<<endl;
    while(1)
	{
		int length = kseq_read(ks1);
		if(length < 0) break;

		max_seq_len = max<int>(max_seq_len, ks1->seq.l);
	
		Sequence seq;
		seq.name = ks1->name.s;
		// seq.comment = ks1->comment.s;
		seq.seq = ks1->seq.s;
        seq.len = length;
		if(seq.seq.size()<50)continue;
        //  cerr<<"11111"<<endl;
	    char * seq_data = seq.seq.data();
		for(int i = 0; i < seq.seq.size(); i++)
		{
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
    int band_width = 20, required_aa1 = 33723;
    
    WorkingBuffer buffer(max_seq_len);  // 创建工作区

    char* seq1 = seqs[1].seq.data();
    char* seq2 = seqs[0].seq.data();
    int len1 = seqs[1].len;
    int len2 = seqs[0].len;
    double t0 = get_time();
	double t1 = get_time();
    // 计算二元对的倒排表
    // ComputeAAP(seq1, len1, buffer);
    
    // cerr << "ComputeAAP time : " << t1 - t0 << " seconds" << endl;
    // 计算最佳带宽和匹配结果
    int best_sum, band_left, band_center, band_right,best_score,tiden_no,alnln;
    int talign_info[5];
    float tiden_pc, distance=0;
    // diag_test_aapn( seq2, len1, len2, buffer, best_sum, band_width, band_left, band_center, band_right, required_aa1);

	diag_no_table(seq1,len1,seq2,len2,buffer, best_sum, band_width, band_left, band_center, band_right, required_aa1);
    int required_aa2 = 31740;
    if ( best_sum < required_aa2 ) exit(0);
    int rc = FAILED_FUNC;
    double t2 = get_time();
    cerr << "diag time : " << t2 - t1 << " seconds" << endl;
    cout << "Best sum: " << best_sum << endl;
    cout << "Band left: " << band_left << ", Band center: " << band_center << ", Band right: " << band_right << endl;
    rc=local_band_align(seq1, seq2, len1, len2, mat,
					best_score, tiden_no, alnln, distance, talign_info,
					band_left, band_center, band_right, buffer);
    double t3 = get_time();
    cerr << "local_band_align time : " << t3 - t2 << " seconds" << endl;
    cout<<"tiden_no  "<<tiden_no<<endl;
    
    return 0;
}
