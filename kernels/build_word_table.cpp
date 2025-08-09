#include "kseq.h"
#include "CLI11.hpp"
#include <zlib.h>

#include <iostream>
#include <string>
#include <cstring>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <sys/time.h>

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

	double t1 = get_time();
	/// encode and build word table
    for (int seq_id = 0; seq_id < (int)seqs.size(); seq_id++) 
	{
		//encode 
        EncodeWords(seqs[seq_id], word_encodes, word_encodes_no, kmer_size);

		int kmer_no = seqs[seq_id].seq.size() - kmer_size + 1;
        for (int j = 0; j < kmer_no; j++) {
            int bucket = word_encodes[j];
            int count  = word_encodes_no[j];
            if (count > 0) {
				//assert(bucket < table_size);
                word_table[bucket].emplace_back(seq_id, count);

            }
        }
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

	gzclose(fp1);	

	return 0;
}
