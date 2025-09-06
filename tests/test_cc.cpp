#include "libcdhit/cdhit.h"
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
#include <unordered_set>

using namespace std;

KSEQ_INIT(gzFile, gzread)

double get_time()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
}
int main(int argc, char* argv[])
{

	gzFile fp1;
	kseq_t *ks1;

	CLI::App app{"test connect component clustering"};
	string filename = "";
	string res_file = "";
	int64_t max_num_seqs = 0;
	int kmer_size = 5;
	double tau = 0.05;   // Jaccard 阈值：按需设置

	auto option_input = app.add_option("-i, --input", filename, "input file name, fasta or gziped fasta formats");
	auto option_max_n = app.add_option("-n, --max_num_seqs", max_num_seqs, "max number of seqs for building word table");
	auto option_k = app.add_option("-k, --kmer_size", kmer_size, "set kmer size, default 5.");
	auto option_tau = app.add_option("-t, --jaccard_thres", tau, "set weighted jaccard threshold, default 0.36.");
	auto option_output = app.add_option("-o, --output", res_file, "output files");

	option_input->required();

	CLI11_PARSE(app, argc, argv);

	fp1 = gzopen(filename.c_str(),"r");
	cerr << "Threshold: " << tau << endl;
	cerr << "Kmer_size: " << kmer_size << endl;

	if(NULL == fp1){
		cerr << "Fail to open file: " << filename << endl;
		return 0;
	}

	ks1 = kseq_init(fp1);


	vector<Sequence_new> seqs;
	int64_t pos = 0;
	int64_t number_seqs = 0;
	int max_seq_len = 0;

	while(1)
	{
		int length = kseq_read(ks1);
		if(length < 0) break;

		int l = ks1->seq.l;
		max_seq_len = std::max(max_seq_len, l);
		seqs.push_back(Sequence_new(number_seqs,strdup(ks1->seq.s)));

		number_seqs++;

		if(max_num_seqs > 0 && number_seqs >= max_num_seqs) break;
	}

	cerr << "number of sequences: " << seqs.size() << endl;
	cerr << "max seq length:" << max_seq_len << endl;
	if(seqs.size() < max_num_seqs) 
		cerr << "No more seqs than required! Total: "<< seqs.size() << "\tRequired: " << max_num_seqs <<endl;

/// init cc clustering buffer
	int nthreads = 2;

	cerr << "Threads: " << nthreads << endl;
	vector<ClusterBuffer> buffs(nthreads, ClusterBuffer(seqs.size(), kmer_size, max_seq_len));
	cerr << "init buffer " << endl;
/// call clustering api
	std::vector<int> parent;
	parent.resize(seqs.size());

	double t1 = get_time();
	int chunk = 1000;
	#pragma omp parallel for num_threads(nthreads) shared(parent, buffs)
	for(int i = 0; i < seqs.size(); i+=chunk)
	{
		int tid = omp_get_thread_num();
		cerr << "tid: " << tid << endl;
		buffs[tid].clear();
		cerr << "after clear()" << endl;
		vector<Sequence_new> seqs_local(seqs.begin() + i, seqs.begin()+i+chunk);
		cluster_sequences_st(buffs[tid], seqs_local, parent, kmer_size, tau);
	}
	double t2 = get_time();
	// 打印结果
	//std::cout << "Parent array:" << std::endl;
	//for (size_t i = 0; i < parent.size(); ++i) {
	//	std::cout << "seq " << i << ": parent = " << parent[i] << std::endl;
	//}

    unordered_set<int> roots;
	roots.reserve(parent.size() * 2);
	for (const auto& pr : parent) roots.insert(pr);
	cerr << "Number of clusters: " << roots.size() << endl;
	cerr << "Clustering time: " << t2 - t1 << " s" << endl;
	cerr << "Avg seq per second: " << (double)number_seqs / (t2 - t1) << endl;

	return 0;
}
