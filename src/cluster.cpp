#include "cluster.h"
//SequenceDB seq_db;
Options options;
void initOptions(){
	options = Options();
    options.SetOptionsInputVector();
    options.Validate();
}

void setOptionsClusterThd(float cluster_thd){
	options.setClusterThd(cluster_thd);
	std::cerr << "Cluster-Threshold of cdhit is: " << options.cluster_thd << std::endl;
    options.Validate();
}

void cluster::cdhit_cluster(std::vector<Sequence_new>& seq, std::vector<int>& parent)
{
	SequenceDB seq_db;
    seq_db = SequenceDB();

    string db_out;

    if (seq.size() == 0) {
        std::cout<<"Empty input sequence"<<std::endl;
        return;
    }


    db_out = options.output;
    InitNAA(MAX_UAA);
    options.NAAN = NAAN_array[options.NAA];
    seq_db.NAAN = NAAN_array[options.NAA];

    seq_db.Readvector(seq, options);

    seq_db.SortDivide(options);

    seq_db.DoClustering(options);

	seq_db.updateParent(parent);
}

void cluster::cdhit_cluster(std::vector<Sequence_new>& seq, std::vector<int>& parent, int neededThread)
{
	SequenceDB seq_db;
    seq_db = SequenceDB();

    string db_out;

    if (seq.size() == 0) {
        std::cout<<"Empty input sequence"<<std::endl;
        return;
    }


    db_out = options.output;
    InitNAA(MAX_UAA);
    options.NAAN = NAAN_array[options.NAA];
    seq_db.NAAN = NAAN_array[options.NAA];

    seq_db.Readvector(seq, options);

    seq_db.SortDivide(options);

	if(neededThread == 1){
    	seq_db.DoClustering(options);
	}else{
    	seq_db.DoClustering(neededThread, options);
	}
	seq_db.updateParent(parent);
}

