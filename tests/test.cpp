#include <iostream>
#include <vector>
#include "libcdhit/cdhit.h"

int main() {

    std::vector<Sequence_new> seqs = {
        Sequence_new(0, "ATCGATCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        Sequence_new(1, "ATCGATCGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
        Sequence_new(2, "GCGTGCGT"),
        Sequence_new(3, "ATCGATCG")
    };

    std::vector<int> parent;
	parent.resize(seqs.size());

    cluster clst;
    clst.cdhit_cluster(seqs, parent);

    std::cout << "Cluster Results (parent array):" << std::endl;
    for (size_t i = 0; i < parent.size(); ++i) {
        std::cout << "Sequence " << i << " --> Cluster Root: " << parent[i] << std::endl;
    }

    return 0;
}

