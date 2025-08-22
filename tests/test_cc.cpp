#include "libcdhit/cdhit.h"
#include <vector>
#include <iostream>

int main() {
	// 示例输入
	std::vector<Sequence_new> test_seqs = {
		Sequence_new(0, "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLLVFAVNSAKSFEDIGTYREQIKRVKDAEEVPMVLVGNKCDLASWNVNNEQHAVKLARELLLNHYQQ"),
		Sequence_new(1, "MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLLVFAVNSAKSFEDIGTYREQIKRVKDAEEVPMVLVGNKCDLASWNVNNEQHAVKLARELLLNHYQQ"),  // 相似
		Sequence_new(2, "ABCDEF")  // 不同
	};

	std::vector<int> parent;
	parent.resize(test_seqs.size());

	cluster_sequences(test_seqs, parent, 5, 0.36,1);

	// 打印结果
	std::cout << "Parent array:" << std::endl;
	for (size_t i = 0; i < parent.size(); ++i) {
		std::cout << "seq " << i << ": parent = " << parent[i] << std::endl;
	}

	return 0;
}
