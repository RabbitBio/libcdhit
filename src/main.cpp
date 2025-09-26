#include <algorithm>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "ips4o.hpp"
#include "CLI11.hpp"
using Clock = std::chrono::high_resolution_clock;

template<class SortFn>
long long operate_once(const std::vector<int>& base, SortFn fn, const char* name){
	std::vector<int> v = base;
	auto t0 = Clock::now();
	fn(v.begin(),v.end());
	auto t1 = Clock::now();
	if(!std::is_sorted(v.begin(),v.end())){
		std::cerr<<"[Error] "<<name<<" result is not sorted!\n";
		std::exit(1);	
	}
	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
	std::cout<<"[Info] "<<name<<" sorted time is [ "<< ms << " ] ms\n";
	return ms;
}

int main(int argc,char** argv){
	size_t N = 10000000;
	int R = 3;
	uint32_t seed = 123;
	bool run_std = false;
	bool run_seq = false;
	bool run_par = false;
	bool descending = false;

	CLI::App app{"IPS4o Test Analyze"};
	app.add_option("-n,--size", N, "元素个数")->check(CLI::PositiveNumber);
    app.add_option("-r,--repeats", R, "重复次数")->check(CLI::PositiveNumber);
    app.add_option("-s,--seed", seed, "随机种子");
    app.add_flag("--std", run_std, "运行 std::sort");
    app.add_flag("--seq", run_seq, "运行 ips4o::sort（顺序）");
#ifdef USE_IPS4O_PARALLEL
    app.add_flag("--par", run_par, "运行 ips4o::parallel::sort（并行）");
#endif
    app.add_flag("--desc,--descending", descending, "降序排序（默认升序）");
	CLI11_PARSE(app,argc,argv);
	if(!(run_std || run_seq
#ifdef USE_IPS4O_PARALLEL
         || run_par
#endif
    )) {
        run_std = true;
        run_seq = true;
    }

	std::mt19937 rng(seed);
	std::uniform_int_distribution<int> dist(std::numeric_limits<int>::min(),
                                            std::numeric_limits<int>::max());
    std::vector<int> base(N);
    for (auto &x : base) x = dist(rng);
	auto cmp_asc = [](const int& a, const int& b){ return a < b; };
    auto cmp_desc = [](const int& a, const int& b){ return a > b; };
    auto cmp = descending ? cmp_desc : cmp_asc;

	auto run_case = [&](const char* name, auto sorter){
		long long best = std::numeric_limits<long long>::max();
		for(int i=0;i<R;i++){
			best = std::min(best,operate_once(base,sorter,name));
		}
		std::cout<<"[Info] " <<name<<" best sort time is "<<best<<" ms .\n";
		std::cout<<"[Info] Avg time is " << (double)best/N*1e6 <<" ns/elem .\n\n";
	};
	if(run_std){
		run_case("std::sort",[&](auto b, auto e){std::sort(b,e,cmp);});
	}
	if(run_seq){
		run_case("ips4o::sort",[&](auto b, auto e){ips4o::sort(b,e,cmp);});
	}
#ifdef USE_IPS4O_PARALLEL
	if(run_par){
		run_case("ips4o::parallel::sort",[&](auto b, auto e){ips4o::parallel::sort(b,e,cmp);});
	}
#endif
	return 0;
}

