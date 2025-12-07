#include <algorithm>
#include <chrono>
#include <cstdint>
#include <functional>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "ips4o.hpp"
#include "CLI11.hpp"

using Clock = std::chrono::high_resolution_clock;

// ----------------------
// 自定义结构体类型
// ----------------------
struct Record {
    uint64_t id;
    uint64_t value;
};

// ----------------------
// 测试函数
// ----------------------
template<class SortFn>
long long operate_once(const std::vector<Record>& base, SortFn fn, const char* name) {
    std::vector<Record> v = base;
    auto t0 = Clock::now();
    fn(v.begin(), v.end());
    auto t1 = Clock::now();

    // 检查是否按 value 排序正确
    if (!std::is_sorted(v.begin(), v.end(), [](const Record& a, const Record& b) {
        return a.value < b.value;
    })) {
        std::cerr << "[Error] " << name << " result is not sorted!\n";
        std::exit(1);
    }

    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << "[Info] " << name << " sorted time is [ " << ms << " ] ms\n";
    return ms;
}

int main(int argc, char** argv) {
    size_t N = 280000000;
    int R = 3;
    uint32_t seed = 123;
    bool run_std = false;
    bool run_seq = false;
    bool run_par = false;
    bool descending = false;
    size_t threads = 1;

    CLI::App app{"IPS4o Struct Sort Benchmark"};
    app.add_option("-n,--size", N, "元素个数")->check(CLI::PositiveNumber);
    app.add_option("-r,--repeats", R, "重复次数")->check(CLI::PositiveNumber);
    app.add_option("-s,--seed", seed, "随机种子");
    app.add_flag("--std", run_std, "运行 std::sort");
    app.add_flag("--seq", run_seq, "运行 ips4o::sort（顺序）");
#ifdef USE_IPS4O_PARALLEL
    app.add_flag("--par", run_par, "运行 ips4o::parallel::sort（并行）");
    app.add_option("-t,--threads", threads, "线程数");
#endif
    app.add_flag("--desc,--descending", descending, "降序排序（默认升序）");
    CLI11_PARSE(app, argc, argv);

#ifdef _OPENMP
    omp_set_dynamic(0);                 // 固定线程数
    omp_set_max_active_levels(1);       // 禁止嵌套并行
    omp_set_num_threads((int)threads);
#endif

    if (!(run_std || run_seq
#ifdef USE_IPS4O_PARALLEL
        || run_par
#endif
    )) {
        run_std = true;
        run_seq = true;
    }

    // 随机生成数据：id 连续递增，value 随机
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<uint64_t> dist(0, std::numeric_limits<uint64_t>::max());
    std::vector<Record> base(N);
    for (size_t i = 0; i < N; ++i) {
        base[i].id = i;
        base[i].value = dist(rng);
    }

    // 比较器：基于 value 排序
    auto cmp_asc  = [](const Record& a, const Record& b) { return a.value < b.value; };
    auto cmp_desc = [](const Record& a, const Record& b) { return a.value > b.value; };
    auto cmp = descending ? cmp_desc : cmp_asc;

    auto run_case = [&](const char* name, auto sorter) {
        long long best = std::numeric_limits<long long>::max();
        for (int i = 0; i < R; i++) {
            best = std::min(best, operate_once(base, sorter, name));
        }
        std::cout << "[Info] " << name << " best sort time is " << best << " ms .\n";
        std::cout << "[Info] Avg time is " << (double)best / N * 1e6 << " ns/elem .\n\n";
    };

    if (run_std) {
        run_case("std::sort", [&](auto b, auto e) {
            if (descending) std::sort(b, e, cmp_desc);
            else            std::sort(b, e, cmp_asc);
        });
    }
    if (run_seq) {
        run_case("ips4o::sort", [&](auto b, auto e) {
            if (descending)
                ips4o::sort(b, e, cmp_desc);
            else
                ips4o::sort(b, e, cmp_asc);
        });
    }
#ifdef USE_IPS4O_PARALLEL
    if (run_par) {
        run_case("ips4o::parallel::sort", [&](auto b, auto e) {
            if (descending)
                ips4o::parallel::sort(b, e, cmp_desc);
            else
                ips4o::parallel::sort(b, e, cmp_asc);
        });
    }
#endif

    return 0;
}
