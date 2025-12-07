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

#include "CLI11.hpp"
#include "ips2ra.hpp"
#include "ips4o.hpp"

using namespace std;
using Clock = std::chrono::high_resolution_clock;

struct Data {
    std::vector<uint64_t> value;
    int id;
    Data(int id, std::vector<uint64_t> value) : id(id), value(value) {}
    Data() {}
};

vector<Data> hash_vec;

void ips4o_sort_single_thread(vector<Data>& hash_vec, int r = 1) {
    auto cmp = [r](const Data& a, const Data& b) {
        for (int i = 0; i < r; i++) {
            if (a.value[i] < b.value[i])
                return true;
            else if (a.value[i] > b.value[i])
                return false;
        }
        return false;
    };
    ips4o::sort(hash_vec.begin(), hash_vec.end(), cmp);
}

void ips4o_sort_multi_thread(vector<Data>& hash_vec, int r = 1) {
    auto cmp = [r](const Data& a, const Data& b) {
        for (int i = 0; i < r; i++) {
            if (a.value[i] < b.value[i])
                return true;
            else if (a.value[i] > b.value[i])
                return false;
        }
        return false;
    };
    ips4o::parallel::sort(hash_vec.begin(), hash_vec.end(), cmp);
}

void ips2ra_sort_single_thread(vector<Data>& hash_vec, const int r = 1) {
    ips2ra::sort(hash_vec.begin(), hash_vec.end(), [](const Data& r) { return r.value[0]; });
}

void ips2ra_sort_multi_thread(vector<Data>& hash_vec, const int r = 1) {
    ips2ra::parallel::sort(hash_vec.begin(), hash_vec.end(), [](const Data& r) { return r.value[0]; });
}

int main(int argc, char** argv) {
    size_t N = 1000000;
    int R = 3;
    uint32_t seed = 123;
    bool run_std = false;
    bool run_seq = false;
    bool run_par = false;
    size_t threads = 1;
    int vector_size = 1;

    CLI::App app{"IPS4o Test Analyze"};
    app.add_option("-n,--size", N, "元素个数")->check(CLI::PositiveNumber);
    app.add_option("-r,--repeats", R, "重复次数")->check(CLI::PositiveNumber);
    app.add_option("-s,--seed", seed, "随机种子");
    app.add_option("-v,--vec_size", vector_size, "抽样数");
    app.add_flag("--std", run_std, "运行 std::sort");
    app.add_flag("--seq", run_seq, "运行 ips4o::sort and ips2ra::sort");
    app.add_flag("--par", run_par, "运行 ips4o::parallel::sort and ips2ra::parallel::sort");
    app.add_option("-t,--threads", threads, "线程数");
    CLI11_PARSE(app, argc, argv);

#ifdef _OPENMP
    omp_set_dynamic(0);           // 固定线程数
    omp_set_max_active_levels(1); // 禁止嵌套并行
    omp_set_num_threads((int) threads);
#endif

    if (!(run_std || run_seq
#ifdef USE_IPS4O_PARALLEL
          || run_par
#endif
          )) {
        run_std = true;
        run_seq = true;
    }

    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<uint64_t> dist(0, std::numeric_limits<uint64_t>::max());

    hash_vec.clear();
    hash_vec.reserve(N);

    for (size_t i = 0; i < N; ++i) {
        std::vector<uint64_t> value(vector_size);
        for (size_t j = 0; j < vector_size; ++j) {
            value[j] = dist(rng);
        }
        hash_vec.emplace_back(i, std::move(value));
    }

    auto run_case = [&](const char* name, auto sorter) {
        uint64_t total_time = 0;
        printf("--------------------------------\n");
        for (int i = 0; i < R; i++) {
            std::vector<Data> v = hash_vec; // 每次复制原始数据

            auto t0 = Clock::now();
            sorter(v, vector_size);
            auto t1 = Clock::now();

            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
            total_time += ms;
            std::cout << "[Info] " << name << " time is " << ms << " ms\n";
        }
        std::cout << "[Info] " << name << " average time is " << total_time / R << " ms\n";
        std::cout << "[Info] " << name << " avg per element: " << (double) total_time / R / N * 1e6 << " ns/elem\n\n";
    };

    // 运行各个排序算法
    if (run_std) {
        run_case("std::sort",
                 [](vector<Data>& v, int vs) { std::sort(v.begin(), v.end(), [](const Data& a, const Data& b) { return a.value < b.value; }); });
    }

    if (run_seq) {
        run_case("ips4o::sort", ips4o_sort_single_thread);
        run_case("ips2ra::sort", ips2ra_sort_single_thread);
    }

#ifdef USE_IPS4O_PARALLEL
    if (run_par) {
        run_case("ips4o::parallel::sort", ips4o_sort_multi_thread);
        run_case("ips2ra::parallel::sort", ips2ra_sort_multi_thread);
    }
#endif

    return 0;
}