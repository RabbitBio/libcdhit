cmake_minimum_required(VERSION 3.16)
project(ips_cli_bench LANGUAGES CXX)

# --- 基本编译设置 ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release CACHE STRING "")

# --- 选项（可按需 -D 开关） ---
option(OPTIMIZE_NATIVE "为本机CPU优化（-O3 -march=native）" ON)
option(ENABLE_PARALLEL "启用并行路径（仅对 ips4o 定义 USE_IPS4O_PARALLEL）" ON)
option(ENABLE_OPENMP   "若可用则链接 OpenMP" OFF)
option(LINK_ATOMIC     "如遇 __atomic_compare_exchange_16 链接错误则开" OFF)

# --- 让 ips4o / ips2ra 当子工程编 ---
# 注意：在 add_subdirectory 之前设置子项目变量更生效
# 如果你不想本机优化，可把下面两个设成 OFF
set(IPS4O_OPTIMIZE_FOR_NATIVE ${OPTIMIZE_NATIVE} CACHE BOOL "" FORCE)
set(IPS4O_USE_MCX16           ${OPTIMIZE_NATIVE} CACHE BOOL "" FORCE)

add_subdirectory(thirdparty/ips4o  ips4o-build)
add_subdirectory(thirdparty/ips2ra ips2ra-build)

# --- 你的可执行程序 ---
# src 中两个测试脚本
add_executable(bench_ips4o src/ips4o.cpp)
add_executable(bench_ips2ra src/ips2ra.cpp)
add_executable(bench_main src/main.cpp)

# 公共包含目录（例如 include/ 里有 CLI.hpp）
target_include_directories(bench_ips4o PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_include_directories(bench_ips2ra PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_include_directories(bench_main PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# 分别链接各自的库（由上面的 add_subdirectory 提供）
target_link_libraries(bench_ips4o PRIVATE ips4o)
target_link_libraries(bench_ips2ra PRIVATE ips2ra)
target_link_libraries(bench_main PRIVATE ips4o ips2ra)

# --- 并行与依赖（可选） ---
if(ENABLE_PARALLEL)
  # 仅对 ips4o 测试启用该宏（若代码需要）
  target_compile_definitions(bench_ips4o PRIVATE USE_PARALLEL)
  target_compile_definitions(bench_ips2ra PRIVATE USE_PARALLEL)
  target_compile_definitions(bench_main PRIVATE USE_PARALLEL)

  # 如系统已安装 TBB，尝试链接（两者都给，找不到也没关系）
  find_package(TBB QUIET)
  if(TBB_FOUND)
    target_link_libraries(bench_ips4o PRIVATE TBB::tbb)
    target_link_libraries(bench_ips2ra PRIVATE TBB::tbb)
    target_link_libraries(bench_main PRIVATE TBB::tbb)
  endif()
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(bench_ips4o PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(bench_ips2ra PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(bench_main PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()

# --- 编译优化（不想改全局 CFLAGS 的就只对目标加） ---
if(OPTIMIZE_NATIVE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(bench_ips4o PRIVATE -O3 -DNDEBUG -march=native)
    target_compile_options(bench_ips2ra PRIVATE -O3 -DNDEBUG -march=native)
    target_compile_options(bench_main PRIVATE -O3 -DNDEBUG -march=native)
  elseif(MSVC)
    target_compile_options(bench_ips4o PRIVATE /O2)
    target_compile_options(bench_ips2ra PRIVATE /O2)
    target_compile_options(bench_main PRIVATE /O2)
  endif()
endif()

# --- 部分平台需要手动链接 libatomic（只在报错时开） ---
if(LINK_ATOMIC)
  if(UNIX AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_link_libraries(bench_ips4o PRIVATE atomic)
    target_link_libraries(bench_ips2ra PRIVATE atomic)
    target_link_libraries(bench_main PRIVATE atomic)
  endif()
endif()
