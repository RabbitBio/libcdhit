cmake_minimum_required(VERSION 3.16)
project(ips_cli_bench LANGUAGES CXX)

# --- 基本编译设置 ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release CACHE STRING "")

# --- 选项（可按需 -D 开关） ---
option(OPTIMIZE_NATIVE "为本机CPU优化（-O3 -march=native）" ON)
option(ENABLE_PARALLEL "启用并行路径（仅对 ips4o 定义 USE_IPS4O_PARALLEL）" ON)
option(LINK_ATOMIC     "如遇 __atomic_compare_exchange_16 链接错误则开" OFF)

# --- 让 ips4o / ips2ra 当子工程编 ---
# 注意：在 add_subdirectory 之前设置子项目变量更生效
# 如果你不想本机优化，可把下面两个设成 OFF
set(IPS4O_OPTIMIZE_FOR_NATIVE ${OPTIMIZE_NATIVE} CACHE BOOL "" FORCE)
set(IPS4O_USE_MCX16           ${OPTIMIZE_NATIVE} CACHE BOOL "" FORCE)

add_subdirectory(thirdparty/ips4o  ips4o-build)
add_subdirectory(thirdparty/ips2ra ips2ra-build)


add_executable(bench_main src/main.cpp)

# 公共包含目录（例如 include/ 里有 CLI.hpp）
target_include_directories(bench_main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# 分别链接各自的库（由上面的 add_subdirectory 提供）
target_link_libraries(bench_main PRIVATE ips4o ips2ra)

# --- 并行与依赖（可选） ---
if(ENABLE_PARALLEL)
  # 仅对 ips4o 测试启用该宏（若代码需要）
  target_compile_definitions(bench_main PRIVATE USE_IPS4O_PARALLEL)
endif()

# --- 编译优化（不想改全局 CFLAGS 的就只对目标加） ---
if(OPTIMIZE_NATIVE)
  target_compile_options(bench_main PRIVATE -O3 -DNDEBUG -march=native)
endif()

# --- 部分平台需要手动链接 libatomic（只在报错时开） ---
if(LINK_ATOMIC)
  if(UNIX AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_link_libraries(bench_main PRIVATE atomic)
  endif()
endif()
